#!/usr/bin/env python3
"""Package and import the complete de-duplicated ModCRE model collection.

Dry-run is the default and operates on an in-memory copy of SQLite.  Apply
mode creates the serving archive atomically, backs up SQLite, applies migration
003 when needed, updates existing public ModCRE rows in place (preserving
foreign-key links), and inserts only missing biological model identities.

The script never edits PWM/motif provenance, model summaries, AF3 assignments,
or the original Baldo coordinate files.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import os
import re
import shutil
import sqlite3
import tarfile
import tempfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, BinaryIO, Iterable


DEFAULT_DATABASE = Path("data/tf_webdb_local_staging.sqlite")
DEFAULT_MANIFEST = Path("outputs/modcre_complete_model_manifest.tsv")
DEFAULT_MANIFEST_QC = Path("outputs/modcre_complete_model_manifest_qc.json")
DEFAULT_MIGRATION = Path("migrations/003_add_structure_model_qc_columns.sql")
DEFAULT_ARCHIVE = Path("/home/patricia/TF_database_Baldo_data/modcre_complete_models.tar.gz")
DEFAULT_PREVIEW = Path("outputs/modcre_complete_model_import_preview.tsv")
DEFAULT_QC_OUTPUT = Path("outputs/modcre_complete_model_import_qc.json")
DEFAULT_BACKUP_DIR = Path("backups")

EXPECTED_MODELS = 58_976
EXPECTED_TFS = 40_296
EXPECTED_DIMERS = 18_680
EXPECTED_TFS_ACCESSIONS = 4_562
EXPECTED_INTERACTION_YES = 58_427
EXPECTED_INTERACTION_NO = 549
EXPECTED_NO_PROTEIN = 4
EXPECTED_EXISTING_ACTIVE = 11_237
EXPECTED_INSERTS_INITIAL = 47_739
EXPECTED_ACTIVE_AF3_ASSIGNMENTS = 114

QC_COLUMNS = {
    "model_kind",
    "protein_chain_count",
    "dna_chain_count",
    "dna_interaction",
    "dna_contact_count",
    "dna_qc_cutoff",
    "dna_qc_status",
    "dna_qc_checked_at",
}
BASE_COLUMNS = {
    "id",
    "source",
    "model_id",
    "tf_id",
    "member_path",
    "archive_path",
    "file_type",
    "status",
    "template_pdb",
    "residue_start",
    "residue_end",
}
PROTECTED_TABLES = (
    "tf",
    "motif_ref",
    "motif_file",
    "motif_structure",
    "tf_primary_annotation",
    "model_summary",
    "structure_model_assignment",
    "structure_confidence_artifact",
    "tf_structure_status",
)
PREVIEW_COLUMNS = (
    "tf_id",
    "model_id",
    "model_kind",
    "coordinate_source",
    "dna_chain_count",
    "dna_interaction",
    "dna_contact_count",
    "existing_structure_file_id",
    "planned_action",
    "serving_member_path",
)

OLD_MODEL_RE = re.compile(
    r"^(TFS|DIMER)_(?:tr|sp)_([^_]+)_[^:]+:(\d+):(\d+)_(.+)$"
)
NEW_MODEL_RE = re.compile(r"^(TFS|DIMER)_([^:]+):(\d+):(\d+)_(.+)$")


class ValidationError(RuntimeError):
    pass


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database", type=Path, default=DEFAULT_DATABASE)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--manifest-qc", type=Path, default=DEFAULT_MANIFEST_QC)
    parser.add_argument("--migration", type=Path, default=DEFAULT_MIGRATION)
    parser.add_argument("--archive", type=Path, default=DEFAULT_ARCHIVE)
    parser.add_argument("--preview-output", type=Path, default=DEFAULT_PREVIEW)
    parser.add_argument("--qc-output", type=Path, default=DEFAULT_QC_OUTPUT)
    parser.add_argument("--backup-dir", type=Path, default=DEFAULT_BACKUP_DIR)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--dry-run", action="store_true", help="Validate using an in-memory DB copy")
    mode.add_argument("--apply", action="store_true", help="Create archive and update staging DB")
    parser.add_argument(
        "--replace-preview",
        action="store_true",
        help="Permit replacing preview/QC output files",
    )
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def table_columns(connection: sqlite3.Connection, table: str) -> list[str]:
    return [str(row[1]) for row in connection.execute(f"PRAGMA table_info({table})")]


def table_exists(connection: sqlite3.Connection, table: str) -> bool:
    return connection.execute(
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name=?", (table,)
    ).fetchone() is not None


def encoded_value(value: Any) -> bytes:
    if value is None:
        return b"N"
    payload = value if isinstance(value, bytes) else str(value).encode("utf-8")
    prefix = b"B" if isinstance(value, bytes) else b"T"
    return prefix + str(len(payload)).encode("ascii") + b":" + payload


def table_snapshot(connection: sqlite3.Connection, table: str) -> tuple[int, str]:
    if not table_exists(connection, table):
        return (0, "ABSENT")
    columns = table_columns(connection, table)
    quoted = ", ".join('"' + column.replace('"', '""') + '"' for column in columns)
    digest = hashlib.sha256()
    count = 0
    for row in connection.execute(f"SELECT {quoted} FROM {table} ORDER BY {quoted}"):
        count += 1
        for value in row:
            digest.update(encoded_value(value))
            digest.update(b"|")
        digest.update(b"\n")
    return count, digest.hexdigest()


def protected_snapshots(connection: sqlite3.Connection) -> dict[str, tuple[int, str]]:
    return {table: table_snapshot(connection, table) for table in PROTECTED_TABLES}


def normalize_model_id(model_id: str) -> str | None:
    match = OLD_MODEL_RE.fullmatch(model_id) or NEW_MODEL_RE.fullmatch(model_id)
    if match is None:
        return None
    kind, accession, start, end, tail = match.groups()
    return f"{kind}_{accession}:{start}:{end}_{tail}"


def read_manifest(path: Path, qc_path: Path) -> tuple[list[dict[str, Any]], str]:
    if not path.is_file():
        raise ValidationError(f"missing manifest: {path}")
    if not qc_path.is_file():
        raise ValidationError(f"missing manifest QC: {qc_path}")
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {
        "tf_id", "model_id", "model_kind", "residue_start", "residue_end",
        "template_pdb", "coordinate_source", "source_path",
        "source_archive_path", "source_member_path", "serving_member_path",
        "file_size", "sha256", "protein_chain_count", "dna_chain_count",
        "dna_interaction", "dna_contact_count", "dna_qc_cutoff", "dna_qc_status",
    }
    if not rows or not required.issubset(rows[0]):
        raise ValidationError("manifest is empty or lacks required columns")
    converted: list[dict[str, Any]] = []
    for row in rows:
        converted.append(
            {
                **row,
                "residue_start": int(row["residue_start"]),
                "residue_end": int(row["residue_end"]),
                "file_size": int(row["file_size"]),
                "protein_chain_count": int(row["protein_chain_count"]),
                "dna_chain_count": int(row["dna_chain_count"]),
                "dna_interaction": int(row["dna_interaction"]),
                "dna_contact_count": int(row["dna_contact_count"] or 0),
                "dna_qc_cutoff": float(row["dna_qc_cutoff"]),
            }
        )
    qc = json.loads(qc_path.read_text())
    if qc.get("all_invariants_passed") is not True:
        raise ValidationError("manifest QC does not report passed invariants")
    checked_at = str(qc.get("generated_at_utc") or "").strip()
    if not checked_at:
        raise ValidationError("manifest QC lacks generated_at_utc")

    if len(converted) != EXPECTED_MODELS:
        raise ValidationError(f"expected {EXPECTED_MODELS} manifest rows, found {len(converted)}")
    if len({row["model_id"] for row in converted}) != EXPECTED_MODELS:
        raise ValidationError("manifest model_id values are not unique")
    if len({row["serving_member_path"] for row in converted}) != EXPECTED_MODELS:
        raise ValidationError("manifest serving_member_path values are not unique")
    if len({row["tf_id"] for row in converted}) != EXPECTED_TFS_ACCESSIONS:
        raise ValidationError("unexpected number of represented TF accessions")
    kinds = Counter(row["model_kind"] for row in converted)
    if kinds != Counter({"TFS": EXPECTED_TFS, "DIMER": EXPECTED_DIMERS}):
        raise ValidationError(f"unexpected model-kind counts: {dict(kinds)}")
    interactions = Counter(row["dna_interaction"] for row in converted)
    if interactions != Counter({1: EXPECTED_INTERACTION_YES, 0: EXPECTED_INTERACTION_NO}):
        raise ValidationError(f"unexpected interaction counts: {dict(interactions)}")
    if sum(row["dna_qc_status"] == "no_protein_atoms" for row in converted) != EXPECTED_NO_PROTEIN:
        raise ValidationError("unexpected DNA-only model count")
    return converted, checked_at


def connect_readonly(path: Path) -> sqlite3.Connection:
    connection = sqlite3.connect(f"file:{path}?mode=ro&immutable=1", uri=True)
    connection.row_factory = sqlite3.Row
    connection.execute("PRAGMA foreign_keys = ON")
    return connection


def in_memory_copy(path: Path) -> sqlite3.Connection:
    source = connect_readonly(path)
    destination = sqlite3.connect(":memory:")
    destination.row_factory = sqlite3.Row
    source.backup(destination)
    source.close()
    destination.execute("PRAGMA foreign_keys = ON")
    return destination


def apply_migration_if_needed(connection: sqlite3.Connection, migration: Path) -> bool:
    columns = set(table_columns(connection, "structure_file"))
    present = columns & QC_COLUMNS
    if present == QC_COLUMNS:
        return False
    if present:
        raise ValidationError(
            "structure_file has a partial migration-003 state: " + ", ".join(sorted(present))
        )
    if not migration.is_file():
        raise ValidationError(f"missing migration: {migration}")
    connection.executescript(migration.read_text())
    after = set(table_columns(connection, "structure_file"))
    if not (BASE_COLUMNS | QC_COLUMNS).issubset(after):
        raise ValidationError("migration 003 did not create all required structure_file columns")
    return True


def desired_values(row: dict[str, Any], archive: Path, checked_at: str) -> tuple[Any, ...]:
    return (
        "modcre",
        row["model_id"],
        row["tf_id"],
        row["serving_member_path"],
        str(archive.resolve()),
        "pdb",
        "active",
        row["template_pdb"],
        row["residue_start"],
        row["residue_end"],
        row["model_kind"],
        row["protein_chain_count"],
        row["dna_chain_count"],
        row["dna_interaction"],
        row["dna_contact_count"],
        row["dna_qc_cutoff"],
        row["dna_qc_status"],
        checked_at,
    )


def existing_values(row: sqlite3.Row) -> tuple[Any, ...]:
    return tuple(
        row[column]
        for column in (
            "source", "model_id", "tf_id", "member_path", "archive_path", "file_type",
            "status", "template_pdb", "residue_start", "residue_end", "model_kind",
            "protein_chain_count", "dna_chain_count", "dna_interaction",
            "dna_contact_count", "dna_qc_cutoff", "dna_qc_status", "dna_qc_checked_at",
        )
    )


def plan_changes(
    connection: sqlite3.Connection,
    rows: list[dict[str, Any]],
    archive: Path,
    checked_at: str,
) -> tuple[list[dict[str, Any]], dict[str, int]]:
    tf_ids = {str(row[0]) for row in connection.execute("SELECT tf_id FROM tf")}
    missing_tf = sorted({str(row["tf_id"]) for row in rows} - tf_ids)
    if missing_tf:
        raise ValidationError("manifest accessions absent from tf: " + ", ".join(missing_tf[:10]))

    existing_rows = connection.execute(
        "SELECT * FROM structure_file WHERE LOWER(source)='modcre' AND file_type='pdb'"
    ).fetchall()
    by_key: dict[str, list[sqlite3.Row]] = defaultdict(list)
    unparsed: list[str] = []
    for existing in existing_rows:
        key = normalize_model_id(str(existing["model_id"]))
        if key is None:
            unparsed.append(str(existing["model_id"]))
        else:
            by_key[key].append(existing)
    if unparsed:
        raise ValidationError("unparsed existing ModCRE model IDs: " + ", ".join(unparsed[:10]))

    target_keys = {str(row["model_id"]) for row in rows}
    active_outside = sorted(
        key
        for key, candidates in by_key.items()
        if key not in target_keys and any(str(item["status"]) == "active" for item in candidates)
    )
    if active_outside:
        raise ValidationError(
            "active ModCRE PDB rows outside complete target: " + ", ".join(active_outside[:10])
        )

    preview: list[dict[str, Any]] = []
    counts = Counter()
    for row in rows:
        candidates = by_key.get(str(row["model_id"]), [])
        active = [item for item in candidates if str(item["status"]) == "active"]
        if len(active) > 1:
            raise ValidationError(f"multiple active rows for model identity: {row['model_id']}")
        keeper = active[0] if active else (min(candidates, key=lambda item: int(item["id"])) if candidates else None)
        desired = desired_values(row, archive, checked_at)
        if keeper is None:
            action = "INSERT"
        elif existing_values(keeper) == desired:
            action = "NO_CHANGE"
        else:
            action = "UPDATE_EXISTING"
        counts[action] += 1
        preview.append(
            {
                **row,
                "existing_structure_file_id": int(keeper["id"]) if keeper is not None else "",
                "planned_action": action,
            }
        )
    if sum(counts.values()) != EXPECTED_MODELS:
        raise ValidationError("change plan does not cover every manifest row")
    return preview, dict(counts)


def execute_changes(
    connection: sqlite3.Connection,
    preview: list[dict[str, Any]],
    archive: Path,
    checked_at: str,
) -> dict[str, int]:
    update_sql = """
        UPDATE structure_file
        SET source=?, model_id=?, tf_id=?, member_path=?, archive_path=?, file_type=?,
            status=?, template_pdb=?, residue_start=?, residue_end=?, model_kind=?,
            protein_chain_count=?, dna_chain_count=?, dna_interaction=?,
            dna_contact_count=?, dna_qc_cutoff=?, dna_qc_status=?, dna_qc_checked_at=?
        WHERE id=?
    """
    insert_sql = """
        INSERT INTO structure_file (
            source, model_id, tf_id, member_path, archive_path, file_type, status,
            template_pdb, residue_start, residue_end, model_kind, protein_chain_count,
            dna_chain_count, dna_interaction, dna_contact_count, dna_qc_cutoff,
            dna_qc_status, dna_qc_checked_at
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """
    counts = Counter()
    with connection:
        for row in preview:
            action = str(row["planned_action"])
            desired = desired_values(row, archive, checked_at)
            if action == "UPDATE_EXISTING":
                connection.execute(update_sql, (*desired, int(row["existing_structure_file_id"])))
            elif action == "INSERT":
                connection.execute(insert_sql, desired)
            elif action != "NO_CHANGE":
                raise ValidationError(f"unknown planned action: {action}")
            counts[action] += 1
    return dict(counts)


def validate_post_state(connection: sqlite3.Connection, archive: Path) -> dict[str, Any]:
    rows = connection.execute(
        """
        SELECT * FROM structure_file
        WHERE LOWER(source)='modcre' AND status='active' AND file_type='pdb'
        """
    ).fetchall()
    if len(rows) != EXPECTED_MODELS:
        raise ValidationError(f"post-import active ModCRE PDB count is {len(rows)}")
    normalized = [normalize_model_id(str(row["model_id"])) for row in rows]
    if None in normalized or len(set(normalized)) != EXPECTED_MODELS:
        raise ValidationError("post-import active ModCRE identities are not unique")
    if {str(row["archive_path"]) for row in rows} != {str(archive.resolve())}:
        raise ValidationError("not every active ModCRE PDB points to the complete archive")
    kinds = Counter(str(row["model_kind"]) for row in rows)
    interactions = Counter(int(row["dna_interaction"]) for row in rows)
    no_protein = sum(str(row["dna_qc_status"]) == "no_protein_atoms" for row in rows)
    if kinds != Counter({"TFS": EXPECTED_TFS, "DIMER": EXPECTED_DIMERS}):
        raise ValidationError(f"post-import model-kind counts differ: {dict(kinds)}")
    if interactions != Counter({1: EXPECTED_INTERACTION_YES, 0: EXPECTED_INTERACTION_NO}):
        raise ValidationError(f"post-import interaction counts differ: {dict(interactions)}")
    if no_protein != EXPECTED_NO_PROTEIN:
        raise ValidationError(f"post-import DNA-only count differs: {no_protein}")
    af3 = 0
    if table_exists(connection, "structure_model_assignment"):
        af3 = int(
            connection.execute(
                "SELECT COUNT(*) FROM structure_model_assignment WHERE LOWER(source)='af3' AND is_active=1"
            ).fetchone()[0]
        )
        if af3 != EXPECTED_ACTIVE_AF3_ASSIGNMENTS:
            raise ValidationError(f"active AF3 assignments changed: {af3}")
    integrity = str(connection.execute("PRAGMA quick_check").fetchone()[0])
    if integrity != "ok":
        raise ValidationError(f"SQLite quick_check failed: {integrity}")
    return {
        "active_modcre_pdb_rows": len(rows),
        "active_modcre_accessions": len({str(row["tf_id"]) for row in rows}),
        "active_tfs_models": kinds["TFS"],
        "active_dimer_models": kinds["DIMER"],
        "dna_interaction_yes": interactions[1],
        "dna_interaction_no_or_unavailable": interactions[0],
        "dna_only_unavailable": no_protein,
        "active_af3_assignments": af3,
        "quick_check": integrity,
    }


class HashingReader:
    def __init__(self, raw: BinaryIO):
        self.raw = raw
        self.digest = hashlib.sha256()
        self.count = 0

    def read(self, size: int = -1) -> bytes:
        data = self.raw.read(size)
        self.digest.update(data)
        self.count += len(data)
        return data

    def hexdigest(self) -> str:
        return self.digest.hexdigest()


def tar_info(member_name: str, size: int) -> tarfile.TarInfo:
    info = tarfile.TarInfo(member_name)
    info.size = size
    info.mode = 0o644
    info.uid = 0
    info.gid = 0
    info.uname = ""
    info.gname = ""
    info.mtime = 0
    return info


def add_stream_verified(
    archive: tarfile.TarFile,
    member_name: str,
    source: BinaryIO,
    expected_size: int,
    expected_sha256: str,
) -> None:
    reader = HashingReader(source)
    archive.addfile(tar_info(member_name, expected_size), reader)
    if reader.count != expected_size or reader.hexdigest() != expected_sha256:
        raise ValidationError(f"source changed since manifest build: {member_name}")


def build_archive(path: Path, rows: list[dict[str, Any]]) -> int:
    if path.exists():
        validate_archive(path, rows, verify_hashes=True)
        return 0
    path.parent.mkdir(parents=True, exist_ok=True)
    raw_rows = [row for row in rows if row["coordinate_source"] != "legacy_archive"]
    legacy_rows = [row for row in rows if row["coordinate_source"] == "legacy_archive"]
    legacy_by_member = {str(row["source_member_path"]): row for row in legacy_rows}
    legacy_archives = {str(row["source_archive_path"]) for row in legacy_rows}
    if len(legacy_archives) != 1:
        raise ValidationError(f"expected one legacy source archive, found {legacy_archives}")

    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    if temporary.exists():
        raise ValidationError(f"temporary archive already exists: {temporary}")
    try:
        with temporary.open("wb") as raw_output:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw_output, mtime=0) as gzip_output:
                with tarfile.open(fileobj=gzip_output, mode="w") as target:
                    for index, row in enumerate(raw_rows, 1):
                        source_path = Path(str(row["source_path"]))
                        if not source_path.is_file():
                            raise ValidationError(f"missing source PDB: {source_path}")
                        with source_path.open("rb") as source:
                            add_stream_verified(
                                target,
                                str(row["serving_member_path"]),
                                source,
                                int(row["file_size"]),
                                str(row["sha256"]),
                            )
                        if index % 5_000 == 0:
                            print(f"Archived raw PDBs: {index}/{len(raw_rows)}", flush=True)
                    with tarfile.open(next(iter(legacy_archives)), "r:gz") as legacy:
                        written = 0
                        for member in legacy:
                            row = legacy_by_member.get(member.name)
                            if row is None:
                                continue
                            source = legacy.extractfile(member)
                            if source is None:
                                raise ValidationError(f"cannot read legacy member: {member.name}")
                            with source:
                                add_stream_verified(
                                    target,
                                    str(row["serving_member_path"]),
                                    source,
                                    int(row["file_size"]),
                                    str(row["sha256"]),
                                )
                            written += 1
                        if written != len(legacy_rows):
                            raise ValidationError(
                                f"wrote {written}/{len(legacy_rows)} selected legacy members"
                            )
        temporary.replace(path)
    except Exception:
        if temporary.exists():
            temporary.unlink()
        raise
    validate_archive(path, rows, verify_hashes=False)
    return len(rows)


def validate_archive(path: Path, rows: list[dict[str, Any]], verify_hashes: bool) -> None:
    expected = {str(row["serving_member_path"]): row for row in rows}
    observed: set[str] = set()
    duplicates: list[str] = []
    with tarfile.open(path, "r:gz") as archive:
        for member in archive:
            if not member.isfile():
                continue
            if member.name in observed:
                duplicates.append(member.name)
            observed.add(member.name)
            row = expected.get(member.name)
            if row is None:
                raise ValidationError(f"unexpected serving archive member: {member.name}")
            if member.size != int(row["file_size"]):
                raise ValidationError(f"serving archive size mismatch: {member.name}")
            if verify_hashes:
                extracted = archive.extractfile(member)
                if extracted is None:
                    raise ValidationError(f"cannot read serving archive member: {member.name}")
                digest = hashlib.sha256()
                with extracted:
                    for chunk in iter(lambda: extracted.read(1024 * 1024), b""):
                        digest.update(chunk)
                if digest.hexdigest() != row["sha256"]:
                    raise ValidationError(f"serving archive hash mismatch: {member.name}")
    if duplicates:
        raise ValidationError(f"duplicate serving archive members: {duplicates[:10]}")
    missing = set(expected) - observed
    if missing:
        raise ValidationError("missing serving archive members: " + ", ".join(sorted(missing)[:10]))


def write_preview(path: Path, rows: list[dict[str, Any]], replace: bool) -> None:
    if path.exists() and not replace:
        raise ValidationError(f"preview already exists: {path}; use --replace-preview")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=PREVIEW_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({column: row.get(column, "") for column in PREVIEW_COLUMNS} for row in rows)


def backup_database(path: Path, backup_dir: Path) -> Path:
    backup_dir.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    backup = backup_dir / f"{path.name}.before_complete_modcre_{stamp}.bak"
    if backup.exists():
        raise ValidationError(f"backup already exists: {backup}")
    shutil.copy2(path, backup)
    if sha256_file(path) != sha256_file(backup):
        backup.unlink(missing_ok=True)
        raise ValidationError("database backup SHA256 does not match source")
    return backup


def main() -> int:
    args = parse_args()
    dry_run = not args.apply
    if not args.database.is_file():
        raise ValidationError(f"missing staging database: {args.database}")
    before_db_hash = sha256_file(args.database)
    rows, checked_at = read_manifest(args.manifest, args.manifest_qc)

    # Always prove the complete migration/import on an in-memory copy first.
    # Apply mode does not open SQLite read-write until the simulation, archive
    # construction and byte-identical backup have all succeeded.
    simulation = in_memory_copy(args.database)
    before_protected = protected_snapshots(simulation)
    migration_applied = apply_migration_if_needed(simulation, args.migration)
    preview, initial_plan = plan_changes(simulation, rows, args.archive, checked_at)
    write_preview(args.preview_output, preview, args.replace_preview)
    simulated_applied = execute_changes(simulation, preview, args.archive, checked_at)
    simulated_post = validate_post_state(simulation, args.archive)
    after_protected = protected_snapshots(simulation)
    if before_protected != after_protected:
        changed = sorted(
            table for table in PROTECTED_TABLES if before_protected[table] != after_protected[table]
        )
        raise ValidationError("protected tables changed: " + ", ".join(changed))

    second_preview, second_plan = plan_changes(simulation, rows, args.archive, checked_at)
    second_changes = second_plan.get("INSERT", 0) + second_plan.get("UPDATE_EXISTING", 0)
    if second_changes != 0:
        raise ValidationError(f"second simulated run would make {second_changes} changes")
    simulation.close()

    backup: Path | None = None
    archive_additions = 0
    applied = simulated_applied
    post = simulated_post
    if args.apply:
        archive_additions = build_archive(args.archive, rows)
        backup = backup_database(args.database, args.backup_dir)
        connection = sqlite3.connect(args.database)
        connection.row_factory = sqlite3.Row
        connection.execute("PRAGMA foreign_keys = ON")
        real_before_protected = protected_snapshots(connection)
        real_migration_applied = apply_migration_if_needed(connection, args.migration)
        if real_migration_applied != migration_applied:
            raise ValidationError("real migration plan differs from in-memory simulation")
        real_preview, real_plan = plan_changes(connection, rows, args.archive, checked_at)
        if real_plan != initial_plan:
            raise ValidationError("real database plan differs from in-memory simulation")
        applied = execute_changes(connection, real_preview, args.archive, checked_at)
        post = validate_post_state(connection, args.archive)
        real_after_protected = protected_snapshots(connection)
        if real_before_protected != real_after_protected:
            changed = sorted(
                table
                for table in PROTECTED_TABLES
                if real_before_protected[table] != real_after_protected[table]
            )
            raise ValidationError("protected tables changed in real import: " + ", ".join(changed))
        connection.commit()
        connection.close()

    after_db_hash = sha256_file(args.database)
    if dry_run and before_db_hash != after_db_hash:
        raise ValidationError("dry run modified the staging database")

    qc = {
        "mode": "DRY_RUN" if dry_run else "APPLY",
        "database": str(args.database.resolve()),
        "database_sha256_before": before_db_hash,
        "database_sha256_after": after_db_hash,
        "database_modified": before_db_hash != after_db_hash,
        "manifest": str(args.manifest.resolve()),
        "manifest_rows": len(rows),
        "target_archive": str(args.archive.resolve()),
        "archive_exists": args.archive.is_file(),
        "archive_members_added": archive_additions,
        "migration_applied": migration_applied,
        "initial_plan": initial_plan,
        "applied_actions": applied,
        "post_state": post,
        "protected_tables_unchanged": before_protected == after_protected,
        "second_run_plan": second_plan,
        "second_run_changes": second_changes,
        "backup": str(backup.resolve()) if backup else "",
        "all_invariants_passed": True,
    }
    if args.qc_output.exists() and not args.replace_preview:
        raise ValidationError(f"QC output already exists: {args.qc_output}; use --replace-preview")
    args.qc_output.parent.mkdir(parents=True, exist_ok=True)
    args.qc_output.write_text(json.dumps(qc, indent=2, sort_keys=True) + "\n")

    print("Mode:", qc["mode"])
    print("Target active ModCRE models:", post["active_modcre_pdb_rows"])
    print("DIMER / TFS:", post["active_dimer_models"], "/", post["active_tfs_models"])
    print("DNA interaction Yes / other:", post["dna_interaction_yes"], "/", post["dna_interaction_no_or_unavailable"])
    print("Initial plan:", initial_plan)
    print("Second-run changes:", second_changes)
    print("Protected tables unchanged: YES")
    print("Database modified:", "NO" if dry_run else "YES")
    print("Archive created:", "NO" if dry_run else ("YES" if archive_additions else "NO (already valid)"))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ValidationError as exc:
        raise SystemExit(f"ERROR: {exc}")
