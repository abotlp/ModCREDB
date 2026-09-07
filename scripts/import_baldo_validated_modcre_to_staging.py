#!/usr/bin/env python3
"""Package and index the frozen, interface-positive Baldo ModCRE model set.

The dry-run path is intentionally read-only for both SQLite and the target
archive.  The apply path is provided for a later, separately authorized run.
It creates one dedicated archive and inserts only missing ``structure_file``
rows; it never updates existing structure rows, PWM tables, model summaries,
or persisted TF counters.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import os
from pathlib import Path
import re
import shutil
import sqlite3
import sys
import tarfile
import tempfile
from typing import Any, Iterable


REPO_ROOT = Path(__file__).resolve().parents[1]
TRUTH_PATH = REPO_ROOT / "outputs/baldo_3690_STRUCTURAL_TRUTH_CURRENT.tsv"
AUDIT_PATH = REPO_ROOT / "outputs/baldo_3479_modcre_interface_audit.tsv"
STAGING_DB = REPO_ROOT / "data/tf_webdb_local_staging.sqlite"
PRODUCTION_DB = REPO_ROOT / "data/tf_webdb.sqlite"
PREVIEW_PATH = REPO_ROOT / "outputs/baldo_3468_modcre_import_preview.tsv"
MODEL_ROOT = Path("/data/sbi/interchange/boliva/patricia/models")
ARCHIVE_PATH = Path(
    "/home/patricia/TF_database_Baldo_data/baldo_validated_modcre_models.tar.gz"
)

EXPECTED_ACCEPTED = 3468
EXPECTED_EMPTY_INTERFACE = 11
EXPECTED_LEGACY_NON_BALDO = 1083
PROTECTED_TABLES = (
    "motif_ref",
    "motif_file",
    "motif_structure",
    "tf_primary_annotation",
)
FILENAME_RE = re.compile(
    r"^(?P<prefix>DIMER|TFS)_(?P<accession>[^:]+):"
    r"(?P<start>[0-9]+):(?P<end>[0-9]+)_(?P<tail>.+)[.]pdb$"
)

PREVIEW_COLUMNS = (
    "accession",
    "tf_id",
    "structural_status",
    "passing_model_path",
    "model_id",
    "member_path",
    "archive_path",
    "residue_start",
    "residue_end",
    "template_pdb",
    "valid_interface_any_model",
    "atom_contacts",
    "protein_interface_residues",
    "dna_interface_residues",
    "file_exists",
    "file_readable",
    "filename_accession_matches",
    "tf_mapping_count",
    "existing_exact_structure_file_rows",
    "model_id_collision",
    "member_path_collision",
    "planned_structure_file_action",
    "planned_archive_member_action",
)


class ValidationError(RuntimeError):
    """Raised when frozen inputs or a target invariant are not satisfied."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Validate or import the frozen 3468-accession Baldo ModCRE set "
            "into the local staging database."
        )
    )
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument(
        "--dry-run",
        action="store_true",
        help="validate and write only the TSV preview; do not alter SQLite or the archive",
    )
    mode.add_argument(
        "--apply",
        action="store_true",
        help="create/reuse the dedicated archive and insert missing staging rows",
    )
    parser.add_argument(
        "--database",
        type=Path,
        default=STAGING_DB,
        help=f"staging SQLite path (must resolve to {STAGING_DB})",
    )
    parser.add_argument(
        "--archive",
        type=Path,
        default=ARCHIVE_PATH,
        help=f"dedicated archive path (must resolve to {ARCHIVE_PATH})",
    )
    parser.add_argument(
        "--preview-output",
        type=Path,
        default=PREVIEW_PATH,
        help=f"dry-run TSV preview path (must resolve to {PREVIEW_PATH})",
    )
    return parser.parse_args()


def require_exact_target(actual: Path, expected: Path, label: str) -> Path:
    actual_resolved = actual.expanduser().resolve()
    expected_resolved = expected.resolve()
    if actual_resolved != expected_resolved:
        raise ValidationError(
            f"{label} must be {expected_resolved}, not {actual_resolved}"
        )
    return actual_resolved


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file():
        raise ValidationError(f"required input is missing: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValidationError(f"TSV has no header: {path}")
        return list(reader.fieldnames), [dict(row) for row in reader]


def require_columns(path: Path, headers: Iterable[str], required: Iterable[str]) -> None:
    missing = sorted(set(required) - set(headers))
    if missing:
        raise ValidationError(f"{path} lacks required columns: {', '.join(missing)}")


def truth_status_column(headers: list[str], rows: list[dict[str, str]]) -> str:
    candidates = [
        name
        for name in ("structural_status", "current_structural_status")
        if name in headers
    ]
    if not candidates:
        raise ValidationError(
            f"{TRUTH_PATH} lacks structural_status/current_structural_status"
        )
    if len(candidates) == 2:
        inconsistent = [
            row.get("accession", "")
            for row in rows
            if row["structural_status"] != row["current_structural_status"]
        ]
        if inconsistent:
            raise ValidationError(
                "truth status columns disagree for: " + ", ".join(inconsistent[:10])
            )
        return "structural_status"
    return candidates[0]


def duplicate_values(values: Iterable[str]) -> list[str]:
    seen: set[str] = set()
    duplicates: set[str] = set()
    for value in values:
        if value in seen:
            duplicates.add(value)
        seen.add(value)
    return sorted(duplicates)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_fingerprint(path: Path) -> tuple[int, int, str]:
    stat_result = path.stat()
    return stat_result.st_size, stat_result.st_mtime_ns, sha256_file(path)


def open_database(path: Path, *, writable: bool) -> sqlite3.Connection:
    if writable:
        connection = sqlite3.connect(path)
    else:
        connection = sqlite3.connect(
            f"file:{path}?mode=ro&immutable=1",
            uri=True,
        )
    connection.row_factory = sqlite3.Row
    connection.execute("PRAGMA foreign_keys = ON")
    return connection


def table_columns(connection: sqlite3.Connection, table: str) -> list[str]:
    return [str(row["name"]) for row in connection.execute(f"PRAGMA table_info({table})")]


def validate_schema(connection: sqlite3.Connection) -> None:
    required: dict[str, set[str]] = {
        "tf": {"tf_id", "active_model_count"},
        "source": {"source"},
        "structure_file": {
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
        },
        "model_summary": {"id", "summary_file_id", "matched_structure_id"},
        "motif_ref": {"id"},
        "motif_file": {"source", "motif_id"},
        "motif_structure": {"motif_ref_id", "structure_file_id"},
        "tf_primary_annotation": {"tf_id"},
    }
    for table, required_columns in required.items():
        columns = set(table_columns(connection, table))
        if not columns:
            raise ValidationError(f"staging database lacks table {table}")
        missing = sorted(required_columns - columns)
        if missing:
            raise ValidationError(
                f"staging table {table} lacks columns: {', '.join(missing)}"
            )
    source_count = connection.execute(
        "SELECT COUNT(*) FROM source WHERE source = 'modcre'"
    ).fetchone()[0]
    if source_count != 1:
        raise ValidationError(
            f"expected exactly one source='modcre' row, found {source_count}"
        )


def encoded_value(value: Any) -> bytes:
    if value is None:
        return b"N"
    if isinstance(value, bytes):
        payload = value
        prefix = b"B"
    else:
        payload = str(value).encode("utf-8")
        prefix = b"T"
    return prefix + str(len(payload)).encode("ascii") + b":" + payload


def table_snapshot(connection: sqlite3.Connection, table: str) -> tuple[int, str]:
    columns = table_columns(connection, table)
    if not columns:
        raise ValidationError(f"cannot snapshot absent table: {table}")
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


def load_selected_records() -> tuple[list[dict[str, Any]], set[str]]:
    truth_headers, truth_rows = read_tsv(TRUTH_PATH)
    require_columns(TRUTH_PATH, truth_headers, ("accession",))
    status_column = truth_status_column(truth_headers, truth_rows)
    truth_duplicates = duplicate_values(row["accession"].strip() for row in truth_rows)
    if truth_duplicates:
        raise ValidationError(
            "duplicate accessions in structural truth: " + ", ".join(truth_duplicates[:10])
        )

    accepted = {
        row["accession"].strip(): row
        for row in truth_rows
        if row[status_column].strip() == "VALID_MODCRE"
    }
    rejected_empty = {
        row["accession"].strip()
        for row in truth_rows
        if row[status_column].strip() == "MODCRE_EMPTY_INTERFACE"
    }
    if len(accepted) != EXPECTED_ACCEPTED:
        raise ValidationError(
            f"expected {EXPECTED_ACCEPTED} VALID_MODCRE accessions, found {len(accepted)}"
        )
    if len(rejected_empty) != EXPECTED_EMPTY_INTERFACE:
        raise ValidationError(
            f"expected {EXPECTED_EMPTY_INTERFACE} MODCRE_EMPTY_INTERFACE accessions, "
            f"found {len(rejected_empty)}"
        )
    if set(accepted) & rejected_empty:
        raise ValidationError("VALID_MODCRE and MODCRE_EMPTY_INTERFACE sets overlap")

    audit_headers, audit_rows = read_tsv(AUDIT_PATH)
    require_columns(
        AUDIT_PATH,
        audit_headers,
        (
            "accession",
            "valid_interface_any_model",
            "passing_model_path",
            "atom_contacts",
            "protein_interface_residues",
            "dna_interface_residues",
        ),
    )
    audit_duplicates = duplicate_values(row["accession"].strip() for row in audit_rows)
    if audit_duplicates:
        raise ValidationError(
            "duplicate accessions in interface audit: " + ", ".join(audit_duplicates[:10])
        )
    audit_by_accession = {row["accession"].strip(): row for row in audit_rows}
    missing_audit = sorted(set(accepted) - set(audit_by_accession))
    if missing_audit:
        raise ValidationError(
            "accepted accessions missing from interface audit: "
            + ", ".join(missing_audit[:10])
        )

    model_root = MODEL_ROOT.resolve()
    records: list[dict[str, Any]] = []
    for accession in sorted(accepted):
        audit = audit_by_accession[accession]
        if audit["valid_interface_any_model"].strip() != "YES":
            raise ValidationError(f"accepted accession lacks passing audit result: {accession}")
        raw_path = audit["passing_model_path"].strip()
        if not raw_path:
            raise ValidationError(f"accepted accession lacks passing_model_path: {accession}")
        model_path = Path(raw_path)
        if not model_path.is_absolute():
            raise ValidationError(f"passing_model_path is not absolute for {accession}: {raw_path}")
        resolved_path = model_path.resolve()
        if resolved_path.parent != model_root:
            raise ValidationError(
                f"selected PDB is not directly under {model_root}: {resolved_path}"
            )
        match = FILENAME_RE.fullmatch(resolved_path.name)
        if match is None:
            raise ValidationError(
                f"selected filename does not match DIMER/TFS convention: {resolved_path.name}"
            )
        embedded_accession = match.group("accession")
        if embedded_accession != accession:
            raise ValidationError(
                f"filename accession mismatch: expected {accession}, found {embedded_accession}"
            )
        residue_start = int(match.group("start"))
        residue_end = int(match.group("end"))
        if residue_start < 1 or residue_end < residue_start:
            raise ValidationError(
                f"invalid residue range for {accession}: {residue_start}-{residue_end}"
            )
        if not resolved_path.is_file():
            raise ValidationError(f"selected PDB does not exist: {resolved_path}")
        try:
            with resolved_path.open("rb") as handle:
                handle.read(1)
        except OSError as exc:
            raise ValidationError(f"selected PDB is not readable: {resolved_path}: {exc}") from exc

        model_id = resolved_path.stem
        template_pdb = match.group("tail").split("_", 1)[0]
        records.append(
            {
                "accession": accession,
                "tf_id": accession,
                "structural_status": "VALID_MODCRE",
                "passing_model_path": str(resolved_path),
                "model_id": model_id,
                "member_path": f"models/{model_id}.pdb",
                "archive_path": str(ARCHIVE_PATH),
                "residue_start": residue_start,
                "residue_end": residue_end,
                "template_pdb": template_pdb,
                "valid_interface_any_model": audit["valid_interface_any_model"].strip(),
                "atom_contacts": audit["atom_contacts"].strip(),
                "protein_interface_residues": audit["protein_interface_residues"].strip(),
                "dna_interface_residues": audit["dna_interface_residues"].strip(),
                "file_exists": "YES",
                "file_readable": "YES",
                "filename_accession_matches": "YES",
            }
        )

    for label, values in (
        ("passing_model_path", (row["passing_model_path"] for row in records)),
        ("model_id", (row["model_id"] for row in records)),
        ("member_path", (row["member_path"] for row in records)),
    ):
        duplicates = duplicate_values(values)
        if duplicates:
            raise ValidationError(
                f"duplicate selected {label} values: " + ", ".join(duplicates[:10])
            )
    if len(records) != EXPECTED_ACCEPTED:
        raise ValidationError(
            f"expected {EXPECTED_ACCEPTED} selected PDBs, found {len(records)}"
        )
    return records, rejected_empty


def desired_structure_values(record: dict[str, Any]) -> tuple[Any, ...]:
    return (
        "modcre",
        record["model_id"],
        record["tf_id"],
        record["member_path"],
        record["archive_path"],
        "pdb",
        "active",
        record["template_pdb"],
        record["residue_start"],
        record["residue_end"],
    )


def existing_structure_values(row: sqlite3.Row) -> tuple[Any, ...]:
    return (
        row["source"],
        row["model_id"],
        row["tf_id"],
        row["member_path"],
        row["archive_path"],
        row["file_type"],
        row["status"],
        row["template_pdb"],
        row["residue_start"],
        row["residue_end"],
    )


def enrich_database_plan(
    connection: sqlite3.Connection,
    records: list[dict[str, Any]],
    rejected_empty: set[str],
) -> dict[str, Any]:
    validate_schema(connection)
    tf_ids = {str(row["tf_id"]) for row in connection.execute("SELECT tf_id FROM tf")}
    unresolved = sorted(row["accession"] for row in records if row["accession"] not in tf_ids)
    for row in records:
        row["tf_mapping_count"] = 1 if row["accession"] in tf_ids else 0
    if unresolved:
        raise ValidationError(
            "accepted accessions without a unique direct tf.tf_id mapping: "
            + ", ".join(unresolved[:10])
        )

    existing_rows = connection.execute("SELECT * FROM structure_file").fetchall()
    by_model_id: dict[str, list[sqlite3.Row]] = {}
    by_member_path: dict[str, list[sqlite3.Row]] = {}
    for existing in existing_rows:
        by_model_id.setdefault(str(existing["model_id"]), []).append(existing)
        by_member_path.setdefault(str(existing["member_path"]), []).append(existing)

    model_collision_count = 0
    member_collision_count = 0
    insert_count = 0
    exact_existing_ids: set[int] = set()
    selected_model_ids = {str(row["model_id"]) for row in records}
    for record in records:
        desired = desired_structure_values(record)
        model_rows = by_model_id.get(str(record["model_id"]), [])
        member_rows = by_member_path.get(str(record["member_path"]), [])
        exact_model_rows = [row for row in model_rows if existing_structure_values(row) == desired]
        exact_member_rows = [row for row in member_rows if existing_structure_values(row) == desired]
        exact_ids = {int(row["id"]) for row in exact_model_rows} & {
            int(row["id"]) for row in exact_member_rows
        }
        model_collision = any(existing_structure_values(row) != desired for row in model_rows)
        member_collision = any(existing_structure_values(row) != desired for row in member_rows)
        if len(exact_ids) > 1:
            model_collision = True
            member_collision = True
        record["existing_exact_structure_file_rows"] = len(exact_ids)
        record["model_id_collision"] = "YES" if model_collision else "NO"
        record["member_path_collision"] = "YES" if member_collision else "NO"
        if model_collision:
            model_collision_count += 1
        if member_collision:
            member_collision_count += 1
        if model_collision or member_collision:
            record["planned_structure_file_action"] = "BLOCKED_COLLISION"
        elif len(exact_ids) == 1:
            exact_existing_ids.update(exact_ids)
            record["planned_structure_file_action"] = "EXISTING_EXACT_NO_CHANGE"
        else:
            insert_count += 1
            record["planned_structure_file_action"] = "INSERT"

    if model_collision_count or member_collision_count:
        raise ValidationError(
            f"incompatible structure_file collisions: model_id={model_collision_count}, "
            f"member_path={member_collision_count}"
        )

    accepted = {str(row["accession"]) for row in records}
    active_modcre_rows = connection.execute(
        """
        SELECT id, tf_id, model_id
        FROM structure_file
        WHERE source = 'modcre' AND status = 'active' AND file_type = 'pdb'
        """
    ).fetchall()
    active_modcre_accessions = {str(row["tf_id"]) for row in active_modcre_rows}
    legacy_non_baldo = active_modcre_accessions - accepted
    if len(legacy_non_baldo) != EXPECTED_LEGACY_NON_BALDO:
        raise ValidationError(
            f"expected {EXPECTED_LEGACY_NON_BALDO} legacy non-Baldo active ModCRE "
            f"accessions, found {len(legacy_non_baldo)}"
        )
    unexpected_selected_rows = [
        row
        for row in active_modcre_rows
        if str(row["tf_id"]) in accepted
        and not (
            int(row["id"]) in exact_existing_ids
            and str(row["model_id"]) in selected_model_ids
        )
    ]
    if unexpected_selected_rows:
        raise ValidationError(
            "accepted accessions have active non-selected ModCRE rows; existing rows "
            "would need adjudication before this insert-only importer can proceed"
        )
    rejected_active = sorted(rejected_empty & active_modcre_accessions)
    if rejected_active:
        raise ValidationError(
            "MODCRE_EMPTY_INTERFACE accessions unexpectedly active in staging: "
            + ", ".join(rejected_active)
        )

    currently_covered = accepted & active_modcre_accessions
    post_covered = currently_covered | {
        str(row["accession"])
        for row in records
        if row["planned_structure_file_action"] in {
            "INSERT",
            "EXISTING_EXACT_NO_CHANGE",
        }
    }
    if len(post_covered) != EXPECTED_ACCEPTED:
        raise ValidationError(
            f"post-import accepted coverage would be {len(post_covered)}, not {EXPECTED_ACCEPTED}"
        )
    return {
        "unresolved_tf_mappings": unresolved,
        "model_id_collisions": model_collision_count,
        "member_path_collisions": member_collision_count,
        "structure_file_inserts": insert_count,
        "newly_covered_accessions": len(post_covered - currently_covered),
        "post_covered_accessions": len(post_covered),
        "legacy_non_baldo_accessions": len(legacy_non_baldo),
        "rejected_empty_active": len(rejected_active),
    }


def archive_file_hashes(path: Path) -> tuple[dict[str, str], list[str]]:
    hashes: dict[str, str] = {}
    duplicates: list[str] = []
    with tarfile.open(path, "r:gz") as archive:
        for member in archive.getmembers():
            if not member.isfile():
                continue
            if member.name in hashes:
                duplicates.append(member.name)
                continue
            extracted = archive.extractfile(member)
            if extracted is None:
                raise ValidationError(f"cannot read archive member: {member.name}")
            digest = hashlib.sha256()
            for block in iter(lambda: extracted.read(1024 * 1024), b""):
                digest.update(block)
            hashes[member.name] = digest.hexdigest()
    return hashes, sorted(set(duplicates))


def plan_archive(records: list[dict[str, Any]], archive_path: Path) -> dict[str, Any]:
    planned_names = {str(row["member_path"]) for row in records}
    if archive_path.exists():
        if not archive_path.is_file():
            raise ValidationError(f"archive target exists but is not a file: {archive_path}")
        archive_hashes, duplicate_members = archive_file_hashes(archive_path)
        if duplicate_members:
            raise ValidationError(
                "dedicated archive has duplicate members: " + ", ".join(duplicate_members[:10])
            )
        if set(archive_hashes) != planned_names:
            missing = sorted(planned_names - set(archive_hashes))
            extra = sorted(set(archive_hashes) - planned_names)
            raise ValidationError(
                f"existing dedicated archive does not match frozen set; "
                f"missing={len(missing)}, extra={len(extra)}"
            )
        by_member = {str(row["member_path"]): row for row in records}
        mismatches = [
            member
            for member, digest in archive_hashes.items()
            if sha256_file(Path(by_member[member]["passing_model_path"])) != digest
        ]
        if mismatches:
            raise ValidationError(
                "existing archive member content differs from selected PDB: "
                + ", ".join(mismatches[:10])
            )
        action = "EXISTING_EXACT_NO_CHANGE"
    else:
        action = "CREATE"
    for record in records:
        record["planned_archive_member_action"] = (
            "EXISTING_EXACT_NO_CHANGE" if action == "EXISTING_EXACT_NO_CHANGE" else "CREATE"
        )
    return {
        "archive_action": action,
        "archive_members_planned": 0 if action == "EXISTING_EXACT_NO_CHANGE" else len(records),
        "archive_created": False,
    }


def render_preview(records: list[dict[str, Any]]) -> str:
    output = io.StringIO(newline="")
    writer = csv.DictWriter(
        output,
        fieldnames=PREVIEW_COLUMNS,
        delimiter="\t",
        lineterminator="\n",
        extrasaction="ignore",
    )
    writer.writeheader()
    writer.writerows(records)
    return output.getvalue()


def write_preview_without_overwrite(path: Path, content: str) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    encoded = content.encode("utf-8")
    if path.exists():
        if path.is_file() and path.read_bytes() == encoded:
            return "EXISTING_IDENTICAL_NO_CHANGE"
        raise ValidationError(f"refusing to overwrite non-identical preview: {path}")
    try:
        descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(encoded)
    except Exception:
        if path.exists():
            path.unlink()
        raise
    return "CREATED"


def create_archive_temp(records: list[dict[str, Any]], archive_path: Path) -> Path:
    archive_path.parent.mkdir(parents=True, exist_ok=True)
    handle = tempfile.NamedTemporaryFile(
        prefix=archive_path.name + ".",
        suffix=".tmp",
        dir=archive_path.parent,
        delete=False,
    )
    temporary_path = Path(handle.name)
    try:
        with handle:
            with gzip.GzipFile(filename="", mode="wb", fileobj=handle, mtime=0) as compressed:
                with tarfile.open(fileobj=compressed, mode="w", format=tarfile.PAX_FORMAT) as archive:
                    for record in records:
                        source_path = Path(record["passing_model_path"])
                        info = tarfile.TarInfo(str(record["member_path"]))
                        info.size = source_path.stat().st_size
                        info.mode = 0o644
                        info.mtime = 0
                        info.uid = 0
                        info.gid = 0
                        info.uname = ""
                        info.gname = ""
                        with source_path.open("rb") as source:
                            archive.addfile(info, source)
        archive_hashes, duplicates = archive_file_hashes(temporary_path)
        if duplicates or set(archive_hashes) != {
            str(row["member_path"]) for row in records
        }:
            raise ValidationError("new archive failed member-list validation")
        return temporary_path
    except Exception:
        temporary_path.unlink(missing_ok=True)
        raise


def apply_import(
    database_path: Path,
    archive_path: Path,
    records: list[dict[str, Any]],
    rejected_empty: set[str],
    protected_before: dict[str, tuple[int, str]],
) -> int:
    archive_created = False
    temporary_archive: Path | None = None
    if not archive_path.exists():
        temporary_archive = create_archive_temp(records, archive_path)

    connection = open_database(database_path, writable=True)
    try:
        connection.execute("BEGIN IMMEDIATE")
        transaction_records = [dict(row) for row in records]
        transaction_plan = enrich_database_plan(
            connection, transaction_records, rejected_empty
        )
        insert_records = [
            row
            for row in transaction_records
            if row["planned_structure_file_action"] == "INSERT"
        ]
        if transaction_plan["structure_file_inserts"] != len(insert_records):
            raise ValidationError("transaction insert plan changed unexpectedly")
        connection.executemany(
            """
            INSERT INTO structure_file (
                source, model_id, tf_id, member_path, archive_path,
                file_type, status, template_pdb, residue_start, residue_end
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            [desired_structure_values(row) for row in insert_records],
        )
        if protected_snapshots(connection) != protected_before:
            raise ValidationError("a protected PWM/provenance table changed in transaction")
        post_records = [dict(row) for row in records]
        post_plan = enrich_database_plan(connection, post_records, rejected_empty)
        if post_plan["structure_file_inserts"] != 0:
            raise ValidationError("post-insert plan is not idempotent")
        if temporary_archive is not None:
            if archive_path.exists():
                raise ValidationError(f"refusing to overwrite archive: {archive_path}")
            os.replace(temporary_archive, archive_path)
            temporary_archive = None
            archive_created = True
        connection.commit()
        return len(insert_records)
    except Exception:
        connection.rollback()
        if archive_created:
            archive_path.unlink(missing_ok=True)
        raise
    finally:
        connection.close()
        if temporary_archive is not None:
            temporary_archive.unlink(missing_ok=True)


def main() -> int:
    args = parse_args()
    database_path = require_exact_target(args.database, STAGING_DB, "database")
    archive_path = require_exact_target(args.archive, ARCHIVE_PATH, "archive")
    preview_path = require_exact_target(args.preview_output, PREVIEW_PATH, "preview output")
    if database_path == PRODUCTION_DB.resolve():
        raise ValidationError("production database is never a valid target")
    if not database_path.is_file():
        raise ValidationError(f"staging database is missing: {database_path}")

    database_before = file_fingerprint(database_path)
    records, rejected_empty = load_selected_records()
    with open_database(database_path, writable=False) as connection:
        protected_before = protected_snapshots(connection)
        database_plan = enrich_database_plan(connection, records, rejected_empty)
    archive_plan = plan_archive(records, archive_path)

    if args.dry_run:
        preview_action = write_preview_without_overwrite(
            preview_path, render_preview(records)
        )
        with open_database(database_path, writable=False) as connection:
            protected_after = protected_snapshots(connection)
        if protected_after != protected_before:
            raise ValidationError("protected PWM/provenance tables changed during dry-run")
        if file_fingerprint(database_path) != database_before:
            raise ValidationError("staging database changed during dry-run")
        if archive_plan["archive_action"] == "CREATE" and archive_path.exists():
            raise ValidationError("dry-run unexpectedly created the dedicated archive")
        print(
            "dry-run: "
            f"preview={preview_action}; accepted={len(records)}; "
            f"selected_pdbs={len({row['passing_model_path'] for row in records})}; "
            f"inserts={database_plan['structure_file_inserts']}; "
            f"archive_members={archive_plan['archive_members_planned']}; "
            f"post_coverage={database_plan['post_covered_accessions']}"
        )
    else:
        inserted = apply_import(
            database_path,
            archive_path,
            records,
            rejected_empty,
            protected_before,
        )
        print(f"apply: inserted_structure_file_rows={inserted}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ValidationError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)
