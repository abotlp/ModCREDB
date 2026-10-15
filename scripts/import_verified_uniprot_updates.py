#!/usr/bin/env python3
"""Apply the verified UniProt name and merged-accession corrections safely.

Dry-run mode copies the database into memory, applies migration 004 when
needed, simulates all row updates twice, and confirms idempotency.  Apply mode
requires an explicit database path and creates a byte-identical backup first.
No TF, PWM, or structural-model table is modified.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import sqlite3
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATABASE = ROOT / "data" / "tf_webdb_local_staging.sqlite"
DEFAULT_MANIFEST = ROOT / "data_sources" / "uniprot_accession_updates_verified_2026-10-15.tsv"
DEFAULT_MIGRATION = ROOT / "migrations" / "004_add_uniprot_merge_target.sql"
DEFAULT_PREVIEW = ROOT / "outputs" / "uniprot_accession_updates_preview.tsv"
DEFAULT_BACKUPS = ROOT / "backups"

EXPECTED_MERGED = 5
EXPECTED_NAME_UPDATES = 6
EXPECTED_DELETED = 258
EXPECTED_TFS = 5417

PROTECTED_TABLES = (
    "tf",
    "tf_primary_annotation",
    "motif_ref",
    "motif_file",
    "motif_structure",
    "structure_file",
    "structure_model_assignment",
    "structure_confidence_artifact",
    "tf_structure_status",
)

PREVIEW_FIELDS = (
    "accession",
    "update_type",
    "database_gene_names_before",
    "database_gene_names_after",
    "database_protein_name_before",
    "database_protein_name_after",
    "merged_into_accession_before",
    "merged_into_accession_after",
    "official_entry_url",
    "uniprot_release",
    "verified_date",
)


class ImportValidationError(RuntimeError):
    """Raised when an expected invariant does not hold."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def connect_read_only(path: Path) -> sqlite3.Connection:
    connection = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    connection.row_factory = sqlite3.Row
    connection.execute("PRAGMA foreign_keys = ON")
    return connection


def table_exists(connection: sqlite3.Connection, table: str) -> bool:
    return connection.execute(
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name=?", (table,)
    ).fetchone() is not None


def column_exists(connection: sqlite3.Connection, table: str, column: str) -> bool:
    return any(row[1] == column for row in connection.execute(f"PRAGMA table_info({table})"))


def database_health(connection: sqlite3.Connection) -> None:
    result = connection.execute("PRAGMA quick_check").fetchone()[0]
    if result != "ok":
        raise ImportValidationError(f"SQLite quick_check failed: {result}")
    violations = connection.execute("PRAGMA foreign_key_check").fetchall()
    if violations:
        raise ImportValidationError(
            f"SQLite foreign_key_check found {len(violations)} violation(s)"
        )


def canonical_value(value: Any) -> Any:
    if isinstance(value, bytes):
        return {"bytes_sha256": hashlib.sha256(value).hexdigest()}
    return value


def table_digest(connection: sqlite3.Connection, table: str) -> tuple[int, str]:
    if not table_exists(connection, table):
        return 0, "TABLE_NOT_PRESENT"
    columns = [row[1] for row in connection.execute(f"PRAGMA table_info({table})")]
    order = ", ".join(f'"{column}"' for column in columns)
    digest = hashlib.sha256()
    count = 0
    for row in connection.execute(f'SELECT * FROM "{table}" ORDER BY {order}'):
        payload = [canonical_value(value) for value in row]
        digest.update(json.dumps(payload, ensure_ascii=False, separators=(",", ":")).encode())
        digest.update(b"\n")
        count += 1
    return count, digest.hexdigest()


def protected_snapshot(connection: sqlite3.Connection) -> dict[str, tuple[int, str]]:
    return {table: table_digest(connection, table) for table in PROTECTED_TABLES}


def load_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != EXPECTED_MERGED + EXPECTED_NAME_UPDATES:
        raise ImportValidationError(f"expected 11 manifest rows, found {len(rows)}")
    accessions = [row["accession"] for row in rows]
    if len(accessions) != len(set(accessions)):
        raise ImportValidationError("manifest contains duplicate accessions")
    counts: dict[str, int] = {}
    for row in rows:
        counts[row["update_type"]] = counts.get(row["update_type"], 0) + 1
        if row["uniprot_release"] != "2026_03" or row["verified_date"] != "2026-10-15":
            raise ImportValidationError(f"unexpected provenance for {row['accession']}")
    if counts != {"MERGED_ACCESSION": 5, "CURRENT_NAME": 6}:
        raise ImportValidationError(f"unexpected manifest update breakdown: {counts}")
    return rows


def migration_up_sql(path: Path) -> str:
    text = path.read_text(encoding="utf-8")
    if "-- migrate:up" not in text or "-- migrate:down" not in text:
        raise ImportValidationError("migration is missing up/down markers")
    return text.split("-- migrate:up", 1)[1].split("-- migrate:down", 1)[0]


def apply_migration_if_needed(connection: sqlite3.Connection, path: Path) -> bool:
    if column_exists(connection, "tf_annotation", "merged_into_accession"):
        return False
    connection.executescript(migration_up_sql(path))
    if not column_exists(connection, "tf_annotation", "merged_into_accession"):
        raise ImportValidationError("migration did not add merged_into_accession")
    return True


def require_expected_source_row(
    connection: sqlite3.Connection, record: dict[str, str]
) -> sqlite3.Row:
    row = connection.execute(
        """
        SELECT tf_id, gene_names, protein_name, merged_into_accession
        FROM tf_annotation WHERE tf_id=? AND uniprot_accession=?
        """,
        (record["accession"], record["accession"]),
    ).fetchone()
    if row is None:
        raise ImportValidationError(f"missing exact tf_annotation row: {record['accession']}")
    expected_gene = record["database_gene_names"] or None
    expected_protein = record["database_protein_name"] or None
    already_applied = (
        record["update_type"] == "CURRENT_NAME"
        and (row["gene_names"] or None) == (record["official_gene_names"] or None)
        and (row["protein_name"] or None) == (record["official_protein_name"] or None)
    ) or (
        record["update_type"] == "MERGED_ACCESSION"
        and row["merged_into_accession"] == record["merged_into_accession"]
    )
    if not already_applied and (
        (row["gene_names"] or None) != expected_gene
        or (row["protein_name"] or None) != expected_protein
    ):
        raise ImportValidationError(
            f"current metadata for {record['accession']} differs from verified manifest: "
            f"gene={row['gene_names']!r}, protein={row['protein_name']!r}"
        )
    return row


def apply_rows(
    connection: sqlite3.Connection, manifest: list[dict[str, str]]
) -> tuple[int, list[dict[str, object]]]:
    changes = 0
    preview: list[dict[str, object]] = []
    for record in manifest:
        before = require_expected_source_row(connection, record)
        after_gene = before["gene_names"]
        after_protein = before["protein_name"]
        after_merged = before["merged_into_accession"]
        if record["update_type"] == "CURRENT_NAME":
            after_gene = record["official_gene_names"] or None
            after_protein = record["official_protein_name"] or None
            cursor = connection.execute(
                """
                UPDATE tf_annotation
                SET gene_names=?, protein_name=?
                WHERE tf_id=?
                  AND (gene_names IS NOT ? OR protein_name IS NOT ?)
                """,
                (after_gene, after_protein, record["accession"], after_gene, after_protein),
            )
        elif record["update_type"] == "MERGED_ACCESSION":
            target = record["merged_into_accession"]
            if connection.execute("SELECT 1 FROM tf WHERE tf_id=?", (target,)).fetchone() is None:
                raise ImportValidationError(
                    f"merged target {target} for {record['accession']} is not a TF row"
                )
            after_merged = target
            cursor = connection.execute(
                """
                UPDATE tf_annotation SET merged_into_accession=?
                WHERE tf_id=? AND merged_into_accession IS NOT ?
                """,
                (target, record["accession"], target),
            )
        else:
            raise ImportValidationError(f"unsupported update type: {record['update_type']}")
        changes += cursor.rowcount
        preview.append(
            {
                "accession": record["accession"],
                "update_type": record["update_type"],
                "database_gene_names_before": before["gene_names"] or "",
                "database_gene_names_after": after_gene or "",
                "database_protein_name_before": before["protein_name"] or "",
                "database_protein_name_after": after_protein or "",
                "merged_into_accession_before": before["merged_into_accession"] or "",
                "merged_into_accession_after": after_merged or "",
                "official_entry_url": record["official_entry_url"],
                "uniprot_release": record["uniprot_release"],
                "verified_date": record["verified_date"],
            }
        )
    return changes, preview


def validate_final_state(connection: sqlite3.Connection, manifest: list[dict[str, str]]) -> None:
    database_health(connection)
    if connection.execute("SELECT COUNT(*) FROM tf").fetchone()[0] != EXPECTED_TFS:
        raise ImportValidationError("TF count changed")
    if connection.execute("SELECT COUNT(*) FROM tf_annotation").fetchone()[0] != EXPECTED_TFS:
        raise ImportValidationError("tf_annotation count changed")
    if connection.execute(
        "SELECT COUNT(*) FROM tf_annotation WHERE lower(trim(protein_name))='deleted'"
    ).fetchone()[0] != EXPECTED_DELETED:
        raise ImportValidationError("deleted-record count changed")
    if connection.execute(
        "SELECT COUNT(*) FROM tf_annotation WHERE lower(trim(protein_name))='merged'"
    ).fetchone()[0] != EXPECTED_MERGED:
        raise ImportValidationError("merged sentinel count changed")
    if connection.execute(
        "SELECT COUNT(*) FROM tf_annotation WHERE merged_into_accession IS NOT NULL"
    ).fetchone()[0] != EXPECTED_MERGED:
        raise ImportValidationError("merged-target count is not five")
    for record in manifest:
        row = connection.execute(
            "SELECT gene_names,protein_name,merged_into_accession FROM tf_annotation WHERE tf_id=?",
            (record["accession"],),
        ).fetchone()
        if record["update_type"] == "CURRENT_NAME":
            if (row[0] or "") != record["official_gene_names"] or row[1] != record["official_protein_name"]:
                raise ImportValidationError(f"name update failed: {record['accession']}")
        elif row[2] != record["merged_into_accession"] or row[1] != "merged":
            raise ImportValidationError(f"merged update failed: {record['accession']}")


def write_preview(path: Path, rows: list[dict[str, object]], overwrite: bool) -> None:
    if path.exists() and not overwrite:
        raise ImportValidationError(f"refusing to overwrite preview: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=PREVIEW_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--dry-run", action="store_true")
    mode.add_argument("--apply", action="store_true")
    parser.add_argument("--database", type=Path, default=DEFAULT_DATABASE)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--migration", type=Path, default=DEFAULT_MIGRATION)
    parser.add_argument("--preview-output", type=Path, default=DEFAULT_PREVIEW)
    parser.add_argument("--backup-dir", type=Path, default=DEFAULT_BACKUPS)
    parser.add_argument("--overwrite-preview", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not args.database.is_file():
        raise ImportValidationError(f"database not found: {args.database}")
    manifest = load_manifest(args.manifest)
    database_sha_before = sha256_file(args.database)
    source = connect_read_only(args.database)
    try:
        database_health(source)
        before_protected = protected_snapshot(source)
        if args.dry_run:
            working = sqlite3.connect(":memory:")
            working.row_factory = sqlite3.Row
            source.backup(working)
    finally:
        source.close()

    if args.dry_run:
        working.execute("PRAGMA foreign_keys = ON")
        backup_path = None
    else:
        timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        args.backup_dir.mkdir(parents=True, exist_ok=True)
        backup_path = args.backup_dir / f"{args.database.name}.before_uniprot_updates_{timestamp}.bak"
        if backup_path.exists():
            raise ImportValidationError(f"backup already exists: {backup_path}")
        shutil.copy2(args.database, backup_path)
        if sha256_file(backup_path) != database_sha_before:
            raise ImportValidationError("backup SHA256 does not match source database")
        working = sqlite3.connect(args.database)
        working.row_factory = sqlite3.Row
        working.execute("PRAGMA foreign_keys = ON")

    try:
        migration_applied = apply_migration_if_needed(working, args.migration)
        first_changes, preview = apply_rows(working, manifest)
        validate_final_state(working, manifest)
        after_first_protected = protected_snapshot(working)
        if before_protected != after_first_protected:
            changed = [table for table in PROTECTED_TABLES if before_protected[table] != after_first_protected[table]]
            raise ImportValidationError(f"protected tables changed: {changed}")
        second_changes, _ = apply_rows(working, manifest)
        validate_final_state(working, manifest)
        if second_changes != 0:
            raise ImportValidationError(f"second run proposed {second_changes} row changes")
        if args.apply:
            working.commit()
    except Exception:
        if args.apply:
            working.rollback()
        raise
    finally:
        working.close()

    write_preview(args.preview_output, preview, args.overwrite_preview)
    database_sha_after = sha256_file(args.database)
    if args.dry_run and database_sha_after != database_sha_before:
        raise ImportValidationError("dry run changed the staging database")

    print("Mode:", "DRY-RUN" if args.dry_run else "APPLY")
    print("Migration required:", "YES" if migration_applied else "NO")
    print("Verified name rows:", EXPECTED_NAME_UPDATES)
    print("Verified merged rows:", EXPECTED_MERGED)
    print("First-run row changes:", first_changes)
    print("Second-run row changes:", second_changes)
    print("Deleted rows unchanged:", EXPECTED_DELETED)
    print("Protected tables unchanged: YES")
    print("Database SHA256 before:", database_sha_before)
    print("Database SHA256 after:", database_sha_after)
    print("Backup:", str(backup_path) if backup_path else "<not created in dry-run>")
    print("Preview:", args.preview_output)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ImportValidationError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
