#!/usr/bin/env python3
"""Finalize the curated status rows for the 108 unresolved Baldo accessions.

Dry-run mode copies the staging database into SQLite memory, applies migration
002 there, performs the status import twice, and verifies that the second pass
is a no-op.  The on-disk staging database is opened read-only in dry-run mode.
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import io
import os
from pathlib import Path
import sqlite3
import sys
from typing import Any, Iterable


REPO_ROOT = Path(__file__).resolve().parents[1]
STAGING_DB = REPO_ROOT / "data/tf_webdb_local_staging.sqlite"
PRODUCTION_DB = REPO_ROOT / "data/tf_webdb.sqlite"
MIGRATION_PATH = REPO_ROOT / "migrations/002_finalize_baldo_unresolved_status.sql"
FINAL_ADJUDICATION_PATH = (
    REPO_ROOT / "outputs/baldo_109_final_structural_adjudication.tsv"
)
UNRESOLVED_REVIEW_PATH = REPO_ROOT / "outputs/baldo_109_unresolved_model_review.tsv"
STRUCTURAL_TRUTH_PATH = REPO_ROOT / "outputs/baldo_3690_STRUCTURAL_TRUTH_CURRENT.tsv"
PREVIEW_PATH = REPO_ROOT / "outputs/baldo_final_108_structure_status_preview.tsv"

SOURCE_UNRESOLVED_STATUS = "PWM_PRESENT_NO_VALID_DNA_BOUND_STRUCTURE"
PUBLIC_STRUCTURAL_STATUS = "PWM_PRESENT_NO_VALIDATED_DNA_BOUND_STRUCTURE"
DISPLAY_RECOMMENDATION = "PWM available; no validated DNA-bound structure"
FINAL_REVIEW_STATUS = "FINAL_BALDO_ADJUDICATED"
SOURCE_PROVENANCE = (
    "outputs/baldo_109_final_structural_adjudication.tsv;"
    "outputs/baldo_109_unresolved_model_review.tsv"
)

EXPECTED_UNRESOLVED = 108
EXPECTED_IDENTICAL_PWM = 60
EXPECTED_HOMOLOGOUS_PWM = 48
EXPECTED_EXISTING_FINAL = 98
EXPECTED_MISSING_FINAL = 10
STALE_ACCESSIONS = {"A0A494C1L3", "H2BNB9", "H7C4N4", "Q709A9"}
NFIB_ACCESSION = "A0A1B0GWI9"

OWN_STATUS_BY_DBD_STATUS = {
    "ACCESSION_LACKS_DBD": "ACCESSION_LACKS_DBD",
    "ACCESSION_CONTAINS_INCOMPLETE_DBD": "ACCESSION_CONTAINS_INCOMPLETE_DBD",
    "COFACTOR_OR_NON_DNA_BINDING_COMPONENT": "NON_DNA_BINDING_COMPONENT",
    "NO_COMPLETE_DBD_DEMONSTRATED": "UNRESOLVED",
}
EXPECTED_OWN_STATUS_COUNTS = {
    "ACCESSION_LACKS_DBD": 76,
    "ACCESSION_CONTAINS_INCOMPLETE_DBD": 10,
    "NON_DNA_BINDING_COMPONENT": 2,
    "UNRESOLVED": 20,
}

PROTECTED_TABLES = (
    "motif_ref",
    "motif_file",
    "motif_structure",
    "tf_primary_annotation",
    "structure_file",
    "structure_model_assignment",
)
STATUS_COLUMNS = (
    "tf_id",
    "uniprot_accession",
    "primary_structural_status",
    "own_accession_structure_status",
    "canonical_reference_status",
    "canonical_reference_accession",
    "canonical_reference_model_path",
    "database_display_recommendation",
    "action_for_baldo",
    "decision_reason",
    "remaining_uncertainty",
    "review_status",
    "source_audit_file",
    "created_at",
    "updated_at",
)
PREVIEW_COLUMNS = (
    "tf_id",
    "uniprot_accession",
    "gene",
    "PWM_status",
    "source_final_structural_status",
    "source_DBD_status",
    "expected_DBD",
    "DBD_evidence_coordinates",
    "DBD_coverage_percent",
    "DBD_identity_percent",
    "accession_coordinates_matching_expected_DBD",
    "pfam_review_category",
    "primary_structural_status",
    "own_accession_structure_status",
    "canonical_reference_status",
    "canonical_reference_accession",
    "canonical_reference_model_path",
    "database_display_recommendation",
    "action_for_baldo",
    "decision_reason",
    "remaining_uncertainty",
    "review_status",
    "source_audit_file",
    "initial_import_action",
)


class ValidationError(RuntimeError):
    """Raised when an input or database invariant is not satisfied."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Finalize the 108 curated Baldo unresolved structure-status rows."
    )
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument(
        "--dry-run",
        action="store_true",
        help="apply migration/import twice to an in-memory database copy only",
    )
    mode.add_argument(
        "--apply",
        action="store_true",
        help="apply migration/import to the local staging database",
    )
    parser.add_argument(
        "--database",
        type=Path,
        default=STAGING_DB,
        help=f"database target (must resolve to {STAGING_DB})",
    )
    parser.add_argument(
        "--preview-output",
        type=Path,
        default=PREVIEW_PATH,
        help=f"dry-run preview (must resolve to {PREVIEW_PATH})",
    )
    return parser.parse_args()


def require_exact_path(actual: Path, expected: Path, label: str) -> Path:
    resolved = actual.expanduser().resolve()
    if resolved != expected.resolve():
        raise ValidationError(f"{label} must be {expected.resolve()}, not {resolved}")
    return resolved


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_fingerprint(path: Path) -> tuple[int, int, str]:
    stat_result = path.stat()
    return stat_result.st_size, stat_result.st_mtime_ns, sha256_file(path)


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


def index_unique(rows: list[dict[str, str]], label: str) -> dict[str, dict[str, str]]:
    indexed: dict[str, dict[str, str]] = {}
    duplicates: list[str] = []
    for row in rows:
        accession = row["accession"].strip()
        if not accession:
            raise ValidationError(f"blank accession in {label}")
        if accession in indexed:
            duplicates.append(accession)
        indexed[accession] = row
    if duplicates:
        raise ValidationError(
            f"duplicate accessions in {label}: {', '.join(sorted(set(duplicates))[:10])}"
        )
    return indexed


def truth_status_column(headers: list[str]) -> str:
    if "current_structural_status" in headers:
        return "current_structural_status"
    if "structural_status" in headers:
        return "structural_status"
    raise ValidationError(
        f"{STRUCTURAL_TRUTH_PATH} lacks current_structural_status/structural_status"
    )


def load_final_records() -> list[dict[str, str]]:
    adjudication_headers, adjudication_rows = read_tsv(FINAL_ADJUDICATION_PATH)
    require_columns(
        FINAL_ADJUDICATION_PATH,
        adjudication_headers,
        (
            "accession",
            "gene",
            "PWM_status",
            "final_structural_status",
            "DBD_status",
            "expected_DBD",
            "DBD_evidence_coordinates",
        ),
    )
    adjudication = index_unique(adjudication_rows, str(FINAL_ADJUDICATION_PATH))
    if len(adjudication) != 109:
        raise ValidationError(f"expected 109 final adjudication rows, found {len(adjudication)}")
    unresolved_ids = {
        accession
        for accession, row in adjudication.items()
        if row["final_structural_status"].strip() == SOURCE_UNRESOLVED_STATUS
    }
    if len(unresolved_ids) != EXPECTED_UNRESOLVED:
        raise ValidationError(
            f"expected {EXPECTED_UNRESOLVED} final unresolved accessions, "
            f"found {len(unresolved_ids)}"
        )
    rescued = {
        accession
        for accession, row in adjudication.items()
        if row["final_structural_status"].strip() == "VALID_NEW_AF3_DBD_MODEL"
    }
    if rescued != {NFIB_ACCESSION} or NFIB_ACCESSION in unresolved_ids:
        raise ValidationError("NFIB rescue/exclusion invariant failed")

    review_headers, review_rows = read_tsv(UNRESOLVED_REVIEW_PATH)
    require_columns(
        UNRESOLVED_REVIEW_PATH,
        review_headers,
        (
            "accession",
            "reason_currently_unresolved",
            "expected_DBD",
            "DBD_coverage_percent",
            "DBD_identity_percent",
            "accession_coordinates_matching_expected_DBD",
            "pfam_review_category",
        ),
    )
    review = index_unique(review_rows, str(UNRESOLVED_REVIEW_PATH))
    if set(review) != set(adjudication):
        raise ValidationError("final adjudication and unresolved-review accession sets differ")

    truth_headers, truth_rows = read_tsv(STRUCTURAL_TRUTH_PATH)
    require_columns(STRUCTURAL_TRUTH_PATH, truth_headers, ("accession",))
    truth = index_unique(truth_rows, str(STRUCTURAL_TRUTH_PATH))
    if len(truth) != 3690:
        raise ValidationError(f"expected 3690 structural-truth rows, found {len(truth)}")
    status_column = truth_status_column(truth_headers)
    valid_modcre = {
        accession
        for accession, row in truth.items()
        if row[status_column].strip() == "VALID_MODCRE"
    }
    valid_af3 = {
        accession
        for accession, row in truth.items()
        if row[status_column].strip() == "VALID_AF3_FULL_LENGTH"
    }
    truth_unresolved_after_rescue = set(truth) - valid_modcre - valid_af3 - rescued
    if truth_unresolved_after_rescue != unresolved_ids:
        raise ValidationError(
            "final adjudication unresolved set does not equal structural truth after NFIB rescue"
        )

    pwm_counts: dict[str, int] = {}
    own_counts: dict[str, int] = {}
    records: list[dict[str, str]] = []
    for accession in sorted(unresolved_ids):
        adjudicated = adjudication[accession]
        reviewed = review[accession]
        pwm_status = adjudicated["PWM_status"].strip()
        pwm_counts[pwm_status] = pwm_counts.get(pwm_status, 0) + 1
        dbd_status = adjudicated["DBD_status"].strip()
        own_status = OWN_STATUS_BY_DBD_STATUS.get(dbd_status)
        if own_status is None:
            raise ValidationError(f"unsupported DBD_status for {accession}: {dbd_status}")
        own_counts[own_status] = own_counts.get(own_status, 0) + 1
        reason = reviewed["reason_currently_unresolved"].strip()
        if not reason:
            raise ValidationError(f"missing reason_currently_unresolved for {accession}")
        records.append(
            {
                "accession": accession,
                "gene": adjudicated["gene"].strip(),
                "PWM_status": pwm_status,
                "source_final_structural_status": adjudicated[
                    "final_structural_status"
                ].strip(),
                "source_DBD_status": dbd_status,
                "expected_DBD": adjudicated["expected_DBD"].strip()
                or reviewed["expected_DBD"].strip(),
                "DBD_evidence_coordinates": adjudicated[
                    "DBD_evidence_coordinates"
                ].strip(),
                "DBD_coverage_percent": reviewed["DBD_coverage_percent"].strip(),
                "DBD_identity_percent": reviewed["DBD_identity_percent"].strip(),
                "accession_coordinates_matching_expected_DBD": reviewed[
                    "accession_coordinates_matching_expected_DBD"
                ].strip(),
                "pfam_review_category": reviewed["pfam_review_category"].strip(),
                "own_accession_structure_status": own_status,
                "decision_reason": reason,
            }
        )
    if pwm_counts != {
        "Identical_PWM": EXPECTED_IDENTICAL_PWM,
        "Homologous_PWM": EXPECTED_HOMOLOGOUS_PWM,
    }:
        raise ValidationError(f"unexpected PWM breakdown: {pwm_counts}")
    if own_counts != EXPECTED_OWN_STATUS_COUNTS:
        raise ValidationError(f"unexpected detailed status breakdown: {own_counts}")
    return records


def open_readonly_database(path: Path) -> sqlite3.Connection:
    connection = sqlite3.connect(f"file:{path}?mode=ro&immutable=1", uri=True)
    connection.row_factory = sqlite3.Row
    connection.execute("PRAGMA query_only = ON")
    return connection


def copy_database_to_memory(source_path: Path) -> sqlite3.Connection:
    source = open_readonly_database(source_path)
    memory = sqlite3.connect(":memory:")
    try:
        source.backup(memory)
    finally:
        source.close()
    memory.row_factory = sqlite3.Row
    memory.execute("PRAGMA foreign_keys = ON")
    memory.execute("PRAGMA temp_store = MEMORY")
    return memory


def table_columns(connection: sqlite3.Connection, table: str) -> list[str]:
    return [str(row["name"]) for row in connection.execute(f"PRAGMA table_info({table})")]


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
        raise ValidationError(f"required protected table is missing: {table}")
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


def foreign_key_snapshot(connection: sqlite3.Connection) -> tuple[tuple[Any, ...], ...]:
    """Capture inherited FK findings so this scoped import can prove it adds none."""
    return tuple(sorted(tuple(row) for row in connection.execute("PRAGMA foreign_key_check")))


def migration_up_sql() -> str:
    text = MIGRATION_PATH.read_text(encoding="utf-8")
    if "-- migrate:up" not in text or "-- migrate:down" not in text:
        raise ValidationError("migration 002 lacks migrate:up/migrate:down markers")
    return text.split("-- migrate:up", 1)[1].split("-- migrate:down", 1)[0]


def status_schema_has_final_value(connection: sqlite3.Connection) -> bool:
    row = connection.execute(
        "SELECT sql FROM sqlite_master WHERE type='table' AND name='tf_structure_status'"
    ).fetchone()
    return row is not None and PUBLIC_STRUCTURAL_STATUS in str(row["sql"])


def validate_status_schema(connection: sqlite3.Connection) -> None:
    columns = table_columns(connection, "tf_structure_status")
    if tuple(columns) != STATUS_COLUMNS:
        raise ValidationError(f"tf_structure_status columns changed: {columns}")
    if not status_schema_has_final_value(connection):
        raise ValidationError("final public status is absent from the CHECK constraint")
    index_names = {
        str(row["name"])
        for row in connection.execute("PRAGMA index_list(tf_structure_status)")
    }
    required_indexes = {
        "idx_tf_structure_status_primary",
        "idx_tf_structure_status_review",
        "idx_tf_structure_status_canonical",
    }
    if not required_indexes <= index_names:
        raise ValidationError("migration did not preserve the three explicit status indexes")
    unique_index_columns = {
        tuple(str(item["name"]) for item in connection.execute(f"PRAGMA index_info({row['name']})"))
        for row in connection.execute("PRAGMA index_list(tf_structure_status)")
        if int(row["unique"]) == 1
    }
    if ("uniprot_accession",) not in unique_index_columns:
        raise ValidationError("migration did not preserve unique uniprot_accession")
    foreign_keys = {
        (str(row["from"]), str(row["table"]), str(row["to"]), str(row["on_delete"]))
        for row in connection.execute("PRAGMA foreign_key_list(tf_structure_status)")
    }
    expected_foreign_keys = {
        ("tf_id", "tf", "tf_id", "CASCADE"),
        ("canonical_reference_accession", "tf", "tf_id", "RESTRICT"),
    }
    if foreign_keys != expected_foreign_keys:
        raise ValidationError(f"migration changed foreign keys: {foreign_keys}")


def apply_migration_if_needed(connection: sqlite3.Connection) -> bool:
    if status_schema_has_final_value(connection):
        validate_status_schema(connection)
        return False
    before_count = connection.execute("SELECT COUNT(*) FROM tf_structure_status").fetchone()[0]
    connection.executescript(migration_up_sql())
    after_count = connection.execute("SELECT COUNT(*) FROM tf_structure_status").fetchone()[0]
    if before_count != after_count:
        raise ValidationError(
            f"migration changed status row count: {before_count} -> {after_count}"
        )
    validate_status_schema(connection)
    return True


def status_rows(connection: sqlite3.Connection) -> dict[str, dict[str, Any]]:
    return {
        str(row["tf_id"]): dict(row)
        for row in connection.execute("SELECT * FROM tf_structure_status")
    }


def validate_initial_state(
    connection: sqlite3.Connection, final_ids: set[str]
) -> dict[str, Any]:
    rows = status_rows(connection)
    current_ids = set(rows)
    if current_ids == final_ids:
        state = "ALREADY_FINAL"
    elif (
        len(rows) == 102
        and len(current_ids & final_ids) == EXPECTED_EXISTING_FINAL
        and len(final_ids - current_ids) == EXPECTED_MISSING_FINAL
        and current_ids - final_ids == STALE_ACCESSIONS
    ):
        state = "PRE_FINAL_IMPORT"
    else:
        raise ValidationError(
            "unexpected tf_structure_status population: "
            f"rows={len(rows)}, represented={len(current_ids & final_ids)}, "
            f"missing={len(final_ids - current_ids)}, extras={sorted(current_ids - final_ids)}"
        )
    return {
        "state": state,
        "rows": rows,
        "represented": len(current_ids & final_ids),
        "missing": len(final_ids - current_ids),
        "extras": current_ids - final_ids,
    }


def now_utc() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def action_for_new_row(own_status: str) -> str:
    if own_status == "UNRESOLVED":
        return "MANUAL_REVIEW_BEFORE_SENDING"
    return "NO_SEND_MARK_DNA_STRUCTURE_UNAVAILABLE"


def desired_row(
    record: dict[str, str], existing: dict[str, Any] | None, timestamp: str
) -> dict[str, Any]:
    accession = record["accession"]
    if existing is None:
        canonical_status = "NOT_APPLICABLE"
        canonical_accession = None
        canonical_model = None
        action = action_for_new_row(record["own_accession_structure_status"])
        remaining_uncertainty = ""
        created_at = timestamp
    else:
        canonical_status = existing["canonical_reference_status"]
        canonical_accession = existing["canonical_reference_accession"]
        canonical_model = existing["canonical_reference_model_path"]
        action = existing["action_for_baldo"]
        remaining_uncertainty = existing["remaining_uncertainty"]
        created_at = existing["created_at"]
    return {
        "tf_id": accession,
        "uniprot_accession": accession,
        "primary_structural_status": PUBLIC_STRUCTURAL_STATUS,
        "own_accession_structure_status": record["own_accession_structure_status"],
        "canonical_reference_status": canonical_status,
        "canonical_reference_accession": canonical_accession,
        "canonical_reference_model_path": canonical_model,
        "database_display_recommendation": DISPLAY_RECOMMENDATION,
        "action_for_baldo": action,
        "decision_reason": record["decision_reason"],
        "remaining_uncertainty": remaining_uncertainty,
        "review_status": FINAL_REVIEW_STATUS,
        "source_audit_file": SOURCE_PROVENANCE,
        "created_at": created_at,
        "updated_at": timestamp,
    }


def substantive_values(row: dict[str, Any]) -> tuple[Any, ...]:
    return tuple(row[column] for column in STATUS_COLUMNS if column != "updated_at")


def import_status_rows(
    connection: sqlite3.Connection,
    records: list[dict[str, str]],
    timestamp: str,
) -> tuple[dict[str, int], dict[str, dict[str, Any]]]:
    final_ids = {record["accession"] for record in records}
    existing_rows = status_rows(connection)
    planned_rows: dict[str, dict[str, Any]] = {}
    counts = {"updated": 0, "inserted": 0, "removed": 0}
    connection.execute("BEGIN IMMEDIATE")
    try:
        for stale in sorted(STALE_ACCESSIONS):
            if stale in existing_rows:
                connection.execute("DELETE FROM tf_structure_status WHERE tf_id = ?", (stale,))
                counts["removed"] += 1
        for record in records:
            accession = record["accession"]
            existing = existing_rows.get(accession)
            desired = desired_row(record, existing, timestamp)
            planned_rows[accession] = desired
            if existing is None:
                placeholders = ", ".join("?" for _ in STATUS_COLUMNS)
                connection.execute(
                    f"INSERT INTO tf_structure_status ({', '.join(STATUS_COLUMNS)}) "
                    f"VALUES ({placeholders})",
                    tuple(desired[column] for column in STATUS_COLUMNS),
                )
                counts["inserted"] += 1
            elif substantive_values(existing) != substantive_values(desired):
                assignments = ", ".join(
                    f"{column} = ?" for column in STATUS_COLUMNS if column != "tf_id"
                )
                connection.execute(
                    f"UPDATE tf_structure_status SET {assignments} WHERE tf_id = ?",
                    tuple(desired[column] for column in STATUS_COLUMNS if column != "tf_id")
                    + (accession,),
                )
                counts["updated"] += 1
        extras = set(status_rows(connection)) - final_ids
        if extras:
            raise ValidationError(f"unexpected status rows remain after import: {sorted(extras)}")
        connection.commit()
    except Exception:
        connection.rollback()
        raise
    return counts, planned_rows


def validate_final_state(
    connection: sqlite3.Connection,
    records: list[dict[str, str]],
    expected_foreign_key_findings: tuple[tuple[Any, ...], ...],
) -> dict[str, int]:
    final_ids = {record["accession"] for record in records}
    rows = status_rows(connection)
    if set(rows) != final_ids or len(rows) != EXPECTED_UNRESOLVED:
        raise ValidationError(
            f"final status population is not exactly 108: rows={len(rows)}, "
            f"missing={len(final_ids - set(rows))}, extras={len(set(rows) - final_ids)}"
        )
    tf_ids = {str(row["tf_id"]) for row in connection.execute("SELECT tf_id FROM tf")}
    missing_tf = sorted(final_ids - tf_ids)
    if missing_tf:
        raise ValidationError(f"final accessions lack tf.tf_id mappings: {missing_tf}")
    if NFIB_ACCESSION in rows:
        raise ValidationError("NFIB was incorrectly assigned the unresolved status")
    own_counts: dict[str, int] = {}
    by_accession = {record["accession"]: record for record in records}
    for accession, row in rows.items():
        if row["uniprot_accession"] != accession:
            raise ValidationError(f"tf_id/accession mismatch for {accession}")
        if row["primary_structural_status"] != PUBLIC_STRUCTURAL_STATUS:
            raise ValidationError(f"incorrect public status for {accession}")
        if row["database_display_recommendation"] != DISPLAY_RECOMMENDATION:
            raise ValidationError(f"incorrect display recommendation for {accession}")
        if row["decision_reason"] != by_accession[accession]["decision_reason"]:
            raise ValidationError(f"decision reason was not preserved for {accession}")
        if row["review_status"] != FINAL_REVIEW_STATUS:
            raise ValidationError(f"incorrect review status for {accession}")
        if row["source_audit_file"] != SOURCE_PROVENANCE:
            raise ValidationError(f"incorrect source provenance for {accession}")
        own_status = str(row["own_accession_structure_status"])
        own_counts[own_status] = own_counts.get(own_status, 0) + 1
    if own_counts != EXPECTED_OWN_STATUS_COUNTS:
        raise ValidationError(f"final detailed status breakdown differs: {own_counts}")
    foreign_key_findings = foreign_key_snapshot(connection)
    if foreign_key_findings != expected_foreign_key_findings:
        raise ValidationError(
            "the import changed the staging database foreign-key findings: "
            f"before={len(expected_foreign_key_findings)}, after={len(foreign_key_findings)}"
        )
    if connection.execute("PRAGMA quick_check").fetchone()[0] != "ok":
        raise ValidationError("in-memory database quick_check failed")
    return own_counts


def render_preview(
    records: list[dict[str, str]],
    final_rows: dict[str, dict[str, Any]],
    initial_actions: dict[str, str],
) -> str:
    output = io.StringIO(newline="")
    writer = csv.DictWriter(
        output,
        fieldnames=PREVIEW_COLUMNS,
        delimiter="\t",
        lineterminator="\n",
        extrasaction="ignore",
    )
    writer.writeheader()
    for source in records:
        accession = source["accession"]
        status = final_rows[accession]
        writer.writerow(
            {
                **source,
                **status,
                "initial_import_action": initial_actions[accession],
            }
        )
    return output.getvalue()


def write_preview_without_overwrite(path: Path, content: str) -> str:
    encoded = content.encode("utf-8")
    if path.exists():
        if path.is_file() and path.read_bytes() == encoded:
            return "EXISTING_IDENTICAL_NO_CHANGE"
        raise ValidationError(f"refusing to overwrite non-identical preview: {path}")
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(encoded)
    except Exception:
        path.unlink(missing_ok=True)
        raise
    return "CREATED"


def dry_run(database_path: Path, preview_path: Path, records: list[dict[str, str]]) -> dict[str, Any]:
    database_before = file_fingerprint(database_path)
    connection = copy_database_to_memory(database_path)
    try:
        protected_before = protected_snapshots(connection)
        foreign_keys_before = foreign_key_snapshot(connection)
        initial = validate_initial_state(
            connection, {record["accession"] for record in records}
        )
        if initial["state"] != "PRE_FINAL_IMPORT":
            raise ValidationError(
                f"this dry-run expected the pre-final 102-row state, found {initial['state']}"
            )
        apply_migration_if_needed(connection)
        timestamp = now_utc()
        first_counts, _ = import_status_rows(connection, records, timestamp)
        own_counts = validate_final_state(connection, records, foreign_keys_before)
        protected_after_first = protected_snapshots(connection)
        if protected_after_first != protected_before:
            raise ValidationError("a protected table changed during first simulated import")
        final_rows = status_rows(connection)
        initial_actions = {
            record["accession"]: (
                "UPDATE" if record["accession"] in initial["rows"] else "INSERT"
            )
            for record in records
        }
        second_counts, _ = import_status_rows(connection, records, timestamp)
        validate_final_state(connection, records, foreign_keys_before)
        if sum(second_counts.values()) != 0:
            raise ValidationError(f"second simulated import was not idempotent: {second_counts}")
        if protected_snapshots(connection) != protected_before:
            raise ValidationError("a protected table changed during second simulated import")
        preview_action = write_preview_without_overwrite(
            preview_path,
            render_preview(records, final_rows, initial_actions),
        )
    finally:
        connection.close()
    if file_fingerprint(database_path) != database_before:
        raise ValidationError("on-disk staging database changed during dry-run")
    return {
        "first_counts": first_counts,
        "second_counts": second_counts,
        "own_counts": own_counts,
        "preview_action": preview_action,
        "protected_unchanged": True,
    }


def apply_to_staging(database_path: Path, records: list[dict[str, str]]) -> dict[str, Any]:
    connection = sqlite3.connect(database_path)
    connection.row_factory = sqlite3.Row
    connection.execute("PRAGMA foreign_keys = ON")
    try:
        protected_before = protected_snapshots(connection)
        foreign_keys_before = foreign_key_snapshot(connection)
        validate_initial_state(connection, {record["accession"] for record in records})
        apply_migration_if_needed(connection)
        counts, _ = import_status_rows(connection, records, now_utc())
        own_counts = validate_final_state(connection, records, foreign_keys_before)
        if protected_snapshots(connection) != protected_before:
            raise ValidationError("a protected table changed during import")
        return {"counts": counts, "own_counts": own_counts}
    finally:
        connection.close()


def main() -> int:
    args = parse_args()
    database_path = require_exact_path(args.database, STAGING_DB, "database")
    preview_path = require_exact_path(args.preview_output, PREVIEW_PATH, "preview output")
    if database_path == PRODUCTION_DB.resolve():
        raise ValidationError("production database is never a valid target")
    if not database_path.is_file():
        raise ValidationError(f"staging database is missing: {database_path}")
    records = load_final_records()
    if args.dry_run:
        result = dry_run(database_path, preview_path, records)
        counts = result["first_counts"]
        print(
            "dry-run: "
            f"preview={result['preview_action']}; updated={counts['updated']}; "
            f"inserted={counts['inserted']}; removed={counts['removed']}; "
            f"second_changes={sum(result['second_counts'].values())}; final_rows={len(records)}"
        )
    else:
        result = apply_to_staging(database_path, records)
        counts = result["counts"]
        print(
            "apply: "
            f"updated={counts['updated']}; inserted={counts['inserted']}; "
            f"removed={counts['removed']}; final_rows={len(records)}"
        )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ValidationError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)
