#!/usr/bin/env python3
"""Import the 102-accession structural adjudication DRAFT into local staging."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sqlite3
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
SOURCE_DATABASE = ROOT / "data" / "tf_webdb.sqlite"
STAGING_DATABASE = ROOT / "data" / "tf_webdb_local_staging.sqlite"
AUDIT_INPUT = ROOT / "outputs" / "af3_failed_102_final_adjudication_DRAFT.tsv"
MIGRATION = ROOT / "migrations" / "001_add_structure_status_tables.sql"
QC_OUTPUT = ROOT / "outputs" / "local_staging_structure_status_import_qc.json"
SUMMARY_OUTPUT = ROOT / "outputs" / "local_staging_structure_status_import_summary.tsv"
DEFAULT_BASELINE = Path("/tmp/tf_webdb_local_staging_baseline.json")
SOURCE_AUDIT_LABEL = "outputs/af3_failed_102_final_adjudication_DRAFT.tsv"

NEW_TABLES = {
    "tf_structure_status",
    "structure_model_assignment",
    "structure_confidence_artifact",
}

ALLOWED_PRIMARY_STATUSES = {
    "VALID_DBD_FRAGMENT_MODEL",
    "DBD_PRESENT_FRAGMENT_FAILED",
    "DBD_PRESENT_FRAGMENT_NOT_SENT",
    "ACCESSION_LACKS_DBD",
    "ACCESSION_CONTAINS_INCOMPLETE_DBD",
    "COFACTOR_OR_NON_DNA_BINDING_COMPONENT",
    "UNCERTAIN_MANUAL_REVIEW",
}

EXPECTED_STATUS_COUNTS = {
    "ACCESSION_LACKS_DBD": 73,
    "ACCESSION_CONTAINS_INCOMPLETE_DBD": 5,
    "COFACTOR_OR_NON_DNA_BINDING_COMPONENT": 2,
    "VALID_DBD_FRAGMENT_MODEL": 1,
    "UNCERTAIN_MANUAL_REVIEW": 21,
}

OWN_STATUS_BY_PRIMARY = {
    "VALID_DBD_FRAGMENT_MODEL": "VALID_DBD_FRAGMENT_MODEL",
    "DBD_PRESENT_FRAGMENT_FAILED": "DBD_PRESENT_NO_VALID_FRAGMENT",
    "DBD_PRESENT_FRAGMENT_NOT_SENT": "DBD_PRESENT_NO_VALID_FRAGMENT",
    "ACCESSION_LACKS_DBD": "ACCESSION_LACKS_DBD",
    "ACCESSION_CONTAINS_INCOMPLETE_DBD": "ACCESSION_CONTAINS_INCOMPLETE_DBD",
    "COFACTOR_OR_NON_DNA_BINDING_COMPONENT": "NON_DNA_BINDING_COMPONENT",
    "UNCERTAIN_MANUAL_REVIEW": "UNRESOLVED",
}

REQUIRED_INPUT_COLUMNS = {
    "failed_accession",
    "reviewed_canonical_accession",
    "canonical_reference_model",
    "canonical_reference_interface_status",
    "primary_structural_status",
    "action_for_baldo",
    "database_display_recommendation",
    "decision_reason",
    "remaining_uncertainty",
    "draft_status",
}


class ImportValidationError(RuntimeError):
    """Raised when the staging import cannot proceed safely."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Validate and import the 102-accession structural DRAFT into "
            "data/tf_webdb_local_staging.sqlite only."
        )
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Validate inputs, mappings, schema, and controls without inserting rows or writing outputs.",
    )
    parser.add_argument(
        "--expected-initial-sha256",
        default="",
        help="SHA256 recorded after the byte-for-byte staging copy and before migration.",
    )
    parser.add_argument(
        "--baseline-json",
        type=Path,
        default=DEFAULT_BASELINE,
        help=(
            "Read-only pre-import baseline captured externally from source and staging "
            f"(default: {DEFAULT_BASELINE})."
        ),
    )
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_fixed_paths() -> None:
    if STAGING_DATABASE.absolute() == SOURCE_DATABASE.absolute():
        raise ImportValidationError("Staging and source database paths resolve to the same file.")
    if not STAGING_DATABASE.is_file():
        raise ImportValidationError(f"Staging database is missing: {STAGING_DATABASE}")
    if not AUDIT_INPUT.is_file():
        raise ImportValidationError(f"Audit input is missing: {AUDIT_INPUT}")
    if not MIGRATION.is_file():
        raise ImportValidationError(f"Migration is missing: {MIGRATION}")


def connect_staging() -> sqlite3.Connection:
    connection = sqlite3.connect(f"file:{STAGING_DATABASE}?mode=rw", uri=True)
    connection.row_factory = sqlite3.Row
    connection.execute("PRAGMA foreign_keys = ON")
    if connection.execute("PRAGMA foreign_keys").fetchone()[0] != 1:
        connection.close()
        raise ImportValidationError("Could not enable SQLite foreign-key enforcement.")
    return connection


def application_tables(connection: sqlite3.Connection) -> list[str]:
    return [
        row[0]
        for row in connection.execute(
            """
            SELECT name
            FROM sqlite_master
            WHERE type = 'table' AND name NOT LIKE 'sqlite_%'
            ORDER BY name
            """
        )
    ]


def table_counts(connection: sqlite3.Connection, tables: list[str]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for table in tables:
        escaped = table.replace('"', '""')
        counts[table] = int(connection.execute(f'SELECT COUNT(*) FROM "{escaped}"').fetchone()[0])
    return counts


def table_schemas(connection: sqlite3.Connection) -> dict[str, str]:
    return {
        row["name"]: row["sql"]
        for row in connection.execute(
            "SELECT name, sql FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%' ORDER BY name"
        )
    }


def foreign_key_signature(connection: sqlite3.Connection) -> list[tuple[Any, ...]]:
    return sorted(tuple(row) for row in connection.execute("PRAGMA foreign_key_check").fetchall())


def read_baseline(path: Path) -> dict[str, Any]:
    try:
        with path.open("r", encoding="utf-8") as handle:
            baseline = json.load(handle)
    except (OSError, json.JSONDecodeError) as exc:
        raise ImportValidationError(f"Could not read baseline JSON {path}: {exc}") from exc
    required = {
        "source_database",
        "source_before",
        "source_after_preimport",
        "source_unchanged_preimport",
        "source_tables",
        "source_counts",
        "source_schemas",
        "source_fk_signature",
        "staging_database",
        "staging_initial_sha256",
    }
    missing = sorted(required - set(baseline))
    if missing:
        raise ImportValidationError(f"Baseline JSON is missing fields: {missing}")
    if baseline["source_database"] != str(SOURCE_DATABASE):
        raise ImportValidationError("Baseline source database path does not match the fixed source path.")
    if baseline["staging_database"] != str(STAGING_DATABASE):
        raise ImportValidationError("Baseline staging database path does not match the fixed staging path.")
    if not baseline["source_unchanged_preimport"]:
        raise ImportValidationError("Baseline reports that the source changed before import.")
    if baseline["source_before"] != baseline["source_after_preimport"]:
        raise ImportValidationError("Source before/pre-import fingerprints differ in baseline.")
    if baseline["staging_initial_sha256"] != baseline["source_before"].get("sha256"):
        raise ImportValidationError("Initial staging SHA256 does not equal the source SHA256 in baseline.")
    return baseline


def read_audit_rows() -> list[dict[str, str]]:
    with AUDIT_INPUT.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ImportValidationError("Audit TSV has no header.")
        missing = sorted(REQUIRED_INPUT_COLUMNS - set(reader.fieldnames))
        if missing:
            raise ImportValidationError(f"Audit TSV is missing required columns: {missing}")
        rows = [dict(row) for row in reader]
    if len(rows) != 102:
        raise ImportValidationError(f"Expected 102 audit rows, found {len(rows)}.")
    accessions = [row["failed_accession"].strip() for row in rows]
    if len(set(accessions)) != 102:
        duplicates = sorted(accession for accession, count in Counter(accessions).items() if count > 1)
        raise ImportValidationError(f"Audit accessions are not unique: {duplicates}")
    if any(not accession for accession in accessions):
        raise ImportValidationError("Audit contains a blank failed_accession.")
    return rows


def require_exact_tf_mapping(
    connection: sqlite3.Connection,
    accession: str,
) -> tuple[str, str]:
    matches = connection.execute(
        """
        SELECT tf.tf_id, ta.uniprot_accession
        FROM tf
        JOIN tf_annotation AS ta ON ta.tf_id = tf.tf_id
        WHERE UPPER(tf.tf_id) = UPPER(?)
        """,
        (accession,),
    ).fetchall()
    if len(matches) != 1:
        raise ImportValidationError(
            f"Accession mapping is missing or ambiguous for {accession}: {len(matches)} matches."
        )
    tf_id = str(matches[0]["tf_id"])
    uniprot_accession = str(matches[0]["uniprot_accession"])
    if tf_id != accession or uniprot_accession != accession:
        raise ImportValidationError(
            f"Accession mapping is not an exact identity for {accession}: "
            f"tf_id={tf_id}, uniprot_accession={uniprot_accession}."
        )
    return tf_id, uniprot_accession


def canonical_status(row: dict[str, str]) -> str:
    accession = row["reviewed_canonical_accession"].strip()
    model_path = row["canonical_reference_model"].strip()
    interface_status = row["canonical_reference_interface_status"].strip()
    if not accession:
        if model_path:
            raise ImportValidationError(
                f"{row['failed_accession']} has a canonical model path but no canonical accession."
            )
        return "NOT_APPLICABLE"
    if not model_path:
        return "SEPARATE_REFERENCE_ACCESSION_ONLY"
    if interface_status == "PASS_NONEMPTY_PROTEIN_DNA_INTERFACE":
        return "SEPARATE_REFERENCE_INTERFACE_VALID"
    return "SEPARATE_REFERENCE_NOT_VALIDATED"


def validate_rows(
    staging: sqlite3.Connection,
    rows: list[dict[str, str]],
) -> list[dict[str, Any]]:
    actual_counts = Counter(row["primary_structural_status"].strip() for row in rows)
    if set(actual_counts) - ALLOWED_PRIMARY_STATUSES:
        raise ImportValidationError(
            f"Unexpected primary statuses: {sorted(set(actual_counts) - ALLOWED_PRIMARY_STATUSES)}"
        )
    if dict(actual_counts) != EXPECTED_STATUS_COUNTS:
        raise ImportValidationError(
            f"Status counts do not match the required counts: {dict(actual_counts)}"
        )

    now = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    prepared: list[dict[str, Any]] = []
    for row in rows:
        accession = row["failed_accession"].strip()
        tf_id, uniprot_accession = require_exact_tf_mapping(staging, accession)
        primary_status = row["primary_structural_status"].strip()
        review_status = row["draft_status"].strip()
        if review_status != "FINAL_ADJUDICATION_DRAFT_NOT_PUBLIC":
            raise ImportValidationError(
                f"{accession} has unexpected review/draft status: {review_status}"
            )
        if primary_status == "UNCERTAIN_MANUAL_REVIEW" and not row["remaining_uncertainty"].strip():
            raise ImportValidationError(
                f"{accession} is UNCERTAIN_MANUAL_REVIEW but remaining_uncertainty is blank."
            )
        canonical_accession = row["reviewed_canonical_accession"].strip()
        if canonical_accession:
            require_exact_tf_mapping(staging, canonical_accession)
            if canonical_accession == accession:
                raise ImportValidationError(
                    f"{accession} incorrectly uses itself as a canonical same-gene reference."
                )
        prepared.append(
            {
                "tf_id": tf_id,
                "uniprot_accession": uniprot_accession,
                "primary_structural_status": primary_status,
                "own_accession_structure_status": OWN_STATUS_BY_PRIMARY[primary_status],
                "canonical_reference_status": canonical_status(row),
                "canonical_reference_accession": canonical_accession or None,
                "canonical_reference_model_path": row["canonical_reference_model"].strip() or None,
                "database_display_recommendation": row["database_display_recommendation"],
                "action_for_baldo": row["action_for_baldo"],
                "decision_reason": row["decision_reason"],
                "remaining_uncertainty": row["remaining_uncertainty"],
                "review_status": review_status,
                "source_audit_file": SOURCE_AUDIT_LABEL,
                "created_at": now,
                "updated_at": now,
            }
        )
    return prepared


def require_migration_tables(staging: sqlite3.Connection) -> None:
    tables = set(application_tables(staging))
    missing = sorted(NEW_TABLES - tables)
    if missing:
        raise ImportValidationError(f"Required migration tables are missing: {missing}")
    columns = {
        row[1]
        for row in staging.execute("PRAGMA table_info(tf_structure_status)").fetchall()
    }
    required = {
        "tf_id", "uniprot_accession", "primary_structural_status",
        "own_accession_structure_status", "canonical_reference_status",
        "canonical_reference_accession", "canonical_reference_model_path",
        "database_display_recommendation", "action_for_baldo", "decision_reason",
        "remaining_uncertainty", "review_status", "source_audit_file",
        "created_at", "updated_at",
    }
    if columns != required:
        raise ImportValidationError(
            f"tf_structure_status columns differ from migration contract: {sorted(columns)}"
        )


def validate_control_rows(staging: sqlite3.Connection) -> dict[str, dict[str, Any]]:
    controls: dict[str, dict[str, Any]] = {}
    for accession in ("E9PN75", "H7C4N4", "B4DHE0", "Q59EF3"):
        row = staging.execute(
            "SELECT * FROM tf_structure_status WHERE tf_id = ?", (accession,)
        ).fetchone()
        if row is None:
            raise ImportValidationError(f"Control row is missing after import: {accession}")
        controls[accession] = dict(row)

    e9 = controls["E9PN75"]
    if e9["primary_structural_status"] != "ACCESSION_LACKS_DBD":
        raise ImportValidationError("E9PN75 status control failed.")
    if e9["canonical_reference_accession"] != "Q9NZC4":
        raise ImportValidationError("E9PN75 canonical reference control failed.")
    if e9["own_accession_structure_status"] != "ACCESSION_LACKS_DBD":
        raise ImportValidationError("E9PN75 own-accession structure control failed.")
    if controls["H7C4N4"]["primary_structural_status"] != "VALID_DBD_FRAGMENT_MODEL":
        raise ImportValidationError("H7C4N4 valid-fragment control failed.")
    for accession in ("B4DHE0", "Q59EF3"):
        if controls[accession]["primary_structural_status"] == "VALID_DBD_FRAGMENT_MODEL":
            raise ImportValidationError(f"{accession} non-DBD fragment control failed.")
    assignment_count = int(staging.execute(
        "SELECT COUNT(*) FROM structure_model_assignment"
    ).fetchone()[0])
    if assignment_count != 0:
        raise ImportValidationError("structure_model_assignment must remain empty.")
    return controls


def atomic_write_text(path: Path, content: str) -> None:
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    try:
        with temporary.open("x", encoding="utf-8", newline="") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def main() -> int:
    args = parse_args()
    try:
        validate_fixed_paths()
        baseline = read_baseline(args.baseline_json)
        if args.expected_initial_sha256 and (
            len(args.expected_initial_sha256) != 64
            or any(char not in "0123456789abcdef" for char in args.expected_initial_sha256)
        ):
            raise ImportValidationError("--expected-initial-sha256 is not a lowercase SHA256 value.")
        baseline_initial_sha = str(baseline["staging_initial_sha256"])
        if args.expected_initial_sha256 and args.expected_initial_sha256 != baseline_initial_sha:
            raise ImportValidationError(
                "--expected-initial-sha256 differs from the externally captured baseline."
            )
        if not args.dry_run:
            existing_outputs = [path for path in (QC_OUTPUT, SUMMARY_OUTPUT) if path.exists()]
            if existing_outputs:
                raise ImportValidationError(
                    f"Refusing to overwrite output files: {[str(path) for path in existing_outputs]}"
                )

        audit_rows = read_audit_rows()
        source_tables = [str(table) for table in baseline["source_tables"]]
        source_counts = {
            str(table): int(count) for table, count in baseline["source_counts"].items()
        }
        source_schemas = {
            str(table): str(schema) for table, schema in baseline["source_schemas"].items()
        }
        source_fk_signature = sorted(
            tuple(row) for row in baseline["source_fk_signature"]
        )

        with connect_staging() as staging:
            require_migration_tables(staging)
            staging_tables = application_tables(staging)
            expected_tables = set(source_tables) | NEW_TABLES
            if set(staging_tables) != expected_tables:
                raise ImportValidationError(
                    "Staging application tables differ from source plus the three migration tables: "
                    f"extra={sorted(set(staging_tables)-expected_tables)}, "
                    f"missing={sorted(expected_tables-set(staging_tables))}"
                )
            staging_schemas = table_schemas(staging)
            for table, source_sql in source_schemas.items():
                if staging_schemas.get(table) != source_sql:
                    raise ImportValidationError(f"Existing table schema changed: {table}")
            for table in NEW_TABLES:
                existing = int(staging.execute(f'SELECT COUNT(*) FROM "{table}"').fetchone()[0])
                if existing != 0:
                    raise ImportValidationError(
                        f"Refusing to import because {table} already has {existing} rows."
                    )
            prepared = validate_rows(staging, audit_rows)

            if args.dry_run:
                print("DRY_RUN_VALIDATED")
                print(f"rows={len(prepared)}")
                print(f"status_counts={dict(Counter(row['primary_structural_status'] for row in prepared))}")
                return 0

            insert_sql = """
                INSERT INTO tf_structure_status (
                    tf_id, uniprot_accession, primary_structural_status,
                    own_accession_structure_status, canonical_reference_status,
                    canonical_reference_accession, canonical_reference_model_path,
                    database_display_recommendation, action_for_baldo,
                    decision_reason, remaining_uncertainty, review_status,
                    source_audit_file, created_at, updated_at
                ) VALUES (
                    :tf_id, :uniprot_accession, :primary_structural_status,
                    :own_accession_structure_status, :canonical_reference_status,
                    :canonical_reference_accession, :canonical_reference_model_path,
                    :database_display_recommendation, :action_for_baldo,
                    :decision_reason, :remaining_uncertainty, :review_status,
                    :source_audit_file, :created_at, :updated_at
                )
            """
            staging.execute("BEGIN IMMEDIATE")
            try:
                staging.executemany(insert_sql, prepared)
                imported_count = int(staging.execute(
                    "SELECT COUNT(*) FROM tf_structure_status"
                ).fetchone()[0])
                distinct_count = int(staging.execute(
                    "SELECT COUNT(DISTINCT tf_id) FROM tf_structure_status"
                ).fetchone()[0])
                if imported_count != 102 or distinct_count != 102:
                    raise ImportValidationError(
                        f"Imported row reconciliation failed: rows={imported_count}, distinct={distinct_count}."
                    )
                imported_status_counts = {
                    row[0]: int(row[1])
                    for row in staging.execute(
                        "SELECT primary_structural_status, COUNT(*) FROM tf_structure_status GROUP BY primary_structural_status"
                    )
                }
                if imported_status_counts != EXPECTED_STATUS_COUNTS:
                    raise ImportValidationError(
                        f"Imported status counts differ: {imported_status_counts}"
                    )
                unresolved_count = int(staging.execute(
                    """
                    SELECT COUNT(*) FROM tf_structure_status
                    WHERE primary_structural_status='UNCERTAIN_MANUAL_REVIEW'
                      AND own_accession_structure_status='UNRESOLVED'
                      AND review_status='FINAL_ADJUDICATION_DRAFT_NOT_PUBLIC'
                    """
                ).fetchone()[0])
                if unresolved_count != 21:
                    raise ImportValidationError(
                        f"Unresolved manual-review preservation failed: {unresolved_count}."
                    )
                controls = validate_control_rows(staging)
                confidence_count = int(staging.execute(
                    "SELECT COUNT(*) FROM structure_confidence_artifact"
                ).fetchone()[0])
                if confidence_count != 0:
                    raise ImportValidationError("structure_confidence_artifact must remain empty.")
                staging.commit()
            except Exception:
                staging.rollback()
                raise

        with connect_staging() as staging:
            quick_check = str(staging.execute("PRAGMA quick_check").fetchone()[0])
            integrity_check = str(staging.execute("PRAGMA integrity_check").fetchone()[0])
            if quick_check != "ok" or integrity_check != "ok":
                raise ImportValidationError(
                    f"Post-import integrity failed: quick={quick_check}, integrity={integrity_check}"
                )
            staging_fk_signature = foreign_key_signature(staging)
            new_table_fk_violations = [
                row for row in staging_fk_signature if row[0] in NEW_TABLES
            ]
            if new_table_fk_violations:
                raise ImportValidationError(
                    f"New tables have foreign-key violations: {new_table_fk_violations[:10]}"
                )
            if staging_fk_signature != source_fk_signature:
                raise ImportValidationError(
                    "Foreign-key violation baseline changed; existing violations must be preserved exactly."
                )
            staging_existing_counts = table_counts(staging, source_tables)
            if staging_existing_counts != source_counts:
                changed = {
                    table: {"source": source_counts[table], "staging": staging_existing_counts[table]}
                    for table in source_tables
                    if source_counts[table] != staging_existing_counts[table]
                }
                raise ImportValidationError(f"Existing table row counts changed: {changed}")
            imported_count = int(staging.execute(
                "SELECT COUNT(*) FROM tf_structure_status"
            ).fetchone()[0])
            assignment_count = int(staging.execute(
                "SELECT COUNT(*) FROM structure_model_assignment"
            ).fetchone()[0])
            confidence_count = int(staging.execute(
                "SELECT COUNT(*) FROM structure_confidence_artifact"
            ).fetchone()[0])
            imported_status_counts = {
                row[0]: int(row[1])
                for row in staging.execute(
                    "SELECT primary_structural_status, COUNT(*) FROM tf_structure_status GROUP BY primary_structural_status"
                )
            }
            duplicate_tf_count = int(staging.execute(
                """
                SELECT COUNT(*) FROM (
                    SELECT tf_id FROM tf_structure_status GROUP BY tf_id HAVING COUNT(*) > 1
                )
                """
            ).fetchone()[0])
            duplicate_uniprot_count = int(staging.execute(
                """
                SELECT COUNT(*) FROM (
                    SELECT uniprot_accession FROM tf_structure_status
                    GROUP BY uniprot_accession HAVING COUNT(*) > 1
                )
                """
            ).fetchone()[0])
            controls = validate_control_rows(staging)

        final_staging_sha = sha256_file(STAGING_DATABASE)
        summary_rows: list[dict[str, Any]] = []
        for status in sorted(EXPECTED_STATUS_COUNTS):
            actual = imported_status_counts.get(status, 0)
            expected = EXPECTED_STATUS_COUNTS[status]
            summary_rows.append(
                {
                    "category": "PRIMARY_STATUS_COUNT",
                    "key": status,
                    "value": actual,
                    "expected": expected,
                    "result": "PASS" if actual == expected else "FAIL",
                    "notes": "DRAFT import count",
                }
            )
        summary_rows.extend(
            [
                {"category":"CONTROL","key":"E9PN75","value":controls["E9PN75"]["primary_structural_status"],"expected":"ACCESSION_LACKS_DBD; canonical=Q9NZC4; no assignment","result":"PASS","notes":"Q9NZC4 is canonical reference evidence only."},
                {"category":"CONTROL","key":"H7C4N4","value":controls["H7C4N4"]["primary_structural_status"],"expected":"VALID_DBD_FRAGMENT_MODEL","result":"PASS","notes":"No structure assignment inserted."},
                {"category":"CONTROL","key":"B4DHE0","value":controls["B4DHE0"]["primary_structural_status"],"expected":"not VALID_DBD_FRAGMENT_MODEL","result":"PASS","notes":"Contacting non-DBD fragment not promoted."},
                {"category":"CONTROL","key":"Q59EF3","value":controls["Q59EF3"]["primary_structural_status"],"expected":"not VALID_DBD_FRAGMENT_MODEL","result":"PASS","notes":"Contacting non-DBD fragment not promoted."},
                {"category":"TABLE_COUNT","key":"tf_structure_status","value":imported_count,"expected":102,"result":"PASS","notes":"One row per failed accession."},
                {"category":"TABLE_COUNT","key":"structure_model_assignment","value":assignment_count,"expected":0,"result":"PASS","notes":"Deliberately not populated."},
                {"category":"TABLE_COUNT","key":"structure_confidence_artifact","value":confidence_count,"expected":0,"result":"PASS","notes":"Deliberately not populated."},
            ]
        )
        fieldnames = ["category", "key", "value", "expected", "result", "notes"]
        from io import StringIO
        stream = StringIO()
        writer = csv.DictWriter(stream, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(summary_rows)
        summary_text = stream.getvalue()

        qc = {
            "audit_input": str(AUDIT_INPUT),
            "audit_input_row_count": len(audit_rows),
            "audit_input_unique_accession_count": len({row["failed_accession"] for row in audit_rows}),
            "controls": {
                "E9PN75": {
                    "result": "PASS",
                    "primary_structural_status": controls["E9PN75"]["primary_structural_status"],
                    "own_accession_structure_status": controls["E9PN75"]["own_accession_structure_status"],
                    "canonical_reference_accession": controls["E9PN75"]["canonical_reference_accession"],
                    "own_model_assignment_count": 0,
                },
                "H7C4N4": {"result":"PASS","primary_structural_status":controls["H7C4N4"]["primary_structural_status"]},
                "B4DHE0": {"result":"PASS","primary_structural_status":controls["B4DHE0"]["primary_structural_status"]},
                "Q59EF3": {"result":"PASS","primary_structural_status":controls["Q59EF3"]["primary_structural_status"]},
            },
            "created_files": [
                str(STAGING_DATABASE), str(MIGRATION), str(Path(__file__).resolve()),
                str(QC_OUTPUT), str(SUMMARY_OUTPUT),
            ],
            "duplicate_tf_id_count": duplicate_tf_count,
            "duplicate_uniprot_accession_count": duplicate_uniprot_count,
            "existing_table_row_counts_changed": False,
            "existing_table_row_counts_source": source_counts,
            "existing_table_row_counts_staging": staging_existing_counts,
            "existing_table_schema_changed": False,
            "foreign_key_check": {
                "source_preexisting_violation_count": len(source_fk_signature),
                "staging_postimport_violation_count": len(staging_fk_signature),
                "signatures_identical": staging_fk_signature == source_fk_signature,
                "new_table_violation_count": len(new_table_fk_violations),
                "result": "PASS_BASELINE_PRESERVED_NO_NEW_VIOLATIONS",
            },
            "imported_status_counts": imported_status_counts,
            "integrity_check": integrity_check,
            "migration_file": str(MIGRATION),
            "migration_tables_created": sorted(NEW_TABLES),
            "network_connection_attempted": False,
            "pwm_annotation_tables_modified": False,
            "quick_check": quick_check,
            "review_status_counts": {"FINAL_ADJUDICATION_DRAFT_NOT_PUBLIC": 102},
            "source_database": str(SOURCE_DATABASE),
            "source_database_after_preimport": baseline["source_after_preimport"],
            "source_database_before": baseline["source_before"],
            "source_database_unchanged_preimport": baseline["source_unchanged_preimport"],
            "source_database_opened_by_importer": False,
            "staging_database": str(STAGING_DATABASE),
            "staging_final_sha256": final_staging_sha,
            "staging_initial_sha256": baseline_initial_sha,
            "structure_confidence_artifact_count": confidence_count,
            "structure_model_assignment_count": assignment_count,
            "tf_structure_status_count": imported_count,
            "uncertain_manual_review_preserved_as_unresolved": 21,
            "validation_status": "PASS",
            "website_code_modified": False,
        }
        atomic_write_text(SUMMARY_OUTPUT, summary_text)
        atomic_write_text(QC_OUTPUT, json.dumps(qc, indent=2, sort_keys=True) + "\n")

        print("IMPORTED_AND_VALIDATED")
        print(f"tf_structure_status={imported_count}")
        print(f"status_counts={imported_status_counts}")
        print(f"source_sha256_from_external_baseline={baseline['source_before']['sha256']}")
        print(f"staging_sha256={final_staging_sha}")
        return 0
    except (ImportValidationError, sqlite3.Error, OSError, csv.Error) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
