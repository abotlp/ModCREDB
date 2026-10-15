#!/usr/bin/env python3
"""Remove only proven stale structural child rows from local staging.

The complete manuscript restoration intentionally did not recreate structural
evidence absent from the frozen integrated TSV.  Older child rows survived for
deleted parent IDs, producing duplicate or parentless relationships.  This
script reconstructs every stale row against the historical pre-removal
database and the frozen integrated TSV, refuses unknown cases, and deletes only
the broken child rows.  It never deletes motif, structure, PWM, or model files.

Dry-run operates on an in-memory SQLite copy.  Apply mode creates and verifies
a timestamped byte-for-byte backup before opening the target database writable.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import sqlite3
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATABASE = ROOT / "data/tf_webdb_local_staging.sqlite"
DEFAULT_HISTORICAL_DATABASE = ROOT / "data/tf_webdb.before_remove_TARDBP.sqlite"
DEFAULT_INTEGRATED_TSV = Path(
    "/home/patricia/TF_database_Baldo_data/"
    "TF_PWM_chart_final_integrated_HOCOMOCO_hierarchical.tsv"
)
DEFAULT_PREVIEW = ROOT / "outputs/stale_structure_foreign_key_cleanup_preview.tsv"
DEFAULT_QC = ROOT / "outputs/stale_structure_foreign_key_cleanup_qc.json"
DEFAULT_BACKUPS = ROOT / "backups"

INTEGRATED_SHA256 = "f71160a2c2d318512be8354302208b55b4d4e7b427178631d98977c4209ed845"
HISTORICAL_SHA256 = "18f0018d1f9babf884025248ec3fd5002a4f238fe9e36bafd7022acc8ca6138e"

EXPECTED_MOTIF_STRUCTURE_STALE = 823
EXPECTED_MODEL_SUMMARY_STALE = 6
EXPECTED_LINK_CLASSES = Counter(
    {
        "replacement_already_present": 55,
        "excluded_unannotated_alphafold": 762,
        "excluded_unannotated_modcre": 1,
        "excluded_not_in_delivered_alphafold_tokens": 5,
    }
)
EXPECTED_NONDELIVERED_MODELS = {
    "A0A0A0N0M3": 4,
    "B4DRW3": 1,
}
EXPECTED_SUMMARIES_BY_TF = {"Q16533": 2, "Q16594": 4}
EXPECTED_FINAL_COUNTS = {"motif_structure": 11685, "model_summary": 7694}

PROTECTED_TABLES = (
    "source",
    "source_release",
    "tf",
    "tf_annotation",
    "tf_primary_annotation",
    "tf_family",
    "motif_ref",
    "motif_file",
    "structure_file",
    "structure_model_assignment",
    "structure_confidence_artifact",
    "tf_structure_status",
)

PREVIEW_COLUMNS = (
    "table_name",
    "row_key",
    "tf_id",
    "source",
    "motif_id",
    "model_id",
    "historical_motif_ref_id",
    "historical_structure_file_id",
    "classification",
    "action",
    "reason",
)


class CleanupError(RuntimeError):
    pass


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--dry-run", action="store_true")
    mode.add_argument("--apply", action="store_true")
    parser.add_argument("--database", type=Path, default=DEFAULT_DATABASE)
    parser.add_argument(
        "--historical-database", type=Path, default=DEFAULT_HISTORICAL_DATABASE
    )
    parser.add_argument("--integrated-tsv", type=Path, default=DEFAULT_INTEGRATED_TSV)
    parser.add_argument("--preview-output", type=Path, default=DEFAULT_PREVIEW)
    parser.add_argument("--qc-output", type=Path, default=DEFAULT_QC)
    parser.add_argument("--backup-dir", type=Path, default=DEFAULT_BACKUPS)
    parser.add_argument("--replace-outputs", action="store_true")
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def connect_ro(path: Path) -> sqlite3.Connection:
    connection = sqlite3.connect(f"file:{path.resolve()}?mode=ro&immutable=1", uri=True)
    connection.row_factory = sqlite3.Row
    return connection


def in_memory_copy(path: Path) -> sqlite3.Connection:
    source = connect_ro(path)
    destination = sqlite3.connect(":memory:")
    destination.row_factory = sqlite3.Row
    source.backup(destination)
    source.close()
    destination.execute("PRAGMA foreign_keys=ON")
    return destination


def table_columns(connection: sqlite3.Connection, table: str) -> list[str]:
    return [str(row[1]) for row in connection.execute(f'PRAGMA table_info("{table}")')]


def encoded(value: Any) -> bytes:
    if value is None:
        return b"N"
    payload = value if isinstance(value, bytes) else str(value).encode("utf-8")
    prefix = b"B" if isinstance(value, bytes) else b"T"
    return prefix + str(len(payload)).encode("ascii") + b":" + payload


def table_snapshot(connection: sqlite3.Connection, table: str) -> tuple[int, str]:
    columns = table_columns(connection, table)
    if not columns:
        return (0, "ABSENT")
    quoted = ",".join('"' + column.replace('"', '""') + '"' for column in columns)
    digest = hashlib.sha256()
    count = 0
    for row in connection.execute(f'SELECT {quoted} FROM "{table}" ORDER BY {quoted}'):
        count += 1
        for value in row:
            digest.update(encoded(value))
            digest.update(b"|")
        digest.update(b"\n")
    return count, digest.hexdigest()


def snapshots(connection: sqlite3.Connection) -> dict[str, tuple[int, str]]:
    return {table: table_snapshot(connection, table) for table in PROTECTED_TABLES}


def read_integrated(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    result = {row["TF_name"].strip(): row for row in rows}
    if len(rows) != 5417 or len(result) != 5417:
        raise CleanupError(
            f"Expected 5,417 unique integrated rows; found {len(rows)} rows/{len(result)} IDs"
        )
    return result


def split_tokens(value: str) -> set[str]:
    return {item.strip() for item in value.split(";") if item.strip()}


def current_equivalent(
    connection: sqlite3.Connection,
    motif: sqlite3.Row,
    structure: sqlite3.Row,
) -> tuple[int | None, int | None, bool]:
    current_motif = connection.execute(
        """
        SELECT id FROM motif_ref
        WHERE tf_id = ? AND source = ? AND motif_id = ?
        ORDER BY id
        """,
        (motif["tf_id"], motif["source"], motif["motif_id"]),
    ).fetchall()
    current_structure = connection.execute(
        """
        SELECT id FROM structure_file
        WHERE tf_id = ? AND source = ? AND model_id = ? AND file_type = ?
        ORDER BY id
        """,
        (
            structure["tf_id"],
            structure["source"],
            structure["model_id"],
            structure["file_type"],
        ),
    ).fetchall()
    if len(current_motif) > 1 or len(current_structure) > 1:
        raise CleanupError(
            "Non-unique replacement candidate for "
            f"{motif['tf_id']} {motif['motif_id']} {structure['model_id']}"
        )
    motif_id = int(current_motif[0][0]) if current_motif else None
    structure_id = int(current_structure[0][0]) if current_structure else None
    linked = False
    if motif_id is not None and structure_id is not None:
        linked = connection.execute(
            "SELECT 1 FROM motif_structure WHERE motif_ref_id=? AND structure_file_id=?",
            (motif_id, structure_id),
        ).fetchone() is not None
    return motif_id, structure_id, linked


def stale_motif_structure_rows(connection: sqlite3.Connection) -> list[sqlite3.Row]:
    return connection.execute(
        """
        SELECT ms.rowid AS child_rowid, ms.motif_ref_id, ms.structure_file_id,
               mr.id AS current_motif_exists, sf.id AS current_structure_exists
        FROM motif_structure AS ms
        LEFT JOIN motif_ref AS mr ON mr.id = ms.motif_ref_id
        LEFT JOIN structure_file AS sf ON sf.id = ms.structure_file_id
        WHERE mr.id IS NULL OR sf.id IS NULL
        ORDER BY ms.rowid
        """
    ).fetchall()


def stale_model_summary_rows(connection: sqlite3.Connection) -> list[sqlite3.Row]:
    return connection.execute(
        """
        SELECT summary.*, sf_summary.id AS current_summary_exists,
               sf_match.id AS current_match_exists
        FROM model_summary AS summary
        LEFT JOIN structure_file AS sf_summary ON sf_summary.id = summary.summary_file_id
        LEFT JOIN structure_file AS sf_match ON sf_match.id = summary.matched_structure_id
        WHERE sf_summary.id IS NULL
           OR (summary.matched_structure_id IS NOT NULL AND sf_match.id IS NULL)
        ORDER BY summary.id
        """
    ).fetchall()


def classify_stale_links(
    connection: sqlite3.Connection,
    historical: sqlite3.Connection,
    integrated: dict[str, dict[str, str]],
) -> tuple[list[dict[str, str]], Counter[str]]:
    preview: list[dict[str, str]] = []
    classes: Counter[str] = Counter()
    for child in stale_motif_structure_rows(connection):
        motif = historical.execute(
            "SELECT * FROM motif_ref WHERE id=?", (child["motif_ref_id"],)
        ).fetchone()
        structure = historical.execute(
            "SELECT * FROM structure_file WHERE id=?", (child["structure_file_id"],)
        ).fetchone()
        if motif is None or structure is None:
            raise CleanupError(
                f"Historical parents unavailable for motif_structure rowid {child['child_rowid']}"
            )
        if motif["tf_id"] != structure["tf_id"]:
            raise CleanupError(
                f"Historical TF mismatch at motif_structure rowid {child['child_rowid']}"
            )
        tf_id = str(motif["tf_id"])
        delivered = integrated.get(tf_id)
        if delivered is None:
            raise CleanupError(f"Historical link TF {tf_id} absent from integrated TSV")

        _, _, linked = current_equivalent(connection, motif, structure)
        if linked:
            classification = "replacement_already_present"
            reason = (
                "Equivalent current motif and structure rows are already linked; "
                "this child row contains superseded historical integer IDs."
            )
        elif delivered["Best_annotation_level"] == "Unannotated":
            classification = f"excluded_unannotated_{str(structure['source']).lower()}"
            reason = (
                "Frozen integrated TSV classifies this TF as Unannotated and contains "
                "no delivered PWM/model token; historical parent rows were intentionally "
                "not restored."
            )
        elif (
            delivered["Best_annotation_level"] == "AlphaFold"
            and str(structure["source"]).lower() == "alphafold"
        ):
            delivered_tokens = split_tokens(delivered["AlphaFold"])
            historical_token = str(motif["original_value"])
            if historical_token not in delivered_tokens:
                classification = "excluded_not_in_delivered_alphafold_tokens"
                reason = (
                    "Historical AlphaFold motif/model is not among the exact AlphaFold "
                    "tokens in the frozen integrated TSV."
                )
            else:
                raise CleanupError(
                    f"Delivered AlphaFold link lacks a current replacement: {tf_id} "
                    f"{historical_token}"
                )
        else:
            raise CleanupError(
                "Unclassified stale motif_structure row: "
                f"TF={tf_id} class={delivered['Best_annotation_level']} "
                f"source={structure['source']} motif={motif['motif_id']}"
            )

        classes[classification] += 1
        preview.append(
            {
                "table_name": "motif_structure",
                "row_key": str(child["child_rowid"]),
                "tf_id": tf_id,
                "source": str(structure["source"]),
                "motif_id": str(motif["motif_id"]),
                "model_id": str(structure["model_id"]),
                "historical_motif_ref_id": str(child["motif_ref_id"]),
                "historical_structure_file_id": str(child["structure_file_id"]),
                "classification": classification,
                "action": "DELETE_BROKEN_CHILD_ROW",
                "reason": reason,
            }
        )
    return preview, classes


def classify_stale_summaries(
    connection: sqlite3.Connection,
    historical: sqlite3.Connection,
    integrated: dict[str, dict[str, str]],
) -> tuple[list[dict[str, str]], Counter[str]]:
    preview: list[dict[str, str]] = []
    by_tf: Counter[str] = Counter()
    for summary in stale_model_summary_rows(connection):
        tf_id = str(summary["tf_id"])
        by_tf[tf_id] += 1
        delivered = integrated.get(tf_id)
        if delivered is None or delivered["Best_annotation_level"] != "Unannotated":
            raise CleanupError(f"Stale model_summary TF {tf_id} is not delivered Unannotated")
        tf = connection.execute(
            "SELECT motif_ref_count, active_model_count FROM tf WHERE tf_id=?", (tf_id,)
        ).fetchone()
        if tf is None or tuple(tf) != (0, 0):
            raise CleanupError(
                f"Expected zero motif/model counts for {tf_id}; found {tuple(tf) if tf else None}"
            )
        historical_summary = historical.execute(
            "SELECT * FROM structure_file WHERE id=?", (summary["summary_file_id"],)
        ).fetchone()
        historical_match = None
        if summary["matched_structure_id"] is not None:
            historical_match = historical.execute(
                "SELECT * FROM structure_file WHERE id=?",
                (summary["matched_structure_id"],),
            ).fetchone()
        if historical_summary is None or historical_match is None:
            raise CleanupError(f"Historical model_summary parents unavailable for row {summary['id']}")
        preview.append(
            {
                "table_name": "model_summary",
                "row_key": str(summary["id"]),
                "tf_id": tf_id,
                "source": str(summary["source"]),
                "motif_id": "",
                "model_id": str(summary["summary_model_id"]),
                "historical_motif_ref_id": "",
                "historical_structure_file_id": (
                    f"summary={summary['summary_file_id']};matched={summary['matched_structure_id']}"
                ),
                "classification": "excluded_unannotated_model_summary",
                "action": "DELETE_BROKEN_CHILD_ROW",
                "reason": (
                    "Both referenced historical structure rows are absent from current "
                    "staging; the frozen integrated TSV classifies this TF as Unannotated "
                    "and its current motif/model counts are zero."
                ),
            }
        )
    return preview, by_tf


def validate_plan(
    connection: sqlite3.Connection,
    link_preview: list[dict[str, str]],
    link_classes: Counter[str],
    summary_preview: list[dict[str, str]],
    summaries_by_tf: Counter[str],
) -> None:
    link_count = len(link_preview)
    summary_count = len(summary_preview)
    if link_count not in {0, EXPECTED_MOTIF_STRUCTURE_STALE}:
        raise CleanupError(
            f"Expected 823 stale motif_structure rows or idempotent 0; found {link_count}"
        )
    if summary_count not in {0, EXPECTED_MODEL_SUMMARY_STALE}:
        raise CleanupError(
            f"Expected 6 stale model_summary rows or idempotent 0; found {summary_count}"
        )
    if bool(link_count) != bool(summary_count):
        raise CleanupError("Partial cleanup state detected; refusing to continue")
    if link_count:
        if link_classes != EXPECTED_LINK_CLASSES:
            raise CleanupError(
                f"Stale-link provenance differs: expected {dict(EXPECTED_LINK_CLASSES)}, "
                f"found {dict(link_classes)}"
            )
        if summaries_by_tf != Counter(EXPECTED_SUMMARIES_BY_TF):
            raise CleanupError(
                f"Stale-summary provenance differs: expected {EXPECTED_SUMMARIES_BY_TF}, "
                f"found {dict(summaries_by_tf)}"
            )
        nondelivered = Counter(
            row["tf_id"]
            for row in link_preview
            if row["classification"] == "excluded_not_in_delivered_alphafold_tokens"
        )
        if dict(nondelivered) != EXPECTED_NONDELIVERED_MODELS:
            raise CleanupError(
                f"Non-delivered AlphaFold models differ: {dict(nondelivered)}"
            )
    violations = connection.execute("PRAGMA foreign_key_check").fetchall()
    expected_violations = 1658 if link_count else 0
    if len(violations) != expected_violations:
        raise CleanupError(
            f"Expected {expected_violations} initial FK violations; found {len(violations)}"
        )


def perform_cleanup(
    connection: sqlite3.Connection,
    link_preview: list[dict[str, str]],
    summary_preview: list[dict[str, str]],
) -> tuple[int, int]:
    motif_deleted = 0
    summary_deleted = 0
    for row in link_preview:
        cursor = connection.execute(
            "DELETE FROM motif_structure WHERE rowid=?", (int(row["row_key"]),)
        )
        motif_deleted += cursor.rowcount
    for row in summary_preview:
        cursor = connection.execute(
            "DELETE FROM model_summary WHERE id=?", (int(row["row_key"]),)
        )
        summary_deleted += cursor.rowcount
    return motif_deleted, summary_deleted


def write_preview(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=PREVIEW_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def make_backup(database: Path, backup_dir: Path) -> tuple[Path, str]:
    backup_dir.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    backup = backup_dir / f"{database.name}.before_stale_fk_cleanup_{stamp}.bak"
    if backup.exists():
        raise CleanupError(f"Backup already exists: {backup}")
    before = sha256_file(database)
    shutil.copy2(database, backup)
    backup_hash = sha256_file(backup)
    if before != backup_hash:
        raise CleanupError("Backup SHA256 does not match source database")
    return backup, before


def main() -> int:
    args = parse_args()
    for path in (args.database, args.historical_database, args.integrated_tsv):
        if not path.is_file():
            raise CleanupError(f"Required input missing: {path}")
    if sha256_file(args.integrated_tsv) != INTEGRATED_SHA256:
        raise CleanupError("Frozen integrated TSV SHA256 mismatch")
    if sha256_file(args.historical_database) != HISTORICAL_SHA256:
        raise CleanupError("Historical database SHA256 mismatch")
    for output in (args.preview_output, args.qc_output):
        if output.exists() and not args.replace_outputs:
            raise CleanupError(f"Output already exists (use --replace-outputs): {output}")

    database_hash_before = sha256_file(args.database)
    backup_path: Path | None = None
    backup_hash: str | None = None
    if args.apply:
        backup_path, backup_hash = make_backup(args.database, args.backup_dir)
        connection = sqlite3.connect(args.database)
        connection.row_factory = sqlite3.Row
        connection.execute("PRAGMA foreign_keys=ON")
    else:
        connection = in_memory_copy(args.database)

    historical = connect_ro(args.historical_database)
    integrated = read_integrated(args.integrated_tsv)
    protected_before = snapshots(connection)
    counts_before = {
        "motif_structure": connection.execute("SELECT COUNT(*) FROM motif_structure").fetchone()[0],
        "model_summary": connection.execute("SELECT COUNT(*) FROM model_summary").fetchone()[0],
    }

    try:
        link_preview, link_classes = classify_stale_links(
            connection, historical, integrated
        )
        summary_preview, summaries_by_tf = classify_stale_summaries(
            connection, historical, integrated
        )
        validate_plan(
            connection,
            link_preview,
            link_classes,
            summary_preview,
            summaries_by_tf,
        )
        connection.execute("BEGIN")
        motif_deleted, summary_deleted = perform_cleanup(
            connection, link_preview, summary_preview
        )

        fk_after = connection.execute("PRAGMA foreign_key_check").fetchall()
        quick_check = connection.execute("PRAGMA quick_check").fetchone()[0]
        protected_after = snapshots(connection)
        counts_after = {
            "motif_structure": connection.execute("SELECT COUNT(*) FROM motif_structure").fetchone()[0],
            "model_summary": connection.execute("SELECT COUNT(*) FROM model_summary").fetchone()[0],
        }
        if fk_after:
            raise CleanupError(f"Foreign-key violations remain after cleanup: {len(fk_after)}")
        if quick_check != "ok":
            raise CleanupError(f"SQLite quick_check failed: {quick_check}")
        if protected_before != protected_after:
            changed = [
                table for table in PROTECTED_TABLES
                if protected_before[table] != protected_after[table]
            ]
            raise CleanupError(f"Protected tables changed: {changed}")
        if counts_after != EXPECTED_FINAL_COUNTS:
            raise CleanupError(
                f"Unexpected final child-table counts: {counts_after}; "
                f"expected {EXPECTED_FINAL_COUNTS}"
            )
        if motif_deleted != len(link_preview) or summary_deleted != len(summary_preview):
            raise CleanupError("Deleted-row count differs from the validated plan")

        if args.apply:
            connection.commit()
        else:
            connection.rollback()
    except Exception:
        if connection.in_transaction:
            connection.rollback()
        raise
    finally:
        historical.close()
        connection.close()

    database_hash_after = sha256_file(args.database)
    if args.dry_run and database_hash_after != database_hash_before:
        raise CleanupError("Dry-run modified the staging database")

    all_preview = link_preview + summary_preview
    write_preview(args.preview_output, all_preview)
    qc = {
        "mode": "APPLY" if args.apply else "DRY-RUN",
        "database": str(args.database.resolve()),
        "database_sha256_before": database_hash_before,
        "database_sha256_after": database_hash_after,
        "backup": str(backup_path.resolve()) if backup_path else None,
        "backup_sha256": backup_hash,
        "backup_matches_pre_apply": bool(backup_hash == database_hash_before) if backup_hash else None,
        "integrated_tsv_sha256": INTEGRATED_SHA256,
        "historical_database_sha256": HISTORICAL_SHA256,
        "counts_before": counts_before,
        "counts_after_simulated_or_applied": counts_after,
        "motif_structure_rows_deleted": motif_deleted,
        "model_summary_rows_deleted": summary_deleted,
        "link_classification": dict(sorted(link_classes.items())),
        "model_summary_by_tf": dict(sorted(summaries_by_tf.items())),
        "foreign_key_violations_after": 0,
        "quick_check": quick_check,
        "protected_tables_unchanged": True,
        "files_or_parent_records_deleted": False,
        "invariants_passed": True,
    }
    write_json(args.qc_output, qc)

    print(f"Mode: {qc['mode']}")
    print(f"Broken motif_structure rows removed: {motif_deleted}")
    print(f"Broken model_summary rows removed: {summary_deleted}")
    print(f"Replacement relationships already preserved: {link_classes.get('replacement_already_present', 0)}")
    print(
        "Intentionally excluded historical relationships removed: "
        f"{motif_deleted - link_classes.get('replacement_already_present', 0)}"
    )
    print("Foreign-key violations after: 0")
    print(f"SQLite quick_check: {quick_check}")
    print("Protected parent/model/PWM tables unchanged: YES")
    print(f"Database SHA256 before: {database_hash_before}")
    print(f"Database SHA256 after: {database_hash_after}")
    print(f"Backup: {backup_path if backup_path else '<not created in dry-run>'}")
    print(f"Preview: {args.preview_output}")
    print(f"QC: {args.qc_output}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except CleanupError as error:
        raise SystemExit(f"ERROR: {error}")
