#!/usr/bin/env python3
"""Reconcile stale TF primary labels from verified current direct-PWM evidence.

The two delivered PWM charts are immutable evidence inputs.  This script opens
them read-only, verifies their frozen SHA256 values, and never rewrites them.

The reconciliation is intentionally narrow:

* identify TFs whose delivered ``Identical_PWM`` cells are empty;
* require a current, usable, non-missing ``motif_ref`` labelled ``identical``;
* promote TFs with either the frozen verified CIS-BP v2 signature or an exact
  TF-accession/motif pair in the frozen official JASPAR API verification
  manifest;
* mark only those exact JASPAR ``motif_ref`` rows as officially verified;
* update ``tf_primary_annotation.best_annotation_level`` and
  ``tf_primary_annotation.best_pwm_or_model``;
* preserve the delivered-row provenance fields
  ``n_nonempty_annotation_columns`` and ``source_table``.

Dry-run mode simulates the update in memory.  Apply mode requires an explicit
database path, creates a byte-identical timestamped backup, writes a row-level
audit TSV, and verifies that protected tables and immutable inputs did not
change.  Rerunning after a successful apply is a no-op.
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
from pathlib import Path
import shutil
import sqlite3
from typing import Iterable


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATABASE = REPO_ROOT / "data/tf_webdb_local_staging.sqlite"
DEFAULT_PREVIEW = REPO_ROOT / "outputs/tf_primary_annotation_reconciliation.tsv"
DEFAULT_BACKUP_DIR = REPO_ROOT / "backups"

ORIGINAL_TSV = Path(
    "/home/patricia/TF_database_Baldo_data/TF_PWM_chart_final.tsv"
)
INTEGRATED_TSV = Path(
    "/home/patricia/TF_database_Baldo_data/"
    "TF_PWM_chart_final_integrated_HOCOMOCO_hierarchical.tsv"
)
JASPAR_VERIFICATION_MANIFEST = (
    REPO_ROOT
    / "data_sources/jaspar_official_direct_mappings_verified_2026-10-08.tsv"
)

FROZEN_INPUT_SHA256 = {
    ORIGINAL_TSV: "1c9b59936d4e27e367ca81289b2c7b0bda7f68f86e2282261406194eed152f81",
    INTEGRATED_TSV: "f71160a2c2d318512be8354302208b55b4d4e7b427178631d98977c4209ed845",
    JASPAR_VERIFICATION_MANIFEST: "f62c353cd99de3b6f1d7b147ce5863a1d692d0dd979411f1d184d4eb00fc2919",
}

EXPECTED_RECONCILIATION_UNIVERSE = 145
EXPECTED_VERIFIED_PROMOTIONS = 145
EXPECTED_PENDING_ONLY = 0
EXPECTED_JASPAR_VERIFIED_TFS = 33
EXPECTED_JASPAR_VERIFIED_LINKS = 48
CONTROL_TF = "Q9Y603"
JASPAR_CONTROL_TF = "Q9BTP6"
JASPAR_VERIFIED_CURATION_STATUS = "verified_official_jaspar_api_uniprot_mapping"

VERIFIED_DIRECT_SIGNATURE = {
    "source": "cisbp",
    "curation_status": "verified_same_version_source_mapping",
    "mapping_type": "direct_cisbp_v2_tf_information",
    "original_column": "CISBP_v2_TF_Information_direct",
}

PROTECTED_TABLES = (
    "tf",
    "tf_family",
    "motif_file",
    "motif_structure",
    "structure_file",
    "structure_model_assignment",
    "structure_confidence_artifact",
    "tf_structure_status",
)

PREVIEW_COLUMNS = (
    "tf_id",
    "gene",
    "old_best_annotation_level",
    "new_best_annotation_level",
    "old_best_pwm_or_model",
    "new_best_pwm_or_model",
    "verified_direct_sources",
    "verified_direct_motif_ids",
    "verified_direct_curation_statuses",
    "verified_direct_bases",
    "official_jaspar_entry_urls",
    "pending_direct_sources",
    "pending_direct_motif_ids",
    "pending_direct_curation_statuses",
    "original_tsv_identical_pwm",
    "integrated_tsv_identical_pwm",
    "original_tsv_sha256",
    "integrated_tsv_sha256",
    "database_sha256_before",
    "action",
    "reason",
)


class ReconciliationError(RuntimeError):
    """Raised when an input or database invariant is not satisfied."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument(
        "--dry-run",
        action="store_true",
        help="simulate updates in memory; the on-disk database remains read-only",
    )
    mode.add_argument(
        "--apply",
        action="store_true",
        help="back up and update the explicitly selected database",
    )
    parser.add_argument(
        "--database",
        type=Path,
        default=DEFAULT_DATABASE,
        help=f"SQLite database (default: {DEFAULT_DATABASE})",
    )
    parser.add_argument(
        "--preview-output",
        type=Path,
        default=DEFAULT_PREVIEW,
        help=f"row-level reconciliation audit TSV (default: {DEFAULT_PREVIEW})",
    )
    parser.add_argument(
        "--backup-dir",
        type=Path,
        default=DEFAULT_BACKUP_DIR,
        help=f"backup directory used by --apply (default: {DEFAULT_BACKUP_DIR})",
    )
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_fingerprint(path: Path) -> tuple[int, int, str]:
    stat = path.stat()
    return stat.st_size, stat.st_mtime_ns, sha256_file(path)


def verify_frozen_inputs() -> dict[Path, tuple[int, int, str]]:
    fingerprints: dict[Path, tuple[int, int, str]] = {}
    for path, expected_sha256 in FROZEN_INPUT_SHA256.items():
        if not path.is_file():
            raise ReconciliationError(f"immutable input missing: {path}")
        fingerprint = file_fingerprint(path)
        if fingerprint[2] != expected_sha256:
            raise ReconciliationError(
                f"immutable input SHA256 mismatch for {path}: "
                f"expected {expected_sha256}, observed {fingerprint[2]}"
            )
        fingerprints[path] = fingerprint
    return fingerprints


def require_inputs_unchanged(
    before: dict[Path, tuple[int, int, str]],
) -> None:
    for path, fingerprint in before.items():
        observed = file_fingerprint(path)
        if observed != fingerprint:
            raise ReconciliationError(
                f"immutable input changed during reconciliation: {path}"
            )


def read_chart(path: Path) -> dict[str, dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            headers = [header.strip() for header in next(reader)]
        except StopIteration as exc:
            raise ReconciliationError(f"empty TSV: {path}") from exc
        if "TF_name" not in headers or "Identical_PWM" not in headers:
            raise ReconciliationError(
                f"{path} lacks TF_name/Identical_PWM columns after header normalization"
            )
        records: dict[str, dict[str, str]] = {}
        for values in reader:
            if not values:
                continue
            values += [""] * (len(headers) - len(values))
            row = {
                headers[index]: values[index].strip()
                for index in range(len(headers))
            }
            tf_id = row.get("TF_name", "")
            if not tf_id:
                continue
            if tf_id in records:
                raise ReconciliationError(f"duplicate TF_name {tf_id} in {path}")
            records[tf_id] = row
    return records


def read_jaspar_verification_manifest(
    path: Path,
) -> dict[tuple[str, str], dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "tf_id",
            "motif_id",
            "jaspar_name",
            "uniprot_ids",
            "species_tax_ids",
            "official_entry_url",
            "verified_date",
            "verification_result",
            "notes",
        }
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ReconciliationError(
                f"JASPAR verification manifest lacks required columns: {path}"
            )
        records: dict[tuple[str, str], dict[str, str]] = {}
        for raw in reader:
            row = {key: str(value or "").strip() for key, value in raw.items()}
            key = (row["tf_id"], row["motif_id"])
            if not all(key) or key in records:
                raise ReconciliationError(
                    f"invalid or duplicate JASPAR verification key: {key}"
                )
            if row["verification_result"] != "EXACT_UNIPROT_MATCH":
                raise ReconciliationError(
                    f"non-exact JASPAR verification result for {key}"
                )
            if row["tf_id"] not in row["uniprot_ids"].split(";"):
                raise ReconciliationError(
                    f"JASPAR manifest does not list {row['tf_id']} for {key}"
                )
            expected_url = (
                "https://jaspar.elixir.no/api/v1/matrix/"
                f"{row['motif_id']}/"
            )
            if row["official_entry_url"] != expected_url:
                raise ReconciliationError(
                    f"unexpected official JASPAR URL for {key}"
                )
            records[key] = row
    if len(records) != EXPECTED_JASPAR_VERIFIED_LINKS:
        raise ReconciliationError(
            f"expected {EXPECTED_JASPAR_VERIFIED_LINKS} verified JASPAR links, "
            f"observed {len(records)}"
        )
    tf_count = len({tf_id for tf_id, _ in records})
    if tf_count != EXPECTED_JASPAR_VERIFIED_TFS:
        raise ReconciliationError(
            f"expected {EXPECTED_JASPAR_VERIFIED_TFS} verified JASPAR TFs, "
            f"observed {tf_count}"
        )
    return records


def jaspar_evidence_note(record: dict[str, str]) -> str:
    note = (
        "Official JASPAR API exact UniProt mapping verified "
        f"{record['verified_date']}: {record['official_entry_url']} lists "
        f"UniProt {record['tf_id']}; JASPAR species taxon(s) "
        f"{record['species_tax_ids']}."
    )
    if record["notes"]:
        note += f" {record['notes']}"
    return note


def connect_read_only(path: Path) -> sqlite3.Connection:
    connection = sqlite3.connect(
        f"file:{path.resolve()}?mode=ro&immutable=1", uri=True
    )
    connection.row_factory = sqlite3.Row
    connection.execute("PRAGMA query_only = ON")
    return connection


def required_tables(connection: sqlite3.Connection) -> None:
    expected = {"tf", "tf_primary_annotation", "motif_ref", "motif_file"}
    observed = {
        row[0]
        for row in connection.execute(
            "SELECT name FROM sqlite_master WHERE type='table'"
        )
    }
    missing = sorted(expected - observed)
    if missing:
        raise ReconciliationError(
            f"database lacks required tables: {', '.join(missing)}"
        )


def table_exists(connection: sqlite3.Connection, table: str) -> bool:
    return (
        connection.execute(
            "SELECT 1 FROM sqlite_master WHERE type='table' AND name=?", (table,)
        ).fetchone()
        is not None
    )


def table_digest(connection: sqlite3.Connection, table: str) -> tuple[int, str]:
    if not table_exists(connection, table):
        return -1, "<TABLE_ABSENT>"
    columns = [
        row[1] for row in connection.execute(f'PRAGMA table_info("{table}")')
    ]
    if not columns:
        return 0, hashlib.sha256(b"").hexdigest()
    quoted = ",".join(f'"{column}"' for column in columns)
    digest = hashlib.sha256()
    count = 0
    for row in connection.execute(f'SELECT {quoted} FROM "{table}" ORDER BY {quoted}'):
        digest.update(repr(tuple(row)).encode("utf-8"))
        digest.update(b"\n")
        count += 1
    return count, digest.hexdigest()


def protected_snapshot(connection: sqlite3.Connection) -> dict[str, tuple[int, str]]:
    return {table: table_digest(connection, table) for table in PROTECTED_TABLES}


def motif_ref_partition_snapshots(
    connection: sqlite3.Connection,
    jaspar_verified: dict[tuple[str, str], dict[str, str]],
) -> tuple[tuple[int, str], tuple[int, str]]:
    columns = [row[1] for row in connection.execute("PRAGMA table_info(motif_ref)")]
    mutable = {"curation_status", "evidence_note"}
    immutable_columns = [column for column in columns if column not in mutable]
    quoted = ",".join(f'"{column}"' for column in columns)
    non_target_digest = hashlib.sha256()
    target_immutable_digest = hashlib.sha256()
    non_target_count = 0
    target_count = 0
    seen: set[tuple[str, str]] = set()
    for row in connection.execute(f"SELECT {quoted} FROM motif_ref ORDER BY id"):
        record = dict(row)
        key = (str(record["tf_id"]), str(record["motif_id"]))
        is_target = str(record["source"]) == "jaspar" and key in jaspar_verified
        if is_target:
            if key in seen:
                raise ReconciliationError(f"duplicate target JASPAR motif_ref row: {key}")
            seen.add(key)
            if (
                str(record["evidence_type"]) != "identical"
                or str(record["original_column"]) != "JASPAR2024_metadata"
                or str(record["mapping_type"]) != "direct_or_identical"
                or int(record["missing_local_file"] or 0) != 0
            ):
                raise ReconciliationError(
                    f"target JASPAR motif_ref metadata is not direct/usable for {key}"
                )
            target_immutable_digest.update(
                repr(tuple(record[column] for column in immutable_columns)).encode("utf-8")
            )
            target_immutable_digest.update(b"\n")
            target_count += 1
        else:
            non_target_digest.update(
                repr(tuple(record[column] for column in columns)).encode("utf-8")
            )
            non_target_digest.update(b"\n")
            non_target_count += 1
    missing = sorted(set(jaspar_verified) - seen)
    if missing:
        raise ReconciliationError(
            f"verified JASPAR links missing exact motif_ref rows: {missing}"
        )
    if target_count != EXPECTED_JASPAR_VERIFIED_LINKS:
        raise ReconciliationError(
            f"expected {EXPECTED_JASPAR_VERIFIED_LINKS} target JASPAR motif_ref rows, "
            f"observed {target_count}"
        )
    return (
        (non_target_count, non_target_digest.hexdigest()),
        (target_count, target_immutable_digest.hexdigest()),
    )


def validate_jaspar_curation(
    connection: sqlite3.Connection,
    jaspar_verified: dict[tuple[str, str], dict[str, str]],
) -> None:
    for key, manifest_row in sorted(jaspar_verified.items()):
        rows = connection.execute(
            """
            SELECT curation_status, evidence_note
            FROM motif_ref
            WHERE tf_id = ? AND source = 'jaspar' AND motif_id = ?
              AND evidence_type = 'identical'
            """,
            key,
        ).fetchall()
        if len(rows) != 1:
            raise ReconciliationError(
                f"expected one exact JASPAR motif_ref row for {key}, observed {len(rows)}"
            )
        if str(rows[0]["curation_status"] or "") != JASPAR_VERIFIED_CURATION_STATUS:
            raise ReconciliationError(f"JASPAR curation status not finalized for {key}")
        if str(rows[0]["evidence_note"] or "") != jaspar_evidence_note(manifest_row):
            raise ReconciliationError(f"JASPAR evidence note not finalized for {key}")


def direct_rows(connection: sqlite3.Connection) -> dict[str, list[dict[str, object]]]:
    by_tf: dict[str, list[dict[str, object]]] = {}
    for row in connection.execute(
        """
        SELECT
            mr.id AS motif_ref_id,
            mr.tf_id,
            mr.source,
            mr.motif_id,
            mr.original_column,
            mr.mapping_type,
            mr.curation_status,
            mr.evidence_note,
            mr.display_priority
        FROM motif_ref AS mr
        JOIN motif_file AS mf
          ON mf.source = mr.source
         AND mf.motif_id = mr.motif_id
        WHERE mr.evidence_type = 'identical'
          AND mr.missing_local_file = 0
          AND mf.matrix_status = 'usable'
        ORDER BY mr.tf_id, mr.display_priority, mr.source, mr.motif_id, mr.id
        """
    ):
        by_tf.setdefault(str(row["tf_id"]), []).append(dict(row))
    return by_tf


def verified_direct_basis(
    row: dict[str, object],
    jaspar_verified: dict[tuple[str, str], dict[str, str]],
) -> str:
    if all(
        str(row.get(key) or "") == value
        for key, value in VERIFIED_DIRECT_SIGNATURE.items()
    ):
        return "CISBP_V2_VERIFIED_DIRECT"
    key = (str(row.get("tf_id") or ""), str(row.get("motif_id") or ""))
    if (
        str(row.get("source") or "") == "jaspar"
        and str(row.get("original_column") or "") == "JASPAR2024_metadata"
        and str(row.get("mapping_type") or "") == "direct_or_identical"
        and key in jaspar_verified
    ):
        return "OFFICIAL_JASPAR_API_EXACT_UNIPROT"
    return ""


def is_verified_direct(
    row: dict[str, object],
    jaspar_verified: dict[tuple[str, str], dict[str, str]],
) -> bool:
    return bool(verified_direct_basis(row, jaspar_verified))


def joined_unique(values: Iterable[object]) -> str:
    return ";".join(dict.fromkeys(str(value) for value in values if str(value)))


def build_plan(
    connection: sqlite3.Connection,
    original: dict[str, dict[str, str]],
    integrated: dict[str, dict[str, str]],
    jaspar_verified: dict[tuple[str, str], dict[str, str]],
    database_sha256: str,
) -> list[dict[str, str]]:
    required_tables(connection)
    primary = {
        str(row["tf_id"]): dict(row)
        for row in connection.execute(
            """
            SELECT tf_id, best_annotation_level, best_pwm_or_model,
                   n_nonempty_annotation_columns, source_table
            FROM tf_primary_annotation
            """
        )
    }
    genes = {
        str(row["tf_id"]): str(row["gene_names"] or "").split()[0]
        if str(row["gene_names"] or "").strip()
        else ""
        for row in connection.execute(
            """
            SELECT tf.tf_id, ta.gene_names
            FROM tf
            LEFT JOIN tf_annotation AS ta ON ta.tf_id = tf.tf_id
            """
        )
    }
    current_direct = direct_rows(connection)
    original_sha = FROZEN_INPUT_SHA256[ORIGINAL_TSV]
    integrated_sha = FROZEN_INPUT_SHA256[INTEGRATED_TSV]
    plan: list[dict[str, str]] = []

    for tf_id, evidence in sorted(current_direct.items()):
        if tf_id not in original or tf_id not in integrated or tf_id not in primary:
            continue
        original_identical = original[tf_id].get("Identical_PWM", "").strip()
        integrated_identical = integrated[tf_id].get("Identical_PWM", "").strip()
        if original_identical or integrated_identical:
            continue

        verified = [
            row for row in evidence if is_verified_direct(row, jaspar_verified)
        ]
        pending = [
            row for row in evidence if not is_verified_direct(row, jaspar_verified)
        ]
        old_level = str(primary[tf_id]["best_annotation_level"] or "")
        old_best = str(primary[tf_id]["best_pwm_or_model"] or "")
        if verified:
            new_level = "Identical_PWM"
            new_best = joined_unique(row["motif_id"] for row in verified)
            if old_level == new_level and old_best == new_best:
                action = "NO_CHANGE_ALREADY_RECONCILED"
                reason = (
                    "verified current direct evidence is already reflected in the primary annotation"
                )
            else:
                action = "UPDATE_TO_IDENTICAL_PWM"
                bases = {
                    verified_direct_basis(row, jaspar_verified) for row in verified
                }
                reason = (
                    "delivered Identical_PWM cells are empty, but current usable motif_ref "
                    "contains verified exact same-TF direct evidence from "
                    + " and ".join(sorted(bases))
                )
        else:
            new_level = old_level
            new_best = old_best
            action = "PENDING_CONFIRMATION_NO_UPDATE"
            reason = (
                "usable identical evidence exists, but none has the frozen verified-direct "
                "signature; retain the current primary annotation pending source confirmation"
            )

        plan.append(
            {
                "tf_id": tf_id,
                "gene": genes.get(tf_id, ""),
                "old_best_annotation_level": old_level,
                "new_best_annotation_level": new_level,
                "old_best_pwm_or_model": old_best,
                "new_best_pwm_or_model": new_best,
                "verified_direct_sources": joined_unique(row["source"] for row in verified),
                "verified_direct_motif_ids": joined_unique(row["motif_id"] for row in verified),
                "verified_direct_curation_statuses": joined_unique(
                    row["curation_status"] for row in verified
                ),
                "verified_direct_bases": joined_unique(
                    verified_direct_basis(row, jaspar_verified) for row in verified
                ),
                "official_jaspar_entry_urls": joined_unique(
                    jaspar_verified[(str(row["tf_id"]), str(row["motif_id"]))][
                        "official_entry_url"
                    ]
                    for row in verified
                    if str(row["source"]) == "jaspar"
                ),
                "pending_direct_sources": joined_unique(row["source"] for row in pending),
                "pending_direct_motif_ids": joined_unique(row["motif_id"] for row in pending),
                "pending_direct_curation_statuses": joined_unique(
                    row["curation_status"] for row in pending
                ),
                "original_tsv_identical_pwm": original_identical,
                "integrated_tsv_identical_pwm": integrated_identical,
                "original_tsv_sha256": original_sha,
                "integrated_tsv_sha256": integrated_sha,
                "database_sha256_before": database_sha256,
                "action": action,
                "reason": reason,
            }
        )

    if len(plan) != EXPECTED_RECONCILIATION_UNIVERSE:
        raise ReconciliationError(
            f"expected {EXPECTED_RECONCILIATION_UNIVERSE} reconciliation TFs, "
            f"observed {len(plan)}"
        )
    verified_count = sum(bool(row["verified_direct_motif_ids"]) for row in plan)
    pending_only_count = sum(
        not row["verified_direct_motif_ids"] and bool(row["pending_direct_motif_ids"])
        for row in plan
    )
    if verified_count != EXPECTED_VERIFIED_PROMOTIONS:
        raise ReconciliationError(
            f"expected {EXPECTED_VERIFIED_PROMOTIONS} verified TFs, observed {verified_count}"
        )
    if pending_only_count != EXPECTED_PENDING_ONLY:
        raise ReconciliationError(
            f"expected {EXPECTED_PENDING_ONLY} pending-only TFs, observed {pending_only_count}"
        )
    control = next((row for row in plan if row["tf_id"] == CONTROL_TF), None)
    if control is None or control["new_best_annotation_level"] != "Identical_PWM":
        raise ReconciliationError(
            f"control {CONTROL_TF} is absent or is not promoted to Identical_PWM"
        )
    expected_control_motifs = {
        "M04712_2.00",
        "M04713_2.00",
        "M04714_2.00",
        "M04715_2.00",
        "MA1708.1",
        "MA1708.2",
    }
    if set(control["verified_direct_motif_ids"].split(";")) != expected_control_motifs:
        raise ReconciliationError(
            f"control {CONTROL_TF} verified motifs do not match the frozen direct set"
        )
    jaspar_control = next(
        (row for row in plan if row["tf_id"] == JASPAR_CONTROL_TF), None
    )
    if (
        jaspar_control is None
        or jaspar_control["new_best_annotation_level"] != "Identical_PWM"
        or set(jaspar_control["verified_direct_motif_ids"].split(";"))
        != {"MA1971.1", "MA1971.2"}
    ):
        raise ReconciliationError(
            f"control {JASPAR_CONTROL_TF} is not resolved to its two official JASPAR motifs"
        )
    return plan


def write_preview(path: Path, plan: list[dict[str, str]]) -> None:
    path = path.expanduser().resolve()
    if path in {input_path.resolve() for input_path in FROZEN_INPUT_SHA256}:
        raise ReconciliationError("preview output may not overwrite an immutable input")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=PREVIEW_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(plan)


def copy_to_memory(source: sqlite3.Connection) -> sqlite3.Connection:
    target = sqlite3.connect(":memory:")
    target.row_factory = sqlite3.Row
    source.backup(target)
    target.execute("PRAGMA foreign_keys = ON")
    return target


def apply_plan(
    connection: sqlite3.Connection,
    plan: list[dict[str, str]],
    jaspar_verified: dict[tuple[str, str], dict[str, str]],
) -> tuple[int, int]:
    updates = [row for row in plan if row["action"] == "UPDATE_TO_IDENTICAL_PWM"]
    motif_ref_updates = 0
    with connection:
        for key, manifest_row in sorted(jaspar_verified.items()):
            current = connection.execute(
                """
                SELECT curation_status, evidence_note
                FROM motif_ref
                WHERE tf_id = ? AND source = 'jaspar' AND motif_id = ?
                  AND evidence_type = 'identical'
                  AND original_column = 'JASPAR2024_metadata'
                  AND mapping_type = 'direct_or_identical'
                  AND missing_local_file = 0
                """,
                key,
            ).fetchall()
            if len(current) != 1:
                raise ReconciliationError(
                    f"expected one updatable JASPAR motif_ref row for {key}, "
                    f"observed {len(current)}"
                )
            expected_note = jaspar_evidence_note(manifest_row)
            if (
                str(current[0]["curation_status"] or "")
                == JASPAR_VERIFIED_CURATION_STATUS
                and str(current[0]["evidence_note"] or "") == expected_note
            ):
                continue
            cursor = connection.execute(
                """
                UPDATE motif_ref
                   SET curation_status = ?, evidence_note = ?
                 WHERE tf_id = ? AND source = 'jaspar' AND motif_id = ?
                   AND evidence_type = 'identical'
                   AND original_column = 'JASPAR2024_metadata'
                   AND mapping_type = 'direct_or_identical'
                   AND missing_local_file = 0
                """,
                (
                    JASPAR_VERIFIED_CURATION_STATUS,
                    expected_note,
                    key[0],
                    key[1],
                ),
            )
            if cursor.rowcount != 1:
                raise ReconciliationError(
                    f"concurrent or unexpected JASPAR motif_ref state for {key}"
                )
            motif_ref_updates += 1
        for row in updates:
            cursor = connection.execute(
                """
                UPDATE tf_primary_annotation
                   SET best_annotation_level = ?,
                       best_pwm_or_model = ?
                 WHERE tf_id = ?
                   AND COALESCE(best_annotation_level, '') = ?
                   AND COALESCE(best_pwm_or_model, '') = ?
                """,
                (
                    row["new_best_annotation_level"],
                    row["new_best_pwm_or_model"],
                    row["tf_id"],
                    row["old_best_annotation_level"],
                    row["old_best_pwm_or_model"],
                ),
            )
            if cursor.rowcount != 1:
                raise ReconciliationError(
                    f"concurrent or unexpected primary row state for {row['tf_id']}"
                )
    return len(updates), motif_ref_updates


def validate_post_state(
    connection: sqlite3.Connection,
    original: dict[str, dict[str, str]],
    integrated: dict[str, dict[str, str]],
    jaspar_verified: dict[tuple[str, str], dict[str, str]],
    database_sha256: str,
) -> list[dict[str, str]]:
    post_plan = build_plan(
        connection,
        original,
        integrated,
        jaspar_verified,
        database_sha256=database_sha256,
    )
    remaining_updates = [
        row for row in post_plan if row["action"] == "UPDATE_TO_IDENTICAL_PWM"
    ]
    if remaining_updates:
        raise ReconciliationError(
            f"post-update validation still finds {len(remaining_updates)} verified promotions"
        )
    if sum(
        row["action"] == "PENDING_CONFIRMATION_NO_UPDATE" for row in post_plan
    ) != EXPECTED_PENDING_ONLY:
        raise ReconciliationError("pending-only case count changed unexpectedly")
    validate_jaspar_curation(connection, jaspar_verified)
    return post_plan


def create_backup(database: Path, backup_dir: Path) -> Path:
    backup_dir.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    backup = backup_dir / f"{database.name}.before_primary_reconciliation_{stamp}.bak"
    if backup.exists():
        raise ReconciliationError(f"refusing to overwrite backup: {backup}")
    shutil.copy2(database, backup)
    if sha256_file(backup) != sha256_file(database):
        raise ReconciliationError("backup SHA256 does not match the pre-update database")
    return backup


def print_summary(
    *,
    mode: str,
    plan: list[dict[str, str]],
    primary_changes: int,
    motif_ref_changes: int,
    preview: Path,
    database_sha_before: str,
    database_sha_after: str,
    backup: Path | None,
) -> None:
    print(f"Mode: {mode}")
    print(f"Reconciliation universe: {len(plan)}")
    print(
        "Verified direct TFs: "
        f"{sum(bool(row['verified_direct_motif_ids']) for row in plan)}"
    )
    print(
        "Pending-only TFs: "
        f"{sum(not row['verified_direct_motif_ids'] for row in plan)}"
    )
    print(f"Primary rows changed: {primary_changes}")
    print(f"Official JASPAR TFs: {EXPECTED_JASPAR_VERIFIED_TFS}")
    print(f"Official JASPAR links: {EXPECTED_JASPAR_VERIFIED_LINKS}")
    print(f"JASPAR motif_ref rows changed: {motif_ref_changes}")
    print(f"Q9Y603 target: Identical_PWM")
    print(f"Q9BTP6 target: Identical_PWM")
    print(f"Original TSVs modified: NO")
    print(f"Preview: {preview.resolve()}")
    print(f"Database SHA256 before: {database_sha_before}")
    print(f"Database SHA256 after: {database_sha_after}")
    print(f"Backup: {backup.resolve() if backup else '<not created in dry-run>'}")


def main() -> None:
    args = parse_args()
    database = args.database.expanduser().resolve()
    preview = args.preview_output.expanduser().resolve()
    backup_dir = args.backup_dir.expanduser().resolve()
    if not database.is_file():
        raise ReconciliationError(f"database not found: {database}")
    if database in {input_path.resolve() for input_path in FROZEN_INPUT_SHA256}:
        raise ReconciliationError("database path resolves to an immutable input")

    frozen_before = verify_frozen_inputs()
    database_fingerprint_before = file_fingerprint(database)
    database_sha_before = database_fingerprint_before[2]
    original = read_chart(ORIGINAL_TSV)
    integrated = read_chart(INTEGRATED_TSV)
    jaspar_verified = read_jaspar_verification_manifest(
        JASPAR_VERIFICATION_MANIFEST
    )

    backup: Path | None = None
    if args.dry_run:
        with connect_read_only(database) as disk_connection:
            protected_before = protected_snapshot(disk_connection)
            motif_ref_partitions_before = motif_ref_partition_snapshots(
                disk_connection, jaspar_verified
            )
            plan = build_plan(
                disk_connection,
                original,
                integrated,
                jaspar_verified,
                database_sha_before,
            )
            memory = copy_to_memory(disk_connection)
        try:
            primary_changes, motif_ref_changes = apply_plan(
                memory, plan, jaspar_verified
            )
            validate_post_state(
                memory,
                original,
                integrated,
                jaspar_verified,
                database_sha_before,
            )
            protected_after = protected_snapshot(memory)
            if protected_after != protected_before:
                changed = sorted(
                    table
                    for table in protected_before
                    if protected_before[table] != protected_after[table]
                )
                raise ReconciliationError(
                    f"protected tables changed in dry-run: {', '.join(changed)}"
                )
            motif_ref_partitions_after = motif_ref_partition_snapshots(
                memory, jaspar_verified
            )
            if motif_ref_partitions_after != motif_ref_partitions_before:
                raise ReconciliationError(
                    "non-target or immutable motif_ref content changed in dry-run"
                )
        finally:
            memory.close()
        if file_fingerprint(database) != database_fingerprint_before:
            raise ReconciliationError("on-disk database changed during dry-run")
        database_sha_after = database_sha_before
        mode = "DRY-RUN"
    else:
        backup = create_backup(database, backup_dir)
        connection = sqlite3.connect(database)
        connection.row_factory = sqlite3.Row
        connection.execute("PRAGMA foreign_keys = ON")
        try:
            protected_before = protected_snapshot(connection)
            motif_ref_partitions_before = motif_ref_partition_snapshots(
                connection, jaspar_verified
            )
            plan = build_plan(
                connection,
                original,
                integrated,
                jaspar_verified,
                database_sha_before,
            )
            primary_changes, motif_ref_changes = apply_plan(
                connection, plan, jaspar_verified
            )
            validate_post_state(
                connection,
                original,
                integrated,
                jaspar_verified,
                database_sha_before,
            )
            protected_after = protected_snapshot(connection)
            if protected_after != protected_before:
                changed = sorted(
                    table
                    for table in protected_before
                    if protected_before[table] != protected_after[table]
                )
                raise ReconciliationError(
                    f"protected tables changed during apply: {', '.join(changed)}"
                )
            motif_ref_partitions_after = motif_ref_partition_snapshots(
                connection, jaspar_verified
            )
            if motif_ref_partitions_after != motif_ref_partitions_before:
                raise ReconciliationError(
                    "non-target or immutable motif_ref content changed during apply"
                )
            if connection.execute("PRAGMA quick_check").fetchone()[0] != "ok":
                raise ReconciliationError("post-update PRAGMA quick_check failed")
        finally:
            connection.close()
        database_sha_after = sha256_file(database)
        mode = "APPLY"

    require_inputs_unchanged(frozen_before)
    write_preview(preview, plan)
    print_summary(
        mode=mode,
        plan=plan,
        primary_changes=primary_changes,
        motif_ref_changes=motif_ref_changes,
        preview=preview,
        database_sha_before=database_sha_before,
        database_sha_after=database_sha_after,
        backup=backup,
    )


if __name__ == "__main__":
    main()
