#!/usr/bin/env python3
"""Reconcile the frozen 114-model AF3 set into the local staging database."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import sqlite3
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATABASE = ROOT / "data/tf_webdb_local_staging.sqlite"
SOURCE_DATABASE_LABEL = ROOT / "data/tf_webdb.sqlite"
BASE_SELECTED_INPUT = ROOT / "outputs/local_af3_model_import_selected.tsv"
FINAL_TRUTH_INPUT = ROOT / "outputs/baldo_3690_STRUCTURAL_TRUTH_CURRENT.tsv"
VALID_INTERFACE_INPUT = ROOT / "outputs/baldo_3690_af3_valid_interface.tsv"
COMPLETE_AF3_AUDIT_INPUT = ROOT / "outputs/baldo_af3_complete_inventory.tsv"
FINAL_ADJUDICATION_INPUT = ROOT / "outputs/baldo_109_final_structural_adjudication.tsv"
FINAL_FRAGMENT_AUDIT_INPUT = ROOT / "outputs/baldo_AF3_unresolved_fragment_model_audit.tsv"
DEFAULT_QC_OUTPUT = ROOT / "outputs/local_af3_staging_import_qc.json"
DEFAULT_SUMMARY_OUTPUT = ROOT / "outputs/local_af3_staging_import_summary.tsv"

PASS_RESULT = "PASS_NONEMPTY_PROTEIN_DNA_INTERFACE"
FINAL_FULL_STATUS = "VALID_AF3_FULL_LENGTH"
FINAL_FRAGMENT_STATUS = "VALID_NEW_AF3_DBD_MODEL"
FRAGMENT_RUN_PATTERN = re.compile(
    r"^(?P<accession>[A-Za-z0-9]+)_(?P<start>\d+)_(?P<end>\d+)$"
)
EXPECTED_ROLE_COUNTS = {"OWN_ACCESSION_MODEL": 113, "DBD_FRAGMENT_MODEL": 1}
EXPECTED_ARTIFACT_COUNTS = {"MMCIF": 114, "PLDDT_JSON": 114, "PAE_JSON": 114}
EXPECTED_FINAL_ACCESSION_COUNT = 114
PROTECTED_PWM_TABLES = ("motif_ref", "motif_file", "motif_structure", "tf_primary_annotation")


class ImportValidationError(RuntimeError):
    """Raised when the staging-only import fails a safety or QC condition."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Reconcile the frozen 113 full-length plus one DBD-fragment AF3 "
            "assignments and their available confidence artifacts into local staging."
        )
    )
    parser.add_argument("--database", type=Path, default=DEFAULT_DATABASE)
    parser.add_argument("--qc-output", type=Path, default=DEFAULT_QC_OUTPUT)
    parser.add_argument("--summary-output", type=Path, default=DEFAULT_SUMMARY_OUTPUT)
    parser.add_argument(
        "--expected-preimport-sha256",
        default="",
        help="Optional lowercase SHA256 required to match the target database before import.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help=(
            "Simulate the reconciliation in an in-memory copy opened from staging "
            "read-only; do not write the database or report files."
        ),
    )
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fingerprint(path: Path) -> dict[str, Any]:
    stat = path.stat()
    return {
        "path": str(path),
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "sha256": sha256_file(path),
    }


def read_tsv(path: Path) -> list[dict[str, str]]:
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None:
                raise ImportValidationError(f"TSV has no header: {path}")
            return [dict(row) for row in reader]
    except (OSError, csv.Error) as exc:
        raise ImportValidationError(f"Could not read {path}: {exc}") from exc


def application_table_counts(connection: sqlite3.Connection) -> dict[str, int]:
    tables = [
        row[0]
        for row in connection.execute(
            "SELECT name FROM sqlite_master "
            "WHERE type='table' AND name NOT LIKE 'sqlite_%' ORDER BY name"
        )
    ]
    return {
        table: int(connection.execute(f'SELECT COUNT(*) FROM "{table}"').fetchone()[0])
        for table in tables
    }


def foreign_key_signature(connection: sqlite3.Connection) -> list[tuple[Any, ...]]:
    return sorted(tuple(row) for row in connection.execute("PRAGMA foreign_key_check"))


def readable_file(path: Path) -> bool:
    if not path.is_file():
        return False
    try:
        with path.open("rb") as handle:
            handle.read(1)
    except OSError:
        return False
    return True


def full_length_ownership_ok(accession: str, path: Path) -> bool:
    sample_dir = path.parent.name
    if re.fullmatch(r"seed-\d+_sample-\d+", sample_dir, flags=re.IGNORECASE):
        run_dir = path.parent.parent
        expected_name = f"{run_dir.name}_{sample_dir}_model.cif"
    else:
        run_dir = path.parent
        expected_name = f"{run_dir.name}_model.cif"
    normalized = re.sub(r"_(dimer|monomer)$", "", run_dir.name, flags=re.IGNORECASE)
    return (
        normalized.casefold() == accession.casefold()
        and path.name.casefold() == expected_name.casefold()
    )


def fragment_parts(accession: str, path: Path) -> tuple[int, int]:
    if path.name != "model.cif" or not re.fullmatch(
        r"seed-\d+_sample-\d+", path.parent.name, flags=re.IGNORECASE
    ):
        raise ImportValidationError(f"Unexpected final fragment path layout: {path}")
    run_name = path.parent.parent.name
    match = FRAGMENT_RUN_PATTERN.fullmatch(run_name)
    if not match or match.group("accession").casefold() != accession.casefold():
        raise ImportValidationError(f"Fragment path/accession mismatch: {accession}/{path}")
    return int(match.group("start")), int(match.group("end"))


def exact_database_mapping_ok(connection: sqlite3.Connection, accession: str) -> bool:
    mappings = connection.execute(
        "SELECT tf.tf_id, ta.uniprot_accession FROM tf "
        "JOIN tf_annotation AS ta ON ta.tf_id=tf.tf_id "
        "WHERE UPPER(tf.tf_id)=UPPER(?)",
        (accession,),
    ).fetchall()
    return (
        len(mappings) == 1
        and mappings[0]["tf_id"] == accession
        and mappings[0]["uniprot_accession"] == accession
    )


def validate_inputs(connection: sqlite3.Connection) -> tuple[list[dict[str, str]], set[str]]:
    truth_rows = read_tsv(FINAL_TRUTH_INPUT)
    full_rows = [row for row in truth_rows if row["current_structural_status"] == FINAL_FULL_STATUS]
    full_accessions = {row["accession"] for row in full_rows}
    if len(full_rows) != 113 or len(full_accessions) != 113:
        raise ImportValidationError(
            f"Frozen truth must contain 113 unique {FINAL_FULL_STATUS} accessions."
        )

    base_selection = read_tsv(BASE_SELECTED_INPUT)
    base_full_rows = [row for row in base_selection if row["candidate_type"] == "FULL_LENGTH_AF3"]
    base_full_by_accession = {row["accession"]: row for row in base_full_rows}
    if len(base_full_rows) != 109 or len(base_full_by_accession) != 109:
        raise ImportValidationError("Existing preselected manifest must contain 109 unique full-length models.")
    if not set(base_full_by_accession).issubset(full_accessions):
        raise ImportValidationError("A preselected full-length accession is outside frozen final truth.")
    prior_fragments = [row for row in base_selection if row["candidate_type"] == "DBD_FRAGMENT_AF3"]
    if (
        len(prior_fragments) != 1
        or prior_fragments[0]["accession"] != "H7C4N4"
        or "H7C4N4_PF03299_2-70" not in prior_fragments[0]["candidate_model_path"]
    ):
        raise ImportValidationError("Prior manifest H7C4N4 fragment control differs.")
    newly_resolved_full = full_accessions - set(base_full_by_accession)
    if newly_resolved_full != {"A0A494C1L3", "H2BNB9", "Q709A9", "H7C4N4"}:
        raise ImportValidationError(
            f"Newly resolved full-length set differs: {sorted(newly_resolved_full)}"
        )

    interface_rows = read_tsv(VALID_INTERFACE_INPUT)
    interface_by_accession = {row["accession"]: row for row in interface_rows}
    if len(interface_by_accession) != len(interface_rows):
        raise ImportValidationError("Valid-interface manifest contains duplicate accessions.")

    complete_audit_rows = read_tsv(COMPLETE_AF3_AUDIT_INPUT)
    complete_audit_by_path = {
        (row["accession"], row["model_path"]): row for row in complete_audit_rows
    }
    if len(complete_audit_by_path) != len(complete_audit_rows):
        raise ImportValidationError("Complete AF3 audit contains duplicate accession/model paths.")

    selected: list[dict[str, str]] = []
    for accession in sorted(full_accessions):
        prior = base_full_by_accession.get(accession)
        if prior is not None:
            if (
                prior["tf_id"] != accession
                or prior["interface_result"] != PASS_RESULT
                or prior["proposed_assignment_role"] != "OWN_ACCESSION_MODEL"
                or prior["proposed_is_active"] != "YES"
                or prior["decision"] != "PROPOSE_IMPORT_ACTIVE"
            ):
                raise ImportValidationError(f"Existing preselected row is invalid: {accession}")
            model_path = prior["candidate_model_path"]
            model_rank = prior["candidate_model_rank"]
            source_audit = str(BASE_SELECTED_INPUT)
        else:
            summary = interface_by_accession.get(accession)
            if summary is None:
                raise ImportValidationError(f"No selected-interface row for {accession}.")
            if (
                summary["valid_interface_any_model"] != "YES"
                or summary["final_return_status"] != "VALID_AF3_PROTEIN_DNA_MODEL"
            ):
                raise ImportValidationError(f"Selected-interface summary is invalid: {accession}")
            model_path = summary["best_model_path"]
            model_rank = summary["best_model_rank"]
            source_audit = str(VALID_INTERFACE_INPUT)

        audit = complete_audit_by_path.get((accession, model_path))
        if audit is None:
            raise ImportValidationError(f"Selected model is absent from complete AF3 audit: {accession}")
        path = Path(model_path)
        if (
            audit["file_exists"] != "YES"
            or audit["readable_cif"] != "YES"
            or audit["protein_present"] != "YES"
            or audit["dna_present"] != "YES"
            or audit["interface_result"] != PASS_RESULT
        ):
            raise ImportValidationError(f"Selected full-length audit evidence is invalid: {accession}")
        evidence = prior if prior is not None else interface_by_accession[accession]
        for key in ("atom_contacts", "protein_interface_residues", "dna_interface_residues"):
            if evidence[key] != audit[key]:
                raise ImportValidationError(f"Interface-count mismatch for {accession}: {key}")
        if not full_length_ownership_ok(accession, path):
            raise ImportValidationError(f"Full-length path/accession mismatch: {accession}/{path}")
        selected.append(
            {
                "accession": accession,
                "tf_id": accession,
                "candidate_type": "FULL_LENGTH_AF3",
                "candidate_model_path": model_path,
                "candidate_model_rank": model_rank,
                "domain_or_fragment_id": "FULL_LENGTH",
                "fragment_start": "",
                "fragment_end": "",
                "reviewed_expected_DBD": "NOT_APPLICABLE",
                "interface_result": PASS_RESULT,
                "atom_contacts": audit["atom_contacts"],
                "protein_interface_residues": audit["protein_interface_residues"],
                "dna_interface_residues": audit["dna_interface_residues"],
                "proposed_assignment_role": "OWN_ACCESSION_MODEL",
                "proposed_is_active": "YES",
                "decision": FINAL_FULL_STATUS,
                "source_audit": source_audit,
            }
        )

    adjudication_rows = read_tsv(FINAL_ADJUDICATION_INPUT)
    fragment_rows = [
        row for row in adjudication_rows
        if row["final_structural_status"] == FINAL_FRAGMENT_STATUS
    ]
    if len(fragment_rows) != 1:
        raise ImportValidationError(
            f"Frozen adjudication must contain one {FINAL_FRAGMENT_STATUS} accession."
        )
    final_fragment = fragment_rows[0]
    accession = final_fragment["accession"]
    model_path = final_fragment["final_model_path"]
    path = Path(model_path)
    start, end = fragment_parts(accession, path)
    if (accession, start, end) != ("A0A1B0GWI9", 62, 166):
        raise ImportValidationError(
            f"Final fragment control is not A0A1B0GWI9 residues 62-166: {accession}/{start}-{end}"
        )
    fragment_audit_rows = [
        row for row in read_tsv(FINAL_FRAGMENT_AUDIT_INPUT)
        if row["accession"] == accession and row["model_path"] == model_path
    ]
    if len(fragment_audit_rows) != 1:
        raise ImportValidationError("Final fragment path is not uniquely identified by its audit.")
    fragment_audit = fragment_audit_rows[0]
    if (
        fragment_audit["protein_present"] != "YES"
        or fragment_audit["dna_present"] != "YES"
        or fragment_audit["interface_result"] != PASS_RESULT
        or fragment_audit["error"]
    ):
        raise ImportValidationError("Final NFIB fragment lacks passing audit evidence.")
    selected.append(
        {
            "accession": accession,
            "tf_id": accession,
            "candidate_type": "DBD_FRAGMENT_AF3",
            "candidate_model_path": model_path,
            "candidate_model_rank": "",
            "domain_or_fragment_id": f"{accession}_{start}-{end}",
            "fragment_start": str(start),
            "fragment_end": str(end),
            "reviewed_expected_DBD": "YES",
            "interface_result": PASS_RESULT,
            "atom_contacts": fragment_audit["atom_contacts"],
            "protein_interface_residues": fragment_audit["protein_interface_residues"],
            "dna_interface_residues": fragment_audit["dna_interface_residues"],
            "proposed_assignment_role": "DBD_FRAGMENT_MODEL",
            "proposed_is_active": "YES",
            "decision": FINAL_FRAGMENT_STATUS,
            "source_audit": str(FINAL_ADJUDICATION_INPUT),
        }
    )

    if len(selected) != EXPECTED_FINAL_ACCESSION_COUNT:
        raise ImportValidationError(f"Expected 114 final selections, found {len(selected)}.")
    if len({row["accession"] for row in selected}) != EXPECTED_FINAL_ACCESSION_COUNT:
        raise ImportValidationError("Final accepted AF3 accessions are not unique.")
    if Counter(row["proposed_assignment_role"] for row in selected) != Counter(EXPECTED_ROLE_COUNTS):
        raise ImportValidationError("Final selection role counts do not equal 113 full-length plus one fragment.")

    for row in selected:
        accession = row["accession"]
        path = Path(row["candidate_model_path"])
        if not exact_database_mapping_ok(connection, accession):
            raise ImportValidationError(f"Missing or ambiguous exact database mapping: {accession}")
        if not path.is_absolute() or not readable_file(path):
            raise ImportValidationError(f"Selected coordinate file is missing/unreadable: {path}")
    return selected, {row["accession"] for row in selected}


def adjacent_confidence_path(model_path: Path) -> Path:
    suffix = "_model.cif"
    if model_path.name.endswith(suffix):
        name = model_path.name[: -len(suffix)] + "_confidences.json"
    elif model_path.name == "model.cif":
        name = "confidences.json"
    else:
        raise ImportValidationError(f"Unexpected selected AF3 filename: {model_path}")
    return model_path.with_name(name)


def build_artifacts(rows: list[dict[str, str]], created_at: str) -> list[dict[str, str]]:
    artifacts: list[dict[str, str]] = []
    digest_cache: dict[Path, str] = {}

    def digest(path: Path) -> str:
        if path not in digest_cache:
            digest_cache[path] = sha256_file(path)
        return digest_cache[path]

    for row in rows:
        model_path = Path(row["candidate_model_path"])
        model_kind = (
            "full-length" if row["candidate_type"] == "FULL_LENGTH_AF3" else "DBD-fragment"
        )
        artifacts.append(
            {
                "model_path": str(model_path),
                "artifact_type": "MMCIF",
                "artifact_path": str(model_path),
                "sha256": digest(model_path),
                "notes": f"Selected active AF3 {model_kind} coordinate file.",
                "created_at": created_at,
            }
        )
        confidence_path = adjacent_confidence_path(model_path)
        if not confidence_path.exists():
            continue
        if not readable_file(confidence_path):
            raise ImportValidationError(f"Confidence JSON is unreadable: {confidence_path}")
        try:
            payload = json.loads(confidence_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise ImportValidationError(f"Invalid confidence JSON {confidence_path}: {exc}") from exc
        if not isinstance(payload, dict):
            raise ImportValidationError(f"Confidence JSON is not an object: {confidence_path}")
        if any(key in payload for key in ("atom_plddts", "plddt", "plddts")):
            artifacts.append(
                {
                    "model_path": str(model_path),
                    "artifact_type": "PLDDT_JSON",
                    "artifact_path": str(confidence_path),
                    "sha256": digest(confidence_path),
                    "notes": "AF3 confidence JSON containing per-atom pLDDT values; not yet displayed.",
                    "created_at": created_at,
                }
            )
        if any(key in payload for key in ("pae", "predicted_aligned_error")):
            artifacts.append(
                {
                    "model_path": str(model_path),
                    "artifact_type": "PAE_JSON",
                    "artifact_path": str(confidence_path),
                    "sha256": digest(confidence_path),
                    "notes": "AF3 confidence JSON containing PAE values; not yet displayed.",
                    "created_at": created_at,
                }
            )
    return artifacts


def write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    if path.exists():
        raise ImportValidationError(f"Refusing to overwrite output: {path}")
    if not path.parent.is_dir():
        raise ImportValidationError(f"Output directory does not exist: {path.parent}")
    fieldnames = ["category", "key", "value", "expected", "result", "notes"]
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    try:
        with temporary.open("x", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def write_json(path: Path, payload: dict[str, Any]) -> None:
    if path.exists():
        raise ImportValidationError(f"Refusing to overwrite output: {path}")
    if not path.parent.is_dir():
        raise ImportValidationError(f"Output directory does not exist: {path.parent}")
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    try:
        with temporary.open("x", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def prepare_assignments(
    selected_rows: list[dict[str, str]], created_at: str
) -> list[dict[str, Any]]:
    assignments: list[dict[str, Any]] = []
    for row in selected_rows:
        notes = {
            "candidate_type": row["candidate_type"],
            "domain_or_fragment_id": row["domain_or_fragment_id"],
            "fragment_start": row["fragment_start"],
            "fragment_end": row["fragment_end"],
            "final_decision": row["decision"],
            "interface_atom_contacts": row["atom_contacts"],
            "protein_interface_residues": row["protein_interface_residues"],
            "dna_interface_residues": row["dna_interface_residues"],
            "source_audit": row["source_audit"],
        }
        assignments.append(
            {
                "tf_id": row["tf_id"],
                "display_accession": row["accession"],
                "model_accession": row["accession"],
                "model_path": row["candidate_model_path"],
                "model_role": row["proposed_assignment_role"],
                "interface_status": row["interface_result"],
                "is_active": 1,
                "source": "AF3",
                "notes": json.dumps(notes, sort_keys=True, separators=(",", ":")),
                "created_at": created_at,
                "updated_at": created_at,
            }
        )
    return assignments


def assignment_key(row: dict[str, Any] | sqlite3.Row) -> tuple[str, str, str]:
    return (row["display_accession"], row["model_path"], row["model_role"])


def reconcile_assignments(
    connection: sqlite3.Connection,
    assignments: list[dict[str, Any]],
    updated_at: str,
) -> dict[str, int]:
    expected_keys = {assignment_key(row) for row in assignments}
    changes = {"inserted": 0, "updated": 0, "deactivated": 0}
    active_rows = connection.execute(
        "SELECT * FROM structure_model_assignment WHERE is_active=1 AND UPPER(source)='AF3'"
    ).fetchall()
    for row in active_rows:
        if assignment_key(row) not in expected_keys:
            connection.execute(
                "UPDATE structure_model_assignment SET is_active=0, updated_at=? WHERE id=?",
                (updated_at, row["id"]),
            )
            changes["deactivated"] += 1

    compared_fields = (
        "tf_id", "display_accession", "model_accession", "model_path", "model_role",
        "interface_status", "is_active", "source",
    )
    for desired in assignments:
        existing = connection.execute(
            "SELECT * FROM structure_model_assignment "
            "WHERE display_accession=? AND model_path=? AND model_role=?",
            assignment_key(desired),
        ).fetchone()
        if existing is None:
            connection.execute(
                """
                INSERT INTO structure_model_assignment (
                    tf_id, display_accession, model_accession, model_path, model_role,
                    interface_status, is_active, source, notes, created_at, updated_at
                ) VALUES (
                    :tf_id, :display_accession, :model_accession, :model_path, :model_role,
                    :interface_status, :is_active, :source, :notes, :created_at, :updated_at
                )
                """,
                desired,
            )
            changes["inserted"] += 1
            continue
        if any(existing[field] != desired[field] for field in compared_fields):
            connection.execute(
                """
                UPDATE structure_model_assignment SET
                    tf_id=:tf_id,
                    model_accession=:model_accession,
                    interface_status=:interface_status,
                    is_active=:is_active,
                    source=:source,
                    notes=:notes,
                    updated_at=:updated_at
                WHERE id=:id
                """,
                {**desired, "id": existing["id"]},
            )
            changes["updated"] += 1
    return changes


def reconcile_artifacts(
    connection: sqlite3.Connection, artifacts: list[dict[str, str]]
) -> dict[str, int]:
    changes = {"inserted": 0, "updated": 0}
    for desired in artifacts:
        existing = connection.execute(
            "SELECT * FROM structure_confidence_artifact "
            "WHERE model_path=? AND artifact_type=? AND artifact_path=?",
            (
                desired["model_path"],
                desired["artifact_type"],
                desired["artifact_path"],
            ),
        ).fetchone()
        if existing is None:
            connection.execute(
                """
                INSERT INTO structure_confidence_artifact (
                    model_path, artifact_type, artifact_path, sha256, notes, created_at
                ) VALUES (
                    :model_path, :artifact_type, :artifact_path, :sha256, :notes, :created_at
                )
                """,
                desired,
            )
            changes["inserted"] += 1
        elif existing["sha256"] != desired["sha256"] or existing["notes"] != desired["notes"]:
            connection.execute(
                "UPDATE structure_confidence_artifact SET sha256=?, notes=? WHERE id=?",
                (desired["sha256"], desired["notes"], existing["id"]),
            )
            changes["updated"] += 1
    return changes


def validate_final_state(
    connection: sqlite3.Connection,
    assignments: list[dict[str, Any]],
    artifacts: list[dict[str, str]],
) -> dict[str, Any]:
    expected_keys = {assignment_key(row) for row in assignments}
    active_rows = connection.execute(
        "SELECT * FROM structure_model_assignment WHERE is_active=1 AND UPPER(source)='AF3'"
    ).fetchall()
    active_keys = {assignment_key(row) for row in active_rows}
    role_counts = dict(Counter(row["model_role"] for row in active_rows))
    unexpected = sorted(active_keys - expected_keys)
    missing = sorted(expected_keys - active_keys)
    if len(active_rows) != EXPECTED_FINAL_ACCESSION_COUNT:
        raise ImportValidationError(f"Active AF3 assignment count differs: {len(active_rows)}")
    if len({row["display_accession"] for row in active_rows}) != EXPECTED_FINAL_ACCESSION_COUNT:
        raise ImportValidationError("Active AF3 assignments are not 114 unique accessions.")
    if role_counts != EXPECTED_ROLE_COUNTS:
        raise ImportValidationError(f"Active AF3 role counts differ: {role_counts}")
    if unexpected or missing:
        raise ImportValidationError(
            f"Active AF3 set differs; unexpected={unexpected}, missing={missing}"
        )

    final_by_accession = {row["display_accession"]: row for row in assignments}
    nfib_expected = final_by_accession["A0A1B0GWI9"]
    nfib_rows = [
        row for row in active_rows
        if row["display_accession"] == "A0A1B0GWI9"
        and assignment_key(row) == assignment_key(nfib_expected)
    ]
    if len(nfib_rows) != 1 or nfib_rows[0]["model_role"] != "DBD_FRAGMENT_MODEL":
        raise ImportValidationError("Active NFIB DBD-fragment invariant failed.")

    h7_expected = final_by_accession["H7C4N4"]
    h7_rows = [
        row for row in active_rows
        if row["display_accession"] == "H7C4N4"
        and assignment_key(row) == assignment_key(h7_expected)
    ]
    old_h7_path = str(
        ROOT
        / "external/baldo_model_inventory/fragments_cif_audit/AF3"
        / "H7C4N4_PF03299_2-70/H7C4N4_PF03299_2-70_model.cif"
    )
    old_h7_active = connection.execute(
        "SELECT COUNT(*) FROM structure_model_assignment "
        "WHERE display_accession='H7C4N4' AND model_path=? AND is_active=1",
        (old_h7_path,),
    ).fetchone()[0]
    if (
        len(h7_rows) != 1
        or h7_rows[0]["model_role"] != "OWN_ACCESSION_MODEL"
        or old_h7_active != 0
    ):
        raise ImportValidationError("H7C4N4 full-length replacement invariant failed.")

    selected_paths = {row["model_path"] for row in assignments}
    artifact_counts = dict(
        connection.execute(
            "SELECT artifact_type, COUNT(*) FROM structure_confidence_artifact "
            f"WHERE model_path IN ({','.join('?' for _ in selected_paths)}) "
            "GROUP BY artifact_type",
            tuple(sorted(selected_paths)),
        )
    )
    if artifact_counts != EXPECTED_ARTIFACT_COUNTS:
        raise ImportValidationError(f"Final AF3 artifact counts differ: {artifact_counts}")
    if len(artifacts) != sum(EXPECTED_ARTIFACT_COUNTS.values()):
        raise ImportValidationError("Prepared final AF3 artifacts are not unique and complete.")
    return {
        "active_assignment_count": len(active_rows),
        "role_counts": role_counts,
        "artifact_counts": artifact_counts,
        "unexpected_active_assignments": unexpected,
        "missing_active_assignments": missing,
        "nfib_ready": True,
        "h7_replacement_ready": True,
        "old_h7_fragment_active": old_h7_active,
    }


def open_working_connection(database_path: Path, dry_run: bool) -> sqlite3.Connection:
    if not dry_run:
        return sqlite3.connect(f"file:{database_path}?mode=rw", uri=True)
    source = sqlite3.connect(f"file:{database_path}?mode=ro&immutable=1", uri=True)
    working = sqlite3.connect(":memory:")
    try:
        source.backup(working)
    finally:
        source.close()
    return working


def main() -> int:
    args = parse_args()
    database_path = args.database.resolve()
    qc_output = args.qc_output.resolve()
    summary_output = args.summary_output.resolve()
    connection: sqlite3.Connection | None = None
    try:
        if database_path == SOURCE_DATABASE_LABEL.resolve():
            raise ImportValidationError("Refusing to open the source database for writing.")
        if database_path != DEFAULT_DATABASE.resolve() and Path("/tmp") not in database_path.parents:
            raise ImportValidationError("Target database must be local staging or a test copy under /tmp.")
        if not database_path.is_file():
            raise ImportValidationError(f"Target database is missing: {database_path}")
        if not args.dry_run and (
            qc_output == summary_output or qc_output.exists() or summary_output.exists()
        ):
            raise ImportValidationError("Output paths collide or already exist; refusing overwrite.")
        if args.expected_preimport_sha256 and not re.fullmatch(r"[0-9a-f]{64}", args.expected_preimport_sha256):
            raise ImportValidationError("--expected-preimport-sha256 must be a lowercase SHA256.")

        before = fingerprint(database_path)
        if args.expected_preimport_sha256 and before["sha256"] != args.expected_preimport_sha256:
            raise ImportValidationError(
                f"Pre-import database SHA256 mismatch: {before['sha256']}"
            )
        connection = open_working_connection(database_path, args.dry_run)
        connection.row_factory = sqlite3.Row
        connection.execute("PRAGMA foreign_keys=ON")
        if connection.execute("PRAGMA foreign_keys").fetchone()[0] != 1:
            raise ImportValidationError("Could not enable foreign-key enforcement.")
        for table in ("tf_structure_status", "structure_model_assignment", "structure_confidence_artifact"):
            if connection.execute(
                "SELECT 1 FROM sqlite_master WHERE type='table' AND name=?", (table,)
            ).fetchone() is None:
                raise ImportValidationError(f"Required staging table is missing: {table}")

        counts_before = application_table_counts(connection)
        protected_before = {table: counts_before.get(table) for table in PROTECTED_PWM_TABLES}
        fk_before = foreign_key_signature(connection)
        selected_rows, final_accessions = validate_inputs(connection)
        if len(final_accessions) != EXPECTED_FINAL_ACCESSION_COUNT:
            raise ImportValidationError("Final accepted accession-set size differs from 114.")
        created_at = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
        assignments = prepare_assignments(selected_rows, created_at)
        artifacts = build_artifacts(selected_rows, created_at)
        if Counter(row["model_role"] for row in assignments) != Counter(EXPECTED_ROLE_COUNTS):
            raise ImportValidationError("Prepared assignment role counts do not match 113/1.")
        if Counter(row["artifact_type"] for row in artifacts) != Counter(EXPECTED_ARTIFACT_COUNTS):
            raise ImportValidationError(
                f"Available artifact counts differ from final local evidence: "
                f"{Counter(row['artifact_type'] for row in artifacts)}"
            )

        connection.execute("BEGIN IMMEDIATE")
        try:
            assignment_changes = reconcile_assignments(connection, assignments, created_at)
            artifact_changes = reconcile_artifacts(connection, artifacts)
            final_state = validate_final_state(connection, assignments, artifacts)
            idempotent_second_pass = None
            if args.dry_run:
                second_assignment_changes = reconcile_assignments(
                    connection, assignments, created_at
                )
                second_artifact_changes = reconcile_artifacts(connection, artifacts)
                validate_final_state(connection, assignments, artifacts)
                idempotent_second_pass = {
                    "assignments": second_assignment_changes,
                    "artifacts": second_artifact_changes,
                }
                if any(second_assignment_changes.values()) or any(second_artifact_changes.values()):
                    raise ImportValidationError("A second reconciliation pass was not idempotent.")
            quick_check = str(connection.execute("PRAGMA quick_check").fetchone()[0])
            integrity_check = str(connection.execute("PRAGMA integrity_check").fetchone()[0])
            if quick_check != "ok" or integrity_check != "ok":
                raise ImportValidationError(
                    f"Post-reconciliation database check failed: quick={quick_check}, integrity={integrity_check}"
                )
            fk_after = foreign_key_signature(connection)
            if fk_after != fk_before:
                raise ImportValidationError("Foreign-key check signature changed after reconciliation.")
            counts_after = application_table_counts(connection)
            for table, count in counts_before.items():
                if table not in {"structure_model_assignment", "structure_confidence_artifact"}:
                    if counts_after[table] != count:
                        raise ImportValidationError(f"Unrelated table row count changed: {table}")
            protected_after = {table: counts_after.get(table) for table in PROTECTED_PWM_TABLES}
            pwm_motif_tables_touched = protected_after != protected_before
            if pwm_motif_tables_touched:
                raise ImportValidationError("A protected PWM/motif table would change.")
            if args.dry_run:
                connection.rollback()
            else:
                connection.commit()
        except Exception:
            connection.rollback()
            raise
        finally:
            connection.close()
            connection = None

        after = fingerprint(database_path)
        if args.dry_run and after != before:
            raise ImportValidationError("Dry-run changed the staging database fingerprint.")

        if not args.dry_run:
            summary_rows: list[dict[str, Any]] = []
            for role, expected in EXPECTED_ROLE_COUNTS.items():
                summary_rows.append({
                    "category": "ASSIGNMENT_ROLE", "key": role,
                    "value": final_state["role_counts"][role], "expected": expected,
                    "result": "PASS", "notes": "Active exact-accession final AF3 assignment.",
                })
            for artifact_type, expected in EXPECTED_ARTIFACT_COUNTS.items():
                summary_rows.append({
                    "category": "CONFIDENCE_ARTIFACT", "key": artifact_type,
                    "value": final_state["artifact_counts"][artifact_type], "expected": expected,
                    "result": "PASS", "notes": "Artifact for a final accepted AF3 model.",
                })
            qc = {
                "validation_status": "PASS",
                "database": str(database_path),
                "database_before": before,
                "database_after": after,
                "quick_check": quick_check,
                "integrity_check": integrity_check,
                "selected_input_count": len(selected_rows),
                "active_assignment_count": final_state["active_assignment_count"],
                "active_assignment_role_counts": final_state["role_counts"],
                "confidence_artifact_type_counts": final_state["artifact_counts"],
                "assignment_changes": assignment_changes,
                "artifact_changes": artifact_changes,
                "pwm_motif_tables_touched": pwm_motif_tables_touched,
                "source_database_opened_by_importer": False,
                "remote_connection_attempted": False,
                "model_files_modified": False,
            }
            write_tsv(summary_output, summary_rows)
            write_json(qc_output, qc)

        print("DRY_RUN_VALIDATED" if args.dry_run else "IMPORTED_AND_VALIDATED")
        print(f"target_active_assignments={final_state['active_assignment_count']}")
        print(f"target_role_counts={final_state['role_counts']}")
        print(f"nfib_ready={'YES' if final_state['nfib_ready'] else 'NO'}")
        print(f"h7_replacement_ready={'YES' if final_state['h7_replacement_ready'] else 'NO'}")
        print(f"old_h7_fragment_active={final_state['old_h7_fragment_active']}")
        print(f"unexpected_active_assignments={len(final_state['unexpected_active_assignments'])}")
        print(f"assignment_changes={assignment_changes}")
        print(f"artifact_changes={artifact_changes}")
        if idempotent_second_pass is not None:
            print(f"idempotent_second_pass={idempotent_second_pass}")
        print(f"pwm_motif_tables_touched={'YES' if pwm_motif_tables_touched else 'NO'}")
        print(f"database_before_sha256={before['sha256']}")
        print(f"database_after_sha256={after['sha256']}")
        return 0
    except (ImportValidationError, OSError, csv.Error, sqlite3.Error, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    finally:
        if connection is not None:
            connection.close()


if __name__ == "__main__":
    raise SystemExit(main())
