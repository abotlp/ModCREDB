#!/usr/bin/env python3
"""Build a read-only AF3/fragment model import preview for local staging."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import sqlite3
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable


ROOT = Path(__file__).resolve().parents[1]
STAGING_DATABASE = ROOT / "data" / "tf_webdb_local_staging.sqlite"
SOURCE_DATABASE_LABEL = ROOT / "data" / "tf_webdb.sqlite"
AF3_ROOT = Path("/data/sbi/interchange/boliva/patricia/AF3")
MODCRE_ROOT = Path("/data/sbi/interchange/boliva/patricia/models")

# The two requested generic filenames do not exist locally. These are their
# concrete, provenance-compatible equivalents.
AF3_INVENTORY = ROOT / "external/baldo_model_inventory/af3_unique_tf_inventory.tsv"
AF3_ELIGIBILITY = ROOT / "external/baldo_model_inventory/af3_candidate_db_eligibility.tsv"
AF3_INTERFACE_AUDIT = ROOT / "external/baldo_model_inventory/af3_interface_audit.tsv"
AF3_INTERFACE_SUMMARY = ROOT / "external/baldo_model_inventory/af3_tf_interface_summary.tsv"
FAILED_FULL_LENGTH = ROOT / "external/baldo_model_inventory/af3_failed_tf_for_fragmentation.tsv"
FRAGMENT_SUMMARY = ROOT / "external/baldo_model_inventory/af3_fragment_interface_summary.tsv"
FRAGMENT_ALL_MODELS = ROOT / "external/baldo_model_inventory/af3_fragment_all_model_interface_audit.tsv"
FINAL_ADJUDICATION = ROOT / "outputs/af3_failed_102_final_adjudication_DRAFT.tsv"
DATABASE_DISPLAY = ROOT / "outputs/af3_failed_102_database_display_DRAFT.tsv"
PFAM_ROLES = ROOT / "data_sources/af3_failed_102_pfam_role_reviewed.tsv"

DEFAULT_PREVIEW = ROOT / "outputs/local_af3_model_import_preview.tsv"
DEFAULT_SELECTED = ROOT / "outputs/local_af3_model_import_selected.tsv"
DEFAULT_REJECTED = ROOT / "outputs/local_af3_model_import_rejected.tsv"
DEFAULT_QC = ROOT / "outputs/local_af3_model_import_preview_qc.json"

FIELDS = [
    "accession",
    "tf_id",
    "gene",
    "candidate_type",
    "candidate_model_rank",
    "candidate_model_path",
    "archive_path_if_applicable",
    "archive_member_if_applicable",
    "model_file_exists",
    "model_file_readable",
    "protein_present",
    "dna_present",
    "interface_result",
    "atom_contacts",
    "protein_interface_residues",
    "dna_interface_residues",
    "domain_or_fragment_id",
    "reviewed_expected_DBD",
    "candidate_biological_role",
    "database_structure_record_exists",
    "database_model_summary_exists",
    "currently_active_for_accession",
    "same_accession_active_model_count",
    "proposed_assignment_role",
    "proposed_is_active",
    "decision",
    "decision_reason",
    "manual_review_required",
]

ALLOWED_CANDIDATE_TYPES = {"FULL_LENGTH_AF3", "DBD_FRAGMENT_AF3"}
ALLOWED_ASSIGNMENT_ROLES = {"OWN_ACCESSION_MODEL", "DBD_FRAGMENT_MODEL", "NONE"}
ALLOWED_DECISIONS = {
    "PROPOSE_IMPORT_ACTIVE",
    "ALREADY_INDEXED_ACTIVE",
    "REJECT_EMPTY_INTERFACE",
    "REJECT_NON_DBD_FRAGMENT",
    "REJECT_ACCESSION_MISMATCH",
    "REJECT_DUPLICATE",
    "REJECT_MISSING_FILE",
    "MANUAL_REVIEW",
}
SELECTED_DECISIONS = {"PROPOSE_IMPORT_ACTIVE", "ALREADY_INDEXED_ACTIVE"}
PASS_RESULT = "PASS_NONEMPTY_PROTEIN_DNA_INTERFACE"
EMPTY_RESULT = "FAIL_EMPTY_PROTEIN_DNA_INTERFACE"
FRAGMENT_PATTERN = re.compile(
    r"^(?P<accession>[A-Za-z0-9]+)_(?P<pfam>PF\d{5})_(?P<start>\d+)-(?P<end>\d+)$"
)
VARIANT_PRIORITY = {"base": 0, "dimer": 1, "monomer": 2}


class PreviewValidationError(RuntimeError):
    """Raised when the read-only preview cannot be validated safely."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build a read-only preview of full-length AF3 and reviewed Pfam-fragment "
            "models for the local staging database. No database rows are changed."
        )
    )
    parser.add_argument("--preview-output", type=Path, default=DEFAULT_PREVIEW)
    parser.add_argument("--selected-output", type=Path, default=DEFAULT_SELECTED)
    parser.add_argument("--rejected-output", type=Path, default=DEFAULT_REJECTED)
    parser.add_argument("--qc-output", type=Path, default=DEFAULT_QC)
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


def connect_staging_read_only() -> sqlite3.Connection:
    connection = sqlite3.connect(
        f"file:{STAGING_DATABASE}?mode=ro&immutable=1",
        uri=True,
    )
    connection.row_factory = sqlite3.Row
    return connection


def read_tsv(path: Path) -> list[dict[str, str]]:
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None:
                raise PreviewValidationError(f"TSV has no header: {path}")
            return [dict(row) for row in reader]
    except (OSError, csv.Error) as exc:
        raise PreviewValidationError(f"Could not read {path}: {exc}") from exc


def require_columns(path: Path, rows: list[dict[str, str]], required: set[str]) -> None:
    if not rows:
        raise PreviewValidationError(f"Input has no data rows: {path}")
    missing = sorted(required - set(rows[0]))
    if missing:
        raise PreviewValidationError(f"{path} is missing columns: {missing}")


def split_semicolon(value: str) -> list[str]:
    return [token.strip() for token in value.split(";") if token.strip()]


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


def file_state(path: Path) -> tuple[str, str]:
    if not path.is_file():
        return "NO", "NO"
    try:
        with path.open("rb") as handle:
            handle.read(1)
    except OSError:
        return "YES", "NO"
    return "YES", "YES"


def absolute_model_path(raw_path: str) -> Path:
    path = Path(raw_path)
    return path if path.is_absolute() else (ROOT / path).resolve()


def primary_gene(gene_names: str) -> str:
    return gene_names.strip().split()[0] if gene_names.strip() else ""


def full_length_ownership_ok(accession: str, row: dict[str, str], path: Path) -> bool:
    directory_name = row["directory_name"].strip()
    if path.parent.name.casefold() != directory_name.casefold():
        return False
    normalized = directory_name
    for suffix in ("_dimer", "_monomer"):
        if normalized.casefold().endswith(suffix):
            normalized = normalized[: -len(suffix)]
            break
    if normalized.casefold() != accession.casefold():
        return False
    return path.name.casefold() == f"{directory_name}_model.cif".casefold()


def fragment_parts(fragment_id: str) -> tuple[str, str, int, int]:
    match = FRAGMENT_PATTERN.fullmatch(fragment_id)
    if not match:
        raise PreviewValidationError(f"Could not parse fragment ID: {fragment_id}")
    return (
        match.group("accession").upper(),
        match.group("pfam"),
        int(match.group("start")),
        int(match.group("end")),
    )


def candidate_structure_matches(
    records: list[dict[str, Any]],
    candidate_path: Path,
    candidate_type: str,
    run_name: str,
) -> list[dict[str, Any]]:
    candidate_string = str(candidate_path).casefold()
    candidate_name = candidate_path.name.casefold()
    candidate_stem = candidate_path.stem.casefold()
    matches: list[dict[str, Any]] = []
    for record in records:
        archive_path = str(record["archive_path"] or "")
        member_path = str(record["member_path"] or "")
        model_id = str(record["model_id"] or "")
        direct_match = candidate_string in {
            archive_path.casefold(),
            member_path.casefold(),
        }
        basename_match = candidate_name == Path(member_path).name.casefold()
        model_id_match = candidate_stem == model_id.casefold()
        archive_member_equivalent = False
        if candidate_type == "FULL_LENGTH_AF3" and str(record["source"]).casefold() == "alphafold":
            prefix = f"{run_name.casefold()}.af3."
            archive_member_equivalent = (
                model_id.casefold().startswith(prefix)
                or Path(member_path).name.casefold().startswith(prefix)
            )
        if direct_match or basename_match or model_id_match or archive_member_equivalent:
            matches.append(record)
    return sorted(matches, key=lambda row: (row["status"] != "active", int(row["id"])))


def database_fields(
    accession: str,
    candidate_path: Path,
    candidate_type: str,
    run_name: str,
    structures_by_tf: dict[str, list[dict[str, Any]]],
    summaries_by_tf: dict[str, list[dict[str, Any]]],
) -> dict[str, Any]:
    records = structures_by_tf.get(accession, [])
    matches = candidate_structure_matches(records, candidate_path, candidate_type, run_name)
    match_ids = {int(row["id"]) for row in matches}
    summaries = summaries_by_tf.get(accession, [])
    summary_matches = [
        row
        for row in summaries
        if int(row["summary_file_id"]) in match_ids
        or (row["matched_structure_id"] is not None and int(row["matched_structure_id"]) in match_ids)
        or str(row["summary_model_id"] or "").casefold() == candidate_path.stem.casefold()
    ]
    active_coordinate_records = [
        row
        for row in records
        if row["status"] == "active" and str(row["file_type"]).casefold() in {"pdb", "cif", "mmcif"}
    ]
    first = matches[0] if matches else None
    return {
        "database_structure_record_exists": "YES" if matches else "NO",
        "database_model_summary_exists": "YES" if summary_matches else "NO",
        "currently_active_for_accession": (
            "YES" if any(row["status"] == "active" for row in matches) else "NO"
        ),
        "same_accession_active_model_count": len(active_coordinate_records),
        "archive_path_if_applicable": str(first["archive_path"]) if first else "",
        "archive_member_if_applicable": str(first["member_path"]) if first else "",
    }


def infer_polymer_presence(result: str, protein_value: str = "", dna_value: str = "") -> tuple[str, str]:
    if protein_value or dna_value:
        return ("YES" if protein_value else "NO", "YES" if dna_value else "NO")
    if result in {PASS_RESULT, EMPTY_RESULT}:
        return "YES", "YES"
    if result == "FAIL_NO_PROTEIN":
        return "NO", "YES"
    if result == "FAIL_NO_DNA":
        return "YES", "NO"
    return "UNCERTAIN", "UNCERTAIN"


def base_row(
    accession: str,
    annotation: dict[str, Any],
    candidate_type: str,
    path: Path,
    result: str,
    atom_contacts: str,
    protein_residues: str,
    dna_residues: str,
    protein_present: str,
    dna_present: str,
    domain_or_fragment_id: str,
    expected_dbd: str,
    biological_role: str,
    db_fields: dict[str, Any],
) -> dict[str, Any]:
    exists, readable = file_state(path)
    row: dict[str, Any] = {
        "accession": accession,
        "tf_id": str(annotation["tf_id"]),
        "gene": primary_gene(str(annotation["gene_names"] or "")),
        "candidate_type": candidate_type,
        "candidate_model_rank": "1",
        "candidate_model_path": str(path),
        "archive_path_if_applicable": db_fields["archive_path_if_applicable"],
        "archive_member_if_applicable": db_fields["archive_member_if_applicable"],
        "model_file_exists": exists,
        "model_file_readable": readable,
        "protein_present": protein_present,
        "dna_present": dna_present,
        "interface_result": result,
        "atom_contacts": atom_contacts,
        "protein_interface_residues": protein_residues,
        "dna_interface_residues": dna_residues,
        "domain_or_fragment_id": domain_or_fragment_id,
        "reviewed_expected_DBD": expected_dbd,
        "candidate_biological_role": biological_role,
        "database_structure_record_exists": db_fields["database_structure_record_exists"],
        "database_model_summary_exists": db_fields["database_model_summary_exists"],
        "currently_active_for_accession": db_fields["currently_active_for_accession"],
        "same_accession_active_model_count": db_fields["same_accession_active_model_count"],
        "proposed_assignment_role": "NONE",
        "proposed_is_active": "NO",
        "decision": "",
        "decision_reason": "",
        "manual_review_required": "NO",
    }
    return row


def choose_primary_decision(row: dict[str, Any], role: str) -> tuple[str, str, str]:
    if row["currently_active_for_accession"] == "YES":
        return (
            "ALREADY_INDEXED_ACTIVE",
            role,
            "YES",
        )
    if row["database_structure_record_exists"] == "YES":
        return "REJECT_DUPLICATE", "NONE", "NO"
    if int(row["same_accession_active_model_count"]) > 0:
        return "REJECT_DUPLICATE", "NONE", "NO"
    return "PROPOSE_IMPORT_ACTIVE", role, "YES"


def write_tsv(path: Path, rows: Iterable[dict[str, Any]]) -> None:
    rows = list(rows)
    if path.exists():
        raise PreviewValidationError(f"Refusing to overwrite existing output: {path}")
    if not path.parent.is_dir():
        raise PreviewValidationError(f"Output directory does not exist: {path.parent}")
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    try:
        with temporary.open("x", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
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
        raise PreviewValidationError(f"Refusing to overwrite existing output: {path}")
    if not path.parent.is_dir():
        raise PreviewValidationError(f"Output directory does not exist: {path.parent}")
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


def main() -> int:
    args = parse_args()
    output_paths = [
        args.preview_output.resolve(),
        args.selected_output.resolve(),
        args.rejected_output.resolve(),
        args.qc_output.resolve(),
    ]
    try:
        if len(set(output_paths)) != 4:
            raise PreviewValidationError("The four output paths must be distinct.")
        existing = [str(path) for path in output_paths if path.exists()]
        if existing:
            raise PreviewValidationError(f"Refusing to overwrite outputs: {existing}")
        required_paths = [
            STAGING_DATABASE,
            AF3_ROOT,
            MODCRE_ROOT,
            AF3_INVENTORY,
            AF3_ELIGIBILITY,
            AF3_INTERFACE_AUDIT,
            AF3_INTERFACE_SUMMARY,
            FAILED_FULL_LENGTH,
            FRAGMENT_SUMMARY,
            FRAGMENT_ALL_MODELS,
            FINAL_ADJUDICATION,
            DATABASE_DISPLAY,
            PFAM_ROLES,
        ]
        missing = [str(path) for path in required_paths if not path.exists()]
        if missing:
            raise PreviewValidationError(f"Missing required local inputs: {missing}")

        staging_before = fingerprint(STAGING_DATABASE)
        inventory_rows = read_tsv(AF3_INVENTORY)
        eligibility_rows = read_tsv(AF3_ELIGIBILITY)
        interface_rows = read_tsv(AF3_INTERFACE_AUDIT)
        interface_summary_rows = read_tsv(AF3_INTERFACE_SUMMARY)
        failed_rows = read_tsv(FAILED_FULL_LENGTH)
        fragment_summary_rows = read_tsv(FRAGMENT_SUMMARY)
        fragment_model_rows = read_tsv(FRAGMENT_ALL_MODELS)
        adjudication_rows = read_tsv(FINAL_ADJUDICATION)
        display_rows = read_tsv(DATABASE_DISPLAY)
        pfam_rows = read_tsv(PFAM_ROLES)

        require_columns(AF3_INVENTORY, inventory_rows, {"tf_id", "db_class", "selected_cifs"})
        require_columns(
            AF3_ELIGIBILITY,
            eligibility_rows,
            {"tf_id", "eligible_for_af3_import", "selected_cifs", "all_structure_count"},
        )
        require_columns(
            AF3_INTERFACE_AUDIT,
            interface_rows,
            {
                "tf_id", "directory_name", "variant", "model_cif", "result",
                "protein_chains", "nucleic_chains", "atom_contacts",
                "protein_interface_residues", "nucleic_interface_residues",
            },
        )
        require_columns(AF3_INTERFACE_SUMMARY, interface_summary_rows, {"tf_id", "n_pass", "n_fail"})
        require_columns(FAILED_FULL_LENGTH, failed_rows, {"tf_id", "results"})
        require_columns(
            FRAGMENT_SUMMARY,
            fragment_summary_rows,
            {"fragment_id", "selected_result", "classification", "best_cif_path"},
        )
        require_columns(
            FRAGMENT_ALL_MODELS,
            fragment_model_rows,
            {
                "fragment_id", "model_kind", "sample_name", "cif_path", "result",
                "atom_contacts", "protein_interface_residues", "dna_interface_residues",
            },
        )
        require_columns(FINAL_ADJUDICATION, adjudication_rows, {"failed_accession", "primary_structural_status"})
        require_columns(DATABASE_DISPLAY, display_rows, {"failed_accession", "draft_only"})
        require_columns(
            PFAM_ROLES,
            pfam_rows,
            {"pfam_id", "is_expected_TF_DBD", "proposed_domain_role", "review_status"},
        )

        inventory_by_tf = {row["tf_id"].upper(): row for row in inventory_rows}
        if len(inventory_by_tf) != len(inventory_rows):
            raise PreviewValidationError("AF3 inventory contains duplicate accessions.")
        eligibility_by_tf = {row["tf_id"].upper(): row for row in eligibility_rows}
        if len(eligibility_by_tf) != len(eligibility_rows):
            raise PreviewValidationError("Eligibility TSV contains duplicate accessions.")
        inventory_absent_ids = {
            accession
            for accession, row in inventory_by_tf.items()
            if row["db_class"] == "absent_from_alphafold_db"
        }
        inventory_existing_ids = {
            accession
            for accession, row in inventory_by_tf.items()
            if row["db_class"] == "existing_active"
        }
        if set(eligibility_by_tf) != inventory_absent_ids:
            raise PreviewValidationError(
                "Eligibility TSV is not exactly the absent-from-AlphaFold subset of the inventory."
            )
        eligible_rows = [row for row in eligibility_rows if row["eligible_for_af3_import"].lower() == "yes"]
        eligible_ids = {row["tf_id"].upper() for row in eligible_rows}
        audit_ids = {row["tf_id"].upper() for row in interface_rows}
        if eligible_ids != audit_ids:
            raise PreviewValidationError(
                "Eligible/accession audit reconciliation failed: "
                f"missing_audit={sorted(eligible_ids-audit_ids)}, unexpected_audit={sorted(audit_ids-eligible_ids)}"
            )
        expected_models = {
            (row["tf_id"].upper(), str(absolute_model_path(path)))
            for row in eligible_rows
            for path in split_semicolon(row["selected_cifs"])
        }
        audited_models = {
            (row["tf_id"].upper(), str(absolute_model_path(row["model_cif"])))
            for row in interface_rows
        }
        if len(audited_models) != len(interface_rows) or expected_models != audited_models:
            raise PreviewValidationError("Eligible selected-CIF and interface-audit model sets differ.")
        interface_summary_by_tf = {
            row["tf_id"].upper(): row for row in interface_summary_rows
        }
        if len(interface_summary_by_tf) != len(interface_summary_rows) or set(interface_summary_by_tf) != audit_ids:
            raise PreviewValidationError("Full-length interface summary accession set differs from the detailed audit.")
        detailed_by_tf: dict[str, list[dict[str, str]]] = defaultdict(list)
        for row in interface_rows:
            detailed_by_tf[row["tf_id"].upper()].append(row)
        for accession, detailed in detailed_by_tf.items():
            summary = interface_summary_by_tf[accession]
            passing_count = sum(row["result"] == PASS_RESULT for row in detailed)
            failing_count = len(detailed) - passing_count
            if (
                int(summary["n_variants"]) != len(detailed)
                or int(summary["n_pass"]) != passing_count
                or int(summary["n_fail"]) != failing_count
            ):
                raise PreviewValidationError(f"Full-length interface summary count mismatch: {accession}")

        failed_ids = {row["failed_accession"].upper() for row in adjudication_rows}
        if len(failed_ids) != 102 or failed_ids != {row["tf_id"].upper() for row in failed_rows}:
            raise PreviewValidationError("The two failed-full-length evidence sets do not reconcile to 102 accessions.")
        if {row["failed_accession"].upper() for row in display_rows} != failed_ids:
            raise PreviewValidationError("Database-display DRAFT does not match the 102 failed accessions.")
        if any(row["draft_only"] != "YES" for row in display_rows):
            raise PreviewValidationError("Database-display evidence is not uniformly marked DRAFT-only.")

        pfam_by_id = {row["pfam_id"]: row for row in pfam_rows}
        if len(pfam_by_id) != len(pfam_rows):
            raise PreviewValidationError("Reviewed Pfam role table contains duplicate Pfam IDs.")
        if any(row["review_status"] != "REVIEWED_FOR_AF3_102_AUDIT" for row in pfam_rows):
            raise PreviewValidationError("Pfam role input contains a non-reviewed row.")

        selected_fragment_audits: dict[str, dict[str, str]] = {}
        all_fragment_audits: dict[str, list[dict[str, str]]] = defaultdict(list)
        for row in fragment_model_rows:
            all_fragment_audits[row["fragment_id"]].append(row)
            if row["model_kind"] == "selected":
                if row["fragment_id"] in selected_fragment_audits:
                    raise PreviewValidationError(f"Multiple selected fragment models: {row['fragment_id']}")
                selected_fragment_audits[row["fragment_id"]] = row
        fragment_summary_by_id = {row["fragment_id"]: row for row in fragment_summary_rows}
        if set(selected_fragment_audits) != set(fragment_summary_by_id):
            raise PreviewValidationError("Fragment selected-model and summary ID sets differ.")
        if any(len(rows) != 6 for rows in all_fragment_audits.values()):
            raise PreviewValidationError("Every fragment must retain exactly six audited model variants.")
        for fragment_id, summary in fragment_summary_by_id.items():
            selected_audit = selected_fragment_audits[fragment_id]
            if selected_audit["result"] != summary["selected_result"]:
                raise PreviewValidationError(f"Fragment selected-result mismatch: {fragment_id}")

        with connect_staging_read_only() as database:
            quick_check_before = str(database.execute("PRAGMA quick_check").fetchone()[0])
            if quick_check_before != "ok":
                raise PreviewValidationError(f"Staging quick_check failed: {quick_check_before}")
            table_counts_before = application_table_counts(database)
            annotations = {
                row["tf_id"].upper(): dict(row)
                for row in database.execute(
                    "SELECT tf.tf_id, tf.active_model_count, ta.uniprot_accession, "
                    "ta.gene_names, ta.protein_name FROM tf "
                    "JOIN tf_annotation AS ta ON ta.tf_id=tf.tf_id"
                )
            }
            structures_by_tf: dict[str, list[dict[str, Any]]] = defaultdict(list)
            for row in database.execute("SELECT * FROM structure_file"):
                if row["tf_id"]:
                    structures_by_tf[str(row["tf_id"]).upper()].append(dict(row))
            summaries_by_tf: dict[str, list[dict[str, Any]]] = defaultdict(list)
            for row in database.execute("SELECT * FROM model_summary"):
                if row["tf_id"]:
                    summaries_by_tf[str(row["tf_id"]).upper()].append(dict(row))
            active_alphafold_records = int(database.execute(
                "SELECT COUNT(*) FROM structure_file WHERE source='alphafold' AND status='active'"
            ).fetchone()[0])
            active_alphafold_accessions = int(database.execute(
                "SELECT COUNT(DISTINCT tf_id) FROM structure_file "
                "WHERE source='alphafold' AND status='active'"
            ).fetchone()[0])
            active_alphafold_ids = {
                str(row[0]).upper()
                for row in database.execute(
                    "SELECT DISTINCT tf_id FROM structure_file "
                    "WHERE source='alphafold' AND status='active'"
                )
            }

        if active_alphafold_ids != inventory_existing_ids:
            raise PreviewValidationError(
                "Inventory existing-active accessions differ from active AlphaFold accessions in staging."
            )

        for accession in eligible_ids:
            annotation = annotations.get(accession)
            if annotation is None or annotation["tf_id"] != accession or annotation["uniprot_accession"] != accession:
                raise PreviewValidationError(f"Missing or ambiguous exact database mapping: {accession}")
        fragment_accessions = {fragment_parts(row["fragment_id"])[0] for row in fragment_summary_rows}
        for accession in fragment_accessions:
            annotation = annotations.get(accession)
            if annotation is None or annotation["tf_id"] != accession or annotation["uniprot_accession"] != accession:
                raise PreviewValidationError(f"Missing exact fragment accession mapping: {accession}")

        preview_rows: list[dict[str, Any]] = []
        audit_by_tf: dict[str, list[dict[str, str]]] = defaultdict(list)
        for row in interface_rows:
            audit_by_tf[row["tf_id"].upper()].append(row)

        full_length_valid_accessions: set[str] = set()
        full_length_empty_only_accessions: set[str] = set()
        full_length_base_fail_alternate_pass: set[str] = set()
        for accession in sorted(audit_by_tf):
            group = audit_by_tf[accession]
            passing = [row for row in group if row["result"] == PASS_RESULT]
            if passing:
                full_length_valid_accessions.add(accession)
            elif all(row["result"] == EMPTY_RESULT for row in group):
                full_length_empty_only_accessions.add(accession)
            base = [row for row in group if row["variant"] == "base"]
            if len(base) != 1:
                raise PreviewValidationError(f"Expected one base full-length variant: {accession}")
            if base[0]["result"] != PASS_RESULT and passing:
                full_length_base_fail_alternate_pass.add(accession)
            preferred = min(
                passing,
                key=lambda row: (VARIANT_PRIORITY.get(row["variant"], 99), row["model_cif"]),
            ) if passing else None

            for audit in sorted(group, key=lambda row: (VARIANT_PRIORITY.get(row["variant"], 99), row["model_cif"])):
                path = absolute_model_path(audit["model_cif"])
                ownership_ok = full_length_ownership_ok(accession, audit, path)
                db_values = database_fields(
                    accession,
                    path,
                    "FULL_LENGTH_AF3",
                    audit["directory_name"],
                    structures_by_tf,
                    summaries_by_tf,
                )
                protein_present, dna_present = infer_polymer_presence(
                    audit["result"], audit["protein_chains"], audit["nucleic_chains"]
                )
                row = base_row(
                    accession,
                    annotations[accession],
                    "FULL_LENGTH_AF3",
                    path,
                    audit["result"],
                    audit["atom_contacts"],
                    audit["protein_interface_residues"],
                    audit["nucleic_interface_residues"],
                    protein_present,
                    dna_present,
                    f"FULL_LENGTH_VARIANT:{audit['variant']}",
                    "NOT_APPLICABLE",
                    "EXACT_ACCESSION_FULL_LENGTH_TF_MODEL",
                    db_values,
                )
                variant_label = audit["variant"]
                if not ownership_ok:
                    row["decision"] = "REJECT_ACCESSION_MISMATCH"
                    row["decision_reason"] = (
                        f"The {variant_label} audit path/directory does not encode exact accession {accession}."
                    )
                elif row["model_file_exists"] != "YES" or row["model_file_readable"] != "YES":
                    row["decision"] = "REJECT_MISSING_FILE"
                    row["decision_reason"] = (
                        f"The exact-accession {variant_label} selected/rank-1 coordinate file is missing or unreadable."
                    )
                elif audit["result"] == EMPTY_RESULT:
                    row["decision"] = "REJECT_EMPTY_INTERFACE"
                    failed_note = (
                        " It is one of the 102 failed full-length accessions and is not proposed for activation."
                        if accession in failed_ids else ""
                    )
                    row["decision_reason"] = (
                        f"The audited exact-accession {variant_label} selected/rank-1 full-length model "
                        "contains protein and DNA but has an empty protein-DNA interface."
                        + failed_note
                    )
                elif audit["result"] != PASS_RESULT:
                    row["decision"] = "MANUAL_REVIEW"
                    row["manual_review_required"] = "YES"
                    row["decision_reason"] = (
                        f"The {variant_label} model audit result is {audit['result']}; no automatic selection is safe."
                    )
                elif audit is not preferred:
                    row["decision"] = "REJECT_DUPLICATE"
                    row["decision_reason"] = (
                        f"This exact-accession {variant_label} selected/rank-1 variant has a nonempty interface, "
                        "but the passing base variant is the single preferred primary model. It is retained only "
                        "as non-primary evidence."
                    )
                else:
                    decision, role, active = choose_primary_decision(row, "OWN_ACCESSION_MODEL")
                    row["decision"] = decision
                    row["proposed_assignment_role"] = role
                    row["proposed_is_active"] = active
                    if decision == "PROPOSE_IMPORT_ACTIVE":
                        row["decision_reason"] = (
                            f"The exact-accession {variant_label} selected/rank-1 AF3 model is readable, contains "
                            f"protein and DNA, and has a nonempty audited interface ({audit['atom_contacts']} atom "
                            "contacts). The accession is database-eligible and has no active coordinate model."
                        )
                    elif decision == "ALREADY_INDEXED_ACTIVE":
                        row["decision_reason"] = (
                            "The preferred passing exact-accession model is already represented and active; no new import is needed."
                        )
                    else:
                        row["decision_reason"] = (
                            "The preferred passing model would duplicate an existing structure record or active accession model."
                        )
                preview_rows.append(row)

        non_dbd_contacting_fragments_rejected = 0
        reviewed_valid_fragment_models = 0
        lower_ranked_fragment_rescues_rejected = 0
        for fragment_id in sorted(fragment_summary_by_id):
            summary = fragment_summary_by_id[fragment_id]
            audit = selected_fragment_audits[fragment_id]
            accession, pfam_id, _, _ = fragment_parts(fragment_id)
            role = pfam_by_id.get(pfam_id)
            if role is None:
                raise PreviewValidationError(f"No reviewed Pfam role for fragment {fragment_id}")
            expected_dbd = role["is_expected_TF_DBD"]
            if expected_dbd not in {"YES", "NO", "UNCERTAIN"}:
                raise PreviewValidationError(f"Invalid reviewed DBD value for {pfam_id}: {expected_dbd}")
            path = absolute_model_path(audit["cif_path"])
            path_fragment = path.parent.name
            ownership_ok = path_fragment == fragment_id and fragment_parts(path_fragment)[0] == accession
            db_values = database_fields(
                accession,
                path,
                "DBD_FRAGMENT_AF3",
                fragment_id,
                structures_by_tf,
                summaries_by_tf,
            )
            protein_present, dna_present = infer_polymer_presence(audit["result"])
            row = base_row(
                accession,
                annotations[accession],
                "DBD_FRAGMENT_AF3",
                path,
                audit["result"],
                audit["atom_contacts"],
                audit["protein_interface_residues"],
                audit["dna_interface_residues"],
                protein_present,
                dna_present,
                fragment_id,
                expected_dbd,
                role["proposed_domain_role"],
                db_values,
            )
            passing_models = [
                variant for variant in all_fragment_audits[fragment_id] if variant["result"] == PASS_RESULT
            ]
            selected_passes = audit["result"] == PASS_RESULT
            lower_passing = [variant for variant in passing_models if variant["model_kind"] != "selected"]

            if not ownership_ok:
                row["decision"] = "REJECT_ACCESSION_MISMATCH"
                row["decision_reason"] = (
                    f"Fragment path does not encode the exact accession/fragment identity {fragment_id}."
                )
            elif row["model_file_exists"] != "YES" or row["model_file_readable"] != "YES":
                row["decision"] = "REJECT_MISSING_FILE"
                row["decision_reason"] = "The selected fragment coordinate file is missing or unreadable."
            elif passing_models and expected_dbd == "NO":
                row["decision"] = "REJECT_NON_DBD_FRAGMENT"
                non_dbd_contacting_fragments_rejected += 1
                if selected_passes:
                    contact_note = (
                        f"the selected/rank-1 fragment contacts DNA ({audit['atom_contacts']} atom contacts)"
                    )
                else:
                    best = max(passing_models, key=lambda variant: int(variant["atom_contacts"] or 0))
                    contact_note = (
                        f"the selected/rank-1 fragment has an empty interface but lower-ranked {best['sample_name']} "
                        f"contacts DNA ({best['atom_contacts']} atom contacts)"
                    )
                    lower_ranked_fragment_rescues_rejected += 1
                row["decision_reason"] = (
                    f"Reviewed {pfam_id} evidence says expected TF DBD=NO; although {contact_note}, "
                    "interface contacts alone do not make this a DBD. No fragment variant is proposed."
                )
            elif selected_passes and expected_dbd == "YES":
                decision, assignment_role, active = choose_primary_decision(row, "DBD_FRAGMENT_MODEL")
                row["decision"] = decision
                row["proposed_assignment_role"] = assignment_role
                row["proposed_is_active"] = active
                reviewed_valid_fragment_models += 1
                if decision == "PROPOSE_IMPORT_ACTIVE":
                    row["decision_reason"] = (
                        f"Reviewed {pfam_id} evidence identifies the expected TF DBD; the exact-accession "
                        f"selected/rank-1 fragment is readable and has a nonempty interface "
                        f"({audit['atom_contacts']} atom contacts)."
                    )
                elif decision == "ALREADY_INDEXED_ACTIVE":
                    row["decision_reason"] = (
                        "The reviewed DBD fragment is already represented and active; no new import is needed."
                    )
                else:
                    row["decision_reason"] = (
                        "The reviewed DBD fragment passes but would duplicate an existing structure or active model."
                    )
            elif passing_models and expected_dbd == "UNCERTAIN":
                row["decision"] = "MANUAL_REVIEW"
                row["manual_review_required"] = "YES"
                row["decision_reason"] = (
                    f"At least one {pfam_id} fragment variant contacts DNA, but reviewed expected-DBD status is UNCERTAIN."
                )
            elif audit["result"] == EMPTY_RESULT:
                row["decision"] = "REJECT_EMPTY_INTERFACE"
                row["decision_reason"] = (
                    f"The exact-accession selected/rank-1 {pfam_id} fragment contains protein and DNA but has "
                    f"an empty interface; reviewed expected TF DBD={expected_dbd}."
                )
            else:
                row["decision"] = "MANUAL_REVIEW"
                row["manual_review_required"] = "YES"
                row["decision_reason"] = (
                    f"Selected fragment audit result {audit['result']} cannot be resolved automatically."
                )
            preview_rows.append(row)

        preview_rows.sort(
            key=lambda row: (
                row["accession"], row["candidate_type"], row["domain_or_fragment_id"], row["candidate_model_path"]
            )
        )
        selected_rows = [row for row in preview_rows if row["decision"] in SELECTED_DECISIONS]
        rejected_rows = [row for row in preview_rows if row["decision"] not in SELECTED_DECISIONS]

        invalid_types = sorted({row["candidate_type"] for row in preview_rows} - ALLOWED_CANDIDATE_TYPES)
        invalid_roles = sorted({row["proposed_assignment_role"] for row in preview_rows} - ALLOWED_ASSIGNMENT_ROLES)
        invalid_decisions = sorted({row["decision"] for row in preview_rows} - ALLOWED_DECISIONS)
        if invalid_types or invalid_roles or invalid_decisions:
            raise PreviewValidationError(
                f"Invalid controlled values: types={invalid_types}, roles={invalid_roles}, decisions={invalid_decisions}"
            )
        selected_accessions = [row["accession"] for row in selected_rows]
        if len(selected_accessions) != len(set(selected_accessions)):
            duplicates = sorted(
                accession for accession, count in Counter(selected_accessions).items() if count > 1
            )
            raise PreviewValidationError(f"Selected accessions are not unique: {duplicates}")
        for row in selected_rows:
            if row["model_file_exists"] != "YES" or row["model_file_readable"] != "YES":
                raise PreviewValidationError(f"Selected model is missing/unreadable: {row['candidate_model_path']}")
            if row["protein_present"] != "YES" or row["dna_present"] != "YES":
                raise PreviewValidationError(f"Selected model lacks protein or DNA: {row['candidate_model_path']}")
            if row["interface_result"] != PASS_RESULT or int(row["atom_contacts"] or 0) <= 0:
                raise PreviewValidationError(f"Selected model lacks a nonempty interface: {row['candidate_model_path']}")
            if row["tf_id"] != row["accession"]:
                raise PreviewValidationError(f"Selected accession mismatch: {row['accession']}")
            if row["candidate_type"] == "DBD_FRAGMENT_AF3" and row["reviewed_expected_DBD"] != "YES":
                raise PreviewValidationError(f"Selected fragment is not a reviewed expected DBD: {row['accession']}")
        selected_failed_full = [
            row for row in selected_rows
            if row["candidate_type"] == "FULL_LENGTH_AF3" and row["accession"] in failed_ids
        ]
        if selected_failed_full:
            raise PreviewValidationError("A failed full-length accession was proposed for activation.")
        if full_length_valid_accessions & failed_ids:
            raise PreviewValidationError("One of the 102 failed accessions unexpectedly has a passing full-length audit.")
        if full_length_empty_only_accessions != failed_ids:
            raise PreviewValidationError("The 102 failed set is not exactly the full-length empty-interface set.")

        controls: dict[str, Any] = {}
        e9_selected = [row for row in selected_rows if row["accession"] == "E9PN75"]
        e9_foreign = [
            row for row in preview_rows
            if row["accession"] == "E9PN75" and "Q9NZC4" in (
                row["candidate_model_path"] + row["archive_member_if_applicable"]
            )
        ]
        if e9_selected or e9_foreign:
            raise PreviewValidationError("E9PN75 control failed: own/foreign model proposal detected.")
        controls["E9PN75"] = {
            "result": "PASS",
            "selected_model_count": 0,
            "Q9NZC4_assigned_as_own_model": False,
        }

        h7 = [
            row for row in preview_rows
            if row["accession"] == "H7C4N4" and row["domain_or_fragment_id"] == "H7C4N4_PF03299_2-70"
        ]
        if len(h7) != 1 or h7[0]["decision"] not in SELECTED_DECISIONS or h7[0]["reviewed_expected_DBD"] != "YES":
            raise PreviewValidationError("H7C4N4/PF03299 positive control failed.")
        controls["H7C4N4"] = {
            "result": "PASS",
            "fragment_id": h7[0]["domain_or_fragment_id"],
            "decision": h7[0]["decision"],
        }

        for accession, fragment_id in (
            ("B4DHE0", "B4DHE0_PF01017_146-183"),
            ("Q59EF3", "Q59EF3_PF17725_160-368"),
        ):
            matches = [
                row for row in preview_rows
                if row["accession"] == accession and row["domain_or_fragment_id"] == fragment_id
            ]
            if len(matches) != 1 or matches[0]["decision"] != "REJECT_NON_DBD_FRAGMENT":
                raise PreviewValidationError(f"{accession} negative fragment control failed.")
            controls[accession] = {
                "result": "PASS",
                "fragment_id": fragment_id,
                "decision": matches[0]["decision"],
            }

        decision_counts = dict(sorted(Counter(row["decision"] for row in preview_rows).items()))
        candidate_type_counts = dict(sorted(Counter(row["candidate_type"] for row in preview_rows).items()))
        manual_review_count = decision_counts.get("MANUAL_REVIEW", 0)
        if any(row["manual_review_required"] != ("YES" if row["decision"] == "MANUAL_REVIEW" else "NO") for row in preview_rows):
            raise PreviewValidationError("manual_review_required does not agree with decision.")
        if len(selected_rows) + len(rejected_rows) != len(preview_rows):
            raise PreviewValidationError("Selected/rejected reconciliation failed.")

        with connect_staging_read_only() as database:
            quick_check_after = str(database.execute("PRAGMA quick_check").fetchone()[0])
            table_counts_after = application_table_counts(database)
        staging_after = fingerprint(STAGING_DATABASE)
        if staging_after != staging_before:
            raise PreviewValidationError("Staging database fingerprint changed during read-only preview.")
        if table_counts_after != table_counts_before:
            raise PreviewValidationError("Staging database row counts changed during read-only preview.")
        if quick_check_after != "ok":
            raise PreviewValidationError(f"Final staging quick_check failed: {quick_check_after}")

        qc = {
            "validation_status": "PASS",
            "read_only": True,
            "source_database_opened_by_script": False,
            "staging_database_open_mode": "mode=ro&immutable=1",
            "staging_database_before": staging_before,
            "staging_database_after": staging_after,
            "staging_database_unchanged": staging_before == staging_after,
            "staging_quick_check_before": quick_check_before,
            "staging_quick_check_after": quick_check_after,
            "staging_table_counts_before": table_counts_before,
            "staging_table_counts_after": table_counts_after,
            "input_path_resolution": {
                "requested_af3_model_inventory.tsv": {
                    "status": "LOCAL_FILENAME_DIFFERS",
                    "resolved_equivalent": str(AF3_INVENTORY),
                },
                "requested_af3_model_audit.tsv": {
                    "status": "LOCAL_FILENAME_DIFFERS",
                    "resolved_equivalent": str(AF3_INTERFACE_AUDIT),
                },
            },
            "inputs": {
                "af3_inventory": str(AF3_INVENTORY),
                "af3_eligibility": str(AF3_ELIGIBILITY),
                "af3_interface_audit": str(AF3_INTERFACE_AUDIT),
                "af3_interface_summary": str(AF3_INTERFACE_SUMMARY),
                "failed_full_length": str(FAILED_FULL_LENGTH),
                "fragment_summary": str(FRAGMENT_SUMMARY),
                "fragment_all_models": str(FRAGMENT_ALL_MODELS),
                "final_adjudication": str(FINAL_ADJUDICATION),
                "database_display": str(DATABASE_DISPLAY),
                "pfam_roles": str(PFAM_ROLES),
            },
            "model_roots": {
                "full_length_af3": str(AF3_ROOT),
                "modcre": str(MODCRE_ROOT),
                "fragment_audit_coordinates": str(
                    ROOT / "external/baldo_model_inventory/fragments_cif_audit/AF3"
                ),
            },
            "coordinate_composition_evidence": "Existing local interface audits; no new interface calculation was run.",
            "archive_member_equivalence_checked": True,
            "total_af3_accessions_inventoried": len(inventory_rows),
            "inventory_db_class_counts": dict(sorted(Counter(row["db_class"] for row in inventory_rows).items())),
            "total_absent_from_existing_alphafold_class": len(eligibility_rows),
            "total_database_eligible_af3_accessions": len(eligible_ids),
            "database_ineligible_unmapped_af3_accessions": len(eligibility_rows) - len(eligible_ids),
            "full_length_candidate_model_count": len(interface_rows),
            "full_length_accessions_with_valid_model": len(full_length_valid_accessions),
            "full_length_accessions_only_empty_interface": len(full_length_empty_only_accessions),
            "full_length_base_fail_alternate_pass_accessions": sorted(full_length_base_fail_alternate_pass),
            "failed_102_active_full_length_proposal_count": len(selected_failed_full),
            "fragment_candidate_count": len(fragment_summary_rows),
            "valid_reviewed_dbd_fragment_model_count": reviewed_valid_fragment_models,
            "non_dbd_contacting_fragments_rejected": non_dbd_contacting_fragments_rejected,
            "lower_ranked_non_dbd_fragment_rescues_rejected": lower_ranked_fragment_rescues_rejected,
            "preview_row_count": len(preview_rows),
            "selected_row_count": len(selected_rows),
            "rejected_row_count": len(rejected_rows),
            "selected_unique_accession_count": len(set(selected_accessions)),
            "decision_counts": decision_counts,
            "candidate_type_counts": candidate_type_counts,
            "manual_review_count": manual_review_count,
            "candidate_models_already_indexed_active": decision_counts.get("ALREADY_INDEXED_ACTIVE", 0),
            "models_proposed_for_new_import": decision_counts.get("PROPOSE_IMPORT_ACTIVE", 0),
            "inventory_accessions_existing_active_need_no_new_activation": active_alphafold_accessions,
            "database_active_alphafold_structure_records": active_alphafold_records,
            "selected_model_files_all_exist_and_readable": True,
            "selected_accessions_unique": True,
            "selected_accession_mismatch_count": 0,
            "selected_nonreviewed_fragment_count": 0,
            "controls": controls,
            "created_files": [str(path) for path in output_paths],
            "database_rows_changed": False,
            "website_code_modified": False,
            "model_files_modified": False,
            "remote_connection_attempted": False,
            "sql_created": False,
            "structure_model_assignment_modified": False,
            "source_database_label_only_not_opened": str(SOURCE_DATABASE_LABEL),
        }

        write_tsv(output_paths[0], preview_rows)
        write_tsv(output_paths[1], selected_rows)
        write_tsv(output_paths[2], rejected_rows)
        write_json(output_paths[3], qc)

        print("VALIDATED_READ_ONLY_PREVIEW")
        print(f"preview_rows={len(preview_rows)}")
        print(f"selected_rows={len(selected_rows)}")
        print(f"rejected_rows={len(rejected_rows)}")
        print(f"decision_counts={decision_counts}")
        print(f"candidate_type_counts={candidate_type_counts}")
        return 0
    except (PreviewValidationError, OSError, csv.Error, sqlite3.Error, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
