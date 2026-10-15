#!/usr/bin/env python3
"""Reconcile the read-only pre-Baldo ModCREDB structure population.

The source and staging SQLite databases are opened with mode=ro and immutable=1.
Only the requested TSV/JSON reports are written, and existing outputs are refused.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import socket
import sqlite3
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_FASTA = Path("/home/patricia/TF_database_Baldo_data/TF_without_model.fasta")
MODEL_ROOT = Path("/data/sbi/interchange/boliva/patricia/models")
AF3_ROOT = Path("/data/sbi/interchange/boliva/patricia/AF3")

MASTER_COLUMNS = [
    "tf_id", "accession", "gene", "protein_name", "PWM_annotation", "PWM_count",
    "PWM_source_or_provenance", "original_active_modcre_model",
    "original_active_modcre_model_count", "original_model_path", "sent_to_baldo",
    "baldo_fasta_path", "baldo_full_length_AF3_found", "baldo_valid_interface_found",
    "baldo_empty_interface_only", "baldo_fragment_found", "expected_DBD",
    "previous_pfam_status", "final_DBD_interpretation", "canonical_accession",
    "reason_no_original_modcre_model", "reason_not_sent_to_baldo", "notes",
]

PFAM_REVIEW_COLUMNS = [
    "tf_id", "accession", "gene", "protein_name", "PWM_annotation", "PWM_count",
    "PWM_source_or_provenance", "reviewed_status", "sequence_length",
    "canonical_accession", "expected_DBD", "previous_pfam_status",
    "current_sqlite_pfam_status", "DBD_coverage_percent", "DBD_identity_percent",
    "failed_coordinates_corresponding_to_DBD", "final_DBD_interpretation",
    "evidence_interpretation", "remaining_uncertainty",
]

RETURN_COLUMNS = [
    "accession", "sent_to_baldo", "full_length_AF3_found", "full_length_model_paths",
    "full_length_variant_count", "valid_full_length_interface_found",
    "empty_interface_only", "fragment_found", "fragment_ids", "fragment_model_paths",
    "valid_DBD_fragment_rescue", "no_returned_model_found", "notes",
]

OUTPUT_NAMES = {
    "master": "pre_baldo_modcre_structure_reconciliation.tsv",
    "no_model": "pre_baldo_no_modcre_model.tsv",
    "sent": "pre_baldo_sent_to_baldo.tsv",
    "not_sent": "pre_baldo_no_model_not_sent.tsv",
    "pfam": "pre_baldo_pwm_without_pfam_review.tsv",
    "returns": "pre_baldo_baldo_return_reconciliation.tsv",
    "qc": "pre_baldo_structure_reconciliation_qc.json",
}

ALLOWED_NOT_SENT_REASONS = {
    "EXTERNAL_PWM_ONLY_NOT_STRUCTURALLY_MODELLED",
    "MISSING_OR_INVALID_SEQUENCE",
    "DUPLICATE_OR_ALTERNATIVE_ACCESSION",
    "CANONICAL_ACCESSION_USED_INSTEAD",
    "NON_DNA_BINDING_COMPONENT",
    "PARTIAL_ACCESSION_LACKS_DBD",
    "INCOMPLETE_DBD",
    "FUSION_OR_SEQUENCE_PROBLEM",
    "ALREADY_HAD_UNINDEXED_VALID_MODCRE_MODEL",
    "UNKNOWN_FILTER_REASON",
}

ALLOWED_DBD_INTERPRETATIONS = {
    "NO_PFAM_ANNOTATION",
    "NO_ACCEPTED_HIT_UNDER_PREVIOUS_THRESHOLDS",
    "PARTIAL_ACCESSION_MISSING_DBD",
    "INCOMPLETE_DBD",
    "TRUE_NO_EXPECTED_DBD",
    "CANONICAL_ONLY_HAS_DBD",
    "FUSION_OR_COMPLEX_CASE",
    "PFAM_MISSED_BUT_SEQUENCE_EVIDENCE_SUPPORTS_DBD",
    "UNRESOLVED",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build the read-only pre-Baldo structure/FASTA/Pfam reconciliation reports."
    )
    parser.add_argument("--database", type=Path, default=ROOT / "data/tf_webdb.sqlite")
    parser.add_argument(
        "--staging-database", type=Path,
        default=ROOT / "data/tf_webdb_local_staging.sqlite",
    )
    parser.add_argument("--baldo-fasta", type=Path, default=DEFAULT_FASTA)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "outputs")
    parser.add_argument(
        "--failed-without-pfam", type=Path,
        default=ROOT / "external/baldo_model_inventory/af3_failed_without_pfam.tsv",
    )
    parser.add_argument(
        "--structural-draft", type=Path,
        default=ROOT / "outputs/af3_failed_102_structural_status_draft.tsv",
    )
    parser.add_argument(
        "--final-adjudication", type=Path,
        default=ROOT / "outputs/af3_failed_102_final_adjudication_DRAFT.tsv",
    )
    parser.add_argument(
        "--af3-inventory", type=Path,
        default=ROOT / "external/baldo_model_inventory/af3_candidate_db_eligibility.tsv",
    )
    parser.add_argument(
        "--af3-interface-summary", type=Path,
        default=ROOT / "external/baldo_model_inventory/af3_tf_interface_summary.tsv",
    )
    parser.add_argument(
        "--domain-inventory", type=Path,
        default=ROOT / "outputs/af3_failed_102_domain_inventory.tsv",
    )
    parser.add_argument(
        "--pwm-structure-audit", type=Path,
        default=ROOT / "outputs/pwm_accession_structure_audit.tsv",
    )
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Iterable[dict[str, Any]], columns: list[str]) -> None:
    with path.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def sqlite_connection(path: Path) -> sqlite3.Connection:
    connection = sqlite3.connect(f"file:{path.resolve()}?mode=ro&immutable=1", uri=True)
    connection.row_factory = sqlite3.Row
    return connection


def database_snapshot(path: Path) -> dict[str, Any]:
    stat = path.stat()
    with sqlite_connection(path) as connection:
        quick_check = connection.execute("PRAGMA quick_check").fetchone()[0]
        tables = [
            row[0] for row in connection.execute(
                "SELECT name FROM sqlite_master WHERE type='table' ORDER BY name"
            )
        ]
        row_counts = {
            table: connection.execute(f'SELECT COUNT(*) FROM "{table}"').fetchone()[0]
            for table in tables
        }
    return {
        "path": str(path.resolve()), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns,
        "sha256": sha256(path), "quick_check": quick_check, "row_counts": row_counts,
        "total_rows": sum(row_counts.values()),
    }


def tree_snapshot(path: Path) -> dict[str, Any]:
    root_stat = path.stat()
    file_count = 0
    total_size = 0
    maximum_mtime_ns = 0
    for directory, _, files in os.walk(path):
        for name in files:
            stat = (Path(directory) / name).stat()
            file_count += 1
            total_size += stat.st_size
            maximum_mtime_ns = max(maximum_mtime_ns, stat.st_mtime_ns)
    return {
        "path": str(path.resolve()), "root_size": root_stat.st_size,
        "root_mtime_ns": root_stat.st_mtime_ns, "file_count": file_count,
        "total_file_size": total_size, "maximum_file_mtime_ns": maximum_mtime_ns,
    }


def website_snapshot() -> dict[str, str]:
    candidates = [ROOT / "app.py"]
    for base in [ROOT / "templates", ROOT / "static"]:
        if base.exists():
            candidates.extend(path for path in base.rglob("*") if path.is_file())
    return {str(path.relative_to(ROOT)): sha256(path) for path in sorted(set(candidates))}


def parse_fasta(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    records: list[dict[str, str]] = []
    errors: list[str] = []
    header = ""
    sequence: list[str] = []

    def finish() -> None:
        if not header:
            return
        accession = header.split("|", 1)[0].split()[0].upper()
        seq = "".join(sequence).upper()
        metadata: dict[str, str] = {}
        prefix = header.split(" ", 1)[0]
        for token in prefix.split("|")[1:]:
            if "=" in token:
                key, value = token.split("=", 1)
                metadata[key] = value
        if not accession:
            errors.append(f"Unparseable header: {header}")
        records.append({
            "accession": accession, "header": header, "sequence": seq,
            "gene": metadata.get("gene", ""), "primary": metadata.get("primary", ""),
            "family": metadata.get("family", ""),
        })

    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if line.startswith(">"):
                finish()
                header = line[1:]
                sequence = []
            elif line:
                if not header:
                    errors.append("Sequence content before first FASTA header")
                sequence.append(line)
    finish()
    return records, errors


def clean(value: Any) -> str:
    return "" if value is None else str(value).strip()


def primary_gene(gene_names: str, fasta_gene: str = "") -> str:
    value = clean(gene_names)
    if value:
        return value.split()[0]
    value = clean(fasta_gene)
    return value.split(",", 1)[0]


def add_check(checks: dict[str, Any], name: str, observed: Any, expected: Any) -> None:
    checks[name] = {"observed": observed, "expected": expected, "passed": observed == expected}


def fusion_expected_dbd(final_row: dict[str, str]) -> str:
    raw = clean(final_row.get("fusion_parent_analysis"))
    if not raw:
        return ""
    try:
        parent_rows = json.loads(raw)
    except json.JSONDecodeError:
        return ""
    return ";".join(sorted({clean(row.get("expected_DBD_type")) for row in parent_rows if clean(row.get("expected_DBD_type"))}))


def dbd_interpretation(final_row: dict[str, str]) -> str:
    status = clean(final_row.get("primary_structural_status"))
    uncertainty = clean(final_row.get("remaining_uncertainty"))
    if status == "ACCESSION_LACKS_DBD":
        return "PARTIAL_ACCESSION_MISSING_DBD"
    if status == "ACCESSION_CONTAINS_INCOMPLETE_DBD":
        return "INCOMPLETE_DBD"
    if status == "COFACTOR_OR_NON_DNA_BINDING_COMPONENT":
        return "TRUE_NO_EXPECTED_DBD"
    if "Fusion" in uncertainty or clean(final_row.get("fusion_parent_analysis")):
        return "FUSION_OR_COMPLEX_CASE"
    if status == "VALID_DBD_FRAGMENT_MODEL" or clean(final_row.get("DBD_disposition")) == "COMPLETE_DBD":
        return "PFAM_MISSED_BUT_SEQUENCE_EVIDENCE_SUPPORTS_DBD"
    return "UNRESOLVED"


def main() -> int:
    args = parse_args()
    if socket.gethostname() != "masada":
        raise SystemExit(f"Refusing to run: hostname is {socket.gethostname()!r}, expected 'masada'")

    inputs = [
        args.database, args.staging_database, args.baldo_fasta, args.failed_without_pfam,
        args.structural_draft, args.final_adjudication, args.af3_inventory,
        args.af3_interface_summary, args.domain_inventory, args.pwm_structure_audit,
    ]
    missing = [str(path) for path in inputs if not path.is_file()]
    if missing:
        raise SystemExit("Missing required input(s): " + ", ".join(missing))

    args.output_dir.mkdir(parents=True, exist_ok=True)
    output_paths = {key: args.output_dir / name for key, name in OUTPUT_NAMES.items()}
    existing = [str(path) for path in output_paths.values() if path.exists()]
    if existing:
        raise SystemExit("Refusing to overwrite existing output(s): " + ", ".join(existing))

    database_before = {
        "source": database_snapshot(args.database),
        "staging": database_snapshot(args.staging_database),
    }
    model_before = {str(path): tree_snapshot(path) for path in [MODEL_ROOT, AF3_ROOT]}
    website_before = website_snapshot()

    fasta_records, fasta_errors = parse_fasta(args.baldo_fasta)
    fasta_counts = Counter(row["accession"] for row in fasta_records)
    fasta_by_accession = {row["accession"]: row for row in fasta_records}
    fasta_ids = set(fasta_by_accession)

    with sqlite_connection(args.database) as connection:
        tf_rows = [dict(row) for row in connection.execute(
            """
            SELECT t.tf_id, t.family_text, t.motif_ref_count,
                   ta.uniprot_accession, ta.gene_names, ta.protein_name,
                   ta.reviewed, ta.sequence_length,
                   tpa.best_annotation_level, tpa.best_pwm_or_model
            FROM tf AS t
            LEFT JOIN tf_annotation AS ta ON ta.tf_id = t.tf_id
            LEFT JOIN tf_primary_annotation AS tpa ON tpa.tf_id = t.tf_id
            ORDER BY t.tf_id
            """
        )]
        tf_ids = {clean(row["tf_id"]).upper() for row in tf_rows}

        active_rows: dict[str, list[dict[str, Any]]] = defaultdict(list)
        active_source_accessions: dict[str, set[str]] = defaultdict(set)
        for row_obj in connection.execute(
            """SELECT source, model_id, tf_id, member_path, archive_path
               FROM structure_file
               WHERE status='active' AND file_type='pdb'
               ORDER BY tf_id, source, model_id, member_path"""
        ):
            row = dict(row_obj)
            accession = clean(row["tf_id"]).upper()
            if accession:
                active_rows[accession].append(row)
                active_source_accessions[clean(row["source"])].add(accession)

        pwm_counts = {
            clean(row["tf_id"]).upper(): int(row["count"])
            for row in connection.execute(
                """
                SELECT mr.tf_id, COUNT(DISTINCT mr.source || '|' || mr.motif_id) AS count
                FROM motif_ref AS mr
                JOIN motif_file AS mf ON mf.source=mr.source AND mf.motif_id=mr.motif_id
                WHERE mr.missing_local_file=0 AND mf.matrix_status='usable'
                GROUP BY mr.tf_id
                """
            )
        }
        provenance: dict[str, list[str]] = defaultdict(list)
        for row in connection.execute(
            """SELECT tf_id, source, evidence_type, COUNT(DISTINCT motif_id) AS count
               FROM motif_ref GROUP BY tf_id, source, evidence_type
               ORDER BY tf_id, source, evidence_type"""
        ):
            accession = clean(row["tf_id"]).upper()
            provenance[accession].append(
                f"{clean(row['source'])}:{clean(row['evidence_type'])}:{int(row['count'])}"
            )
        pfam_by_accession: dict[str, list[str]] = defaultdict(list)
        for row in connection.execute(
            """SELECT tf_id, pfam_id, pfam_name, start, end
               FROM tf_pfam_annotation ORDER BY tf_id, start, end, pfam_id"""
        ):
            pfam_by_accession[clean(row["tf_id"]).upper()].append(
                f"{clean(row['pfam_id'])}:{row['start']}-{row['end']}:{clean(row['pfam_name'])}"
            )

    structural_rows = {row["failed_accession"].upper(): row for row in read_tsv(args.structural_draft)}
    final_rows = {row["failed_accession"].upper(): row for row in read_tsv(args.final_adjudication)}
    no_pfam_source = read_tsv(args.failed_without_pfam)
    no_pfam_ids = {row["tf_id"].upper() for row in no_pfam_source}
    pwm_audit = {row["tf_id"].upper(): row for row in read_tsv(args.pwm_structure_audit)}

    af3_inventory = {
        row["tf_id"].upper(): row for row in read_tsv(args.af3_inventory)
        if clean(row.get("eligible_for_af3_import")).lower() == "yes"
    }
    interface_summary = {
        row["tf_id"].upper(): row for row in read_tsv(args.af3_interface_summary)
    }
    fragment_rows: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in read_tsv(args.domain_inventory):
        if clean(row.get("fragment_id")) and int(clean(row.get("fragment_models_checked")) or 0) > 0:
            fragment_rows[row["failed_accession"].upper()].append(row)

    full_length_ids = fasta_ids & set(af3_inventory)
    valid_full_ids = {
        accession for accession in full_length_ids
        if int(clean(interface_summary.get(accession, {}).get("n_pass")) or 0) > 0
    }
    empty_full_ids = {
        accession for accession in full_length_ids
        if int(clean(interface_summary.get(accession, {}).get("n_pass")) or 0) == 0
    }
    fragment_ids = fasta_ids & set(fragment_rows)
    valid_fragment_ids = {
        accession for accession in fragment_ids
        if clean(final_rows.get(accession, {}).get("valid_existing_fragment_model")) == "YES"
    }

    master_rows: list[dict[str, Any]] = []
    for source in tf_rows:
        accession = clean(source["tf_id"]).upper()
        fasta_row = fasta_by_accession.get(accession, {})
        active = active_rows.get(accession, [])
        structural = structural_rows.get(accession, {})
        final = final_rows.get(accession, {})
        expected_dbd = clean(final.get("expected_DBD_type")) or clean(structural.get("canonical_reference_DBD_context"))
        if not expected_dbd:
            expected_dbd = fusion_expected_dbd(final)
        previous_pfam = "NO_ACCEPTED_HIT_UNDER_PREVIOUS_THRESHOLDS" if accession in no_pfam_ids else ""
        interpretation = dbd_interpretation(final) if accession in no_pfam_ids else ""

        if active:
            reason_no_original = ""
        elif accession == "P01100":
            reason_no_original = "ALREADY_HAD_UNINDEXED_VALID_MODCRE_MODEL"
        else:
            reason_no_original = "EXTERNAL_PWM_ONLY_NOT_STRUCTURALLY_MODELLED"

        reason_not_sent = ""
        if not active and accession not in fasta_ids:
            reason_not_sent = "MISSING_OR_INVALID_SEQUENCE"

        notes: list[str] = []
        if active:
            notes.append(
                "Pre-Baldo ModCREDB structure population uses the website SQL definition: "
                "status=active and file_type=pdb across indexed structure sources."
            )
        else:
            notes.append(
                "Original SQLite has no active indexed PDB and no ModCRE/AlphaFold motif or structure row; "
                "PWM evidence is external-resource provenance."
            )
        if reason_not_sent:
            notes.append(
                "Record has external CIS-BP/JASPAR PWM evidence but protein_name=deleted/merged, no local "
                "UniProt sequence annotation/length, and no sequence in the Baldo FASTA."
            )
        if accession == "P01100":
            p01100 = pwm_audit.get(accession, {})
            notes.append(
                "Control: original website 3D=0 reflects a database-indexing gap, not model-generation "
                f"failure; {clean(p01100.get('valid_unindexed_structure_count'))} exact-P01100 local "
                f"ModCRE models were validated, best={clean(p01100.get('best_valid_model_path'))}."
            )
        if accession == "E9PN75":
            notes.append(
                "Control: Q9NZC4 is canonical same-gene reference only; alignment review finds E9PN75 "
                "lacks the expected PF00178 ETS DBD."
            )
        if accession == "H7C4N4":
            notes.append(
                "Control: full-length AF3 is empty-interface; the reviewed PF03299 fragment is the valid DBD rescue."
            )
        if accession == "B4DHE0":
            notes.append("Control: contacting PF01017 fragment is non-DBD and is not a valid DBD structure.")
        if accession == "Q59EF3":
            notes.append("Control: contacting PF17725 YAP-binding fragment is not the PF01285 TEA DBD.")

        original_paths = [
            f"{clean(row['source'])}|{clean(row['archive_path'])}::{clean(row['member_path'])}"
            for row in active
        ]
        master_rows.append({
            "tf_id": accession,
            "accession": clean(source.get("uniprot_accession")) or accession,
            "gene": primary_gene(clean(source.get("gene_names")), clean(fasta_row.get("gene"))),
            "protein_name": clean(source.get("protein_name")),
            "PWM_annotation": clean(source.get("best_annotation_level")) or clean(fasta_row.get("primary")),
            "PWM_count": pwm_counts.get(accession, 0),
            "PWM_source_or_provenance": ";".join(provenance.get(accession, [])),
            "original_active_modcre_model": "YES" if active else "NO",
            "original_active_modcre_model_count": len(active),
            "original_model_path": ";".join(original_paths),
            "sent_to_baldo": "YES" if accession in fasta_ids else "NO",
            "baldo_fasta_path": str(args.baldo_fasta.resolve()) if accession in fasta_ids else "",
            "baldo_full_length_AF3_found": "YES" if accession in full_length_ids else "NO",
            "baldo_valid_interface_found": "YES" if accession in valid_full_ids else "NO",
            "baldo_empty_interface_only": "YES" if accession in empty_full_ids else "NO",
            "baldo_fragment_found": "YES" if accession in fragment_ids else "NO",
            "expected_DBD": expected_dbd,
            "previous_pfam_status": previous_pfam,
            "final_DBD_interpretation": interpretation,
            "canonical_accession": clean(final.get("reviewed_canonical_accession")) or clean(structural.get("canonical_reference_accession")),
            "reason_no_original_modcre_model": reason_no_original,
            "reason_not_sent_to_baldo": reason_not_sent,
            "notes": " ".join(notes),
        })

    master_by_accession = {row["accession"]: row for row in master_rows}
    no_model_rows = [row for row in master_rows if row["original_active_modcre_model"] == "NO"]
    sent_rows = [row for row in master_rows if row["sent_to_baldo"] == "YES"]
    not_sent_rows = [
        row for row in master_rows
        if row["original_active_modcre_model"] == "NO" and row["sent_to_baldo"] == "NO"
    ]

    pfam_review_rows: list[dict[str, Any]] = []
    for source in sorted(no_pfam_source, key=lambda row: row["tf_id"]):
        accession = source["tf_id"].upper()
        master = master_by_accession[accession]
        structural = structural_rows[accession]
        final = final_rows[accession]
        expected_dbd = master["expected_DBD"]
        interpretation = master["final_DBD_interpretation"]
        evidence = (
            f"Previous failed-accession Pfam38 fragment scan: {clean(source.get('reason'))}; "
            f"canonical/reference alignment review: status={clean(final.get('primary_structural_status'))}, "
            f"DBD disposition={clean(final.get('DBD_disposition'))}, coverage="
            f"{clean(structural.get('canonical_DBD_coverage_percent')) or clean(final.get('DBD_coverage_percent'))}%, "
            f"identity={clean(structural.get('canonical_DBD_identity_percent')) or clean(final.get('DBD_identity_percent'))}%. "
            f"{clean(final.get('decision_reason'))}"
        )
        pfam_review_rows.append({
            "tf_id": accession, "accession": accession, "gene": master["gene"],
            "protein_name": master["protein_name"], "PWM_annotation": master["PWM_annotation"],
            "PWM_count": master["PWM_count"],
            "PWM_source_or_provenance": master["PWM_source_or_provenance"],
            "reviewed_status": clean(structural.get("reviewed_status")),
            "sequence_length": clean(structural.get("failed_sequence_length")),
            "canonical_accession": master["canonical_accession"], "expected_DBD": expected_dbd,
            "previous_pfam_status": "NO_ACCEPTED_HIT_UNDER_PREVIOUS_THRESHOLDS",
            "current_sqlite_pfam_status": "NO_PFAM_ANNOTATION" if not pfam_by_accession.get(accession) else ";".join(pfam_by_accession[accession]),
            "DBD_coverage_percent": clean(structural.get("canonical_DBD_coverage_percent")) or clean(final.get("DBD_coverage_percent")),
            "DBD_identity_percent": clean(structural.get("canonical_DBD_identity_percent")) or clean(final.get("DBD_identity_percent")),
            "failed_coordinates_corresponding_to_DBD": clean(structural.get("failed_coordinates_corresponding_to_DBD")) or clean(final.get("failed_coordinates_corresponding_to_DBD")),
            "final_DBD_interpretation": interpretation,
            "evidence_interpretation": evidence,
            "remaining_uncertainty": clean(final.get("remaining_uncertainty")),
        })

    return_rows: list[dict[str, Any]] = []
    for accession in sorted(fasta_ids):
        inventory = af3_inventory.get(accession, {})
        summary = interface_summary.get(accession, {})
        fragments = fragment_rows.get(accession, [])
        full = accession in full_length_ids
        valid = accession in valid_full_ids
        empty = accession in empty_full_ids
        fragment = accession in fragment_ids
        valid_fragment = accession in valid_fragment_ids
        notes = []
        if full:
            notes.append(
                f"Existing local AF3 interface audit: n_pass={clean(summary.get('n_pass'))}, "
                f"n_fail={clean(summary.get('n_fail'))}."
            )
        if fragment:
            notes.append("Fragment return is secondary and is not assumed to be a DBD without reviewed domain-role evidence.")
        if accession == "H7C4N4":
            notes.append("PF03299 fragment is the one reviewed valid DBD rescue; full-length AF3 interface is empty.")
        if accession == "E9PN75":
            notes.append("PF02198 is a SAM/Pointed region, not the missing PF00178 ETS DBD; Q9NZC4 remains reference-only.")
        if accession == "B4DHE0":
            notes.append("Contacting PF01017 fragment is non-DBD and does not rescue the accession structure.")
        if accession == "Q59EF3":
            notes.append("Contacting PF17725 fragment is not the PF01285 TEA DBD and does not rescue the accession structure.")
        return_rows.append({
            "accession": accession, "sent_to_baldo": "YES",
            "full_length_AF3_found": "YES" if full else "NO",
            "full_length_model_paths": clean(inventory.get("selected_cifs")),
            "full_length_variant_count": clean(summary.get("n_variants")) if full else "0",
            "valid_full_length_interface_found": "YES" if valid else "NO",
            "empty_interface_only": "YES" if empty else "NO",
            "fragment_found": "YES" if fragment else "NO",
            "fragment_ids": ";".join(clean(row.get("fragment_id")) for row in fragments),
            "fragment_model_paths": ";".join(clean(row.get("fragment_best_model_path")) for row in fragments if clean(row.get("fragment_best_model_path"))),
            "valid_DBD_fragment_rescue": "YES" if valid_fragment else "NO",
            "no_returned_model_found": "YES" if not full and not fragment else "NO",
            "notes": " ".join(notes),
        })

    checks: dict[str, Any] = {}
    active_ids = set(active_rows)
    no_active_ids = tf_ids - active_ids
    add_check(checks, "hostname_is_masada", socket.gethostname(), "masada")
    add_check(checks, "master_rows_unique", (len(master_rows), len(master_by_accession)), (5267, 5267))
    add_check(checks, "original_active_indexed_population", len(active_ids), 1344)
    add_check(checks, "original_no_active_indexed_population", len(no_active_ids), 3923)
    add_check(checks, "fasta_records_unique", (len(fasta_records), len(fasta_ids)), (3690, 3690))
    add_check(checks, "fasta_parser_errors", fasta_errors, [])
    add_check(checks, "fasta_duplicate_accessions", sorted(k for k, v in fasta_counts.items() if v > 1), [])
    add_check(checks, "fasta_unmapped_accessions", sorted(fasta_ids - tf_ids), [])
    add_check(checks, "original_model_accessions_sent", len(fasta_ids & active_ids), 0)
    add_check(checks, "no_model_accessions_sent", len(fasta_ids & no_active_ids), 3690)
    add_check(checks, "no_model_accessions_not_sent", len(no_active_ids - fasta_ids), 233)
    add_check(checks, "not_sent_reasons_allowed", sorted({row["reason_not_sent_to_baldo"] for row in not_sent_rows}), ["MISSING_OR_INVALID_SEQUENCE"])
    add_check(checks, "not_sent_deleted_or_merged_without_sequence", dict(sorted(Counter(
        row["protein_name"] for row in not_sent_rows
        if row["reason_not_sent_to_baldo"] == "MISSING_OR_INVALID_SEQUENCE"
    ).items())), {"deleted": 230, "merged": 3})
    add_check(checks, "previous_no_pfam_cohort", (len(pfam_review_rows), len(no_pfam_ids)), (45, 45))
    add_check(checks, "dbd_interpretations_allowed", sorted({row["final_DBD_interpretation"] for row in pfam_review_rows} - ALLOWED_DBD_INTERPRETATIONS), [])
    add_check(checks, "full_length_af3_returned", len(full_length_ids), 211)
    add_check(checks, "valid_full_length_interfaces", len(valid_full_ids), 109)
    add_check(checks, "empty_interface_only_full_length", len(empty_full_ids), 102)
    add_check(checks, "fragment_returned_accessions", len(fragment_ids), 57)
    add_check(checks, "valid_dbd_fragment_rescues", sorted(valid_fragment_ids), ["H7C4N4"])
    add_check(checks, "no_baldo_model_found", sum(row["no_returned_model_found"] == "YES" for row in return_rows), 3479)
    add_check(checks, "all_full_length_returns_have_interface_summary", sorted(full_length_ids - set(interface_summary)), [])

    controls = {
        "P01100": {
            "sent_to_baldo": master_by_accession["P01100"]["sent_to_baldo"],
            "original_indexed_active_model": master_by_accession["P01100"]["original_active_modcre_model"],
            "reason": master_by_accession["P01100"]["reason_no_original_modcre_model"],
            "valid_unindexed_models": pwm_audit["P01100"]["valid_unindexed_structure_count"],
            "best_model": pwm_audit["P01100"]["best_valid_model_path"],
            "model_mtime_ns": Path(pwm_audit["P01100"]["best_valid_model_path"]).stat().st_mtime_ns,
            "earliest_returned_af3_mtime_ns": min(
                Path(path).stat().st_mtime_ns
                for row in af3_inventory.values()
                for path in clean(row.get("selected_cifs")).split(";") if path and Path(path).is_file()
            ),
        },
        "E9PN75": {
            "own_expected_ETS_DBD": "ABSENT_BY_CANONICAL_ALIGNMENT_REVIEW",
            "canonical_reference": master_by_accession["E9PN75"]["canonical_accession"],
            "canonical_reference_counted_as_own": False,
            "full_length_empty_interface": "E9PN75" in empty_full_ids,
        },
        "H7C4N4": {
            "sent_to_baldo": "H7C4N4" in fasta_ids,
            "full_length_empty_interface": "H7C4N4" in empty_full_ids,
            "PF03299_fragment_found": "H7C4N4" in fragment_ids,
            "valid_DBD_fragment_rescue": "H7C4N4" in valid_fragment_ids,
        },
        "B4DHE0": {
            "full_length_empty_interface": "B4DHE0" in empty_full_ids,
            "contacting_fragment": "PF01017",
            "fragment_counted_as_valid_DBD": "B4DHE0" in valid_fragment_ids,
        },
        "Q59EF3": {
            "full_length_empty_interface": "Q59EF3" in empty_full_ids,
            "contacting_fragment": "PF17725",
            "expected_DBD": master_by_accession["Q59EF3"]["expected_DBD"],
            "fragment_counted_as_valid_DBD": "Q59EF3" in valid_fragment_ids,
        },
    }
    controls_pass = (
        controls["P01100"]["reason"] == "ALREADY_HAD_UNINDEXED_VALID_MODCRE_MODEL"
        and int(controls["P01100"]["valid_unindexed_models"]) == 20
        and controls["P01100"]["model_mtime_ns"] < controls["P01100"]["earliest_returned_af3_mtime_ns"]
        and controls["E9PN75"]["canonical_reference"] == "Q9NZC4"
        and not controls["E9PN75"]["canonical_reference_counted_as_own"]
        and controls["E9PN75"]["own_expected_ETS_DBD"] == "ABSENT_BY_CANONICAL_ALIGNMENT_REVIEW"
        and controls["H7C4N4"]["sent_to_baldo"]
        and controls["H7C4N4"]["full_length_empty_interface"]
        and controls["H7C4N4"]["PF03299_fragment_found"]
        and controls["H7C4N4"]["valid_DBD_fragment_rescue"]
        and controls["B4DHE0"]["contacting_fragment"] == "PF01017"
        and not controls["B4DHE0"]["fragment_counted_as_valid_DBD"]
        and controls["Q59EF3"]["contacting_fragment"] == "PF17725"
        and "PF01285" in controls["Q59EF3"]["expected_DBD"]
        and not controls["Q59EF3"]["fragment_counted_as_valid_DBD"]
    )
    add_check(checks, "all_controls_pass", controls_pass, True)

    if not all(check["passed"] for check in checks.values()):
        failed = [name for name, check in checks.items() if not check["passed"]]
        raise SystemExit("QC failed before output creation: " + ", ".join(failed))

    write_tsv(output_paths["master"], master_rows, MASTER_COLUMNS)
    write_tsv(output_paths["no_model"], no_model_rows, MASTER_COLUMNS)
    write_tsv(output_paths["sent"], sent_rows, MASTER_COLUMNS)
    write_tsv(output_paths["not_sent"], not_sent_rows, MASTER_COLUMNS)
    write_tsv(output_paths["pfam"], pfam_review_rows, PFAM_REVIEW_COLUMNS)
    write_tsv(output_paths["returns"], return_rows, RETURN_COLUMNS)

    database_after = {
        "source": database_snapshot(args.database),
        "staging": database_snapshot(args.staging_database),
    }
    model_after = {str(path): tree_snapshot(path) for path in [MODEL_ROOT, AF3_ROOT]}
    website_after = website_snapshot()
    add_check(checks, "databases_unchanged", database_after, database_before)
    add_check(checks, "model_directories_unchanged", model_after, model_before)
    add_check(checks, "website_files_unchanged", website_after, website_before)
    add_check(checks, "no_remote_connection", False, False)

    qc = {
        "validation_status": "PASS" if all(check["passed"] for check in checks.values()) else "FAIL",
        "read_only": True,
        "remote_connection_attempted": False,
        "database_open_mode": "mode=ro&immutable=1",
        "definition": {
            "original_active_modcre_model": (
                "Any structure_file row in the original ModCREDB database with status='active' and "
                "file_type='pdb', matching the website 3D-count logic. This is the 1,344-accession "
                "union across source='modcre' and legacy source='alphafold'; source is retained in paths."
            ),
            "narrow_source_counts": {
                source: len(accessions) for source, accessions in sorted(active_source_accessions.items())
            },
            "pwm_count": (
                "COUNT(DISTINCT source|motif_id) joined to motif_file where missing_local_file=0 "
                "and matrix_status='usable', matching website scan-ready PWM logic."
            ),
            "previous_no_pfam_cohort": str(args.failed_without_pfam.resolve()),
        },
        "counts": {
            "total_original_accessions": len(tf_ids),
            "original_active_indexed_model_accessions": len(active_ids),
            "original_no_active_indexed_model_accessions": len(no_active_ids),
            "fasta_records": len(fasta_records), "fasta_unique_accessions": len(fasta_ids),
            "fasta_duplicate_accessions": len(fasta_records) - len(fasta_ids),
            "no_model_sent": len(no_active_ids & fasta_ids),
            "no_model_not_sent": len(no_active_ids - fasta_ids),
            "original_model_sent": len(active_ids & fasta_ids),
            "fasta_not_in_database": len(fasta_ids - tf_ids),
            "previous_no_pfam": len(no_pfam_ids),
            "baldo_full_length_af3": len(full_length_ids),
            "baldo_valid_full_length_interface": len(valid_full_ids),
            "baldo_empty_interface_only": len(empty_full_ids),
            "baldo_fragment_accessions": len(fragment_ids),
            "baldo_valid_dbd_fragment_rescue": len(valid_fragment_ids),
            "baldo_no_model_found": sum(row["no_returned_model_found"] == "YES" for row in return_rows),
        },
        "reason_no_original_modcre_model_counts": dict(sorted(Counter(
            row["reason_no_original_modcre_model"] or "ORIGINAL_ACTIVE_INDEXED_MODEL"
            for row in master_rows
        ).items())),
        "reason_not_sent_to_baldo_counts": dict(sorted(Counter(
            row["reason_not_sent_to_baldo"] for row in not_sent_rows
        ).items())),
        "final_DBD_interpretation_counts": dict(sorted(Counter(
            row["final_DBD_interpretation"] for row in pfam_review_rows
        ).items())),
        "controls": controls,
        "checks": checks,
        "database_before": database_before,
        "database_after": database_after,
        "model_directories_before": model_before,
        "model_directories_after": model_after,
        "website_before": website_before,
        "website_after": website_after,
        "fasta": {
            "path": str(args.baldo_fasta.resolve()), "sha256": sha256(args.baldo_fasta),
            "identical_repository_copy": (
                sha256(ROOT / "outputs/TF_without_model.fasta") == sha256(args.baldo_fasta)
                if (ROOT / "outputs/TF_without_model.fasta").is_file() else False
            ),
            "unmapped_headers": sorted(fasta_ids - tf_ids),
        },
        "created_files": [str(output_paths[key].resolve()) for key in ["master", "no_model", "sent", "not_sent", "pfam", "returns", "qc"]],
    }
    with output_paths["qc"].open("x", encoding="utf-8") as handle:
        json.dump(qc, handle, indent=2, sort_keys=True)
        handle.write("\n")

    if qc["validation_status"] != "PASS":
        raise SystemExit("Post-write protected-state QC failed")
    print(json.dumps({
        "validation_status": qc["validation_status"], "counts": qc["counts"],
        "reason_no_original_modcre_model_counts": qc["reason_no_original_modcre_model_counts"],
        "reason_not_sent_to_baldo_counts": qc["reason_not_sent_to_baldo_counts"],
        "final_DBD_interpretation_counts": qc["final_DBD_interpretation_counts"],
        "controls": qc["controls"], "created_files": qc["created_files"],
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
