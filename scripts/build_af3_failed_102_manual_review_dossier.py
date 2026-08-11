#!/usr/bin/env python3
"""Build a read-only manual-review dossier for the 36 unresolved AF3 cases."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sqlite3
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from classify_af3_failed_102_structural_status import (
    CANONICAL_DBD_PFAMS,
    canonical_domains,
    choose_model_path,
    db_signature,
    fasta_records,
    file_sha256,
    local_alignment_metrics,
    parse_model_interval,
    reconcile_pdb,
)


ROOT = Path(__file__).resolve().parents[1]
DECISIONS = {
    "LIKELY_ACCESSION_LACKS_DBD",
    "LIKELY_DBD_PRESENT",
    "PARTIAL_DBD_PRESENT",
    "FUSION_COMPONENT_REQUIRES_REVIEW",
    "COFACTOR_OR_NON_DNA_BINDING_RECORD",
    "GENE_MAPPING_UNRESOLVED",
    "INSUFFICIENT_EVIDENCE",
}
ACTIONS = {
    "LIKELY_NO_SEND_ACCESSION_LACKS_DBD",
    "LIKELY_SEND_NEW_DBD_FRAGMENT",
    "LIKELY_SEND_CORRECTED_DBD_FRAGMENT",
    "LIKELY_NO_SEND_CANONICAL_REFERENCE_ONLY",
    "REQUIRES_MANUAL_REVIEW",
}
FUSION_CONTROLS = {"A0A0A7M2K0", "H2BNB9", "B1NY96", "D1LUZ8"}
UNRESOLVED_CONTROLS = {"A8K549", "B4DTN3", "Q68D60"}

DOSSIER_COLUMNS = [
    "failed_accession", "resolved_primary_gene", "gene_resolution_status",
    "gene_resolution_evidence", "failed_protein_description",
    "failed_sequence_length", "failed_sequence_sha256",
    "reviewed_canonical_accession", "reviewed_canonical_gene_names",
    "reviewed_canonical_sequence_length", "reviewed_canonical_selection_reason",
    "model_bearing_accession", "model_bearing_is_reviewed",
    "model_bearing_sequence_length", "model_bearing_selection_reason",
    "model_path", "model_interval", "existing_interface_result",
    "canonical_DBD_type", "canonical_DBD_pfam_ids",
    "canonical_DBD_start", "canonical_DBD_end",
    "canonical_DBD_length", "canonical_DBD_coordinate_source",
    "failed_to_canonical_alignment_columns",
    "failed_to_canonical_aligned_residues",
    "failed_to_canonical_identity_percent",
    "canonical_DBD_aligned_residues", "canonical_DBD_coverage_percent",
    "canonical_DBD_identity_percent",
    "failed_coordinates_corresponding_to_DBD",
    "failed_matched_canonical_span", "DBD_sequence_relationship",
    "aligned_non_DBD_residues", "aligned_non_DBD_identity_percent",
    "pdb_protein_chain", "pdb_protein_sequence",
    "pdb_protein_sequence_length", "pdb_protein_sequence_sha256",
    "pdb_to_model_accession_aligned_residues",
    "pdb_to_model_accession_identity_percent",
    "pdb_to_model_accession_coverage_percent",
    "pdb_model_accession_reconciliation_status",
    "pdb_to_canonical_DBD_aligned_residues",
    "pdb_to_canonical_DBD_coverage_percent",
    "pdb_to_canonical_DBD_identity_percent",
    "pdb_represents_canonical_DBD", "canonical_and_model_accessions_separate",
    "suggested_decision", "suggested_baldo_action",
    "decision_evidence", "manual_review_notes",
]

ALIGNMENT_COLUMNS = [
    "failed_accession", "reviewed_canonical_accession",
    "model_bearing_accession", "canonical_DBD_start", "canonical_DBD_end",
    "failed_alignment", "canonical_alignment",
    "canonical_DBD_alignment_mask", "reference_to_failed_coordinate_map",
    "failed_matched_canonical_positions", "DBD_sequence_relationship",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build the 36-accession read-only manual-review dossier.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ext = ROOT / "external/baldo_model_inventory"
    parser.add_argument("--manual-review", type=Path, default=ROOT / "outputs/af3_failed_102_manual_review.tsv")
    parser.add_argument("--structural-draft", type=Path, default=ROOT / "outputs/af3_failed_102_structural_status_draft.tsv")
    parser.add_argument("--reference-alignments", type=Path, default=ROOT / "outputs/af3_failed_102_reference_alignments.tsv")
    parser.add_argument("--domain-inventory", type=Path, default=ROOT / "outputs/af3_failed_102_domain_inventory.tsv")
    parser.add_argument("--reviewed-pfam-map", type=Path, default=ROOT / "data_sources/af3_failed_102_pfam_role_reviewed.tsv")
    parser.add_argument("--failed-fasta", type=Path, default=ext / "af3_failed_102_full_length.fasta")
    parser.add_argument("--same-gene-candidates", type=Path, default=ext / "af3_failed_same_gene_modcre_candidates_v2.tsv")
    parser.add_argument("--same-gene-audit", type=Path, default=ext / "af3_failed_same_gene_modcre_interface_audit.tsv")
    parser.add_argument("--database", type=Path, default=ROOT / "data/tf_webdb.sqlite")
    parser.add_argument("--all-tf-fasta", type=Path, default=Path("/home/patricia/TF_database_Baldo_data/TF_without_model.fasta"))
    parser.add_argument("--models-root", type=Path, default=Path("/data/sbi/interchange/boliva/patricia/models"))
    parser.add_argument("--output-dossier", type=Path, default=ROOT / "outputs/af3_failed_102_manual_review_dossier.tsv")
    parser.add_argument("--output-alignments", type=Path, default=ROOT / "outputs/af3_failed_102_manual_review_alignments.tsv")
    parser.add_argument("--output-qc", type=Path, default=ROOT / "outputs/af3_failed_102_manual_review_qc.json")
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def sequence_sha256(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def global_alignment(failed: str, canonical: str) -> dict[str, Any]:
    """Needleman-Wunsch alignment retaining aligned strings and coordinate map."""
    n, m = len(failed), len(canonical)
    gap = -2
    previous = [j * gap for j in range(m + 1)]
    trace = [bytearray(m + 1) for _ in range(n + 1)]
    for j in range(1, m + 1):
        trace[0][j] = 3
    for i in range(1, n + 1):
        trace[i][0] = 2
        current = [i * gap] + [0] * m
        aa = failed[i - 1]
        for j in range(1, m + 1):
            diagonal = previous[j - 1] + (2 if aa == canonical[j - 1] else -1)
            up = previous[j] + gap
            left = current[j - 1] + gap
            best = max(diagonal, up, left)
            current[j] = best
            trace[i][j] = 1 if best == diagonal else (2 if best == up else 3)
        previous = current
    i, j = n, m
    fa: list[str] = []
    ca: list[str] = []
    pairs: list[tuple[int | None, int | None, bool]] = []
    while i or j:
        direction = trace[i][j]
        if direction == 1:
            fa.append(failed[i - 1])
            ca.append(canonical[j - 1])
            pairs.append((i, j, failed[i - 1] == canonical[j - 1]))
            i -= 1
            j -= 1
        elif direction == 2:
            fa.append(failed[i - 1])
            ca.append("-")
            pairs.append((i, None, False))
            i -= 1
        else:
            fa.append("-")
            ca.append(canonical[j - 1])
            pairs.append((None, j, False))
            j -= 1
    fa.reverse()
    ca.reverse()
    pairs.reverse()
    mapping = {r: (f, match) for f, r, match in pairs if f and r}
    aligned = len(mapping)
    matches = sum(match for _, match in mapping.values())
    return {
        "failed_alignment": "".join(fa),
        "canonical_alignment": "".join(ca),
        "pairs": pairs,
        "mapping": mapping,
        "columns": len(pairs),
        "aligned": aligned,
        "matches": matches,
        "identity": 100.0 * matches / aligned if aligned else 0.0,
    }


def exact_gene_tokens(gene_names: str) -> set[str]:
    return {
        token.upper()
        for token in re.split(r"[\s,;/]+", gene_names)
        if token
    }


def select_reviewed_canonical(
    connection: sqlite3.Connection, gene: str
) -> tuple[str, str, int] | None:
    if not gene:
        return None
    candidates = connection.execute(
        """SELECT uniprot_accession, gene_names, sequence_length
           FROM tf_annotation WHERE reviewed=1
           ORDER BY annotation_score DESC, sequence_length DESC, uniprot_accession"""
    ).fetchall()
    gene_upper = gene.upper()
    exact = [
        (accession, genes or "", int(length))
        for accession, genes, length in candidates
        if gene_upper in exact_gene_tokens(genes or "")
    ]
    return exact[0] if exact else None


def canonical_dbd(
    domains: list[tuple[str, int, int, str]]
) -> tuple[list[tuple[str, int, int, str]], tuple[int, int] | None]:
    dbds = [domain for domain in domains if domain[0] in CANONICAL_DBD_PFAMS]
    if not dbds:
        return [], None
    return dbds, (min(d[1] for d in dbds), max(d[2] for d in dbds))


def topology_and_metrics(
    alignment: dict[str, Any], interval: tuple[int, int] | None,
    canonical_length: int,
) -> dict[str, Any]:
    if not interval:
        return {
            "length": 0, "covered": 0, "matches": 0, "coverage": 0.0,
            "identity": 0.0, "failed_coords": "", "matched_span": "",
            "relationship": "NO_CANONICAL_DBD_COORDINATES",
            "non_dbd_aligned": 0, "non_dbd_identity": 0.0,
        }
    start, end = interval
    mapped = [
        (position, alignment["mapping"][position])
        for position in range(start, end + 1)
        if position in alignment["mapping"]
    ]
    covered = len(mapped)
    matches = sum(value[1] for _, value in mapped)
    dbd_length = end - start + 1
    coverage = 100.0 * covered / dbd_length
    identity = 100.0 * matches / covered if covered else 0.0
    failed_positions = [value[0] for _, value in mapped]
    exact_reference_positions = sorted(
        position for position, (_, match) in alignment["mapping"].items() if match
    )
    matched_span = (
        f"{exact_reference_positions[0]}-{exact_reference_positions[-1]}"
        if exact_reference_positions else ""
    )
    before = sum(position < start for position in exact_reference_positions)
    inside = sum(start <= position <= end for position in exact_reference_positions)
    after = sum(position > end for position in exact_reference_positions)
    if covered >= 30 and coverage >= 80 and identity >= 80:
        relationship = "CONTAINS_DBD_COMPLETELY"
    elif identity >= 80 and covered >= 10 and 10 < coverage < 80:
        relationship = "PARTIALLY_CONTAINS_DBD"
    elif coverage <= 10:
        if before >= 20 and not after:
            relationship = "ENDS_BEFORE_DBD"
        elif after >= 20 and not before:
            relationship = "STARTS_AFTER_DBD"
        elif before >= 20 and after >= 20 and inside < 10:
            relationship = "INTERNAL_DBD_DELETION"
        else:
            relationship = "LACKS_DBD_UNRESOLVED_TOPOLOGY"
    elif coverage >= 80 and identity < 80:
        relationship = "COVERS_DBD_COORDINATE_SPAN_LOW_IDENTITY"
    else:
        relationship = "PARTIAL_OR_DIVERGENT_DBD_UNRESOLVED"
    non_dbd = [
        value for position, value in alignment["mapping"].items()
        if position < start or position > end
    ]
    non_matches = sum(value[1] for value in non_dbd)
    return {
        "length": dbd_length, "covered": covered, "matches": matches,
        "coverage": coverage, "identity": identity,
        "failed_coords": (
            f"{min(failed_positions)}-{max(failed_positions)}"
            if failed_positions else ""
        ),
        "matched_span": matched_span, "relationship": relationship,
        "non_dbd_aligned": len(non_dbd),
        "non_dbd_identity": (
            100.0 * non_matches / len(non_dbd) if non_dbd else 0.0
        ),
    }


def add_check(
    checks: dict[str, Any], name: str, passed: bool,
    observed: Any, expected: Any,
) -> None:
    checks[name] = {
        "passed": bool(passed), "observed": observed, "expected": expected,
    }


def main() -> int:
    args = parse_args()
    inputs = [
        args.manual_review, args.structural_draft, args.reference_alignments,
        args.domain_inventory, args.reviewed_pfam_map, args.failed_fasta,
        args.same_gene_candidates, args.same_gene_audit, args.database,
        args.all_tf_fasta,
    ]
    missing = [str(path) for path in inputs + [args.models_root] if not path.exists()]
    if missing:
        raise SystemExit("Missing required input(s): " + ", ".join(missing))
    outputs = [args.output_dossier, args.output_alignments, args.output_qc]
    existing = [str(path) for path in outputs if path.exists()]
    if existing:
        raise SystemExit("Refusing to overwrite existing output(s): " + ", ".join(existing))
    for path in outputs:
        path.parent.mkdir(parents=True, exist_ok=True)

    before = db_signature(args.database)
    input_hashes = {str(path): file_sha256(path) for path in inputs}
    manual = read_tsv(args.manual_review)
    structural = {r["failed_accession"]: r for r in read_tsv(args.structural_draft)}
    read_tsv(args.reference_alignments)
    read_tsv(args.domain_inventory)
    read_tsv(args.reviewed_pfam_map)
    candidates = read_tsv(args.same_gene_candidates)
    audits = read_tsv(args.same_gene_audit)
    failed_records = fasta_records(args.failed_fasta)
    all_records = fasta_records(args.all_tf_fasta)
    candidate_by_failed: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in candidates:
        candidate_by_failed[row["failed_tf_id"]].append(row)
    audit_by_accession = {r["best_model_accession"]: r for r in audits}

    if len(manual) != 36 or len({r["failed_accession"] for r in manual}) != 36:
        raise SystemExit("Manual-review input is not exactly 36 unique accessions")

    connection = sqlite3.connect(f"file:{args.database.resolve()}?mode=ro", uri=True)
    connection.execute("PRAGMA query_only=ON")
    dossier: list[dict[str, Any]] = []
    alignments: list[dict[str, Any]] = []
    selected_pdbs: set[Path] = set()
    try:
        for source in manual:
            accession = source["failed_accession"]
            failed_sequence = failed_records[accession]["sequence"]
            gene = source["gene_names"].strip()
            is_fusion = (
                accession in FUSION_CONTROLS
                or "fusion" in source["protein_description"].lower()
                or "fusion" in source["gene_mapping_source"].lower()
            )
            unresolved = accession in UNRESOLVED_CONTROLS
            if is_fusion:
                gene_status = "FUSION_MULTIPLE_COMPONENTS"
                primary_gene = ""
                gene_evidence = source["gene_mapping_source"]
            elif unresolved or not gene:
                gene_status = "UNRESOLVED"
                primary_gene = ""
                gene_evidence = source["gene_mapping_source"]
            else:
                gene_status = "RESOLVED_SINGLE_GENE"
                primary_gene = gene.split(";")[0].split()[0]
                gene_evidence = source["gene_mapping_source"]

            canonical_accession = ""
            canonical_genes = ""
            canonical_length: int | str = ""
            canonical_reason = ""
            canonical_sequence = ""
            canonical_meta = None
            canonical_domains_list: list[tuple[str, int, int, str]] = []
            if primary_gene:
                selected = select_reviewed_canonical(connection, primary_gene)
                if selected:
                    canonical_accession, canonical_genes, canonical_length = selected
                    canonical_reason = (
                        f"REVIEWED_UNIPROT_EXACT_GENE_TOKEN:{primary_gene}"
                    )
                    canonical_sequence = all_records.get(
                        canonical_accession, {}
                    ).get("sequence", "")
                    canonical_meta, canonical_domains_list = canonical_domains(
                        connection, canonical_accession
                    )

            model_accession = source["canonical_reference_accession"]
            model_path_text = source["canonical_reference_model_path"]
            model_reason = source["canonical_reference_selection_reason"]
            if not model_accession and primary_gene:
                modeled = [
                    row for row in candidate_by_failed[accession]
                    if int(row.get("total_pdb_count") or 0) > 0
                    and row.get("example_pdb")
                ]
                if modeled:
                    chosen_candidate = sorted(
                        modeled,
                        key=lambda row: (
                            -int(row.get("total_pdb_count") or 0),
                            row["candidate_tf_id"],
                        ),
                    )[0]
                    model_accession = chosen_candidate["candidate_tf_id"]
                    model_path_text = chosen_candidate["example_pdb"]
                    model_reason = "SAME_GENE_CANDIDATE_WITH_EXISTING_MODCRE_MODEL"
            model_meta = None
            model_domains: list[tuple[str, int, int, str]] = []
            if model_accession:
                model_meta, model_domains = canonical_domains(
                    connection, model_accession
                )
            if not model_accession and canonical_accession:
                paths = sorted(args.models_root.glob(f"TFS_{canonical_accession}:*.pdb"))
                chosen_path = choose_model_path(paths, canonical_domains_list)
                if chosen_path:
                    model_accession = canonical_accession
                    model_path_text = str(chosen_path)
                    model_reason = "REVIEWED_CANONICAL_ACCESSION_HAS_LOCAL_MODCRE_MODEL"
                    model_meta, model_domains = canonical_domains(
                        connection, model_accession
                    )
            if model_accession and not model_path_text:
                paths = sorted(args.models_root.glob(f"TFS_{model_accession}:*.pdb"))
                chosen_path = choose_model_path(paths, model_domains)
                if chosen_path:
                    model_path_text = str(chosen_path)

            dbd_domains, dbd_interval = canonical_dbd(canonical_domains_list)
            dbd_type = ";".join(
                f"{pfam}:{name}" for pfam, _, _, name in dbd_domains
            )
            dbd_pfams = ";".join(dict.fromkeys(d[0] for d in dbd_domains))
            alignment = None
            metrics = topology_and_metrics(
                {"mapping": {}}, None, int(canonical_length or 0)
            )
            if canonical_accession and canonical_sequence and dbd_interval:
                alignment = global_alignment(failed_sequence, canonical_sequence)
                metrics = topology_and_metrics(
                    alignment, dbd_interval, len(canonical_sequence)
                )

            model_sequence = all_records.get(model_accession, {}).get(
                "sequence", ""
            )
            model_interval = (
                parse_model_interval(model_path_text) if model_path_text else None
            )
            reconciliation = {
                "chain": "", "sequence": "", "aligned": 0, "identity": 0.0,
                "coverage": 0.0, "status": "NO_MODEL_PDB", "note": "",
            }
            pdb_dbd = {
                "aligned_target": 0, "target_coverage_percent": 0.0,
                "identity_percent": 0.0,
            }
            if model_path_text:
                model_path = Path(model_path_text)
                if not model_path.exists():
                    raise SystemExit(f"Selected PDB does not exist: {model_path}")
                if not model_sequence:
                    raise SystemExit(
                        f"Complete model-bearing sequence missing: {model_accession}"
                    )
                selected_pdbs.add(model_path)
                reconciliation = reconcile_pdb(
                    model_path, model_sequence, model_interval
                )
                if reconciliation["sequence"] and dbd_interval and canonical_sequence:
                    canonical_dbd_sequence = canonical_sequence[
                        dbd_interval[0] - 1:dbd_interval[1]
                    ]
                    pdb_dbd = local_alignment_metrics(
                        reconciliation["sequence"], canonical_dbd_sequence
                    )
            pdb_represents = "NO_MODEL_PDB"
            if reconciliation["sequence"] and dbd_interval:
                if (
                    pdb_dbd["aligned_target"] >= 30
                    and pdb_dbd["target_coverage_percent"] >= 80
                    and pdb_dbd["identity_percent"] >= 80
                ):
                    pdb_represents = "YES"
                elif pdb_dbd["aligned_target"]:
                    pdb_represents = "PARTIAL_OR_DIVERGENT"
                else:
                    pdb_represents = "NO"
            elif reconciliation["sequence"]:
                pdb_represents = "CANONICAL_DBD_UNDEFINED"

            high_present = (
                metrics["covered"] >= 30
                and metrics["coverage"] >= 80
                and metrics["identity"] >= 80
            )
            high_absent = (
                metrics["coverage"] <= 10
                and metrics["non_dbd_aligned"] >= 40
                and metrics["non_dbd_identity"] >= 80
            )
            if is_fusion:
                decision = "FUSION_COMPONENT_REQUIRES_REVIEW"
                action = "REQUIRES_MANUAL_REVIEW"
                evidence = "Fusion components cannot be collapsed to one canonical accession."
            elif unresolved:
                decision = "GENE_MAPPING_UNRESOLVED"
                action = "REQUIRES_MANUAL_REVIEW"
                evidence = "Primary gene is unresolved; canonical selection is intentionally blank."
            elif not canonical_accession or not dbd_interval:
                decision = "INSUFFICIENT_EVIDENCE"
                action = "REQUIRES_MANUAL_REVIEW"
                evidence = "No reviewed canonical accession with supported DBD coordinates was resolved."
            elif high_present:
                decision = "LIKELY_DBD_PRESENT"
                action = "LIKELY_SEND_NEW_DBD_FRAGMENT"
                evidence = "Canonical DBD coverage, identity, and aligned-residue thresholds all pass."
            elif (
                metrics["relationship"] == "PARTIALLY_CONTAINS_DBD"
            ):
                decision = "PARTIAL_DBD_PRESENT"
                action = "REQUIRES_MANUAL_REVIEW"
                evidence = "A high-identity portion of the canonical DBD is present, but coverage is below 80%."
            elif high_absent:
                decision = "LIKELY_ACCESSION_LACKS_DBD"
                audit = audit_by_accession.get(model_accession)
                action = (
                    "LIKELY_NO_SEND_CANONICAL_REFERENCE_ONLY"
                    if audit and audit.get("result") == "PASS_NONEMPTY_PROTEIN_DNA_INTERFACE"
                    else "LIKELY_NO_SEND_ACCESSION_LACKS_DBD"
                )
                evidence = "DBD coverage is at most 10%; at least 40 non-DBD residues align at at least 80% identity."
            else:
                decision = "INSUFFICIENT_EVIDENCE"
                action = "REQUIRES_MANUAL_REVIEW"
                evidence = "Coverage/identity does not satisfy a conservative automatic rule."

            audit = audit_by_accession.get(model_accession, {})
            row = {
                "failed_accession": accession,
                "resolved_primary_gene": primary_gene,
                "gene_resolution_status": gene_status,
                "gene_resolution_evidence": gene_evidence,
                "failed_protein_description": source["protein_description"],
                "failed_sequence_length": len(failed_sequence),
                "failed_sequence_sha256": sequence_sha256(failed_sequence),
                "reviewed_canonical_accession": canonical_accession,
                "reviewed_canonical_gene_names": canonical_genes,
                "reviewed_canonical_sequence_length": canonical_length,
                "reviewed_canonical_selection_reason": canonical_reason,
                "model_bearing_accession": model_accession,
                "model_bearing_is_reviewed": (
                    "YES" if model_meta and model_meta["reviewed"]
                    else "NO" if model_meta else ""
                ),
                "model_bearing_sequence_length": len(model_sequence) if model_sequence else "",
                "model_bearing_selection_reason": model_reason,
                "model_path": model_path_text,
                "model_interval": (
                    f"{model_interval[0]}-{model_interval[1]}"
                    if model_interval else ""
                ),
                "existing_interface_result": audit.get(
                    "result", source["canonical_reference_interface_result"]
                ),
                "canonical_DBD_type": dbd_type,
                "canonical_DBD_pfam_ids": dbd_pfams,
                "canonical_DBD_start": dbd_interval[0] if dbd_interval else "",
                "canonical_DBD_end": dbd_interval[1] if dbd_interval else "",
                "canonical_DBD_length": metrics["length"],
                "canonical_DBD_coordinate_source": (
                    "SQLITE_REVIEWED_CANONICAL_PFAM_ARCHITECTURE"
                    if dbd_interval else ""
                ),
                "failed_to_canonical_alignment_columns": (
                    alignment["columns"] if alignment else ""
                ),
                "failed_to_canonical_aligned_residues": (
                    alignment["aligned"] if alignment else ""
                ),
                "failed_to_canonical_identity_percent": (
                    f"{alignment['identity']:.2f}" if alignment else ""
                ),
                "canonical_DBD_aligned_residues": metrics["covered"],
                "canonical_DBD_coverage_percent": f"{metrics['coverage']:.2f}",
                "canonical_DBD_identity_percent": f"{metrics['identity']:.2f}",
                "failed_coordinates_corresponding_to_DBD": metrics["failed_coords"],
                "failed_matched_canonical_span": metrics["matched_span"],
                "DBD_sequence_relationship": metrics["relationship"],
                "aligned_non_DBD_residues": metrics["non_dbd_aligned"],
                "aligned_non_DBD_identity_percent": f"{metrics['non_dbd_identity']:.2f}",
                "pdb_protein_chain": reconciliation["chain"],
                "pdb_protein_sequence": reconciliation["sequence"],
                "pdb_protein_sequence_length": len(reconciliation["sequence"]),
                "pdb_protein_sequence_sha256": (
                    sequence_sha256(reconciliation["sequence"])
                    if reconciliation["sequence"] else ""
                ),
                "pdb_to_model_accession_aligned_residues": reconciliation["aligned"],
                "pdb_to_model_accession_identity_percent": f"{reconciliation['identity']:.2f}",
                "pdb_to_model_accession_coverage_percent": f"{reconciliation['coverage']:.2f}",
                "pdb_model_accession_reconciliation_status": reconciliation["status"],
                "pdb_to_canonical_DBD_aligned_residues": pdb_dbd["aligned_target"],
                "pdb_to_canonical_DBD_coverage_percent": f"{pdb_dbd['target_coverage_percent']:.2f}",
                "pdb_to_canonical_DBD_identity_percent": f"{pdb_dbd['identity_percent']:.2f}",
                "pdb_represents_canonical_DBD": pdb_represents,
                "canonical_and_model_accessions_separate": (
                    "YES" if canonical_accession and model_accession
                    and canonical_accession != model_accession
                    else "NO" if canonical_accession and model_accession
                    else "NOT_APPLICABLE"
                ),
                "suggested_decision": decision,
                "suggested_baldo_action": action,
                "decision_evidence": evidence,
                "manual_review_notes": "",
            }
            dossier.append(row)

            if alignment and dbd_interval:
                mask = "".join(
                    "^" if rpos and dbd_interval[0] <= rpos <= dbd_interval[1]
                    else " "
                    for _, rpos, _ in alignment["pairs"]
                )
                mapping_text = ";".join(
                    f"{r}:{value[0]}" for r, value in sorted(
                        alignment["mapping"].items()
                    ) if dbd_interval[0] <= r <= dbd_interval[1]
                )
                matched_positions = ";".join(
                    str(position) for position, (_, match)
                    in sorted(alignment["mapping"].items()) if match
                )
                alignments.append({
                    "failed_accession": accession,
                    "reviewed_canonical_accession": canonical_accession,
                    "model_bearing_accession": model_accession,
                    "canonical_DBD_start": dbd_interval[0],
                    "canonical_DBD_end": dbd_interval[1],
                    "failed_alignment": alignment["failed_alignment"],
                    "canonical_alignment": alignment["canonical_alignment"],
                    "canonical_DBD_alignment_mask": mask,
                    "reference_to_failed_coordinate_map": mapping_text,
                    "failed_matched_canonical_positions": matched_positions,
                    "DBD_sequence_relationship": metrics["relationship"],
                })
    finally:
        connection.close()

    after_analysis = db_signature(args.database)
    by_id = {r["failed_accession"]: r for r in dossier}
    checks: dict[str, Any] = {}
    add_check(checks, "exactly_36_unique_accessions",
              len(dossier) == len({r["failed_accession"] for r in dossier}) == 36,
              {"rows": len(dossier), "unique": len({r["failed_accession"] for r in dossier})},
              {"rows": 36, "unique": 36})
    add_check(checks, "one_dossier_row_per_input_accession",
              {r["failed_accession"] for r in dossier} == {r["failed_accession"] for r in manual},
              sorted(r["failed_accession"] for r in dossier),
              sorted(r["failed_accession"] for r in manual))
    add_check(checks, "canonical_and_model_accession_columns_separate",
              all("reviewed_canonical_accession" in r and "model_bearing_accession" in r for r in dossier),
              True, True)
    forbidden = {
        "primary_structural_status", "database_display_recommendation",
        "final_structural_status", "final_database_status",
    } & set(DOSSIER_COLUMNS)
    add_check(checks, "no_final_structural_status_assigned", not forbidden, sorted(forbidden), [])
    add_check(checks, "only_allowed_suggestions",
              all(r["suggested_decision"] in DECISIONS and r["suggested_baldo_action"] in ACTIONS for r in dossier),
              True, True)
    invalid_send = [
        r["failed_accession"] for r in dossier
        if r["suggested_baldo_action"] in {
            "LIKELY_SEND_NEW_DBD_FRAGMENT",
            "LIKELY_SEND_CORRECTED_DBD_FRAGMENT",
        } and not (
            float(r["canonical_DBD_coverage_percent"]) >= 80
            and float(r["canonical_DBD_identity_percent"]) >= 80
            and int(r["canonical_DBD_aligned_residues"]) >= 30
        )
    ]
    add_check(checks, "no_send_suggestion_below_80_80_30", not invalid_send, invalid_send, [])
    add_check(checks, "database_unchanged_after_analysis", before == after_analysis, after_analysis, before)

    add_check(checks, "control_B4DNX4_EGR1_topology_determined",
              by_id["B4DNX4"]["reviewed_canonical_accession"] == "P18146"
              and by_id["B4DNX4"]["DBD_sequence_relationship"] in {
                  "ENDS_BEFORE_DBD", "STARTS_AFTER_DBD",
                  "INTERNAL_DBD_DELETION", "LACKS_DBD_UNRESOLVED_TOPOLOGY",
              },
              {k: by_id["B4DNX4"][k] for k in [
                  "reviewed_canonical_accession", "canonical_DBD_type",
                  "DBD_sequence_relationship", "failed_matched_canonical_span",
              ]}, "EGR1 canonical C2H2 topology explicitly determined")
    for accession in ["Q59EF3", "E9PKB7"]:
        add_check(checks, f"control_{accession}_TEAD1_PF17725_not_DBD",
                  by_id[accession]["reviewed_canonical_accession"] == "P28347"
                  and "PF01285" in by_id[accession]["canonical_DBD_pfam_ids"]
                  and "PF17725" not in by_id[accession]["canonical_DBD_pfam_ids"],
                  {k: by_id[accession][k] for k in [
                      "reviewed_canonical_accession", "model_bearing_accession",
                      "canonical_DBD_pfam_ids",
                  ]}, {"canonical": "P28347", "DBD": "PF01285", "not_DBD": "PF17725"})
    for accession in ["Q5VVR5", "X5D8U0", "H0YNG1"]:
        add_check(checks, f"control_{accession}_high_coverage_low_identity_visible",
                  bool(by_id[accession]["canonical_DBD_coverage_percent"])
                  and bool(by_id[accession]["canonical_DBD_identity_percent"])
                  and by_id[accession]["suggested_baldo_action"] != "LIKELY_SEND_NEW_DBD_FRAGMENT"
                  if float(by_id[accession]["canonical_DBD_identity_percent"]) < 80 else True,
                  {k: by_id[accession][k] for k in [
                      "canonical_DBD_coverage_percent",
                      "canonical_DBD_identity_percent",
                      "DBD_sequence_relationship", "suggested_baldo_action",
                  ]}, "low identity cannot trigger likely send")
    for accession in ["A0A804HLH1", "Q6L9M1"]:
        add_check(checks, f"control_{accession}_reviewed_PPARG",
                  by_id[accession]["reviewed_canonical_accession"] == "P37231",
                  by_id[accession]["reviewed_canonical_accession"], "P37231")
    add_check(checks, "control_all_fusions_flagged",
              all(by_id[a]["suggested_decision"] == "FUSION_COMPONENT_REQUIRES_REVIEW" for a in FUSION_CONTROLS),
              {a: by_id[a]["suggested_decision"] for a in sorted(FUSION_CONTROLS)},
              "all fusion controls flagged")
    add_check(checks, "control_all_uncertain_genes_unresolved",
              all(by_id[a]["suggested_decision"] == "GENE_MAPPING_UNRESOLVED" for a in UNRESOLVED_CONTROLS),
              {a: by_id[a]["suggested_decision"] for a in sorted(UNRESOLVED_CONTROLS)},
              "all unresolved-gene controls flagged")

    failed_checks = [name for name, item in checks.items() if not item["passed"]]
    if failed_checks:
        print("FAILED_QC_CHECKS")
        for name in failed_checks:
            print(name, json.dumps(checks[name], sort_keys=True), sep="\t")
        raise SystemExit("QC failed before output write")

    write_tsv(args.output_dossier, dossier, DOSSIER_COLUMNS)
    write_tsv(args.output_alignments, alignments, ALIGNMENT_COLUMNS)
    after_write = db_signature(args.database)
    add_check(checks, "database_unchanged_after_output_write", before == after_write, after_write, before)
    output_hashes = {
        str(path): file_sha256(path)
        for path in [args.output_dossier, args.output_alignments]
    }
    qc = {
        "all_checks_passed": all(item["passed"] for item in checks.values()),
        "read_only_dossier": True,
        "counts": {
            "suggested_decision": dict(sorted(Counter(r["suggested_decision"] for r in dossier).items())),
            "suggested_baldo_action": dict(sorted(Counter(r["suggested_baldo_action"] for r in dossier).items())),
            "alignments": len(alignments),
        },
        "likely_send_accessions": sorted(
            r["failed_accession"] for r in dossier
            if r["suggested_baldo_action"] in {
                "LIKELY_SEND_NEW_DBD_FRAGMENT",
                "LIKELY_SEND_CORRECTED_DBD_FRAGMENT",
            }
        ),
        "fusion_accessions": sorted(
            r["failed_accession"] for r in dossier
            if r["suggested_decision"] == "FUSION_COMPONENT_REQUIRES_REVIEW"
        ),
        "unresolved_gene_accessions": sorted(
            r["failed_accession"] for r in dossier
            if r["suggested_decision"] == "GENE_MAPPING_UNRESOLVED"
        ),
        "database_signature_before": before,
        "database_signature_after": after_write,
        "input_hashes": input_hashes,
        "selected_pdb_hashes": {
            str(path): file_sha256(path) for path in sorted(selected_pdbs)
        },
        "output_hashes_excluding_qc_self_hash": output_hashes,
        "checks": checks,
    }
    with args.output_qc.open("x", encoding="utf-8") as handle:
        json.dump(qc, handle, indent=2, sort_keys=True)
        handle.write("\n")
    if not qc["all_checks_passed"]:
        raise SystemExit("Final QC failed")
    print(f"Wrote {args.output_dossier}")
    print(f"Wrote {args.output_alignments}")
    print(f"Wrote {args.output_qc}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
