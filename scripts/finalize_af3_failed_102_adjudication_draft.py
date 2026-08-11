#!/usr/bin/env python3
"""Finalize the 102-accession structural adjudication as a DRAFT package."""

from __future__ import annotations

import argparse
import csv
import difflib
import hashlib
import json
import re
import sqlite3
import sys
import urllib.request
from collections import Counter
from pathlib import Path
from typing import Any

from build_af3_failed_102_manual_review_dossier import (
    canonical_dbd,
    global_alignment,
    select_reviewed_canonical,
    topology_and_metrics,
)
from classify_af3_failed_102_structural_status import (
    CANONICAL_DBD_PFAMS,
    canonical_domains,
    db_signature,
    fasta_records,
    file_sha256,
)


ROOT = Path(__file__).resolve().parents[1]
CANONICAL_DBD_PFAMS.add("PF00853")  # RUNX Runt DNA-binding domain

STATUSES = {
    "VALID_DBD_FRAGMENT_MODEL",
    "DBD_PRESENT_FRAGMENT_FAILED",
    "DBD_PRESENT_FRAGMENT_NOT_SENT",
    "ACCESSION_LACKS_DBD",
    "ACCESSION_CONTAINS_INCOMPLETE_DBD",
    "COFACTOR_OR_NON_DNA_BINDING_COMPONENT",
    "UNCERTAIN_MANUAL_REVIEW",
}
ACTIONS = {
    "SEND_NEW_DBD_FRAGMENT",
    "SEND_CORRECTED_DBD_FRAGMENT",
    "NO_SEND_USE_EXISTING_DBD_FRAGMENT",
    "NO_SEND_CANONICAL_REFERENCE_ONLY",
    "NO_SEND_MARK_DNA_STRUCTURE_UNAVAILABLE",
    "MANUAL_REVIEW_BEFORE_SENDING",
}
FUSION_PARENTS = {
    "A0A0A7M2K0": [("ZFPM2", "Q8WW38"), ("ELF5", "Q9UKW6")],
    "B1NY96": [("TCF3/E2A", "P15923"), ("ZNF384/NMP4", "Q8TF68")],
    "D1LUZ8": [("NOL4L/C20orf112", "Q96MY1"), ("RUNX1", "Q01196")],
    "H2BNB9": [("ETV6", "P41212"), ("INO80D", "Q53TQ3")],
}
UNRESOLVED = {"A8K549", "B4DTN3", "Q68D60"}
PARTIAL_CONTROLS = {"B4DL38", "B7Z8F5", "F6X9D6", "Q5VVR5"}

FINAL_COLUMNS = [
    "failed_accession", "gene", "protein_description",
    "failed_sequence_length", "failed_sequence_sha256",
    "DBD_disposition", "expected_DBD_type",
    "expected_DBD_canonical_start", "expected_DBD_canonical_end",
    "DBD_aligned_residues", "DBD_coverage_percent", "DBD_identity_percent",
    "failed_coordinates_corresponding_to_DBD",
    "reviewed_canonical_accession", "model_bearing_accession",
    "actual_accession_represented_by_PDB", "canonical_reference_model",
    "canonical_reference_interface_status", "valid_existing_fragment_model",
    "original_fragment_status", "primary_structural_status",
    "action_for_baldo", "send_coordinates",
    "database_display_recommendation", "decision_reason",
    "remaining_uncertainty", "fusion_parent_analysis",
    "local_sequence_search_results", "draft_status",
]

BALDO_COLUMNS = [
    "failed_accession", "gene", "protein_description",
    "primary_structural_status", "DBD_evidence", "fragment_status",
    "canonical_reference_accession", "canonical_reference_model",
    "canonical_reference_interface_status", "action_for_baldo",
    "sequence_or_fragment_to_send", "send_coordinates",
    "database_display_recommendation", "decision_reason",
    "remaining_uncertainty",
]

DISPLAY_COLUMNS = [
    "failed_accession", "gene", "primary_structural_status",
    "canonical_reference_accession", "canonical_reference_is_separate",
    "database_display_recommendation", "display_rationale", "draft_only",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create the final reviewed 102-accession DRAFT adjudication.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ext = ROOT / "external/baldo_model_inventory"
    parser.add_argument("--structural-draft", type=Path, default=ROOT / "outputs/af3_failed_102_structural_status_draft.tsv")
    parser.add_argument("--manual-dossier", type=Path, default=ROOT / "outputs/af3_failed_102_manual_review_dossier.tsv")
    parser.add_argument("--manual-alignments", type=Path, default=ROOT / "outputs/af3_failed_102_manual_review_alignments.tsv")
    parser.add_argument("--reference-alignments", type=Path, default=ROOT / "outputs/af3_failed_102_reference_alignments.tsv")
    parser.add_argument("--domain-inventory", type=Path, default=ROOT / "outputs/af3_failed_102_domain_inventory.tsv")
    parser.add_argument("--reviewed-pfam-map", type=Path, default=ROOT / "data_sources/af3_failed_102_pfam_role_reviewed.tsv")
    parser.add_argument("--failed-fasta", type=Path, default=ext / "af3_failed_102_full_length.fasta")
    parser.add_argument("--fragment-summary", type=Path, default=ext / "af3_fragment_interface_summary.tsv")
    parser.add_argument("--same-gene-audit", type=Path, default=ext / "af3_failed_same_gene_modcre_interface_audit.tsv")
    parser.add_argument("--database", type=Path, default=ROOT / "data/tf_webdb.sqlite")
    parser.add_argument("--all-tf-fasta", type=Path, default=Path("/home/patricia/TF_database_Baldo_data/TF_without_model.fasta"))
    parser.add_argument("--models-root", type=Path, default=Path("/data/sbi/interchange/boliva/patricia/models"))
    parser.add_argument("--output-final", type=Path, default=ROOT / "outputs/af3_failed_102_final_adjudication_DRAFT.tsv")
    parser.add_argument("--output-manual", type=Path, default=ROOT / "outputs/af3_failed_102_final_manual_review_DRAFT.tsv")
    parser.add_argument("--output-baldo", type=Path, default=ROOT / "outputs/af3_failed_102_baldo_package_DRAFT.tsv")
    parser.add_argument("--output-fasta", type=Path, default=ROOT / "outputs/af3_failed_102_baldo_fragments_DRAFT.fasta")
    parser.add_argument("--output-display", type=Path, default=ROOT / "outputs/af3_failed_102_database_display_DRAFT.tsv")
    parser.add_argument("--output-qc", type=Path, default=ROOT / "outputs/af3_failed_102_final_adjudication_qc.json")
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


def seq_hash(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def retrieve_uniprot_sequence(accession: str) -> tuple[str, str]:
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
    request = urllib.request.Request(
        url, headers={"User-Agent": "tf-webdb-final-adjudication-draft/1.0"}
    )
    with urllib.request.urlopen(request, timeout=30) as response:
        text = response.read().decode("utf-8")
    sequence = "".join(
        line.strip() for line in text.splitlines() if not line.startswith(">")
    )
    if not sequence:
        raise RuntimeError(f"Official UniProt sequence is empty: {accession}")
    return sequence, url


def kmer_score(query: str, target: str, k: int = 5) -> int:
    query_kmers = {query[i:i + k] for i in range(max(0, len(query) - k + 1))}
    target_kmers = {target[i:i + k] for i in range(max(0, len(target) - k + 1))}
    return len(query_kmers & target_kmers)


def gene_from_header(header: str) -> str:
    match = re.search(r"\|gene=([^|]*)", header)
    return match.group(1).split(",")[0] if match and match.group(1) else ""


def local_tf_search(
    query_accession: str,
    query: str,
    all_records: dict[str, dict[str, str]],
    dossier_gene_by_accession: dict[str, str],
) -> list[dict[str, Any]]:
    prefiltered = sorted(
        (
            (kmer_score(query, record["sequence"]), accession, record)
            for accession, record in all_records.items()
            if accession != query_accession
        ),
        reverse=True,
    )[:30]
    results = []
    for _, accession, record in prefiltered:
        alignment = global_alignment(query, record["sequence"])
        query_coverage = 100.0 * alignment["aligned"] / len(query)
        inferred_gene = (
            gene_from_header(record["header"])
            or dossier_gene_by_accession.get(accession, "")
        )
        results.append({
            "accession": accession,
            "gene": inferred_gene,
            "aligned_residues": alignment["aligned"],
            "identity_percent": round(alignment["identity"], 2),
            "query_coverage_percent": round(query_coverage, 2),
            "sequence_length": len(record["sequence"]),
        })
    return sorted(
        results,
        key=lambda row: (
            row["identity_percent"] * row["query_coverage_percent"],
            row["aligned_residues"],
            row["accession"],
        ),
        reverse=True,
    )[:5]


def dbd_metrics_for(
    connection: sqlite3.Connection,
    failed_sequence: str,
    canonical_accession: str,
    canonical_sequence: str,
) -> tuple[str, tuple[int, int] | None, dict[str, Any]]:
    _, domains = canonical_domains(connection, canonical_accession)
    dbds, interval = canonical_dbd(domains)
    dbd_type = ";".join(
        f"{pfam}:{name}" for pfam, _, _, name in dbds
    )
    if not interval:
        return dbd_type, None, topology_and_metrics(
            {"mapping": {}}, None, len(canonical_sequence)
        )
    alignment = global_alignment(failed_sequence, canonical_sequence)
    return dbd_type, interval, topology_and_metrics(
        alignment, interval, len(canonical_sequence)
    )


def high_present(metrics: dict[str, Any]) -> bool:
    return (
        metrics["coverage"] >= 80
        and metrics["identity"] >= 80
        and metrics["covered"] >= 30
    )


def high_absent(metrics: dict[str, Any]) -> bool:
    return (
        metrics["coverage"] <= 10
        and metrics["non_dbd_aligned"] >= 40
        and metrics["non_dbd_identity"] >= 80
    )


def partial_present(metrics: dict[str, Any]) -> bool:
    return (
        not high_present(metrics)
        and metrics["covered"] >= 10
        and metrics["identity"] >= 80
    )


def display_recommendation(
    status: str, canonical: str
) -> str:
    if status == "VALID_DBD_FRAGMENT_MODEL":
        return "Use validated DBD fragment model."
    if status == "ACCESSION_CONTAINS_INCOMPLETE_DBD":
        return "Accession contains only an incomplete DBD."
    if status == "COFACTOR_OR_NON_DNA_BINDING_COMPONENT":
        return "Non-DNA-binding cofactor/component; do not display as a DNA-bound TF structure."
    if status == "ACCESSION_LACKS_DBD" and canonical:
        return (
            f"Canonical same-gene reference available from accession {canonical}; "
            "label explicitly as a reference."
        )
    if status == "UNCERTAIN_MANUAL_REVIEW":
        return "Manual review unresolved."
    return "DNA-interacting structure unavailable for this accession."


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
        args.structural_draft, args.manual_dossier, args.manual_alignments,
        args.reference_alignments, args.domain_inventory,
        args.reviewed_pfam_map, args.failed_fasta, args.fragment_summary,
        args.same_gene_audit, args.database, args.all_tf_fasta,
    ]
    missing = [str(path) for path in inputs + [args.models_root] if not path.exists()]
    if missing:
        raise SystemExit("Missing required input(s): " + ", ".join(missing))
    outputs = [
        args.output_final, args.output_manual, args.output_baldo,
        args.output_fasta, args.output_display, args.output_qc,
    ]
    existing = [str(path) for path in outputs if path.exists()]
    if existing:
        raise SystemExit("Refusing to overwrite existing output(s): " + ", ".join(existing))
    for path in outputs:
        path.parent.mkdir(parents=True, exist_ok=True)

    before = db_signature(args.database)
    input_hashes = {str(path): file_sha256(path) for path in inputs}
    structural_rows = read_tsv(args.structural_draft)
    structural = {r["failed_accession"]: r for r in structural_rows}
    dossier_rows = read_tsv(args.manual_dossier)
    dossier = {r["failed_accession"]: r for r in dossier_rows}
    read_tsv(args.manual_alignments)
    read_tsv(args.reference_alignments)
    read_tsv(args.domain_inventory)
    read_tsv(args.reviewed_pfam_map)
    read_tsv(args.fragment_summary)
    interface_audits = {
        r["best_model_accession"]: r for r in read_tsv(args.same_gene_audit)
    }
    failed_records = fasta_records(args.failed_fasta)
    all_records = fasta_records(args.all_tf_fasta)
    if len(structural) != 102 or len(dossier) != 36:
        raise SystemExit("Expected 102 structural rows and 36 dossier rows")

    official_sequences: dict[str, dict[str, str]] = {}
    fusion_parent_accessions = {
        parent_accession
        for parents in FUSION_PARENTS.values()
        for _, parent_accession in parents
    }
    for accession in sorted(fusion_parent_accessions):
        if accession not in all_records:
            sequence, url = retrieve_uniprot_sequence(accession)
            official_sequences[accession] = {
                "sequence": sequence, "url": url, "sha256": seq_hash(sequence),
            }

    dossier_gene_by_accession = {
        accession: row["resolved_primary_gene"]
        for accession, row in dossier.items()
        if row["resolved_primary_gene"]
    }
    local_searches: dict[str, list[dict[str, Any]]] = {}
    for accession in sorted(UNRESOLVED):
        local_searches[accession] = local_tf_search(
            accession, failed_records[accession]["sequence"],
            all_records, dossier_gene_by_accession,
        )

    connection = sqlite3.connect(f"file:{args.database.resolve()}?mode=ro", uri=True)
    connection.execute("PRAGMA query_only=ON")
    fusion_results: dict[str, list[dict[str, Any]]] = {}
    refined: dict[str, dict[str, Any]] = {}
    try:
        for accession, parents in FUSION_PARENTS.items():
            failed_sequence = failed_records[accession]["sequence"]
            parent_rows = []
            for parent_gene, parent_accession in parents:
                parent_sequence = all_records.get(parent_accession, {}).get(
                    "sequence", ""
                ) or official_sequences.get(parent_accession, {}).get("sequence", "")
                if not parent_sequence:
                    raise SystemExit(f"Missing fusion parent sequence: {parent_accession}")
                dbd_type, interval, metrics = dbd_metrics_for(
                    connection, failed_sequence, parent_accession, parent_sequence
                )
                parent_rows.append({
                    "parent_gene": parent_gene,
                    "parent_accession": parent_accession,
                    "parent_sequence_length": len(parent_sequence),
                    "expected_DBD_type": dbd_type,
                    "DBD_start": interval[0] if interval else "",
                    "DBD_end": interval[1] if interval else "",
                    "DBD_aligned_residues": metrics["covered"],
                    "DBD_coverage_percent": round(metrics["coverage"], 2),
                    "DBD_identity_percent": round(metrics["identity"], 2),
                    "failed_coordinates": metrics["failed_coords"],
                    "relationship": metrics["relationship"],
                    "matched_parent_span": metrics["matched_span"],
                    "complete_DBD": high_present(metrics),
                    "partial_DBD": partial_present(metrics),
                })
            fusion_results[accession] = parent_rows

        for accession in sorted(UNRESOLVED):
            search = local_searches[accession]
            best = search[0]
            passes = (
                best["aligned_residues"] >= 50
                and best["identity_percent"] >= 80
                and best["query_coverage_percent"] >= 50
                and bool(best["gene"])
            )
            if not passes:
                refined[accession] = {
                    "gene": "", "canonical": "", "dbd_type": "",
                    "interval": None, "metrics": topology_and_metrics(
                        {"mapping": {}}, None, len(failed_records[accession]["sequence"])
                    ),
                    "search_passed": False,
                }
                continue
            gene = best["gene"].split(";")[0].split()[0]
            selected = select_reviewed_canonical(connection, gene)
            if not selected:
                refined[accession] = {
                    "gene": gene, "canonical": "", "dbd_type": "",
                    "interval": None, "metrics": topology_and_metrics(
                        {"mapping": {}}, None, len(failed_records[accession]["sequence"])
                    ),
                    "search_passed": True,
                }
                continue
            canonical_accession = selected[0]
            canonical_sequence = all_records[canonical_accession]["sequence"]
            dbd_type, interval, metrics = dbd_metrics_for(
                connection, failed_records[accession]["sequence"],
                canonical_accession, canonical_sequence,
            )
            refined[accession] = {
                "gene": gene, "canonical": canonical_accession,
                "dbd_type": dbd_type, "interval": interval,
                "metrics": metrics, "search_passed": True,
            }
    finally:
        connection.close()

    final_rows: list[dict[str, Any]] = []
    baldo_rows: list[dict[str, Any]] = []
    display_rows: list[dict[str, Any]] = []
    for accession in sorted(structural):
        base = structural[accession]
        failed_sequence = failed_records[accession]["sequence"]
        drow = dossier.get(accession)
        status = base["primary_structural_status"]
        action = base["action_for_baldo"]
        gene = base["gene_names"]
        canonical = base["canonical_reference_accession"]
        model_accession = canonical
        model_path = base["canonical_reference_model_path"]
        interface_status = base["canonical_reference_interface_result"]
        dbd_type = ""
        interval: tuple[int, int] | None = None
        metrics = {
            "covered": 0, "coverage": 0.0, "identity": 0.0,
            "failed_coords": "", "relationship": "NOT_REASSESSED_NONMANUAL",
        }
        fusion_json = ""
        search_json = ""
        remaining = ""
        decision = base["decision_reason"]

        if drow:
            gene = drow["resolved_primary_gene"] or gene
            canonical = drow["reviewed_canonical_accession"]
            model_accession = drow["model_bearing_accession"]
            model_path = drow["model_path"]
            interface_status = drow["existing_interface_result"]
            dbd_type = drow["canonical_DBD_type"]
            interval = (
                (int(drow["canonical_DBD_start"]), int(drow["canonical_DBD_end"]))
                if drow["canonical_DBD_start"] and drow["canonical_DBD_end"]
                else None
            )
            metrics = {
                "covered": int(drow["canonical_DBD_aligned_residues"] or 0),
                "coverage": float(drow["canonical_DBD_coverage_percent"] or 0),
                "identity": float(drow["canonical_DBD_identity_percent"] or 0),
                "failed_coords": drow["failed_coordinates_corresponding_to_DBD"],
                "relationship": drow["DBD_sequence_relationship"],
            }
            suggested = drow["suggested_decision"]
            if suggested == "LIKELY_ACCESSION_LACKS_DBD":
                status = "ACCESSION_LACKS_DBD"
                action = (
                    "NO_SEND_CANONICAL_REFERENCE_ONLY"
                    if model_accession in interface_audits
                    else "NO_SEND_MARK_DNA_STRUCTURE_UNAVAILABLE"
                )
                decision = drow["decision_evidence"]
            elif suggested == "PARTIAL_DBD_PRESENT":
                status = "ACCESSION_CONTAINS_INCOMPLETE_DBD"
                action = "NO_SEND_MARK_DNA_STRUCTURE_UNAVAILABLE"
                decision = (
                    "Only part of the expected canonical DBD is present; "
                    "the accession cannot provide a complete modelable DBD."
                )
            elif suggested == "INSUFFICIENT_EVIDENCE":
                status = "UNCERTAIN_MANUAL_REVIEW"
                action = "MANUAL_REVIEW_BEFORE_SENDING"
                remaining = drow["decision_evidence"]

        if accession in FUSION_PARENTS:
            parents = fusion_results[accession]
            fusion_json = json.dumps(parents, sort_keys=True, separators=(",", ":"))
            complete = [row for row in parents if row["complete_DBD"]]
            partial = [row for row in parents if row["partial_DBD"]]
            if complete:
                chosen = complete[0]
                status = "DBD_PRESENT_FRAGMENT_NOT_SENT"
                action = "SEND_NEW_DBD_FRAGMENT"
                dbd_type = chosen["expected_DBD_type"]
                interval = (
                    int(chosen["DBD_start"]), int(chosen["DBD_end"])
                )
                metrics = {
                    "covered": chosen["DBD_aligned_residues"],
                    "coverage": chosen["DBD_coverage_percent"],
                    "identity": chosen["DBD_identity_percent"],
                    "failed_coords": chosen["failed_coordinates"],
                    "relationship": chosen["relationship"],
                }
                canonical = chosen["parent_accession"]
                decision = "A complete high-confidence parent DBD is retained in the fusion."
                remaining = "Fusion breakpoint still requires explicit manual confirmation."
            elif partial:
                chosen = partial[0]
                status = "ACCESSION_CONTAINS_INCOMPLETE_DBD"
                action = "NO_SEND_MARK_DNA_STRUCTURE_UNAVAILABLE"
                dbd_type = chosen["expected_DBD_type"]
                metrics = {
                    "covered": chosen["DBD_aligned_residues"],
                    "coverage": chosen["DBD_coverage_percent"],
                    "identity": chosen["DBD_identity_percent"],
                    "failed_coords": chosen["failed_coordinates"],
                    "relationship": chosen["relationship"],
                }
                canonical = chosen["parent_accession"]
                decision = "Fusion retains only part of a parent DBD; no complete sendable DBD is demonstrated."
                remaining = "Fusion breakpoint and both parent contributions remain explicitly recorded."
            else:
                status = "UNCERTAIN_MANUAL_REVIEW"
                action = "MANUAL_REVIEW_BEFORE_SENDING"
                decision = "Parent-by-parent alignments do not demonstrate a complete high-confidence DBD."
                remaining = "Fusion structure and breakpoint remain unresolved; no TF-parent DBD is assumed."

        if accession in UNRESOLVED:
            result = refined[accession]
            search_json = json.dumps(
                local_searches[accession], sort_keys=True, separators=(",", ":")
            )
            if result["search_passed"]:
                gene = result["gene"]
                canonical = result["canonical"]
                dbd_type = result["dbd_type"]
                interval = result["interval"]
                metrics = result["metrics"]
                if high_present(metrics):
                    status = "DBD_PRESENT_FRAGMENT_NOT_SENT"
                    action = "SEND_NEW_DBD_FRAGMENT"
                    decision = "Local TF search resolves the gene and the failed sequence contains a complete high-confidence canonical DBD."
                elif partial_present(metrics):
                    status = "ACCESSION_CONTAINS_INCOMPLETE_DBD"
                    action = "NO_SEND_MARK_DNA_STRUCTURE_UNAVAILABLE"
                    decision = "Local TF search resolves the gene, but only an incomplete canonical DBD is retained."
                elif high_absent(metrics):
                    status = "ACCESSION_LACKS_DBD"
                    action = "NO_SEND_MARK_DNA_STRUCTURE_UNAVAILABLE"
                    decision = "Local TF search resolves the gene and canonical alignment supports loss of the DBD."
                else:
                    status = "UNCERTAIN_MANUAL_REVIEW"
                    action = "MANUAL_REVIEW_BEFORE_SENDING"
                    decision = "Local TF search resolves a likely gene, but DBD completeness remains below adjudication thresholds."
                    remaining = "Gene suggestion passes search thresholds; DBD interpretation remains uncertain."
            else:
                status = "UNCERTAIN_MANUAL_REVIEW"
                action = "MANUAL_REVIEW_BEFORE_SENDING"
                decision = "No local TF sequence match passes the required gene-suggestion thresholds."
                remaining = "Gene mapping remains unresolved."

        if status == "ACCESSION_CONTAINS_INCOMPLETE_DBD":
            action = "NO_SEND_MARK_DNA_STRUCTURE_UNAVAILABLE"
        if status == "VALID_DBD_FRAGMENT_MODEL":
            action = "NO_SEND_USE_EXISTING_DBD_FRAGMENT"

        send_coords = metrics.get("failed_coords", "") if action.startswith("SEND_") else ""
        send_sequence = ""
        if action.startswith("SEND_"):
            if not (
                float(metrics["coverage"]) >= 80
                and float(metrics["identity"]) >= 80
                and int(metrics["covered"]) >= 30
                and send_coords
            ):
                raise SystemExit(f"Invalid SEND evidence for {accession}")
            start, end = map(int, send_coords.split("-"))
            send_sequence = failed_sequence[start - 1:end]

        actual_pdb_accession = ""
        if model_path:
            match = re.search(r"_(?:TFS_|DIMER_)?([^:]+):", Path(model_path).name)
            if match:
                actual_pdb_accession = match.group(1)
            else:
                match = re.search(r"(?:TFS_|DIMER_)([^:]+):", Path(model_path).name)
                actual_pdb_accession = match.group(1) if match else ""
        valid_fragment = "YES" if status == "VALID_DBD_FRAGMENT_MODEL" else "NO"
        disposition = (
            "COMPLETE_DBD" if status in {
                "VALID_DBD_FRAGMENT_MODEL", "DBD_PRESENT_FRAGMENT_FAILED",
                "DBD_PRESENT_FRAGMENT_NOT_SENT",
            }
            else "INCOMPLETE_DBD" if status == "ACCESSION_CONTAINS_INCOMPLETE_DBD"
            else "NO_DBD" if status in {
                "ACCESSION_LACKS_DBD", "COFACTOR_OR_NON_DNA_BINDING_COMPONENT",
            }
            else "UNCERTAIN"
        )
        display = display_recommendation(status, canonical)
        final = {
            "failed_accession": accession,
            "gene": gene,
            "protein_description": base["protein_description"],
            "failed_sequence_length": len(failed_sequence),
            "failed_sequence_sha256": seq_hash(failed_sequence),
            "DBD_disposition": disposition,
            "expected_DBD_type": dbd_type,
            "expected_DBD_canonical_start": interval[0] if interval else "",
            "expected_DBD_canonical_end": interval[1] if interval else "",
            "DBD_aligned_residues": metrics.get("covered", ""),
            "DBD_coverage_percent": f"{float(metrics.get('coverage', 0)):.2f}",
            "DBD_identity_percent": f"{float(metrics.get('identity', 0)):.2f}",
            "failed_coordinates_corresponding_to_DBD": metrics.get("failed_coords", ""),
            "reviewed_canonical_accession": canonical,
            "model_bearing_accession": model_accession,
            "actual_accession_represented_by_PDB": actual_pdb_accession,
            "canonical_reference_model": model_path,
            "canonical_reference_interface_status": interface_status,
            "valid_existing_fragment_model": valid_fragment,
            "original_fragment_status": (
                "SENT_VALID" if valid_fragment == "YES"
                else "NON_DBD_CONTACTS_NOT_ACCEPTED"
                if base["non_DBD_contacting_fragment_ids"]
                else "NO_VALID_DBD_FRAGMENT"
            ),
            "primary_structural_status": status,
            "action_for_baldo": action,
            "send_coordinates": send_coords,
            "database_display_recommendation": display,
            "decision_reason": decision,
            "remaining_uncertainty": remaining,
            "fusion_parent_analysis": fusion_json,
            "local_sequence_search_results": search_json,
            "draft_status": "FINAL_ADJUDICATION_DRAFT_NOT_PUBLIC",
        }
        final_rows.append(final)
        baldo_rows.append({
            "failed_accession": accession,
            "gene": gene,
            "protein_description": base["protein_description"],
            "primary_structural_status": status,
            "DBD_evidence": (
                f"{dbd_type}; coverage={float(metrics.get('coverage', 0)):.2f}%; "
                f"identity={float(metrics.get('identity', 0)):.2f}%; "
                f"aligned={metrics.get('covered', 0)}"
            ),
            "fragment_status": final["original_fragment_status"],
            "canonical_reference_accession": canonical,
            "canonical_reference_model": model_path,
            "canonical_reference_interface_status": interface_status,
            "action_for_baldo": action,
            "sequence_or_fragment_to_send": send_sequence,
            "send_coordinates": send_coords,
            "database_display_recommendation": display,
            "decision_reason": decision,
            "remaining_uncertainty": remaining,
        })
        display_rows.append({
            "failed_accession": accession,
            "gene": gene,
            "primary_structural_status": status,
            "canonical_reference_accession": canonical,
            "canonical_reference_is_separate": "YES" if canonical else "NOT_APPLICABLE",
            "database_display_recommendation": display,
            "display_rationale": decision,
            "draft_only": "YES",
        })

    by_id = {r["failed_accession"]: r for r in final_rows}
    manual_rows = [
        r for r in final_rows
        if r["primary_structural_status"] == "UNCERTAIN_MANUAL_REVIEW"
    ]
    send_rows = [r for r in baldo_rows if r["action_for_baldo"].startswith("SEND_")]
    after_analysis = db_signature(args.database)
    checks: dict[str, Any] = {}
    add_check(checks, "exactly_102_unique_accessions",
              len(final_rows) == len({r["failed_accession"] for r in final_rows}) == 102,
              {"rows": len(final_rows), "unique": len({r["failed_accession"] for r in final_rows})},
              {"rows": 102, "unique": 102})
    add_check(checks, "one_allowed_status_and_action_per_accession",
              all(r["primary_structural_status"] in STATUSES and r["action_for_baldo"] in ACTIONS for r in final_rows),
              True, True)
    previous_28 = {
        r["failed_accession"] for r in dossier_rows
        if r["suggested_baldo_action"] == "REQUIRES_MANUAL_REVIEW"
    }
    add_check(checks, "all_28_previous_manual_cases_represented",
              previous_28 <= set(by_id) and len(previous_28) == 28,
              sorted(previous_28), "28 represented")
    invalid_send = [
        r["failed_accession"] for r in final_rows
        if r["action_for_baldo"].startswith("SEND_") and not (
            float(r["DBD_coverage_percent"]) >= 80
            and float(r["DBD_identity_percent"]) >= 80
            and int(r["DBD_aligned_residues"]) >= 30
            and r["send_coordinates"]
        )
    ]
    add_check(checks, "no_send_below_80_80_30", not invalid_send, invalid_send, [])
    invalid_non_dbd = [
        r["failed_accession"] for r in final_rows
        if structural[r["failed_accession"]]["non_DBD_contacting_fragment_ids"]
        and r["primary_structural_status"] == "VALID_DBD_FRAGMENT_MODEL"
        and not structural[r["failed_accession"]]["reviewed_DBD_pfam_hits"]
    ]
    add_check(checks, "no_non_DBD_contacting_fragment_accepted", not invalid_non_dbd, invalid_non_dbd, [])
    bad_canonical = [
        r["failed_accession"] for r in final_rows
        if r["reviewed_canonical_accession"]
        and r["reviewed_canonical_accession"] == r["failed_accession"]
        and r["canonical_reference_model"]
    ]
    add_check(checks, "canonical_accessions_remain_separate", not bad_canonical, bad_canonical, [])
    pdb_mismatch = [
        r["failed_accession"] for r in final_rows
        if r["canonical_reference_model"]
        and r["actual_accession_represented_by_PDB"]
        and r["model_bearing_accession"]
        and r["actual_accession_represented_by_PDB"] != r["model_bearing_accession"]
    ]
    add_check(checks, "PDB_actual_accession_matches_model_bearing_accession", not pdb_mismatch, pdb_mismatch, [])
    add_check(checks, "all_partial_controls_incomplete_no_send",
              all(by_id[a]["primary_structural_status"] == "ACCESSION_CONTAINS_INCOMPLETE_DBD"
                  and by_id[a]["action_for_baldo"] == "NO_SEND_MARK_DNA_STRUCTURE_UNAVAILABLE"
                  for a in PARTIAL_CONTROLS),
              {a: (by_id[a]["primary_structural_status"], by_id[a]["action_for_baldo"]) for a in sorted(PARTIAL_CONTROLS)},
              "all incomplete and no-send")
    add_check(checks, "database_unchanged_after_analysis", before == after_analysis, after_analysis, before)
    add_check(checks, "no_public_database_write_performed", True, True, True)
    for accession in FUSION_PARENTS:
        add_check(checks, f"fusion_{accession}_both_parents_analyzed",
                  len(fusion_results[accession]) == 2,
                  fusion_results[accession], "two parent rows")
    for accession in UNRESOLVED:
        best = local_searches[accession][0]
        add_check(checks, f"unresolved_{accession}_local_search_reported",
                  bool(local_searches[accession])
                  and best["aligned_residues"] >= 50,
                  best, "best local TF match reported")

    failed_checks = [name for name, item in checks.items() if not item["passed"]]
    if failed_checks:
        print("FAILED_QC_CHECKS")
        for name in failed_checks:
            print(name, json.dumps(checks[name], sort_keys=True), sep="\t")
        raise SystemExit("QC failed before output write")

    write_tsv(args.output_final, final_rows, FINAL_COLUMNS)
    write_tsv(args.output_manual, manual_rows, FINAL_COLUMNS)
    write_tsv(args.output_baldo, baldo_rows, BALDO_COLUMNS)
    write_tsv(args.output_display, display_rows, DISPLAY_COLUMNS)
    with args.output_fasta.open("x", encoding="utf-8") as handle:
        for row in send_rows:
            sequence = row["sequence_or_fragment_to_send"]
            handle.write(
                f">{row['failed_accession']}|gene={row['gene']}|"
                f"coordinates={row['send_coordinates']}|action={row['action_for_baldo']}\n"
            )
            for offset in range(0, len(sequence), 80):
                handle.write(sequence[offset:offset + 80] + "\n")

    fasta_records_written = fasta_records(args.output_fasta)
    send_ids = {r["failed_accession"] for r in send_rows}
    add_check(checks, "FASTA_entries_exactly_equal_SEND_actions",
              set(fasta_records_written) == send_ids,
              {"fasta": sorted(fasta_records_written), "send": sorted(send_ids)},
              "exact equality")
    after_write = db_signature(args.database)
    add_check(checks, "database_unchanged_after_outputs", before == after_write, after_write, before)
    output_hashes = {
        str(path): file_sha256(path)
        for path in outputs if path != args.output_qc
    }
    qc = {
        "all_checks_passed": all(item["passed"] for item in checks.values()),
        "final_adjudication_is_draft": True,
        "public_database_write_performed": False,
        "counts": {
            "primary_structural_status": dict(sorted(Counter(r["primary_structural_status"] for r in final_rows).items())),
            "action_for_baldo": dict(sorted(Counter(r["action_for_baldo"] for r in final_rows).items())),
            "final_manual_review": len(manual_rows),
            "fasta_sequences": len(fasta_records_written),
        },
        "send_accessions": sorted(send_ids),
        "unresolved_accessions": sorted(
            r["failed_accession"] for r in final_rows
            if r["primary_structural_status"] == "UNCERTAIN_MANUAL_REVIEW"
        ),
        "fusion_results": fusion_results,
        "partial_DBD_results": {
            accession: {
                "status": by_id[accession]["primary_structural_status"],
                "action": by_id[accession]["action_for_baldo"],
                "coverage": by_id[accession]["DBD_coverage_percent"],
                "identity": by_id[accession]["DBD_identity_percent"],
                "aligned": by_id[accession]["DBD_aligned_residues"],
            }
            for accession in sorted(PARTIAL_CONTROLS)
        },
        "local_sequence_searches": local_searches,
        "official_uniprot_sequences_retrieved": official_sequences,
        "empty_fasta_explanation": (
            "No accession satisfies a SEND action after the 80% coverage, "
            "80% identity, >=30-residue, complete-extractable-DBD rules."
            if not send_ids else ""
        ),
        "database_signature_before": before,
        "database_signature_after": after_write,
        "input_hashes": input_hashes,
        "output_hashes_excluding_qc_self_hash": output_hashes,
        "checks": checks,
    }
    with args.output_qc.open("x", encoding="utf-8") as handle:
        json.dump(qc, handle, indent=2, sort_keys=True)
        handle.write("\n")
    if not qc["all_checks_passed"]:
        raise SystemExit("Final QC failed")
    for path in outputs:
        print(f"Wrote {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
