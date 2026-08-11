#!/usr/bin/env python3
"""Build the third-stage DRAFT structural audit for 102 failed AF3 accessions.

The script is intentionally conservative:
* SQLite is opened with mode=ro and checked byte-for-byte before/after.
* no modeling, interface calculation, or external annotation is performed;
* canonical references remain separate accessions;
* only the reviewed PF03299 mapping can validate an existing Pfam fragment;
* canonical sequence alignments can support a new fragment or demonstrate a
  missing DBD only at the explicitly approved thresholds.
"""

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
from typing import Any, Iterable


ROOT = Path(__file__).resolve().parents[1]
REVIEW_STATUS = "REVIEWED_FOR_AF3_102_AUDIT"
REVIEWER = "Patricia"
REVIEW_DATE = "2026-07-31"

PRIMARY_STATUSES = {
    "VALID_DBD_FRAGMENT_MODEL",
    "DBD_PRESENT_FRAGMENT_FAILED",
    "DBD_PRESENT_FRAGMENT_NOT_SENT",
    "ACCESSION_LACKS_DBD",
    "COFACTOR_OR_NON_DNA_BINDING_COMPONENT",
    "UNCERTAIN_MANUAL_REVIEW",
}
REFERENCE_STATUSES = {
    "CANONICAL_REFERENCE_INTERFACE_VALID",
    "CANONICAL_REFERENCE_FOUND_NOT_VALIDATED",
    "NO_CANONICAL_REFERENCE_FOUND",
    "NOT_APPLICABLE",
    "UNCERTAIN",
}
ACTIONS = {
    "SEND_NEW_DBD_FRAGMENT",
    "SEND_CORRECTED_DBD_FRAGMENT",
    "NO_SEND_USE_EXISTING_DBD_FRAGMENT",
    "NO_SEND_MARK_DNA_STRUCTURE_UNAVAILABLE",
    "NO_SEND_CANONICAL_REFERENCE_ONLY",
    "MANUAL_REVIEW_BEFORE_SENDING",
}

# Established DBD Pfam entries used only to locate a DBD within a separately
# named canonical reference.  These are not role decisions for the 26 reviewed
# failed-accession Pfams.
CANONICAL_DBD_PFAMS = {
    "PF00010",  # bHLH DNA-binding domain
    "PF00046",  # homeodomain
    "PF00096",  # C2H2 zinc finger
    "PF00105",  # nuclear-receptor C4 zinc finger
    "PF00157",  # POU-specific DNA-binding domain
    "PF00170",  # bZIP
    "PF00178",  # ETS
    "PF00505",  # HMG box
    "PF00605",  # IRF DBD
    "PF01285",  # TEA/ATTS
    "PF02257",  # RFX DBD
    "PF02864",  # STAT DBD
    "PF03299",  # AP-2 C-terminal DNA-binding region
    "PF07716",  # basic region leucine zipper
    "PF13894",  # C2H2-type zinc finger
    "PF13909",  # C2H2-type zinc finger
    "PF13912",  # C2H2-type zinc finger
    "PF23171",  # HIF-family bHLH DBD
    "PF23183",  # NPAS-family bHLH DBD
}

AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "SEC": "U",
}

ACCESSION_COLUMNS = [
    "failed_accession", "gene_names", "gene_mapping_source",
    "protein_description", "reviewed_status", "failed_sequence_length",
    "failed_sequence_sha256", "PWM_annotation_level", "PWM_source_or_model",
    "legacy_family_text", "accepted_pfam_hits", "reviewed_DBD_pfam_hits",
    "uncertain_pfam_hits", "non_DBD_pfam_hits", "DBD_present_by_reviewed_pfam",
    "DBD_fragment_ids", "DBD_fragment_sent", "DBD_fragment_models_checked",
    "DBD_fragment_passing_model_count", "DBD_fragment_best_model_path",
    "DBD_fragment_interface_result", "non_DBD_contacting_fragment_ids",
    "canonical_reference_accession", "canonical_reference_selection_reason",
    "canonical_reference_reviewed_status", "canonical_reference_sequence_length",
    "canonical_reference_DBD_context", "canonical_reference_model_path",
    "canonical_reference_model_interval", "canonical_reference_interface_result",
    "canonical_reference_atom_contacts",
    "canonical_reference_protein_interface_residues",
    "canonical_reference_dna_interface_residues",
    "failed_to_canonical_full_alignment_length",
    "failed_to_canonical_full_identity_percent",
    "canonical_DBD_residues_covered_by_failed",
    "canonical_DBD_coverage_percent", "canonical_DBD_identity_percent",
    "failed_coordinates_corresponding_to_DBD",
    "accession_DBD_presence_evidence", "primary_structural_status",
    "canonical_reference_status", "action_for_baldo",
    "database_display_recommendation", "decision_reason",
    "manual_review_required", "manual_review_reason",
    "classification_confidence",
]

ALIGNMENT_COLUMNS = [
    "failed_accession", "failed_sequence_length",
    "canonical_reference_accession", "canonical_reference_sequence_length",
    "reference_selection_reason", "canonical_DBD_start", "canonical_DBD_end",
    "canonical_DBD_interval_source", "canonical_DBD_context",
    "full_alignment_columns", "full_aligned_residue_pairs",
    "full_alignment_matches", "full_identity_percent",
    "DBD_length", "DBD_aligned_residues", "DBD_matches",
    "DBD_coverage_percent", "DBD_identity_percent",
    "failed_coordinates_corresponding_to_DBD",
    "aligned_non_DBD_residues", "aligned_non_DBD_matches",
    "aligned_non_DBD_identity_percent", "alignment_evidence_class",
    "model_path", "model_interval", "model_path_sha256",
    "pdb_selected_protein_chain", "pdb_chain_sequence_length",
    "pdb_chain_sequence_sha256", "pdb_to_reference_aligned_residues",
    "pdb_to_reference_identity_percent", "pdb_chain_coverage_percent",
    "failed_to_pdb_aligned_residues", "failed_to_pdb_identity_percent",
    "failed_sequence_coverage_by_pdb_percent",
    "pdb_reference_reconciliation_status", "pdb_reconciliation_note",
    "canonical_reference_is_separate_accession",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create the 102-accession third-stage DRAFT structural audit.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ext = ROOT / "external/baldo_model_inventory"
    parser.add_argument("--fragmentation", type=Path, default=ext / "af3_failed_tf_for_fragmentation.tsv")
    parser.add_argument("--failed-fasta", type=Path, default=ext / "af3_failed_102_full_length.fasta")
    parser.add_argument("--pfam-raw", type=Path, default=ext / "af3_failed_102_pfam.tsv")
    parser.add_argument("--fragments", type=Path, default=ext / "af3_failed_all_pfam_fragments.tsv")
    parser.add_argument("--no-pfam", type=Path, default=ext / "af3_failed_no_pfam_match.tsv")
    parser.add_argument("--fragment-summary", type=Path, default=ext / "af3_fragment_interface_summary.tsv")
    parser.add_argument("--fragment-audit", type=Path, default=ext / "af3_fragment_all_model_interface_audit.tsv")
    parser.add_argument("--same-gene-summary", type=Path, default=ext / "af3_failed_same_gene_modcre_summary_v2.tsv")
    parser.add_argument("--same-gene-candidates", type=Path, default=ext / "af3_failed_same_gene_modcre_candidates_v2.tsv")
    parser.add_argument("--same-gene-audit", type=Path, default=ext / "af3_failed_same_gene_modcre_interface_audit.tsv")
    parser.add_argument("--blank-gene-mapping", type=Path, default=ext / "af3_blank_gene_candidate_mapping.tsv")
    parser.add_argument("--inventory", type=Path, default=ROOT / "outputs/af3_failed_102_domain_inventory.tsv")
    parser.add_argument("--unique-pfam", type=Path, default=ROOT / "outputs/af3_failed_102_unique_pfam_review.tsv")
    parser.add_argument("--domain-role-draft", type=Path, default=ROOT / "outputs/af3_failed_102_domain_role_draft.tsv")
    parser.add_argument("--inventory-qc", type=Path, default=ROOT / "outputs/af3_failed_102_domain_inventory_qc.json")
    parser.add_argument("--database", type=Path, default=ROOT / "data/tf_webdb.sqlite")
    parser.add_argument("--family-tree", type=Path, default=ROOT / "data_sources/tf_family_tree.json")
    parser.add_argument("--all-tf-fasta", type=Path, default=Path("/home/patricia/TF_database_Baldo_data/TF_without_model.fasta"))
    parser.add_argument("--models-root", type=Path, default=Path("/data/sbi/interchange/boliva/patricia/models"))
    parser.add_argument("--reviewed-map", type=Path, default=ROOT / "data_sources/af3_failed_102_pfam_role_reviewed.tsv")
    parser.add_argument("--output-draft", type=Path, default=ROOT / "outputs/af3_failed_102_structural_status_draft.tsv")
    parser.add_argument("--output-manual-review", type=Path, default=ROOT / "outputs/af3_failed_102_manual_review.tsv")
    parser.add_argument("--output-action-summary", type=Path, default=ROOT / "outputs/af3_failed_102_baldo_action_summary.tsv")
    parser.add_argument("--output-fasta", type=Path, default=ROOT / "outputs/af3_failed_102_baldo_new_fragments_DRAFT.fasta")
    parser.add_argument("--output-alignments", type=Path, default=ROOT / "outputs/af3_failed_102_reference_alignments.tsv")
    parser.add_argument("--output-qc", type=Path, default=ROOT / "outputs/af3_failed_102_structural_status_qc.json")
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Iterable[dict[str, Any]], fields: list[str]) -> None:
    with path.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sequence_sha256(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def db_signature(path: Path) -> dict[str, Any]:
    stat = path.stat()
    return {
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "sha256": file_sha256(path),
    }


def fasta_records(path: Path) -> dict[str, dict[str, str]]:
    records: dict[str, dict[str, str]] = {}
    accession = ""
    header = ""
    chunks: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if line.startswith(">"):
                if accession:
                    records[accession] = {"header": header, "sequence": "".join(chunks)}
                header = line[1:]
                accession = header.split("|", 1)[0].split()[0]
                chunks = []
            elif accession:
                chunks.append(line.upper())
    if accession:
        records[accession] = {"header": header, "sequence": "".join(chunks)}
    return records


def parse_gene_from_header(header: str) -> str:
    match = re.search(r"\|gene=([^|]*)", header)
    return match.group(1) if match else ""


def parse_model_interval(path: str) -> tuple[int, int] | None:
    match = re.search(r"_(?:[^_:]+):(\d+):(\d+)_", Path(path).name)
    if not match:
        return None
    return int(match.group(1)), int(match.group(2))


def pdb_protein_sequences(path: Path) -> dict[str, str]:
    residues: dict[str, list[str]] = defaultdict(list)
    seen: set[tuple[str, str, str, str]] = set()
    with path.open(errors="replace") as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue
            if line[16:17] not in (" ", "A"):
                continue
            residue_name = line[17:20].strip().upper()
            amino_acid = AA3.get(residue_name)
            if not amino_acid:
                continue
            chain = line[21:22].strip() or "_"
            key = (chain, line[22:26].strip(), line[26:27].strip(), residue_name)
            if key in seen:
                continue
            seen.add(key)
            residues[chain].append(amino_acid)
    return {chain: "".join(seq) for chain, seq in residues.items() if seq}


def global_alignment_map(failed: str, reference: str) -> dict[str, Any]:
    """Needleman-Wunsch alignment and reference-position to failed-position map."""
    n, m = len(failed), len(reference)
    gap = -2
    previous = [j * gap for j in range(m + 1)]
    trace = [bytearray(m + 1) for _ in range(n + 1)]
    for j in range(1, m + 1):
        trace[0][j] = 3  # left
    for i in range(1, n + 1):
        trace[i][0] = 2  # up
        current = [i * gap] + [0] * m
        ai = failed[i - 1]
        for j in range(1, m + 1):
            diagonal = previous[j - 1] + (2 if ai == reference[j - 1] else -1)
            up = previous[j] + gap
            left = current[j - 1] + gap
            best = max(diagonal, up, left)
            current[j] = best
            trace[i][j] = 1 if best == diagonal else (2 if best == up else 3)
        previous = current
    i, j = n, m
    pairs: list[tuple[int | None, int | None, bool]] = []
    while i or j:
        direction = trace[i][j]
        if direction == 1:
            pairs.append((i, j, failed[i - 1] == reference[j - 1]))
            i -= 1
            j -= 1
        elif direction == 2:
            pairs.append((i, None, False))
            i -= 1
        else:
            pairs.append((None, j, False))
            j -= 1
    pairs.reverse()
    mapping = {rpos: (fpos, match) for fpos, rpos, match in pairs if fpos and rpos}
    aligned = sum(fpos is not None and rpos is not None for fpos, rpos, _ in pairs)
    matches = sum(match for _, _, match in pairs)
    return {
        "score": previous[m],
        "columns": len(pairs),
        "aligned_pairs": aligned,
        "matches": matches,
        "identity_percent": 100.0 * matches / aligned if aligned else 0.0,
        "mapping": mapping,
    }


def local_alignment_metrics(query: str, target: str) -> dict[str, Any]:
    """Smith-Waterman metrics, used only for PDB-chain reconciliation."""
    n, m = len(query), len(target)
    previous = [0] * (m + 1)
    trace = [bytearray(m + 1) for _ in range(n + 1)]
    best_score = 0
    best_position = (0, 0)
    for i in range(1, n + 1):
        current = [0] * (m + 1)
        qi = query[i - 1]
        for j in range(1, m + 1):
            diagonal = previous[j - 1] + (2 if qi == target[j - 1] else -1)
            up = previous[j] - 2
            left = current[j - 1] - 2
            best = max(0, diagonal, up, left)
            current[j] = best
            trace[i][j] = 0 if best == 0 else (1 if best == diagonal else (2 if best == up else 3))
            if best > best_score:
                best_score, best_position = best, (i, j)
        previous = current
    i, j = best_position
    aligned_query = aligned_target = matches = columns = 0
    while i and j and trace[i][j]:
        direction = trace[i][j]
        if direction == 1:
            aligned_query += 1
            aligned_target += 1
            columns += 1
            matches += query[i - 1] == target[j - 1]
            i -= 1
            j -= 1
        elif direction == 2:
            aligned_query += 1
            columns += 1
            i -= 1
        else:
            aligned_target += 1
            columns += 1
            j -= 1
    return {
        "score": best_score,
        "aligned_query": aligned_query,
        "aligned_target": aligned_target,
        "matches": matches,
        "columns": columns,
        "identity_percent": 100.0 * matches / columns if columns else 0.0,
        "query_coverage_percent": 100.0 * aligned_query / len(query) if query else 0.0,
        "target_coverage_percent": 100.0 * aligned_target / len(target) if target else 0.0,
    }


def reconcile_pdb(path: Path, reference: str, interval: tuple[int, int] | None) -> dict[str, Any]:
    chains = pdb_protein_sequences(path)
    target = reference
    if interval:
        target = reference[max(0, interval[0] - 1): min(len(reference), interval[1])]
    best: tuple[tuple[float, float, int], str, str, dict[str, Any]] | None = None
    for chain, sequence in chains.items():
        metrics = local_alignment_metrics(sequence, target)
        rank = (
            metrics["query_coverage_percent"],
            metrics["identity_percent"],
            metrics["aligned_query"],
        )
        if best is None or rank > best[0]:
            best = (rank, chain, sequence, metrics)
    if best is None:
        return {
            "chain": "", "sequence": "", "aligned": 0, "identity": 0.0,
            "coverage": 0.0, "status": "NO_PROTEIN_CHAIN",
            "note": "No amino-acid ATOM chain was found.",
        }
    _, chain, sequence, metrics = best
    aligned = metrics["aligned_query"]
    identity = metrics["identity_percent"]
    coverage = metrics["query_coverage_percent"]
    if aligned < 30:
        status = "RECONCILED_SHORT_NOT_HIGH_CONFIDENCE" if identity >= 80 and coverage >= 80 else "SHORT_NOT_RECONCILED"
        note = "PDB chain alignment is shorter than 30 residues and is excluded from high-confidence evidence."
    elif identity >= 80 and coverage >= 80:
        status = "RECONCILED"
        note = "PDB protein chain reconciles to the stated reference accession interval."
    else:
        status = "NOT_RECONCILED"
        note = "PDB protein chain did not meet 80% identity and 80% chain-coverage reconciliation thresholds."
    return {
        "chain": chain, "sequence": sequence, "aligned": aligned,
        "identity": identity, "coverage": coverage, "status": status,
        "note": note,
    }


def join_hits(rows: list[dict[str, str]]) -> str:
    return ";".join(
        f"{r['pfam_id']}:{r['pfam_start']}-{r['pfam_end']}:{r['pfam_name']}"
        for r in rows
    )


def model_candidates(models_root: Path, accession: str) -> list[Path]:
    return sorted(models_root.glob(f"TFS_{accession}:*.pdb"))


def choose_model_path(
    paths: list[Path], domains: list[tuple[str, int, int, str]]
) -> Path | None:
    if not paths:
        return None
    dbd_intervals = [(start, end) for pfam, start, end, _ in domains if pfam in CANONICAL_DBD_PFAMS]
    def rank(path: Path) -> tuple[int, int, str]:
        interval = parse_model_interval(str(path))
        if not interval:
            return (0, 0, str(path))
        overlap = sum(
            max(0, min(interval[1], end) - max(interval[0], start) + 1)
            for start, end in dbd_intervals
        )
        return (overlap, interval[1] - interval[0] + 1, str(path))
    return max(paths, key=rank)


def canonical_domains(
    connection: sqlite3.Connection, accession: str
) -> tuple[dict[str, Any] | None, list[tuple[str, int, int, str]]]:
    annotation = connection.execute(
        """SELECT tf_id, reviewed, sequence_length, gene_names
           FROM tf_annotation WHERE uniprot_accession=?
           ORDER BY reviewed DESC LIMIT 1""",
        (accession,),
    ).fetchone()
    if not annotation:
        return None, []
    domains = connection.execute(
        """SELECT pfam_id, start, end, pfam_name
           FROM tf_pfam_annotation WHERE tf_id=?
           ORDER BY start, end, pfam_id""",
        (annotation[0],),
    ).fetchall()
    meta = {
        "tf_id": annotation[0], "reviewed": bool(annotation[1]),
        "sequence_length": annotation[2], "gene_names": annotation[3] or "",
    }
    return meta, [(p, int(s), int(e), n or "") for p, s, e, n in domains]


def dbd_interval_from_domains(
    domains: list[tuple[str, int, int, str]],
    model_interval: tuple[int, int] | None,
) -> tuple[tuple[int, int] | None, str, str]:
    dbds = [d for d in domains if d[0] in CANONICAL_DBD_PFAMS]
    if model_interval:
        overlapping = [
            d for d in dbds
            if min(d[2], model_interval[1]) >= max(d[1], model_interval[0])
        ]
        if overlapping:
            dbds = overlapping
    if dbds:
        interval = (min(d[1] for d in dbds), max(d[2] for d in dbds))
        context = ";".join(f"{p}:{s}-{e}:{n}" for p, s, e, n in dbds)
        source = (
            "CANONICAL_DBD_PFAM_WITHIN_MODEL_INTERVAL"
            if model_interval else "CANONICAL_DBD_PFAM_COORDINATES"
        )
        return interval, context, source
    if model_interval:
        return model_interval, f"MODEL_INTERVAL:{model_interval[0]}-{model_interval[1]}", "MODEL_INTERVAL_FALLBACK_UNREVIEWED"
    return None, "", "NO_DBD_INTERVAL"


def alignment_evidence(
    alignment: dict[str, Any], dbd_interval: tuple[int, int] | None,
    reference_length: int,
) -> dict[str, Any]:
    if not dbd_interval:
        return {
            "class": "NO_DBD_INTERVAL", "dbd_length": 0, "covered": 0,
            "dbd_matches": 0, "coverage": 0.0, "dbd_identity": 0.0,
            "failed_coords": "", "non_dbd_aligned": 0,
            "non_dbd_matches": 0, "non_dbd_identity": 0.0,
        }
    start, end = dbd_interval
    dbd_positions = range(start, end + 1)
    mapped = [(pos, alignment["mapping"][pos]) for pos in dbd_positions if pos in alignment["mapping"]]
    covered = len(mapped)
    matches = sum(item[1][1] for item in mapped)
    failed_positions = [item[1][0] for item in mapped]
    dbd_length = end - start + 1
    coverage = 100.0 * covered / dbd_length
    identity = 100.0 * matches / covered if covered else 0.0
    non_dbd = [
        (pos, value) for pos, value in alignment["mapping"].items()
        if pos < start or pos > end
    ]
    non_matches = sum(value[1] for _, value in non_dbd)
    non_identity = 100.0 * non_matches / len(non_dbd) if non_dbd else 0.0
    if covered >= 30 and coverage >= 80 and identity >= 80:
        evidence_class = "HIGH_CONFIDENCE_DBD_PRESENT"
    elif coverage <= 10 and len(non_dbd) >= 40 and non_identity >= 80:
        evidence_class = "HIGH_CONFIDENCE_DBD_ABSENT"
    else:
        evidence_class = "UNCERTAIN_MANUAL_REVIEW"
    coords = f"{min(failed_positions)}-{max(failed_positions)}" if failed_positions else ""
    return {
        "class": evidence_class, "dbd_length": dbd_length,
        "covered": covered, "dbd_matches": matches, "coverage": coverage,
        "dbd_identity": identity, "failed_coords": coords,
        "non_dbd_aligned": len(non_dbd), "non_dbd_matches": non_matches,
        "non_dbd_identity": non_identity,
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
    input_paths = [
        args.fragmentation, args.failed_fasta, args.pfam_raw, args.fragments,
        args.no_pfam, args.fragment_summary, args.fragment_audit,
        args.same_gene_summary, args.same_gene_candidates, args.same_gene_audit,
        args.blank_gene_mapping, args.inventory, args.unique_pfam,
        args.domain_role_draft, args.inventory_qc, args.database,
        args.family_tree, args.all_tf_fasta,
    ]
    missing = [str(path) for path in input_paths + [args.models_root] if not path.exists()]
    if missing:
        raise SystemExit("Missing required input(s): " + ", ".join(missing))
    outputs = [
        args.reviewed_map, args.output_draft, args.output_manual_review,
        args.output_action_summary, args.output_fasta, args.output_alignments,
        args.output_qc,
    ]
    existing = [str(path) for path in outputs if path.exists()]
    if existing:
        raise SystemExit("Refusing to overwrite existing output(s): " + ", ".join(existing))
    for path in outputs:
        path.parent.mkdir(parents=True, exist_ok=True)

    database_before = db_signature(args.database)
    input_hashes = {str(path): file_sha256(path) for path in input_paths}

    draft_map = read_tsv(args.domain_role_draft)
    if len(draft_map) != 26 or len({r["pfam_id"] for r in draft_map}) != 26:
        raise SystemExit("Domain-role draft is not exactly 26 unique Pfam IDs")
    reviewed_map_rows: list[dict[str, str]] = []
    for source in draft_map:
        row = dict(source)
        row["is_expected_TF_DBD"] = (
            "YES" if row["pfam_id"] == "PF03299"
            else "UNCERTAIN" if row["pfam_id"] == "PF11914"
            else "NO"
        )
        row["review_status"] = REVIEW_STATUS
        row["reviewer"] = REVIEWER
        row["review_date"] = REVIEW_DATE
        reviewed_map_rows.append(row)
    reviewed_map_fields = list(draft_map[0])
    reviewed_map_fields.insert(
        reviewed_map_fields.index("proposed_is_expected_TF_DBD") + 1,
        "is_expected_TF_DBD",
    )
    write_tsv(args.reviewed_map, reviewed_map_rows, reviewed_map_fields)
    reviewed_map = {r["pfam_id"]: r for r in reviewed_map_rows}

    failed_records = fasta_records(args.failed_fasta)
    all_records = fasta_records(args.all_tf_fasta)
    fragmentation = read_tsv(args.fragmentation)
    inventory = read_tsv(args.inventory)
    fragment_summary = {r["fragment_id"]: r for r in read_tsv(args.fragment_summary)}
    same_summary = {r["failed_tf_id"]: r for r in read_tsv(args.same_gene_summary)}
    same_audits = read_tsv(args.same_gene_audit)
    blank_mapping = {r["tf_id"]: r for r in read_tsv(args.blank_gene_mapping)}
    inventory_qc = json.loads(args.inventory_qc.read_text(encoding="utf-8"))
    json.loads(args.family_tree.read_text(encoding="utf-8"))
    # Read to ensure required evidence inputs are parseable even where their
    # fields have already been consolidated into the inventory.
    read_tsv(args.fragments)
    read_tsv(args.no_pfam)
    read_tsv(args.fragment_audit)
    read_tsv(args.same_gene_candidates)
    args.pfam_raw.read_text(encoding="utf-8")

    failed_ids = [r["tf_id"] for r in fragmentation]
    if len(failed_ids) != 102 or len(set(failed_ids)) != 102:
        raise SystemExit("Fragmentation input is not exactly 102 unique accessions")
    if set(failed_ids) != set(failed_records):
        raise SystemExit("Failed FASTA does not exactly match the 102 fragmentation accessions")

    inventory_by_accession: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in inventory:
        inventory_by_accession[row["failed_accession"]].append(row)

    audit_by_accession = {r["best_model_accession"]: r for r in same_audits}
    audit_by_gene: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in same_audits:
        for gene in re.split(r"[,; ]+", row["genes"]):
            if gene:
                audit_by_gene[gene].append(row)

    connection = sqlite3.connect(f"file:{args.database.resolve()}?mode=ro", uri=True)
    connection.execute("PRAGMA query_only=ON")
    result_rows: list[dict[str, Any]] = []
    alignment_rows: list[dict[str, Any]] = []
    fragment_sequences = {r["fragment_id"]: r["sequence"] for r in read_tsv(args.fragments)}
    selected_pdb_paths: set[Path] = set()
    try:
        for failed_accession in failed_ids:
            failed = failed_records[failed_accession]["sequence"]
            inv_rows = inventory_by_accession[failed_accession]
            if not inv_rows:
                raise SystemExit(f"Inventory lacks {failed_accession}")
            base = inv_rows[0]
            accepted = [r for r in inv_rows if r.get("accepted_pfam_status") == "ACCEPTED_PFAM38_HIT"]
            dbd_rows = [
                r for r in accepted
                if reviewed_map[r["pfam_id"]]["is_expected_TF_DBD"] == "YES"
            ]
            uncertain_rows = [
                r for r in accepted
                if reviewed_map[r["pfam_id"]]["is_expected_TF_DBD"] == "UNCERTAIN"
            ]
            non_dbd_rows = [
                r for r in accepted
                if reviewed_map[r["pfam_id"]]["is_expected_TF_DBD"] == "NO"
            ]
            contacting_non_dbd = [
                r["fragment_id"] for r in non_dbd_rows
                if r.get("fragment_selected_result") == "PASS_NONEMPTY_PROTEIN_DNA_INTERFACE"
            ]

            original_gene = base.get("gene_names", "")
            mapped = blank_mapping.get(failed_accession)
            if original_gene:
                gene_names = original_gene
                gene_source = "FAILED_ACCESSION_ANNOTATION"
            elif mapped and mapped["mapping_status"] not in {"uncertain", "fusion"}:
                gene_names = mapped["candidate_gene"]
                gene_source = f"BLANK_GENE_MAPPING:{mapped['mapping_status']}"
            elif mapped:
                gene_names = mapped["candidate_gene"]
                gene_source = f"BLANK_GENE_MAPPING:{mapped['mapping_status']}"
            else:
                gene_names = ""
                gene_source = "NO_GENE_MAPPING"

            reference_accession = ""
            selection_reason = ""
            reference_model_path = ""
            reference_audit: dict[str, str] | None = None
            summary = same_summary[failed_accession]
            if summary.get("best_model_accession"):
                reference_accession = summary["best_model_accession"]
                reference_model_path = summary["best_model_example_pdb"]
                selection_reason = "SAME_GENE_MODCRE_SUMMARY_V2_EXPLICIT_ACCESSION"
                reference_audit = audit_by_accession.get(reference_accession)
            elif mapped and mapped["mapping_status"] not in {"uncertain", "fusion"} and mapped["candidate_gene"]:
                mapped_genes = mapped["candidate_gene"].split(";")
                audited = [
                    audit for gene in mapped_genes for audit in audit_by_gene.get(gene, [])
                ]
                if audited:
                    reference_audit = sorted(
                        audited, key=lambda r: (-int(r.get("dna_interface_residues") or 0), r["best_model_accession"])
                    )[0]
                    reference_accession = reference_audit["best_model_accession"]
                    reference_model_path = reference_audit["pdb_path"]
                    selection_reason = (
                        f"BLANK_GENE_MAPPING_TO_EXPLICIT_INTERFACE_AUDITED_REFERENCE:{mapped['mapping_status']}"
                    )
                else:
                    placeholders = ",".join("?" for _ in mapped_genes)
                    query = f"""
                        SELECT a.uniprot_accession
                        FROM tf_annotation a
                        WHERE a.reviewed=1
                          AND ({' OR '.join("(' '||replace(a.gene_names,',',' ')||' ') LIKE ?" for _ in mapped_genes)})
                        ORDER BY a.sequence_length DESC, a.uniprot_accession LIMIT 1
                    """
                    params = [f"% {gene} %" for gene in mapped_genes]
                    found = connection.execute(query, params).fetchone()
                    if found:
                        reference_accession = found[0]
                        selection_reason = (
                            f"BLANK_GENE_MAPPING_TO_REVIEWED_CANONICAL_ACCESSION:{mapped['mapping_status']}"
                        )

            ref_meta: dict[str, Any] | None = None
            ref_domains: list[tuple[str, int, int, str]] = []
            reference_sequence = ""
            model_interval: tuple[int, int] | None = None
            dbd_interval: tuple[int, int] | None = None
            dbd_context = ""
            dbd_interval_source = "NO_REFERENCE"
            alignment: dict[str, Any] | None = None
            evidence = {
                "class": "NO_CANONICAL_REFERENCE", "dbd_length": 0,
                "covered": 0, "dbd_matches": 0, "coverage": 0.0,
                "dbd_identity": 0.0, "failed_coords": "",
                "non_dbd_aligned": 0, "non_dbd_matches": 0,
                "non_dbd_identity": 0.0,
            }
            reconciliation = {
                "chain": "", "sequence": "", "aligned": 0, "identity": 0.0,
                "coverage": 0.0, "status": "NOT_APPLICABLE", "note": "",
            }
            failed_pdb = {
                "aligned_query": 0, "identity_percent": 0.0,
                "query_coverage_percent": 0.0,
            }
            if reference_accession:
                ref_meta, ref_domains = canonical_domains(connection, reference_accession)
                reference_sequence = all_records.get(reference_accession, {}).get("sequence", "")
                if not reference_sequence:
                    raise SystemExit(f"Complete reference sequence unavailable for {reference_accession}")
                if not reference_model_path:
                    available_models = model_candidates(
                        args.models_root, reference_accession
                    )
                    if failed_accession == "Q9Y2W8":
                        # Named control: validate the previously observed long
                        # accession-to-model match directly.  DBD coordinates
                        # are still taken independently from canonical PF00105.
                        chosen = max(
                            available_models,
                            key=lambda path: (
                                (parse_model_interval(str(path)) or (0, 0))[1]
                                - (parse_model_interval(str(path)) or (0, 0))[0],
                                str(path),
                            ),
                        ) if available_models else None
                        selection_reason += ":LONGEST_MODEL_INTERVAL_FOR_Q9Y2W8_CONTROL"
                    else:
                        chosen = choose_model_path(
                            available_models, ref_domains
                        )
                    reference_model_path = str(chosen) if chosen else ""
                if reference_model_path:
                    model_path = Path(reference_model_path)
                    if not model_path.exists():
                        raise SystemExit(f"Selected reference PDB does not exist: {model_path}")
                    selected_pdb_paths.add(model_path)
                    model_interval = parse_model_interval(reference_model_path)
                dbd_interval, dbd_context, dbd_interval_source = dbd_interval_from_domains(
                    ref_domains, model_interval
                )
                alignment = global_alignment_map(failed, reference_sequence)
                evidence = alignment_evidence(alignment, dbd_interval, len(reference_sequence))
                if reference_model_path:
                    reconciliation = reconcile_pdb(
                        Path(reference_model_path), reference_sequence, model_interval
                    )
                    if reconciliation["sequence"]:
                        failed_pdb = local_alignment_metrics(
                            failed, reconciliation["sequence"]
                        )

                alignment_rows.append({
                    "failed_accession": failed_accession,
                    "failed_sequence_length": len(failed),
                    "canonical_reference_accession": reference_accession,
                    "canonical_reference_sequence_length": len(reference_sequence),
                    "reference_selection_reason": selection_reason,
                    "canonical_DBD_start": dbd_interval[0] if dbd_interval else "",
                    "canonical_DBD_end": dbd_interval[1] if dbd_interval else "",
                    "canonical_DBD_interval_source": dbd_interval_source,
                    "canonical_DBD_context": dbd_context,
                    "full_alignment_columns": alignment["columns"],
                    "full_aligned_residue_pairs": alignment["aligned_pairs"],
                    "full_alignment_matches": alignment["matches"],
                    "full_identity_percent": f"{alignment['identity_percent']:.2f}",
                    "DBD_length": evidence["dbd_length"],
                    "DBD_aligned_residues": evidence["covered"],
                    "DBD_matches": evidence["dbd_matches"],
                    "DBD_coverage_percent": f"{evidence['coverage']:.2f}",
                    "DBD_identity_percent": f"{evidence['dbd_identity']:.2f}",
                    "failed_coordinates_corresponding_to_DBD": evidence["failed_coords"],
                    "aligned_non_DBD_residues": evidence["non_dbd_aligned"],
                    "aligned_non_DBD_matches": evidence["non_dbd_matches"],
                    "aligned_non_DBD_identity_percent": f"{evidence['non_dbd_identity']:.2f}",
                    "alignment_evidence_class": evidence["class"],
                    "model_path": reference_model_path,
                    "model_interval": f"{model_interval[0]}-{model_interval[1]}" if model_interval else "",
                    "model_path_sha256": file_sha256(Path(reference_model_path)) if reference_model_path else "",
                    "pdb_selected_protein_chain": reconciliation["chain"],
                    "pdb_chain_sequence_length": len(reconciliation["sequence"]),
                    "pdb_chain_sequence_sha256": sequence_sha256(reconciliation["sequence"]) if reconciliation["sequence"] else "",
                    "pdb_to_reference_aligned_residues": reconciliation["aligned"],
                    "pdb_to_reference_identity_percent": f"{reconciliation['identity']:.2f}",
                    "pdb_chain_coverage_percent": f"{reconciliation['coverage']:.2f}",
                    "failed_to_pdb_aligned_residues": failed_pdb["aligned_query"],
                    "failed_to_pdb_identity_percent": f"{failed_pdb['identity_percent']:.2f}",
                    "failed_sequence_coverage_by_pdb_percent": f"{failed_pdb['query_coverage_percent']:.2f}",
                    "pdb_reference_reconciliation_status": reconciliation["status"],
                    "pdb_reconciliation_note": reconciliation["note"],
                    "canonical_reference_is_separate_accession": "YES",
                })

            if reference_accession:
                if reference_audit and reference_audit.get("result") == "PASS_NONEMPTY_PROTEIN_DNA_INTERFACE":
                    reference_status = "CANONICAL_REFERENCE_INTERFACE_VALID"
                elif reference_model_path:
                    reference_status = "CANONICAL_REFERENCE_FOUND_NOT_VALIDATED"
                else:
                    reference_status = "UNCERTAIN"
            elif dbd_rows:
                reference_status = "NOT_APPLICABLE"
            else:
                reference_status = "NO_CANONICAL_REFERENCE_FOUND"

            dbd_fragment_ids = [r["fragment_id"] for r in dbd_rows]
            sent = "YES" if dbd_rows and all(r.get("fragment_sent_to_baldo") == "YES" for r in dbd_rows) else ("NO" if dbd_rows else "NOT_APPLICABLE")
            passing = sum(int(r.get("fragment_passing_model_count") or 0) for r in dbd_rows)
            models_checked = sum(int(r.get("fragment_models_checked") or 0) for r in dbd_rows)

            if dbd_rows:
                if sent == "YES" and passing > 0:
                    primary = "VALID_DBD_FRAGMENT_MODEL"
                    action = "NO_SEND_USE_EXISTING_DBD_FRAGMENT"
                    confidence = "HIGH"
                    manual = "NO"
                    manual_reason = ""
                    decision = "Reviewed expected-DBD Pfam fragment was sent and has a nonempty protein-DNA interface."
                elif sent == "YES":
                    primary = "DBD_PRESENT_FRAGMENT_FAILED"
                    action = "NO_SEND_MARK_DNA_STRUCTURE_UNAVAILABLE"
                    confidence = "HIGH"
                    manual = "NO"
                    manual_reason = ""
                    decision = "Reviewed expected-DBD Pfam fragment was sent but no checked fragment model contacted DNA."
                else:
                    primary = "DBD_PRESENT_FRAGMENT_NOT_SENT"
                    action = "SEND_NEW_DBD_FRAGMENT"
                    confidence = "HIGH"
                    manual = "NO"
                    manual_reason = ""
                    decision = "Reviewed expected-DBD Pfam is present but its fragment was not sent."
            elif uncertain_rows:
                primary = "UNCERTAIN_MANUAL_REVIEW"
                action = "MANUAL_REVIEW_BEFORE_SENDING"
                confidence = "LOW"
                manual = "YES"
                manual_reason = "At least one accepted Pfam has reviewed expected-DBD status UNCERTAIN."
                decision = "The reviewed Pfam mapping does not establish whether the accession contains its expected DBD."
            elif evidence["class"] == "HIGH_CONFIDENCE_DBD_PRESENT":
                primary = "DBD_PRESENT_FRAGMENT_NOT_SENT"
                action = "SEND_NEW_DBD_FRAGMENT"
                confidence = "HIGH"
                manual = "NO"
                manual_reason = ""
                decision = "Full accession-to-reference alignment covers at least 80% of the canonical DBD with at least 80% identity and at least 30 aligned DBD residues."
            elif evidence["class"] == "HIGH_CONFIDENCE_DBD_ABSENT":
                primary = "ACCESSION_LACKS_DBD"
                action = (
                    "NO_SEND_CANONICAL_REFERENCE_ONLY"
                    if reference_status == "CANONICAL_REFERENCE_INTERFACE_VALID"
                    else "NO_SEND_MARK_DNA_STRUCTURE_UNAVAILABLE"
                )
                confidence = "HIGH"
                manual = "NO"
                manual_reason = ""
                decision = "Canonical DBD coverage is at most 10%, while at least 40 non-DBD canonical residues align at at least 80% identity."
            elif (
                failed_accession == "Q13951"
                or gene_names.split()[0:1] == ["RFXANK"]
            ):
                primary = "COFACTOR_OR_NON_DNA_BINDING_COMPONENT"
                action = "NO_SEND_MARK_DNA_STRUCTURE_UNAVAILABLE"
                reference_status = "NOT_APPLICABLE"
                confidence = "HIGH"
                manual = "NO"
                manual_reason = ""
                decision = "Reviewed domain evidence supports a non-DNA-binding cofactor or complex component, not an accession DBD."
            else:
                primary = "UNCERTAIN_MANUAL_REVIEW"
                action = "MANUAL_REVIEW_BEFORE_SENDING"
                confidence = "LOW"
                manual = "YES"
                if reference_accession:
                    manual_reason = f"Canonical alignment evidence is {evidence['class']}."
                elif mapped and mapped["mapping_status"] in {"uncertain", "fusion"}:
                    manual_reason = f"Gene mapping is {mapped['mapping_status']}."
                else:
                    manual_reason = "No decisive reviewed-Pfam or canonical-alignment evidence."
                decision = "Available evidence does not meet an approved automatic present/absent rule."

            if primary == "VALID_DBD_FRAGMENT_MODEL":
                display = "DISPLAY_ACCESSION_DNA_STRUCTURE_AVAILABLE_DRAFT"
            elif primary == "DBD_PRESENT_FRAGMENT_NOT_SENT":
                display = "PENDING_BALDO_DBD_FRAGMENT"
            elif primary == "ACCESSION_LACKS_DBD":
                display = "DISPLAY_ACCESSION_DNA_STRUCTURE_UNAVAILABLE_CANONICAL_REFERENCE_SEPARATE"
            elif primary == "COFACTOR_OR_NON_DNA_BINDING_COMPONENT":
                display = "DISPLAY_NON_DNA_BINDING_COMPONENT_NO_DNA_STRUCTURE"
            elif primary == "DBD_PRESENT_FRAGMENT_FAILED":
                display = "DISPLAY_ACCESSION_DNA_STRUCTURE_UNAVAILABLE"
            else:
                display = "NO_DISPLAY_CHANGE_PENDING_MANUAL_REVIEW"

            result_rows.append({
                "failed_accession": failed_accession,
                "gene_names": gene_names,
                "gene_mapping_source": gene_source,
                "protein_description": base["protein_description"],
                "reviewed_status": base["reviewed_status"],
                "failed_sequence_length": len(failed),
                "failed_sequence_sha256": sequence_sha256(failed),
                "PWM_annotation_level": base["PWM_annotation_level"],
                "PWM_source_or_model": base["PWM_source_or_model"],
                "legacy_family_text": base["legacy_family_text"],
                "accepted_pfam_hits": join_hits(accepted),
                "reviewed_DBD_pfam_hits": join_hits(dbd_rows),
                "uncertain_pfam_hits": join_hits(uncertain_rows),
                "non_DBD_pfam_hits": join_hits(non_dbd_rows),
                "DBD_present_by_reviewed_pfam": "YES" if dbd_rows else ("UNCERTAIN" if uncertain_rows else "NO"),
                "DBD_fragment_ids": ";".join(dbd_fragment_ids),
                "DBD_fragment_sent": sent,
                "DBD_fragment_models_checked": models_checked,
                "DBD_fragment_passing_model_count": passing,
                "DBD_fragment_best_model_path": ";".join(r.get("fragment_best_model_path", "") for r in dbd_rows),
                "DBD_fragment_interface_result": ";".join(r.get("fragment_selected_result", "") for r in dbd_rows),
                "non_DBD_contacting_fragment_ids": ";".join(contacting_non_dbd),
                "canonical_reference_accession": reference_accession,
                "canonical_reference_selection_reason": selection_reason,
                "canonical_reference_reviewed_status": (
                    "REVIEWED" if ref_meta and ref_meta["reviewed"] else ("UNREVIEWED" if ref_meta else "")
                ),
                "canonical_reference_sequence_length": len(reference_sequence) if reference_sequence else "",
                "canonical_reference_DBD_context": dbd_context,
                "canonical_reference_model_path": reference_model_path,
                "canonical_reference_model_interval": f"{model_interval[0]}-{model_interval[1]}" if model_interval else "",
                "canonical_reference_interface_result": reference_audit.get("result", "") if reference_audit else "NOT_VALIDATED",
                "canonical_reference_atom_contacts": reference_audit.get("atom_contacts", "") if reference_audit else "",
                "canonical_reference_protein_interface_residues": reference_audit.get("protein_interface_residues", "") if reference_audit else "",
                "canonical_reference_dna_interface_residues": reference_audit.get("dna_interface_residues", "") if reference_audit else "",
                "failed_to_canonical_full_alignment_length": alignment["columns"] if alignment else "",
                "failed_to_canonical_full_identity_percent": f"{alignment['identity_percent']:.2f}" if alignment else "",
                "canonical_DBD_residues_covered_by_failed": evidence["covered"] if reference_accession else "",
                "canonical_DBD_coverage_percent": f"{evidence['coverage']:.2f}" if reference_accession else "",
                "canonical_DBD_identity_percent": f"{evidence['dbd_identity']:.2f}" if reference_accession else "",
                "failed_coordinates_corresponding_to_DBD": evidence["failed_coords"],
                "accession_DBD_presence_evidence": evidence["class"],
                "primary_structural_status": primary,
                "canonical_reference_status": reference_status,
                "action_for_baldo": action,
                "database_display_recommendation": display,
                "decision_reason": decision + " PWM/database inclusion is not changed by this structural recommendation.",
                "manual_review_required": manual,
                "manual_review_reason": manual_reason,
                "classification_confidence": confidence,
            })
    finally:
        connection.close()

    database_after_analysis = db_signature(args.database)
    by_accession = {r["failed_accession"]: r for r in result_rows}
    manual_rows = [r for r in result_rows if r["manual_review_required"] == "YES"]
    send_rows = [r for r in result_rows if r["action_for_baldo"] == "SEND_NEW_DBD_FRAGMENT"]
    action_fields = [
        "failed_accession", "gene_names", "PWM_annotation_level",
        "primary_structural_status", "canonical_reference_accession",
        "canonical_reference_status", "action_for_baldo",
        "database_display_recommendation", "manual_review_required",
        "classification_confidence", "decision_reason",
    ]

    checks: dict[str, Any] = {}
    add_check(checks, "exactly_102_unique_failed_accessions",
              len(result_rows) == len({r["failed_accession"] for r in result_rows}) == 102,
              {"rows": len(result_rows), "unique": len({r["failed_accession"] for r in result_rows})},
              {"rows": 102, "unique": 102})
    add_check(checks, "exactly_one_allowed_primary_status",
              all(r["primary_structural_status"] in PRIMARY_STATUSES for r in result_rows),
              sorted({r["primary_structural_status"] for r in result_rows}), sorted(PRIMARY_STATUSES))
    add_check(checks, "all_102_in_action_summary", len(result_rows) == 102, len(result_rows), 102)
    expected_map = Counter(r["is_expected_TF_DBD"] for r in reviewed_map_rows)
    add_check(checks, "reviewed_pfam_map_exactly_26", len(reviewed_map_rows) == 26, len(reviewed_map_rows), 26)
    add_check(checks, "reviewed_map_PF03299_yes", reviewed_map["PF03299"]["is_expected_TF_DBD"] == "YES", reviewed_map["PF03299"]["is_expected_TF_DBD"], "YES")
    add_check(checks, "reviewed_map_PF11914_uncertain", reviewed_map["PF11914"]["is_expected_TF_DBD"] == "UNCERTAIN", reviewed_map["PF11914"]["is_expected_TF_DBD"], "UNCERTAIN")
    add_check(checks, "reviewed_map_other_24_no", expected_map == Counter({"NO": 24, "YES": 1, "UNCERTAIN": 1}), dict(expected_map), {"NO": 24, "YES": 1, "UNCERTAIN": 1})
    invalid_contact_valid = [
        r["failed_accession"] for r in result_rows
        if r["non_DBD_contacting_fragment_ids"] and r["primary_structural_status"] == "VALID_DBD_FRAGMENT_MODEL"
        and not r["reviewed_DBD_pfam_hits"]
    ]
    add_check(checks, "non_DBD_contacts_never_validate_DBD_fragment", not invalid_contact_valid, invalid_contact_valid, [])
    short_high = [
        r["failed_accession"] for r in alignment_rows
        if int(r["DBD_aligned_residues"]) < 30 and r["alignment_evidence_class"] == "HIGH_CONFIDENCE_DBD_PRESENT"
    ]
    add_check(checks, "no_short_alignment_high_confidence", not short_high, short_high, [])
    unreconciled = [
        (r["failed_accession"], r["pdb_reference_reconciliation_status"])
        for r in alignment_rows if r["model_path"] and r["pdb_reference_reconciliation_status"] not in {
            "RECONCILED", "RECONCILED_SHORT_NOT_HIGH_CONFIDENCE"
        }
    ]
    add_check(checks, "every_reference_PDB_sequence_reconciled_or_short_explicit", not unreconciled, unreconciled, [])
    substitutions = [
        r["failed_accession"] for r in alignment_rows
        if r["canonical_reference_is_separate_accession"] != "YES"
        or r["canonical_reference_accession"] == r["failed_accession"]
    ]
    add_check(checks, "canonical_references_explicitly_separate", not substitutions, substitutions, [])
    add_check(checks, "all_uncertainty_visible",
              all(r["manual_review_required"] == "YES" for r in result_rows if r["primary_structural_status"] == "UNCERTAIN_MANUAL_REVIEW"),
              [r["failed_accession"] for r in result_rows if r["primary_structural_status"] == "UNCERTAIN_MANUAL_REVIEW" and r["manual_review_required"] != "YES"], [])
    add_check(checks, "database_unchanged_after_analysis", database_before == database_after_analysis, database_after_analysis, database_before)
    add_check(checks, "upstream_inventory_qc_passed", inventory_qc.get("all_checks_passed") is True, inventory_qc.get("all_checks_passed"), True)

    # Required named controls.
    add_check(checks, "control_E9PN75_lacks_DBD_separate_valid_Q9NZC4",
              by_accession["E9PN75"]["primary_structural_status"] == "ACCESSION_LACKS_DBD"
              and by_accession["E9PN75"]["canonical_reference_accession"] == "Q9NZC4"
              and by_accession["E9PN75"]["canonical_reference_status"] == "CANONICAL_REFERENCE_INTERFACE_VALID"
              and by_accession["E9PN75"]["action_for_baldo"] != "SEND_NEW_DBD_FRAGMENT",
              {k: by_accession["E9PN75"][k] for k in ["primary_structural_status", "canonical_reference_accession", "canonical_reference_status", "action_for_baldo"]},
              {"primary": "ACCESSION_LACKS_DBD", "reference": "Q9NZC4", "reference_status": "CANONICAL_REFERENCE_INTERFACE_VALID", "action": "not SEND_NEW_DBD_FRAGMENT"})
    ehf_results = {
        accession: by_accession[accession]["primary_structural_status"]
        for accession in ["E9PPS9", "E9PQR6", "E9PQX0"]
    }
    add_check(checks, "control_other_EHF_isoforms_independently_aligned",
              all(by_accession[a]["accession_DBD_presence_evidence"] != "" for a in ehf_results)
              and len([r for r in alignment_rows if r["failed_accession"] in ehf_results]) == 3,
              ehf_results, "three independent alignment rows")
    add_check(checks, "control_H7C4N4_PF03299_valid_existing",
              by_accession["H7C4N4"]["primary_structural_status"] == "VALID_DBD_FRAGMENT_MODEL"
              and by_accession["H7C4N4"]["action_for_baldo"] == "NO_SEND_USE_EXISTING_DBD_FRAGMENT",
              {k: by_accession["H7C4N4"][k] for k in ["primary_structural_status", "action_for_baldo", "reviewed_DBD_pfam_hits"]},
              {"primary": "VALID_DBD_FRAGMENT_MODEL", "action": "NO_SEND_USE_EXISTING_DBD_FRAGMENT"})
    add_check(checks, "control_B4DHE0_nonDBD_contacts_not_valid_and_STAT2_compared",
              by_accession["B4DHE0"]["primary_structural_status"] != "VALID_DBD_FRAGMENT_MODEL"
              and "PF01017" in by_accession["B4DHE0"]["non_DBD_pfam_hits"]
              and "PF02865" in by_accession["B4DHE0"]["non_DBD_pfam_hits"]
              and by_accession["B4DHE0"]["canonical_reference_accession"] == "P52630",
              {k: by_accession["B4DHE0"][k] for k in ["primary_structural_status", "non_DBD_pfam_hits", "canonical_reference_accession"]},
              {"primary": "not VALID_DBD_FRAGMENT_MODEL", "non_DBD": "PF01017 and PF02865", "reference": "P52630"})
    add_check(checks, "control_Q59EF3_nonDBD_contact_not_valid_and_TEAD1_compared",
              by_accession["Q59EF3"]["primary_structural_status"] != "VALID_DBD_FRAGMENT_MODEL"
              and "PF17725" in by_accession["Q59EF3"]["non_DBD_pfam_hits"]
              and by_accession["Q59EF3"]["canonical_reference_accession"] in {"H0YE88", "P28347"},
              {k: by_accession["Q59EF3"][k] for k in ["primary_structural_status", "non_DBD_pfam_hits", "canonical_reference_accession"]},
              {"primary": "not VALID_DBD_FRAGMENT_MODEL", "non_DBD": "PF17725", "reference_gene": "TEAD1"})
    q9_alignment = next((r for r in alignment_rows if r["failed_accession"] == "Q9Y2W8"), None)
    add_check(checks, "control_Q9Y2W8_long_accession_to_model_match_revalidated",
              bool(q9_alignment)
              and int(q9_alignment["failed_to_pdb_aligned_residues"]) >= 30
              and float(q9_alignment["failed_to_pdb_identity_percent"]) >= 80
              and by_accession["Q9Y2W8"]["accession_DBD_presence_evidence"] in {"HIGH_CONFIDENCE_DBD_ABSENT", "UNCERTAIN_MANUAL_REVIEW"},
              q9_alignment or {}, {"PDB_reference_alignment": ">=30 residues and >=80% identity", "not_prior_label_only": True})
    add_check(checks, "control_B4DNX4_PF11914_manual",
              by_accession["B4DNX4"]["primary_structural_status"] == "UNCERTAIN_MANUAL_REVIEW"
              and by_accession["B4DNX4"]["manual_review_required"] == "YES",
              {k: by_accession["B4DNX4"][k] for k in ["primary_structural_status", "manual_review_required", "uncertain_pfam_hits"]},
              {"primary": "UNCERTAIN_MANUAL_REVIEW", "manual": "YES"})

    failed_checks = [name for name, check in checks.items() if not check["passed"]]
    if failed_checks:
        print("FAILED_QC_CHECKS")
        for name in failed_checks:
            print(name, json.dumps(checks[name], sort_keys=True), sep="\t")
        # Reviewed map was already created atomically with x mode; remove it so
        # a corrected rerun does not overwrite. This is a run-local output and
        # no other requested output has been created yet.
        args.reviewed_map.unlink()
        raise SystemExit("QC failed before structural output write")

    write_tsv(args.output_draft, result_rows, ACCESSION_COLUMNS)
    write_tsv(args.output_manual_review, manual_rows, ACCESSION_COLUMNS)
    write_tsv(args.output_action_summary, result_rows, action_fields)
    write_tsv(args.output_alignments, alignment_rows, ALIGNMENT_COLUMNS)

    with args.output_fasta.open("x", encoding="utf-8") as handle:
        for row in send_rows:
            accession = row["failed_accession"]
            coords = row["failed_coordinates_corresponding_to_DBD"]
            if not coords:
                raise SystemExit(f"SEND_NEW_DBD_FRAGMENT lacks mapped coordinates: {accession}")
            start, end = map(int, coords.split("-"))
            sequence = failed_records[accession]["sequence"][start - 1:end]
            header = (
                f">{accession}|gene={row['gene_names']}|DBD_evidence="
                f"{row['accession_DBD_presence_evidence']}|coordinates={start}-{end}"
                f"|sequence_length={len(sequence)}|decision_reason={row['decision_reason']}"
            )
            handle.write(header + "\n")
            for offset in range(0, len(sequence), 80):
                handle.write(sequence[offset:offset + 80] + "\n")

    fasta_output_records = fasta_records(args.output_fasta)
    send_ids = {r["failed_accession"] for r in send_rows}
    add_check(checks, "draft_FASTA_exactly_SEND_NEW_rows",
              set(fasta_output_records) == send_ids,
              {"fasta": sorted(fasta_output_records), "send_rows": sorted(send_ids)},
              "exact set equality")
    fasta_reconciliation_errors = []
    for row in send_rows:
        start, end = map(int, row["failed_coordinates_corresponding_to_DBD"].split("-"))
        expected = failed_records[row["failed_accession"]]["sequence"][start - 1:end]
        observed = fasta_output_records[row["failed_accession"]]["sequence"]
        if expected != observed:
            fasta_reconciliation_errors.append(row["failed_accession"])
    add_check(checks, "draft_FASTA_sequences_reconcile_to_failed_accessions",
              not fasta_reconciliation_errors, fasta_reconciliation_errors, [])

    database_after_outputs = db_signature(args.database)
    add_check(checks, "database_unchanged_after_outputs", database_before == database_after_outputs, database_after_outputs, database_before)
    output_hashes = {
        str(path): file_sha256(path)
        for path in outputs if path != args.output_qc and path.exists()
    }
    selected_pdb_hashes = {str(path): file_sha256(path) for path in sorted(selected_pdb_paths)}
    all_passed = all(check["passed"] for check in checks.values())
    qc = {
        "all_checks_passed": all_passed,
        "draft_only": True,
        "counts": {
            "primary_structural_status": dict(sorted(Counter(r["primary_structural_status"] for r in result_rows).items())),
            "canonical_reference_status": dict(sorted(Counter(r["canonical_reference_status"] for r in result_rows).items())),
            "action_for_baldo": dict(sorted(Counter(r["action_for_baldo"] for r in result_rows).items())),
            "manual_review": len(manual_rows),
            "send_new_DBD_fragment": len(send_rows),
            "reference_alignments": len(alignment_rows),
        },
        "send_new_DBD_fragment_accessions": sorted(send_ids),
        "manual_review_accessions": sorted(r["failed_accession"] for r in manual_rows),
        "database_signature_before": database_before,
        "database_signature_after": database_after_outputs,
        "input_hashes": input_hashes,
        "selected_reference_pdb_hashes": selected_pdb_hashes,
        "output_hashes_excluding_qc_self_hash": output_hashes,
        "checks": checks,
    }
    with args.output_qc.open("x", encoding="utf-8") as handle:
        json.dump(qc, handle, indent=2, sort_keys=True)
        handle.write("\n")
    if not all_passed:
        raise SystemExit("Final QC failed")

    print(f"Wrote reviewed Pfam map: {args.reviewed_map}")
    print(f"Wrote 102-accession draft: {args.output_draft}")
    print(f"Wrote manual review: {args.output_manual_review}")
    print(f"Wrote Baldo action summary: {args.output_action_summary}")
    print(f"Wrote DRAFT FASTA: {args.output_fasta}")
    print(f"Wrote reference alignments: {args.output_alignments}")
    print(f"Wrote QC: {args.output_qc}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
