#!/usr/bin/env python3
"""Build a read-only evidence inventory for the 102 failed AF3 accessions.

This script inventories existing evidence only. It does not classify domains,
assign structural status, run modelling tools, or modify SQLite.
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
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable


ROOT = Path(__file__).resolve().parents[1]
BALDO = ROOT / "external" / "baldo_model_inventory"

DEFAULT_FAILED = BALDO / "af3_failed_tf_for_fragmentation.tsv"
DEFAULT_FAILED_FASTA = BALDO / "af3_failed_102_full_length.fasta"
DEFAULT_PFAM = BALDO / "af3_failed_102_pfam.tsv"
DEFAULT_FRAGMENTS = BALDO / "af3_failed_all_pfam_fragments.tsv"
DEFAULT_NO_PFAM = BALDO / "af3_failed_no_pfam_match.tsv"
DEFAULT_FRAGMENT_SUMMARY = BALDO / "af3_fragment_interface_summary.tsv"
DEFAULT_SAME_GENE = BALDO / "af3_failed_same_gene_modcre_candidates_v2.tsv"
DEFAULT_REFERENCE_AUDIT = (
    BALDO / "af3_failed_same_gene_modcre_interface_audit.tsv"
)
DEFAULT_DB = ROOT / "data" / "tf_webdb.sqlite"
DEFAULT_SOURCE_FASTA = Path(
    "/home/patricia/TF_database_Baldo_data/TF_without_model.fasta"
)
DEFAULT_MODELS = Path("/data/sbi/interchange/boliva/patricia/models")
DEFAULT_FAMILY_TREE = ROOT / "data_sources" / "tf_family_tree.json"

DEFAULT_INVENTORY_OUT = (
    ROOT / "outputs" / "af3_failed_102_domain_inventory.tsv"
)
DEFAULT_PFAM_REVIEW_OUT = (
    ROOT / "outputs" / "af3_failed_102_unique_pfam_review.tsv"
)
DEFAULT_QC_OUT = (
    ROOT / "outputs" / "af3_failed_102_domain_inventory_qc.json"
)

EXPECTED_FAILED = 102
EXPECTED_PFAM_ROWS = 67
EXPECTED_WITH_PFAM = 57
EXPECTED_WITHOUT_PFAM = 45
EXPECTED_FRAGMENT_MODELS = 402
EXPECTED_CONTACTING_FRAGMENTS = 9

INVENTORY_FIELDS = [
    "failed_accession",
    "failed_sequence_length",
    "failed_sequence_sha256",
    "gene_names",
    "protein_description",
    "reviewed_status",
    "PWM_annotation_level",
    "PWM_source_or_model",
    "legacy_family_text",
    "accepted_pfam_status",
    "pfam_id",
    "pfam_name",
    "pfam_start",
    "pfam_end",
    "pfam_fragment_length",
    "pfam_analysis",
    "pfam_evalue",
    "pfam_match_status",
    "pfam_scan_date",
    "interpro_id",
    "interpro_name",
    "pfam_raw_row",
    "fragment_id",
    "fragment_sent_to_baldo",
    "fragment_models_checked",
    "fragment_selected_result",
    "fragment_passing_model_count",
    "fragment_best_model_kind",
    "fragment_best_model_path",
    "fragment_atom_contacts",
    "fragment_protein_interface_residues",
    "fragment_dna_interface_residues",
    "fragment_contact_metric_source",
    "same_gene_candidate_accessions",
    "reviewed_same_gene_candidate_accessions",
    "longer_same_gene_candidate_accessions",
    "same_gene_candidate_pfam_domains",
    "candidate_modcre_model_paths",
    "candidate_modcre_model_intervals",
    "candidate_interface_audit_results",
    "domain_role_review",
    "is_expected_DBD_review",
    "domain_role_evidence",
    "manual_review_notes",
]

PFAM_REVIEW_FIELDS = [
    "pfam_id",
    "pfam_name",
    "failed_accession_count",
    "failed_accessions",
    "fragment_count",
    "contacting_fragment_count",
    "legacy_family_tokens",
    "same_gene_candidate_context",
    "domain_role_review",
    "is_TF_DBD_review",
    "review_evidence",
    "reviewer",
    "review_date",
    "review_status",
]

PFAM_TSV_COLUMNS = [
    "sequence_id",
    "sequence_md5",
    "sequence_length",
    "analysis",
    "pfam_id",
    "pfam_name",
    "start",
    "end",
    "evalue",
    "match_status",
    "scan_date",
    "interpro_id",
    "interpro_name",
    "go_terms",
    "pathways",
]

FORBIDDEN_FINAL_VALUES = {
    "VALID_DBD_FRAGMENT_MODEL",
    "DBD_PRESENT_FRAGMENT_FAILED",
    "DBD_PRESENT_FRAGMENT_NOT_SENT",
    "ACCESSION_LACKS_DBD",
    "CANONICAL_REFERENCE_AVAILABLE",
    "COFACTOR_OR_NON_DNA_BINDING_COMPONENT",
    "UNCERTAIN_MANUAL_REVIEW",
    "SAME_GENE_ALTERNATIVE_HAS_MODCRE_MODEL",
    "EXACT_ACCESSION_HAS_MODCRE_MODEL",
    "NO_SAME_GENE_MODCRE_MODEL_FOUND",
}

MODEL_INTERVAL_RE = re.compile(
    r"^(?:TFS|DIMER)_[^:]+:(\d+):(\d+)_"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--failed", type=Path, default=DEFAULT_FAILED)
    parser.add_argument("--failed-fasta", type=Path, default=DEFAULT_FAILED_FASTA)
    parser.add_argument("--pfam", type=Path, default=DEFAULT_PFAM)
    parser.add_argument("--fragments", type=Path, default=DEFAULT_FRAGMENTS)
    parser.add_argument("--no-pfam", type=Path, default=DEFAULT_NO_PFAM)
    parser.add_argument(
        "--fragment-summary", type=Path, default=DEFAULT_FRAGMENT_SUMMARY
    )
    parser.add_argument("--same-gene", type=Path, default=DEFAULT_SAME_GENE)
    parser.add_argument(
        "--reference-interface-audit",
        type=Path,
        default=DEFAULT_REFERENCE_AUDIT,
    )
    parser.add_argument("--db", type=Path, default=DEFAULT_DB)
    parser.add_argument("--source-fasta", type=Path, default=DEFAULT_SOURCE_FASTA)
    parser.add_argument("--models-dir", type=Path, default=DEFAULT_MODELS)
    parser.add_argument("--family-tree", type=Path, default=DEFAULT_FAMILY_TREE)
    parser.add_argument(
        "--inventory-out", type=Path, default=DEFAULT_INVENTORY_OUT
    )
    parser.add_argument(
        "--pfam-review-out", type=Path, default=DEFAULT_PFAM_REVIEW_OUT
    )
    parser.add_argument("--qc-out", type=Path, default=DEFAULT_QC_OUT)
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Allow replacement of the three requested output files.",
    )
    return parser.parse_args()


def clean(value: Any) -> str:
    return "" if value is None else str(value).strip()


def to_int(value: Any, label: str) -> int:
    try:
        return int(clean(value))
    except ValueError as exc:
        raise ValueError(f"Invalid integer for {label}: {value!r}") from exc


def stable_json(value: Any) -> str:
    if value in ({}, [], None, ""):
        return ""
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def hash_file(path: Path) -> dict[str, Any]:
    digest = hashlib.sha256()
    size = 0
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
            size += len(block)
    return {
        "path": str(path),
        "kind": "file",
        "size_bytes": size,
        "sha256": digest.hexdigest(),
    }


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8", errors="replace") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_pfam_tsv(path: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.rstrip("\n")
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != len(PFAM_TSV_COLUMNS):
                raise ValueError(
                    f"{path}:{line_number}: expected "
                    f"{len(PFAM_TSV_COLUMNS)} columns, found {len(fields)}"
                )
            row = dict(zip(PFAM_TSV_COLUMNS, fields))
            row["raw_row"] = line
            row["accession"] = row["sequence_id"].split("|", 1)[0].upper()
            rows.append(row)
    return rows


def read_fasta(path: Path) -> tuple[dict[str, str], dict[str, dict[str, str]]]:
    sequences: dict[str, str] = {}
    metadata: dict[str, dict[str, str]] = {}
    accession = ""
    chunks: list[str] = []

    def save() -> None:
        if not accession:
            return
        if accession in sequences:
            raise ValueError(f"Duplicate FASTA accession in {path}: {accession}")
        sequences[accession] = "".join(chunks)

    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                save()
                token, _, description = line[1:].partition(" ")
                parts = token.split("|")
                accession = parts[0].upper()
                chunks = []
                info = {"description": description, "raw_header": line[1:]}
                for part in parts[1:]:
                    if "=" in part:
                        key, value = part.split("=", 1)
                        info[key] = value
                metadata[accession] = info
            else:
                chunks.append(line)
    save()
    return sequences, metadata


def connect_read_only(path: Path) -> sqlite3.Connection:
    connection = sqlite3.connect(f"file:{path.resolve()}?mode=ro", uri=True)
    connection.row_factory = sqlite3.Row
    return connection


def load_database_evidence(
    db_path: Path, accessions: set[str]
) -> tuple[
    dict[str, dict[str, Any]],
    dict[str, list[dict[str, Any]]],
]:
    metadata: dict[str, dict[str, Any]] = {}
    pfam: dict[str, list[dict[str, Any]]] = defaultdict(list)
    with connect_read_only(db_path) as connection:
        quick_check = connection.execute("PRAGMA quick_check").fetchone()[0]
        if quick_check != "ok":
            raise RuntimeError(f"SQLite quick_check failed: {quick_check}")

        for accession in sorted(accessions):
            row = connection.execute(
                """
                SELECT
                    tf.tf_id,
                    tf.family_text,
                    ta.gene_names,
                    ta.protein_name,
                    ta.reviewed,
                    ta.sequence_length,
                    tpa.best_annotation_level,
                    tpa.best_pwm_or_model
                FROM tf
                LEFT JOIN tf_annotation AS ta ON ta.tf_id = tf.tf_id
                LEFT JOIN tf_primary_annotation AS tpa ON tpa.tf_id = tf.tf_id
                WHERE tf.tf_id = ?
                """,
                (accession,),
            ).fetchone()
            if row is not None:
                metadata[accession] = dict(row)

        if accessions:
            placeholders = ",".join("?" for _ in accessions)
            query = f"""
                SELECT
                    tf_id, pfam_id, pfam_name, pfam_type,
                    interpro_id, interpro_name, start, end,
                    source, source_release, assignment_method
                FROM tf_pfam_annotation
                WHERE tf_id IN ({placeholders})
                ORDER BY tf_id, COALESCE(start, 999999),
                         COALESCE(end, 999999), pfam_id
            """
            for row in connection.execute(query, sorted(accessions)):
                pfam[clean(row["tf_id"]).upper()].append(dict(row))
    return metadata, pfam


def split_family_tokens(text: str) -> list[str]:
    return sorted({item.strip() for item in text.split(",") if item.strip()})


def normalize_fragment_id(fragment_id: str) -> str:
    return fragment_id.replace("|", "_")


def model_interval(path: Path) -> str:
    match = MODEL_INTERVAL_RE.match(path.name)
    return f"{match.group(1)}:{match.group(2)}" if match else ""


def inventory_candidate_models(
    models_dir: Path, candidate_accessions: Iterable[str]
) -> tuple[
    dict[str, list[str]],
    dict[str, list[str]],
    list[Path],
]:
    paths_by_accession: dict[str, list[str]] = {}
    intervals_by_accession: dict[str, list[str]] = {}
    all_paths: list[Path] = []
    for accession in sorted(set(candidate_accessions)):
        paths = sorted(
            list(models_dir.glob(f"TFS_{accession}:*.pdb"))
            + list(models_dir.glob(f"DIMER_{accession}:*.pdb"))
        )
        if not paths:
            continue
        paths_by_accession[accession] = [str(path) for path in paths]
        intervals_by_accession[accession] = sorted(
            {interval for path in paths if (interval := model_interval(path))}
        )
        all_paths.extend(paths)
    return paths_by_accession, intervals_by_accession, all_paths


def hash_model_manifest(models_dir: Path, model_paths: Iterable[Path]) -> dict[str, Any]:
    entries = []
    total_size = 0
    for path in sorted(set(model_paths)):
        info = hash_file(path)
        total_size += int(info["size_bytes"])
        entries.append(
            "\t".join(
                [
                    str(path.relative_to(models_dir)),
                    str(info["size_bytes"]),
                    str(info["sha256"]),
                ]
            )
        )
    manifest = "\n".join(entries)
    if manifest:
        manifest += "\n"
    return {
        "path": str(models_dir),
        "kind": "directory_manifest_of_candidate_pdb_files",
        "file_count": len(entries),
        "size_bytes": total_size,
        "sha256": sha256_bytes(manifest.encode("utf-8")),
        "manifest_definition": (
            "SHA256 of sorted tab-separated relative_path, size_bytes, "
            "file_sha256 lines for candidate ModCRE PDB files used"
        ),
    }


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            extrasaction="raise",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def check(
    checks: dict[str, dict[str, Any]],
    name: str,
    passed: bool,
    *,
    observed: Any,
    expected: Any,
    detail: str = "",
) -> None:
    checks[name] = {
        "passed": bool(passed),
        "observed": observed,
        "expected": expected,
        "detail": detail,
    }


def main() -> int:
    args = parse_args()

    regular_inputs = {
        "failed_accessions": args.failed,
        "failed_fasta": args.failed_fasta,
        "pfam38_tsv": args.pfam,
        "sent_fragments": args.fragments,
        "no_pfam": args.no_pfam,
        "fragment_interface_summary": args.fragment_summary,
        "same_gene_candidates": args.same_gene,
        "reference_interface_audit": args.reference_interface_audit,
        "sqlite_database": args.db,
        "source_fasta": args.source_fasta,
        "family_tree": args.family_tree,
    }
    for label, path in regular_inputs.items():
        if not path.is_file():
            raise SystemExit(f"Missing required input {label}: {path}")
    if not args.models_dir.is_dir():
        raise SystemExit(f"Missing models directory: {args.models_dir}")

    output_paths = [args.inventory_out, args.pfam_review_out, args.qc_out]
    if len({path.resolve() for path in output_paths}) != len(output_paths):
        raise SystemExit("Output paths must be distinct")
    existing = [path for path in output_paths if path.exists()]
    if existing and not args.overwrite:
        raise SystemExit(
            "Refusing to overwrite existing output(s): "
            + ", ".join(str(path) for path in existing)
        )

    failed_rows = read_tsv(args.failed)
    failed_ids = [clean(row.get("tf_id")).upper() for row in failed_rows]
    failed_set = set(failed_ids)
    failed_sequences, failed_fasta_meta = read_fasta(args.failed_fasta)
    source_sequences, source_fasta_meta = read_fasta(args.source_fasta)
    pfam_rows = read_pfam_tsv(args.pfam)
    fragment_rows = read_tsv(args.fragments)
    no_pfam_rows = read_tsv(args.no_pfam)
    fragment_summary_rows = read_tsv(args.fragment_summary)
    candidate_rows = read_tsv(args.same_gene)
    reference_audit_rows = read_tsv(args.reference_interface_audit)

    with args.family_tree.open(encoding="utf-8") as handle:
        family_tree = json.load(handle)
    if not isinstance(family_tree, list):
        raise ValueError("Family tree must contain a JSON list")

    all_candidate_ids = {
        clean(row.get("candidate_tf_id")).upper()
        for row in candidate_rows
        if clean(row.get("candidate_tf_id"))
    }
    database_ids = failed_set | all_candidate_ids
    db_metadata, db_pfam = load_database_evidence(args.db, database_ids)

    pfam_by_key: dict[tuple[str, str, int, int], dict[str, str]] = {}
    for row in pfam_rows:
        key = (
            row["accession"],
            clean(row["pfam_id"]),
            to_int(row["start"], "Pfam start"),
            to_int(row["end"], "Pfam end"),
        )
        if key in pfam_by_key:
            raise ValueError(f"Duplicate accepted Pfam key: {key}")
        pfam_by_key[key] = row

    fragments_by_key: dict[tuple[str, str, int, int], dict[str, str]] = {}
    fragments_by_accession: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in fragment_rows:
        accession = clean(row.get("tf_id")).upper()
        key = (
            accession,
            clean(row.get("pfam_id")),
            to_int(row.get("start"), "fragment start"),
            to_int(row.get("end"), "fragment end"),
        )
        if key in fragments_by_key:
            raise ValueError(f"Duplicate fragment key: {key}")
        fragments_by_key[key] = row
        fragments_by_accession[accession].append(row)

    interface_by_id: dict[str, dict[str, str]] = {}
    for row in fragment_summary_rows:
        fragment_id = clean(row.get("fragment_id"))
        if fragment_id in interface_by_id:
            raise ValueError(f"Duplicate interface summary ID: {fragment_id}")
        interface_by_id[fragment_id] = row

    candidates_by_failed: dict[str, list[dict[str, str]]] = defaultdict(list)
    candidate_pair_counts: Counter[tuple[str, str]] = Counter()
    for row in candidate_rows:
        failed_id = clean(row.get("failed_tf_id")).upper()
        candidate_id = clean(row.get("candidate_tf_id")).upper()
        candidates_by_failed[failed_id].append(row)
        candidate_pair_counts[(failed_id, candidate_id)] += 1

    audits_by_failed: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in reference_audit_rows:
        for failed_id in clean(row.get("failed_tf_ids")).split(","):
            failed_id = failed_id.strip().upper()
            if failed_id:
                audits_by_failed[failed_id].append(row)

    alternative_candidate_ids = {
        clean(row.get("candidate_tf_id")).upper()
        for row in candidate_rows
        if clean(row.get("candidate_tf_id"))
        and clean(row.get("candidate_tf_id")).upper()
        != clean(row.get("failed_tf_id")).upper()
    }
    (
        candidate_model_paths,
        candidate_model_intervals,
        candidate_model_files,
    ) = inventory_candidate_models(args.models_dir, alternative_candidate_ids)

    ambiguous_joins: list[str] = []
    duplicate_pairs = [
        f"{failed_id}->{candidate_id} ({count} rows)"
        for (failed_id, candidate_id), count in sorted(candidate_pair_counts.items())
        if count != 1
    ]
    if duplicate_pairs:
        ambiguous_joins.append(
            "Duplicate failed/candidate pairs: " + "; ".join(duplicate_pairs)
        )

    missing_source_evidence: list[str] = []
    missing_db_failed = sorted(failed_set - db_metadata.keys())
    if missing_db_failed:
        missing_source_evidence.append(
            "Failed accessions missing database metadata: "
            + ",".join(missing_db_failed)
        )
    missing_source_sequences = sorted(failed_set - source_sequences.keys())
    if missing_source_sequences:
        missing_source_evidence.append(
            "Failed accessions missing source FASTA: "
            + ",".join(missing_source_sequences)
        )

    expected_model_counts: dict[str, int] = defaultdict(int)
    for row in candidate_rows:
        failed_id = clean(row.get("failed_tf_id")).upper()
        candidate_id = clean(row.get("candidate_tf_id")).upper()
        if candidate_id == failed_id:
            continue
        expected_model_counts[candidate_id] = max(
            expected_model_counts[candidate_id],
            to_int(row.get("total_pdb_count") or 0, "total_pdb_count"),
        )
    missing_candidate_models = sorted(
        accession
        for accession, count in expected_model_counts.items()
        if count > 0 and accession not in candidate_model_paths
    )
    if missing_candidate_models:
        missing_source_evidence.append(
            "Candidate accessions report PDBs but no matching files were found: "
            + ",".join(missing_candidate_models)
        )

    candidate_context: dict[str, dict[str, Any]] = {}
    for failed_id in sorted(failed_set):
        failed_length = len(failed_sequences.get(failed_id, ""))
        alternatives = []
        for row in candidates_by_failed.get(failed_id, []):
            candidate_id = clean(row.get("candidate_tf_id")).upper()
            if not candidate_id or candidate_id == failed_id:
                continue
            candidate_meta = db_metadata.get(candidate_id, {})
            candidate_length_text = clean(row.get("candidate_length"))
            candidate_length = (
                to_int(candidate_length_text, "candidate_length")
                if candidate_length_text
                else candidate_meta.get("sequence_length")
            )
            alternatives.append(
                {
                    "accession": candidate_id,
                    "candidate_length": candidate_length,
                    "gene_names": clean(
                        candidate_meta.get("gene_names")
                        or row.get("candidate_gene_names")
                    ),
                    "reviewed_status": (
                        "REVIEWED"
                        if candidate_meta.get("reviewed") == 1
                        else (
                            "UNREVIEWED"
                            if candidate_meta.get("reviewed") == 0
                            else ""
                        )
                    ),
                    "is_longer_than_failed": (
                        bool(candidate_length and candidate_length > failed_length)
                    ),
                    "pfam_domains_raw": clean(row.get("pfam_domains")),
                    "reported_total_pdb_count": to_int(
                        row.get("total_pdb_count") or 0,
                        "reported_total_pdb_count",
                    ),
                    "reported_best_model_domain": clean(
                        row.get("best_model_domain")
                    ),
                    "reported_example_pdb": clean(row.get("example_pdb")),
                }
            )
        alternatives.sort(key=lambda item: item["accession"])
        candidate_ids = [item["accession"] for item in alternatives]
        reviewed_ids = [
            item["accession"]
            for item in alternatives
            if item["reviewed_status"] == "REVIEWED"
        ]
        longer_ids = [
            item["accession"]
            for item in alternatives
            if item["is_longer_than_failed"]
        ]
        pfam_context = {
            item["accession"]: (
                [
                    domain
                    for domain in item["pfam_domains_raw"].split("|")
                    if domain
                ]
            )
            for item in alternatives
            if item["pfam_domains_raw"]
        }
        path_context = {
            accession: candidate_model_paths[accession]
            for accession in candidate_ids
            if accession in candidate_model_paths
        }
        interval_context = {
            accession: candidate_model_intervals[accession]
            for accession in candidate_ids
            if accession in candidate_model_intervals
        }
        audit_context = []
        for row in audits_by_failed.get(failed_id, []):
            audit_context.append(
                {
                    "reference_accession": clean(
                        row.get("best_model_accession")
                    ).upper(),
                    "pdb_path": clean(row.get("pdb_path")),
                    "reported_domain": clean(row.get("best_model_domain")),
                    "result": clean(row.get("result")),
                    "atom_contacts": clean(row.get("atom_contacts")),
                    "protein_interface_residues": clean(
                        row.get("protein_interface_residues")
                    ),
                    "dna_interface_residues": clean(
                        row.get("dna_interface_residues")
                    ),
                }
            )
        audit_context.sort(
            key=lambda item: (
                item["reference_accession"],
                item["pdb_path"],
            )
        )
        candidate_context[failed_id] = {
            "alternatives": alternatives,
            "candidate_ids": candidate_ids,
            "reviewed_ids": reviewed_ids,
            "longer_ids": longer_ids,
            "pfam_context": pfam_context,
            "path_context": path_context,
            "interval_context": interval_context,
            "audit_context": audit_context,
        }

    blank_gene_without_alternatives = sorted(
        failed_id
        for failed_id in failed_set
        if not clean(db_metadata.get(failed_id, {}).get("gene_names"))
        and not candidate_context[failed_id]["candidate_ids"]
    )
    if blank_gene_without_alternatives:
        missing_source_evidence.append(
            "Blank-gene failed accessions have no alternative same-gene "
            "evidence in the declared v2 candidate input: "
            + ",".join(blank_gene_without_alternatives)
        )

    failed_without_reference_interface_audit = sorted(
        failed_set - audits_by_failed.keys()
    )
    if failed_without_reference_interface_audit:
        missing_source_evidence.append(
            "Failed accessions have no row in the declared same-gene "
            "reference interface audit: "
            + ",".join(failed_without_reference_interface_audit)
        )

    inventory_rows: list[dict[str, Any]] = []
    sequence_errors: list[str] = []
    interface_join_errors: list[str] = []
    source_sequence_mismatches: list[str] = []

    for failed_id in sorted(failed_set):
        sequence = failed_sequences.get(failed_id, "")
        fasta_info = failed_fasta_meta.get(failed_id, {})
        db_info = db_metadata.get(failed_id, {})
        source_sequence = source_sequences.get(failed_id)
        if source_sequence is not None and source_sequence != sequence:
            source_sequence_mismatches.append(failed_id)

        common = {
            "failed_accession": failed_id,
            "failed_sequence_length": len(sequence),
            "failed_sequence_sha256": sha256_bytes(sequence.encode("ascii")),
            "gene_names": clean(
                db_info.get("gene_names") or fasta_info.get("gene")
            ),
            "protein_description": clean(
                db_info.get("protein_name") or fasta_info.get("description")
            ),
            "reviewed_status": (
                "REVIEWED"
                if db_info.get("reviewed") == 1
                else ("UNREVIEWED" if db_info.get("reviewed") == 0 else "")
            ),
            "PWM_annotation_level": clean(
                db_info.get("best_annotation_level") or fasta_info.get("primary")
            ),
            "PWM_source_or_model": clean(db_info.get("best_pwm_or_model")),
            "legacy_family_text": clean(
                db_info.get("family_text") or fasta_info.get("family")
            ),
            "same_gene_candidate_accessions": ";".join(
                candidate_context[failed_id]["candidate_ids"]
            ),
            "reviewed_same_gene_candidate_accessions": ";".join(
                candidate_context[failed_id]["reviewed_ids"]
            ),
            "longer_same_gene_candidate_accessions": ";".join(
                candidate_context[failed_id]["longer_ids"]
            ),
            "same_gene_candidate_pfam_domains": stable_json(
                candidate_context[failed_id]["pfam_context"]
            ),
            "candidate_modcre_model_paths": stable_json(
                candidate_context[failed_id]["path_context"]
            ),
            "candidate_modcre_model_intervals": stable_json(
                candidate_context[failed_id]["interval_context"]
            ),
            "candidate_interface_audit_results": stable_json(
                candidate_context[failed_id]["audit_context"]
            ),
            "domain_role_review": "UNREVIEWED",
            "is_expected_DBD_review": "UNREVIEWED",
            "domain_role_evidence": "",
            "manual_review_notes": "",
        }

        accession_fragments = sorted(
            fragments_by_accession.get(failed_id, []),
            key=lambda row: (
                to_int(row.get("start"), "fragment sort start"),
                to_int(row.get("end"), "fragment sort end"),
                clean(row.get("pfam_id")),
            ),
        )
        if not accession_fragments:
            inventory_rows.append(
                {
                    **common,
                    "accepted_pfam_status": "NO_ACCEPTED_PFAM38_HIT",
                    "fragment_sent_to_baldo": "NO",
                }
            )
            continue

        for fragment in accession_fragments:
            start = to_int(fragment.get("start"), "fragment start")
            end = to_int(fragment.get("end"), "fragment end")
            pfam_id = clean(fragment.get("pfam_id"))
            key = (failed_id, pfam_id, start, end)
            raw_pfam = pfam_by_key.get(key)
            expected_fragment = sequence[start - 1 : end]
            if not (
                1 <= start <= end <= len(sequence)
                and expected_fragment == clean(fragment.get("sequence"))
                and len(expected_fragment)
                == to_int(fragment.get("fragment_length"), "fragment length")
                == end - start + 1
            ):
                sequence_errors.append(clean(fragment.get("fragment_id")))

            normalized_id = normalize_fragment_id(
                clean(fragment.get("fragment_id"))
            )
            interface = interface_by_id.get(normalized_id)
            if interface is None:
                interface_join_errors.append(clean(fragment.get("fragment_id")))
                interface = {}

            passing_count = clean(interface.get("passing_model_count"))
            has_passing = bool(passing_count and int(passing_count) > 0)
            if has_passing:
                atom_contacts = clean(interface.get("best_atom_contacts"))
                protein_residues = clean(
                    interface.get("best_protein_interface_residues")
                )
                dna_residues = clean(
                    interface.get("best_dna_interface_residues")
                )
                metric_source = "BEST_PASSING_MODEL"
            else:
                atom_contacts = clean(interface.get("selected_atom_contacts"))
                protein_residues = clean(
                    interface.get("selected_protein_interface_residues")
                )
                dna_residues = clean(
                    interface.get("selected_dna_interface_residues")
                )
                metric_source = "SELECTED_MODEL"

            inventory_rows.append(
                {
                    **common,
                    "accepted_pfam_status": "ACCEPTED_PFAM38_HIT",
                    "pfam_id": pfam_id,
                    "pfam_name": clean(fragment.get("pfam_name")),
                    "pfam_start": start,
                    "pfam_end": end,
                    "pfam_fragment_length": end - start + 1,
                    "pfam_analysis": clean(
                        raw_pfam.get("analysis") if raw_pfam else ""
                    ),
                    "pfam_evalue": clean(
                        raw_pfam.get("evalue") if raw_pfam else ""
                    ),
                    "pfam_match_status": clean(
                        raw_pfam.get("match_status") if raw_pfam else ""
                    ),
                    "pfam_scan_date": clean(
                        raw_pfam.get("scan_date") if raw_pfam else ""
                    ),
                    "interpro_id": clean(
                        raw_pfam.get("interpro_id") if raw_pfam else ""
                    ),
                    "interpro_name": clean(
                        raw_pfam.get("interpro_name") if raw_pfam else ""
                    ),
                    "pfam_raw_row": clean(
                        raw_pfam.get("raw_row") if raw_pfam else ""
                    ),
                    "fragment_id": clean(fragment.get("fragment_id")),
                    "fragment_sent_to_baldo": "YES",
                    "fragment_models_checked": clean(
                        interface.get("models_checked")
                    ),
                    "fragment_selected_result": clean(
                        interface.get("selected_result")
                    ),
                    "fragment_passing_model_count": passing_count,
                    "fragment_best_model_kind": clean(
                        interface.get("best_model_kind")
                    ),
                    "fragment_best_model_path": clean(
                        interface.get("best_cif_path")
                    ),
                    "fragment_atom_contacts": atom_contacts,
                    "fragment_protein_interface_residues": protein_residues,
                    "fragment_dna_interface_residues": dna_residues,
                    "fragment_contact_metric_source": metric_source,
                }
            )

    pfam_inventory_rows = [
        row
        for row in inventory_rows
        if row["accepted_pfam_status"] == "ACCEPTED_PFAM38_HIT"
    ]
    pfam_review_rows: list[dict[str, Any]] = []
    for pfam_id in sorted({row["pfam_id"] for row in pfam_inventory_rows}):
        represented = [
            row for row in pfam_inventory_rows if row["pfam_id"] == pfam_id
        ]
        failed_accessions = sorted(
            {clean(row["failed_accession"]) for row in represented}
        )
        family_tokens = sorted(
            {
                token
                for row in represented
                for token in split_family_tokens(clean(row["legacy_family_text"]))
            }
        )
        context = []
        for failed_id in failed_accessions:
            c = candidate_context[failed_id]
            context.append(
                {
                    "failed_accession": failed_id,
                    "same_gene_candidate_accessions": c["candidate_ids"],
                    "reviewed_same_gene_candidate_accessions": c["reviewed_ids"],
                    "longer_same_gene_candidate_accessions": c["longer_ids"],
                }
            )
        pfam_review_rows.append(
            {
                "pfam_id": pfam_id,
                "pfam_name": sorted(
                    {
                        clean(row["pfam_name"])
                        for row in represented
                        if clean(row["pfam_name"])
                    }
                )[0],
                "failed_accession_count": len(failed_accessions),
                "failed_accessions": ";".join(failed_accessions),
                "fragment_count": len(represented),
                "contacting_fragment_count": sum(
                    int(clean(row["fragment_passing_model_count"]) or 0) > 0
                    for row in represented
                ),
                "legacy_family_tokens": ";".join(family_tokens),
                "same_gene_candidate_context": stable_json(context),
                "domain_role_review": "UNREVIEWED",
                "is_TF_DBD_review": "UNREVIEWED",
                "review_evidence": "",
                "reviewer": "",
                "review_date": "",
                "review_status": "UNREVIEWED",
            }
        )

    no_pfam_ids = {
        clean(row.get("tf_id")).upper() for row in no_pfam_rows
    }
    pfam_accessions = {
        clean(row["failed_accession"]) for row in pfam_inventory_rows
    }
    inventory_no_pfam_ids = {
        clean(row["failed_accession"])
        for row in inventory_rows
        if row["accepted_pfam_status"] == "NO_ACCEPTED_PFAM38_HIT"
    }
    raw_pfam_keys = set(pfam_by_key)
    inventory_pfam_keys = {
        (
            clean(row["failed_accession"]),
            clean(row["pfam_id"]),
            int(row["pfam_start"]),
            int(row["pfam_end"]),
        )
        for row in pfam_inventory_rows
    }
    total_fragment_models = sum(
        to_int(row.get("models_checked"), "models_checked")
        for row in fragment_summary_rows
    )
    contacting_fragment_ids = {
        clean(row.get("fragment_id"))
        for row in fragment_summary_rows
        if to_int(row.get("passing_model_count") or 0, "passing_model_count") > 0
    }

    checks: dict[str, dict[str, Any]] = {}
    check(
        checks,
        "exactly_102_unique_failed_accessions",
        len(failed_rows) == len(failed_set) == EXPECTED_FAILED,
        observed={"rows": len(failed_rows), "unique": len(failed_set)},
        expected={"rows": EXPECTED_FAILED, "unique": EXPECTED_FAILED},
    )
    check(
        checks,
        "failed_fasta_exact_accession_set",
        set(failed_sequences) == failed_set,
        observed=len(failed_sequences),
        expected=EXPECTED_FAILED,
    )
    check(
        checks,
        "failed_sequences_match_source_fasta",
        not source_sequence_mismatches and not missing_source_sequences,
        observed={
            "mismatches": source_sequence_mismatches,
            "missing": missing_source_sequences,
        },
        expected={"mismatches": [], "missing": []},
    )
    check(
        checks,
        "exactly_67_accepted_pfam_rows",
        len(pfam_rows) == len(raw_pfam_keys) == EXPECTED_PFAM_ROWS,
        observed={"rows": len(pfam_rows), "unique_keys": len(raw_pfam_keys)},
        expected=EXPECTED_PFAM_ROWS,
    )
    check(
        checks,
        "exactly_57_accessions_with_accepted_pfam",
        len(pfam_accessions) == EXPECTED_WITH_PFAM,
        observed=len(pfam_accessions),
        expected=EXPECTED_WITH_PFAM,
    )
    check(
        checks,
        "exactly_45_accessions_without_accepted_pfam",
        len(no_pfam_ids)
        == len(inventory_no_pfam_ids)
        == EXPECTED_WITHOUT_PFAM
        and no_pfam_ids == inventory_no_pfam_ids == failed_set - pfam_accessions,
        observed={
            "source_no_pfam": len(no_pfam_ids),
            "inventory_no_pfam": len(inventory_no_pfam_ids),
        },
        expected=EXPECTED_WITHOUT_PFAM,
    )
    check(
        checks,
        "every_accepted_pfam_row_in_inventory",
        raw_pfam_keys == inventory_pfam_keys,
        observed={
            "inventory_keys": len(inventory_pfam_keys),
            "missing": sorted(raw_pfam_keys - inventory_pfam_keys),
            "unexpected": sorted(inventory_pfam_keys - raw_pfam_keys),
        },
        expected={"inventory_keys": EXPECTED_PFAM_ROWS, "missing": [], "unexpected": []},
    )
    check(
        checks,
        "fragment_coordinates_and_sequences_reconcile",
        not sequence_errors,
        observed=sequence_errors,
        expected=[],
    )
    check(
        checks,
        "every_fragment_maps_to_interface_summary",
        not interface_join_errors
        and len(fragment_summary_rows) == EXPECTED_PFAM_ROWS,
        observed={
            "missing": interface_join_errors,
            "summary_rows": len(fragment_summary_rows),
        },
        expected={"missing": [], "summary_rows": EXPECTED_PFAM_ROWS},
    )
    check(
        checks,
        "exactly_402_fragment_models_represented",
        total_fragment_models == EXPECTED_FRAGMENT_MODELS,
        observed=total_fragment_models,
        expected=EXPECTED_FRAGMENT_MODELS,
        detail=(
            "Derived by summing models_checked in the declared "
            "af3_fragment_interface_summary.tsv input."
        ),
    )
    check(
        checks,
        "exactly_nine_contacting_fragments",
        len(contacting_fragment_ids) == EXPECTED_CONTACTING_FRAGMENTS,
        observed={
            "count": len(contacting_fragment_ids),
            "fragment_ids": sorted(contacting_fragment_ids),
        },
        expected=EXPECTED_CONTACTING_FRAGMENTS,
    )
    review_values_ok = all(
        row["domain_role_review"] == "UNREVIEWED"
        and row["is_expected_DBD_review"] == "UNREVIEWED"
        for row in inventory_rows
    ) and all(
        row["domain_role_review"] == "UNREVIEWED"
        and row["is_TF_DBD_review"] == "UNREVIEWED"
        and row["review_status"] == "UNREVIEWED"
        for row in pfam_review_rows
    )
    serialized_rows = json.dumps(
        {"inventory": inventory_rows, "pfam_review": pfam_review_rows},
        sort_keys=True,
    )
    forbidden_found = sorted(
        value for value in FORBIDDEN_FINAL_VALUES if value in serialized_rows
    )
    check(
        checks,
        "no_final_dbd_or_rescue_classification_assigned",
        review_values_ok and not forbidden_found,
        observed={
            "review_fields_are_unreviewed": review_values_ok,
            "forbidden_values_found": forbidden_found,
        },
        expected={
            "review_fields_are_unreviewed": True,
            "forbidden_values_found": [],
        },
    )

    e9_pfams = {
        row["pfam_id"]
        for row in inventory_rows
        if row["failed_accession"] == "E9PN75" and row["pfam_id"]
    }
    check(
        checks,
        "control_E9PN75_has_PF02198_not_PF00178",
        "PF02198" in e9_pfams and "PF00178" not in e9_pfams,
        observed=sorted(e9_pfams),
        expected={"contains": "PF02198", "does_not_contain": "PF00178"},
    )
    q9_failed = any(
        row["failed_accession"] == "Q9NZC4" for row in inventory_rows
    )
    q9_evidence_rows = [
        row
        for row in inventory_rows
        if "Q9NZC4" in row["same_gene_candidate_accessions"].split(";")
    ]
    check(
        checks,
        "control_Q9NZC4_is_separate_reference_evidence_only",
        not q9_failed
        and bool(q9_evidence_rows)
        and all(row["failed_accession"] != "Q9NZC4" for row in q9_evidence_rows),
        observed={
            "is_failed_accession": q9_failed,
            "reference_evidence_row_count": len(q9_evidence_rows),
            "failed_accessions_with_reference": sorted(
                {row["failed_accession"] for row in q9_evidence_rows}
            ),
        },
        expected={"is_failed_accession": False, "reference_evidence_present": True},
    )
    for accession, pfam_id in [
        ("B4DHE0", "PF01017"),
        ("Q59EF3", "PF17725"),
    ]:
        rows = [
            row
            for row in inventory_rows
            if row["failed_accession"] == accession
            and row["pfam_id"] == pfam_id
        ]
        passed = (
            len(rows) == 1
            and rows[0]["domain_role_review"] == "UNREVIEWED"
            and rows[0]["is_expected_DBD_review"] == "UNREVIEWED"
            and int(rows[0]["fragment_passing_model_count"] or 0) > 0
        )
        check(
            checks,
            f"control_{accession}_{pfam_id}_contact_remains_unreviewed",
            passed,
            observed=rows,
            expected={
                "row_count": 1,
                "passing_model_count": ">0",
                "domain_role_review": "UNREVIEWED",
                "is_expected_DBD_review": "UNREVIEWED",
            },
        )

    input_provenance = {
        label: hash_file(path) for label, path in regular_inputs.items()
    }
    input_provenance["models_directory"] = hash_model_manifest(
        args.models_dir, candidate_model_files
    )

    if any(not item["passed"] for item in checks.values()):
        failed_checks = [
            name for name, item in checks.items() if not item["passed"]
        ]
        print(
            "QC failed before outputs were written: " + ", ".join(failed_checks),
            file=sys.stderr,
        )
        return 1

    write_tsv(args.inventory_out, INVENTORY_FIELDS, inventory_rows)
    write_tsv(
        args.pfam_review_out,
        PFAM_REVIEW_FIELDS,
        pfam_review_rows,
    )

    qc = {
        "schema_version": "1.0",
        "purpose": (
            "Read-only evidence inventory; no domain-role, DBD, rescue, "
            "or final structural classification is assigned."
        ),
        "generated_at_utc": datetime.now(timezone.utc)
        .replace(microsecond=0)
        .isoformat(),
        "script": str(Path(__file__).resolve()),
        "sqlite_open_mode": "ro",
        "inputs": input_provenance,
        "outputs": {
            "domain_inventory": {
                "path": str(args.inventory_out),
                "row_count": len(inventory_rows),
            },
            "unique_pfam_review": {
                "path": str(args.pfam_review_out),
                "row_count": len(pfam_review_rows),
            },
            "qc": {"path": str(args.qc_out)},
        },
        "counts": {
            "failed_accessions": len(failed_set),
            "domain_inventory_rows": len(inventory_rows),
            "accepted_pfam_rows": len(pfam_inventory_rows),
            "accessions_with_accepted_pfam": len(pfam_accessions),
            "accessions_without_accepted_pfam": len(inventory_no_pfam_ids),
            "unique_accepted_pfam_ids": len(pfam_review_rows),
            "fragment_models_represented": total_fragment_models,
            "contacting_fragments": len(contacting_fragment_ids),
            "same_gene_candidate_source_rows": len(candidate_rows),
            "alternative_same_gene_candidate_accessions": len(
                alternative_candidate_ids
            ),
            "candidate_modcre_pdb_files_inventoried": len(
                set(candidate_model_files)
            ),
            "reference_interface_audit_source_rows": len(reference_audit_rows),
        },
        "checks": checks,
        "all_checks_passed": all(item["passed"] for item in checks.values()),
        "ambiguous_joins": ambiguous_joins,
        "missing_source_evidence": missing_source_evidence,
        "notes": [
            (
                "The 402-model count is derived from the declared fragment "
                "summary input; interface calculations were not rerun."
            ),
            (
                "Candidate model paths and intervals are inventory evidence "
                "from existing filenames and are not accession substitutions."
            ),
            (
                "The family-tree JSON is provenance-hashed and structurally "
                "validated, but its labels are not used to infer domain roles."
            ),
        ],
    }
    args.qc_out.parent.mkdir(parents=True, exist_ok=True)
    args.qc_out.write_text(
        json.dumps(qc, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    print(f"Domain inventory rows: {len(inventory_rows)}")
    print(f"Unique accepted Pfam IDs: {len(pfam_review_rows)}")
    print(f"QC checks passed: {sum(x['passed'] for x in checks.values())}")
    print(f"QC checks failed: {sum(not x['passed'] for x in checks.values())}")
    print(f"Inventory: {args.inventory_out}")
    print(f"Pfam review: {args.pfam_review_out}")
    print(f"QC: {args.qc_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
