#!/usr/bin/env python3
import csv
import sqlite3
from collections import defaultdict
from pathlib import Path

ROOT = Path("external/baldo_model_inventory")
FAILED = ROOT / "af3_failed_tf_for_fragmentation.tsv"
FASTA = Path("/home/patricia/TF_database_Baldo_data/TF_without_model.fasta")
DB = Path("data/tf_webdb.sqlite")
MODELS = Path("/data/sbi/interchange/boliva/patricia/models")

OUT_SUMMARY = ROOT / "af3_failed_same_gene_modcre_summary.tsv"
OUT_CANDIDATES = ROOT / "af3_failed_same_gene_modcre_candidates.tsv"


def read_fasta(path):
    records = {}
    by_gene = defaultdict(list)
    accession = None
    meta = {}
    desc = ""
    chunks = []

    def save():
        if accession is None:
            return
        seq = "".join(chunks)
        gene = meta.get("gene", "").strip()
        records[accession] = {
            "tf_id": accession,
            "gene": gene,
            "gene_key": gene.upper(),
            "primary": meta.get("primary", "").strip(),
            "family": meta.get("family", "").strip(),
            "description": desc,
            "length": len(seq),
        }
        if gene:
            by_gene[gene.upper()].append(accession)

    with path.open() as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                save()
                token, _, desc = line[1:].partition(" ")
                parts = token.split("|")
                accession = parts[0].upper()
                meta = {}
                chunks = []
                for part in parts[1:]:
                    if "=" in part:
                        key, value = part.split("=", 1)
                        meta[key] = value
            else:
                chunks.append(line)
    save()
    return records, by_gene


def load_pfam(path):
    out = defaultdict(set)
    con = sqlite3.connect(path)
    try:
        for tf_id, pfam_id, pfam_name, start, end in con.execute(
            "SELECT tf_id,pfam_id,pfam_name,start,end FROM tf_pfam_annotation"
        ):
            if tf_id:
                out[str(tf_id).upper()].add(
                    f"{pfam_id}:{start}-{end}:{pfam_name or ''}"
                )
    finally:
        con.close()
    return out


def ffloat(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def parse_summary(path, model_type):
    rows = []
    if not path.is_file():
        return rows
    with path.open(errors="replace") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split(";")
            if len(fields) < 13:
                continue
            rows.append({
                "model_type": model_type,
                "template": fields[2],
                "domain": fields[10],
                "identity": ffloat(fields[11]),
                "similarity": ffloat(fields[12]),
            })
    return rows


def inspect_models(accession):
    rows = (
        parse_summary(MODELS / f"TFS_{accession}_model.summary.txt", "TFS")
        + parse_summary(MODELS / f"DIMER_{accession}_model.summary.txt", "DIMER")
    )
    tfs_pdbs = sorted(MODELS.glob(f"TFS_{accession}:*.pdb"))
    dimer_pdbs = sorted(MODELS.glob(f"DIMER_{accession}:*.pdb"))

    def score(row):
        return (
            row["identity"] if row["identity"] is not None else float("-inf"),
            row["similarity"] if row["similarity"] is not None else float("-inf"),
        )

    best = max(rows, key=score) if rows else None
    all_pdbs = tfs_pdbs + dimer_pdbs

    return {
        "tfs_summary_rows": sum(r["model_type"] == "TFS" for r in rows),
        "dimer_summary_rows": sum(r["model_type"] == "DIMER" for r in rows),
        "tfs_pdb_count": len(tfs_pdbs),
        "dimer_pdb_count": len(dimer_pdbs),
        "total_pdb_count": len(all_pdbs),
        "best": best,
        "example_pdb": str(all_pdbs[0]) if all_pdbs else "",
    }


def main():
    for path in (FAILED, FASTA, DB, MODELS):
        if not path.exists():
            raise SystemExit(f"Missing required input: {path}")

    with FAILED.open() as fh:
        failed_ids = sorted({
            row["tf_id"].upper()
            for row in csv.DictReader(fh, delimiter="\t")
        })

    records, by_gene = read_fasta(FASTA)
    pfam = load_pfam(DB)

    missing = sorted(set(failed_ids) - set(records))
    if missing:
        raise RuntimeError(
            "Failed accessions missing from FASTA: " + ",".join(missing)
        )

    model_cache = {}

    def model_info(accession):
        if accession not in model_cache:
            model_cache[accession] = inspect_models(accession)
        return model_cache[accession]

    candidate_rows = []
    summary_rows = []

    for failed_id in failed_ids:
        failed = records[failed_id]
        gene = failed["gene"]

        if gene:
            candidates = sorted(set(by_gene[gene.upper()]))
        else:
            candidates = [failed_id]

        details = []

        for candidate_id in candidates:
            candidate = records[candidate_id]
            info = model_info(candidate_id)
            best = info["best"]

            row = {
                "failed_tf_id": failed_id,
                "gene": gene,
                "failed_length": failed["length"],
                "failed_primary": failed["primary"],
                "failed_family": failed["family"],
                "candidate_tf_id": candidate_id,
                "candidate_length": candidate["length"],
                "candidate_primary": candidate["primary"],
                "candidate_family": candidate["family"],
                "is_failed_accession": int(candidate_id == failed_id),
                "is_longer_than_failed": int(
                    candidate["length"] > failed["length"]
                ),
                "pfam_domains": "|".join(sorted(pfam.get(candidate_id, set()))),
                "tfs_summary_rows": info["tfs_summary_rows"],
                "dimer_summary_rows": info["dimer_summary_rows"],
                "tfs_pdb_count": info["tfs_pdb_count"],
                "dimer_pdb_count": info["dimer_pdb_count"],
                "total_pdb_count": info["total_pdb_count"],
                "best_model_type": best["model_type"] if best else "",
                "best_template": best["template"] if best else "",
                "best_model_domain": best["domain"] if best else "",
                "best_identity_percent": best["identity"] if best else "",
                "best_similarity_percent": best["similarity"] if best else "",
                "example_pdb": info["example_pdb"],
            }
            candidate_rows.append(row)
            details.append(row)

        exact = next(
            row for row in details
            if row["candidate_tf_id"] == failed_id
        )
        alternatives = [
            row for row in details
            if row["candidate_tf_id"] != failed_id
        ]
        alternatives_with_models = [
            row for row in alternatives
            if row["total_pdb_count"] > 0
        ]
        all_with_models = [
            row for row in details
            if row["total_pdb_count"] > 0
        ]

        def candidate_score(row):
            ident = row["best_identity_percent"]
            sim = row["best_similarity_percent"]
            return (
                float(ident) if ident != "" else float("-inf"),
                float(sim) if sim != "" else float("-inf"),
                row["total_pdb_count"],
                row["candidate_length"],
            )

        best_model = (
            max(all_with_models, key=candidate_score)
            if all_with_models else None
        )

        if not gene:
            classification = "BLANK_GENE_REQUIRES_MANUAL_MAPPING"
        elif exact["total_pdb_count"] > 0:
            classification = "EXACT_ACCESSION_HAS_MODCRE_MODEL"
        elif alternatives_with_models:
            classification = "SAME_GENE_ALTERNATIVE_HAS_MODCRE_MODEL"
        else:
            classification = "NO_SAME_GENE_MODCRE_MODEL_FOUND"

        summary_rows.append({
            "failed_tf_id": failed_id,
            "gene": gene,
            "failed_length": failed["length"],
            "failed_primary": failed["primary"],
            "failed_family": failed["family"],
            "same_gene_candidate_count": len(details),
            "same_gene_alternative_count": len(alternatives),
            "longer_same_gene_alternative_count": sum(
                row["is_longer_than_failed"] for row in alternatives
            ),
            "exact_accession_modcre_pdb_count": exact["total_pdb_count"],
            "alternative_accessions_with_modcre": "|".join(
                row["candidate_tf_id"]
                for row in alternatives_with_models
            ),
            "alternative_modcre_pdb_count": sum(
                row["total_pdb_count"]
                for row in alternatives_with_models
            ),
            "best_model_accession": (
                best_model["candidate_tf_id"] if best_model else ""
            ),
            "best_model_type": (
                best_model["best_model_type"] if best_model else ""
            ),
            "best_template": (
                best_model["best_template"] if best_model else ""
            ),
            "best_model_domain": (
                best_model["best_model_domain"] if best_model else ""
            ),
            "best_identity_percent": (
                best_model["best_identity_percent"] if best_model else ""
            ),
            "best_similarity_percent": (
                best_model["best_similarity_percent"] if best_model else ""
            ),
            "best_model_example_pdb": (
                best_model["example_pdb"] if best_model else ""
            ),
            "classification": classification,
        })

    with OUT_CANDIDATES.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=candidate_rows[0].keys(), delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(candidate_rows)

    with OUT_SUMMARY.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=summary_rows[0].keys(), delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(summary_rows)

    counts = defaultdict(int)
    for row in summary_rows:
        counts[row["classification"]] += 1

    explicit = sum(bool(records[x]["gene"]) for x in failed_ids)
    blank = len(failed_ids) - explicit

    print("Failed accessions:", len(failed_ids))
    print("Explicit-gene accessions:", explicit)
    print("Blank-gene accessions:", blank)
    print("Candidate rows written:", len(candidate_rows))
    print("Summary rows written:", len(summary_rows))
    print("Classifications:")
    for key in sorted(counts):
        print(f"  {key}: {counts[key]}")
    print("Summary:", OUT_SUMMARY)
    print("Candidates:", OUT_CANDIDATES)


if __name__ == "__main__":
    main()
