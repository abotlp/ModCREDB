#!/usr/bin/env python3

import csv
import re
import sqlite3
from pathlib import Path

ROOT = Path("external/baldo_model_inventory")
MAPPING = ROOT / "af3_blank_gene_candidate_mapping.tsv"
OUT = ROOT / "af3_blank_gene_database_matches.tsv"
DB = Path("data/tf_webdb.sqlite")


def gene_tokens(text):
    return {
        token
        for token in re.split(
            r"[\s,;/]+",
            (text or "").upper(),
        )
        if token
    }


with MAPPING.open() as handle:
    mappings = list(
        csv.DictReader(handle, delimiter="\t")
    )

if len(mappings) != 20:
    raise RuntimeError(
        f"Expected 20 mappings, found {len(mappings)}"
    )

connection = sqlite3.connect(DB)
connection.row_factory = sqlite3.Row

annotations = connection.execute("""
    SELECT
        tf_id,
        gene_names,
        protein_name,
        sequence_length
    FROM tf_annotation
    WHERE tf_id IS NOT NULL
""").fetchall()

connection.close()

output_rows = []

for mapping in mappings:
    tf_id = mapping["tf_id"].upper()
    candidate_field = mapping["candidate_gene"].strip()

    candidate_genes = [
        gene.upper()
        for gene in candidate_field.split(";")
        if gene.strip()
    ]

    matches = []

    for row in annotations:
        tokens = gene_tokens(row["gene_names"])

        matched_genes = [
            gene
            for gene in candidate_genes
            if gene in tokens
        ]

        if matched_genes:
            matches.append({
                "failed_tf_id": tf_id,
                "candidate_gene": candidate_field,
                "mapping_status": mapping["mapping_status"],
                "matched_gene": ",".join(matched_genes),
                "candidate_accession": row["tf_id"],
                "sequence_length": row["sequence_length"],
                "gene_names": row["gene_names"] or "",
                "protein_name": row["protein_name"] or "",
            })

    output_rows.extend(matches)

    longest = sorted(
        matches,
        key=lambda row: (
            -(row["sequence_length"] or 0),
            row["candidate_accession"],
        ),
    )[:3]

    print(
        tf_id,
        f"candidate={candidate_field or 'NONE'}",
        f"status={mapping['mapping_status']}",
        f"database_matches={len(matches)}",
        sep="\t",
    )

    for row in longest:
        print(
            "  ",
            row["candidate_accession"],
            f"length={row['sequence_length']}",
            f"genes={row['gene_names']}",
            sep="\t",
        )

fields = [
    "failed_tf_id",
    "candidate_gene",
    "mapping_status",
    "matched_gene",
    "candidate_accession",
    "sequence_length",
    "gene_names",
    "protein_name",
]

with OUT.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=fields,
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(output_rows)

print()
print("Mapping rows checked:", len(mappings))
print("Database candidate rows written:", len(output_rows))
print("Output:", OUT)
