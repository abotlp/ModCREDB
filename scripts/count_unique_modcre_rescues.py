#!/usr/bin/env python3

import csv
from collections import Counter
from pathlib import Path

path = Path(
    "external/baldo_model_inventory/"
    "af3_failed_same_gene_modcre_summary_v2.tsv"
)

if not path.is_file():
    raise SystemExit(f"Missing input: {path}")

with path.open() as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))

rescued = [
    row
    for row in rows
    if row["classification"]
    == "SAME_GENE_ALTERNATIVE_HAS_MODCRE_MODEL"
]

genes = [row["gene"] for row in rescued]
accessions = [row["best_model_accession"] for row in rescued]
pdbs = [row["best_model_example_pdb"] for row in rescued]

missing_pdb = [
    row["failed_tf_id"]
    for row in rescued
    if not row["best_model_example_pdb"]
]

print("Candidate rescued failed accessions:", len(rescued))
print("Unique genes represented:", len(set(genes)))
print("Unique best model accessions:", len(set(accessions)))
print("Unique best PDB files:", len(set(pdbs)))
print("Rows with missing best PDB:", len(missing_pdb))

if missing_pdb:
    print("\nMissing best PDB rows:")
    for tf_id in missing_pdb:
        print(tf_id)

print("\nGenes represented by multiple failed accessions:")
for gene, count in sorted(
    Counter(genes).items(),
    key=lambda item: (-item[1], item[0]),
):
    if count > 1:
        print(gene, count)
