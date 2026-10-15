#!/usr/bin/env python3

import csv
from collections import defaultdict
from pathlib import Path

ROOT = Path("external/baldo_model_inventory")
MATCHES = ROOT / "af3_blank_gene_database_matches.tsv"
MAPPING = ROOT / "af3_blank_gene_candidate_mapping.tsv"
MODELS = Path("/data/sbi/interchange/boliva/patricia/models")

OUT = ROOT / "af3_blank_gene_modcre_model_inventory.tsv"

with MATCHES.open() as handle:
    matches = list(csv.DictReader(handle, delimiter="\t"))

with MAPPING.open() as handle:
    mappings = list(csv.DictReader(handle, delimiter="\t"))

rows = []
summary = defaultdict(lambda: {
    "candidate_accessions": set(),
    "accessions_with_pdb": set(),
    "pdb_count": 0,
})

for row in matches:
    failed_id = row["failed_tf_id"]
    accession = row["candidate_accession"].upper()

    tfs_pdbs = sorted(MODELS.glob(f"TFS_{accession}:*.pdb"))
    dimer_pdbs = sorted(MODELS.glob(f"DIMER_{accession}:*.pdb"))
    all_pdbs = tfs_pdbs + dimer_pdbs

    rows.append({
        **row,
        "tfs_pdb_count": len(tfs_pdbs),
        "dimer_pdb_count": len(dimer_pdbs),
        "total_pdb_count": len(all_pdbs),
        "example_pdb": str(all_pdbs[0]) if all_pdbs else "",
    })

    summary[failed_id]["candidate_accessions"].add(accession)

    if all_pdbs:
        summary[failed_id]["accessions_with_pdb"].add(accession)
        summary[failed_id]["pdb_count"] += len(all_pdbs)

if rows:
    with OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=rows[0].keys(),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)

print("Blank-gene mappings:", len(mappings))
print("Database candidate rows:", len(matches))
print()

for mapping in mappings:
    tf_id = mapping["tf_id"]
    info = summary[tf_id]

    print(
        tf_id,
        f"candidate_gene={mapping['candidate_gene'] or 'NONE'}",
        f"candidate_accessions={len(info['candidate_accessions'])}",
        f"accessions_with_ModCRE_PDB={len(info['accessions_with_pdb'])}",
        f"total_ModCRE_PDBs={info['pdb_count']}",
        sep="\t",
    )

    if info["accessions_with_pdb"]:
        print(
            "  model_accessions=" +
            ",".join(sorted(info["accessions_with_pdb"]))
        )

print()
print("Output:", OUT)
