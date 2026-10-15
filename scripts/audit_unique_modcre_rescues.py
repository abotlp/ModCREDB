#!/usr/bin/env python3

import csv
import re
import subprocess
from collections import Counter
from pathlib import Path

ROOT = Path("external/baldo_model_inventory")

SUMMARY = ROOT / "af3_failed_same_gene_modcre_summary_v2.tsv"
AUDITOR = Path("scripts/audit_one_af3_interface.py")
AUDIT_PYTHON = Path(".conda-af3-audit/bin/python")

OUT = ROOT / "af3_failed_same_gene_modcre_interface_audit.tsv"


def extract_integer(text, label):
    match = re.search(
        rf"^{re.escape(label)}:\s*(\d+)\s*$",
        text,
        flags=re.MULTILINE,
    )
    return int(match.group(1)) if match else ""


def extract_result(text):
    match = re.search(
        r"^RESULT:\s*(\S+)\s*$",
        text,
        flags=re.MULTILINE,
    )
    return match.group(1) if match else ""


for required in (SUMMARY, AUDITOR, AUDIT_PYTHON):
    if not required.is_file():
        raise SystemExit(f"Missing required file: {required}")

with SUMMARY.open() as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))

rescued = [
    row
    for row in rows
    if row["classification"]
    == "SAME_GENE_ALTERNATIVE_HAS_MODCRE_MODEL"
]

groups = {}

for row in rescued:
    pdb_path = row["best_model_example_pdb"].strip()

    if not pdb_path:
        raise RuntimeError(
            f"Missing PDB path for {row['failed_tf_id']}"
        )

    group = groups.setdefault(
        pdb_path,
        {
            "best_model_accession": row["best_model_accession"],
            "best_model_type": row["best_model_type"],
            "best_template": row["best_template"],
            "best_model_domain": row["best_model_domain"],
            "bestbest_template": row["best_template"],
            "best_model_domain": row["best_model_domain"],
            "best_identity_percent": row["best_identity_percent"],
            "best_similarity_percent": row["best_similarity_percent"],
            "genes": set(),
            "failed_tf_ids": set(),
        },
    )

    group["genes"].add(row["gene"])
    group["failed_tf_ids"].add(row["failed_tf_id"])

if len(rescued) != 74:
    raise RuntimeError(
        f"Expected 74 rescued accession rows, found {len(rescued)}"
    )

if len(groups) != 37:
    raise RuntimeError(
        f"Expected 37 unique PDB files, found {len(groups)}"
    )

audit_rows = []

for index, (pdb_path, group) in enumerate(
    sorted(groups.items()),
    start=1,
):
    pdb = Path(pdb_path)

    print(
        f"[{index}/37] Auditing "
        f"{group['best_model_accession']} "
        f"{pdb.name}"
    )

    if not pdb.is_file():
        audit_rows.append({
            "best_model_accession": group["best_model_accession"],
            "genes": ",".join(sorted(group["genes"])),
            "failed_tf_ids": ",".join(
                sorted(group["failed_tf_ids"])
            ),
            "failed_accession_count": len(
                group["failed_tf_ids"]
            ),
            "best_model_type": group["best_model_type"],
            "best_template": group["best_template"],
            "best_model_domain": group["best_model_domain"],
            "best_identity_percent": group["best_identity_percent"],
            "best_similarity_percent": group[
                "best_similarity_percent"
            ],
            "pdb_path": pdb_path,
            "atom_contacts": "",
            "protein_interface_residues": "",
            "dna_interface_residues": "",
            "result": "ERROR_MISSING_PDB",
            "returncode": "",
            "stderr": "",
        })
        continue

    completed = subprocess.run(
        [
            str(AUDIT_PYTHON),
            str(AUDITOR),
            str(pdb),
        ],
        text=True,
        capture_output=True,
    )

    result = extract_result(completed.stdout)

    if completed.returncode != 0:
        result = "ERROR_AUDIT_FAILED"
    elif not result:
        result = "ERROR_RESULT_NOT_PARSED"

    audit_rows.append({
        "best_model_accession": group["best_model_accession"],
        "genes": ",".join(sorted(group["genes"])),
        "failed_tf_ids": ",".join(
            sorted(group["failed_tf_ids"])
        ),
        "failed_accession_count": len(
            group["failed_tf_ids"]
        ),
        "best_model_type": group["best_model_type"],
        "best_template": group["best_template"],
        "best_model_domain": group["best_model_domain"],
        "best_identity_percent": group["best_identity_percent"],
        "best_similarity_percent": group[
            "best_similarity_percent"
        ],
        "pdb_path": pdb_path,
        "atom_contacts": extract_integer(
            completed.stdout,
            "ATOM_CONTACTS",
        ),
        "protein_interface_residues": extract_integer(
            completed.stdout,
            "PROTEIN_INTERFACE_RESIDUES",
        ),
        "dna_interface_residues": extract_integer(
            completed.stdout,
            "DNA_INTERFACE_RESIDUES",
        ),
        "result": result,
        "returncode": completed.returncode,
        "stderr": " ".join(
            completed.stderr.strip().splitlines()
        ),
    })

with OUT.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=audit_rows[0].keys(),
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(audit_rows)

counts = Counter(row["result"] for row in audit_rows)

print()
print("Rescued accession rows:", len(rescued))
print("Unique PDB files:", len(groups))
print("Audit rows written:", len(audit_rows))
print("Results:")

for result, count in sorted(counts.items()):
    print(f"  {result}: {count}")

print(
    "Audit errors:",
    sum(
        count
        for result, count in counts.items()
        if result.startswith("ERROR")
    ),
)

print("Output:", OUT)
