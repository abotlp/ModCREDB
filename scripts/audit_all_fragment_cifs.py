#!/usr/bin/env python3

import csv
import re
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(
    "external/baldo_model_inventory/fragments_cif_audit/AF3"
)
PYTHON = Path(".conda-af3-audit/bin/python")
AUDITOR = Path("scripts/audit_one_af3_interface.py")

OUT_MODELS = Path(
    "external/baldo_model_inventory/"
    "af3_fragment_all_model_interface_audit.tsv"
)
OUT_SUMMARY = Path(
    "external/baldo_model_inventory/"
    "af3_fragment_interface_summary.tsv"
)


def parse_int(text, label):
    match = re.search(
        rf"^{re.escape(label)}:\s*(\d+)\s*$",
        text,
        re.MULTILINE,
    )
    return int(match.group(1)) if match else ""


def parse_result(text):
    match = re.search(
        r"^RESULT:\s*(\S+)\s*$",
        text,
        re.MULTILINE,
    )
    return match.group(1) if match else ""


for path in (ROOT, PYTHON, AUDITOR):
    if not path.exists():
        raise SystemExit(f"Missing required path: {path}")

cif_files = sorted(ROOT.rglob("*_model.cif"))

if len(cif_files) != 402:
    raise SystemExit(
        f"Expected 402 CIF files, found {len(cif_files)}"
    )

model_rows = []

for index, cif_path in enumerate(cif_files, start=1):
    relative = cif_path.relative_to(ROOT)
    parts = relative.parts

    if len(parts) == 2:
        fragment_id = parts[0]
        model_kind = "selected"
        sample_name = "selected"
    elif len(parts) == 3:
        fragment_id = parts[0]
        model_kind = "sample"
        sample_name = parts[1]
    else:
        raise RuntimeError(
            f"Unexpected path structure: {relative}"
        )

    completed = subprocess.run(
        [str(PYTHON), str(AUDITOR), str(cif_path)],
        text=True,
        capture_output=True,
    )

    result = parse_result(completed.stdout)

    if completed.returncode != 0:
        result = "ERROR_AUDIT_FAILED"
    elif not result:
        result = "ERROR_RESULT_NOT_PARSED"

    model_rows.append({
        "fragment_id": fragment_id,
        "model_kind": model_kind,
        "sample_name": sample_name,
        "cif_path": str(cif_path),
        "atom_contacts": parse_int(
            completed.stdout,
            "ATOM_CONTACTS",
        ),
        "protein_interface_residues": parse_int(
            completed.stdout,
            "PROTEIN_INTERFACE_RESIDUES",
        ),
        "dna_interface_residues": parse_int(
            completed.stdout,
            "DNA_INTERFACE_RESIDUES",
        ),
        "result": result,
        "returncode": completed.returncode,
        "stderr": " ".join(
            completed.stderr.strip().splitlines()
        ),
    })

    print(
        f"[{index}/402] {fragment_id} {sample_name}: {result}"
    )

with OUT_MODELS.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=model_rows[0].keys(),
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(model_rows)

by_fragment = defaultdict(list)

for row in model_rows:
    by_fragment[row["fragment_id"]].append(row)

if len(by_fragment) != 67:
    raise RuntimeError(
        f"Expected 67 fragments, found {len(by_fragment)}"
    )

summary_rows = []

for fragment_id in sorted(by_fragment):
    rows = by_fragment[fragment_id]

    if len(rows) != 6:
        raise RuntimeError(
            f"{fragment_id}: expected 6 models, found {len(rows)}"
        )

    selected = [
        row for row in rows
        if row["model_kind"] == "selected"
    ]

    if len(selected) != 1:
        raise RuntimeError(
            f"{fragment_id}: expected one selected model"
        )

    selected = selected[0]

    passing = [
        row for row in rows
        if row["result"]
        == "PASS_NONEMPTY_PROTEIN_DNA_INTERFACE"
    ]

    passing_samples = [
        row for row in passing
        if row["model_kind"] == "sample"
    ]

    errors = [
        row for row in rows
        if row["result"].startswith("ERROR")
    ]

    best = None

    if passing:
        best = max(
            passing,
            key=lambda row: (
                int(row["atom_contacts"] or -1),
                int(row["protein_interface_residues"] or -1),
                int(row["dna_interface_residues"] or -1),
            ),
        )

    selected_passes = (
        selected["result"]
        == "PASS_NONEMPTY_PROTEIN_DNA_INTERFACE"
    )

    if errors:
        classification = "AUDIT_ERROR"
    elif selected_passes:
        classification = "SELECTED_MODEL_CONTACTS_DNA"
    elif passing_samples:
        classification = "LOWER_RANKED_SAMPLE_RESCUE"
    else:
        classification = "NO_MODEL_CONTACTS_DNA"

    summary_rows.append({
        "fragment_id": fragment_id,
        "models_checked": len(rows),
        "selected_result": selected["result"],
        "selected_atom_contacts": selected["atom_contacts"],
        "selected_protein_interface_residues":
            selected["protein_interface_residues"],
        "selected_dna_interface_residues":
            selected["dna_interface_residues"],
        "passing_model_count": len(passing),
        "passing_sample_count": len(passing_samples),
        "best_model_kind":
            best["model_kind"] if best else "",
        "best_sample_name":
            best["sample_name"] if best else "",
        "best_atom_contacts":
            best["atom_contacts"] if best else "",
        "best_protein_interface_residues":
            best["protein_interface_residues"] if best else "",
        "best_dna_interface_residues":
            best["dna_interface_residues"] if best else "",
        "best_cif_path":
            best["cif_path"] if best else "",
        "classification": classification,
    })

with OUT_SUMMARY.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=summary_rows[0].keys(),
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(summary_rows)

model_counts = Counter(
    row["result"] for row in model_rows
)
fragment_counts = Counter(
    row["classification"] for row in summary_rows
)

errors = sum(
    count
    for result, count in model_counts.items()
    if result.startswith("ERROR")
)

print()
print("CIF models checked:", len(model_rows))
print("Fragments summarized:", len(summary_rows))
print(
    "Models per fragment:",
    dict(sorted(Counter(
        len(rows) for rows in by_fragment.values()
    ).items())),
)
print("Audit errors:", errors)

print("\nModel results:")
for result, count in sorted(model_counts.items()):
    print(f"  {result}: {count}")

print("\nFragment classifications:")
for classification, count in sorted(fragment_counts.items()):
    print(f"  {classification}: {count}")

print("\nDetailed output:", OUT_MODELS)
print("Fragment summary:", OUT_SUMMARY)
