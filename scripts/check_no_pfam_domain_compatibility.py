#!/usr/bin/env python3

import csv
import re
import subprocess
from collections import defaultdict
from pathlib import Path

ROOT = Path("external/baldo_model_inventory")
INPUT = ROOT / "af3_no_pfam_unresolved_10.tsv"
MODELS = Path("/data/sbi/interchange/boliva/patricia/models")

AUDIT_PYTHON = Path(".conda-af3-audit/bin/python")
AUDITOR = Path("scripts/audit_one_af3_interface.py")

OUT_ALL = ROOT / "af3_no_pfam_candidate_modcre_compatibility_all.tsv"
OUT_BEST = ROOT / "af3_no_pfam_candidate_modcre_compatibility_best.tsv"

CASES = {
    "A8MTU8": ("SREBF1", "P36956"),
    "B3KU95": ("HOXA2", "O43364"),
    "B4E2G3": ("CREB3L4", "Q8TEY5"),
    "B7Z5Z3": ("SOX5", "P35711"),
    "Q9Y2W8": ("ESR1", "P03372"),
}

AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D",
    "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
    "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
    "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M",
}


def parse_result(text):
    match = re.search(
        r"^RESULT:\s*(\S+)\s*$",
        text,
        flags=re.MULTILINE,
    )
    return match.group(1) if match else ""


def parse_int(text, label):
    match = re.search(
        rf"^{re.escape(label)}:\s*(\d+)\s*$",
        text,
        flags=re.MULTILINE,
    )
    return int(match.group(1)) if match else 0


def pdb_protein_sequences(path):
    residues = defaultdict(list)
    seen = set()

    with path.open(errors="replace") as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue

            altloc = line[16:17]
            if altloc not in (" ", "A"):
                continue

            residue_name = line[17:20].strip().upper()
            amino_acid = AA3.get(residue_name)

            if not amino_acid:
                continue

            chain = line[21:22].strip() or "_"
            residue_number = line[22:26].strip()
            insertion_code = line[26:27].strip()

            key = (
                chain,
                residue_number,
                insertion_code,
                residue_name,
            )

            if key in seen:
                continue

            seen.add(key)
            residues[chain].append(amino_acid)

    return {
        chain: "".join(sequence)
        for chain, sequence in residues.items()
        if sequence
    }


def local_alignment(sequence_a, sequence_b):
    rows = len(sequence_a) + 1
    columns = len(sequence_b) + 1

    score = [
        [0] * columns
        for _ in range(rows)
    ]
    pointer = [
        [0] * columns
        for _ in range(rows)
    ]

    best_score = 0
    best_position = (0, 0)

    for i in range(1, rows):
        for j in range(1, columns):
            diagonal = (
                score[i - 1][j - 1]
                + (
                    2
                    if sequence_a[i - 1] == sequence_b[j - 1]
                    else -1
                )
            )
            up = score[i - 1][j] - 2
            left = score[i][j - 1] - 2

            best = max(0, diagonal, up, left)
            score[i][j] = best

            if best == 0:
                pointer[i][j] = 0
            elif best == diagonal:
                pointer[i][j] = 1
            elif best == up:
                pointer[i][j] = 2
            else:
                pointer[i][j] = 3

            if best > best_score:
                best_score = best
                best_position = (i, j)

    i, j = best_position
    matches = 0
    aligned_a = 0
    aligned_b = 0
    aligned_columns = 0

    while i > 0 and j > 0 and score[i][j] > 0:
        direction = pointer[i][j]

        if direction == 1:
            aligned_a += 1
            aligned_b += 1
            aligned_columns += 1

            if sequence_a[i - 1] == sequence_b[j - 1]:
                matches += 1

            i -= 1
            j -= 1

        elif direction == 2:
            aligned_a += 1
            aligned_columns += 1
            i -= 1

        elif direction == 3:
            aligned_b += 1
            aligned_columns += 1
            j -= 1

        else:
            break

    identity = (
        100.0 * matches / aligned_columns
        if aligned_columns
        else 0.0
    )

    model_coverage = (
        100.0 * aligned_b / len(sequence_b)
        if sequence_b
        else 0.0
    )

    failed_coverage = (
        100.0 * aligned_a / len(sequence_a)
        if sequence_a
        else 0.0
    )

    return {
        "alignment_score": best_score,
        "matches": matches,
        "aligned_columns": aligned_columns,
        "identity_percent": identity,
        "model_chain_coverage_percent": model_coverage,
        "failed_sequence_coverage_percent": failed_coverage,
    }


with INPUT.open() as handle:
    sequence_rows = {
        row["tf_id"]: row
        for row in csv.DictReader(handle, delimiter="\t")
    }

missing = [
    tf_id
    for tf_id in CASES
    if tf_id not in sequence_rows
]

if missing:
    raise SystemExit(
        "Missing failed sequences: " + ",".join(missing)
    )

all_rows = []

for failed_tf_id, (
    candidate_gene,
    representative_accession,
) in CASES.items():

    failed_sequence = sequence_rows[failed_tf_id]["sequence"]

    pdb_files = sorted(
        MODELS.glob(
            f"TFS_{representative_accession}:*.pdb"
        )
    )

    print(
        failed_tf_id,
        candidate_gene,
        representative_accession,
        f"TFS_models={len(pdb_files)}",
        sep="\t",
    )

    for pdb_path in pdb_files:
        completed = subprocess.run(
            [
                str(AUDIT_PYTHON),
                str(AUDITOR),
                str(pdb_path),
            ],
            text=True,
            capture_output=True,
        )

        result = parse_result(completed.stdout)

        if completed.returncode != 0:
            result = "ERROR_AUDIT_FAILED"
        elif not result:
            result = "ERROR_RESULT_NOT_PARSED"

        chain_sequences = pdb_protein_sequences(pdb_path)

        best_chain = ""
        best_alignment = None
        best_chain_length = 0

        for chain, model_sequence in chain_sequences.items():
            alignment = local_alignment(
                failed_sequence,
                model_sequence,
            )

            ranking = (
                alignment["model_chain_coverage_percent"],
                alignment["identity_percent"],
                alignment["matches"],
            )

            if (
                best_alignment is None
                or ranking > (
                    best_alignment[
                        "model_chain_coverage_percent"
                    ],
                    best_alignment["identity_percent"],
                    best_alignment["matches"],
                )
            ):
                best_chain = chain
                best_alignment = alignment
                best_chain_length = len(model_sequence)

        if best_alignment is None:
            best_alignment = {
                "alignment_score": 0,
                "matches": 0,
                "aligned_columns": 0,
                "identity_percent": 0.0,
                "model_chain_coverage_percent": 0.0,
                "failed_sequence_coverage_percent": 0.0,
            }

        all_rows.append({
            "failed_tf_id": failed_tf_id,
            "candidate_gene": candidate_gene,
            "representative_accession":
                representative_accession,
            "pdb_path": str(pdb_path),
            "interface_result": result,
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
            "best_protein_chain": best_chain,
            "model_chain_length": best_chain_length,
            "alignment_score":
                best_alignment["alignment_score"],
            "alignment_matches":
                best_alignment["matches"],
            "alignment_columns":
                best_alignment["aligned_columns"],
            "identity_percent": round(
                best_alignment["identity_percent"],
                2,
            ),
            "model_chain_coverage_percent": round(
                best_alignment[
                    "model_chain_coverage_percent"
                ],
                2,
            ),
            "failed_sequence_coverage_percent": round(
                best_alignment[
                    "failed_sequence_coverage_percent"
                ],
                2,
            ),
        })

fields = [
    "failed_tf_id",
    "candidate_gene",
    "representative_accession",
    "pdb_path",
    "interface_result",
    "atom_contacts",
    "protein_interface_residues",
    "dna_interface_residues",
    "best_protein_chain",
    "model_chain_length",
    "alignment_score",
    "alignment_matches",
    "alignment_columns",
    "identity_percent",
    "model_chain_coverage_percent",
    "failed_sequence_coverage_percent",
]

with OUT_ALL.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=fields,
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(all_rows)

best_rows = []

for failed_tf_id in CASES:
    candidates = [
        row
        for row in all_rows
        if row["failed_tf_id"] == failed_tf_id
        and row["interface_result"]
        == "PASS_NONEMPTY_PROTEIN_DNA_INTERFACE"
    ]

    candidates.sort(
        key=lambda row: (
            float(row["model_chain_coverage_percent"]),
            float(row["identity_percent"]),
            int(row["atom_contacts"]),
        ),
        reverse=True,
    )

    if candidates:
        best_rows.append(candidates[0])

with OUT_BEST.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=fields,
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(best_rows)

audit_errors = sum(
    row["interface_result"].startswith("ERROR")
    for row in all_rows
)

print()
print("Cases checked:", len(CASES))
print("Models checked:", len(all_rows))
print("Audit errors:", audit_errors)
print(
    "Cases with a DNA-contacting candidate:",
    len(best_rows),
)
print()

for row in best_rows:
    print(
        row["failed_tf_id"],
        f"gene={row['candidate_gene']}",
        f"representative={row['representative_accession']}",
        f"identity={row['identity_percent']}%",
        f"model_coverage="
        f"{row['model_chain_coverage_percent']}%",
        f"failed_sequence_coverage="
        f"{row['failed_sequence_coverage_percent']}%",
        f"contacts={row['atom_contacts']}",
        f"model={Path(row['pdb_path']).name}",
        sep="\t",
    )

print()
print("All models:", OUT_ALL)
print("Best candidates:", OUT_BEST)
