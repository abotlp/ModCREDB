#!/usr/bin/env python3

import csv
import textwrap
from pathlib import Path

ROOT = Path("external/baldo_model_inventory")
FAILED = ROOT / "af3_failed_tf_for_fragmentation.tsv"
FASTA = ROOT / "af3_failed_102_full_length.fasta"
IPR = ROOT / "af3_failed_102_pfam.tsv"

OUT_TSV = ROOT / "af3_failed_all_pfam_fragments.tsv"
OUT_FASTA = ROOT / "af3_failed_all_pfam_fragments.fasta"
OUT_NO_HIT = ROOT / "af3_failed_no_pfam_match.tsv"


def read_fasta(path: Path):
    sequences = {}
    metadata = {}
    accession = None
    chunks = []

    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue

            if line.startswith(">"):
                if accession is not None:
                    sequences[accession] = "".join(chunks)

                token, _, description = line[1:].partition(" ")
                fields = token.split("|")
                accession = fields[0].upper()
                chunks = []

                info = {
                    "gene": "",
                    "primary": "",
                    "description": description,
                }

                for field in fields[1:]:
                    if "=" in field:
                        key, value = field.split("=", 1)
                        info[key] = value

                metadata[accession] = info
            else:
                chunks.append(line)

    if accession is not None:
        sequences[accession] = "".join(chunks)

    return sequences, metadata


def main() -> None:
    for required in (FAILED, FASTA, IPR):
        if not required.is_file():
            raise SystemExit(f"Missing required input: {required}")

    with FAILED.open() as handle:
        failed_ids = {
            row["tf_id"].upper()
            for row in csv.DictReader(handle, delimiter="\t")
        }

    sequences, metadata = read_fasta(FASTA)

    missing_sequences = sorted(failed_ids - sequences.keys())
    if missing_sequences:
        raise RuntimeError(
            "Missing sequences: " + ",".join(missing_sequences)
        )

    hits = []
    seen = set()
    invalid = []
    unexpected = set()
    raw_rows = 0

    with IPR.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue

            columns = line.rstrip("\n").split("\t")
            if len(columns) < 8:
                invalid.append((line_number, "too_few_columns"))
                continue

            tf_id = columns[0].split("|", 1)[0].upper()
            analysis = columns[3]
            pfam_id = columns[4]
            pfam_name = columns[5]

            try:
                start = int(columns[6])
                end = int(columns[7])
            except ValueError:
                invalid.append((line_number, "invalid_coordinates"))
                continue

            raw_rows += 1

            if tf_id not in failed_ids:
                unexpected.add(tf_id)
                continue

            sequence = sequences[tf_id]

            if not (1 <= start <= end <= len(sequence)):
                invalid.append(
                    (
                        line_number,
                        f"coordinates_outside_sequence:{tf_id}:{start}-{end}:len={len(sequence)}",
                    )
                )
                continue

            key = (tf_id, pfam_id, start, end)
            if key in seen:
                continue
            seen.add(key)

            fragment = sequence[start - 1:end]
            expected_length = end - start + 1

            if len(fragment) != expected_length:
                raise RuntimeError(
                    f"{tf_id} {pfam_id} {start}-{end}: "
                    f"length mismatch ({len(fragment)} != {expected_length})"
                )

            info = metadata.get(tf_id, {})

            hits.append({
                "fragment_id": f"{tf_id}|{pfam_id}|{start}-{end}",
                "tf_id": tf_id,
                "gene": info.get("gene", ""),
                "primary": info.get("primary", ""),
                "protein_length": len(sequence),
                "analysis": analysis,
                "pfam_id": pfam_id,
                "pfam_name": pfam_name,
                "start": start,
                "end": end,
                "fragment_length": len(fragment),
                "sequence": fragment,
                "reason": "empty_full_length_AF3_protein_DNA_interface",
            })

    if invalid:
        raise RuntimeError(f"Invalid InterProScan rows: {invalid}")

    if unexpected:
        raise RuntimeError(
            "Unexpected TF IDs: " + ",".join(sorted(unexpected))
        )

    if not hits:
        raise RuntimeError("No valid Pfam hits were parsed")

    hits.sort(
        key=lambda row: (
            row["tf_id"],
            int(row["start"]),
            int(row["end"]),
            row["pfam_id"],
        )
    )

    tf_with_hits = {row["tf_id"] for row in hits}
    tf_without_hits = sorted(failed_ids - tf_with_hits)

    if len(tf_with_hits) + len(tf_without_hits) != len(failed_ids):
        raise RuntimeError("Failed-TF reconciliation did not close")

    with OUT_TSV.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=hits[0].keys(),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(hits)

    with OUT_FASTA.open("w") as handle:
        for row in hits:
            handle.write(
                f">{row['fragment_id']}"
                f"|gene={row['gene']}"
                f"|primary={row['primary']}"
                f"|pfam_name={row['pfam_name']}"
                f"|fragment_length={row['fragment_length']}"
                f"|reason=empty_full_length_AF3_interface\n"
            )
            handle.write(
                "\n".join(textwrap.wrap(row["sequence"], 80)) + "\n"
            )

    fields = [
        "tf_id",
        "gene",
        "primary",
        "protein_length",
        "description",
        "reason",
    ]

    with OUT_NO_HIT.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
        )
        writer.writeheader()

        for tf_id in tf_without_hits:
            info = metadata.get(tf_id, {})
            writer.writerow({
                "tf_id": tf_id,
                "gene": info.get("gene", ""),
                "primary": info.get("primary", ""),
                "protein_length": len(sequences[tf_id]),
                "description": info.get("description", ""),
                "reason": "no_Pfam_38.0_match",
            })

    print("Original failed TFs:", len(failed_ids))
    print("Raw InterProScan rows:", raw_rows)
    print("Unique Pfam fragments:", len(hits))
    print("TFs with fragments:", len(tf_with_hits))
    print("TFs without Pfam hit:", len(tf_without_hits))
    print(
        "Reconciliation:",
        len(tf_with_hits) + len(tf_without_hits),
        "of",
        len(failed_ids),
    )
    print("Invalid rows:", len(invalid))
    print("Unexpected TF IDs:", len(unexpected))
    print("Fragment TSV:", OUT_TSV)
    print("Fragment FASTA:", OUT_FASTA)
    print("No-Pfam list:", OUT_NO_HIT)


if __name__ == "__main__":
    main()
