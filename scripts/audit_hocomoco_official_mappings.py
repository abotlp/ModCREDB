#!/usr/bin/env python3
"""Confirm HOCOMOCO v11 TF-motif links against official local source files.

This audit is read-only with respect to SQLite. It compares ModCREDB's
HOCOMOCO motif_ref links with the official HOCOMOCO v11 annotation/MEME files
and with the integrated hierarchical chart, then writes TSV and Markdown
reports for review.
"""

from __future__ import annotations

import argparse
import csv
import re
import sqlite3
from collections import Counter
from pathlib import Path


HOCOMOCO_TOKEN = re.compile(r"\.H11MO\.")
EVIDENCE_COLUMNS = {
    "Identical_PWM": "identical",
    "Homologous_PWM": "homologous",
    "Relatively_Homologous_PWM": "relative_homologous",
    "ModCRE": "modcre",
    "AlphaFold": "alphafold",
}
OUTPUT_COLUMNS = (
    "db_tf_id",
    "db_entry_name",
    "db_gene_names",
    "hocomoco_motif_id",
    "official_transcription_factor",
    "official_uniprot_id",
    "official_uniprot_ac",
    "official_model_release",
    "official_data_source",
    "official_model_in_annotation",
    "official_model_in_meme",
    "db_motif_file_exists",
    "matrix_status",
    "chart_evidence_column",
    "db_evidence_type",
    "db_mapping_type",
    "uniprot_accession_match",
    "entry_name_match",
    "gene_token_match",
    "mapping_result",
    "notes",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", type=Path, required=True)
    parser.add_argument("--hocomoco-annotation", type=Path, required=True)
    parser.add_argument("--hocomoco-meme", type=Path, required=True)
    parser.add_argument("--hierarchical-tsv", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return [
            {(key or "").strip(): (value or "").strip() for key, value in row.items()}
            for row in csv.DictReader(handle, delimiter="\t")
        ]


def meme_models(path: Path) -> set[str]:
    return {
        match.group(1)
        for line in path.read_text(encoding="utf-8", errors="replace").splitlines()
        if (match := re.match(r"^MOTIF\s+(\S+)", line))
    }


def split_tokens(value: str) -> list[str]:
    return [token.strip() for token in re.split(r"[;,]", value) if token.strip()]


def chart_columns_for_motif(row: dict[str, str], motif_id: str) -> list[str]:
    return [
        column
        for column in EVIDENCE_COLUMNS
        if motif_id in split_tokens(row.get(column, ""))
    ]


def gene_tokens(value: str | None) -> set[str]:
    return {token.upper() for token in re.split(r"[\s;,/]+", value or "") if token}


def write_summary(
    path: Path,
    rows: list[dict[str, str]],
    annotation_count: int,
    meme_count: int,
) -> None:
    results = Counter(row["mapping_result"] for row in rows)
    matrix_statuses = Counter(row["matrix_status"] for row in rows)
    controls = [
        row
        for row in rows
        if (row["db_tf_id"], row["hocomoco_motif_id"])
        in {
            ("P04637", "P53_HUMAN.H11MO.0.A"),
            ("P37231", "PPARG_HUMAN.H11MO.0.A"),
        }
    ]
    mismatch_rows = [row for row in rows if row["mapping_result"] != "exact_uniprot_confirmed"]
    lines = [
        "# HOCOMOCO Official Mapping Audit",
        "",
        "## Scope",
        "",
        "This is a read-only comparison of the promoted ModCREDB staging database, "
        "the local official HOCOMOCO v11 human CORE annotation table, the matching "
        "MEME motif file, and the final integrated hierarchical chart. No SQLite rows "
        "were added, changed, or removed.",
        "",
        "## Summary",
        "",
        f"- Official HOCOMOCO annotation rows: {annotation_count}",
        f"- Official HOCOMOCO MEME motifs: {meme_count}",
        f"- ModCREDB HOCOMOCO links audited: {len(rows)}",
        f"- Exact official UniProt-accession confirmations: {results['exact_uniprot_confirmed']}",
        f"- Non-confirming/mismatch rows: {len(mismatch_rows)}",
        "",
        "## Matrix Status",
        "",
        "| Matrix status | Links |",
        "| --- | ---: |",
        *[f"| {status} | {count} |" for status, count in sorted(matrix_statuses.items())],
        "",
        "## Mapping Results",
        "",
        "| Result | Links |",
        "| --- | ---: |",
        *[f"| {status} | {count} |" for status, count in sorted(results.items())],
        "",
        "`exact_uniprot_confirmed` means the official HOCOMOCO `UniProt AC` equals "
        "the ModCREDB TF accession. The chart-column and evidence-tier checks are "
        "also included in this result.",
        "",
        "## Controls",
        "",
        "| DB TF | HOCOMOCO model | Official UniProt AC | Chart column | Result |",
        "| --- | --- | --- | --- | --- |",
        *[
            f"| {row['db_tf_id']} | {row['hocomoco_motif_id']} | "
            f"{row['official_uniprot_ac']} | {row['chart_evidence_column']} | "
            f"{row['mapping_result']} |"
            for row in controls
        ],
        "",
        "## Interpretation",
        "",
    ]
    if mismatch_rows:
        lines.extend(
            [
                "The rows below require review; they were not changed by this audit.",
                "",
                "| DB TF | Motif | Result | Notes |",
                "| --- | --- | --- | --- |",
                *[
                    f"| {row['db_tf_id']} | {row['hocomoco_motif_id']} | "
                    f"{row['mapping_result']} | {row['notes']} |"
                    for row in mismatch_rows[:50]
                ],
            ]
        )
    else:
        lines.append(
            "All audited HOCOMOCO links are confirmed by the official HOCOMOCO "
            "UniProt accession, matching local MEME matrix, final-chart evidence column, "
            "and stored evidence tier. No database repair is needed."
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    for path in (args.db, args.hocomoco_annotation, args.hocomoco_meme, args.hierarchical_tsv):
        if not path.is_file():
            raise SystemExit(f"Required file not found: {path}")

    official_rows = read_tsv(args.hocomoco_annotation)
    official = {row["Model"]: row for row in official_rows if row.get("Model")}
    official_meme = meme_models(args.hocomoco_meme)
    chart = {row["TF_name"]: row for row in read_tsv(args.hierarchical_tsv) if row.get("TF_name")}

    connection = sqlite3.connect(f"file:{args.db.resolve()}?mode=ro", uri=True)
    connection.row_factory = sqlite3.Row
    try:
        db_rows = connection.execute(
            """
            SELECT
                mr.tf_id,
                mr.motif_id,
                mr.evidence_type,
                mr.mapping_type,
                ta.entry_name,
                ta.gene_names,
                mf.matrix_status,
                mf.motif_id AS db_motif_id
            FROM motif_ref AS mr
            LEFT JOIN motif_file AS mf
              ON mf.source = mr.source AND mf.motif_id = mr.motif_id
            LEFT JOIN tf_annotation AS ta ON ta.tf_id = mr.tf_id
            WHERE mr.source = 'hocomoco'
            ORDER BY mr.tf_id, mr.motif_id
            """
        ).fetchall()
    finally:
        connection.close()

    output_rows: list[dict[str, str]] = []
    for db_row in db_rows:
        motif_id = db_row["motif_id"]
        annotation = official.get(motif_id)
        chart_row = chart.get(db_row["tf_id"])
        chart_columns = chart_columns_for_motif(chart_row or {}, motif_id)
        expected_evidence = {
            EVIDENCE_COLUMNS[column] for column in chart_columns
        }
        official_ac = (annotation or {}).get("UniProt AC", "")
        official_entry = (annotation or {}).get("UniProt ID", "")
        official_tf = (annotation or {}).get("Transcription factor", "")
        accession_match = bool(annotation) and official_ac == db_row["tf_id"]
        entry_match = bool(annotation) and official_entry == (db_row["entry_name"] or "")
        gene_match = bool(annotation) and official_tf.upper() in gene_tokens(db_row["gene_names"])
        notes: list[str] = []
        if not annotation:
            result = "missing_official_annotation"
            notes.append("Model ID is absent from the official annotation table.")
        elif motif_id not in official_meme:
            result = "missing_official_meme"
            notes.append("Model ID is absent from the official MEME file.")
        elif not db_row["db_motif_id"]:
            result = "missing_db_motif_file"
            notes.append("No matching HOCOMOCO motif_file row exists in the database.")
        elif db_row["matrix_status"] != "usable":
            result = "matrix_not_usable"
            notes.append(f"Database matrix_status is {db_row['matrix_status']}.")
        elif not chart_columns:
            result = "chart_mismatch"
            notes.append("Motif is absent from every evidence column for this TF in the integrated chart.")
        elif db_row["evidence_type"] not in expected_evidence:
            result = "evidence_mismatch"
            notes.append(
                f"Chart expects {', '.join(sorted(expected_evidence))}; database stores {db_row['evidence_type']}."
            )
        elif not accession_match:
            result = "uniprot_accession_mismatch"
            notes.append(f"Official UniProt AC is {official_ac or 'missing'}.")
        else:
            result = "exact_uniprot_confirmed"
            if not entry_match:
                notes.append("UniProt accession agrees; entry-name spelling differs or is unavailable locally.")
            if not gene_match:
                notes.append("UniProt accession agrees; HOCOMOCO TF label is not an exact local gene-name token.")
        output_rows.append(
            {
                "db_tf_id": db_row["tf_id"],
                "db_entry_name": db_row["entry_name"] or "",
                "db_gene_names": db_row["gene_names"] or "",
                "hocomoco_motif_id": motif_id,
                "official_transcription_factor": official_tf,
                "official_uniprot_id": official_entry,
                "official_uniprot_ac": official_ac,
                "official_model_release": (annotation or {}).get("Model release", ""),
                "official_data_source": (annotation or {}).get("Data source", ""),
                "official_model_in_annotation": "yes" if annotation else "no",
                "official_model_in_meme": "yes" if motif_id in official_meme else "no",
                "db_motif_file_exists": "yes" if db_row["db_motif_id"] else "no",
                "matrix_status": db_row["matrix_status"] or "missing",
                "chart_evidence_column": ";".join(chart_columns),
                "db_evidence_type": db_row["evidence_type"],
                "db_mapping_type": db_row["mapping_type"],
                "uniprot_accession_match": "yes" if accession_match else "no",
                "entry_name_match": "yes" if entry_match else "no",
                "gene_token_match": "yes" if gene_match else "no",
                "mapping_result": result,
                "notes": " ".join(notes),
            }
        )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)
    write_summary(args.summary, output_rows, len(official), len(official_meme))

    results = Counter(row["mapping_result"] for row in output_rows)
    print(f"Official annotation models: {len(official)}")
    print(f"Official MEME motifs: {len(official_meme)}")
    print(f"ModCREDB HOCOMOCO links audited: {len(output_rows)}")
    for result, count in sorted(results.items()):
        print(f"{result}: {count}")
    print(f"Wrote TSV: {args.out}")
    print(f"Wrote summary: {args.summary}")


if __name__ == "__main__":
    main()
