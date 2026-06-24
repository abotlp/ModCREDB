#!/usr/bin/env python3
"""Audit JASPAR motif-file coverage and TF-to-motif links without changing a DB.

The supplied JASPAR archive contains MEME matrices. This script determines
whether it also provides enough metadata to audit direct JASPAR UniProt links,
then compares the current SQLite motif_ref rows to two externally verified
positive controls. It is intentionally read-only with respect to SQLite.
"""

from __future__ import annotations

import argparse
import csv
import re
import sqlite3
import tarfile
from collections import Counter
from pathlib import Path


# These two controls are intentionally explicit. They are not inferred from the
# local MEME archive: MA0106.3 -> P04637 is verified on the public JASPAR page,
# and MA0065.1 -> P37231 is the known local positive control.
EXTERNAL_CONTROL_MAPPINGS = {
    "MA0106.3": {"name": "TP53", "uniprot_id": "P04637"},
    "MA0065.1": {"name": "PPARG", "uniprot_id": "P37231"},
}
MOTIF_LINE = re.compile(r"^MOTIF\s+(\S+)(?:\s+(.+))?$", re.MULTILINE)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", type=Path, required=True)
    parser.add_argument("--jaspar-archive", type=Path, required=True)
    parser.add_argument("--output-tsv", type=Path, required=True)
    parser.add_argument("--output-md", type=Path, required=True)
    return parser.parse_args()


def connect_read_only(path: Path) -> sqlite3.Connection:
    uri = f"file:{path.resolve()}?mode=ro"
    connection = sqlite3.connect(uri, uri=True)
    connection.row_factory = sqlite3.Row
    return connection


def archive_inventory(path: Path) -> tuple[list[str], list[str]]:
    with tarfile.open(path, "r:gz") as archive:
        regular_members = [member.name for member in archive.getmembers() if member.isfile()]
    meme_members = [name for name in regular_members if name.lower().endswith(".meme")]
    metadata_members = [name for name in regular_members if not name.lower().endswith(".meme")]
    return meme_members, metadata_members


def motif_name(content: object, fallback: str) -> str:
    match = MOTIF_LINE.search(str(content or ""))
    if not match:
        return ""
    return (match.group(2) or fallback).strip()


def joined(values: list[object]) -> str:
    return ";".join(sorted({str(value) for value in values if value not in (None, "")}))


def markdown_cell(value: object) -> str:
    return str(value or "").replace("|", "\\|")


def main() -> None:
    args = parse_args()
    if not args.db.is_file():
        raise SystemExit(f"Database not found: {args.db}")
    if not args.jaspar_archive.is_file():
        raise SystemExit(f"JASPAR archive not found: {args.jaspar_archive}")

    meme_members, metadata_members = archive_inventory(args.jaspar_archive)

    with connect_read_only(args.db) as connection:
        tf_ids = {
            row["tf_id"]
            for row in connection.execute("SELECT tf_id FROM tf WHERE tf_id IN (?, ?)", ("P04637", "P37231"))
        }
        raw_rows = connection.execute(
            """
            SELECT
                mf.motif_id,
                mf.content,
                mf.consensus,
                mf.width,
                mf.matrix_status,
                COUNT(mr.id) AS motif_ref_count,
                GROUP_CONCAT(DISTINCT mr.evidence_type) AS evidence_types,
                GROUP_CONCAT(DISTINCT mr.mapping_type) AS mapping_types,
                GROUP_CONCAT(DISTINCT mr.tf_id) AS linked_tf_ids
            FROM motif_file AS mf
            LEFT JOIN motif_ref AS mr
              ON mr.source = mf.source AND mr.motif_id = mf.motif_id
            WHERE mf.source = 'jaspar'
            GROUP BY mf.source, mf.motif_id
            ORDER BY mf.motif_id
            """
        ).fetchall()

    rows: list[dict[str, object]] = []
    for raw in raw_rows:
        motif_id = str(raw["motif_id"])
        linked_tf_ids = str(raw["linked_tf_ids"] or "").split(",") if raw["linked_tf_ids"] else []
        control = EXTERNAL_CONTROL_MAPPINGS.get(motif_id)
        expected_tf = control["uniprot_id"] if control else ""
        link_exists = bool(raw["motif_ref_count"])
        control_link_present = bool(expected_tf and expected_tf in linked_tf_ids)
        if control:
            if expected_tf not in tf_ids:
                status = "missing_tf_record"
            elif not link_exists:
                status = "missing_motif_ref_link"
            elif not control_link_present:
                status = "missing_motif_ref_link"
            elif "identical" not in str(raw["evidence_types"] or "").split(","):
                status = "conflicting_evidence_type"
            else:
                status = "ok_link_present"
        else:
            status = "metadata_unavailable"
        rows.append(
            {
                "motif_id": motif_id,
                "jaspar_name": control["name"] if control else motif_name(raw["content"], motif_id),
                "jaspar_uniprot_id": expected_tf,
                "tf_exists_in_modcredb": "yes" if expected_tf in tf_ids else ("no" if expected_tf else "not_assessed"),
                "motif_file_exists": "yes",
                "motif_ref_exists": "yes" if link_exists else "no",
                "db_evidence_type": joined(str(raw["evidence_types"] or "").split(",")),
                "db_mapping_type": joined(str(raw["mapping_types"] or "").split(",")),
                "matrix_status": str(raw["matrix_status"] or "unknown"),
                "status": status,
                "consensus": str(raw["consensus"] or ""),
                "width": raw["width"] if raw["width"] is not None else "",
            }
        )

    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    args.output_md.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "motif_id",
        "jaspar_name",
        "jaspar_uniprot_id",
        "tf_exists_in_modcredb",
        "motif_file_exists",
        "motif_ref_exists",
        "db_evidence_type",
        "db_mapping_type",
        "matrix_status",
        "status",
    ]
    with args.output_tsv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    total = len(rows)
    linked = sum(row["motif_ref_exists"] == "yes" for row in rows)
    zero_link = total - linked
    statuses = Counter(str(row["status"]) for row in rows)
    zero_link_rows = [row for row in rows if row["motif_ref_exists"] == "no"]
    control_rows = {row["motif_id"]: row for row in rows if row["motif_id"] in EXTERNAL_CONTROL_MAPPINGS}

    lines = [
        "# JASPAR UniProt Link Audit",
        "",
        "## Executive Summary",
        "",
        f"- JASPAR motif files in the staging DB: **{total:,}**.",
        f"- Motif files with at least one existing `motif_ref` link: **{linked:,}**.",
        f"- Motif files with zero existing `motif_ref` links: **{zero_link:,}**.",
        f"- Local archive MEME members: **{len(meme_members):,}**.",
        f"- Local archive metadata members: **{len(metadata_members):,}**.",
        "",
        "The supplied `jaspar.tar.gz` contains MEME motif matrices only; it does not contain a JASPAR profile metadata table mapping motif IDs to UniProt accessions. Therefore, this audit can measure current DB link coverage, but cannot infer a complete expected JASPAR-to-UniProt mapping from the local files alone.",
        "",
        "The local JASPAR archive appears to contain motif matrices but not enough metadata to infer all direct TF-to-UniProt links. Missing links such as `MA0106.3 -> P04637` cannot be repaired from MEME files alone without an external JASPAR metadata table, API, or release download.",
        "",
        "A zero-link motif is not automatically an import error: the local JASPAR archive can include non-human profiles or profiles whose source UniProt record is absent from the provided TF sequence-record set. However, an externally verified control such as `MA0106.3 -> P04637` demonstrates that at least one direct JASPAR link is missing.",
        "",
        "## Archive Inventory",
        "",
        f"- Archive label: `{args.jaspar_archive.name}`.",
        "- Archive motif directory: `jaspar/jaspar_2024/PWMS/`.",
        "- Relevant metadata members: none found.",
        "- The MEME header stores a motif name, for example `MOTIF MA0106.3 TP53`, but does not store a UniProt accession. A motif name alone is not sufficient for reliable TF-record mapping.",
        "",
        "## Positive and Negative Controls",
        "",
        "| Motif | Expected external mapping | Current DB result | Status |",
        "|---|---|---|---|",
    ]
    for motif_id, expected in EXTERNAL_CONTROL_MAPPINGS.items():
        row = control_rows.get(motif_id)
        if row is None:
            current = "motif file missing"
            status = "missing_motif_file"
        else:
            current = (
                f"motif_ref={row['motif_ref_exists']}; evidence={row['db_evidence_type'] or 'none'}; "
                f"mapping={row['db_mapping_type'] or 'none'}"
            )
            status = str(row["status"])
        lines.append(
            f"| `{motif_id}` | `{expected['name']} -> {expected['uniprot_id']}` | {markdown_cell(current)} | `{status}` |"
        )

    lines.extend(
        [
            "",
            "## First 20 JASPAR Motifs With Zero Current TF Links",
            "",
            "| Motif ID | Name from MEME | Consensus | Width | Matrix status |",
            "|---|---|---|---:|---|",
        ]
    )
    for row in zero_link_rows[:20]:
        lines.append(
            f"| `{row['motif_id']}` | {markdown_cell(row['jaspar_name'])} | `{row['consensus']}` | {row['width']} | `{row['matrix_status']}` |"
        )

    lines.extend(
        [
            "",
            "## Status Counts",
            "",
            "| Status | Count |",
            "|---|---:|",
        ]
    )
    for status, count in sorted(statuses.items()):
        lines.append(f"| `{status}` | {count:,} |")

    lines.extend(
        [
            "",
            "## Recommendation",
            "",
            "Create a reproducible JASPAR metadata importer from a specified JASPAR release metadata download or API response. It should retain release provenance, profile name, species, collection, source UniProt accession(s), and a documented rule for matching those accessions to ModCREDB TF records. Do not repair these links with manual SQLite edits.",
            "",
            "Before applying that future importer, use a metadata source matching the `jaspar_2024` archive release. The live JASPAR website can change releases over time, so it should not be silently treated as the provenance source for the local 2024 bundle.",
            "",
            "## Reproducibility",
            "",
            "Generated by `scripts/audit_jaspar_uniprot_links.py` in read-only mode. The TSV contains one row for every JASPAR motif file currently in the staging DB. Only two rows use externally verified control mappings; all other rows are marked `metadata_unavailable` because no local UniProt metadata table exists in the archive.",
        ]
    )
    args.output_md.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"JASPAR motif files: {total}")
    print(f"With motif_ref links: {linked}")
    print(f"With zero motif_ref links: {zero_link}")
    print(f"Archive MEME members: {len(meme_members)}")
    print(f"Archive metadata members: {len(metadata_members)}")
    print(f"Wrote: {args.output_tsv}")
    print(f"Wrote: {args.output_md}")


if __name__ == "__main__":
    main()
