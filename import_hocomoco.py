#!/usr/bin/env python3
"""Add HOCOMOCO v11 motifs to a staging copy of the TF web database.

This script intentionally copies an existing SQLite database first. It adds
HOCOMOCO as extra evidence and optionally records the hierarchical/best
annotation table without deleting any existing JASPAR, CIS-BP, ModCRE, or
AlphaFold evidence.
"""

from __future__ import annotations

import argparse
import csv
import shutil
import sqlite3
from pathlib import Path

from import_db import ensure_motif_ref_semantics_columns, load_source_releases, motif_ref_semantics, parse_meme_qc, split_list_field


APP_DIR = Path(__file__).resolve().parent
DEFAULT_INPUT_DB = APP_DIR / "data" / "tf_webdb.sqlite"
DEFAULT_OUTPUT_DB = APP_DIR / "data" / "tf_webdb_hocomoco_staging.sqlite"

HOCOMOCO_SOURCE = "hocomoco"
HOCOMOCO_EVIDENCE_COLUMNS = (
    ("Identical_PWM", "identical"),
    ("Homologous_PWM", "homologous"),
    ("Relatively_Homologous_PWM", "relative_homologous"),
    ("ModCRE", "modcre"),
    ("AlphaFold", "alphafold"),
)
PRIMARY_LEVEL_TO_EVIDENCE = {
    "Identical_PWM": "identical",
    "Homologous_PWM": "homologous",
    "Relatively_Homologous_PWM": "relative_homologous",
    "ModCRE": "modcre",
    "AlphaFold": "alphafold",
    "AlphaFold_ModCRE": "alphafold",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows: list[dict[str, str]] = []
        for row in reader:
            rows.append({(key or "").strip(): (value or "").strip() for key, value in row.items()})
        return rows


def split_meme_file(path: Path) -> tuple[str, dict[str, str]]:
    text = path.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()
    motif_indices = [idx for idx, line in enumerate(lines) if line.startswith("MOTIF ")]
    if not motif_indices:
        raise ValueError(f"No MOTIF records found in {path}")

    header = "\n".join(lines[: motif_indices[0]]).rstrip()
    motifs: dict[str, str] = {}
    for pos, start in enumerate(motif_indices):
        end = motif_indices[pos + 1] if pos + 1 < len(motif_indices) else len(lines)
        motif_lines = lines[start:end]
        parts = motif_lines[0].split()
        if len(parts) < 2:
            continue
        motif_id = parts[1]
        content = f"{header}\n\n" + "\n".join(motif_lines).strip() + "\n"
        motifs[motif_id] = content
    return header, motifs


def hocomoco_tokens(row: dict[str, str]) -> list[tuple[str, str | None, str]]:
    """Return HOCOMOCO motif IDs together with their evidence provenance.

    Evidence columns are authoritative. Best_PWM_or_model is a summary field,
    so it is only used when a HOCOMOCO token is absent from every evidence
    column and Best_annotation_level unambiguously supplies its tier.
    """
    tokens: list[tuple[str, str | None, str]] = []
    seen_evidence_tokens: set[tuple[str, str, str]] = set()
    evidence_motif_ids: set[str] = set()
    for column, evidence_type in HOCOMOCO_EVIDENCE_COLUMNS:
        for token in split_list_field(row.get(column, "")):
            token = token.strip()
            if ".H11MO." not in token:
                continue
            key = (token, evidence_type, column)
            if key not in seen_evidence_tokens:
                tokens.append(key)
                seen_evidence_tokens.add(key)
            evidence_motif_ids.add(token)

    best_evidence_type = PRIMARY_LEVEL_TO_EVIDENCE.get(row.get("Best_annotation_level", ""))
    for token in split_list_field(row.get("Best_PWM_or_model", "")):
        token = token.strip()
        if ".H11MO." not in token or token in evidence_motif_ids:
            continue
        tokens.append((token, best_evidence_type, "Best_PWM_or_model"))
    return tokens


def ensure_primary_annotation_table(conn: sqlite3.Connection) -> None:
    conn.execute(
        """
        CREATE TABLE IF NOT EXISTS tf_primary_annotation (
            tf_id TEXT PRIMARY KEY REFERENCES tf(tf_id) ON DELETE CASCADE,
            best_annotation_level TEXT,
            best_pwm_or_model TEXT,
            n_nonempty_annotation_columns INTEGER,
            source_table TEXT NOT NULL
        )
        """
    )


def import_hocomoco(
    input_db: Path,
    output_db: Path,
    meme_path: Path,
    annotation_path: Path | None,
    hierarchical_tsv: Path,
) -> None:
    if output_db.exists():
        output_db.unlink()
    output_db.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(input_db, output_db)

    _, motif_contents = split_meme_file(meme_path)
    rows = read_tsv(hierarchical_tsv)
    hocomoco_ids = sorted({motif_id for row in rows for motif_id, _, _ in hocomoco_tokens(row)})

    conn = sqlite3.connect(output_db)
    try:
        conn.execute("PRAGMA foreign_keys = ON")
        conn.execute(
            "INSERT OR IGNORE INTO source (source, label, description) VALUES (?, ?, ?)",
            (HOCOMOCO_SOURCE, "HOCOMOCO", "HOCOMOCO v11 CORE human mononucleotide motifs"),
        )
        ensure_motif_ref_semantics_columns(conn)
        load_source_releases(conn)

        imported_motifs = 0
        for motif_id, content in motif_contents.items():
            parsed = parse_meme_qc(content)
            conn.execute(
                """
                INSERT OR REPLACE INTO motif_file
                    (source, motif_id, member_path, archive_path, content, width, nsites, consensus, matrix_json,
                     matrix_status, matrix_row_count, matrix_expected_width, matrix_row_sum_min,
                     matrix_row_sum_max, matrix_warning)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    HOCOMOCO_SOURCE,
                    motif_id,
                    motif_id,
                    str(meme_path),
                    content,
                    parsed["width"],
                    parsed["nsites"],
                    parsed["consensus"],
                    parsed["matrix_json"],
                    parsed["matrix_status"],
                    parsed["matrix_row_count"],
                    parsed["matrix_expected_width"],
                    parsed["matrix_row_sum_min"],
                    parsed["matrix_row_sum_max"],
                    parsed["matrix_warning"],
                ),
            )
            imported_motifs += 1

        ensure_primary_annotation_table(conn)
        new_tfs = 0
        hocomoco_refs = 0
        missing_refs = 0

        for row in rows:
            tf_id = row["TF_name"]
            family_text = row.get("TF_family", "")
            existed = conn.execute("SELECT 1 FROM tf WHERE tf_id = ?", (tf_id,)).fetchone() is not None
            conn.execute(
                "INSERT OR IGNORE INTO tf (tf_id, family_text) VALUES (?, ?)",
                (tf_id, family_text),
            )
            if not existed:
                new_tfs += 1
            if family_text:
                conn.execute(
                    "UPDATE tf SET family_text = CASE WHEN family_text = '' THEN ? ELSE family_text END WHERE tf_id = ?",
                    (family_text, tf_id),
                )
            for family in split_list_field(family_text.replace(",", ";")):
                conn.execute(
                    "INSERT OR IGNORE INTO tf_family (tf_id, family) VALUES (?, ?)",
                    (tf_id, family),
                )

            conn.execute(
                """
                INSERT OR REPLACE INTO tf_primary_annotation
                    (tf_id, best_annotation_level, best_pwm_or_model,
                     n_nonempty_annotation_columns, source_table)
                VALUES (?, ?, ?, ?, ?)
                """,
                (
                    tf_id,
                    row.get("Best_annotation_level", ""),
                    row.get("Best_PWM_or_model", ""),
                    int(row.get("N_nonempty_annotation_columns") or 0),
                    str(hierarchical_tsv),
                ),
            )

            for motif_id, evidence_type, original_column in hocomoco_tokens(row):
                if evidence_type is None:
                    conn.execute(
                        """
                        INSERT INTO import_issue
                            (severity, category, message, tf_id, source, motif_id)
                        VALUES (?, ?, ?, ?, ?, ?)
                        """,
                        (
                            "warning",
                            "unresolved_hocomoco_evidence_tier",
                            "Best_PWM_or_model references a HOCOMOCO motif, but Best_annotation_level does not identify an evidence tier.",
                            tf_id,
                            HOCOMOCO_SOURCE,
                            motif_id,
                        ),
                    )
                    continue
                missing = 0 if motif_id in motif_contents else 1
                existing_ref = conn.execute(
                    """
                    SELECT 1
                    FROM motif_ref
                    WHERE tf_id = ? AND evidence_type = ? AND source = ? AND motif_id = ?
                    """,
                    (tf_id, evidence_type, HOCOMOCO_SOURCE, motif_id),
                ).fetchone()
                semantics = motif_ref_semantics(evidence_type, HOCOMOCO_SOURCE, original_column)
                if existing_ref is None:
                    conn.execute(
                        """
                        INSERT INTO motif_ref
                            (tf_id, evidence_type, source, motif_id, original_value,
                             identity_percent, missing_local_file, original_column, mapping_type,
                             curation_status, evidence_note, display_priority)
                        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                        """,
                        (
                            tf_id,
                            evidence_type,
                            HOCOMOCO_SOURCE,
                            motif_id,
                            motif_id,
                            None,
                            missing,
                            semantics["original_column"],
                            semantics["mapping_type"],
                            semantics["curation_status"],
                            semantics["evidence_note"],
                            semantics["display_priority"],
                        ),
                    )
                    hocomoco_refs += 1
                else:
                    conn.execute(
                        """
                        UPDATE motif_ref
                           SET original_column = ?,
                               mapping_type = ?,
                               curation_status = ?,
                               evidence_note = ?,
                               display_priority = ?
                         WHERE tf_id = ? AND evidence_type = ? AND source = ? AND motif_id = ?
                        """,
                        (
                            semantics["original_column"],
                            semantics["mapping_type"],
                            semantics["curation_status"],
                            semantics["evidence_note"],
                            semantics["display_priority"],
                            tf_id,
                            evidence_type,
                            HOCOMOCO_SOURCE,
                            motif_id,
                        ),
                    )
                if missing:
                    missing_refs += 1
                    conn.execute(
                        """
                        INSERT INTO import_issue
                            (severity, category, message, tf_id, source, motif_id)
                        VALUES (?, ?, ?, ?, ?, ?)
                        """,
                        (
                            "warning",
                            "missing_hocomoco_motif_file",
                            "Hierarchical table references a HOCOMOCO motif absent from the HOCOMOCO MEME file.",
                            tf_id,
                            HOCOMOCO_SOURCE,
                            motif_id,
                        ),
                    )

        conn.execute(
            """
            UPDATE tf
               SET motif_ref_count = (
                       SELECT COUNT(*) FROM motif_ref WHERE motif_ref.tf_id = tf.tf_id
                   )
            """
        )
        conn.execute(
            "INSERT OR REPLACE INTO metadata (key, value) VALUES (?, ?)",
            ("hocomoco_meme_file", str(meme_path)),
        )
        if annotation_path:
            conn.execute(
                "INSERT OR REPLACE INTO metadata (key, value) VALUES (?, ?)",
                ("hocomoco_annotation_file", str(annotation_path)),
            )
        conn.execute(
            "INSERT OR REPLACE INTO metadata (key, value) VALUES (?, ?)",
            ("hocomoco_hierarchical_tsv", str(hierarchical_tsv)),
        )
        conn.commit()

        print(f"Staging DB: {output_db}")
        print(f"HOCOMOCO motifs in MEME file: {imported_motifs}")
        print(f"HOCOMOCO motif IDs referenced by TSV: {len(hocomoco_ids)}")
        print(f"HOCOMOCO motif refs added: {hocomoco_refs}")
        print(f"New TF rows added: {new_tfs}")
        print(f"Missing HOCOMOCO motif refs: {missing_refs}")
    finally:
        conn.close()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Import HOCOMOCO v11 into a staging TF web DB.")
    parser.add_argument("--input-db", type=Path, default=DEFAULT_INPUT_DB)
    parser.add_argument("--output-db", type=Path, default=DEFAULT_OUTPUT_DB)
    parser.add_argument("--hocomoco-meme", type=Path, required=True)
    parser.add_argument("--hocomoco-annotation", type=Path)
    parser.add_argument("--hierarchical-tsv", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    import_hocomoco(
        args.input_db,
        args.output_db,
        args.hocomoco_meme,
        args.hocomoco_annotation,
        args.hierarchical_tsv,
    )


if __name__ == "__main__":
    main()
