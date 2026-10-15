#!/usr/bin/env python3
"""Parse ModCRE model summary files into the existing SQLite database."""

from __future__ import annotations

import argparse
import sqlite3
import tarfile
from pathlib import Path

from import_db import insert_model_summary_rows, summary_base_id, update_model_summary_links


APP_DIR = Path(__file__).resolve().parent
DEFAULT_DB = APP_DIR / "data" / "tf_webdb.sqlite"

MODEL_SUMMARY_SCHEMA = """
CREATE TABLE IF NOT EXISTS model_summary (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    summary_file_id INTEGER NOT NULL REFERENCES structure_file(id) ON DELETE CASCADE,
    matched_structure_id INTEGER REFERENCES structure_file(id) ON DELETE SET NULL,
    source TEXT NOT NULL REFERENCES source(source),
    status TEXT NOT NULL,
    tf_id TEXT,
    summary_model_id TEXT NOT NULL,
    model_rank INTEGER,
    n TEXT,
    template_pdb TEXT,
    n_tails TEXT,
    c_tails TEXT,
    protein_chain TEXT,
    dna_chain TEXT,
    identities TEXT,
    coverage TEXT,
    template_by_rmsd TEXT,
    domain TEXT,
    identity_percent REAL,
    similarity_percent REAL
);

CREATE INDEX IF NOT EXISTS idx_model_summary_tf ON model_summary(tf_id);
CREATE INDEX IF NOT EXISTS idx_model_summary_template ON model_summary(template_pdb);
CREATE INDEX IF NOT EXISTS idx_model_summary_matched_structure ON model_summary(matched_structure_id);
"""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Index model summary rows from models.tar.gz.")
    parser.add_argument("--db", type=Path, default=DEFAULT_DB)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    conn = sqlite3.connect(args.db)
    try:
        conn.executescript(MODEL_SUMMARY_SCHEMA)
        conn.execute("DELETE FROM model_summary")
        summary_files = {
            row[1]: {
                "id": row[0],
                "source": row[2],
                "status": row[3],
                "tf_id": row[4],
                "model_id": row[5],
                "archive_path": row[6],
            }
            for row in conn.execute(
                """
                SELECT id, member_path, source, status, tf_id, model_id, archive_path
                FROM structure_file
                WHERE file_type = 'summary'
                """
            )
        }
        archive_paths = sorted({Path(row["archive_path"]) for row in summary_files.values()})
        inserted = 0
        for archive_path in archive_paths:
            wanted = {
                member_path: info
                for member_path, info in summary_files.items()
                if Path(info["archive_path"]) == archive_path
            }
            with tarfile.open(archive_path, "r:gz") as archive:
                for member in archive:
                    if not member.isfile() or member.name not in wanted:
                        continue
                    extracted = archive.extractfile(member)
                    if extracted is None:
                        continue
                    info = wanted[member.name]
                    content = extracted.read().decode("utf-8", errors="replace")
                    inserted += insert_model_summary_rows(
                        conn,
                        info["id"],
                        info["source"],
                        info["status"],
                        info["tf_id"],
                        summary_base_id(info["model_id"]),
                        content,
                    )
        update_model_summary_links(conn)
        linked = conn.execute(
            "SELECT COUNT(*) FROM model_summary WHERE matched_structure_id IS NOT NULL"
        ).fetchone()[0]
        conn.execute(
            "INSERT OR REPLACE INTO metadata (key, value) VALUES (?, ?)",
            ("model_summary_rows", str(inserted)),
        )
        conn.execute(
            "INSERT OR REPLACE INTO metadata (key, value) VALUES (?, ?)",
            ("model_summary_linked_rows", str(linked)),
        )
        conn.commit()
        print(f"Inserted {inserted} model summary rows")
        print(f"{linked} summary rows are linked to model files")
    finally:
        conn.close()


if __name__ == "__main__":
    main()
