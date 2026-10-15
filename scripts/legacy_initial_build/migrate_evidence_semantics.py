#!/usr/bin/env python3
"""Backfill explicit TF-motif evidence semantics for existing ModCREDB databases."""

from __future__ import annotations

import argparse
import sqlite3
from pathlib import Path

from import_db import DEFAULT_DB, ensure_motif_ref_semantics_columns, motif_ref_semantics


def migrate(db_path: Path) -> dict[str, int]:
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    try:
        conn.execute("PRAGMA foreign_keys = ON")
        ensure_motif_ref_semantics_columns(conn)
        rows = conn.execute("SELECT id, evidence_type, source FROM motif_ref ORDER BY id").fetchall()
        mapping_counts: dict[str, int] = {}
        for row in rows:
            semantics = motif_ref_semantics(row["evidence_type"], row["source"])
            conn.execute(
                """
                UPDATE motif_ref
                   SET original_column = ?,
                       mapping_type = ?,
                       curation_status = ?,
                       evidence_note = ?,
                       display_priority = ?
                 WHERE id = ?
                """,
                (
                    semantics["original_column"],
                    semantics["mapping_type"],
                    semantics["curation_status"],
                    semantics["evidence_note"],
                    semantics["display_priority"],
                    row["id"],
                ),
            )
            mapping = str(semantics["mapping_type"])
            mapping_counts[mapping] = mapping_counts.get(mapping, 0) + 1
        conn.commit()
        return mapping_counts
    finally:
        conn.close()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Backfill motif_ref evidence mapping and curation fields.")
    parser.add_argument("--db", type=Path, default=DEFAULT_DB)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    counts = migrate(args.db)
    print(f"Backfilled evidence semantics for {args.db}")
    for mapping_type, count in sorted(counts.items()):
        print(f"{mapping_type}	{count}")


if __name__ == "__main__":
    main()
