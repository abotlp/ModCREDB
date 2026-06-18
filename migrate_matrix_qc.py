#!/usr/bin/env python3
"""Classify MEME/PWM matrix QC fields for existing ModCREDB SQLite databases."""

from __future__ import annotations

import argparse
import sqlite3
from pathlib import Path

from import_db import DEFAULT_DB, parse_meme_qc

MATRIX_QC_COLUMNS = {
    "matrix_status": "TEXT NOT NULL DEFAULT 'unknown'",
    "matrix_row_count": "INTEGER",
    "matrix_expected_width": "INTEGER",
    "matrix_row_sum_min": "REAL",
    "matrix_row_sum_max": "REAL",
    "matrix_warning": "TEXT",
}


def existing_columns(conn: sqlite3.Connection, table_name: str) -> set[str]:
    return {row[1] for row in conn.execute(f"PRAGMA table_info({table_name})")}


def ensure_matrix_qc_columns(conn: sqlite3.Connection) -> None:
    columns = existing_columns(conn, "motif_file")
    for column, definition in MATRIX_QC_COLUMNS.items():
        if column not in columns:
            conn.execute(f"ALTER TABLE motif_file ADD COLUMN {column} {definition}")


def migrate(db_path: Path) -> dict[str, int]:
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    try:
        conn.execute("PRAGMA foreign_keys = ON")
        ensure_matrix_qc_columns(conn)
        status_counts: dict[str, int] = {}
        rows = conn.execute("SELECT source, motif_id, content FROM motif_file ORDER BY source, motif_id").fetchall()
        for row in rows:
            parsed = parse_meme_qc(row["content"] or "")
            conn.execute(
                """
                UPDATE motif_file
                   SET width = ?,
                       nsites = ?,
                       consensus = ?,
                       matrix_json = ?,
                       matrix_status = ?,
                       matrix_row_count = ?,
                       matrix_expected_width = ?,
                       matrix_row_sum_min = ?,
                       matrix_row_sum_max = ?,
                       matrix_warning = ?
                 WHERE source = ? AND motif_id = ?
                """,
                (
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
                    row["source"],
                    row["motif_id"],
                ),
            )
            status = str(parsed["matrix_status"])
            status_counts[status] = status_counts.get(status, 0) + 1
        conn.commit()
        return status_counts
    finally:
        conn.close()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Classify motif_file matrix QC status fields.")
    parser.add_argument("--db", type=Path, default=DEFAULT_DB)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    counts = migrate(args.db)
    print(f"Classified matrix QC for {args.db}")
    for status, count in sorted(counts.items()):
        print(f"{status}	{count}")


if __name__ == "__main__":
    main()
