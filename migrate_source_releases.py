#!/usr/bin/env python3
"""Add or refresh source-release metadata in an existing ModCREDB SQLite database."""

from __future__ import annotations

import argparse
import sqlite3
from pathlib import Path

from import_db import DEFAULT_DB, DEFAULT_SOURCE_RELEASES, load_source_releases


def migrate(db_path: Path, config_path: Path) -> int:
    conn = sqlite3.connect(db_path)
    try:
        conn.execute("PRAGMA foreign_keys = ON")
        loaded = load_source_releases(conn, config_path)
        conn.commit()
        return loaded
    finally:
        conn.close()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Load ModCREDB source-release metadata into SQLite.")
    parser.add_argument("--db", type=Path, default=DEFAULT_DB)
    parser.add_argument("--config", type=Path, default=DEFAULT_SOURCE_RELEASES)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    loaded = migrate(args.db, args.config)
    print(f"Loaded {loaded} source_release rows into {args.db}")


if __name__ == "__main__":
    main()
