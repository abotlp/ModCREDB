#!/usr/bin/env python3
"""Fetch UniProt metadata for TF accessions stored in the SQLite database."""

from __future__ import annotations

import argparse
import csv
import sqlite3
import time
from datetime import datetime, timezone
from io import StringIO
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen


APP_DIR = Path(__file__).resolve().parent
DEFAULT_DB = APP_DIR / "data" / "tf_webdb.sqlite"
UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
FIELDS = [
    "accession",
    "id",
    "gene_names",
    "protein_name",
    "organism_name",
    "organism_id",
    "reviewed",
    "length",
    "annotation_score",
]


CREATE_TABLE_SQL = """
CREATE TABLE IF NOT EXISTS tf_annotation (
    tf_id TEXT PRIMARY KEY REFERENCES tf(tf_id) ON DELETE CASCADE,
    uniprot_accession TEXT NOT NULL,
    entry_name TEXT,
    gene_names TEXT,
    protein_name TEXT,
    organism_name TEXT,
    organism_id INTEGER,
    reviewed INTEGER,
    sequence_length INTEGER,
    annotation_score REAL,
    uniprot_url TEXT,
    fetched_at TEXT NOT NULL
);

CREATE INDEX IF NOT EXISTS idx_tf_annotation_gene_names ON tf_annotation(gene_names);
CREATE INDEX IF NOT EXISTS idx_tf_annotation_organism_name ON tf_annotation(organism_name);
"""


def chunks(items: list[str], size: int) -> list[list[str]]:
    return [items[index : index + size] for index in range(0, len(items), size)]


def fetch_chunk(accessions: list[str], timeout: int) -> list[dict[str, str]]:
    query = " OR ".join(f"accession:{accession}" for accession in accessions)
    params = urlencode(
        {
            "query": query,
            "format": "tsv",
            "fields": ",".join(FIELDS),
            "size": len(accessions),
        }
    )
    request = Request(
        f"{UNIPROT_SEARCH_URL}?{params}",
        headers={"User-Agent": "tf-webdb-local-enrichment/0.1"},
    )
    with urlopen(request, timeout=timeout) as response:
        text = response.read().decode("utf-8")
    if not text.strip():
        return []
    reader = csv.DictReader(StringIO(text), delimiter="\t")
    return list(reader)


def get_accessions(conn: sqlite3.Connection, missing_only: bool, limit: int | None) -> list[str]:
    where = ""
    if missing_only:
        where = "WHERE NOT EXISTS (SELECT 1 FROM tf_annotation WHERE tf_annotation.tf_id = tf.tf_id)"
    sql = f"SELECT tf_id FROM tf {where} ORDER BY tf_id"
    if limit:
        sql += " LIMIT ?"
        rows = conn.execute(sql, (limit,)).fetchall()
    else:
        rows = conn.execute(sql).fetchall()
    return [row[0] for row in rows]


def parse_int(value: str | None) -> int | None:
    if not value:
        return None
    try:
        return int(value)
    except ValueError:
        return None


def parse_float(value: str | None) -> float | None:
    if not value:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def insert_annotations(conn: sqlite3.Connection, rows: list[dict[str, str]], fetched_at: str) -> int:
    inserted = 0
    for row in rows:
        accession = row.get("Entry") or row.get("accession")
        if not accession:
            continue
        reviewed_value = (row.get("Reviewed") or "").strip().lower()
        reviewed = 1 if reviewed_value == "reviewed" else 0 if reviewed_value else None
        conn.execute(
            """
            INSERT OR REPLACE INTO tf_annotation
                (tf_id, uniprot_accession, entry_name, gene_names, protein_name,
                 organism_name, organism_id, reviewed, sequence_length,
                 annotation_score, uniprot_url, fetched_at)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                accession,
                accession,
                row.get("Entry Name") or None,
                row.get("Gene Names") or None,
                row.get("Protein names") or None,
                row.get("Organism") or None,
                parse_int(row.get("Organism (ID)")),
                reviewed,
                parse_int(row.get("Length")),
                parse_float(row.get("Annotation")),
                f"https://www.uniprot.org/uniprotkb/{accession}/entry",
                fetched_at,
            ),
        )
        inserted += 1
    return inserted


def main() -> None:
    parser = argparse.ArgumentParser(description="Enrich TF rows with UniProt metadata.")
    parser.add_argument("--db", type=Path, default=DEFAULT_DB)
    parser.add_argument("--chunk-size", type=int, default=50)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--timeout", type=int, default=30)
    parser.add_argument("--sleep", type=float, default=0.2)
    parser.add_argument(
        "--all",
        action="store_true",
        help="Refresh every TF annotation instead of only missing annotations.",
    )
    args = parser.parse_args()

    conn = sqlite3.connect(args.db)
    try:
        conn.executescript(CREATE_TABLE_SQL)
        conn.commit()
        accessions = get_accessions(conn, missing_only=not args.all, limit=args.limit)
        print(f"Fetching UniProt annotations for {len(accessions)} TFs")
        fetched_at = datetime.now(timezone.utc).isoformat(timespec="seconds")
        done = 0
        for group in chunks(accessions, args.chunk_size):
            try:
                rows = fetch_chunk(group, timeout=args.timeout)
            except (HTTPError, URLError, TimeoutError) as exc:
                print(f"WARNING: failed chunk {group[0]}..{group[-1]}: {exc}")
                continue
            inserted = insert_annotations(conn, rows, fetched_at)
            conn.commit()
            done += len(group)
            print(f"{done}/{len(accessions)} TFs checked; {inserted} annotations returned in last chunk")
            if args.sleep:
                time.sleep(args.sleep)
        total = conn.execute("SELECT COUNT(*) FROM tf_annotation").fetchone()[0]
        conn.execute(
            "INSERT OR REPLACE INTO metadata (key, value) VALUES (?, ?)",
            ("uniprot_annotations", str(total)),
        )
        conn.execute(
            "INSERT OR REPLACE INTO metadata (key, value) VALUES (?, ?)",
            ("uniprot_fetched_at", fetched_at),
        )
        conn.commit()
        print(f"Stored {total} UniProt annotations in {args.db}")
    finally:
        conn.close()


if __name__ == "__main__":
    main()
