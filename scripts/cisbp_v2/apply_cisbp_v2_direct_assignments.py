#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import shutil
import sqlite3
from datetime import datetime
from pathlib import Path


REQUIRED_COLUMNS = [
    "tf_id",
    "evidence_type",
    "source",
    "motif_id",
    "original_value",
    "identity_percent",
    "missing_local_file",
    "original_column",
    "mapping_type",
    "curation_status",
    "evidence_note",
    "display_priority",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Apply verified direct CIS-BP v2 Homo sapiens TF-to-motif assignments to motif_ref."
    )
    parser.add_argument("--db", default="data/tf_webdb.sqlite", help="SQLite database path")
    parser.add_argument(
        "--candidates",
        default="data/cisbp_v2_direct_human_motif_ref_insert_candidates.tsv",
        help="Verified candidate TSV",
    )
    parser.add_argument("--apply", action="store_true", help="Actually modify the DB. Default is dry-run.")
    parser.add_argument("--backup", action="store_true", help="Create timestamped DB backup before applying.")
    return parser.parse_args()


def read_candidates(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = [c for c in REQUIRED_COLUMNS if c not in (reader.fieldnames or [])]
        if missing:
            raise SystemExit(f"Candidate TSV missing columns: {missing}")
        rows = list(reader)

    seen = set()
    duplicates = []
    for row in rows:
        key = (row["tf_id"], row["source"], row["motif_id"])
        if key in seen:
            duplicates.append(key)
        seen.add(key)

    if duplicates:
        raise SystemExit(f"Duplicate tf_id/source/motif_id keys found; first duplicate: {duplicates[0]}")

    return rows


def main() -> None:
    args = parse_args()
    db = Path(args.db)
    candidates = Path(args.candidates)

    if not db.exists():
        raise SystemExit(f"Missing DB: {db}")
    if not candidates.exists():
        raise SystemExit(f"Missing candidate TSV: {candidates}")

    rows = read_candidates(candidates)

    con = sqlite3.connect(str(db))
    cur = con.cursor()

    existing_refs = set(cur.execute("SELECT tf_id, source, motif_id FROM motif_ref").fetchall())
    valid_tfs = set(r[0] for r in cur.execute("SELECT tf_id FROM tf"))
    valid_motifs = set(
        r[0] for r in cur.execute("SELECT motif_id FROM motif_file WHERE source='cisbp'")
    )

    insert_rows = []
    skipped_existing = 0
    skipped_missing_tf = 0
    skipped_missing_motif = 0

    for row in rows:
        key = (row["tf_id"], row["source"], row["motif_id"])

        if key in existing_refs:
            skipped_existing += 1
            continue
        if row["tf_id"] not in valid_tfs:
            skipped_missing_tf += 1
            continue
        if row["motif_id"] not in valid_motifs:
            skipped_missing_motif += 1
            continue

        insert_rows.append(
            (
                row["tf_id"],
                row["evidence_type"],
                row["source"],
                row["motif_id"],
                row["original_value"],
                None if row["identity_percent"] == "" else float(row["identity_percent"]),
                int(row["missing_local_file"]),
                row["original_column"],
                row["mapping_type"],
                row["curation_status"],
                row["evidence_note"],
                int(row["display_priority"]),
            )
        )

    print(f"Candidate rows: {len(rows)}")
    print(f"Already present: {skipped_existing}")
    print(f"Skipped missing TF: {skipped_missing_tf}")
    print(f"Skipped missing CIS-BP motif_file: {skipped_missing_motif}")
    print(f"Rows to insert: {len(insert_rows)}")
    print(f"Mode: {'APPLY' if args.apply else 'DRY-RUN'}")

    if not args.apply:
        return

    if args.backup:
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        backup = db.with_name(f"{db.name}.before_cisbp_v2_direct_{stamp}.bak")
        shutil.copy2(db, backup)
        print(f"Backup written: {backup}")

    with con:
        cur.executemany(
            """
            INSERT INTO motif_ref (
                tf_id,
                evidence_type,
                source,
                motif_id,
                original_value,
                identity_percent,
                missing_local_file,
                original_column,
                mapping_type,
                curation_status,
                evidence_note,
                display_priority
            )
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            insert_rows,
        )

        cur.execute(
            """
            UPDATE tf
            SET motif_ref_count = (
                SELECT COUNT(*)
                FROM motif_ref
                WHERE motif_ref.tf_id = tf.tf_id
                  AND motif_ref.missing_local_file = 0
            )
            """
        )

    print("Applied successfully.")


if __name__ == "__main__":
    main()
