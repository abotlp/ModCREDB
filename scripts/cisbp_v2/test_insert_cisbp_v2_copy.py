from pathlib import Path
import csv
import sqlite3
from collections import Counter

db = Path("outputs/cisbp_v2/tf_webdb.cisbp_v2_test.sqlite")
srcdb = Path("data/tf_webdb.sqlite")
cand = Path("outputs/cisbp_v2/cisbp_v2_direct_human_motif_ref_insert_candidates.tsv")

if not srcdb.exists():
    raise SystemExit(f"Missing source DB: {srcdb}")
if not cand.exists():
    raise SystemExit(f"Missing candidates: {cand}")

db.write_bytes(srcdb.read_bytes())

con = sqlite3.connect(str(db))
cur = con.cursor()

rows = []
with cand.open() as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    rows.extend(reader)

keys = [(r["tf_id"], r["source"], r["motif_id"]) for r in rows]
dups = [k for k, n in Counter(keys).items() if n > 1]
print("Candidate rows:", len(rows))
print("Duplicate candidate keys:", len(dups))
if dups:
    print("First duplicate:", dups[0])
    raise SystemExit("STOP: duplicate candidate keys")

existing = set(cur.execute("SELECT tf_id, source, motif_id FROM motif_ref").fetchall())

insert_rows = []
skipped_existing = 0
for r in rows:
    key = (r["tf_id"], r["source"], r["motif_id"])
    if key in existing:
        skipped_existing += 1
        continue
    insert_rows.append((
        r["tf_id"],
        r["evidence_type"],
        r["source"],
        r["motif_id"],
        r["original_value"],
        None if r["identity_percent"] == "" else float(r["identity_percent"]),
        int(r["missing_local_file"]),
        r["original_column"],
        r["mapping_type"],
        r["curation_status"],
        r["evidence_note"],
        int(r["display_priority"]),
    ))

print("Rows to insert:", len(insert_rows))
print("Skipped existing:", skipped_existing)

with con:
    cur.executemany("""
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
    """, insert_rows)

    cur.execute("""
        UPDATE tf
        SET motif_ref_count = (
            SELECT COUNT(*)
            FROM motif_ref
            WHERE motif_ref.tf_id = tf.tf_id
              AND motif_ref.missing_local_file = 0
        )
    """)

print("\nP04637 motif_ref rows in TEST DB:")
for row in cur.execute("""
    SELECT tf_id, evidence_type, source, motif_id, original_column,
           mapping_type, curation_status, display_priority
    FROM motif_ref
    WHERE tf_id = 'P04637'
    ORDER BY display_priority, source, motif_id
"""):
    print("\t".join("" if x is None else str(x) for x in row))

print("\nP04637 source counts in TEST DB:")
for row in cur.execute("""
    SELECT source, COUNT(*)
    FROM motif_ref
    WHERE tf_id = 'P04637'
      AND missing_local_file = 0
    GROUP BY source
    ORDER BY source
"""):
    print(row[0], row[1])

print("\nP04637 stored tf.motif_ref_count:")
print(cur.execute("SELECT motif_ref_count FROM tf WHERE tf_id='P04637'").fetchone()[0])

print("\nMissing Transfac-licensed TP53 motifs should remain absent:")
for motif in ["M11199_2.00", "M11200_2.00", "M11201_2.00", "M11202_2.00", "M11203_2.00"]:
    n = cur.execute("""
        SELECT COUNT(*) FROM motif_ref
        WHERE tf_id='P04637' AND source='cisbp' AND motif_id=?
    """, (motif,)).fetchone()[0]
    print(motif, n)

print("\nTEST DB written:", db)
