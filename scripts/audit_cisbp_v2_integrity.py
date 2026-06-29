#!/usr/bin/env python3
from pathlib import Path
import csv
import sqlite3
import sys

DB = Path("data/tf_webdb.sqlite")
LINKS = Path("data/external_motif_links.tsv")
MAP = Path("data/cisbp_v2/cisbp_v2_direct_motif_to_tf_report.tsv")

EXPECTED_TP53 = {
    "M09337_2.00": "T311040_2.00",
    "M09621_2.00": "T311040_2.00",
    "M11197_2.00": "T311040_2.00",
    "M11198_2.00": "T311040_2.00",
}

errors = []

def fail(msg):
    errors.append(msg)

for p in [DB, LINKS, MAP]:
    if not p.exists():
        fail(f"Missing required file: {p}")

if errors:
    print("\n".join(errors))
    sys.exit(1)

con = sqlite3.connect(DB)
cur = con.cursor()

# 1. Raw PWM counts.
motif_file_counts = dict(cur.execute("""
    SELECT source, COUNT(*)
    FROM motif_file
    GROUP BY source
""").fetchall())

print("motif_file counts:", motif_file_counts)

if motif_file_counts.get("cisbp") != 10035:
    fail(f"Unexpected CIS-BP motif_file count: {motif_file_counts.get('cisbp')}")

# 2. Verified direct CIS-BP assignment count.
verified_direct = cur.execute("""
    SELECT COUNT(*)
    FROM motif_ref
    WHERE source='cisbp'
      AND mapping_type='direct_cisbp_v2_tf_information'
      AND curation_status='verified_same_version_source_mapping'
""").fetchone()[0]

verified_direct_tfs = cur.execute("""
    SELECT COUNT(DISTINCT tf_id)
    FROM motif_ref
    WHERE source='cisbp'
      AND mapping_type='direct_cisbp_v2_tf_information'
      AND curation_status='verified_same_version_source_mapping'
""").fetchone()[0]

print("verified direct CIS-BP v2 motif_ref rows:", verified_direct)
print("TFs with verified direct CIS-BP v2 motif_ref rows:", verified_direct_tfs)

if verified_direct != 3943:
    fail(f"Unexpected verified direct CIS-BP motif_ref count: {verified_direct}")

# 3. Every verified direct motif_ref has local PWM.
missing_pwm = cur.execute("""
    SELECT mr.tf_id, mr.motif_id
    FROM motif_ref mr
    LEFT JOIN motif_file mf
      ON mf.source='cisbp'
     AND mf.motif_id=mr.motif_id
    WHERE mr.source='cisbp'
      AND mr.mapping_type='direct_cisbp_v2_tf_information'
      AND mr.curation_status='verified_same_version_source_mapping'
      AND mf.motif_id IS NULL
    LIMIT 20
""").fetchall()

if missing_pwm:
    fail(f"Verified CIS-BP motif_ref rows missing local PWM: {missing_pwm[:5]}")

# 4. TP53 control.
tp53_rows = dict(cur.execute("""
    SELECT motif_id, mapping_type
    FROM motif_ref
    WHERE tf_id='P04637'
      AND source='cisbp'
""").fetchall())

print("P04637 CIS-BP rows:", tp53_rows)

if set(tp53_rows) != set(EXPECTED_TP53):
    fail(f"P04637 CIS-BP motifs wrong: observed={sorted(tp53_rows)} expected={sorted(EXPECTED_TP53)}")

tp53_count = cur.execute("SELECT motif_ref_count FROM tf WHERE tf_id='P04637'").fetchone()[0]
print("P04637 motif_ref_count:", tp53_count)

if tp53_count != 8:
    fail(f"P04637 motif_ref_count should be 8, observed {tp53_count}")

# 5. Load link map.
with MAP.open(newline="") as fh:
    map_rows = list(csv.DictReader(fh, delimiter="\t"))

map_by_motif = {r["motif_id"]: r for r in map_rows}
print("CIS-BP v2 direct motif-to-report map rows:", len(map_rows))

if len(map_rows) != len(map_by_motif):
    fail("Duplicate motif_id values in CIS-BP motif-to-report map")

# 6. Load external links.
with LINKS.open(newline="") as fh:
    link_rows = list(csv.DictReader(fh, delimiter="\t"))

cisbp_links = [r for r in link_rows if r["source"] == "cisbp"]
links_by_motif = {r["motif_id"]: r for r in cisbp_links}

print("external_motif_links CIS-BP rows:", len(cisbp_links))

if len(cisbp_links) != 10035:
    fail(f"Unexpected CIS-BP external link count: {len(cisbp_links)}")

tf_report_links = 0
fallback_links = 0
bad_v3 = []
bad_map_links = []
missing_link_rows = []

for motif_id, m in map_by_motif.items():
    link = links_by_motif.get(motif_id)
    if link is None:
        missing_link_rows.append(motif_id)
        continue

    expected_url = m["url"]
    observed_url = link["url"]

    if observed_url != expected_url:
        bad_map_links.append((motif_id, observed_url, expected_url))

for r in cisbp_links:
    url = r.get("url", "")
    if "TFnewreport.php?searchTF=" in url:
        tf_report_links += 1
        if "_2.00" not in url or "_3.10" in url:
            bad_v3.append((r["motif_id"], url))
    else:
        fallback_links += 1

print("CIS-BP TF-report links:", tf_report_links)
print("CIS-BP fallback links:", fallback_links)

if missing_link_rows:
    fail(f"Motifs in v2 map missing from external_motif_links: {missing_link_rows[:5]}")

if bad_map_links:
    fail(f"Mapped CIS-BP links do not match v2 report map: {bad_map_links[:5]}")

if bad_v3:
    fail(f"Bad CIS-BP v3/non-v2 links found: {bad_v3[:5]}")

if tf_report_links != len(map_by_motif):
    fail(f"TF-report link count {tf_report_links} does not equal map rows {len(map_by_motif)}")

# 7. TP53 link control.
for motif_id, tf_id in EXPECTED_TP53.items():
    link = links_by_motif.get(motif_id)
    expected_url = f"https://cisbp.ccbr.utoronto.ca/TFnewreport.php?searchTF={tf_id}"
    if not link:
        fail(f"Missing external link row for {motif_id}")
    elif link["url"] != expected_url:
        fail(f"Wrong TP53 link for {motif_id}: {link['url']} expected {expected_url}")

if errors:
    print("\nFAIL")
    for e in errors:
        print(" -", e)
    sys.exit(1)

print("\nPASS: CIS-BP v2 assignments and links are internally consistent.")
