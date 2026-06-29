#!/usr/bin/env python3
from pathlib import Path
import csv
import sqlite3
import sys
from collections import Counter, defaultdict

DB = Path("data/tf_webdb.sqlite")
LINKS = Path("data/external_motif_links.tsv")
CISBP_MAP = Path("data/cisbp_v2/cisbp_v2_direct_motif_to_tf_report.tsv")

errors = []

def fail(msg):
    errors.append(msg)

for p in [DB, LINKS, CISBP_MAP]:
    if not p.exists():
        fail(f"Missing required file: {p}")

if errors:
    print("\n".join(errors))
    sys.exit(1)

con = sqlite3.connect(DB)
cur = con.cursor()

# Load external links.
with LINKS.open(newline="") as fh:
    links = list(csv.DictReader(fh, delimiter="\t"))

required_link_cols = {"source", "motif_id", "url"}
missing = required_link_cols - set(links[0].keys())
if missing:
    fail(f"external_motif_links.tsv missing columns: {missing}")

link_key_counts = Counter((r["source"], r["motif_id"]) for r in links)
dups = [k for k, n in link_key_counts.items() if n > 1]
if dups:
    fail(f"Duplicate external link keys found, first examples: {dups[:10]}")

links_by_key = {(r["source"], r["motif_id"]): r for r in links}

# motif_file coverage.
motif_file_rows = list(cur.execute("SELECT source, motif_id FROM motif_file"))
motif_file_keys = set(motif_file_rows)
motif_file_counts = Counter(source for source, motif in motif_file_rows)

print("motif_file counts:")
for k, v in sorted(motif_file_counts.items()):
    print(f"  {k}: {v}")

missing_external = sorted(k for k in motif_file_keys if k not in links_by_key)
extra_external = sorted(k for k in links_by_key if k not in motif_file_keys)

if missing_external:
    fail(f"motif_file rows missing external links: {missing_external[:20]} total={len(missing_external)}")

if extra_external:
    fail(f"external links without motif_file rows: {extra_external[:20]} total={len(extra_external)}")

# motif_ref must point to existing motif_file unless marked missing.
bad_refs = list(cur.execute("""
    SELECT tf_id, source, motif_id
    FROM motif_ref
    WHERE missing_local_file = 0
      AND NOT EXISTS (
        SELECT 1
        FROM motif_file mf
        WHERE mf.source = motif_ref.source
          AND mf.motif_id = motif_ref.motif_id
      )
    LIMIT 50
"""))

if bad_refs:
    fail(f"motif_ref rows point to missing motif_file rows: {bad_refs[:10]}")

# CisBP v2 map integrity.
with CISBP_MAP.open(newline="") as fh:
    cisbp_map_rows = list(csv.DictReader(fh, delimiter="\t"))

cisbp_map_by_motif = {r["motif_id"]: r for r in cisbp_map_rows}

if len(cisbp_map_by_motif) != len(cisbp_map_rows):
    fail("Duplicate motif IDs in cisbp_v2_direct_motif_to_tf_report.tsv")

# Every mapped CIS-BP motif must have exact v2 TF-report URL.
bad_cisbp_mapped_links = []
for motif_id, m in cisbp_map_by_motif.items():
    link = links_by_key.get(("cisbp", motif_id))
    expected_url = m["url"]
    if link is None:
        bad_cisbp_mapped_links.append((motif_id, "missing external link", expected_url))
    elif link["url"] != expected_url:
        bad_cisbp_mapped_links.append((motif_id, link["url"], expected_url))

if bad_cisbp_mapped_links:
    fail(f"CIS-BP mapped motifs with wrong links: {bad_cisbp_mapped_links[:10]} total={len(bad_cisbp_mapped_links)}")

# No CIS-BP URL may point to v3 or versionless guessed TF report.
cisbp_links = [r for r in links if r["source"] == "cisbp"]
cisbp_report = []
cisbp_fallback = []
bad_cisbp_version = []

for r in cisbp_links:
    url = r["url"]
    if "TFnewreport.php?searchTF=" in url:
        cisbp_report.append(r)
        if "_2.00" not in url or "_3.10" in url:
            bad_cisbp_version.append((r["motif_id"], url))
    else:
        cisbp_fallback.append(r)

if bad_cisbp_version:
    fail(f"Bad CIS-BP non-v2/v3 TF report links: {bad_cisbp_version[:10]} total={len(bad_cisbp_version)}")

if len(cisbp_report) != len(cisbp_map_by_motif):
    fail(
        f"CIS-BP TF-report link count does not match map: "
        f"report_links={len(cisbp_report)} map={len(cisbp_map_by_motif)}"
    )

# Every verified direct CIS-BP motif_ref must have a v2 TF-report external link.
verified_direct_refs = list(cur.execute("""
    SELECT tf_id, motif_id
    FROM motif_ref
    WHERE source='cisbp'
      AND mapping_type='direct_cisbp_v2_tf_information'
      AND curation_status='verified_same_version_source_mapping'
      AND missing_local_file=0
"""))

bad_verified_ref_links = []
for tf_id, motif_id in verified_direct_refs:
    link = links_by_key.get(("cisbp", motif_id))
    if link is None or "TFnewreport.php?searchTF=" not in link["url"] or "_2.00" not in link["url"]:
        bad_verified_ref_links.append((tf_id, motif_id, None if link is None else link["url"]))

if bad_verified_ref_links:
    fail(f"Verified CIS-BP motif_refs lacking v2 report links: {bad_verified_ref_links[:10]} total={len(bad_verified_ref_links)}")

# Known sentinels.
sentinels = {
    ("cisbp", "M09337_2.00"): "https://cisbp.ccbr.utoronto.ca/TFnewreport.php?searchTF=T311040_2.00",
    ("cisbp", "M09621_2.00"): "https://cisbp.ccbr.utoronto.ca/TFnewreport.php?searchTF=T311040_2.00",
    ("cisbp", "M11197_2.00"): "https://cisbp.ccbr.utoronto.ca/TFnewreport.php?searchTF=T311040_2.00",
    ("cisbp", "M11198_2.00"): "https://cisbp.ccbr.utoronto.ca/TFnewreport.php?searchTF=T311040_2.00",
    ("jaspar", "MA0106.1"): "https://jaspar.elixir.no/matrix/MA0106.1/",
    ("hocomoco", "P53_HUMAN.H11MO.0.A"): "https://hocomoco14.autosome.org/motif/P53_HUMAN.H14MO.0.A",
}

for key, expected in sentinels.items():
    link = links_by_key.get(key)
    if link is None:
        fail(f"Missing sentinel link: {key}")
    elif link["url"] != expected:
        fail(f"Wrong sentinel link for {key}: observed={link['url']} expected={expected}")

# HOCOMOCO unresolved must not be guessed.
bad_hocomoco = [
    r for r in links
    if r["source"] == "hocomoco"
    and "hocomoco14.autosome.org/motif/" in r["url"]
    and ".H14MO." not in r["url"]
]
if bad_hocomoco:
    fail(f"Suspicious HOCOMOCO links without H14MO target: {bad_hocomoco[:10]}")

print()
print("external_motif_links counts:")
for source, count in sorted(Counter(r["source"] for r in links).items()):
    print(f"  {source}: {count}")

print()
print("CIS-BP link split:")
print(f"  TF-report v2 links: {len(cisbp_report)}")
print(f"  fallback links: {len(cisbp_fallback)}")
print(f"  v2 map rows: {len(cisbp_map_by_motif)}")
print(f"  verified direct motif_ref rows: {len(verified_direct_refs)}")

print()
print("motif_ref rows by source:")
for source, count in cur.execute("""
    SELECT source, COUNT(*)
    FROM motif_ref
    WHERE missing_local_file=0
    GROUP BY source
    ORDER BY source
"""):
    print(f"  {source}: {count}")

if errors:
    print("\nFAIL")
    for e in errors:
        print(" -", e)
    sys.exit(1)

print("\nPASS: external motif links and motif_ref/motif_file joins are internally consistent.")
