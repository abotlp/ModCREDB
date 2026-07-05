#!/usr/bin/env python3
from pathlib import Path
import csv
import sqlite3

DB = Path("data/tf_webdb.sqlite")
LINKS = Path("data/external_motif_links.tsv")
DIRECT_MAP = Path("data/cisbp_v2/cisbp_v2_direct_motif_to_tf_report.tsv")

def read_tsv(path):
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

def get_cisbp_motifs_from_db():
    if not DB.exists():
        raise SystemExit(f"Missing SQLite DB: {DB}")

    con = sqlite3.connect(DB)
    cur = con.cursor()

    tables = {r[0] for r in cur.execute("SELECT name FROM sqlite_master WHERE type='table'")}
    if "motif_file" not in tables:
        raise SystemExit("SQLite DB does not contain table motif_file")

    cols = [r[1] for r in cur.execute("PRAGMA table_info(motif_file)")]

    if "source" not in cols or "motif_id" not in cols:
        raise SystemExit(f"motif_file table columns are unexpected: {cols}")

    motifs = {
        r[0]
        for r in cur.execute(
            "SELECT DISTINCT motif_id FROM motif_file WHERE source = ?",
            ("cisbp",),
        )
    }
    con.close()
    return motifs

link_rows = read_tsv(LINKS)
direct_rows = read_tsv(DIRECT_MAP)

cisbp_motifs = get_cisbp_motifs_from_db()

cisbp_links = [
    r for r in link_rows
    if r.get("source") == "cisbp"
]

direct_by_motif = {}
for r in direct_rows:
    direct_by_motif.setdefault(r["motif_id"], r)

link_by_motif = {}
dups = []
for r in cisbp_links:
    k = r["motif_id"]
    if k in link_by_motif:
        dups.append(k)
    link_by_motif[k] = r

specific = []
fallback = []
bad_specific = []
bad_fallback = []
bad_url = []

for r in cisbp_links:
    motif = r["motif_id"]
    url = r.get("url", "")
    mapped = r.get("mapped_id", "")

    is_specific = "TFnewreport.php?searchTF=" in url

    if is_specific:
        specific.append(r)
        expected = direct_by_motif.get(motif)
        if not expected:
            bad_specific.append((motif, mapped, "specific but motif not in direct-human map"))
        else:
            expected_tf = expected["cisbp_tf_id"]
            if mapped != expected_tf or not url.endswith(expected_tf):
                bad_specific.append((motif, mapped, f"expected {expected_tf}"))
    else:
        fallback.append(r)
        if motif in direct_by_motif:
            bad_fallback.append((motif, mapped, direct_by_motif[motif]["cisbp_tf_id"]))

    if "_3." in url or "_3_" in url or "3_00" in url:
        bad_url.append((motif, url))

missing_links = sorted(cisbp_motifs - set(link_by_motif))
extra_links = sorted(set(link_by_motif) - cisbp_motifs)

sentinels = {}
for motif in ["M11192_2.00", "M09337_2.00", "M08359_2.00"]:
    sentinels[motif] = link_by_motif.get(motif)

print("=== Direct-only CIS-BP external link audit ===")
print("CIS-BP motif_file DB rows:        ", len(cisbp_motifs))
print("CIS-BP external link rows:        ", len(cisbp_links))
print("Specific direct-human links:      ", len(specific))
print("Fallback links:                   ", len(fallback))
print("Direct-human map rows:            ", len(direct_by_motif))
print()
print("Duplicate CIS-BP link keys:       ", len(dups))
print("motif_file rows missing links:    ", len(missing_links))
print("link rows without motif_file rows:", len(extra_links))
print("Bad specific direct links:        ", len(bad_specific))
print("Bad fallback rows:                ", len(bad_fallback))
print("Bad v3/version URLs:              ", len(bad_url))

print("\nSentinels:")
for motif, r in sentinels.items():
    if not r:
        print(motif, "MISSING")
        continue
    print(motif)
    print("  mapped_id:", r.get("mapped_id"))
    print("  url:", r.get("url"))
    print("  note:", r.get("note"))

failed = False
checks = [
    ("duplicate keys", dups),
    ("missing links", missing_links),
    ("extra links", extra_links),
    ("bad specific links", bad_specific),
    ("bad fallback rows", bad_fallback),
    ("bad version URLs", bad_url),
]
for name, vals in checks:
    if vals:
        failed = True
        print(f"\nFAIL examples for {name}:")
        for x in vals[:20]:
            print(x)

expect = {
    "M11192_2.00": ("M11192_2.00", "https://cisbp.ccbr.utoronto.ca/"),
    "M09337_2.00": ("T311040_2.00", "https://cisbp.ccbr.utoronto.ca/TFnewreport.php?searchTF=T311040_2.00"),
    "M08359_2.00": ("T095207_2.00", "https://cisbp.ccbr.utoronto.ca/TFnewreport.php?searchTF=T095207_2.00"),
}
for motif, (mapped, url) in expect.items():
    r = sentinels.get(motif)
    if not r or r.get("mapped_id") != mapped or r.get("url") != url:
        failed = True
        print("\nFAIL sentinel:", motif, "expected", mapped, url, "observed", r)

if failed:
    raise SystemExit("FAIL: direct-only audit did not pass.")

print("\nPASS: direct-only CIS-BP links are internally consistent.")
