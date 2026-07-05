#!/usr/bin/env python3
from pathlib import Path
import csv
from collections import Counter

LINKS = Path("data/external_motif_links.tsv")
DIRECT_MAP = Path("data/cisbp_v2/cisbp_v2_direct_motif_to_tf_report.tsv")
OUT = Path("data/cisbp_v2/cisbp_external_links_direct_only_PREVIEW.tsv")
CMP = Path("data/cisbp_v2/cisbp_external_links_direct_only_COMPARE.tsv")

FALLBACK_URL = "https://cisbp.ccbr.utoronto.ca/"
FALLBACK_NOTE = "CisBP v2 motif identifier; search this motif ID in CisBP using VM/version 2.00."

def read_tsv(path):
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

def write_tsv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

links = read_tsv(LINKS)
direct_rows = read_tsv(DIRECT_MAP)

direct_by_motif = {}
for r in direct_rows:
    motif = r["motif_id"]
    # direct map should already be direct human only.
    direct_by_motif.setdefault(motif, r)

out = []
compare = []

for r in links:
    r2 = dict(r)

    if r.get("source") == "cisbp":
        motif = r["motif_id"]
        old_url = r.get("url", "")
        old_mapped = r.get("mapped_id", "")
        old_note = r.get("note", "")

        if motif in direct_by_motif:
            d = direct_by_motif[motif]
            tfid = d["cisbp_tf_id"]
            r2["mapped_id"] = tfid
            r2["label"] = "CIS-BP v2 direct source"
            r2["url"] = f"https://cisbp.ccbr.utoronto.ca/TFnewreport.php?searchTF={tfid}"
            r2["note"] = (
                f"CIS-BP v2 direct human source mapping from {d.get('metadata_file','CIS-BP metadata')}: "
                f"{motif} -> {tfid} ({d.get('tf_name','')}, {d.get('tf_species','')}, TF_Status={d.get('tf_status','D')})."
            )
            new_class = "direct_human"
        else:
            r2["mapped_id"] = motif
            r2["label"] = "CIS-BP v2"
            r2["url"] = FALLBACK_URL
            r2["note"] = FALLBACK_NOTE
            new_class = "fallback_no_direct_source"

        compare.append({
            "motif_id": motif,
            "old_mapped_id": old_mapped,
            "old_url": old_url,
            "old_note": old_note,
            "new_mapped_id": r2["mapped_id"],
            "new_url": r2["url"],
            "new_note": r2["note"],
            "new_class": new_class,
            "changed": str((old_url != r2["url"]) or (old_mapped != r2["mapped_id"])),
        })

    out.append(r2)

fields = list(links[0].keys())
write_tsv(OUT, out, fields)
write_tsv(CMP, compare, [
    "motif_id", "old_mapped_id", "old_url", "old_note",
    "new_mapped_id", "new_url", "new_note", "new_class", "changed"
])

counts = Counter(r["new_class"] for r in compare)
changed = sum(1 for r in compare if r["changed"] == "True")

print("Wrote preview:", OUT)
print("Wrote compare:", CMP)
print("Counts:")
for k, v in sorted(counts.items()):
    print(k, v)
print("Changed CIS-BP rows:", changed)

print("\nSentinels:")
for motif in ["M11192_2.00", "M09337_2.00", "M08359_2.00"]:
    rows = [r for r in compare if r["motif_id"] == motif]
    for r in rows:
        print(motif)
        print("  old:", r["old_mapped_id"], r["old_url"])
        print("  new:", r["new_mapped_id"], r["new_url"])
        print("  class:", r["new_class"])
