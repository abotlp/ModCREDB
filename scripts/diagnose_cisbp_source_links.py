#!/usr/bin/env python3
from pathlib import Path
import csv
import json
import argparse
from collections import defaultdict, Counter

LINKS = Path("data/external_motif_links.tsv")

DEFAULT_ROOTS = [
    Path("/tmp/cisbp_hs_check"),
    Path("/home/patricia/tf_webdb/external/cisbp"),
]

SENTINELS = {
    "M11192_2.00": "expect direct non-human source T305275_2.00 if all-species metadata is available",
    "M09337_2.00": "TP53 direct human source should remain T311040_2.00",
    "M08359_2.00": "known working CIS-BP example",
}


def read_tsv(path):
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def pct(n, d):
    return 0.0 if d == 0 else round(100.0 * n / d, 2)


def find_metadata_files(roots):
    out = []
    seen = set()
    names = {
        "TF_Information.txt",
        "TF_Information_all_motifs.txt",
        "TF_Information_all_motifs_plus.txt",
    }
    for root in roots:
        if not root.exists():
            continue
        try:
            iterator = root.rglob("*")
            for p in iterator:
                try:
                    if p.is_file() and p.name in names:
                        rp = p.resolve()
                        if rp not in seen:
                            out.append(p)
                            seen.add(rp)
                except OSError as e:
                    print(f"WARNING: skipping unreadable path {p}: {e}")
                    continue
        except OSError as e:
            print(f"WARNING: skipping unreadable root {root}: {e}")
            continue
    return sorted(out)


def classify(species, status):
    if status == "D" and species == "Homo_sapiens":
        return "direct_human", 0
    if status == "D" and species != "Homo_sapiens":
        return "direct_nonhuman", 1
    if status != "D" and species == "Homo_sapiens":
        return "related_or_inferred_human", 2
    return "related_or_inferred_nonhuman", 3


def parse_metadata(meta_files, target_motifs):
    hits_by_motif = defaultdict(list)

    for p in meta_files:
        with p.open(errors="replace") as fh:
            header = fh.readline().rstrip("\n").split("\t")
            for line_no, line in enumerate(fh, start=2):
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 14:
                    continue

                tf_id = parts[0]
                motif_id = parts[3]
                if motif_id not in target_motifs:
                    continue

                tf_external_id = parts[5] if len(parts) > 5 else ""
                tf_name = parts[6] if len(parts) > 6 else ""
                species = parts[7] if len(parts) > 7 else ""
                status = parts[8] if len(parts) > 8 else ""
                source_id = parts[13] if len(parts) > 13 else ""
                pmid = parts[19] if len(parts) > 19 else ""

                link_class, priority = classify(species, status)

                hits_by_motif[motif_id].append({
                    "motif_id": motif_id,
                    "cisbp_tf_id": tf_id,
                    "tf_name": tf_name,
                    "tf_external_id": tf_external_id,
                    "tf_species": species,
                    "tf_status": status,
                    "motif_source_id": source_id,
                    "pmid": pmid,
                    "metadata_file": str(p),
                    "metadata_line": str(line_no),
                    "link_class": link_class,
                    "priority": str(priority),
                    "url": f"https://cisbp.ccbr.utoronto.ca/TFnewreport.php?searchTF={tf_id}",
                })

    return hits_by_motif


def choose_best(hits):
    return sorted(
        hits,
        key=lambda r: (
            int(r["priority"]),
            r["tf_species"] != "Homo_sapiens",
            r["tf_name"].lower(),
            r["cisbp_tf_id"],
            r["metadata_file"],
            int(r["metadata_line"]),
        )
    )[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", action="append", default=[], help="Additional root directory containing CIS-BP TF_Information files")
    ap.add_argument("--outdir", default="data/cisbp_v2/source_link_diagnostic")
    args = ap.parse_args()

    roots = DEFAULT_ROOTS + [Path(x) for x in args.root]
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if not LINKS.exists():
        raise SystemExit(f"Missing {LINKS}")

    links = read_tsv(LINKS)
    cisbp_links = [r for r in links if r.get("source") == "cisbp"]
    target_motifs = {r["motif_id"] for r in cisbp_links}

    meta_files = find_metadata_files(roots)
    print("Metadata files found:")
    for p in meta_files:
        print(" ", p)
    print()

    hits_by_motif = parse_metadata(meta_files, target_motifs)

    all_hits = []
    proposed = []
    compare = []

    current_by_motif = {r["motif_id"]: r for r in cisbp_links}

    for motif_id in sorted(target_motifs):
        hits = hits_by_motif.get(motif_id, [])
        current = current_by_motif[motif_id]
        current_url = current.get("url", "")
        current_tf = ""
        if "TFnewreport.php?searchTF=" in current_url:
            current_tf = current_url.rsplit("searchTF=", 1)[-1]

        if hits:
            best = choose_best(hits)
            proposed.append(best)
            proposed_tf = best["cisbp_tf_id"]
            proposed_url = best["url"]
            proposed_class = best["link_class"]
        else:
            best = None
            proposed_tf = ""
            proposed_url = ""
            proposed_class = "no_metadata"

        for h in hits:
            h2 = dict(h)
            h2["n_hits_for_motif"] = str(len(hits))
            all_hits.append(h2)

        current_is_specific = "TFnewreport.php?searchTF=" in current_url
        proposed_is_specific = bool(proposed_url)

        change_type = "unchanged"
        if not proposed_is_specific and current_is_specific:
            change_type = "current_specific_but_no_metadata_found"
        elif proposed_is_specific and not current_is_specific:
            change_type = "fallback_to_specific"
        elif proposed_is_specific and current_is_specific and proposed_tf != current_tf:
            change_type = "specific_to_better_specific"
        elif not proposed_is_specific and not current_is_specific:
            change_type = "fallback_stays_fallback"

        compare.append({
            "motif_id": motif_id,
            "current_tf": current_tf,
            "current_url": current_url,
            "proposed_tf": proposed_tf,
            "proposed_url": proposed_url,
            "proposed_class": proposed_class,
            "change_type": change_type,
            "n_metadata_hits": str(len(hits)),
        })

    def write_tsv(path, rows, fields):
        with path.open("w", newline="") as fh:
            w = csv.DictWriter(fh, delimiter="\t", fieldnames=fields)
            w.writeheader()
            w.writerows(rows)

    all_fields = [
        "motif_id", "cisbp_tf_id", "tf_name", "tf_external_id", "tf_species",
        "tf_status", "motif_source_id", "pmid", "metadata_file", "metadata_line",
        "link_class", "priority", "n_hits_for_motif", "url",
    ]
    prop_fields = [
        "motif_id", "cisbp_tf_id", "tf_name", "tf_external_id", "tf_species",
        "tf_status", "motif_source_id", "pmid", "metadata_file", "metadata_line",
        "link_class", "priority", "url",
    ]
    cmp_fields = [
        "motif_id", "current_tf", "current_url", "proposed_tf",
        "proposed_url", "proposed_class", "change_type", "n_metadata_hits",
    ]

    write_tsv(outdir / "cisbp_source_link_all_hits.tsv", all_hits, all_fields)
    write_tsv(outdir / "cisbp_source_link_proposed_best.tsv", proposed, prop_fields)
    write_tsv(outdir / "cisbp_source_link_compare_current.tsv", compare, cmp_fields)

    counts = Counter(r["change_type"] for r in compare)
    class_counts = Counter(r["link_class"] for r in proposed)
    current_specific = sum(1 for r in cisbp_links if "TFnewreport.php?searchTF=" in r.get("url", ""))
    current_fallback = len(cisbp_links) - current_specific

    summary = {
        "cisbp_total": len(cisbp_links),
        "current_specific": current_specific,
        "current_fallback": current_fallback,
        "metadata_files": [str(p) for p in meta_files],
        "proposed_specific": len(proposed),
        "proposed_fallback": len(cisbp_links) - len(proposed),
        "change_type_counts": dict(counts),
        "proposed_class_counts": dict(class_counts),
    }

    with (outdir / "summary.json").open("w") as fh:
        json.dump(summary, fh, indent=2, sort_keys=True)

    print("=== Summary ===")
    print(json.dumps(summary, indent=2, sort_keys=True))

    print()
    print("=== Sentinel motifs ===")
    cmp_by_motif = {r["motif_id"]: r for r in compare}
    hits_by_motif_flat = defaultdict(list)
    for h in all_hits:
        hits_by_motif_flat[h["motif_id"]].append(h)

    for motif, expectation in SENTINELS.items():
        print()
        print(motif, "-", expectation)
        print("COMPARE:", cmp_by_motif.get(motif))
        for h in hits_by_motif_flat.get(motif, [])[:20]:
            print(
                "HIT",
                h["cisbp_tf_id"],
                h["tf_name"],
                h["tf_species"],
                "status=" + h["tf_status"],
                h["link_class"],
                h["metadata_file"],
                h["metadata_line"],
                sep="\t",
            )

    print()
    print("Wrote:")
    print(" ", outdir / "cisbp_source_link_all_hits.tsv")
    print(" ", outdir / "cisbp_source_link_proposed_best.tsv")
    print(" ", outdir / "cisbp_source_link_compare_current.tsv")
    print(" ", outdir / "summary.json")


if __name__ == "__main__":
    main()
