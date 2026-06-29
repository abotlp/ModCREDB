#!/usr/bin/env python3
"""Build and audit external motif source links for ModCREDB.

Rules:
- JASPAR motif IDs link directly to the JASPAR matrix page.
- HOCOMOCO v11 IDs are converted to the HOCOMOCO v14 URL ID convention.
- CisBP links are only emitted when a real M-motif to T-record mapping file is supplied.
  The script does not invent CisBP T IDs from M IDs.
"""

from __future__ import annotations

import argparse
import csv
import re
import sqlite3
import urllib.error
import urllib.request
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from urllib.parse import quote

SOURCE_ORDER = {"jaspar": 0, "cisbp": 1, "hocomoco": 2}


def hocomoco_v11_to_v14(motif_id: str) -> str:
    match = re.fullmatch(r"(.+?)(?:_HUMAN)?\.H11MO\.([^.]+)\.([A-Za-z])", str(motif_id or "").strip())
    if not match:
        return str(motif_id or "").strip()
    gene, subtype, quality = match.groups()
    return f"{gene}.H14CORE.{subtype}.P.{quality.upper()}"


def load_cisbp_map(path: str | None) -> dict[str, str]:
    if not path:
        return {}
    p = Path(path)
    if not p.exists():
        raise SystemExit(f"CisBP mapping file does not exist: {p}")
    with p.open(newline="") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        dialect = csv.Sniffer().sniff(sample, delimiters="\t,") if sample.strip() else csv.excel_tab
        rows = list(csv.reader(handle, dialect))
    if not rows:
        return {}

    header = [cell.strip().lower() for cell in rows[0]]
    has_header = any(name in header for name in ["motif_id", "cisbp_motif_id", "cisbp_tf_id", "tf_id", "source_id"])
    start = 1 if has_header else 0
    if has_header:
        def col(*names: str) -> int | None:
            for name in names:
                if name in header:
                    return header.index(name)
            return None
        motif_col = col("motif_id", "cisbp_motif_id", "motif")
        tf_col = col("cisbp_tf_id", "tf_id", "source_id", "cisbp_id")
        if motif_col is None or tf_col is None:
            raise SystemExit("CisBP map header must include motif_id and cisbp_tf_id/tf_id/source_id")
    else:
        motif_col, tf_col = 0, 1

    mapping: dict[str, str] = {}
    for row in rows[start:]:
        if len(row) <= max(motif_col, tf_col):
            continue
        motif_id = row[motif_col].strip()
        tf_id = row[tf_col].strip()
        if motif_id and tf_id:
            mapping[motif_id] = tf_id
    return mapping


def read_motifs(db_path: str) -> list[tuple[str, str]]:
    conn = sqlite3.connect(db_path)
    rows = conn.execute(
        """
        SELECT DISTINCT source, motif_id
        FROM motif_file
        WHERE source IN ('jaspar', 'cisbp', 'hocomoco')
        ORDER BY source, motif_id
        """
    ).fetchall()
    conn.close()
    return [(str(source), str(motif_id)) for source, motif_id in rows]


def source_link(source: str, motif_id: str, cisbp_map: dict[str, str]) -> dict[str, str] | None:
    source = source.lower()
    if source == "jaspar":
        return {
            "source": source,
            "motif_id": motif_id,
            "mapped_id": motif_id,
            "label": "Open in JASPAR",
            "url": f"https://jaspar.elixir.no/matrix/{quote(motif_id, safe='')}/",
            "note": "JASPAR matrix page",
        }
    if source == "hocomoco":
        mapped = hocomoco_v11_to_v14(motif_id)
        return {
            "source": source,
            "motif_id": motif_id,
            "mapped_id": mapped,
            "label": "Open in HOCOMOCO",
            "url": f"https://hocomoco14.autosome.org/motif/{quote(mapped, safe='')}",
            "note": "HOCOMOCO v14 motif page",
        }
    if source == "cisbp":
        mapped = cisbp_map.get(motif_id)
        if not mapped:
            return None
        return {
            "source": source,
            "motif_id": motif_id,
            "mapped_id": mapped,
            "label": "Open in CisBP",
            "url": f"https://cisbp.ccbr.utoronto.ca/TFnewreport.php?searchTF={quote(mapped, safe='')}",
            "note": "CisBP TF report page",
        }
    return None


def build_links(db_path: str, cisbp_map: dict[str, str]) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    links: list[dict[str, str]] = []
    missing: list[dict[str, str]] = []
    for source, motif_id in read_motifs(db_path):
        link = source_link(source, motif_id, cisbp_map)
        if link:
            links.append(link)
        elif source == "cisbp":
            missing.append({
                "source": source,
                "motif_id": motif_id,
                "mapped_id": "",
                "label": "",
                "url": "",
                "note": "Missing exact CisBP T..._2.00 mapping",
            })
    return links, missing


def write_link_map(path: str, rows: list[dict[str, str]]) -> None:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    fields = ["source", "motif_id", "mapped_id", "label", "url", "note"]
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def check_url(row: dict[str, str], timeout: float) -> dict[str, str]:
    request = urllib.request.Request(row["url"], headers={"User-Agent": "ModCREDB-link-audit/1.0"})
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            code = int(getattr(response, "status", 200))
    except urllib.error.HTTPError as exc:
        return {**row, "status": "HTTP_ERROR", "http_code": str(exc.code), "error": str(exc)}
    except Exception as exc:
        return {**row, "status": "ERROR", "http_code": "", "error": repr(exc)}
    return {**row, "status": "OK" if code < 400 else "HTTP_ERROR", "http_code": str(code), "error": ""}


def audit_links(rows: list[dict[str, str]], timeout: float, workers: int) -> list[dict[str, str]]:
    results: list[dict[str, str]] = []
    with ThreadPoolExecutor(max_workers=max(1, workers)) as pool:
        futures = [pool.submit(check_url, row, timeout) for row in rows]
        for future in as_completed(futures):
            results.append(future.result())
    return sorted(results, key=lambda row: (SOURCE_ORDER.get(row["source"], 99), row["motif_id"]))


def write_audit(path: str, rows: list[dict[str, str]]) -> None:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    fields = ["source", "motif_id", "mapped_id", "label", "url", "note", "status", "http_code", "error"]
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def print_counts(title: str, rows: list[dict[str, str]], status_field: str | None = None) -> None:
    counts: dict[str, Counter[str]] = defaultdict(Counter)
    for row in rows:
        key = row[status_field] if status_field else "COUNT"
        counts[row["source"]][key] += 1
    print(title)
    for source in sorted(counts, key=lambda s: SOURCE_ORDER.get(s, 99)):
        total = sum(counts[source].values())
        detail = ", ".join(f"{k}={v}" for k, v in sorted(counts[source].items()))
        print(f"{source}\ttotal={total}\t{detail}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--db", default="data/tf_webdb.sqlite")
    parser.add_argument("--cisbp-map", default="")
    parser.add_argument("--write-map", default="data/external_motif_links.tsv")
    parser.add_argument("--write-missing-cisbp", default="outputs/cisbp_missing_external_ids.tsv")
    parser.add_argument("--out", default="outputs/external_motif_link_audit.tsv")
    parser.add_argument("--timeout", type=float, default=8.0)
    parser.add_argument("--workers", type=int, default=12)
    parser.add_argument("--no-http", action="store_true")
    args = parser.parse_args()

    cisbp_map = load_cisbp_map(args.cisbp_map or None)
    links, missing = build_links(args.db, cisbp_map)
    write_link_map(args.write_map, links)
    write_link_map(args.write_missing_cisbp, missing)
    print_counts("Generated external links", links)
    print(f"Missing CisBP exact mappings: {len(missing)}")

    if not args.no_http:
        audited = audit_links(links, args.timeout, args.workers)
        write_audit(args.out, audited)
        print_counts("HTTP audit", audited, "status")
        bad = [row for row in audited if row["status"] != "OK"][:25]
        if bad:
            print("First non-OK examples:")
            for row in bad:
                print(f"{row['source']}\t{row['motif_id']}\t{row['mapped_id']}\t{row['status']}\t{row['url']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
