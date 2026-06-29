#!/usr/bin/env python3
"""Audit external motif source links for ModCREDB.

This checks distinct source/motif IDs in the SQLite database and reports how many
external links can be generated and how many respond over HTTP.

CisBP is special: ModCREDB stores M... motif IDs, while the CisBP report URL often
expects T... TF IDs. Pass --cisbp-map with a two-column TSV/CSV when available:

    motif_id    cisbp_tf_id
    M08889_2.00 T094907_2.00
"""

from __future__ import annotations

import argparse
import csv
import re
import sqlite3
import sys
import urllib.error
import urllib.request
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from urllib.parse import quote


SOURCE_ORDER = {"jaspar": 0, "cisbp": 1, "hocomoco": 2}


def hocomoco_v11_to_v14(motif_id: str) -> str:
    """Map HOCOMOCO v11-style IDs to the current v14 URL naming pattern.

    Example: FEZF1_HUMAN.H11MO.0.C -> FEZF1.H14CORE.0.P.C
    """
    match = re.match(r"^(.+?)(?:_HUMAN)?\.H11MO\.([^.]+)\.([A-Za-z])$", motif_id)
    if not match:
        return motif_id
    gene, subtype, quality = match.groups()
    return f"{gene}.H14CORE.{subtype}.P.{quality}"


def load_cisbp_map(path: str | None) -> dict[str, str]:
    if not path:
        return {}
    p = Path(path)
    if not p.exists():
        raise SystemExit(f"CisBP mapping file does not exist: {p}")
    mapping: dict[str, str] = {}
    with p.open(newline="") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        dialect = csv.Sniffer().sniff(sample, delimiters="\t,") if sample.strip() else csv.excel_tab
        rows = list(csv.reader(handle, dialect))
    if not rows:
        return mapping

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
            raise SystemExit("CisBP mapping header must include motif_id and cisbp_tf_id/tf_id/source_id columns")
    else:
        motif_col, tf_col = 0, 1

    for row in rows[start:]:
        if len(row) <= max(motif_col, tf_col):
            continue
        motif_id = row[motif_col].strip()
        tf_id = row[tf_col].strip()
        if motif_id and tf_id:
            mapping[motif_id] = tf_id
    return mapping


def source_url(source: str, motif_id: str, cisbp_map: dict[str, str]) -> tuple[str, str, str]:
    source = source.lower()
    if source == "jaspar":
        return "jaspar", motif_id, f"https://jaspar.elixir.no/matrix/{quote(motif_id, safe='')}/"
    if source == "hocomoco":
        mapped = hocomoco_v11_to_v14(motif_id)
        return "hocomoco", mapped, f"https://hocomoco14.autosome.org/motif/{quote(mapped, safe='')}"
    if source == "cisbp":
        mapped = cisbp_map.get(motif_id, motif_id)
        return "cisbp", mapped, f"https://cisbp.ccbr.utoronto.ca/TFnewreport.php?searchTF={quote(mapped, safe='')}"
    raise ValueError(source)


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


def check_url(item: tuple[str, str, str, str], timeout: float) -> dict[str, object]:
    source, motif_id, mapped_id, url = item
    if source == "cisbp" and motif_id.startswith("M") and mapped_id == motif_id:
        return {
            "source": source,
            "motif_id": motif_id,
            "mapped_id": mapped_id,
            "url": url,
            "status": "NEEDS_CISBP_T_ID_MAP",
            "http_code": "",
            "error": "CisBP report URL usually expects T..._2.00, not M..._2.00",
        }
    request = urllib.request.Request(url, headers={"User-Agent": "ModCREDB-link-audit/1.0"})
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            code = int(getattr(response, "status", 200))
            body = response.read(8192).decode("utf-8", "ignore").lower()
    except urllib.error.HTTPError as exc:
        return {"source": source, "motif_id": motif_id, "mapped_id": mapped_id, "url": url, "status": "HTTP_ERROR", "http_code": exc.code, "error": str(exc)}
    except Exception as exc:  # noqa: BLE001
        return {"source": source, "motif_id": motif_id, "mapped_id": mapped_id, "url": url, "status": "ERROR", "http_code": "", "error": repr(exc)}

    mapped_lower = mapped_id.lower()
    if code >= 400:
        status = "HTTP_ERROR"
    elif source == "jaspar" and mapped_lower not in body:
        status = "CHECK_CONTENT"
    elif source == "hocomoco" and mapped_lower not in body and mapped_lower.split(".")[0] not in body:
        status = "CHECK_CONTENT"
    elif source == "cisbp" and mapped_id.startswith("T") and mapped_lower not in body:
        status = "CHECK_CONTENT"
    else:
        status = "OK"
    return {"source": source, "motif_id": motif_id, "mapped_id": mapped_id, "url": url, "status": status, "http_code": code, "error": ""}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--db", default="data/tf_webdb.sqlite")
    parser.add_argument("--cisbp-map", default="")
    parser.add_argument("--timeout", type=float, default=8.0)
    parser.add_argument("--workers", type=int, default=12)
    parser.add_argument("--limit", type=int, default=0, help="0 means all motifs")
    parser.add_argument("--out", default="outputs/external_motif_link_audit.tsv")
    args = parser.parse_args()

    cisbp_map = load_cisbp_map(args.cisbp_map or None)
    motifs = read_motifs(args.db)
    if args.limit:
        motifs = motifs[: args.limit]

    items = []
    for source, motif_id in motifs:
        _, mapped_id, url = source_url(source, motif_id, cisbp_map)
        items.append((source, motif_id, mapped_id, url))

    results = []
    with ThreadPoolExecutor(max_workers=max(1, args.workers)) as pool:
        futures = [pool.submit(check_url, item, args.timeout) for item in items]
        for future in as_completed(futures):
            results.append(future.result())

    results.sort(key=lambda row: (SOURCE_ORDER.get(str(row["source"]), 99), str(row["motif_id"])))
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["source", "motif_id", "mapped_id", "status", "http_code", "url", "error"], delimiter="\t")
        writer.writeheader()
        writer.writerows(results)

    by_source: dict[str, Counter[str]] = defaultdict(Counter)
    for row in results:
        by_source[str(row["source"])][str(row["status"])] += 1

    print(f"Audit TSV: {out}")
    print(f"CisBP map entries loaded: {len(cisbp_map)}")
    for source in sorted(by_source, key=lambda s: SOURCE_ORDER.get(s, 99)):
        counts = by_source[source]
        total = sum(counts.values())
        parts = ", ".join(f"{key}={value}" for key, value in sorted(counts.items()))
        print(f"{source}\ttotal={total}\t{parts}")

    bad = [row for row in results if row["status"] != "OK"][:25]
    if bad:
        print("\nFirst non-OK examples:")
        for row in bad:
            print(f"{row['source']}\t{row['motif_id']}\t{row['mapped_id']}\t{row['status']}\t{row['url']}\t{row['error']}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
