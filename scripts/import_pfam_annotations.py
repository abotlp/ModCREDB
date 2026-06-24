#!/usr/bin/env python3
"""Import real PFAM/InterPro annotations into a copied ModCREDB SQLite DB.

This creates a new tf_pfam_annotation layer. It never relabels or overwrites the
legacy chart-derived TF_family/tf_family values.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
import sqlite3
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import Request, urlopen

INTERPRO_URL = "https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/{accession}/"
PFAM_RE = re.compile(r"^PF\d{5}$")
IPR_RE = re.compile(r"^IPR\d{6}$")
COLUMNS = [
    "tf_id", "uniprot_accession", "gene", "protein_name", "pfam_id", "pfam_name",
    "pfam_type", "interpro_id", "interpro_name", "start", "end", "source",
    "source_release", "source_url", "action", "notes",
]
SCHEMA_SQL = """
CREATE TABLE IF NOT EXISTS tf_pfam_annotation (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    tf_id TEXT NOT NULL,
    uniprot_accession TEXT NOT NULL,
    pfam_id TEXT NOT NULL,
    pfam_name TEXT,
    pfam_type TEXT,
    interpro_id TEXT,
    interpro_name TEXT,
    start INTEGER,
    end INTEGER,
    source TEXT NOT NULL DEFAULT 'interpro_pfam',
    source_release TEXT,
    source_url TEXT,
    fetched_at TEXT,
    assignment_method TEXT,
    raw_json TEXT,
    UNIQUE(tf_id, pfam_id, start, end, source_release)
);
CREATE INDEX IF NOT EXISTS idx_tf_pfam_tf_id ON tf_pfam_annotation(tf_id);
CREATE INDEX IF NOT EXISTS idx_tf_pfam_pfam_id ON tf_pfam_annotation(pfam_id);
CREATE INDEX IF NOT EXISTS idx_tf_pfam_name ON tf_pfam_annotation(pfam_name);
CREATE INDEX IF NOT EXISTS idx_tf_pfam_interpro_id ON tf_pfam_annotation(interpro_id);
"""


def now_utc() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--db", type=Path, required=True)
    mode = p.add_mutually_exclusive_group(required=True)
    mode.add_argument("--dry-run", action="store_true")
    mode.add_argument("--apply", action="store_true")
    p.add_argument("--output-db", type=Path, help="Required with --apply; must not exist.")
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--report", type=Path)
    p.add_argument("--cache", type=Path, default=Path("data/cache/interpro_pfam"))
    p.add_argument("--source-release", default="InterPro/Pfam API current at fetch time")
    p.add_argument("--limit", type=int, help="Optional small preflight limit.")
    p.add_argument("--sleep-seconds", type=float, default=0.1)
    p.add_argument("--refresh-cache", action="store_true")
    a = p.parse_args()
    if a.apply and not a.output_db:
        p.error("--apply requires --output-db")
    if a.output_db and a.output_db.exists():
        p.error(f"--output-db already exists: {a.output_db}")
    return a


def connect_ro(db: Path) -> sqlite3.Connection:
    c = sqlite3.connect(f"file:{db.resolve()}?mode=ro", uri=True)
    c.row_factory = sqlite3.Row
    return c


def table_exists(c: sqlite3.Connection, name: str) -> bool:
    return c.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name=?", (name,)).fetchone() is not None


def tf_records(c: sqlite3.Connection, limit: int | None) -> list[dict[str, str]]:
    sql = """
    SELECT tf.tf_id, COALESCE(ta.gene_names, '') AS gene_names,
           COALESCE(ta.protein_name, '') AS protein_name
    FROM tf LEFT JOIN tf_annotation AS ta ON ta.tf_id = tf.tf_id
    ORDER BY tf.tf_id
    """
    rows = c.execute(sql + (" LIMIT ?" if limit else ""), ((limit,) if limit else ())).fetchall()
    out = []
    for r in rows:
        tf_id = str(r["tf_id"] or "").strip()
        if tf_id:
            out.append({
                "tf_id": tf_id,
                "uniprot_accession": tf_id,
                "gene": str(r["gene_names"] or "").split()[0] if r["gene_names"] else "",
                "protein_name": str(r["protein_name"] or ""),
            })
    return out


def existing_keys(c: sqlite3.Connection) -> set[tuple[str, str, str, str, str]]:
    if not table_exists(c, "tf_pfam_annotation"):
        return set()
    return {
        (str(r["tf_id"]), str(r["pfam_id"]), str(r["start"] or ""), str(r["end"] or ""), str(r["source_release"] or ""))
        for r in c.execute("SELECT tf_id, pfam_id, start, end, source_release FROM tf_pfam_annotation")
    }


def cache_file(cache: Path, accession: str) -> Path:
    return cache / (re.sub(r"[^A-Za-z0-9_.-]", "_", accession) + ".json")


def fetch_results(accession: str, cache: Path, refresh: bool, sleep: float) -> tuple[list[dict], str, str]:
    cache.mkdir(parents=True, exist_ok=True)
    path = cache_file(cache, accession)
    url = INTERPRO_URL.format(accession=quote(accession, safe=""))
    if path.exists() and not refresh:
        payload = json.loads(path.read_text(encoding="utf-8"))
        return payload.get("results", []), url, "cached"
    results: list[dict] = []
    notes = "fetched"
    next_url = url
    try:
        while next_url:
            req = Request(next_url, headers={"User-Agent": "ModCREDB-PFAM-import/1.0"})
            with urlopen(req, timeout=60) as resp:
                payload = json.loads(resp.read().decode("utf-8"))
            results.extend([x for x in payload.get("results", []) if isinstance(x, dict)])
            next_url = payload.get("next") or ""
    except HTTPError as e:
        notes = f"http_error_{e.code}"
    except URLError as e:
        notes = f"url_error_{e.reason}"
    except json.JSONDecodeError as e:
        notes = f"json_error_{e.msg}"
    path.write_text(json.dumps({"url": url, "fetched_at": now_utc(), "notes": notes, "results": results}, indent=2), encoding="utf-8")
    if notes == "fetched" and sleep > 0:
        time.sleep(sleep)
    return results, url, notes


def find_ipr(obj) -> tuple[str, str]:
    if isinstance(obj, dict):
        acc = str(obj.get("accession") or obj.get("id") or "")
        if IPR_RE.match(acc):
            return acc, str(obj.get("name") or obj.get("short_name") or "")
        for v in obj.values():
            found = find_ipr(v)
            if found[0]:
                return found
    if isinstance(obj, list):
        for v in obj:
            found = find_ipr(v)
            if found[0]:
                return found
    return "", ""


def fragments(obj) -> list[tuple[int | None, int | None]]:
    out = []
    def walk(x):
        if isinstance(x, dict):
            if isinstance(x.get("fragments"), list):
                for f in x["fragments"]:
                    if isinstance(f, dict):
                        out.append((to_int(f.get("start") or f.get("from")), to_int(f.get("end") or f.get("to"))))
            elif any(k in x for k in ("start", "from", "end", "to")):
                out.append((to_int(x.get("start") or x.get("from")), to_int(x.get("end") or x.get("to"))))
            for v in x.values():
                walk(v)
        elif isinstance(x, list):
            for v in x:
                walk(v)
    walk(obj)
    seen = set()
    uniq = []
    for pair in out or [(None, None)]:
        if pair not in seen:
            seen.add(pair)
            uniq.append(pair)
    return uniq


def to_int(x) -> int | None:
    try:
        return int(x) if x not in (None, "") else None
    except (TypeError, ValueError):
        return None


def parse_assignments(results: list[dict]) -> list[dict]:
    out = []
    for r in results:
        meta = r.get("metadata") if isinstance(r.get("metadata"), dict) else r
        pfam_id = str(meta.get("accession") or meta.get("id") or "") if isinstance(meta, dict) else ""
        if not PFAM_RE.match(pfam_id):
            continue
        ipr_id, ipr_name = find_ipr(r)
        for start, end in fragments(r):
            out.append({
                "pfam_id": pfam_id,
                "pfam_name": str(meta.get("name") or meta.get("short_name") or ""),
                "pfam_type": str(meta.get("type") or meta.get("entry_type") or ""),
                "interpro_id": ipr_id,
                "interpro_name": ipr_name,
                "start": start,
                "end": end,
            })
    seen = set(); uniq = []
    for a in out:
        key = (a["pfam_id"], a["start"], a["end"], a["interpro_id"])
        if key not in seen:
            seen.add(key); uniq.append(a)
    return uniq


def build_rows(args, records, existing):
    rows = []
    stats = Counter()
    no_pfam = []
    for rec in records:
        stats["tf_records"] += 1
        results, url, notes = fetch_results(rec["uniprot_accession"], args.cache, args.refresh_cache, args.sleep_seconds)
        stats[notes] += 1
        anns = parse_assignments(results)
        if anns:
            stats["accessions_with_pfam"] += 1
        else:
            no_pfam.append(rec["uniprot_accession"])
        for ann in anns:
            key = (rec["tf_id"], ann["pfam_id"], str(ann["start"] or ""), str(ann["end"] or ""), args.source_release)
            rows.append({
                **rec, **ann,
                "source": "interpro_pfam",
                "source_release": args.source_release,
                "source_url": url,
                "fetched_at": now_utc(),
                "assignment_method": "InterPro API PFAM matches for UniProt accession",
                "raw_json": "",
                "action": "already_present" if key in existing else "add_assignment",
                "notes": notes,
            })
            stats["assignments"] += 1
    return rows, stats, no_pfam


def write_tsv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as h:
        w = csv.DictWriter(h, fieldnames=COLUMNS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({c: "" if r.get(c) is None else r.get(c) for c in COLUMNS})


def apply_rows(db: Path, rows: list[dict]) -> tuple[int, int]:
    inserted = skipped = 0
    with sqlite3.connect(db) as c:
        c.executescript(SCHEMA_SQL)
        for r in rows:
            if r["action"] != "add_assignment":
                skipped += 1; continue
            before = c.total_changes
            c.execute("""
                INSERT OR IGNORE INTO tf_pfam_annotation
                (tf_id, uniprot_accession, pfam_id, pfam_name, pfam_type, interpro_id,
                 interpro_name, start, end, source, source_release, source_url, fetched_at,
                 assignment_method, raw_json)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, (r["tf_id"], r["uniprot_accession"], r["pfam_id"], r["pfam_name"], r["pfam_type"],
                  r["interpro_id"], r["interpro_name"], r["start"], r["end"], r["source"],
                  r["source_release"], r["source_url"], r["fetched_at"], r["assignment_method"], r["raw_json"]))
            inserted += 1 if c.total_changes > before else 0
            skipped += 0 if c.total_changes > before else 1
        c.commit()
    return inserted, skipped


def write_report(path: Path, rows: list[dict], stats: Counter, no_pfam: list[str], inserted=None, skipped=None) -> None:
    pfams = Counter(r["pfam_id"] for r in rows)
    names = {r["pfam_id"]: r.get("pfam_name", "") for r in rows}
    tf_by_pfam = defaultdict(set)
    for r in rows:
        tf_by_pfam[r["pfam_id"]].add(r["tf_id"])
    lines = [
        "# PFAM annotation report", "",
        "Real PFAM/InterPro annotations are stored separately from legacy TF family text.", "",
        f"- TF records queried: {stats['tf_records']}",
        f"- Accessions with PFAM: {stats['accessions_with_pfam']}",
        f"- PFAM assignment rows: {stats['assignments']}",
        f"- Unique PFAM IDs: {len(pfams)}",
    ]
    if inserted is not None:
        lines += [f"- Rows inserted: {inserted}", f"- Existing/skipped rows: {skipped}"]
    lines += ["", "## Top PFAM IDs", "", "| PFAM ID | Name | TF records | Rows |", "|---|---|---:|---:|"]
    for pfam_id, count in pfams.most_common(25):
        lines.append(f"| {pfam_id} | {names.get(pfam_id, '')} | {len(tf_by_pfam[pfam_id])} | {count} |")
    lines += ["", f"Accessions without PFAM assignment: {len(no_pfam)}"]
    if no_pfam:
        lines.append("First 50: " + ", ".join(no_pfam[:50]))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    with connect_ro(args.db) as c:
        integrity = c.execute("PRAGMA integrity_check").fetchone()[0]
        if integrity != "ok":
            raise SystemExit(f"integrity_check failed: {integrity}")
        records = tf_records(c, args.limit)
        existing = existing_keys(c)
    rows, stats, no_pfam = build_rows(args, records, existing)
    write_tsv(args.out, rows)
    inserted = skipped = None
    if args.apply:
        assert args.output_db is not None
        args.output_db.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(args.db, args.output_db)
        inserted, skipped = apply_rows(args.output_db, rows)
    if args.report:
        write_report(args.report, rows, stats, no_pfam, inserted, skipped)
    print(f"TF records queried: {stats['tf_records']}")
    print(f"Accessions with PFAM: {stats['accessions_with_pfam']}")
    print(f"PFAM assignment rows: {stats['assignments']}")
    if inserted is not None:
        print(f"Rows inserted into {args.output_db}: {inserted}; skipped: {skipped}")
    print(f"Wrote TSV: {args.out}")
    if args.report:
        print(f"Wrote report: {args.report}")


if __name__ == "__main__":
    main()
