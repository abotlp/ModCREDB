#!/usr/bin/env python3
"""Audit ModCREDB TF identity metadata against the official UniProt REST API.

The SQLite database is always opened read-only.  This script never updates it.
Official responses may be cached outside the repository so an interrupted audit
can be resumed without repeatedly querying UniProt.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sqlite3
import sys
import time
from collections import Counter
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATABASE = ROOT / "data" / "tf_webdb_local_staging.sqlite"
DEFAULT_CACHE = Path("/tmp/modcredb_uniprot_identity_audit_cache")
DEFAULT_OUTPUT = Path("/tmp/modcredb_uniprot_identity_audit.tsv")
DEFAULT_QC = Path("/tmp/modcredb_uniprot_identity_audit_qc.json")
UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
USER_AGENT = "ModCREDB-UniProt-read-only-audit/1.0"


OUTPUT_FIELDS = [
    "requested_accession",
    "mapping_status",
    "official_primary_accession",
    "official_entry_type",
    "official_inactive_reason",
    "official_deleted_reason",
    "official_merge_targets",
    "official_uniparc_id",
    "database_entry_name",
    "official_entry_name",
    "database_gene_names",
    "official_primary_gene",
    "gene_comparison",
    "database_protein_name",
    "official_primary_protein_name",
    "protein_name_comparison",
    "database_organism_id",
    "official_organism_id",
    "database_reviewed",
    "official_reviewed",
    "database_sequence_length",
    "official_sequence_length",
    "organism_comparison",
    "reviewed_comparison",
    "length_comparison",
    "database_uniprot_url",
    "expected_uniprot_url",
    "url_comparison",
]


class AuditError(RuntimeError):
    """Raised when the audit cannot be completed safely."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def chunks(items: list[str], size: int) -> list[list[str]]:
    return [items[index : index + size] for index in range(0, len(items), size)]


def normalize_text(value: object) -> str:
    return re.sub(r"\s+", " ", str(value or "")).strip()


def normalize_name(value: object) -> str:
    return re.sub(r"[^A-Z0-9]+", " ", normalize_text(value).upper()).strip()


def parse_reviewed(entry_type: object) -> int | None:
    value = normalize_text(entry_type).casefold()
    if "unreviewed" in value or "trembl" in value:
        return 0
    if "reviewed" in value or "swiss-prot" in value:
        return 1
    return None


def extract_primary_gene(record: dict[str, Any]) -> str:
    for gene in record.get("genes") or []:
        value = ((gene.get("geneName") or {}).get("value"))
        if value:
            return normalize_text(value)
    return ""


def nested_name(block: object) -> str:
    if not isinstance(block, dict):
        return ""
    full_name = block.get("fullName")
    if isinstance(full_name, dict):
        return normalize_text(full_name.get("value"))
    return normalize_text(full_name)


def extract_primary_protein(record: dict[str, Any]) -> str:
    description = record.get("proteinDescription") or {}
    value = nested_name(description.get("recommendedName"))
    if value:
        return value
    for key in ("submissionNames", "alternativeNames"):
        for block in description.get(key) or []:
            value = nested_name(block)
            if value:
                return value
    return ""


def extract_organism_id(record: dict[str, Any]) -> int | None:
    value = (record.get("organism") or {}).get("taxonId")
    try:
        return int(value) if value is not None else None
    except (TypeError, ValueError):
        return None


def extract_length(record: dict[str, Any]) -> int | None:
    value = (record.get("sequence") or {}).get("length")
    try:
        return int(value) if value is not None else None
    except (TypeError, ValueError):
        return None


def compare_scalar(database_value: object, official_value: object) -> str:
    if database_value is None and official_value is None:
        return "BOTH_MISSING"
    if database_value is None:
        return "DATABASE_MISSING"
    if official_value is None:
        return "OFFICIAL_MISSING"
    return "MATCH" if database_value == official_value else "DIFFERENT"


def compare_gene(database_gene_names: object, official_gene: object) -> str:
    database = normalize_text(database_gene_names)
    official = normalize_text(official_gene)
    if not database and not official:
        return "BOTH_MISSING"
    if not database:
        return "DATABASE_MISSING"
    if not official:
        return "OFFICIAL_MISSING"
    if database.casefold() == official.casefold():
        return "MATCH"
    first_database_gene = re.split(r"[\s,;/]+", database)[0]
    return "MATCH" if first_database_gene.casefold() == official.casefold() else "DIFFERENT"


def compare_protein(database_name: object, official_name: object) -> str:
    database = normalize_name(database_name)
    official = normalize_name(official_name)
    if not database and not official:
        return "BOTH_MISSING"
    if not database:
        return "DATABASE_MISSING"
    if not official:
        return "OFFICIAL_MISSING"
    if database == official:
        return "MATCH"
    if database.startswith(official + " ") or official.startswith(database + " "):
        return "MATCH_WITH_ALIASES"
    return "DIFFERENT"


def fetch_chunk(
    accessions: list[str], cache_dir: Path, timeout: int, offline: bool
) -> tuple[dict[str, Any], str, str]:
    key = hashlib.sha256("\n".join(accessions).encode("utf-8")).hexdigest()
    cache_path = cache_dir / f"{key}.json"
    metadata_path = cache_dir / f"{key}.headers.json"
    if cache_path.is_file():
        payload = json.loads(cache_path.read_text(encoding="utf-8"))
        metadata = (
            json.loads(metadata_path.read_text(encoding="utf-8"))
            if metadata_path.is_file()
            else {}
        )
        return payload, str(metadata.get("release", "")), str(metadata.get("release_date", ""))
    if offline:
        raise AuditError(f"missing cached UniProt response for {accessions[0]}..{accessions[-1]}")

    query = " OR ".join(f"accession:{accession}" for accession in accessions)
    url = f"{UNIPROT_SEARCH_URL}?{urlencode({'query': query, 'format': 'json', 'size': 500})}"
    request = Request(url, headers={"User-Agent": USER_AGENT, "Accept": "application/json"})
    try:
        with urlopen(request, timeout=timeout) as response:
            payload = json.loads(response.read().decode("utf-8"))
            release = response.headers.get("X-UniProt-Release", "")
            release_date = response.headers.get("X-UniProt-Release-Date", "")
    except (HTTPError, URLError, TimeoutError, json.JSONDecodeError) as exc:
        raise AuditError(
            f"official UniProt request failed for {accessions[0]}..{accessions[-1]}: {exc}"
        ) from exc
    cache_path.write_text(json.dumps(payload, sort_keys=True), encoding="utf-8")
    metadata_path.write_text(
        json.dumps({"release": release, "release_date": release_date}, indent=2) + "\n",
        encoding="utf-8",
    )
    return payload, release, release_date


def map_chunk_records(
    requested: list[str], records: list[dict[str, Any]]
) -> dict[str, tuple[str, dict[str, Any] | None]]:
    candidates: dict[str, list[tuple[str, dict[str, Any]]]] = {item: [] for item in requested}
    requested_set = set(requested)
    for record in records:
        primary = normalize_text(record.get("primaryAccession"))
        if primary in requested_set:
            candidates[primary].append(("EXACT_PRIMARY_ACCESSION", record))
        for secondary in record.get("secondaryAccessions") or []:
            secondary = normalize_text(secondary)
            if secondary in requested_set:
                candidates[secondary].append(("SECONDARY_ACCESSION", record))

    mapped: dict[str, tuple[str, dict[str, Any] | None]] = {}
    for accession, matches in candidates.items():
        if not matches:
            mapped[accession] = ("NOT_RETURNED", None)
        elif len(matches) == 1:
            mapped[accession] = matches[0]
        else:
            mapped[accession] = ("AMBIGUOUS_MULTIPLE_RECORDS", None)
    return mapped


def load_database_rows(database: Path) -> list[dict[str, Any]]:
    connection = sqlite3.connect(f"file:{database}?mode=ro", uri=True)
    connection.row_factory = sqlite3.Row
    try:
        if connection.execute("PRAGMA quick_check").fetchone()[0] != "ok":
            raise AuditError("SQLite quick_check failed")
        foreign_keys = connection.execute("PRAGMA foreign_key_check").fetchall()
        if foreign_keys:
            raise AuditError(f"SQLite foreign_key_check found {len(foreign_keys)} violations")
        rows = connection.execute(
            """
            SELECT tf.tf_id, ta.*
            FROM tf
            LEFT JOIN tf_annotation AS ta ON ta.tf_id = tf.tf_id
            ORDER BY tf.tf_id
            """
        ).fetchall()
        return [dict(row) for row in rows]
    finally:
        connection.close()


def build_output_row(
    database_row: dict[str, Any], mapping_status: str, record: dict[str, Any] | None
) -> dict[str, object]:
    requested = str(database_row["tf_id"])
    official_accession = normalize_text((record or {}).get("primaryAccession"))
    official_gene = extract_primary_gene(record or {})
    official_protein = extract_primary_protein(record or {})
    official_organism = extract_organism_id(record or {})
    official_reviewed = parse_reviewed((record or {}).get("entryType"))
    official_length = extract_length(record or {})
    expected_url = (
        f"https://www.uniprot.org/uniprotkb/{official_accession}/entry"
        if official_accession
        else ""
    )
    database_url = normalize_text(database_row.get("uniprot_url"))
    inactive_reason = (record or {}).get("inactiveReason") or {}
    if mapping_status == "EXACT_PRIMARY_ACCESSION" and normalize_text(
        (record or {}).get("entryType")
    ) == "Inactive":
        mapping_status = "EXACT_INACTIVE_ACCESSION"
    return {
        "requested_accession": requested,
        "mapping_status": mapping_status,
        "official_primary_accession": official_accession,
        "official_entry_type": normalize_text((record or {}).get("entryType")),
        "official_inactive_reason": normalize_text(inactive_reason.get("inactiveReasonType")),
        "official_deleted_reason": normalize_text(inactive_reason.get("deletedReason")),
        "official_merge_targets": ";".join(inactive_reason.get("mergeDemergeTo") or []),
        "official_uniparc_id": normalize_text(
            ((record or {}).get("extraAttributes") or {}).get("uniParcId")
        ),
        "database_entry_name": normalize_text(database_row.get("entry_name")),
        "official_entry_name": normalize_text((record or {}).get("uniProtkbId")),
        "database_gene_names": normalize_text(database_row.get("gene_names")),
        "official_primary_gene": official_gene,
        "gene_comparison": compare_gene(database_row.get("gene_names"), official_gene),
        "database_protein_name": normalize_text(database_row.get("protein_name")),
        "official_primary_protein_name": official_protein,
        "protein_name_comparison": compare_protein(database_row.get("protein_name"), official_protein),
        "database_organism_id": database_row.get("organism_id"),
        "official_organism_id": official_organism,
        "database_reviewed": database_row.get("reviewed"),
        "official_reviewed": official_reviewed,
        "database_sequence_length": database_row.get("sequence_length"),
        "official_sequence_length": official_length,
        "organism_comparison": compare_scalar(database_row.get("organism_id"), official_organism),
        "reviewed_comparison": compare_scalar(database_row.get("reviewed"), official_reviewed),
        "length_comparison": compare_scalar(database_row.get("sequence_length"), official_length),
        "database_uniprot_url": database_url,
        "expected_uniprot_url": expected_url,
        "url_comparison": "MATCH" if database_url == expected_url else "DIFFERENT",
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database", type=Path, default=DEFAULT_DATABASE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--qc-output", type=Path, default=DEFAULT_QC)
    parser.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--chunk-size", type=int, default=50)
    parser.add_argument("--timeout", type=int, default=60)
    parser.add_argument("--sleep", type=float, default=0.1)
    parser.add_argument("--offline", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not args.database.is_file():
        raise AuditError(f"database not found: {args.database}")
    if args.chunk_size < 1 or args.chunk_size > 100:
        raise AuditError("--chunk-size must be between 1 and 100")
    for output in (args.output, args.qc_output):
        if output.exists() and not args.overwrite:
            raise AuditError(f"refusing to overwrite existing output: {output}")
        output.parent.mkdir(parents=True, exist_ok=True)
    args.cache_dir.mkdir(parents=True, exist_ok=True)

    database_sha256_before = sha256_file(args.database)
    database_rows = load_database_rows(args.database)
    if len(database_rows) != 5417:
        raise AuditError(f"expected 5,417 TF rows, found {len(database_rows)}")
    accessions = [str(row["tf_id"]) for row in database_rows]
    if len(accessions) != len(set(accessions)):
        raise AuditError("duplicate tf_id values found")

    mappings: dict[str, tuple[str, dict[str, Any] | None]] = {}
    releases: set[str] = set()
    release_dates: set[str] = set()
    groups = chunks(accessions, args.chunk_size)
    for number, group in enumerate(groups, start=1):
        payload, release, release_date = fetch_chunk(
            group, args.cache_dir, args.timeout, args.offline
        )
        if release:
            releases.add(release)
        if release_date:
            release_dates.add(release_date)
        mappings.update(map_chunk_records(group, payload.get("results") or []))
        print(f"UniProt chunks: {number}/{len(groups)}", file=sys.stderr)
        if args.sleep and not args.offline:
            time.sleep(args.sleep)

    output_rows = [
        build_output_row(row, *mappings[str(row["tf_id"])])
        for row in database_rows
    ]
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)

    mapping_counts = Counter(str(row["mapping_status"]) for row in output_rows)
    gene_counts = Counter(str(row["gene_comparison"]) for row in output_rows)
    protein_counts = Counter(str(row["protein_name_comparison"]) for row in output_rows)
    organism_counts = Counter(str(row["organism_comparison"]) for row in output_rows)
    reviewed_counts = Counter(str(row["reviewed_comparison"]) for row in output_rows)
    length_counts = Counter(str(row["length_comparison"]) for row in output_rows)
    database_sha256_after = sha256_file(args.database)
    qc = {
        "database": str(args.database.resolve()),
        "database_sha256_before": database_sha256_before,
        "database_sha256_after": database_sha256_after,
        "database_unchanged": database_sha256_before == database_sha256_after,
        "tf_rows": len(output_rows),
        "uniprot_releases": sorted(releases),
        "uniprot_release_dates": sorted(release_dates),
        "mapping_status_counts": dict(sorted(mapping_counts.items())),
        "gene_comparison_counts": dict(sorted(gene_counts.items())),
        "protein_name_comparison_counts": dict(sorted(protein_counts.items())),
        "organism_comparison_counts": dict(sorted(organism_counts.items())),
        "reviewed_comparison_counts": dict(sorted(reviewed_counts.items())),
        "length_comparison_counts": dict(sorted(length_counts.items())),
    }
    args.qc_output.write_text(json.dumps(qc, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print("TF rows:", len(output_rows))
    print("UniProt release(s):", ", ".join(sorted(releases)) or "not recorded")
    print("Mapping status:", dict(sorted(mapping_counts.items())))
    print("Gene comparison:", dict(sorted(gene_counts.items())))
    print("Protein comparison:", dict(sorted(protein_counts.items())))
    print("Organism comparison:", dict(sorted(organism_counts.items())))
    print("Reviewed comparison:", dict(sorted(reviewed_counts.items())))
    print("Length comparison:", dict(sorted(length_counts.items())))
    print("Database unchanged:", "YES" if database_sha256_before == database_sha256_after else "NO")
    print("Output:", args.output)
    print("QC:", args.qc_output)

    if database_sha256_before != database_sha256_after:
        raise AuditError("database SHA256 changed during read-only audit")
    if len(output_rows) != 5417:
        raise AuditError("output row count is not 5,417")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except AuditError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
