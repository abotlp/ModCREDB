#!/usr/bin/env python3
"""Audit or apply JASPAR metadata TF-to-motif links reproducibly.

The dry-run mode reads a ModCREDB SQLite database without changing it and
writes reviewable link candidates. Apply mode first copies the input database
to a new output path, then adds only exact JASPAR metadata links to that copy.
Existing evidence is never removed or overwritten.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import re
import shutil
import sqlite3
from collections import Counter
from pathlib import Path
from urllib.parse import urlparse
from urllib.request import Request, urlopen


SQL_TABLES = {"MATRIX", "MATRIX_PROTEIN", "MATRIX_SPECIES", "TAX"}
UNIPROT_SPLIT = re.compile(r"[;,|]\s*|\s+\|\s+")
PROJECT_ROOT = Path(__file__).resolve().parents[1]
SOURCE_RELEASE_COLUMNS = (
    "source",
    "release_name",
    "collection",
    "species_scope",
    "motif_or_model_type",
    "source_url",
    "download_url",
    "citation",
    "license_note",
    "downloaded_at",
    "local_file_label",
    "checksum_sha256",
    "confirmation_status",
    "notes",
)
JASPAR_LINK_NOTE = "Direct JASPAR 2024 metadata UniProt mapping from the official JASPAR SQL metadata dump."


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--db",
        type=Path,
        required=True,
        help="Existing database to audit or copy before applying links.",
    )
    parser.add_argument(
        "--metadata",
        required=True,
        help="Local metadata path or URL. Supports JASPAR2024 SQL(.gz) and TSV.",
    )
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument(
        "--dry-run",
        action="store_true",
        help="Read the database only and write reviewable candidates.",
    )
    mode.add_argument(
        "--apply",
        action="store_true",
        help="Copy --db to --output-db and add only verified direct JASPAR links to the copy.",
    )
    parser.add_argument(
        "--output-db",
        type=Path,
        help="New candidate database path. Required with --apply and must not already exist.",
    )
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument(
        "--report",
        type=Path,
        help="Optional Markdown apply/dry-run summary. The report contains no local absolute paths.",
    )
    return parser.parse_args()


def read_metadata_bytes(location: str) -> tuple[bytes, str]:
    parsed = urlparse(location)
    if parsed.scheme in {"http", "https"}:
        request = Request(location, headers={"User-Agent": "ModCREDB-JASPAR-audit/1.0"})
        with urlopen(request, timeout=60) as response:
            return response.read(), location
    path = Path(location)
    if not path.is_file():
        raise SystemExit(f"Metadata file not found: {path}")
    return path.read_bytes(), str(path)


def maybe_decompress_gzip(data: bytes, label: str) -> bytes:
    if data[:2] == b"\x1f\x8b" or label.lower().endswith(".gz"):
        return gzip.decompress(data)
    return data


def connect_read_only(path: Path) -> sqlite3.Connection:
    connection = sqlite3.connect(f"file:{path.resolve()}?mode=ro", uri=True)
    connection.row_factory = sqlite3.Row
    return connection


def metadata_from_sql(sql_bytes: bytes) -> list[dict[str, object]]:
    source = sqlite3.connect(":memory:")
    source.row_factory = sqlite3.Row
    try:
        source.executescript(sql_bytes.decode("utf-8"))
        tables = {
            row["name"]
            for row in source.execute("SELECT name FROM sqlite_master WHERE type = 'table'")
        }
        missing_tables = SQL_TABLES - tables
        if missing_tables:
            raise ValueError(f"Not a supported JASPAR SQL dump; missing tables: {sorted(missing_tables)}")
        rows = source.execute(
            """
            SELECT
                m.BASE_ID || '.' || m.VERSION AS motif_id,
                m.NAME AS jaspar_name,
                m.COLLECTION AS jaspar_collection,
                GROUP_CONCAT(DISTINCT tax.SPECIES) AS jaspar_species,
                GROUP_CONCAT(DISTINCT protein.ACC) AS jaspar_uniprot_ids
            FROM MATRIX AS m
            LEFT JOIN MATRIX_PROTEIN AS protein ON protein.ID = m.ID
            LEFT JOIN MATRIX_SPECIES AS species ON species.ID = m.ID
            LEFT JOIN TAX AS tax ON tax.TAX_ID = species.TAX_ID
            GROUP BY m.ID, m.BASE_ID, m.VERSION, m.NAME, m.COLLECTION
            ORDER BY motif_id
            """
        ).fetchall()
    finally:
        source.close()
    return [dict(row) for row in rows]


def normalized_header(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", value.lower())


def find_header(headers: list[str], candidates: tuple[str, ...], required: bool = False) -> str | None:
    normalized = {normalized_header(header): header for header in headers}
    for candidate in candidates:
        if candidate in normalized:
            return normalized[candidate]
    if required:
        raise ValueError(f"Metadata TSV lacks required column. Tried: {candidates}")
    return None


def metadata_from_tsv(data: bytes) -> list[dict[str, object]]:
    text = data.decode("utf-8-sig")
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    if not reader.fieldnames:
        raise ValueError("Metadata TSV has no header")
    headers = list(reader.fieldnames)
    motif_column = find_header(headers, ("matrixid", "motifid", "matrix"), required=True)
    name_column = find_header(headers, ("name", "tfname", "transcriptionfactor"))
    collection_column = find_header(headers, ("collection",))
    species_column = find_header(headers, ("species", "taxon", "taxonomicgroup"))
    uniprot_column = find_header(headers, ("uniprotid", "uniprotids", "uniprot", "uniprotaccession"))
    if uniprot_column is None:
        raise ValueError("Metadata TSV has no recognizable UniProt accession column")
    return [
        {
            "motif_id": str(row.get(motif_column) or "").strip(),
            "jaspar_name": str(row.get(name_column) or "").strip() if name_column else "",
            "jaspar_collection": str(row.get(collection_column) or "").strip() if collection_column else "",
            "jaspar_species": str(row.get(species_column) or "").strip() if species_column else "",
            "jaspar_uniprot_ids": str(row.get(uniprot_column) or "").strip(),
        }
        for row in reader
        if str(row.get(motif_column) or "").strip()
    ]


def load_metadata(location: str) -> tuple[list[dict[str, object]], str, str]:
    raw_bytes, source_label = read_metadata_bytes(location)
    checksum = hashlib.sha256(raw_bytes).hexdigest()
    data = maybe_decompress_gzip(raw_bytes, source_label)
    stripped = data.lstrip()
    if stripped.startswith(b"CREATE TABLE") or b"CREATE TABLE `MATRIX`" in data[:100000]:
        return metadata_from_sql(data), "JASPAR SQL dump", checksum
    return metadata_from_tsv(data), "JASPAR metadata TSV", checksum


def split_uniprot_ids(value: object) -> list[str]:
    return sorted({token.strip() for token in UNIPROT_SPLIT.split(str(value or "")) if token.strip()})


def db_inventory(
    connection: sqlite3.Connection,
) -> tuple[set[str], set[str], dict[tuple[str, str], list[sqlite3.Row]]]:
    tf_ids = {row["tf_id"] for row in connection.execute("SELECT tf_id FROM tf")}
    motif_ids = {
        row["motif_id"]
        for row in connection.execute("SELECT motif_id FROM motif_file WHERE source = 'jaspar'")
    }
    existing_links: dict[tuple[str, str], list[sqlite3.Row]] = {}
    for row in connection.execute(
        """
        SELECT motif_id, tf_id, evidence_type, mapping_type, original_column
        FROM motif_ref
        WHERE source = 'jaspar'
        """
    ):
        existing_links.setdefault((row["motif_id"], row["tf_id"]), []).append(row)
    return tf_ids, motif_ids, existing_links


def has_direct_link(rows: list[sqlite3.Row]) -> bool:
    """Return true when the pair already has direct/identical evidence.

    A pre-existing homology link must not suppress a direct JASPAR metadata
    link. Both are useful provenance and are intentionally preserved.
    """

    return any(
        row["evidence_type"] == "identical"
        or row["mapping_type"] == "direct_or_identical"
        for row in rows
    )


def build_candidate_rows(
    metadata_rows: list[dict[str, object]],
    tf_ids: set[str],
    motif_ids: set[str],
    existing_links: dict[tuple[str, str], list[sqlite3.Row]],
) -> tuple[list[dict[str, str]], int]:
    output_rows: list[dict[str, str]] = []
    pair_count = 0
    for metadata in metadata_rows:
        motif_id = str(metadata.get("motif_id") or "").strip()
        if not motif_id:
            continue
        uniprot_ids = split_uniprot_ids(metadata.get("jaspar_uniprot_ids"))
        common = {
            "motif_id": motif_id,
            "jaspar_name": str(metadata.get("jaspar_name") or ""),
            "jaspar_collection": str(metadata.get("jaspar_collection") or ""),
            "jaspar_species": str(metadata.get("jaspar_species") or ""),
        }
        if not uniprot_ids:
            output_rows.append(
                {
                    **common,
                    "jaspar_uniprot_id": "",
                    "motif_file_exists": "yes" if motif_id in motif_ids else "no",
                    "tf_exists": "not_assessed",
                    "existing_motif_ref": "no",
                    "proposed_evidence_type": "",
                    "proposed_mapping_type": "",
                    "action": "skip_no_uniprot_metadata",
                }
            )
            continue
        for uniprot_id in uniprot_ids:
            pair_count += 1
            motif_exists = motif_id in motif_ids
            tf_exists = uniprot_id in tf_ids
            pair_links = existing_links.get((motif_id, uniprot_id), [])
            link_exists = bool(pair_links)
            direct_link_exists = has_direct_link(pair_links)
            if not motif_exists:
                action = "missing_motif_file"
            elif not tf_exists:
                action = "missing_tf"
            elif direct_link_exists:
                action = "already_present"
            else:
                action = "add_link"
            output_rows.append(
                {
                    **common,
                    "jaspar_uniprot_id": uniprot_id,
                    "motif_file_exists": "yes" if motif_exists else "no",
                    "tf_exists": "yes" if tf_exists else "no",
                    "existing_motif_ref": "yes" if link_exists else "no",
                    "proposed_evidence_type": "identical" if action == "add_link" else "",
                    "proposed_mapping_type": "direct_or_identical" if action == "add_link" else "",
                    "action": action,
                }
            )
    output_rows.sort(key=lambda row: (row["motif_id"], row["jaspar_uniprot_id"]))
    return output_rows, pair_count


def write_candidates(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "motif_id",
        "jaspar_name",
        "jaspar_collection",
        "jaspar_species",
        "jaspar_uniprot_id",
        "motif_file_exists",
        "tf_exists",
        "existing_motif_ref",
        "proposed_evidence_type",
        "proposed_mapping_type",
        "action",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def sync_jaspar_source_release(connection: sqlite3.Connection) -> None:
    """Replace only the candidate DB's JASPAR provenance row from project config."""

    config_path = PROJECT_ROOT / "config" / "source_releases.tsv"
    with config_path.open(newline="", encoding="utf-8") as handle:
        jaspar_rows = [
            row for row in csv.DictReader(handle, delimiter="\t") if row.get("source") == "jaspar"
        ]
    if len(jaspar_rows) != 1:
        raise ValueError(f"Expected one JASPAR source-release row in {config_path.name}")
    row = jaspar_rows[0]
    values = tuple((row.get(column) or "").strip() for column in SOURCE_RELEASE_COLUMNS)
    connection.execute("DELETE FROM source_release WHERE source = 'jaspar'")
    connection.execute(
        """
        INSERT INTO source_release
            (source, release_name, collection, species_scope, motif_or_model_type,
             source_url, download_url, citation, license_note, downloaded_at,
             local_file_label, checksum_sha256, confirmation_status, notes)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        values,
    )


def apply_candidates(
    input_db: Path,
    output_db: Path,
    rows: list[dict[str, str]],
) -> int:
    if output_db.exists():
        raise FileExistsError(f"Refusing to overwrite existing candidate database: {output_db}")
    output_db.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(input_db, output_db)
    connection = sqlite3.connect(output_db)
    try:
        connection.execute("PRAGMA foreign_keys = ON")
        inserts = [row for row in rows if row["action"] == "add_link"]
        connection.executemany(
            """
            INSERT INTO motif_ref
                (tf_id, evidence_type, source, motif_id, original_value,
                 identity_percent, missing_local_file, original_column,
                 mapping_type, curation_status, evidence_note, display_priority)
            VALUES (?, 'identical', 'jaspar', ?, ?, NULL, 0,
                    'JASPAR2024_metadata', 'direct_or_identical',
                    'pending_confirmation', ?, 10)
            """,
            [
                (
                    row["jaspar_uniprot_id"],
                    row["motif_id"],
                    row["motif_id"],
                    JASPAR_LINK_NOTE,
                )
                for row in inserts
            ],
        )
        sync_jaspar_source_release(connection)
        connection.commit()
        return len(inserts)
    except Exception:
        connection.rollback()
        raise
    finally:
        connection.close()


def write_report(
    path: Path,
    *,
    mode: str,
    metadata_kind: str,
    metadata_source: str,
    checksum: str,
    pair_count: int,
    actions: Counter[str],
    added_links: int,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    output_db_label = "a separate candidate database" if mode == "apply" else "no database output"
    path.write_text(
        "\n".join(
            [
                "# JASPAR 2024 Metadata Link Import",
                "",
                f"Mode: `{mode}`.",
                f"Metadata format: {metadata_kind}.",
                f"Metadata source: {metadata_source}.",
                f"Metadata SHA256: `{checksum}`.",
                f"Result: {output_db_label}.",
                "",
                "## Candidate Counts",
                "",
                f"- Metadata motif-UniProt pairs: {pair_count}",
                f"- Existing direct links: {actions['already_present']}",
                f"- Direct links added/proposed: {actions['add_link']}",
                f"- Missing TF records: {actions['missing_tf']}",
                f"- Missing local motif files: {actions['missing_motif_file']}",
                f"- Metadata rows without UniProt IDs: {actions['skip_no_uniprot_metadata']}",
                f"- Applied rows: {added_links}",
                "",
                "## Semantics",
                "",
                "Inserted rows preserve existing evidence and use `evidence_type=identical`, "
                "`mapping_type=direct_or_identical`, `original_column=JASPAR2024_metadata`, and "
                "`curation_status=pending_confirmation`.",
                "",
                "A pre-existing homology link does not block a direct JASPAR metadata link; both "
                "evidence records are retained. No missing TFs or missing motif files are created.",
                "",
                "## Controls",
                "",
                "- `MA0106.3` / `P04637` is expected to receive a direct JASPAR metadata link.",
                "- `MA0065.1` / `P37231` is expected to remain an already-present direct link.",
                "",
            ]
        ),
        encoding="utf-8",
    )


def main() -> None:
    args = parse_args()
    if not args.db.is_file():
        raise SystemExit(f"Database not found: {args.db}")
    if args.apply and args.output_db is None:
        raise SystemExit("--output-db is required with --apply")
    if args.dry_run and args.output_db is not None:
        raise SystemExit("--output-db is only valid with --apply")

    metadata_rows, metadata_kind, checksum = load_metadata(args.metadata)
    with connect_read_only(args.db) as connection:
        tf_ids, motif_ids, existing_links = db_inventory(connection)
    output_rows, pair_count = build_candidate_rows(
        metadata_rows, tf_ids, motif_ids, existing_links
    )
    write_candidates(args.out, output_rows)

    actions = Counter(row["action"] for row in output_rows)
    motif_exists = sum(row["motif_file_exists"] == "yes" for row in output_rows)
    tf_exists = sum(row["tf_exists"] == "yes" for row in output_rows)
    added_links = 0
    if args.apply:
        assert args.output_db is not None
        added_links = apply_candidates(args.db, args.output_db, output_rows)
    if args.report:
        write_report(
            args.report,
            mode="apply" if args.apply else "dry-run",
            metadata_kind=metadata_kind,
            metadata_source=args.metadata,
            checksum=checksum,
            pair_count=pair_count,
            actions=actions,
            added_links=added_links,
        )
    print(f"Metadata source: {args.metadata}")
    print(f"Metadata format: {metadata_kind}")
    print(f"Metadata SHA256: {checksum}")
    print(f"Metadata rows parsed: {len(metadata_rows)}")
    print(f"Motif-UniProt pairs parsed: {pair_count}")
    print(f"Pairs where motif_file exists: {motif_exists}")
    print(f"Pairs where TF exists: {tf_exists}")
    print(f"Already-present direct links: {actions['already_present']}")
    print(f"New direct links proposed: {actions['add_link']}")
    print(f"Missing TF records: {actions['missing_tf']}")
    print(f"Missing motif files: {actions['missing_motif_file']}")
    print(f"No-UniProt metadata rows skipped: {actions['skip_no_uniprot_metadata']}")
    print(f"Wrote review candidates: {args.out}")
    if args.apply:
        print(f"Applied direct links to new candidate DB: {added_links}")
        print(f"Candidate database: {args.output_db}")
    if args.report:
        print(f"Wrote report: {args.report}")


if __name__ == "__main__":
    main()
