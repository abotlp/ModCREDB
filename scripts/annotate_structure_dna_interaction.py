#!/usr/bin/env python3
"""Annotate active PDB structure records with protein-DNA interaction QC.

The QC rule matches the Baldo review diagnostic used for ModCREDB:
a model is DNA-interacting when at least one protein/nucleic-acid non-hydrogen
atom pair is within the configured cutoff (default 4.5 A).

The script updates per-structure derived QC columns on ``structure_file``.
It first reads an extracted model from ``--model-cache``; if that file is not
present, it falls back to the row's indexed ``archive_path`` + ``member_path``.
This is required for the curated Baldo ModCRE set, which is intentionally kept
inside ``baldo_validated_modcre_models.tar.gz`` rather than duplicated on disk.

The script is intentionally standard-library-only so it can run in the
production Python environment without installing packages.
"""

from __future__ import annotations

import argparse
import io
import math
import sqlite3
import tarfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable


PROTEIN_RESIDUES = frozenset(
    {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
        "ASX", "GLX", "HYP", "MSE", "PYL", "SEC", "XLE",
    }
)

NUCLEIC_RESIDUES = frozenset(
    {
        "DA", "DC", "DG", "DT", "DI",
        "A", "C", "G", "T", "U", "I",
        "ADE", "CYT", "GUA", "THY", "URI",
    }
)

QC_COLUMNS = {
    "dna_interaction": "INTEGER",
    "dna_contact_count": "INTEGER",
    "dna_qc_cutoff": "REAL",
    "dna_qc_status": "TEXT",
    "dna_qc_checked_at": "TEXT",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Annotate active structure_file PDB rows with DNA-interaction QC."
    )
    parser.add_argument("--db", required=True, type=Path, help="SQLite database to update")
    parser.add_argument(
        "--model-cache",
        required=True,
        type=Path,
        help="Root containing extracted structure_file member_path files",
    )
    parser.add_argument("--cutoff", type=float, default=4.5, help="Contact cutoff in Angstrom")
    parser.add_argument(
        "--apply",
        action="store_true",
        help="Write QC columns and values. Without this flag the run is read-only.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Recompute selected rows even when they already have QC at the same cutoff.",
    )
    parser.add_argument(
        "--only-unresolved",
        action="store_true",
        help="Select only rows whose current DNA QC status is absent or not 'ok'.",
    )
    parser.add_argument(
        "--ids",
        nargs="*",
        type=int,
        help="Optional structure_file IDs to restrict the run.",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=250,
        help="Print progress every N processed rows.",
    )
    return parser.parse_args()


def ensure_qc_columns(conn: sqlite3.Connection) -> None:
    existing = {row[1] for row in conn.execute("PRAGMA table_info(structure_file)")}
    for name, sql_type in QC_COLUMNS.items():
        if name not in existing:
            conn.execute(f"ALTER TABLE structure_file ADD COLUMN {name} {sql_type}")


def is_hydrogen(line: str, atom_name: str) -> bool:
    element = line[76:78].strip().upper() if len(line) >= 78 else ""
    if element in {"H", "D"}:
        return True
    if not element:
        stripped = atom_name.lstrip("0123456789").upper()
        return stripped.startswith(("H", "D"))
    return False


def parse_pdb_lines(
    lines: Iterable[str],
) -> tuple[list[tuple[float, float, float]], list[tuple[float, float, float]]]:
    protein: list[tuple[float, float, float]] = []
    nucleic: list[tuple[float, float, float]] = []

    for line in lines:
        if not (line.startswith("ATOM  ") or line.startswith("HETATM")):
            continue
        if len(line) < 54:
            continue

        residue = line[17:20].strip().upper()
        if residue not in PROTEIN_RESIDUES and residue not in NUCLEIC_RESIDUES:
            continue

        atom_name = line[12:16].strip().upper()
        if is_hydrogen(line, atom_name):
            continue

        try:
            xyz = (
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54]),
            )
        except ValueError:
            continue

        if residue in PROTEIN_RESIDUES:
            protein.append(xyz)
        else:
            nucleic.append(xyz)

    return protein, nucleic


def parse_pdb_path(
    path: Path,
) -> tuple[list[tuple[float, float, float]], list[tuple[float, float, float]]]:
    with path.open("rt", errors="replace") as handle:
        return parse_pdb_lines(handle)


def parse_pdb_archive_member(
    archive: tarfile.TarFile,
    member_path: str,
) -> tuple[list[tuple[float, float, float]], list[tuple[float, float, float]]]:
    try:
        member = archive.getmember(member_path)
    except KeyError as exc:
        raise FileNotFoundError(member_path) from exc
    raw = archive.extractfile(member)
    if raw is None:
        raise OSError(f"Could not read archive member: {member_path}")
    with raw:
        with io.TextIOWrapper(raw, encoding="ascii", errors="replace") as handle:
            return parse_pdb_lines(handle)


def contact_count(
    protein: list[tuple[float, float, float]],
    nucleic: list[tuple[float, float, float]],
    cutoff: float,
) -> int:
    """Count exact atom pairs within cutoff using a cubic spatial hash."""
    if not protein or not nucleic:
        return 0

    cutoff2 = cutoff * cutoff
    cell = cutoff
    grid: dict[tuple[int, int, int], list[tuple[float, float, float]]] = defaultdict(list)

    for x, y, z in nucleic:
        key = (math.floor(x / cell), math.floor(y / cell), math.floor(z / cell))
        grid[key].append((x, y, z))

    contacts = 0
    for x, y, z in protein:
        gx, gy, gz = (
            math.floor(x / cell),
            math.floor(y / cell),
            math.floor(z / cell),
        )
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for nx, ny, nz in grid.get((gx + dx, gy + dy, gz + dz), ()):
                        ddx = x - nx
                        ddy = y - ny
                        ddz = z - nz
                        if ddx * ddx + ddy * ddy + ddz * ddz <= cutoff2:
                            contacts += 1
    return contacts


def classify(
    row: sqlite3.Row,
    model_cache: Path,
    cutoff: float,
    archive_cache: dict[Path, tarfile.TarFile],
) -> tuple[int | None, int | None, str]:
    member_path = str(row["member_path"])
    local_path = model_cache / member_path

    try:
        if local_path.is_file():
            protein, nucleic = parse_pdb_path(local_path)
        else:
            archive_text = str(row["archive_path"] or "").strip()
            if not archive_text:
                return None, None, "missing_file"
            archive_path = Path(archive_text)
            if not archive_path.is_file():
                return None, None, "missing_file"
            archive = archive_cache.get(archive_path)
            if archive is None:
                archive = tarfile.open(archive_path, "r:gz")
                archive_cache[archive_path] = archive
            try:
                protein, nucleic = parse_pdb_archive_member(archive, member_path)
            except FileNotFoundError:
                return None, None, "missing_file"
    except (OSError, tarfile.TarError):
        return None, None, "read_error"
    except Exception:
        return None, None, "parse_error"

    if not protein:
        return None, 0, "no_protein_atoms"
    if not nucleic:
        return 0, 0, "no_nucleic_atoms"

    contacts = contact_count(protein, nucleic, cutoff)
    return (1 if contacts > 0 else 0), contacts, "ok"


def fetch_rows(
    conn: sqlite3.Connection,
    ids: list[int] | None,
    only_unresolved: bool,
) -> list[sqlite3.Row]:
    sql = """
        SELECT id, source, tf_id, model_id, member_path, archive_path,
               dna_interaction, dna_contact_count, dna_qc_cutoff, dna_qc_status
        FROM structure_file
        WHERE status = 'active' AND file_type = 'pdb'
    """
    params: list[object] = []
    if ids:
        sql += " AND id IN (%s)" % ",".join("?" for _ in ids)
        params.extend(ids)
    if only_unresolved:
        sql += " AND (dna_qc_status IS NULL OR dna_qc_status <> 'ok')"
    sql += " ORDER BY id"
    return conn.execute(sql, params).fetchall()


def main() -> int:
    args = parse_args()
    if args.cutoff <= 0:
        raise SystemExit("--cutoff must be > 0")
    if not args.db.is_file():
        raise SystemExit(f"Database not found: {args.db}")
    if not args.model_cache.is_dir():
        raise SystemExit(f"Model cache not found: {args.model_cache}")

    conn = sqlite3.connect(args.db)
    conn.row_factory = sqlite3.Row
    archive_cache: dict[Path, tarfile.TarFile] = {}

    try:
        if args.apply:
            conn.execute("BEGIN IMMEDIATE")
            ensure_qc_columns(conn)
        else:
            existing = {row[1] for row in conn.execute("PRAGMA table_info(structure_file)")}
            missing = [name for name in QC_COLUMNS if name not in existing]
            if missing:
                raise SystemExit(
                    "Dry-run requires existing QC columns. Run on a disposable DB copy with --apply first. "
                    f"Missing: {', '.join(missing)}"
                )

        rows = fetch_rows(conn, args.ids, args.only_unresolved)
        print(f"active PDB rows selected: {len(rows)}")
        print(f"contact rule: protein/nucleic non-H atom pair <= {args.cutoff:.2f} A")
        print(f"mode: {'APPLY' if args.apply else 'READ-ONLY'}")
        if args.only_unresolved:
            print("selection: unresolved QC rows only")

        summary = Counter()
        source_summary: dict[str, Counter] = defaultdict(Counter)
        processed = 0
        skipped = 0
        checked_at = datetime.now(timezone.utc).isoformat(timespec="seconds")

        for row in rows:
            if (
                not args.force
                and row["dna_qc_status"]
                and row["dna_qc_cutoff"] is not None
                and abs(float(row["dna_qc_cutoff"]) - args.cutoff) < 1e-9
            ):
                skipped += 1
                continue

            interaction, contacts, status = classify(
                row,
                args.model_cache,
                args.cutoff,
                archive_cache,
            )
            processed += 1

            label = "yes" if interaction == 1 else "no" if interaction == 0 else "unknown"
            summary[(status, label)] += 1
            source_summary[str(row["source"])][label] += 1
            source_summary[str(row["source"])][f"status:{status}"] += 1

            if args.apply:
                conn.execute(
                    """
                    UPDATE structure_file
                    SET dna_interaction = ?,
                        dna_contact_count = ?,
                        dna_qc_cutoff = ?,
                        dna_qc_status = ?,
                        dna_qc_checked_at = ?
                    WHERE id = ?
                    """,
                    (interaction, contacts, args.cutoff, status, checked_at, row["id"]),
                )

            if args.progress_every > 0 and processed % args.progress_every == 0:
                print(f"processed {processed}/{len(rows)}")

        if args.apply:
            conn.commit()

        print("\n=== QC SUMMARY ===")
        print(f"processed: {processed}")
        print(f"skipped existing: {skipped}")
        for (status, label), count in sorted(summary.items()):
            print(f"{status:20s} {label:7s} {count}")

        print("\n=== BY SOURCE ===")
        for source in sorted(source_summary):
            counts = source_summary[source]
            print(
                f"{source}: yes={counts['yes']} no={counts['no']} unknown={counts['unknown']}"
            )
            statuses = sorted(
                (key.split(':', 1)[1], value)
                for key, value in counts.items()
                if key.startswith("status:")
            )
            for status, count in statuses:
                print(f"  {status}: {count}")

        if args.apply:
            print("\n=== BALDO CHECK CASES ===")
            checks = conn.execute(
                """
                SELECT id, tf_id, source, model_id, dna_interaction,
                       dna_contact_count, dna_qc_status, dna_qc_cutoff
                FROM structure_file
                WHERE id IN (3381, 5846) OR tf_id = 'B3KN18'
                ORDER BY tf_id, id
                """
            ).fetchall()
            for check_row in checks:
                print(dict(check_row))

        return 0
    except Exception:
        if args.apply:
            conn.rollback()
        raise
    finally:
        for archive in archive_cache.values():
            archive.close()
        conn.close()


if __name__ == "__main__":
    raise SystemExit(main())
