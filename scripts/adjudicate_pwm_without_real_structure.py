#!/usr/bin/env python3
"""Adjudicate PWM-bearing accessions without a confirmed own-accession DNA structure.

This program is deliberately read-only with respect to SQLite and all model
inventories.  It writes only the requested TSV/JSON review products and refuses
to overwrite any existing output.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import sqlite3
import tarfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
MODELS_ROOT = Path("/data/sbi/interchange/boliva/patricia/models")
AF3_ROOT = Path("/data/sbi/interchange/boliva/patricia/AF3")
ALL_TF_FASTA = Path("/home/patricia/TF_database_Baldo_data/TF_without_model.fasta")
INTERFACE_CUTOFF = 4.5

OUTPUT_COLUMNS = [
    "uniprot_accession",
    "tf_id",
    "gene",
    "protein_name",
    "PWM_annotation_level",
    "PWM_count",
    "previous_structure_category",
    "current_web_3D_count",
    "exact_accession_model_files_found",
    "valid_exact_accession_models_found",
    "canonical_reference_accession",
    "canonical_reference_model_found",
    "expected_DBD",
    "DBD_present",
    "DBD_complete",
    "best_DBD_coordinates",
    "best_model_path",
    "best_model_source",
    "protein_present",
    "dna_present",
    "interface_valid",
    "interface_residue_count",
    "primary_reason",
    "detailed_explanation",
    "needs_baldo_review",
    "recommended_action",
]

PRIMARY_REASONS = {
    "VALID_MODEL_FOUND_BUT_NOT_INDEXED",
    "NO_MODEL_FILE_FOUND",
    "ONLY_EMPTY_INTERFACE_MODELS",
    "ONLY_WRONG_ACCESSION_OR_CANONICAL_REFERENCE",
    "ACCESSION_SEQUENCE_LACKS_DBD",
    "ACCESSION_CONTAINS_INCOMPLETE_DBD",
    "NON_DNA_BINDING_COMPONENT_OR_COFACTOR",
    "NON_DBD_FRAGMENT_ONLY",
    "FUSION_OR_BREAKPOINT_UNRESOLVED",
    "MODEL_EXISTS_BUT_INTERFACE_INVALID",
    "DBD_PRESENT_BUT_NO_VALID_DNA_BOUND_MODEL",
    "UNRESOLVED_MANUAL_REVIEW",
}

RECOMMENDED_ACTIONS = {
    "IMPORT_VALID_MODEL",
    "SEND_COMPLETE_DBD_FOR_MODELLING",
    "DO_NOT_SEND_ACCESSION_LACKS_DBD",
    "DO_NOT_SEND_INCOMPLETE_DBD",
    "DO_NOT_SEND_NON_DNA_BINDING_COMPONENT",
    "DO_NOT_SEND_NON_DBD_FRAGMENT",
    "CHECK_FUSION_MANUALLY",
    "MANUAL_REVIEW",
    "NO_ACTION_ALREADY_RESOLVED",
}

OUTPUT_NAMES = {
    "adjudication": "pwm_without_real_structure_adjudication.tsv",
    "true_gaps": "pwm_true_structure_gaps.tsv",
    "valid_found": "pwm_valid_models_found_during_review.tsv",
    "baldo": "pwm_structure_gaps_for_baldo.tsv",
    "manual": "pwm_without_real_structure_manual_review.tsv",
    "qc": "pwm_without_real_structure_adjudication_qc.json",
}

FUSION_UNRESOLVED = {"A0A0A7M2K0", "B1NY96", "H2BNB9"}
NON_DBD_FRAGMENT_CONTROLS = {"E9PKB7", "Q59EF3"}

# These eleven accessions have exact-accession local ModCRE coordinate files
# but no accepted interface.  Each mapping is supported by a reviewed local
# UniProt record plus global sequence alignment; three otherwise gene-less
# records use an exact reviewed protein-name match documented in the output.
SEQUENCE_REFERENCE_MAP: dict[str, tuple[str, str]] = {
    "A0A0K0TMT5": ("P01100", "PF00170"),
    "A0A0K0TMT7": ("P01100", "PF00170"),
    "A0A1B0GWI9": ("O00712", "PF03165"),
    "A0A224A4N4": ("P11473", "PF00105"),
    "A0A994J6X9": ("P35711", "PF00505"),
    "B3KVM5": ("O75469", "PF00105"),
    "B4DS34": ("Q16665", "PF23171"),
    "C9JWL3": ("Q9NZI7", "PF04516"),
    "E9PK33": ("Q96BA8", "PF00170"),
    "F8W133": ("P35638", "PF07716"),
    "Q6LDR0": ("P10827", "PF00105"),
}

# These accessions have an accepted full expected DBD interval directly on the
# accession.  The cited Pfam names in SQLite are explicit DNA-binding families;
# no inference is made for other Pfams or for records lacking this evidence.
DIRECT_COMPLETE_DBD: dict[str, str] = {
    "E7EQY4": "PF00249",
    "E7ESY4": "PF00249",
    "E7EV10": "PF00249",
    "E9PCS8": "PF00249",
    "F6RRE2": "PF00249",
    "Q59FW1": "PF00249",
    "Q6ZNB4": "PF02257",
    "Q9BTC8": "PF00249",
}

AA_RESIDUES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
    "TYR", "VAL", "ASX", "GLX", "SEC", "PYL", "MSE",
}
NUCLEIC_RESIDUES = {
    "A", "C", "G", "T", "U", "I", "DA", "DC", "DG", "DT", "DU", "DI",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build the read-only whole-database PWM structure-gap adjudication.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--database", type=Path, default=ROOT / "data/tf_webdb_local_staging.sqlite")
    parser.add_argument("--production-database", type=Path, default=ROOT / "data/tf_webdb.sqlite")
    parser.add_argument("--pwm-audit", type=Path, default=ROOT / "outputs/pwm_accession_structure_audit.tsv")
    parser.add_argument("--pwm-without", type=Path, default=ROOT / "outputs/pwm_without_valid_structure.tsv")
    parser.add_argument("--pwm-manual", type=Path, default=ROOT / "outputs/pwm_structure_manual_review.tsv")
    parser.add_argument("--pwm-visible", type=Path, default=ROOT / "outputs/pwm_with_valid_visible_structure.tsv")
    parser.add_argument("--pwm-qc", type=Path, default=ROOT / "outputs/pwm_accession_structure_audit_qc.json")
    parser.add_argument("--final-adjudication", type=Path, default=ROOT / "outputs/af3_failed_102_final_adjudication_DRAFT.tsv")
    parser.add_argument("--import-preview", type=Path, default=ROOT / "outputs/local_af3_model_import_preview.tsv")
    parser.add_argument("--af3-vs-database", type=Path, default=ROOT / "external/baldo_model_inventory/af3_vs_database.tsv")
    parser.add_argument("--reviewed-pfam", type=Path, default=ROOT / "data_sources/af3_failed_102_pfam_role_reviewed.tsv")
    parser.add_argument("--all-tf-fasta", type=Path, default=ALL_TF_FASTA)
    parser.add_argument("--models-root", type=Path, default=MODELS_ROOT)
    parser.add_argument("--af3-root", type=Path, default=AF3_ROOT)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "outputs")
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=OUTPUT_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def database_signature(path: Path) -> dict[str, Any]:
    stat = path.stat()
    connection = sqlite3.connect(f"file:{path.resolve()}?mode=ro", uri=True)
    connection.execute("PRAGMA query_only=ON")
    try:
        quick_check = connection.execute("PRAGMA quick_check").fetchone()[0]
        row_counts = {
            table: int(connection.execute(f'SELECT COUNT(*) FROM "{table}"').fetchone()[0])
            for (table,) in connection.execute(
                "SELECT name FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%' ORDER BY name"
            )
        }
    finally:
        connection.close()
    return {
        "path": str(path.resolve()),
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "sha256": file_sha256(path),
        "quick_check": quick_check,
        "row_counts": row_counts,
    }


def file_signature(path: Path, include_hash: bool = False) -> dict[str, Any]:
    stat = path.stat()
    result: dict[str, Any] = {
        "path": str(path.resolve()), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns,
    }
    if include_hash:
        result["sha256"] = file_sha256(path)
    return result


def fasta_records(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    accession = ""
    chunks: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if accession:
                    records[accession] = "".join(chunks)
                accession = line[1:].split("|", 1)[0].split()[0].upper()
                chunks = []
            else:
                chunks.append(line)
    if accession:
        records[accession] = "".join(chunks)
    return records


def global_alignment(target: str, canonical: str) -> dict[str, Any]:
    """Needleman-Wunsch alignment using the established project scoring."""
    n, m = len(target), len(canonical)
    gap = -2
    previous = [j * gap for j in range(m + 1)]
    trace = [bytearray(m + 1) for _ in range(n + 1)]
    for j in range(1, m + 1):
        trace[0][j] = 3
    for i in range(1, n + 1):
        trace[i][0] = 2
        current = [i * gap] + [0] * m
        aa = target[i - 1]
        for j in range(1, m + 1):
            diagonal = previous[j - 1] + (2 if aa == canonical[j - 1] else -1)
            up = previous[j] + gap
            left = current[j - 1] + gap
            best = max(diagonal, up, left)
            current[j] = best
            trace[i][j] = 1 if best == diagonal else (2 if best == up else 3)
        previous = current
    i, j = n, m
    pairs: list[tuple[int | None, int | None, bool]] = []
    while i or j:
        direction = trace[i][j]
        if direction == 1:
            pairs.append((i, j, target[i - 1] == canonical[j - 1]))
            i -= 1
            j -= 1
        elif direction == 2:
            pairs.append((i, None, False))
            i -= 1
        else:
            pairs.append((None, j, False))
            j -= 1
    pairs.reverse()
    mapping = {canonical_pos: (target_pos, match) for target_pos, canonical_pos, match in pairs if target_pos and canonical_pos}
    aligned = len(mapping)
    matches = sum(match for _, match in mapping.values())
    return {
        "mapping": mapping,
        "aligned": aligned,
        "identity": 100.0 * matches / aligned if aligned else 0.0,
    }


def dbd_alignment_metrics(alignment: dict[str, Any], start: int, end: int) -> dict[str, Any]:
    mapped = [
        (position, alignment["mapping"][position])
        for position in range(start, end + 1)
        if position in alignment["mapping"]
    ]
    covered = len(mapped)
    matches = sum(value[1] for _, value in mapped)
    length = end - start + 1
    target_positions = [value[0] for _, value in mapped]
    non_dbd = [
        value for position, value in alignment["mapping"].items()
        if position < start or position > end
    ]
    non_matches = sum(value[1] for value in non_dbd)
    return {
        "length": length,
        "covered": covered,
        "coverage": 100.0 * covered / length,
        "identity": 100.0 * matches / covered if covered else 0.0,
        "target_coordinates": f"{min(target_positions)}-{max(target_positions)}" if target_positions else "",
        "non_dbd_aligned": len(non_dbd),
        "non_dbd_identity": 100.0 * non_matches / len(non_dbd) if non_dbd else 0.0,
    }


def model_interval(path: str) -> tuple[int, int] | None:
    match = re.search(r":(\d+):(\d+)_", Path(path.split("::")[-1]).name)
    return (int(match.group(1)), int(match.group(2))) if match else None


def intervals_overlap(first: tuple[int, int] | None, second: tuple[int, int] | None) -> bool:
    if not first or not second:
        return False
    return max(first[0], second[0]) <= min(first[1], second[1])


def inspect_pdb_text(text: str) -> dict[str, Any]:
    protein_atoms: list[tuple[float, float, float, str, str]] = []
    dna_atoms: list[tuple[float, float, float, str, str]] = []
    for line in text.splitlines():
        if not (line.startswith("ATOM  ") or line.startswith("HETATM")):
            continue
        if len(line) < 54:
            continue
        residue = line[17:20].strip().upper()
        try:
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError:
            continue
        chain = line[21:22].strip()
        resid = line[22:27].strip()
        atom = (xyz[0], xyz[1], xyz[2], chain, resid)
        if residue in AA_RESIDUES:
            protein_atoms.append(atom)
        elif residue in NUCLEIC_RESIDUES:
            dna_atoms.append(atom)
    result = {
        "readable": "YES",
        "protein_present": "YES" if protein_atoms else "NO",
        "dna_present": "YES" if dna_atoms else "NO",
        "atom_contacts": 0,
        "protein_interface_residues": 0,
        "result": "FAIL_NO_PROTEIN" if not protein_atoms else ("FAIL_NO_DNA" if not dna_atoms else "FAIL_EMPTY_PROTEIN_DNA_INTERFACE"),
        "error": "",
    }
    if not protein_atoms or not dna_atoms:
        return result
    cell = INTERFACE_CUTOFF
    grid: dict[tuple[int, int, int], list[tuple[float, float, float, str, str]]] = defaultdict(list)
    for atom in dna_atoms:
        grid[(math.floor(atom[0] / cell), math.floor(atom[1] / cell), math.floor(atom[2] / cell))].append(atom)
    cutoff_sq = INTERFACE_CUTOFF * INTERFACE_CUTOFF
    contacts = 0
    interface_residues: set[tuple[str, str]] = set()
    for x, y, z, chain, resid in protein_atoms:
        key = (math.floor(x / cell), math.floor(y / cell), math.floor(z / cell))
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for nx, ny, nz, _, _ in grid.get((key[0] + dx, key[1] + dy, key[2] + dz), []):
                        if (x - nx) ** 2 + (y - ny) ** 2 + (z - nz) ** 2 <= cutoff_sq:
                            contacts += 1
                            interface_residues.add((chain, resid))
    result["atom_contacts"] = contacts
    result["protein_interface_residues"] = len(interface_residues)
    result["result"] = "PASS_NONEMPTY_PROTEIN_DNA_INTERFACE" if contacts else "FAIL_EMPTY_PROTEIN_DNA_INTERFACE"
    return result


def inspect_pdb_path(path: Path) -> dict[str, Any]:
    try:
        return inspect_pdb_text(path.read_text(encoding="utf-8", errors="replace"))
    except Exception as exc:
        return {
            "readable": "NO", "protein_present": "NO", "dna_present": "NO",
            "atom_contacts": 0, "protein_interface_residues": 0,
            "result": "FAIL_UNREADABLE", "error": f"{type(exc).__name__}: {exc}",
        }


def exact_gene_tokens(value: str) -> set[str]:
    return {token.upper() for token in re.split(r"[\s,;/]+", value or "") if token}


def main() -> int:
    args = parse_args()
    inputs = [
        args.database, args.production_database, args.pwm_audit, args.pwm_without,
        args.pwm_manual, args.pwm_visible, args.pwm_qc, args.final_adjudication,
        args.import_preview, args.af3_vs_database, args.reviewed_pfam,
        args.all_tf_fasta, args.models_root, args.af3_root,
    ]
    missing = [str(path) for path in inputs if not path.exists()]
    if missing:
        raise SystemExit("Missing required input(s): " + ", ".join(missing))
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs = {key: output_dir / name for key, name in OUTPUT_NAMES.items()}
    existing = [str(path) for path in outputs.values() if path.exists()]
    if existing:
        raise SystemExit("Refusing to overwrite existing output(s): " + ", ".join(existing))

    database_before = {
        "staging": database_signature(args.database),
        "production": database_signature(args.production_database),
    }
    protected_files = [
        ROOT / "app.py", ROOT / "templates/model.html", ROOT / "templates/tf_clean_models.html",
        Path("/home/patricia/TF_database_Baldo_data/models.tar.gz"),
    ]
    protected_before = {
        str(path): file_signature(path, include_hash=path.suffix in {".py", ".html"})
        for path in protected_files
    }
    protected_dirs = [args.models_root.resolve(), args.af3_root.resolve()]
    protected_dirs_before = {str(path): file_signature(path) for path in protected_dirs}

    pwm_rows = read_tsv(args.pwm_audit)
    pwm_by_id = {row["tf_id"]: row for row in pwm_rows}
    no_valid = {row["tf_id"] for row in pwm_rows if row["final_structure_category"] == "NO_VALID_STRUCTURE"}
    canonical_only = {row["tf_id"] for row in pwm_rows if row["final_structure_category"] == "CANONICAL_REFERENCE_ONLY"}
    manual_group = {row["tf_id"] for row in pwm_rows if row["final_structure_category"] == "MANUAL_REVIEW"}
    invalid_visible = {
        row["tf_id"] for row in pwm_rows
        if int(row["current_web_3D_count"]) > 0 and int(row["valid_indexed_structure_count"]) == 0
    }
    population = no_valid | canonical_only | manual_group | invalid_visible
    group_sets = {
        "NO_VALID_STRUCTURE": no_valid,
        "CANONICAL_REFERENCE_ONLY": canonical_only,
        "MANUAL_REVIEW": manual_group,
        "INVALID_VISIBLE": invalid_visible,
    }

    final_by_id = {row["failed_accession"]: row for row in read_tsv(args.final_adjudication)}
    preview_by_id: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in read_tsv(args.import_preview):
        if row["accession"] in population:
            preview_by_id[row["accession"]].append(row)
    af3_inventory_by_id: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in read_tsv(args.af3_vs_database):
        if row["tf_id"] in population:
            af3_inventory_by_id[row["tf_id"]].append(row)
    reviewed_roles = {row["pfam_id"]: row for row in read_tsv(args.reviewed_pfam)}
    sequences = fasta_records(args.all_tf_fasta)

    connection = sqlite3.connect(f"file:{args.database.resolve()}?mode=ro", uri=True)
    connection.row_factory = sqlite3.Row
    connection.execute("PRAGMA query_only=ON")
    annotation = {row["tf_id"]: dict(row) for row in connection.execute("SELECT * FROM tf_annotation")}
    family = {row["tf_id"]: row["family_text"] for row in connection.execute("SELECT tf_id,family_text FROM tf")}
    pfams: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in connection.execute("SELECT * FROM tf_pfam_annotation ORDER BY tf_id,start,end"):
        pfams[row["tf_id"]].append(dict(row))
    structure_rows: dict[str, list[dict[str, Any]]] = defaultdict(list)
    requested_archive_members: dict[str, dict[str, list[dict[str, Any]]]] = defaultdict(lambda: defaultdict(list))
    for row_obj in connection.execute("SELECT * FROM structure_file WHERE file_type='pdb'"):
        row = dict(row_obj)
        tf_id = row["tf_id"]
        if tf_id in population:
            structure_rows[tf_id].append(row)
            requested_archive_members[row["archive_path"]][row["member_path"]].append(row)
    assignments: dict[str, list[dict[str, Any]]] = defaultdict(list)
    if connection.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name='structure_model_assignment'").fetchone():
        for row in connection.execute("SELECT * FROM structure_model_assignment"):
            if row["display_accession"] in population:
                assignments[row["display_accession"]].append(dict(row))
    connection.close()

    candidates: dict[str, dict[str, dict[str, Any]]] = defaultdict(dict)

    # Directly validate every archive member referenced by the review population.
    for archive_path, members in requested_archive_members.items():
        archive = Path(archive_path)
        seen: set[str] = set()
        if archive.is_file() and os.access(archive, os.R_OK):
            with tarfile.open(archive, "r:gz") as handle:
                for member in handle:
                    if member.name not in members:
                        continue
                    seen.add(member.name)
                    extracted = handle.extractfile(member)
                    scan = (
                        inspect_pdb_text(extracted.read().decode("utf-8", errors="replace"))
                        if extracted is not None
                        else {
                            "readable": "NO", "protein_present": "NO", "dna_present": "NO",
                            "atom_contacts": 0, "protein_interface_residues": 0,
                            "result": "FAIL_UNREADABLE", "error": "archive member not readable",
                        }
                    )
                    for db_row in members[member.name]:
                        accession = db_row["tf_id"]
                        key = f"{archive_path}::{member.name}"
                        existing_candidate = candidates[accession].get(key)
                        statuses = {db_row["status"]}
                        if existing_candidate:
                            statuses.update(existing_candidate["database_statuses"])
                        candidates[accession][key] = {
                            "path": key,
                            "source": f"{str(db_row['source']).upper()}_INDEXED_ARCHIVE",
                            "database_statuses": sorted(statuses),
                            "protein_present": scan["protein_present"],
                            "dna_present": scan["dna_present"],
                            "atom_contacts": int(scan["atom_contacts"]),
                            "interface_residues": int(scan["protein_interface_residues"]),
                            "interface_result": scan["result"],
                            "readable": scan["readable"],
                            "biological_role": "UNKNOWN",
                            "model_interval": (db_row["residue_start"], db_row["residue_end"])
                            if db_row["residue_start"] is not None and db_row["residue_end"] is not None
                            else model_interval(member.name),
                        }
        for member_name, db_rows in members.items():
            if member_name in seen:
                continue
            for db_row in db_rows:
                accession = db_row["tf_id"]
                key = f"{archive_path}::{member_name}"
                candidates[accession][key] = {
                    "path": key, "source": f"{str(db_row['source']).upper()}_INDEXED_ARCHIVE",
                    "database_statuses": [db_row["status"]], "protein_present": "NO",
                    "dna_present": "NO", "atom_contacts": 0, "interface_residues": 0,
                    "interface_result": "FAIL_MISSING_FILE", "readable": "NO",
                    "biological_role": "UNKNOWN", "model_interval": model_interval(member_name),
                }

    # Validate exact-accession local ModCRE PDB files.
    for accession in sorted(population):
        for path in sorted(args.models_root.glob(f"TFS_{accession}:*.pdb")) + sorted(args.models_root.glob(f"DIMER_{accession}:*.pdb")):
            scan = inspect_pdb_path(path)
            candidates[accession][str(path)] = {
                "path": str(path), "source": "MODCRE_LOCAL_INVENTORY", "database_statuses": [],
                "protein_present": scan["protein_present"], "dna_present": scan["dna_present"],
                "atom_contacts": int(scan["atom_contacts"]),
                "interface_residues": int(scan["protein_interface_residues"]),
                "interface_result": scan["result"], "readable": scan["readable"],
                "biological_role": "UNKNOWN", "model_interval": model_interval(str(path)),
            }

    # Use the existing project audit for AF3 CIF candidates and reviewed fragments.
    for accession, rows in preview_by_id.items():
        for row in rows:
            path = row["candidate_model_path"]
            biological_role = row["candidate_biological_role"] or "UNKNOWN"
            if row["decision"] == "REJECT_NON_DBD_FRAGMENT" or row["reviewed_expected_DBD"] == "NO":
                biological_role = "NON_DBD_FRAGMENT"
            elif row["reviewed_expected_DBD"] == "YES":
                biological_role = "EXPECTED_DBD_FRAGMENT"
            candidates[accession][path] = {
                "path": path, "source": "AF3_REVIEWED_CANDIDATE",
                "database_statuses": [], "protein_present": row["protein_present"],
                "dna_present": row["dna_present"], "atom_contacts": int(row["atom_contacts"] or 0),
                "interface_residues": int(row["protein_interface_residues"] or 0),
                "interface_result": row["interface_result"],
                "readable": row["model_file_readable"], "biological_role": biological_role,
                "model_interval": None, "decision": row["decision"],
                "domain_or_fragment_id": row["domain_or_fragment_id"],
            }

    # Include direct assignments if any population member has one; none are assumed valid from status alone.
    for accession, rows in assignments.items():
        for row in rows:
            path = row["model_path"]
            if path in candidates[accession]:
                continue
            candidates[accession][path] = {
                "path": path, "source": f"{str(row['source']).upper()}_DIRECT_ASSIGNMENT",
                "database_statuses": ["active" if int(row["is_active"]) else "inactive"],
                "protein_present": "YES" if Path(path).is_file() else "NO",
                "dna_present": "YES" if str(row["interface_status"] or "").startswith("PASS") else "UNCERTAIN",
                "atom_contacts": 1 if str(row["interface_status"] or "").startswith("PASS") else 0,
                "interface_residues": 0,
                "interface_result": row["interface_status"] or "UNVALIDATED",
                "readable": "YES" if Path(path).is_file() and os.access(path, os.R_OK) else "NO",
                "biological_role": row["model_role"], "model_interval": model_interval(path),
            }

    # Build reviewed canonical selection by exact primary gene token.
    reviewed_by_gene: dict[str, list[str]] = defaultdict(list)
    for accession, row in annotation.items():
        if int(row.get("reviewed") or 0) != 1:
            continue
        for token in exact_gene_tokens(row.get("gene_names") or ""):
            reviewed_by_gene[token].append(accession)

    sequence_evidence: dict[str, dict[str, Any]] = {}
    for accession, (canonical, dbd_pfam) in SEQUENCE_REFERENCE_MAP.items():
        canonical_domains = [domain for domain in pfams[canonical] if domain["pfam_id"] == dbd_pfam]
        if accession not in sequences or canonical not in sequences or len(canonical_domains) != 1:
            raise RuntimeError(f"Missing unique local sequence/DBD evidence for {accession} -> {canonical} {dbd_pfam}")
        domain = canonical_domains[0]
        alignment = global_alignment(sequences[accession], sequences[canonical])
        metrics = dbd_alignment_metrics(alignment, int(domain["start"]), int(domain["end"]))
        if metrics["covered"] >= 30 and metrics["coverage"] >= 80 and metrics["identity"] >= 80:
            disposition = "COMPLETE"
        elif metrics["covered"] >= 10 and 10 < metrics["coverage"] < 80 and metrics["identity"] >= 80:
            disposition = "INCOMPLETE"
        elif (
            metrics["coverage"] <= 10
            and metrics["non_dbd_aligned"] >= 40
            and metrics["non_dbd_identity"] >= 80
        ):
            disposition = "LACKS"
        else:
            disposition = "UNCERTAIN"
        sequence_evidence[accession] = {
            "canonical": canonical,
            "domain": domain,
            "metrics": metrics,
            "alignment": alignment,
            "disposition": disposition,
        }

    direct_dbd: dict[str, dict[str, Any]] = {}
    for accession, pfam_id in DIRECT_COMPLETE_DBD.items():
        matches = [domain for domain in pfams[accession] if domain["pfam_id"] == pfam_id]
        if len(matches) != 1:
            raise RuntimeError(f"Expected exactly one direct complete DBD {pfam_id} for {accession}")
        direct_dbd[accession] = matches[0]

    def canonical_model_found(accession: str) -> str:
        if not accession or accession not in pwm_by_id:
            return "NO"
        row = pwm_by_id[accession]
        return "YES" if int(row["valid_indexed_structure_count"]) + int(row["valid_unindexed_structure_count"]) > 0 else "NO"

    output_rows: list[dict[str, Any]] = []
    candidate_debug: dict[str, list[dict[str, Any]]] = {}
    for accession in sorted(population):
        source = pwm_by_id[accession]
        ann = annotation.get(accession, {})
        gene = (ann.get("gene_names") or "").split()[0] if ann.get("gene_names") else ""
        final = final_by_id.get(accession)
        accession_candidates = list(candidates.get(accession, {}).values())
        candidate_debug[accession] = accession_candidates

        canonical_accession = ""
        canonical_found = "NO"
        expected_dbd = "UNRESOLVED_FROM_LOCAL_EVIDENCE"
        dbd_present = "UNCERTAIN"
        dbd_complete = "UNCERTAIN"
        best_dbd_coordinates = ""
        evidence_kind = "UNRESOLVED"
        explanation_parts: list[str] = []
        expected_interval: tuple[int, int] | None = None

        if final:
            canonical_accession = final["reviewed_canonical_accession"]
            canonical_found = (
                "YES" if final["canonical_reference_model"] and final["canonical_reference_interface_status"].startswith("PASS")
                else canonical_model_found(canonical_accession)
            )
            expected_dbd = final["expected_DBD_type"] or "UNRESOLVED_FROM_FINAL_AF3_102_REVIEW"
            best_dbd_coordinates = final["failed_coordinates_corresponding_to_DBD"]
            if best_dbd_coordinates and re.fullmatch(r"\d+-\d+", best_dbd_coordinates):
                start, end = map(int, best_dbd_coordinates.split("-"))
                expected_interval = (start, end)
            if final["primary_structural_status"] == "ACCESSION_LACKS_DBD":
                dbd_present, dbd_complete, evidence_kind = "NO", "NO", "LACKS"
            elif final["primary_structural_status"] == "ACCESSION_CONTAINS_INCOMPLETE_DBD":
                dbd_present, dbd_complete, evidence_kind = "YES", "NO", "INCOMPLETE"
            elif final["primary_structural_status"] == "COFACTOR_OR_NON_DNA_BINDING_COMPONENT":
                dbd_present, dbd_complete, evidence_kind = "NO", "NO", "COFACTOR"
            elif accession in FUSION_UNRESOLVED:
                evidence_kind = "FUSION"
            elif accession in NON_DBD_FRAGMENT_CONTROLS:
                evidence_kind = "NON_DBD_FRAGMENT"
            else:
                evidence_kind = "UNRESOLVED"
            explanation_parts.append(
                "AF3-102 final draft evidence: "
                f"DBD_disposition={final['DBD_disposition']}; primary_structural_status={final['primary_structural_status']}; "
                f"DBD coverage={final['DBD_coverage_percent'] or 'not resolved'}%; identity={final['DBD_identity_percent'] or 'not resolved'}%. "
                f"{final['decision_reason']}"
            )
            if final["remaining_uncertainty"]:
                explanation_parts.append("Remaining uncertainty: " + final["remaining_uncertainty"])
        elif accession in sequence_evidence:
            evidence = sequence_evidence[accession]
            domain = evidence["domain"]
            metrics = evidence["metrics"]
            canonical_accession = evidence["canonical"]
            canonical_found = canonical_model_found(canonical_accession)
            expected_dbd = f"{domain['pfam_id']}:{domain['pfam_name']}"
            best_dbd_coordinates = metrics["target_coordinates"]
            expected_interval = (
                tuple(map(int, best_dbd_coordinates.split("-")))
                if best_dbd_coordinates else None
            )
            evidence_kind = evidence["disposition"]
            if evidence_kind == "COMPLETE":
                dbd_present, dbd_complete = "YES", "YES"
            elif evidence_kind == "INCOMPLETE":
                dbd_present, dbd_complete = "YES", "NO"
            elif evidence_kind == "LACKS":
                dbd_present, dbd_complete = "NO", "NO"
            explanation_parts.append(
                f"Global alignment to reviewed canonical {canonical_accession} maps {metrics['covered']}/{metrics['length']} expected DBD residues "
                f"({metrics['coverage']:.2f}% coverage, {metrics['identity']:.2f}% identity); "
                f"{metrics['non_dbd_aligned']} non-DBD residues align at {metrics['non_dbd_identity']:.2f}% identity. "
                f"The conservative sequence disposition is {evidence_kind}."
            )
            if accession in {"B3KVM5", "B4DS34", "Q6LDR0"}:
                explanation_parts.append(
                    "The reviewed canonical was selected by exact reviewed protein identity/description because the accession gene field is blank; this choice is recorded explicitly and no low-identity automatic absence call is made."
                )
        elif accession in direct_dbd:
            domain = direct_dbd[accession]
            expected_dbd = f"{domain['pfam_id']}:{domain['pfam_name']}"
            best_dbd_coordinates = f"{domain['start']}-{domain['end']}"
            expected_interval = (int(domain["start"]), int(domain["end"]))
            dbd_present, dbd_complete, evidence_kind = "YES", "YES", "COMPLETE"
            explanation_parts.append(
                f"The accession has an accepted full {domain['pfam_id']} {domain['pfam_name']} interval at {domain['start']}-{domain['end']}; "
                "this directly supports a complete expected DBD without inferring sequence absence from model failure."
            )
            if gene:
                reviewed = [item for item in reviewed_by_gene.get(gene.upper(), []) if item != accession]
                if reviewed:
                    canonical_accession = sorted(reviewed)[0]
                    canonical_found = canonical_model_found(canonical_accession)
        else:
            if gene:
                reviewed = [item for item in reviewed_by_gene.get(gene.upper(), []) if item != accession]
                if reviewed:
                    canonical_accession = sorted(reviewed)[0]
                    canonical_found = canonical_model_found(canonical_accession)

        # Assign biological relevance without allowing contacts alone to create a DBD model.
        for candidate in accession_candidates:
            if candidate["biological_role"] == "NON_DBD_FRAGMENT":
                continue
            if candidate["biological_role"] == "EXPECTED_DBD_FRAGMENT":
                continue
            if evidence_kind == "COMPLETE" and intervals_overlap(candidate.get("model_interval"), expected_interval):
                candidate["biological_role"] = "EXPECTED_DBD_INTERVAL"
            elif candidate["source"] == "AF3_REVIEWED_CANDIDATE" and candidate["biological_role"] == "EXACT_ACCESSION_FULL_LENGTH_TF_MODEL":
                candidate["biological_role"] = "FULL_LENGTH_WITH_EXPECTED_DBD" if evidence_kind == "COMPLETE" else "FULL_LENGTH_DBD_UNRESOLVED"

        valid_candidates = [
            candidate for candidate in accession_candidates
            if candidate["readable"] == "YES"
            and candidate["protein_present"] == "YES"
            and candidate["dna_present"] == "YES"
            and int(candidate["atom_contacts"]) > 0
            and candidate["biological_role"] in {
                "EXPECTED_DBD_FRAGMENT", "EXPECTED_DBD_INTERVAL", "FULL_LENGTH_WITH_EXPECTED_DBD",
            }
        ]
        best_candidate = None
        if accession_candidates:
            best_candidate = sorted(
                accession_candidates,
                key=lambda item: (
                    0 if item in valid_candidates else 1,
                    -int(item["atom_contacts"]),
                    0 if item["source"] == "MODCRE_LOCAL_INVENTORY" else 1,
                    item["path"],
                ),
            )[0]

        if valid_candidates:
            primary_reason = "VALID_MODEL_FOUND_BUT_NOT_INDEXED"
            recommended_action = "IMPORT_VALID_MODEL"
            needs_review = "NO"
            explanation_parts.append(
                f"Careful review found {len(valid_candidates)} valid exact-accession DBD model(s) with protein, DNA, and a nonempty interface; the best is not currently indexed as valid."
            )
        elif evidence_kind == "LACKS":
            primary_reason = "ACCESSION_SEQUENCE_LACKS_DBD"
            recommended_action = "DO_NOT_SEND_ACCESSION_LACKS_DBD"
            needs_review = "NO"
        elif evidence_kind == "INCOMPLETE":
            primary_reason = "ACCESSION_CONTAINS_INCOMPLETE_DBD"
            recommended_action = "DO_NOT_SEND_INCOMPLETE_DBD"
            needs_review = "YES" if accession == "D1LUZ8" else "NO"
        elif evidence_kind == "COFACTOR":
            primary_reason = "NON_DNA_BINDING_COMPONENT_OR_COFACTOR"
            recommended_action = "DO_NOT_SEND_NON_DNA_BINDING_COMPONENT"
            needs_review = "NO"
        elif evidence_kind == "FUSION":
            primary_reason = "FUSION_OR_BREAKPOINT_UNRESOLVED"
            recommended_action = "CHECK_FUSION_MANUALLY"
            needs_review = "YES"
        elif evidence_kind == "NON_DBD_FRAGMENT":
            primary_reason = "NON_DBD_FRAGMENT_ONLY"
            recommended_action = "DO_NOT_SEND_NON_DBD_FRAGMENT"
            needs_review = "YES"
        elif evidence_kind == "COMPLETE":
            primary_reason = "DBD_PRESENT_BUT_NO_VALID_DNA_BOUND_MODEL"
            recommended_action = "SEND_COMPLETE_DBD_FOR_MODELLING"
            needs_review = "NO"
        elif evidence_kind == "UNRESOLVED" and (final or accession in sequence_evidence):
            # Preserve the explicit accession/sequence uncertainty from the
            # prior adjudication even when every attempted model is empty.
            primary_reason = "UNRESOLVED_MANUAL_REVIEW"
            recommended_action = "MANUAL_REVIEW"
            needs_review = "YES"
        elif not accession_candidates:
            primary_reason = "NO_MODEL_FILE_FOUND"
            recommended_action = "MANUAL_REVIEW"
            needs_review = "YES"
        elif all(
            candidate["readable"] == "YES"
            and candidate["protein_present"] == "YES"
            and candidate["dna_present"] == "YES"
            and int(candidate["atom_contacts"]) == 0
            for candidate in accession_candidates
        ):
            primary_reason = "ONLY_EMPTY_INTERFACE_MODELS"
            recommended_action = "MANUAL_REVIEW"
            needs_review = "YES"
        elif any(candidate["biological_role"] == "NON_DBD_FRAGMENT" and int(candidate["atom_contacts"]) > 0 for candidate in accession_candidates):
            primary_reason = "NON_DBD_FRAGMENT_ONLY"
            recommended_action = "DO_NOT_SEND_NON_DBD_FRAGMENT"
            needs_review = "YES"
        elif any(int(candidate["atom_contacts"]) > 0 for candidate in accession_candidates):
            # A geometric interface without defensible expected-DBD relevance
            # is not an invalid interface and is never promoted to a structure.
            primary_reason = "UNRESOLVED_MANUAL_REVIEW"
            recommended_action = "MANUAL_REVIEW"
            needs_review = "YES"
        else:
            primary_reason = "MODEL_EXISTS_BUT_INTERFACE_INVALID"
            recommended_action = "MANUAL_REVIEW"
            needs_review = "YES"

        if accession_candidates:
            result_counts = Counter(candidate["interface_result"] for candidate in accession_candidates)
            explanation_parts.append(
                f"Exact-accession coordinate evidence: {len(accession_candidates)} unique audited model file(s); "
                + ", ".join(f"{key}={value}" for key, value in sorted(result_counts.items()))
                + "."
            )
        else:
            explanation_parts.append(
                "No exact-accession coordinate candidate is present in structure_file/direct assignments, the local ModCRE PDB inventory, or the audited AF3 candidate inventory. Model absence is not treated as evidence that the sequence lacks a DBD."
            )
        if canonical_accession:
            explanation_parts.append(
                f"Reviewed/same-gene canonical accession {canonical_accession} is kept separate; canonical_reference_model_found={canonical_found} and it is not counted as this accession's structure."
            )
        if accession in af3_inventory_by_id:
            inventory_classes = sorted({row["db_class"] for row in af3_inventory_by_id[accession]})
            explanation_parts.append("AF3 inventory class(es): " + ",".join(inventory_classes) + ".")

        row = {
            "uniprot_accession": ann.get("uniprot_accession") or accession,
            "tf_id": accession,
            "gene": gene or (final["gene"] if final else ""),
            "protein_name": ann.get("protein_name") or (final["protein_description"] if final else ""),
            "PWM_annotation_level": source["PWM_annotation_level"],
            "PWM_count": source["PWM_count"],
            "previous_structure_category": source["final_structure_category"],
            "current_web_3D_count": source["current_web_3D_count"],
            "exact_accession_model_files_found": len(accession_candidates),
            "valid_exact_accession_models_found": len(valid_candidates),
            "canonical_reference_accession": canonical_accession,
            "canonical_reference_model_found": canonical_found,
            "expected_DBD": expected_dbd,
            "DBD_present": dbd_present,
            "DBD_complete": dbd_complete,
            "best_DBD_coordinates": best_dbd_coordinates,
            "best_model_path": best_candidate["path"] if best_candidate else "",
            "best_model_source": best_candidate["source"] if best_candidate else "",
            "protein_present": best_candidate["protein_present"] if best_candidate else "",
            "dna_present": best_candidate["dna_present"] if best_candidate else "",
            "interface_valid": "YES" if best_candidate and int(best_candidate["atom_contacts"]) > 0 else ("NO" if best_candidate else ""),
            "interface_residue_count": best_candidate["interface_residues"] if best_candidate else "",
            "primary_reason": primary_reason,
            "detailed_explanation": " ".join(part.strip() for part in explanation_parts if part.strip()),
            "needs_baldo_review": needs_review,
            "recommended_action": recommended_action,
        }
        output_rows.append(row)

    output_rows.sort(key=lambda row: row["uniprot_accession"])
    valid_found_rows = [row for row in output_rows if row["primary_reason"] == "VALID_MODEL_FOUND_BUT_NOT_INDEXED"]
    true_gap_rows = [row for row in output_rows if row["primary_reason"] != "VALID_MODEL_FOUND_BUT_NOT_INDEXED"]
    action_priority = {
        "SEND_COMPLETE_DBD_FOR_MODELLING": 0,
        "MANUAL_REVIEW": 1,
        "CHECK_FUSION_MANUALLY": 2,
        "DO_NOT_SEND_INCOMPLETE_DBD": 3,
        "DO_NOT_SEND_ACCESSION_LACKS_DBD": 4,
        "DO_NOT_SEND_NON_DNA_BINDING_COMPONENT": 5,
        "DO_NOT_SEND_NON_DBD_FRAGMENT": 6,
        "IMPORT_VALID_MODEL": 7,
        "NO_ACTION_ALREADY_RESOLVED": 8,
    }
    baldo_rows = sorted(
        true_gap_rows,
        key=lambda row: (action_priority[row["recommended_action"]], row["uniprot_accession"]),
    )
    manual_rows = [row for row in output_rows if row["needs_baldo_review"] == "YES"]
    reason_counts = Counter(row["primary_reason"] for row in output_rows)
    action_counts = Counter(row["recommended_action"] for row in output_rows)

    checks: dict[str, dict[str, Any]] = {}

    def check(name: str, passed: bool, observed: Any, expected: Any) -> None:
        checks[name] = {"passed": bool(passed), "observed": observed, "expected": expected}

    check("input_group_counts", {key: len(value) for key, value in group_sets.items()} == {
        "NO_VALID_STRUCTURE": 273, "CANONICAL_REFERENCE_ONLY": 78,
        "MANUAL_REVIEW": 21, "INVALID_VISIBLE": 27,
    }, {key: len(value) for key, value in group_sets.items()}, {
        "NO_VALID_STRUCTURE": 273, "CANONICAL_REFERENCE_ONLY": 78,
        "MANUAL_REVIEW": 21, "INVALID_VISIBLE": 27,
    })
    check("invalid_visible_overlap", invalid_visible <= no_valid and len(invalid_visible & no_valid) == 27, len(invalid_visible & no_valid), 27)
    check("unique_population", len(population) == 372, len(population), 372)
    check("one_output_row_per_accession", len(output_rows) == len({row["tf_id"] for row in output_rows}) == 372, len(output_rows), 372)
    check("output_population_exact", {row["tf_id"] for row in output_rows} == population, len({row["tf_id"] for row in output_rows} ^ population), 0)
    check("allowed_reasons", set(reason_counts) <= PRIMARY_REASONS, sorted(reason_counts), sorted(PRIMARY_REASONS))
    check("allowed_actions", set(action_counts) <= RECOMMENDED_ACTIONS, sorted(action_counts), sorted(RECOMMENDED_ACTIONS))
    check("valid_found_models_fully_valid", all(
        int(row["valid_exact_accession_models_found"]) > 0
        and row["protein_present"] == "YES" and row["dna_present"] == "YES"
        and row["interface_valid"] == "YES" and row["best_model_path"]
        for row in valid_found_rows
    ), len(valid_found_rows), "all valid-found rows pass exact file/protein/DNA/interface/DBD checks")
    check(
        "valid_found_exact_accession_paths",
        all(
            row["uniprot_accession"].upper() in Path(row["best_model_path"].split("::")[-1]).name.upper()
            for row in valid_found_rows
        ),
        [row["uniprot_accession"] for row in valid_found_rows],
        "every best model filename contains the exact accession",
    )
    check(
        "newly_resolved_valid_models",
        {row["uniprot_accession"] for row in valid_found_rows} == {"E7ESY4", "Q59FW1"},
        sorted(row["uniprot_accession"] for row in valid_found_rows),
        ["E7ESY4", "Q59FW1"],
    )
    check("no_canonical_counted_as_own", all(
        not row["canonical_reference_accession"]
        or row["canonical_reference_accession"] != row["uniprot_accession"]
        for row in output_rows
    ), sum(row["canonical_reference_accession"] == row["uniprot_accession"] for row in output_rows), 0)
    check("non_dbd_fragments_not_valid", all(
        not any(
            candidate["biological_role"] == "NON_DBD_FRAGMENT" and candidate in [valid for valid in candidate_debug[row["tf_id"]] if valid["readable"] == "YES" and valid["protein_present"] == "YES" and valid["dna_present"] == "YES" and int(valid["atom_contacts"]) > 0 and valid["biological_role"] in {"EXPECTED_DBD_FRAGMENT", "EXPECTED_DBD_INTERVAL", "FULL_LENGTH_WITH_EXPECTED_DBD"}]
            for candidate in candidate_debug[row["tf_id"]]
        )
        for row in output_rows
    ), 0, 0)
    failed_102_population = population & set(final_by_id)
    check("failed_102_count", len(failed_102_population) == 101, len(failed_102_population), 101)
    check("empty_full_length_af3_not_rescued", all(
        not (
            row["candidate_type"] == "FULL_LENGTH_AF3"
            and row["interface_result"] == "FAIL_EMPTY_PROTEIN_DNA_INTERFACE"
            and any(output["tf_id"] == accession and output["primary_reason"] == "VALID_MODEL_FOUND_BUT_NOT_INDEXED" for output in output_rows)
        )
        for accession in failed_102_population for row in preview_by_id[accession]
    ), 0, 0)
    check("send_rows_have_complete_dbd", all(
        row["DBD_present"] == "YES" and row["DBD_complete"] == "YES"
        and row["best_DBD_coordinates"] and int(row["valid_exact_accession_models_found"]) == 0
        for row in output_rows if row["recommended_action"] == "SEND_COMPLETE_DBD_FOR_MODELLING"
    ), action_counts["SEND_COMPLETE_DBD_FOR_MODELLING"], "all have complete mapped DBD and no valid own model")
    check("true_gaps_exclude_valid_found", len(true_gap_rows) + len(valid_found_rows) == 372 and not ({row["tf_id"] for row in true_gap_rows} & {row["tf_id"] for row in valid_found_rows}), (len(true_gap_rows), len(valid_found_rows)), "partition of 372")
    check("baldo_contains_only_true_gaps", {row["tf_id"] for row in baldo_rows} == {row["tf_id"] for row in true_gap_rows}, len(baldo_rows), len(true_gap_rows))
    check("manual_file_matches_flag", {row["tf_id"] for row in manual_rows} == {row["tf_id"] for row in output_rows if row["needs_baldo_review"] == "YES"}, len(manual_rows), "all needs_baldo_review=YES")

    by_id = {row["tf_id"]: row for row in output_rows}
    controls = {
        "P01100": {
            "pass": "P01100" not in population
            and pwm_by_id["P01100"]["final_structure_category"] == "VALID_STRUCTURE_NOT_INDEXED"
            and int(pwm_by_id["P01100"]["valid_unindexed_structure_count"]) > 0
            and Path(pwm_by_id["P01100"]["best_valid_model_path"]).is_file(),
            "result": "Excluded: valid exact-P01100 local ModCRE model exists and is not indexed.",
        },
        "H7C4N4": {
            "pass": "H7C4N4" not in population
            and pwm_by_id["H7C4N4"]["final_structure_category"] == "VALID_STRUCTURE_VISIBLE"
            and "PF03299" in pwm_by_id["H7C4N4"]["best_valid_model_path"],
            "result": "Excluded: active valid exact-accession PF03299 DBD fragment is visible.",
        },
        "E9PN75": {
            "pass": by_id["E9PN75"]["valid_exact_accession_models_found"] == 0
            and by_id["E9PN75"]["canonical_reference_accession"] == "Q9NZC4"
            and by_id["E9PN75"]["primary_reason"] == "ACCESSION_SEQUENCE_LACKS_DBD",
            "result": by_id["E9PN75"],
        },
        "B4DHE0": {
            "pass": by_id["B4DHE0"]["valid_exact_accession_models_found"] == 0
            and by_id["B4DHE0"]["primary_reason"] == "ACCESSION_SEQUENCE_LACKS_DBD"
            and any(candidate["biological_role"] == "NON_DBD_FRAGMENT" for candidate in candidate_debug["B4DHE0"]),
            "result": by_id["B4DHE0"],
        },
        "Q59EF3": {
            "pass": by_id["Q59EF3"]["valid_exact_accession_models_found"] == 0
            and by_id["Q59EF3"]["primary_reason"] == "NON_DBD_FRAGMENT_ONLY"
            and any(candidate["biological_role"] == "NON_DBD_FRAGMENT" for candidate in candidate_debug["Q59EF3"]),
            "result": by_id["Q59EF3"],
        },
    }
    check("all_controls_pass", all(item["pass"] for item in controls.values()), {key: item["pass"] for key, item in controls.items()}, {key: True for key in controls})

    failed_checks = [name for name, value in checks.items() if not value["passed"]]
    if failed_checks:
        raise RuntimeError("QC failed before output creation: " + ", ".join(failed_checks))

    write_tsv(outputs["adjudication"], output_rows)
    write_tsv(outputs["true_gaps"], true_gap_rows)
    write_tsv(outputs["valid_found"], valid_found_rows)
    write_tsv(outputs["baldo"], baldo_rows)
    write_tsv(outputs["manual"], manual_rows)

    database_after = {
        "staging": database_signature(args.database),
        "production": database_signature(args.production_database),
    }
    protected_after = {
        str(path): file_signature(path, include_hash=path.suffix in {".py", ".html"})
        for path in protected_files
    }
    protected_dirs_after = {str(path): file_signature(path) for path in protected_dirs}
    check("databases_unchanged", database_before == database_after, database_after, database_before)
    check("database_quick_checks_ok", all(item["quick_check"] == "ok" for item in database_after.values()), {key: value["quick_check"] for key, value in database_after.items()}, {"staging": "ok", "production": "ok"})
    check("website_and_model_archive_files_unchanged", protected_before == protected_after, protected_after, protected_before)
    check("model_directories_unchanged", protected_dirs_before == protected_dirs_after, protected_dirs_after, protected_dirs_before)
    failed_checks = [name for name, value in checks.items() if not value["passed"]]
    if failed_checks:
        raise RuntimeError("Post-output QC failed: " + ", ".join(failed_checks))

    overlap = {
        f"{first}__{second}": len(group_sets[first] & group_sets[second])
        for index, first in enumerate(group_sets)
        for second in list(group_sets)[index + 1:]
    }
    qc = {
        "status": "PASS",
        "hostname_required": "masada",
        "read_only": True,
        "remote_connection_occurred": False,
        "population": {
            "group_counts": {key: len(value) for key, value in group_sets.items()},
            "pairwise_overlap_counts": overlap,
            "sum_before_deduplication": sum(len(value) for value in group_sets.values()),
            "unique_accessions_reviewed": len(population),
        },
        "scientific_rules": {
            "valid_structure": "Exact UniProt accession; readable model; protein and DNA; >=1 protein-DNA atom contact at 4.5 A; expected DBD/full-length complete-DBD relevance; canonical references and non-DBD contacts excluded.",
            "sequence_complete_threshold": "At least 30 canonical DBD residues, >=80% DBD coverage and >=80% DBD identity.",
            "sequence_incomplete_threshold": "At least 10 DBD residues, >10% and <80% DBD coverage, and >=80% identity.",
            "sequence_absence_threshold": "<=10% DBD coverage plus at least 40 non-DBD aligned residues at >=80% identity.",
            "no_model_caution": "No model file never implies DBD absence; such cases remain unresolved unless separate sequence/domain evidence is decisive.",
        },
        "counts": {
            "by_primary_reason": dict(sorted(reason_counts.items())),
            "by_recommended_action": dict(sorted(action_counts.items())),
            "valid_models_found_during_review": len(valid_found_rows),
            "true_structure_gaps": len(true_gap_rows),
            "structure_gaps_for_baldo": len(baldo_rows),
            "manual_review_rows": len(manual_rows),
        },
        "controls": controls,
        "checks": checks,
        "database_before": database_before,
        "database_after": database_after,
        "protected_files_before": protected_before,
        "protected_files_after": protected_after,
        "protected_directories_before": protected_dirs_before,
        "protected_directories_after": protected_dirs_after,
        "input_hashes": {str(path): file_sha256(path) for path in inputs if path.is_file()},
        "outputs": {
            key: {
                "path": str(path),
                "rows": {
                    "adjudication": len(output_rows), "true_gaps": len(true_gap_rows),
                    "valid_found": len(valid_found_rows), "baldo": len(baldo_rows),
                    "manual": len(manual_rows),
                }.get(key),
                "sha256": file_sha256(path) if key != "qc" else None,
            }
            for key, path in outputs.items()
        },
    }
    with outputs["qc"].open("x", encoding="utf-8") as handle:
        json.dump(qc, handle, indent=2, sort_keys=True)
        handle.write("\n")

    print(json.dumps({
        "status": "PASS",
        "unique_accessions_reviewed": len(population),
        "overlap": overlap,
        "by_primary_reason": dict(sorted(reason_counts.items())),
        "by_recommended_action": dict(sorted(action_counts.items())),
        "true_structure_gaps": len(true_gap_rows),
        "controls": {key: item["pass"] for key, item in controls.items()},
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
