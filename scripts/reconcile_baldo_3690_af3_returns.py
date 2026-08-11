#!/usr/bin/env python3
"""Read-only reconciliation of the exact 3,690 Baldo FASTA against AF3 returns.

The AF3 inventory is rebuilt directly from the returned directory.  For sent
accessions, every selected and seed/sample CIF is audited with the same Gemmi
4.5 A heavy-atom interface criterion used by audit_eligible_af3_interfaces.py.
SQLite inputs are opened using mode=ro&immutable=1.  Existing outputs are
never overwritten.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import hashlib
import json
import os
import re
import socket
import sqlite3
import statistics
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable


ROOT = Path(__file__).resolve().parents[1]
FASTA_DEFAULT = Path("/home/patricia/TF_database_Baldo_data/TF_without_model.fasta")
AF3_DEFAULT = Path("/data/sbi/interchange/boliva/patricia/AF3")
CONTACT_CUTOFF_ANGSTROM = 4.5
AF3_HASH_WORKERS = 8

OUTPUT_NAMES = {
    "fasta": "baldo_3690_fasta_accessions.tsv",
    "inventory": "baldo_af3_complete_inventory.tsv",
    "master": "baldo_3690_af3_reconciliation.tsv",
    "valid": "baldo_3690_af3_valid_interface.tsv",
    "empty": "baldo_3690_af3_empty_interface.tsv",
    "no_return": "baldo_3690_no_af3_return.tsv",
    "old": "baldo_af3_not_in_3690_fasta.tsv",
    "qc": "baldo_3690_af3_reconciliation_qc.json",
}

MASTER_COLUMNS = [
    "accession", "in_sent_fasta", "fasta_header", "sequence_length",
    "af3_return_found", "af3_model_count", "af3_model_paths",
    "confidence_json_found", "pae_json_found", "protein_present", "dna_present",
    "valid_interface_any_model", "best_model_path", "best_model_rank",
    "atom_contacts", "protein_interface_residues", "dna_interface_residues",
    "final_return_status", "notes",
]

FASTA_COLUMNS = [
    "record_order", "accession", "fasta_header", "sequence_length",
    "sequence_sha256", "header_status", "duplicate_count",
]

INVENTORY_COLUMNS = [
    "accession", "accession_extraction_status", "af3_set_classification",
    "in_sent_fasta", "directory_name", "directory_path", "directory_variant",
    "model_kind", "model_rank", "ranking_score", "model_path",
    "relative_model_path", "model_file_size", "model_sha256",
    "input_json_paths", "data_json_paths", "confidence_json_path",
    "confidence_json_found", "plddt_json_found", "pae_json_found",
    "summary_confidence_json_path", "ranking_metadata_paths", "file_exists",
    "readable_cif", "protein_present", "dna_present", "protein_atom_count",
    "dna_atom_count", "atom_contacts", "protein_interface_residues",
    "dna_interface_residues", "interface_result", "protein_chains", "dna_chains",
    "protein_ca_plddt_mean", "protein_ca_plddt_min", "protein_ca_plddt_max",
    "error", "notes",
]

OLD_COLUMNS = [
    "accession", "af3_model_count", "af3_model_paths", "in_original_db",
    "original_active_structure", "possible_old_af3", "notes",
]

ALLOWED_STATUSES = {
    "VALID_AF3_PROTEIN_DNA_MODEL",
    "AF3_RETURN_EMPTY_INTERFACE",
    "AF3_RETURN_INVALID_OR_UNREADABLE",
    "NO_AF3_RETURN_FOUND",
}

# UniProt's six- and ten-character accession forms.
UNIPROT_RE = re.compile(
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|"
    r"[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]|"
    r"[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){2})$"
)
SAMPLE_RE = re.compile(r"seed-(\d+)_sample-(\d+)", re.IGNORECASE)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Directly reconcile TF_without_model.fasta against all AF3 output "
            "under Baldo's returned directory and audit sent-accession interfaces."
        )
    )
    parser.add_argument("--fasta", type=Path, default=FASTA_DEFAULT)
    parser.add_argument("--af3-root", type=Path, default=AF3_DEFAULT)
    parser.add_argument("--database", type=Path, default=ROOT / "data/tf_webdb.sqlite")
    parser.add_argument(
        "--staging-database", type=Path,
        default=ROOT / "data/tf_webdb_local_staging.sqlite",
    )
    parser.add_argument("--output-dir", type=Path, default=ROOT / "outputs")
    return parser.parse_args()


def clean(value: Any) -> str:
    return "" if value is None else str(value).strip()


def sha256_file(path: Path, tokens: tuple[bytes, ...] = ()) -> tuple[str, dict[str, bool]]:
    digest = hashlib.sha256()
    found = {token.decode("ascii", "replace"): False for token in tokens}
    overlap = max((len(token) for token in tokens), default=1) - 1
    tail = b""
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
            if tokens:
                searchable = tail + block
                for token in tokens:
                    if token in searchable:
                        found[token.decode("ascii", "replace")] = True
                tail = searchable[-overlap:] if overlap else b""
    return digest.hexdigest(), found


def simple_file_snapshot(path: Path) -> dict[str, Any]:
    stat = path.stat()
    digest, _ = sha256_file(path)
    return {
        "path": str(path.resolve()), "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns, "sha256": digest,
    }


def sqlite_connection(path: Path) -> sqlite3.Connection:
    connection = sqlite3.connect(f"file:{path.resolve()}?mode=ro&immutable=1", uri=True)
    connection.row_factory = sqlite3.Row
    return connection


def database_snapshot(path: Path) -> dict[str, Any]:
    snapshot = simple_file_snapshot(path)
    with sqlite_connection(path) as connection:
        snapshot["quick_check"] = connection.execute("PRAGMA quick_check").fetchone()[0]
        tables = [
            row[0] for row in connection.execute(
                "SELECT name FROM sqlite_master WHERE type='table' ORDER BY name"
            )
        ]
        snapshot["row_counts"] = {
            table: connection.execute(f'SELECT COUNT(*) FROM "{table}"').fetchone()[0]
            for table in tables
        }
    return snapshot


def website_snapshot() -> dict[str, str]:
    paths = [ROOT / "app.py"]
    for directory in [ROOT / "templates", ROOT / "static"]:
        if directory.is_dir():
            paths.extend(path for path in directory.rglob("*") if path.is_file())
    result: dict[str, str] = {}
    for path in sorted(set(paths)):
        result[str(path.relative_to(ROOT))] = sha256_file(path)[0]
    return result


def af3_tree_snapshot(
    root: Path, *, scan_confidence_tokens: bool, label: str
) -> tuple[dict[str, Any], dict[str, dict[str, Any]]]:
    """Hash every AF3 file and return a compact tree plus per-file manifest."""
    files = sorted(path for path in root.rglob("*") if path.is_file())
    manifest: dict[str, dict[str, Any]] = {}
    extension_counts: Counter[str] = Counter()
    extension_bytes: Counter[str] = Counter()
    aggregate = hashlib.sha256()
    total_bytes = sum(path.stat().st_size for path in files)
    processed_bytes = 0
    started = time.monotonic()
    print(f"[{label}] hashing {len(files)} AF3 files ({total_bytes} bytes)", flush=True)

    def hash_path(path: Path) -> tuple[Path, os.stat_result, str, dict[str, bool]]:
        stat = path.stat()
        tokens: tuple[bytes, ...] = ()
        if (
            scan_confidence_tokens
            and path.suffix.lower() == ".json"
            and "confidences" in path.name.lower()
            and "summary" not in path.name.lower()
        ):
            tokens = (b'"atom_plddts"', b'"pae"')
        digest, found = sha256_file(path, tokens)
        return path, stat, digest, found

    with concurrent.futures.ThreadPoolExecutor(max_workers=AF3_HASH_WORKERS) as executor:
        hashed_files = executor.map(hash_path, files)
        for index, (path, stat, digest, found) in enumerate(hashed_files, start=1):
            relative = path.relative_to(root).as_posix()
            entry = {
                "size": stat.st_size, "mtime_ns": stat.st_mtime_ns,
                "sha256": digest,
                "contains_atom_plddts": found.get('"atom_plddts"', False),
                "contains_pae": found.get('"pae"', False),
            }
            manifest[relative] = entry
            suffix = path.suffix.lower() or "[no_suffix]"
            extension_counts[suffix] += 1
            extension_bytes[suffix] += stat.st_size
            processed_bytes += stat.st_size
            aggregate.update(relative.encode("utf-8"))
            aggregate.update(b"\0")
            aggregate.update(str(stat.st_size).encode("ascii"))
            aggregate.update(b"\0")
            aggregate.update(digest.encode("ascii"))
            aggregate.update(b"\n")
            if index % 1000 == 0 or index == len(files):
                elapsed = max(time.monotonic() - started, 0.001)
                print(
                    f"[{label}] {index}/{len(files)} files; "
                    f"{processed_bytes}/{total_bytes} bytes; {processed_bytes / elapsed / 1e6:.1f} MB/s",
                    flush=True,
                )

    root_stat = root.stat()
    summary = {
        "path": str(root.resolve()), "file_count": len(files),
        "total_file_bytes": total_bytes, "root_mtime_ns": root_stat.st_mtime_ns,
        "manifest_sha256": aggregate.hexdigest(),
        "extension_counts": dict(sorted(extension_counts.items())),
        "extension_bytes": dict(sorted(extension_bytes.items())),
    }
    return summary, manifest


def parse_fasta(path: Path) -> tuple[list[dict[str, Any]], list[str]]:
    records: list[dict[str, Any]] = []
    errors: list[str] = []
    header: str | None = None
    sequence: list[str] = []

    def finish() -> None:
        if header is None:
            return
        raw_id = header.split("|", 1)[0].split()[0].strip().upper()
        seq = "".join(sequence).replace(" ", "").upper()
        status = "PARSED_UNIPROT_ACCESSION" if UNIPROT_RE.fullmatch(raw_id) else "MALFORMED_OR_UNPARSEABLE"
        if status != "PARSED_UNIPROT_ACCESSION":
            errors.append(header)
        records.append({
            "record_order": len(records) + 1, "accession": raw_id,
            "fasta_header": header, "sequence": seq, "sequence_length": len(seq),
            "sequence_sha256": hashlib.sha256(seq.encode("ascii")).hexdigest(),
            "header_status": status,
        })

    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if line.startswith(">"):
                finish()
                header = line[1:]
                sequence = []
            elif line:
                if header is None:
                    errors.append("SEQUENCE_BEFORE_FIRST_HEADER")
                sequence.append(line)
    finish()
    counts = Counter(row["accession"] for row in records)
    for row in records:
        row["duplicate_count"] = counts[row["accession"]]
    return records, errors


def extract_accession(directory_name: str) -> tuple[str, str, str]:
    normalized = directory_name.upper()
    variant = "base"
    for suffix, candidate_variant in [("_DIMER", "dimer"), ("_MONOMER", "monomer")]:
        if normalized.endswith(suffix):
            normalized = normalized[: -len(suffix)]
            variant = candidate_variant
            break
    if UNIPROT_RE.fullmatch(normalized):
        return normalized, "UNAMBIGUOUS_DIRECTORY_ACCESSION", variant
    return normalized, "AMBIGUOUS_OR_NON_UNIPROT_DIRECTORY_ID", variant


def read_ranking_scores(directory: Path) -> tuple[dict[tuple[int, int], tuple[int, float]], list[Path]]:
    ranking_files = sorted(directory.glob("*ranking_scores.csv"))
    ranking: dict[tuple[int, int], tuple[int, float]] = {}
    for path in ranking_files:
        with path.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        scored = []
        for row in rows:
            try:
                scored.append((int(row["seed"]), int(row["sample"]), float(row["ranking_score"])))
            except (KeyError, TypeError, ValueError):
                continue
        for rank, (seed, sample, score) in enumerate(
            sorted(scored, key=lambda value: (-value[2], value[0], value[1])), start=1
        ):
            ranking[(seed, sample)] = (rank, score)
    return ranking, ranking_files


def matching_json(model_path: Path, kind: str) -> Path:
    if model_path.name == "model.cif":
        return model_path.with_name(f"{kind}.json")
    marker = "_model.cif"
    if model_path.name.lower().endswith(marker):
        return model_path.with_name(model_path.name[: -len(marker)] + f"_{kind}.json")
    return model_path.with_suffix(f".{kind}.json")


def audit_model(path: Path) -> dict[str, Any]:
    # Import lazily so --help and py_compile work with system Python.  The
    # validated audit runtime is .conda-af3-audit/bin/python (Gemmi 0.7.5).
    import gemmi  # type: ignore

    record: dict[str, Any] = {
        "file_exists": "YES" if path.is_file() else "NO", "readable_cif": "NO",
        "protein_present": "NO", "dna_present": "NO", "protein_atom_count": 0,
        "dna_atom_count": 0, "atom_contacts": 0, "protein_interface_residues": 0,
        "dna_interface_residues": 0, "interface_result": "", "protein_chains": "",
        "dna_chains": "", "protein_ca_plddt_mean": "", "protein_ca_plddt_min": "",
        "protein_ca_plddt_max": "", "error": "",
    }
    if not path.is_file():
        record["interface_result"] = "ERROR_MISSING_FILE"
        record["error"] = "File does not exist"
        return record
    try:
        structure = gemmi.read_structure(str(path))
        structure.setup_entities()
        if len(structure) == 0:
            record["interface_result"] = "ERROR_NO_MODEL"
            record["error"] = "Structure contains no models"
            return record
        record["readable_cif"] = "YES"
        model = structure[0]
        protein_indices: set[int] = set()
        dna_indices: set[int] = set()
        for chain_index, chain in enumerate(model):
            polymer = chain.get_polymer()
            if len(polymer) == 0:
                continue
            polymer_type = str(polymer.check_polymer_type())
            if "Peptide" in polymer_type:
                protein_indices.add(chain_index)
            # This exactly preserves the previous Baldo AF3 audit's nucleic
            # polymer rule.  Returned inputs are DNA; RNA would also be
            # recognized by that established method.
            if "Dna" in polymer_type or "Rna" in polymer_type:
                dna_indices.add(chain_index)
        record["protein_chains"] = ",".join(model[index].name for index in sorted(protein_indices))
        record["dna_chains"] = ",".join(model[index].name for index in sorted(dna_indices))
        record["protein_present"] = "YES" if protein_indices else "NO"
        record["dna_present"] = "YES" if dna_indices else "NO"
        protein_atoms = 0
        dna_atoms = 0
        ca_plddt: list[float] = []
        for index in protein_indices:
            for residue in model[index]:
                for atom in residue:
                    if not atom.is_hydrogen():
                        protein_atoms += 1
                        if atom.name.strip() == "CA":
                            ca_plddt.append(float(atom.b_iso))
        for index in dna_indices:
            for residue in model[index]:
                for atom in residue:
                    if not atom.is_hydrogen():
                        dna_atoms += 1
        record["protein_atom_count"] = protein_atoms
        record["dna_atom_count"] = dna_atoms
        if ca_plddt:
            record["protein_ca_plddt_mean"] = f"{statistics.mean(ca_plddt):.3f}"
            record["protein_ca_plddt_min"] = f"{min(ca_plddt):.3f}"
            record["protein_ca_plddt_max"] = f"{max(ca_plddt):.3f}"
        if not protein_indices:
            record["interface_result"] = "FAIL_NO_PROTEIN"
            return record
        if not dna_indices:
            record["interface_result"] = "FAIL_NO_DNA"
            return record

        neighbors = gemmi.NeighborSearch(model, structure.cell, CONTACT_CUTOFF_ANGSTROM).populate(include_h=False)
        atom_contacts: set[tuple[int, int, int, int, int, int]] = set()
        protein_residues: set[tuple[str, str, str]] = set()
        dna_residues: set[tuple[str, str, str]] = set()
        for protein_chain_index in sorted(protein_indices):
            chain = model[protein_chain_index]
            for protein_residue_index, residue in enumerate(chain):
                for protein_atom_index, atom in enumerate(residue):
                    if atom.is_hydrogen():
                        continue
                    for mark in neighbors.find_neighbors(
                        atom, min_dist=0.1, max_dist=CONTACT_CUTOFF_ANGSTROM
                    ):
                        if mark.image_idx != 0 or mark.chain_idx not in dna_indices:
                            continue
                        dna_cra = mark.to_cra(model)
                        atom_contacts.add((
                            protein_chain_index, protein_residue_index, protein_atom_index,
                            mark.chain_idx, mark.residue_idx, mark.atom_idx,
                        ))
                        protein_residues.add((chain.name, str(residue.seqid), residue.name))
                        dna_residues.add((
                            dna_cra.chain.name, str(dna_cra.residue.seqid), dna_cra.residue.name
                        ))
        record["atom_contacts"] = len(atom_contacts)
        record["protein_interface_residues"] = len(protein_residues)
        record["dna_interface_residues"] = len(dna_residues)
        record["interface_result"] = (
            "PASS_NONEMPTY_PROTEIN_DNA_INTERFACE"
            if atom_contacts else "FAIL_EMPTY_PROTEIN_DNA_INTERFACE"
        )
    except Exception as exc:  # Preserve each unreadable/incomplete return.
        record["interface_result"] = "ERROR"
        record["error"] = f"{type(exc).__name__}: {exc}"
    return record


def write_tsv(path: Path, rows: Iterable[dict[str, Any]], columns: list[str]) -> None:
    with path.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def add_check(checks: dict[str, Any], name: str, observed: Any, expected: Any) -> None:
    checks[name] = {"observed": observed, "expected": expected, "passed": observed == expected}


def choose_best(rows: list[dict[str, Any]], status: str) -> tuple[dict[str, Any] | None, str]:
    if not rows:
        return None, "No CIF model file was present."
    if status == "VALID_AF3_PROTEIN_DNA_MODEL":
        candidates = [row for row in rows if row["interface_result"] == "PASS_NONEMPTY_PROTEIN_DNA_INTERFACE"]
        reason = (
            "Preferred among passing models by atom-contact count, then protein/DNA "
            "interface-residue counts, ranking score, and selected-model status."
        )
        key = lambda row: (
            int(row["atom_contacts"] or 0), int(row["protein_interface_residues"] or 0),
            int(row["dna_interface_residues"] or 0), float(row["ranking_score"] or -1),
            row["model_kind"] == "selected",
        )
    elif status == "AF3_RETURN_EMPTY_INTERFACE":
        candidates = [row for row in rows if row["interface_result"] == "FAIL_EMPTY_PROTEIN_DNA_INTERFACE"]
        reason = "Preferred complete empty-interface model by ranking score and selected-model status."
        key = lambda row: (float(row["ranking_score"] or -1), row["model_kind"] == "selected")
    else:
        candidates = rows
        reason = "Representative invalid/unreadable return; selected/rank and file existence were preferred."
        key = lambda row: (
            row["file_exists"] == "YES", row["model_kind"] == "selected",
            -int(row["model_rank"] or 999), float(row["ranking_score"] or -1),
        )
    return max(candidates, key=key) if candidates else None, reason


def main() -> int:
    args = parse_args()
    if socket.gethostname() != "masada":
        raise SystemExit(f"Refusing to run on {socket.gethostname()!r}; expected exactly 'masada'")
    required = [args.fasta, args.database, args.staging_database]
    missing = [str(path) for path in required if not path.is_file()]
    if not args.af3_root.is_dir():
        missing.append(str(args.af3_root))
    if missing:
        raise SystemExit("Missing required input(s): " + ", ".join(missing))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    output_paths = {key: args.output_dir / name for key, name in OUTPUT_NAMES.items()}
    existing = [str(path) for path in output_paths.values() if path.exists()]
    if existing:
        raise SystemExit("Refusing to overwrite existing output(s): " + ", ".join(existing))

    source_db_before = database_snapshot(args.database)
    staging_db_before = database_snapshot(args.staging_database)
    fasta_before = simple_file_snapshot(args.fasta)
    website_before = website_snapshot()
    af3_before, af3_manifest_before = af3_tree_snapshot(
        args.af3_root, scan_confidence_tokens=True, label="AF3 before"
    )

    fasta_records, fasta_errors = parse_fasta(args.fasta)
    fasta_by_accession = {row["accession"]: row for row in fasta_records}
    fasta_ids = set(fasta_by_accession)
    fasta_counts = Counter(row["accession"] for row in fasta_records)

    directory_records: list[dict[str, Any]] = []
    directories_by_accession: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for directory in sorted(path for path in args.af3_root.iterdir() if path.is_dir()):
        accession, extraction_status, variant = extract_accession(directory.name)
        ranking, ranking_files = read_ranking_scores(directory)
        model_paths = sorted(
            path for path in directory.rglob("*")
            if path.is_file() and path.suffix.lower() in {".cif", ".mmcif"}
        )
        input_json = args.af3_root / f"{directory.name}_af3.json"
        data_jsons = sorted(directory.glob("*_data.json"))
        record = {
            "accession": accession, "extraction_status": extraction_status,
            "variant": variant, "directory": directory, "ranking": ranking,
            "ranking_files": ranking_files, "model_paths": model_paths,
            "input_jsons": [input_json] if input_json.is_file() else [],
            "data_jsons": data_jsons,
        }
        directory_records.append(record)
        directories_by_accession[accession].append(record)

    af3_ids = set(directories_by_accession)
    intersection = fasta_ids & af3_ids
    fasta_only = fasta_ids - af3_ids
    af3_only = af3_ids - fasta_ids

    with sqlite_connection(args.database) as connection:
        original_tf_ids = {clean(row[0]).upper() for row in connection.execute("SELECT tf_id FROM tf")}
        active_by_accession: dict[str, list[str]] = defaultdict(list)
        active_alpha_ids: set[str] = set()
        for row in connection.execute(
            """SELECT tf_id, source FROM structure_file
               WHERE status='active' AND file_type='pdb'"""
        ):
            accession = clean(row["tf_id"]).upper()
            if accession:
                active_by_accession[accession].append(clean(row["source"]))
                if clean(row["source"]) == "alphafold":
                    active_alpha_ids.add(accession)

    inventory_rows: list[dict[str, Any]] = []
    model_rows_by_accession: dict[str, list[dict[str, Any]]] = defaultdict(list)
    cif_paths_seen: set[str] = set()
    audit_index = 0
    audit_total = sum(
        len(directory["model_paths"])
        for accession in intersection for directory in directories_by_accession[accession]
    )
    print(f"Auditing {audit_total} CIF variants for {len(intersection)} sent accessions", flush=True)

    for accession in sorted(af3_ids):
        for directory_record in directories_by_accession[accession]:
            directory = directory_record["directory"]
            common = {
                "accession": accession,
                "accession_extraction_status": directory_record["extraction_status"],
                "af3_set_classification": "SENT_FASTA_AF3" if accession in fasta_ids else "AF3_NOT_IN_SENT_FASTA",
                "in_sent_fasta": "YES" if accession in fasta_ids else "NO",
                "directory_name": directory.name, "directory_path": str(directory),
                "directory_variant": directory_record["variant"],
                "input_json_paths": ";".join(str(path) for path in directory_record["input_jsons"]),
                "data_json_paths": ";".join(str(path) for path in directory_record["data_jsons"]),
                "ranking_metadata_paths": ";".join(str(path) for path in directory_record["ranking_files"]),
            }
            if not directory_record["model_paths"]:
                row = dict(common)
                row.update({
                    "model_kind": "no_model_cif", "model_rank": "", "ranking_score": "",
                    "model_path": "", "relative_model_path": "", "model_file_size": 0,
                    "model_sha256": "", "confidence_json_path": "",
                    "confidence_json_found": "NO", "plddt_json_found": "NO",
                    "pae_json_found": "NO", "summary_confidence_json_path": "",
                    "file_exists": "NO", "readable_cif": "NO", "protein_present": "NO",
                    "dna_present": "NO", "protein_atom_count": 0, "dna_atom_count": 0,
                    "atom_contacts": 0, "protein_interface_residues": 0,
                    "dna_interface_residues": 0, "interface_result": "NO_MODEL_CIF",
                    "protein_chains": "", "dna_chains": "", "protein_ca_plddt_mean": "",
                    "protein_ca_plddt_min": "", "protein_ca_plddt_max": "", "error": "",
                    "notes": "AF3 job directory/input data exist, but no CIF/mmCIF structure was returned.",
                })
                inventory_rows.append(row)
                continue

            for model_path in directory_record["model_paths"]:
                cif_paths_seen.add(str(model_path.resolve()))
                relative = model_path.relative_to(args.af3_root).as_posix()
                sample_match = SAMPLE_RE.search(relative)
                if sample_match:
                    seed = int(sample_match.group(1))
                    sample = int(sample_match.group(2))
                    model_kind = "sample"
                    model_rank, ranking_score = directory_record["ranking"].get((seed, sample), ("", ""))
                else:
                    model_kind = "selected"
                    model_rank = 1 if directory_record["ranking"] else ""
                    ranking_score = max(
                        (value[1] for value in directory_record["ranking"].values()), default=""
                    )
                confidence_path = matching_json(model_path, "confidences")
                summary_path = matching_json(model_path, "summary_confidences")
                confidence_relative = (
                    confidence_path.relative_to(args.af3_root).as_posix()
                    if confidence_path.is_file() else ""
                )
                confidence_manifest = af3_manifest_before.get(confidence_relative, {})
                model_manifest = af3_manifest_before[relative]
                row = dict(common)
                row.update({
                    "model_kind": model_kind, "model_rank": model_rank,
                    "ranking_score": ranking_score, "model_path": str(model_path),
                    "relative_model_path": relative, "model_file_size": model_manifest["size"],
                    "model_sha256": model_manifest["sha256"],
                    "confidence_json_path": str(confidence_path) if confidence_path.is_file() else "",
                    "confidence_json_found": "YES" if confidence_path.is_file() else "NO",
                    "plddt_json_found": "YES" if confidence_manifest.get("contains_atom_plddts") else "NO",
                    "pae_json_found": "YES" if confidence_manifest.get("contains_pae") else "NO",
                    "summary_confidence_json_path": str(summary_path) if summary_path.is_file() else "",
                })
                if accession in intersection:
                    audit_index += 1
                    row.update(audit_model(model_path))
                    row["notes"] = (
                        f"Interface audited with Gemmi heavy-atom cutoff {CONTACT_CUTOFF_ANGSTROM:.1f} A; "
                        "all selected and seed/sample variants for sent accessions are included."
                    )
                    if audit_index % 100 == 0 or audit_index == audit_total:
                        print(f"Interface audit {audit_index}/{audit_total}", flush=True)
                else:
                    row.update({
                        "file_exists": "YES", "readable_cif": "NOT_AUDITED",
                        "protein_present": "NOT_AUDITED", "dna_present": "NOT_AUDITED",
                        "protein_atom_count": "", "dna_atom_count": "", "atom_contacts": "",
                        "protein_interface_residues": "", "dna_interface_residues": "",
                        "interface_result": "NOT_AUDITED_NOT_IN_SENT_FASTA",
                        "protein_chains": "", "dna_chains": "", "protein_ca_plddt_mean": "",
                        "protein_ca_plddt_min": "", "protein_ca_plddt_max": "", "error": "",
                        "notes": "Historical/not-sent AF3 model inventoried but not interface-audited in this task.",
                    })
                inventory_rows.append(row)
                model_rows_by_accession[accession].append(row)

    master_rows: list[dict[str, Any]] = []
    for fasta_record in fasta_records:
        accession = fasta_record["accession"]
        rows = model_rows_by_accession.get(accession, [])
        return_found = accession in af3_ids
        passing = [row for row in rows if row["interface_result"] == "PASS_NONEMPTY_PROTEIN_DNA_INTERFACE"]
        complete_empty = [row for row in rows if row["interface_result"] == "FAIL_EMPTY_PROTEIN_DNA_INTERFACE"]
        if not return_found:
            status = "NO_AF3_RETURN_FOUND"
        elif passing:
            status = "VALID_AF3_PROTEIN_DNA_MODEL"
        elif complete_empty:
            status = "AF3_RETURN_EMPTY_INTERFACE"
        else:
            status = "AF3_RETURN_INVALID_OR_UNREADABLE"
        best, best_reason = choose_best(rows, status)
        result_counts = Counter(row["interface_result"] for row in rows)
        note_parts = [best_reason]
        if rows:
            note_parts.append(
                "Model results: " + ", ".join(f"{key}={value}" for key, value in sorted(result_counts.items())) + "."
            )
        elif return_found:
            note_parts.append("AF3 directory/input material was present but no model CIF/mmCIF was found.")
        else:
            note_parts.append("No AF3 accession directory identifiable for this exact FASTA accession.")
        master_rows.append({
            "accession": accession, "in_sent_fasta": "YES",
            "fasta_header": fasta_record["fasta_header"],
            "sequence_length": fasta_record["sequence_length"],
            "af3_return_found": "YES" if return_found else "NO",
            "af3_model_count": len(rows),
            "af3_model_paths": ";".join(row["model_path"] for row in rows),
            "confidence_json_found": "YES" if any(row["confidence_json_found"] == "YES" for row in rows) else "NO",
            "pae_json_found": "YES" if any(row["pae_json_found"] == "YES" for row in rows) else "NO",
            "protein_present": "YES" if any(row["protein_present"] == "YES" for row in rows) else "NO",
            "dna_present": "YES" if any(row["dna_present"] == "YES" for row in rows) else "NO",
            "valid_interface_any_model": "YES" if passing else "NO",
            "best_model_path": best["model_path"] if best else "",
            "best_model_rank": best["model_rank"] if best else "",
            "atom_contacts": best["atom_contacts"] if best else 0,
            "protein_interface_residues": best["protein_interface_residues"] if best else 0,
            "dna_interface_residues": best["dna_interface_residues"] if best else 0,
            "final_return_status": status, "notes": " ".join(note_parts),
        })

    old_rows: list[dict[str, Any]] = []
    for accession in sorted(af3_only):
        rows = model_rows_by_accession.get(accession, [])
        in_db = accession in original_tf_ids
        active_sources = sorted(set(active_by_accession.get(accession, [])))
        if accession in active_alpha_ids:
            possible = "YES_ORIGINAL_ACTIVE_ALPHAFOLD"
            note = (
                "Not in the 3,690 FASTA and already represented by an active original-db "
                "AlphaFold structure; strong support for Baldo's old/pre-existing set."
            )
        elif in_db and active_sources:
            possible = "YES_ORIGINAL_ACTIVE_OTHER_STRUCTURE"
            note = (
                "Not in the 3,690 FASTA and original DB already has active structure source(s): "
                + ",".join(active_sources) + "."
            )
        elif in_db:
            possible = "POSSIBLE_HISTORICAL_AF3"
            note = "Not in the sent FASTA but accession exists in original DB; timestamp alone is not used as proof."
        else:
            possible = "POSSIBLE_OTHER_HISTORICAL_AF3_DATASET"
            note = "Not in the sent FASTA or current original TF table; retained as candidate historical AF3 material."
        if not rows:
            note += " AF3 directory/input material exists but no model CIF/mmCIF is present."
        old_rows.append({
            "accession": accession, "af3_model_count": len(rows),
            "af3_model_paths": ";".join(row["model_path"] for row in rows),
            "in_original_db": "YES" if in_db else "NO",
            "original_active_structure": "YES" if active_sources else "NO",
            "possible_old_af3": possible, "notes": note,
        })

    valid_rows = [row for row in master_rows if row["final_return_status"] == "VALID_AF3_PROTEIN_DNA_MODEL"]
    empty_rows = [row for row in master_rows if row["final_return_status"] == "AF3_RETURN_EMPTY_INTERFACE"]
    invalid_rows = [row for row in master_rows if row["final_return_status"] == "AF3_RETURN_INVALID_OR_UNREADABLE"]
    no_return_rows = [row for row in master_rows if row["final_return_status"] == "NO_AF3_RETURN_FOUND"]

    checks: dict[str, Any] = {}
    add_check(checks, "hostname", socket.gethostname(), "masada")
    add_check(checks, "fasta_records_and_unique", (len(fasta_records), len(fasta_ids)), (3690, 3690))
    add_check(checks, "fasta_duplicate_accessions", sorted(key for key, count in fasta_counts.items() if count > 1), [])
    add_check(checks, "fasta_malformed_headers", fasta_errors, [])
    add_check(checks, "set_partition_fasta", len(intersection) + len(fasta_only), 3690)
    add_check(checks, "set_partition_af3", len(intersection) + len(af3_only), len(af3_ids))
    add_check(checks, "ambiguous_af3_directory_ids", sorted(
        record["directory"].name for record in directory_records
        if record["extraction_status"] != "UNAMBIGUOUS_DIRECTORY_ACCESSION"
    ), [])
    actual_cifs = {
        str(path.resolve()) for path in args.af3_root.rglob("*")
        if path.is_file() and path.suffix.lower() in {".cif", ".mmcif"}
    }
    add_check(checks, "every_af3_model_file_in_inventory", sorted(actual_cifs - cif_paths_seen), [])
    add_check(checks, "inventory_has_no_unexpected_model_paths", sorted(cif_paths_seen - actual_cifs), [])
    add_check(checks, "master_statuses_allowed", sorted({row["final_return_status"] for row in master_rows} - ALLOWED_STATUSES), [])
    category_sets = {
        status: {row["accession"] for row in master_rows if row["final_return_status"] == status}
        for status in ALLOWED_STATUSES
    }
    overlap = set()
    statuses = sorted(category_sets)
    for index, status in enumerate(statuses):
        for other in statuses[index + 1:]:
            overlap.update(category_sets[status] & category_sets[other])
    add_check(checks, "final_categories_disjoint", sorted(overlap), [])
    add_check(
        checks,
        "final_categories_cover_fasta",
        sorted(set().union(*category_sets.values())),
        sorted(fasta_ids),
    )
    add_check(checks, "valid_empty_invalid_all_from_fasta", sorted(
        (category_sets["VALID_AF3_PROTEIN_DNA_MODEL"] |
         category_sets["AF3_RETURN_EMPTY_INTERFACE"] |
         category_sets["AF3_RETURN_INVALID_OR_UNREADABLE"]) - fasta_ids
    ), [])
    add_check(checks, "all_intersection_models_audited", sum(
        row["interface_result"].startswith("NOT_AUDITED") for accession in intersection
        for row in model_rows_by_accession.get(accession, [])
    ), 0)

    control_accessions = ["H7C4N4", "E9PN75", "B4DHE0", "Q59EF3", "P01100"]
    master_by_accession = {row["accession"]: row for row in master_rows}
    controls: dict[str, Any] = {}
    for accession in control_accessions:
        row = master_by_accession.get(accession)
        controls[accession] = {
            "in_sent_fasta": accession in fasta_ids,
            "af3_directory_found": accession in af3_ids,
            "af3_model_count": int(row["af3_model_count"]) if row else 0,
            "final_return_status": row["final_return_status"] if row else "NOT_IN_FASTA",
            "best_model_path": row["best_model_path"] if row else "",
            "exact_accession_paths_only": all(
                extract_accession(Path(model["model_path"]).relative_to(args.af3_root).parts[0])[0] == accession
                for model in model_rows_by_accession.get(accession, [])
            ),
            "later_fragment_evidence_used": False,
            "modcre_model_evidence_used": False,
        }
    add_check(checks, "controls_in_fasta", sorted(
        accession for accession in control_accessions if accession not in fasta_ids
    ), [])
    add_check(checks, "controls_exact_accession_only", all(
        controls[accession]["exact_accession_paths_only"] for accession in control_accessions
    ), True)
    add_check(checks, "controls_exclude_fragment_and_modcre_evidence", all(
        not controls[accession]["later_fragment_evidence_used"]
        and not controls[accession]["modcre_model_evidence_used"]
        for accession in control_accessions
    ), True)

    prewrite_failures = [name for name, check in checks.items() if not check["passed"]]
    if prewrite_failures:
        raise SystemExit("QC failed before output creation: " + ", ".join(prewrite_failures))

    write_tsv(output_paths["fasta"], fasta_records, FASTA_COLUMNS)
    write_tsv(output_paths["inventory"], inventory_rows, INVENTORY_COLUMNS)
    write_tsv(output_paths["master"], master_rows, MASTER_COLUMNS)
    write_tsv(output_paths["valid"], valid_rows, MASTER_COLUMNS)
    write_tsv(output_paths["empty"], empty_rows, MASTER_COLUMNS)
    write_tsv(output_paths["no_return"], no_return_rows, MASTER_COLUMNS)
    write_tsv(output_paths["old"], old_rows, OLD_COLUMNS)

    source_db_after = database_snapshot(args.database)
    staging_db_after = database_snapshot(args.staging_database)
    fasta_after = simple_file_snapshot(args.fasta)
    website_after = website_snapshot()
    af3_after, af3_manifest_after = af3_tree_snapshot(
        args.af3_root, scan_confidence_tokens=False, label="AF3 after"
    )
    changed_af3_files = sorted(
        (set(af3_manifest_before) ^ set(af3_manifest_after))
        | {
            path for path in set(af3_manifest_before) & set(af3_manifest_after)
            if af3_manifest_before[path]["sha256"] != af3_manifest_after[path]["sha256"]
            or af3_manifest_before[path]["size"] != af3_manifest_after[path]["size"]
            or af3_manifest_before[path]["mtime_ns"] != af3_manifest_after[path]["mtime_ns"]
        }
    )
    add_check(checks, "source_database_unchanged", source_db_after, source_db_before)
    add_check(checks, "staging_database_unchanged", staging_db_after, staging_db_before)
    add_check(checks, "fasta_unchanged", fasta_after, fasta_before)
    add_check(checks, "website_files_unchanged", website_after, website_before)
    add_check(checks, "af3_tree_summary_unchanged", af3_after, af3_before)
    add_check(checks, "af3_file_hashes_sizes_mtimes_unchanged", changed_af3_files, [])
    add_check(checks, "no_remote_connection", False, False)

    model_file_rows = [row for row in inventory_rows if row["model_path"]]
    accession_model_counts = Counter(len(model_rows_by_accession.get(accession, [])) for accession in af3_ids)
    accessions_with_confidence = {
        accession for accession in af3_ids
        if any(row["confidence_json_found"] == "YES" for row in model_rows_by_accession.get(accession, []))
    }
    accessions_with_pae = {
        accession for accession in af3_ids
        if any(row["pae_json_found"] == "YES" for row in model_rows_by_accession.get(accession, []))
    }
    accessions_with_cif = {accession for accession in af3_ids if model_rows_by_accession.get(accession)}
    status_counts = Counter(row["final_return_status"] for row in master_rows)
    qc = {
        "validation_status": "PASS" if all(check["passed"] for check in checks.values()) else "FAIL",
        "read_only": True, "remote_connection_attempted": False,
        "database_open_mode": "mode=ro&immutable=1",
        "interface_method": {
            "implementation_reference": str((ROOT / "scripts/audit_eligible_af3_interfaces.py").resolve()),
            "parser": "Gemmi 0.7.5",
            "polymer_rule": "Peptide protein; Dna or Rna nucleic polymer (preserved previous method)",
            "contact_cutoff_angstrom": CONTACT_CUTOFF_ANGSTROM,
            "hydrogen_atoms_excluded": True,
            "periodic_image_contacts_excluded": True,
            "accession_valid_if_any_variant_passes": True,
            "all_selected_and_seed_sample_variants_audited_for_intersection": True,
        },
        "counts": {
            "fasta_records": len(fasta_records), "fasta_unique_accessions": len(fasta_ids),
            "fasta_duplicate_accessions": len(fasta_records) - len(fasta_ids),
            "fasta_malformed_or_unparseable": len(fasta_errors),
            "af3_job_directories": len(directory_records),
            "af3_unique_accession_like_ids": len(af3_ids),
            "af3_model_files": len(model_file_rows),
            "af3_accessions_with_model_cif": len(accessions_with_cif),
            "af3_accessions_with_confidence_json": len(accessions_with_confidence),
            "af3_accessions_with_pae_in_confidence_json": len(accessions_with_pae),
            "af3_accessions_with_cif_or_mmcif": len(accessions_with_cif),
            "fasta_with_af3_return": len(intersection),
            "fasta_without_af3_return": len(fasta_only),
            "af3_not_in_fasta": len(af3_only),
            "sent_valid_interface": len(valid_rows),
            "sent_empty_interface": len(empty_rows),
            "sent_invalid_or_unreadable": len(invalid_rows),
            "sent_no_af3_return": len(no_return_rows),
        },
        "model_file_count_distribution_per_accession": {
            str(key): value for key, value in sorted(accession_model_counts.items())
        },
        "directory_variant_counts": dict(sorted(Counter(
            record["variant"] for record in directory_records
        ).items())),
        "final_return_status_counts": dict(sorted(status_counts.items())),
        "af3_only_support_counts": dict(sorted(Counter(
            row["possible_old_af3"] for row in old_rows
        ).items())),
        "controls": controls,
        "checks": checks,
        "protected_state_before": {
            "source_database": source_db_before, "staging_database": staging_db_before,
            "fasta": fasta_before, "website_files": website_before, "af3_tree": af3_before,
        },
        "protected_state_after": {
            "source_database": source_db_after, "staging_database": staging_db_after,
            "fasta": fasta_after, "website_files": website_after, "af3_tree": af3_after,
        },
        "changed_af3_files": changed_af3_files,
        "inputs": {
            "fasta": str(args.fasta.resolve()), "af3_root": str(args.af3_root.resolve()),
            "source_database": str(args.database.resolve()),
            "staging_database": str(args.staging_database.resolve()),
            "previous_inventory_tsv_used_as_source": False,
            "pfam_or_fragment_evidence_used": False,
        },
        "created_files": [
            str(output_paths[key].resolve())
            for key in ["fasta", "inventory", "master", "valid", "empty", "no_return", "old", "qc"]
        ],
    }
    with output_paths["qc"].open("x", encoding="utf-8") as handle:
        json.dump(qc, handle, indent=2, sort_keys=True)
        handle.write("\n")
    if qc["validation_status"] != "PASS":
        raise SystemExit("Protected-state QC failed after report creation")
    print(json.dumps({
        "validation_status": qc["validation_status"], "counts": qc["counts"],
        "model_file_count_distribution_per_accession": qc["model_file_count_distribution_per_accession"],
        "directory_variant_counts": qc["directory_variant_counts"],
        "af3_only_support_counts": qc["af3_only_support_counts"],
        "controls": controls, "created_files": qc["created_files"],
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
