#!/usr/bin/env python3
"""Build the complete, de-duplicated ModCRE coordinate manifest.

The three input collections are treated as immutable evidence:

* the legacy ``models.tar.gz`` archive;
* Baldo's original 46,892-PDB directory;
* Baldo's corrected 15,957-PDB DNA rerun directory.

Files are de-duplicated by biological model identity (accession, TFS/DIMER,
residue interval, template, template chain and rank).  A corrected-rerun file
is selected for an exact identity when available, otherwise the original file
is selected, followed by the legacy archive.  This avoids showing a broken and
corrected coordinate version as two independent biological models while still
retaining every distinct model identity.

The script is read-only with respect to all source files and SQLite.  It writes
only the requested TSV/JSON outputs and never creates the serving archive.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import sqlite3
import tarfile
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATABASE = Path("data/tf_webdb_local_staging.sqlite")
DEFAULT_LEGACY_ARCHIVE = Path("/home/patricia/TF_database_Baldo_data/models.tar.gz")
DEFAULT_ORIGINAL_DIR = Path("/data/sbi/interchange/boliva/patricia/models")
DEFAULT_RERUN_DIR = Path("/data/sbi/interchange/boliva/patricia/models_failed_in_DNA")
DEFAULT_INTERFACE_AUDIT = ROOT / "outputs/all_baldo_model_interface_audit.tsv"
DEFAULT_OUTPUT = Path("outputs/modcre_complete_model_manifest.tsv")
DEFAULT_REJECTED = Path("outputs/modcre_complete_model_manifest_rejected.tsv")
DEFAULT_QC = Path("outputs/modcre_complete_model_manifest_qc.json")

EXPECTED_LEGACY_MODELS = 7_776
EXPECTED_ORIGINAL_MODELS = 46_892
EXPECTED_RERUN_MODELS = 15_957
EXPECTED_ORIGINAL_RERUN_OVERLAP = 7_851
EXPECTED_LEGACY_RERUN_OVERLAP = 3_791
EXPECTED_DISTINCT_MODELS = 58_983
EXPECTED_MAPPED_MODELS = 58_976
EXPECTED_REJECTED_MODELS = 7
EXPECTED_MAPPED_NO_PROTEIN_MODELS = 4
CONTACT_CUTOFF = 4.5

OLD_MODEL_RE = re.compile(
    r"^(?P<kind>TFS|DIMER)_(?:tr|sp)_(?P<accession>[^_]+)_[^:]+:"
    r"(?P<start>\d+):(?P<end>\d+)_(?P<tail>.+)$"
)
NEW_MODEL_RE = re.compile(
    r"^(?P<kind>TFS|DIMER)_(?P<accession>[^:]+):"
    r"(?P<start>\d+):(?P<end>\d+)_(?P<tail>.+)$"
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()

PROTEIN_RESIDUES = frozenset(
    {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
        "ASX", "GLX", "HYP", "MSE", "PYL", "SEC", "XLE",
    }
)
DNA_RESIDUES = frozenset(
    {
        "DA", "DC", "DG", "DT", "DI",
        "A", "C", "G", "T", "U", "I",
        "ADE", "CYT", "GUA", "THY", "URI",
    }
)

MANIFEST_COLUMNS = (
    "tf_id",
    "model_id",
    "model_kind",
    "residue_start",
    "residue_end",
    "template_pdb",
    "template_chain",
    "model_rank",
    "coordinate_source",
    "selection_reason",
    "source_path",
    "source_archive_path",
    "source_member_path",
    "serving_member_path",
    "file_size",
    "sha256",
    "protein_chain_count",
    "protein_chains",
    "dna_chain_count",
    "dna_chains",
    "dna_interaction",
    "dna_contact_count",
    "dna_qc_cutoff",
    "dna_qc_status",
    "display_priority",
)

REJECTED_COLUMNS = MANIFEST_COLUMNS + ("rejection_reason",)


class ValidationError(RuntimeError):
    pass


@dataclass(frozen=True)
class Candidate:
    model_id: str
    model_kind: str
    tf_id: str
    residue_start: int
    residue_end: int
    template_pdb: str
    template_chain: str
    model_rank: str
    coordinate_source: str
    source_path: str
    source_archive_path: str
    source_member_path: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database", type=Path, default=DEFAULT_DATABASE)
    parser.add_argument("--legacy-archive", type=Path, default=DEFAULT_LEGACY_ARCHIVE)
    parser.add_argument("--original-dir", type=Path, default=DEFAULT_ORIGINAL_DIR)
    parser.add_argument("--rerun-dir", type=Path, default=DEFAULT_RERUN_DIR)
    parser.add_argument(
        "--interface-audit",
        type=Path,
        default=DEFAULT_INTERFACE_AUDIT,
        help="Existing complete per-PDB original/rerun interface audit",
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--rejected-output", type=Path, default=DEFAULT_REJECTED)
    parser.add_argument("--qc-output", type=Path, default=DEFAULT_QC)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--force", action="store_true", help="Permit replacing output files")
    return parser.parse_args()


def ensure_output_available(path: Path, force: bool) -> None:
    if path.exists() and not force:
        raise ValidationError(f"output already exists: {path}; use --force to replace it")
    path.parent.mkdir(parents=True, exist_ok=True)


def normalize_model_name(name: str) -> tuple[str, re.Match[str]] | None:
    stem = Path(name).stem
    match = OLD_MODEL_RE.fullmatch(stem) or NEW_MODEL_RE.fullmatch(stem)
    if match is None:
        return None
    start = int(match.group("start"))
    end = int(match.group("end"))
    if start < 1 or end < start:
        raise ValidationError(f"invalid residue interval in model name: {name}")
    model_id = (
        f"{match.group('kind')}_{match.group('accession')}:"
        f"{start}:{end}_{match.group('tail')}"
    )
    return model_id, match


def candidate_from_name(
    name: str,
    coordinate_source: str,
    *,
    source_path: str = "",
    source_archive_path: str = "",
    source_member_path: str = "",
) -> Candidate | None:
    parsed = normalize_model_name(name)
    if parsed is None:
        return None
    model_id, match = parsed
    tail_parts = match.group("tail").rsplit("_", 2)
    if len(tail_parts) != 3:
        raise ValidationError(f"cannot parse template/chain/rank from {name}")
    template_pdb, template_chain, model_rank = tail_parts
    return Candidate(
        model_id=model_id,
        model_kind=match.group("kind"),
        tf_id=match.group("accession"),
        residue_start=int(match.group("start")),
        residue_end=int(match.group("end")),
        template_pdb=template_pdb,
        template_chain=template_chain,
        model_rank=model_rank,
        coordinate_source=coordinate_source,
        source_path=source_path,
        source_archive_path=source_archive_path,
        source_member_path=source_member_path,
    )


def collect_directory(path: Path, label: str) -> dict[str, Candidate]:
    if not path.is_dir():
        raise ValidationError(f"missing source directory: {path}")
    result: dict[str, Candidate] = {}
    unparsed: list[str] = []
    for pdb_path in sorted(path.glob("*.pdb")):
        candidate = candidate_from_name(
            pdb_path.name,
            label,
            source_path=str(pdb_path.resolve()),
        )
        if candidate is None:
            unparsed.append(pdb_path.name)
            continue
        if candidate.model_id in result:
            raise ValidationError(f"duplicate normalized model in {path}: {candidate.model_id}")
        result[candidate.model_id] = candidate
    if unparsed:
        raise ValidationError(f"unparsed PDB names in {path}: {', '.join(unparsed[:10])}")
    return result


def collect_legacy(path: Path) -> dict[str, Candidate]:
    if not path.is_file():
        raise ValidationError(f"missing legacy archive: {path}")
    result: dict[str, Candidate] = {}
    with tarfile.open(path, "r:gz") as archive:
        for member in archive.getmembers():
            # The legacy bundle deliberately contains a separate
            # ``models_failed/`` copy for records rejected by the historical
            # import.  Those files often share a biological model identity
            # with the public ``models/`` member but are not an additional
            # model.  Only the serving collection is part of this inventory.
            if (
                not member.isfile()
                or not member.name.startswith("models/")
                or not member.name.casefold().endswith(".pdb")
            ):
                continue
            candidate = candidate_from_name(
                member.name,
                "legacy_archive",
                source_archive_path=str(path.resolve()),
                source_member_path=member.name,
            )
            if candidate is None:
                continue
            if candidate.model_id in result:
                raise ValidationError(
                    f"duplicate normalized model in legacy archive: {candidate.model_id}"
                )
            result[candidate.model_id] = candidate
    return result


def load_interface_audit(path: Path) -> dict[str, dict[str, str]]:
    if not path.is_file():
        raise ValidationError(f"missing complete interface audit: {path}")
    result: dict[str, dict[str, str]] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "batch", "accession", "model_id", "model_path", "model_type",
            "residue_start", "residue_end", "dna_interaction",
            "dna_contact_count", "dna_qc_status",
        }
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValidationError(f"interface audit lacks required columns: {path}")
        for row in reader:
            resolved = str(Path(row["model_path"]).resolve())
            if resolved in result:
                raise ValidationError(f"duplicate model_path in interface audit: {resolved}")
            result[resolved] = row
    return result


def is_hydrogen(line: str, atom_name: str) -> bool:
    element = line[76:78].strip().upper() if len(line) >= 78 else ""
    if element in {"H", "D"}:
        return True
    if not element:
        return atom_name.lstrip("0123456789").upper().startswith(("H", "D"))
    return False


def analyze_pdb_bytes(content: bytes, calculate_contacts: bool) -> dict[str, object]:
    protein_chains: set[str] = set()
    dna_chains: set[str] = set()
    protein_atoms: list[tuple[float, float, float]] = []
    dna_atoms: list[tuple[float, float, float]] = []
    for raw_line in content.splitlines():
        line = raw_line.decode("ascii", errors="replace")
        if not (line.startswith("ATOM  ") or line.startswith("HETATM")) or len(line) < 54:
            continue
        residue = line[17:20].strip().upper()
        if residue not in PROTEIN_RESIDUES and residue not in DNA_RESIDUES:
            continue
        chain = line[21:22].strip() or "_"
        atom_name = line[12:16].strip().upper()
        if is_hydrogen(line, atom_name):
            continue
        try:
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError:
            continue
        if residue in PROTEIN_RESIDUES:
            protein_chains.add(chain)
            if calculate_contacts:
                protein_atoms.append(xyz)
        else:
            dna_chains.add(chain)
            if calculate_contacts:
                dna_atoms.append(xyz)

    contacts: int | None = None
    status = "ok"
    if not protein_chains:
        status = "no_protein_atoms"
    elif not dna_chains:
        status = "no_nucleic_atoms"
    elif calculate_contacts:
        cutoff2 = CONTACT_CUTOFF * CONTACT_CUTOFF
        grid: dict[tuple[int, int, int], list[tuple[float, float, float]]] = defaultdict(list)
        for x, y, z in dna_atoms:
            key = (
                math.floor(x / CONTACT_CUTOFF),
                math.floor(y / CONTACT_CUTOFF),
                math.floor(z / CONTACT_CUTOFF),
            )
            grid[key].append((x, y, z))
        contacts = 0
        for x, y, z in protein_atoms:
            gx, gy, gz = (
                math.floor(x / CONTACT_CUTOFF),
                math.floor(y / CONTACT_CUTOFF),
                math.floor(z / CONTACT_CUTOFF),
            )
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    for dz in (-1, 0, 1):
                        for nx, ny, nz in grid.get((gx + dx, gy + dy, gz + dz), ()):
                            if (x - nx) ** 2 + (y - ny) ** 2 + (z - nz) ** 2 <= cutoff2:
                                contacts += 1
    return {
        "file_size": len(content),
        "sha256": hashlib.sha256(content).hexdigest(),
        "protein_chain_count": len(protein_chains),
        "protein_chains": ",".join(sorted(protein_chains)),
        "dna_chain_count": len(dna_chains),
        "dna_chains": ",".join(sorted(dna_chains)),
        "dna_contact_count": contacts,
        "dna_qc_status": status,
    }


def analyze_path(item: tuple[str, str]) -> tuple[str, dict[str, object]]:
    model_id, raw_path = item
    path = Path(raw_path)
    return model_id, analyze_pdb_bytes(path.read_bytes(), calculate_contacts=False)


def display_priority(row: dict[str, object]) -> int:
    interaction = int(row["dna_interaction"])
    complete_dna = int(row["dna_chain_count"]) >= 2
    dimer = row["model_kind"] == "DIMER"
    if interaction and complete_dna:
        return 0 if dimer else 1
    if interaction:
        return 2 if dimer else 3
    return 4 if dimer else 5


def write_tsv(path: Path, columns: Iterable[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(columns), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def database_tf_ids(path: Path) -> set[str]:
    if not path.is_file():
        raise ValidationError(f"missing staging database: {path}")
    connection = sqlite3.connect(f"file:{path}?mode=ro&immutable=1", uri=True)
    try:
        return {str(row[0]) for row in connection.execute("SELECT tf_id FROM tf")}
    finally:
        connection.close()


def main() -> int:
    args = parse_args()
    if args.workers < 1:
        raise ValidationError("--workers must be at least 1")
    for output in (args.output, args.rejected_output, args.qc_output):
        ensure_output_available(output, args.force)

    legacy = collect_legacy(args.legacy_archive)
    original = collect_directory(args.original_dir, "original")
    rerun = collect_directory(args.rerun_dir, "corrected_rerun")
    if len(legacy) != EXPECTED_LEGACY_MODELS:
        raise ValidationError(f"expected {EXPECTED_LEGACY_MODELS} legacy models, found {len(legacy)}")
    if len(original) != EXPECTED_ORIGINAL_MODELS:
        raise ValidationError(f"expected {EXPECTED_ORIGINAL_MODELS} original models, found {len(original)}")
    if len(rerun) != EXPECTED_RERUN_MODELS:
        raise ValidationError(f"expected {EXPECTED_RERUN_MODELS} rerun models, found {len(rerun)}")
    if len(set(original) & set(rerun)) != EXPECTED_ORIGINAL_RERUN_OVERLAP:
        raise ValidationError("unexpected original/rerun identity overlap")
    if len(set(legacy) & set(rerun)) != EXPECTED_LEGACY_RERUN_OVERLAP:
        raise ValidationError("unexpected legacy/rerun identity overlap")

    selected: dict[str, Candidate] = dict(legacy)
    selected.update(original)
    selected.update(rerun)
    if len(selected) != EXPECTED_DISTINCT_MODELS:
        raise ValidationError(
            f"expected {EXPECTED_DISTINCT_MODELS} distinct models, found {len(selected)}"
        )

    interface_audit = load_interface_audit(args.interface_audit)
    raw_selected = {
        model_id: candidate
        for model_id, candidate in selected.items()
        if candidate.coordinate_source in {"original", "corrected_rerun"}
    }
    raw_items = [(model_id, candidate.source_path) for model_id, candidate in raw_selected.items()]
    raw_analysis: dict[str, dict[str, object]] = {}
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        for index, (model_id, analysis) in enumerate(executor.map(analyze_path, raw_items, chunksize=40), 1):
            raw_analysis[model_id] = analysis
            if index % 5_000 == 0:
                print(f"Analyzed raw PDBs: {index}/{len(raw_items)}", flush=True)

    legacy_selected = {
        model_id: candidate
        for model_id, candidate in selected.items()
        if candidate.coordinate_source == "legacy_archive"
    }
    legacy_analysis: dict[str, dict[str, object]] = {}
    with tarfile.open(args.legacy_archive, "r:gz") as archive:
        wanted = {candidate.source_member_path: model_id for model_id, candidate in legacy_selected.items()}
        for member in archive:
            model_id = wanted.get(member.name)
            if model_id is None:
                continue
            extracted = archive.extractfile(member)
            if extracted is None:
                raise ValidationError(f"cannot read legacy archive member: {member.name}")
            with extracted:
                legacy_analysis[model_id] = analyze_pdb_bytes(
                    extracted.read(), calculate_contacts=True
                )
    if set(legacy_analysis) != set(legacy_selected):
        missing = sorted(set(legacy_selected) - set(legacy_analysis))
        raise ValidationError(f"selected legacy members not read: {', '.join(missing[:10])}")

    tf_ids = database_tf_ids(args.database)
    accepted_rows: list[dict[str, object]] = []
    rejected_rows: list[dict[str, object]] = []
    source_counts = Counter()
    for model_id in sorted(selected):
        candidate = selected[model_id]
        analysis = (
            legacy_analysis[model_id]
            if candidate.coordinate_source == "legacy_archive"
            else raw_analysis[model_id]
        )
        if candidate.coordinate_source == "legacy_archive":
            contact_count = analysis["dna_contact_count"]
            dna_status = str(analysis["dna_qc_status"])
            interaction = 1 if dna_status == "ok" and int(contact_count or 0) > 0 else 0
        else:
            audit = interface_audit.get(str(Path(candidate.source_path).resolve()))
            if audit is None:
                raise ValidationError(f"selected raw model absent from interface audit: {candidate.source_path}")
            if audit["accession"] != candidate.tf_id or audit["model_id"] != candidate.model_id:
                raise ValidationError(f"interface-audit identity mismatch: {candidate.source_path}")
            contact_count = int(audit["dna_contact_count"] or 0)
            interaction = int(audit["dna_interaction"] or 0)
            audit_status = audit["dna_qc_status"]
            dna_status = (
                "ok"
                if audit_status in {
                    "PASS_NONEMPTY_PROTEIN_DNA_INTERFACE",
                    "EMPTY_PROTEIN_DNA_INTERFACE",
                }
                else audit_status.casefold()
            )
        selection_reason = {
            "corrected_rerun": "CORRECTED_EXACT_VERSION_SELECTED",
            "original": "ORIGINAL_MODEL_NO_CORRECTED_EXACT_VERSION",
            "legacy_archive": "LEGACY_MODEL_NO_CORRECTED_EXACT_VERSION",
        }[candidate.coordinate_source]
        row: dict[str, object] = {
            "tf_id": candidate.tf_id,
            "model_id": candidate.model_id,
            "model_kind": candidate.model_kind,
            "residue_start": candidate.residue_start,
            "residue_end": candidate.residue_end,
            "template_pdb": candidate.template_pdb,
            "template_chain": candidate.template_chain,
            "model_rank": candidate.model_rank,
            "coordinate_source": candidate.coordinate_source,
            "selection_reason": selection_reason,
            "source_path": candidate.source_path,
            "source_archive_path": candidate.source_archive_path,
            "source_member_path": candidate.source_member_path,
            "serving_member_path": f"models/{candidate.model_id}.pdb",
            "file_size": analysis["file_size"],
            "sha256": analysis["sha256"],
            "protein_chain_count": analysis["protein_chain_count"],
            "protein_chains": analysis["protein_chains"],
            "dna_chain_count": analysis["dna_chain_count"],
            "dna_chains": analysis["dna_chains"],
            "dna_interaction": interaction,
            "dna_contact_count": contact_count,
            "dna_qc_cutoff": CONTACT_CUTOFF,
            "dna_qc_status": dna_status,
        }
        row["display_priority"] = display_priority(row)
        source_counts[candidate.coordinate_source] += 1
        if candidate.tf_id in tf_ids:
            accepted_rows.append(row)
        else:
            rejected_rows.append({**row, "rejection_reason": "ACCESSION_NOT_IN_TF_TABLE"})

    if len(accepted_rows) != EXPECTED_MAPPED_MODELS:
        raise ValidationError(
            f"expected {EXPECTED_MAPPED_MODELS} mapped models, found {len(accepted_rows)}"
        )
    if len(rejected_rows) != EXPECTED_REJECTED_MODELS:
        raise ValidationError(
            f"expected {EXPECTED_REJECTED_MODELS} rejected models, found {len(rejected_rows)}"
        )
    if len({row["model_id"] for row in accepted_rows}) != len(accepted_rows):
        raise ValidationError("accepted manifest contains duplicate model_id values")
    if len({row["serving_member_path"] for row in accepted_rows}) != len(accepted_rows):
        raise ValidationError("accepted manifest contains duplicate serving_member_path values")
    parse_state_counts = Counter(str(row["dna_qc_status"]) for row in accepted_rows)
    unexpected_states = set(parse_state_counts) - {"ok", "no_protein_atoms"}
    if unexpected_states:
        raise ValidationError(f"unexpected selected PDB parse states: {dict(parse_state_counts)}")
    if parse_state_counts["no_protein_atoms"] != EXPECTED_MAPPED_NO_PROTEIN_MODELS:
        raise ValidationError(
            "expected "
            f"{EXPECTED_MAPPED_NO_PROTEIN_MODELS} mapped DNA-only legacy PDBs, found "
            f"{parse_state_counts['no_protein_atoms']}"
        )

    accepted_rows.sort(
        key=lambda row: (
            str(row["tf_id"]),
            int(row["display_priority"]),
            int(row["residue_start"]),
            str(row["model_id"]),
        )
    )
    rejected_rows.sort(key=lambda row: (str(row["tf_id"]), str(row["model_id"])))
    write_tsv(args.output, MANIFEST_COLUMNS, accepted_rows)
    write_tsv(args.rejected_output, REJECTED_COLUMNS, rejected_rows)

    qc = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "database": str(args.database.resolve()),
        "legacy_archive": str(args.legacy_archive.resolve()),
        "original_directory": str(args.original_dir.resolve()),
        "corrected_rerun_directory": str(args.rerun_dir.resolve()),
        "interface_audit": str(args.interface_audit.resolve()),
        "interface_audit_sha256": sha256_file(args.interface_audit),
        "physical_parseable_counts": {
            "legacy": len(legacy),
            "original": len(original),
            "corrected_rerun": len(rerun),
        },
        "identity_overlaps": {
            "original_corrected_rerun": len(set(original) & set(rerun)),
            "legacy_corrected_rerun": len(set(legacy) & set(rerun)),
            "legacy_original": len(set(legacy) & set(original)),
        },
        "distinct_model_identities": len(selected),
        "mapped_models": len(accepted_rows),
        "rejected_models": len(rejected_rows),
        "rejected_accessions": sorted({str(row["tf_id"]) for row in rejected_rows}),
        "selected_coordinate_sources": dict(sorted(source_counts.items())),
        "selected_model_kinds": dict(
            sorted(Counter(str(row["model_kind"]) for row in accepted_rows).items())
        ),
        "dna_chain_counts": dict(
            sorted(Counter(str(row["dna_chain_count"]) for row in accepted_rows).items())
        ),
        "dna_interaction_counts": dict(
            sorted(Counter(str(row["dna_interaction"]) for row in accepted_rows).items())
        ),
        "dna_qc_status_counts": dict(sorted(parse_state_counts.items())),
        "display_priority_counts": dict(
            sorted(Counter(str(row["display_priority"]) for row in accepted_rows).items())
        ),
        "all_invariants_passed": True,
    }
    args.qc_output.write_text(json.dumps(qc, indent=2, sort_keys=True) + "\n")
    print(f"Distinct biological models: {len(selected)}")
    print(f"Mapped public models: {len(accepted_rows)}")
    print(f"Unmapped preserved models: {len(rejected_rows)}")
    print(f"Coordinate sources: {dict(sorted(source_counts.items()))}")
    print(f"Manifest: {args.output}")
    print(f"QC: {args.qc_output}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ValidationError as exc:
        raise SystemExit(f"ERROR: {exc}")
