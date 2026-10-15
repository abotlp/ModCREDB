#!/usr/bin/env python3
"""Audit protein-DNA contacts for the two Baldo ModCRE coordinate batches.

This is the reproducible producer for the per-PDB evidence consumed by
``build_complete_modcre_model_manifest.py``.  It does not modify coordinates or
SQLite.  A model is interface-positive when at least one non-hydrogen protein
atom and one non-hydrogen DNA atom are within the requested cutoff.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import multiprocessing as mp
import re
from collections import Counter
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ORIGINAL = Path("/data/sbi/interchange/boliva/patricia/models")
DEFAULT_CORRECTED = Path("/data/sbi/interchange/boliva/patricia/models_failed_in_DNA")
DEFAULT_OUTPUT = REPO_ROOT / "outputs/all_baldo_model_interface_audit.tsv"
MODEL_RE = re.compile(r"^(TFS|DIMER)_([^:]+):(\d+):(\d+)_(.+)\.pdb$")
FIELDS = (
    "batch",
    "accession",
    "model_id",
    "model_path",
    "model_type",
    "residue_start",
    "residue_end",
    "dna_interaction",
    "dna_contact_count",
    "dna_qc_status",
)


def load_dna_qc_module() -> Any:
    source = REPO_ROOT / "scripts/annotate_structure_dna_interaction.py"
    spec = importlib.util.spec_from_file_location("modcredb_dna_qc", source)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load DNA-QC implementation: {source}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


DNA_QC = load_dna_qc_module()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Audit protein-DNA contacts in original and corrected ModCRE PDB batches."
    )
    parser.add_argument("--original-dir", type=Path, default=DEFAULT_ORIGINAL)
    parser.add_argument("--corrected-dir", type=Path, default=DEFAULT_CORRECTED)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--cutoff", type=float, default=4.5)
    parser.add_argument("--workers", type=int, default=min(8, mp.cpu_count()))
    return parser.parse_args()


def audit(item: tuple[str, str, float]) -> tuple[object, ...]:
    batch, raw_path, cutoff = item
    path = Path(raw_path)
    match = MODEL_RE.fullmatch(path.name)
    if match is None:
        return (batch, "", path.stem, str(path), "", "", "", "", "", "UNPARSED_NAME")
    kind, accession, start, end, _ = match.groups()
    try:
        protein, dna = DNA_QC.parse_pdb_path(path)
        if not protein:
            return (batch, accession, path.stem, str(path), kind, start, end, 0, 0, "NO_PROTEIN_ATOMS")
        if not dna:
            return (batch, accession, path.stem, str(path), kind, start, end, 0, 0, "NO_DNA_ATOMS")
        contacts = DNA_QC.contact_count(protein, dna, cutoff)
        return (
            batch,
            accession,
            path.stem,
            str(path),
            kind,
            start,
            end,
            1 if contacts else 0,
            contacts,
            "PASS_NONEMPTY_PROTEIN_DNA_INTERFACE" if contacts else "EMPTY_PROTEIN_DNA_INTERFACE",
        )
    except Exception as exc:  # Preserve the failing path and exact parser error.
        return (
            batch,
            accession,
            path.stem,
            str(path),
            kind,
            start,
            end,
            "",
            "",
            f"ERROR:{type(exc).__name__}:{exc}",
        )


def main() -> int:
    args = parse_args()
    if args.cutoff <= 0:
        raise SystemExit("--cutoff must be positive")
    if args.workers <= 0:
        raise SystemExit("--workers must be positive")
    if args.output.exists():
        raise SystemExit(f"output exists; refusing to overwrite: {args.output}")

    items: list[tuple[str, str, float]] = []
    for batch, root in (("original", args.original_dir), ("corrected_rerun", args.corrected_dir)):
        if not root.is_dir():
            raise SystemExit(f"coordinate directory not found: {root}")
        paths = sorted(root.glob("*.pdb"))
        print(f"{batch}: {len(paths)} PDBs", flush=True)
        items.extend((batch, str(path.resolve()), args.cutoff) for path in paths)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    counts: Counter[tuple[str, str]] = Counter()
    with args.output.open("x", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(FIELDS)
        with mp.Pool(processes=min(args.workers, mp.cpu_count())) as pool:
            for index, row in enumerate(pool.imap_unordered(audit, items, chunksize=20), start=1):
                writer.writerow(row)
                counts[(str(row[0]), str(row[-1]).split(":", 1)[0])] += 1
                if index % 2000 == 0:
                    print(f"processed {index}/{len(items)}", flush=True)

    print("SUMMARY", flush=True)
    for (batch, status), count in sorted(counts.items()):
        print(f"{batch}\t{status}\t{count}", flush=True)
    print(f"output\t{args.output.resolve()}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
