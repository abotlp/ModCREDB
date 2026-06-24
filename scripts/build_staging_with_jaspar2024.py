#!/usr/bin/env python3
"""Build a fresh ModCREDB staging database with HOCOMOCO and JASPAR 2024 links.

The workflow is deliberately non-destructive. It creates all intermediate
SQLite databases in a temporary directory and only writes the requested
output database after every import step succeeds.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_JASPAR_METADATA = "https://jaspar2024.elixir.no/download/database/JASPAR2024.sql.gz"


def run(command: list[str]) -> None:
    print("+", " ".join(command), flush=True)
    subprocess.run(command, check=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--hocomoco-meme", type=Path, required=True)
    parser.add_argument("--hocomoco-annotation", type=Path, required=True)
    parser.add_argument("--hierarchical-tsv", type=Path, required=True)
    parser.add_argument("--output-db", type=Path, required=True)
    parser.add_argument("--audit-dir", type=Path, required=True)
    parser.add_argument("--jaspar-metadata", default=DEFAULT_JASPAR_METADATA)
    parser.add_argument("--uniprot-chunk-size", type=int, default=50)
    parser.add_argument("--uniprot-timeout", type=int, default=30)
    parser.add_argument("--uniprot-sleep", type=float, default=0.2)
    parser.add_argument(
        "--skip-uniprot",
        action="store_true",
        help="Do not fetch UniProt annotations. Intended only for offline parser checks.",
    )
    return parser.parse_args()


def require_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise SystemExit(f"{label} not found: {path}")


def main() -> None:
    args = parse_args()
    for path, label in (
        (args.data_dir / "TF_PWM_chart_final.tsv", "Base TF chart"),
        (args.data_dir / "jaspar.tar.gz", "JASPAR archive"),
        (args.data_dir / "cisbp.tar.gz", "CIS-BP archive"),
        (args.data_dir / "pwms.tar.gz", "Predicted-PWM archive"),
        (args.data_dir / "models.tar.gz", "Model archive"),
        (args.hocomoco_meme, "HOCOMOCO MEME file"),
        (args.hocomoco_annotation, "HOCOMOCO annotation file"),
        (args.hierarchical_tsv, "Integrated hierarchical chart"),
    ):
        require_file(path, label)

    output_db = args.output_db.resolve()
    if output_db.exists():
        raise SystemExit(f"Refusing to overwrite existing output database: {output_db}")
    args.audit_dir.mkdir(parents=True, exist_ok=True)
    output_db.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="modcredb-rebuild-") as temp_dir_text:
        temp_dir = Path(temp_dir_text)
        base_db = temp_dir / "base.sqlite"
        hocomoco_db = temp_dir / "hocomoco.sqlite"
        final_db = temp_dir / "final.sqlite"
        candidate_tsv = temp_dir / "jaspar_metadata_link_candidates.tsv"
        apply_report = temp_dir / "jaspar_metadata_link_apply.md"

        run(
            [
                sys.executable,
                str(PROJECT_ROOT / "import_db.py"),
                "--data-dir",
                str(args.data_dir.resolve()),
                "--db",
                str(base_db),
            ]
        )
        run(
            [
                sys.executable,
                str(PROJECT_ROOT / "import_hocomoco.py"),
                "--input-db",
                str(base_db),
                "--output-db",
                str(hocomoco_db),
                "--hocomoco-meme",
                str(args.hocomoco_meme.resolve()),
                "--hocomoco-annotation",
                str(args.hocomoco_annotation.resolve()),
                "--hierarchical-tsv",
                str(args.hierarchical_tsv.resolve()),
            ]
        )
        if not args.skip_uniprot:
            run(
                [
                    sys.executable,
                    str(PROJECT_ROOT / "enrich_uniprot.py"),
                    "--db",
                    str(hocomoco_db),
                    "--chunk-size",
                    str(args.uniprot_chunk_size),
                    "--timeout",
                    str(args.uniprot_timeout),
                    "--sleep",
                    str(args.uniprot_sleep),
                ]
            )
        run(
            [
                sys.executable,
                str(PROJECT_ROOT / "scripts" / "import_jaspar_metadata_links.py"),
                "--db",
                str(hocomoco_db),
                "--metadata",
                args.jaspar_metadata,
                "--apply",
                "--output-db",
                str(final_db),
                "--out",
                str(candidate_tsv),
                "--report",
                str(apply_report),
            ]
        )

        shutil.copy2(final_db, output_db)
        shutil.copy2(candidate_tsv, args.audit_dir / candidate_tsv.name)
        shutil.copy2(apply_report, args.audit_dir / apply_report.name)

    print(f"Built fresh staging database: {output_db}")
    print(f"JASPAR audit files: {args.audit_dir}")


if __name__ == "__main__":
    main()
