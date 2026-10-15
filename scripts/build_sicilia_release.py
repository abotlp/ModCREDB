#!/usr/bin/env python3
"""Build a checksummed, non-deploying ModCREDB release for Sicilia.

The release contains the runtime source, a copy of the final SQLite database
with AF3 paths translated to the established production-tree convention, and
only the archive-backed cache members absent from the previous production
database.  AF3 coordinate/confidence files are not duplicated in the release;
an exact SHA256 manifest is emitted so deployment must verify the copies that
already exist below the production tree.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import py_compile
import shutil
import sqlite3
import subprocess
import tarfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATABASE = ROOT / "data/tf_webdb_local_staging.sqlite"
DEFAULT_RESTORATION_PREVIEW = ROOT / "outputs/manuscript_5417_restoration_preview.tsv"
DEFAULT_CACHE = Path("/home/patricia/TF_database_Baldo_data/modcredb_model_cache")
DEFAULT_OUTPUT = ROOT / "deploy/modcredb_sicilia_release_20261015"
DEFAULT_ARCHIVE = ROOT / "deploy/modcredb_sicilia_release_20261015.tar.gz"
DEFAULT_PRODUCTION_ROOT = Path("/home/patricia/modcredb_candidate")

EXPECTED = {
    "tf": 5417,
    "primary": 5417,
    "motif_ref": 117372,
    "active_pdb": 62505,
    "active_af3": 114,
    "final_unresolved": 108,
}


class ReleaseError(RuntimeError):
    pass


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build (but do not deploy) the checksummed Sicilia runtime release."
    )
    parser.add_argument("--database", type=Path, default=DEFAULT_DATABASE)
    parser.add_argument(
        "--restoration-preview", type=Path, default=DEFAULT_RESTORATION_PREVIEW
    )
    parser.add_argument("--model-cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--archive", type=Path, default=DEFAULT_ARCHIVE)
    parser.add_argument("--production-root", type=Path, default=DEFAULT_PRODUCTION_ROOT)
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def connect_ro(path: Path) -> sqlite3.Connection:
    if not path.is_file():
        raise ReleaseError(f"required database is missing: {path}")
    connection = sqlite3.connect(f"file:{path.resolve()}?mode=ro", uri=True)
    connection.row_factory = sqlite3.Row
    return connection


def active_members(path: Path) -> dict[str, sqlite3.Row]:
    with connect_ro(path) as connection:
        rows = connection.execute(
            """
            SELECT id, tf_id, source, model_id, member_path
            FROM structure_file
            WHERE status='active' AND file_type='pdb'
            ORDER BY member_path
            """
        ).fetchall()
    result: dict[str, sqlite3.Row] = {}
    for row in rows:
        member = PurePosixPath(str(row["member_path"]))
        if member.is_absolute() or any(part in {"", ".", ".."} for part in member.parts):
            raise ReleaseError(f"unsafe cache member path: {member}")
        if member.as_posix() in result:
            raise ReleaseError(f"duplicate active cache member path: {member}")
        result[member.as_posix()] = row
    return result


def restored_alphafold_accessions(path: Path) -> dict[str, int]:
    import csv

    if not path.is_file():
        raise ReleaseError(f"restoration preview is missing: {path}")
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    selected = {
        str(row["tf_id"]): int(row["active_models_planned"])
        for row in rows
        if int(row.get("active_models_planned") or 0) > 0
    }
    if len(selected) != 10 or sum(selected.values()) != 55:
        raise ReleaseError(
            "restoration preview does not identify the frozen 10-accession/55-model "
            f"AlphaFold restoration: {selected}"
        )
    return selected


def validate_database(connection: sqlite3.Connection) -> dict[str, Any]:
    if connection.execute("PRAGMA quick_check").fetchone()[0] != "ok":
        raise ReleaseError("SQLite quick_check failed")
    foreign_keys = connection.execute("PRAGMA foreign_key_check").fetchall()
    if foreign_keys:
        raise ReleaseError(f"SQLite foreign-key violations: {len(foreign_keys)}")
    observed = {
        "tf": connection.execute("SELECT COUNT(*) FROM tf").fetchone()[0],
        "primary": connection.execute(
            "SELECT COUNT(*) FROM tf_primary_annotation"
        ).fetchone()[0],
        "motif_ref": connection.execute("SELECT COUNT(*) FROM motif_ref").fetchone()[0],
        "active_pdb": connection.execute(
            "SELECT COUNT(*) FROM structure_file WHERE status='active' AND file_type='pdb'"
        ).fetchone()[0],
        "active_af3": connection.execute(
            """SELECT COUNT(*) FROM structure_model_assignment
               WHERE is_active=1 AND UPPER(source)='AF3'"""
        ).fetchone()[0],
        "final_unresolved": connection.execute(
            """SELECT COUNT(*) FROM tf_structure_status
               WHERE primary_structural_status=
                 'PWM_PRESENT_NO_VALIDATED_DNA_BOUND_STRUCTURE'"""
        ).fetchone()[0],
    }
    if observed != EXPECTED:
        raise ReleaseError(f"release database counts differ: {observed} != {EXPECTED}")
    return observed


def target_af3_dir(production_root: Path, accession: str, assignment_id: int) -> Path:
    return (
        production_root
        / "external/baldo_model_inventory/fragments_cif_audit/AF3/production_active"
        / f"{accession}__assignment_{assignment_id}"
    )


def rewrite_af3_paths(database: Path, production_root: Path, manifest: Path) -> int:
    connection = sqlite3.connect(database)
    connection.row_factory = sqlite3.Row
    connection.execute("PRAGMA foreign_keys=ON")
    manifest_rows: list[dict[str, Any]] = []
    source_hashes: dict[str, str] = {}
    try:
        connection.execute("BEGIN IMMEDIATE")
        assignments = connection.execute(
            """SELECT * FROM structure_model_assignment
               WHERE is_active=1 AND UPPER(source)='AF3' ORDER BY id"""
        ).fetchall()
        if len(assignments) != EXPECTED["active_af3"]:
            raise ReleaseError(f"expected 114 AF3 assignments, found {len(assignments)}")
        for assignment in assignments:
            assignment_id = int(assignment["id"])
            accession = str(assignment["display_accession"])
            old_model = Path(str(assignment["model_path"]))
            if not old_model.is_file():
                raise ReleaseError(f"AF3 model source is missing: {old_model}")
            destination_dir = target_af3_dir(production_root, accession, assignment_id)
            new_model = destination_dir / old_model.name
            artifacts = connection.execute(
                """SELECT * FROM structure_confidence_artifact
                   WHERE model_path=? ORDER BY id""",
                (str(old_model),),
            ).fetchall()
            if Counter(str(row["artifact_type"]) for row in artifacts) != Counter(
                {"MMCIF": 1, "PLDDT_JSON": 1, "PAE_JSON": 1}
            ):
                raise ReleaseError(f"AF3 artifact set is incomplete for assignment {assignment_id}")
            connection.execute(
                "UPDATE structure_model_assignment SET model_path=? WHERE id=?",
                (str(new_model), assignment_id),
            )
            for artifact in artifacts:
                old_artifact = Path(str(artifact["artifact_path"]))
                if not old_artifact.is_file():
                    raise ReleaseError(f"AF3 artifact source is missing: {old_artifact}")
                new_artifact = destination_dir / old_artifact.name
                observed_hash = source_hashes.setdefault(
                    str(old_artifact), sha256_file(old_artifact)
                )
                expected_hash = str(artifact["sha256"] or "").lower()
                if observed_hash != expected_hash:
                    raise ReleaseError(f"AF3 artifact SHA mismatch: {old_artifact}")
                connection.execute(
                    """UPDATE structure_confidence_artifact
                       SET model_path=?, artifact_path=? WHERE id=?""",
                    (str(new_model), str(new_artifact), int(artifact["id"])),
                )
                manifest_rows.append(
                    {
                        "assignment_id": assignment_id,
                        "accession": accession,
                        "artifact_type": str(artifact["artifact_type"]),
                        "target_path": str(new_artifact),
                        "sha256": expected_hash,
                        "size_bytes": old_artifact.stat().st_size,
                    }
                )
        validate_database(connection)
        connection.commit()
    except BaseException:
        connection.rollback()
        raise
    finally:
        connection.close()

    import csv

    with manifest.open("x", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "assignment_id",
                "accession",
                "artifact_type",
                "target_path",
                "sha256",
                "size_bytes",
            ),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(manifest_rows)
    return len(manifest_rows)


def git_commit() -> str:
    completed = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=True,
    )
    return completed.stdout.strip()


def write_checksums(output_dir: Path) -> None:
    lines = []
    for path in sorted(output_dir.rglob("*")):
        if path.is_file() and path.name != "SHA256SUMS":
            lines.append(f"{sha256_file(path)}  {path.relative_to(output_dir).as_posix()}")
    (output_dir / "SHA256SUMS").write_text("\n".join(lines) + "\n")


def build_release(args: argparse.Namespace) -> dict[str, Any]:
    database = args.database.resolve()
    restoration_preview = args.restoration_preview.resolve()
    cache = args.model_cache.resolve()
    output = args.output_dir.resolve()
    archive = args.archive.resolve()
    if output.exists() or archive.exists():
        raise ReleaseError(f"release output already exists: {output} or {archive}")
    if not cache.is_dir():
        raise ReleaseError(f"model cache is missing: {cache}")
    source_hash_before = sha256_file(database)
    with connect_ro(database) as connection:
        counts = validate_database(connection)

    current_members = active_members(database)
    restored_accessions = restored_alphafold_accessions(restoration_preview)
    cache_delta = sorted(
        member
        for member, row in current_members.items()
        if str(row["tf_id"]) in restored_accessions
        and str(row["source"]).casefold() == "alphafold"
    )
    observed_by_accession = Counter(
        str(current_members[member]["tf_id"]) for member in cache_delta
    )
    if len(cache_delta) != 55 or dict(observed_by_accession) != restored_accessions:
        raise ReleaseError(
            "database members do not match the frozen restored AlphaFold set: "
            f"files={len(cache_delta)}, accessions={dict(observed_by_accession)}"
        )

    output.mkdir(parents=True)
    try:
        shutil.copy2(ROOT / "app.py", output / "app.py")
        shutil.copy2(ROOT / "requirements.txt", output / "requirements.txt")
        shutil.copytree(ROOT / "templates", output / "templates")
        shutil.copytree(ROOT / "static", output / "static")
        (output / "data").mkdir()
        shutil.copy2(ROOT / "data/external_motif_links.tsv", output / "data/external_motif_links.tsv")
        next_database = output / "data/tf_webdb.next.sqlite"
        shutil.copy2(database, next_database)
        artifact_count = rewrite_af3_paths(
            next_database,
            args.production_root,
            output / "AF3_EXPECTED_SHA256.tsv",
        )
        delta_root = output / "cache_delta"
        delta_bytes = 0
        for member in cache_delta:
            source = cache.joinpath(*PurePosixPath(member).parts)
            if not source.is_file():
                raise ReleaseError(f"cache delta source is missing: {source}")
            destination = delta_root.joinpath(*PurePosixPath(member).parts)
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, destination)
            delta_bytes += source.stat().st_size

        py_compile.compile(str(output / "app.py"), doraise=True)
        shutil.rmtree(output / "__pycache__", ignore_errors=True)
        metadata = {
            "built_at_utc": datetime.now(timezone.utc).isoformat(),
            "git_commit": git_commit(),
            "source_database": str(database),
            "source_database_sha256": source_hash_before,
            "release_database_sha256": sha256_file(next_database),
            "restoration_preview": str(restoration_preview),
            "restoration_preview_sha256": sha256_file(restoration_preview),
            "production_root": str(args.production_root),
            "counts": counts,
            "cache_delta_files": len(cache_delta),
            "cache_delta_bytes": delta_bytes,
            "af3_artifact_manifest_rows": artifact_count,
            "database_or_source_data_modified": False,
        }
        (output / "RELEASE_MANIFEST.json").write_text(
            json.dumps(metadata, indent=2, sort_keys=True) + "\n"
        )
        write_checksums(output)
        archive.parent.mkdir(parents=True, exist_ok=True)
        with tarfile.open(archive, "w:gz") as tar:
            tar.add(output, arcname=output.name)
    except BaseException:
        shutil.rmtree(output, ignore_errors=True)
        archive.unlink(missing_ok=True)
        raise
    if sha256_file(database) != source_hash_before:
        raise ReleaseError("source database changed while building release")
    metadata["release_archive"] = str(archive)
    metadata["release_archive_sha256"] = sha256_file(archive)
    return metadata


def main() -> int:
    args = parse_args()
    metadata = build_release(args)
    print("Release directory:", args.output_dir.resolve())
    print("Release archive:", metadata["release_archive"])
    print("Release archive SHA256:", metadata["release_archive_sha256"])
    print("Git commit:", metadata["git_commit"])
    print("TF rows:", metadata["counts"]["tf"])
    print("Active PDB rows:", metadata["counts"]["active_pdb"])
    print("Active AF3 assignments:", metadata["counts"]["active_af3"])
    print("Cache delta files:", metadata["cache_delta_files"])
    print("AF3 manifest rows:", metadata["af3_artifact_manifest_rows"])
    print("Deployed: NO")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ReleaseError as exc:
        raise SystemExit(f"ERROR: {exc}")
