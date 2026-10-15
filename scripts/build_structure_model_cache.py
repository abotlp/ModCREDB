#!/usr/bin/env python3
"""Build the filesystem cache used to serve active PDB models quickly.

The SQLite database remains the source of the exact active member list and the
tar.gz archives remain the recoverable coordinate source.  The cache mirrors
``member_path`` below one root so ``app.py`` can use it through
``TF_WEBDB_MODEL_CACHE_DIR`` without changing database records.

Dry-run is the default.  Apply mode refuses to overwrite a non-matching cache;
an already complete cache is validated and reported as an idempotent no-op.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import sqlite3
import tarfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath
from typing import BinaryIO


APP_DIR = Path(__file__).resolve().parents[1]
DEFAULT_DATABASE = APP_DIR / "data/tf_webdb_local_staging.sqlite"
DEFAULT_CACHE = Path("/home/patricia/TF_database_Baldo_data/modcredb_model_cache")
DEFAULT_QC = APP_DIR / "outputs/modcredb_model_cache_qc.json"


class CacheError(RuntimeError):
    pass


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build or validate the extracted cache for active structure_file PDBs."
    )
    parser.add_argument("--database", type=Path, default=DEFAULT_DATABASE)
    parser.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--qc-output", type=Path, default=DEFAULT_QC)
    parser.add_argument(
        "--resume-temp",
        type=Path,
        help="Validate and atomically activate an interrupted temporary cache",
    )
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--dry-run", action="store_true", help="Validate and plan only (default)")
    mode.add_argument("--apply", action="store_true", help="Create the cache atomically")
    parser.add_argument(
        "--sync-existing",
        action="store_true",
        help=(
            "With --apply, add only indexed files missing from an existing cache. "
            "Refuses caches containing extra files or wrong-size files."
        ),
    )
    parser.add_argument(
        "--replace-qc",
        action="store_true",
        help="Replace only the generated QC JSON; never replaces the cache",
    )
    return parser.parse_args()


def safe_member_path(value: str) -> PurePosixPath:
    member = PurePosixPath(value)
    if member.is_absolute() or not member.parts or any(part in {"", ".", ".."} for part in member.parts):
        raise CacheError(f"unsafe archive member path: {value!r}")
    return member


def connect_readonly(path: Path) -> sqlite3.Connection:
    if not path.is_file():
        raise CacheError(f"database not found: {path}")
    connection = sqlite3.connect(f"file:{path}?mode=ro&immutable=1", uri=True)
    connection.row_factory = sqlite3.Row
    return connection


def active_plan(path: Path) -> list[dict[str, object]]:
    connection = connect_readonly(path)
    quick_check = str(connection.execute("PRAGMA quick_check").fetchone()[0])
    if quick_check != "ok":
        raise CacheError(f"database quick_check failed: {quick_check}")
    rows = connection.execute(
        """
        SELECT id, tf_id, source, archive_path, member_path
        FROM structure_file
        WHERE status = 'active' AND file_type = 'pdb'
        ORDER BY archive_path, member_path, id
        """
    ).fetchall()
    connection.close()
    if not rows:
        raise CacheError("no active PDB rows found")

    by_output: dict[str, tuple[str, str]] = {}
    records: list[dict[str, object]] = []
    for row in rows:
        archive = str(Path(str(row["archive_path"])).resolve())
        member = safe_member_path(str(row["member_path"])).as_posix()
        previous = by_output.get(member)
        if previous is not None and previous != (archive, member):
            raise CacheError(
                f"cache member collision for {member}: {previous[0]} and {archive}"
            )
        by_output[member] = (archive, member)
        records.append(
            {
                "structure_file_id": int(row["id"]),
                "tf_id": str(row["tf_id"]),
                "source": str(row["source"]),
                "archive_path": archive,
                "member_path": member,
            }
        )
    if len(by_output) != len(records):
        raise CacheError("multiple active database rows point to the same cache member")
    return records


def resolve_archive_members(records: list[dict[str, object]]) -> tuple[dict[str, dict[str, int]], int]:
    wanted: dict[str, set[str]] = defaultdict(set)
    for row in records:
        wanted[str(row["archive_path"])].add(str(row["member_path"]))

    resolved: dict[str, dict[str, int]] = {}
    total_bytes = 0
    for archive_name, members in sorted(wanted.items()):
        archive_path = Path(archive_name)
        if not archive_path.is_file():
            raise CacheError(f"archive not found: {archive_path}")
        found: dict[str, int] = {}
        with tarfile.open(archive_path, "r:gz") as archive:
            for info in archive:
                if info.name in members:
                    if not info.isfile():
                        raise CacheError(f"active archive member is not a file: {info.name}")
                    if info.name in found:
                        raise CacheError(f"duplicate archive member: {archive_path}:{info.name}")
                    found[info.name] = int(info.size)
        missing = sorted(members - set(found))
        if missing:
            raise CacheError(
                f"{len(missing)} active members missing from {archive_path}: {missing[:5]}"
            )
        resolved[archive_name] = found
        total_bytes += sum(found.values())
    return resolved, total_bytes


def sha256_stream_and_copy(source: BinaryIO, destination: Path) -> str:
    digest = hashlib.sha256()
    with destination.open("wb") as output:
        while True:
            chunk = source.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
            output.write(chunk)
    return digest.hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def cache_mismatches(
    cache: Path,
    resolved: dict[str, dict[str, int]],
) -> tuple[list[str], list[str], list[str]]:
    expected = {
        member: size
        for archive_members in resolved.values()
        for member, size in archive_members.items()
    }
    if not cache.is_dir():
        return sorted(expected), [], []
    actual = {
        path.relative_to(cache).as_posix(): path.stat().st_size
        for path in cache.rglob("*")
        if path.is_file()
    }
    missing = sorted(set(expected) - set(actual))
    extra = sorted(set(actual) - set(expected))
    wrong_size = sorted(
        member for member in set(expected) & set(actual) if expected[member] != actual[member]
    )
    return missing, extra, wrong_size


def validate_cache(cache: Path, resolved: dict[str, dict[str, int]], hashes: dict[str, str] | None) -> None:
    missing, extra, wrong_size = cache_mismatches(cache, resolved)
    if missing or extra or wrong_size:
        raise CacheError(
            "cache mismatch: "
            f"missing={missing[:5]}, extra={extra[:5]}, wrong_size={wrong_size[:5]}"
        )
    if hashes is not None:
        for index, (member, expected_hash) in enumerate(sorted(hashes.items()), start=1):
            if sha256_file(cache / PurePosixPath(member)) != expected_hash:
                raise CacheError(f"cached content hash mismatch: {member}")
            if index % 5000 == 0:
                print(f"Validated cached PDBs: {index}/{len(hashes)}", flush=True)


def sync_missing_cache_files(
    cache: Path,
    resolved: dict[str, dict[str, int]],
    missing: list[str],
) -> tuple[int, int]:
    """Add exact missing archive members without replacing existing cache files.

    Each member is first written beside the cache in a temporary directory and
    then moved atomically into place.  An interrupted run is therefore safe to
    repeat.  Extra or wrong-size cache files are rejected by the caller.
    """
    wanted = set(missing)
    if not wanted:
        return 0, 0
    temporary = cache.with_name(f"{cache.name}.sync.{os.getpid()}")
    if temporary.exists():
        raise CacheError(f"temporary sync path already exists: {temporary}")
    temporary.mkdir(parents=True)
    created = 0
    copied_bytes = 0
    try:
        for archive_name, members in sorted(resolved.items()):
            remaining = wanted & set(members)
            if not remaining:
                continue
            with tarfile.open(archive_name, "r:gz") as archive:
                for info in archive:
                    if info.name not in remaining:
                        continue
                    source = archive.extractfile(info)
                    if source is None:
                        raise CacheError(f"cannot read archive member: {archive_name}:{info.name}")
                    staged = temporary.joinpath(*PurePosixPath(info.name).parts)
                    staged.parent.mkdir(parents=True, exist_ok=True)
                    sha256_stream_and_copy(source, staged)
                    if staged.stat().st_size != members[info.name]:
                        raise CacheError(f"extracted size mismatch: {info.name}")
                    destination = cache.joinpath(*PurePosixPath(info.name).parts)
                    destination.parent.mkdir(parents=True, exist_ok=True)
                    if destination.exists():
                        if destination.stat().st_size != members[info.name]:
                            raise CacheError(f"cache file appeared with wrong size: {info.name}")
                        staged.unlink()
                    else:
                        os.replace(staged, destination)
                        created += 1
                        copied_bytes += int(info.size)
                    remaining.remove(info.name)
            if remaining:
                raise CacheError(
                    f"missing members were not found during sync in {archive_name}: "
                    f"{sorted(remaining)[:5]}"
                )
        validate_cache(cache, resolved, hashes=None)
    finally:
        shutil.rmtree(temporary, ignore_errors=True)
    return created, copied_bytes


def archive_member_hashes(resolved: dict[str, dict[str, int]]) -> dict[str, str]:
    hashes: dict[str, str] = {}
    expected_count = sum(len(members) for members in resolved.values())
    for archive_name, members in sorted(resolved.items()):
        remaining = set(members)
        with tarfile.open(archive_name, "r:gz") as archive:
            for info in archive:
                if info.name not in remaining:
                    continue
                source = archive.extractfile(info)
                if source is None:
                    raise CacheError(f"cannot read archive member: {archive_name}:{info.name}")
                digest = hashlib.sha256()
                for chunk in iter(lambda: source.read(1024 * 1024), b""):
                    digest.update(chunk)
                hashes[info.name] = digest.hexdigest()
                remaining.remove(info.name)
                if len(hashes) % 5000 == 0:
                    print(f"Hashed source PDBs: {len(hashes)}/{expected_count}", flush=True)
        if remaining:
            raise CacheError(f"archive changed during source hashing: {archive_name}")
    if len(hashes) != expected_count:
        raise CacheError(f"source hash count mismatch: {len(hashes)} != {expected_count}")
    return hashes


def validate_resume_path(cache: Path, temporary: Path) -> Path:
    cache = cache.resolve()
    temporary = temporary.resolve()
    if temporary.parent != cache.parent or not temporary.name.startswith(f"{cache.name}.tmp."):
        raise CacheError(
            f"resume path must be an interrupted {cache.name}.tmp.* directory beside the cache"
        )
    if not temporary.is_dir():
        raise CacheError(f"resume directory not found: {temporary}")
    if cache.exists():
        raise CacheError(f"cannot resume because final cache already exists: {cache}")
    return temporary


def build_cache(
    cache: Path,
    resolved: dict[str, dict[str, int]],
    resume_temp: Path | None = None,
) -> tuple[int, int]:
    if cache.exists():
        if resume_temp is not None:
            raise CacheError("--resume-temp cannot be used when the final cache already exists")
        if not cache.is_dir():
            raise CacheError(f"cache target exists and is not a directory: {cache}")
        validate_cache(cache, resolved, hashes=None)
        return 0, 0

    if resume_temp is not None:
        temporary = validate_resume_path(cache, resume_temp)
        hashes = archive_member_hashes(resolved)
        validate_cache(temporary, resolved, hashes=hashes)
        temporary.rename(cache)
        return 0, 0

    temporary = cache.with_name(f"{cache.name}.tmp.{os.getpid()}")
    if temporary.exists():
        raise CacheError(f"temporary cache path already exists: {temporary}")
    temporary.mkdir(parents=True)
    hashes: dict[str, str] = {}
    created = 0
    copied_bytes = 0
    try:
        for archive_name, members in sorted(resolved.items()):
            remaining = set(members)
            with tarfile.open(archive_name, "r:gz") as archive:
                for info in archive:
                    if info.name not in remaining:
                        continue
                    source = archive.extractfile(info)
                    if source is None:
                        raise CacheError(f"cannot read archive member: {archive_name}:{info.name}")
                    destination = temporary.joinpath(*PurePosixPath(info.name).parts)
                    destination.parent.mkdir(parents=True, exist_ok=True)
                    hashes[info.name] = sha256_stream_and_copy(source, destination)
                    if destination.stat().st_size != members[info.name]:
                        raise CacheError(f"extracted size mismatch: {info.name}")
                    copied_bytes += destination.stat().st_size
                    created += 1
                    remaining.remove(info.name)
                    if created % 5000 == 0:
                        print(f"Cached PDBs: {created}", flush=True)
            if remaining:
                raise CacheError(f"archive changed during cache build: {archive_name}")
        validate_cache(temporary, resolved, hashes=hashes)
        temporary.rename(cache)
    except BaseException:
        shutil.rmtree(temporary, ignore_errors=True)
        raise
    return created, copied_bytes


def main() -> int:
    args = parse_args()
    if args.resume_temp is not None and not args.apply:
        raise CacheError("--resume-temp requires --apply")
    if args.sync_existing and not args.apply:
        raise CacheError("--sync-existing requires --apply")
    records = active_plan(args.database)
    resolved, total_bytes = resolve_archive_members(records)
    archive_count = len(resolved)
    source_counts = Counter(str(row["source"]).casefold() for row in records)
    usage = shutil.disk_usage(args.cache_dir.parent)
    required_with_margin = int(total_bytes * 1.05)
    if args.apply and not args.cache_dir.exists() and usage.free < required_with_margin:
        raise CacheError(
            f"insufficient free space: need at least {required_with_margin} bytes, have {usage.free}"
        )

    created = 0
    copied_bytes = 0
    missing_before, extra_before, wrong_size_before = cache_mismatches(
        args.cache_dir, resolved
    )
    expected_sizes = {
        member: size
        for archive_members in resolved.values()
        for member, size in archive_members.items()
    }
    missing_bytes = sum(expected_sizes[member] for member in missing_before)
    if (
        args.apply
        and args.cache_dir.exists()
        and args.sync_existing
        and usage.free < int(missing_bytes * 1.05)
    ):
        raise CacheError(
            "insufficient free space for cache sync: "
            f"need at least {int(missing_bytes * 1.05)} bytes, have {usage.free}"
        )
    if args.apply:
        if args.cache_dir.exists() and args.sync_existing:
            if extra_before or wrong_size_before:
                raise CacheError(
                    "refusing to sync unsafe existing cache: "
                    f"extra={extra_before[:5]}, wrong_size={wrong_size_before[:5]}"
                )
            created, copied_bytes = sync_missing_cache_files(
                args.cache_dir, resolved, missing_before
            )
        else:
            created, copied_bytes = build_cache(args.cache_dir, resolved, args.resume_temp)

    missing_after, extra_after, wrong_size_after = cache_mismatches(
        args.cache_dir, resolved
    )

    qc = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "mode": "APPLY" if args.apply else "DRY_RUN",
        "database": str(args.database.resolve()),
        "cache_dir": str(args.cache_dir.resolve()),
        "active_pdb_rows": len(records),
        "unique_cache_members": sum(len(value) for value in resolved.values()),
        "archives": archive_count,
        "source_counts": dict(sorted(source_counts.items())),
        "planned_uncompressed_bytes": total_bytes,
        "free_bytes_before": usage.free,
        "files_created": created,
        "bytes_copied": copied_bytes,
        "missing_files_before": len(missing_before),
        "missing_bytes_before": missing_bytes,
        "extra_files_before": len(extra_before),
        "wrong_size_files_before": len(wrong_size_before),
        "missing_files_after": len(missing_after),
        "extra_files_after": len(extra_after),
        "wrong_size_files_after": len(wrong_size_after),
        "resumed_temporary_cache": str(args.resume_temp.resolve()) if args.resume_temp else "",
        "cache_exists": args.cache_dir.is_dir(),
        "all_invariants_passed": not (missing_after or extra_after or wrong_size_after),
    }
    if args.qc_output.exists() and not args.replace_qc:
        raise CacheError(f"QC output exists: {args.qc_output}; use --replace-qc")
    args.qc_output.parent.mkdir(parents=True, exist_ok=True)
    args.qc_output.write_text(json.dumps(qc, indent=2, sort_keys=True) + "\n")

    print("Mode:", qc["mode"])
    print("Active PDB rows:", qc["active_pdb_rows"])
    print("Archive sources:", archive_count)
    print("Uncompressed cache bytes:", total_bytes)
    print("Source counts:", dict(sorted(source_counts.items())))
    print("Files created:", created)
    print("Missing cache files before:", len(missing_before))
    print("Missing cache files after:", len(missing_after))
    print("Cache ready:", "YES" if args.cache_dir.is_dir() else "NO (dry run)")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except CacheError as exc:
        raise SystemExit(f"ERROR: {exc}")
