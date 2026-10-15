#!/usr/bin/env python3
"""Restore the complete 5,417-row manuscript TF universe to local staging.

The delivered TSV is immutable membership/provenance evidence.  Most omitted
records are restored from the historical pre-TARDBP-removal SQLite snapshot.
Q96T58, the only delivered accession absent from that snapshot's ``tf`` table,
uses a frozen official UniProt annotation recorded in this script.

The 11 delivered direct-PWM rows referenced withdrawn JASPAR records.  Their
replacement evidence and sequence-transfer decisions are frozen in
``data_sources/manuscript_missing_annotated_restoration_2026-10-14.tsv``.
The 10 TARDBP records recover only the exact 55 AlphaFold motif/model tokens in
the delivered TSV, not additional historical candidates.  Each recovered PDB
is re-audited with the standard 4.5 A protein-DNA heavy-atom rule.

Dry-run is the default safe workflow and operates on an in-memory SQLite copy.
Apply mode first creates and verifies a timestamped byte-for-byte backup.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import shutil
import sqlite3
import tarfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATABASE = ROOT / "data/tf_webdb_local_staging.sqlite"
HISTORICAL_DATABASE = ROOT / "data/tf_webdb.before_remove_TARDBP.sqlite"
INTEGRATED_TSV = Path(
    "/home/patricia/TF_database_Baldo_data/"
    "TF_PWM_chart_final_integrated_HOCOMOCO_hierarchical.tsv"
)
DECISION_MANIFEST = (
    ROOT / "data_sources/manuscript_missing_annotated_restoration_2026-10-14.tsv"
)
DEFAULT_MODEL_CACHE = Path(
    "/home/patricia/TF_database_Baldo_data/modcredb_model_cache"
)
DEFAULT_PREVIEW = ROOT / "outputs/manuscript_5417_restoration_preview.tsv"
DEFAULT_QC = ROOT / "outputs/manuscript_5417_restoration_qc.json"
DEFAULT_BACKUPS = ROOT / "backups"

INTEGRATED_SHA256 = "f71160a2c2d318512be8354302208b55b4d4e7b427178631d98977c4209ed845"
EXPECTED_TFS = 5417
EXPECTED_INITIAL_MISSING = 150
EXPECTED_MISSING_BREAKDOWN = Counter(
    {"Unannotated": 129, "Identical_PWM": 11, "AlphaFold": 10}
)
EXPECTED_TARDBP_TFS = 10
EXPECTED_TARDBP_MODELS = 55
EXPECTED_TARDBP_INTERFACE_YES = 54
EXPECTED_TARDBP_INTERFACE_NO = 1
EXPECTED_CORRECTED_LEVELS = Counter({"Identical_PWM": 2, "Homologous_PWM": 9})
CONTACT_CUTOFF = 4.5

TARDBP_ACCESSIONS = frozenset(
    {
        "A0A087WV68", "A0A087WXQ5", "A0A087WYY0", "A0A087X260",
        "A0A0A0N0L3", "A0A0A0N0M3", "B1AKP7", "B4DRW3", "G3V162",
        "Q13148",
    }
)

PROTEIN_RESIDUES = frozenset(
    "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR "
    "TRP TYR VAL ASX GLX HYP MSE PYL SEC XLE".split()
)
NUCLEIC_RESIDUES = frozenset(
    "DA DC DG DT DI A C G T U I ADE CYT GUA THY URI".split()
)

PREVIEW_COLUMNS = (
    "tf_id", "delivered_primary_class", "restored_primary_class", "gene",
    "record_action", "motif_refs_planned", "active_models_planned", "reason",
)

PROTECTED_TABLES = (
    "source", "source_release", "motif_file", "model_summary",
    "structure_model_assignment", "structure_confidence_artifact",
    "tf_structure_status",
)

Q96T58_ANNOTATION = {
    "uniprot_accession": "Q96T58",
    "entry_name": "MINT_HUMAN",
    "gene_names": "SPEN KIAA0929 MINT SHARP",
    "protein_name": "Msx2-interacting protein",
    "organism_name": "Homo sapiens (Human)",
    "organism_id": 9606,
    "reviewed": 1,
    "sequence_length": 3664,
    "annotation_score": 5.0,
    "uniprot_url": "https://www.uniprot.org/uniprotkb/Q96T58/entry",
    "fetched_at": "2026-10-14T00:00:00+00:00",
}


class RestorationError(RuntimeError):
    pass


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--dry-run", action="store_true")
    mode.add_argument("--apply", action="store_true")
    parser.add_argument("--database", type=Path, default=DEFAULT_DATABASE)
    parser.add_argument("--historical-database", type=Path, default=HISTORICAL_DATABASE)
    parser.add_argument("--integrated-tsv", type=Path, default=INTEGRATED_TSV)
    parser.add_argument("--decision-manifest", type=Path, default=DECISION_MANIFEST)
    parser.add_argument("--model-cache", type=Path, default=DEFAULT_MODEL_CACHE)
    parser.add_argument("--preview-output", type=Path, default=DEFAULT_PREVIEW)
    parser.add_argument("--qc-output", type=Path, default=DEFAULT_QC)
    parser.add_argument("--backup-dir", type=Path, default=DEFAULT_BACKUPS)
    parser.add_argument("--replace-outputs", action="store_true")
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def connect_ro(path: Path) -> sqlite3.Connection:
    connection = sqlite3.connect(f"file:{path}?mode=ro&immutable=1", uri=True)
    connection.row_factory = sqlite3.Row
    return connection


def in_memory_copy(path: Path) -> sqlite3.Connection:
    source = connect_ro(path)
    destination = sqlite3.connect(":memory:")
    destination.row_factory = sqlite3.Row
    source.backup(destination)
    source.close()
    destination.execute("PRAGMA foreign_keys=ON")
    return destination


def table_columns(connection: sqlite3.Connection, table: str) -> list[str]:
    return [str(row[1]) for row in connection.execute(f"PRAGMA table_info({table})")]


def encoded(value: Any) -> bytes:
    if value is None:
        return b"N"
    payload = value if isinstance(value, bytes) else str(value).encode("utf-8")
    return (b"B" if isinstance(value, bytes) else b"T") + str(len(payload)).encode() + b":" + payload


def table_snapshot(connection: sqlite3.Connection, table: str) -> tuple[int, str]:
    columns = table_columns(connection, table)
    if not columns:
        return (0, "ABSENT")
    quoted = ",".join('"' + item.replace('"', '""') + '"' for item in columns)
    digest = hashlib.sha256()
    count = 0
    for row in connection.execute(f"SELECT {quoted} FROM {table} ORDER BY {quoted}"):
        count += 1
        for value in row:
            digest.update(encoded(value)); digest.update(b"|")
        digest.update(b"\n")
    return count, digest.hexdigest()


def snapshots(connection: sqlite3.Connection) -> dict[str, tuple[int, str]]:
    return {table: table_snapshot(connection, table) for table in PROTECTED_TABLES}


def is_hydrogen(line: str, atom_name: str) -> bool:
    element = line[76:78].strip().upper() if len(line) >= 78 else ""
    if element in {"H", "D"}:
        return True
    return not element and atom_name.lstrip("0123456789").upper().startswith(("H", "D"))


def audit_pdb_lines(lines: Iterable[str]) -> dict[str, Any]:
    protein: list[tuple[float, float, float]] = []
    dna: list[tuple[float, float, float]] = []
    protein_chains: set[str] = set()
    dna_chains: set[str] = set()
    for line in lines:
        if not line.startswith(("ATOM  ", "HETATM")) or len(line) < 54:
            continue
        if line[16:17] not in {" ", "A"}:
            continue
        residue = line[17:20].strip().upper()
        if residue not in PROTEIN_RESIDUES and residue not in NUCLEIC_RESIDUES:
            continue
        atom_name = line[12:16].strip().upper()
        if is_hydrogen(line, atom_name):
            continue
        try:
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError:
            continue
        chain = line[21:22].strip() or "_"
        if residue in PROTEIN_RESIDUES:
            protein.append(xyz); protein_chains.add(chain)
        else:
            dna.append(xyz); dna_chains.add(chain)

    if not protein:
        return {"protein_chain_count": 0, "dna_chain_count": len(dna_chains),
                "dna_interaction": None, "dna_contact_count": 0,
                "dna_qc_status": "no_protein_atoms"}
    if not dna:
        return {"protein_chain_count": len(protein_chains), "dna_chain_count": 0,
                "dna_interaction": 0, "dna_contact_count": 0,
                "dna_qc_status": "no_nucleic_atoms"}

    grid: dict[tuple[int, int, int], list[tuple[float, float, float]]] = defaultdict(list)
    for x, y, z in dna:
        grid[(math.floor(x / CONTACT_CUTOFF), math.floor(y / CONTACT_CUTOFF),
              math.floor(z / CONTACT_CUTOFF))].append((x, y, z))
    cutoff2 = CONTACT_CUTOFF * CONTACT_CUTOFF
    contacts = 0
    for x, y, z in protein:
        gx, gy, gz = (math.floor(x / CONTACT_CUTOFF), math.floor(y / CONTACT_CUTOFF),
                      math.floor(z / CONTACT_CUTOFF))
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for nx, ny, nz in grid.get((gx + dx, gy + dy, gz + dz), ()):
                        if (x - nx) ** 2 + (y - ny) ** 2 + (z - nz) ** 2 <= cutoff2:
                            contacts += 1
    return {"protein_chain_count": len(protein_chains),
            "dna_chain_count": len(dna_chains),
            "dna_interaction": 1 if contacts else 0,
            "dna_contact_count": contacts, "dna_qc_status": "ok"}


def audit_member(
    member_path: str,
    archive_path: str,
    model_cache: Path,
    archives: dict[Path, tarfile.TarFile],
) -> dict[str, Any]:
    extracted = model_cache / member_path
    if extracted.is_file():
        with extracted.open(errors="replace") as handle:
            return audit_pdb_lines(handle)
    archive_file = Path(archive_path)
    if not archive_file.is_file():
        raise RestorationError(f"missing coordinate archive: {archive_file}")
    archive = archives.get(archive_file)
    if archive is None:
        archive = tarfile.open(archive_file, "r:gz")
        archives[archive_file] = archive
    try:
        raw = archive.extractfile(member_path)
    except KeyError as exc:
        raise RestorationError(f"missing coordinate member: {member_path}") from exc
    if raw is None:
        raise RestorationError(f"unreadable coordinate member: {member_path}")
    with raw, io.TextIOWrapper(raw, encoding="ascii", errors="replace") as handle:
        return audit_pdb_lines(handle)


def validate_inputs(args: argparse.Namespace) -> tuple[list[dict[str, str]], dict[str, dict[str, str]]]:
    for path in (args.database, args.historical_database, args.integrated_tsv,
                 args.decision_manifest):
        if not path.is_file():
            raise RestorationError(f"missing input: {path}")
    if sha256_file(args.integrated_tsv) != INTEGRATED_SHA256:
        raise RestorationError("immutable integrated TSV SHA256 mismatch")
    delivered = read_tsv(args.integrated_tsv)
    if len(delivered) != EXPECTED_TFS or len({r["TF_name"] for r in delivered}) != EXPECTED_TFS:
        raise RestorationError("integrated TSV is not the expected 5,417 unique TF universe")
    decisions_list = read_tsv(args.decision_manifest)
    decisions = {row["tf_id"]: row for row in decisions_list}
    if len(decisions) != 11 or len(decisions_list) != 11:
        raise RestorationError("decision manifest must contain 11 unique TFs")
    if Counter(row["corrected_primary_class"] for row in decisions_list) != EXPECTED_CORRECTED_LEVELS:
        raise RestorationError("unexpected corrected-level breakdown in decision manifest")
    delivered_by_id = {row["TF_name"]: row for row in delivered}
    for accession, decision in decisions.items():
        if delivered_by_id.get(accession, {}).get("Best_annotation_level") != "Identical_PWM":
            raise RestorationError(f"decision accession is not delivered Identical_PWM: {accession}")
        if delivered_by_id[accession]["Best_PWM_or_model"] != decision["withdrawn_original_motif"]:
            raise RestorationError(f"withdrawn motif mismatch: {accession}")
    delivered_tar = {r["TF_name"] for r in delivered if r["TF_name"] in TARDBP_ACCESSIONS}
    if delivered_tar != TARDBP_ACCESSIONS:
        raise RestorationError("delivered TARDBP set mismatch")
    return delivered, decisions


def row_dict(row: sqlite3.Row) -> dict[str, Any]:
    return {key: row[key] for key in row.keys()}


def insert_tf_records(
    connection: sqlite3.Connection,
    historical: sqlite3.Connection,
    delivered: list[dict[str, str]],
    integrated_path: Path,
) -> tuple[set[str], Counter[str]]:
    current = {str(row[0]) for row in connection.execute("SELECT tf_id FROM tf")}
    missing_rows = [row for row in delivered if row["TF_name"] not in current]
    breakdown = Counter(row["Best_annotation_level"] or "Unannotated" for row in missing_rows)
    if len(missing_rows) not in {0, EXPECTED_INITIAL_MISSING}:
        raise RestorationError(f"unexpected pre-restoration missing count: {len(missing_rows)}")
    if missing_rows and breakdown != EXPECTED_MISSING_BREAKDOWN:
        raise RestorationError(f"unexpected missing breakdown: {dict(breakdown)}")

    inserted: set[str] = set()
    for delivered_row in missing_rows:
        accession = delivered_row["TF_name"]
        old_tf = historical.execute("SELECT * FROM tf WHERE tf_id=?", (accession,)).fetchone()
        family_text = delivered_row["TF_family"]
        if old_tf is not None and str(old_tf["family_text"]) != family_text:
            raise RestorationError(f"family text differs from delivered TSV: {accession}")
        connection.execute(
            "INSERT INTO tf(tf_id,family_text,motif_ref_count,active_model_count) VALUES(?,?,0,0)",
            (accession, family_text),
        )
        old_annotation = historical.execute(
            "SELECT * FROM tf_annotation WHERE tf_id=?", (accession,)
        ).fetchone()
        if old_annotation is None:
            if accession != "Q96T58":
                raise RestorationError(f"historical annotation missing unexpectedly: {accession}")
            annotation = {"tf_id": accession, **Q96T58_ANNOTATION}
        else:
            annotation = row_dict(old_annotation)
        columns = table_columns(connection, "tf_annotation")
        connection.execute(
            f"INSERT INTO tf_annotation({','.join(columns)}) VALUES({','.join('?' for _ in columns)})",
            tuple(annotation[column] for column in columns),
        )
        primary = connection.execute(
            "SELECT 1 FROM tf_primary_annotation WHERE tf_id=?", (accession,)
        ).fetchone()
        if primary is None:
            connection.execute(
                """INSERT INTO tf_primary_annotation
                   (tf_id,best_annotation_level,best_pwm_or_model,
                    n_nonempty_annotation_columns,source_table)
                   VALUES(?,?,?,?,?)""",
                (accession, delivered_row["Best_annotation_level"],
                 delivered_row["Best_PWM_or_model"],
                 int(delivered_row["N_nonempty_annotation_columns"] or 0),
                 str(integrated_path.resolve())),
            )
        inserted.add(accession)
    return inserted, breakdown


def motif_natural_id(
    connection: sqlite3.Connection, tf_id: str, source: str, motif_id: str
) -> int | None:
    rows = connection.execute(
        "SELECT id FROM motif_ref WHERE tf_id=? AND source=? AND motif_id=? ORDER BY id",
        (tf_id, source, motif_id),
    ).fetchall()
    if len(rows) > 1:
        raise RestorationError(f"duplicate motif natural key: {tf_id}/{source}/{motif_id}")
    return int(rows[0][0]) if rows else None


def structure_natural_id(
    connection: sqlite3.Connection, tf_id: str, source: str, model_id: str
) -> int | None:
    rows = connection.execute(
        """SELECT id FROM structure_file
           WHERE tf_id=? AND LOWER(source)=LOWER(?) AND model_id=? AND file_type='pdb'
           ORDER BY id""",
        (tf_id, source, model_id),
    ).fetchall()
    if len(rows) > 1:
        raise RestorationError(f"duplicate structure natural key: {tf_id}/{source}/{model_id}")
    return int(rows[0][0]) if rows else None


def restore_corrected_pwm_evidence(
    connection: sqlite3.Connection,
    decisions: dict[str, dict[str, str]],
    decision_manifest: Path,
) -> tuple[int, int]:
    inserted = 0
    primary_updates = 0
    for accession, decision in sorted(decisions.items()):
        source = decision["replacement_source"]
        canonical = decision["canonical_accession"]
        motif_ids = [item for item in decision["replacement_motif_ids"].split(";") if item]
        level = decision["corrected_primary_class"]
        evidence_type = "identical" if level == "Identical_PWM" else "homologous"
        identity = None if evidence_type == "identical" else float(decision["identity_percent"])
        for motif_id in motif_ids:
            if connection.execute(
                "SELECT 1 FROM motif_file WHERE source=? AND motif_id=?", (source, motif_id)
            ).fetchone() is None:
                raise RestorationError(f"replacement motif file absent: {source}/{motif_id}")
            canonical_row = connection.execute(
                """SELECT * FROM motif_ref WHERE tf_id=? AND source=? AND motif_id=?
                   ORDER BY id LIMIT 1""", (canonical, source, motif_id)
            ).fetchone()
            if canonical_row is None:
                raise RestorationError(f"canonical motif evidence absent: {canonical}/{source}/{motif_id}")
            if motif_natural_id(connection, accession, source, motif_id) is not None:
                continue
            relationship = decision["sequence_relationship"]
            mapping = (
                "identical_sequence_to_reviewed_canonical"
                if evidence_type == "identical" else "same_gene_homologous_sequence_transfer"
            )
            curation = "verified_same_gene_sequence_transfer"
            note = (
                f"Transferred from current direct {source.upper()} evidence for reviewed "
                f"{decision['gene']} accession {canonical}; local sequence audit: "
                f"{relationship}, {decision['aligned_residues']}/{decision['canonical_length']} "
                f"aligned residues, {decision['identity_percent']}% identity; "
                f"{decision['dbd_evidence']}. The delivered {decision['withdrawn_original_motif']} "
                "assignment was not reactivated because that JASPAR matrix is withdrawn."
            )
            connection.execute(
                """INSERT INTO motif_ref
                   (tf_id,evidence_type,source,motif_id,original_value,identity_percent,
                    missing_local_file,original_column,mapping_type,curation_status,
                    evidence_note,display_priority)
                   VALUES(?,?,?,?,?,?,?,?,?,?,?,?)""",
                (accession, evidence_type, source, motif_id, motif_id, identity, 0,
                 "manuscript_missing_record_restoration", mapping, curation, note,
                 10 if source == "jaspar" else 15),
            )
            inserted += 1
        desired_pwm = ";".join(motif_ids)
        current = connection.execute(
            "SELECT best_annotation_level,best_pwm_or_model FROM tf_primary_annotation WHERE tf_id=?",
            (accession,),
        ).fetchone()
        if current is None:
            raise RestorationError(f"primary row missing after TF restoration: {accession}")
        if (current[0], current[1]) != (level, desired_pwm):
            connection.execute(
                """UPDATE tf_primary_annotation
                   SET best_annotation_level=?, best_pwm_or_model=?,
                       n_nonempty_annotation_columns=1, source_table=?
                   WHERE tf_id=?""",
                (level, desired_pwm, str(decision_manifest.resolve()), accession),
            )
            primary_updates += 1
    return inserted, primary_updates


def restore_tardbp(
    connection: sqlite3.Connection,
    historical: sqlite3.Connection,
    delivered_by_id: dict[str, dict[str, str]],
    model_cache: Path,
    checked_at: str,
) -> tuple[int, int, int, Counter[int]]:
    motif_inserts = structure_inserts = link_inserts = 0
    interactions: Counter[int] = Counter()
    archives: dict[Path, tarfile.TarFile] = {}
    try:
        desired: list[tuple[str, str]] = []
        for accession in sorted(TARDBP_ACCESSIONS):
            row = delivered_by_id[accession]
            if row["Best_annotation_level"] != "AlphaFold":
                raise RestorationError(f"TARDBP row is not delivered AlphaFold: {accession}")
            tokens = [token.removesuffix(".meme") for token in row["AlphaFold"].split(";") if token]
            desired.extend((accession, token) for token in tokens)
        if len(desired) != EXPECTED_TARDBP_MODELS or len(set(desired)) != EXPECTED_TARDBP_MODELS:
            raise RestorationError("delivered TARDBP token set is not the expected 55 unique pairs")

        for accession, model_id in desired:
            old_motif = historical.execute(
                """SELECT * FROM motif_ref
                   WHERE tf_id=? AND source='alphafold' AND motif_id=?""",
                (accession, model_id),
            ).fetchall()
            old_structure = historical.execute(
                """SELECT * FROM structure_file
                   WHERE tf_id=? AND source='alphafold' AND model_id=? AND file_type='pdb'""",
                (accession, model_id),
            ).fetchall()
            if len(old_motif) != 1 or len(old_structure) != 1:
                raise RestorationError(f"historical TARDBP pair not unique: {accession}/{model_id}")
            motif_id = motif_natural_id(connection, accession, "alphafold", model_id)
            if motif_id is None:
                old = old_motif[0]
                columns = [c for c in table_columns(connection, "motif_ref") if c != "id"]
                connection.execute(
                    f"INSERT INTO motif_ref({','.join(columns)}) VALUES({','.join('?' for _ in columns)})",
                    tuple(old[c] for c in columns),
                )
                motif_id = int(connection.execute("SELECT last_insert_rowid()").fetchone()[0])
                motif_inserts += 1

            old = old_structure[0]
            qc = audit_member(str(old["member_path"]), str(old["archive_path"]),
                              model_cache, archives)
            if qc["dna_qc_status"] != "ok":
                raise RestorationError(f"TARDBP PDB is not assessable: {accession}/{model_id}")
            interactions[int(qc["dna_interaction"])] += 1
            structure_id = structure_natural_id(connection, accession, "alphafold", model_id)
            expected = {
                "source": old["source"], "model_id": old["model_id"], "tf_id": old["tf_id"],
                "member_path": old["member_path"], "archive_path": old["archive_path"],
                "file_type": old["file_type"], "status": "active",
                "template_pdb": old["template_pdb"], "residue_start": old["residue_start"],
                "residue_end": old["residue_end"], "model_kind": None,
                "protein_chain_count": qc["protein_chain_count"],
                "dna_chain_count": qc["dna_chain_count"],
                "dna_interaction": qc["dna_interaction"],
                "dna_contact_count": qc["dna_contact_count"],
                "dna_qc_cutoff": CONTACT_CUTOFF, "dna_qc_status": qc["dna_qc_status"],
                "dna_qc_checked_at": checked_at,
            }
            if structure_id is None:
                columns = [c for c in table_columns(connection, "structure_file") if c != "id"]
                connection.execute(
                    f"INSERT INTO structure_file({','.join(columns)}) VALUES({','.join('?' for _ in columns)})",
                    tuple(expected.get(c) for c in columns),
                )
                structure_id = int(connection.execute("SELECT last_insert_rowid()").fetchone()[0])
                structure_inserts += 1
            else:
                existing = connection.execute(
                    "SELECT * FROM structure_file WHERE id=?", (structure_id,)
                ).fetchone()
                compare = {key: existing[key] for key in expected}
                expected_compare = dict(expected)
                # A prior successful run has a valid earlier audit timestamp.
                expected_compare["dna_qc_checked_at"] = compare["dna_qc_checked_at"]
                if compare != expected_compare:
                    raise RestorationError(f"existing TARDBP model differs from frozen row: {model_id}")
            exists = connection.execute(
                "SELECT 1 FROM motif_structure WHERE motif_ref_id=? AND structure_file_id=?",
                (motif_id, structure_id),
            ).fetchone()
            if exists is None:
                connection.execute(
                    "INSERT INTO motif_structure(motif_ref_id,structure_file_id) VALUES(?,?)",
                    (motif_id, structure_id),
                )
                link_inserts += 1
    finally:
        for archive in archives.values():
            archive.close()
    return motif_inserts, structure_inserts, link_inserts, interactions


def update_counts(connection: sqlite3.Connection, accessions: set[str]) -> int:
    changed = 0
    for accession in sorted(accessions):
        motif_count = int(connection.execute(
            "SELECT COUNT(*) FROM motif_ref WHERE tf_id=?", (accession,)
        ).fetchone()[0])
        model_count = int(connection.execute(
            """SELECT COUNT(*) FROM structure_file
               WHERE tf_id=? AND status='active' AND file_type='pdb'""", (accession,)
        ).fetchone()[0])
        current = connection.execute(
            "SELECT motif_ref_count,active_model_count FROM tf WHERE tf_id=?", (accession,)
        ).fetchone()
        if current is None:
            raise RestorationError(f"TF row missing while updating counts: {accession}")
        if tuple(current) != (motif_count, model_count):
            connection.execute(
                "UPDATE tf SET motif_ref_count=?,active_model_count=? WHERE tf_id=?",
                (motif_count, model_count, accession),
            )
            changed += 1
    return changed


def validate_final(
    connection: sqlite3.Connection,
    delivered: list[dict[str, str]],
    decisions: dict[str, dict[str, str]],
) -> dict[str, Any]:
    delivered_ids = {row["TF_name"] for row in delivered}
    final_ids = {str(row[0]) for row in connection.execute("SELECT tf_id FROM tf")}
    if final_ids != delivered_ids:
        raise RestorationError(
            f"final TF universe mismatch: missing={len(delivered_ids-final_ids)}, extra={len(final_ids-delivered_ids)}"
        )
    if connection.execute("SELECT COUNT(*) FROM tf_annotation").fetchone()[0] != EXPECTED_TFS:
        raise RestorationError("final tf_annotation count is not 5,417")
    levels = Counter()
    for accession in decisions:
        row = connection.execute(
            "SELECT best_annotation_level FROM tf_primary_annotation WHERE tf_id=?", (accession,)
        ).fetchone()
        if row is None:
            raise RestorationError(f"missing corrected primary row: {accession}")
        levels[str(row[0])] += 1
    if levels != EXPECTED_CORRECTED_LEVELS:
        raise RestorationError(f"corrected primary breakdown mismatch: {dict(levels)}")

    placeholders = ",".join("?" for _ in TARDBP_ACCESSIONS)
    models = connection.execute(
        f"""SELECT COUNT(*),SUM(dna_interaction=1),SUM(dna_interaction=0)
            FROM structure_file WHERE tf_id IN ({placeholders})
              AND source='alphafold' AND file_type='pdb' AND status='active'""",
        tuple(sorted(TARDBP_ACCESSIONS)),
    ).fetchone()
    model_counts = tuple(int(value or 0) for value in models)
    if model_counts != (EXPECTED_TARDBP_MODELS, EXPECTED_TARDBP_INTERFACE_YES,
                        EXPECTED_TARDBP_INTERFACE_NO):
        raise RestorationError(f"TARDBP model/QC counts mismatch: {model_counts}")
    motif_count = int(connection.execute(
        f"""SELECT COUNT(*) FROM motif_ref WHERE tf_id IN ({placeholders})
            AND source='alphafold'""", tuple(sorted(TARDBP_ACCESSIONS))
    ).fetchone()[0])
    link_count = int(connection.execute(
        f"""SELECT COUNT(*) FROM motif_structure ms
            JOIN motif_ref mr ON mr.id=ms.motif_ref_id
            WHERE mr.tf_id IN ({placeholders}) AND mr.source='alphafold'""",
        tuple(sorted(TARDBP_ACCESSIONS)),
    ).fetchone()[0])
    if (motif_count, link_count) != (EXPECTED_TARDBP_MODELS, EXPECTED_TARDBP_MODELS):
        raise RestorationError(f"TARDBP motif/link counts mismatch: {motif_count}/{link_count}")
    fk = list(connection.execute("PRAGMA foreign_key_check"))
    fk_tables = Counter(str(row[0]) for row in fk)
    if fk_tables.get("tf_family", 0) or fk_tables.get("tf_primary_annotation", 0):
        raise RestorationError(f"TF-parent foreign key violations remain: {dict(fk_tables)}")
    return {
        "final_tf_rows": len(final_ids), "final_tf_annotation_rows": EXPECTED_TFS,
        "corrected_primary_levels": dict(levels), "tardbp_motif_refs": motif_count,
        "tardbp_active_models": model_counts[0], "tardbp_interface_yes": model_counts[1],
        "tardbp_interface_no": model_counts[2], "tardbp_links": link_count,
        "remaining_foreign_key_violations": len(fk),
        "remaining_foreign_key_violation_tables": dict(fk_tables),
    }


def write_preview(
    path: Path,
    delivered: list[dict[str, str]],
    decisions: dict[str, dict[str, str]],
    initial_missing: set[str],
    replace: bool,
) -> None:
    if path.exists() and not replace:
        raise RestorationError(f"preview exists; use --replace-outputs: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    mode = "w" if replace else "x"
    with path.open(mode, newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=PREVIEW_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in delivered:
            accession = row["TF_name"]
            if accession not in initial_missing:
                continue
            decision = decisions.get(accession)
            if decision:
                corrected = decision["corrected_primary_class"]
                motifs = len([x for x in decision["replacement_motif_ids"].split(";") if x])
                reason = (
                    f"Restore omitted manuscript TF; replace withdrawn {decision['withdrawn_original_motif']} "
                    f"with current {decision['replacement_source']} evidence from {decision['canonical_accession']}."
                )
                models = 0
            elif accession in TARDBP_ACCESSIONS:
                corrected = row["Best_annotation_level"]
                motifs = len([x for x in row["AlphaFold"].split(";") if x])
                models = motifs
                reason = "Restore exact delivered AlphaFold PWM/model pairs and current DNA-interface QC."
            else:
                corrected = row["Best_annotation_level"] or "Unannotated"
                motifs = models = 0
                reason = "Restore delivered manuscript-universe record without inventing annotation evidence."
            writer.writerow({
                "tf_id": accession, "delivered_primary_class": row["Best_annotation_level"] or "Unannotated",
                "restored_primary_class": corrected, "gene": decision["gene"] if decision else "",
                "record_action": "INSERT" if accession in initial_missing else "UNCHANGED",
                "motif_refs_planned": motifs, "active_models_planned": models, "reason": reason,
            })


def backup_database(database: Path, backup_dir: Path) -> tuple[Path, str]:
    backup_dir.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    backup = backup_dir / f"{database.name}.before_manuscript_5417_restore_{timestamp}.bak"
    if backup.exists():
        raise RestorationError(f"backup already exists: {backup}")
    before = sha256_file(database)
    shutil.copy2(database, backup)
    if sha256_file(backup) != before:
        backup.unlink(missing_ok=True)
        raise RestorationError("backup SHA256 mismatch")
    return backup, before


def execute_restoration(
    connection: sqlite3.Connection,
    historical: sqlite3.Connection,
    delivered: list[dict[str, str]],
    decisions: dict[str, dict[str, str]],
    model_cache: Path,
    integrated_path: Path,
    decision_manifest: Path,
) -> dict[str, Any]:
    initial_ids = {str(row[0]) for row in connection.execute("SELECT tf_id FROM tf")}
    delivered_ids = {row["TF_name"] for row in delivered}
    initial_missing = delivered_ids - initial_ids
    protected_before = snapshots(connection)
    checked_at = datetime.now(timezone.utc).isoformat(timespec="seconds")
    inserted, initial_breakdown = insert_tf_records(
        connection, historical, delivered, integrated_path
    )
    motif_evidence_inserts, primary_updates = restore_corrected_pwm_evidence(
        connection, decisions, decision_manifest
    )
    tar_motifs, tar_models, tar_links, tar_interactions = restore_tardbp(
        connection, historical, {r["TF_name"]: r for r in delivered}, model_cache, checked_at
    )
    # Only the restored records are in scope.  Existing TF counters are left
    # untouched because structural/PWM reconciliation outside this 150-row
    # restoration is a separate operation.
    count_updates = update_counts(connection, inserted)
    final = validate_final(connection, delivered, decisions)
    protected_after = snapshots(connection)
    if protected_before != protected_after:
        changed = [table for table in PROTECTED_TABLES if protected_before[table] != protected_after[table]]
        raise RestorationError(f"protected tables changed: {', '.join(changed)}")
    return {
        "initial_tf_rows": len(initial_ids), "initial_missing": len(initial_missing),
        "initial_missing_breakdown": dict(initial_breakdown),
        "tf_rows_inserted": len(inserted), "corrected_motif_refs_inserted": motif_evidence_inserts,
        "primary_rows_updated": primary_updates, "tardbp_motif_refs_inserted": tar_motifs,
        "tardbp_structure_rows_inserted": tar_models, "tardbp_links_inserted": tar_links,
        "tardbp_run_interface_counts": {str(k): v for k, v in sorted(tar_interactions.items())},
        "tf_count_rows_updated": count_updates, "protected_tables_unchanged": True,
        **final,
    }


def main() -> int:
    args = parse_args()
    delivered, decisions = validate_inputs(args)
    for output in (args.preview_output, args.qc_output):
        if output.exists() and not args.replace_outputs:
            raise RestorationError(f"output exists; use --replace-outputs: {output}")
    integrated_before = sha256_file(args.integrated_tsv)
    disk_before = sha256_file(args.database)
    preview_ids: set[str]
    with connect_ro(args.database) as current:
        current_ids = {str(row[0]) for row in current.execute("SELECT tf_id FROM tf")}
        preview_ids = {r["TF_name"] for r in delivered} - current_ids
    write_preview(args.preview_output, delivered, decisions, preview_ids, args.replace_outputs)

    backup: Path | None = None
    historical = connect_ro(args.historical_database)
    if args.dry_run:
        connection = in_memory_copy(args.database)
    else:
        backup, backed_hash = backup_database(args.database, args.backup_dir)
        if backed_hash != disk_before:
            raise RestorationError("database changed before apply")
        connection = sqlite3.connect(args.database)
        connection.row_factory = sqlite3.Row
        connection.execute("PRAGMA foreign_keys=ON")
        connection.execute("BEGIN IMMEDIATE")
    try:
        result = execute_restoration(
            connection,
            historical,
            delivered,
            decisions,
            args.model_cache,
            args.integrated_tsv,
            args.decision_manifest,
        )
        if args.apply:
            connection.commit()
        else:
            connection.rollback()
    except Exception:
        if args.apply:
            connection.rollback()
        raise
    finally:
        connection.close(); historical.close()

    if sha256_file(args.integrated_tsv) != integrated_before:
        raise RestorationError("immutable integrated TSV changed")
    disk_after = sha256_file(args.database)
    if args.dry_run and disk_after != disk_before:
        raise RestorationError("dry-run modified the staging database")
    result.update({
        "mode": "APPLY" if args.apply else "DRY-RUN",
        "database_sha256_before": disk_before, "database_sha256_after": disk_after,
        "backup": str(backup) if backup else None,
        "integrated_tsv_sha256": integrated_before,
        "integrated_tsv_unchanged": True,
        "all_invariants_passed": True,
    })
    args.qc_output.parent.mkdir(parents=True, exist_ok=True)
    args.qc_output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")

    print(f"Mode: {result['mode']}")
    print(f"Initial TF rows: {result['initial_tf_rows']}")
    print(f"TF rows inserted: {result['tf_rows_inserted']}")
    print(f"Final TF rows: {result['final_tf_rows']}")
    print(f"Corrected PWM motif refs inserted: {result['corrected_motif_refs_inserted']}")
    print(f"TARDBP motif/model/link inserts: {result['tardbp_motif_refs_inserted']}/{result['tardbp_structure_rows_inserted']}/{result['tardbp_links_inserted']}")
    print(f"TARDBP interface yes/no: {result['tardbp_interface_yes']}/{result['tardbp_interface_no']}")
    print(f"Protected tables unchanged: {'YES' if result['protected_tables_unchanged'] else 'NO'}")
    print(f"Database SHA256 before: {disk_before}")
    print(f"Database SHA256 after:  {disk_after}")
    print(f"Backup: {backup if backup else '<not created in dry-run>'}")
    print("All invariants passed: YES")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
