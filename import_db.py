#!/usr/bin/env python3
"""Build the local TF web database from the TSV and archive files.

The importer keeps the small MEME motif files in SQLite for fast display.
The large model archive is indexed by member path only, so the first version
does not duplicate nearly a gigabyte of structure files.
"""

from __future__ import annotations

import argparse
import csv
import os
import json
import re
import sqlite3
import tarfile
from pathlib import Path


APP_DIR = Path(__file__).resolve().parent
DEFAULT_DATA_DIR = Path(os.environ.get("TF_WEBDB_DATA_DIR", APP_DIR / "data" / "raw"))
DEFAULT_DB = APP_DIR / "data" / "tf_webdb.sqlite"
DEFAULT_SOURCE_RELEASES = APP_DIR / "config" / "source_releases.tsv"

BASES = ("A", "C", "G", "T")
SOURCE_LABELS = {
    "jaspar": "JASPAR",
    "cisbp": "CisBP",
    "hocomoco": "HOCOMOCO",
    "modcre": "ModCRE",
    "alphafold": "AlphaFold",
    "uniprot": "UniProt",
}

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


def table_exists(conn: sqlite3.Connection, table_name: str) -> bool:
    return (
        conn.execute(
            "SELECT 1 FROM sqlite_master WHERE type = 'table' AND name = ?",
            (table_name,),
        ).fetchone()
        is not None
    )


def ensure_source_release_table(conn: sqlite3.Connection) -> None:
    conn.execute(
        """
        CREATE TABLE IF NOT EXISTS source_release (
            source TEXT NOT NULL REFERENCES source(source),
            release_name TEXT NOT NULL,
            collection TEXT,
            species_scope TEXT,
            motif_or_model_type TEXT,
            source_url TEXT,
            download_url TEXT,
            citation TEXT,
            license_note TEXT,
            downloaded_at TEXT,
            local_file_label TEXT,
            checksum_sha256 TEXT,
            confirmation_status TEXT NOT NULL,
            notes TEXT,
            PRIMARY KEY (source, release_name, collection)
        )
        """
    )
    conn.execute("CREATE INDEX IF NOT EXISTS idx_source_release_source ON source_release(source)")


def load_source_releases(
    conn: sqlite3.Connection,
    config_path: Path = DEFAULT_SOURCE_RELEASES,
) -> int:
    if not config_path.exists():
        return 0
    ensure_source_release_table(conn)
    metadata = {}
    if table_exists(conn, "metadata"):
        metadata = dict(conn.execute("SELECT key, value FROM metadata").fetchall())

    loaded = 0
    with config_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            clean = {col: (row.get(col) or "").strip() for col in SOURCE_RELEASE_COLUMNS}
            source = clean["source"]
            if not source:
                continue
            clean["release_name"] = clean["release_name"] or "pending confirmation"
            clean["collection"] = clean["collection"] or "pending confirmation"
            clean["confirmation_status"] = clean["confirmation_status"] or "pending confirmation"
            if source == "uniprot" and metadata.get("uniprot_fetched_at"):
                clean["downloaded_at"] = metadata["uniprot_fetched_at"]
            conn.execute(
                "INSERT OR IGNORE INTO source (source, label, description) VALUES (?, ?, ?)",
                (source, SOURCE_LABELS.get(source, source), f"{SOURCE_LABELS.get(source, source)} source metadata"),
            )
            conn.execute(
                """
                INSERT OR REPLACE INTO source_release
                    (source, release_name, collection, species_scope, motif_or_model_type,
                     source_url, download_url, citation, license_note, downloaded_at,
                     local_file_label, checksum_sha256, confirmation_status, notes)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                tuple(clean[col] for col in SOURCE_RELEASE_COLUMNS),
            )
            loaded += 1
    return loaded


SCHEMA = """
PRAGMA foreign_keys = ON;

DROP TABLE IF EXISTS motif_structure;
DROP TABLE IF EXISTS model_summary;
DROP TABLE IF EXISTS motif_ref;
DROP TABLE IF EXISTS tf_family;
DROP TABLE IF EXISTS motif_file;
DROP TABLE IF EXISTS structure_file;
DROP TABLE IF EXISTS tf_annotation;
DROP TABLE IF EXISTS tf;
DROP TABLE IF EXISTS source_release;
DROP TABLE IF EXISTS source;
DROP TABLE IF EXISTS import_issue;
DROP TABLE IF EXISTS metadata;

CREATE TABLE metadata (
    key TEXT PRIMARY KEY,
    value TEXT NOT NULL
);

CREATE TABLE source (
    source TEXT PRIMARY KEY,
    label TEXT NOT NULL,
    description TEXT NOT NULL
);

CREATE TABLE source_release (
    source TEXT NOT NULL REFERENCES source(source),
    release_name TEXT NOT NULL,
    collection TEXT,
    species_scope TEXT,
    motif_or_model_type TEXT,
    source_url TEXT,
    download_url TEXT,
    citation TEXT,
    license_note TEXT,
    downloaded_at TEXT,
    local_file_label TEXT,
    checksum_sha256 TEXT,
    confirmation_status TEXT NOT NULL,
    notes TEXT,
    PRIMARY KEY (source, release_name, collection)
);

CREATE TABLE tf (
    tf_id TEXT PRIMARY KEY,
    family_text TEXT NOT NULL,
    motif_ref_count INTEGER NOT NULL DEFAULT 0,
    active_model_count INTEGER NOT NULL DEFAULT 0
);

CREATE TABLE tf_annotation (
    tf_id TEXT PRIMARY KEY REFERENCES tf(tf_id) ON DELETE CASCADE,
    uniprot_accession TEXT NOT NULL,
    entry_name TEXT,
    gene_names TEXT,
    protein_name TEXT,
    organism_name TEXT,
    organism_id INTEGER,
    reviewed INTEGER,
    sequence_length INTEGER,
    annotation_score REAL,
    uniprot_url TEXT,
    fetched_at TEXT NOT NULL
);

CREATE TABLE tf_family (
    tf_id TEXT NOT NULL REFERENCES tf(tf_id) ON DELETE CASCADE,
    family TEXT NOT NULL,
    PRIMARY KEY (tf_id, family)
);

CREATE TABLE motif_file (
    source TEXT NOT NULL REFERENCES source(source),
    motif_id TEXT NOT NULL,
    member_path TEXT NOT NULL,
    archive_path TEXT NOT NULL,
    content TEXT NOT NULL,
    width INTEGER,
    nsites TEXT,
    consensus TEXT,
    matrix_json TEXT,
    PRIMARY KEY (source, motif_id)
);

CREATE TABLE motif_ref (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    tf_id TEXT NOT NULL REFERENCES tf(tf_id) ON DELETE CASCADE,
    evidence_type TEXT NOT NULL,
    source TEXT NOT NULL REFERENCES source(source),
    motif_id TEXT NOT NULL,
    original_value TEXT NOT NULL,
    identity_percent REAL,
    missing_local_file INTEGER NOT NULL DEFAULT 0
);

CREATE TABLE structure_file (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    source TEXT NOT NULL REFERENCES source(source),
    model_id TEXT NOT NULL,
    tf_id TEXT,
    member_path TEXT NOT NULL,
    archive_path TEXT NOT NULL,
    file_type TEXT NOT NULL,
    status TEXT NOT NULL,
    template_pdb TEXT,
    residue_start INTEGER,
    residue_end INTEGER
);

CREATE TABLE model_summary (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    summary_file_id INTEGER NOT NULL REFERENCES structure_file(id) ON DELETE CASCADE,
    matched_structure_id INTEGER REFERENCES structure_file(id) ON DELETE SET NULL,
    source TEXT NOT NULL REFERENCES source(source),
    status TEXT NOT NULL,
    tf_id TEXT,
    summary_model_id TEXT NOT NULL,
    model_rank INTEGER,
    n TEXT,
    template_pdb TEXT,
    n_tails TEXT,
    c_tails TEXT,
    protein_chain TEXT,
    dna_chain TEXT,
    identities TEXT,
    coverage TEXT,
    template_by_rmsd TEXT,
    domain TEXT,
    identity_percent REAL,
    similarity_percent REAL
);

CREATE TABLE motif_structure (
    motif_ref_id INTEGER NOT NULL REFERENCES motif_ref(id) ON DELETE CASCADE,
    structure_file_id INTEGER NOT NULL REFERENCES structure_file(id) ON DELETE CASCADE,
    PRIMARY KEY (motif_ref_id, structure_file_id)
);

CREATE TABLE import_issue (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    severity TEXT NOT NULL,
    category TEXT NOT NULL,
    message TEXT NOT NULL,
    tf_id TEXT,
    source TEXT,
    motif_id TEXT
);

CREATE INDEX idx_source_release_source ON source_release(source);
CREATE INDEX idx_tf_family_family ON tf_family(family);
CREATE INDEX idx_tf_annotation_gene_names ON tf_annotation(gene_names);
CREATE INDEX idx_tf_annotation_organism_name ON tf_annotation(organism_name);
CREATE INDEX idx_motif_ref_tf ON motif_ref(tf_id);
CREATE INDEX idx_motif_ref_motif ON motif_ref(motif_id);
CREATE INDEX idx_motif_ref_evidence ON motif_ref(evidence_type);
CREATE INDEX idx_structure_file_tf ON structure_file(tf_id);
CREATE INDEX idx_structure_file_model ON structure_file(source, model_id);
CREATE INDEX idx_structure_file_status ON structure_file(status);
CREATE INDEX idx_model_summary_tf ON model_summary(tf_id);
CREATE INDEX idx_model_summary_template ON model_summary(template_pdb);
CREATE INDEX idx_model_summary_matched_structure ON model_summary(matched_structure_id);
"""


def strip_known_extension(filename: str) -> str:
    if filename.endswith(".summary.txt"):
        return filename[: -len(".txt")]
    for suffix in (".meme", ".pdb", ".pir"):
        if filename.endswith(suffix):
            return filename[: -len(suffix)]
    return filename


def normalize_motif_token(token: str) -> tuple[str, float | None]:
    token = token.strip()
    match = re.match(r"^(.*?)\s*\((\d+(?:\.\d+)?)%\)\s*$", token)
    if match:
        return strip_known_extension(match.group(1).strip()), float(match.group(2))
    return strip_known_extension(token), None


def source_for_motif_id(motif_id: str, evidence_type: str) -> str:
    if evidence_type == "modcre":
        return "modcre"
    if evidence_type == "alphafold":
        return "alphafold"
    if motif_id.startswith("MA"):
        return "jaspar"
    if motif_id.startswith("M"):
        return "cisbp"
    return "unknown"


def parse_meme(content: str) -> tuple[int | None, str | None, str | None, str | None]:
    width = None
    nsites = None
    matrix: list[list[float]] = []
    lines = content.splitlines()

    matrix_start = None
    for index, line in enumerate(lines):
        if "letter-probability matrix:" in line:
            matrix_start = index + 1
            width_match = re.search(r"\bw=\s*(\d+)", line)
            nsites_match = re.search(r"\bnsites=\s*([^\s]+)", line)
            if width_match:
                width = int(width_match.group(1))
            if nsites_match:
                nsites = nsites_match.group(1)
            break

    if matrix_start is not None:
        for line in lines[matrix_start:]:
            parts = line.split()
            if len(parts) < 4:
                if matrix:
                    break
                continue
            try:
                row = [float(value) for value in parts[:4]]
            except ValueError:
                if matrix:
                    break
                continue
            matrix.append(row)
            if width is not None and len(matrix) >= width:
                break

    if width is None and matrix:
        width = len(matrix)

    consensus = None
    matrix_json = None
    if matrix:
        consensus = "".join(BASES[max(range(4), key=lambda i: row[i])] for row in matrix)
        matrix_json = json.dumps(matrix, separators=(",", ":"))

    return width, nsites, consensus, matrix_json


def motif_source_from_member(member_path: str) -> str | None:
    if not member_path.endswith(".meme"):
        return None
    if member_path.startswith("jaspar/"):
        return "jaspar"
    if member_path.startswith("pbm/"):
        return "cisbp"
    if member_path.startswith("pwms_modcre_all/"):
        return "modcre"
    if member_path.startswith("pwms_af3_all/"):
        return "alphafold"
    return None


def index_motif_archive(
    conn: sqlite3.Connection,
    archive_path: Path,
    expected_source: str | None = None,
) -> set[tuple[str, str]]:
    present: set[tuple[str, str]] = set()
    with tarfile.open(archive_path, "r:gz") as archive:
        for member in archive:
            if not member.isfile():
                continue
            source = motif_source_from_member(member.name)
            if source is None:
                continue
            if expected_source and source != expected_source:
                continue
            extracted = archive.extractfile(member)
            if extracted is None:
                continue
            content = extracted.read().decode("utf-8", errors="replace")
            motif_id = strip_known_extension(Path(member.name).name)
            width, nsites, consensus, matrix_json = parse_meme(content)
            conn.execute(
                """
                INSERT OR REPLACE INTO motif_file
                    (source, motif_id, member_path, archive_path, content, width, nsites, consensus, matrix_json)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    source,
                    motif_id,
                    member.name,
                    str(archive_path),
                    content,
                    width,
                    nsites,
                    consensus,
                    matrix_json,
                ),
            )
            present.add((source, motif_id))
    return present


def infer_model_source(member_path: str) -> tuple[str | None, str | None]:
    top_dir = member_path.split("/", 1)[0]
    status = "failed" if "failed" in top_dir.lower() else "active"
    if top_dir.startswith("MODELS_AF3"):
        return "alphafold", status
    if top_dir.startswith("models"):
        return "modcre", status
    return None, None


def parse_model_metadata(model_id: str, source: str) -> tuple[str | None, str | None, int | None, int | None]:
    tf_id = None
    template_pdb = None
    residue_start = None
    residue_end = None

    if source == "alphafold":
        tf_id = model_id.split(".", 1)[0].replace("_dimer", "").replace("_monomer", "").upper()
    else:
        match = re.match(r"^(?:TFS|DIMER)_(?:sp|tr)_([A-Z0-9]+)_", model_id)
        if match:
            tf_id = match.group(1)
        template_match = re.search(r"_([0-9A-Za-z]{4})_[A-Za-z0-9]+_\d+$", model_id)
        if template_match:
            template_pdb = template_match.group(1).lower()
        range_match = re.search(r":(\d+):(\d+)_", model_id)
        if range_match:
            residue_start = int(range_match.group(1))
            residue_end = int(range_match.group(2))

    return tf_id, template_pdb, residue_start, residue_end


def parse_int(value: str | None) -> int | None:
    if value is None or value == "" or value == "NA":
        return None
    try:
        return int(value)
    except ValueError:
        return None


def parse_float(value: str | None) -> float | None:
    if value is None or value == "" or value == "NA":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def summary_base_id(model_id: str) -> str:
    return model_id.removesuffix("_model.summary")


def split_multi_value(value: str | None) -> list[str]:
    if not value or value == "NA":
        return []
    return [item.strip() for item in value.split(",") if item.strip()]


def candidate_model_ids(base_id: str, row: dict[str, str]) -> list[str]:
    template = (row.get("template") or "").strip().lower()
    model_number = (row.get("n") or "").strip()
    n_tails = split_multi_value(row.get("N-tails"))
    c_tails = split_multi_value(row.get("C-tails"))
    protein_chains = split_multi_value(row.get("protein-chain"))
    if not template or not model_number or not n_tails or not c_tails or not protein_chains:
        return []
    count = max(len(n_tails), len(c_tails), len(protein_chains))
    candidates = []
    for index in range(count):
        n_tail = n_tails[index] if index < len(n_tails) else n_tails[-1]
        c_tail = c_tails[index] if index < len(c_tails) else c_tails[-1]
        chain = protein_chains[index] if index < len(protein_chains) else protein_chains[-1]
        candidates.append(f"{base_id}:{n_tail}:{c_tail}_{template}_{chain}_{model_number}")
    return candidates


def insert_model_summary_rows(
    conn: sqlite3.Connection,
    summary_file_id: int,
    source: str,
    status: str,
    tf_id: str | None,
    summary_model_id: str,
    content: str,
) -> int:
    lines = [line for line in content.splitlines() if line.strip()]
    if len(lines) < 2:
        return 0
    header = lines[0].lstrip("#").split(";")
    inserted = 0
    for line in lines[1:]:
        values = line.split(";")
        values += [""] * (len(header) - len(values))
        row = {header[index]: values[index].strip() for index in range(len(header))}
        matched_id = None
        for candidate in candidate_model_ids(summary_model_id, row):
            match = conn.execute(
                """
                SELECT id FROM structure_file
                WHERE source = ? AND status = ? AND model_id = ? AND file_type IN ('pdb', 'pir')
                ORDER BY CASE file_type WHEN 'pdb' THEN 0 ELSE 1 END
                LIMIT 1
                """,
                (source, status, candidate),
            ).fetchone()
            if match:
                matched_id = match[0]
                break
        conn.execute(
            """
            INSERT INTO model_summary
                (summary_file_id, matched_structure_id, source, status, tf_id,
                 summary_model_id, model_rank, n, template_pdb, n_tails,
                 c_tails, protein_chain, dna_chain, identities, coverage,
                 template_by_rmsd, domain, identity_percent, similarity_percent)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                summary_file_id,
                matched_id,
                source,
                status,
                tf_id,
                summary_model_id,
                parse_int(row.get("model")),
                row.get("n") or None,
                (row.get("template") or "").lower() or None,
                row.get("N-tails") or None,
                row.get("C-tails") or None,
                row.get("protein-chain") or None,
                row.get("DNA-chain") or None,
                row.get("#identities") or None,
                row.get("#coverage") or None,
                row.get("template_by_RMSD") or None,
                row.get("domain") or None,
                parse_float(row.get("%_identity")),
                parse_float(row.get("%_similarity")),
            ),
        )
        inserted += 1
    return inserted


def update_model_summary_links(conn: sqlite3.Connection) -> int:
    rows = conn.execute(
        """
        SELECT id, source, status, summary_model_id, n, template_pdb, n_tails, c_tails, protein_chain
        FROM model_summary
        WHERE matched_structure_id IS NULL
        """
    ).fetchall()
    linked = 0
    for row in rows:
        summary_row = {
            "n": row[4] or "",
            "template": row[5] or "",
            "N-tails": row[6] or "",
            "C-tails": row[7] or "",
            "protein-chain": row[8] or "",
        }
        for candidate in candidate_model_ids(row[3], summary_row):
            match = conn.execute(
                """
                SELECT id FROM structure_file
                WHERE source = ? AND status = ? AND model_id = ? AND file_type IN ('pdb', 'pir')
                ORDER BY CASE file_type WHEN 'pdb' THEN 0 ELSE 1 END
                LIMIT 1
                """,
                (row[1], row[2], candidate),
            ).fetchone()
            if match:
                conn.execute(
                    "UPDATE model_summary SET matched_structure_id = ? WHERE id = ?",
                    (match[0], row[0]),
                )
                linked += 1
                break
    return linked


def index_model_archive(conn: sqlite3.Connection, archive_path: Path) -> int:
    count = 0
    with tarfile.open(archive_path, "r:gz") as archive:
        for member in archive:
            if not member.isfile():
                continue
            name = Path(member.name).name
            if name.endswith(".summary.txt"):
                file_type = "summary"
            elif name.endswith(".pdb"):
                file_type = "pdb"
            elif name.endswith(".pir"):
                file_type = "pir"
            else:
                continue
            source, status = infer_model_source(member.name)
            if source is None or status is None:
                continue
            model_id = strip_known_extension(name)
            tf_id, template_pdb, residue_start, residue_end = parse_model_metadata(model_id, source)
            cursor = conn.execute(
                """
                INSERT INTO structure_file
                    (source, model_id, tf_id, member_path, archive_path, file_type, status,
                     template_pdb, residue_start, residue_end)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    source,
                    model_id,
                    tf_id,
                    member.name,
                    str(archive_path),
                    file_type,
                    status,
                    template_pdb,
                    residue_start,
                    residue_end,
                ),
            )
            if file_type == "summary":
                extracted = archive.extractfile(member)
                if extracted is not None:
                    content = extracted.read().decode("utf-8", errors="replace")
                    insert_model_summary_rows(
                        conn,
                        cursor.lastrowid,
                        source,
                        status,
                        tf_id,
                        summary_base_id(model_id),
                        content,
                    )
            count += 1
    return count


def read_chart(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        headers = [header.strip() for header in next(reader)]
        rows = []
        for values in reader:
            if not values:
                continue
            values = values + [""] * (len(headers) - len(values))
            rows.append({headers[i]: values[i].strip() for i in range(len(headers))})
    return rows


def split_list_field(value: str) -> list[str]:
    if not value.strip():
        return []
    return [item.strip() for item in value.split(";") if item.strip()]


def insert_chart_rows(
    conn: sqlite3.Connection,
    rows: list[dict[str, str]],
    present_motifs: set[tuple[str, str]],
) -> None:
    evidence_columns = [
        ("Identical_PWM", "identical"),
        ("Homologous_PWM", "homologous"),
        ("Relatively_Homologous_PWM", "relative_homologous"),
        ("ModCRE", "modcre"),
        ("AlphaFold", "alphafold"),
    ]

    for row in rows:
        tf_id = row["TF_name"]
        family_text = row.get("TF_family", "")
        conn.execute(
            "INSERT OR REPLACE INTO tf (tf_id, family_text) VALUES (?, ?)",
            (tf_id, family_text),
        )
        for family in split_list_field(family_text.replace(",", ";")):
            conn.execute(
                "INSERT OR IGNORE INTO tf_family (tf_id, family) VALUES (?, ?)",
                (tf_id, family),
            )

        for column, evidence_type in evidence_columns:
            for raw_token in split_list_field(row.get(column, "")):
                motif_id, identity_percent = normalize_motif_token(raw_token)
                source = source_for_motif_id(motif_id, evidence_type)
                missing = 0 if (source, motif_id) in present_motifs else 1
                conn.execute(
                    """
                    INSERT INTO motif_ref
                        (tf_id, evidence_type, source, motif_id, original_value,
                         identity_percent, missing_local_file)
                    VALUES (?, ?, ?, ?, ?, ?, ?)
                    """,
                    (
                        tf_id,
                        evidence_type,
                        source,
                        motif_id,
                        raw_token,
                        identity_percent,
                        missing,
                    ),
                )
                if missing:
                    conn.execute(
                        """
                        INSERT INTO import_issue
                            (severity, category, message, tf_id, source, motif_id)
                        VALUES (?, ?, ?, ?, ?, ?)
                        """,
                        (
                            "warning",
                            "missing_motif_file",
                            "Chart references a motif that is not present in the local motif archives.",
                            tf_id,
                            source,
                            motif_id,
                        ),
                    )


def finalize_links_and_counts(conn: sqlite3.Connection) -> None:
    conn.execute(
        """
        INSERT OR IGNORE INTO motif_structure (motif_ref_id, structure_file_id)
        SELECT mr.id, sf.id
        FROM motif_ref AS mr
        JOIN structure_file AS sf
          ON sf.source = mr.source
         AND sf.model_id = mr.motif_id
         AND sf.file_type = 'pdb'
        WHERE mr.source IN ('modcre', 'alphafold')
        """
    )
    conn.execute(
        """
        UPDATE tf
           SET motif_ref_count = (
                   SELECT COUNT(*) FROM motif_ref WHERE motif_ref.tf_id = tf.tf_id
               ),
               active_model_count = (
                   SELECT COUNT(*)
                   FROM structure_file
                   WHERE structure_file.tf_id = tf.tf_id
                     AND structure_file.status = 'active'
                     AND structure_file.file_type = 'pdb'
               )
        """
    )


def build_database(data_dir: Path, db_path: Path, skip_model_index: bool = False) -> None:
    data_dir = data_dir.resolve()
    db_path.parent.mkdir(parents=True, exist_ok=True)
    if db_path.exists():
        db_path.unlink()

    chart_path = data_dir / "TF_PWM_chart_final.tsv"
    archives = {
        "jaspar": data_dir / "jaspar.tar.gz",
        "cisbp": data_dir / "cisbp.tar.gz",
        "pwms": data_dir / "pwms.tar.gz",
        "models": data_dir / "models.tar.gz",
    }

    for path in [chart_path, *archives.values()]:
        if not path.exists():
            raise FileNotFoundError(path)

    conn = sqlite3.connect(db_path)
    try:
        conn.executescript(SCHEMA)
        for source, label in SOURCE_LABELS.items():
            conn.execute(
                "INSERT INTO source (source, label, description) VALUES (?, ?, ?)",
                (source, label, f"{label} motif/model evidence"),
            )
        load_source_releases(conn)

        present_motifs: set[tuple[str, str]] = set()
        present_motifs |= index_motif_archive(conn, archives["jaspar"], "jaspar")
        present_motifs |= index_motif_archive(conn, archives["cisbp"], "cisbp")
        present_motifs |= index_motif_archive(conn, archives["pwms"])

        rows = read_chart(chart_path)
        insert_chart_rows(conn, rows, present_motifs)

        model_count = 0
        if not skip_model_index:
            model_count = index_model_archive(conn, archives["models"])
            update_model_summary_links(conn)

        finalize_links_and_counts(conn)

        stats = {
            "chart_rows": str(len(rows)),
            "motif_files": str(len(present_motifs)),
            "model_files_indexed": str(model_count),
            "model_summary_rows": str(conn.execute("SELECT COUNT(*) FROM model_summary").fetchone()[0]),
            "source_data_dir": str(data_dir),
        }
        for key, value in stats.items():
            conn.execute("INSERT INTO metadata (key, value) VALUES (?, ?)", (key, value))
        conn.commit()
    finally:
        conn.close()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build the local TF web database.")
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--db", type=Path, default=DEFAULT_DB)
    parser.add_argument(
        "--skip-model-index",
        action="store_true",
        help="Import motifs and TF rows only. Useful for quick UI testing.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    build_database(args.data_dir, args.db, args.skip_model_index)
    print(f"Built {args.db}")


if __name__ == "__main__":
    main()
