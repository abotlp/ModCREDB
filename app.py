#!/usr/bin/env python3
"""Tiny local web interface for the TF database.

This intentionally uses the Python standard library plus Jinja2, because the
Masada environment already has those pieces available.
"""

from __future__ import annotations

import argparse
import csv
import html
import io
import json
import math
import mimetypes
import os
import re
import shutil
import sqlite3
import subprocess
import tarfile
import tempfile
import traceback
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, quote, unquote, urlparse

from jinja2 import Environment, FileSystemLoader, select_autoescape


APP_DIR = Path(__file__).resolve().parent
DEFAULT_DB = APP_DIR / "data" / "tf_webdb.sqlite"
DEFAULT_MAX_POST_BYTES = int(os.environ.get("TF_WEBDB_MAX_POST_BYTES", "1000000"))
EVIDENCE_LABELS = {
    "identical": "Identical PWM",
    "homologous": "Homologous PWM",
    "relative_homologous": "Relatively Homologous PWM",
    "modcre": "ModCRE",
    "alphafold": "AlphaFold",
}
EVIDENCE_INFO = {
    "identical": {
        "trust": "Highest",
        "description": "Known PWM assigned directly to the same TF accession or an apparently identical TF entry.",
        "use": "Use as the primary motif when available.",
    },
    "homologous": {
        "trust": "High",
        "description": "PWM transferred from a close homolog. In this chart the reported identities are 70.0-99.9%.",
        "use": "Good candidate motif, but keep the source TF and identity percentage visible.",
    },
    "relative_homologous": {
        "trust": "Medium",
        "description": "PWM transferred from a more distant homolog. In this chart the reported identities are 50.0-69.9%.",
        "use": "Useful when no closer motif exists; should be interpreted cautiously.",
    },
    "modcre": {
        "trust": "Predicted",
        "description": "Structure-based PWM predicted from ModCRE TF-DNA models and template complexes.",
        "use": "Useful as a predicted specificity model; inspect template, residue range, identity, and coverage.",
    },
    "alphafold": {
        "trust": "Predicted",
        "description": "PWM predicted from AlphaFold-derived TF-DNA models.",
        "use": "Useful exploratory evidence; needs structural and biological validation.",
    },
}
TF_SCAN_DEFAULT_EVIDENCE = ["identical", "homologous", "modcre", "alphafold"]
TF_SCAN_MAX_MOTIFS = 500
SOURCE_LABELS = {
    "jaspar": "JASPAR",
    "cisbp": "CisBP",
    "hocomoco": "HOCOMOCO",
    "modcre": "ModCRE",
    "alphafold": "AlphaFold",
}
SOURCE_HOME_URLS = {
    "jaspar": "https://jaspar2024.elixir.no/",
    "cisbp": "https://cisbp.ccbr.utoronto.ca/",
    "hocomoco": "https://hocomoco11.autosome.org/",
}
BASE_COLORS = {
    "A": "#23845b",
    "C": "#2f65d9",
    "G": "#d29612",
    "T": "#c73f43",
}
BASE_INDEX = {"A": 0, "C": 1, "G": 2, "T": 3}
REVCOMP = str.maketrans("ACGT", "TGCA")


def connect(db_path: Path) -> sqlite3.Connection:
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    return conn


def dict_row(row: sqlite3.Row | None) -> dict | None:
    return dict(row) if row is not None else None


def rows_with_percent(rows: list[sqlite3.Row], value_key: str = "count") -> list[dict[str, object]]:
    dicts = [dict(row) for row in rows]
    max_value = max((int(row.get(value_key) or 0) for row in dicts), default=0)
    for row in dicts:
        value = int(row.get(value_key) or 0)
        row["percent"] = 0 if max_value == 0 else round((value / max_value) * 100, 1)
    return dicts


def safe_count(conn: sqlite3.Connection, sql: str) -> int:
    try:
        return conn.execute(sql).fetchone()[0]
    except sqlite3.OperationalError:
        return 0


def motif_external_links(source: str, motif_id: str) -> list[dict[str, str]]:
    if source == "jaspar" and motif_id.startswith("MA"):
        return [
            {
                "label": "JASPAR matrix",
                "url": f"https://jaspar2024.elixir.no/matrix/{quote(motif_id, safe='')}/",
                "note": "Direct matrix page",
            }
        ]
    if source == "cisbp" and motif_id.startswith("M"):
        return [
            {
                "label": "CIS-BP database",
                "url": SOURCE_HOME_URLS["cisbp"],
                "note": "Search this motif ID in CIS-BP",
            }
        ]
    if source == "hocomoco":
        return [
            {
                "label": "HOCOMOCO v11",
                "url": SOURCE_HOME_URLS["hocomoco"],
                "note": "Search this motif ID in HOCOMOCO v11",
            }
        ]
    return []


def pdb_url(pdb_id: str | None) -> str:
    if not pdb_id:
        return ""
    return f"https://www.rcsb.org/structure/{quote(str(pdb_id).upper(), safe='')}"


def clean_dna_sequence(text: str) -> tuple[str, int]:
    sequence_parts = []
    ignored = 0
    for line in text.splitlines():
        if line.startswith(">"):
            continue
        for char in line.upper():
            if char in "ACGTN":
                sequence_parts.append(char)
            elif char.strip():
                ignored += 1
    return "".join(sequence_parts), ignored


def reverse_complement(sequence: str) -> str:
    return sequence.translate(REVCOMP)[::-1]


def parse_motif_specs(text: str) -> list[tuple[str, str]]:
    specs = []
    for token in re.split(r"[\n,]+", text):
        token = token.strip()
        if not token:
            continue
        if "|" in token:
            source, motif_id = token.split("|", 1)
        elif "\t" in token:
            source, motif_id = token.split("\t", 1)
        elif " " in token:
            source, motif_id = token.split(None, 1)
        elif ":" in token and token.split(":", 1)[0] in SOURCE_LABELS:
            source, motif_id = token.split(":", 1)
        else:
            source, motif_id = "", token
        specs.append((source.strip().lower(), motif_id.strip()))
    return specs


def motif_spec(source: str, motif_id: str) -> str:
    return f"{source}|{motif_id}"


def search_scan_motifs(
    conn: sqlite3.Connection,
    query: str = "",
    source: str = "",
    limit: int = 15,
) -> list[sqlite3.Row]:
    where = ["mf.matrix_json IS NOT NULL"]
    args: list[object] = []
    if source:
        where.append("mf.source = ?")
        args.append(source)
    if query:
        like = f"%{query}%"
        where.append(
            """
            (mf.motif_id LIKE ? OR mf.consensus LIKE ? OR mr.tf_id LIKE ?
             OR ta.gene_names LIKE ? OR ta.protein_name LIKE ? OR ta.organism_name LIKE ?)
            """
        )
        args.extend([like, like, like, like, like, like])
    where_sql = " AND ".join(where)
    return conn.execute(
        f"""
        SELECT mf.source, mf.motif_id, mf.width, mf.consensus,
               COUNT(DISTINCT mr.tf_id) AS tf_count,
               GROUP_CONCAT(
                   DISTINCT COALESCE(NULLIF(ta.gene_names, ''), mr.tf_id)
               ) AS example_tfs
        FROM motif_file AS mf
        LEFT JOIN motif_ref AS mr
          ON mr.source = mf.source
         AND mr.motif_id = mf.motif_id
        LEFT JOIN tf_annotation AS ta ON ta.tf_id = mr.tf_id
        WHERE {where_sql}
        GROUP BY mf.source, mf.motif_id
        ORDER BY
          CASE mf.source
            WHEN 'jaspar' THEN 0
            WHEN 'cisbp' THEN 1
            WHEN 'modcre' THEN 2
            WHEN 'alphafold' THEN 3
            ELSE 4
          END,
          tf_count DESC,
          mf.motif_id
        LIMIT ?
        """,
        [*args, limit],
    ).fetchall()


def normalize_tf_id(tf_id: str) -> str:
    return tf_id.strip().upper()


def selected_evidence_from_params(params: dict[str, list[str]], tf_id: str) -> list[str]:
    raw_values = params.get("tf_evidence", [])
    selected = [value for value in raw_values if value in EVIDENCE_LABELS]
    if not selected and tf_id:
        return list(TF_SCAN_DEFAULT_EVIDENCE)
    return selected


def parse_region_value(value: str) -> tuple[int, int] | None:
    match = re.fullmatch(r"\s*(\d+)\s*-\s*(\d+)\s*", value or "")
    if not match:
        return None
    start = int(match.group(1))
    end = int(match.group(2))
    if start <= 0 or end < start:
        return None
    return start, end


def load_tf_scan_regions(conn: sqlite3.Connection, tf_id: str) -> list[dict[str, object]]:
    rows = conn.execute(
        """
        SELECT DISTINCT mr.source, mr.motif_id, sf.residue_start, sf.residue_end
        FROM motif_ref AS mr
        JOIN motif_structure AS ms ON ms.motif_ref_id = mr.id
        JOIN structure_file AS sf ON sf.id = ms.structure_file_id
        JOIN motif_file AS mf
          ON mf.source = mr.source
         AND mf.motif_id = mr.motif_id
        WHERE UPPER(mr.tf_id) = UPPER(?)
          AND sf.status = 'active'
          AND sf.residue_start IS NOT NULL
          AND sf.residue_end IS NOT NULL
          AND mf.matrix_json IS NOT NULL
        ORDER BY sf.residue_start, sf.residue_end
        """,
        (tf_id,),
    ).fetchall()
    intervals = sorted(
        {
            (int(row["residue_start"]), int(row["residue_end"]))
            for row in rows
            if row["residue_start"] is not None and row["residue_end"] is not None
        }
    )
    merged: list[list[int]] = []
    for start, end in intervals:
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)

    regions: list[dict[str, object]] = []
    for start, end in merged:
        motif_keys = {
            (row["source"], row["motif_id"])
            for row in rows
            if overlaps(start, end, int(row["residue_start"]), int(row["residue_end"]))
        }
        value = f"{start}-{end}"
        regions.append(
            {
                "start": start,
                "end": end,
                "value": value,
                "label": f"Region {value}",
                "motif_count": len(motif_keys),
            }
        )
    return regions


def load_tf_scan_motifs(
    conn: sqlite3.Connection,
    tf_id: str,
    evidence_types: list[str],
    region: tuple[int, int] | None = None,
    limit: int = TF_SCAN_MAX_MOTIFS,
) -> tuple[dict | None, list[sqlite3.Row], dict[str, object]]:
    tf = dict_row(
        conn.execute(
            """
            SELECT tf.tf_id, tf.family_text, ta.gene_names, ta.protein_name,
                   ta.organism_name, ta.sequence_length
            FROM tf
            LEFT JOIN tf_annotation AS ta ON ta.tf_id = tf.tf_id
            WHERE UPPER(tf.tf_id) = UPPER(?)
            """,
            (tf_id,),
        ).fetchone()
    )
    if not tf:
        return None, [], {}

    selected = [evidence for evidence in evidence_types if evidence in EVIDENCE_LABELS]
    if not selected:
        selected = list(TF_SCAN_DEFAULT_EVIDENCE)
    placeholders = ",".join("?" for _ in selected)
    args: list[object] = [tf["tf_id"], *selected]
    region_clause = ""
    if region:
        region_clause = """
              AND EXISTS (
                  SELECT 1
                  FROM motif_structure AS ms_region
                  JOIN structure_file AS sf_region
                    ON sf_region.id = ms_region.structure_file_id
                  WHERE ms_region.motif_ref_id = mr.id
                    AND sf_region.status = 'active'
                    AND sf_region.residue_start IS NOT NULL
                    AND sf_region.residue_end IS NOT NULL
                    AND sf_region.residue_start <= ?
                    AND sf_region.residue_end >= ?
              )
        """
        args.extend([region[1], region[0]])

    total_ready = conn.execute(
        f"""
        SELECT COUNT(*) FROM (
            SELECT mr.source, mr.motif_id
            FROM motif_ref AS mr
            JOIN motif_file AS mf
              ON mf.source = mr.source
             AND mf.motif_id = mr.motif_id
            WHERE mr.tf_id = ?
              AND mr.evidence_type IN ({placeholders})
              AND mf.matrix_json IS NOT NULL
              {region_clause}
            GROUP BY mr.source, mr.motif_id
        )
        """,
        args,
    ).fetchone()[0]
    skipped_no_matrix = conn.execute(
        f"""
        SELECT COUNT(*) FROM (
            SELECT mr.source, mr.motif_id
            FROM motif_ref AS mr
            LEFT JOIN motif_file AS mf
              ON mf.source = mr.source
             AND mf.motif_id = mr.motif_id
            WHERE mr.tf_id = ?
              AND mr.evidence_type IN ({placeholders})
              AND (mf.matrix_json IS NULL)
              {region_clause}
            GROUP BY mr.source, mr.motif_id
        )
        """,
        args,
    ).fetchone()[0]
    evidence_counts = conn.execute(
        f"""
        SELECT mr.evidence_type, COUNT(DISTINCT mr.source || '|' || mr.motif_id) AS count
        FROM motif_ref AS mr
        JOIN motif_file AS mf
          ON mf.source = mr.source
         AND mf.motif_id = mr.motif_id
        WHERE mr.tf_id = ?
          AND mr.evidence_type IN ({placeholders})
          AND mf.matrix_json IS NOT NULL
          {region_clause}
        GROUP BY mr.evidence_type
        ORDER BY
          CASE mr.evidence_type
            WHEN 'identical' THEN 0
            WHEN 'homologous' THEN 1
            WHEN 'relative_homologous' THEN 2
            WHEN 'modcre' THEN 3
            WHEN 'alphafold' THEN 4
            ELSE 5
          END
        """,
        args,
    ).fetchall()
    rows = conn.execute(
        f"""
        SELECT mf.source, mf.motif_id, mf.width, mf.nsites, mf.consensus, mf.matrix_json,
               GROUP_CONCAT(DISTINCT mr.evidence_type) AS evidence_types,
               COUNT(*) AS link_count,
               MIN(
                 CASE mr.evidence_type
                   WHEN 'identical' THEN 0
                   WHEN 'homologous' THEN 1
                   WHEN 'relative_homologous' THEN 2
                   WHEN 'modcre' THEN 3
                   WHEN 'alphafold' THEN 4
                   ELSE 5
                 END
               ) AS evidence_rank
        FROM motif_ref AS mr
        JOIN motif_file AS mf
          ON mf.source = mr.source
         AND mf.motif_id = mr.motif_id
        WHERE mr.tf_id = ?
          AND mr.evidence_type IN ({placeholders})
          AND mf.matrix_json IS NOT NULL
          {region_clause}
        GROUP BY mf.source, mf.motif_id
        ORDER BY evidence_rank, mf.source, mf.motif_id
        LIMIT ?
        """,
        [*args, limit],
    ).fetchall()
    summary = {
        "tf": tf,
        "selected_evidence": selected,
        "used_count": len(rows),
        "total_ready": total_ready,
        "skipped_no_matrix": skipped_no_matrix,
        "limited": total_ready > limit,
        "limit": limit,
        "evidence_counts": [dict(row) for row in evidence_counts],
        "selected_region": f"{region[0]}-{region[1]}" if region else "",
    }
    return tf, rows, summary


def parse_int_list(value: object) -> list[int]:
    if value is None:
        return []
    return [int(match) for match in re.findall(r"\d+", str(value))]


def summary_interval(row: sqlite3.Row) -> tuple[int, int] | None:
    starts = parse_int_list(row["n_tails"])
    ends = parse_int_list(row["c_tails"])
    if not starts or not ends:
        return None
    start = min(starts)
    end = max(ends)
    if start <= 0 or end < start:
        return None
    return start, end


def overlaps(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    return a_start <= b_end and b_start <= a_end


def row_interval(row: sqlite3.Row, start_key: str, end_key: str) -> tuple[int, int] | None:
    start = row[start_key]
    end = row[end_key]
    if start is None or end is None:
        return None
    start_int = int(start)
    end_int = int(end)
    if start_int <= 0 or end_int < start_int:
        return None
    return start_int, end_int


def build_region_groups(
    tf: dict,
    motif_rows: list[sqlite3.Row],
    active_models: list[sqlite3.Row],
    model_summaries: list[sqlite3.Row],
) -> list[dict[str, object]]:
    intervals: list[tuple[int, int]] = []
    for row in motif_rows:
        interval = row_interval(row, "region_start", "region_end")
        if interval:
            intervals.append(interval)
    for row in active_models:
        interval = row_interval(row, "residue_start", "residue_end")
        if interval:
            intervals.append(interval)
    for row in model_summaries:
        interval = summary_interval(row)
        if interval:
            intervals.append(interval)

    if not intervals:
        return []

    intervals.sort()
    merged: list[list[int]] = []
    for start, end in intervals:
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)

    sequence_length = int(tf.get("sequence_length") or 0)
    groups: list[dict[str, object]] = []
    for start, end in merged:
        motifs = [
            row
            for row in motif_rows
            if (interval := row_interval(row, "region_start", "region_end"))
            and overlaps(start, end, interval[0], interval[1])
        ]
        models = [
            row
            for row in active_models
            if (interval := row_interval(row, "residue_start", "residue_end"))
            and overlaps(start, end, interval[0], interval[1])
        ]
        summaries = [
            row
            for row in model_summaries
            if (interval := summary_interval(row)) and overlaps(start, end, interval[0], interval[1])
        ]
        templates = sorted(
            {
                str(value).upper()
                for value in [
                    *(row["template_pdb"] for row in models),
                    *(row["template_pdb"] for row in summaries),
                ]
                if value
            }
        )
        sources = sorted({SOURCE_LABELS.get(row["source"], row["source"]) for row in motifs})
        track_left = 0.0
        track_width = 0.0
        if sequence_length > 0:
            track_left = max(0.0, min(100.0, ((start - 1) / sequence_length) * 100.0))
            track_width = max(0.5, min(100.0 - track_left, ((end - start + 1) / sequence_length) * 100.0))
        groups.append(
            {
                "start": start,
                "end": end,
                "width": end - start + 1,
                "motifs": motifs,
                "models": models,
                "summaries": summaries,
                "templates": templates[:8],
                "template_count": len(templates),
                "sources": sources,
                "track_left": round(track_left, 3),
                "track_width": round(track_width, 3),
            }
        )
    return groups


def score_pwm_window(matrix: list[list[float]], window: str) -> tuple[float, float] | None:
    score = 0.0
    min_score = 0.0
    max_score = 0.0
    for row, base in zip(matrix, window):
        if base not in BASE_INDEX:
            return None
        log_odds = [math.log2(max(probability, 1e-6) / 0.25) for probability in row]
        score += log_odds[BASE_INDEX[base]]
        min_score += min(log_odds)
        max_score += max(log_odds)
    relative_score = 1.0 if max_score == min_score else (score - min_score) / (max_score - min_score)
    return score, relative_score


def scan_pwm(
    motif: sqlite3.Row,
    sequence: str,
    threshold: float,
    max_hits: int,
) -> list[dict[str, object]]:
    matrix = json.loads(motif["matrix_json"])
    width = len(matrix)
    hits = []
    if width == 0 or len(sequence) < width:
        return hits
    for start in range(0, len(sequence) - width + 1):
        forward = sequence[start : start + width]
        for strand, window in (("+", forward), ("-", reverse_complement(forward))):
            scored = score_pwm_window(matrix, window)
            if scored is None:
                continue
            score, relative_score = scored
            if relative_score < threshold:
                continue
            hits.append(
                {
                    "source": motif["source"],
                    "motif_id": motif["motif_id"],
                    "width": width,
                    "start": start + 1,
                    "end": start + width,
                    "strand": strand,
                    "sequence": forward if strand == "+" else window,
                    "score": score,
                    "relative_score": relative_score,
                    "pvalue": None,
                    "qvalue": None,
                    "engine": "local",
                }
            )
    hits.sort(key=lambda hit: (-float(hit["relative_score"]), -float(hit["score"]), int(hit["start"])))
    return hits[:max_hits]


def fimo_safe_id(source: str, motif_id: str) -> str:
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", motif_id).strip("_")
    return f"{source}__{safe}"


def build_combined_meme(motifs: list[sqlite3.Row]) -> tuple[str, dict[str, tuple[str, str]]]:
    lines = [
        "MEME version 4",
        "",
        "ALPHABET= ACGT",
        "",
        "strands: + -",
        "",
        "Background letter frequencies:",
        "A 0.25 C 0.25 G 0.25 T 0.25",
        "",
    ]
    id_map: dict[str, tuple[str, str]] = {}
    for motif in motifs:
        matrix = json.loads(motif["matrix_json"])
        safe_id = fimo_safe_id(motif["source"], motif["motif_id"])
        id_map[safe_id] = (motif["source"], motif["motif_id"])
        nsites = motif["nsites"] or "20"
        lines.append(f"MOTIF {safe_id} {motif['source']}|{motif['motif_id']}")
        lines.append(f"letter-probability matrix: alength= 4 w= {len(matrix)} nsites= {nsites} E= 0")
        for row in matrix:
            lines.append(" ".join(f"{float(value):.6f}" for value in row))
        lines.append("")
    return "\n".join(lines), id_map


def parse_float_for_sort(value: str | None, default: float = 0.0) -> float:
    if value is None or value == "":
        return default
    try:
        return float(value)
    except ValueError:
        return default


def csv_data_uri(fieldnames: list[str], rows: list[dict[str, object]]) -> str:
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=fieldnames, extrasaction="ignore")
    writer.writeheader()
    for row in rows:
        writer.writerow({field: row.get(field, "") for field in fieldnames})
    return "data:text/csv;charset=utf-8," + quote(buffer.getvalue(), safe="")


def build_hits_csv_uri(hits: list[dict[str, object]]) -> str:
    fields = [
        "source",
        "motif_id",
        "start",
        "end",
        "strand",
        "matched_sequence",
        "score",
        "pvalue",
        "qvalue",
    ]
    rows = [
        {
            "source": hit.get("source", ""),
            "motif_id": hit.get("motif_id", ""),
            "start": hit.get("start", ""),
            "end": hit.get("end", ""),
            "strand": hit.get("strand", ""),
            "matched_sequence": hit.get("sequence", ""),
            "score": hit.get("score", ""),
            "pvalue": hit.get("pvalue", ""),
            "qvalue": hit.get("qvalue", ""),
        }
        for hit in hits
    ]
    return csv_data_uri(fields, rows)


def build_profile_csv_uri(profile_summary: dict[str, object] | None) -> str:
    if not profile_summary:
        return ""
    rows = profile_summary.get("profile_rows")
    if not isinstance(rows, list):
        return ""
    fields = [
        "bin",
        "start",
        "end",
        "hit_count",
        "best_neg_log10_pvalue",
        "best_source",
        "best_motif",
        "best_pvalue",
    ]
    return csv_data_uri(fields, rows)


def build_fimo_profile_svg(
    hits: list[dict[str, object]],
    sequence_length: int,
    pvalue_threshold: float,
    max_bins: int = 900,
) -> tuple[str, dict[str, object] | None]:
    if not hits or sequence_length <= 0:
        return "", None

    bins = max(1, min(sequence_length, max_bins))
    best_signal = [0.0] * bins
    best_motif = [""] * bins
    best_source = [""] * bins
    best_pvalue = [""] * bins
    hit_counts = [0] * bins
    source_counts: dict[str, int] = {}
    profiled_hits = 0
    for hit in hits:
        start = hit.get("start")
        end = hit.get("end")
        pvalue = hit.get("pvalue")
        if not isinstance(start, int) or not isinstance(end, int):
            continue
        if not isinstance(pvalue, float) or pvalue <= 0:
            continue
        signal = -math.log10(max(pvalue, 1e-300))
        first_bin = max(0, min(bins - 1, int(((start - 1) / sequence_length) * bins)))
        last_bin = max(0, min(bins - 1, int(((end - 1) / sequence_length) * bins)))
        source = str(hit.get("source") or "")
        motif_id = str(hit.get("motif_id") or "")
        source_counts[source] = source_counts.get(source, 0) + 1
        for idx in range(first_bin, last_bin + 1):
            if signal > best_signal[idx]:
                best_signal[idx] = signal
                best_source[idx] = source
                best_motif[idx] = motif_id
                best_pvalue[idx] = f"{pvalue:.6g}"
            hit_counts[idx] += 1
        profiled_hits += 1

    max_signal = max(best_signal, default=0.0)
    max_count = max(hit_counts, default=0)
    if max_signal <= 0 or profiled_hits == 0:
        return "", None

    profile_rows: list[dict[str, object]] = []
    peak_rows: list[dict[str, object]] = []
    covered_bins = 0
    for idx in range(bins):
        start_pos = int((idx / bins) * sequence_length) + 1
        end_pos = max(start_pos, int(((idx + 1) / bins) * sequence_length))
        row = {
            "bin": idx + 1,
            "start": start_pos,
            "end": end_pos,
            "hit_count": hit_counts[idx],
            "best_neg_log10_pvalue": round(best_signal[idx], 4),
            "best_source": best_source[idx],
            "best_motif": best_motif[idx],
            "best_pvalue": best_pvalue[idx],
        }
        profile_rows.append(row)
        if hit_counts[idx] > 0:
            covered_bins += 1
            peak_rows.append(row)
    peak_rows = sorted(
        peak_rows,
        key=lambda row: (
            -float(row["best_neg_log10_pvalue"]),
            -int(row["hit_count"]),
            int(row["start"]),
        ),
    )[:10]
    source_count_rows = [
        {
            "source": source,
            "label": SOURCE_LABELS.get(source, source or "Unknown"),
            "count": count,
        }
        for source, count in sorted(source_counts.items(), key=lambda item: (-item[1], item[0]))
    ]

    width = 980
    height = 280
    left = 58
    right = 22
    top = 24
    bottom = 44
    plot_width = width - left - right
    plot_height = height - top - bottom
    baseline = top + plot_height

    def x_at(idx: int) -> float:
        if bins == 1:
            return left + plot_width / 2
        return left + (idx / (bins - 1)) * plot_width

    def y_at(value: float) -> float:
        return baseline - (value / max_signal) * plot_height

    count_rects = []
    bar_width = max(1.0, plot_width / bins)
    for idx, count in enumerate(hit_counts):
        if count <= 0:
            continue
        x = left + (idx / bins) * plot_width
        bar_height = (count / max_count) * plot_height if max_count else 0
        start_pos = int((idx / bins) * sequence_length) + 1
        end_pos = int(((idx + 1) / bins) * sequence_length)
        count_rects.append(
            f'<rect x="{x:.2f}" y="{baseline - bar_height:.2f}" width="{bar_width:.2f}" '
            f'height="{bar_height:.2f}"><title>Positions {start_pos}-{max(start_pos, end_pos)}; '
            f'{count} FIMO hit(s); best -log10(p) {best_signal[idx]:.2f}</title></rect>'
        )

    points = " ".join(f"{x_at(idx):.2f},{y_at(value):.2f}" for idx, value in enumerate(best_signal))
    half_signal = max_signal / 2
    threshold_signal = -math.log10(max(pvalue_threshold, 1e-300))
    threshold_y = y_at(min(threshold_signal, max_signal))
    threshold_line = ""
    if threshold_signal <= max_signal:
        threshold_line = (
            f'<line class="profile-threshold" x1="{left}" y1="{threshold_y:.2f}" '
            f'x2="{left + plot_width}" y2="{threshold_y:.2f}" />'
            f'<text class="axis-label" x="{left + plot_width - 4}" y="{threshold_y - 5:.2f}" '
            f'text-anchor="end">threshold</text>'
        )
    svg = f"""
<svg class="fimo-profile-svg" viewBox="0 0 {width} {height}" role="img" aria-label="FIMO binding profile">
  <title>FIMO profile: best negative log10 p-value per DNA position. Pale bars show how many FIMO hits cover each position.</title>
  <rect class="plot-bg" x="{left}" y="{top}" width="{plot_width}" height="{plot_height}" />
  <line class="axis" x1="{left}" y1="{baseline}" x2="{left + plot_width}" y2="{baseline}" />
  <line class="axis" x1="{left}" y1="{top}" x2="{left}" y2="{baseline}" />
  <line class="grid" x1="{left}" y1="{y_at(half_signal):.2f}" x2="{left + plot_width}" y2="{y_at(half_signal):.2f}" />
  <line class="grid" x1="{left}" y1="{top}" x2="{left + plot_width}" y2="{top}" />
  {threshold_line}
  <g class="profile-count-bars">
    {''.join(count_rects)}
  </g>
  <polyline class="profile-line" points="{points}" />
  <text class="axis-label" x="{left}" y="{height - 12}">1</text>
  <text class="axis-label" x="{left + plot_width}" y="{height - 12}" text-anchor="end">{sequence_length:,} bp</text>
  <text class="axis-label" x="12" y="{top + 4}" transform="rotate(-90 12,{top + 4})">best -log10(p)</text>
  <text class="axis-label" x="{left - 8}" y="{baseline + 4}" text-anchor="end">0</text>
  <text class="axis-label" x="{left - 8}" y="{y_at(half_signal) + 4:.2f}" text-anchor="end">{half_signal:.1f}</text>
  <text class="axis-label" x="{left - 8}" y="{top + 4}" text-anchor="end">{max_signal:.1f}</text>
</svg>
""".strip()
    summary = {
        "profiled_hits": profiled_hits,
        "bins": bins,
        "bin_size": max(1, math.ceil(sequence_length / bins)),
        "max_signal": max_signal,
        "max_count": max_count,
        "covered_bins": covered_bins,
        "coverage_pct": (covered_bins / bins) * 100,
        "source_counts": source_count_rows,
        "peak_rows": peak_rows,
        "profile_rows": profile_rows,
    }
    return svg, summary


def run_fimo_scan(
    motifs: list[sqlite3.Row],
    sequence: str,
    pvalue_threshold: float,
    max_hits: int,
    timeout: int = 60,
) -> tuple[list[dict[str, object]], list[str], str, dict[str, object] | None]:
    fimo_path = shutil.which("fimo")
    if not fimo_path:
        return [], ["FIMO is not installed or not available in PATH on this machine."], "", None

    motif_text, id_map = build_combined_meme(motifs)
    fasta_text = f">query\n{sequence}\n"
    errors: list[str] = []
    with tempfile.TemporaryDirectory(prefix="tf_webdb_fimo_") as tmpdir:
        tmpdir_path = Path(tmpdir)
        motif_path = tmpdir_path / "motifs.meme"
        sequence_path = tmpdir_path / "sequence.fa"
        output_dir = tmpdir_path / "fimo_out"
        motif_path.write_text(motif_text, encoding="utf-8")
        sequence_path.write_text(fasta_text, encoding="utf-8")
        try:
            completed = subprocess.run(
                [
                    fimo_path,
                    "--verbosity",
                    "1",
                    "--thresh",
                    str(pvalue_threshold),
                    "--oc",
                    str(output_dir),
                    str(motif_path),
                    str(sequence_path),
                ],
                check=False,
                capture_output=True,
                text=True,
                timeout=timeout,
            )
        except subprocess.TimeoutExpired:
            return [], [f"FIMO timed out after {timeout} seconds."], "", None

        if completed.returncode != 0:
            message = completed.stderr.strip() or completed.stdout.strip() or "FIMO failed without an error message."
            return [], [message], "", None

        result_path = output_dir / "fimo.tsv"
        if not result_path.exists():
            legacy_path = output_dir / "fimo.txt"
            if legacy_path.exists():
                result_path = legacy_path
            else:
                return [], ["FIMO finished, but no fimo.tsv result file was created."], "", None

        data_lines = [
            line
            for line in result_path.read_text(encoding="utf-8", errors="replace").splitlines()
            if line.strip() and not line.startswith("#")
        ]
        if not data_lines:
            return [], errors, "", None

        reader = csv.DictReader(data_lines, delimiter="\t")
        hits: list[dict[str, object]] = []
        for row in reader:
            raw_motif_id = row.get("motif_id") or row.get("motif_id ") or ""
            source, motif_id = id_map.get(raw_motif_id, ("", raw_motif_id))
            start = row.get("start") or row.get("start ") or ""
            stop = row.get("stop") or row.get("end") or ""
            score = row.get("score") or ""
            pvalue = row.get("p-value") or row.get("pvalue") or ""
            qvalue = row.get("q-value") or row.get("qvalue") or ""
            matched_sequence = row.get("matched_sequence") or row.get("matched sequence") or ""
            pvalue_value = parse_float_for_sort(pvalue, default=float("inf")) if pvalue else None
            qvalue_value = parse_float_for_sort(qvalue, default=float("nan")) if qvalue else None
            hits.append(
                {
                    "source": source,
                    "motif_id": motif_id,
                    "start": int(start) if str(start).isdigit() else start,
                    "end": int(stop) if str(stop).isdigit() else stop,
                    "strand": row.get("strand") or "",
                    "sequence": matched_sequence,
                    "score": parse_float_for_sort(score),
                    "relative_score": None,
                    "pvalue": pvalue_value,
                    "qvalue": qvalue_value,
                    "engine": "fimo",
                }
            )
    hits.sort(
        key=lambda hit: (
            parse_float_for_sort(str(hit.get("pvalue")), default=float("inf")),
            -parse_float_for_sort(str(hit.get("score"))),
            str(hit.get("motif_id")),
            int(hit["start"]) if isinstance(hit["start"], int) else 0,
        )
    )
    profile_svg, profile_summary = build_fimo_profile_svg(hits, len(sequence), pvalue_threshold)
    return hits[:max_hits], errors, profile_svg, profile_summary


def render_pwm_svg(matrix_json: str | None, compact: bool = False) -> str:
    if not matrix_json:
        return ""
    matrix = json.loads(matrix_json)
    if not matrix:
        return ""

    cell_w = 22 if compact else 34
    plot_h = 58 if compact else 118
    top = 8 if compact else 16
    bottom = 18 if compact else 32
    left = 8 if compact else 42
    right = 8 if compact else 16
    max_bits = 2.0
    pixels_per_bit = plot_h / max_bits
    font_size = 18 if compact else 28
    width = left * 2 + cell_w * len(matrix)
    width = left + right + cell_w * len(matrix)
    height = top + plot_h + bottom
    pieces = [
        f'<svg class="pwm-svg sequence-logo{" compact" if compact else ""}" viewBox="0 0 {width} {height}" '
        'role="img" aria-label="DNA sequence logo">'
    ]

    if not compact:
        axis_x = left - 10
        y_bottom = top + plot_h
        pieces.append(f'<line x1="{axis_x}" y1="{top}" x2="{axis_x}" y2="{y_bottom}" class="logo-axis" />')
        for bits in (0, 1, 2):
            y = y_bottom - bits * pixels_per_bit
            pieces.append(f'<line x1="{axis_x - 3}" y1="{y:.1f}" x2="{axis_x}" y2="{y:.1f}" class="logo-axis" />')
            pieces.append(
                f'<text x="{axis_x - 6}" y="{y + 3:.1f}" text-anchor="end" class="logo-axis-label">{bits}</text>'
            )
        pieces.append(f'<text x="{axis_x - 22}" y="{top + 10}" class="logo-axis-label">bits</text>')

    for pos, row in enumerate(matrix):
        x = left + pos * cell_w
        y_bottom = top + plot_h
        values = [max(0.0, float(value)) for value in row[:4]]
        total = sum(values)
        if total <= 0:
            continue
        probabilities = [value / total for value in values]
        entropy = -sum(probability * math.log2(probability) for probability in probabilities if probability > 0)
        information = max(0.0, max_bits - entropy)
        y_cursor = y_bottom
        for base, probability in sorted(zip(("A", "C", "G", "T"), probabilities), key=lambda item: item[1]):
            letter_h = max(0.0, probability * information * pixels_per_bit)
            if letter_h <= 1.0:
                continue
            y_cursor -= letter_h
            scale_y = max(letter_h / font_size, 0.01)
            pieces.append(
                f'<g transform="translate({x + cell_w / 2:.1f} {y_cursor + letter_h:.1f}) scale(1 {scale_y:.4f})">'
                f'<text class="logo-letter" text-anchor="middle" dominant-baseline="text-after-edge" '
                f'font-size="{font_size}" fill="{BASE_COLORS[base]}">{base}'
                f'<title>Position {pos + 1}, {base}: probability {probability:.3f}, stack {information:.2f} bits</title>'
                '</text></g>'
            )
        if not compact:
            pieces.append(
                f'<text x="{x + cell_w / 2:.1f}" y="{height - 5}" text-anchor="middle" '
                'class="logo-position-label">'
                f"{pos + 1}</text>"
            )
    pieces.append("</svg>")
    return "".join(pieces)


def evidence_sort_key(row: sqlite3.Row) -> tuple[int, str]:
    order = {
        "identical": 0,
        "homologous": 1,
        "relative_homologous": 2,
        "modcre": 3,
        "alphafold": 4,
    }
    return (order.get(row["evidence_type"], 99), row["motif_id"])


class TFWebApp:
    def __init__(self, db_path: Path, enable_debug: bool = False, show_errors: bool = False):
        self.db_path = db_path
        self.enable_debug = enable_debug
        self.show_errors = show_errors
        self.templates = Environment(
            loader=FileSystemLoader(APP_DIR / "templates"),
            autoescape=select_autoescape(("html", "xml")),
        )
        self.templates.filters["pwm_svg"] = render_pwm_svg
        self.templates.filters["quote"] = lambda value: quote(str(value), safe="")
        self.templates.globals["motif_external_links"] = motif_external_links
        self.templates.globals["pdb_url"] = pdb_url
        self.templates.globals["motif_spec"] = motif_spec

    def stats(self, conn: sqlite3.Connection) -> dict[str, int | str]:
        return {
            "tf_count": conn.execute("SELECT COUNT(*) FROM tf").fetchone()[0],
            "annotated_tf_count": safe_count(conn, "SELECT COUNT(*) FROM tf_annotation"),
            "motif_ref_count": conn.execute("SELECT COUNT(*) FROM motif_ref").fetchone()[0],
            "motif_file_count": conn.execute("SELECT COUNT(*) FROM motif_file").fetchone()[0],
            "active_model_count": conn.execute(
                "SELECT COUNT(*) FROM structure_file WHERE status = 'active' AND file_type = 'pdb'"
            ).fetchone()[0],
            "model_summary_count": safe_count(conn, "SELECT COUNT(*) FROM model_summary"),
            "missing_motif_count": conn.execute(
                "SELECT COUNT(*) FROM motif_ref WHERE missing_local_file = 1"
            ).fetchone()[0],
        }

    def render(self, template_name: str, **context: object) -> bytes:
        template = self.templates.get_template(template_name)
        context.setdefault("evidence_labels", EVIDENCE_LABELS)
        context.setdefault("evidence_info", EVIDENCE_INFO)
        context.setdefault("source_labels", SOURCE_LABELS)
        context.setdefault("debug_enabled", self.enable_debug)
        html_text = template.render(**context)
        return html_text.encode("utf-8")

    def index(self) -> tuple[bytes, str, int]:
        with connect(self.db_path) as conn:
            stats = self.stats(conn)
            examples = conn.execute(
                """
                SELECT tf.tf_id, tf.family_text, tf.motif_ref_count, tf.active_model_count,
                       ta.gene_names, ta.protein_name, ta.organism_name, ta.reviewed
                FROM tf
                LEFT JOIN tf_annotation AS ta ON ta.tf_id = tf.tf_id
                WHERE motif_ref_count > 0
                ORDER BY tf.active_model_count DESC, tf.motif_ref_count DESC, tf.tf_id
                LIMIT 8
                """
            ).fetchall()
        return self.render("index.html", stats=stats, examples=examples), "text/html", 200

    def search(self, params: dict[str, list[str]]) -> tuple[bytes, str, int]:
        q = params.get("q", [""])[0].strip()
        evidence = params.get("evidence", [""])[0].strip()
        source = params.get("source", [""])[0].strip()
        limit = 60

        where = []
        args: list[object] = []
        joins = """
        LEFT JOIN motif_ref AS mr ON mr.tf_id = tf.tf_id
        LEFT JOIN tf_annotation AS ta ON ta.tf_id = tf.tf_id
        """
        if q:
            like = f"%{q}%"
            where.append(
                """
                (tf.tf_id LIKE ? OR tf.family_text LIKE ? OR mr.motif_id LIKE ?
                 OR ta.gene_names LIKE ? OR ta.protein_name LIKE ? OR ta.organism_name LIKE ?)
                """
            )
            args.extend([like, like, like, like, like, like])
        if evidence:
            where.append("mr.evidence_type = ?")
            args.append(evidence)
        if source:
            where.append("mr.source = ?")
            args.append(source)
        where_sql = "WHERE " + " AND ".join(where) if where else ""

        with connect(self.db_path) as conn:
            rows = conn.execute(
                f"""
                SELECT DISTINCT tf.tf_id, tf.family_text, tf.motif_ref_count, tf.active_model_count,
                       ta.gene_names, ta.protein_name, ta.organism_name, ta.reviewed
                FROM tf
                {joins}
                {where_sql}
                ORDER BY tf.active_model_count DESC, tf.motif_ref_count DESC, tf.tf_id
                LIMIT ?
                """,
                [*args, limit],
            ).fetchall()
            total = conn.execute(
                f"""
                SELECT COUNT(*) FROM (
                    SELECT DISTINCT tf.tf_id
                    FROM tf
                    {joins}
                    {where_sql}
                )
                """,
                args,
            ).fetchone()[0]
        return (
            self.render(
                "search.html",
                q=q,
                evidence=evidence,
                source=source,
                rows=rows,
                total=total,
                limit=limit,
            ),
            "text/html",
            200,
        )

    def tf_detail(self, tf_id: str) -> tuple[bytes, str, int]:
        tf_id = unquote(tf_id)
        with connect(self.db_path) as conn:
            tf = dict_row(
                conn.execute(
                    """
                    SELECT tf.*, ta.entry_name, ta.gene_names, ta.protein_name,
                           ta.organism_name, ta.organism_id, ta.reviewed,
                           ta.sequence_length, ta.annotation_score, ta.uniprot_url,
                           ta.fetched_at
                    FROM tf
                    LEFT JOIN tf_annotation AS ta ON ta.tf_id = tf.tf_id
                    WHERE tf.tf_id = ?
                    """,
                    (tf_id,),
                ).fetchone()
            )
            if not tf:
                return self.not_found(f"Unknown TF: {html.escape(tf_id)}")
            families = conn.execute(
                "SELECT family FROM tf_family WHERE tf_id = ? ORDER BY family", (tf_id,)
            ).fetchall()
            motif_rows = conn.execute(
                """
                SELECT mr.*, mf.consensus, mf.width, mf.nsites, mf.matrix_json,
                       region.region_start, region.region_end, region.region_model_count
                FROM motif_ref AS mr
                LEFT JOIN motif_file AS mf
                  ON mf.source = mr.source
                 AND mf.motif_id = mr.motif_id
                LEFT JOIN (
                    SELECT ms.motif_ref_id,
                           MIN(sf.residue_start) AS region_start,
                           MAX(sf.residue_end) AS region_end,
                           COUNT(DISTINCT sf.id) AS region_model_count
                    FROM motif_structure AS ms
                    JOIN structure_file AS sf ON sf.id = ms.structure_file_id
                    WHERE sf.status = 'active'
                      AND sf.residue_start IS NOT NULL
                      AND sf.residue_end IS NOT NULL
                    GROUP BY ms.motif_ref_id
                ) AS region ON region.motif_ref_id = mr.id
                WHERE mr.tf_id = ?
                ORDER BY mr.evidence_type, mr.identity_percent DESC, mr.motif_id
                """,
                (tf_id,),
            ).fetchall()
            active_models = conn.execute(
                """
                SELECT * FROM structure_file
                WHERE tf_id = ? AND status = 'active' AND file_type = 'pdb'
                ORDER BY source, model_id
                LIMIT 100
                """,
                (tf_id,),
            ).fetchall()
            model_summaries = conn.execute(
                """
                SELECT ms.*, sf.id AS matched_file_id, sf.model_id AS matched_model_id,
                       sf.file_type AS matched_file_type
                FROM model_summary AS ms
                LEFT JOIN structure_file AS sf ON sf.id = ms.matched_structure_id
                WHERE ms.tf_id = ? AND ms.status = 'active'
                ORDER BY ms.identity_percent DESC, ms.similarity_percent DESC, ms.model_rank
                LIMIT 12
                """,
                (tf_id,),
            ).fetchall()

        grouped: dict[str, list[sqlite3.Row]] = {key: [] for key in EVIDENCE_LABELS}
        for row in sorted(motif_rows, key=evidence_sort_key):
            grouped.setdefault(row["evidence_type"], []).append(row)
        region_groups = build_region_groups(tf, list(motif_rows), list(active_models), list(model_summaries))

        return (
            self.render(
                "tf.html",
                tf=tf,
                families=families,
                grouped=grouped,
                region_groups=region_groups,
                active_models=active_models,
                model_summaries=model_summaries,
                visible_limit=100,
            ),
            "text/html",
            200,
        )

    def motif_detail(self, params: dict[str, list[str]]) -> tuple[bytes, str, int]:
        source = params.get("source", [""])[0].strip()
        motif_id = params.get("id", [""])[0].strip()
        if not source or not motif_id:
            return self.not_found("Motif source and id are required.")
        with connect(self.db_path) as conn:
            motif = dict_row(
                conn.execute(
                    "SELECT * FROM motif_file WHERE source = ? AND motif_id = ?",
                    (source, motif_id),
                ).fetchone()
            )
            refs = conn.execute(
                """
                SELECT mr.*, tf.family_text
                FROM motif_ref AS mr
                JOIN tf ON tf.tf_id = mr.tf_id
                WHERE mr.source = ? AND mr.motif_id = ?
                ORDER BY mr.evidence_type, mr.identity_percent DESC, mr.tf_id
                LIMIT 100
                """,
                (source, motif_id),
            ).fetchall()
            structures = conn.execute(
                """
                SELECT DISTINCT sf.*
                FROM motif_ref AS mr
                JOIN motif_structure AS ms ON ms.motif_ref_id = mr.id
                JOIN structure_file AS sf ON sf.id = ms.structure_file_id
                WHERE mr.source = ? AND mr.motif_id = ? AND sf.status = 'active'
                ORDER BY sf.file_type, sf.model_id
                """,
                (source, motif_id),
            ).fetchall()
        if not motif and not refs:
            return self.not_found(f"Unknown motif: {html.escape(source)} / {html.escape(motif_id)}")
        return (
            self.render("motif.html", motif=motif, refs=refs, structures=structures, source=source, motif_id=motif_id),
            "text/html",
            200,
        )

    def scan(self, params: dict[str, list[str]], method: str = "GET") -> tuple[bytes, str, int]:
        source = params.get("source", [""])[0].strip()
        motif_id = params.get("id", [""])[0].strip()
        default_motifs = motif_spec(source, motif_id) if source and motif_id else ""
        motifs_text = params.get("motifs", [default_motifs])[0].strip()
        tf_id = normalize_tf_id(params.get("tf_id", [""])[0])
        tf_evidence = selected_evidence_from_params(params, tf_id)
        tf_region_raw = params.get("tf_region", [""])[0].strip()
        tf_region = parse_region_value(tf_region_raw)
        sequence_text = params.get("sequence", [""])[0]
        motif_q = params.get("motif_q", [""])[0].strip()
        motif_source = params.get("motif_source", [""])[0].strip().lower()
        if motif_source not in SOURCE_LABELS:
            motif_source = ""
        engine = "fimo"
        pvalue_raw = params.get("pvalue", ["1e-4"])[0]
        max_hits_raw = params.get("max_hits", ["200"])[0]
        errors: list[str] = []

        try:
            pvalue = float(pvalue_raw)
        except ValueError:
            pvalue = 1e-4
            errors.append("FIMO p-value threshold was not numeric; using 1e-4.")
        pvalue = min(max(pvalue, 1e-12), 1.0)

        try:
            max_hits = int(max_hits_raw)
        except ValueError:
            max_hits = 200
            errors.append("Max hits was not numeric; using 200.")
        max_hits = min(max(max_hits, 1), 1000)

        motifs: list[sqlite3.Row] = []
        seen_motifs: set[tuple[str, str]] = set()
        tf_scan_summary: dict[str, object] | None = None
        tf_scan_regions: list[dict[str, object]] = []

        def add_motif(row: sqlite3.Row) -> None:
            key = (row["source"], row["motif_id"])
            if key not in seen_motifs:
                motifs.append(row)
                seen_motifs.add(key)

        if tf_id:
            with connect(self.db_path) as conn:
                tf_scan_regions = load_tf_scan_regions(conn, tf_id)
                available_region_values = {str(region["value"]) for region in tf_scan_regions}
                if tf_region_raw and tf_region_raw not in available_region_values:
                    errors.append(f"DNA-binding region is not available for this TF: {tf_region_raw}")
                    tf_region = None
                    tf_region_raw = ""
                tf, tf_motifs, summary = load_tf_scan_motifs(conn, tf_id, tf_evidence, region=tf_region)
            if not tf:
                errors.append(f"TF accession not found: {tf_id}")
            else:
                tf_scan_summary = summary
                for row in tf_motifs:
                    add_motif(row)
                if summary.get("limited"):
                    errors.append(
                        f"TF motif set was limited to {summary['limit']} motifs. "
                        "Use stricter evidence filters for a smaller scan."
                    )
                if not tf_motifs:
                    if tf_region:
                        errors.append(
                            "No FIMO-ready motifs were found for this TF with the selected evidence and region filters."
                        )
                    else:
                        errors.append("No FIMO-ready motifs were found for this TF with the selected evidence filters.")

        requested_specs = parse_motif_specs(motifs_text)
        if requested_specs:
            with connect(self.db_path) as conn:
                for requested_source, requested_id in requested_specs:
                    if requested_source:
                        row = conn.execute(
                            """
                            SELECT source, motif_id, width, nsites, consensus, matrix_json
                            FROM motif_file
                            WHERE source = ? AND motif_id = ? AND matrix_json IS NOT NULL
                            """,
                            (requested_source, requested_id),
                        ).fetchone()
                    else:
                        row = conn.execute(
                            """
                            SELECT source, motif_id, width, nsites, consensus, matrix_json
                            FROM motif_file
                            WHERE motif_id = ? AND matrix_json IS NOT NULL
                            ORDER BY source
                            LIMIT 1
                            """,
                            (requested_id,),
                        ).fetchone()
                    if row:
                        add_motif(row)
                    else:
                        label = f"{requested_source}|{requested_id}" if requested_source else requested_id
                        errors.append(f"Motif not found or has no matrix: {label}")

        with connect(self.db_path) as conn:
            motif_results = search_scan_motifs(conn, motif_q, motif_source)

        sequence, ignored_count = clean_dna_sequence(sequence_text)
        if ignored_count:
            errors.append(f"Ignored {ignored_count} non-DNA characters.")

        hits: list[dict[str, object]] = []
        profile_svg = ""
        profile_summary: dict[str, object] | None = None
        hits_csv_uri = ""
        profile_csv_uri = ""
        if method == "POST":
            if not motifs:
                errors.append("Choose at least one motif.")
            if not sequence:
                errors.append("Paste a DNA sequence.")
            if len(sequence) > 100000:
                errors.append("Sequence is longer than 100,000 bases; please use a shorter sequence for this interactive scanner.")
            if motifs and sequence and len(sequence) <= 100000:
                hits, fimo_errors, profile_svg, profile_summary = run_fimo_scan(motifs, sequence, pvalue, max_hits)
                errors.extend(fimo_errors)
                if hits:
                    hits_csv_uri = build_hits_csv_uri(hits)
                if profile_summary:
                    profile_csv_uri = build_profile_csv_uri(profile_summary)

        return (
            self.render(
                "scan.html",
                motifs_text=motifs_text,
                tf_id=tf_id,
                tf_evidence=tf_evidence,
                tf_region=tf_region_raw,
                tf_scan_regions=tf_scan_regions,
                tf_scan_summary=tf_scan_summary,
                tf_scan_default_evidence=TF_SCAN_DEFAULT_EVIDENCE,
                tf_scan_max_motifs=TF_SCAN_MAX_MOTIFS,
                sequence_text=sequence_text,
                motif_q=motif_q,
                motif_source=motif_source,
                motif_results=motif_results,
                engine=engine,
                pvalue=pvalue,
                max_hits=max_hits,
                motifs=motifs,
                hits=hits,
                profile_svg=profile_svg,
                profile_summary=profile_summary,
                hits_csv_uri=hits_csv_uri,
                profile_csv_uri=profile_csv_uri,
                errors=errors,
                sequence_length=len(sequence),
                method=method,
            ),
            "text/html",
            200,
        )

    def model_summaries(self, tf_id: str) -> tuple[bytes, str, int]:
        tf_id = unquote(tf_id)
        with connect(self.db_path) as conn:
            tf = dict_row(
                conn.execute(
                    """
                    SELECT tf.*, ta.gene_names, ta.protein_name, ta.organism_name, ta.uniprot_url
                    FROM tf
                    LEFT JOIN tf_annotation AS ta ON ta.tf_id = tf.tf_id
                    WHERE tf.tf_id = ?
                    """,
                    (tf_id,),
                ).fetchone()
            )
            if not tf:
                return self.not_found(f"Unknown TF: {html.escape(tf_id)}")
            rows = conn.execute(
                """
                SELECT ms.*, sf.id AS matched_file_id, sf.model_id AS matched_model_id,
                       sf.file_type AS matched_file_type
                FROM model_summary AS ms
                LEFT JOIN structure_file AS sf ON sf.id = ms.matched_structure_id
                WHERE ms.tf_id = ? AND ms.status = 'active'
                ORDER BY ms.identity_percent DESC, ms.similarity_percent DESC, ms.model_rank
                LIMIT 500
                """,
                (tf_id,),
            ).fetchall()
            total = conn.execute(
                "SELECT COUNT(*) FROM model_summary WHERE tf_id = ? AND status = 'active'",
                (tf_id,),
            ).fetchone()[0]
        return (
            self.render("model_summaries.html", tf=tf, rows=rows, total=total, limit=500),
            "text/html",
            200,
        )

    def debug(self) -> tuple[bytes, str, int]:
        with connect(self.db_path) as conn:
            issues = conn.execute(
                """
                SELECT category, source, motif_id, COUNT(*) AS count,
                       MIN(tf_id) AS example_tf_id
                FROM import_issue
                GROUP BY category, source, motif_id
                ORDER BY count DESC, source, motif_id
                LIMIT 200
                """
            ).fetchall()
            failed_counts = conn.execute(
                """
                SELECT source, file_type, COUNT(*) AS count,
                       MIN(member_path) AS example_member_path
                FROM structure_file
                WHERE status = 'failed'
                GROUP BY source, file_type
                ORDER BY source, file_type
                """
            ).fetchall()
            metadata = conn.execute("SELECT key, value FROM metadata ORDER BY key").fetchall()
        return self.render("debug.html", issues=issues, failed_counts=failed_counts, metadata=metadata), "text/html", 200

    def docs(self) -> tuple[bytes, str, int]:
        with connect(self.db_path) as conn:
            stats = self.stats(conn)
            stats["unique_missing_motif_count"] = safe_count(
                conn,
                """
                SELECT COUNT(*) FROM (
                    SELECT DISTINCT source, motif_id
                    FROM motif_ref
                    WHERE missing_local_file = 1
                )
                """,
            )
            stats["failed_model_file_count"] = safe_count(
                conn,
                "SELECT COUNT(*) FROM structure_file WHERE status = 'failed'",
            )
            stats["fimo_ready_motif_count"] = safe_count(
                conn,
                "SELECT COUNT(*) FROM motif_file WHERE matrix_json IS NOT NULL",
            )
            stats["tf_with_usable_pwm_count"] = safe_count(
                conn,
                """
                SELECT COUNT(*) FROM (
                    SELECT DISTINCT mr.tf_id
                    FROM motif_ref AS mr
                    JOIN motif_file AS mf
                      ON mf.source = mr.source
                     AND mf.motif_id = mr.motif_id
                    WHERE mf.matrix_json IS NOT NULL
                )
                """,
            )
            stats["tf_with_active_model_count"] = safe_count(
                conn,
                """
                SELECT COUNT(*) FROM (
                    SELECT DISTINCT tf_id
                    FROM structure_file
                    WHERE status = 'active'
                      AND file_type = 'pdb'
                      AND tf_id IS NOT NULL
                )
                """,
            )
            motif_files_by_source = rows_with_percent(
                conn.execute(
                    """
                    SELECT source, COUNT(*) AS count,
                           SUM(CASE WHEN matrix_json IS NOT NULL THEN 1 ELSE 0 END) AS with_matrix,
                           SUM(CASE WHEN matrix_json IS NULL THEN 1 ELSE 0 END) AS without_matrix,
                           SUM(CASE WHEN width = 0 THEN 1 ELSE 0 END) AS width_zero
                    FROM motif_file
                    GROUP BY source
                    ORDER BY count DESC, source
                    """
                ).fetchall()
            )
            motif_links_by_evidence = rows_with_percent(
                conn.execute(
                    """
                    SELECT evidence_type, COUNT(*) AS count,
                           SUM(CASE WHEN missing_local_file = 1 THEN 1 ELSE 0 END) AS missing_count
                    FROM motif_ref
                    GROUP BY evidence_type
                    ORDER BY count DESC, evidence_type
                    """
                ).fetchall()
            )
            motif_links_by_source = rows_with_percent(
                conn.execute(
                    """
                    SELECT source, COUNT(*) AS count,
                           SUM(CASE WHEN missing_local_file = 1 THEN 1 ELSE 0 END) AS missing_count
                    FROM motif_ref
                    GROUP BY source
                    ORDER BY count DESC, source
                    """
                ).fetchall()
            )
            model_files_by_group = rows_with_percent(
                conn.execute(
                    """
                    SELECT source, status, file_type, COUNT(*) AS count
                    FROM structure_file
                    GROUP BY source, status, file_type
                    ORDER BY source, status, file_type
                    """
                ).fetchall()
            )
            model_summary_by_status = rows_with_percent(
                conn.execute(
                    """
                    SELECT status, COUNT(*) AS count,
                           SUM(CASE WHEN matched_structure_id IS NOT NULL THEN 1 ELSE 0 END) AS linked_count
                    FROM model_summary
                    GROUP BY status
                    ORDER BY status
                    """
                ).fetchall()
            )
            annotation_status = rows_with_percent(
                conn.execute(
                    """
                    SELECT CASE reviewed WHEN 1 THEN 'Reviewed Swiss-Prot'
                                         ELSE 'Unreviewed TrEMBL' END AS label,
                           COUNT(*) AS count
                    FROM tf_annotation
                    GROUP BY reviewed
                    ORDER BY reviewed DESC
                    """
                ).fetchall()
            )
            top_families = rows_with_percent(
                conn.execute(
                    """
                    SELECT family AS label, COUNT(*) AS count
                    FROM tf_family
                    GROUP BY family
                    ORDER BY count DESC, family
                    LIMIT 12
                    """
                ).fetchall()
            )
            missing_motifs = conn.execute(
                """
                SELECT source, motif_id, COUNT(*) AS count, MIN(tf_id) AS example_tf_id
                FROM motif_ref
                WHERE missing_local_file = 1
                GROUP BY source, motif_id
                ORDER BY count DESC, source, motif_id
                """
            ).fetchall()
            metadata = conn.execute("SELECT key, value FROM metadata ORDER BY key").fetchall()
        return (
            self.render(
                "docs.html",
                stats=stats,
                motif_files_by_source=motif_files_by_source,
                motif_links_by_evidence=motif_links_by_evidence,
                motif_links_by_source=motif_links_by_source,
                model_files_by_group=model_files_by_group,
                model_summary_by_status=model_summary_by_status,
                annotation_status=annotation_status,
                top_families=top_families,
                missing_motifs=missing_motifs,
                metadata=metadata,
            ),
            "text/html",
            200,
        )

    def evidence(self) -> tuple[bytes, str, int]:
        return self.render("evidence.html"), "text/html", 200

    def download_motif(self, params: dict[str, list[str]]) -> tuple[bytes, str, int, dict[str, str]]:
        source = params.get("source", [""])[0].strip()
        motif_id = params.get("id", [""])[0].strip()
        with connect(self.db_path) as conn:
            row = conn.execute(
                "SELECT content FROM motif_file WHERE source = ? AND motif_id = ?",
                (source, motif_id),
            ).fetchone()
        if not row:
            body, content_type, status = self.not_found("Motif file is not available locally.")
            return body, content_type, status, {}
        filename = f"{motif_id}.meme"
        return (
            row["content"].encode("utf-8"),
            "text/plain; charset=utf-8",
            200,
            {"Content-Disposition": f'attachment; filename="{filename}"'},
        )

    def get_model_file(self, file_id: str) -> tuple[sqlite3.Row | None, bytes | None, str | None]:
        if not file_id.isdigit():
            return None, None, None
        with connect(self.db_path) as conn:
            row = conn.execute("SELECT * FROM structure_file WHERE id = ?", (int(file_id),)).fetchone()
        if not row:
            return None, None, None
        with tarfile.open(row["archive_path"], "r:gz") as archive:
            extracted = archive.extractfile(row["member_path"])
            if extracted is None:
                return row, None, None
            content = extracted.read()
        return row, content, Path(row["member_path"]).name

    def model_viewer(self, params: dict[str, list[str]]) -> tuple[bytes, str, int]:
        file_id = params.get("id", [""])[0].strip()
        if not file_id.isdigit():
            return self.not_found("Model id is required.")
        with connect(self.db_path) as conn:
            model = dict_row(
                conn.execute(
                    """
                    SELECT sf.*, ta.gene_names, ta.protein_name, ta.organism_name
                    FROM structure_file AS sf
                    LEFT JOIN tf_annotation AS ta ON ta.tf_id = sf.tf_id
                    WHERE sf.id = ?
                    """,
                    (int(file_id),),
                ).fetchone()
            )
            if not model:
                return self.not_found("Model file is not indexed.")
            summary = conn.execute(
                """
                SELECT * FROM model_summary
                WHERE matched_structure_id = ?
                ORDER BY identity_percent DESC, similarity_percent DESC, model_rank
                LIMIT 1
                """,
                (int(file_id),),
            ).fetchone()
        if model["file_type"] != "pdb":
            return self.not_found("Only PDB files can be viewed in 3D.")
        return self.render("model.html", model=model, summary=summary), "text/html", 200

    def model_data(self, params: dict[str, list[str]]) -> tuple[bytes, str, int, dict[str, str]]:
        file_id = params.get("id", [""])[0].strip()
        row, content, filename = self.get_model_file(file_id)
        if not row:
            body, content_type, status = self.not_found("Model file is not indexed.")
            return body, content_type, status, {}
        if content is None or filename is None:
            body, content_type, status = self.not_found("Model file could not be extracted.")
            return body, content_type, status, {}
        return content, "chemical/x-pdb" if filename.endswith(".pdb") else "text/plain", 200, {
            "Cache-Control": "no-store"
        }

    def download_model(self, params: dict[str, list[str]]) -> tuple[bytes, str, int, dict[str, str]]:
        file_id = params.get("id", [""])[0].strip()
        row, content, filename = self.get_model_file(file_id)
        if not row:
            body, content_type, status = self.not_found("Model file is not indexed.")
            return body, content_type, status, {}
        if content is None or filename is None:
            body, content_type, status = self.not_found("Model file could not be extracted.")
            return body, content_type, status, {}
        return (
            content,
            "chemical/x-pdb" if filename.endswith(".pdb") else "text/plain",
            200,
            {"Content-Disposition": f'attachment; filename="{filename}"'},
        )

    def static_file(self, path: str) -> tuple[bytes, str, int]:
        rel_path = unquote(path).removeprefix("/static/")
        file_path = (APP_DIR / "static" / rel_path).resolve()
        static_root = (APP_DIR / "static").resolve()
        if static_root not in file_path.parents and file_path != static_root:
            return self.not_found("Invalid static path.")
        if not file_path.exists() or not file_path.is_file():
            return self.not_found("Static file not found.")
        return file_path.read_bytes(), mimetypes.guess_type(file_path.name)[0] or "application/octet-stream", 200

    def not_found(self, message: str) -> tuple[bytes, str, int]:
        return self.render("error.html", title="Not found", message=message), "text/html", 404

    def request_too_large(self, max_bytes: int) -> tuple[bytes, str, int]:
        message = f"Request too large. Maximum accepted POST body is {max_bytes:,} bytes."
        return self.render("error.html", title="Request too large", message=message), "text/html", 413

    def server_error(self, exc: Exception) -> tuple[bytes, str, int]:
        traceback.print_exc()
        message = html.escape(str(exc)) if self.show_errors else "Internal server error. Please contact the database maintainer if this persists."
        return self.render("error.html", title="Server error", message=message), "text/html", 500


def make_handler(app: TFWebApp):
    class Handler(BaseHTTPRequestHandler):
        def send_body(
            self,
            body: bytes,
            content_type: str,
            status: int,
            headers: dict[str, str] | None = None,
        ) -> None:
            self.send_response(status)
            self.send_header("Content-Type", content_type)
            self.send_header("Content-Length", str(len(body)))
            for key, value in (headers or {}).items():
                self.send_header(key, value)
            self.end_headers()
            try:
                self.wfile.write(body)
            except (BrokenPipeError, ConnectionResetError):
                print(f"{self.address_string()} - client closed connection while sending response")

        def do_GET(self) -> None:
            parsed = urlparse(self.path)
            params = parse_qs(parsed.query)
            headers: dict[str, str] = {}
            try:
                if parsed.path == "/":
                    body, content_type, status = app.index()
                elif parsed.path == "/search":
                    body, content_type, status = app.search(params)
                elif parsed.path.startswith("/tf/"):
                    body, content_type, status = app.tf_detail(parsed.path.removeprefix("/tf/"))
                elif parsed.path == "/motif":
                    body, content_type, status = app.motif_detail(params)
                elif parsed.path == "/scan":
                    body, content_type, status = app.scan(params)
                elif parsed.path.startswith("/model-summaries/"):
                    body, content_type, status = app.model_summaries(
                        parsed.path.removeprefix("/model-summaries/")
                    )
                elif parsed.path == "/model":
                    body, content_type, status = app.model_viewer(params)
                elif parsed.path == "/model-data":
                    body, content_type, status, headers = app.model_data(params)
                elif parsed.path == "/debug":
                    if app.enable_debug:
                        body, content_type, status = app.debug()
                    else:
                        body, content_type, status = app.not_found("Page not found.")
                elif parsed.path == "/docs":
                    body, content_type, status = app.docs()
                elif parsed.path == "/evidence":
                    body, content_type, status = app.evidence()
                elif parsed.path == "/download/motif":
                    body, content_type, status, headers = app.download_motif(params)
                elif parsed.path == "/download/model":
                    body, content_type, status, headers = app.download_model(params)
                elif parsed.path.startswith("/static/"):
                    body, content_type, status = app.static_file(parsed.path)
                else:
                    body, content_type, status = app.not_found("Page not found.")
            except Exception as exc:  # pragma: no cover - defensive response
                body, content_type, status = app.server_error(exc)

            self.send_body(body, content_type, status, headers)

        def do_POST(self) -> None:
            parsed = urlparse(self.path)
            headers: dict[str, str] = {}
            try:
                content_length = int(self.headers.get("Content-Length", "0") or "0")
            except ValueError:
                body = app.render("error.html", title="Bad request", message="Invalid Content-Length header.")
                self.send_body(body, "text/html", 400, headers)
                return
            if content_length > self.server.max_post_bytes:
                body, content_type, status = app.request_too_large(self.server.max_post_bytes)
                self.send_body(body, content_type, status, headers)
                return
            raw_body = self.rfile.read(content_length).decode("utf-8", errors="replace")
            params = parse_qs(raw_body)
            try:
                if parsed.path == "/scan":
                    body, content_type, status = app.scan(params, method="POST")
                else:
                    body, content_type, status = app.not_found("Page not found.")
            except Exception as exc:  # pragma: no cover - defensive response
                body, content_type, status = app.server_error(exc)

            self.send_body(body, content_type, status, headers)

        def log_message(self, fmt: str, *args: object) -> None:
            print(f"{self.address_string()} - {fmt % args}")

    return Handler


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the local TF web database.")
    parser.add_argument("--db", type=Path, default=DEFAULT_DB)
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8090)
    parser.add_argument("--max-post-bytes", type=int, default=DEFAULT_MAX_POST_BYTES)
    parser.add_argument("--enable-debug", action="store_true", help="Expose the internal /debug page. Keep disabled in production.")
    parser.add_argument("--show-errors", action="store_true", help="Show exception messages in HTTP 500 responses. Keep disabled in production.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not args.db.exists():
        raise SystemExit(f"Database does not exist: {args.db}. Run import_db.py first.")
    app = TFWebApp(args.db, enable_debug=args.enable_debug, show_errors=args.show_errors)
    server = ThreadingHTTPServer((args.host, args.port), make_handler(app))
    server.max_post_bytes = args.max_post_bytes
    print(f"TF web database running at http://{args.host}:{args.port}/")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nStopping.")


if __name__ == "__main__":
    main()
