#!/usr/bin/env python3
"""Build homepage TF-family tree JSON from the ModCREDB SQLite database.

The family tree is a display/navigation layer. It maps clean family labels to
existing family/domain tokens in tf_family.family and tf.family_text, computes
counts from the SQLite database, and writes a static JSON file consumed by the
homepage JavaScript.
"""

from __future__ import annotations

import argparse
import json
import sqlite3
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DB = ROOT / "data" / "tf_webdb.sqlite"
DEFAULT_TAXONOMY = ROOT / "data_sources" / "tf_family_tree.json"
DEFAULT_OUT = ROOT / "static" / "family_tree_data.generated.json"

EVIDENCE_TO_PRIMARY = {
    "identical": "Identical_PWM",
    "homologous": "Homologous_PWM",
    "relative_homologous": "Relatively_Homologous_PWM",
    "modcre": "ModCRE",
    "alphafold": "AlphaFold",
}


def connect(db_path: Path) -> sqlite3.Connection:
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    return conn


def table_exists(conn: sqlite3.Connection, table_name: str) -> bool:
    return conn.execute(
        "SELECT 1 FROM sqlite_master WHERE type = 'table' AND name = ?",
        (table_name,),
    ).fetchone() is not None


def table_columns(conn: sqlite3.Connection, table_name: str) -> set[str]:
    if not table_exists(conn, table_name):
        return set()
    return {str(row["name"]) for row in conn.execute(f"PRAGMA table_info({table_name})").fetchall()}


def placeholders(values: list[str]) -> str:
    return ", ".join("?" for _ in values)


def family_tf_ids(conn: sqlite3.Connection, tokens: list[str]) -> set[str]:
    clean = sorted({token.strip() for token in tokens if token and token.strip()})
    if not clean:
        return set()

    ids: set[str] = set()
    if table_exists(conn, "tf_family"):
        rows = conn.execute(
            f"SELECT DISTINCT tf_id FROM tf_family WHERE family IN ({placeholders(clean)})",
            clean,
        ).fetchall()
        ids.update(str(row["tf_id"]) for row in rows)

    clauses = []
    args: list[str] = []
    for token in clean:
        clauses.append("family_text LIKE ?")
        args.append(f"%{token}%")
    if clauses:
        rows = conn.execute(
            f"SELECT DISTINCT tf_id FROM tf WHERE {' OR '.join(clauses)}",
            args,
        ).fetchall()
        ids.update(str(row["tf_id"]) for row in rows)
    return ids


def count_primary_levels(conn: sqlite3.Connection, ids: set[str]) -> dict[str, int]:
    counts = {
        "known_count": 0,
        "homologous_count": 0,
        "relative_homologous_count": 0,
        "modcre_count": 0,
        "alphafold_count": 0,
        "unannotated_count": 0,
        "predicted_low_count": 0,
    }
    if not ids:
        return counts

    ordered = sorted(ids)
    if table_exists(conn, "tf_primary_annotation"):
        rows = conn.execute(
            f"""
            SELECT best_annotation_level, COUNT(*) AS count
            FROM tf_primary_annotation
            WHERE tf_id IN ({placeholders(ordered)})
            GROUP BY best_annotation_level
            """,
            ordered,
        ).fetchall()
        by_level = {str(row["best_annotation_level"]): int(row["count"] or 0) for row in rows}
    elif table_exists(conn, "motif_ref"):
        rows = conn.execute(
            f"""
            SELECT tf_id, evidence_type
            FROM motif_ref
            WHERE tf_id IN ({placeholders(ordered)})
            """,
            ordered,
        ).fetchall()
        priority = {"Identical_PWM": 0, "Homologous_PWM": 1, "Relatively_Homologous_PWM": 2, "ModCRE": 3, "AlphaFold": 4}
        best: dict[str, str] = {}
        for row in rows:
            level = EVIDENCE_TO_PRIMARY.get(str(row["evidence_type"]), "Unannotated")
            tf_id = str(row["tf_id"])
            if tf_id not in best or priority.get(level, 99) < priority.get(best[tf_id], 99):
                best[tf_id] = level
        by_level: dict[str, int] = {}
        for tf_id in ordered:
            level = best.get(tf_id, "Unannotated")
            by_level[level] = by_level.get(level, 0) + 1
    else:
        by_level = {"Unannotated": len(ids)}

    counts["known_count"] = by_level.get("Identical_PWM", 0)
    counts["homologous_count"] = by_level.get("Homologous_PWM", 0)
    counts["relative_homologous_count"] = by_level.get("Relatively_Homologous_PWM", 0)
    counts["modcre_count"] = by_level.get("ModCRE", 0)
    counts["alphafold_count"] = by_level.get("AlphaFold", 0)
    counts["unannotated_count"] = by_level.get("Unannotated", 0)
    counts["predicted_low_count"] = counts["modcre_count"] + counts["alphafold_count"]
    return counts


def generated_pwm_tf_count(conn: sqlite3.Connection, ids: set[str]) -> int:
    if not ids or not table_exists(conn, "motif_ref") or not table_exists(conn, "motif_file"):
        return 0
    ordered = sorted(ids)
    motif_ref_columns = table_columns(conn, "motif_ref")
    filters = [f"mr.tf_id IN ({placeholders(ordered)})", "mf.matrix_status = 'usable'"]
    if "missing_local_file" in motif_ref_columns:
        filters.append("mr.missing_local_file = 0")
    return int(conn.execute(
        f"""
        SELECT COUNT(*) FROM (
            SELECT DISTINCT mr.tf_id
            FROM motif_ref AS mr
            JOIN motif_file AS mf ON mf.source = mr.source AND mf.motif_id = mr.motif_id
            WHERE {' AND '.join(filters)}
        )
        """,
        ordered,
    ).fetchone()[0] or 0)


def model_counts(conn: sqlite3.Connection, ids: set[str]) -> dict[str, int]:
    counts = {"model_tf_count": 0, "monomer_model_tf_count": 0, "dimer_model_tf_count": 0}
    if not ids or not table_exists(conn, "structure_file"):
        return counts
    columns = table_columns(conn, "structure_file")
    if "tf_id" not in columns:
        return counts

    path_column = next((c for c in ("path", "file_path", "local_path", "model_path", "pdb_path", "structure_path", "source_path", "relative_path", "filename") if c in columns), None)
    ordered = sorted(ids)
    where = [f"tf_id IN ({placeholders(ordered)})"]
    if "status" in columns:
        where.append("status = 'active'")
    if "file_type" in columns:
        where.append("file_type = 'pdb'")

    path_expr = f"{path_column} AS model_path" if path_column else "'' AS model_path"
    rows = conn.execute(
        f"""
        SELECT tf_id, {path_expr}
        FROM structure_file
        WHERE {' AND '.join(where)}
        """,
        ordered,
    ).fetchall()

    model_ids = {str(row["tf_id"]) for row in rows}
    monomer_ids = {str(row["tf_id"]) for row in rows if "monomer" in str(row["model_path"] or "").lower()}
    dimer_ids = {str(row["tf_id"]) for row in rows if "dimer" in str(row["model_path"] or "").lower()}
    counts["model_tf_count"] = len(model_ids)
    counts["monomer_model_tf_count"] = len(monomer_ids)
    counts["dimer_model_tf_count"] = len(dimer_ids)
    return counts


def descendants(nodes_by_parent: dict[str | None, list[dict[str, Any]]], family_id: str) -> list[dict[str, Any]]:
    stack = list(nodes_by_parent.get(family_id, []))
    out: list[dict[str, Any]] = []
    while stack:
        node = stack.pop()
        out.append(node)
        stack.extend(nodes_by_parent.get(str(node["family_id"]), []))
    return out


def build_tree(conn: sqlite3.Connection, taxonomy: list[dict[str, Any]]) -> list[dict[str, Any]]:
    nodes_by_parent: dict[str | None, list[dict[str, Any]]] = {}
    for node in taxonomy:
        nodes_by_parent.setdefault(node.get("parent_id"), []).append(node)

    ids_by_family: dict[str, set[str]] = {}
    for node in taxonomy:
        family_id = str(node["family_id"])
        ids_by_family[family_id] = family_tf_ids(conn, list(node.get("tokens") or []))

    for node in sorted(taxonomy, key=lambda item: len(str(item["family_id"])), reverse=True):
        family_id = str(node["family_id"])
        if ids_by_family[family_id]:
            continue
        inherited: set[str] = set()
        for child in descendants(nodes_by_parent, family_id):
            inherited.update(ids_by_family.get(str(child["family_id"]), set()))
        ids_by_family[family_id] = inherited

    all_tf_ids = {str(row["tf_id"]) for row in conn.execute("SELECT tf_id FROM tf").fetchall()}
    if "root" in ids_by_family:
        ids_by_family["root"] = all_tf_ids

    out: list[dict[str, Any]] = []
    for node in taxonomy:
        family_id = str(node["family_id"])
        ids = ids_by_family.get(family_id, set())
        enriched = dict(node)
        enriched["tf_count"] = len(ids)
        enriched.update(count_primary_levels(conn, ids))
        enriched["generated_pwm_tf_count"] = generated_pwm_tf_count(conn, ids)
        enriched.update(model_counts(conn, ids))
        query = str(enriched.get("search_query") or "")
        enriched["open_url"] = "/search" if family_id == "root" else f"/search?q={query}" if query else "/search"
        out.append(enriched)
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description="Build ModCREDB homepage family-tree data from SQLite.")
    parser.add_argument("--db", type=Path, default=DEFAULT_DB)
    parser.add_argument("--taxonomy", type=Path, default=DEFAULT_TAXONOMY)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()

    if not args.db.exists():
        raise SystemExit(f"Database not found: {args.db}")
    if not args.taxonomy.exists():
        raise SystemExit(f"Taxonomy JSON not found: {args.taxonomy}")

    taxonomy = json.loads(args.taxonomy.read_text(encoding="utf-8"))
    with connect(args.db) as conn:
        tree = build_tree(conn, taxonomy)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({"nodes": tree}, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    populated = sum(1 for node in tree if int(node.get("tf_count") or 0) > 0)
    print(f"Wrote {args.out}")
    print(f"Family nodes: {len(tree)}; populated nodes: {populated}")


if __name__ == "__main__":
    main()
