#!/usr/bin/env python3
"""Local patch helper for build_family_tree_data.py.

Replaces the model_counts function so the builder does not assume the model path
column is named path.
"""

from pathlib import Path

p = Path(__file__).resolve().parents[1] / "scripts" / "build_family_tree_data.py"
s = p.read_text()
start = s.index("def model_counts(")
end = s.index("\ndef descendants(", start)
new = '''def model_counts(conn: sqlite3.Connection, ids: set[str]) -> dict[str, int]:
    counts = {"model_tf_count": 0, "monomer_model_tf_count": 0, "dimer_model_tf_count": 0}
    if not ids or not table_exists(conn, "structure_file"):
        return counts
    columns = {str(row["name"]) for row in conn.execute("PRAGMA table_info(structure_file)").fetchall()}
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

'''
p.write_text(s[:start] + new + s[end:])
print(f"Patched {p}")
