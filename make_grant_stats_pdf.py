#!/usr/bin/env python3
"""Generate a small grant-facing PDF summary from the TF web database.

This uses only the Python standard library so it works on Masada without
installing reportlab/pandoc/weasyprint.
"""

from __future__ import annotations

import datetime as dt
import os
import sqlite3
import textwrap
from pathlib import Path


APP_DIR = Path(__file__).resolve().parent
DB_PATH = APP_DIR / "data" / "tf_webdb.sqlite"
OUT_DIR = Path(os.environ.get("TF_WEBDB_REPORT_DIR", APP_DIR / "reports"))
PDF_PATH = OUT_DIR / "TF_database_grant_stats_summary.pdf"
MD_PATH = OUT_DIR / "TF_database_grant_stats_summary.md"

SOURCE_LABELS = {
    "jaspar": "JASPAR",
    "cisbp": "CIS-BP",
    "hocomoco": "HOCOMOCO",
    "modcre": "ModCRE",
    "alphafold": "AlphaFold/AF3",
}

EVIDENCE_LABELS = {
    "identical": "Identical PWM",
    "homologous": "Homologous PWM",
    "relative_homologous": "Relatively homologous PWM",
    "modcre": "ModCRE predicted PWM",
    "alphafold": "AlphaFold/AF3 predicted PWM",
}


def connect() -> sqlite3.Connection:
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    return conn


def scalar(conn: sqlite3.Connection, query: str) -> int:
    return int(conn.execute(query).fetchone()[0])


def rows(conn: sqlite3.Connection, query: str) -> list[dict[str, object]]:
    return [dict(row) for row in conn.execute(query).fetchall()]


def pct(value: int, total: int) -> str:
    if total == 0:
        return "0.0%"
    percent = (value / total) * 100
    if 0 < percent < 0.1:
        return f"{percent:.2f}%"
    return f"{percent:.1f}%"


def collect_stats() -> dict[str, object]:
    with connect() as conn:
        stats = {
            "tf_count": scalar(conn, "SELECT COUNT(*) FROM tf"),
            "annotated_tf_count": scalar(conn, "SELECT COUNT(*) FROM tf_annotation"),
            "motif_ref_count": scalar(conn, "SELECT COUNT(*) FROM motif_ref"),
            "motif_file_count": scalar(conn, "SELECT COUNT(*) FROM motif_file"),
            "fimo_ready_motif_count": scalar(conn, "SELECT COUNT(*) FROM motif_file WHERE matrix_json IS NOT NULL"),
            "active_pdb_count": scalar(
                conn,
                "SELECT COUNT(*) FROM structure_file WHERE status = 'active' AND file_type = 'pdb'",
            ),
            "active_model_summary_count": scalar(
                conn,
                "SELECT COUNT(*) FROM model_summary WHERE status = 'active'",
            ),
            "active_linked_model_summary_count": scalar(
                conn,
                "SELECT COUNT(*) FROM model_summary WHERE status = 'active' AND matched_structure_id IS NOT NULL",
            ),
            "model_summary_count": scalar(conn, "SELECT COUNT(*) FROM model_summary"),
            "missing_unique_count": scalar(
                conn,
                """
                SELECT COUNT(*) FROM (
                    SELECT DISTINCT source, motif_id
                    FROM motif_ref
                    WHERE missing_local_file = 1
                )
                """,
            ),
            "missing_link_count": scalar(conn, "SELECT COUNT(*) FROM motif_ref WHERE missing_local_file = 1"),
            "failed_file_count": scalar(conn, "SELECT COUNT(*) FROM structure_file WHERE status = 'failed'"),
        }
        stats["evidence_rows"] = rows(
            conn,
            """
            SELECT evidence_type, COUNT(*) AS count,
                   SUM(CASE WHEN missing_local_file = 1 THEN 1 ELSE 0 END) AS missing
            FROM motif_ref
            GROUP BY evidence_type
            ORDER BY count DESC
            """,
        )
        stats["source_rows"] = rows(
            conn,
            """
            SELECT source, COUNT(*) AS count,
                   SUM(CASE WHEN missing_local_file = 1 THEN 1 ELSE 0 END) AS missing
            FROM motif_ref
            GROUP BY source
            ORDER BY count DESC
            """,
        )
        stats["motif_file_rows"] = rows(
            conn,
            """
            SELECT source, COUNT(*) AS count,
                   SUM(CASE WHEN matrix_json IS NOT NULL THEN 1 ELSE 0 END) AS matrix_count
            FROM motif_file
            GROUP BY source
            ORDER BY count DESC
            """,
        )
        stats["active_model_rows"] = rows(
            conn,
            """
            SELECT source, COUNT(*) AS count
            FROM structure_file
            WHERE status = 'active' AND file_type = 'pdb'
            GROUP BY source
            ORDER BY count DESC
            """,
        )
        stats["summary_rows"] = rows(
            conn,
            """
            SELECT status, COUNT(*) AS count,
                   SUM(CASE WHEN matched_structure_id IS NOT NULL THEN 1 ELSE 0 END) AS linked_count
            FROM model_summary
            GROUP BY status
            ORDER BY status
            """,
        )
        stats["top_families"] = rows(
            conn,
            """
            SELECT family, COUNT(*) AS count
            FROM tf_family
            GROUP BY family
            ORDER BY count DESC, family
            LIMIT 8
            """,
        )
        stats["missing_motifs"] = rows(
            conn,
            """
            SELECT source, motif_id, COUNT(*) AS count
            FROM motif_ref
            WHERE missing_local_file = 1
            GROUP BY source, motif_id
            ORDER BY count DESC
            """,
        )
    return stats


def pdf_escape(text: object) -> str:
    value = str(text)
    value = value.replace("\\", "\\\\").replace("(", "\\(").replace(")", "\\)")
    return value.encode("latin-1", errors="replace").decode("latin-1")


class PDF:
    def __init__(self) -> None:
        self.pages: list[list[str]] = []
        self.current: list[str] = []
        self.width = 595
        self.height = 842
        self.margin = 54
        self.y = self.height - self.margin

    def new_page(self) -> None:
        if self.current:
            self.pages.append(self.current)
        self.current = []
        self.y = self.height - self.margin

    def ensure(self, needed: float) -> None:
        if self.y - needed < self.margin:
            self.new_page()

    def text(self, text: str, x: float | None = None, size: int = 10, font: str = "F1", leading: int | None = None) -> None:
        if x is None:
            x = self.margin
        if leading is None:
            leading = size + 5
        safe = pdf_escape(text)
        self.current.append(f"BT /{font} {size} Tf {x:.1f} {self.y:.1f} Td ({safe}) Tj ET")
        self.y -= leading

    def paragraph(self, text: str, size: int = 10, font: str = "F1", width_chars: int = 88, space: int = 8) -> None:
        lines = textwrap.wrap(text, width=width_chars)
        self.ensure(len(lines) * (size + 4) + space)
        for line in lines:
            self.text(line, size=size, font=font, leading=size + 4)
        self.y -= space

    def heading(self, text: str, level: int = 1) -> None:
        size = 18 if level == 1 else 13
        space = 12 if level == 1 else 8
        self.ensure(size + 24)
        self.text(text, size=size, font="F2", leading=size + 8)
        self.y -= space

    def rule(self) -> None:
        self.current.append(f"{self.margin:.1f} {self.y:.1f} m {self.width - self.margin:.1f} {self.y:.1f} l S")
        self.y -= 16

    def table(self, headers: list[str], rows_: list[list[object]], widths: list[float], size: int = 9) -> None:
        row_h = size + 9
        self.ensure((len(rows_) + 2) * row_h)
        x_positions = [self.margin]
        for w in widths[:-1]:
            x_positions.append(x_positions[-1] + w)
        self.current.append(f"0.93 0.95 0.98 rg {self.margin:.1f} {self.y - 4:.1f} {sum(widths):.1f} {row_h:.1f} re f 0 g")
        header_y = self.y
        for x, header in zip(x_positions, headers):
            self.current.append(f"BT /F2 {size} Tf {x + 4:.1f} {header_y:.1f} Td ({pdf_escape(header)}) Tj ET")
        self.y -= row_h
        for row in rows_:
            for x, value in zip(x_positions, row):
                text = str(value)
                max_chars = max(8, int((widths[x_positions.index(x)] - 8) / (size * 0.5)))
                if len(text) > max_chars:
                    text = text[: max_chars - 1] + "."
                self.current.append(f"BT /F1 {size} Tf {x + 4:.1f} {self.y:.1f} Td ({pdf_escape(text)}) Tj ET")
            self.y -= row_h
        self.y -= 8

    def bar_table(self, title: str, items: list[tuple[str, int]], total: int | None = None) -> None:
        self.heading(title, 2)
        if not items:
            return
        max_value = max(value for _, value in items) or 1
        label_x = self.margin
        value_x = self.margin + 230
        bar_x = self.margin + 300
        bar_w = 190
        row_h = 18
        self.ensure(len(items) * row_h + 20)
        for label, value in items:
            self.current.append(f"BT /F1 9 Tf {label_x:.1f} {self.y:.1f} Td ({pdf_escape(label[:42])}) Tj ET")
            value_text = f"{value:,}"
            if total:
                value_text += f" ({pct(value, total)})"
            self.current.append(f"BT /F1 9 Tf {value_x:.1f} {self.y:.1f} Td ({pdf_escape(value_text)}) Tj ET")
            filled = (value / max_value) * bar_w
            self.current.append(f"0.93 0.95 0.98 rg {bar_x:.1f} {self.y - 3:.1f} {bar_w:.1f} 7 re f")
            self.current.append(f"0.14 0.39 0.67 rg {bar_x:.1f} {self.y - 3:.1f} {filled:.1f} 7 re f 0 g")
            self.y -= row_h
        self.y -= 8

    def save(self, path: Path) -> None:
        if self.current:
            self.pages.append(self.current)
        objects: list[bytes] = []
        objects.append(b"<< /Type /Catalog /Pages 2 0 R >>")
        page_refs = " ".join(f"{3 + i * 2} 0 R" for i in range(len(self.pages)))
        objects.append(f"<< /Type /Pages /Kids [{page_refs}] /Count {len(self.pages)} >>".encode("latin-1"))
        font_object_id = 3 + len(self.pages) * 2
        for index, commands in enumerate(self.pages):
            content_id = 4 + index * 2
            page = (
                f"<< /Type /Page /Parent 2 0 R /MediaBox [0 0 {self.width} {self.height}] "
                f"/Resources << /Font << /F1 {font_object_id} 0 R /F2 {font_object_id + 1} 0 R /F3 {font_object_id + 2} 0 R >> >> "
                f"/Contents {content_id} 0 R >>"
            )
            stream = "\n".join(commands).encode("latin-1", errors="replace")
            content = b"<< /Length " + str(len(stream)).encode("ascii") + b" >>\nstream\n" + stream + b"\nendstream"
            objects.append(page.encode("latin-1"))
            objects.append(content)
        objects.append(b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>")
        objects.append(b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica-Bold >>")
        objects.append(b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica-Oblique >>")

        output = bytearray(b"%PDF-1.4\n")
        offsets = [0]
        for number, obj in enumerate(objects, start=1):
            offsets.append(len(output))
            output.extend(f"{number} 0 obj\n".encode("ascii"))
            output.extend(obj)
            output.extend(b"\nendobj\n")
        xref_pos = len(output)
        output.extend(f"xref\n0 {len(objects) + 1}\n".encode("ascii"))
        output.extend(b"0000000000 65535 f \n")
        for offset in offsets[1:]:
            output.extend(f"{offset:010d} 00000 n \n".encode("ascii"))
        output.extend(
            f"trailer\n<< /Size {len(objects) + 1} /Root 1 0 R >>\nstartxref\n{xref_pos}\n%%EOF\n".encode("ascii")
        )
        path.write_bytes(output)


def build_markdown(stats: dict[str, object]) -> str:
    evidence_rows = stats["evidence_rows"]
    motif_file_rows = stats["motif_file_rows"]
    active_model_rows = stats["active_model_rows"]
    missing_motifs = stats["missing_motifs"]
    lines = [
        "# TF Motif and Model Database - Grant Statistics Summary",
        "",
        f"Generated: {dt.date.today().isoformat()}",
        "",
        "Prepared from Baldo-provided TF/PWM/model files. ModCRE is a structure-based TF-DNA modeling approach developed by Baldo Oliva.",
        "",
        "## Key numbers",
        "",
        f"- TF entries: {stats['tf_count']:,}",
        f"- UniProt-annotated TFs: {stats['annotated_tf_count']:,} ({pct(stats['annotated_tf_count'], stats['tf_count'])})",
        f"- TF-motif evidence links: {stats['motif_ref_count']:,}",
        f"- Local motif files: {stats['motif_file_count']:,}",
        f"- FIMO-ready parsed motif matrices: {stats['fimo_ready_motif_count']:,}",
        f"- Active PDB models: {stats['active_pdb_count']:,}",
        f"- Active model summary rows: {stats['active_model_summary_count']:,}",
        "",
        "## Motif evidence links",
        "",
    ]
    for row in evidence_rows:
        lines.append(
            f"- {EVIDENCE_LABELS.get(row['evidence_type'], row['evidence_type'])}: {row['count']:,}"
        )
    lines.extend(["", "## Local motif files", ""])
    for row in motif_file_rows:
        label = SOURCE_LABELS.get(row["source"], row["source"])
        lines.append(f"- {label}: {row['count']:,} files, {row['matrix_count']:,} parsed matrices")
    lines.extend(["", "## Active model files", ""])
    for row in active_model_rows:
        label = SOURCE_LABELS.get(row["source"], row["source"])
        lines.append(f"- {label}: {row['count']:,} active PDB models")
    lines.extend(
        [
            "",
            "## Curation note",
            "",
            f"- Missing unique motif IDs: {stats['missing_unique_count']:,}",
            f"- Affected motif links: {stats['missing_link_count']:,} of {stats['motif_ref_count']:,} ({pct(stats['missing_link_count'], stats['motif_ref_count'])})",
        ]
    )
    for row in missing_motifs:
        lines.append(f"  - {SOURCE_LABELS.get(row['source'], row['source'])} {row['motif_id']}: {row['count']} links")
    return "\n".join(lines) + "\n"


def build_pdf(stats: dict[str, object]) -> PDF:
    pdf = PDF()
    pdf.heading("TF Motif and Model Database", 1)
    pdf.text("Grant Statistics Summary", size=14, font="F2", leading=20)
    pdf.text(f"Generated: {dt.date.today().isoformat()}", size=9, font="F3", leading=16)
    pdf.rule()
    pdf.paragraph(
        "Prepared from Baldo-provided transcription factor PWM and model files. "
        "ModCRE is a structure-based TF-DNA modeling approach developed by Baldo Oliva. "
        "This summary focuses on grant-relevant database coverage and omits internal debugging detail."
    )

    key_rows = [
        ["TF entries", f"{stats['tf_count']:,}"],
        [
            "UniProt-annotated TFs",
            f"{stats['annotated_tf_count']:,} ({pct(stats['annotated_tf_count'], stats['tf_count'])})",
        ],
        ["TF-motif evidence links", f"{stats['motif_ref_count']:,}"],
        ["Local motif files", f"{stats['motif_file_count']:,}"],
        ["FIMO-ready motif matrices", f"{stats['fimo_ready_motif_count']:,}"],
        ["Active PDB models", f"{stats['active_pdb_count']:,}"],
        ["Active model summary rows", f"{stats['active_model_summary_count']:,}"],
        [
            "Summary rows linked to models",
            f"{stats['active_linked_model_summary_count']:,} "
            f"({pct(stats['active_linked_model_summary_count'], stats['active_model_summary_count'])})",
        ],
    ]
    pdf.heading("Core Database Coverage", 2)
    pdf.table(["Metric", "Value"], key_rows, [310, 160], size=10)

    evidence_items = [
        (EVIDENCE_LABELS.get(row["evidence_type"], row["evidence_type"]), int(row["count"]))
        for row in stats["evidence_rows"]
    ]
    pdf.bar_table("TF-Motif Evidence Links", evidence_items, total=int(stats["motif_ref_count"]))

    motif_file_rows = []
    for row in stats["motif_file_rows"]:
        motif_file_rows.append(
            [
                SOURCE_LABELS.get(row["source"], row["source"]),
                f"{int(row['count']):,}",
                f"{int(row['matrix_count'] or 0):,}",
            ]
        )
    pdf.heading("Local Motif File Coverage", 2)
    pdf.table(["Source", "Files", "Parsed matrices"], motif_file_rows, [190, 120, 160], size=10)

    active_model_items = [
        (SOURCE_LABELS.get(row["source"], row["source"]), int(row["count"]))
        for row in stats["active_model_rows"]
    ]
    pdf.bar_table("Active TF-DNA Model Coverage", active_model_items, total=int(stats["active_pdb_count"]))

    pdf.new_page()
    top_families = [(row["family"], int(row["count"])) for row in stats["top_families"]]
    pdf.bar_table("Top TF Families / Domains", top_families)

    pdf.heading("Curation Note", 2)
    pdf.paragraph(
        f"The imported chart references {stats['missing_unique_count']:,} unique motif IDs not present in the local motif archives. "
        f"These affect {stats['missing_link_count']:,} of {stats['motif_ref_count']:,} TF-motif links "
        f"({pct(stats['missing_link_count'], stats['motif_ref_count'])}). "
        "The affected links are preserved for traceability and can be resolved by adding the exact MEME file or documenting a version remap."
    )
    missing_rows = []
    for row in stats["missing_motifs"]:
        missing_rows.append([SOURCE_LABELS.get(row["source"], row["source"]), row["motif_id"], row["count"]])
    pdf.table(["Source", "Motif ID", "Links"], missing_rows, [160, 180, 90], size=10)

    pdf.heading("Grant-Relevant Interpretation", 2)
    pdf.paragraph(
        "The resource integrates curated and predicted TF binding specificity data with structural model evidence. "
        "It combines public motif resources (JASPAR and CIS-BP), ModCRE predicted PWMs and TF-DNA models, "
        "AlphaFold/AF3-derived predictions, UniProt annotation, sequence-logo visualization, and FIMO-based DNA scanning."
    )
    pdf.paragraph(
        "For proposal language, the strongest headline numbers are the 5,384 TF entries, 112,693 TF-motif evidence links, "
        "23,638 FIMO-ready motif matrices, and 12,072 active PDB model files."
    )
    return pdf


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    stats = collect_stats()
    MD_PATH.write_text(build_markdown(stats), encoding="utf-8")
    pdf = build_pdf(stats)
    pdf.save(PDF_PATH)
    print(PDF_PATH)
    print(MD_PATH)


if __name__ == "__main__":
    main()
