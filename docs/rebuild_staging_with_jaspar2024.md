# Rebuilding ModCREDB Staging Data

This workflow creates a fresh staging database without changing the active
staging database or the raw archives. It preserves every supplied evidence
record and then adds direct JASPAR 2024 UniProt mappings from the official
JASPAR metadata dump.

## Import Order

1. Import the original Baldo chart and all supplied motif/model archives.
2. Add HOCOMOCO v11 motifs and the hierarchical primary-annotation table.
3. Enrich the complete TF set from UniProt. This happens after HOCOMOCO
   because HOCOMOCO adds 33 TF records absent from the original chart.
4. Add direct JASPAR 2024 `motif_id -> UniProt accession` links. Existing
   JASPAR, CIS-BP, ModCRE, AlphaFold3-assisted, and HOCOMOCO evidence is never
   removed. A direct JASPAR link can coexist with an older homologous link.

## Command

Set these paths for the machine where raw data are stored:

```bash
RAW_DATA_DIR=/path/to/TF_database_Baldo_data
HOCOMOCO_DIR="$RAW_DATA_DIR/external_sources/HOCOMOCO_v11_core_HUMAN_mono"

cd /path/to/tf_webdb

python3 scripts/build_staging_with_jaspar2024.py \
  --data-dir "$RAW_DATA_DIR" \
  --hocomoco-meme "$HOCOMOCO_DIR/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme" \
  --hocomoco-annotation "$HOCOMOCO_DIR/HOCOMOCOv11_core_annotation_HUMAN_mono.tsv" \
  --hierarchical-tsv "$RAW_DATA_DIR/TF_PWM_chart_final_integrated_HOCOMOCO_hierarchical.tsv" \
  --output-db data/tf_webdb_hocomoco_jaspar2024_rebuild.sqlite \
  --audit-dir docs/audits/jaspar2024_rebuild
```

The output database path must not already exist. The script writes all
intermediates into a temporary directory; it only creates the requested output
database after all four stages finish successfully.

## Expected JASPAR Controls

- `MA0106.3 -> P04637` (TP53) is a direct JASPAR 2024 metadata link.
- `MA0065.1 -> P37231` (PPARG) is already present in the supplied chart and
  must remain present.

The JASPAR metadata importer writes direct rows with:

```text
evidence_type     identical
mapping_type      direct_or_identical
original_column   JASPAR2024_metadata
curation_status   pending_confirmation
```

The exact metadata release and checksum are recorded in `source_release` from
`config/source_releases.tsv`. The JASPAR 2024 SQL dump used for this workflow
has SHA256 `33f00db0f06adaa5fff497827a0f63551da4c30e5c6dd601c08ac625cecc60c4`.

## Validation

After a rebuild, compare it with the current staging database:

```bash
python3 - <<'PY'
import sqlite3

for db in (
    "data/tf_webdb_hocomoco_staging.sqlite",
    "data/tf_webdb_hocomoco_jaspar2024_rebuild.sqlite",
):
    con = sqlite3.connect(db)
    print(db, con.execute("PRAGMA integrity_check").fetchone()[0])
    for table in ("tf", "motif_file", "motif_ref", "structure_file", "model_summary"):
        print(table, con.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0])
    print("TP53 MA0106.3", con.execute(
        "SELECT COUNT(*) FROM motif_ref WHERE tf_id='P04637' AND source='jaspar' AND motif_id='MA0106.3'"
    ).fetchone()[0])
    con.close()
PY
```

Do not replace `data/tf_webdb_hocomoco_staging.sqlite` until the comparison and
browser checks pass.
