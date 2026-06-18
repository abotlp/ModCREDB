# TF Motif Web Database

Local first-version web database for transcription factor PWM/model files.

The code is kept separate from the large private data files. Build/import commands expect a local data directory containing the source TSV and archives.

The first version is intentionally small and easy to move later. It uses:

- SQLite for the local database.
- Python standard library HTTP server for local browsing.
- Jinja2 templates for HTML.
- MEME motif files stored inside SQLite for fast motif pages and downloads.
- Structure/model archive paths indexed in SQLite, without unpacking the large `models.tar.gz`.

## Data Files

SQLite databases, raw archives, extracted PDB/model files, and UniProt-enriched
production data are not committed to GitHub. Use `.gitignore` as the source of
truth for excluded local artifacts.

For local import, provide a directory containing:

```text
TF_PWM_chart_final.tsv
jaspar.tar.gz
cisbp.tar.gz
pwms.tar.gz
models.tar.gz
```

HOCOMOCO can be added to a staging database with `import_hocomoco.py` after
downloading the HOCOMOCO v11 CORE human mononucleotide MEME file.

## Build The Database

From this directory:

```bash
python3 import_db.py --data-dir /path/to/source_data
```

That creates:

```text
data/tf_webdb.sqlite
```

For a quicker UI-only test that skips the large model archive index:

```bash
python3 import_db.py --skip-model-index
```

## Add UniProt Annotation

After the first import, enrich TF accessions with gene names, protein names,
organisms, reviewed status, taxonomy IDs, and sequence lengths:

```bash
python3 enrich_uniprot.py
```

For a small test run:

```bash
python3 enrich_uniprot.py --limit 20
```

The enrichment script only fetches missing annotations by default. Use
`--all` when you intentionally want to refresh everything.

## Run Locally

```bash
python3 app.py --host 127.0.0.1 --port 8090 --db data/tf_webdb.sqlite
```

For local curation only, expose the internal debug page and exception details with:

```bash
python3 app.py --host 127.0.0.1 --port 8090 --db data/tf_webdb.sqlite --enable-debug --show-errors
```

Keep `--enable-debug` and `--show-errors` disabled in production.

Then open:

```text
http://127.0.0.1:8090/
```

## FIMO Scanner

The `/scan` page uses MEME/FIMO. Install MEME Suite or load the MEME module,
then make sure the `fimo` command is available in `PATH` before starting the
web app:

```bash
command -v fimo
fimo --version
```

## Current Pages

- `/` home page with database stats.
- `/search` searchable TF list.
- `/tf/<TF_ID>` TF detail page with motif evidence grouped by evidence type.
- `/motif?source=<source>&id=<motif_id>` motif page with PWM visualization, MEME download, and external source links when available.
- `/model-summaries/<TF_ID>` model summary rows parsed from ModCRE `.summary.txt` files.
- `/model?id=<MODEL_FILE_ID>` interactive 3D viewer for indexed PDB model files.
- `/scan` FIMO scanner for pasted DNA/FASTA sequences, with motif search/selection and example input.
- `/docs` documentation landing page with database statistics and import coverage.
- `/evidence` explanation of evidence levels and interpretation.
- `/debug` optional local-only curation page for missing motif files and failed model counts, enabled with `--enable-debug`.

## V1 Design Choices

- Failed models are indexed as `status = failed`, but hidden from normal TF pages.
- The importer trims the TSV header ` AlphaFold` to `AlphaFold`.
- Missing local motifs are kept as motif references and flagged in `import_issue`.
- UniProt annotations are optional. Pages still work before enrichment, but
  search becomes more useful after enrichment because users can search gene,
  protein, and organism names.
- ModCRE summary rows are stored in `model_summary`, with links to PDB/PIR
  model files when the model filename can be reconstructed.
- The app can download active PDB model files from the archive, but this can be slow because `models.tar.gz` is compressed.
- The 3D viewer currently uses NGL from a public CDN and streams PDB files
  from `models.tar.gz`. This is fine for local testing, but the public server
  should use a local NGL copy and unpacked/object-stored PDB files for speed.
- The scanner writes the selected motifs as a temporary MEME file, runs
  MEME/FIMO against the pasted sequence, and reports FIMO p-values/q-values.
  The web app process must be started in an environment where `fimo` is
  available in `PATH`.

## Next Useful Additions

1. Vendor the NGL JavaScript locally instead of loading it from a CDN.
2. Add export formats for scanner results and a queued job mode for longer sequences.
3. Replace the tiny local server with Flask/FastAPI when deploying publicly.
4. Unpack or object-store model files on the university server for fast download and visualization.


## Source-release metadata

Public provenance metadata is stored in `config/source_releases.tsv` and loaded into the `source_release` SQLite table by `import_db.py` or `migrate_source_releases.py`. Unknown release, license, citation, or workflow details are deliberately marked `pending confirmation` until PI/Baldo confirmation.


## Matrix QC metadata

Motif scan/logo readiness is stored on `motif_file` using `matrix_status` and related QC columns. New imports classify these fields automatically. Existing SQLite databases can be updated without rebuilding by running `python3 migrate_matrix_qc.py --db path/to/database.sqlite`. Only `matrix_status = 'usable'` motifs are selected for FIMO scanning.
