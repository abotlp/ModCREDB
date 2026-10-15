# ModCREDB

ModCREDB is a human transcription-factor database combining PWM evidence,
protein–DNA structural models, domain/family annotations, and a FIMO DNA
scanner. The web application is a small Python/Jinja2 service backed by
SQLite.

The complete implementation and release workflow are documented in
[`docs/MODCREDB_TECHNICAL_GUIDE.md`](docs/MODCREDB_TECHNICAL_GUIDE.md).

## Local release database

The publication-final working database is intentionally not stored in Git:

```text
data/tf_webdb_local_staging.sqlite
```

Its corresponding model cache is also local:

```text
/home/patricia/TF_database_Baldo_data/modcredb_model_cache
```

Set `TF_WEBDB_MODEL_CACHE_DIR` when starting the application so active PDBs
are served from the extracted cache rather than repeatedly decompressed from
their archives.

```bash
cd /home/patricia/tf_webdb
export TF_WEBDB_MODEL_CACHE_DIR=/home/patricia/TF_database_Baldo_data/modcredb_model_cache
python3 app.py --host 127.0.0.1 --port 8094 \
  --db data/tf_webdb_local_staging.sqlite
```

Open `http://127.0.0.1:8094/` locally, or use an SSH tunnel from another
machine.

## Runtime dependencies

```bash
python3 -m pip install -r requirements.txt
```

The `/scan` route additionally requires MEME Suite/FIMO in `PATH`. Structural
audit scripts that parse mmCIF use the separately documented audit dependency:

```bash
python3 -m pip install -r requirements-audit.txt
```

## Important data rules

- Original Baldo/manuscript TSVs and model archives are immutable source data.
- Primary PWM annotation follows the evidence hierarchy documented in the
  technical guide; later verified direct evidence may legitimately promote a
  TF without deleting its supporting predicted motifs/models.
- A model is not considered DNA-interacting merely because protein and DNA are
  present. The stored interface QC uses explicit protein–DNA atom contacts.
- The 108 curated accessions with PWM evidence but no validated DNA-bound
  structure are represented explicitly; the label is never inferred merely
  from a zero model count.
- Public PDB links use the active stable `structure_file.model_id`:
  `/model?model=<model_id>`. Numeric row IDs are internal and deliberately not
  accepted as public model identifiers after a database rebuild.

## Safety

Database-changing scripts default to dry-run where implemented and validate
their frozen inputs. Read each script's `--help`, create a verified backup,
and validate on a staging copy before applying it. Production databases,
large model archives, extracted caches, backups, and generated reports are
excluded from Git.
