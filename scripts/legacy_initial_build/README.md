# Legacy initial database build

These files are exact copies of the original build/enrichment scripts from
the parent of commit `6749dec` (`Clean ModCREDB serving repository`). They are
retained so the first TSV/archive-to-SQLite workflow remains visible in the
current source tree.

They describe the historical initial build; they do **not** incorporate the
later schema migrations, structural adjudications, complete-model import,
direct-evidence reconciliation, or UniProt corrections. Do not run them
against the final staging or production database. Use the current scripts and
migrations in the repository root for those controlled later stages.

Restored historical files:

- `import_db.py`
- `enrich_uniprot.py`
- `import_hocomoco.py`
- `index_model_summaries.py`
- `migrate_evidence_semantics.py`
- `migrate_matrix_qc.py`
- `migrate_source_releases.py`
