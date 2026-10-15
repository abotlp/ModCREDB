# ModCREDB technical guide

## Purpose and biological model

ModCREDB connects human transcription-factor accessions to PWM evidence and to
protein–DNA structural evidence. These are related but independent axes:

- **PWM provenance** records whether a motif is direct/identical, transferred
  from a close or more distant homolog, predicted by ModCRE, or associated with
  an AlphaFold-assisted workflow.
- **Structural provenance** records the exact coordinate model, its source, the
  modeled protein interval, DNA content, and protein–DNA interface QC.

A direct PWM does not make a structural model valid, and an atom contact alone
does not prove that a modeled fragment is the biologically expected DNA-binding
domain.

## Immutable source data

The delivered source data live outside this repository. The original TSVs
must not be edited in place:

```text
/data/sbi/interchange/boliva/patricia/TF_PWM_chart_final.tsv
/home/patricia/TF_database_Baldo_data/TF_PWM_chart_final_integrated_HOCOMOCO_hierarchical.tsv
/data/sbi/interchange/boliva/patricia/models
/data/sbi/interchange/boliva/patricia/models_failed_in_DNA
/data/sbi/interchange/boliva/patricia/AF3
/data/sbi/interchange/boliva/patricia/AF3_unresolved
```

The repository stores scripts, migrations, small verified mapping manifests,
and selected audit truth tables. SQLite databases, coordinate archives, and
extracted caches remain local because of size and data-provenance constraints.

## Current Masada staging release

The local release database is:

```text
/home/patricia/tf_webdb/data/tf_webdb_local_staging.sqlite
```

The validated 2026-10-15 snapshot contains:

- 5,417 TF accessions and 5,417 primary-annotation rows;
- 117,372 motif evidence links and 26,330 local motif matrices;
- 62,505 active archive-backed PDB records: 58,976 ModCRE and 3,529
  AlphaFold-source PDBs;
- 114 active final AF3 assignments stored separately in
  `structure_model_assignment`;
- 108 curated accessions explicitly marked
  `PWM_PRESENT_NO_VALIDATED_DNA_BOUND_STRUCTURE`.

The primary classes in this release are 2,296 `Identical_PWM`, 1,730
`Homologous_PWM`, 631 `Relatively_Homologous_PWM`, 390 `ModCRE`, 241
`AlphaFold`, and 129 `Unannotated`. These counts differ from the immutable
delivered TSV because later direct JASPAR/CisBP evidence was verified and
reconciled. The original classifications remain recoverable from the source
TSVs and audit provenance.

## Model storage and serving

`structure_file` stores archive-backed PDB records. Its `archive_path` and
`member_path` identify the authoritative recoverable coordinate. For normal
web serving, `app.py` first checks:

```text
TF_WEBDB_MODEL_CACHE_DIR/<member_path>
```

and falls back to the archive when the cache is unavailable. Build or verify
the cache with:

```bash
python3 scripts/build_structure_model_cache.py --dry-run \
  --database data/tf_webdb_local_staging.sqlite \
  --cache-dir /home/patricia/TF_database_Baldo_data/modcredb_model_cache \
  --qc-output /tmp/modcredb_model_cache_qc.json
```

After a controlled database import adds records, add only missing indexed
members to an otherwise valid cache with `--apply --sync-existing`. This mode
refuses extra or wrong-size files and is idempotent.

The selected contemporary AF3 mmCIF files are represented in
`structure_model_assignment`; their mmCIF, pLDDT JSON, and PAE JSON metadata
are recorded in `structure_confidence_artifact`. They are not conflated with
legacy archive-backed AlphaFold PDB rows.

Public archive-backed model links use `/model?model=<model_id>`. SQLite numeric
row IDs may change during rebuilds and are therefore rejected as obsolete
public links. AF3 assignment routes remain `/model?assignment=<id>` because
they resolve through the separately curated assignment table.

## Structural QC

Migration 003 adds model type, protein/DNA chain counts, contact count,
interaction status, cutoff, and audit status to `structure_file`.
`scripts/annotate_structure_dna_interaction.py` records contacts between
non-hydrogen protein and nucleic-acid atoms at a 4.5 Å cutoff. The UI reports:

- **Yes**: protein and DNA are present and at least one contact exists;
- **No**: protein and DNA are present but no contact was found;
- **Unavailable**: no protein component exists in the coordinate file;
- **Not assessed**: no completed QC result is stored.

The complete ModCRE build preserves distinct model outputs and orders verified
DNA-contacting two-DNA-chain models first, with DIMER before TFS where other
evidence is equal. DIMER/TFS describes protein arrangement; it is not a DNA
chain count.

## PWM evidence and scan behavior

The public labels map to the stored evidence hierarchy:

- `identical` → Known;
- `homologous` → Nearest Neighbor (>70%);
- `relative_homologous` → Nearest Neighbor (50–70%);
- `modcre` and `alphafold` → predicted evidence.

TF-level scanning selects all usable matrices assigned to the exact TF across
those five evidence classes. A selected protein region further restricts the
set to motifs explicitly linked to an active structure overlapping that
interval. It is not a genomic-region selector. The global scanner has a
separate 5,000-motif default ceiling. FIMO must be installed and visible in
the web process `PATH`.

## UniProt identity handling

Deleted historical accessions remain in the database because their PWM/model
evidence is part of the delivered scientific record. They are labelled as
deleted rather than silently discarded. Verified merged records store their
replacement in `tf_annotation.merged_into_accession` (migration 004), and the
TF page links to the current UniProt entry. Verified current-name corrections
are imported by `scripts/import_verified_uniprot_updates.py` from the dated
manifest in `data_sources/`.

## Reproducibility map

The exact original TSV/archive import and early enrichment utilities removed
from the serving checkout in commit `6749dec` are preserved under
`scripts/legacy_initial_build/`. They document the initial database build but
must not be run against the final schema.

The main controlled stages are:

1. audit the two ModCRE batches with
   `audit_complete_modcre_interfaces.py`;
2. construct the complete de-duplicated model manifest with
   `build_complete_modcre_model_manifest.py`;
3. import/package it with
   `import_complete_modcre_models_to_staging.py`;
4. import the final AF3 set with
   `import_local_af3_models_to_staging.py`;
5. import the final 108 status records with
   `import_baldo_final_108_structure_status.py`;
6. restore the complete 5,417-accession delivered universe with
   `restore_manuscript_5417_to_staging.py`;
7. reconcile verified direct motif evidence with
   `reconcile_tf_primary_annotations.py`;
8. apply verified UniProt updates with
   `import_verified_uniprot_updates.py`;
9. remove only proven stale structural foreign-key rows with
   `cleanup_stale_structure_foreign_keys.py`;
10. build/validate the extracted runtime cache with
    `build_structure_model_cache.py`.

Intermediate audit scripts are retained because they document how the frozen
truth tables were produced. Scripts with write modes must always be run on a
staging copy first; generated previews and QC reports should be reviewed before
deployment.

## Release validation checklist

Before packaging or deployment:

1. verify `PRAGMA quick_check` and `PRAGMA foreign_key_check`;
2. record the staging database SHA256;
3. validate exact TF, motif, structure, AF3 assignment, and final-108 counts;
4. confirm the model cache has one readable file per active PDB record;
5. run Python compilation and `git diff --check`;
6. test a Known-motif TF, a model-rich TF, a ModCRE PDB, an AlphaFold PDB, a
   contemporary AF3 mmCIF, a final-108 banner case, and a normal model-less TF;
7. test selected-motif, TF-wide, and structure-region FIMO scans;
8. create a checksummed release bundle and deploy it only after a verified
   production backup.
