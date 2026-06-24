# ModCREDB Technical Guide

This is the practical reference for the current ModCREDB prototype: what it stores, how its data are imported, what users can search, and what the results mean.

## 1. What ModCREDB Is

ModCREDB is a searchable database of **human transcription-factor UniProt sequence records**, their DNA-binding motif evidence, and available structural-model evidence.

It is an evidence database. It indexes precomputed information from public motif resources and from the supplied ModCRE/ModCRElib and AlphaFold3-assisted prediction datasets. It does **not** run a new ModCRElib or AlphaFold3 modelling job when a user opens a page.

The central scientific rule is:

```text
One TF sequence record has one primary annotation label.
All motif and model evidence links remain stored and visible.
```

Therefore, the primary annotation is a compact summary, not a replacement for the complete evidence record.

## 2. What A User Can Do

| User task | Example input | Result |
| --- | --- | --- |
| Find a TF record | `P04637`, `TP53`, `CTCF`, `C2H2` | One or more UniProt sequence records, their annotations, motifs, and models. |
| Find a motif | `MA0106.3`, `PPARG_HUMAN.H11MO.0.A` | Motif matrix, logo, linked TF records, provenance, and model links where available. |
| Browse a TF | `/tf/P04637` | Primary annotation, FIMO-ready motifs, all preserved evidence links, and active structural models. |
| Scan DNA | DNA sequence plus motif(s) or TF accession | FIMO motif hits with scores, p-values, q-values, and a positional hit summary. |
| Inspect a model | A model linked from a TF or motif page | Model metadata, template, residue interval, summary information, PDB download, and 3D viewer where available. |

Search results are **UniProt sequence records**, not one row per gene. A gene symbol can legitimately return several reviewed or unreviewed UniProt records.

## 3. Main Pages And Their Meaning

| Route | Purpose |
| --- | --- |
| `/` | Entry page with search and examples. |
| `/search?q=...` | Searches accessions, gene symbols, motif IDs, protein names, and TF-family text. |
| `/tf/<accession>` | One TF sequence record and all retained motif/model evidence. |
| `/motif?source=<source>&id=<motif_id>` | One motif record and its source, matrix QC, TF links, and model links. |
| `/scan` | FIMO DNA-sequence scan. |
| `/model?id=<id>` | Structural-model detail and 3D view. |
| `/model-summaries/<accession>` | Model-summary rows for a TF record. |
| `/evidence` | Explanation of evidence layers. |
| `/docs` | Short public documentation and statistics. |

## 4. Evidence Vocabulary

Several concepts must stay separate. They answer different questions.

| Concept | Meaning |
| --- | --- |
| `source` | Where a motif or model originated, for example JASPAR, CIS-BP, HOCOMOCO, ModCRE, or AlphaFold3-assisted ModCRE. |
| `evidence_type` | Why a TF record is linked to a motif: `identical`, `homologous`, `relative_homologous`, `modcre`, or `alphafold`. |
| `mapping_type` | How the link was made, for example `direct_or_identical`, `close_homolog`, `distant_homolog`, `structure_predicted`, or `af3_structure_predicted`. |
| primary annotation | One final hierarchical label for a TF record. It is stored separately from all links. |
| `matrix_status` | Whether a motif contains a usable DNA probability matrix for logo display and FIMO. |
| model status | Whether a structure file is active/public evidence or failed/debug material. |

The public display label for raw database value `AlphaFold` is **AlphaFold3-assisted ModCRE**. The stored database value is retained unchanged for reproducibility.

### Evidence Layers

The current evidence hierarchy is presented to users as:

1. Direct PWM: a direct motif assignment from JASPAR, CIS-BP, or HOCOMOCO.
2. Close homologous PWM: transferred from a close homolog.
3. Distant homologous candidate: transferred from a more distant homolog.
4. ModCRE predicted motif: predicted from a precomputed protein-DNA structural model.
5. AlphaFold3-assisted ModCRE: predicted using an AlphaFold3-assisted structural model.
6. Unannotated: no primary motif/model was selected in the supplied hierarchy.

The hierarchy is a useful single-label summary. It must never be interpreted as deleting lower-level or alternative evidence stored in the database.

## 5. Input Files

Raw archives and original supplied files are deliberately kept outside the Git repository. The app database is built from controlled copies of the following inputs.

| Input | Content | Role in ModCREDB |
| --- | --- | --- |
| `TF_PWM_chart_final.tsv` | Original supplied annotation chart with seven evidence columns. | Base TF records, families, and original motif/model references. |
| `TF_PWM_chart_final_integrated_HOCOMOCO_hierarchical.tsv` | Final supplied chart with the original evidence columns plus hierarchy fields. | Loads `tf_primary_annotation`; it does not erase the original evidence links. |
| `jaspar.tar.gz` | JASPAR MEME motif files. | Creates local JASPAR motif records and matrices. |
| `cisbp.tar.gz` | CIS-BP MEME motif files. | Creates local CIS-BP motif records and matrices. |
| `pwms.tar.gz` | ModCRE and AlphaFold3-assisted MEME/PWM files. | Creates predicted motif records, including no-matrix files retained for provenance. |
| `models.tar.gz` | Structural PDB files, summaries, and failed/debug records. | Creates active model, summary, and debug-file indexes. |
| HOCOMOCO v11 human CORE mononucleotide MEME and annotation files | Public HOCOMOCO release files. | Creates HOCOMOCO motif records and exact motif IDs. |
| JASPAR 2024 metadata database | Official JASPAR metadata mapping source. | Reproducibly proposes or applies exact direct JASPAR motif-to-UniProt links. |

The final hierarchical chart has these columns:

```text
TF_name
TF_family
Identical_PWM
Homologous_PWM
Relatively_Homologous_PWM
ModCRE
AlphaFold
Best_annotation_level
Best_PWM_or_model
N_nonempty_annotation_columns
```

The first seven columns preserve the evidence supplied in the original chart. The final three columns describe the selected primary annotation summary.

## 6. Database Files

The active app database is a local SQLite file selected when starting the server. The current development/staging database is normally named:

```text
data/tf_webdb_hocomoco_staging.sqlite
```

SQLite databases and large raw archives are intentionally ignored by Git. They should not be uploaded as part of a public source-code repository without a separate data-sharing decision.

The full SQL schema snapshot for review is stored in:

```text
docs/review/modcredb_schema.sql
```

### Core Tables

| Table | What one row represents |
| --- | --- |
| `tf` | One UniProt TF sequence record. |
| `tf_annotation` | UniProt-derived annotation fields such as gene, protein name, organism, review status, and sequence length. |
| `tf_family` | A TF-family label attached to a TF record. |
| `tf_primary_annotation` | One final hierarchical annotation summary for a TF record. |
| `motif_file` | One locally available motif file and its parsed matrix/QC state. |
| `motif_ref` | One preserved TF-to-motif evidence link. This is many-to-many. |
| `structure_file` | One indexed model file, including active and failed/debug files. |
| `model_summary` | One parsed structure-summary row. |
| `motif_structure` | An exact motif-to-structure relationship when the supplied files establish one. |
| `source_release` | Source/release/provenance metadata for JASPAR, CIS-BP, HOCOMOCO, ModCRE, AlphaFold3-assisted ModCRE, and UniProt. |
| `import_issue` | References that could not be resolved during import, such as missing local motif files. |
| `metadata` | Import-time metadata. It is operational metadata and may contain non-public local paths, so it must be redacted before sharing a database. |

### Important Relationships

```text
tf --< motif_ref >-- motif_file
tf --< structure_file
motif_file --< motif_structure >-- structure_file
tf --1 tf_primary_annotation
tf --< tf_annotation
tf --< tf_family
```

One TF can have many motifs because it may have several external annotations, homologous candidates, DNA-binding regions, templates, conformations, or model sources. One motif can link to multiple TF records.

## 7. Matrix Quality Control And FIMO

`motif_file.matrix_status` is the authoritative public classification for whether a current DNA motif can be rendered as a logo or supplied to FIMO.

| Matrix status | Meaning | Logo/FIMO |
| --- | --- | --- |
| `usable` | A parsable A/C/G/T probability matrix has a consistent width and usable rows. | Yes |
| `width_zero_no_matrix` | The MEME file explicitly reports `w=0`; no usable probability matrix was produced. | No |
| `no_parsed_matrix` | A file exists but no usable rows were parsed. | No |
| `malformed_matrix` | Matrix rows are invalid, non-numeric, or have the wrong shape. | No |
| `width_mismatch` | Declared matrix width and parsed row count disagree. | No |
| `unsupported_alphabet` | The file is not usable by the current DNA scanner. | No |
| `missing_local_file` | The evidence link exists but the local motif file is absent. | No |
| `unknown` | Legacy or not-yet-classified import state. | Not assumed usable |

No-matrix records are not necessarily failed TFs. A TF may have a structural model and an empty motif matrix for one model, while still having another usable motif from the same or a different source.

The scan page runs **MEME/FIMO** only with distinct motif records whose `matrix_status` is `usable`. It reports FIMO hits, scores, p-values, q-values, and a positional summary. This is a motif scan, not a full ModCRE structural binding profile or a new structural simulation.

## 8. Structures And Models

`structure_file` indexes the supplied model data. The normal public pages use active PDB model files. Failed/debug files are retained for quality-control provenance but are not presented as public model evidence by default.

`model_summary` provides available information such as a template PDB/chain, residue interval, alignment identity, coverage, and model metadata. A summary row is useful supporting information; it is not automatically proof that every expected chain or interface relationship is biologically validated.

Some TFs have multiple DNA-binding regions. In those cases, motifs and models should remain separated by their residue intervals and model IDs. ModCREDB must not merge them into one artificial PWM.

## 9. Import And Rebuild Workflow

The ordinary reproducible staging rebuild is documented in:

```text
docs/rebuild_staging_with_jaspar2024.md
```

The wrapper script is:

```text
scripts/build_staging_with_jaspar2024.py
```

Its intended order is:

1. Import base chart, motif archives, and model archive with `import_db.py`.
2. Import HOCOMOCO motifs and primary hierarchy with `import_hocomoco.py`.
3. Add UniProt annotations with `enrich_uniprot.py`.
4. Apply official JASPAR 2024 exact motif-to-UniProt metadata links with `scripts/import_jaspar_metadata_links.py`.
5. Validate the new staging database before promoting it.

Always build a new SQLite file first. Compare counts, test representative pages, and only then decide whether to replace the currently served database.

### What Not To Do

- Do not edit individual SQLite rows by hand to fix a biological mapping.
- Do not infer direct JASPAR links from motif names alone.
- Do not overwrite the existing staging database before validating a newly rebuilt one.
- Do not treat all no-matrix records as TF failures.
- Do not collapse all evidence links to a single best record.

## 10. Running The App Locally

The prototype uses Python's standard-library HTTP server, Jinja2 templates, and SQLite. It is not a Flask or Django application.

FIMO must be available on `PATH` for scan functionality. On the local cluster environment it has been used through a MEME module.

```bash
module load MEME/5.5.7
cd /path/to/tf_webdb
python3 app.py --host 127.0.0.1 --port 8198 --db data/tf_webdb_hocomoco_staging.sqlite
```

The 3D viewer uses NGL JavaScript. For public deployment, decide whether the existing external CDN dependency is acceptable or whether a vetted local copy should be hosted with the application.

## 11. Routine Validation

Before a rebuild, deployment, or code review, run the relevant checks:

```bash
python3 -m py_compile app.py import_db.py import_hocomoco.py \
  enrich_uniprot.py scripts/build_staging_with_jaspar2024.py \
  scripts/import_jaspar_metadata_links.py

sqlite3 data/tf_webdb_hocomoco_staging.sqlite 'PRAGMA integrity_check;'
git diff --check
```

Useful representative records include:

| Example | Why it is useful |
| --- | --- |
| `P04637` / TP53 | Direct JASPAR metadata-link control with `MA0106.3`. |
| `P37231` / PPARG | Direct multi-source control including JASPAR `MA0065.1` and HOCOMOCO. |
| `P49711` / CTCF | Direct motif and ModCRE model comparison control. |
| `A0A087WX29` | Active models with no FIMO-ready PWM; the TF scan action must be unavailable. |
| `A0A0S2Z4K5` | Relatively homologous evidence and a long candidate list. |

## 12. Known Limits And Questions Still Open

- The exact ModCRE/ModCRElib version, script settings, template library, and prediction parameters used for the supplied precomputed dataset need PI/Baldo confirmation.
- The exact AlphaFold3-assisted generation workflow and causes of some no-matrix outputs remain to be confirmed.
- Failed/debug model IDs can overlap active model IDs. They are audit material, not automatically evidence of a biological failure.
- HOCOMOCO v11 source release is known, but the exact mapping semantics used by the supplied hierarchical chart should remain labelled as pending confirmation until confirmed by the data provider.
- CIS-BP does not yet have the same release-pinned direct UniProt metadata import workflow as JASPAR.
- The local database includes operational import metadata that may contain local paths. Share only a redacted sample database or a cleaned deployment database.

## 13. Further Documentation

| File or directory | Use |
| --- | --- |
| `README.md` | Quick start and repository overview. |
| `docs/review/README.md` | Review-bundle description. |
| `docs/review/modcredb_schema.sql` | Full SQLite schema snapshot. |
| `docs/rebuild_staging_with_jaspar2024.md` | Reproducible rebuild instructions. |
| `docs/audits/` | JASPAR, HOCOMOCO, matrix, source-mapping, and rebuild audits. |
| `tests/fixtures/tf_webdb_sample.sqlite` | Small redacted review fixture, not the full private database. |

## 14. Safe Review Rule

When asking another person or an AI to review ModCREDB, share source code, this guide, schema exports, audit reports, and the small review fixture. Do **not** share the full private SQLite database, raw model archive, or operational metadata containing local filesystem paths unless there is an explicit data-sharing decision.
