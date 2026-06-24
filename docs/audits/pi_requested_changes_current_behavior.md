# PI-requested terminology and search audit — current local behavior

Audited: 2026-06-24  
Scope: read-only inspection of the current local application and configured SQLite database. No application code, database, schema, or import data was changed. The two files in this audit are the only intended additions.

## Executive result

The configured local application database is [`data/tf_webdb.sqlite`](../../data/tf_webdb.sqlite) (`/home/patricia/tf_webdb/data/tf_webdb.sqlite`). It exists and `PRAGMA integrity_check` returned `ok`.

The app has no running local server process at audit time, so this is the configured default database rather than a runtime-observed `--db` override. `app.py` defaults to this path and supports an explicit `--db` argument; no `.env`, startup script, or database environment-variable override was found. A deployment service that invokes `app.py --db <other-path>` would need to be checked separately.

PFAM search is **not currently available as verified PFAM search**. The database has no PFAM table, accession, name, source, or provenance column. The normalized `tf_family` values come directly from the chart's `TF_family` text. Many values resemble PFAM labels (for example `zf-C2H2`, `Homeobox`, and `Fork_head`), but there are no `PFxxxxx` accessions and no recorded PFAM source. Publicly renaming this UI to **PFAM Families (by similarity)** would therefore be unsupported by the current data.

Most literal terminology substitutions are display-only changes. The exceptions are PFAM terminology/search, the AlphaFold “Full sequence” region behavior, and complete monomer/dimer grouping: those require implementation logic and, for scientifically complete classification, a data/provenance decision.

## 1. Active configured database

| Check | Result | Evidence |
| --- | --- | --- |
| Configured default path | `/home/patricia/tf_webdb/data/tf_webdb.sqlite` | `app.py:34` (`DEFAULT_DB`); `app.py:2898` (`--db`, defaulting to `DEFAULT_DB`) |
| Exists | Yes | Opened successfully with SQLite URI `mode=ro` |
| Runtime override | Possible with `--db`; not observed because no app process was running | `app.py:2896-2910`; `README.md:82,88` |
| Environment override | No DB-path environment override found | No local `.env`; only unrelated `TF_WEBDB_MAX_POST_BYTES` and model-cache environment settings are read in `app.py:35-40` |
| `PRAGMA integrity_check` | `ok` | Read-only SQLite query |

### Table list and row counts

| Table | Rows | Role / relevant schema summary |
| --- | ---: | --- |
| `tf` | 5,384 | One UniProt TF record: `tf_id`, chart-derived `family_text`, cached motif/model counts. |
| `tf_annotation` | 5,384 | UniProt annotation: accession, entry/gene/protein/organism, review status, sequence length, URL, fetch time. |
| `tf_family` | 36,636 | Normalized split values of the chart `TF_family` field: `(tf_id, family)` primary key. |
| `motif_file` | 25,929 | Local motif file and parsed PWM/QC fields including `matrix_status` and `matrix_json`. |
| `motif_ref` | 112,693 | TF-to-motif links: evidence type, source, motif ID, identity, mapping/curation fields. |
| `motif_structure` | 12,453 | Exact link table from a motif reference to an indexed structure. |
| `structure_file` | 18,223 | Indexed files: source, model ID, archive/member path, type, status, template, residue interval. |
| `model_summary` | 7,700 | Parsed model-summary rows: source/status, model rank, template, tails/chains, coverage/identity/similarity. |
| `source` | 6 | Source labels/descriptions: `alphafold`, `cisbp`, `hocomoco`, `jaspar`, `modcre`, `uniprot`. |
| `source_release` | 6 | Source release/provenance metadata. |
| `metadata` | 8 | Import counters and local source-data metadata. |
| `import_issue` | 31 | Import/audit issues. |
| `sqlite_sequence` | internal | SQLite AUTOINCREMENT bookkeeping. |

`tf_primary_annotation` is **not present** in this active database, although the application handles it as optional. The UI falls back to its evidence/model-derived primary-status logic (`app.py:257-303`, `app.py:581-618`).

The metadata reports 7,700 model-summary rows, 18,223 indexed model files, 25,929 motif files, and 5,384 UniProt annotations. The UniProt annotation fetch timestamp is `2026-05-26T13:47:13+00:00`.

## 2. PFAM/family audit

### Schema and provenance results

| Question | Result | Evidence |
| --- | --- | --- |
| Dedicated PFAM table? | No | Complete table list has no PFAM table. |
| PFAM ID column, such as `PF00096`? | No | No column name contains `pfam` or `interpro`; zero `tf_family.family` values match a `PF` + five-digit accession pattern. |
| PFAM name column? | No | Only `tf.family_text`, `tf_family.family`, and `model_summary.domain` contain family/domain-like terms. `model_summary.domain` is model-summary metadata, not a TF PFAM assignment. |
| Is the existing family field from old `TF_family`? | Yes | `import_db.py:885-895` reads `row["TF_family"]`, stores it in `tf.family_text`, and splits it into `tf_family.family`. |
| Are values proven PFAM assignments? | No | Values look PFAM-like but have neither accession nor provenance. Do not infer a PFAM assignment from a label alone. |
| Does current search include the existing family text? | Yes, via `tf.family_text` | `app.py:1877-1885`; the route does not join or search the normalized `tf_family` table directly. |

The most frequent `tf_family.family` values include `zf-C2H2_4` (1,646), `zf-C2H2` (1,640), `zf-H2C2_2` (1,602), `zf-C2H2_jaz` (1,502), `Homeobox` (751), `Fork_head` (174), `bZIP_1` (166), and `Hormone_recep` (303). These are compatible with PFAM-style labels but are not sufficiently sourced to make a public PFAM claim.

There are 2,391 exact matches among a deliberately narrow set of names (`zf-C2H2`, `Homeobox`, `Forkhead`, `bZIP`, `nuclear receptor`); that is not PFAM evidence. `Homeobox` is stored in family text (751 TF records). `Forkhead` is not stored as that exact family string (`Fork_head` is stored for 179 records); the `Forkhead` free-text search hits are therefore attributable to the other route fields, especially protein annotations, rather than verified PFAM-family search.

### Recommendation

Keep public labels as **Families** or **TF families** until a source with explicit PFAM accessions and provenance is imported. A proper PFAM implementation needs an assignment relation such as `(tf_id, pfam_id, pfam_name, source, release, method/provenance)`, an importer or controlled mapping workflow, and explicit search support. It should not relabel the legacy chart field by assumption.

## 3. Search behavior

### What `/search` currently searches

The TF-result query joins `tf`, `motif_ref`, and `tf_annotation` and matches a non-empty `q` against:

* `tf.tf_id` — UniProt sequence-record accession;
* `tf.family_text` — the raw chart `TF_family` text;
* `motif_ref.motif_id`;
* `tf_annotation.gene_names`, `protein_name`, and `organism_name`.

The separate motif-result query matches `motif_file.motif_id`, `motif_file.source`, and the linked TF gene/protein text. The form also has exact `source` and `evidence` filters (`templates/search.html:51-65`); `source` filtering is applied to `motif_ref.source` for TF results and `motif_file.source` for motif results.

| Search capability | Current result | Notes |
| --- | --- | --- |
| UniProt accession | Yes | `tf.tf_id LIKE ?`; the requested `P04637` has no row in this active DB, so that control returns zero results. |
| Gene symbol | Yes | `gene_names LIKE ?`; a token-aware gene-summary view is also produced for symbol-like inputs. |
| Motif ID | Yes | `motif_ref.motif_id` in TF results and `motif_file.motif_id` in motif results. |
| Source | Yes, as a form filter; motif free-text also matches source | Use `?source=jaspar`, etc. Free-text source is not a main TF-result predicate. |
| TF family | Yes, raw chart text only | `tf.family_text LIKE ?`; not the normalized table and not PFAM-aware. |
| PFAM ID/name | No verified PFAM field exists | `PF00096` returns nothing. PFAM-like legacy text happens to be free-text searchable. |

### Requested query checks

No server was running, so these are read-only executions of the route's SQL predicates rather than HTTP requests. TF and motif counts below are unbounded query matches before the UI's result limits (60 TFs, 100 motif files).

| Query | TF matches | Motif-file matches | Interpretation |
| --- | ---: | ---: | --- |
| `P04637` | 0 | 0 | The accession is absent from both `tf` and `tf_annotation` in this active DB. This does not disprove accession-search support. |
| `TP53` | 58 | 6 | Gene/protein text is searchable; result is sequence-record level. |
| `MA0106.3` | 0 | 1 | Direct motif-file search returns `jaspar|MA0106.3`. |
| `P53_HUMAN.H11MO.0.A` | 0 | 0 | This HOCOMOCO motif is absent in this active DB. |
| `C2H2` | 1,674 | 310 | Legacy family-text search works. |
| `PF00096` | 0 | 0 | No PFAM accession data is present. |
| `Homeobox` | 796 | 1,321 | Legacy family text contributes matching records. |
| `Forkhead` | 157 | 169 | Matches other searchable text; legacy stored family spelling is `Fork_head`, not `Forkhead`. |

The source filter has current TF/motif-file coverage for `alphafold` (560 / 4,296), `cisbp` (3,199 / 10,035), `jaspar` (3,627 / 4,279), and `modcre` (1,082 / 7,319). The UI also exposes HOCOMOCO and UniProt filters, but this configured database currently has zero motif-file and TF-link results for those two sources.

### Ranking

TF rows are loaded distinctly, status-enriched, then sorted in Python (`app.py:370-382`, `app.py:1905-1913`) by:

1. exact accession;
2. exact gene-name token;
3. reviewed before unreviewed;
4. primary-evidence priority (identical, homologous, relative-homologous, ModCRE, AlphaFold, unannotated);
5. scan capability;
6. more scan-ready PWMs;
7. more active models;
8. accession.

Motif rows are SQL-sorted (`app.py:418-545`) by exact motif ID, preferred gene-summary TF link, exact gene link, evidence rank, usable matrix, mapping rank, reviewed linked TF, linked-TF count, source priority (JASPAR, CIS-BP, HOCOMOCO, ModCRE, AlphaFold), then source and motif ID.

## 4. Public UI terminology inventory

`static/style.css` contains class names only; no requested public labels are emitted from static assets. The following inventory covers the template text and the route/configuration values that render it. Stored database codes (for example `evidence_type`, `matrix_status`, and `source`) should remain unchanged unless an explicit migration is approved.

| PI request | Current public location/text | Recommended public label/action | Change class |
| --- | --- | --- | --- |
| `evidence` → `prediction` | Broad public use in `templates/index.html`, `search.html`, `tf.html`, `motif.html`, `docs.html`, and `evidence.html`; server values in `app.py` feed labels | Change only user-facing wording where the PI means *prediction*. Keep DB/API field names and the `evidence` URL/filter parameter stable. | UI-only once wording scope is confirmed; broad wording is a decision. |
| `Identical PWM` → `Known` | `app.py:42`, rendered in search, scan, TF, motif, docs, evidence-guide templates | Change `EVIDENCE_LABELS['identical']` to `Known`; update explanatory docs. | UI-only. |
| `Homologous PWM` → `Nearest Neighbor (>70%)` | `app.py:43` and the same rendered label locations | Change the label and related docs. | UI-only label, but retain current identity semantics. |
| `Relatively Homologous PWM` → `Nearest Neighbor (70% - 40%)` | `app.py:44`; current description says 50.0–69.9% (`app.py:62-64`) | Do not silently publish `70% - 40%`: it overlaps the prior label at 70 and contradicts the documented current 50–69.9% interval. Confirm whether the intended lower bound is 40 or 50 and whether stored links need reclassification. | Needs PI decision; may need data/import work. |
| `Families` → `PFAM Families (by similarity)` | `templates/tf.html:83-85`, `templates/motif.html:227`, `templates/search.html:118`, plus home/docs wording | Do **not** use PFAM label with current data. Keep `Families`/`TF families` pending PFAM assignments and provenance. | DB/import/search/provenance work required. |
| `Scan this region` → `Scan` | `templates/tf.html:103` | Change link text to `Scan`. | UI-only. |
| `# active PDB models` → `# 3D models` | Current variants: `active models` (`templates/tf.html:19`, `docs.html:60`), `active PDB model files` (`templates/docs.html:22`, `tf.html:261,263,276`) | Replace public count labels with `3D models`; retain `status='active' AND file_type='pdb'` predicates. | UI-only. |
| `w=0 / no matrix` → `Missing MEME` | `app.py:132`; shown via the matrix-status map in TF, motif, search, scan, and docs pages | Change the display map for `width_zero_no_matrix` to `Missing MEME`; retain the raw status code for QC. | UI-only. |
| `FIMO Ready` → `Generated PWM` | Current spelling is `FIMO-ready` in `app.py:131` and many templates, e.g. `tf.html:168`, `search.html:25`, `motif.html:12` | Change display badges/text only if PI intends it. `FIMO-ready` means scanner-usable, whereas `Generated PWM` may be read as provenance, so confirm desired scientific meaning. | UI-only after terminology decision. |
| `active models` → `3D models` | `templates/tf.html:19,36,112`, `search.html:78,122`, `docs.html:21,60`, and public prose | Replace display text, not active-status filtering. | UI-only. |
| `Motif evidence` → `Motif Prediction` | Primary heading is `templates/tf.html:128-133`; additional prose in index/search/docs/motif/evidence-guide templates | Change this heading to `Motif Prediction`; audit surrounding prose under the broad evidence→prediction decision. | UI-only for the heading; wording decision for all prose. |
| `PREDICTED = LOW` | `app.py:65,70` supplies `Predicted` trust badges rendered at `templates/evidence.html:11-19`; docs also use predicted wording | If PI means the badge value, render `Low` (or `LOW`) for ModCRE/AlphaFold prediction cards. Confirm capitalization and whether this is a confidence category rather than a general prose replacement. | UI-only after decision. |
| `Active model files` → `3D Models` | `templates/tf.html:263`; nearby docs say `active model files` at `docs.html:401` | Change the summary/display label to `3D Models`; retain debug/internal wording where it deliberately describes files. | UI-only. |
| `TF Motif DB` → `ModCRE DB` | No exact `TF Motif DB` string appears in public routes/templates. `make_grant_stats_pdf.py:313,360` has the non-runtime report title `TF Motif and Model Database`. | No live UI change found. Rename the report title only if PI includes generated grant reports. | UI/report-only. |
| Remove the quoted `w=0` explanatory text | The exact quoted two-sentence string is not present. Equivalent public text is in `app.py:149`, `templates/motif.html:188-199`, and `templates/docs.html:378-380`. | Remove/rewrite those display explanations while retaining a short non-diagnostic `Missing MEME` status. | UI-only. |
| Model summaries: remove rank, duplicates, residues, model | Embedded summary table: `templates/tf.html:280-289`; full page: `templates/model_summaries.html:23-69` | Remove specified display columns; define duplicate key before adding display deduplication. Do not delete DB rows. | UI/query-only after duplicate-key decision. |
| AlphaFold `Full sequence` region | No current special case; `app.py:1194-1276` requires numeric intervals. | Add an explicit AlphaFold full-sequence group and matching scan filter; see section 6. | Route/template logic required; no schema change for display-only version. |
| Split 3D models into DIMERS and MONOMERS | Current `templates/tf.html:258-296` renders one combined active-PDB list. | Classify only explicit `_dimer`/`_monomer` paths at display time; preserve unclassified models as a third group pending provenance. | UI-only partial; data/import work for complete classification. |

### Model-summary column distinction

The embedded TF-page summary table displays `Rank`, `Residues`, and `Model` (`templates/tf.html:283-288`). The full `/model-summaries/<TF>` table displays `Rank` and `Model file` but does not have a `Residues` column (`templates/model_summaries.html:23-69`). The separate structural-model and motif matched-structure tables also have `Model`/`Residues` columns (`templates/tf.html:266-269`, `templates/motif.html:260-280`); they are not model-summary displays. The PI should confirm whether the removal applies only to model summaries, as stated, or also to those distinct structural tables.

## 5. Model-summary display audit

### Storage and current display

`model_summary` is the storage table. It has 7,700 total rows, 5,988 active rows, and all 7,700 currently link to a `structure_file`. Active rows are all `source='modcre'`; AlphaFold has indexed PDB files but no rows in `model_summary`.

The TF detail page selects at most 12 active summary rows, ordered by identity percentage, similarity percentage, then model rank (`app.py:2049-2060`). The full summary page selects at most 500 rows using the same ordering (`app.py:2424-2439`). The columns rendered are listed in the inventory above.

### Duplicates

No duplicates were found under any of the following read-only checks:

* exact equality of every stored summary content field other than `id`;
* active `(tf_id, source, status, summary_model_id, model_rank)`;
* active `(tf_id, source, status, matched_structure_id)`.

Therefore the requested “remove duplicates” has no current duplicate data to remove under these defensible keys. It can be implemented at display/query level (for example, a `ROW_NUMBER()` partition or Python de-duplication), without deleting database rows, **after the PI chooses the intended identity key and winner rule**. A silent key such as template alone would collapse scientifically distinct model summaries.

### Source and oligomeric state

Model source is available: `structure_file.source` and `model_summary.source` (`modcre` or `alphafold`). No explicit oligomeric-state database column exists. `model_summary.protein_chain` is free text (1,536 active rows contain comma-separated chains), which is not by itself a safe biological monomer/dimer classifier.

The archive path/model ID does explicitly mark some structures, however:

| Active PDB path classification | Count |
| --- | ---: |
| Explicit `dimer` path/model ID | 1,941 |
| Explicit `monomer` path/model ID | 402 |
| Unclassified by path/model ID | 9,729 |

There are 67 TFs with both explicitly marked monomer and dimer active PDB models. This supports a conservative display-only split for explicitly labelled models. It does **not** support assigning the 9,729 unclassified models to either category without confirmation or a data-level classification workflow.

## 6. AlphaFold/AF3 DNA-binding-region audit

### Current behavior

The TF-page DNA-binding regions are built from numeric intervals in motif-region aggregates, active `structure_file` records, and parsed model summaries (`app.py:1194-1276`). The scan-page region selector uses only active, usable-PWM motif-to-structure links with non-null `residue_start`/`residue_end` (`app.py:938-989`). Selecting a region then applies an overlap condition to that numeric interval (`app.py:1018-1034`).

There is no AlphaFold-specific branch. The source identifier is `structure_file.source = 'alphafold'`, with `source_labels['alphafold']` rendered as `AlphaFold3-assisted ModCRE` (`app.py:77-84`).

The configured DB has 4,296 active AlphaFold PDB files for 561 TFs, and **all 4,296 have null residue intervals**. It has 3,360 motif-to-AlphaFold-structure links, so the data exists but cannot enter the current numeric-region UI. Of the 561 AlphaFold-model TFs, 560 have a positive UniProt sequence length; one has no usable sequence length for a safe full-sequence interval.

### Exact implementation required for “Full sequence”

This is not currently possible. A display-only implementation needs:

1. In `build_region_groups` (`app.py:1194-1276`), collect active AlphaFold PDB models with a known `tf.sequence_length` and create a separate group flagged, for example, `is_full_sequence=True`, `start=1`, `end=sequence_length`, and `models=<AlphaFold models>`. Do not pretend null residue intervals are measured DNA-binding intervals.
2. In `templates/tf.html:89-125`, render that flagged group as **Full sequence** rather than `Region 1-N`, and render the AlphaFold model IDs/source/action links within the group so the results are actually visible there (the current panels show counts only).
3. If the scan selector is also meant to support the new region, extend `load_tf_scan_regions`, `parse_region_value`, and `load_tf_scan_motifs` (`app.py:927-1147`) with a distinct sentinel such as `full_sequence`, not a fake numeric interval. Its SQL must include usable motifs linked through `motif_structure` to active `structure_file.source='alphafold'` records rather than using the existing residue-overlap predicate.
4. In `templates/scan.html:35-50`, render `Full sequence` and ensure the TF-page scan link passes the sentinel. For the one TF without a sequence length, omit this group or label it as a data-quality exception.

No database-schema change is necessary for that conservative UI/route implementation. Persisting an explicit region type or predicted coverage interval would be a separate data-model decision.

## 7. Implementation sequence

1. **PI terminology decisions:** confirm the intended meaning of `PREDICTED = LOW`; whether `Generated PWM` deliberately replaces the scanner-availability meaning of `FIMO-ready`; and whether the 40% lower bound is correct for the relative-neighbor label.
2. **Safe UI terminology patch:** centralize display mappings in `app.py`, update the listed templates/docs, retain raw codes and URL/query parameters, and remove the requested no-matrix explanatory prose.
3. **Model-summary UI patch:** remove only the specified summary-table columns; add no deduplication until the PI supplies a semantic duplicate key. Do not delete rows.
4. **3D-model grouping:** implement explicit path-based `DIMERS`, `MONOMERS`, and `Unclassified 3D Models`; surface the 67 TFs that have both. Obtain source confirmation/import classification before hiding or recategorizing unclassified models.
5. **AlphaFold Full sequence:** implement the flagged full-sequence group and, if requested, a matching scan sentinel using the exact AlphaFold motif-structure links.
6. **PFAM:** obtain/import release-pinned PFAM assignments with IDs, names, sources, and method/provenance; then add schema/search/UI support and only then use `PFAM Families (by similarity)`.

## 8. Worktree checks

`git status --short` before this audit already showed user changes outside the two new audit files:

```text
 M README.md
 M app.py
 M config/source_releases.tsv
 M docs/audits/hocomoco_evidence_type_audit.md
 M static/style.css
 M templates/index.html
 M templates/scan.html
 M templates/search.html
 M templates/tf.html
?? docs/MODCREDB_TECHNICAL_GUIDE.md
?? docs/audits/cisbp_metadata_link_preflight.md
?? docs/audits/hocomoco_official_mapping_audit.md
?? docs/audits/hocomoco_official_mapping_audit.tsv
?? docs/audits/jaspar2024_rebuild_validation/
?? docs/audits/jaspar_metadata_link_apply.md
?? docs/audits/jaspar_metadata_link_candidates.tsv
?? docs/audits/jaspar_metadata_link_dry_run.md
?? docs/audits/jaspar_uniprot_link_audit.md
?? docs/audits/jaspar_uniprot_link_audit.tsv
?? docs/rebuild_staging_with_jaspar2024.md
?? scripts/
```

The final status and whitespace check are recorded after the TSV is added.
