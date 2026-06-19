# ModCREDB Scientific Design Review Bundle

This folder contains safe, reviewable exports generated from the local HOCOMOCO staging SQLite database. The purpose is to review schema design, evidence modelling, motif/model records, FIMO scanner interpretation, and page organization before further coding.

The exports intentionally preserve all evidence types. They do not collapse TFs to a single best annotation. HOCOMOCO is reported as a separate public motif source with mapping semantics marked as unconfirmed.

Generated artifacts include:

- final HOCOMOCO hierarchical chart provenance (`final_hierarchical_chart_provenance.md`)
- schema SQL
- table/source/evidence inventories
- source-release provenance inventory
- matrix QC inventory (`matrix_qc_inventory.tsv`) and representative status examples (`matrix_qc_examples.tsv`)
- model and model-summary inventories
- TF/motif/model case-study reports
- FIMO scanner review notes
- redaction notes
- recommended public information architecture
- a small representative sample SQLite DB under `tests/fixtures/`

The sample DB preserves the real table layout but contains only representative TFs and redacted path fields.

PR3 evidence-semantics exports add:

- `mapping_type_counts.tsv`
- `curation_status_counts.tsv`
- `source_by_mapping_type_counts.tsv`

These files describe how TF-motif links are interpreted without collapsing evidence to one best annotation.
