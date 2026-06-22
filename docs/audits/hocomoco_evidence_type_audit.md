# HOCOMOCO Evidence-Type Audit

## Scope and safeguards

This is a read-only audit of the current ModCREDB staging database, the
integrated HOCOMOCO hierarchical chart, and the importer source code. It does
not modify application code, SQLite schema or rows, raw data, import scripts,
or FIMO behavior.

Database audited:

\`\`\`text
data/tf_webdb_hocomoco_staging.sqlite
\`\`\`

Chart audited:

\`\`\`text
/home/patricia/TF_database_Baldo_data/TF_PWM_chart_final_integrated_HOCOMOCO_hierarchical.tsv
\`\`\`

## Executive summary

**There is no current HOCOMOCO evidence-tier mismatch in the staging database.**

- The database contains **401** HOCOMOCO \`motif_ref\` rows.
- All 401 have \`evidence_type='identical'\` and mapping type
  \`public_database_mapping_unconfirmed\`.
- All 401 belong to TFs whose primary annotation is \`Identical_PWM\`.
- Comparing every HOCOMOCO link with the integrated TSV found all 401 motif
  tokens in that TF's \`Identical_PWM\` column. There are zero links in
  \`Homologous_PWM\`, \`Relatively_Homologous_PWM\`, \`ModCRE\`, or
  \`AlphaFold\`; zero Best-only assignments; zero unmatched links; and zero
  evidence-type mismatches.

However, the current **importer code has a real future-rebuild risk**:
\`import_hocomoco.py\` scans HOCOMOCO IDs from several evidence columns but
inserts every one with the hard-coded evidence type \`identical\`. It also stores
the generic original column \`HOCOMOCO\`, not the source TSV column. A future
integrated chart with a HOCOMOCO motif in \`Homologous_PWM\` or
\`Relatively_Homologous_PWM\` would therefore be mislabeled as direct/identical
evidence.

**Classification:** importer code bug risk, but the current database is
unaffected.

**Recommendation:** make a code-only importer fix before any future rebuild;
do not rebuild or edit the current staging database for this issue.

## Code inspection

### Hard-coded HOCOMOCO evidence type

\`import_hocomoco.py\` explicitly sets HOCOMOCO references to \`identical\`:

\`\`\`python
# import_hocomoco.py:25-26
HOCOMOCO_SOURCE = "hocomoco"
HOCOMOCO_EVIDENCE = "identical"
\`\`\`

### HOCOMOCO tokens are collected from multiple hierarchy columns

The importer scans all of the following fields, without retaining the field
that supplied each token:

\`\`\`python
# import_hocomoco.py:59-73
def hocomoco_tokens(row: dict[str, str]) -> list[str]:
    tokens: list[str] = []
    for column in (
        "Identical_PWM",
        "Homologous_PWM",
        "Relatively_Homologous_PWM",
        "ModCRE",
        "AlphaFold",
        "Best_PWM_or_model",
    ):
        for token in split_list_field(row.get(column, "")):
            token = token.strip()
            if ".H11MO." in token and token not in tokens:
                tokens.append(token)
    return tokens
\`\`\`

The importer uses this function both to find all motif IDs
(\`import_hocomoco.py:102-104\`) and to add TF-to-motif links
(\`import_hocomoco.py:189-199\`). At insertion, it unconditionally uses the
constant rather than an evidence tier derived from the column:

\`\`\`python
# import_hocomoco.py:189-211
for motif_id in hocomoco_tokens(row):
    ...
    existing_ref = conn.execute(..., (
        tf_id, HOCOMOCO_EVIDENCE, HOCOMOCO_SOURCE, motif_id
    )).fetchone()
    semantics = motif_ref_semantics(
        HOCOMOCO_EVIDENCE, HOCOMOCO_SOURCE, "HOCOMOCO"
    )
    ...
    VALUES (..., HOCOMOCO_EVIDENCE, HOCOMOCO_SOURCE, motif_id, ...)
\`\`\`

Therefore:

| Question | Code answer |
| --- | --- |
| Is \`HOCOMOCO_EVIDENCE\` hard-coded to \`identical\`? | Yes, at line 26. |
| Which fields are scanned? | \`Identical_PWM\`, \`Homologous_PWM\`, \`Relatively_Homologous_PWM\`, \`ModCRE\`, \`AlphaFold\`, and \`Best_PWM_or_model\`. |
| Is the source TSV column retained for each HOCOMOCO token? | No. The helper returns only motif IDs. |
| Is evidence inferred from \`Best_annotation_level\`? | No. It is copied into \`tf_primary_annotation\` at lines 173-186, but not used to set HOCOMOCO \`motif_ref.evidence_type\`. |
| Can a future HOCOMOCO token from a homologous/candidate column be inserted as identical? | Yes. |

### Contrast with the base chart importer

The generic chart importer preserves the TSV evidence column correctly for
JASPAR, CIS-BP, ModCRE, and AlphaFold tokens:

\`\`\`python
# import_db.py:877-903
evidence_columns = [
    ("Identical_PWM", "identical"),
    ("Homologous_PWM", "homologous"),
    ("Relatively_Homologous_PWM", "relative_homologous"),
    ("ModCRE", "modcre"),
    ("AlphaFold", "alphafold"),
]
...
for column, evidence_type in evidence_columns:
    for raw_token in split_list_field(row.get(column, "")):
        motif_id, identity_percent = normalize_motif_token(raw_token)
        source = source_for_motif_id(motif_id, evidence_type)
        semantics = motif_ref_semantics(evidence_type, source, column)
\`\`\`

\`import_db.py:113-125\` applies HOCOMOCO-specific mapping semantics whenever
the source is HOCOMOCO. That does not itself overwrite the supplied
\`evidence_type\`, but the separate HOCOMOCO staging importer currently does so
through its hard-coded constant.

## Current database evidence distribution

The current staging DB has one HOCOMOCO evidence/mapping combination:

| Source | Evidence type | Mapping type | Links |
| --- | --- | --- | ---: |
| HOCOMOCO | identical | public_database_mapping_unconfirmed | 401 |

The same 401 links grouped by their TF primary annotation are:

| Primary annotation | HOCOMOCO evidence type | Mapping type | Links |
| --- | --- | --- | ---: |
| Identical_PWM | identical | public_database_mapping_unconfirmed | 401 |

The explicit suspicious-row query returned zero rows:

\`\`\`text
source = hocomoco
AND evidence_type = identical
AND primary annotation is not Identical_PWM

result: 0
\`\`\`

## TSV comparison for every HOCOMOCO motif link

For every \`motif_ref\` row with \`source='hocomoco'\`, the audit loaded the TSV
row with the same \`TF_name\`/UniProt accession and inspected semicolon-separated
tokens in these fields:

| TSV field containing HOCOMOCO token | Expected evidence type |
| --- | --- |
| \`Identical_PWM\` | identical |
| \`Homologous_PWM\` | homologous |
| \`Relatively_Homologous_PWM\` | relative_homologous |
| \`ModCRE\` | unexpected HOCOMOCO placement; flag |
| \`AlphaFold\` | unexpected HOCOMOCO placement; flag |
| \`Best_PWM_or_model\` only | infer from \`Best_annotation_level\` only if unambiguous |
| No field | unmatched; flag |

Results:

| Comparison category | Links |
| --- | ---: |
| Total HOCOMOCO \`motif_ref\` rows | 401 |
| Matched \`Identical_PWM\` | 401 |
| Matched \`Homologous_PWM\` | 0 |
| Matched \`Relatively_Homologous_PWM\` | 0 |
| Matched \`ModCRE\` | 0 |
| Matched \`AlphaFold\` | 0 |
| Only \`Best_PWM_or_model\` | 0 |
| Unmatched in the same TF TSV row | 0 |
| Expected evidence type equals DB evidence type | 401 |
| Expected evidence type differs from DB evidence type | 0 |

There are no mismatch examples to list: the requested first 100 mismatch table
is empty.

Examples from successful direct mappings:

| TF | HOCOMOCO motif | TSV column | DB evidence type |
| --- | --- | --- | --- |
| A0PJY2 | \`FEZF1_HUMAN.H11MO.0.C\` | \`Identical_PWM\` | identical |
| O00327 | \`BMAL1_HUMAN.H11MO.0.A\` | \`Identical_PWM\` | identical |
| O00470 | \`MEIS1_HUMAN.H11MO.0.B\` | \`Identical_PWM\` | identical |
| P04637 | \`P53_HUMAN.H11MO.0.A\` | \`Identical_PWM\` | identical |
| P37231 | \`PPARG_HUMAN.H11MO.0.A\` | \`Identical_PWM\` | identical |

## PPARG and TP53 checks

### PPARG/P37231

The TSV correctly records both public motifs in \`Identical_PWM\`:

\`\`\`text
MA0065.1;PPARG_HUMAN.H11MO.0.A
\`\`\`

The database preserves:

| Source | Motif | DB evidence type | Mapping type |
| --- | --- | --- | --- |
| JASPAR | \`MA0065.1\` | identical | direct_or_identical |
| HOCOMOCO | \`PPARG_HUMAN.H11MO.0.A\` | identical | public_database_mapping_unconfirmed |

Thus the HOCOMOCO \`identical\` label is correct for P37231.

### TP53/P04637

The integrated TSV records only:

\`\`\`text
P53_HUMAN.H11MO.0.A
\`\`\`

in \`Identical_PWM\`. The database gives this HOCOMOCO link evidence type
\`identical\` and mapping type \`public_database_mapping_unconfirmed\`. That is
correct for the current chart and is consistent with the separate TP53/PPARG
source-mapping audit.

## Classification and recommendation

**Classification: importer code bug risk, but current DB unaffected.**

The current input data happened to place all HOCOMOCO tokens in
\`Identical_PWM\`, so the importer's hard-coded label produced correct current
rows. The code nevertheless violates the evidence-tier preservation rule if a
future TSV contains HOCOMOCO tokens in other evidence columns.

**Recommendation: code fix required before future rebuilds only.**

The safe future implementation should make the HOCOMOCO token extractor return
both \`motif_id\` and its source TSV column/evidence type. It should preserve the
real column in \`original_column\`, use the corresponding evidence tier, and
deduplicate only within the appropriate TF/source/motif/evidence combination.

No current database rebuild, manual curation, or user-facing warning is
required from this audit alone.

## Repository state

Only this audit report was added for this task. No code, database, schema,
importer, raw data, or scanner file was modified. No commit was made.
