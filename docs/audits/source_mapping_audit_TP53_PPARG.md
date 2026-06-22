# Source Mapping Audit: TP53/P04637 and PPARG/P37231

## Scope and safeguards

This is a read-only audit of the ModCREDB staging database and the supplied
raw archives. It does not modify application code, importers, SQLite schema,
SQLite rows, or raw data. The question is whether the public motif links for
two reviewed human TF sequence records are consistent with the supplied chart
and motif archives.

Database audited:

```text
data/tf_webdb_hocomoco_staging.sqlite
```

Raw input folder audited:

```text
/home/patricia/TF_database_Baldo_data
```

## Executive summary

The SQLite database is faithfully representing the integrated hierarchical
chart for both test cases.

- **PPARG/P37231 is a positive control.** The chart lists both JASPAR
  `MA0065.1` and HOCOMOCO `PPARG_HUMAN.H11MO.0.A`; the database preserves both
  as direct/identical supporting evidence.
- **TP53/P04637 is HOCOMOCO-only because the integrated chart makes it
  HOCOMOCO-only.** `P04637` is absent from the original chart and was added by
  the HOCOMOCO integration with only `P53_HUMAN.H11MO.0.A`.
- JASPAR `MA0106.1`, `MA0106.2`, and `MA0106.3` are present as valid local
  JASPAR motif files. Their raw MEME headers identify them as TP53 motifs.
  However, only `MA0106.1` has chart-derived links, and these are to numerous
  mostly unreviewed TP53-like sequence records, not to reviewed `P04637`.
- `MA0106.3` has no `motif_ref` rows because no supplied chart entry references
  it. ModCREDB correctly does not create a TF--motif link merely because a
  motif's title contains `TP53`.

**Classification:** this is primarily an **expected original/integrated chart
behavior with a missing JASPAR mapping/curation decision for P04637**, not a
matrix parser or motif-link importer bug. The source-release table also needs
a later metadata update because the raw archive itself identifies JASPAR 2024
and CIS-BP 2019, while their database release fields still say pending
confirmation.

## Database structure and source metadata

Relevant tables in the staging database:

| Table | Role in this audit |
| --- | --- |
| `tf` | One TF-like sequence record per accession. |
| `tf_annotation` | UniProt annotation: entry name, gene names, review status, length. |
| `tf_primary_annotation` | One final hierarchical label per TF record. |
| `motif_file` | Imported local MEME motif record and parsed probability matrix. |
| `motif_ref` | All TF-to-motif supporting-evidence links. |
| `source` / `source_release` | Source labels and release/provenance metadata. |
| `metadata` | Import-level provenance and counts. |

The source metadata is useful but incomplete for the supplied archives:

| Source | Database `source_release` | Raw input evidence | Audit status |
| --- | --- | --- | --- |
| JASPAR | `pending confirmation` | `jaspar/jaspar_2024/PWMS/*.meme` | Raw archive supports JASPAR 2024; DB metadata is stale/pending. |
| CIS-BP | `pending confirmation` | `pbm/CisBP_2019/pwms/*_2.00.meme` | Raw archive supports CIS-BP 2019 and motif version 2.00 naming; DB metadata is stale/pending. |
| HOCOMOCO | HOCOMOCO v11, CORE, human | Separate downloaded HOCOMOCO v11 human CORE mononucleotide source | Source release confirmed; exact TF-mapping semantics remain pending confirmation. |
| ModCRE | ModCRE/ModCRElib-derived precomputed dataset | Supplied ModCRE PWM/model archive | Exact scripts, version, and settings are pending PI/Baldo confirmation. |
| AlphaFold | AF3-derived precomputed dataset | Supplied AF3 PWM/model archive | Exact generation workflow is pending PI/Baldo confirmation. |
| UniProt | Fetch timestamp stored | Imported annotations | Exact UniProt release snapshot remains pending confirmation. |

The staging database metadata records 5,384 original chart rows and HOCOMOCO
integration paths. The `hocomoco_hierarchical_tsv` metadata path still points
to an earlier Downloads location; this is provenance housekeeping only and
does not affect the current TP53/PPARG links.

## Raw archive inventory

Archive member counts and top-level layouts were read without extracting the
large archives:

| Archive | Members | Apparent layout / interpretation |
| --- | ---: | --- |
| `jaspar.tar.gz` | 4,280 | `jaspar/jaspar_2024/PWMS/<motif>.meme` |
| `cisbp.tar.gz` | 10,035 | `pbm/CisBP_2019/pwms/<motif>_2.00.meme` |
| `pwms.tar.gz` | 11,615 | `pwms_modcre_all` (7,319) and `pwms_af3_all` (4,296) |
| `models.tar.gz` | 18,230 | `models` (11,156), `models_failed` (2,366), `MODELS_AF3` (4,297), and `MODELS_AF3_failed` (411) |

These are **archive-member counts**, not all PDB-model counts. The model
directories also contain summary/PIR records and directory entries.

JASPAR member names encode motif IDs, not a UniProt accession. The raw MEME
headers provide names such as `MOTIF MA0106.3 TP53`. CIS-BP member names are
database motif IDs such as `M03438_2.00` and `M11197_2.00`. HOCOMOCO is not
inside `pwms.tar.gz`; it was supplied separately as a downloaded HOCOMOCO v11
MEME collection.

## Integrated chart audit

The canonical integrated chart has 5,418 lines: one header plus 5,417 TF
rows. Its relevant columns are `TF_name`, `TF_family`, the five evidence
columns, `Best_annotation_level`, `Best_PWM_or_model`, and
`N_nonempty_annotation_columns`.

| TF | Present in original chart? | Integrated `Identical_PWM` | Other evidence columns | Primary result |
| --- | --- | --- | --- | --- |
| P04637 (TP53) | No | `P53_HUMAN.H11MO.0.A` | Empty | `Identical_PWM`; one nonempty annotation column |
| P37231 (PPARG) | Yes | `MA0065.1;PPARG_HUMAN.H11MO.0.A` | Empty | `Identical_PWM`; one nonempty annotation column |

This distinction explains the database result. P04637 was a HOCOMOCO-added
record in the integrated chart. P37231 already had a JASPAR assignment, and
the integration appended the HOCOMOCO direct motif to the same evidence layer.

## TP53/P04637 audit

### TF record and primary annotation

`P04637` exists in the database as reviewed UniProt entry `P53_HUMAN`, with
gene names `TP53 P53`, sequence length 393, and family
`p53-related factors{6.3.1}`. Its primary annotation exactly matches the
integrated chart:

```text
best_annotation_level: Identical_PWM
best_pwm_or_model:     P53_HUMAN.H11MO.0.A
n_nonempty_annotation_columns: 1
```

It has one supporting motif link:

| Source | Motif | Evidence | Mapping type | Matrix status | Consensus |
| --- | --- | --- | --- | --- | --- |
| HOCOMOCO | `P53_HUMAN.H11MO.0.A` | identical | `public_database_mapping_unconfirmed` | usable | `AGACATGCCCAGACATGCCC` |

The mapping type intentionally remains conservative: the source is a confirmed
HOCOMOCO v11 collection, while the precise mapping policy is marked pending
confirmation rather than silently claiming a stronger external identifier
mapping than was documented.

### JASPAR MA0106 records

All three MA0106 versions are present as usable JASPAR motif files:

| Motif | Width | Link count | Reviewed links | TP53-gene links | P04637 link? |
| --- | ---: | ---: | ---: | ---: | --- |
| `MA0106.1` | 20 | 54 | 0 | 45 | No |
| `MA0106.2` | 15 | 0 | 0 | 0 | No |
| `MA0106.3` | 18 | 0 | 0 | 0 | No |

The `MA0106.1` links are dominated by unreviewed TP53-like records, with
chart-derived evidence/mapping values such as homology transfers. This is
consistent with the resource being **sequence-record-level**, not a canonical
one-gene-one-row resource.

Raw MEME content for both `MA0106.1` and `MA0106.3` identifies the motif as
TP53. The SQLite `motif_file.content` is byte-for-byte equivalent to the raw
archive member, and the stored `matrix_json` matches the MEME
letter-probability matrix row-by-row. For example, MA0106.1 position 1 is:

```text
raw MEME / SQLite matrix_json: A=0.294118, C=0.470588, G=0.117647, T=0.117647
row sum: 1.000000
```

For MA0106.3 position 1, the same check gives:

```text
raw MEME / SQLite matrix_json: A=0.433264, C=0.059557, G=0.376924, T=0.130255
row sum: 1.000000
```

Therefore, MA0106.3 being unlinked is **not** caused by a matrix parsing or
normalization problem. It is a valid local motif file without a chart-derived
TF link.

## PPARG/P37231 positive control

`P37231` is reviewed UniProt entry `PPARG_HUMAN` (gene names `PPARG NR1C3`,
length 505). Its primary annotation is:

```text
best_annotation_level: Identical_PWM
best_pwm_or_model:     MA0065.1;PPARG_HUMAN.H11MO.0.A
n_nonempty_annotation_columns: 1
```

The database retains both direct supporting records from the integrated chart:

| Source | Motif | Evidence | Mapping type | Matrix status | Consensus |
| --- | --- | --- | --- | --- | --- |
| JASPAR | `MA0065.1` | identical | `direct_or_identical` | usable | `CCAGGGGTCAAAGGTCATCG` |
| HOCOMOCO | `PPARG_HUMAN.H11MO.0.A` | identical | `public_database_mapping_unconfirmed` | usable | `AAGTGGGGCAAAGGTCA` |

This is the expected direct multi-source pattern. It demonstrates that the
importer can preserve JASPAR and HOCOMOCO together for a reviewed canonical
record when both are present in the integrated chart.

## Broader direct-source sanity checks

Among reviewed TF sequence records whose primary annotation is
`Identical_PWM`, the database contains direct supporting links from all three
public motif sources:

| Source | Distinct reviewed TF records with direct/identical links | Links |
| --- | ---: | ---: |
| JASPAR | 531 | 531 |
| CIS-BP | 254 | 254 |
| HOCOMOCO | 401 | 401 |

P04637 is not an isolated HOCOMOCO-only record. Other reviewed HOCOMOCO direct
examples include BMAL1/O00327, BACH1/O14867, TP73/O15350, and CLOCK/O15516.

JASPAR includes a broad local archive rather than only motifs referenced by the
chart:

```text
JASPAR motif_file records: 4,279
JASPAR motif_file records with zero motif_ref links: 3,563
```

An unlinked local JASPAR motif is therefore a normal state in the current
archive-plus-chart import design. It is retained for provenance and possible
future curation but is not automatically attached to a TF.

## Interpretation and recommended next action

### What is correct

1. The integrated chart and database agree for P04637 and P37231.
2. The JASPAR raw MEME probability matrices are preserved correctly in SQLite.
3. The importer preserves multiple supporting records when the chart supplies
   them; PPARG is the positive control.
4. The importer appropriately avoids inferring a TF link from a motif title.

### What requires follow-up

1. **JASPAR mapping curation for TP53/P04637:** Ask the PI/Baldo whether
   reviewed P04637 should explicitly receive `MA0106.1`, `MA0106.3`, another
   current JASPAR TP53 motif/version, or remain HOCOMOCO-only by design.
2. **Source metadata update:** Record JASPAR 2024 and CIS-BP 2019/2.00 in the
   source-release metadata after confirming these archive labels are the
   intended release provenance.
3. **HOCOMOCO mapping provenance:** Retain the current conservative
   `public_database_mapping_unconfirmed` label until the mapping procedure is
   confirmed.

### Recommended action now

**No database change and no importer change.** The correct next action is an
auditable manual-curation decision for TP53/P04637, backed by authoritative
JASPAR TF-to-motif metadata or written PI/Baldo confirmation. If approved,
add the selected TP53 JASPAR mapping explicitly with a provenance note. Do not
infer it simply from `TP53` appearing in a MEME title.

## Limitations

- This audit compares supplied archives, the supplied/integrated chart, and
  the local SQLite database. It does not assert which JASPAR TP53 version is
  scientifically preferred.
- The raw JASPAR archive identifies its directory as `jaspar_2024`, but the
  database `source_release` row has not yet been updated from pending
  confirmation.
- ModCRE and AF3 model-generation settings were outside the scope of this
  source-mapping audit.
