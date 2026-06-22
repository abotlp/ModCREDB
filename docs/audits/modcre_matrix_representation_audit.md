# ModCRE Matrix Representation Audit

## Scope and safeguards

This is a read-only audit of the ModCREDB staging SQLite database, the
supplied PWM/model archives, and the extracted local files. It does not modify
application code, SQLite schema or rows, import scripts, FIMO behavior,
`matrix_json`, or raw data.

Database audited:

```text
data/tf_webdb_hocomoco_staging.sqlite
```

Raw archives audited:

```text
/home/patricia/TF_database_Baldo_data/pwms.tar.gz
/home/patricia/TF_database_Baldo_data/models.tar.gz
```

## Executive summary

1. **The local ModCRE import is correct for the verified local control.** For
   `TFS_tr_A0A024R0Y4_A0A024R0Y4_HUMAN:211:243_6r2v_A_1`, the raw archive MEME
   content, the `motif_file.content` field, and `motif_file.matrix_json` agree
   row-by-row. The local file is explicitly a MEME `letter-probability matrix`.
2. **The requested CTCF/P49711 live-webserver model is not present in the
   supplied local ModCRE/AF3 datasets.** P49711 has a usable HOCOMOCO motif,
   but has zero local ModCRE/AF3 motif links and zero local structure files.
   No local record combines P49711/CTCF, template `2wbu`, and residues
   `299:373`.
3. The complete 4-by-4 value block copied from the webserver screenshot was
   not found in any local `matrix_json`, local MEME content, or extracted MEME
   file, either in the listed orientation or transposed.
4. Therefore the live CTCF webserver matrix cannot currently be classified as
   the same or a transformed version of a local ModCREDB matrix. It is most
   conservatively classified as **not present in local data / possibly a
   different job, dataset version, or settings**. The representation used by
   that live page remains unresolved until its actual downloadable output is
   compared.

## Important interpretation of the screenshot numbers

The copied screenshot block was:

```text
0.106000  0.316000  0.504000  0.074000
0.001000  0.861000  0.043000  0.095000
0.001000  0.829000  0.106000  0.064000
0.022000  0.777000  0.064000  0.137000
```

Each listed row sums to exactly 1.0:

```text
1.000000, 1.000000, 1.000000, 1.000000
```

The columns sum to `0.130`, `2.783`, `0.717`, and `0.370`. Thus the conclusion
depends on the visual orientation of the live page:

- If the four listed vectors are **motif positions**, with values ordered
  A/C/G/T, they are ordinary probability rows.
- If they are four **base rows** across positions, then they are not a standard
  position-normalized MEME probability display.

The copied numbers alone do not establish that the webserver uses a scoring or
weight matrix. The live page layout or downloadable PWM/MEME file is needed to
resolve orientation and representation.

## CTCF/P49711 findings

### TF record

| Field | Local database value |
| --- | --- |
| TF accession | `P49711` |
| UniProt entry | `CTCF_HUMAN` |
| Gene names | `CTCF` |
| Protein | Transcriptional repressor CTCF (11-zinc finger protein) |
| Reviewed | Yes |
| Family | `More than 3 adjacent zinc finger factors{2.3.3}` |
| Primary annotation | `Identical_PWM` |
| Primary motif/model | `CTCF_HUMAN.H11MO.0.A` |

### Local motif and model evidence

P49711 has exactly one local `motif_ref` record:

| Source | Motif | Evidence | Mapping type | Consensus | Width | Matrix status |
| --- | --- | --- | --- | --- | --- | --- |
| HOCOMOCO | `CTCF_HUMAN.H11MO.0.A` | identical | `public_database_mapping_unconfirmed` | `TGGCCACCAGGGGGCGCCA` | 19 | usable |

Its stored row-sum range is effectively 1.0 (`0.9999999999999999` to `1.0`).
The HOCOMOCO raw MEME file begins with a normal letter-probability matrix:

```text
letter-probability matrix: alength= 4 w= 19 nsites= 500
0.094  0.360  0.126  0.420
0.152  0.168  0.522  0.158
0.258  0.084  0.554  0.104
0.068  0.856  0.036  0.040
```

There are no P49711 records in `motif_ref` with source `modcre` or
`alphafold`, and no `structure_file` records for P49711. Consequently, there
is no local CTCF ModCRE model, no linked ModCRE PWM, and no local member path
to print for the live-webserver test case.

### Search for the webserver model identifiers

Searches of SQLite, `pwms.tar.gz`, `models.tar.gz`, and `extracted/` found no
member containing the combination P49711/CTCF with `2wbu` and `299:373`.

The supplied local archive does contain three active ModCRE models using
template `2wbu`, but they are different TF records and regions:

| TF | Model ID | Region |
| --- | --- | --- |
| `Q3SY56` | `TFS_sp_Q3SY56_SP6_HUMAN:254:336_2wbu_A_1` | 254:336 |
| `Q96IQ9` | `TFS_sp_Q96IQ9_ZN414_HUMAN:111:195_2wbu_A_1` | 111:195 |
| `B3KXP2` | `TFS_tr_B3KXP2_B3KXP2_HUMAN:254:336_2wbu_A_1` | 254:336 |

Template identity alone is therefore not enough to identify the same TF model
or the same predicted PWM.

## Screenshot-value search

The six individual numbers (`0.106000`, `0.316000`, `0.504000`, `0.861000`,
`0.829000`, and `0.777000`) occur separately in many unrelated local motif
files. That is expected for three-decimal probability values and is not
evidence of a match.

More specific searches were negative:

| Search | Result |
| --- | --- |
| Complete four-row screenshot block in any SQLite `matrix_json` | 0 matches |
| Same block after a 4-by-4 transpose | 0 matches |
| Complete listed sequence in `motif_file.content` | 0 matches |
| Complete transposed sequence in `motif_file.content` | 0 matches |
| Complete listed sequence in extracted `*.meme` files | 0 matches |
| Complete transposed sequence in extracted `*.meme` files | 0 matches |

There is no matching local CTCF model output file to search for a separate
raw PWM/weight-matrix representation.

## Verified local ModCRE control: A0A024R0Y4

Tested motif/model:

```text
source:   modcre
motif:    TFS_tr_A0A024R0Y4_A0A024R0Y4_HUMAN:211:243_6r2v_A_1
model:    TFS_tr_A0A024R0Y4_A0A024R0Y4_HUMAN:211:243_6r2v_A_1
template: 6r2v (chain A in the model identifier)
region:   211:243
```

The motif file is present at:

```text
pwms_modcre_all/TFS_tr_A0A024R0Y4_A0A024R0Y4_HUMAN:211:243_6r2v_A_1.meme
```

The matching active PDB is present at:

```text
models/TFS_tr_A0A024R0Y4_A0A024R0Y4_HUMAN:211:243_6r2v_A_1.pdb
```

An identically named failed/debug PDB also exists under `models_failed/`. That
duplicate status is archive provenance and does not alter the PWM comparison.

The raw MEME header is:

```text
letter-probability matrix: alength= 4 w= 10 nsites= 20 E= 0
```

Raw MEME rows, in A/C/G/T order, are:

| Position | A | C | G | T | Sum |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 0.221 | 0.152 | 0.346 | 0.281 | 1.000 |
| 2 | 0.368 | 0.009 | 0.286 | 0.338 | 1.001 |
| 3 | 0.411 | 0.004 | 0.004 | 0.580 | 0.999 |
| 4 | 0.017 | 0.004 | 0.004 | 0.974 | 0.999 |
| 5 | 0.160 | 0.069 | 0.043 | 0.727 | 0.999 |
| 6 | 0.385 | 0.338 | 0.022 | 0.255 | 1.000 |
| 7 | 0.065 | 0.563 | 0.004 | 0.368 | 1.000 |
| 8 | 0.004 | 0.983 | 0.009 | 0.004 | 1.000 |
| 9 | 0.004 | 0.957 | 0.004 | 0.035 | 1.000 |
| 10 | 0.009 | 0.545 | 0.074 | 0.372 | 1.000 |

Checks performed:

| Check | Result |
| --- | --- |
| Raw archive MEME content equals SQLite `motif_file.content` | Yes, byte-for-byte |
| Parsed raw MEME numeric rows equal SQLite `matrix_json` | Yes |
| Matrix status | usable |
| Stored row-sum range | 0.999 to 1.001 |
| Additional `letter-probability matrix` sections in MEME file | No; exactly one |
| `weight`, `log-odds`, or `raw PWM` text in MEME file | No |
| Additional local raw PWM/weight-matrix file for this model | None found; the related archive entries are active and failed PDBs |

This is strong evidence that ModCREDB displays the supplied local ModCRE MEME
probability matrix correctly for a model that is actually present in the local
dataset.

## Classification

| Question | Classification |
| --- | --- |
| Is the local ModCRE MEME-to-SQLite import broken? | No evidence of an importer bug. The verified local control matches exactly. |
| Is the live CTCF/P49711 webserver matrix the same local record? | No. The corresponding CTCF/2wbu/299:373 record is absent locally. |
| Is the live matrix a transpose or reverse complement of a local matrix? | No exact local block match in either tested orientation; reverse-complement testing is not meaningful without a candidate local record. |
| Is the live page definitely showing a scoring/weight matrix? | Unresolved. The copied rows each sum to 1.0; display orientation must be verified from the live output. |
| Best current explanation | Different model/run/version/settings or a record not included in the supplied precomputed dataset. |

## Recommended next action

Do not change ModCREDB matrix display yet. Obtain one of the following from the
same live CTCF webserver experiment:

1. its downloadable MEME/PWM file, or
2. the full plain-text matrix with the A/C/G/T and position headers intact,
   plus its job identifier and tool/version/settings.

Then compare that raw downloaded file with a local record only after confirming
the same TF accession, model identifier, template chain, residue interval, and
generation workflow. If the live job is supposed to have been part of the
supplied ModCRElib dataset, its absence should be raised as a data-provenance
question rather than repaired by transforming the UI matrix.

## Repository state

Only this audit report was added for this task. No code, database, schema,
importer, raw archive, or scanner file was modified. No commit was made.
