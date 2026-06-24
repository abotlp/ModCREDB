# HOCOMOCO Official Mapping Audit

## Scope

This is a read-only comparison of the promoted ModCREDB staging database, the local official HOCOMOCO v11 human CORE annotation table, the matching MEME motif file, and the final integrated hierarchical chart. No SQLite rows were added, changed, or removed.

## Summary

- Official HOCOMOCO annotation rows: 401
- Official HOCOMOCO MEME motifs: 401
- ModCREDB HOCOMOCO links audited: 401
- Exact official UniProt-accession confirmations: 401
- Non-confirming/mismatch rows: 0

## Matrix Status

| Matrix status | Links |
| --- | ---: |
| usable | 401 |

## Mapping Results

| Result | Links |
| --- | ---: |
| exact_uniprot_confirmed | 401 |

`exact_uniprot_confirmed` means the official HOCOMOCO `UniProt AC` equals the ModCREDB TF accession. The chart-column and evidence-tier checks are also included in this result.

## Controls

| DB TF | HOCOMOCO model | Official UniProt AC | Chart column | Result |
| --- | --- | --- | --- | --- |
| P04637 | P53_HUMAN.H11MO.0.A | P04637 | Identical_PWM | exact_uniprot_confirmed |
| P37231 | PPARG_HUMAN.H11MO.0.A | P37231 | Identical_PWM | exact_uniprot_confirmed |

## Interpretation

All audited HOCOMOCO links are confirmed by the official HOCOMOCO UniProt accession, matching local MEME matrix, final-chart evidence column, and stored evidence tier. No database repair is needed.
