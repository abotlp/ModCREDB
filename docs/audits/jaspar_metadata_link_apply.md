# JASPAR 2024 Metadata Link Import

Mode: `apply`.
Metadata format: JASPAR SQL dump.
Metadata source: https://jaspar2024.elixir.no/download/database/JASPAR2024.sql.gz.
Metadata SHA256: `33f00db0f06adaa5fff497827a0f63551da4c30e5c6dd601c08ac625cecc60c4`.
Result: a separate candidate database.

## Candidate Counts

- Metadata motif-UniProt pairs: 5292
- Existing direct links: 531
- Direct links added/proposed: 1087
- Missing TF records: 2602
- Missing local motif files: 1072
- Metadata rows without UniProt IDs: 299
- Applied rows: 1087

## Semantics

Inserted rows preserve existing evidence and use `evidence_type=identical`, `mapping_type=direct_or_identical`, `original_column=JASPAR2024_metadata`, and `curation_status=pending_confirmation`.

A pre-existing homology link does not block a direct JASPAR metadata link; both evidence records are retained. No missing TFs or missing motif files are created.

## Controls

- `MA0106.3` / `P04637` is expected to receive a direct JASPAR metadata link.
- `MA0065.1` / `P37231` is expected to remain an already-present direct link.
