# CIS-BP Metadata Link Preflight

## Scope

This is a read-only preflight for potential CIS-BP TF-to-motif metadata links.
It does not modify SQLite, importer code, raw data, or any evidence links.

## Local Archive

The supplied archive is labelled as a CIS-BP 2019 PWM collection. Its 10,035
members all follow this pattern:

```text
pbm/CisBP_2019/pwms/<motif_id>.meme
```

Examples:

```text
M00001_2.00.meme
M04524_2.00.meme
```

The archive contains MEME motif matrices only. It contains no TSV, CSV, JSON,
SQL, or other metadata member that maps a motif ID to a protein, gene, taxon,
or UniProt accession.

The official CIS-BP update log identifies **Build 2.00** as the release used
for the Lambert et al. 2019 follow-up publication. Together with the local
`CisBP_2019` directory label and the `_2.00` motif IDs, this identifies the
local archive as a CIS-BP Build-2.00/2019 snapshot. The individual legacy PWM
files are still reachable under the official `data/2_00` path.

## Current ModCREDB Coverage

The validated staging database contains:

| Measure | Count |
| --- | ---: |
| CIS-BP motif files | 10,035 |
| Usable DNA matrices | 10,035 |
| Motif files with at least one existing TF evidence link | 1,315 |
| Motif files with zero current TF evidence links | 8,720 |
| Existing CIS-BP evidence-link rows | 78,959 |
| Direct/identical evidence rows | 762 |
| Close-homolog evidence rows | 3,706 |
| Distant-homolog evidence rows | 74,491 |

The 8,720 zero-link motif files are **not automatically missing links**. The
collection includes motifs from many species and proteins outside the supplied
human TF record set. Without a version-matched metadata table, a zero-link
motif cannot safely be assigned to a TF by motif name or consensus sequence.

## Official CIS-BP Site

The current CIS-BP bulk-download page exposes a current CIS-BP 3.10 SQL dump:

```text
https://cisbp.ccbr.utoronto.ca/data/3_10/DataFiles/SQLDumps/SQLArchive_cisbp_3_10.zip
```

It is useful as a future reference, but it is not confirmed to match the local
2019 PWM archive. Applying current CIS-BP 3.10 metadata to the older local
motif snapshot could introduce mismatched or historical motif-ID mappings.

## Decision

Do **not** build or apply a CIS-BP direct-link importer yet.

The best next input is a version-matched CIS-BP Build-2.00/2019 metadata
export, ideally the source `TF_Information` table/file that was downloaded
with the supplied PWM archive. It must contain at least:

```text
motif ID -> protein or UniProt accession -> species/taxon
```

Once that file is available, use the same safe process as JASPAR:

1. write a dry-run metadata-link audit;
2. inspect known positive controls and proposed counts;
3. apply only into a new candidate SQLite database;
4. validate before promotion.

## Alternative Retrieval Routes

The missing table may be recoverable without asking Baldo for it. Safe routes,
in order, are:

1. Retrieve the historical Build-2.00 `TF_Information` export from the CIS-BP
   maintainers or an archived CIS-BP release.
2. Locate a trusted academic mirror that explicitly preserves the Build-2.00
   `motif ID -> protein/UniProt` table and verify a checksum and sample rows
   against the local motif IDs.
3. Use current CIS-BP Build 3.10 only for a **read-only cross-release audit**:
   identify exact unchanged `Mxxxxx_2.00` IDs, then manually review any
   proposed UniProt mappings before creating an importer. Build 3.10 must not
   be applied automatically to the Build-2.00 archive.

Do not infer a direct link from a motif name, consensus, TF family, or a
similarity search. Those methods can help investigate candidates, but they do
not establish a reproducible direct CIS-BP mapping.

## Build-3.10 Cross-Release Result

This fallback was tested with the official CIS-BP Build-3.10
`TF_Information_all_motifs` download. The current file contains **zero** exact
matches to the 10,035 local `Mxxxxx_2.00` motif IDs. CIS-BP Build 3.10 has
renumbered or replaced the relevant motif records, so it cannot provide an
exact-ID crosswalk for this local Build-2.00 archive.

Therefore the only defensible paths are recovery of Build-2.00 metadata or
leaving the existing chart-derived CIS-BP links unchanged.

## Why This Is Different From JASPAR

JASPAR was safe to repair because we confirmed the local archive was JASPAR
2024 and obtained the matching official JASPAR 2024 SQL metadata dump. For
CIS-BP, the local archive and currently downloadable metadata are visibly from
different releases. Preserving the current supplied CIS-BP evidence is safer
than creating unverified direct links.
