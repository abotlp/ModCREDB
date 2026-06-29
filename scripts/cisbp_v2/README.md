# CIS-BP v2 direct Homo sapiens motif assignment

Problem fixed:
- CIS-BP v2 PWM files were present in `motif_file`.
- Direct CIS-BP TF-to-motif assignments were missing from `motif_ref` for many reviewed human UniProt TFs.
- Example: TP53/P04637 had JASPAR and HOCOMOCO motifs, but lacked direct CIS-BP motifs from CIS-BP TF report `T311040_2.00`.

Primary source:
- CIS-BP v2 Homo sapiens bulk download.
- Required table: `TF_Information.txt`.

Important parsing rule:
- `TF_Information.txt` contains duplicate `DBID` column names.
- Parse by column index, not by `csv.DictReader`.
- Column 5 is TF external ID, e.g. Ensembl gene ID.
- Column 13 is motif/source identifier, e.g. HOCOMOCO/TRANSFAC source motif name.

Insertion filter:
- `TF_Species == Homo_sapiens`
- `TF_Status == D`
- selected ModCREDB TF is a unique reviewed UniProt entry
- local CIS-BP PWM exists in `motif_file`
- motif_ref row does not already exist

Do not use `TF_Information_all_motifs_plus.txt` for Known/direct rows. It includes related/inferred motifs.
