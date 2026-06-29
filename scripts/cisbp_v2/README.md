# CIS-BP v2 direct motif assignments

This workflow fixes missing direct CIS-BP v2 TF-to-motif assignments in `motif_ref`.

Issue:
- CIS-BP PWM files existed in `motif_file`.
- Some TF pages lacked direct CIS-BP motifs because `motif_ref` was not populated from CIS-BP v2 `TF_Information.txt`.
- Example: TP53/P04637 had JASPAR and HOCOMOCO motifs but lacked direct CIS-BP motifs from CIS-BP TF report `T311040_2.00`.

Versioned patch:
- `data/cisbp_v2_direct_human_motif_ref_insert_candidates.tsv`
- `scripts/cisbp_v2/apply_cisbp_v2_direct_assignments.py`

Insertion criteria:
- direct CIS-BP v2 mapping from `TF_Information.txt`
- `TF_Species == Homo_sapiens`
- `TF_Status == D`
- selected ModCREDB TF is a unique reviewed UniProt entry
- local CIS-BP PWM exists in `motif_file`
- exact `tf_id/source/motif_id` row does not already exist

Do not use `TF_Information_all_motifs_plus.txt` for direct/known rows; it includes related/inferred motifs.
