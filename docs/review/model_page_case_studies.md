# Model Page Case Studies

These examples are for reviewing how active, failed/debug, summary-linked, and chain-mismatch model records should be explained to users. Failed/debug examples should not be presented as public biological failure counts.

| case | tf_id | source | status | file_type | model_id | template_pdb | residues | linked_summary_rows | summary_model_id | summary_protein_chain | actual_protein_chains | summary_dna_chain | actual_dna_chains | interpretation |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| active ModCRE PDB model | A0A024R0Y4 | modcre | active | pdb | TFS_tr_A0A024R0Y4_A0A024R0Y4_HUMAN:211:243_6r2v_A_1 | 6r2v | 211-243 | 1 |  |  |  |  |  |  |
| active AF3 PDB model | A0A024R0Y4 | alphafold | active | pdb | a0a024r0y4.af3.0 |  |  | 0 |  |  |  |  |  |  |
| model summary linked to active PDB | A0A024R0Y4 | modcre | active | pdb | TFS_tr_A0A024R0Y4_A0A024R0Y4_HUMAN:211:243_6r2v_A_1 | 6r2v | 211-243 | 1 |  |  |  |  |  |  |
| model summary not linked to active PDB | A0A024R0Y4 | modcre | failed | pdb | TFS_tr_A0A024R0Y4_A0A024R0Y4_HUMAN:211:243_6r2v_A_1 | 6r2v | 211-243 | 1 |  |  |  |  |  |  |
| failed/debug duplicate example | A0A075B7G4 | alphafold | active:1; failed:1 | pdb | a0a075b7g4_dimer.af3.2 |  |  |  |  |  |  |  |  |  |
| chain mismatch example | Q14582 | modcre |  |  | TFS_sp_Q14582_MAD4_HUMAN:55:111_1an4_A_1 |  |  |  | TFS_sp_Q14582_MAD4_HUMAN | A,B | A | C,D | C,D | summary/template chain labels do not map directly onto active PDB chain IDs |
