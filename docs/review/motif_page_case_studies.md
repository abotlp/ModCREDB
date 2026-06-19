# Motif Page Case Studies

These examples cover usable and no-matrix motifs from the main motif/model sources. HOCOMOCO is represented through `mapping_type=public_database_mapping_unconfirmed` and `curation_status=pending_confirmation` until exact TF mapping semantics are confirmed.

| case | tf_id | source | evidence_type | mapping_type | curation_status | motif_id | width | nsites | consensus | matrix_status | evidence_note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| JASPAR example | A0A024R2I8 | jaspar | identical | direct_or_identical | pending_confirmation | MA1574.1 | 15 | 8627 | GGGGTCAAAGGTCAT | usable | Exact meaning of Identical_PWM needs PI/Baldo confirmation. |
| CIS-BP example | A0A024R3C6 | cisbp | relative_homologous | distant_homolog | pending_confirmation | M00143_2.00 | 10 | 20 | CCCCCCCACG | usable | Relative homology threshold/method needs PI/Baldo confirmation. |
| HOCOMOCO example | A0PJY2 | hocomoco | identical | public_database_mapping_unconfirmed | pending_confirmation | FEZF1_HUMAN.H11MO.0.C | 12 | 500 | GCTGCTCTTTTT | usable | HOCOMOCO v11 human CORE mononucleotide motif imported as public motif evidence; exact TF mapping semantics require PI/Baldo confirmation. |
| ModCRE usable-matrix example | A0A024R0Y4 | modcre | modcre | structure_predicted | imported | TFS_tr_A0A024R0Y4_A0A024R0Y4_HUMAN:211:243_6r2v_A_1 | 10 | 20 | GATTTACCCC | usable | Structure-derived predicted PWM/model evidence from provided ModCRE/ModCRElib-related dataset; exact generation settings pending confirmation. |
| ModCRE no-matrix/w=0 example | A0A1W2PPK0 | modcre | modcre | structure_predicted | imported | DIMER_sp_A0A1W2PPK0_CPHL2_HUMAN:22:78_1puf_B_1 | 0 | 20 |  | width_zero_no_matrix | Structure-derived predicted PWM/model evidence from provided ModCRE/ModCRElib-related dataset; exact generation settings pending confirmation. |
| AF3 usable-matrix example | A0A024R0Y4 | alphafold | alphafold | af3_structure_predicted | imported | a0a024r0y4.af3.0 | 14 | 20 | TAATCACGGCGCCC | usable | AF3-derived predicted PWM/model evidence from provided dataset; exact generation workflow pending confirmation. |
| AF3 no-matrix/w=0 example | A0A087WX29 | alphafold | alphafold | af3_structure_predicted | imported | a0a087wx29.af3.0 | 0 | 20 |  | width_zero_no_matrix | AF3-derived predicted PWM/model evidence from provided dataset; exact generation workflow pending confirmation. |
| ModCRE no-parsed-matrix example | Q92966 | modcre | modcre | structure_predicted | imported | TFS_sp_Q92966_SNPC3_HUMAN:21:411_8iue_3_1 | 2 | 20 |  | no_parsed_matrix | Structure-derived predicted PWM/model evidence from provided ModCRE/ModCRElib-related dataset; exact generation settings pending confirmation. |
