# TF Page Case Studies

These examples are selected from the staging database to review page organization and evidence preservation. HOCOMOCO is shown as a separate source; exact evidence semantics should be confirmed before public labelling.

| case | tf_id | gene | sources | evidence_types | usable_motifs | no_matrix_motifs | active_models | active_regions | family_excerpt |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| TF with JASPAR + CIS-BP + HOCOMOCO evidence | O00327 | BMAL1 | cisbp:3, hocomoco:1, jaspar:1 | homologous:3, public_database_motif_unconfirmed:1, relative_homologous:1 | 5 | 0 | 0 | 0 | HLH,PAS,PAS_2,PAS_3,PAS_9,PAS_4 |
| TF with direct/public motif only | A0A024R2I8 | NR1A2 | jaspar:1 | identical:1 | 1 | 0 | 0 | 0 | Hormone_recep,zf-C4 |
| TF with homologous motif only | A0A068JFF8 |  | jaspar:1 | homologous:1 | 1 | 0 | 0 | 0 |  |
| TF with relatively homologous motif evidence | A0A024R3C6 | ZBTB16 | alphafold:6, cisbp:4, jaspar:4, modcre:10 | alphafold:6, modcre:10, relative_homologous:8 | 23 | 1 | 16 | 8 | zf-C2H2_4,zf-C2H2,Endonuclease_7,Prok-RING_1,zf-C2H2_jaz,zf-met,zf-H2C2_2,DZR,BT |
| TF with ModCRE predicted PWM + active model | A0A1W2PPF3 | DUXB | alphafold:6, jaspar:1, modcre:19 | alphafold:6, modcre:19, relative_homologous:1 | 26 | 0 | 37 | 8 | CENP-B_N,MerR_1,HTH_3,Phage_CII,HTH_23,Homeobox,HTH_Tnp_1,HTH_26,Homeobox_KN |
| TF with AF3-derived predicted PWM + active model | A0A024R0Y4 | TADA2A | alphafold:6, modcre:3 | alphafold:6, modcre:3 | 9 | 0 | 9 | 3 | zf-C2H2_4,RPA_C,ZZ,Myb_DNA-bind_6,SWIRM,Myb_DNA-binding |
| TF with model but no usable PWM | A0A087WX29 | TARDBP | alphafold:6 | alphafold:6 | 0 | 6 | 6 | 0 | RRM_3,RRM_6,RRM_1 |
| TF with no-matrix/w=0 motif | A0A024R3C6 | ZBTB16 | alphafold:6, cisbp:4, jaspar:4, modcre:10 | alphafold:6, modcre:10, relative_homologous:8 | 23 | 1 | 16 | 8 | zf-C2H2_4,zf-C2H2,Endonuclease_7,Prok-RING_1,zf-C2H2_jaz,zf-met,zf-H2C2_2,DZR,BT |
| TF with multiple protein regions/domains | A0A024R0Y4 | TADA2A | alphafold:6, modcre:3 | alphafold:6, modcre:3 | 9 | 0 | 9 | 3 | zf-C2H2_4,RPA_C,ZZ,Myb_DNA-bind_6,SWIRM,Myb_DNA-binding |
| TF with active/failed duplicate model ID | A0A075B7G4 | ZNF595 | alphafold:6, cisbp:67, jaspar:5, modcre:5 | alphafold:6, modcre:5, relative_homologous:72 | 82 | 1 | 24 | 5 | zf-C2H2_4,C1_4,zf-C2H2,DNA_RNApol_7kD,C1_3,zf-C2H2_jaz,zf-trcl,Ogr_Delta,KRAB,XP |
| TF with summary-chain mismatch | Q14582 | MXD4 | cisbp:1, jaspar:1, modcre:12 | modcre:12, relative_homologous:2 | 14 | 0 | 13 | 5 | HLH |
