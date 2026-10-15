BEGIN IMMEDIATE;

ALTER TABLE structure_file ADD COLUMN model_kind TEXT
    CHECK (model_kind IS NULL OR model_kind IN ('DIMER', 'TFS'));
ALTER TABLE structure_file ADD COLUMN protein_chain_count INTEGER
    CHECK (protein_chain_count IS NULL OR protein_chain_count >= 0);
ALTER TABLE structure_file ADD COLUMN dna_chain_count INTEGER
    CHECK (dna_chain_count IS NULL OR dna_chain_count >= 0);
ALTER TABLE structure_file ADD COLUMN dna_interaction INTEGER
    CHECK (dna_interaction IS NULL OR dna_interaction IN (0, 1));
ALTER TABLE structure_file ADD COLUMN dna_contact_count INTEGER
    CHECK (dna_contact_count IS NULL OR dna_contact_count >= 0);
ALTER TABLE structure_file ADD COLUMN dna_qc_cutoff REAL
    CHECK (dna_qc_cutoff IS NULL OR dna_qc_cutoff > 0);
ALTER TABLE structure_file ADD COLUMN dna_qc_status TEXT;
ALTER TABLE structure_file ADD COLUMN dna_qc_checked_at TEXT;

CREATE INDEX idx_structure_file_active_tf_quality
    ON structure_file(
        tf_id,
        status,
        file_type,
        source,
        dna_interaction,
        dna_chain_count,
        model_kind
    );

COMMIT;
