CREATE TABLE import_issue (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    severity TEXT NOT NULL,
    category TEXT NOT NULL,
    message TEXT NOT NULL,
    tf_id TEXT,
    source TEXT,
    motif_id TEXT
);

CREATE TABLE metadata (
    key TEXT PRIMARY KEY,
    value TEXT NOT NULL
);

CREATE TABLE model_summary (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    summary_file_id INTEGER NOT NULL REFERENCES structure_file(id) ON DELETE CASCADE,
    matched_structure_id INTEGER REFERENCES structure_file(id) ON DELETE SET NULL,
    source TEXT NOT NULL REFERENCES source(source),
    status TEXT NOT NULL,
    tf_id TEXT,
    summary_model_id TEXT NOT NULL,
    model_rank INTEGER,
    n TEXT,
    template_pdb TEXT,
    n_tails TEXT,
    c_tails TEXT,
    protein_chain TEXT,
    dna_chain TEXT,
    identities TEXT,
    coverage TEXT,
    template_by_rmsd TEXT,
    domain TEXT,
    identity_percent REAL,
    similarity_percent REAL
);

CREATE TABLE motif_file (
    source TEXT NOT NULL REFERENCES source(source),
    motif_id TEXT NOT NULL,
    member_path TEXT NOT NULL,
    archive_path TEXT NOT NULL,
    content TEXT NOT NULL,
    width INTEGER,
    nsites TEXT,
    consensus TEXT,
    matrix_json TEXT,
    PRIMARY KEY (source, motif_id)
);

CREATE TABLE motif_ref (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    tf_id TEXT NOT NULL REFERENCES tf(tf_id) ON DELETE CASCADE,
    evidence_type TEXT NOT NULL,
    source TEXT NOT NULL REFERENCES source(source),
    motif_id TEXT NOT NULL,
    original_value TEXT NOT NULL,
    identity_percent REAL,
    missing_local_file INTEGER NOT NULL DEFAULT 0
);

CREATE TABLE motif_structure (
    motif_ref_id INTEGER NOT NULL REFERENCES motif_ref(id) ON DELETE CASCADE,
    structure_file_id INTEGER NOT NULL REFERENCES structure_file(id) ON DELETE CASCADE,
    PRIMARY KEY (motif_ref_id, structure_file_id)
);

CREATE TABLE source (
    source TEXT PRIMARY KEY,
    label TEXT NOT NULL,
    description TEXT NOT NULL
);

CREATE TABLE structure_file (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    source TEXT NOT NULL REFERENCES source(source),
    model_id TEXT NOT NULL,
    tf_id TEXT,
    member_path TEXT NOT NULL,
    archive_path TEXT NOT NULL,
    file_type TEXT NOT NULL,
    status TEXT NOT NULL,
    template_pdb TEXT,
    residue_start INTEGER,
    residue_end INTEGER
);

CREATE TABLE tf (
    tf_id TEXT PRIMARY KEY,
    family_text TEXT NOT NULL,
    motif_ref_count INTEGER NOT NULL DEFAULT 0,
    active_model_count INTEGER NOT NULL DEFAULT 0
);

CREATE TABLE tf_annotation (
    tf_id TEXT PRIMARY KEY REFERENCES tf(tf_id) ON DELETE CASCADE,
    uniprot_accession TEXT NOT NULL,
    entry_name TEXT,
    gene_names TEXT,
    protein_name TEXT,
    organism_name TEXT,
    organism_id INTEGER,
    reviewed INTEGER,
    sequence_length INTEGER,
    annotation_score REAL,
    uniprot_url TEXT,
    fetched_at TEXT NOT NULL
);

CREATE TABLE tf_family (
    tf_id TEXT NOT NULL REFERENCES tf(tf_id) ON DELETE CASCADE,
    family TEXT NOT NULL,
    PRIMARY KEY (tf_id, family)
);

CREATE TABLE tf_primary_annotation (
            tf_id TEXT PRIMARY KEY REFERENCES tf(tf_id) ON DELETE CASCADE,
            best_annotation_level TEXT,
            best_pwm_or_model TEXT,
            n_nonempty_annotation_columns INTEGER,
            source_table TEXT NOT NULL
        );

CREATE INDEX idx_model_summary_matched_structure ON model_summary(matched_structure_id);

CREATE INDEX idx_model_summary_template ON model_summary(template_pdb);

CREATE INDEX idx_model_summary_tf ON model_summary(tf_id);

CREATE INDEX idx_motif_ref_evidence ON motif_ref(evidence_type);

CREATE INDEX idx_motif_ref_motif ON motif_ref(motif_id);

CREATE INDEX idx_motif_ref_tf ON motif_ref(tf_id);

CREATE INDEX idx_structure_file_model ON structure_file(source, model_id);

CREATE INDEX idx_structure_file_status ON structure_file(status);

CREATE INDEX idx_structure_file_tf ON structure_file(tf_id);

CREATE INDEX idx_tf_annotation_gene_names ON tf_annotation(gene_names);

CREATE INDEX idx_tf_annotation_organism_name ON tf_annotation(organism_name);

CREATE INDEX idx_tf_family_family ON tf_family(family);
