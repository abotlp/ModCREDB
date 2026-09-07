-- ModCREDB local staging migration 001.
-- Apply only to a verified copy of the database, never to data/tf_webdb.sqlite.

-- migrate:up
PRAGMA foreign_keys = ON;
BEGIN IMMEDIATE;

CREATE TABLE tf_structure_status (
    tf_id TEXT PRIMARY KEY REFERENCES tf(tf_id) ON DELETE CASCADE,
    uniprot_accession TEXT NOT NULL UNIQUE,
    primary_structural_status TEXT NOT NULL CHECK (
        primary_structural_status IN (
            'VALID_DBD_FRAGMENT_MODEL',
            'DBD_PRESENT_FRAGMENT_FAILED',
            'DBD_PRESENT_FRAGMENT_NOT_SENT',
            'ACCESSION_LACKS_DBD',
            'ACCESSION_CONTAINS_INCOMPLETE_DBD',
            'COFACTOR_OR_NON_DNA_BINDING_COMPONENT',
            'UNCERTAIN_MANUAL_REVIEW'
        )
    ),
    own_accession_structure_status TEXT NOT NULL CHECK (
        own_accession_structure_status IN (
            'VALID_DBD_FRAGMENT_MODEL',
            'DBD_PRESENT_NO_VALID_FRAGMENT',
            'ACCESSION_LACKS_DBD',
            'ACCESSION_CONTAINS_INCOMPLETE_DBD',
            'NON_DNA_BINDING_COMPONENT',
            'UNRESOLVED'
        )
    ),
    canonical_reference_status TEXT NOT NULL CHECK (
        canonical_reference_status IN (
            'NOT_APPLICABLE',
            'SEPARATE_REFERENCE_ACCESSION_ONLY',
            'SEPARATE_REFERENCE_INTERFACE_VALID',
            'SEPARATE_REFERENCE_NOT_VALIDATED'
        )
    ),
    canonical_reference_accession TEXT REFERENCES tf(tf_id) ON DELETE RESTRICT,
    canonical_reference_model_path TEXT,
    database_display_recommendation TEXT NOT NULL,
    action_for_baldo TEXT NOT NULL,
    decision_reason TEXT NOT NULL,
    remaining_uncertainty TEXT NOT NULL DEFAULT '',
    review_status TEXT NOT NULL,
    source_audit_file TEXT NOT NULL,
    created_at TEXT NOT NULL,
    updated_at TEXT NOT NULL,
    CHECK (uniprot_accession = tf_id),
    CHECK (canonical_reference_accession IS NULL OR canonical_reference_accession <> tf_id),
    CHECK (
        (canonical_reference_status = 'NOT_APPLICABLE'
         AND canonical_reference_accession IS NULL
         AND canonical_reference_model_path IS NULL)
        OR
        (canonical_reference_status <> 'NOT_APPLICABLE'
         AND canonical_reference_accession IS NOT NULL)
    )
);

CREATE INDEX idx_tf_structure_status_primary
    ON tf_structure_status(primary_structural_status);
CREATE INDEX idx_tf_structure_status_review
    ON tf_structure_status(review_status);
CREATE INDEX idx_tf_structure_status_canonical
    ON tf_structure_status(canonical_reference_accession);

CREATE TABLE structure_model_assignment (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    tf_id TEXT NOT NULL REFERENCES tf(tf_id) ON DELETE CASCADE,
    display_accession TEXT NOT NULL REFERENCES tf(tf_id) ON DELETE CASCADE,
    model_accession TEXT NOT NULL REFERENCES tf(tf_id) ON DELETE RESTRICT,
    model_path TEXT NOT NULL,
    model_role TEXT NOT NULL CHECK (
        model_role IN (
            'OWN_ACCESSION_MODEL',
            'DBD_FRAGMENT_MODEL',
            'CANONICAL_SAME_GENE_REFERENCE'
        )
    ),
    interface_status TEXT,
    is_active INTEGER NOT NULL DEFAULT 0 CHECK (is_active IN (0, 1)),
    source TEXT NOT NULL,
    notes TEXT,
    created_at TEXT NOT NULL,
    updated_at TEXT NOT NULL,
    CHECK (tf_id = display_accession),
    CHECK (
        (model_role = 'CANONICAL_SAME_GENE_REFERENCE' AND display_accession <> model_accession)
        OR
        (model_role IN ('OWN_ACCESSION_MODEL', 'DBD_FRAGMENT_MODEL')
         AND display_accession = model_accession)
    ),
    UNIQUE (display_accession, model_path, model_role)
);

CREATE INDEX idx_structure_model_assignment_display
    ON structure_model_assignment(display_accession, is_active);
CREATE INDEX idx_structure_model_assignment_model
    ON structure_model_assignment(model_accession, model_role);

CREATE TABLE structure_confidence_artifact (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    model_path TEXT NOT NULL,
    artifact_type TEXT NOT NULL CHECK (
        artifact_type IN ('MMCIF', 'PLDDT_JSON', 'PAE_JSON', 'PLDDT_PLOT')
    ),
    artifact_path TEXT NOT NULL,
    sha256 TEXT NOT NULL CHECK (length(sha256) = 64),
    notes TEXT,
    created_at TEXT NOT NULL,
    UNIQUE (model_path, artifact_type, artifact_path)
);

CREATE INDEX idx_structure_confidence_artifact_model
    ON structure_confidence_artifact(model_path, artifact_type);

COMMIT;

-- migrate:down
PRAGMA foreign_keys = ON;
BEGIN IMMEDIATE;

DROP TABLE IF EXISTS structure_confidence_artifact;
DROP TABLE IF EXISTS structure_model_assignment;
DROP TABLE IF EXISTS tf_structure_status;

COMMIT;
