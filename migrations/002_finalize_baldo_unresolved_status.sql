-- Add the finalized public status for the 108 curated Baldo accessions.
-- Apply only to a verified staging copy; data/tf_webdb.sqlite is not a target.

-- migrate:up
PRAGMA foreign_keys = ON;
BEGIN IMMEDIATE;

CREATE TABLE tf_structure_status__new (
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
            'UNCERTAIN_MANUAL_REVIEW',
            'PWM_PRESENT_NO_VALIDATED_DNA_BOUND_STRUCTURE'
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

INSERT INTO tf_structure_status__new (
    tf_id,
    uniprot_accession,
    primary_structural_status,
    own_accession_structure_status,
    canonical_reference_status,
    canonical_reference_accession,
    canonical_reference_model_path,
    database_display_recommendation,
    action_for_baldo,
    decision_reason,
    remaining_uncertainty,
    review_status,
    source_audit_file,
    created_at,
    updated_at
)
SELECT
    tf_id,
    uniprot_accession,
    primary_structural_status,
    own_accession_structure_status,
    canonical_reference_status,
    canonical_reference_accession,
    canonical_reference_model_path,
    database_display_recommendation,
    action_for_baldo,
    decision_reason,
    remaining_uncertainty,
    review_status,
    source_audit_file,
    created_at,
    updated_at
FROM tf_structure_status;

DROP TABLE tf_structure_status;
ALTER TABLE tf_structure_status__new RENAME TO tf_structure_status;

CREATE INDEX idx_tf_structure_status_primary
    ON tf_structure_status(primary_structural_status);
CREATE INDEX idx_tf_structure_status_review
    ON tf_structure_status(review_status);
CREATE INDEX idx_tf_structure_status_canonical
    ON tf_structure_status(canonical_reference_accession);

COMMIT;

-- migrate:down
PRAGMA foreign_keys = ON;
BEGIN IMMEDIATE;

CREATE TABLE tf_structure_status__old (
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

INSERT INTO tf_structure_status__old (
    tf_id,
    uniprot_accession,
    primary_structural_status,
    own_accession_structure_status,
    canonical_reference_status,
    canonical_reference_accession,
    canonical_reference_model_path,
    database_display_recommendation,
    action_for_baldo,
    decision_reason,
    remaining_uncertainty,
    review_status,
    source_audit_file,
    created_at,
    updated_at
)
SELECT
    tf_id,
    uniprot_accession,
    primary_structural_status,
    own_accession_structure_status,
    canonical_reference_status,
    canonical_reference_accession,
    canonical_reference_model_path,
    database_display_recommendation,
    action_for_baldo,
    decision_reason,
    remaining_uncertainty,
    review_status,
    source_audit_file,
    created_at,
    updated_at
FROM tf_structure_status;

DROP TABLE tf_structure_status;
ALTER TABLE tf_structure_status__old RENAME TO tf_structure_status;

CREATE INDEX idx_tf_structure_status_primary
    ON tf_structure_status(primary_structural_status);
CREATE INDEX idx_tf_structure_status_review
    ON tf_structure_status(review_status);
CREATE INDEX idx_tf_structure_status_canonical
    ON tf_structure_status(canonical_reference_accession);

COMMIT;
