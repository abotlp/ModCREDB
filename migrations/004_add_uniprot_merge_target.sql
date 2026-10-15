-- Record the current UniProt accession for historical entries that were merged.
-- Apply only through the controlled importer after validating all target TF rows.

-- migrate:up
PRAGMA foreign_keys = ON;
BEGIN IMMEDIATE;

ALTER TABLE tf_annotation ADD COLUMN merged_into_accession TEXT
    REFERENCES tf(tf_id) ON DELETE RESTRICT;

CREATE INDEX idx_tf_annotation_merged_into
    ON tf_annotation(merged_into_accession);

COMMIT;

-- migrate:down
PRAGMA foreign_keys = ON;
BEGIN IMMEDIATE;

DROP INDEX IF EXISTS idx_tf_annotation_merged_into;
ALTER TABLE tf_annotation DROP COLUMN merged_into_accession;

COMMIT;
