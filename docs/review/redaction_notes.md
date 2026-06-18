# Redaction Notes

The review bundle is safe to commit publicly.

Redactions applied:

- Full private SQLite databases are not included.
- Raw private archives are not included.
- Absolute local filesystem paths are excluded from TSV/Markdown reports.
- In the sample SQLite DB, `archive_path`, metadata path values, and hierarchy source-table paths are replaced with redacted placeholders.
- The sample DB keeps motif MEME text and matrices because these are small and needed for schema/evidence review.
- PDB contents are not embedded in the sample DB.
- HOCOMOCO links are labelled in review outputs as `public_database_motif_unconfirmed` until exact mapping semantics are confirmed.
