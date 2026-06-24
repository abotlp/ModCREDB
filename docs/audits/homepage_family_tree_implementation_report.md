# Homepage family-tree implementation report

## Summary

This branch adds a tree-first homepage implementation for ModCREDB. The design keeps the existing direct search/results workflow and replaces flat family cards with a HOCOMOCO-inspired, ModCREDB-specific DNA-binding family tree.

The tree is intended to be interactive:

- Click a family/branch node to update the selected-family panel on the right.
- Hover a node to show a tooltip with counts.
- Click **Open family** to go to filtered search results.
- Circle size reflects TF sequence-record count.
- Inner circle size reflects generated-PWM TF coverage when generated DB counts are available.

## Files added

- `templates/index_family_tree.html`
  - Replacement homepage template.
  - Keeps top search and prediction/source filters.
  - Adds the interactive family-tree section.
  - Keeps prediction/source/data-availability browse sections.

- `static/home_family_tree.css`
  - Homepage family-tree styles using the existing ModCREDB CSS variables.

- `static/family_tree.js`
  - Vanilla JavaScript renderer for the SVG family tree.
  - Reads generated DB-backed data from `/static/family_tree_data.generated.json` when present.
  - Falls back to embedded homepage JSON if generated data is absent.

- `data_sources/tf_family_tree.json`
  - HOCOMOCO-inspired display/navigation taxonomy.
  - Maps clean family labels to existing `tf_family.family` / `tf.family_text` tokens.
  - This is a display layer, not a schema migration.

- `scripts/build_family_tree_data.py`
  - Builds DB-backed `static/family_tree_data.generated.json` from SQLite.
  - Uses the taxonomy file and existing DB tables.

- `scripts/apply_homepage_family_tree.py`
  - Local helper to install the replacement homepage template safely.
  - Backs up `templates/index.html` and copies `templates/index_family_tree.html` into place.
  - Optionally patches the top navigation to add Browse and rename Scan DNA to Scan.

## Files intentionally not directly overwritten in this branch

The GitHub connector allowed adding files, but direct replacement of existing tracked templates was blocked. Therefore this branch includes an apply helper instead of directly replacing:

- `templates/index.html`
- `templates/base.html`

To apply locally:

```bash
cd /home/patricia/tf_webdb
python3 scripts/apply_homepage_family_tree.py
python3 scripts/build_family_tree_data.py --db data/tf_webdb.sqlite
```

Then start the app:

```bash
python3 app.py --db data/tf_webdb.sqlite --port 8198
```

## Tables used by `build_family_tree_data.py`

The builder reads:

- `tf`
- `tf_family`
- `tf_primary_annotation` if available
- `motif_ref`
- `motif_file`
- `structure_file`

It does not modify the database.

## Count logic

For each family node:

1. Find TF IDs matching taxonomy tokens.
   - Exact matches in `tf_family.family`.
   - Fallback substring matches in `tf.family_text` for legacy comma-separated family/domain strings.
2. If a node is a branch and has no direct tokens, use the union of all descendant-family TF IDs.
3. Count:
   - TF sequence records.
   - Known / Homologous / Relative homologous / ModCRE / AlphaFold primary levels.
   - Generated PWM TFs using `motif_ref` + `motif_file.matrix_status = 'usable'`.
   - TFs with active PDB models from `structure_file`.
   - Monomer/dimer model TFs from explicit path text containing `monomer` or `dimer`.

## Current limitations

- `/search?family_id=...` is not implemented yet because this first branch avoids editing `app.py` directly.
- Family links currently use `/search?q=<family token>` or generated `open_url` values.
- Some taxonomy tokens are approximate and should be reviewed against real `tf_family.family` values after the generated JSON is inspected.
- The selected-family ChimeraX figure is still a placeholder SVG. Real ChimeraX renders should be added later.
- Generated JSON should not be committed until we decide whether DB-derived static artifacts belong in version control.

## Suggested validation commands

```bash
cd /home/patricia/tf_webdb
python3 -m py_compile scripts/build_family_tree_data.py scripts/apply_homepage_family_tree.py
python3 scripts/apply_homepage_family_tree.py --dry-run
python3 scripts/build_family_tree_data.py --db data/tf_webdb.sqlite
python3 scripts/apply_homepage_family_tree.py
python3 app.py --db data/tf_webdb.sqlite --port 8198
```

Manual checks:

- Open `/` and verify the family tree appears.
- Click C2H2 and verify the right panel updates.
- Click a branch such as Zinc-coordinating domains and verify branch summary appears.
- Click Open family and verify search results open.
- Test `/search?q=TP53`.
- Test `/search?source=jaspar`.
- Test `/scan`.

## Next implementation step

After reviewing the homepage visually, the next code step should be extending `app.py` so `/search?family_id=<id>` uses the same taxonomy mapping as the homepage tree instead of relying on raw `q` searches.
