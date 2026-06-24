# Homepage family-tree implementation report

This clean branch adds a tree-first homepage implementation without patch-helper scripts.

## Added files

- `templates/index_family_tree.html`
- `static/family_tree.js`
- `static/home_family_tree.css`
- `data_sources/tf_family_tree.json`
- `scripts/build_family_tree_data.py`
- `scripts/apply_homepage_family_tree.py`

## Important difference from the previous working branch

This branch does not include `scripts/fix_family_tree_builder_schema.py` or `docs/audits/homepage_tree_schema_fix_note.md`. The schema fix is implemented directly in `scripts/build_family_tree_data.py`.

## Local application

```bash
cd /home/patricia/tf_webdb
git fetch origin homepage-tree-clean
git switch -c homepage-tree-clean --track origin/homepage-tree-clean
python3 -m py_compile scripts/build_family_tree_data.py scripts/apply_homepage_family_tree.py
python3 scripts/build_family_tree_data.py --db data/tf_webdb.sqlite
python3 scripts/apply_homepage_family_tree.py
python3 app.py --db data/tf_webdb.sqlite --port 8198
```

## Builder schema safety

`build_family_tree_data.py` now inspects `PRAGMA table_info(structure_file)` and does not assume a column named `path`. It chooses one of these path-like columns when available:

- `path`
- `file_path`
- `local_path`
- `model_path`
- `pdb_path`
- `structure_path`
- `source_path`
- `relative_path`
- `filename`

If no path-like column exists, it still counts TFs with active models but leaves monomer/dimer counts as zero.

## Next step

After visual review, add true `/search?family_id=<id>` handling in `app.py`. Current Open-family links use `/search?q=<family token>`.
