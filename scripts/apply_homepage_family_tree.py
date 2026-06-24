#!/usr/bin/env python3
"""Install the tree-first homepage template in a local ModCREDB checkout."""

from __future__ import annotations

import argparse
import shutil
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
INDEX = ROOT / "templates" / "index.html"
NEW_INDEX = ROOT / "templates" / "index_family_tree.html"
BASE = ROOT / "templates" / "base.html"

OLD_NAV = '''      <nav>
        <a href="/search">Search</a>
        <a href="/scan">Scan DNA</a>
        <a href="/docs">Docs</a>
        <a href="/evidence">Evidence</a>
        {% if debug_enabled %}<a href="/debug">Debug</a>{% endif %}
      </nav>'''

NEW_NAV = '''      <nav>
        <a href="/#browse">Browse</a>
        <a href="/search">Search</a>
        <a href="/scan">Scan</a>
        <a href="/docs">Docs</a>
        <a href="/evidence">Evidence</a>
        {% if debug_enabled %}<a href="/debug">Debug</a>{% endif %}
      </nav>'''


def backup(path: Path) -> Path:
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out = path.with_suffix(path.suffix + f".family_tree_backup_{stamp}")
    shutil.copy2(path, out)
    return out


def apply_index(dry_run: bool) -> None:
    if not NEW_INDEX.exists():
        raise SystemExit(f"Missing new homepage template: {NEW_INDEX}")
    if dry_run:
        print(f"Would copy {NEW_INDEX} -> {INDEX}")
        return
    backup_path = backup(INDEX)
    shutil.copy2(NEW_INDEX, INDEX)
    print(f"Backed up {INDEX} to {backup_path}")
    print(f"Installed {NEW_INDEX} as {INDEX}")


def apply_nav(dry_run: bool) -> None:
    text = BASE.read_text(encoding="utf-8")
    if NEW_NAV in text:
        print("Navigation already contains Browse link.")
        return
    if OLD_NAV not in text:
        print("Navigation block did not match expected old template; leaving base.html unchanged.")
        return
    if dry_run:
        print("Would update templates/base.html navigation.")
        return
    backup_path = backup(BASE)
    BASE.write_text(text.replace(OLD_NAV, NEW_NAV), encoding="utf-8")
    print(f"Backed up {BASE} to {backup_path}")
    print("Updated navigation in templates/base.html")


def main() -> None:
    parser = argparse.ArgumentParser(description="Install the ModCREDB tree-first homepage template.")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--skip-nav", action="store_true")
    args = parser.parse_args()
    apply_index(args.dry_run)
    if not args.skip_nav:
        apply_nav(args.dry_run)
    print("Next:")
    print("  python3 scripts/build_family_tree_data.py --db data/tf_webdb.sqlite")
    print("  python3 app.py --db data/tf_webdb.sqlite --port 8198")


if __name__ == "__main__":
    main()
