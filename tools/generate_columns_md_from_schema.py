#!/usr/bin/env python3
"""Generate a Markdown column reference from docs/schema.json (Task 11B).

This generator produces a *generated* documentation artifact so that column
metadata does not drift over time.

Output:
- docs/columns_schema.md (generated)

Example:
  python tools/generate_columns_md_from_schema.py
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Tuple


from molprop_toolkit.schema import load_schema


def md_escape(s: str) -> str:
    return (s or "").replace("|", "\\|")


def main() -> None:
    ap = argparse.ArgumentParser(description="Generate docs/columns_schema.md from docs/schema.json")
    ap.add_argument("--schema", default=None)
    ap.add_argument("--out", default="docs/columns_schema.md")
    args = ap.parse_args()

    schema = load_schema(args.schema)
    cols: Dict[str, Any] = schema.get("columns", {}) or {}
    patterns: List[Dict[str, Any]] = schema.get("patterns", []) or []

    # Group columns by first tag.
    by_tag: Dict[str, List[Tuple[str, Dict[str, Any]]]] = defaultdict(list)
    for name, meta in cols.items():
        tags = meta.get("tags", []) or ["(untagged)"]
        tag = str(tags[0])
        by_tag[tag].append((name, meta))

    for tag in by_tag:
        by_tag[tag].sort(key=lambda x: x[0])

    order = sorted(by_tag.keys())

    lines: List[str] = []
    lines.append("# Column reference (generated)")
    lines.append("")
    lines.append("This file is generated from `docs/schema.json`. Do not edit by hand; update the schema and re-run the generator.")
    lines.append("")
    lines.append(f"Generated: `{datetime.now().isoformat(timespec='seconds')}`")
    lines.append("")

    for tag in order:
        lines.append(f"## {tag}")
        lines.append("")
        lines.append("| column | type | units | required | description |")
        lines.append("| --- | --- | --- | --- | --- |")
        for name, meta in by_tag[tag]:
            t = str(meta.get("type", "auto"))
            u = meta.get("units")
            req = "yes" if meta.get("required") else "no"
            desc = str(meta.get("description", ""))
            lines.append(
                f"| `{md_escape(name)}` | `{md_escape(t)}` | `{md_escape(str(u) if u is not None else '')}` | {req} | {md_escape(desc)} |"
            )
        lines.append("")

    lines.append("## Pattern-based column families")
    lines.append("")
    lines.append("Patterns define metadata for whole families of columns, typically produced with shared prefixes.")
    lines.append("")
    lines.append("| regex | type | tags | description |")
    lines.append("| --- | --- | --- | --- |")
    for p in patterns:
        rx = str(p.get("match", ""))
        t = str(p.get("type", "auto"))
        tags = ", ".join(p.get("tags", []) or [])
        desc = str(p.get("description", ""))
        lines.append(f"| `{md_escape(rx)}` | `{md_escape(t)}` | {md_escape(tags)} | {md_escape(desc)} |")

    outp = Path(args.out)
    outp.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {outp}")


if __name__ == "__main__":
    main()

