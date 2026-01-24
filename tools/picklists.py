#!/usr/bin/env python3
"""Picklist generator CLI.

Generates operational decision filters from a results CSV.

Outputs:
- one CSV per picklist
- optional HTML index + per-picklist HTML pages
- optional copying of structure sketches (so HTML is portable)

Examples:
  python tools/picklists.py results.csv --list-builtins
  python tools/picklists.py results.csv --outdir picklists/run1
  python tools/picklists.py results.csv --picklists top_cns_followup,oral_low_ddi_risk

With a user config (JSON/YAML):
  python tools/picklists.py results.csv --config my_picklists.json --html

With structure sketches:
  python tools/picklists.py results.csv --html --sketch-dir sketches_compare --structure-views calc
"""

from __future__ import annotations

import argparse
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

# Allow running directly without installing the package.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from molprop_toolkit.core import detect_id_column, read_csv
from molprop_toolkit.picklists.builtins import BUILTIN_PICKLISTS
from molprop_toolkit.picklists.engine import (
    PicklistDefinition,
    build_picklists,
    load_picklist_config,
    registry_category_columns,
    select_output_columns,
)
from molprop_toolkit.picklists.html import StructureAssetSpec, copy_structure_assets, write_index_html, write_picklist_html


def _parse_csv_list(s: Optional[str]) -> List[str]:
    if not s:
        return []
    return [x.strip() for x in s.split(",") if x.strip()]


def _resolve_definitions(
    *,
    picklist_keys: List[str],
    config_defs: List[PicklistDefinition],
) -> List[PicklistDefinition]:
    # Merge built-ins and config-based picklists into a single name space.
    all_defs: Dict[str, PicklistDefinition] = dict(BUILTIN_PICKLISTS)
    for d in config_defs:
        all_defs[d.key] = d

    if not picklist_keys:
        # Default: all built-ins + all config definitions.
        return list(all_defs.values())

    missing = [k for k in picklist_keys if k not in all_defs]
    if missing:
        raise SystemExit(f"Unknown picklists: {missing}. Use --list-builtins and/or check your config.")

    return [all_defs[k] for k in picklist_keys]


def main() -> None:
    ap = argparse.ArgumentParser(description="MolProp Toolkit picklist generator")
    ap.add_argument("input", help="Results CSV")
    ap.add_argument("--list-builtins", action="store_true", help="List built-in picklists and exit")
    ap.add_argument(
        "--picklists",
        help="Comma-separated picklist keys to run (default: all built-ins + config picklists)",
        default=None,
    )
    ap.add_argument("--config", help="Picklist config (JSON/YAML)", default=None)
    ap.add_argument(
        "--outdir",
        help="Output directory (default: picklists/<csvstem>_picklists_<timestamp>)",
        default=None,
    )
    ap.add_argument("--html", action="store_true", help="Write an index.html and per-picklist HTML pages")
    ap.add_argument(
        "--sketch-dir",
        help="Directory produced by tools/sketch_structures.py (used to copy structure assets into output)",
        default=None,
    )
    ap.add_argument(
        "--structure-views",
        help="Comma-separated views to embed when using --sketch-dir (compare mode: input,prepared,calc_base,protomer,calc). Default: calc",
        default="calc",
    )
    ap.add_argument(
        "--structure-formats",
        help="Comma-separated formats to copy/embed (svg,png). Default: svg",
        default="svg",
    )
    ap.add_argument(
        "--columns",
        choices=["minimal", "category", "all"],
        default="minimal",
        help="Which columns to include in each picklist CSV/HTML table",
    )

    args = ap.parse_args()

    if args.list_builtins:
        print("Built-in picklists:\n")
        for k, d in BUILTIN_PICKLISTS.items():
            logic = getattr(d, "logic", "all")
            print(f"- {k}: {d.name} (logic={logic}, filters={len(d.filters)})")
        return

    csv_path = args.input
    stem = Path(csv_path).stem
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    outdir = Path(args.outdir or (Path("picklists") / f"{stem}_picklists_{ts}"))
    outdir.mkdir(parents=True, exist_ok=True)

    config_defs: List[PicklistDefinition] = []
    if args.config:
        config_defs = load_picklist_config(args.config)

    picklist_keys = _parse_csv_list(args.picklists)
    defs = _resolve_definitions(picklist_keys=picklist_keys, config_defs=config_defs)

    df_full = read_csv(csv_path)
    id_col = detect_id_column(df_full)

    # Run picklists
    _, _, results = build_picklists(csv_path, defs)

    # Optional structure asset setup
    asset_spec: Optional[StructureAssetSpec] = None
    if args.sketch_dir:
        asset_spec = StructureAssetSpec(
            sketch_dir=Path(args.sketch_dir),
            views=_parse_csv_list(args.structure_views) or ["calc"],
            formats=_parse_csv_list(args.structure_formats) or ["svg"],
        )

    index_entries = []

    for res in results:
        d = res.definition

        # Choose output columns
        if args.columns == "all":
            df_out = res.df.copy()
        elif args.columns == "category":
            cols = [id_col]
            # union categories
            for cat in d.include_categories:
                for c in registry_category_columns(df_full, cat):
                    if c not in cols:
                        cols.append(c)
            # also add filter/sort cols
            for f in d.filters:
                if f.column in df_full.columns and f.column not in cols:
                    cols.append(f.column)
            for s in d.sort:
                if s.column in df_full.columns and s.column not in cols:
                    cols.append(s.column)
            df_out = res.df[[c for c in cols if c in res.df.columns]].copy()
        else:
            cols = select_output_columns(df_full, id_col, d)
            df_out = res.df[[c for c in cols if c in res.df.columns]].copy()

        # Write CSV
        csv_out = outdir / f"{d.key}.csv"
        csv_out.write_text(df_out.to_csv(index=False), encoding="utf-8")

        # HTML per picklist
        html_href = ""
        struct_map = None
        if args.html:
            html_out = outdir / f"{d.key}.html"

            if asset_spec is not None and id_col in df_out.columns:
                cids = df_out[id_col].astype(str).tolist()
                struct_map = copy_structure_assets(cids=cids, asset=asset_spec, outdir=outdir)

            subtitle = f"{d.name} • rows: {len(df_out)} • source: {Path(csv_path).name}"
            write_picklist_html(
                outpath=html_out,
                title=f"Picklist: {d.name}",
                subtitle=subtitle,
                df=df_out,
                id_col=id_col,
                warnings=res.warnings,
                structure_map=struct_map,
            )
            html_href = html_out.name

        index_entries.append((d.key, d.name, int(len(df_out)), html_href or csv_out.name))

    if args.html:
        write_index_html(
            outpath=outdir / "index.html",
            title="MolProp Toolkit Picklists",
            subtitle=f"Source: {Path(csv_path).name} • Generated: {datetime.now().isoformat(timespec='seconds')}",
            entries=index_entries,
        )

    print(f"Picklists written to: {outdir}")
    for k, name, count, href in index_entries:
        print(f"- {k}: {count} rows → {href}")


if __name__ == "__main__":
    main()

