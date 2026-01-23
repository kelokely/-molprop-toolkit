#!/usr/bin/env python3
"""══════════════════════════════════════════════════════════════════════════════
MASTER MOLECULAR PROPERTY ANALYZER
═══════════════════════════════════════════════════════════════════════════════

This is the umbrella CLI for analyzing a results CSV produced by the MolProp
calculators (or any compatible table). It is designed to keep the compound
identifier visible in outputs and to print a narrative description of each
category before computing summaries.

Task 10 moves shared logic (column detection and the category registry) into an
importable package namespace so other tools (reports, picklists, schema) can
reuse it.
"""

from __future__ import annotations

import argparse
import sys
from datetime import datetime
from pathlib import Path
from typing import List, Optional

# Allow running this file directly (without installing the package).
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd

from molprop_toolkit.core import detect_id_column, read_csv
from molprop_toolkit.core.registry import CATEGORY_SPECS


class MasterAnalyzer:
    def __init__(self, csv_path: str):
        self.csv_path = csv_path
        self.df = read_csv(csv_path)
        self.id_col = detect_id_column(self.df)

    def list_categories(self) -> None:
        print("\nAvailable categories:\n")
        for k, spec in CATEGORY_SPECS.items():
            cols = self.category_columns(k)
            print(f"  {k:18s} | {spec.name} | columns available: {len(cols)}")

    def category_columns(self, category: str) -> List[str]:
        if category not in CATEGORY_SPECS:
            return []
        info = CATEGORY_SPECS[category]

        cols: List[str] = []
        cols.extend([c for c in info.columns if c in self.df.columns])

        if info.prefix:
            for c in self.df.columns:
                if c.startswith(info.prefix) and c not in cols:
                    cols.append(c)

        return cols

    def analyze_category(self, category: str, top_n: int = 20) -> pd.DataFrame:
        if category not in CATEGORY_SPECS:
            raise ValueError(f"Unknown category: {category}")

        info = CATEGORY_SPECS[category]
        cols = self.category_columns(category)
        if not cols:
            raise ValueError(f"No columns found for category: {category}")

        print("\n" + "=" * 80)
        print(info.name)
        print("=" * 80)
        print(info.description.strip() + "\n")

        out_cols = [self.id_col] + cols
        out_cols = [c for c in out_cols if c in self.df.columns]
        df_cat = self.df[out_cols].copy()

        key = info.key_column
        if key and key in df_cat.columns:
            print("-" * 80)
            print("Summary")
            print("-" * 80)
            if pd.api.types.is_numeric_dtype(df_cat[key]):
                vals = df_cat[key].dropna()
                if len(vals) > 0:
                    print(
                        f"{key}: mean={vals.mean():.3f}, median={vals.median():.3f}, "
                        f"min={vals.min():.3f}, max={vals.max():.3f}, std={vals.std():.3f}"
                    )
            else:
                vc = df_cat[key].value_counts(dropna=False)
                for k2, v2 in vc.items():
                    print(f"{key}={k2}: {v2} ({100*v2/len(df_cat):.1f}%)")

        print("\n" + "-" * 80)
        print("Compound-level preview")
        print("-" * 80)

        preview_cols = [self.id_col]
        if key and key in df_cat.columns:
            preview_cols.append(key)
        for c in cols:
            if c not in preview_cols:
                preview_cols.append(c)
            if len(preview_cols) >= 8:
                break

        print(df_cat[preview_cols].head(top_n).to_string(index=False))
        return df_cat

    def report_text(self) -> str:
        lines: List[str] = []
        lines.append("=" * 80)
        lines.append("MolProp Toolkit report")
        lines.append("=" * 80)
        lines.append(f"dataset: {self.csv_path}")
        lines.append(f"compounds: {len(self.df)}")
        lines.append(f"columns: {len(self.df.columns)}")
        lines.append(f"generated: {datetime.now().isoformat(timespec='seconds')}")

        for cat, info in CATEGORY_SPECS.items():
            cols = self.category_columns(cat)
            if not cols:
                continue

            key = info.key_column
            lines.append("\n" + "-" * 80)
            lines.append(f"{info.name} ({cat})")
            lines.append("-" * 80)

            if key and key in self.df.columns:
                if pd.api.types.is_numeric_dtype(self.df[key]):
                    vals = self.df[key].dropna()
                    if len(vals) > 0:
                        lines.append(
                            f"{key}: mean={vals.mean():.3f}, median={vals.median():.3f}, "
                            f"min={vals.min():.3f}, max={vals.max():.3f}"
                        )
                        top = self.df[[self.id_col, key]].nlargest(min(5, len(self.df)), key)
                        lines.append("top compounds:")
                        for _, r in top.iterrows():
                            lines.append(f"  {r[self.id_col]} | {key}={r[key]}")
                else:
                    vc = self.df[key].value_counts(dropna=False)
                    for k2, v2 in vc.items():
                        lines.append(f"{key}={k2}: {v2} ({100*v2/len(self.df):.1f}%)")

        return "\n".join(lines)

    def write_report(self, output_path: str) -> None:
        Path(output_path).write_text(self.report_text(), encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser(description="MolProp Toolkit master analyzer")
    ap.add_argument("input", help="Results CSV")
    ap.add_argument("--list", action="store_true", help="List categories")
    ap.add_argument("--category", "-c", help="Analyze a single category")
    ap.add_argument("--top", type=int, default=20, help="Preview rows")
    ap.add_argument("--report", action="store_true", help="Generate an overall report")
    ap.add_argument("-o", "--output", help="Output path for report")
    args = ap.parse_args()

    analyzer = MasterAnalyzer(args.input)

    if args.list:
        analyzer.list_categories()
        return

    if args.report:
        text = analyzer.report_text()
        if args.output:
            analyzer.write_report(args.output)
        print(text)
        return

    if args.category:
        analyzer.analyze_category(args.category, top_n=args.top)
        return

    print("Use --list to see categories, --category <name> to analyze one, or --report to generate a report.")


if __name__ == "__main__":
    main()

