#!/usr/bin/env python3
"""══════════════════════════════════════════════════════════════════════════════
METABOLISM ANALYZER
═══════════════════════════════════════════════════════════════════════════════

This analyzer summarizes metabolic soft spots and stability class estimates.
It looks for columns prefixed with `Met_` and keeps compound identifiers
visible in all outputs.

The signals are heuristic; they are most useful for ranking within a chemical
series and for prioritizing metabolic stability experiments.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Allow running this file directly (without installing the package).
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd

from molprop_toolkit.core import detect_id_column, read_table, write_csv


DESCRIPTION = """Metabolism and soft spots

The Met_ columns typically represent counts of structural motifs that are
commonly oxidized or dealkylated, as well as conjugation handles that may be
glucuronidated or sulfated. A higher total soft spot count often correlates
with faster metabolism, especially for lipophilic compounds.

The `Met_Stability` class is a coarse category derived from a combination of
soft spots and lipophilicity. It is not a substitute for microsomal or
hepatocyte stability.
"""


def analyze_metabolism(csv_path: str, output_path: str | None = None, top_n: int = 20) -> pd.DataFrame:
    print("=" * 80)
    print("METABOLISM ANALYSIS")
    print("=" * 80)
    print(DESCRIPTION.strip() + "\n")

    df = read_table(csv_path)
    id_col = detect_id_column(df)

    met_cols = [c for c in df.columns if c.startswith("Met_")]
    cols = [id_col]
    if "SMILES" in df.columns:
        cols.append("SMILES")
    for c in ["MolWt", "LogP", "TPSA"]:
        if c in df.columns:
            cols.append(c)
    cols.extend([c for c in met_cols if c not in cols])
    out = df[cols].copy()

    if "Met_Stability" in out.columns:
        print("-" * 80)
        print("Stability distribution")
        print("-" * 80)
        vc = out["Met_Stability"].value_counts(dropna=False)
        for k, v in vc.items():
            print(f"{k}: {v} ({100*v/len(out):.1f}%)")

    if "Met_Total_Sites" in out.columns:
        s = out["Met_Total_Sites"].dropna()
        if len(s) > 0:
            print("\n" + "-" * 80)
            print("Total site summary")
            print("-" * 80)
            print(f"mean={s.mean():.2f}, median={s.median():.2f}, min={s.min():.0f}, max={s.max():.0f}")

    print("\n" + "-" * 80)
    print("Most common soft spot columns")
    print("-" * 80)

    counts = {}
    for c in met_cols:
        if c in ("Met_Stability", "Met_Total_Sites"):
            continue
        if pd.api.types.is_numeric_dtype(out[c]):
            nflag = (out[c] > 0).sum()
        else:
            nflag = (out[c] == True).sum()
        if nflag > 0:
            counts[c] = int(nflag)

    for name, nflag in sorted(counts.items(), key=lambda kv: -kv[1])[:25]:
        print(f"{name}: {nflag} ({100*nflag/len(out):.1f}%)")

    print("\n" + "-" * 80)
    print("Compounds with lowest Met_Total_Sites")
    print("-" * 80)

    if "Met_Total_Sites" in out.columns:
        show_cols = [c for c in [id_col, "Met_Stability", "Met_Total_Sites", "LogP", "MolWt"] if c in out.columns]
        print(out.sort_values("Met_Total_Sites", ascending=True)[show_cols].head(top_n).to_string(index=False))

    if output_path:
        write_csv(out, output_path)
        print(f"\nExported: {output_path}")

    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="Metabolism analyzer")
    ap.add_argument("input", help="Results CSV")
    ap.add_argument("-o", "--output", help="Export CSV")
    ap.add_argument("--top", type=int, default=20, help="Rows to print")
    args = ap.parse_args()

    analyze_metabolism(args.input, args.output, args.top)


if __name__ == "__main__":
    main()

