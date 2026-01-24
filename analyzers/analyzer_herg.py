#!/usr/bin/env python3
"""══════════════════════════════════════════════════════════════════════════════
hERG ANALYZER
═══════════════════════════════════════════════════════════════════════════════

This analyzer summarizes hERG risk columns and prints compound-level tables.
The heuristics are meant for triage and prioritization of experimental testing.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Allow running this file directly (without installing the package).
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd

from molprop_toolkit.core import detect_id_column, read_table, write_csv


DESCRIPTION = """hERG liability

The hERG potassium channel is a major cardiac safety liability. Block of hERG
can lead to QT prolongation and torsades de pointes risk. The hERG outputs in
this toolkit are heuristic signals based on common structure–property
correlations, especially the presence of a basic nitrogen together with
lipophilicity and aromatic/hydrophobic features.

This report is intended to help you decide which compounds should be prioritized
for early hERG assays and where medicinal chemistry mitigation might be needed.
"""


def analyze_herg(csv_path: str, output_path: str | None = None, top_n: int = 20) -> pd.DataFrame:
    print("=" * 80)
    print("hERG RISK ANALYSIS")
    print("=" * 80)
    print(DESCRIPTION.strip() + "\n")

    df = read_table(csv_path)
    id_col = detect_id_column(df)

    herg_cols = [c for c in df.columns if c.startswith("hERG_")]
    cols = [id_col]
    if "SMILES" in df.columns:
        cols.append("SMILES")
    for c in ["MolWt", "LogP", "TPSA", "pKa_Basic", "AromaticRings"]:
        if c in df.columns:
            cols.append(c)
    cols.extend([c for c in herg_cols if c not in cols])
    out = df[cols].copy()

    if "hERG_Risk" in out.columns:
        print("-" * 80)
        print("Risk distribution")
        print("-" * 80)
        vc = out["hERG_Risk"].value_counts(dropna=False)
        for k, v in vc.items():
            print(f"{k}: {v} ({100*v/len(out):.1f}%)")

    if "hERG_Alerts" in out.columns:
        s = out["hERG_Alerts"].dropna()
        if len(s) > 0:
            print("\n" + "-" * 80)
            print("Alert count summary")
            print("-" * 80)
            print(f"mean={s.mean():.2f}, median={s.median():.2f}, min={s.min():.0f}, max={s.max():.0f}")

    print("\n" + "-" * 80)
    print("High risk compounds")
    print("-" * 80)

    if "hERG_Risk" in out.columns:
        high = out[out["hERG_Risk"].astype(str).str.lower() == "high"].copy()
        print(f"high risk count: {len(high)}")
        if len(high) > 0:
            show_cols = [c for c in [id_col, "hERG_Risk", "hERG_Alerts", "MolWt", "LogP", "pKa_Basic"] if c in out.columns]
            print(high[show_cols].head(top_n).to_string(index=False))

    if output_path:
        write_csv(out, output_path)
        print(f"\nExported: {output_path}")

    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="hERG analyzer")
    ap.add_argument("input", help="Results CSV")
    ap.add_argument("-o", "--output", help="Export CSV")
    ap.add_argument("--top", type=int, default=20, help="Rows to print")
    args = ap.parse_args()

    analyze_herg(args.input, args.output, args.top)


if __name__ == "__main__":
    main()

