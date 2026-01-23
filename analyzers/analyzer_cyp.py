#!/usr/bin/env python3
"""══════════════════════════════════════════════════════════════════════════════
CYP INHIBITION ANALYZER
═══════════════════════════════════════════════════════════════════════════════

This analyzer summarizes CYP inhibition risk flags and keeps compound names
visible in all outputs.

The signals are heuristic and are not substitutes for in vitro CYP inhibition
panels or time-dependent inhibition (TDI) studies.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Allow running this file directly (without installing the package).
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd

from molprop_toolkit.core import detect_id_column, read_csv, write_csv


DESCRIPTION = """CYP inhibition and DDI risk

Cytochrome P450 enzymes metabolize a large fraction of small-molecule drugs.
Inhibition can increase exposure of co-administered drugs and create clinically
significant drug–drug interactions.

This analyzer looks for typical columns such as CYP3A4_Inhibitor,
CYP2D6_Inhibitor, CYP2C9_Inhibitor, CYP_MBI_Risk, and CYP_Total_Inhibition.
When the total count is present, compounds inhibiting multiple isoforms are
highlighted as higher risk.

If a mechanism-based inactivation flag is present, it should trigger early TDI
experiments because irreversible inhibition has different clinical behavior
than reversible inhibition.
"""


def analyze_cyp(csv_path: str, output_path: str | None = None, top_n: int = 20) -> pd.DataFrame:
    print("=" * 80)
    print("CYP INHIBITION ANALYSIS")
    print("=" * 80)
    print(DESCRIPTION.strip() + "\n")

    df = read_csv(csv_path)
    id_col = detect_id_column(df)

    cols = [id_col]
    if "SMILES" in df.columns:
        cols.append("SMILES")

    for c in [
        "CYP3A4_Inhibitor",
        "CYP2D6_Inhibitor",
        "CYP2C9_Inhibitor",
        "CYP_MBI_Risk",
        "CYP_Total_Inhibition",
        "MolWt",
        "LogP",
        "pKa_Basic",
    ]:
        if c in df.columns and c not in cols:
            cols.append(c)

    out = df[cols].copy()

    print("-" * 80)
    print("Inhibitor rates")
    print("-" * 80)

    n = len(out)
    for c in ["CYP3A4_Inhibitor", "CYP2D6_Inhibitor", "CYP2C9_Inhibitor", "CYP_MBI_Risk"]:
        if c in out.columns:
            hits = (out[c] == True).sum()
            print(f"{c}: {hits}/{n} ({100*hits/n:.1f}%)")

    if "CYP_Total_Inhibition" in out.columns:
        print("\n" + "-" * 80)
        print("Total CYP inhibition distribution")
        print("-" * 80)
        vc = out["CYP_Total_Inhibition"].value_counts(dropna=False).sort_index()
        for k, v in vc.items():
            print(f"{k}: {v} ({100*v/n:.1f}%)")

        multi = out[out["CYP_Total_Inhibition"].fillna(0) >= 2]
        print("\n" + "-" * 80)
        print("Multi-isoform inhibitors (>=2)")
        print("-" * 80)
        print(f"count: {len(multi)}")
        if len(multi) > 0:
            show_cols = [c for c in [id_col, "CYP_Total_Inhibition", "CYP_MBI_Risk", "CYP3A4_Inhibitor", "CYP2D6_Inhibitor", "CYP2C9_Inhibitor"] if c in out.columns]
            print(multi[show_cols].head(top_n).to_string(index=False))

    if output_path:
        write_csv(out, output_path)
        print(f"\nExported: {output_path}")

    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="CYP inhibition analyzer")
    ap.add_argument("input", help="Results CSV")
    ap.add_argument("-o", "--output", help="Export CSV")
    ap.add_argument("--top", type=int, default=20, help="Rows to print")
    args = ap.parse_args()

    analyze_cyp(args.input, args.output, args.top)


if __name__ == "__main__":
    main()

