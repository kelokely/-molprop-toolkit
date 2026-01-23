#!/usr/bin/env python3
"""══════════════════════════════════════════════════════════════════════════════
CNS MPO ANALYZER
═══════════════════════════════════════════════════════════════════════════════

This script analyzes CNS MPO outputs in a results CSV and prints a report that
keeps compound names visible in all summaries.

The score interpretation implemented here follows the widely used Pfizer CNS
MPO framing. The toolkit’s CNS MPO is a heuristic intended for triage and
ranking rather than a replacement for experimental BBB data.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Allow running this file directly (without installing the package).
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd

from molprop_toolkit.core import detect_id_column, read_csv, write_csv


DESCRIPTION = """CNS Multiparameter Optimization (MPO)

CNS MPO is a multi-parameter score designed to describe a region of property
space where CNS drugs tend to succeed. The score is composed of six component
scores that are mapped to desirability functions and summed to a total score
from 0 to 6.

A common working interpretation is that scores at or above 4 indicate a
reasonably CNS-like balance, scores at or above 5 indicate a strong CNS-like
balance, and scores below 3 indicate a poor balance for CNS penetration.

The component breakdown is particularly useful because it often points to the
dominant reason a compound is failing, such as excessive TPSA, too many HBDs,
excessive basicity, or excessive lipophilicity.
"""


def analyze_cns(csv_path: str, output_path: str | None = None, top_n: int = 20) -> pd.DataFrame:
    print("=" * 80)
    print("CNS MPO ANALYSIS")
    print("=" * 80)
    print(DESCRIPTION.strip() + "\n")

    df = read_csv(csv_path)
    id_col = detect_id_column(df)

    wanted = [
        id_col,
        "SMILES",
        "CNS_MPO",
        "CNS_MPO_Class",
        "BBB_Likely",
        "MPO_MW",
        "MPO_LogP",
        "MPO_LogD",
        "MPO_TPSA",
        "MPO_HBD",
        "MPO_pKa",
        "cLogD_7.4",
        "pKa_Basic",
        "pKa_Acidic",
        "MolWt",
        "LogP",
        "TPSA",
        "HBD",
    ]
    cols = [c for c in wanted if c in df.columns]
    out = df[cols].copy()

    if "CNS_MPO" in out.columns:
        s = out["CNS_MPO"].dropna()
        if len(s) > 0:
            print("-" * 80)
            print("Score summary")
            print("-" * 80)
            print(f"count={len(s)}, mean={s.mean():.2f}, median={s.median():.2f}, min={s.min():.2f}, max={s.max():.2f}, std={s.std():.2f}")

    if "CNS_MPO_Class" in out.columns:
        print("\n" + "-" * 80)
        print("Category distribution")
        print("-" * 80)
        vc = out["CNS_MPO_Class"].value_counts(dropna=False)
        for k, v in vc.items():
            print(f"{k}: {v} ({100*v/len(out):.1f}%)")

    if "BBB_Likely" in out.columns:
        print("\n" + "-" * 80)
        print("BBB likelihood")
        print("-" * 80)
        yes = (out["BBB_Likely"] == True).sum()
        print(f"BBB_Likely=True: {yes}/{len(out)} ({100*yes/len(out):.1f}%)")

    print("\n" + "-" * 80)
    print("Top candidates by CNS_MPO")
    print("-" * 80)
    if "CNS_MPO" in out.columns:
        show_cols = [c for c in [id_col, "CNS_MPO", "CNS_MPO_Class", "BBB_Likely", "MolWt", "LogP", "TPSA", "HBD"] if c in out.columns]
        print(out.sort_values("CNS_MPO", ascending=False)[show_cols].head(top_n).to_string(index=False))
    else:
        print(out.head(top_n).to_string(index=False))

    if output_path:
        write_csv(out, output_path)
        print(f"\nExported: {output_path}")

    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="CNS MPO analyzer")
    ap.add_argument("input", help="Results CSV")
    ap.add_argument("-o", "--output", help="Export CSV")
    ap.add_argument("--top", type=int, default=20, help="Rows to print")
    args = ap.parse_args()

    analyze_cns(args.input, args.output, args.top)


if __name__ == "__main__":
    main()

