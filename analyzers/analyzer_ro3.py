#!/usr/bin/env python3
"""Rule of 3 analyzer.

This script focuses on fragment-likeness Rule of 3 columns and prints pass
rates and violation distributions with compound names.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Allow running this file directly (without installing the package).
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd

from molprop_toolkit.core import detect_id_column, read_table, write_csv

def analyze_ro3(csv_path: str, output_path: str | None = None, top_n: int = 20) -> pd.DataFrame:
    df = read_table(csv_path)
    id_col = detect_id_column(df)

    cols = [id_col]
    if "SMILES" in df.columns:
        cols.append("SMILES")
    for c in ["RO3_Pass", "RO3_Violations", "RO3_MW", "RO3_LogP", "RO3_HBD", "RO3_HBA", "RO3_RotBonds", "RO3_TPSA", "MolWt", "LogP", "TPSA"]:
        if c in df.columns and c not in cols:
            cols.append(c)

    out = df[cols].copy()

    print("=" * 80)
    print("RULE OF 3 ANALYSIS")
    print("=" * 80)

    if "RO3_Pass" in out.columns:
        passed = (out["RO3_Pass"] == True).sum()
        print(f"RO3_Pass: {passed}/{len(out)} ({100*passed/len(out):.1f}%)")

    if "RO3_Violations" in out.columns:
        vc = out["RO3_Violations"].value_counts(dropna=False).sort_index()
        print("Violations distribution:")
        for k, v in vc.items():
            print(f"{k}: {v} ({100*v/len(out):.1f}%)")

        print("\nWorst compounds by violations:")
        show_cols = [c for c in [id_col, "RO3_Violations", "MolWt", "LogP", "TPSA"] if c in out.columns]
        print(out.sort_values("RO3_Violations", ascending=False)[show_cols].head(top_n).to_string(index=False))

    if output_path:
        write_csv(out, output_path)
        print(f"\nExported: {output_path}")

    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="Rule of 3 analyzer")
    ap.add_argument("input", help="Results CSV")
    ap.add_argument("-o", "--output", help="Export CSV")
    ap.add_argument("--top", type=int, default=20)
    args = ap.parse_args()

    analyze_ro3(args.input, args.output, args.top)


if __name__ == "__main__":
    main()

