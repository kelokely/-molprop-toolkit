#!/usr/bin/env python3
"""Rule of 5 analyzer.

This script focuses on Lipinski Rule of 5 columns and prints pass rates and the
largest violation counts with compound names.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Allow running this file directly (without installing the package).
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd

from molprop_toolkit.core import detect_id_column, read_table, write_csv

def analyze_ro5(csv_path: str, output_path: str | None = None, top_n: int = 20) -> pd.DataFrame:
    df = read_table(csv_path)
    id_col = detect_id_column(df)

    cols = [id_col]
    if "SMILES" in df.columns:
        cols.append("SMILES")
    for c in ["RO5_Pass", "RO5_Violations", "RO5_MW", "RO5_LogP", "RO5_HBD", "RO5_HBA", "MolWt", "LogP", "HBD", "HBA"]:
        if c in df.columns and c not in cols:
            cols.append(c)

    out = df[cols].copy()

    print("=" * 80)
    print("RULE OF 5 ANALYSIS")
    print("=" * 80)

    if "RO5_Pass" in out.columns:
        passed = (out["RO5_Pass"] == True).sum()
        print(f"RO5_Pass: {passed}/{len(out)} ({100*passed/len(out):.1f}%)")

    if "RO5_Violations" in out.columns:
        vc = out["RO5_Violations"].value_counts(dropna=False).sort_index()
        print("Violations distribution:")
        for k, v in vc.items():
            print(f"{k}: {v} ({100*v/len(out):.1f}%)")

        print("\nWorst compounds by violations:")
        show_cols = [c for c in [id_col, "RO5_Violations", "MolWt", "LogP", "HBD", "HBA"] if c in out.columns]
        print(out.sort_values("RO5_Violations", ascending=False)[show_cols].head(top_n).to_string(index=False))

    if output_path:
        write_csv(out, output_path)
        print(f"\nExported: {output_path}")

    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="Rule of 5 analyzer")
    ap.add_argument("input", help="Results CSV")
    ap.add_argument("-o", "--output", help="Export CSV")
    ap.add_argument("--top", type=int, default=20)
    args = ap.parse_args()

    analyze_ro5(args.input, args.output, args.top)


if __name__ == "__main__":
    main()

