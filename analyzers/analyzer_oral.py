#!/usr/bin/env python3
"""══════════════════════════════════════════════════════════════════════════════
ORAL BIOAVAILABILITY ANALYZER
═══════════════════════════════════════════════════════════════════════════════

This analyzer summarizes oral bioavailability related rule outputs such as
Veber, Egan, Ghose, Muegge, and a combined Bioavail_Score.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Allow running this file directly (without installing the package).
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd

from molprop_toolkit.core import detect_id_column, read_csv, write_csv


DESCRIPTION = """Oral bioavailability filters

Oral bioavailability is influenced by solubility, permeability, stability, and
first-pass metabolism. Several rule sets are commonly used as fast, interpretable
filters. Veber focuses on TPSA and rotatable bonds, Egan uses a TPSA–LogP
boundary, and Ghose and Muegge are broader drug-likeness filters.

Bioavail_Score, when present, is typically the count of passed filters and is
useful for ranking.
"""


def analyze_oral(csv_path: str, output_path: str | None = None, top_n: int = 20) -> pd.DataFrame:
    print("=" * 80)
    print("ORAL BIOAVAILABILITY ANALYSIS")
    print("=" * 80)
    print(DESCRIPTION.strip() + "\n")

    df = read_csv(csv_path)
    id_col = detect_id_column(df)

    cols = [id_col]
    if "SMILES" in df.columns:
        cols.append("SMILES")

    for c in [
        "Bioavail_Score",
        "Veber_Pass",
        "Egan_Pass",
        "Ghose_Pass",
        "Muegge_Pass",
        "RO5_Pass",
        "MolWt",
        "LogP",
        "TPSA",
        "RotBonds",
        "HBD",
        "HBA",
    ]:
        if c in df.columns and c not in cols:
            cols.append(c)

    out = df[cols].copy()

    if "Bioavail_Score" in out.columns:
        s = out["Bioavail_Score"].dropna()
        if len(s) > 0:
            print("-" * 80)
            print("Bioavail_Score summary")
            print("-" * 80)
            print(f"mean={s.mean():.2f}, median={s.median():.2f}, min={s.min():.0f}, max={s.max():.0f}")

            print("\n" + "-" * 80)
            print("Top compounds by Bioavail_Score")
            print("-" * 80)
            show_cols = [c for c in [id_col, "Bioavail_Score", "Veber_Pass", "Egan_Pass", "Ghose_Pass", "Muegge_Pass"] if c in out.columns]
            print(out.sort_values("Bioavail_Score", ascending=False)[show_cols].head(top_n).to_string(index=False))

    for rule in ["Veber_Pass", "Egan_Pass", "Ghose_Pass", "Muegge_Pass", "RO5_Pass"]:
        if rule in out.columns:
            passed = (out[rule] == True).sum()
            print(f"{rule}: {passed}/{len(out)} ({100*passed/len(out):.1f}%)")

    if output_path:
        write_csv(out, output_path)
        print(f"\nExported: {output_path}")

    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="Oral bioavailability analyzer")
    ap.add_argument("input", help="Results CSV")
    ap.add_argument("-o", "--output", help="Export CSV")
    ap.add_argument("--top", type=int, default=20, help="Rows to print")
    args = ap.parse_args()

    analyze_oral(args.input, args.output, args.top)


if __name__ == "__main__":
    main()

