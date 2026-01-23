#!/usr/bin/env python3
"""══════════════════════════════════════════════════════════════════════════════
MEDCHEM FLAGS ANALYZER
═══════════════════════════════════════════════════════════════════════════════

This analyzer summarizes common medicinal chemistry structural flags such as
PAINS and BRENK and reports flagged compounds with their names.

A key design intent is to keep compound identifiers visible and to provide a
short narrative explaining what each filter is intended to catch.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Allow running this file directly (without installing the package).
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd

from molprop_toolkit.core import detect_id_column, read_csv, write_csv


DESCRIPTION = """MedChem structural flags

PAINS patterns are often used to identify frequent assay interferers and
promiscuous binders. BRENK and related filters collect fragments that
medicinal chemists frequently avoid due to instability, reactivity, or
historical risk. These flags are not absolute, but they should trigger
orthogonal validation and careful assay design.

The summary column `MedChem_Alerts` is typically the count of filter hits among
PAINS, BRENK, NIH, and ZINC.
"""


def analyze_medchem(csv_path: str, output_path: str | None = None, top_n: int = 20) -> pd.DataFrame:
    print("=" * 80)
    print("MEDCHEM FLAGS ANALYSIS")
    print("=" * 80)
    print(DESCRIPTION.strip() + "\n")

    df = read_csv(csv_path)
    id_col = detect_id_column(df)

    cols = [id_col]
    if "SMILES" in df.columns:
        cols.append("SMILES")

    for c in ["PAINS", "PAINS_Type", "BRENK", "NIH", "ZINC", "MedChem_Alerts"]:
        if c in df.columns:
            cols.append(c)

    out = df[cols].copy()

    print("-" * 80)
    print("Flag rates")
    print("-" * 80)

    n = len(out)
    for c in ["PAINS", "BRENK", "NIH", "ZINC"]:
        if c in out.columns:
            flagged = (out[c] == True).sum()
            print(f"{c}: {flagged}/{n} ({100*flagged/n:.1f}%)")

    if "MedChem_Alerts" in out.columns:
        vc = out["MedChem_Alerts"].value_counts(dropna=False).sort_index()
        print("\nTotal alerts distribution:")
        for k, v in vc.items():
            print(f"{k}: {v} ({100*v/n:.1f}%)")

    if "PAINS" in out.columns:
        pains = out[out["PAINS"] == True]
        print("\n" + "-" * 80)
        print("PAINS flagged compounds")
        print("-" * 80)
        print(f"count: {len(pains)}")
        if len(pains) > 0:
            show = [c for c in [id_col, "PAINS_Type", "SMILES"] if c in pains.columns]
            print(pains[show].head(top_n).to_string(index=False))

    if "MedChem_Alerts" in out.columns:
        clean = out[out["MedChem_Alerts"].fillna(0) == 0]
        print("\n" + "-" * 80)
        print("Compounds with zero medchem flags")
        print("-" * 80)
        print(f"count: {len(clean)}")
        print(clean[[id_col]].head(top_n).to_string(index=False))

    if output_path:
        write_csv(out, output_path)
        print(f"\nExported: {output_path}")

    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="MedChem flags analyzer")
    ap.add_argument("input", help="Results CSV")
    ap.add_argument("-o", "--output", help="Export CSV")
    ap.add_argument("--top", type=int, default=20, help="Rows to print")
    args = ap.parse_args()

    analyze_medchem(args.input, args.output, args.top)


if __name__ == "__main__":
    main()

