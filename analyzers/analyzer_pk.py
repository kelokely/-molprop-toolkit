#!/usr/bin/env python3
"""══════════════════════════════════════════════════════════════════════════════
PK ANALYZER
═══════════════════════════════════════════════════════════════════════════════

This analyzer summarizes heuristic PK-related columns such as PPB, Vd,
clearance risk, half-life class, and bioavailability class.

These are coarse indicators and are not substitutes for measured PK.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Allow running this file directly (without installing the package).
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd

from molprop_toolkit.core import detect_id_column, read_table, write_csv


DESCRIPTION = """Pharmacokinetics heuristics

This report is a triage-oriented view of PK-related columns. It summarizes
protein binding class, distribution class, clearance risk, half-life class, and
bioavailability class when present.

The intent is to highlight compounds that are likely to suffer from very low
bioavailability, rapid clearance, or extreme protein binding, so that those
issues can be addressed early by redesign or by targeted experiments.
"""


def analyze_pk(csv_path: str, output_path: str | None = None, top_n: int = 20) -> pd.DataFrame:
    print("=" * 80)
    print("PK ANALYSIS")
    print("=" * 80)
    print(DESCRIPTION.strip() + "\n")

    df = read_table(csv_path)
    id_col = detect_id_column(df)

    cols = [id_col]
    if "SMILES" in df.columns:
        cols.append("SMILES")

    for c in [
        "PPB_Est_%",
        "PPB_Class",
        "Vd_Est_L_kg",
        "Vd_Class",
        "CL_Risk",
        "CL_SoftSpots",
        "T_half_Est",
        "F_Score",
        "F_Est_%",
        "F_Class",
        "MolWt",
        "LogP",
        "Met_Stability",
        "Perm_GI_Absorption",
    ]:
        if c in df.columns and c not in cols:
            cols.append(c)

    out = df[cols].copy()

    def dist(col: str) -> None:
        if col in out.columns:
            print("\n" + "-" * 80)
            print(f"{col} distribution")
            print("-" * 80)
            vc = out[col].value_counts(dropna=False)
            for k, v in vc.items():
                print(f"{k}: {v} ({100*v/len(out):.1f}%)")

    dist("PPB_Class")
    dist("Vd_Class")
    dist("CL_Risk")
    dist("T_half_Est")
    dist("F_Class")

    print("\n" + "-" * 80)
    print("Compounds with low predicted bioavailability")
    print("-" * 80)

    if "F_Class" in out.columns:
        low = out[out["F_Class"].astype(str).str.lower() == "low"].copy()
        print(f"count: {len(low)}")
        show_cols = [c for c in [id_col, "F_Class", "F_Est_%", "CL_Risk", "Met_Stability", "Perm_GI_Absorption"] if c in low.columns]
        if len(low) > 0:
            print(low[show_cols].head(top_n).to_string(index=False))

    if output_path:
        write_csv(out, output_path)
        print(f"\nExported: {output_path}")

    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="PK analyzer")
    ap.add_argument("input", help="Results CSV")
    ap.add_argument("-o", "--output", help="Export CSV")
    ap.add_argument("--top", type=int, default=20, help="Rows to print")
    args = ap.parse_args()

    analyze_pk(args.input, args.output, args.top)


if __name__ == "__main__":
    main()

