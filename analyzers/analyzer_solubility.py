#!/usr/bin/env python3
"""══════════════════════════════════════════════════════════════════════════════
SOLUBILITY ANALYZER
═══════════════════════════════════════════════════════════════════════════════

This analyzer summarizes predicted aqueous solubility columns prefixed with
`Sol_`. It is designed for ranking and triage rather than absolute prediction.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Allow running this file directly (without installing the package).
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd

from molprop_toolkit.core import detect_id_column, read_csv, write_csv


DESCRIPTION = """Solubility

Aqueous solubility is a frequent bottleneck for oral exposure. Simple logS
estimators such as ESOL can be useful for ranking within a series. A logS below
about -4 is often treated as a warning sign for oral programs unless there is a
clear formulation plan.

This analyzer summarizes common solubility columns such as Sol_ESOL_LogS and
Sol_ESOL_Class when present and lists the least soluble compounds.
"""


def analyze_solubility(csv_path: str, output_path: str | None = None, top_n: int = 20) -> pd.DataFrame:
    print("=" * 80)
    print("SOLUBILITY ANALYSIS")
    print("=" * 80)
    print(DESCRIPTION.strip() + "\n")

    df = read_csv(csv_path)
    id_col = detect_id_column(df)

    sol_cols = [c for c in df.columns if c.startswith("Sol_")]
    cols = [id_col]
    if "SMILES" in df.columns:
        cols.append("SMILES")
    for c in ["MolWt", "LogP", "TPSA", "AromaticRings", "FractionCSP3"]:
        if c in df.columns:
            cols.append(c)
    cols.extend([c for c in sol_cols if c not in cols])

    out = df[cols].copy()

    if "Sol_ESOL_Class" in out.columns:
        print("-" * 80)
        print("Class distribution")
        print("-" * 80)
        vc = out["Sol_ESOL_Class"].value_counts(dropna=False)
        for k, v in vc.items():
            print(f"{k}: {v} ({100*v/len(out):.1f}%)")

    if "Sol_ESOL_LogS" in out.columns:
        s = out["Sol_ESOL_LogS"].dropna()
        if len(s) > 0:
            print("\n" + "-" * 80)
            print("ESOL logS summary")
            print("-" * 80)
            print(f"mean={s.mean():.2f}, median={s.median():.2f}, min={s.min():.2f}, max={s.max():.2f}")

            print("\n" + "-" * 80)
            print("Least soluble compounds")
            print("-" * 80)
            show_cols = [c for c in [id_col, "Sol_ESOL_LogS", "Sol_ESOL_Class", "LogP", "MolWt"] if c in out.columns]
            print(out.sort_values("Sol_ESOL_LogS", ascending=True)[show_cols].head(top_n).to_string(index=False))

    if output_path:
        write_csv(out, output_path)
        print(f"\nExported: {output_path}")

    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="Solubility analyzer")
    ap.add_argument("input", help="Results CSV")
    ap.add_argument("-o", "--output", help="Export CSV")
    ap.add_argument("--top", type=int, default=20, help="Rows to print")
    args = ap.parse_args()

    analyze_solubility(args.input, args.output, args.top)


if __name__ == "__main__":
    main()

