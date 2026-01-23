#!/usr/bin/env python3
"""══════════════════════════════════════════════════════════════════════════════
DEVELOPABILITY ANALYZER (Task 8)
═══════════════════════════════════════════════════════════════════════════════

Summarizes simple, interpretable developability indices and risk burdens.

Expected columns (when computed):
- Dev_Score, Dev_Class
- Dev_GSK_4_400_Pass, Dev_Pfizer_3_75_Pass, Dev_GoldenTriangle_Pass
- Dev_Risk_Burden, Dev_Risk_Class

These columns are produced by calculators (v4/v5) via calculators/developability.py.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Allow running this file directly (without installing the package).
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd

from molprop_toolkit.core import detect_id_column, read_csv, write_csv


DESCRIPTION = """Developability indices

These are pragmatic, rule-style filters and a composite risk/score intended to
help triage compounds that are likely to be harder to develop due to property
space (size/lipophilicity/polarity) and obvious red flags (tox/hERG/etc.).

The composite Dev_Score is heuristic and should be tuned for modality and
project context.
"""


def analyze_developability(csv_path: str, output_path: str | None = None, top_n: int = 20) -> pd.DataFrame:
    print("=" * 80)
    print("DEVELOPABILITY ANALYSIS")
    print("=" * 80)
    print(DESCRIPTION.strip() + "\n")

    df = read_csv(csv_path)
    id_col = detect_id_column(df)

    dev_cols = [c for c in df.columns if c.startswith("Dev_")]
    cols = [id_col]
    if "SMILES" in df.columns:
        cols.append("SMILES")

    # Include context columns when present
    for c in ["CNS_MPO", "Bioavail_Score", "RO5_Violations", "Tox_Risk", "hERG_Risk", "CYP_Total_Inhibition", "CL_Risk", "Sol_ESOL_LogS"]:
        if c in df.columns and c not in cols:
            cols.append(c)

    cols.extend([c for c in dev_cols if c not in cols])

    out = df[cols].copy()

    # Pass-rate summary
    for rule in ["Dev_GSK_4_400_Pass", "Dev_Pfizer_3_75_Pass", "Dev_GoldenTriangle_Pass"]:
        if rule in out.columns:
            passed = int((out[rule] == True).sum())
            print(f"{rule}: {passed}/{len(out)} ({100*passed/len(out):.1f}%)")

    if "Dev_Risk_Class" in out.columns:
        print("\nRisk class distribution:")
        print(out["Dev_Risk_Class"].value_counts(dropna=False).to_string())

    if "Dev_Class" in out.columns:
        print("\nDevelopability class distribution:")
        print(out["Dev_Class"].value_counts(dropna=False).to_string())

    if "Dev_Score" in out.columns:
        s = pd.to_numeric(out["Dev_Score"], errors="coerce").dropna()
        if len(s) > 0:
            print("\nDev_Score summary:")
            print(f"mean={s.mean():.2f}, median={s.median():.2f}, min={s.min():.0f}, max={s.max():.0f}")

            show = [c for c in [id_col, "Dev_Score", "Dev_Class", "Dev_Risk_Burden", "Dev_Risk_Class"] if c in out.columns]
            print("\nTop compounds by Dev_Score (then lowest risk burden):")
            sort_cols = [c for c in ["Dev_Score", "Dev_Risk_Burden"] if c in out.columns]
            asc = [False, True][: len(sort_cols)]
            print(out.sort_values(sort_cols, ascending=asc)[show].head(top_n).to_string(index=False))

    if output_path:
        write_csv(out, output_path)
        print(f"\nExported: {output_path}")

    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="Developability analyzer (Task 8)")
    ap.add_argument("input", help="Results CSV")
    ap.add_argument("-o", "--output", help="Export CSV")
    ap.add_argument("--top", type=int, default=20, help="Rows to print")
    args = ap.parse_args()

    analyze_developability(args.input, args.output, args.top)


if __name__ == "__main__":
    main()

