#!/usr/bin/env python3
"""══════════════════════════════════════════════════════════════════════════════
LEAD OPTIMIZATION METRICS ANALYZER
═══════════════════════════════════════════════════════════════════════════════

This analyzer summarizes lead optimization metrics such as LE, LLE and LELP.
These metrics are most meaningful when potency data are included in the table.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Allow running this file directly (without installing the package).
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd

from molprop_toolkit.core import detect_id_column, read_table, write_csv


DESCRIPTION = """Lead optimization metrics

Ligand efficiency (LE) and lipophilic ligand efficiency (LLE) are commonly used
for hit-to-lead and lead optimization to avoid increasing molecular size and
lipophilicity without gaining potency in an efficient manner.

If your table includes Lead_LE and Lead_LLE, this analyzer will summarize their
distributions and list the compounds with the best apparent efficiency.
"""


def analyze_lead(csv_path: str, output_path: str | None = None, top_n: int = 20) -> pd.DataFrame:
    print("=" * 80)
    print("LEAD OPTIMIZATION ANALYSIS")
    print("=" * 80)
    print(DESCRIPTION.strip() + "\n")

    df = read_table(csv_path)
    id_col = detect_id_column(df)

    lead_cols = [c for c in df.columns if c.startswith("Lead_")]
    cols = [id_col]
    if "SMILES" in df.columns:
        cols.append("SMILES")
    for c in ["MolWt", "LogP", "TPSA", "QED"]:
        if c in df.columns:
            cols.append(c)
    cols.extend([c for c in lead_cols if c not in cols])

    out = df[cols].copy()

    if "Lead_LE" in out.columns:
        s = out["Lead_LE"].dropna()
        if len(s) > 0:
            print("-" * 80)
            print("LE summary")
            print("-" * 80)
            print(f"mean={s.mean():.3f}, median={s.median():.3f}, min={s.min():.3f}, max={s.max():.3f}")

    if "Lead_LLE" in out.columns:
        s = out["Lead_LLE"].dropna()
        if len(s) > 0:
            print("\n" + "-" * 80)
            print("LLE summary")
            print("-" * 80)
            print(f"mean={s.mean():.2f}, median={s.median():.2f}, min={s.min():.2f}, max={s.max():.2f}")

    print("\n" + "-" * 80)
    print("Top compounds by LLE then LE")
    print("-" * 80)

    if "Lead_LLE" in out.columns:
        sort_cols = ["Lead_LLE"]
        if "Lead_LE" in out.columns:
            sort_cols.append("Lead_LE")
        show_cols = [c for c in [id_col, "Lead_LLE", "Lead_LE", "Lead_LELP", "MolWt", "LogP"] if c in out.columns]
        print(out.sort_values(sort_cols, ascending=False)[show_cols].head(top_n).to_string(index=False))
    else:
        print(out.head(top_n).to_string(index=False))

    if output_path:
        write_csv(out, output_path)
        print(f"\nExported: {output_path}")

    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="Lead optimization analyzer")
    ap.add_argument("input", help="Results CSV")
    ap.add_argument("-o", "--output", help="Export CSV")
    ap.add_argument("--top", type=int, default=20, help="Rows to print")
    args = ap.parse_args()

    analyze_lead(args.input, args.output, args.top)


if __name__ == "__main__":
    main()

