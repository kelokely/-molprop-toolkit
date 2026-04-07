#!/usr/bin/env python3
"""══════════════════════════════════════════════════════════════════════════════
QED ANALYZER
═══════════════════════════════════════════════════════════════════════════════

This analyzer summarizes QED (Quantitative Estimate of Drug-likeness) outputs.
QED is a continuous score between 0 and 1 that can be used for ranking.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Allow running this file directly (without installing the package).
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd

from molprop_toolkit.core import detect_id_column, read_table, write_csv


DESCRIPTION = """QED drug-likeness

QED is a continuous drug-likeness score derived from desirability functions fit
to distributions of known drugs. It provides a ranking-friendly alternative to
binary rule filters.

Typical informal interpretation uses thresholds around 0.7 (high), 0.5 (medium)
and 0.3 (low), though the right target depends on modality and project.
"""


def analyze_qed(csv_path: str, output_path: str | None = None, top_n: int = 20) -> pd.DataFrame:
    print("=" * 80)
    print("QED ANALYSIS")
    print("=" * 80)
    print(DESCRIPTION.strip() + "\n")

    df = read_table(csv_path)
    id_col = detect_id_column(df)

    qed_cols = [c for c in df.columns if c.startswith("QED")]
    cols = [id_col]
    if "SMILES" in df.columns:
        cols.append("SMILES")
    for c in ["MolWt", "LogP", "TPSA", "HBD", "HBA", "RotBonds", "AromaticRings"]:
        if c in df.columns:
            cols.append(c)
    cols.extend([c for c in qed_cols if c not in cols])

    out = df[cols].copy()

    if "QED" in out.columns:
        s = out["QED"].dropna()
        if len(s) > 0:
            print("-" * 80)
            print("Score summary")
            print("-" * 80)
            print(f"count={len(s)}, mean={s.mean():.3f}, median={s.median():.3f}, min={s.min():.3f}, max={s.max():.3f}, std={s.std():.3f}")

            print("\n" + "-" * 80)
            print("Top compounds by QED")
            print("-" * 80)
            show_cols = [c for c in [id_col, "QED", "MolWt", "LogP", "TPSA"] if c in out.columns]
            print(out.sort_values("QED", ascending=False)[show_cols].head(top_n).to_string(index=False))

            print("\n" + "-" * 80)
            print("Low QED compounds (<0.4)")
            print("-" * 80)
            low = out[out["QED"] < 0.4]
            print(f"count: {len(low)}")
            if len(low) > 0:
                print(low[show_cols].head(top_n).to_string(index=False))

    if output_path:
        write_csv(out, output_path)
        print(f"\nExported: {output_path}")

    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="QED analyzer")
    ap.add_argument("input", help="Results CSV")
    ap.add_argument("-o", "--output", help="Export CSV")
    ap.add_argument("--top", type=int, default=20, help="Rows to print")
    args = ap.parse_args()

    analyze_qed(args.input, args.output, args.top)


if __name__ == "__main__":
    main()

