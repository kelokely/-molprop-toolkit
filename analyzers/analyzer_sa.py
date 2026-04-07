#!/usr/bin/env python3
"""══════════════════════════════════════════════════════════════════════════════
SYNTHETIC ACCESSIBILITY + COMPLEXITY ANALYZER (Task 6)
═══════════════════════════════════════════════════════════════════════════════

This analyzer summarizes:
- SA_Score (Ertl & Schuffenhauer, 1–10; lower is easier)
- a small bundle of interpretable complexity descriptors

This is intended for triage and prioritization. SA_Score is a heuristic: it is
useful for ranking within a library/series, but it is not a substitute for
retrosynthesis or chemist review.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

# Allow running this file directly (without installing the package).
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from molprop_toolkit.core import detect_id_column, read_table, write_csv


DESCRIPTION = """Synthetic accessibility and complexity

SA_Score is a widely used heuristic that combines fragment contributions and
molecular complexity penalties. Values are scaled to 1–10.

A common informal interpretation is:
- 1–3: easy / common medicinal chemistry
- 3–6: moderate
- 6–10: hard (special chemistry / route risk)

The accompanying Complexity_* columns are intended to help interpret *why* a
molecule scores as difficult (spiro/bridgeheads, stereocenters, aromaticity,
Bertz complexity, etc.).
"""


def analyze_sa(csv_path: str, output_path: str | None = None, top_n: int = 20) -> pd.DataFrame:
    print("=" * 80)
    print("SYNTHETIC ACCESSIBILITY + COMPLEXITY")
    print("=" * 80)
    print(DESCRIPTION.strip() + "\n")

    df = read_table(csv_path)
    id_col = detect_id_column(df)

    # Discover relevant columns.
    sa_cols = [c for c in df.columns if c.startswith("SA_") or c in ("SA_Score", "SA_Class")]
    cx_cols = [c for c in df.columns if c.startswith("Complexity_")]

    cols = [id_col]
    if "SMILES" in df.columns:
        cols.append("SMILES")

    for c in ["SA_Score", "SA_Class", "MolWt", "LogP", "TPSA", "FractionCSP3", "AromaticRings", "Rings"]:
        if c in df.columns and c not in cols:
            cols.append(c)

    for c in sa_cols:
        if c not in cols:
            cols.append(c)

    for c in sorted(cx_cols):
        if c not in cols:
            cols.append(c)

    out = df[cols].copy()

    if "SA_Score" in out.columns:
        s = pd.to_numeric(out["SA_Score"], errors="coerce").dropna()
        if len(s) > 0:
            print("-" * 80)
            print("SA_Score summary")
            print("-" * 80)
            print(
                f"count={len(s)}, mean={s.mean():.3f}, median={s.median():.3f}, min={s.min():.3f}, max={s.max():.3f}, std={s.std():.3f}"
            )

            show_cols = [c for c in [id_col, "SA_Score", "SA_Class", "MolWt", "LogP", "FractionCSP3"] if c in out.columns]

            print("\n" + "-" * 80)
            print("Hardest-to-synthesize (highest SA_Score)")
            print("-" * 80)
            print(out.sort_values("SA_Score", ascending=False)[show_cols].head(top_n).to_string(index=False))

            print("\n" + "-" * 80)
            print("Easiest-to-synthesize (lowest SA_Score)")
            print("-" * 80)
            print(out.sort_values("SA_Score", ascending=True)[show_cols].head(top_n).to_string(index=False))

            print("\n" + "-" * 80)
            print("Hard set (SA_Score >= 6)")
            print("-" * 80)
            hard = out[pd.to_numeric(out["SA_Score"], errors="coerce") >= 6]
            print(f"count: {len(hard)}")
            if len(hard) > 0:
                print(hard[show_cols].head(top_n).to_string(index=False))

    if output_path:
        write_csv(out, output_path)
        print(f"\nExported: {output_path}")

    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="Synthetic accessibility + complexity analyzer")
    ap.add_argument("input", help="Results CSV")
    ap.add_argument("-o", "--output", help="Export CSV")
    ap.add_argument("--top", type=int, default=20, help="Rows to print")
    args = ap.parse_args()

    analyze_sa(args.input, args.output, args.top)


if __name__ == "__main__":
    main()

