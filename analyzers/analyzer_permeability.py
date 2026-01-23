#!/usr/bin/env python3
"""══════════════════════════════════════════════════════════════════════════════
PERMEABILITY ANALYZER
═══════════════════════════════════════════════════════════════════════════════

This analyzer summarizes permeability and absorption related columns prefixed
with `Perm_` and prints compound-level lists.

The outputs are coarse estimates intended for ranking and triage.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Allow running this file directly (without installing the package).
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd

from molprop_toolkit.core import detect_id_column, read_csv, write_csv


DESCRIPTION = """Permeability

Permeability is central to oral absorption and, for CNS projects, BBB
penetration. This analyzer summarizes the GI absorption label, the Caco-2 and
PAMPA log-permeability estimates when present, and a P-gp substrate likelihood
heuristic.

Because transporter effects and ionization can dominate in cells, these are not
intended as assay replacements. They are designed to highlight compounds that
are likely permeability-limited or efflux-limited.
"""


def analyze_permeability(csv_path: str, output_path: str | None = None, top_n: int = 20) -> pd.DataFrame:
    print("=" * 80)
    print("PERMEABILITY ANALYSIS")
    print("=" * 80)
    print(DESCRIPTION.strip() + "\n")

    df = read_csv(csv_path)
    id_col = detect_id_column(df)

    perm_cols = [c for c in df.columns if c.startswith("Perm_")]
    cols = [id_col]
    if "SMILES" in df.columns:
        cols.append("SMILES")
    for c in ["MolWt", "LogP", "TPSA", "HBD", "HBA", "BBB_Likely"]:
        if c in df.columns:
            cols.append(c)
    cols.extend([c for c in perm_cols if c not in cols])

    out = df[cols].copy()

    if "Perm_GI_Absorption" in out.columns:
        print("-" * 80)
        print("GI absorption distribution")
        print("-" * 80)
        vc = out["Perm_GI_Absorption"].value_counts(dropna=False)
        for k, v in vc.items():
            print(f"{k}: {v} ({100*v/len(out):.1f}%)")

    if "Perm_Pgp_Likelihood" in out.columns:
        print("\n" + "-" * 80)
        print("P-gp likelihood distribution")
        print("-" * 80)
        vc = out["Perm_Pgp_Likelihood"].value_counts(dropna=False)
        for k, v in vc.items():
            print(f"{k}: {v} ({100*v/len(out):.1f}%)")

    if "Perm_Caco2_LogPapp" in out.columns:
        s = out["Perm_Caco2_LogPapp"].dropna()
        if len(s) > 0:
            print("\n" + "-" * 80)
            print("Caco-2 LogPapp summary")
            print("-" * 80)
            print(f"mean={s.mean():.2f}, median={s.median():.2f}, min={s.min():.2f}, max={s.max():.2f}")

    print("\n" + "-" * 80)
    print("High absorption and low P-gp risk (screen)")
    print("-" * 80)

    if "Perm_GI_Absorption" in out.columns:
        subset = out[out["Perm_GI_Absorption"].astype(str).str.lower() == "high"].copy()
        if "Perm_Pgp_Likelihood" in subset.columns:
            subset = subset[subset["Perm_Pgp_Likelihood"].astype(str).str.lower().isin(["low", "medium"])]
        show_cols = [c for c in [id_col, "Perm_GI_Absorption", "Perm_Pgp_Likelihood", "TPSA", "LogP", "MolWt", "BBB_Likely"] if c in subset.columns]
        print(f"count: {len(subset)}")
        if len(subset) > 0:
            print(subset[show_cols].head(top_n).to_string(index=False))

    if output_path:
        write_csv(out, output_path)
        print(f"\nExported: {output_path}")

    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="Permeability analyzer")
    ap.add_argument("input", help="Results CSV")
    ap.add_argument("-o", "--output", help="Export CSV")
    ap.add_argument("--top", type=int, default=20, help="Rows to print")
    args = ap.parse_args()

    analyze_permeability(args.input, args.output, args.top)


if __name__ == "__main__":
    main()

