#!/usr/bin/env python3
"""══════════════════════════════════════════════════════════════════════════════
TOXICITY ALERT ANALYZER
═══════════════════════════════════════════════════════════════════════════════

This analyzer reads a results CSV and summarizes toxicity structural alert
columns. It expects summary columns such as `Tox_Risk` and `Tox_Total_Alerts`
when available, and will also include every column that begins with `Tox_`.

The outputs are structural alert signals rather than validated toxicity
predictions. They should be interpreted as prompts for follow-up experiments.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Allow running this file directly (without installing the package).
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd

from molprop_toolkit.core import detect_id_column, read_csv, write_csv


DESCRIPTION = """Toxicity structural alerts

This report summarizes substructure alerts that are commonly associated with
reactivity, mutagenicity risk, reactive metabolite formation, or historical
liability motifs. A higher alert count does not guarantee toxicity. The intent
is triage: identify compounds that require early experimental follow-up or
redesign.

When `Tox_Risk` is present, it is typically derived from the number of alert
families triggered. When `Tox_Total_Alerts` is present, it is a count of how
many alert types were triggered at least once.
"""


def analyze_toxicity(csv_path: str, output_path: str | None = None, top_n: int = 20) -> pd.DataFrame:
    print("=" * 80)
    print("TOXICITY ALERT ANALYSIS")
    print("=" * 80)
    print(DESCRIPTION.strip() + "\n")

    df = read_csv(csv_path)
    id_col = detect_id_column(df)

    tox_cols = [c for c in df.columns if c.startswith("Tox_")]
    cols = [id_col]
    if "SMILES" in df.columns:
        cols.append("SMILES")
    cols.extend([c for c in tox_cols if c not in cols])

    out = df[cols].copy()

    if "Tox_Risk" in out.columns:
        print("-" * 80)
        print("Risk distribution")
        print("-" * 80)
        vc = out["Tox_Risk"].value_counts(dropna=False)
        for k, v in vc.items():
            print(f"{k}: {v} ({100*v/len(out):.1f}%)")

    if "Tox_Total_Alerts" in out.columns:
        s = out["Tox_Total_Alerts"].dropna()
        if len(s) > 0:
            print("\n" + "-" * 80)
            print("Total alerts summary")
            print("-" * 80)
            print(f"mean={s.mean():.2f}, median={s.median():.2f}, min={s.min():.0f}, max={s.max():.0f}")

    print("\n" + "-" * 80)
    print("Most common alert types")
    print("-" * 80)

    alert_counts = {}
    for c in tox_cols:
        if c in ("Tox_Risk", "Tox_Total_Alerts"):
            continue
        series = out[c]
        if pd.api.types.is_numeric_dtype(series):
            nflag = (series > 0).sum()
        else:
            nflag = (series == True).sum()
        if nflag > 0:
            alert_counts[c] = int(nflag)

    for name, nflag in sorted(alert_counts.items(), key=lambda kv: -kv[1])[:25]:
        print(f"{name}: {nflag} ({100*nflag/len(out):.1f}%)")

    print("\n" + "-" * 80)
    print("Highest risk compounds")
    print("-" * 80)

    if "Tox_Risk" in out.columns:
        high = out[out["Tox_Risk"].astype(str).str.lower() == "high"].copy()
        print(f"high risk count: {len(high)}")
        if len(high) > 0:
            base_cols = [c for c in [id_col, "Tox_Risk", "Tox_Total_Alerts"] if c in out.columns]
            print(high[base_cols].head(top_n).to_string(index=False))

    if output_path:
        write_csv(out, output_path)
        print(f"\nExported: {output_path}")

    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="Toxicity alert analyzer")
    ap.add_argument("input", help="Results CSV")
    ap.add_argument("-o", "--output", help="Export CSV")
    ap.add_argument("--top", type=int, default=20, help="Rows to print")
    args = ap.parse_args()

    analyze_toxicity(args.input, args.output, args.top)


if __name__ == "__main__":
    main()

