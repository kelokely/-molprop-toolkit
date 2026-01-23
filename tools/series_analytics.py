#!/usr/bin/env python3
"""Series analytics runner (Task 13).

Convenience wrapper that runs:
- Murcko scaffolds (Task 13A)
- Similarity clustering (Task 13B)

and writes an enriched CSV with:
- Scaffold_* columns
- Cluster_* columns

This lets downstream tooling (master analyzer, reports, picklists) consume the
series annotations directly.
"""

from __future__ import annotations

import argparse
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd

# Allow running directly without installing.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from molprop_toolkit.core import detect_id_column, read_csv

from analyzers.analyzer_clustering import analyze_clustering
from analyzers.analyzer_scaffolds import analyze_scaffolds


def main() -> None:
    ap = argparse.ArgumentParser(description="MolProp Toolkit series analytics (Task 13)")
    ap.add_argument("input", help="Results CSV")
    ap.add_argument("--outdir", default=None, help="Output directory")
    ap.add_argument("--smiles-col", default=None, help="Override SMILES column")
    ap.add_argument("--similarity", type=float, default=0.7, help="Tanimoto similarity threshold for clustering")
    ap.add_argument("--no-images", action="store_true", help="Disable scaffold/cluster SVG generation")
    ap.add_argument("--enriched-name", default=None, help="Name for enriched CSV (default: <stem>_series.csv)")
    args = ap.parse_args()

    in_path = Path(args.input)
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    outdir = Path(args.outdir) if args.outdir else Path("series") / f"series_{in_path.stem}_{ts}"
    outdir.mkdir(parents=True, exist_ok=True)

    top_images = 0 if args.no_images else 50

    # Run scaffold analysis
    scaf_dir = outdir / "scaffolds"
    analyze_scaffolds(args.input, outdir=scaf_dir, smiles_col=args.smiles_col, top_images=top_images)

    # Run clustering
    clus_dir = outdir / "clusters"
    analyze_clustering(args.input, outdir=clus_dir, smiles_col=args.smiles_col, similarity=args.similarity, top_images=top_images)

    # Build enriched CSV by merging on Compound_ID
    df = read_csv(args.input)
    id_col = detect_id_column(df)

    df_scaf = pd.read_csv(scaf_dir / "scaffold_assignments.csv")
    df_clus = pd.read_csv(clus_dir / "cluster_assignments.csv")

    keep_scaf = [id_col, "Scaffold_ID", "Scaffold_SMILES", "Scaffold_Size"]
    keep_clus = [id_col, "Cluster_ID", "Cluster_Size", "Cluster_Representative_ID", "Cluster_MedoidMeanSim"]

    df_enriched = df.merge(df_scaf[[c for c in keep_scaf if c in df_scaf.columns]], on=id_col, how="left")
    df_enriched = df_enriched.merge(df_clus[[c for c in keep_clus if c in df_clus.columns]], on=id_col, how="left")

    out_name = args.enriched_name or f"{in_path.stem}_series.csv"
    out_csv = outdir / out_name
    df_enriched.to_csv(out_csv, index=False)

    print(f"\nEnriched CSV written to: {out_csv}")
    print("Tip: run reports/analyzers on the enriched CSV to include series columns.")


if __name__ == "__main__":
    main()

