#!/usr/bin/env python3
"""Series analytics runner.

Convenience wrapper that runs Murcko scaffolds and similarity clustering and then writes an enriched CSV with `Scaffold_*`
and `Cluster_*` columns so downstream tooling (analyzers, reports, picklists) can consume series annotations directly.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd

# Allow running directly without installing.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

<<<<<<< HEAD
from molprop_toolkit.core import detect_id_column, read_table, write_run_metadata, write_table
=======
from molprop_toolkit.core import detect_id_column, read_table
>>>>>>> 5d70153 (Core: prefer Calc_Canonical_SMILES + add Parquet table IO)

from analyzers.analyzer_clustering import analyze_clustering
from analyzers.analyzer_scaffolds import analyze_scaffolds


def _sha256_file(path: Path, *, max_bytes: int = 200 * 1024 * 1024) -> str | None:
    try:
        st = path.stat()
        if st.st_size > max_bytes:
            return None
        h = hashlib.sha256()
        with path.open("rb") as f:
            for chunk in iter(lambda: f.read(1024 * 1024), b""):
                h.update(chunk)
        return h.hexdigest()
    except Exception:
        return None


def main() -> None:
    ap = argparse.ArgumentParser(description="MolProp Toolkit series analytics")
    ap.add_argument("input", help="Results table (CSV/TSV/Parquet)")
    ap.add_argument("--outdir", default=None, help="Output directory")
    ap.add_argument("--smiles-col", default=None, help="Override SMILES column")
    ap.add_argument("--similarity", type=float, default=0.7, help="Similarity threshold for clustering")

    ap.add_argument(
        "--fp",
        default="morgan",
        choices=["morgan", "rdkit", "maccs", "atom_pair", "torsion", "usr", "usrcat"],
        help="Fingerprint kind for clustering",
    )
    ap.add_argument("--fp-radius", type=int, default=2, help="Morgan radius")
    ap.add_argument("--fp-nbits", type=int, default=2048, help="Fingerprint size (where applicable)")
    ap.add_argument("--fp-chirality", action="store_true", help="Use chirality in Morgan fingerprints")
    ap.add_argument("--fp-features", action="store_true", help="Use feature-invariants in Morgan fingerprints")
    ap.add_argument("--fp-3d-num-confs", type=int, default=10, help="Conformers for 3D fingerprints (usr/usrcat)")
    ap.add_argument("--fp-3d-minimize", choices=["none", "uff", "mmff"], default="none", help="Optional minimization for 3D fingerprints")

    ap.add_argument("--no-images", action="store_true", help="Disable scaffold/cluster SVG generation")
    ap.add_argument(
        "--enriched-name",
        default=None,
        help="Name for enriched table (default: <stem>_series.csv). Use .parquet to write Parquet.",
    )
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
    analyze_clustering(
        args.input,
        outdir=clus_dir,
        smiles_col=args.smiles_col,
        similarity=args.similarity,
        fp=args.fp,
        radius=args.fp_radius,
        nbits=args.fp_nbits,
        use_chirality=args.fp_chirality,
        use_features=args.fp_features,
        num_confs=args.fp_3d_num_confs,
        minimize=args.fp_3d_minimize,
        top_images=top_images,
    )

<<<<<<< HEAD
    # Build enriched table by merging on Compound_ID
=======
    # Build enriched CSV by merging on Compound_ID
>>>>>>> 5d70153 (Core: prefer Calc_Canonical_SMILES + add Parquet table IO)
    df = read_table(args.input)
    id_col = detect_id_column(df)

    df_scaf = pd.read_csv(scaf_dir / "scaffold_assignments.csv")
    df_clus = pd.read_csv(clus_dir / "cluster_assignments.csv")

    keep_scaf = [id_col, "Scaffold_ID", "Scaffold_SMILES", "Scaffold_Size"]
    keep_clus = [id_col, "Cluster_ID", "Cluster_Size", "Cluster_Representative_ID", "Cluster_MedoidMeanSim"]

    df_enriched = df.merge(df_scaf[[c for c in keep_scaf if c in df_scaf.columns]], on=id_col, how="left")
    df_enriched = df_enriched.merge(df_clus[[c for c in keep_clus if c in df_clus.columns]], on=id_col, how="left")

    out_name = args.enriched_name or f"{in_path.stem}_series.csv"
    out_table = outdir / out_name
    write_table(df_enriched, str(out_table))

<<<<<<< HEAD
    # Task 3: sidecar metadata for the enriched table.
    try:
        # pick best smiles col (or user override) to record what series analytics used.
        smiles_used = args.smiles_col
        if not smiles_used:
            from molprop_toolkit.core.columns import detect_best_smiles_column

            smiles_used = detect_best_smiles_column(df.columns)

        write_run_metadata(
            tool="molprop-series",
            output_table_path=out_table,
            input_path=in_path,
            smiles_column=smiles_used,
            parameters={
                "similarity": float(args.similarity),
                "fp": str(args.fp),
                "fp_radius": int(args.fp_radius),
                "fp_nbits": int(args.fp_nbits),
                "fp_chirality": bool(args.fp_chirality),
                "fp_features": bool(args.fp_features),
                "fp_3d_num_confs": int(args.fp_3d_num_confs),
                "fp_3d_minimize": str(args.fp_3d_minimize),
                "no_images": bool(args.no_images),
            },
            artifacts={
                "series_metadata": str((outdir / "series_metadata.json").resolve()),
                "scaffold_assignments": str((scaf_dir / "scaffold_assignments.csv").resolve()),
                "cluster_assignments": str((clus_dir / "cluster_assignments.csv").resolve()),
            },
        )
    except Exception:
        pass

=======
>>>>>>> 5d70153 (Core: prefer Calc_Canonical_SMILES + add Parquet table IO)
    # Series-level reproducibility metadata (ties scaffold + cluster artifacts together)
    series_meta = {
        "tool": "molprop-series",
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "command": {"argv": sys.argv, "cwd": os.getcwd()},
        "input": {
            "csv_name": str(in_path.name),
            "csv_path": str(in_path.resolve()),
            "csv_sha256": _sha256_file(in_path),
            "csv_mtime": datetime.fromtimestamp(in_path.stat().st_mtime).isoformat(timespec="seconds") if in_path.exists() else None,
            "csv_size_bytes": int(in_path.stat().st_size) if in_path.exists() else None,
            "smiles_column": args.smiles_col,
        },
        "parameters": {
            "similarity": float(args.similarity),
            "fp": str(args.fp),
            "fp_radius": int(args.fp_radius),
            "fp_nbits": int(args.fp_nbits),
            "fp_chirality": bool(args.fp_chirality),
            "fp_features": bool(args.fp_features),
            "fp_3d_num_confs": int(args.fp_3d_num_confs),
            "fp_3d_minimize": str(args.fp_3d_minimize),
            "no_images": bool(args.no_images),
        },
        "artifacts": {
            "outdir": str(outdir.resolve()),
<<<<<<< HEAD
            "enriched_csv": str(out_table.resolve()),
=======
            "enriched_csv": str(out_csv.resolve()),
>>>>>>> 5d70153 (Core: prefer Calc_Canonical_SMILES + add Parquet table IO)
            "scaffolds_dir": str(scaf_dir.resolve()),
            "clusters_dir": str(clus_dir.resolve()),
            "scaffold_assignments": str((scaf_dir / "scaffold_assignments.csv").resolve()),
            "scaffold_summary": str((scaf_dir / "scaffold_summary.csv").resolve()),
            "cluster_assignments": str((clus_dir / "cluster_assignments.csv").resolve()),
            "cluster_summary": str((clus_dir / "cluster_summary.csv").resolve()),
            "cluster_metadata": str((clus_dir / "cluster_metadata.json").resolve()),
        },
    }
    (outdir / "series_metadata.json").write_text(json.dumps(series_meta, indent=2) + "\n", encoding="utf-8")

<<<<<<< HEAD
    print(f"\nEnriched table written to: {out_table}")
=======
    print(f"\nEnriched CSV written to: {out_csv}")
>>>>>>> 5d70153 (Core: prefer Calc_Canonical_SMILES + add Parquet table IO)
    print(f"Series metadata written to: {outdir / 'series_metadata.json'}")
    print("Tip: run reports/analyzers on the enriched CSV to include series columns.")


if __name__ == "__main__":
    main()

