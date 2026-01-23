#!/usr/bin/env python3
"""══════════════════════════════════════════════════════════════════════════════
SIMILARITY CLUSTERING ANALYZER (Task 13B)
═══════════════════════════════════════════════════════════════════════════════

Performs Morgan fingerprint similarity clustering (Butina) and writes:
- cluster_assignments.csv
- cluster_summary.csv
- optional cluster_representatives/ (SVG depictions for cluster medoids)

The main output expectation for Task 13B is a "Cluster_ID" assignment column.
This script produces that assignment and a few helpful companion columns.
"""

from __future__ import annotations

import argparse
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import pandas as pd

# Allow running directly without installing.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from molprop_toolkit.core import detect_best_smiles_column, detect_id_column, read_csv

try:
    from rdkit import Chem, DataStructs
    from rdkit.Chem import rdFingerprintGenerator
    from rdkit.ML.Cluster import Butina
except Exception as e:  # pragma: no cover
    raise SystemExit(
        "RDKit is required for clustering. Install via conda-forge: conda install -c conda-forge rdkit"
    ) from e

from calculators.depict import mol_from_smiles, sanitize_filename, save_svg


def _morgan_fps(mols: Sequence[Chem.Mol], radius: int = 2, nbits: int = 2048):
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=nbits)
    return [gen.GetFingerprint(m) for m in mols]


def _butina_clusters(fps, cutoff: float) -> List[Tuple[int, ...]]:
    # cutoff is a distance cutoff; if you want similarity threshold t, use cutoff = 1 - t
    dists: List[float] = []
    for i in range(1, len(fps)):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        dists.extend([1.0 - x for x in sims])
    clusters = Butina.ClusterData(dists, len(fps), cutoff, isDistData=True)
    # clusters are tuples of indices
    return [tuple(c) for c in clusters]


def _cluster_medoid(indices: Sequence[int], fps) -> Tuple[int, float]:
    # medoid = compound with highest mean similarity to others in the cluster
    if len(indices) == 1:
        return int(indices[0]), 1.0

    best_i = int(indices[0])
    best_mean = -1.0
    for i in indices:
        sims = []
        for j in indices:
            if i == j:
                continue
            sims.append(DataStructs.TanimotoSimilarity(fps[i], fps[j]))
        mean = float(sum(sims) / len(sims)) if sims else 0.0
        if mean > best_mean:
            best_mean = mean
            best_i = int(i)
    return best_i, best_mean


def analyze_clustering(
    csv_path: str,
    outdir: str | Path,
    smiles_col: Optional[str] = None,
    similarity: float = 0.7,
    radius: int = 2,
    nbits: int = 2048,
    top_images: int = 50,
    image_size: Tuple[int, int] = (260, 190),
) -> Dict[str, Path]:
    if similarity <= 0 or similarity >= 1:
        raise SystemExit("--similarity must be between 0 and 1 (e.g., 0.7)")

    df = read_csv(csv_path)
    id_col = detect_id_column(df)
    smi_col = smiles_col or detect_best_smiles_column(df.columns)
    if not smi_col or smi_col not in df.columns:
        raise SystemExit("No SMILES column found for clustering")

    outdir_p = Path(outdir)
    outdir_p.mkdir(parents=True, exist_ok=True)

    # Build mol list and keep mapping back to rows.
    mols: List[Chem.Mol] = []
    row_idx: List[int] = []
    for i, smi in enumerate(df[smi_col].astype(str).tolist()):
        m = Chem.MolFromSmiles(smi)
        if m is None:
            continue
        mols.append(m)
        row_idx.append(i)

    if not mols:
        raise SystemExit("No valid molecules found for clustering")

    fps = _morgan_fps(mols, radius=radius, nbits=nbits)
    cutoff = 1.0 - float(similarity)
    clusters = _butina_clusters(fps, cutoff=cutoff)

    # Assign cluster IDs (1..K) by descending size.
    clusters_sorted = sorted(clusters, key=lambda c: len(c), reverse=True)
    cluster_id_for_local: Dict[int, int] = {}
    cluster_size_for_local: Dict[int, int] = {}
    cluster_rep_for_local: Dict[int, str] = {}
    cluster_rep_sim_for_local: Dict[int, float] = {}

    for cid, members in enumerate(clusters_sorted, 1):
        for mi in members:
            cluster_id_for_local[int(mi)] = cid
            cluster_size_for_local[int(mi)] = int(len(members))

        medoid_i, medoid_mean = _cluster_medoid(members, fps)
        rep_global_row = row_idx[int(medoid_i)]
        rep_id = str(df.iloc[rep_global_row][id_col])

        for mi in members:
            cluster_rep_for_local[int(mi)] = rep_id
            cluster_rep_sim_for_local[int(mi)] = float(medoid_mean)

    # Create assignment table aligned to original df
    df_assign = df[[id_col, smi_col]].copy()

    # Use nullable dtypes so we can set ints/floats without pandas dtype issues.
    df_assign["Cluster_ID"] = pd.Series([pd.NA] * len(df_assign), dtype="Int64")
    df_assign["Cluster_Size"] = pd.Series([pd.NA] * len(df_assign), dtype="Int64")
    df_assign["Cluster_Representative_ID"] = pd.Series(["" for _ in range(len(df_assign))], dtype="string")
    df_assign["Cluster_MedoidMeanSim"] = pd.Series([pd.NA] * len(df_assign), dtype="Float64")

    for local_i, global_i in enumerate(row_idx):
        cid = cluster_id_for_local.get(local_i)
        if cid is None:
            continue
        df_assign.at[global_i, "Cluster_ID"] = int(cid)
        df_assign.at[global_i, "Cluster_Size"] = int(cluster_size_for_local.get(local_i, 0))
        df_assign.at[global_i, "Cluster_Representative_ID"] = str(cluster_rep_for_local.get(local_i, ""))
        df_assign.at[global_i, "Cluster_MedoidMeanSim"] = float(cluster_rep_sim_for_local.get(local_i, 0.0))

    # Add a few context numeric columns
    context_num = [c for c in ["CNS_MPO", "MolWt", "LogP", "TPSA", "Sol_ESOL_LogS", "QED", "SA_Score"] if c in df.columns]
    for c in context_num:
        df_assign[c] = df[c]

    p_assign = outdir_p / "cluster_assignments.csv"
    df_assign.to_csv(p_assign, index=False)

    # Cluster summary
    g = df_assign[df_assign["Cluster_ID"].notna()].groupby("Cluster_ID")
    rows: List[Dict[str, object]] = []
    for cid, sub in g:
        rep = sub["Cluster_Representative_ID"].iloc[0] if len(sub) else ""
        row: Dict[str, object] = {
            "Cluster_ID": int(cid),
            "Count": int(len(sub)),
            "Representative": rep,
            "Example_Compounds": ", ".join(sub[id_col].astype(str).head(5).tolist()),
        }
        for c in context_num:
            vals = pd.to_numeric(sub[c], errors="coerce").dropna()
            if vals.empty:
                row[f"{c}_mean"] = None
                row[f"{c}_median"] = None
            else:
                row[f"{c}_mean"] = float(vals.mean())
                row[f"{c}_median"] = float(vals.median())
        rows.append(row)

    df_summary = pd.DataFrame(rows)
    if not df_summary.empty and "Count" in df_summary.columns:
        df_summary = df_summary.sort_values("Count", ascending=False)
    p_summary = outdir_p / "cluster_summary.csv"
    df_summary.to_csv(p_summary, index=False)

    # Optional cluster representative images
    if top_images and top_images > 0 and len(df_summary) > 0:
        img_dir = outdir_p / "cluster_representatives"
        img_dir.mkdir(parents=True, exist_ok=True)
        reps = df_summary.head(int(top_images))["Representative"].astype(str).tolist()
        rep_set = set(reps)
        # map representative id -> smiles
        rep_smiles: Dict[str, str] = {}
        for _, r in df_assign.iterrows():
            rid = str(r.get("Cluster_Representative_ID") or "")
            if rid and rid in rep_set and rid not in rep_smiles:
                rep_smiles[rid] = str(r.get(smi_col) or "")

        for rid, smi in rep_smiles.items():
            m = mol_from_smiles(smi)
            if m is None:
                continue
            fname = sanitize_filename(rid) + ".svg"
            save_svg(img_dir / fname, m, size=image_size, legend=rid)

    print(f"Clustering artifacts written to: {outdir_p}")
    return {"assignments": p_assign, "summary": p_summary}


def main() -> None:
    ap = argparse.ArgumentParser(description="Similarity clustering (Butina) analyzer (Task 13B)")
    ap.add_argument("input", help="Results CSV")
    ap.add_argument("--outdir", default=None, help="Output directory")
    ap.add_argument("--smiles-col", default=None, help="Override SMILES column")
    ap.add_argument("--similarity", type=float, default=0.7, help="Tanimoto similarity threshold for clustering (default 0.7)")
    ap.add_argument("--radius", type=int, default=2, help="Morgan radius")
    ap.add_argument("--nbits", type=int, default=2048, help="Fingerprint size")
    ap.add_argument("--top-images", type=int, default=50, help="Render SVG depictions for top N cluster representatives (0 disables)")
    args = ap.parse_args()

    in_path = Path(args.input)
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    outdir = Path(args.outdir) if args.outdir else Path("series") / f"clusters_{in_path.stem}_{ts}"

    analyze_clustering(
        args.input,
        outdir=outdir,
        smiles_col=args.smiles_col,
        similarity=args.similarity,
        radius=args.radius,
        nbits=args.nbits,
        top_images=args.top_images,
    )


if __name__ == "__main__":
    main()

