#!/usr/bin/env python3
"""Similarity clustering analyzer.

Performs fingerprint-based similarity clustering (Butina) and writes:
- cluster_assignments.csv
- cluster_summary.csv
- optional cluster_representatives/ (SVG depictions for cluster medoids)

Supported fingerprints:
- 2D (bit-vector, Tanimoto): Morgan (ECFP-style), RDKit topological, MACCS keys, atom-pair (hashed), torsion (hashed)
- 3D (vector, distance-based): USR, USRCAT

For 2D bit-vector fingerprints, a similarity threshold `t` is converted to a Butina distance cutoff `1 - t`.
For USR/USRCAT vectors, we define similarity as `1/(1+euclidean_distance)` and convert the threshold to a distance cutoff
`(1/t) - 1`.

The main output is a per-compound `Cluster_ID` column plus a few companion fields:
- Cluster_Size
- Cluster_Representative_ID (medoid)
- Cluster_MedoidMeanSim (mean similarity of the medoid to cluster members)
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys
from datetime import datetime
from importlib import metadata as importlib_metadata
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

# Allow running directly without installing.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from molprop_toolkit.core import detect_best_smiles_column, detect_id_column, read_table
from molprop_toolkit.fingerprints import FingerprintSpec, FingerprintKind, fingerprint_from_mol
from molprop_toolkit.schema import load_schema

try:
    from rdkit import Chem, DataStructs
    from rdkit.ML.Cluster import Butina
except Exception as e:  # pragma: no cover
    raise SystemExit(
        "RDKit is required for clustering. Install via conda-forge: conda install -c conda-forge rdkit"
    ) from e

from calculators.depict import mol_from_smiles, sanitize_filename, save_svg


def _sha256_file(path: Path, *, max_bytes: int = 200 * 1024 * 1024) -> Optional[str]:
    """Return sha256 hex digest of a file, or None if file is too large or unreadable."""

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


def _pairwise_distances_tanimoto(fps) -> List[float]:
    """Return condensed distance list (len=n*(n-1)/2) where dist=1-tanimoto."""

    dists: List[float] = []
    for i in range(1, len(fps)):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        dists.extend([1.0 - float(x) for x in sims])
    return dists


def _pairwise_distances_euclidean(x: np.ndarray) -> List[float]:
    """Return condensed Euclidean distance list for vectors."""

    dists: List[float] = []
    for i in range(1, x.shape[0]):
        diff = x[:i, :] - x[i, :]
        di = np.sqrt(np.sum(diff * diff, axis=1))
        dists.extend([float(v) for v in di.tolist()])
    return dists


def _butina_cluster_from_distances(dists: Sequence[float], n: int, cutoff: float) -> List[Tuple[int, ...]]:
    clusters = Butina.ClusterData(list(dists), n, float(cutoff), isDistData=True)
    return [tuple(c) for c in clusters]


def _usr_similarity_from_distance(dist: float) -> float:
    return 1.0 / (1.0 + float(dist))


def _cluster_medoid_bit(indices: Sequence[int], fps) -> Tuple[int, float]:
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
            sims.append(float(DataStructs.TanimotoSimilarity(fps[int(i)], fps[int(j)])))
        mean = float(sum(sims) / len(sims)) if sims else 0.0
        if mean > best_mean:
            best_mean = mean
            best_i = int(i)
    return best_i, best_mean


def _cluster_medoid_vec(indices: Sequence[int], x: np.ndarray) -> Tuple[int, float]:
    if len(indices) == 1:
        return int(indices[0]), 1.0

    idx = np.array([int(i) for i in indices], dtype=int)
    best_i = int(idx[0])
    best_mean = -1.0

    for i in idx:
        others = idx[idx != i]
        if len(others) == 0:
            continue
        diff = x[others, :] - x[i, :]
        dist = np.sqrt(np.sum(diff * diff, axis=1))
        sim = 1.0 / (1.0 + dist)
        mean = float(np.mean(sim))
        if mean > best_mean:
            best_mean = mean
            best_i = int(i)

    return best_i, best_mean


def analyze_clustering(
    csv_path: str,
    outdir: str | Path,
    *,
    smiles_col: Optional[str] = None,
    similarity: float = 0.7,
    fp: FingerprintKind = "morgan",
    radius: int = 2,
    nbits: int = 2048,
    use_chirality: bool = False,
    use_features: bool = False,
    num_confs: int = 10,
    minimize: str = "none",
    top_images: int = 50,
    image_size: Tuple[int, int] = (260, 190),
) -> Dict[str, Path]:
    if similarity <= 0 or similarity >= 1:
        raise SystemExit("--similarity must be between 0 and 1 (e.g., 0.7)")

    df = read_table(csv_path)
    id_col = detect_id_column(df)
    smi_col = smiles_col or detect_best_smiles_column(df.columns)
    if not smi_col or smi_col not in df.columns:
        raise SystemExit("No SMILES column found for clustering")

    outdir_p = Path(outdir)
    outdir_p.mkdir(parents=True, exist_ok=True)

    spec = FingerprintSpec(
        kind=fp,
        radius=int(radius),
        nbits=int(nbits),
        use_chirality=bool(use_chirality),
        use_features=bool(use_features),
        num_confs=int(num_confs),
        minimize=minimize,  # type: ignore[arg-type]
    )

    # Build fingerprints and keep mapping back to original df rows.
    fps: List[object] = []
    row_idx: List[int] = []
    invalid_smiles = 0
    fp_failures = 0

    for i, smi in enumerate(df[smi_col].astype(str).tolist()):
        m = Chem.MolFromSmiles(smi)
        if m is None:
            invalid_smiles += 1
            continue
        try:
            fp_val = fingerprint_from_mol(m, spec)
        except Exception:
            fp_failures += 1
            continue
        fps.append(fp_val)
        row_idx.append(i)

    if not fps:
        raise SystemExit("No valid molecules/fingerprints found for clustering")

    # Cluster
    is_vec = fp in ("usr", "usrcat")
    if is_vec:
        x = np.asarray(fps, dtype=float)
        dists = _pairwise_distances_euclidean(x)
        cutoff = (1.0 / float(similarity)) - 1.0
        clusters = _butina_cluster_from_distances(dists, n=len(fps), cutoff=cutoff)
    else:
        dists = _pairwise_distances_tanimoto(fps)
        cutoff = 1.0 - float(similarity)
        clusters = _butina_cluster_from_distances(dists, n=len(fps), cutoff=cutoff)

    # Reproducibility metadata
    clusters_sorted_preview = sorted(clusters, key=lambda c: len(c), reverse=True)

    # Versions / environment
    rdkit_version = None
    try:  # pragma: no cover
        import rdkit  # type: ignore

        rdkit_version = getattr(rdkit, "__version__", None)
    except Exception:
        rdkit_version = None

    pkg_version = None
    try:  # pragma: no cover
        pkg_version = importlib_metadata.version("molprop-toolkit")
    except Exception:
        pkg_version = None

    schema_version = None
    try:  # pragma: no cover
        schema_version = str(load_schema().get("schema_version"))
    except Exception:
        schema_version = None

    git_sha = None
    try:  # pragma: no cover
        repo_root = Path(__file__).resolve().parents[1]
        git_sha = (
            subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], cwd=repo_root)
            .decode("utf-8")
            .strip()
        ) or None
    except Exception:
        git_sha = None

    distance_metric = "euclidean" if is_vec else "1 - tanimoto"
    similarity_definition = "1/(1+euclidean_distance)" if is_vec else "tanimoto"

    meta = {
        "tool": "analyzer_clustering",
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "molprop_toolkit_version": pkg_version,
        "schema_version": schema_version,
        "rdkit_version": rdkit_version,
        "git_commit": git_sha,
        "runtime": {
            "python_version": sys.version,
            "python_executable": sys.executable,
            "platform": platform.platform(),
            "machine": platform.machine(),
            "numpy_version": getattr(np, "__version__", None),
            "pandas_version": getattr(pd, "__version__", None),
        },
        "command": {
            "argv": sys.argv,
            "cwd": os.getcwd(),
        },
        "input": {
            "csv_name": str(Path(csv_path).name),
            "csv_path": str(Path(csv_path).resolve()),
            "csv_sha256": _sha256_file(Path(csv_path)),
            "csv_mtime": datetime.fromtimestamp(Path(csv_path).stat().st_mtime).isoformat(timespec="seconds")
            if Path(csv_path).exists()
            else None,
            "csv_size_bytes": int(Path(csv_path).stat().st_size) if Path(csv_path).exists() else None,
            "smiles_column": str(smi_col),
            "id_column": str(id_col),
        },
        "parameters": {
            "similarity_threshold": float(similarity),
            "butina_distance_cutoff": float(cutoff),
            "distance_metric": distance_metric,
            "similarity_definition": similarity_definition,
            "fingerprint": {
                "kind": spec.kind,
                "radius": int(spec.radius),
                "nbits": int(spec.nbits),
                "use_chirality": bool(spec.use_chirality),
                "use_features": bool(spec.use_features),
                "num_confs": int(spec.num_confs),
                "minimize": str(spec.minimize),
            },
        },
        "counts": {
            "rows_in_input": int(len(df)),
            "molecules_clustered": int(len(fps)),
            "invalid_smiles": int(invalid_smiles),
            "fingerprint_failures": int(fp_failures),
        },
        "clusters": {
            "count": int(len(clusters)),
            "sizes_top10": [int(len(c)) for c in clusters_sorted_preview[:10]],
        },
    }

    (outdir_p / "cluster_metadata.json").write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")

    # Assign cluster IDs (1..K) by descending size.
    clusters_sorted = clusters_sorted_preview
    cluster_id_for_local: Dict[int, int] = {}
    cluster_size_for_local: Dict[int, int] = {}
    cluster_rep_for_local: Dict[int, str] = {}
    cluster_rep_sim_for_local: Dict[int, float] = {}

    for cid, members in enumerate(clusters_sorted, 1):
        for mi in members:
            cluster_id_for_local[int(mi)] = cid
            cluster_size_for_local[int(mi)] = int(len(members))

        if is_vec:
            medoid_i, medoid_mean = _cluster_medoid_vec(members, x)  # type: ignore[arg-type]
        else:
            medoid_i, medoid_mean = _cluster_medoid_bit(members, fps)

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
    ap = argparse.ArgumentParser(description="Similarity clustering (Butina) analyzer")
    ap.add_argument("input", help="Results CSV")
    ap.add_argument("--outdir", default=None, help="Output directory")
    ap.add_argument("--smiles-col", default=None, help="Override SMILES column")
    ap.add_argument("--similarity", type=float, default=0.7, help="Similarity threshold for clustering (default 0.7)")

    ap.add_argument(
        "--fp",
        default="morgan",
        choices=["morgan", "rdkit", "maccs", "atom_pair", "torsion", "usr", "usrcat"],
        help="Fingerprint kind for clustering",
    )

    # Keep older flags as aliases.
    ap.add_argument("--fp-radius", "--radius", dest="fp_radius", type=int, default=2, help="Morgan radius")
    ap.add_argument("--fp-nbits", "--nbits", dest="fp_nbits", type=int, default=2048, help="Fingerprint size (where applicable)")
    ap.add_argument("--fp-chirality", action="store_true", help="Use chirality in Morgan fingerprints")
    ap.add_argument("--fp-features", action="store_true", help="Use feature-invariants in Morgan fingerprints")
    ap.add_argument("--fp-3d-num-confs", type=int, default=10, help="Conformers for 3D fingerprints (usr/usrcat)")
    ap.add_argument("--fp-3d-minimize", choices=["none", "uff", "mmff"], default="none", help="Optional minimization for 3D fingerprints")

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
        fp=args.fp,
        radius=args.fp_radius,
        nbits=args.fp_nbits,
        use_chirality=args.fp_chirality,
        use_features=args.fp_features,
        num_confs=args.fp_3d_num_confs,
        minimize=args.fp_3d_minimize,
        top_images=args.top_images,
    )


if __name__ == "__main__":
    main()

