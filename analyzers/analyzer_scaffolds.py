#!/usr/bin/env python3
"""══════════════════════════════════════════════════════════════════════════════
MURCKO SCAFFOLD ANALYZER (Task 13A)
═══════════════════════════════════════════════════════════════════════════════

Generates Murcko scaffold assignments and scaffold summary tables.

Outputs (written to --outdir):
- scaffold_assignments.csv  (compound-level)
- scaffold_summary.csv      (per scaffold)
- scaffold_outliers.csv     (best/worst per scaffold for chosen key metrics)
- optional scaffold_images/ (SVG depictions for top scaffolds)

Notes:
- This is intended for series triage and project communication.
- Scaffolds are computed from the "best" SMILES column (Calc_Canonical_SMILES
  preferred) using RDKit's Murcko scaffold definition.
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
    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold
except Exception as e:  # pragma: no cover
    raise SystemExit(
        "RDKit is required for scaffold analysis. Install via conda-forge: conda install -c conda-forge rdkit"
    ) from e

from calculators.depict import mol_from_smiles, sanitize_filename, save_svg


def _murcko_scaffold_smiles(smiles: str) -> str:
    if not smiles or not isinstance(smiles, str):
        return ""
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return ""
    try:
        scaf = MurckoScaffold.GetScaffoldForMol(m)
        if scaf is None or scaf.GetNumAtoms() == 0:
            return ""
        return Chem.MolToSmiles(scaf, canonical=True, isomericSmiles=False)
    except Exception:
        return ""


def _make_scaffold_ids(scaf_smiles: pd.Series) -> Tuple[pd.Series, Dict[str, str]]:
    counts = scaf_smiles.value_counts(dropna=False)
    # assign IDs by descending count then lexicographic for stability
    uniq = [s for s in counts.index.tolist() if isinstance(s, str) and s]
    mapping: Dict[str, str] = {s: f"SCAF_{i:04d}" for i, s in enumerate(uniq, 1)}
    ids = scaf_smiles.map(lambda s: mapping.get(s, ""))
    return ids, mapping


def _numeric_cols(df: pd.DataFrame, candidates: Sequence[str]) -> List[str]:
    out: List[str] = []
    for c in candidates:
        if c in df.columns and pd.api.types.is_numeric_dtype(df[c]):
            out.append(c)
    return out


def analyze_scaffolds(
    csv_path: str,
    outdir: str | Path,
    smiles_col: Optional[str] = None,
    top_images: int = 50,
    image_size: Tuple[int, int] = (260, 190),
) -> Dict[str, Path]:
    df = read_csv(csv_path)
    id_col = detect_id_column(df)
    smi_col = smiles_col or detect_best_smiles_column(df.columns)
    if not smi_col or smi_col not in df.columns:
        raise SystemExit("No SMILES column found. Expected one of: Calc_Canonical_SMILES, Calc_Base_SMILES, Canonical_SMILES, Input_Canonical_SMILES, SMILES")

    outdir_p = Path(outdir)
    outdir_p.mkdir(parents=True, exist_ok=True)

    # Compute scaffold
    scaff_smiles = df[smi_col].astype(str).map(_murcko_scaffold_smiles)
    df_out = df[[id_col, smi_col]].copy()
    df_out["Scaffold_SMILES"] = scaff_smiles

    scaffold_id, mapping = _make_scaffold_ids(df_out["Scaffold_SMILES"])
    df_out["Scaffold_ID"] = scaffold_id

    scaf_sizes = df_out["Scaffold_ID"].value_counts()
    df_out["Scaffold_Size"] = df_out["Scaffold_ID"].map(lambda x: int(scaf_sizes.get(x, 0)))

    # Add a few useful numeric context columns when present
    context_num = _numeric_cols(df, ["CNS_MPO", "MolWt", "LogP", "TPSA", "Sol_ESOL_LogS", "QED", "SA_Score"])
    for c in context_num:
        df_out[c] = df[c]

    # scaffold assignments
    p_assign = outdir_p / "scaffold_assignments.csv"
    df_out.to_csv(p_assign, index=False)

    # scaffold summary
    agg_cols = context_num
    g = df_out[df_out["Scaffold_ID"] != ""].groupby("Scaffold_ID", dropna=False)

    summary_rows: List[Dict[str, object]] = []
    for scaf_id, sub in g:
        scaf_smi = sub["Scaffold_SMILES"].iloc[0] if len(sub) else ""
        row: Dict[str, object] = {
            "Scaffold_ID": scaf_id,
            "Scaffold_SMILES": scaf_smi,
            "Count": int(len(sub)),
            "Example_Compounds": ", ".join(sub[id_col].astype(str).head(5).tolist()),
        }
        for c in agg_cols:
            vals = pd.to_numeric(sub[c], errors="coerce").dropna()
            if len(vals) == 0:
                row[f"{c}_mean"] = None
                row[f"{c}_median"] = None
            else:
                row[f"{c}_mean"] = float(vals.mean())
                row[f"{c}_median"] = float(vals.median())
        summary_rows.append(row)

    df_summary = pd.DataFrame(summary_rows)
    if "Count" in df_summary.columns:
        df_summary = df_summary.sort_values("Count", ascending=False)

    p_summary = outdir_p / "scaffold_summary.csv"
    df_summary.to_csv(p_summary, index=False)

    # outliers per scaffold (best/worst by selected metrics)
    metrics = [c for c in ["CNS_MPO", "QED", "SA_Score"] if c in df_out.columns]
    outlier_rows: List[Dict[str, object]] = []
    for scaf_id, sub in g:
        if len(sub) < 2:
            continue
        for m in metrics:
            v = pd.to_numeric(sub[m], errors="coerce")
            if v.dropna().empty:
                continue
            sub2 = sub.copy()
            sub2[m] = v
            best = sub2.sort_values(m, ascending=False).iloc[0]
            worst = sub2.sort_values(m, ascending=True).iloc[0]
            outlier_rows.append(
                {
                    "Scaffold_ID": scaf_id,
                    "Metric": m,
                    "Best_Compound": str(best[id_col]),
                    "Best_Value": best[m],
                    "Worst_Compound": str(worst[id_col]),
                    "Worst_Value": worst[m],
                }
            )

    df_outliers = pd.DataFrame(outlier_rows)
    p_outliers = outdir_p / "scaffold_outliers.csv"
    df_outliers.to_csv(p_outliers, index=False)

    # Optional scaffold images for top scaffolds
    if top_images and top_images > 0 and len(df_summary) > 0:
        img_dir = outdir_p / "scaffold_images"
        img_dir.mkdir(parents=True, exist_ok=True)
        for _, r in df_summary.head(int(top_images)).iterrows():
            scaf_id = str(r.get("Scaffold_ID") or "")
            smi = str(r.get("Scaffold_SMILES") or "")
            if not scaf_id or not smi:
                continue
            m = mol_from_smiles(smi)
            if m is None:
                continue
            fname = sanitize_filename(scaf_id) + ".svg"
            save_svg(img_dir / fname, m, size=image_size, legend=scaf_id)

    print(f"Scaffold artifacts written to: {outdir_p}")
    return {"assignments": p_assign, "summary": p_summary, "outliers": p_outliers}


def main() -> None:
    ap = argparse.ArgumentParser(description="Murcko scaffold analyzer (Task 13A)")
    ap.add_argument("input", help="Results CSV")
    ap.add_argument("--outdir", default=None, help="Output directory")
    ap.add_argument("--smiles-col", default=None, help="Override SMILES column")
    ap.add_argument("--top-images", type=int, default=50, help="Render SVG depictions for top N scaffolds (0 disables)")
    args = ap.parse_args()

    in_path = Path(args.input)
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    outdir = Path(args.outdir) if args.outdir else Path("series") / f"scaffolds_{in_path.stem}_{ts}"

    analyze_scaffolds(args.input, outdir=outdir, smiles_col=args.smiles_col, top_images=args.top_images)


if __name__ == "__main__":
    main()

