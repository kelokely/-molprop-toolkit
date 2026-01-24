#!/usr/bin/env python3
"""molprop-featurize: scalable feature generation for ML, SAR, and clustering.

This tool consumes either:
- a MolProp results CSV (recommended), or
- a SMILES file (name<tab>smiles / smiles<tab>name / whitespace-separated).

It produces a reproducible feature pack directory with:
- sparse hashed fingerprint + fragment + functional-group matrix (CSR NPZ),
- optional interpretable FG/fragment summary CSV,
- deep learning exports: packed-bit fingerprints, SMILES export, and graph tensors,
- metadata JSON recording every parameter and input hash.

By default, the structure basis is auto-detected and prefers Calc_Canonical_SMILES when present.

References:
- RDKit BRICS: https://www.rdkit.org/docs/source/rdkit.Chem.BRICS.html
- RDKit RECAP: https://www.rdkit.org/docs/source/rdkit.Chem.Recap.html
- RDKit fingerprint generators: https://www.rdkit.org/docs/source/rdkit.Chem.rdFingerprintGenerator.html
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from dataclasses import asdict
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from molprop_toolkit.core.columns import detect_best_smiles_column, detect_id_column
from molprop_toolkit.core.io import read_table
from molprop_toolkit.features.fragments import FragmentSpec, fragment_inventory
from molprop_toolkit.features.functional_groups import compute_functional_group_features
from molprop_toolkit.features.fingerprints_ml import MLFPConfig, compute_default_sparse_blocks, morgan_bitvect, pack_bitvect
from molprop_toolkit.features.graphs import GraphSpec, graphs_from_smiles_list
from molprop_toolkit.features.sparse_matrix import CSRBuilder, save_csr_npz


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_smiles_file(path: str) -> List[Tuple[str, str]]:
    items: List[Tuple[str, str]] = []
    with open(path, "r", encoding="utf-8") as f:
        for i, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t") if "\t" in line else line.split(",") if "," in line else line.split()
            if len(parts) == 1:
                smiles, name = parts[0], f"cmpd_{i}"
            else:
                a, b = parts[0], parts[1]
                looks_smiles_a = any(ch in a for ch in "[]=#@()")
                looks_smiles_b = any(ch in b for ch in "[]=#@()")
                if looks_smiles_a and not looks_smiles_b:
                    smiles, name = a, b
                elif looks_smiles_b and not looks_smiles_a:
                    name, smiles = a, b
                else:
                    name, smiles = a, b
            items.append((name, smiles))
    return items


def _require_rdkit():
    try:
        from rdkit import Chem  # type: ignore
        import rdkit  # type: ignore

        return Chem, rdkit
    except Exception as e:
        raise RuntimeError(
            "RDKit is required for molprop-featurize. Install via conda-forge (recommended)."
        ) from e


def select_numeric_columns(df: pd.DataFrame, explicit: Optional[str] = None) -> List[str]:
    if explicit:
        cols = [c.strip() for c in explicit.split(",") if c.strip()]
        return [c for c in cols if c in df.columns]

    # Sensible defaults: core physchem + a few common triage columns if present.
    preferred = [
        "MolWt",
        "LogP",
        "TPSA",
        "HBD",
        "HBA",
        "RotBonds",
        "Rings",
        "AromaticRings",
        "HeavyAtoms",
        "FractionCSP3",
        "QED",
        "cLogD_7.4",
        "pKa_Basic",
        "pKa_Acidic",
        "CNS_MPO",
        "Dev_Score",
    ]

    out = [c for c in preferred if c in df.columns]

    # Add any numeric Dev_* columns if present.
    for c in df.columns:
        if c.startswith("Dev_") and c not in out:
            out.append(c)

    return out


def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(prog="molprop-featurize", description="Generate ML-ready feature packs from MolProp tables")
    ap.add_argument("input", help="Input results CSV (recommended) or SMILES file")
    ap.add_argument("-o", "--outdir", default=None, help="Output directory (default: features/<stem>_features_<timestamp>)")

    ap.add_argument(
        "--smiles-col",
        default=None,
        help="SMILES column to use from a results table (default: auto; prefers Calc_Canonical_SMILES when present)",
    )
    ap.add_argument("--id-col", default=None, help="ID column to use (default: auto-detect)")

    ap.add_argument(
        "--fragments",
        choices=["none", "brics", "recap", "both"],
        default="both",
        help="Fragment vocab to generate",
    )
    ap.add_argument("--max-frags-per-mol", type=int, default=128, help="Max fragments per molecule per method")

    ap.add_argument("--morgan-radius", type=int, default=2)
    ap.add_argument("--morgan-nbits", type=int, default=4096)
    ap.add_argument("--ap-nbits", type=int, default=2048)
    ap.add_argument("--torsion-nbits", type=int, default=2048)
    ap.add_argument("--morgan-use-chirality", action="store_true", default=True)

    # Deep learning exports
    ap.add_argument("--dl", choices=["none", "fp", "smiles", "graphs", "all"], default="all")
    ap.add_argument("--dl-fp-nbits", type=int, default=2048, help="Bit fingerprint size for DL export")
    ap.add_argument(
        "--dl-fp-format",
        choices=["packed", "float32"],
        default="packed",
        help="DL fingerprint export format. packed saves uint8 packed bits; float32 saves dense 0/1 array (large).",
    )
    ap.add_argument("--graphs-include-hs", action="store_true", help="Include explicit hydrogens in graph export")

    ap.add_argument(
        "--numeric-cols",
        default=None,
        help="Comma-separated numeric columns from the input CSV to export as X_numeric (default: a reasonable set)",
    )

    ap.add_argument("--write-interpretable", action="store_true", default=True, help="Write interpretable FG/fragment summary")
    ap.add_argument("--no-interpretable", dest="write_interpretable", action="store_false")

    args = ap.parse_args(list(argv) if argv is not None else None)

    in_path = str(args.input)
    in_stem = Path(in_path).stem

    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    outdir = args.outdir or str(Path("features") / f"{in_stem}_features_{ts}")
    Path(outdir).mkdir(parents=True, exist_ok=True)

    # Record input hash for reproducibility.
    input_sha = sha256_file(in_path) if Path(in_path).is_file() else ""

    # Load data
    rows: List[Dict[str, object]] = []
    if in_path.lower().endswith((".csv", ".tsv")):
        df = read_table(in_path)

        smiles_col = args.smiles_col
        if not smiles_col or smiles_col not in df.columns:
            guess = detect_best_smiles_column(df.columns)
            if guess:
                smiles_col = guess

        id_col = args.id_col or detect_id_column(df.columns)
        if id_col not in df.columns:
            id_col = df.columns[0]

        # Keep original order; drop rows with missing SMILES.
        df2 = df.copy()
        df2[smiles_col] = df2[smiles_col].astype(str)
        df2[id_col] = df2[id_col].astype(str)

        for _, r in df2.iterrows():
            smi = str(r.get(smiles_col, "") or "").strip()
            cid = str(r.get(id_col, "") or "").strip()
            if not cid:
                cid = ""
            if not smi or smi.lower() in ("nan", "none"):
                continue
            rows.append({"Compound_ID": cid, "SMILES": smi, "_row": r.to_dict()})

        numeric_cols = select_numeric_columns(df2, explicit=args.numeric_cols)
    else:
        # SMILES input. Use RDKit preparation only to canonicalize input.
        items = parse_smiles_file(in_path)
        numeric_cols = []
        for cid, smi in items:
            rows.append({"Compound_ID": cid, "SMILES": smi, "_row": {}})

    if not rows:
        raise SystemExit("No valid structures found")

    Chem, rdkit = _require_rdkit()

    # Pass 1: build fragment vocab (BRICS + RECAP) and store per-molecule fragment lists for reuse.
    frag_spec = FragmentSpec(method=args.fragments, max_per_mol=int(args.max_frags_per_mol))

    brics_per_mol: List[List[str]] = []
    recap_per_mol: List[List[str]] = []

    brics_freq: Dict[str, int] = {}
    recap_freq: Dict[str, int] = {}

    for i, it in enumerate(rows):
        smi = str(it["SMILES"])
        m = Chem.MolFromSmiles(smi)
        if m is None:
            brics_per_mol.append([])
            recap_per_mol.append([])
            continue

        inv = fragment_inventory(m, frag_spec)
        br = inv.get("BRICS", [])
        rc = inv.get("RECAP", [])

        brics_per_mol.append(br)
        recap_per_mol.append(rc)

        for s in br:
            brics_freq[s] = brics_freq.get(s, 0) + 1
        for s in rc:
            recap_freq[s] = recap_freq.get(s, 0) + 1

        if (i + 1) % 5000 == 0:
            print(f"[molprop-featurize] pass1 fragments: {i+1}/{len(rows)}", file=sys.stderr)

    # Build fragment index maps.
    # We keep all unique fragments by default; callers can downstream-filter by frequency.
    brics_vocab = sorted(brics_freq.keys())
    recap_vocab = sorted(recap_freq.keys())

    # Fixed hashed fingerprint blocks.
    cfg = MLFPConfig(
        morgan_radius=int(args.morgan_radius),
        morgan_nbits=int(args.morgan_nbits),
        morgan_use_chirality=bool(args.morgan_use_chirality),
        atom_pair_nbits=int(args.ap_nbits),
        torsion_nbits=int(args.torsion_nbits),
        bit_nbits=int(args.dl_fp_nbits),
    )

    # FG block names are stable.
    # Build once from a single molecule so ordering is deterministic.
    fg_names: List[str] = []
    try:
        # Use a tiny molecule to initialize.
        m0 = Chem.MolFromSmiles("CC")
        fg_names = sorted(compute_functional_group_features(m0).keys()) if m0 is not None else []
    except Exception:
        fg_names = []

    # Assign feature blocks and offsets.
    blocks = []
    offset = 0

    blocks.append({"name": "FP_Morgan", "kind": "hashed_count", "size": int(cfg.morgan_nbits), "offset": offset})
    offset += int(cfg.morgan_nbits)

    blocks.append({"name": "FP_AtomPair", "kind": "hashed_count", "size": int(cfg.atom_pair_nbits), "offset": offset})
    offset += int(cfg.atom_pair_nbits)

    blocks.append({"name": "FP_Torsion", "kind": "hashed_count", "size": int(cfg.torsion_nbits), "offset": offset})
    offset += int(cfg.torsion_nbits)

    blocks.append({"name": "FG", "kind": "named_numeric", "size": len(fg_names), "offset": offset})
    fg_offset = offset
    offset += len(fg_names)

    blocks.append({"name": "Frag_BRICS", "kind": "named_presence", "size": len(brics_vocab), "offset": offset})
    brics_offset = offset
    offset += len(brics_vocab)

    blocks.append({"name": "Frag_RECAP", "kind": "named_presence", "size": len(recap_vocab), "offset": offset})
    recap_offset = offset
    offset += len(recap_vocab)

    n_features = offset

    brics_index = {s: i for i, s in enumerate(brics_vocab)}
    recap_index = {s: i for i, s in enumerate(recap_vocab)}

    # Prepare DL fingerprint storage
    dl_mode = (args.dl or "none").lower()
    do_dl = dl_mode != "none"
    do_dl_fp = do_dl and dl_mode in ("fp", "all")
    do_dl_smiles = do_dl and dl_mode in ("smiles", "all")
    do_dl_graphs = do_dl and dl_mode in ("graphs", "all")

    dl_fp_format = str(args.dl_fp_format)
    dl_bit_nbits = int(args.dl_fp_nbits)
    dl_nbytes = (dl_bit_nbits + 7) // 8

    dl_fp_packed = np.zeros((len(rows), dl_nbytes), dtype=np.uint8) if (do_dl_fp and dl_fp_format == "packed") else None
    dl_fp_dense = np.zeros((len(rows), dl_bit_nbits), dtype=np.float32) if (do_dl_fp and dl_fp_format == "float32") else None

    # Numeric features (dense)
    X_num = None
    if numeric_cols:
        X_num = np.zeros((len(rows), len(numeric_cols)), dtype=np.float32)

    # Interpretable output table
    interp_rows: List[Dict[str, object]] = []

    # Build CSR sparse matrix
    csr = CSRBuilder(n_rows=len(rows), n_cols=n_features)

    smiles_for_graphs: List[str] = []

    for i, it in enumerate(rows):
        cid = str(it.get("Compound_ID", ""))
        smi = str(it.get("SMILES", ""))

        m = Chem.MolFromSmiles(smi)
        if m is None:
            csr.add_row({})
            if X_num is not None:
                X_num[i, :] = np.nan
            if dl_fp_packed is not None:
                dl_fp_packed[i, :] = 0
            if dl_fp_dense is not None:
                dl_fp_dense[i, :] = 0
            smiles_for_graphs.append("")
            continue

        # Sparse blocks: hashed counts
        blocks_sparse = compute_default_sparse_blocks(m, cfg)

        row_feats: Dict[int, float] = {}

        # Morgan
        for k, v in blocks_sparse.get("Morgan", {}).items():
            row_feats[int(k) + blocks[0]["offset"]] = float(v)
        # AtomPair
        for k, v in blocks_sparse.get("AtomPair", {}).items():
            row_feats[int(k) + blocks[1]["offset"]] = float(v)
        # Torsion
        for k, v in blocks_sparse.get("Torsion", {}).items():
            row_feats[int(k) + blocks[2]["offset"]] = float(v)

        # FG features
        fg = compute_functional_group_features(m) if fg_names else {}
        for j, name in enumerate(fg_names):
            val = fg.get(name)
            if val is None:
                continue
            try:
                fv = float(val)
            except Exception:
                continue
            if fv != 0.0 and not np.isnan(fv):
                row_feats[fg_offset + j] = fv

        # Fragments as presence features
        for s in brics_per_mol[i]:
            idx = brics_index.get(s)
            if idx is not None:
                row_feats[brics_offset + int(idx)] = 1.0
        for s in recap_per_mol[i]:
            idx = recap_index.get(s)
            if idx is not None:
                row_feats[recap_offset + int(idx)] = 1.0

        csr.add_row(row_feats)

        # Numeric features
        if X_num is not None and numeric_cols:
            r = it.get("_row", {}) or {}
            for j, c in enumerate(numeric_cols):
                try:
                    X_num[i, j] = float(r.get(c))
                except Exception:
                    X_num[i, j] = np.nan

        # DL fingerprint
        if do_dl_fp:
            fp = morgan_bitvect(m, radius=int(args.morgan_radius), nbits=int(dl_bit_nbits), use_chirality=bool(args.morgan_use_chirality))
            if dl_fp_packed is not None:
                b = pack_bitvect(fp)
                # b length should be dl_nbytes
                bb = b if len(b) == dl_nbytes else b[:dl_nbytes].ljust(dl_nbytes, b"\x00")
                dl_fp_packed[i, :] = np.frombuffer(bb, dtype=np.uint8)
            if dl_fp_dense is not None:
                # Dense 0/1 array (large). Keep simple.
                s_bits = fp.ToBitString()
                dl_fp_dense[i, :] = np.fromiter((1.0 if ch == "1" else 0.0 for ch in s_bits), dtype=np.float32, count=dl_bit_nbits)

        # Graph smiles
        smiles_for_graphs.append(smi if do_dl_graphs else "")

        # Interpretable summary row
        if args.write_interpretable:
            ir: Dict[str, object] = {
                "Compound_ID": cid,
                "Canonical_SMILES": smi,
                "Frag_BRICS_Count": len(brics_per_mol[i]),
                "Frag_RECAP_Count": len(recap_per_mol[i]),
                "Frag_BRICS_TopK": ";".join(brics_per_mol[i][:10]),
                "Frag_RECAP_TopK": ";".join(recap_per_mol[i][:10]),
            }
            for name in fg_names:
                if name in fg:
                    ir[name] = fg.get(name)
            interp_rows.append(ir)

        if (i + 1) % 5000 == 0:
            print(f"[molprop-featurize] pass2 features: {i+1}/{len(rows)}", file=sys.stderr)

    # Save CSR matrix
    data, indices, indptr, shape = csr.to_arrays()
    save_csr_npz(str(Path(outdir) / "X_sparse_counts_csr.npz"), data=data, indices=indices, indptr=indptr, shape=shape)

    # Save numeric block
    if X_num is not None and numeric_cols:
        np.save(str(Path(outdir) / "X_numeric.npy"), X_num)
        with open(Path(outdir) / "X_numeric_columns.json", "w", encoding="utf-8") as f:
            json.dump({"columns": numeric_cols}, f, indent=2)

    # Save ID map
    ids_df = pd.DataFrame(
        {
            "row_index": list(range(len(rows))),
            "Compound_ID": [str(r.get("Compound_ID", "")) for r in rows],
            "Canonical_SMILES": [str(r.get("SMILES", "")) for r in rows],
        }
    )
    ids_df.to_csv(Path(outdir) / "ids.csv", index=False)

    # Save fragment vocabs + frequencies
    with open(Path(outdir) / "fragments_brics.json", "w", encoding="utf-8") as f:
        json.dump({"vocab": brics_vocab, "freq": brics_freq}, f)
    with open(Path(outdir) / "fragments_recap.json", "w", encoding="utf-8") as f:
        json.dump({"vocab": recap_vocab, "freq": recap_freq}, f)

    # Save interpretable summary
    if args.write_interpretable and interp_rows:
        pd.DataFrame(interp_rows).to_csv(Path(outdir) / "features_interpretable.csv", index=False)

    # DL exports
    if do_dl_smiles:
        ids_df[["Compound_ID", "Canonical_SMILES"]].to_csv(Path(outdir) / "dl_smiles.csv", index=False)

    if do_dl_fp:
        if dl_fp_packed is not None:
            np.save(Path(outdir) / "dl_fp_packed_uint8.npy", dl_fp_packed)
        if dl_fp_dense is not None:
            np.save(Path(outdir) / "dl_fp_dense_float32.npy", dl_fp_dense)

    if do_dl_graphs:
        gspec = GraphSpec(include_hs=bool(args.graphs_include_hs))
        graphs = graphs_from_smiles_list(smiles_for_graphs, gspec)
        np.savez_compressed(
            Path(outdir) / "dl_graphs.npz",
            atom_features=graphs["atom_features"],
            edge_index=graphs["edge_index"],
            edge_features=graphs["edge_features"],
            mol_atom_ptr=graphs["mol_atom_ptr"],
            mol_edge_ptr=graphs["mol_edge_ptr"],
        )

    # Metadata
    meta = {
        "created_at": ts,
        "input": {
            "path": os.path.abspath(in_path),
            "sha256": input_sha,
            "size_bytes": os.path.getsize(in_path) if os.path.exists(in_path) else None,
        },
        "structure_basis": "auto (prefers Calc_Canonical_SMILES)",
        "smiles_col_requested": args.smiles_col,
        "fragments": asdict(frag_spec),
        "fingerprints": asdict(cfg),
        "blocks": blocks,
        "fg_names": fg_names,
        "numeric_cols": numeric_cols,
        "dl": {
            "mode": dl_mode,
            "fp_nbits": dl_bit_nbits,
            "fp_format": dl_fp_format,
            "graphs_include_hs": bool(args.graphs_include_hs),
        },
        "versions": {
            "python": sys.version.split()[0],
            "rdkit": getattr(rdkit, "__version__", "unknown"),
            "numpy": np.__version__,
            "pandas": pd.__version__,
        },
    }

    with open(Path(outdir) / "features_metadata.json", "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2)

    print(str(Path(outdir).resolve()))


if __name__ == "__main__":
    main()

