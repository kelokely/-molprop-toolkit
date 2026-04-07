#!/usr/bin/env python3
"""molprop-featurize: scalable feature generation for ML, SAR, and clustering.

This tool consumes either:
- a MolProp results table (CSV/TSV/Parquet; recommended), or
- a SMILES file (name<tab>smiles / smiles<tab>name / whitespace-separated).

It produces a reproducible feature pack directory with:
- sparse hashed fingerprint + fragment + functional-group matrix (CSR NPZ),
- optional interpretable FG/fragment summary table (CSV/Parquet),
- optional deep-learning exports: packed-bit fingerprints, SMILES export, and graph tensors,
- metadata JSON recording parameters, input hash, and the selected structure-of-record column.

By default, the structure basis is auto-detected and prefers Calc_Canonical_SMILES when present.

High-end ML workflows often need *stable columns across runs*. The `--spec` option lets you freeze
fragment vocabularies and block layout to avoid train/test feature mismatch.

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
from typing import Any, Dict, List, Optional, Sequence, Tuple

# Allow running directly without installing the package.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd

from molprop_toolkit.core.columns import detect_best_smiles_column, detect_id_column
from molprop_toolkit.core.io import read_table, write_table
from molprop_toolkit.features.fragments import FragmentSpec, fragment_inventory
from molprop_toolkit.features.functional_groups import compute_functional_group_features
from molprop_toolkit.features.fingerprints_ml import MLFPConfig, compute_default_sparse_blocks, morgan_bitvect, pack_bitvect
from molprop_toolkit.features.graphs import GraphSpec, graphs_from_smiles_list
from molprop_toolkit.features.sparse_matrix import CSRBuilder, save_csr_npz


TABLE_EXTS = {".csv", ".tsv", ".tab", ".parquet", ".pq"}


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _json_load(path: str | Path) -> Dict[str, Any]:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def _json_dump(obj: Any, path: str | Path) -> None:
    Path(path).write_text(json.dumps(obj, indent=2) + "\n", encoding="utf-8")


def _looks_like_table(path: str) -> bool:
    return Path(path).suffix.lower() in TABLE_EXTS


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


def _build_blocks(*, cfg: MLFPConfig, fg_names: List[str], brics_vocab: List[str], recap_vocab: List[str]) -> Tuple[List[Dict[str, Any]], int, int, int]:
    """Return (blocks, fg_offset, brics_offset, recap_offset)."""

    blocks: List[Dict[str, Any]] = []
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

    return blocks, fg_offset, brics_offset, recap_offset


def _load_spec(path: str) -> Dict[str, Any]:
    spec = _json_load(path)
    if not isinstance(spec, dict) or "spec_version" not in spec:
        raise SystemExit(f"Invalid spec file (missing spec_version): {path}")
    return spec


def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(prog="molprop-featurize", description="Generate ML-ready feature packs from MolProp tables")
    ap.add_argument("input", help="Input results table (CSV/TSV/Parquet; recommended) or SMILES file")
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
        help="Comma-separated numeric columns from the input table to export as X_numeric (default: a reasonable set)",
    )

    ap.add_argument("--write-interpretable", action="store_true", default=True, help="Write interpretable FG/fragment summary")
    ap.add_argument("--no-interpretable", dest="write_interpretable", action="store_false")

    # Cross-run stability
    ap.add_argument(
        "--spec",
        default=None,
        help=(
            "Optional feature spec JSON to freeze fragment vocab + block layout across runs. "
            "If provided, featurize will reuse vocab/blocks from the spec and ignore conflicting flags."
        ),
    )
    ap.add_argument(
        "--no-write-spec",
        dest="write_spec",
        action="store_false",
        default=True,
        help="Do not write features_spec.json to the output directory.",
    )

    ap.add_argument(
        "--id-map-format",
        choices=["csv", "parquet", "both"],
        default="csv",
        help="Format for row→ID map output (ids.csv / ids.parquet). Parquet requires pyarrow.",
    )
    ap.add_argument(
        "--interpretable-format",
        choices=["csv", "parquet", "both"],
        default="csv",
        help="Format for interpretable summary output. Parquet requires pyarrow.",
    )

    args = ap.parse_args(list(argv) if argv is not None else None)

    # Optional: load a frozen spec to keep fragment vocab / block layout stable across runs.
    spec_path = str(args.spec) if args.spec else None
    spec: Optional[Dict[str, Any]] = _load_spec(spec_path) if spec_path else None
    spec_numeric_cols: Optional[List[str]] = list(spec.get("numeric_cols", [])) if spec is not None else None

    in_path = str(args.input)
    in_stem = Path(in_path).stem

    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    outdir = args.outdir or str(Path("features") / f"{in_stem}_features_{ts}")
    Path(outdir).mkdir(parents=True, exist_ok=True)

    # Record input hash for reproducibility.
    input_sha = sha256_file(in_path) if Path(in_path).is_file() else ""

    # Load data
    rows: List[Dict[str, object]] = []
    table_input = _looks_like_table(in_path)
    table_format: Optional[str] = None

    smiles_col_selected: Optional[str] = None
    id_col_selected: Optional[str] = None

    numeric_cols: List[str] = []

    if table_input:
        table_format = Path(in_path).suffix.lower().lstrip(".")
        df = read_table(in_path)

        smiles_col = args.smiles_col
        if not smiles_col or smiles_col not in df.columns:
            guess = detect_best_smiles_column(df.columns)
            if guess:
                smiles_col = guess

        if not smiles_col or smiles_col not in df.columns:
            raise SystemExit(
                "Could not determine a SMILES column. Provide --smiles-col explicitly. "
                "(Expected one of Calc_Canonical_SMILES / Calc_Base_SMILES / Canonical_SMILES / SMILES.)"
            )

        id_col = args.id_col or detect_id_column(df.columns)
        if not id_col or id_col not in df.columns:
            id_col = str(df.columns[0])

        smiles_col_selected = str(smiles_col)
        id_col_selected = str(id_col)

        # Keep original order; drop rows with missing SMILES.
        df2 = df.copy()
        df2[smiles_col_selected] = df2[smiles_col_selected].astype(str)
        df2[id_col_selected] = df2[id_col_selected].astype(str)

        # Numeric columns can be frozen via --spec to keep train/test feature layouts stable.
        if spec_numeric_cols is not None:
            numeric_cols = [c for c in spec_numeric_cols if c in df2.columns]
        else:
            numeric_cols = select_numeric_columns(df2, explicit=args.numeric_cols)

        for _, r in df2.iterrows():
            smi = str(r.get(smiles_col_selected, "") or "").strip()
            cid = str(r.get(id_col_selected, "") or "").strip()
            if not smi or smi.lower() in ("nan", "none"):
                continue
            num_row = {c: r.get(c) for c in numeric_cols} if numeric_cols else {}
            rows.append({"Compound_ID": cid, "SMILES": smi, "_row": num_row})
    else:
        # SMILES input.
        items = parse_smiles_file(in_path)
        for cid, smi in items:
            rows.append({"Compound_ID": cid, "SMILES": smi, "_row": {}})

    if not rows:
        raise SystemExit("No valid structures found")

    # RDKit required for all featurization paths.
    Chem, rdkit = _require_rdkit()

    # Resolve config (from spec if provided)
    if spec is not None:
        frag_spec = FragmentSpec(**spec.get("fragments", {}))
        cfg = MLFPConfig(**spec.get("fingerprints", {}))
        fg_names = list(spec.get("fg_names", []))
        brics_vocab = list(spec.get("fragment_vocab", {}).get("BRICS", []))
        recap_vocab = list(spec.get("fragment_vocab", {}).get("RECAP", []))
    else:
        frag_spec = FragmentSpec(method=args.fragments, max_per_mol=int(args.max_frags_per_mol))
        cfg = MLFPConfig(
            morgan_radius=int(args.morgan_radius),
            morgan_nbits=int(args.morgan_nbits),
            morgan_use_chirality=bool(args.morgan_use_chirality),
            atom_pair_nbits=int(args.ap_nbits),
            torsion_nbits=int(args.torsion_nbits),
            bit_nbits=int(args.dl_fp_nbits),
        )

        # FG block names are stable. Build once so ordering is deterministic.
        fg_names = []
        try:
            m0 = Chem.MolFromSmiles("CC")
            fg_names = sorted(compute_functional_group_features(m0).keys()) if m0 is not None else []
        except Exception:
            fg_names = []

        brics_vocab = []
        recap_vocab = []

    # Pass 1: compute fragments per molecule. If no spec, this pass also defines the vocab.
    brics_per_mol: List[List[str]] = []
    recap_per_mol: List[List[str]] = []

    brics_freq: Dict[str, int] = {}
    recap_freq: Dict[str, int] = {}

    unknown_brics = 0
    unknown_recap = 0

    brics_allowed = set(brics_vocab) if spec is not None else None
    recap_allowed = set(recap_vocab) if spec is not None else None

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
            if brics_allowed is not None and s not in brics_allowed:
                unknown_brics += 1
                continue
            brics_freq[s] = brics_freq.get(s, 0) + 1
        for s in rc:
            if recap_allowed is not None and s not in recap_allowed:
                unknown_recap += 1
                continue
            recap_freq[s] = recap_freq.get(s, 0) + 1

        if (i + 1) % 5000 == 0:
            print(f"[molprop-featurize] pass1 fragments: {i+1}/{len(rows)}", file=sys.stderr)

    if spec is None:
        # Keep all unique fragments by default; callers can downstream-filter by frequency.
        brics_vocab = sorted(brics_freq.keys())
        recap_vocab = sorted(recap_freq.keys())

    # Assign feature blocks and offsets.
    blocks, fg_offset, brics_offset, recap_offset = _build_blocks(
        cfg=cfg, fg_names=fg_names, brics_vocab=brics_vocab, recap_vocab=recap_vocab
    )
    n_features = blocks[-1]["offset"] + blocks[-1]["size"] if blocks else 0

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
    csr = CSRBuilder(n_rows=len(rows), n_cols=int(n_features))

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
            row_feats[int(k) + int(blocks[0]["offset"])] = float(v)
        # AtomPair
        for k, v in blocks_sparse.get("AtomPair", {}).items():
            row_feats[int(k) + int(blocks[1]["offset"])] = float(v)
        # Torsion
        for k, v in blocks_sparse.get("Torsion", {}).items():
            row_feats[int(k) + int(blocks[2]["offset"])] = float(v)

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
                row_feats[int(fg_offset) + j] = fv

        # Fragments as presence features (respect spec vocab if provided)
        for s in brics_per_mol[i]:
            idx = brics_index.get(s)
            if idx is not None:
                row_feats[int(brics_offset) + int(idx)] = 1.0
        for s in recap_per_mol[i]:
            idx = recap_index.get(s)
            if idx is not None:
                row_feats[int(recap_offset) + int(idx)] = 1.0

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
            fp = morgan_bitvect(
                m,
                radius=int(cfg.morgan_radius),
                nbits=int(dl_bit_nbits),
                use_chirality=bool(cfg.morgan_use_chirality),
            )
            if dl_fp_packed is not None:
                b = pack_bitvect(fp)
                bb = b if len(b) == dl_nbytes else b[:dl_nbytes].ljust(dl_nbytes, b"\x00")
                dl_fp_packed[i, :] = np.frombuffer(bb, dtype=np.uint8)
            if dl_fp_dense is not None:
                s_bits = fp.ToBitString()
                dl_fp_dense[i, :] = np.fromiter((1.0 if ch == "1" else 0.0 for ch in s_bits), dtype=np.float32, count=dl_bit_nbits)

        # Graph smiles
        smiles_for_graphs.append(smi if do_dl_graphs else "")

        # Interpretable summary row
        if args.write_interpretable:
            ir: Dict[str, object] = {
                "Compound_ID": cid,
                "Canonical_SMILES": smi,
                "SMILES_Source_Column": smiles_col_selected or "SMILES",
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
        _json_dump({"columns": numeric_cols}, Path(outdir) / "X_numeric_columns.json")

    # Save ID map
    ids_df = pd.DataFrame(
        {
            "row_index": list(range(len(rows))),
            "Compound_ID": [str(r.get("Compound_ID", "")) for r in rows],
            "Canonical_SMILES": [str(r.get("SMILES", "")) for r in rows],
            "SMILES_Source_Column": [smiles_col_selected or "SMILES"] * len(rows),
        }
    )
    if args.id_map_format in ("csv", "both"):
        ids_df.to_csv(Path(outdir) / "ids.csv", index=False)
    if args.id_map_format in ("parquet", "both"):
        write_table(ids_df, str(Path(outdir) / "ids.parquet"))

    # Save fragment vocabs + frequencies (freq values reflect the current dataset; vocab can be frozen via --spec).
    _json_dump(
        {"vocab": brics_vocab, "freq": brics_freq, "unknown_dropped": unknown_brics if spec is not None else 0},
        Path(outdir) / "fragments_brics.json",
    )
    _json_dump(
        {"vocab": recap_vocab, "freq": recap_freq, "unknown_dropped": unknown_recap if spec is not None else 0},
        Path(outdir) / "fragments_recap.json",
    )

    # Save interpretable summary
    if args.write_interpretable and interp_rows:
        interp_df = pd.DataFrame(interp_rows)
        if args.interpretable_format in ("csv", "both"):
            interp_df.to_csv(Path(outdir) / "features_interpretable.csv", index=False)
        if args.interpretable_format in ("parquet", "both"):
            write_table(interp_df, str(Path(outdir) / "features_interpretable.parquet"))

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

    # Write a reusable spec (fragment vocab + block layout) unless disabled.
    spec_written_path = None
    if args.write_spec:
        spec_out = {
            "spec_version": 1,
            "created_at": ts,
            "fragments": asdict(frag_spec),
            "fingerprints": asdict(cfg),
            "fg_names": fg_names,
            "numeric_cols": numeric_cols,
            "blocks": blocks,
            "fragment_vocab": {"BRICS": brics_vocab, "RECAP": recap_vocab},
        }
        spec_written_path = str(Path(outdir) / "features_spec.json")
        _json_dump(spec_out, spec_written_path)

    # Metadata
    meta = {
        "created_at": ts,
        "input": {
            "path": os.path.abspath(in_path),
            "sha256": input_sha,
            "size_bytes": os.path.getsize(in_path) if os.path.exists(in_path) else None,
            "is_table": bool(table_input),
            "format": table_format,
        },
        "structure_basis": {
            "policy": "auto (prefers Calc_Canonical_SMILES)",
            "smiles_col_requested": args.smiles_col,
            "smiles_col_selected": smiles_col_selected,
            "id_col_selected": id_col_selected,
        },
        "spec": {
            "used": spec_path,
            "written": spec_written_path,
        },
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
        "outputs": {
            "id_map_format": args.id_map_format,
            "interpretable_format": args.interpretable_format,
        },
        "versions": {
            "python": sys.version.split()[0],
            "rdkit": getattr(rdkit, "__version__", "unknown"),
            "numpy": np.__version__,
            "pandas": pd.__version__,
        },
    }

    _json_dump(meta, Path(outdir) / "features_metadata.json")

    print(str(Path(outdir).resolve()))


if __name__ == "__main__":
    main()

