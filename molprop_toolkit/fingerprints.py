"""Compatibility helpers for legacy fingerprint imports.

The similarity refactor moved the main implementation under
``molprop_toolkit.similarity``. Several top-level scripts still import the old
``molprop_toolkit.fingerprints`` module, so this shim keeps that public surface
working while the rest of the codebase finishes consolidating around the new
package layout.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal, TypeAlias

import numpy as np
import pandas as pd

from molprop_toolkit.core import detect_best_smiles_column, detect_id_column, read_table
from molprop_toolkit.similarity.fingerprints import get_fingerprint
from molprop_toolkit.similarity.metrics import bulk_similarity
from molprop_toolkit.similarity.search import _load_smiles_file

FingerprintKind: TypeAlias = Literal[
    "morgan",
    "rdkit",
    "maccs",
    "atom_pair",
    "torsion",
    "usr",
    "usrcat",
]


@dataclass(frozen=True)
class FingerprintSpec:
    kind: FingerprintKind = "morgan"
    radius: int = 2
    nbits: int = 2048
    use_chirality: bool = False
    use_features: bool = False
    num_confs: int = 10
    minimize: str = "none"


@dataclass(frozen=True)
class ReferenceSmiles:
    label: str
    smiles: str


def _normalize_kind(kind: str) -> str:
    normalized = str(kind or "morgan").strip().lower()
    if normalized == "atompair":
        return "atom_pair"
    return normalized


def _similarity_column_suffix(kind: str) -> str:
    normalized = _normalize_kind(kind)
    if normalized == "atom_pair":
        return "atom_pair"
    return normalized


def _backend_kind(kind: str) -> str:
    normalized = _normalize_kind(kind)
    if normalized == "atom_pair":
        return "atompair"
    return normalized


def _bitvect_kwargs(spec: FingerprintSpec) -> dict[str, object]:
    kind = _backend_kind(spec.kind)
    if kind == "morgan":
        return {
            "radius": int(spec.radius),
            "nBits": int(spec.nbits),
            "useChirality": bool(spec.use_chirality),
            "useFeatures": bool(spec.use_features),
        }
    if kind == "rdkit":
        return {"fpSize": int(spec.nbits)}
    if kind in {"atompair", "torsion"}:
        return {"nBits": int(spec.nbits)}
    return {}


def _best_conf_id(mol, spec: FingerprintSpec) -> tuple[object | None, int | None]:
    from rdkit import Chem
    from rdkit.Chem import AllChem

    prepared = Chem.AddHs(Chem.Mol(mol))

    params = AllChem.ETKDGv3()
    params.randomSeed = 0
    conf_ids = list(AllChem.EmbedMultipleConfs(prepared, numConfs=int(spec.num_confs), params=params))
    if not conf_ids:
        return None, None

    minimize = str(spec.minimize or "none").lower()
    if minimize == "mmff":
        props = AllChem.MMFFGetMoleculeProperties(prepared, mmffVariant="MMFF94s")
        if props is None:
            minimize = "uff"

    if minimize == "mmff":
        results = AllChem.MMFFOptimizeMoleculeConfs(prepared, numThreads=0, mmffVariant="MMFF94s")
        energies = {int(conf_id): float(energy) for conf_id, (_, energy) in zip(conf_ids, results)}
        best_conf = min(energies, key=energies.get)
    elif minimize == "uff":
        results = AllChem.UFFOptimizeMoleculeConfs(prepared, numThreads=0)
        energies = {int(conf_id): float(energy) for conf_id, (_, energy) in zip(conf_ids, results)}
        best_conf = min(energies, key=energies.get)
    else:
        best_conf = int(conf_ids[0])

    return Chem.RemoveHs(prepared), best_conf


def fingerprint_from_mol(mol, spec: FingerprintSpec):
    """Compute a fingerprint from an RDKit molecule using the legacy spec API."""

    kind = _normalize_kind(spec.kind)
    if kind in {"usr", "usrcat"}:
        from rdkit.Chem import rdMolDescriptors

        prepared, conf_id = _best_conf_id(mol, spec)
        if prepared is None or conf_id is None:
            return None
        if kind == "usr":
            return np.asarray(rdMolDescriptors.GetUSR(prepared, confId=int(conf_id)), dtype=float)
        return np.asarray(rdMolDescriptors.GetUSRCAT(prepared, confId=int(conf_id)), dtype=float)

    return get_fingerprint(
        mol,
        fp_type=_backend_kind(kind),
        **_bitvect_kwargs(spec),
    )


def load_reference_smiles(
    *,
    df: pd.DataFrame,
    id_col: str,
    smiles_col: str,
    ref_ids: list[str] | None = None,
    ref_smiles: list[str] | None = None,
    ref_file: str | None = None,
) -> list[ReferenceSmiles]:
    """Resolve reference compounds from IDs, SMILES strings, or a file."""

    refs: list[ReferenceSmiles] = []
    ref_ids = list(ref_ids or [])
    ref_smiles = list(ref_smiles or [])

    if ref_ids:
        missing = [rid for rid in ref_ids if rid not in set(df[id_col].astype(str))]
        if missing:
            raise SystemExit(f"Reference IDs not found in {id_col}: {missing}")
        for rid in ref_ids:
            row = df[df[id_col].astype(str) == str(rid)].iloc[0]
            refs.append(ReferenceSmiles(label=str(rid), smiles=str(row[smiles_col])))

    for idx, smiles in enumerate(ref_smiles, start=1):
        refs.append(ReferenceSmiles(label=f"SMILES_{idx}", smiles=str(smiles)))

    if ref_file:
        ref_path = Path(ref_file)
        if ref_path.suffix.lower() in {".csv", ".tsv", ".parquet"}:
            ref_df = read_table(str(ref_path))
        else:
            ref_df = _load_smiles_file(str(ref_path))
        ref_smiles_col = detect_best_smiles_column(ref_df)
        if not ref_smiles_col:
            raise SystemExit(f"Could not detect a SMILES column in reference file: {ref_file}")
        ref_id_col = detect_id_column(ref_df)
        for _, row in ref_df.iterrows():
            refs.append(
                ReferenceSmiles(
                    label=str(row.get(ref_id_col) or f"Ref_{len(refs) + 1}"),
                    smiles=str(row[ref_smiles_col]),
                )
            )

    if not refs:
        raise SystemExit("No reference compounds were provided.")

    return refs


def _similarity_to_references(query_fp, reference_fps: list[object], *, use_euclidean: bool) -> np.ndarray:
    if use_euclidean:
        if query_fp is None:
            return np.full(len(reference_fps), np.nan, dtype=float)
        values = []
        for ref_fp in reference_fps:
            if ref_fp is None:
                values.append(np.nan)
                continue
            distance = float(np.linalg.norm(np.asarray(query_fp, dtype=float) - np.asarray(ref_fp, dtype=float)))
            values.append(1.0 / (1.0 + distance))
        return np.asarray(values, dtype=float)

    return bulk_similarity(query_fp, reference_fps, metric="tanimoto")


def compute_similarity_columns(
    df: pd.DataFrame,
    *,
    refs: list[ReferenceSmiles],
    spec: FingerprintSpec,
    smiles_col: str,
    id_col: str,
    column_prefix: str = "Sim",
) -> pd.DataFrame:
    """Append similarity-to-reference columns used by reports and picklists."""

    from rdkit import Chem

    out = df.copy()
    column_name = f"{column_prefix}_{_similarity_column_suffix(spec.kind)}"
    best_ref_name = f"{column_name}_BestRef"

    reference_fps: list[object] = []
    reference_labels: list[str] = []
    for ref in refs:
        mol = Chem.MolFromSmiles(str(ref.smiles))
        fp = fingerprint_from_mol(mol, spec) if mol is not None else None
        if fp is None:
            continue
        reference_labels.append(str(ref.label))
        reference_fps.append(fp)

    if not reference_fps:
        raise SystemExit("Could not compute fingerprints for any reference compounds.")

    use_euclidean = _normalize_kind(spec.kind) in {"usr", "usrcat"}
    similarities: list[float] = []
    best_refs: list[str] = []

    for smiles in out[smiles_col].astype(str).tolist():
        mol = Chem.MolFromSmiles(smiles)
        fp = fingerprint_from_mol(mol, spec) if mol is not None else None
        sims = _similarity_to_references(fp, reference_fps, use_euclidean=use_euclidean)
        if np.isnan(sims).all():
            similarities.append(np.nan)
            best_refs.append("")
            continue
        best_idx = int(np.nanargmax(sims))
        similarities.append(float(sims[best_idx]))
        best_refs.append(reference_labels[best_idx])

    out[column_name] = similarities
    out[best_ref_name] = best_refs
    return out


__all__ = [
    "FingerprintKind",
    "FingerprintSpec",
    "ReferenceSmiles",
    "compute_similarity_columns",
    "fingerprint_from_mol",
    "load_reference_smiles",
]
