#!/usr/bin/env python3
"""Stereochemistry handling utilities (Task 4).

The toolkit already reports stereochemistry completeness during preparation
(`Prep_StereoStatus`, `Prep_UnassignedChiral`, etc.). This Task 4 module adds an
explicit, controllable *stereo handling step* in the calculator pipeline.

Supported modes:
- keep: Preserve whatever stereochemistry is present (default).
- strip: Remove stereochemistry (treat as unspecified / racemic).
- enumerate: Enumerate stereoisomers for *unassigned* stereocenters/bonds and
  select a representative stereoisomer.

Notes:
- Many 2D descriptors used in this toolkit are stereo-invariant; nevertheless,
  explicitly tracking stereo choices improves traceability and supports future
  3D/structure-dependent extensions.
- Enumeration can explode combinatorially; always keep max_isomers bounded.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Tuple

from rdkit import Chem

try:
    from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
except Exception:  # pragma: no cover
    EnumerateStereoisomers = None  # type: ignore
    StereoEnumerationOptions = None  # type: ignore


@dataclass
class StereoConfig:
    mode: str = "keep"  # keep | strip | enumerate
    max_isomers: int = 32
    top_k: int = 5
    select: str = "canonical"  # canonical | first


@dataclass
class StereoReport:
    mode: str
    count: int
    selected_smiles: str
    topk_smiles: str
    changed: bool
    notes: str = ""


def _to_isomeric_smiles(m: Chem.Mol) -> str:
    return Chem.MolToSmiles(m, canonical=True, isomericSmiles=True)


def strip_stereo(mol: Chem.Mol) -> Chem.Mol:
    m = Chem.Mol(mol)
    try:
        Chem.RemoveStereochemistry(m)
    except Exception:
        pass
    return m


def enumerate_stereoisomers(mol: Chem.Mol, max_isomers: int) -> List[Chem.Mol]:
    if EnumerateStereoisomers is None or StereoEnumerationOptions is None:
        return [mol]

    try:
        opts = StereoEnumerationOptions(
            tryEmbedding=False,
            unique=True,
            onlyUnassigned=True,
            maxIsomers=int(max_isomers),
        )
        isos = list(EnumerateStereoisomers(mol, options=opts))
        if not isos:
            return [mol]
        return isos
    except Exception:
        return [mol]


def select_stereoisomer(mols: List[Chem.Mol], strategy: str = "canonical") -> Chem.Mol:
    if not mols:
        raise ValueError("No stereoisomers provided")

    strategy = (strategy or "first").lower()
    if strategy == "first":
        return mols[0]

    # canonical: choose the lexicographically smallest canonical isomeric SMILES
    smiles = [(_to_isomeric_smiles(m), m) for m in mols]
    smiles.sort(key=lambda x: x[0])
    return smiles[0][1]


def apply_stereo(mol: Chem.Mol, cfg: Optional[StereoConfig]) -> Tuple[Chem.Mol, StereoReport]:
    """Apply stereo mode to a molecule and return (mol_out, report)."""

    if cfg is None:
        cfg = StereoConfig()

    in_smiles = _to_isomeric_smiles(mol)
    mode = (cfg.mode or "keep").lower()

    if mode == "strip":
        m = strip_stereo(mol)
        out_smiles = _to_isomeric_smiles(m)
        rep = StereoReport(
            mode="strip",
            count=1,
            selected_smiles=out_smiles,
            topk_smiles=out_smiles,
            changed=(out_smiles != in_smiles),
            notes="",
        )
        return m, rep

    if mode == "enumerate":
        isos = enumerate_stereoisomers(mol, max_isomers=cfg.max_isomers)
        sel = select_stereoisomer(isos, strategy=cfg.select)
        iso_smiles = [_to_isomeric_smiles(m) for m in isos]
        topk = ";".join(iso_smiles[: max(1, int(cfg.top_k))])
        sel_smi = _to_isomeric_smiles(sel)
        rep = StereoReport(
            mode="enumerate",
            count=len(isos),
            selected_smiles=sel_smi,
            topk_smiles=topk,
            changed=(sel_smi != in_smiles),
            notes="",
        )
        return Chem.Mol(sel), rep

    # keep
    out_smiles = in_smiles
    rep = StereoReport(
        mode="keep",
        count=1,
        selected_smiles=out_smiles,
        topk_smiles=out_smiles,
        changed=False,
        notes="",
    )
    return mol, rep

