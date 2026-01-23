#!/usr/bin/env python3
"""Tautomer enumeration and selection utilities.

This is Task 3: expanded tautomer handling.

We support three practical modes:

1) prep-canonical
   The preparation pipeline (prep.py) canonicalizes the tautomer using RDKit's
   TautomerEnumerator. This is fast and deterministic, but you only see the
   resulting tautomer, not the size of the tautomer set.

2) none
   Tautomer canonicalization is disabled.

3) enumerate
   Enumerate all tautomers (up to a max) from the prepared parent structure and
   select a representative tautomer (RDKit canonical tautomer by default).
   This mode exposes useful metadata: tautomer count and top-k tautomer SMILES.

Important notes:
- RDKit tautomer enumeration is rule-based and may not match experimental
  tautomer distributions.
- The ordering of enumerated tautomers can vary across RDKit versions.
  We attempt to score and sort using ScoreTautomer when available.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Tuple

from rdkit import Chem

try:
    from rdkit.Chem.MolStandardize import rdMolStandardize
except Exception:  # pragma: no cover
    rdMolStandardize = None  # type: ignore


@dataclass
class TautomerConfig:
    mode: str = "prep-canonical"  # prep-canonical | none | enumerate
    max_variants: int = 64
    top_k: int = 5


@dataclass
class TautomerReport:
    mode: str
    count: int
    selected_smiles: str
    topk_smiles: str
    changed: bool
    notes: str = ""


def _mol_to_can_smiles(m: Chem.Mol) -> str:
    return Chem.MolToSmiles(m, canonical=True)


def enumerate_tautomers(mol: Chem.Mol, max_variants: int = 64) -> List[Chem.Mol]:
    if rdMolStandardize is None:
        return [mol]

    te = rdMolStandardize.TautomerEnumerator()

    # RDKit has an optional MaxTransforms/MaxTautomers depending on version.
    try:
        params = te.GetParameters()
        if hasattr(params, "maxTautomers"):
            params.maxTautomers = int(max_variants)
            te.SetParameters(params)
    except Exception:
        pass

    try:
        tautomers = list(te.Enumerate(mol))
        if not tautomers:
            return [mol]
        # Hard cap if RDKit ignores max.
        return tautomers[:max_variants]
    except Exception:
        return [mol]


def _score_tautomers(te, mols: List[Chem.Mol]) -> Optional[List[float]]:
    try:
        score_fn = getattr(te, "ScoreTautomer", None)
        if score_fn is None:
            return None
        return [float(score_fn(m)) for m in mols]
    except Exception:
        return None


def select_representative_tautomer(
    mol: Chem.Mol,
    max_variants: int = 64,
    top_k: int = 5,
) -> Tuple[Chem.Mol, TautomerReport]:
    """Enumerate tautomers and select a representative tautomer.

    Selection strategy: RDKit canonical tautomer (TautomerEnumerator.Canonicalize).
    We additionally report the top-k tautomer SMILES (sorted by ScoreTautomer when
    available; otherwise enumerator order).
    """

    in_smiles = _mol_to_can_smiles(mol)

    if rdMolStandardize is None:
        rep = TautomerReport(
            mode="enumerate",
            count=1,
            selected_smiles=in_smiles,
            topk_smiles=in_smiles,
            changed=False,
            notes="MolStandardize unavailable; tautomer enumeration disabled",
        )
        return mol, rep

    te = rdMolStandardize.TautomerEnumerator()

    # Try to apply max.
    try:
        params = te.GetParameters()
        if hasattr(params, "maxTautomers"):
            params.maxTautomers = int(max_variants)
            te.SetParameters(params)
    except Exception:
        pass

    tautomers = enumerate_tautomers(mol, max_variants=max_variants)
    selected = mol
    try:
        selected = te.Canonicalize(mol)
    except Exception:
        selected = mol

    # Prepare top-k reporting
    smiles_list = [_mol_to_can_smiles(m) for m in tautomers]

    scores = _score_tautomers(te, tautomers)
    if scores is not None:
        ranked = sorted(zip(smiles_list, scores), key=lambda x: x[1], reverse=True)
        smiles_list = [s for s, _ in ranked]

    topk = ";".join(smiles_list[: max(1, int(top_k))])

    sel_smiles = _mol_to_can_smiles(selected)
    changed = sel_smiles != in_smiles

    rep = TautomerReport(
        mode="enumerate",
        count=len(tautomers),
        selected_smiles=sel_smiles,
        topk_smiles=topk,
        changed=changed,
        notes="",
    )

    return selected, rep

