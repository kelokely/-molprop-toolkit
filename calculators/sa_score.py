#!/usr/bin/env python3
"""Synthetic accessibility + complexity metrics (Task 6).

This module provides:
- SA_Score: the Ertl & Schuffenhauer synthetic accessibility score (1–10)
- simple complexity descriptors (BertzCT, BalabanJ, etc.)

The SA score implementation is adapted from RDKit's Contrib/SA_Score
(sascorer.py) and uses the accompanying fragment score file.

Reference:
  Peter Ertl and Ansgar Schuffenhauer,
  "Estimation of Synthetic Accessibility Score of Drug-like Molecules based on
  Molecular Complexity and Fragment Contributions", Journal of Cheminformatics
  1:8 (2009). http://www.jcheminf.com/content/1/1/8

License:
  The upstream RDKit Contrib implementation is BSD-3-Clause (Novartis).
"""

from __future__ import annotations

import gzip
import math
import os.path as op
import pickle
from functools import lru_cache
from typing import Dict, Optional

from rdkit import Chem
from rdkit.Chem import Descriptors, rdFingerprintGenerator, rdMolDescriptors


# Morgan fingerprint generator for the SA score fragment contributions.
_mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2)


@lru_cache(maxsize=1)
def _read_fragment_scores() -> Dict[int, float]:
    """Load fragment contributions for SA score.

    The fragment score file is shipped in this repo at:
      calculators/data/fpscores.pkl.gz
    """

    path = op.join(op.dirname(__file__), "data", "fpscores.pkl.gz")
    try:
        data = pickle.load(gzip.open(path))
    except FileNotFoundError as e:  # pragma: no cover
        raise RuntimeError(
            f"Missing SA fragment score file: {path}. "
            "Ensure calculators/data/fpscores.pkl.gz is present in the repo."
        ) from e

    out: Dict[int, float] = {}
    for row in data:
        score = float(row[0])
        for j in range(1, len(row)):
            out[int(row[j])] = score
    return out


def calc_sa_score(mol: Chem.Mol) -> Optional[float]:
    """Calculate the SA score (1–10, lower is easier).

    Returns None if mol is empty.
    """

    if mol is None or mol.GetNumAtoms() == 0:
        return None

    fscores = _read_fragment_scores()

    # 1) fragment score (average fragment contributions weighted by counts)
    sfp = _mfpgen.GetSparseCountFingerprint(mol)
    nze = sfp.GetNonzeroElements()

    score1 = 0.0
    nf = 0
    for fid, count in nze.items():
        nf += count
        score1 += fscores.get(int(fid), -4.0) * count

    score1 /= float(nf) if nf else 1.0

    # 2) feature penalties
    n_atoms = mol.GetNumAtoms()
    n_chiral = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))

    ri = mol.GetRingInfo()
    n_spiro = rdMolDescriptors.CalcNumSpiroAtoms(mol)
    n_bridge = rdMolDescriptors.CalcNumBridgeheadAtoms(mol)

    n_macro = 0
    for ring in ri.AtomRings():
        if len(ring) > 8:
            n_macro += 1

    size_penalty = n_atoms**1.005 - n_atoms
    stereo_penalty = math.log10(n_chiral + 1)
    spiro_penalty = math.log10(n_spiro + 1)
    bridge_penalty = math.log10(n_bridge + 1)

    # differs slightly from the paper; matches RDKit contrib
    macro_penalty = math.log10(2) if n_macro > 0 else 0.0

    score2 = 0.0 - size_penalty - stereo_penalty - spiro_penalty - bridge_penalty - macro_penalty

    # 3) fingerprint density correction (symmetry)
    score3 = 0.0
    num_bits = len(nze)
    if n_atoms > num_bits and num_bits > 0:
        score3 = math.log(float(n_atoms) / float(num_bits)) * 0.5

    sascore = score1 + score2 + score3

    # transform raw value into [1, 10]
    min_v = -4.0
    max_v = 2.5
    sascore = 11.0 - (sascore - min_v + 1.0) / (max_v - min_v) * 9.0

    # smooth the 10-end
    if sascore > 8.0:
        sascore = 8.0 + math.log(sascore + 1.0 - 9.0)

    if sascore > 10.0:
        sascore = 10.0
    elif sascore < 1.0:
        sascore = 1.0

    return float(sascore)


def sa_class(sa_score: Optional[float]) -> Optional[str]:
    if sa_score is None:
        return None
    if sa_score <= 3.0:
        return "Easy"
    if sa_score <= 6.0:
        return "Moderate"
    return "Hard"


def calc_sa_and_complexity(mol: Chem.Mol) -> Dict[str, object]:
    """Compute SA_Score plus a small, interpretable complexity bundle."""

    out: Dict[str, object] = {}

    sa = calc_sa_score(mol)
    out["SA_Score"] = round(sa, 3) if sa is not None else None
    out["SA_Class"] = sa_class(sa)

    # Complexity metrics: keep small and interpretable.
    try:
        out["Complexity_BertzCT"] = round(float(Descriptors.BertzCT(mol)), 3)
    except Exception:
        out["Complexity_BertzCT"] = None

    try:
        out["Complexity_BalabanJ"] = round(float(Descriptors.BalabanJ(mol)), 3)
    except Exception:
        out["Complexity_BalabanJ"] = None

    try:
        out["Complexity_HallKierAlpha"] = round(float(Descriptors.HallKierAlpha(mol)), 3)
    except Exception:
        out["Complexity_HallKierAlpha"] = None

    try:
        out["Complexity_NumSpiroAtoms"] = int(rdMolDescriptors.CalcNumSpiroAtoms(mol))
    except Exception:
        out["Complexity_NumSpiroAtoms"] = None

    try:
        out["Complexity_NumBridgeheadAtoms"] = int(rdMolDescriptors.CalcNumBridgeheadAtoms(mol))
    except Exception:
        out["Complexity_NumBridgeheadAtoms"] = None

    try:
        out["Complexity_NumStereoCenters"] = int(rdMolDescriptors.CalcNumAtomStereoCenters(mol))
        out["Complexity_NumUnspecifiedStereoCenters"] = int(rdMolDescriptors.CalcNumUnspecifiedAtomStereoCenters(mol))
    except Exception:
        out["Complexity_NumStereoCenters"] = None
        out["Complexity_NumUnspecifiedStereoCenters"] = None

    try:
        out["Complexity_NumHeteroAtoms"] = int(rdMolDescriptors.CalcNumHeteroatoms(mol))
    except Exception:
        out["Complexity_NumHeteroAtoms"] = None

    # Aromatic proportion (atoms) is often a useful proxy for "flatness".
    try:
        n_heavy = int(mol.GetNumHeavyAtoms()) or 1
        n_arom = sum(1 for a in mol.GetAtoms() if a.GetIsAromatic())
        out["Complexity_AromaticAtomFrac"] = round(float(n_arom) / float(n_heavy), 3)
    except Exception:
        out["Complexity_AromaticAtomFrac"] = None

    return out

