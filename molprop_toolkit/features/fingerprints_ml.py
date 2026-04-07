"""Fingerprint helpers for ML exports.

This module provides:
- a config dataclass describing sizes and options,
- hashed *count* fingerprints (sparse dicts),
- bit-vector Morgan fingerprints (for DL exports),
- a consistent packing helper.

RDKit is imported lazily so importing molprop_toolkit does not require RDKit.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict

import numpy as np


@dataclass(frozen=True)
class MLFPConfig:
    morgan_radius: int = 2
    morgan_nbits: int = 4096
    morgan_use_chirality: bool = True
    morgan_use_features: bool = False

    atom_pair_nbits: int = 2048
    torsion_nbits: int = 2048

    # DL export bit size (not used for sparse counts)
    bit_nbits: int = 2048


def _fold_counts(counts: Dict[int, int], nbits: int) -> Dict[int, int]:
    out: Dict[int, int] = {}
    for k, v in counts.items():
        kk = int(k) % int(nbits)
        out[kk] = int(out.get(kk, 0) + int(v))
    return out


def compute_default_sparse_blocks(
    mol: Any, cfg: MLFPConfig
) -> Dict[str, Dict[int, int]]:
    """Return sparse hashed-count blocks for a molecule.

    Output is a dict with keys: Morgan, AtomPair, Torsion.
    Each value is a dict: {column_index: count} with indices already folded to [0, nbits).
    """

    try:
        from rdkit.Chem import AllChem, rdMolDescriptors  # type: ignore

        morgan_raw = AllChem.GetMorganFingerprint(
            mol,
            int(cfg.morgan_radius),
            useChirality=bool(cfg.morgan_use_chirality),
            useFeatures=bool(cfg.morgan_use_features),
        )
        morgan_counts = morgan_raw.GetNonzeroElements()  # type: ignore[attr-defined]

        ap_raw = rdMolDescriptors.GetHashedAtomPairFingerprint(
            mol, nBits=int(cfg.atom_pair_nbits)
        )
        ap_counts = ap_raw.GetNonzeroElements()  # type: ignore[attr-defined]

        tor_raw = rdMolDescriptors.GetHashedTopologicalTorsionFingerprint(
            mol, nBits=int(cfg.torsion_nbits)
        )
        tor_counts = tor_raw.GetNonzeroElements()  # type: ignore[attr-defined]

        return {
            "Morgan": _fold_counts(
                {int(k): int(v) for k, v in morgan_counts.items()},
                int(cfg.morgan_nbits),
            ),
            "AtomPair": _fold_counts(
                {int(k): int(v) for k, v in ap_counts.items()}, int(cfg.atom_pair_nbits)
            ),
            "Torsion": _fold_counts(
                {int(k): int(v) for k, v in tor_counts.items()}, int(cfg.torsion_nbits)
            ),
        }
    except Exception:
        return {"Morgan": {}, "AtomPair": {}, "Torsion": {}}


def morgan_bitvect(
    mol: Any, *, radius: int, nbits: int, use_chirality: bool = True
) -> Any:
    """Return an RDKit ExplicitBitVect (Morgan)."""

    from rdkit.Chem import AllChem  # type: ignore

    return AllChem.GetMorganFingerprintAsBitVect(
        mol, int(radius), nBits=int(nbits), useChirality=bool(use_chirality)
    )


def pack_bitvect(bitvect: Any) -> bytes:
    """Pack an RDKit bit-vector into bytes (np.packbits over a 0/1 array)."""

    from rdkit import DataStructs  # type: ignore

    nbits = int(bitvect.GetNumBits())
    arr = np.zeros((nbits,), dtype=np.uint8)
    DataStructs.ConvertToNumpyArray(bitvect, arr)
    packed = np.packbits(arr)
    return packed.tobytes()
