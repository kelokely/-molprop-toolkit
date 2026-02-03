"""Graph tensor export helpers.

This is a lightweight, framework-agnostic graph export intended for batch GNN training.

The output is a concatenation of per-molecule graphs:
- atom_features: (N_atoms, F)
- edge_index: (2, N_edges) with directed edges
- edge_features: (N_edges, E)
- mol_atom_ptr: (n_mols+1,) offsets into atom_features
- mol_edge_ptr: (n_mols+1,) offsets into edge arrays

RDKit is imported lazily.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List

import numpy as np


@dataclass(frozen=True)
class GraphSpec:
    include_hs: bool = False


def _atom_features(atom: Any) -> List[float]:
    # Compact numeric feature vector (stable length).
    z = float(atom.GetAtomicNum())
    deg = float(atom.GetTotalDegree())
    charge = float(atom.GetFormalCharge())
    aromatic = 1.0 if atom.GetIsAromatic() else 0.0
    hcount = float(atom.GetTotalNumHs())
    mass = float(atom.GetMass())
    return [z, deg, charge, aromatic, hcount, mass]


def _bond_features(bond: Any) -> List[float]:
    bt = str(bond.GetBondType())
    # One-hot-ish bond type
    single = 1.0 if bt == "SINGLE" else 0.0
    double = 1.0 if bt == "DOUBLE" else 0.0
    triple = 1.0 if bt == "TRIPLE" else 0.0
    aromatic = 1.0 if bt == "AROMATIC" else 0.0
    conj = 1.0 if bond.GetIsConjugated() else 0.0
    ring = 1.0 if bond.IsInRing() else 0.0
    return [single, double, triple, aromatic, conj, ring]


def graphs_from_smiles_list(
    smiles: List[str], spec: GraphSpec
) -> Dict[str, np.ndarray]:
    try:
        from rdkit import Chem  # type: ignore
    except Exception as e:  # pragma: no cover
        raise RuntimeError("RDKit is required for graph exports") from e

    atom_feat: List[List[float]] = []
    edge_index_src: List[int] = []
    edge_index_dst: List[int] = []
    edge_feat: List[List[float]] = []

    mol_atom_ptr = [0]
    mol_edge_ptr = [0]

    atom_offset = 0
    edge_offset = 0

    for smi in smiles:
        if not smi:
            mol_atom_ptr.append(atom_offset)
            mol_edge_ptr.append(edge_offset)
            continue

        m = Chem.MolFromSmiles(smi)
        if m is None:
            mol_atom_ptr.append(atom_offset)
            mol_edge_ptr.append(edge_offset)
            continue

        if spec.include_hs:
            m = Chem.AddHs(m)

        # Atoms
        for a in m.GetAtoms():
            atom_feat.append(_atom_features(a))
        n_atoms = m.GetNumAtoms()

        # Bonds (directed)
        for b in m.GetBonds():
            i = int(b.GetBeginAtomIdx()) + atom_offset
            j = int(b.GetEndAtomIdx()) + atom_offset
            bf = _bond_features(b)
            edge_index_src.extend([i, j])
            edge_index_dst.extend([j, i])
            edge_feat.extend([bf, bf])

        n_edges = len(m.GetBonds()) * 2

        atom_offset += n_atoms
        edge_offset += n_edges
        mol_atom_ptr.append(atom_offset)
        mol_edge_ptr.append(edge_offset)

    atom_features = np.asarray(atom_feat, dtype=np.float32)
    edge_index = np.asarray([edge_index_src, edge_index_dst], dtype=np.int64)
    edge_features = np.asarray(edge_feat, dtype=np.float32)

    return {
        "atom_features": atom_features,
        "edge_index": edge_index,
        "edge_features": edge_features,
        "mol_atom_ptr": np.asarray(mol_atom_ptr, dtype=np.int64),
        "mol_edge_ptr": np.asarray(mol_edge_ptr, dtype=np.int64),
    }
