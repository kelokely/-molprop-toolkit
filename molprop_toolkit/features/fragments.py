"""Fragment inventory utilities (BRICS / RECAP).

The featurizer uses fragments primarily as *presence* features and for interpretable summaries.

Implementation notes:
- We use RDKit's BRICS and RECAP decompositions.
- RDKit imports are lazy so importing molprop_toolkit does not require RDKit.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List


@dataclass(frozen=True)
class FragmentSpec:
    """Fragment extraction configuration."""

    method: str = "both"  # one of: none|brics|recap|both
    max_per_mol: int = 128


def fragment_inventory(mol: Any, spec: FragmentSpec) -> Dict[str, List[str]]:
    """Return a fragment inventory for a molecule.

    Output keys are always present: "BRICS" and "RECAP".
    Values are lists of SMILES strings, truncated to spec.max_per_mol.
    """

    method = (spec.method or "both").lower()
    if method not in ("none", "brics", "recap", "both"):
        raise ValueError(f"Unknown fragment method: {spec.method}")

    out: Dict[str, List[str]] = {"BRICS": [], "RECAP": []}

    if method in ("brics", "both"):
        try:
            from rdkit.Chem import BRICS  # type: ignore

            frags = sorted(BRICS.BRICSDecompose(mol))
            out["BRICS"] = frags[: int(spec.max_per_mol)]
        except Exception:
            out["BRICS"] = []

    if method in ("recap", "both"):
        try:
            from rdkit.Chem import Recap  # type: ignore

            # RecapDecompose returns a tree; children keys are SMILES.
            tree = Recap.RecapDecompose(mol)
            if tree is not None:
                out["RECAP"] = sorted(list(tree.GetAllChildren().keys()))[: int(spec.max_per_mol)]
        except Exception:
            out["RECAP"] = []

    return out

