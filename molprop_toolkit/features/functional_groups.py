"""Functional group and environment features.

These features are intended to be:
- deterministic,
- stable in their set of keys across runs,
- interpretable (counts/flags that can be read in SAR review).

RDKit is imported lazily so the overall toolkit can be installed without RDKit.
"""

from __future__ import annotations

from typing import Any, Dict

# A small, pragmatic functional-group panel. Keys are stable and always returned.
SMARTS = {
    "FG_Amide": "C(=O)N",
    "FG_Urea": "NC(=O)N",
    "FG_Carbamate": "OC(=O)N",
    "FG_Ester": "C(=O)O",
    "FG_CarboxylicAcid": "C(=O)[O;H,-]",
    "FG_Sulfonamide": "S(=O)(=O)N",
    "FG_Sulfone": "S(=O)(=O)",
    "FG_Alcohol": "[OX2H]",
    "FG_Phenol": "c[OX2H]",
    "FG_Amine": "[NX3;H2,H1;!$(NC=O)]",
    "FG_Aniline": "c[NX3;H2,H1]",
    "FG_Nitrile": "C#N",
    "FG_Nitro": "[N+](=O)[O-]",
    "FG_Halogen": "[F,Cl,Br,I]",
    "FG_Thiol": "[SX2H]",
    "FG_Thioether": "[SX2][#6]",
    "FG_Ether": "[OD2]([#6])[#6]",
    "FG_AromaticRing": "a1aaaaa1",
}


def compute_functional_group_features(mol: Any) -> Dict[str, float]:
    """Compute a stable set of functional-group and environment features.

    Returns a dict with all keys in SMARTS plus a small set of environment descriptors.
    Values are numeric (counts or continuous). Missing values are represented as 0.0.
    """

    # Base output with zeros (ensures stable key set).
    out: Dict[str, float] = dict.fromkeys(SMARTS.keys(), 0.0)

    try:
        from rdkit import Chem  # type: ignore
        from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors  # type: ignore

        # Compile patterns once per call (simple; could be cached later).
        patterns = {k: Chem.MolFromSmarts(v) for k, v in SMARTS.items()}

        for key, patt in patterns.items():
            if patt is None:
                continue
            try:
                out[key] = float(len(mol.GetSubstructMatches(patt)))
            except Exception:
                out[key] = 0.0

        # Environment features (also stable keys)
        out["Env_HeavyAtomCount"] = float(mol.GetNumHeavyAtoms())
        out["Env_RotatableBonds"] = float(Lipinski.NumRotatableBonds(mol))
        out["Env_Rings"] = float(rdMolDescriptors.CalcNumRings(mol))
        out["Env_AromaticRings"] = float(rdMolDescriptors.CalcNumAromaticRings(mol))
        out["Env_FractionCSP3"] = float(Descriptors.FractionCSP3(mol))

        # HBD/HBA are highly used in medchem; include them explicitly.
        out["Env_HBD"] = float(Lipinski.NumHDonors(mol))
        out["Env_HBA"] = float(Lipinski.NumHAcceptors(mol))

    except Exception:
        # RDKit missing or failed; return the stable zero-filled dict.
        out["Env_HeavyAtomCount"] = 0.0
        out["Env_RotatableBonds"] = 0.0
        out["Env_Rings"] = 0.0
        out["Env_AromaticRings"] = 0.0
        out["Env_FractionCSP3"] = 0.0
        out["Env_HBD"] = 0.0
        out["Env_HBA"] = 0.0

    return out
