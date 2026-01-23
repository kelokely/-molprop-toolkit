#!/usr/bin/env python3
"""3D conformer generation + 3D descriptors (Task 7).

This module is intentionally self-contained so calculators can optionally compute
3D shape descriptors without cluttering the main v5 calculator.

Implementation notes:
- Uses RDKit ETKDG to embed multiple conformers.
- Optionally minimizes with MMFF94s (preferred) or UFF.
- Selects the lowest-energy conformer as the representative for 3D descriptors.

All outputs are returned as a dict of `3D_*` columns.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors3D, rdMolDescriptors


@dataclass(frozen=True)
class ThreeDConfig:
    enabled: bool = False
    num_confs: int = 10
    prune_rms: float = 0.5
    max_attempts: int = 20
    random_seed: int = 0
    minimize: str = "mmff"  # mmff|uff|none


def _pick_lowest_energy_conf(energies: Dict[int, float]) -> Optional[int]:
    if not energies:
        return None
    return min(energies.items(), key=lambda kv: kv[1])[0]


def compute_3d_descriptors(mol: Chem.Mol, cfg: ThreeDConfig) -> Dict[str, object]:
    """Compute 3D descriptor columns for a molecule.

    Returns a dict with metadata columns plus a set of RDKit 3D shape
    descriptors. If embedding fails, returns a dict with status fields and
    descriptor values set to None.
    """

    out: Dict[str, object] = {
        "3D_Mode": "off" if not cfg.enabled else "etkdg",
        "3D_NumConfs": None,
        "3D_SelectedConfId": None,
        "3D_Minimizer": cfg.minimize,
        "3D_Converged": None,
        "3D_BestEnergy": None,
        "3D_Asphericity": None,
        "3D_Eccentricity": None,
        "3D_InertialShapeFactor": None,
        "3D_NPR1": None,
        "3D_NPR2": None,
        "3D_PMI1": None,
        "3D_PMI2": None,
        "3D_PMI3": None,
        "3D_RadiusOfGyration": None,
        "3D_SpherocityIndex": None,
        "3D_PBF": None,
    }

    if not cfg.enabled:
        return out

    if mol is None:
        out["3D_Mode"] = "fail"
        return out

    m = Chem.AddHs(Chem.Mol(mol))

    params = AllChem.ETKDGv3()
    params.randomSeed = int(cfg.random_seed)
    params.pruneRmsThresh = float(cfg.prune_rms)
    params.maxAttempts = int(cfg.max_attempts)

    conf_ids = list(AllChem.EmbedMultipleConfs(m, numConfs=int(cfg.num_confs), params=params))
    out["3D_NumConfs"] = int(len(conf_ids))

    if not conf_ids:
        out["3D_Mode"] = "fail"
        out["3D_Converged"] = False
        return out

    energies_by_conf: Dict[int, float] = {}
    converged = True

    minimize = (cfg.minimize or "none").lower()
    if minimize == "mmff":
        props = AllChem.MMFFGetMoleculeProperties(m, mmffVariant="MMFF94s")
        if props is None:
            minimize = "uff"

    if minimize == "mmff":
        results = AllChem.MMFFOptimizeMoleculeConfs(m, numThreads=0, mmffVariant="MMFF94s")
        for conf_id, (status, energy) in zip(conf_ids, results):
            if status != 0:
                converged = False
            energies_by_conf[int(conf_id)] = float(energy)
    elif minimize == "uff":
        results = AllChem.UFFOptimizeMoleculeConfs(m, numThreads=0)
        for conf_id, (status, energy) in zip(conf_ids, results):
            if status != 0:
                converged = False
            energies_by_conf[int(conf_id)] = float(energy)
    else:
        converged = True

    out["3D_Converged"] = bool(converged)

    best_conf = _pick_lowest_energy_conf(energies_by_conf) if energies_by_conf else int(conf_ids[0])
    out["3D_SelectedConfId"] = int(best_conf)
    out["3D_BestEnergy"] = float(energies_by_conf.get(int(best_conf))) if energies_by_conf else None

    # Compute descriptors on heavy-atom molecule; keep conformer coordinates.
    m2 = Chem.RemoveHs(m)

    def _safe_float(v) -> Optional[float]:
        try:
            if v is None:
                return None
            return float(v)
        except Exception:
            return None

    def _try(fn) -> Optional[float]:
        try:
            return _safe_float(fn())
        except Exception:
            return None

    out["3D_Asphericity"] = _try(lambda: Descriptors3D.Asphericity(m2, confId=int(best_conf)))
    out["3D_Eccentricity"] = _try(lambda: Descriptors3D.Eccentricity(m2, confId=int(best_conf)))
    out["3D_InertialShapeFactor"] = _try(lambda: Descriptors3D.InertialShapeFactor(m2, confId=int(best_conf)))
    out["3D_NPR1"] = _try(lambda: Descriptors3D.NPR1(m2, confId=int(best_conf)))
    out["3D_NPR2"] = _try(lambda: Descriptors3D.NPR2(m2, confId=int(best_conf)))
    out["3D_PMI1"] = _try(lambda: Descriptors3D.PMI1(m2, confId=int(best_conf)))
    out["3D_PMI2"] = _try(lambda: Descriptors3D.PMI2(m2, confId=int(best_conf)))
    out["3D_PMI3"] = _try(lambda: Descriptors3D.PMI3(m2, confId=int(best_conf)))
    out["3D_RadiusOfGyration"] = _try(lambda: Descriptors3D.RadiusOfGyration(m2, confId=int(best_conf)))
    out["3D_SpherocityIndex"] = _try(lambda: Descriptors3D.SpherocityIndex(m2, confId=int(best_conf)))

    # PBF lives in rdMolDescriptors in older RDKit builds.
    out["3D_PBF"] = _try(lambda: rdMolDescriptors.CalcPBF(m2, confId=int(best_conf)))

    if all(out.get(k) is None for k in [
        "3D_Asphericity",
        "3D_Eccentricity",
        "3D_InertialShapeFactor",
        "3D_NPR1",
        "3D_NPR2",
        "3D_PMI1",
        "3D_PMI2",
        "3D_PMI3",
        "3D_RadiusOfGyration",
        "3D_SpherocityIndex",
        "3D_PBF",
    ]):
        out["3D_Mode"] = "fail"

    return out

