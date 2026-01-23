#!/usr/bin/env python3
"""══════════════════════════════════════════════════════════════════════════════
MolProp Toolkit Calculator (v5, extended)
═══════════════════════════════════════════════════════════════════════════════

This calculator extends v4 with additional, still-heuristic property families
that are often useful for early triage: simple solubility estimates, simple
permeability estimates, and coarse PK heuristics. The goal is compatibility
with the analyzers in this repository rather than perfect physical accuracy.

Usage examples:

  python calculators/mpo_v5.py input.smi -o results.csv --fill-missing 0
  python calculators/mpo_v5.py input.smi -o results.csv --pIC50 7.0

Disclaimer: The solubility/permeability/PK fields are approximate screening
models intended for ranking and qualitative flags.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, QED, rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

# Structure preparation (standardization, tautomer canonicalization, stereo reporting)
# Task 15A: allow running both as a script (local imports) and as an installed package (relative imports).
try:  # package import
    from .prep import PrepConfig, prepare_mol
    from .ionization import ionization_features
    from .tautomers import TautomerConfig, select_representative_tautomer
    from .stereo import StereoConfig, apply_stereo
    from .protomers import ProtomerConfig, protomerize_smiles
    from .sa_score import calc_sa_and_complexity
    from .geom3d import ThreeDConfig, compute_3d_descriptors
    from .developability import developability_features
except Exception:  # script import
    from prep import PrepConfig, prepare_mol
    from ionization import ionization_features
    from tautomers import TautomerConfig, select_representative_tautomer
    from stereo import StereoConfig, apply_stereo
    from protomers import ProtomerConfig, protomerize_smiles
    from sa_score import calc_sa_and_complexity
    from geom3d import ThreeDConfig, compute_3d_descriptors
    from developability import developability_features


# Reuse the same pattern sets as v4 for consistency.
TOX_ALERTS: Dict[str, str] = {
    "acyl_halide": "[CX3](=[OX1])[F,Cl,Br,I]",
    "sulfonyl_halide": "[SX4](=[OX1])(=[OX1])[F,Cl,Br,I]",
    "acid_anhydride": "[CX3](=[OX1])[OX2][CX3](=[OX1])",
    "michael_acceptor": "[CX3]=[CX3][CX3](=[OX1])",
    "epoxide": "C1OC1",
    "aziridine": "C1NC1",
    "isocyanate": "[NX2]=[CX2]=[OX1]",
    "nitro_aromatic": "[c,n][N+](=O)[O-]",
    "n_nitroso": "[NX3][NX2]=O",
    "azide": "[NX1]~[NX2]~[NX1]",
    "hydrazine": "[NX3][NX3]",
    "thiourea": "[NX3][CX3](=[SX1])[NX3]",
    "quinone": "O=C1C=CC(=O)C=C1",
    "aniline": "c1ccccc1[NX3;H2,H1]",
    "catechol": "Oc1ccccc1O",
    "furan": "c1ccoc1",
    "thiophene": "c1ccsc1",
}

MET_SOFTSPOTS: Dict[str, str] = {
    "benzylic_carbon": "[cX3][CX4;H2,H3]",
    "allylic_carbon": "[CX3]=[CX3][CX4;H2,H3]",
    "n_methyl": "[NX3][CH3]",
    "n_ethyl": "[NX3][CH2][CH3]",
    "aromatic_methoxy": "c[OX2][CH3]",
    "ester": "[CX3](=O)[OX2][#6]",
    "amide": "[CX3](=O)[NX3]",
    "phenol": "c[OX2H]",
    "carboxylic_acid": "[CX3](=O)[OX2H1]",
}

HERG_FEATURES: Dict[str, str] = {
    "basic_n": "[NX3;H2,H1,H0;!$(NC=O)]",
    "tertiary_amine": "[NX3;H0;!$(NC=O)]([#6])([#6])[#6]",
    "aromatic": "a1aaaaa1",
}

CYP_FEATURES: Dict[str, str] = {
    "imidazole": "c1cnc[nH]1",
    "triazole": "c1ncn[nH]1",
    "pyridine": "n1ccccc1",
    "furan": "c1ccoc1",
    "thiophene": "c1ccsc1",
    "acetylene": "C#C",
    "methylenedioxy": "c1cc2OCOc2cc1",
}

PKA_PATTERNS: List[Tuple[str, str, float, str]] = [
    ("carboxylic_acid", "[CX3](=O)[OX2H1]", 4.5, "acid"),
    ("phenol", "c[OX2H]", 10.0, "acid"),
    ("sulfonic_acid", "[SX4](=O)(=O)[OX2H1]", -1.0, "acid"),
    ("primary_amine", "[NX3;H2;!$(NC=O)]", 10.5, "base"),
    ("secondary_amine", "[NX3;H1;!$(NC=O)]", 10.0, "base"),
    ("tertiary_amine", "[NX3;H0;!$(NC=O)]", 9.0, "base"),
    ("pyridine", "n1ccccc1", 5.2, "base"),
]


def _mol_from_smarts(smarts: str) -> Optional[Chem.Mol]:
    try:
        return Chem.MolFromSmarts(smarts)
    except Exception:
        return None


def has_smarts(mol: Chem.Mol, smarts: str) -> bool:
    patt = _mol_from_smarts(smarts)
    return bool(patt and mol.HasSubstructMatch(patt))


def count_smarts(mol: Chem.Mol, smarts: str) -> int:
    patt = _mol_from_smarts(smarts)
    if not patt:
        return 0
    return len(mol.GetSubstructMatches(patt))


def monotonic(val: float, low: float, high: float) -> float:
    if val <= low:
        return 1.0
    if val >= high:
        return 0.0
    return 1.0 - (val - low) / (high - low)


def estimate_pka(mol: Chem.Mol) -> Tuple[Optional[float], Optional[float], str]:
    acidic = []
    basic = []
    for name, smarts, pka, typ in PKA_PATTERNS:
        if has_smarts(mol, smarts):
            if typ == "acid":
                acidic.append((name, pka))
            else:
                basic.append((name, pka))
    acidic_pka = min([p for _, p in acidic], default=None)
    basic_pka = max([p for _, p in basic], default=None)
    if acidic_pka is not None and basic_pka is not None:
        return basic_pka, acidic_pka, "zwitterion"
    if basic_pka is not None:
        return basic_pka, acidic_pka, "basic"
    if acidic_pka is not None:
        return basic_pka, acidic_pka, "acidic"
    return None, None, "neutral"


def calc_clogd_74(logp: float, basic_pka: Optional[float], acidic_pka: Optional[float]) -> float:
    if basic_pka is not None and basic_pka > 7.4:
        return logp - math.log10(1 + 10 ** (basic_pka - 7.4))
    if acidic_pka is not None and acidic_pka < 7.4:
        return logp - math.log10(1 + 10 ** (7.4 - acidic_pka))
    return logp


def calc_cns_mpo(mw: float, logp: float, clogd: float, tpsa: float, hbd: float, basic_pka: float) -> Dict[str, float | str | bool]:
    mpo_mw = monotonic(mw, 360, 500)
    mpo_logp = monotonic(logp, 3, 5)
    mpo_logd = monotonic(clogd, 2, 4)
    mpo_tpsa = monotonic(tpsa, 40, 90)
    mpo_hbd = monotonic(hbd, 0.5, 3.5)
    mpo_pka = monotonic(basic_pka, 8, 10)
    total = mpo_mw + mpo_logp + mpo_logd + mpo_tpsa + mpo_hbd + mpo_pka

    if total >= 5:
        cls = "Excellent"
    elif total >= 4:
        cls = "Good"
    elif total >= 3:
        cls = "Moderate"
    else:
        cls = "Poor"

    bbb = (total >= 4) and (tpsa < 90) and (mw < 450) and (hbd <= 3)

    return {
        "CNS_MPO": round(total, 2),
        "CNS_MPO_Class": cls,
        "BBB_Likely": bool(bbb),
        "MPO_MW": round(mpo_mw, 3),
        "MPO_LogP": round(mpo_logp, 3),
        "MPO_LogD": round(mpo_logd, 3),
        "MPO_TPSA": round(mpo_tpsa, 3),
        "MPO_HBD": round(mpo_hbd, 3),
        "MPO_pKa": round(mpo_pka, 3),
    }


def calc_ro5(mw: float, logp: float, hbd: float, hba: float) -> Dict[str, int | bool]:
    v = int(mw > 500) + int(logp > 5) + int(hbd > 5) + int(hba > 10)
    return {
        "RO5_Pass": bool(v <= 1),
        "RO5_Violations": v,
        "RO5_MW": bool(mw <= 500),
        "RO5_LogP": bool(logp <= 5),
        "RO5_HBD": bool(hbd <= 5),
        "RO5_HBA": bool(hba <= 10),
    }


def calc_ro3(mw: float, logp: float, hbd: float, hba: float, rotb: float, tpsa: float) -> Dict[str, int | bool]:
    v = int(mw > 300) + int(logp > 3) + int(hbd > 3) + int(hba > 3) + int(rotb > 3)
    return {
        "RO3_Pass": bool(v == 0),
        "RO3_Violations": v,
        "RO3_MW": bool(mw <= 300),
        "RO3_LogP": bool(logp <= 3),
        "RO3_HBD": bool(hbd <= 3),
        "RO3_HBA": bool(hba <= 3),
        "RO3_RotBonds": bool(rotb <= 3),
        "RO3_TPSA": bool(tpsa <= 60),
    }


def calc_oral_filters(mw: float, logp: float, tpsa: float, rotb: float, hba: float, hbd: float, heavy: float, mr: float, rings: float) -> Dict[str, bool | int]:
    veber = (tpsa <= 140) and (rotb <= 10)
    egan = (tpsa <= 132) and (logp <= 5.88)
    ghose = (160 <= mw <= 480) and (-0.4 <= logp <= 5.6) and (40 <= mr <= 130) and (20 <= heavy <= 130)
    muegge = (200 <= mw <= 600) and (-2 <= logp <= 5) and (tpsa <= 150) and (rings <= 7) and (hba <= 10) and (hbd <= 5) and (rotb <= 15)
    score = int(veber) + int(egan) + int(ghose) + int(muegge)
    return {
        "Veber_Pass": bool(veber),
        "Egan_Pass": bool(egan),
        "Ghose_Pass": bool(ghose),
        "Muegge_Pass": bool(muegge),
        "Bioavail_Score": score,
    }


def calc_medchem_flags(mol: Chem.Mol) -> Dict[str, bool | str | int | None]:
    out: Dict[str, bool | str | int | None] = {}

    def _catalog_hit(catalog_name: FilterCatalogParams.FilterCatalogs) -> Tuple[Optional[bool], Optional[str]]:
        try:
            params = FilterCatalogParams()
            params.AddCatalog(catalog_name)
            cat = FilterCatalog(params)
            m = cat.GetFirstMatch(mol)
            if m is None:
                return False, None
            return True, m.GetDescription()
        except Exception:
            return None, None

    pains, pains_desc = _catalog_hit(FilterCatalogParams.FilterCatalogs.PAINS)
    brenk, _ = _catalog_hit(FilterCatalogParams.FilterCatalogs.BRENK)
    nih, _ = _catalog_hit(FilterCatalogParams.FilterCatalogs.NIH)
    zinc, _ = _catalog_hit(FilterCatalogParams.FilterCatalogs.ZINC)

    out["PAINS"] = pains
    out["PAINS_Type"] = pains_desc
    out["BRENK"] = brenk
    out["NIH"] = nih
    out["ZINC"] = zinc

    total = 0
    for k in ["PAINS", "BRENK", "NIH", "ZINC"]:
        if out.get(k) is True:
            total += 1
    out["MedChem_Alerts"] = total

    return out


def calc_toxicity(mol: Chem.Mol) -> Dict[str, int | str]:
    out: Dict[str, int | str] = {}
    families = 0
    for name, smarts in TOX_ALERTS.items():
        c = count_smarts(mol, smarts)
        out[f"Tox_{name}"] = c
        if c > 0:
            families += 1
    out["Tox_Total_Alerts"] = families
    if families >= 3:
        out["Tox_Risk"] = "High"
    elif families >= 1:
        out["Tox_Risk"] = "Medium"
    else:
        out["Tox_Risk"] = "Low"
    return out


def calc_metabolism(mol: Chem.Mol, logp: float) -> Dict[str, int | bool | str]:
    out: Dict[str, int | bool | str] = {}
    total_sites = 0
    for name, smarts in MET_SOFTSPOTS.items():
        c = count_smarts(mol, smarts)
        out[f"Met_{name}"] = c
        total_sites += c
    out["Met_Total_Sites"] = total_sites

    out["Met_Glucuronidation"] = bool(has_smarts(mol, MET_SOFTSPOTS["phenol"]) or has_smarts(mol, MET_SOFTSPOTS["carboxylic_acid"]))
    out["Met_Sulfation"] = bool(has_smarts(mol, MET_SOFTSPOTS["phenol"]))
    out["Met_Acetylation"] = bool(has_smarts(mol, "[NX3;H2;!$(NC=O)]"))
    out["Met_Methylation"] = bool(has_smarts(mol, "c[OX2H]"))

    if total_sites > 5 and logp > 3:
        out["Met_Stability"] = "Low"
    elif total_sites > 3 or logp > 4:
        out["Met_Stability"] = "Medium"
    else:
        out["Met_Stability"] = "High"

    return out


def calc_herg(mol: Chem.Mol, mw: float, logp: float) -> Dict[str, bool | int | str]:
    out: Dict[str, bool | int | str] = {}
    basic_n = has_smarts(mol, HERG_FEATURES["basic_n"])
    tertiary = has_smarts(mol, HERG_FEATURES["tertiary_amine"])
    aromatic = has_smarts(mol, HERG_FEATURES["aromatic"])

    out["hERG_basic_nitrogen"] = bool(basic_n)
    out["hERG_tertiary_amine"] = bool(tertiary)
    out["hERG_aromatic_ring"] = bool(aromatic)

    alerts = int(basic_n) + int(tertiary) + int(aromatic)
    out["hERG_Alerts"] = alerts

    if basic_n and aromatic and logp > 3.5 and mw > 350:
        out["hERG_Risk"] = "High"
    elif basic_n and logp > 3:
        out["hERG_Risk"] = "Medium"
    else:
        out["hERG_Risk"] = "Low"

    return out


def calc_cyp(mol: Chem.Mol, logp: float) -> Dict[str, bool | int]:
    out: Dict[str, bool | int] = {}

    imid = has_smarts(mol, CYP_FEATURES["imidazole"])
    triaz = has_smarts(mol, CYP_FEATURES["triazole"])
    basic_n = has_smarts(mol, "[NX3;H0,H1,H2;!$(NC=O)]")
    aromatic = has_smarts(mol, "a1aaaaa1")

    out["CYP3A4_Inhibitor"] = bool(imid or triaz)
    out["CYP2D6_Inhibitor"] = bool(basic_n and aromatic and logp > 2)
    out["CYP2C9_Inhibitor"] = bool(has_smarts(mol, "[CX3](=O)[OX2H1]"))

    mbi = has_smarts(mol, CYP_FEATURES["acetylene"]) or has_smarts(mol, CYP_FEATURES["furan"]) or has_smarts(mol, CYP_FEATURES["thiophene"]) or has_smarts(mol, CYP_FEATURES["methylenedioxy"])
    out["CYP_MBI_Risk"] = bool(mbi)

    total = int(out["CYP3A4_Inhibitor"]) + int(out["CYP2D6_Inhibitor"]) + int(out["CYP2C9_Inhibitor"]) + int(out["CYP_MBI_Risk"])
    out["CYP_Total_Inhibition"] = total

    return out


# ---------- v5 additions: solubility, permeability, PK, lead metrics ----------

def calc_solubility_esol(mw: float, logp: float, rotb: float, arom_rings: float) -> Dict[str, object]:
    # ESOL (Delaney): LogS = 0.16 - 0.63*LogP - 0.0062*MW + 0.066*RB - 0.74*AP
    logS = 0.16 - 0.63 * logp - 0.0062 * mw + 0.066 * rotb - 0.74 * arom_rings
    mgL = (10 ** logS) * mw * 1000
    if logS >= 0:
        cls = "Highly Soluble"
    elif logS >= -2:
        cls = "Soluble"
    elif logS >= -4:
        cls = "Moderately Soluble"
    elif logS >= -6:
        cls = "Poorly Soluble"
    else:
        cls = "Insoluble"
    return {
        "Sol_ESOL_LogS": round(logS, 3),
        "Sol_ESOL_mg_L": round(mgL, 2),
        "Sol_ESOL_Class": cls,
    }


def calc_permeability(mw: float, logp: float, tpsa: float, hbd: float, hba: float) -> Dict[str, object]:
    # Very simple screening estimators.
    log_papp_caco2 = -5.47 + 0.58 * logp - 0.01 * tpsa - 0.1 * hbd
    caco2_class = "High" if log_papp_caco2 > -5.15 else "Low"

    log_pe_pampa = -4.8 + 0.45 * logp - 0.008 * tpsa
    pampa_class = "High" if log_pe_pampa > -5.0 else "Low"

    gi = "High" if (tpsa <= 140 and logp <= 5.5) else "Low"

    log_kp = -2.7 + 0.71 * logp - 0.0061 * mw

    # P-gp likelihood heuristic
    pgp_score = 0
    if mw > 400:
        pgp_score += 1
    if hbd > 1:
        pgp_score += 1
    if logp > 2:
        pgp_score += 1
    if hba > 4:
        pgp_score += 1

    if pgp_score >= 3:
        pgp = "High"
    elif pgp_score == 2:
        pgp = "Medium"
    else:
        pgp = "Low"

    return {
        "Perm_Caco2_LogPapp": round(log_papp_caco2, 3),
        "Perm_Caco2_Class": caco2_class,
        "Perm_PAMPA_LogPe": round(log_pe_pampa, 3),
        "Perm_PAMPA_Class": pampa_class,
        "Perm_GI_Absorption": gi,
        "Perm_Skin_LogKp": round(log_kp, 3),
        "Perm_Pgp_Likelihood": pgp,
    }


def calc_pk_heuristics(mw: float, logp: float, tpsa: float, met_sites: int) -> Dict[str, object]:
    # Protein binding heuristic driven by LogP.
    ppb = min(99.9, max(0.0, 85 + 3 * logp))
    if ppb >= 95:
        ppb_cls = "High"
    elif ppb >= 80:
        ppb_cls = "Medium"
    else:
        ppb_cls = "Low"

    # Vd heuristic.
    vd = 0.5 + 0.2 * max(0.0, logp)
    if vd < 0.7:
        vd_cls = "Low (plasma confined)"
    elif vd < 2.0:
        vd_cls = "Medium"
    else:
        vd_cls = "High (tissue distributed)"

    # Clearance risk heuristic.
    risk_score = 0
    if logp > 3:
        risk_score += 2
    elif logp > 2:
        risk_score += 1
    if met_sites > 5:
        risk_score += 2
    elif met_sites > 3:
        risk_score += 1
    if mw > 500:
        risk_score += 1

    if risk_score >= 4:
        cl = "High"
    elif risk_score >= 2:
        cl = "Medium"
    else:
        cl = "Low"

    if cl == "High":
        t_half = "Short (<4h)"
    elif cl == "Medium":
        t_half = "Medium (4-12h)"
    else:
        t_half = "Long (>12h)"

    # Oral bioavailability class (coarse): penalize high TPSA, high MW, high logP.
    f_score = 0
    f_score += int(mw <= 500)
    f_score += int(tpsa <= 140)
    f_score += int(logp <= 5)
    f_score += int(tpsa <= 90)
    if f_score >= 4:
        f_cls, f_est = "High", ">70"
    elif f_score >= 3:
        f_cls, f_est = "Moderate", "50-70"
    elif f_score >= 2:
        f_cls, f_est = "Low-Moderate", "30-50"
    else:
        f_cls, f_est = "Low", "<30"

    return {
        "PPB_Est_%": round(ppb, 1),
        "PPB_FreeEst_%": round(100 - ppb, 1),
        "PPB_Class": ppb_cls,
        "Vd_Est_L_kg": round(vd, 2),
        "Vd_Class": vd_cls,
        "CL_SoftSpots": int(met_sites),
        "CL_Risk": cl,
        "T_half_Est": t_half,
        "F_Score": int(f_score),
        "F_Est_%": f_est,
        "F_Class": f_cls,
    }


def calc_lead_metrics(pIC50: float, heavy_atoms: float, logp: float, mw: float, tpsa: float) -> Dict[str, object]:
    # These are only meaningful if pIC50 corresponds to measured potency.
    if heavy_atoms <= 0:
        heavy_atoms = 1
    le = pIC50 / heavy_atoms
    lle = pIC50 - logp
    lelp = logp / le if le > 0 else None
    bei = pIC50 / (mw / 1000) if mw > 0 else None
    sei = pIC50 / (tpsa / 100) if tpsa > 0 else None

    return {
        "Lead_LE": round(le, 3),
        "Lead_LLE": round(lle, 3),
        "Lead_LELP": round(lelp, 3) if lelp is not None else None,
        "Lead_BEI": round(bei, 3) if bei is not None else None,
        "Lead_SEI": round(sei, 3) if sei is not None else None,
        "Lead_LE_Quality": "Good" if le >= 0.3 else "Poor",
        "Lead_LLE_Quality": "Excellent" if lle >= 7 else "Good" if lle >= 5 else "Poor",
        "Lead_LELP_Quality": "Good" if (lelp is not None and lelp < 10) else "Poor",
    }


def safe_num(x):
    try:
        if x is None:
            return None
        fx = float(x)
        if math.isnan(fx) or math.isinf(fx):
            return None
        return fx
    except Exception:
        return None


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

            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                continue
            items.append((name, smiles))
    return items


def compute(
    smiles: str,
    pIC50: float,
    prep_cfg: Optional[PrepConfig] = None,
    ph: float = 7.4,
    ionization_mode: str = "heuristic",
    tautomer_cfg: Optional[TautomerConfig] = None,
    stereo_cfg: Optional[StereoConfig] = None,
    protomer_cfg: Optional[ProtomerConfig] = None,
    three_d_cfg: Optional[ThreeDConfig] = None,
) -> Optional[Dict[str, object]]:
    mol_in = Chem.MolFromSmiles(smiles)
    if mol_in is None:
        return None

    input_can = Chem.MolToSmiles(mol_in, canonical=True)

    mol = mol_in
    prep_report = None
    if prep_cfg is not None:
        mol, prep_report = prepare_mol(mol_in, prep_cfg)

    prep_can = Chem.MolToSmiles(mol, canonical=True)

    calc_mol = mol
    taut_rep = None
    if tautomer_cfg is not None and tautomer_cfg.mode == "enumerate":
        calc_mol, taut_rep = select_representative_tautomer(
            mol,
            max_variants=tautomer_cfg.max_variants,
            top_k=tautomer_cfg.top_k,
        )

    # Task 4: stereo handling (keep/strip/enumerate)
    stereo_rep = None
    if stereo_cfg is not None:
        calc_mol, stereo_rep = apply_stereo(calc_mol, stereo_cfg)

    calc_base_smiles = Chem.MolToSmiles(calc_mol, canonical=True, isomericSmiles=True)

    # pKa estimation on the tautomer+stereo-selected structure (pre-protomer)
    basic_pka, acidic_pka, pka_type = estimate_pka(calc_mol)

    prot_rep = None
    if protomer_cfg is not None and protomer_cfg.mode.lower() == "dimorphite":
        target_charge = 0
        try:
            try:
                from .ionization import estimate_net_charge
            except Exception:
                from ionization import estimate_net_charge

            target_charge = int(round(estimate_net_charge(basic_pka, acidic_pka, ph)))
        except Exception:
            target_charge = 0

        prot_mol, prot_rep = protomerize_smiles(
            calc_base_smiles,
            protomer_cfg,
            target_charge=target_charge,
            top_k=5,
        )

        if protomer_cfg.calc_on_protomer and prot_mol is not None:
            calc_mol = prot_mol

    can = Chem.MolToSmiles(calc_mol, canonical=True, isomericSmiles=True)

    mw = float(Descriptors.MolWt(calc_mol))
    logp = float(Crippen.MolLogP(calc_mol))
    tpsa = float(rdMolDescriptors.CalcTPSA(calc_mol))
    hbd = float(rdMolDescriptors.CalcNumHBD(calc_mol))
    hba = float(rdMolDescriptors.CalcNumHBA(calc_mol))
    rotb = float(Descriptors.NumRotatableBonds(calc_mol))
    rings = float(rdMolDescriptors.CalcNumRings(calc_mol))
    arom = float(rdMolDescriptors.CalcNumAromaticRings(calc_mol))
    heavy = float(calc_mol.GetNumHeavyAtoms())
    fracsp3 = float(rdMolDescriptors.CalcFractionCSP3(calc_mol))
    mr = float(Descriptors.MolMR(calc_mol))

    basic_pka_for_mpo = basic_pka if basic_pka is not None else 7.0
    clogd = calc_clogd_74(logp, basic_pka, acidic_pka)

    out: Dict[str, object] = {}

    out["Input_Canonical_SMILES"] = input_can
    out["Canonical_SMILES"] = prep_can
    out["Calc_Canonical_SMILES"] = Chem.MolToSmiles(calc_mol, canonical=True, isomericSmiles=True)
    out["Calc_Base_SMILES"] = calc_base_smiles

    # Task 4 stereo metadata
    if stereo_cfg is not None:
        out["Stereo_Mode"] = stereo_cfg.mode
        if stereo_rep is not None:
            out["Stereo_IsomerCount"] = stereo_rep.count
            out["Stereo_Selected_SMILES"] = stereo_rep.selected_smiles
            out["Stereo_TopK_SMILES"] = stereo_rep.topk_smiles
            out["Stereo_Changed"] = stereo_rep.changed
            out["Stereo_Notes"] = stereo_rep.notes
        else:
            out["Stereo_IsomerCount"] = 1
            out["Stereo_Selected_SMILES"] = calc_base_smiles
            out["Stereo_TopK_SMILES"] = calc_base_smiles
            out["Stereo_Changed"] = False
            out["Stereo_Notes"] = ""
    if protomer_cfg is not None:
        out["Protomer_Mode"] = protomer_cfg.mode
        out["Protomer_pH_Min"] = getattr(protomer_cfg, "ph_min", None)
        out["Protomer_pH_Max"] = getattr(protomer_cfg, "ph_max", None)
        out["Protomer_Precision"] = getattr(protomer_cfg, "precision", None)
        out["Protomer_MaxVariants"] = getattr(protomer_cfg, "max_variants", None)
        out["Protomer_CalcOnProtomer"] = bool(getattr(protomer_cfg, "calc_on_protomer", False))

        if prot_rep is not None:
            out["Protomer_Count"] = prot_rep.count
            out["Protomer_Selected_SMILES"] = prot_rep.selected_smiles
            out["Protomer_Selected_Charge"] = prot_rep.selected_charge
            out["Protomer_Target_Charge"] = prot_rep.target_charge
            out["Protomer_TopK_SMILES"] = prot_rep.topk_smiles
            out["Protomer_Notes"] = prot_rep.notes
        else:
            out["Protomer_Count"] = 1
            out["Protomer_Selected_SMILES"] = calc_base_smiles
            out["Protomer_Selected_Charge"] = None
            out["Protomer_Target_Charge"] = None
            out["Protomer_TopK_SMILES"] = calc_base_smiles
            out["Protomer_Notes"] = ""

    if tautomer_cfg is not None:
        out["Tautomer_Mode"] = tautomer_cfg.mode
        if tautomer_cfg.mode == "enumerate" and taut_rep is not None:
            out["Tautomer_Count"] = taut_rep.count
            out["Tautomer_Selected_SMILES"] = taut_rep.selected_smiles
            out["Tautomer_TopK_SMILES"] = taut_rep.topk_smiles
            out["Tautomer_Changed"] = taut_rep.changed
        else:
            out["Tautomer_Count"] = 1
            out["Tautomer_Selected_SMILES"] = prep_can
            out["Tautomer_TopK_SMILES"] = prep_can
            out["Tautomer_Changed"] = False

    if prep_report is not None:
        out["Prep_Changed"] = prep_report.changed
        out["Prep_SaltsStripped"] = prep_report.salts_stripped
        out["Prep_Uncharged"] = prep_report.uncharged
        out["Prep_Reionized"] = prep_report.reionized
        out["Prep_TautomerChanged"] = prep_report.tautomer_changed
        out["Prep_StereoStatus"] = prep_report.stereo_status
        out["Prep_UnassignedChiral"] = prep_report.unassigned_chiral_centers
        out["Prep_UnassignedStereoBonds"] = prep_report.unassigned_stereo_bonds
        out["Prep_Notes"] = prep_report.notes
    out["MolWt"] = round(mw, 3)
    out["LogP"] = round(logp, 3)
    out["TPSA"] = round(tpsa, 3)
    out["HBD"] = int(hbd)
    out["HBA"] = int(hba)
    out["RotBonds"] = int(rotb)
    out["Rings"] = int(rings)
    out["AromaticRings"] = int(arom)
    out["HeavyAtoms"] = int(heavy)
    out["FractionCSP3"] = round(fracsp3, 3)

    # Task 6: synthetic accessibility + complexity metrics
    out.update(calc_sa_and_complexity(calc_mol))

    # Task 7: optional 3D conformer-based descriptors
    if three_d_cfg is not None and getattr(three_d_cfg, "enabled", False):
        try:
            out.update(compute_3d_descriptors(calc_mol, three_d_cfg))
        except Exception:
            # Keep calculator robust; schema validator can catch missing columns.
            out.update(compute_3d_descriptors(calc_mol, ThreeDConfig(enabled=False)))

    out["pKa_Basic"] = safe_num(basic_pka)
    out["pKa_Acidic"] = safe_num(acidic_pka)
    out["pKa_Type"] = pka_type
    out["cLogD_7.4"] = round(clogd, 3)

    out["Ion_Mode"] = ionization_mode

    if ionization_mode and ionization_mode.lower() in ("heuristic", "dimorphite"):
        out.update(ionization_features(logp, basic_pka, acidic_pka, ph=ph))

    out.update(calc_cns_mpo(mw, logp, clogd, tpsa, hbd, basic_pka_for_mpo))
    out.update(calc_oral_filters(mw, logp, tpsa, rotb, hba, hbd, heavy, mr, rings))
    out.update(calc_ro5(mw, logp, hbd, hba))
    out.update(calc_ro3(mw, logp, hbd, hba, rotb, tpsa))

    out.update(calc_medchem_flags(calc_mol))
    out.update(calc_toxicity(calc_mol))
    out.update(calc_metabolism(calc_mol, logp))
    out.update(calc_herg(calc_mol, mw, logp))
    out.update(calc_cyp(calc_mol, logp))

    out.update(calc_solubility_esol(mw, logp, rotb, arom))
    out.update(calc_permeability(mw, logp, tpsa, hbd, hba))

    met_sites = int(out.get("Met_Total_Sites", 0) or 0)
    out.update(calc_pk_heuristics(mw, logp, tpsa, met_sites))
    out.update(calc_lead_metrics(pIC50, heavy, logp, mw, tpsa))

    try:
        out["QED"] = round(QED.qed(calc_mol), 4)
    except Exception:
        out["QED"] = None

    # Task 8: developability indices / filters
    try:
        out.update(developability_features(out))
    except Exception:
        # Never fail the calculator due to a summary feature.
        pass

    return out


def write_csv(rows: List[Dict[str, object]], out_path: str, fill_missing: Optional[float]) -> None:
    if not rows:
        raise SystemExit("No rows to write")

    preferred = [
        "Compound_ID",
        "SMILES",
        "Input_Canonical_SMILES",
        "Canonical_SMILES",
        "Prep_Changed",
        "Prep_SaltsStripped",
        "Prep_Uncharged",
        "Prep_Reionized",
        "Prep_TautomerChanged",
        "Prep_StereoStatus",
        "Prep_UnassignedChiral",
        "Prep_UnassignedStereoBonds",
        "Prep_Notes",
        "Ion_Mode",
        "Ion_pH",
        "Ion_State",
        "Ion_NetCharge_Est",
        "Ion_NetCharge_Round",
        "Ion_FracProtonatedBase",
        "Ion_FracDeprotonatedAcid",
        "Ion_LogD_pH_Est",
        "Tautomer_Mode",
        "Tautomer_Count",
        "Tautomer_Selected_SMILES",
        "Tautomer_TopK_SMILES",
        "Tautomer_Changed",
        "Stereo_Mode",
        "Stereo_IsomerCount",
        "Stereo_Selected_SMILES",
        "Stereo_TopK_SMILES",
        "Stereo_Changed",
        "Stereo_Notes",
        "Calc_Base_SMILES",
        "Protomer_Mode",
        "Protomer_pH_Min",
        "Protomer_pH_Max",
        "Protomer_Precision",
        "Protomer_MaxVariants",
        "Protomer_Count",
        "Protomer_Selected_SMILES",
        "Protomer_Selected_Charge",
        "Protomer_Target_Charge",
        "Protomer_TopK_SMILES",
        "Protomer_CalcOnProtomer",
        "Protomer_Notes",
        "Calc_Canonical_SMILES",
        "CNS_MPO",
        "CNS_MPO_Class",
        "BBB_Likely",
        "Bioavail_Score",
        "Veber_Pass",
        "Egan_Pass",
        "Ghose_Pass",
        "Muegge_Pass",
        "RO5_Pass",
        "RO5_Violations",
        "RO3_Pass",
        "RO3_Violations",
        "PAINS",
        "PAINS_Type",
        "BRENK",
        "NIH",
        "ZINC",
        "MedChem_Alerts",
        "Tox_Risk",
        "Tox_Total_Alerts",
        "Met_Stability",
        "Met_Total_Sites",
        "hERG_Risk",
        "CYP_Total_Inhibition",
        "Sol_ESOL_LogS",
        "Sol_ESOL_Class",
        "Perm_GI_Absorption",
        "Perm_Pgp_Likelihood",
        "PPB_Class",
        "CL_Risk",
        "T_half_Est",
        "F_Class",
        "MolWt",
        "LogP",
        "cLogD_7.4",
        "TPSA",
        "HBD",
        "HBA",
        "RotBonds",
        "Rings",
        "AromaticRings",
        "FractionCSP3",
        "SA_Score",
        "SA_Class",
        "Complexity_BertzCT",
        "Complexity_BalabanJ",
        "Complexity_HallKierAlpha",
        "Complexity_NumSpiroAtoms",
        "Complexity_NumBridgeheadAtoms",
        "Complexity_NumStereoCenters",
        "Complexity_NumUnspecifiedStereoCenters",
        "Complexity_NumHeteroAtoms",
        "Complexity_AromaticAtomFrac",
        "pKa_Basic",
        "pKa_Acidic",
        "QED",
        # Task 8: developability indices (optional but cheap)
        "Dev_Score",
        "Dev_Class",
        "Dev_Risk_Class",
        "Dev_Risk_Burden",
        "Dev_BasePoints",
        "Dev_GSK_4_400_Pass",
        "Dev_Pfizer_3_75_Pass",
        "Dev_GoldenTriangle_Pass",
        "Dev_Solubility_Risk",
        "Dev_Flags",

        # Task 7: 3D shape descriptors (optional)
        "3D_Mode",
        "3D_NumConfs",
        "3D_SelectedConfId",
        "3D_Minimizer",
        "3D_Converged",
        "3D_BestEnergy",
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
    ]

    keys = set()
    for r in rows:
        keys.update(r.keys())

    header: List[str] = []
    for c in preferred:
        if c in keys and c not in header:
            header.append(c)

    for c in sorted(keys):
        if c not in header:
            header.append(c)

    if fill_missing is not None:
        for r in rows:
            for c in header:
                v = r.get(c)
                if v is None:
                    r[c] = fill_missing
                elif isinstance(v, float) and (math.isnan(v) or math.isinf(v)):
                    r[c] = fill_missing

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=header)
        w.writeheader()
        w.writerows(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description="MolProp Toolkit calculator (v5)")
    ap.add_argument("input", help="Input SMILES file")
    ap.add_argument("-o", "--output", default=None, help="Output CSV")
    ap.add_argument("--fill-missing", type=float, default=None, help="Fill missing values with constant, e.g. 0")
    ap.add_argument("--pIC50", type=float, default=5.0, help="Optional potency (pIC50) used for lead metrics")

    # Task 7: 3D descriptor mode (optional)
    ap.add_argument("--3d", dest="enable_3d", action="store_true", help="Compute ETKDG conformers + 3D shape descriptors (v5 only)")
    ap.add_argument("--3d-num-confs", dest="d3_num_confs", type=int, default=10, help="Number of conformers to embed")
    ap.add_argument("--3d-prune-rms", dest="d3_prune_rms", type=float, default=0.5, help="RMS prune threshold (Angstrom) for conformer embedding")
    ap.add_argument("--3d-max-attempts", dest="d3_max_attempts", type=int, default=20, help="Max attempts per conformer during embedding")
    ap.add_argument("--3d-seed", dest="d3_seed", type=int, default=0, help="Random seed for ETKDG (0 = deterministic in RDKit)")
    ap.add_argument(
        "--3d-minimize",
        dest="d3_minimize",
        choices=["mmff", "uff", "none"],
        default="mmff",
        help="Minimizer for conformers before computing 3D descriptors",
    )

    # Structure preparation controls
    ap.add_argument("--no-prep", action="store_true", help="Disable RDKit structure preparation")
    ap.add_argument("--no-tautomer", action="store_true", help="(Deprecated) Disable tautomer canonicalization during preparation")
    ap.add_argument("--keep-charges", action="store_true", help="Disable uncharging during preparation")

    # Task 3: tautomer handling
    ap.add_argument(
        "--tautomer-mode",
        choices=["prep-canonical", "none", "enumerate"],
        default="prep-canonical",
        help="Tautomer handling: prep-canonical uses preparation canonical tautomer; enumerate reports tautomer set and selects representative; none disables tautomer processing.",
    )
    ap.add_argument("--tautomer-max", type=int, default=64, help="Max tautomer variants to enumerate (when --tautomer-mode enumerate)")
    ap.add_argument("--tautomer-topk", type=int, default=5, help="How many tautomer SMILES to include in Tautomer_TopK_SMILES")

    # Task 4: stereochemistry handling
    ap.add_argument(
        "--stereo-mode",
        choices=["keep", "strip", "enumerate"],
        default="keep",
        help="Stereo handling: keep preserves stereo; strip removes stereo; enumerate assigns unassigned stereo by enumerating stereoisomers.",
    )
    ap.add_argument("--stereo-max", type=int, default=32, help="Max stereoisomers to enumerate (when --stereo-mode enumerate)")
    ap.add_argument("--stereo-topk", type=int, default=5, help="How many stereoisomer SMILES to include in Stereo_TopK_SMILES")
    ap.add_argument(
        "--stereo-select",
        choices=["canonical", "first"],
        default="canonical",
        help="How to select representative stereoisomer when enumerating",
    )

    # Ionization / pH controls
    ap.add_argument("--ph", type=float, default=7.4, help="pH for heuristic ionization features (default: 7.4)")
    ap.add_argument(
        "--ionization",
        choices=["none", "heuristic", "dimorphite"],
        default="heuristic",
        help="Ionization handling mode. 'heuristic' adds Ion_* columns; 'dimorphite' enumerates protomers via Dimorphite-DL; 'none' disables.",
    )
    ap.add_argument("--ph-min", type=float, default=None, help="Dimorphite-DL minimum pH (default: ph-0.5)")
    ap.add_argument("--ph-max", type=float, default=None, help="Dimorphite-DL maximum pH (default: ph+0.5)")
    ap.add_argument("--dimorphite-precision", type=float, default=1.0, help="Dimorphite-DL precision (default 1.0)")
    ap.add_argument("--dimorphite-max-variants", type=int, default=64, help="Max protomers generated by Dimorphite-DL")
    ap.add_argument(
        "--protomer-select",
        choices=["closest-charge", "neutral", "first"],
        default="closest-charge",
        help="How to select representative protomer from enumeration",
    )
    ap.add_argument(
        "--calc-on-protomer",
        action="store_true",
        help="If set and --ionization dimorphite, compute descriptors on selected protomer (otherwise only record protomer SMILES)",
    )

    args = ap.parse_args()

    out_path = args.output or (Path(args.input).stem + "_results.csv")

    items = parse_smiles_file(args.input)
    if not items:
        raise SystemExit("No valid SMILES found")

    tautomer_mode = args.tautomer_mode
    if args.no_tautomer and args.tautomer_mode == "prep-canonical":
        tautomer_mode = "none"

    prep_cfg: Optional[PrepConfig]
    if args.no_prep:
        prep_cfg = None
    else:
        prep_cfg = PrepConfig(
            canonicalize_tautomer=(tautomer_mode == "prep-canonical"),
            uncharge=not args.keep_charges,
        )

    stereo_cfg: Optional[StereoConfig] = StereoConfig(
        mode=args.stereo_mode,
        max_isomers=args.stereo_max,
        top_k=args.stereo_topk,
        select=args.stereo_select,
    )

    if args.ionization == "dimorphite":
        ph_min = args.ph_min if args.ph_min is not None else args.ph - 0.5
        ph_max = args.ph_max if args.ph_max is not None else args.ph + 0.5
        prot_cfg: Optional[ProtomerConfig] = ProtomerConfig(
            mode="dimorphite",
            ph_min=ph_min,
            ph_max=ph_max,
            precision=args.dimorphite_precision,
            max_variants=args.dimorphite_max_variants,
            select=args.protomer_select,
            calc_on_protomer=args.calc_on_protomer,
        )
    else:
        prot_cfg = ProtomerConfig(mode="none")

    taut_cfg: Optional[TautomerConfig]
    if tautomer_mode == "none":
        taut_cfg = TautomerConfig(mode="none")
    elif tautomer_mode == "enumerate":
        taut_cfg = TautomerConfig(mode="enumerate", max_variants=args.tautomer_max, top_k=args.tautomer_topk)
    else:
        taut_cfg = TautomerConfig(mode="prep-canonical")

    three_d_cfg: Optional[ThreeDConfig] = None
    if args.enable_3d:
        three_d_cfg = ThreeDConfig(
            enabled=True,
            num_confs=int(args.d3_num_confs),
            prune_rms=float(args.d3_prune_rms),
            max_attempts=int(args.d3_max_attempts),
            random_seed=int(args.d3_seed),
            minimize=str(args.d3_minimize),
        )

    rows: List[Dict[str, object]] = []
    for name, smiles in items:
        props = compute(
            smiles,
            args.pIC50,
            prep_cfg,
            ph=args.ph,
            ionization_mode=args.ionization,
            tautomer_cfg=taut_cfg,
            stereo_cfg=stereo_cfg,
            protomer_cfg=prot_cfg,
            three_d_cfg=three_d_cfg,
        )
        if props is None:
            continue
        props["Compound_ID"] = name
        props["SMILES"] = smiles
        rows.append(props)

    write_csv(rows, out_path, args.fill_missing)
    print(f"Saved: {out_path} (compounds={len(rows)}, columns={len(rows[0])})")


if __name__ == "__main__":
    main()

