#!/usr/bin/env python3
"""Developability indices and risk summaries (Task 8).

This module adds a few commonly used, interpretable property-space filters plus
an aggregated risk/score output intended for triage.

Design goals:
- Minimal extra dependencies (pure Python).
- Works for both v4 and v5 calculator outputs: it only uses columns that are
  present in the provided feature dict.

Notes on rules implemented:
- GSK 4/400: MW <= 400 and LogP <= 4.
- Pfizer 3/75: LogP <= 3 and TPSA <= 75.
- Golden Triangle (coarse): 200 <= MW <= 500 and 1 <= cLogD_7.4 <= 3.

The composite Dev_Risk_Burden and Dev_Score are intentionally heuristic.
"""

from __future__ import annotations

from typing import Dict, List, Optional


def _as_float(x) -> Optional[float]:
    try:
        if x is None:
            return None
        return float(x)
    except Exception:
        return None


def _as_int(x) -> Optional[int]:
    try:
        if x is None:
            return None
        return int(x)
    except Exception:
        return None


def _risk_points_from_class(v: object) -> int:
    if v is None:
        return 0
    s = str(v).strip().lower()
    if s == "high":
        return 2
    if s == "medium":
        return 1
    return 0


def developability_features(props: Dict[str, object]) -> Dict[str, object]:
    mw = _as_float(props.get("MolWt"))
    logp = _as_float(props.get("LogP"))
    tpsa = _as_float(props.get("TPSA"))
    clogd = _as_float(props.get("cLogD_7.4"))

    qed = _as_float(props.get("QED"))
    sa = _as_float(props.get("SA_Score"))

    # Rule-style filters
    gsk_4_400 = bool(mw is not None and logp is not None and mw <= 400 and logp <= 4)
    pfizer_3_75 = bool(tpsa is not None and logp is not None and tpsa <= 75 and logp <= 3)
    golden_triangle = bool(mw is not None and clogd is not None and 200 <= mw <= 500 and 1 <= clogd <= 3)

    # Solubility risk (if available)
    sol_logS = _as_float(props.get("Sol_ESOL_LogS"))
    sol_risk: Optional[str]
    if sol_logS is None:
        sol_risk = None
    elif sol_logS <= -5:
        sol_risk = "High"
    elif sol_logS <= -3:
        sol_risk = "Medium"
    else:
        sol_risk = "Low"

    # Aggregate risk burden (only counts what is present)
    burden = 0
    flags: List[str] = []

    tox = props.get("Tox_Risk")
    herg = props.get("hERG_Risk")
    cl_risk = props.get("CL_Risk")

    burden += _risk_points_from_class(tox)
    if str(tox).lower() == "high":
        flags.append("tox_high")

    burden += _risk_points_from_class(herg)
    if str(herg).lower() == "high":
        flags.append("herg_high")

    if sol_risk is not None:
        burden += _risk_points_from_class(sol_risk)
        if sol_risk.lower() == "high":
            flags.append("solubility_high_risk")

    if cl_risk is not None:
        burden += _risk_points_from_class(cl_risk)
        if str(cl_risk).lower() == "high":
            flags.append("clearance_high_risk")

    cyp_total = _as_int(props.get("CYP_Total_Inhibition"))
    if cyp_total is not None:
        if cyp_total >= 3:
            burden += 2
            flags.append("cyp_burden_high")
        elif cyp_total == 2:
            burden += 1
            flags.append("cyp_burden_medium")

    ro5_viol = _as_int(props.get("RO5_Violations"))
    if ro5_viol is not None:
        if ro5_viol >= 2:
            burden += 2
            flags.append("ro5_2plus")
        elif ro5_viol == 1:
            burden += 1

    medchem_alerts = _as_int(props.get("MedChem_Alerts"))
    if medchem_alerts is not None:
        if medchem_alerts >= 2:
            burden += 2
            flags.append("medchem_alerts_2plus")
        elif medchem_alerts == 1:
            burden += 1

    # Positive points (only counts what is present)
    base_points = 0
    base_points += int(gsk_4_400)
    base_points += int(pfizer_3_75)
    base_points += int(golden_triangle)
    if qed is not None:
        base_points += int(qed >= 0.5)
    if sa is not None:
        base_points += int(sa <= 3.0)

    bioavail = _as_int(props.get("Bioavail_Score"))
    if bioavail is not None:
        base_points += int(bioavail >= 3)

    score = max(0, base_points - burden)

    if burden >= 6:
        risk_class = "High"
    elif burden >= 3:
        risk_class = "Medium"
    else:
        risk_class = "Low"

    if score >= 4 and risk_class == "Low":
        dev_class = "High"
    elif score >= 2 and risk_class in ("Low", "Medium"):
        dev_class = "Medium"
    else:
        dev_class = "Low"

    return {
        "Dev_GSK_4_400_Pass": gsk_4_400,
        "Dev_Pfizer_3_75_Pass": pfizer_3_75,
        "Dev_GoldenTriangle_Pass": golden_triangle,
        "Dev_Solubility_Risk": sol_risk,
        "Dev_BasePoints": int(base_points),
        "Dev_Risk_Burden": int(burden),
        "Dev_Risk_Class": risk_class,
        "Dev_Score": int(score),
        "Dev_Class": dev_class,
        "Dev_Flags": ",".join(flags),
    }

