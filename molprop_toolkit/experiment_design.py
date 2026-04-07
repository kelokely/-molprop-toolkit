"""Lightweight experiment-design prompts for decision products.

These rules are intentionally heuristic. The goal is not to replace project
planning, but to make the generated outputs answer a more actionable question:
"what is the smallest next test or analog set that would resolve the current
uncertainty?"
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable

import pandas as pd

from molprop_toolkit.decision_presets import get_stage_preset


@dataclass(frozen=True)
class ExperimentPlan:
    hypothesis: str
    minimal_test_set: str
    minimal_analog_set: str


def _first_present(row: pd.Series, names: Iterable[str]) -> Any:
    for name in names:
        if name in row.index:
            value = row.get(name)
            if value is not None and not pd.isna(value):
                return value
    return None


def _numeric(row: pd.Series, *names: str) -> float | None:
    value = _first_present(row, names)
    if value is None:
        return None
    numeric = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    if pd.isna(numeric):
        return None
    return float(numeric)


def _text(row: pd.Series, *names: str) -> str | None:
    value = _first_present(row, names)
    if value is None:
        return None
    return str(value).strip()


def _truthy(row: pd.Series, *names: str) -> bool:
    value = _first_present(row, names)
    if value is None:
        return False
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    return text in {"true", "1", "yes", "y"}


def _risk_level(row: pd.Series, *names: str) -> str | None:
    text = _text(row, *names)
    if not text:
        return None
    lowered = text.lower()
    if lowered in {"high", "medium", "low"}:
        return lowered
    return None


def _append_unique(items: list[str], value: str) -> None:
    if value and value not in items:
        items.append(value)


def suggest_experiment_plan(
    row: pd.Series,
    *,
    queue: str,
    stage: str | None = None,
) -> ExperimentPlan:
    potency_missing = bool(
        "Measured_pIC50" not in row.index or pd.isna(row.get("Measured_pIC50"))
    )
    tox_risk = _risk_level(row, "Tox_Risk")
    herg_risk = _risk_level(row, "hERG_Risk")
    cyp_total = _numeric(row, "CYP_Total_Inhibition")
    cyp_mbi = _truthy(row, "CYP_MBI_Risk")
    medchem_alerts = _numeric(row, "MedChem_Alerts")
    logp = _numeric(row, "LogP")
    logs = _numeric(row, "Sol_Consensus_LogS")
    retro_rank = _numeric(row, "Retro_Rerank_Rank")
    retro_stock = _numeric(row, "Retro_Stock_Coverage")
    retro_steps = _numeric(row, "Retro_StepCount")
    retro_unsourced = _numeric(row, "Retro_NotInStockCount")
    permeability = _numeric(
        row,
        "Caco2_AtoB",
        "PAMPA_Pe",
        "Perm_Caco2_AtoB",
        "Perm_PAMPA",
    )

    route_risk = bool(
        _truthy(row, "Retro_RecordMatched") is False and "Retro_RecordMatched" in row.index
    )
    route_risk = route_risk or bool(
        _truthy(row, "Retro_IsSolved") is False and "Retro_IsSolved" in row.index
    )
    route_risk = route_risk or (
        retro_rank is not None and retro_rank > 20
    )
    route_risk = route_risk or (
        retro_stock is not None and retro_stock < 0.5
    )
    route_risk = route_risk or (
        retro_steps is not None and retro_steps > 6
    )
    route_risk = route_risk or (
        retro_unsourced is not None and retro_unsourced >= 2
    )

    liability_mode = any(
        (
            tox_risk in {"medium", "high"},
            herg_risk in {"medium", "high"},
            cyp_total is not None and cyp_total >= 2,
            cyp_mbi,
            medchem_alerts is not None and medchem_alerts >= 1,
        )
    )

    stage_preset = get_stage_preset(stage)
    stage_label = stage_preset.label if stage_preset is not None else None

    if queue == "make_next":
        hypothesis = (
            "Current balance looks strong enough to justify near-term synthesis execution."
        )
    elif queue == "test_next":
        hypothesis = (
            "A small test set should decide whether this chemistry deserves a make-next round."
        )
    elif queue == "de_risk_next":
        hypothesis = (
            "The next cycle should resolve whether the current liability is intrinsic to the scaffold or rescueable with a targeted chemistry move."
        )
    elif queue == "backup_series":
        hypothesis = (
            "A small orthogonal follow-up set can show whether this chemistry is strong enough to keep alive as a backup path."
        )
    elif queue == "hold":
        hypothesis = (
            "The best next step is to avoid broad expansion until a clearer rescue hypothesis emerges."
        )
    else:
        hypothesis = (
            "Use a small, discriminating next-step set to decide whether this queue should advance or stall."
        )

    if potency_missing:
        hypothesis = (
            "The first discriminating question is whether measured potency actually supports the current computed balance."
        )
    elif route_risk and queue in {"make_next", "backup_series"}:
        hypothesis = (
            "Separate chemistry quality from route burden before committing more synthesis bandwidth."
        )
    elif liability_mode and queue in {"de_risk_next", "hold"}:
        hypothesis = (
            "Use the next cycle to distinguish a rescueable liability from a scaffold-level dead end."
        )

    tests: list[str] = []
    if potency_missing:
        _append_unique(tests, "primary potency confirmation")
    if logs is None and queue in {"test_next", "make_next", "backup_series"}:
        _append_unique(tests, "kinetic/equilibrium solubility panel")
    if logs is not None and logs <= -5:
        _append_unique(tests, "kinetic/equilibrium solubility panel")
    if permeability is None and queue in {"test_next", "backup_series", "de_risk_next"}:
        _append_unique(tests, "Caco-2 or PAMPA permeability check")
    if herg_risk in {"medium", "high"} or (logp is not None and logp >= 3.5 and queue == "de_risk_next"):
        _append_unique(tests, "hERG counterscreen")
    if cyp_total is not None and cyp_total >= 2:
        _append_unique(tests, "CYP inhibition panel")
    if cyp_mbi:
        _append_unique(tests, "time-dependent CYP inhibition follow-up")
    if tox_risk in {"medium", "high"} or (medchem_alerts is not None and medchem_alerts >= 1):
        _append_unique(tests, "orthogonal artifact / alert counterscreens")
    if route_risk:
        _append_unique(tests, "route review with stock check")
    if not tests:
        _append_unique(tests, "focused potency + developability confirmation set")

    if route_risk:
        analog_set = (
            "Do not launch a broad analog round yet. Identify 1-2 route-supported representatives or a shorter route before expansion."
        )
    elif queue == "de_risk_next":
        if herg_risk in {"medium", "high"} or (logp is not None and logp >= 3.5):
            analog_set = (
                "Make 2-3 lower-logP or lower-basicity analogs around the same scaffold to test whether the liability is lipophilicity-driven."
            )
        elif cyp_total is not None and cyp_total >= 2 or cyp_mbi:
            analog_set = (
                "Make 2-3 analogs that soften the likely metabolic handle or reduce exposed lipophilicity before scaling the series."
            )
        else:
            analog_set = (
                "Make 2-3 alert-removal or polarity-tuning analogs that preserve the core potency motif."
            )
    elif queue == "test_next":
        analog_set = (
            "Hold a broad analog round for now; reserve 2-3 follow-up analogs for whichever test result best resolves the current question."
        )
    elif queue == "backup_series":
        analog_set = (
            "Make 2-3 orthogonal representatives that keep the backup path real without mirroring the lead-series liabilities."
        )
    elif queue == "hold":
        analog_set = (
            "Do not start a new analog round until stronger measured evidence or a clearer rescue hypothesis appears."
        )
    else:
        analog_set = (
            "Make 2-4 close analogs that preserve route support while probing one substituent axis at a time."
        )

    if stage_label == "Hit expansion" and queue in {"make_next", "backup_series"}:
        analog_set = (
            "Make 3-4 breadth-preserving analogs so the campaign still learns across chemistry space rather than only tightening one local series."
        )
    elif stage_label == "Lead optimization" and queue == "make_next":
        analog_set = (
            "Make 2-3 potency-normalized analogs that preserve route support while testing one balance-driving substituent change."
        )
    elif stage_label == "Liability de-risking":
        analog_set = (
            "Keep the analog set narrow: 2-3 rescue analogs tied directly to the current liability hypothesis."
        )

    return ExperimentPlan(
        hypothesis=hypothesis,
        minimal_test_set="; ".join(tests[:4]),
        minimal_analog_set=analog_set,
    )
