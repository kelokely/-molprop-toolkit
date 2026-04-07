"""Stage-aware decision presets for report and picklist workflows.

These presets are intentionally lightweight:
- they describe the decision stage in human terms,
- they surface the columns that matter most for that stage,
- and they map the stage to a sensible default subset of built-in picklists.

They do not replace project-specific rules. The goal is to make the public
workflow more campaign-aware without hard-coding a heavy scoring model.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class StagePreset:
    key: str
    label: str
    summary: str
    focus_points: tuple[str, ...]
    recommended_picklists: tuple[str, ...]
    priority_columns: tuple[str, ...]


STAGE_PRESETS: dict[str, StagePreset] = {
    "hit-expansion": StagePreset(
        key="hit-expansion",
        label="Hit expansion",
        summary=(
            "Bias toward tractable compounds that broaden the chemistry space without"
            " walking immediately into obvious liabilities."
        ),
        focus_points=(
            "Prefer compounds that stay easy to make and easy to iterate.",
            "Use the stage to screen out obvious risk before expensive assay expansion.",
            "Keep some breadth so the campaign can still branch into backup ideas.",
        ),
        recommended_picklists=(
            "test_next",
            "de_risk_next",
            "backup_series",
        ),
        priority_columns=(
            "QED",
            "SA_Score",
            "CNS_MPO",
            "Bioavail_Score",
            "RO5_Violations",
            "MedChem_Alerts",
        ),
    ),
    "lead-optimization": StagePreset(
        key="lead-optimization",
        label="Lead optimization",
        summary=(
            "Favor potency-normalized balance and execution-ready compounds that can"
            " advance the current lead hypothesis without ignoring route support."
        ),
        focus_points=(
            "Watch potency-adjusted balance before adding more chemistry complexity.",
            "Use route support and stock coverage to keep synthesis bandwidth focused.",
            "Treat liabilities as tradeoffs to manage, not just reasons to discard ideas.",
        ),
        recommended_picklists=(
            "make_next",
            "de_risk_next",
            "hold",
        ),
        priority_columns=(
            "Measured_pIC50",
            "Lead_LLE",
            "CNS_MPO",
            "Retro_Rerank_Rank",
            "Retro_Stock_Coverage",
            "Retro_StepCount",
        ),
    ),
    "backup-series": StagePreset(
        key="backup-series",
        label="Backup series management",
        summary=(
            "Keep orthogonal series alive by favoring cleaner liability profiles and"
            " compounds that remain realistic to synthesize if the lead series stalls."
        ),
        focus_points=(
            "Look for orthogonal chemotypes that avoid current series failure modes.",
            "Prefer compounds with enough route support to keep a backup path real.",
            "Balance broad coverage against avoidable complexity and risk accumulation.",
        ),
        recommended_picklists=(
            "backup_series",
            "make_next",
            "hold",
        ),
        priority_columns=(
            "Lead_LLE",
            "Bioavail_Score",
            "Retro_Rerank_Rank",
            "Retro_Stock_Coverage",
            "CYP_Total_Inhibition",
            "Tox_Risk",
        ),
    ),
    "liability-de-risking": StagePreset(
        key="liability-de-risking",
        label="Liability de-risking",
        summary=(
            "Focus on compounds and series that can most quickly answer whether a key"
            " liability is intrinsic, rescueable, or best abandoned."
        ),
        focus_points=(
            "Separate chemistry worth rescuing from chemistry that is only adding risk.",
            "Use route burden to avoid spending synthesis bandwidth on weak rescue ideas.",
            "Push the workflow toward discriminating assays and counterscreens.",
        ),
        recommended_picklists=(
            "de_risk_next",
            "hold",
        ),
        priority_columns=(
            "Tox_Total_Alerts",
            "MedChem_Alerts",
            "CYP_Total_Inhibition",
            "Retro_StepCount",
            "Retro_NotInStockCount",
            "Lead_LLE",
        ),
    ),
}


STAGE_CHOICES: tuple[str, ...] = tuple(STAGE_PRESETS.keys())


def get_stage_preset(stage: str | None) -> StagePreset | None:
    if not stage:
        return None
    return STAGE_PRESETS[str(stage)]
