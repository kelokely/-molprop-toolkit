"""Built-in picklists.

These are intended as useful defaults that work across many projects.
They are heuristic and should be tuned with project context.

All filters use a "missing policy" of warn_skip unless the column is critical.
"""

from __future__ import annotations

from typing import Dict

from .engine import FilterRule, PicklistDefinition, SortRule


BUILTIN_PICKLISTS: Dict[str, PicklistDefinition] = {
    "top_cns_followup": PicklistDefinition(
        key="top_cns_followup",
        name="Top CNS follow-up",
        description=(
            "Prioritize compounds with strong CNS MPO balance and no obvious red flags. "
            "This is intended as a short list for follow-up assays rather than a hard gate."
        ),
        filters=(
            FilterRule("CNS_MPO", ">=", 4.5, missing="warn_skip"),
            FilterRule("RO5_Pass", "==", True, missing="warn_skip"),
            FilterRule("Tox_Risk", "notin", ["High"], missing="warn_skip"),
            FilterRule("hERG_Risk", "notin", ["High"], missing="warn_skip"),
        ),
        sort=(SortRule("CNS_MPO", ascending=False),),
        limit=50,
        include_categories=("cns_mpo", "rule_of_5", "toxicity", "herg"),
    ),
    "oral_low_ddi_risk": PicklistDefinition(
        key="oral_low_ddi_risk",
        name="Oral low-DDI risk",
        description=(
            "Oral-leaning filter emphasizing rule-based oral properties plus low CYP inhibition burden. "
            "Useful for triaging toward cleaner DDI profiles early."
        ),
        filters=(
            FilterRule("Bioavail_Score", ">=", 3, missing="warn_skip"),
            FilterRule("RO5_Violations", "<=", 1, missing="warn_skip"),
            FilterRule("CYP_Total_Inhibition", "<=", 1, missing="warn_skip"),
            FilterRule("CYP_MBI_Risk", "notin", [True, "True", "Yes", "Y"], missing="warn_skip"),
            FilterRule("Tox_Risk", "notin", ["High"], missing="warn_skip"),
        ),
        sort=(
            SortRule("Bioavail_Score", ascending=False),
            SortRule("CYP_Total_Inhibition", ascending=True),
        ),
        limit=100,
        include_categories=("oral_bioavailability", "rule_of_5", "cyp", "toxicity"),
    ),
    "high_risk_deprioritize": PicklistDefinition(
        key="high_risk_deprioritize",
        name="High risk – deprioritize",
        description=(
            "OR-style deprioritization list: any strong red flag (toxicity, hERG, DDI, etc.) "
            "will include the compound. Intended for quick review and cleanup triage, not as a final decision."
        ),
        logic="any",
        filters=(
            FilterRule("Tox_Risk", "in", ["High"], missing="warn_skip"),
            FilterRule("hERG_Risk", "in", ["High"], missing="warn_skip"),
            FilterRule("CYP_Total_Inhibition", ">=", 3, missing="warn_skip"),
            FilterRule("CYP_MBI_Risk", "in", [True, "True", "Yes", "Y"], missing="warn_skip"),
            FilterRule("MedChem_Alerts", ">=", 2, missing="warn_skip"),
            FilterRule("RO5_Violations", ">=", 2, missing="warn_skip"),
            FilterRule("Sol_Consensus_LogS", "<=", -5, missing="warn_skip"),
        ),
        sort=(
            SortRule("Tox_Total_Alerts", ascending=False),
            SortRule("MedChem_Alerts", ascending=False),
        ),
        limit=0,
        include_categories=("toxicity", "herg", "cyp", "medchem_flags", "rule_of_5", "solubility"),
    ),
}


