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
        queue="test_next",
        next_action="Prioritize for CNS-focused follow-up assays and confirmatory profiling.",
        rationale="Balanced CNS-leaning shortlist with enough signal quality to justify additional testing.",
        evidence_columns=("CNS_MPO", "Lead_LLE", "Tox_Risk", "hERG_Risk"),
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
        queue="backup_series",
        next_action="Keep as an oral-leaning or backup-series option while DDI risk remains low.",
        rationale="Cleaner oral profile that can support broader campaign branching without obvious DDI burden.",
        evidence_columns=(
            "Bioavail_Score",
            "RO5_Violations",
            "CYP_Total_Inhibition",
            "Tox_Risk",
        ),
        filters=(
            FilterRule("Bioavail_Score", ">=", 3, missing="warn_skip"),
            FilterRule("RO5_Violations", "<=", 1, missing="warn_skip"),
            FilterRule("CYP_Total_Inhibition", "<=", 1, missing="warn_skip"),
            FilterRule(
                "CYP_MBI_Risk", "notin", [True, "True", "Yes", "Y"], missing="warn_skip"
            ),
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
        queue="hold",
        next_action="Hold expansion until the liabilities are resolved or the chemistry case materially improves.",
        rationale="Liability-heavy compound that should not stay in the front of the active queue without new evidence.",
        evidence_columns=(
            "Tox_Risk",
            "hERG_Risk",
            "CYP_Total_Inhibition",
            "MedChem_Alerts",
            "Sol_Consensus_LogS",
        ),
        logic="any",
        filters=(
            FilterRule("Tox_Risk", "in", ["High"], missing="warn_skip"),
            FilterRule("hERG_Risk", "in", ["High"], missing="warn_skip"),
            FilterRule("CYP_Total_Inhibition", ">=", 3, missing="warn_skip"),
            FilterRule(
                "CYP_MBI_Risk", "in", [True, "True", "Yes", "Y"], missing="warn_skip"
            ),
            FilterRule("MedChem_Alerts", ">=", 2, missing="warn_skip"),
            FilterRule("RO5_Violations", ">=", 2, missing="warn_skip"),
            FilterRule("Sol_Consensus_LogS", "<=", -5, missing="warn_skip"),
        ),
        sort=(
            SortRule("Tox_Total_Alerts", ascending=False),
            SortRule("MedChem_Alerts", ascending=False),
        ),
        limit=0,
        include_categories=(
            "toxicity",
            "herg",
            "cyp",
            "medchem_flags",
            "rule_of_5",
            "solubility",
        ),
    ),
    "make_next_synthesis_ready": PicklistDefinition(
        key="make_next_synthesis_ready",
        name="Make next – synthesis ready",
        description=(
            "Shortlist compounds that already look route-supported in the retro-enriched table. "
            "This is intended for quick make-next prioritization when route support and tractable burden matter."
        ),
        queue="make_next",
        next_action="Advance to synthesis planning or immediate make-next review.",
        rationale="Synthesis-ready compound with route support already in hand.",
        evidence_columns=(
            "Retro_Rerank_Rank",
            "Retro_Stock_Coverage",
            "Retro_StepCount",
            "Lead_LLE",
            "CNS_MPO",
        ),
        filters=(
            FilterRule("Retro_RecordMatched", "==", True, missing="warn_skip"),
            FilterRule("Retro_IsSolved", "==", True, missing="warn_skip"),
            FilterRule("Retro_Rerank_Rank", "<=", 20, missing="warn_skip"),
            FilterRule("Retro_Stock_Coverage", ">=", 0.5, missing="warn_skip"),
            FilterRule("Retro_StepCount", "<=", 6, missing="warn_skip"),
        ),
        sort=(
            SortRule("Retro_Rerank_Rank", ascending=True),
            SortRule("Retro_StepCount", ascending=True),
            SortRule("Retro_Stock_Coverage", ascending=False),
        ),
        diversity_columns=("Scaffold_ID", "Series_ID", "Cluster_ID"),
        max_per_diversity_group=2,
        limit=30,
        output_columns=(
            "Retro_Rerank_Rank",
            "Retro_Rerank_Score",
            "Retro_Stock_Coverage",
            "Retro_StepCount",
            "Retro_TopScore",
            "Lead_LLE",
            "CNS_MPO",
        ),
    ),
    "route_risk_deprioritize": PicklistDefinition(
        key="route_risk_deprioritize",
        name="Route risk – deprioritize",
        description=(
            "Catch compounds that look difficult to support synthetically so teams can separate route risk from chemistry quality."
        ),
        queue="de_risk_next",
        next_action="Use targeted route review or route-aware counterscreening before spending more chemistry bandwidth.",
        rationale="Promising chemistry is still exposed to route risk that should be resolved explicitly.",
        evidence_columns=(
            "Retro_IsSolved",
            "Retro_Rerank_Rank",
            "Retro_Stock_Coverage",
            "Retro_StepCount",
            "Retro_NotInStockCount",
        ),
        logic="any",
        filters=(
            FilterRule("Retro_RecordMatched", "==", False, missing="warn_skip"),
            FilterRule("Retro_IsSolved", "==", False, missing="warn_skip"),
            FilterRule("Retro_Stock_Coverage", "<", 0.25, missing="warn_skip"),
            FilterRule("Retro_StepCount", ">=", 8, missing="warn_skip"),
            FilterRule("Retro_NotInStockCount", ">=", 3, missing="warn_skip"),
        ),
        sort=(
            SortRule("Retro_NotInStockCount", ascending=False),
            SortRule("Retro_StepCount", ascending=False),
            SortRule("Retro_Rerank_Score", ascending=True),
        ),
        limit=0,
        output_columns=(
            "Retro_IsSolved",
            "Retro_Rerank_Rank",
            "Retro_Rerank_Score",
            "Retro_Stock_Coverage",
            "Retro_StepCount",
            "Retro_NotInStockCount",
            "Retro_HTMLPage",
        ),
    ),
    "make_next": PicklistDefinition(
        key="make_next",
        name="Make next",
        description=(
            "Balanced make-next queue for compounds that look potent enough, route-supported enough, and clean enough to warrant near-term synthesis attention."
        ),
        queue="make_next",
        next_action="Advance to synthesis planning and nominate for the next make cycle.",
        rationale="Balanced compound that looks ready to move from analysis into active synthesis planning.",
        evidence_columns=(
            "Measured_pIC50",
            "Lead_LLE",
            "Retro_Rerank_Rank",
            "Retro_Stock_Coverage",
            "Tox_Risk",
            "hERG_Risk",
        ),
        filters=(
            FilterRule("Lead_LLE", ">=", 5.0, missing="warn_skip"),
            FilterRule("Tox_Risk", "notin", ["High"], missing="warn_skip"),
            FilterRule("hERG_Risk", "notin", ["High"], missing="warn_skip"),
            FilterRule("Retro_Rerank_Rank", "<=", 20, missing="warn_skip"),
            FilterRule("Retro_Stock_Coverage", ">=", 0.5, missing="warn_skip"),
        ),
        sort=(
            SortRule("Measured_pIC50", ascending=False),
            SortRule("Lead_LLE", ascending=False),
            SortRule("Retro_Rerank_Rank", ascending=True),
        ),
        diversity_columns=("Scaffold_ID", "Series_ID", "Cluster_ID"),
        max_per_diversity_group=2,
        limit=24,
        output_columns=(
            "Measured_pIC50",
            "Lead_LLE",
            "Retro_Rerank_Rank",
            "Retro_Stock_Coverage",
            "Tox_Risk",
            "hERG_Risk",
        ),
    ),
    "test_next": PicklistDefinition(
        key="test_next",
        name="Test next",
        description=(
            "Queue promising compounds that still need discriminating experimental data before they can compete for synthesis-heavy effort."
        ),
        queue="test_next",
        next_action="Prioritize for the next discriminating assay or counterscreen set.",
        rationale="Promising chemistry is still measurement-limited and should be clarified experimentally before deeper expansion.",
        evidence_columns=(
            "Measured_pIC50",
            "Lead_LLE",
            "CNS_MPO",
            "Bioavail_Score",
            "Tox_Risk",
        ),
        filters=(
            FilterRule("Measured_pIC50", "isnull", missing="false"),
            FilterRule("Lead_LLE", ">=", 4.5, missing="warn_skip"),
            FilterRule("Tox_Risk", "notin", ["High"], missing="warn_skip"),
            FilterRule("hERG_Risk", "notin", ["High"], missing="warn_skip"),
        ),
        sort=(
            SortRule("Lead_LLE", ascending=False),
            SortRule("CNS_MPO", ascending=False),
            SortRule("Bioavail_Score", ascending=False),
        ),
        diversity_columns=("Scaffold_ID", "Series_ID", "Cluster_ID"),
        max_per_diversity_group=2,
        limit=24,
        output_columns=("Measured_pIC50", "Lead_LLE", "CNS_MPO", "Bioavail_Score", "Tox_Risk"),
    ),
    "de_risk_next": PicklistDefinition(
        key="de_risk_next",
        name="De-risk next",
        description=(
            "Operational queue for compounds that still look worth learning from, but now need counterscreens or liability-focused experiments."
        ),
        queue="de_risk_next",
        next_action="Run targeted counterscreens or liability-focused assays before more chemistry investment.",
        rationale="Compound carries a meaningful liability signal that should be resolved explicitly rather than ignored.",
        evidence_columns=(
            "Measured_pIC50",
            "Lead_LLE",
            "Tox_Risk",
            "hERG_Risk",
            "CYP_Total_Inhibition",
            "MedChem_Alerts",
        ),
        logic="any",
        filters=(
            FilterRule("Tox_Risk", "in", ["Medium", "High"], missing="warn_skip"),
            FilterRule("hERG_Risk", "in", ["Medium", "High"], missing="warn_skip"),
            FilterRule("CYP_Total_Inhibition", ">=", 2, missing="warn_skip"),
            FilterRule("MedChem_Alerts", ">=", 1, missing="warn_skip"),
            FilterRule("Retro_NotInStockCount", ">=", 2, missing="warn_skip"),
        ),
        sort=(
            SortRule("Lead_LLE", ascending=False),
            SortRule("Measured_pIC50", ascending=False),
            SortRule("Tox_Total_Alerts", ascending=False),
        ),
        limit=40,
        output_columns=(
            "Measured_pIC50",
            "Lead_LLE",
            "Tox_Risk",
            "hERG_Risk",
            "CYP_Total_Inhibition",
            "MedChem_Alerts",
            "Retro_NotInStockCount",
        ),
    ),
    "hold": PicklistDefinition(
        key="hold",
        name="Hold",
        description=(
            "Queue for compounds that still look interesting enough to remember, but are not strong enough for immediate make-next or test-next action."
        ),
        queue="hold",
        next_action="Keep visible, but wait for stronger evidence or a clearer strategic reason before moving it forward.",
        rationale="Worth keeping in campaign memory, but not yet compelling enough for the active front of the queue.",
        evidence_columns=(
            "Lead_LLE",
            "Tox_Risk",
            "Retro_Rerank_Rank",
            "Retro_Stock_Coverage",
            "Measured_pIC50",
        ),
        filters=(
            FilterRule("Lead_LLE", ">=", 4.0, missing="warn_skip"),
            FilterRule("Tox_Risk", "notin", ["High"], missing="warn_skip"),
            FilterRule("Retro_Rerank_Rank", ">", 20, missing="warn_skip"),
        ),
        sort=(
            SortRule("Lead_LLE", ascending=False),
            SortRule("Measured_pIC50", ascending=False),
        ),
        limit=40,
        output_columns=(
            "Measured_pIC50",
            "Lead_LLE",
            "Tox_Risk",
            "Retro_Rerank_Rank",
            "Retro_Stock_Coverage",
        ),
    ),
    "backup_series": PicklistDefinition(
        key="backup_series",
        name="Backup series",
        description=(
            "Maintain a practical backup path by keeping cleaner, route-supported compounds from orthogonal chemistry in view."
        ),
        queue="backup_series",
        next_action="Keep this chemotype active as a backup path and protect enough synthesis bandwidth for it.",
        rationale="Cleaner compound with enough tractability to justify keeping the backup series alive.",
        evidence_columns=(
            "Lead_LLE",
            "Bioavail_Score",
            "CYP_Total_Inhibition",
            "Tox_Risk",
            "Retro_Rerank_Rank",
            "Retro_Stock_Coverage",
        ),
        filters=(
            FilterRule("Lead_LLE", ">=", 4.5, missing="warn_skip"),
            FilterRule("CYP_Total_Inhibition", "<=", 1, missing="warn_skip"),
            FilterRule("Tox_Risk", "notin", ["High"], missing="warn_skip"),
            FilterRule("Retro_Rerank_Rank", "<=", 30, missing="warn_skip"),
        ),
        sort=(
            SortRule("Lead_LLE", ascending=False),
            SortRule("Bioavail_Score", ascending=False),
            SortRule("Retro_Rerank_Rank", ascending=True),
        ),
        diversity_columns=("Series_ID", "Scaffold_ID", "Cluster_ID"),
        max_per_diversity_group=1,
        limit=30,
        output_columns=(
            "Measured_pIC50",
            "Lead_LLE",
            "Bioavail_Score",
            "CYP_Total_Inhibition",
            "Tox_Risk",
            "Retro_Rerank_Rank",
            "Retro_Stock_Coverage",
        ),
    ),
}
