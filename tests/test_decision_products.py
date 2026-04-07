from __future__ import annotations

import pandas as pd

from molprop_toolkit.picklists.builtins import BUILTIN_PICKLISTS
from molprop_toolkit.picklists.engine import run_picklist
from molprop_toolkit.reporting import ReportBuilder


def test_make_next_picklist_adds_decision_columns() -> None:
    df = pd.DataFrame(
        {
            "Compound_ID": ["CMPD_001", "CMPD_002"],
            "Measured_pIC50": [7.1, 6.0],
            "Measured_pIC50_Provenance": ["measured", "measured"],
            "Lead_LLE": [5.8, 4.0],
            "Lead_LLE_Provenance": ["measured", "heuristic"],
            "Tox_Risk": ["Low", "High"],
            "hERG_Risk": ["Low", "Medium"],
            "Retro_Rerank_Rank": [4, 35],
            "Retro_Stock_Coverage": [0.9, 0.1],
        }
    )

    result = run_picklist(
        df,
        BUILTIN_PICKLISTS["make_next"],
        stage="lead-optimization",
    )

    assert list(result.df["Compound_ID"]) == ["CMPD_001"]
    row = result.df.iloc[0]
    assert row["Decision_Queue"] == "make_next"
    assert row["Decision_Confidence"] == "high"
    assert "Advance to synthesis planning" in row["Decision_Next_Action"]
    assert "Lead_LLE=5.8 meets >= 5.0" in row["Decision_Rationale"]
    assert "measured: Measured_pIC50, Lead_LLE" in row["Decision_Evidence_Summary"]
    assert "near-term synthesis execution" in row["Decision_Hypothesis"]
    assert "kinetic/equilibrium solubility panel" in row["Decision_Minimal_Test_Set"]
    assert "potency-normalized analogs" in row["Decision_Minimal_Analog_Set"]


def test_report_builder_includes_decision_queue_preview(tmp_path) -> None:
    input_path = tmp_path / "campaign.csv"
    outdir = tmp_path / "report"

    pd.DataFrame(
        {
            "Compound_ID": ["CMPD_001", "CMPD_002", "CMPD_003"],
            "Measured_pIC50": [7.2, 6.4, 6.8],
            "Measured_pIC50_Provenance": ["measured", "measured", "measured"],
            "Lead_LLE": [5.9, 4.4, 5.1],
            "Lead_LLE_Provenance": ["measured", "measured", "heuristic"],
            "Tox_Risk": ["Low", "Medium", "Low"],
            "hERG_Risk": ["Low", "Medium", "Low"],
            "Retro_Rerank_Rank": [3, 28, 45],
            "Retro_Stock_Coverage": [0.9, 0.6, 0.2],
            "Retro_StepCount": [3, 6, 8],
            "Retro_NotInStockCount": [0, 1, 4],
        }
    ).to_csv(input_path, index=False)

    builder = ReportBuilder(
        csv_path=str(input_path),
        top_n=2,
        stage="lead-optimization",
    )
    artifacts = builder.build(
        outdir=outdir,
        write_markdown=True,
        write_html=True,
        make_plots=False,
        title="Decision product report",
    )

    report_md = artifacts.markdown_path.read_text(encoding="utf-8")
    report_html = artifacts.html_path.read_text(encoding="utf-8")

    assert "## Decision products" in report_md
    assert "Operational queues" in report_md
    assert "Experiment prompts" in report_md
    assert "Queue: Make next" in report_md
    assert "CMPD_001" in report_md
    assert "Decision products" in report_html
    assert "Experiment prompts" in report_html
    assert "Queue: Make next" in report_html


def test_de_risk_picklist_adds_targeted_experiment_prompts() -> None:
    df = pd.DataFrame(
        {
            "Compound_ID": ["CMPD_101"],
            "Measured_pIC50": [6.5],
            "Measured_pIC50_Provenance": ["measured"],
            "Lead_LLE": [4.9],
            "Lead_LLE_Provenance": ["measured"],
            "Tox_Risk": ["Medium"],
            "hERG_Risk": ["Medium"],
            "CYP_Total_Inhibition": [3],
            "LogP": [4.2],
            "MedChem_Alerts": [1],
            "Retro_RecordMatched": [True],
            "Retro_IsSolved": [True],
            "Retro_Rerank_Rank": [8],
            "Retro_Stock_Coverage": [0.8],
        }
    )

    result = run_picklist(
        df,
        BUILTIN_PICKLISTS["de_risk_next"],
        stage="liability-de-risking",
    )

    row = result.df.iloc[0]
    assert "rescueable" in row["Decision_Hypothesis"]
    assert "hERG counterscreen" in row["Decision_Minimal_Test_Set"]
    assert "CYP inhibition panel" in row["Decision_Minimal_Test_Set"]
    assert "orthogonal artifact / alert counterscreens" in row["Decision_Minimal_Test_Set"]
    assert "rescue analogs" in row["Decision_Minimal_Analog_Set"]


def test_make_next_picklist_spreads_across_scaffolds_when_available() -> None:
    df = pd.DataFrame(
        {
            "Compound_ID": ["CMPD_A1", "CMPD_A2", "CMPD_A3", "CMPD_B1", "CMPD_C1"],
            "Scaffold_ID": ["SCF_A", "SCF_A", "SCF_A", "SCF_B", "SCF_C"],
            "Measured_pIC50": [8.2, 7.9, 7.7, 7.2, 7.1],
            "Measured_pIC50_Provenance": ["measured"] * 5,
            "Lead_LLE": [6.1, 5.9, 5.7, 5.6, 5.5],
            "Lead_LLE_Provenance": ["measured"] * 5,
            "Tox_Risk": ["Low"] * 5,
            "hERG_Risk": ["Low"] * 5,
            "Retro_Rerank_Rank": [1, 2, 3, 4, 5],
            "Retro_Stock_Coverage": [0.95, 0.9, 0.85, 0.8, 0.75],
        }
    )

    result = run_picklist(
        df,
        BUILTIN_PICKLISTS["make_next"],
        stage="lead-optimization",
    )

    assert list(result.df["Compound_ID"]) == ["CMPD_A1", "CMPD_A2", "CMPD_B1", "CMPD_C1"]
    assert result.df["Scaffold_ID"].value_counts().max() == 2
