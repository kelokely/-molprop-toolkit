from __future__ import annotations

import json

import pandas as pd

from molprop_toolkit.tools.compare_cli import main as compare_main


def test_compare_cli_writes_summary_bundle(tmp_path) -> None:
    reference_path = tmp_path / "round4.csv"
    candidate_path = tmp_path / "round5.csv"
    output_dir = tmp_path / "compare"

    reference_df = pd.DataFrame(
        {
            "Compound_ID": ["CMPD_001", "CMPD_002", "CMPD_003"],
            "Scaffold_ID": ["SCF_A", "SCF_A", "SCF_B"],
            "Lead_LLE": [4.2, 4.8, 5.1],
            "Measured_pIC50": [6.1, 6.4, 6.0],
            "Sol_Consensus_LogS": [-3.0, -3.1, -3.4],
            "CNS_MPO": [4.2, None, 4.8],
            "Retro_IsSolved": [True, True, False],
        }
    )
    candidate_df = pd.DataFrame(
        {
            "Compound_ID": ["CMPD_002", "CMPD_003", "CMPD_004"],
            "Scaffold_ID": ["SCF_A", "SCF_B", "SCF_C"],
            "Lead_LLE": [5.4, 5.3, 6.0],
            "Measured_pIC50": [7.0, 6.8, 7.0],
            "Sol_Consensus_LogS": [-3.6, -3.8, -2.9],
            "CNS_MPO": [None, 5.1, 5.5],
            "Retro_IsSolved": [True, False, True],
        }
    )
    reference_df.to_csv(reference_path, index=False)
    candidate_df.to_csv(candidate_path, index=False)

    compare_main([str(reference_path), str(candidate_path), "-o", str(output_dir)])

    summary = json.loads((output_dir / "summary.json").read_text())
    endpoint_deltas = pd.read_csv(output_dir / "endpoint_deltas.csv")
    scaffold_deltas = pd.read_csv(output_dir / "scaffold_id_deltas.csv")
    round_review = json.loads((output_dir / "round_review.json").read_text())
    round_review_md = (output_dir / "round_review.md").read_text(encoding="utf-8")
    metadata = json.loads(
        (output_dir / "endpoint_deltas.metadata.json").read_text(encoding="utf-8")
    )

    assert summary["reference_rows"] == 3
    assert summary["candidate_rows"] == 3
    assert summary["shared_ids"] == 2
    assert summary["candidate_only_ids"] == 1
    assert summary["reference_only_ids"] == 1
    assert summary["group_column"] == "Scaffold_ID"

    lle_row = endpoint_deltas.loc[endpoint_deltas["endpoint"] == "Lead_LLE"].iloc[0]
    assert round(float(lle_row["reference_mean"]), 4) == 4.7
    assert round(float(lle_row["candidate_mean"]), 4) == 5.5667
    assert round(float(lle_row["paired_mean_delta"]), 4) == 0.4

    scf_a = scaffold_deltas.loc[scaffold_deltas["Scaffold_ID"] == "SCF_A"].iloc[0]
    scf_c = scaffold_deltas.loc[scaffold_deltas["Scaffold_ID"] == "SCF_C"].iloc[0]
    assert int(scf_a["reference_count"]) == 2
    assert int(scf_a["candidate_count"]) == 1
    assert int(scf_c["count_delta"]) == 1

    assert round_review["positive_shifts"][0]["endpoint"] == "Measured_pIC50"
    assert round_review["negative_shifts"][0]["endpoint"] == "Sol_Consensus_LogS"
    assert round_review["low_coverage_endpoints"][0]["endpoint"] == "CNS_MPO"
    assert "Largest positive paired shifts" in round_review_md
    assert "Open questions" in round_review_md

    assert metadata["tool"] == "molscope-compare"
    assert metadata["artifacts"]["endpoint_summary"].endswith("endpoint_deltas.csv")
    assert metadata["artifacts"]["round_review_json"].endswith("round_review.json")
