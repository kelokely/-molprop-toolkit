from __future__ import annotations

import json

import pandas as pd

from molprop_toolkit.tools.sar_cli import main as sar_main


def test_sar_cli_writes_scaffold_bundle(tmp_path) -> None:
    input_path = tmp_path / "results_series.csv"
    output_dir = tmp_path / "sar"

    df = pd.DataFrame(
        {
            "Compound_ID": ["CMPD_001", "CMPD_002", "CMPD_003", "CMPD_004"],
            "Scaffold_ID": ["SCF_A", "SCF_A", "SCF_B", "SCF_B"],
            "Scaffold_SMILES": ["c1ccccc1", "c1ccccc1", "C1CCCCC1", "C1CCCCC1"],
            "Measured_pIC50": [6.4, 7.1, 5.8, 6.2],
            "Lead_LLE": [4.0, 5.1, 3.8, 4.4],
            "CNS_MPO": [4.2, 4.8, 3.9, 4.1],
        }
    )
    df.to_csv(input_path, index=False)

    sar_main([str(input_path), "-o", str(output_dir)])

    series_summary = pd.read_csv(output_dir / "series_summary.csv")
    top_compounds = pd.read_csv(output_dir / "top_compounds.csv")
    metric_summary = pd.read_csv(output_dir / "metric_summary.csv")
    summary = json.loads((output_dir / "summary.json").read_text())
    metadata = json.loads(
        (output_dir / "series_summary.metadata.json").read_text(encoding="utf-8")
    )

    assert summary["group_column"] == "Scaffold_ID"
    assert summary["priority_metric"] == "Measured_pIC50"
    assert summary["groups"] == 2

    scf_a = series_summary.loc[series_summary["Scaffold_ID"] == "SCF_A"].iloc[0]
    assert int(scf_a["Count"]) == 2
    assert round(float(scf_a["Measured_pIC50_mean"]), 4) == 6.75
    assert scf_a["Measured_pIC50_Top_Compound"] == "CMPD_002"

    top_a = top_compounds.loc[top_compounds["Scaffold_ID"] == "SCF_A"].iloc[0]
    assert top_a["Compound_ID"] == "CMPD_002"
    assert round(float(top_a["Measured_pIC50"]), 4) == 7.1

    assert set(metric_summary["metric"]) >= {"Measured_pIC50", "Lead_LLE", "CNS_MPO"}
    assert metadata["tool"] == "molscope-sar"
    assert metadata["artifacts"]["top_compounds"].endswith("top_compounds.csv")
