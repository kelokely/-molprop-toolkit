from __future__ import annotations

import json

import pandas as pd

from molprop_toolkit.tools.integrate_cli import main as integrate_main


def test_integrate_cli_merges_assays_and_recomputes_lead_metrics(tmp_path) -> None:
    results_path = tmp_path / "results.csv"
    assays_path = tmp_path / "assays.csv"
    output_path = tmp_path / "integrated.csv"

    results_df = pd.DataFrame(
        {
            "Compound_ID": ["CMPD_001", "CMPD_002", "CMPD_003"],
            "SMILES": ["CCO", "c1ccccc1", "CCN"],
            "MolWt": [46.07, 78.11, 45.08],
            "LogP": [-0.1, 2.1, 1.0],
            "TPSA": [20.2, 0.0, 12.0],
            "HeavyAtoms": [3, 6, 3],
            "Lead_LE": [0.22, 0.25, 0.21],
            "Lead_LLE": [1.2, 1.5, 1.1],
        }
    )
    assay_df = pd.DataFrame(
        {
            "Compound_ID": ["CMPD_001", "CMPD_003"],
            "pIC50": [6.2, 7.1],
            "Solubility_uM": [18.5, 4.2],
        }
    )
    results_df.to_csv(results_path, index=False)
    assay_df.to_csv(assays_path, index=False)

    integrate_main(
        [
            str(results_path),
            str(assays_path),
            "-o",
            str(output_path),
            "--potency-column",
            "pIC50",
        ]
    )

    integrated = pd.read_csv(output_path)
    metadata = json.loads(output_path.with_name("integrated.metadata.json").read_text())

    assert output_path.exists()
    assert list(integrated["Measured_Record_Matched"]) == [True, False, True]
    assert "Measured_pIC50" in integrated.columns
    assert "Measured_Solubility_uM" in integrated.columns
    assert "Measured_pIC50_Provenance" in integrated.columns
    assert "Lead_LLE_Provenance" in integrated.columns
    assert "Heuristic_Lead_LLE" in integrated.columns

    cmpd1 = integrated.loc[integrated["Compound_ID"] == "CMPD_001"].iloc[0]
    cmpd2 = integrated.loc[integrated["Compound_ID"] == "CMPD_002"].iloc[0]

    assert round(float(cmpd1["Measured_pIC50"]), 3) == 6.2
    assert round(float(cmpd1["Lead_LLE"]), 3) == 6.3
    assert cmpd1["Lead_LLE_Provenance"] == "measured"
    assert round(float(cmpd2["Lead_LLE"]), 3) == 1.5
    assert cmpd2["Lead_LLE_Provenance"] == "heuristic"
    assert metadata["tool"] == "molscope-integrate"
    assert metadata["artifacts"]["matched_rows"] == 2
