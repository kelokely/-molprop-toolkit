from __future__ import annotations

import json

import pandas as pd
import pytest

pytest.importorskip("rdkit", reason="molscope mmp smoke test requires RDKit")

from molprop_toolkit.tools.mmp_cli import main as mmp_main


def test_mmp_cli_writes_pair_and_transform_tables(tmp_path) -> None:
    input_path = tmp_path / "results.csv"
    output_dir = tmp_path / "mmp"

    pd.DataFrame(
        {
            "Compound_ID": ["CMPD_001", "CMPD_002", "CMPD_003"],
            "SMILES": ["CCOc1ccccc1", "CCNc1ccccc1", "CCSc1ccccc1"],
            "Scaffold_ID": ["SCF_A", "SCF_A", "SCF_A"],
            "Lead_LLE": [4.6, 4.1, 5.0],
            "Measured_pIC50": [6.8, 6.1, 7.0],
        }
    ).to_csv(input_path, index=False)

    mmp_main(
        [
            str(input_path),
            "-o",
            str(output_dir),
            "--metrics",
            "Lead_LLE,Measured_pIC50",
            "--priority-metric",
            "Lead_LLE",
        ]
    )

    pair_table = pd.read_csv(output_dir / "pair_table.csv")
    transform_summary = pd.read_csv(output_dir / "transform_summary.csv")
    summary = json.loads((output_dir / "summary.json").read_text())
    metadata = json.loads(
        (output_dir / "pair_table.metadata.json").read_text(encoding="utf-8")
    )

    assert not pair_table.empty
    assert not transform_summary.empty
    assert summary["pair_rows"] > 0
    assert summary["transforms"] > 0
    assert summary["priority_metric"] == "Lead_LLE"

    expected_transform = "N([*])[*]>>O([*])[*]"
    assert expected_transform in set(transform_summary["transform"])

    transform_row = transform_summary.loc[
        transform_summary["transform"] == expected_transform
    ].iloc[0]
    assert int(transform_row["pair_count"]) >= 1
    assert round(float(transform_row["Lead_LLE_mean_delta"]), 4) == 0.5

    pair_row = pair_table.loc[pair_table["transform"] == expected_transform].iloc[0]
    assert round(float(pair_row["Lead_LLE_delta"]), 4) == 0.5
    assert pair_row["left_compound_id"] == "CMPD_002"
    assert pair_row["right_compound_id"] == "CMPD_001"

    assert metadata["tool"] == "molscope-mmp"
    assert metadata["artifacts"]["transform_summary"].endswith("transform_summary.csv")
