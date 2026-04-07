from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from molprop_toolkit.tools.learnings_cli import main


def _write(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def test_learnings_cli_builds_campaign_learning_cards(tmp_path) -> None:
    compare_dir = tmp_path / "compare"
    sar_dir = tmp_path / "sar"
    mmp_dir = tmp_path / "mmp"
    picklists_dir = tmp_path / "picklists"
    outdir = tmp_path / "learnings"

    _write(
        compare_dir / "summary.json",
        json.dumps(
            {
                "reference_name": "round4",
                "candidate_name": "round5",
                "group_column": "Scaffold_ID",
            },
            indent=2,
        ),
    )
    pd.DataFrame(
        [
            {"endpoint": "Lead_LLE", "paired_mean_delta": 0.4},
            {"endpoint": "Tox_Risk_Index", "paired_mean_delta": -0.2},
        ]
    ).to_csv(compare_dir / "endpoint_deltas.csv", index=False)
    pd.DataFrame(
        [{"Scaffold_ID": "SCF_C", "reference_count": 1, "candidate_count": 3, "count_delta": 2}]
    ).to_csv(compare_dir / "scaffold_id_deltas.csv", index=False)

    _write(
        sar_dir / "summary.json",
        json.dumps(
            {"group_column": "Scaffold_ID", "priority_metric": "Measured_pIC50"},
            indent=2,
        ),
    )
    pd.DataFrame(
        [
            {
                "Scaffold_ID": "SCF_A",
                "Count": 4,
                "Measured_pIC50_Top_Value": 7.1,
                "Measured_pIC50_Top_Compound": "CMPD_002",
            }
        ]
    ).to_csv(sar_dir / "series_summary.csv", index=False)

    _write(
        mmp_dir / "summary.json",
        json.dumps({"priority_metric": "Lead_LLE"}, indent=2),
    )
    pd.DataFrame(
        [
            {"transform": "A>>B", "pair_count": 3, "Lead_LLE_mean_delta": 0.5},
            {"transform": "B>>C", "pair_count": 2, "Lead_LLE_mean_delta": -0.4},
        ]
    ).to_csv(mmp_dir / "transform_summary.csv", index=False)

    _write(
        picklists_dir / "picklists_run.metadata.json",
        json.dumps({"stage": "lead-optimization"}, indent=2),
    )
    pd.DataFrame(
        [
            {
                "Compound_ID": "CMPD_010",
                "Decision_Queue": "make_next",
                "Decision_Rationale": "Balanced queue leader.",
                "Decision_Minimal_Test_Set": "kinetic/equilibrium solubility panel",
                "Decision_Minimal_Analog_Set": "Make 2-3 potency-normalized analogs",
            }
        ]
    ).to_csv(picklists_dir / "make_next.csv", index=False)

    main(
        [
            "--compare-dir",
            str(compare_dir),
            "--sar-dir",
            str(sar_dir),
            "--mmp-dir",
            str(mmp_dir),
            "--picklists-dir",
            str(picklists_dir),
            "-o",
            str(outdir),
            "--title",
            "Campaign Learnings",
        ]
    )

    cards = pd.read_csv(outdir / "learning_cards.csv")
    summary = json.loads((outdir / "summary.json").read_text(encoding="utf-8"))
    html = (outdir / "index.html").read_text(encoding="utf-8")
    metadata = json.loads((outdir / "learning_cards.metadata.json").read_text(encoding="utf-8"))

    assert not cards.empty
    assert set(cards["source"]) >= {"compare", "sar", "mmp", "picklists"}
    assert "Lead_LLE improved in round5" in set(cards["title"])
    assert "SCF_A carries the strongest Measured_pIC50 signal" in set(cards["title"])
    assert "Transform A>>B tends to improve Lead_LLE" in set(cards["title"])
    assert "CMPD_010 leads the make_next queue" in set(cards["title"])
    assert summary["stage"] == "lead-optimization"
    assert summary["learning_cards"] >= 4
    assert "Campaign Learnings" in html
    assert "learning trapped inside separate artifacts" in html
    assert metadata["tool"] == "molscope-learnings"

