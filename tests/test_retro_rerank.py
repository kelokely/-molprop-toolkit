from __future__ import annotations

import gzip
import json
from pathlib import Path

import pandas as pd

import tools.retrosynthesis as retro


def _write_fake_aizynth_output(path: Path) -> None:
    df = pd.DataFrame(
        {
            "target": ["CCO", "CCN", "CCS"],
            "is_solved": [True, True, False],
            "top_score": [0.72, 0.88, 0.12],
            "number_of_routes": [3, 4, 0],
            "number_of_solved_routes": [2, 3, 0],
            "search_time": [12.0, 11.0, 8.0],
            "max_transforms": [4, 5, 5],
            "number_of_steps": [2, 5, 7],
            "number_of_precursors": [2, 4, 5],
            "precursors_in_stock": [2, 1, 0],
            "precursors_not_in_stock": [0, 3, 5],
            "trees": [[], [], []],
        }
    )
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(df.to_json(orient="table"))


def test_retro_reranking_and_enriched_table(tmp_path, monkeypatch) -> None:
    input_path = tmp_path / "results.csv"
    outdir = tmp_path / "retrosynthesis"
    pd.DataFrame(
        {
            "Compound_ID": ["CMPD_001", "CMPD_002", "CMPD_003"],
            "SMILES": ["CCO", "CCN", "CCS"],
            "Lead_LLE": [4.8, 5.1, 3.2],
        }
    ).to_csv(input_path, index=False)

    def fake_run(
        config: str,
        smiles_file: Path,
        output_file: Path,
        nproc: int,
        extra_args: list[str],
    ) -> None:
        assert config == "fake-config.yml"
        assert smiles_file.exists()
        assert nproc == 1
        assert extra_args == []
        _write_fake_aizynth_output(output_file)

    monkeypatch.setattr(retro, "_run_aizynthcli", fake_run)

    retro.main(
        [
            str(input_path),
            "--config",
            "fake-config.yml",
            "--outdir",
            str(outdir),
            "--no-images",
            "--rerank-policy",
            "stock-first",
            "--write-retro-table",
        ]
    )

    summary = pd.read_csv(outdir / "routes_summary.csv")
    enriched = pd.read_csv(outdir / "retro_enriched.csv")
    summary_metadata = json.loads(
        (outdir / "routes_summary.metadata.json").read_text(encoding="utf-8")
    )
    enriched_metadata = json.loads(
        (outdir / "retro_enriched.metadata.json").read_text(encoding="utf-8")
    )
    summary_md = (outdir / "summary.md").read_text(encoding="utf-8")

    assert list(summary["Compound_ID"]) == ["CMPD_001", "CMPD_002", "CMPD_003"]
    assert summary.loc[0, "Retro_Rerank_Policy"] == "stock-first"
    assert int(summary.loc[0, "Retro_Rerank_Rank"]) == 1
    assert float(summary.loc[0, "Retro_Stock_Coverage"]) == 1.0
    assert {"Retro_Rerank_Score", "Retro_NotInStock_Score"}.issubset(summary.columns)

    assert {"Retro_Target", "Retro_Rerank_Score", "Retro_Rerank_Rank"}.issubset(
        enriched.columns
    )
    assert enriched["Retro_RecordMatched"].tolist() == [True, True, True]

    assert summary_metadata["tool"] == "molscope-retro"
    assert summary_metadata["parameters"]["rerank_policy"] == "stock-first"
    assert summary_metadata["artifacts"]["retro_table"] == "retro_enriched.csv"
    assert enriched_metadata["tool"] == "molscope-retro"
    assert enriched_metadata["artifacts"]["routes_summary"] == "routes_summary.csv"
    assert "rerank policy: `stock-first`" in summary_md
