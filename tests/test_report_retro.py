from __future__ import annotations

import pandas as pd

from molprop_toolkit.reporting import ReportBuilder


def test_report_builder_includes_retro_section(tmp_path) -> None:
    input_path = tmp_path / "retro_enriched.csv"
    outdir = tmp_path / "report"

    pd.DataFrame(
        {
            "Compound_ID": ["CMPD_001", "CMPD_002", "CMPD_003"],
            "Retro_RecordMatched": [True, True, False],
            "Retro_IsSolved": [True, True, False],
            "Retro_Rerank_Rank": [1, 2, 3],
            "Retro_Rerank_Score": [0.92, 0.61, 0.10],
            "Retro_Stock_Coverage": [1.0, 0.5, 0.0],
            "Retro_StepCount": [2, 5, 8],
            "Retro_TopScore": [0.88, 0.73, 0.12],
            "Lead_LLE": [5.1, 4.4, 3.2],
        }
    ).to_csv(input_path, index=False)

    builder = ReportBuilder(csv_path=str(input_path), top_n=2)
    artifacts = builder.build(
        outdir=outdir,
        write_markdown=True,
        write_html=True,
        make_plots=False,
        title="Retro-aware report",
    )

    report_md = artifacts.markdown_path.read_text(encoding="utf-8")
    report_html = artifacts.html_path.read_text(encoding="utf-8")

    assert "## Retrosynthesis feasibility" in report_md
    assert "Top retro-prioritized compounds" in report_md
    assert "CMPD_001" in report_md
    assert "Retrosynthesis feasibility" in report_html
    assert "CMPD_001" in report_html
