from __future__ import annotations

from pathlib import Path

import molprop_toolkit.cli as cli
from molprop_toolkit.reporting import ReportBuilder

FIXTURE_RESULTS = Path(__file__).resolve().parent / "data" / "fixture_results.csv"


def test_report_builder_includes_stage_section(tmp_path) -> None:
    outdir = tmp_path / "report"

    builder = ReportBuilder(
        csv_path=str(FIXTURE_RESULTS),
        top_n=2,
        stage="hit-expansion",
    )
    artifacts = builder.build(
        outdir=outdir,
        write_markdown=True,
        write_html=True,
        make_plots=False,
        title="Stage-aware report",
    )

    report_md = artifacts.markdown_path.read_text(encoding="utf-8")
    report_html = artifacts.html_path.read_text(encoding="utf-8")

    assert "## Decision stage: Hit expansion" in report_md
    assert "recommended_picklists" in report_md
    assert "test_next" in report_md
    assert "## Decision products" in report_md
    assert "Operational queues" in report_md
    assert "Decision stage: Hit expansion" in report_html
    assert "Decision products" in report_html


def test_picklists_stage_selects_recommended_subset(tmp_path) -> None:
    outdir = tmp_path / "picklists"

    cli.main(
        [
            "picklists",
            str(FIXTURE_RESULTS),
            "--outdir",
            str(outdir),
            "--stage",
            "lead-optimization",
        ]
    )

    csv_names = sorted(path.name for path in outdir.glob("*.csv"))
    assert csv_names == [
        "de_risk_next.csv",
        "hold.csv",
        "make_next.csv",
    ]

    run_meta = (outdir / "picklists_run.metadata.json").read_text(encoding="utf-8")
    assert '"stage": "lead-optimization"' in run_meta
