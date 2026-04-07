from __future__ import annotations

import json
from pathlib import Path

from molprop_toolkit.tools.portal_cli import main


def _write(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def test_portal_cli_builds_linked_campaign_workspace(tmp_path) -> None:
    report_dir = tmp_path / "report"
    picklists_dir = tmp_path / "picklists"
    compare_dir = tmp_path / "compare"
    learnings_dir = tmp_path / "learnings"
    dashboard_dir = tmp_path / "dashboard"
    sar_dir = tmp_path / "sar"
    mmp_dir = tmp_path / "mmp"
    retro_dir = tmp_path / "retro"
    outdir = tmp_path / "portal"

    _write(report_dir / "report.html", "<html><body>Report</body></html>")
    _write(report_dir / "report.md", "# Report\n")

    _write(picklists_dir / "index.html", "<html><body>Picklists</body></html>")
    _write(
        picklists_dir / "picklists_run.metadata.json",
        json.dumps({"stage": "lead-optimization"}, indent=2),
    )
    _write(picklists_dir / "make_next.csv", "Compound_ID\nCMPD_001\n")

    _write(
        compare_dir / "summary.json",
        json.dumps(
            {
                "reference_name": "round4",
                "candidate_name": "round5",
                "shared_ids": 18,
                "group_column": "Scaffold_ID",
            },
            indent=2,
        ),
    )
    _write(compare_dir / "summary.md", "# Compare Summary\n")
    _write(
        compare_dir / "round_review.json",
        json.dumps(
            {
                "positive_shifts": [{"endpoint": "Lead_LLE", "paired_mean_delta": 0.4}],
                "negative_shifts": [{"endpoint": "Sol_Consensus_LogS", "paired_mean_delta": -0.3}],
                "low_coverage_endpoints": [{"endpoint": "CNS_MPO", "shared_non_null": 1}],
            },
            indent=2,
        ),
    )
    _write(compare_dir / "round_review.md", "# Round Review\n")
    _write(compare_dir / "endpoint_deltas.csv", "endpoint,mean_delta\nLead_LLE,0.4\n")

    _write(
        learnings_dir / "summary.json",
        json.dumps(
            {"learning_cards": 5, "sources": ["compare", "sar", "mmp", "picklists"], "stage": "lead-optimization"},
            indent=2,
        ),
    )
    _write(learnings_dir / "summary.md", "# Campaign Learnings\n")
    _write(learnings_dir / "index.html", "<html><body>Learnings</body></html>")
    _write(
        learnings_dir / "learning_cards.csv",
        "source,title\ncompare,Lead_LLE improved in round5\n",
    )

    _write(
        dashboard_dir / "summary.json",
        json.dumps(
            {
                "snapshots": 2,
                "total_learning_cards": 7,
                "latest_snapshot": "round5",
                "latest_stage": "lead-optimization",
            },
            indent=2,
        ),
    )
    _write(dashboard_dir / "summary.md", "# Campaign Dashboard\n")
    _write(dashboard_dir / "index.html", "<html><body>Dashboard</body></html>")
    _write(
        dashboard_dir / "snapshot_registry.csv",
        "snapshot_id,snapshot_label\nround4,round4\nround5,round5\n",
    )
    _write(
        dashboard_dir / "knowledge_registry.csv",
        "snapshot_id,title\nround5,Lead_LLE improved in round5\n",
    )

    _write(
        sar_dir / "summary.json",
        json.dumps({"groups": 4, "rows": 21, "priority_metric": "Measured_pIC50"}, indent=2),
    )
    _write(sar_dir / "summary.md", "# SAR Summary\n")
    _write(sar_dir / "series_summary.csv", "Scaffold_ID,Count\nS1,8\n")

    _write(
        mmp_dir / "summary.json",
        json.dumps({"pair_rows": 12, "transforms": 5, "priority_metric": "Lead_LLE"}, indent=2),
    )
    _write(mmp_dir / "summary.md", "# MMP Summary\n")
    _write(mmp_dir / "transform_summary.csv", "transform,pair_count\nA>>B,3\n")

    _write(retro_dir / "index.html", "<html><body>Retro</body></html>")
    _write(retro_dir / "routes_summary.csv", "Compound_ID,target\nCMPD_001,target1\nCMPD_002,target2\n")
    _write(retro_dir / "retro_enriched.csv", "Compound_ID,Retro_IsSolved\nCMPD_001,True\n")

    main(
        [
            "--report-dir",
            str(report_dir),
            "--picklists-dir",
            str(picklists_dir),
            "--compare-dir",
            str(compare_dir),
            "--learnings-dir",
            str(learnings_dir),
            "--dashboard-dir",
            str(dashboard_dir),
            "--sar-dir",
            str(sar_dir),
            "--mmp-dir",
            str(mmp_dir),
            "--retro-dir",
            str(retro_dir),
            "-o",
            str(outdir),
            "--title",
            "Campaign Workspace",
        ]
    )

    html = (outdir / "index.html").read_text(encoding="utf-8")
    metadata = json.loads((outdir / "index.metadata.json").read_text(encoding="utf-8"))

    assert "Campaign Workspace" in html
    assert "report/report.html" in html
    assert "picklists/index.html" in html
    assert "compare/summary.md" in html
    assert "compare/round_review.md" in html
    assert "learnings/index.html" in html
    assert "dashboard/index.html" in html
    assert "sar/series_summary.csv" in html
    assert "mmp/transform_summary.csv" in html
    assert "retro/index.html" in html
    assert "round4" in html
    assert "lead-optimization" in html
    assert metadata["tool"] == "molscope-portal"
    assert set(metadata["parameters"]["sections"]) == {
        "report",
        "picklists",
        "compare",
        "learnings",
        "dashboard",
        "sar",
        "mmp",
        "retro",
    }
