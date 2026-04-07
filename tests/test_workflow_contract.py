from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
EXAMPLE_INPUT = REPO_ROOT / "examples" / "example.smi"


def _runner() -> list[str]:
    console = shutil.which("molscope")
    if console:
        return [console]
    return [sys.executable, "-m", "molprop_toolkit.cli"]


def _run_cli(*args: str) -> subprocess.CompletedProcess[str]:
    completed = subprocess.run(
        [*_runner(), *args],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stdout + "\n" + completed.stderr
    return completed


def test_quickstart_workflow_contract_from_example_data(tmp_path) -> None:
    pytest.importorskip("rdkit")

    results_path = tmp_path / "results.csv"
    report_dir = tmp_path / "report"
    picklists_dir = tmp_path / "picklists"

    calc = _run_cli("calc", str(EXAMPLE_INPUT), "-o", str(results_path))
    assert "Saved:" in calc.stdout
    assert results_path.exists()

    results_df = pd.read_csv(results_path)
    assert len(results_df) == 5
    assert "Compound_ID" in results_df.columns
    assert "Lead_LLE" in results_df.columns

    analyze = _run_cli("analyze", str(results_path), "--list")
    assert "Available categories:" in analyze.stdout
    assert "cns_mpo" in analyze.stdout
    assert "developability" in analyze.stdout

    _run_cli(
        "report",
        str(results_path),
        "--outdir",
        str(report_dir),
        "--no-plots",
        "--title",
        "Workflow Contract Test",
    )
    assert (report_dir / "report.md").exists()
    assert (report_dir / "report.html").exists()
    assert "Workflow Contract Test" in (report_dir / "report.html").read_text(
        encoding="utf-8"
    )

    _run_cli("picklists", str(results_path), "--outdir", str(picklists_dir), "--html")
    assert (picklists_dir / "index.html").exists()
    assert any(picklists_dir.glob("*.csv"))
    html_pages = [path for path in picklists_dir.glob("*.html") if path.name != "index.html"]
    assert html_pages
    index_html = (picklists_dir / "index.html").read_text(encoding="utf-8")
    assert "make_next" in index_html
