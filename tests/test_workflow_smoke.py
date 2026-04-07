from __future__ import annotations

import importlib
from pathlib import Path

import pytest

import molprop_toolkit.cli as cli

FIXTURE_RESULTS = Path(__file__).resolve().parent / "data" / "fixture_results.csv"


def test_molscope_report_smoke_builds_artifacts(tmp_path) -> None:
    outdir = tmp_path / "report"

    cli.main(
        [
            "report",
            str(FIXTURE_RESULTS),
            "--outdir",
            str(outdir),
            "--no-plots",
            "--title",
            "Workflow Smoke Test",
        ]
    )

    assert (outdir / "report.md").exists()
    assert (outdir / "report.html").exists()


def test_molscope_picklists_smoke_builds_outputs(tmp_path) -> None:
    outdir = tmp_path / "picklists"

    cli.main(
        [
            "picklists",
            str(FIXTURE_RESULTS),
            "--outdir",
            str(outdir),
            "--html",
        ]
    )

    assert (outdir / "index.html").exists()
    assert any(outdir.glob("*.csv"))


def test_legacy_fingerprint_compatibility_module_imports() -> None:
    module = importlib.import_module("molprop_toolkit.fingerprints")
    assert hasattr(module, "FingerprintSpec")
    assert hasattr(module, "compute_similarity_columns")


def test_clustering_analyzer_imports_with_compatibility_layer() -> None:
    pytest.importorskip("rdkit")

    module = importlib.import_module("analyzers.analyzer_clustering")
    assert hasattr(module, "analyze_clustering")
