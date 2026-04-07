from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]


def test_campaign_alpha_validation_bundle_builds_end_to_end(tmp_path) -> None:
    pytest.importorskip("rdkit")

    outdir = tmp_path / "campaign_alpha_demo"
    script_path = REPO_ROOT / "examples" / "build_campaign_alpha_demo.py"

    completed = subprocess.run(
        [sys.executable, str(script_path), "-o", str(outdir)],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stdout + "\n" + completed.stderr

    summary = json.loads((outdir / "validation_summary.json").read_text(encoding="utf-8"))
    round_review = json.loads((outdir / "compare" / "round_review.json").read_text(encoding="utf-8"))
    make_next = pd.read_csv(outdir / "picklists" / "make_next.csv")
    learnings = pd.read_csv(outdir / "learnings" / "learning_cards.csv")

    assert (outdir / "portal" / "index.html").exists()
    assert (outdir / "dashboard" / "index.html").exists()
    assert (outdir / "report" / "report.html").exists()
    assert summary["shared_ids"] == 4
    assert summary["group_column"] == "Scaffold_ID"
    assert summary["learning_cards"] >= 5
    assert summary["top_negative_shift"]["endpoint"] == "Measured_Solubility_uM"
    assert summary["low_coverage_endpoints"][0]["endpoint"] == "Measured_Caco2_AtoB"
    assert any(item["group"] == "SCF_C" for item in round_review["expanded_groups"])
    assert set(make_next["Compound_ID"]) >= {"CMPD_A4", "CMPD_C2"}
    assert "compare" in set(learnings["source"])
    assert "picklists" in set(learnings["source"])
