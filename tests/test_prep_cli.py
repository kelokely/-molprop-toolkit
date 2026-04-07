from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

pytest.importorskip("rdkit", reason="molscope prep smoke test requires RDKit")


def test_prep_cli_smoke(tmp_path: Path) -> None:
    fixture = Path(__file__).resolve().parents[1] / "examples" / "example.smi"
    output_path = tmp_path / "prepared.csv"
    summary_path = tmp_path / "prepared-summary.json"

    cmd = [
        sys.executable,
        "-m",
        "molprop_toolkit.tools.prep_cli",
        str(fixture),
        "-o",
        str(output_path),
        "--summary",
        str(summary_path),
        "--preset",
        "triage",
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0, res.stdout + "\n" + res.stderr

    assert output_path.exists()
    assert summary_path.exists()

    df = pd.read_csv(output_path)
    assert not df.empty
    assert "Calc_Canonical_SMILES" in df.columns
    assert "Prep_Changed" in df.columns
    assert "Tautomer_Mode" in df.columns
    assert "Stereo_Mode" in df.columns

    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    assert summary["preset"] == "triage"
    assert int(summary["prepared_count"]) == len(df)
    assert "stereo_status_counts" in summary
