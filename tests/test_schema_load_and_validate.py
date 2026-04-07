from __future__ import annotations

import subprocess
import sys
from pathlib import Path

from molprop_toolkit.schema import load_schema


def test_load_schema_from_bundled_resource() -> None:
    schema = load_schema(None)
    assert "schema_version" in schema
    assert "columns" in schema


def test_validate_csv_schema_cli_smoke() -> None:
    fixture = Path(__file__).resolve().parent / "data" / "fixture_results.csv"
    assert fixture.exists()

    # Run as module so it works in both source and installed layouts.
    cmd = [
        sys.executable,
        "-m",
        "tools.validate_csv_schema",
        str(fixture),
        "--require-category",
        "cns_mpo",
        "--require-category",
        "qed",
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0, res.stdout + "\n" + res.stderr
    assert "Result: PASS" in res.stdout
