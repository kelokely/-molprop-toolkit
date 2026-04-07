from __future__ import annotations

import json

import pandas as pd
import pytest

pytest.importorskip("rdkit", reason="molscope search tests require RDKit")

from molprop_toolkit.tools.search_cli import main as search_main


def test_search_cli_exact_smiles(tmp_path) -> None:
    input_path = tmp_path / "results.csv"
    output_dir = tmp_path / "search"

    pd.DataFrame(
        {
            "Compound_ID": ["CMPD_001", "CMPD_002", "CMPD_003"],
            "SMILES": ["CCOc1ccccc1", "CCNc1ccccc1", "CCSc1ccccc1"],
            "Lead_LLE": [4.6, 4.1, 5.0],
        }
    ).to_csv(input_path, index=False)

    search_main(
        [
            str(input_path),
            "-o",
            str(output_dir),
            "--exact-smiles",
            "CCOc1ccccc1",
            "--query-name",
            "phenoxy",
        ]
    )

    results = pd.read_csv(output_dir / "search_results.csv")
    summary = json.loads((output_dir / "summary.json").read_text())
    metadata = json.loads(
        (output_dir / "search_results.metadata.json").read_text(encoding="utf-8")
    )

    assert list(results["Compound_ID"]) == ["CMPD_001"]
    assert summary["query_type"] == "exact"
    assert summary["matched_rows"] == 1
    assert summary["query_name"] == "phenoxy"
    assert metadata["tool"] == "molscope-search"


def test_search_cli_smarts_and_motif_panel(tmp_path) -> None:
    input_path = tmp_path / "results.csv"
    smarts_outdir = tmp_path / "smarts"
    motif_outdir = tmp_path / "motifs"
    motif_file = tmp_path / "motifs.tsv"

    pd.DataFrame(
        {
            "Compound_ID": ["CMPD_001", "CMPD_002", "CMPD_003"],
            "SMILES": ["CCOc1ccccc1", "CCNc1ccccc1", "CCSc1ccccc1"],
            "Scaffold_ID": ["SCF_A", "SCF_A", "SCF_A"],
        }
    ).to_csv(input_path, index=False)
    motif_file.write_text(
        "name\tsmarts\n"
        "aniline_like\t[#7]-c1ccccc1\n"
        "phenoxy_like\t[#8]-c1ccccc1\n",
        encoding="utf-8",
    )

    search_main(
        [
            str(input_path),
            "-o",
            str(smarts_outdir),
            "--smarts",
            "[#7]-c1ccccc1",
            "--query-name",
            "aniline_like",
        ]
    )
    smarts_results = pd.read_csv(smarts_outdir / "search_results.csv")
    smarts_summary = json.loads((smarts_outdir / "summary.json").read_text())

    assert list(smarts_results["Compound_ID"]) == ["CMPD_002"]
    assert int(smarts_results.iloc[0]["Search_Match_Count"]) >= 1
    assert smarts_summary["query_type"] == "smarts"
    assert smarts_summary["matched_rows"] == 1

    search_main(
        [
            str(input_path),
            "-o",
            str(motif_outdir),
            "--motif-file",
            str(motif_file),
        ]
    )
    motif_results = pd.read_csv(motif_outdir / "search_results.csv")
    motif_hits = pd.read_csv(motif_outdir / "motif_hits.csv")
    motif_summary = pd.read_csv(motif_outdir / "motif_summary.csv")
    motif_json = json.loads((motif_outdir / "summary.json").read_text())

    assert set(motif_results["Compound_ID"]) == {"CMPD_001", "CMPD_002"}
    assert set(motif_hits["motif_name"]) == {"aniline_like", "phenoxy_like"}
    assert set(motif_summary["motif_name"]) == {"aniline_like", "phenoxy_like"}
    assert motif_json["query_type"] == "motif_panel"
    assert motif_json["matched_rows"] == 2
