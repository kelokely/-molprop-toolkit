#!/usr/bin/env python3
"""Build a reproducible validation bundle from the synthetic Campaign Alpha pack."""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path
from typing import Sequence

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))


PACK_ROOT = REPO_ROOT / "examples" / "campaign_alpha"


def _run(argv: Sequence[str]) -> None:
    from molprop_toolkit.cli import main as molscope_main

    molscope_main(list(argv))


def _write_summary(output_dir: Path) -> None:
    compare_dir = output_dir / "compare"
    learnings_dir = output_dir / "learnings"
    dashboard_dir = output_dir / "dashboard"

    compare_summary = json.loads((compare_dir / "summary.json").read_text(encoding="utf-8"))
    round_review = json.loads((compare_dir / "round_review.json").read_text(encoding="utf-8"))
    learnings_summary = json.loads((learnings_dir / "summary.json").read_text(encoding="utf-8"))
    dashboard_summary = json.loads((dashboard_dir / "summary.json").read_text(encoding="utf-8"))

    summary = {
        "campaign": "Campaign Alpha",
        "reference_round": compare_summary.get("reference_name"),
        "candidate_round": compare_summary.get("candidate_name"),
        "shared_ids": compare_summary.get("shared_ids"),
        "group_column": compare_summary.get("group_column"),
        "top_positive_shift": (round_review.get("positive_shifts") or [{}])[0],
        "top_negative_shift": (round_review.get("negative_shifts") or [{}])[0],
        "low_coverage_endpoints": round_review.get("low_coverage_endpoints") or [],
        "learning_cards": learnings_summary.get("learning_cards"),
        "dashboard_snapshots": dashboard_summary.get("snapshots"),
    }
    (output_dir / "validation_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    lines = [
        "# Campaign Alpha Validation Summary",
        "",
        f"- shared IDs: {summary['shared_ids']}",
        f"- top positive shift: {summary['top_positive_shift'].get('endpoint', '(none)')}",
        f"- top negative shift: {summary['top_negative_shift'].get('endpoint', '(none)')}",
        f"- low-coverage endpoints: {', '.join(item.get('endpoint', '(missing)') for item in summary['low_coverage_endpoints']) or '(none)'}",
        f"- learning cards: {summary['learning_cards']}",
        f"- dashboard snapshots: {summary['dashboard_snapshots']}",
        "",
        "This synthetic pack is meant to validate workflow coherence, comparative interpretation, and artifact handoff rather than claim broad scientific calibration.",
        "",
    ]
    (output_dir / "validation_summary.md").write_text(
        "\n".join(lines),
        encoding="utf-8",
    )


def build_validation_bundle(output_dir: Path) -> None:
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    round4_results = PACK_ROOT / "round4_results.csv"
    round5_results = PACK_ROOT / "round5_results.csv"
    round4_assays = PACK_ROOT / "round4_assays.csv"
    round5_assays = PACK_ROOT / "round5_assays.csv"

    integrated_round4 = output_dir / "round4_integrated.csv"
    integrated_round5 = output_dir / "round5_integrated.csv"

    _run(
        [
            "integrate",
            str(round4_results),
            str(round4_assays),
            "-o",
            str(integrated_round4),
            "--potency-column",
            "pIC50",
        ]
    )
    _run(
        [
            "integrate",
            str(round5_results),
            str(round5_assays),
            "-o",
            str(integrated_round5),
            "--potency-column",
            "pIC50",
        ]
    )

    _run(
        [
            "compare",
            str(integrated_round4),
            str(integrated_round5),
            "-o",
            str(output_dir / "compare"),
            "--reference-name",
            "round4",
            "--candidate-name",
            "round5",
        ]
    )
    _run(
        [
            "report",
            str(integrated_round5),
            "--outdir",
            str(output_dir / "report"),
            "--no-plots",
            "--stage",
            "lead-optimization",
            "--title",
            "Campaign Alpha Round 5",
        ]
    )
    _run(
        [
            "picklists",
            str(integrated_round5),
            "--outdir",
            str(output_dir / "picklists"),
            "--html",
            "--stage",
            "lead-optimization",
        ]
    )
    _run(
        [
            "sar",
            str(integrated_round5),
            "-o",
            str(output_dir / "sar"),
            "--top-compounds-per-group",
            "3",
        ]
    )
    _run(
        [
            "mmp",
            str(integrated_round5),
            "-o",
            str(output_dir / "mmp"),
            "--max-pairs-per-constant",
            "8",
        ]
    )
    _run(
        [
            "learnings",
            "--compare-dir",
            str(output_dir / "compare"),
            "--sar-dir",
            str(output_dir / "sar"),
            "--mmp-dir",
            str(output_dir / "mmp"),
            "--picklists-dir",
            str(output_dir / "picklists"),
            "-o",
            str(output_dir / "learnings"),
            "--title",
            "Campaign Alpha Learnings",
        ]
    )
    _run(
        [
            "dashboard",
            "--learnings-dir",
            str(output_dir / "learnings"),
            "-o",
            str(output_dir / "dashboard"),
            "--title",
            "Campaign Alpha Dashboard",
        ]
    )
    _run(
        [
            "portal",
            "--report-dir",
            str(output_dir / "report"),
            "--picklists-dir",
            str(output_dir / "picklists"),
            "--compare-dir",
            str(output_dir / "compare"),
            "--learnings-dir",
            str(output_dir / "learnings"),
            "--dashboard-dir",
            str(output_dir / "dashboard"),
            "--sar-dir",
            str(output_dir / "sar"),
            "--mmp-dir",
            str(output_dir / "mmp"),
            "-o",
            str(output_dir / "portal"),
            "--title",
            "Campaign Alpha Validation Portal",
        ]
    )

    _write_summary(output_dir)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build the synthetic Campaign Alpha validation bundle."
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        default=str(REPO_ROOT / "validation" / "campaign_alpha_demo"),
        help="Output directory for the generated validation bundle.",
    )
    args = parser.parse_args()

    output_dir = Path(args.output_dir).resolve()
    build_validation_bundle(output_dir)
    print(f"[ok] wrote validation bundle to {output_dir}")


if __name__ == "__main__":
    main()
