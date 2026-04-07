from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from molprop_toolkit.tools.dashboard_cli import main


def _write(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def _write_learnings_bundle(
    root: Path,
    *,
    name: str,
    stage: str,
    cards: list[dict[str, str]],
) -> Path:
    bundle = root / name
    bundle.mkdir(parents=True, exist_ok=True)
    _write(
        bundle / "summary.json",
        json.dumps(
            {
                "title": name,
                "stage": stage,
                "sources": sorted({card["source"] for card in cards}),
                "learning_cards": len(cards),
            },
            indent=2,
        ),
    )
    pd.DataFrame(cards).to_csv(bundle / "learning_cards.csv", index=False)
    _write(bundle / "summary.md", f"# {name}\n")
    _write(bundle / "index.html", f"<html><body>{name}</body></html>")
    return bundle


def test_dashboard_cli_persists_campaign_memory_across_snapshots(tmp_path) -> None:
    round4 = _write_learnings_bundle(
        tmp_path,
        name="round4",
        stage="hit-expansion",
        cards=[
            {
                "source": "compare",
                "learning_type": "improvement",
                "title": "Lead_LLE balance improved",
                "evidence": "Lead_LLE moved up versus baseline.",
                "recommendation": "Keep this endpoint in the next design loop.",
                "priority": "promote",
            },
            {
                "source": "picklists",
                "learning_type": "make_next",
                "title": "CMPD_004 leads the make_next queue",
                "evidence": "Balanced queue leader.",
                "recommendation": "Make two analogs around the same core.",
                "priority": "action",
            },
        ],
    )
    round5 = _write_learnings_bundle(
        tmp_path,
        name="round5",
        stage="lead-optimization",
        cards=[
            {
                "source": "compare",
                "learning_type": "improvement",
                "title": "Lead_LLE balance improved",
                "evidence": "Lead_LLE improved again in round5.",
                "recommendation": "Promote the branch into potency-normalized optimization.",
                "priority": "promote",
            },
            {
                "source": "mmp",
                "learning_type": "transform_gain",
                "title": "Transform A>>B tends to improve Lead_LLE",
                "evidence": "Positive mean delta across matched pairs.",
                "recommendation": "Repeat the transform in one more constant environment.",
                "priority": "follow_up",
            },
        ],
    )

    outdir = tmp_path / "dashboard"
    main(
        [
            "--learnings-dir",
            str(round4),
            "-o",
            str(outdir),
            "--title",
            "Campaign Dashboard",
        ]
    )
    main(
        [
            "--learnings-dir",
            str(round5),
            "-o",
            str(outdir),
            "--title",
            "Campaign Dashboard",
        ]
    )

    summary = json.loads((outdir / "summary.json").read_text(encoding="utf-8"))
    snapshots = pd.read_csv(outdir / "snapshot_registry.csv")
    registry = pd.read_csv(outdir / "knowledge_registry.csv")
    recurring = pd.read_csv(outdir / "recurring_lessons.csv")
    active = pd.read_csv(outdir / "active_learning_cards.csv")
    html = (outdir / "index.html").read_text(encoding="utf-8")
    metadata = json.loads((outdir / "index.metadata.json").read_text(encoding="utf-8"))

    assert summary["snapshots"] == 2
    assert summary["total_learning_cards"] == 4
    assert summary["latest_snapshot"] == "round5"
    assert summary["latest_stage"] == "lead-optimization"
    assert set(snapshots["snapshot_id"]) == {"round4", "round5"}
    assert "Lead_LLE balance improved" in set(registry["title"])
    assert recurring.loc[recurring["title"] == "Lead_LLE balance improved", "snapshot_count"].iat[0] == 2
    assert set(active["snapshot_id"]) == {"round5"}
    assert "Campaign Dashboard" in html
    assert "snapshots/round5/index.html" in html
    assert metadata["tool"] == "molscope-dashboard"
