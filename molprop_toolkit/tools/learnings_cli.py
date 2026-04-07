#!/usr/bin/env python3
"""Build a campaign-learnings bundle from existing MolScope artifacts."""

from __future__ import annotations

import argparse
import json
from datetime import datetime
from pathlib import Path
from typing import Any, Optional, Sequence

import pandas as pd

from molprop_toolkit.core.metadata import write_run_metadata


def _default_output_dir() -> str:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return str(Path("learnings") / f"campaign_learnings_{timestamp}")


def _read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def _read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path)


def _as_text(value: Any, default: str = "") -> str:
    if value is None or pd.isna(value):
        return default
    return str(value)


def _html_escape(text: Any) -> str:
    return (
        str(text)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def _card(
    *,
    source: str,
    learning_type: str,
    title: str,
    evidence: str,
    recommendation: str,
    priority: str,
) -> dict[str, str]:
    return {
        "source": source,
        "learning_type": learning_type,
        "title": title,
        "evidence": evidence,
        "recommendation": recommendation,
        "priority": priority,
    }


def _compare_cards(path: Path) -> list[dict[str, str]]:
    summary = _read_json(path / "summary.json")
    endpoint = _read_csv(path / "endpoint_deltas.csv")
    group_column = _as_text(summary.get("group_column"))
    group_df = (
        _read_csv(path / f"{group_column.lower()}_deltas.csv")
        if group_column
        else pd.DataFrame()
    )

    reference_name = _as_text(summary.get("reference_name"), "reference")
    candidate_name = _as_text(summary.get("candidate_name"), "candidate")
    cards: list[dict[str, str]] = []

    if not endpoint.empty:
        delta_column = (
            "paired_mean_delta" if "paired_mean_delta" in endpoint.columns else "mean_delta"
        )
        ranked = endpoint.copy()
        ranked[delta_column] = pd.to_numeric(ranked[delta_column], errors="coerce")
        positive = ranked.loc[ranked[delta_column] > 0].sort_values(
            delta_column, ascending=False
        )
        negative = ranked.loc[ranked[delta_column] < 0].sort_values(
            delta_column, ascending=True
        )

        if not positive.empty:
            row = positive.iloc[0]
            cards.append(
                _card(
                    source="compare",
                    learning_type="improvement",
                    title=f"{row['endpoint']} improved in {candidate_name}",
                    evidence=(
                        f"Paired mean delta {row[delta_column]:.4g} versus {reference_name}"
                    ),
                    recommendation=(
                        "Carry this endpoint into the next round review and check whether the gain is scaffold-specific or campaign-wide."
                    ),
                    priority="promote",
                )
            )

        if not negative.empty:
            row = negative.iloc[0]
            cards.append(
                _card(
                    source="compare",
                    learning_type="regression",
                    title=f"{row['endpoint']} regressed in {candidate_name}",
                    evidence=(
                        f"Paired mean delta {row[delta_column]:.4g} versus {reference_name}"
                    ),
                    recommendation=(
                        "Treat this as a caution signal and inspect the related scaffold or transform views before repeating the same design move."
                    ),
                    priority="caution",
                )
            )

    if not group_df.empty and group_column in group_df.columns:
        growth = group_df.copy()
        growth["count_delta"] = pd.to_numeric(growth["count_delta"], errors="coerce")
        positive = growth.loc[growth["count_delta"] > 0].sort_values(
            "count_delta", ascending=False
        )
        if not positive.empty:
            row = positive.iloc[0]
            cards.append(
                _card(
                    source="compare",
                    learning_type="series_shift",
                    title=f"{row[group_column]} gained share in {candidate_name}",
                    evidence=(
                        f"{group_column} count changed {int(row['reference_count'])} -> {int(row['candidate_count'])}"
                    ),
                    recommendation=(
                        "Check whether the growing chemistry branch is also carrying the endpoint improvements before giving it more synthesis bandwidth."
                    ),
                    priority="follow_up",
                )
            )

    return cards


def _sar_cards(path: Path) -> list[dict[str, str]]:
    summary = _read_json(path / "summary.json")
    series = _read_csv(path / "series_summary.csv")
    if series.empty:
        return []

    group_column = _as_text(summary.get("group_column"), series.columns[0])
    priority_metric = _as_text(summary.get("priority_metric"), "priority_metric")
    top_value_column = f"{priority_metric}_Top_Value"
    top_compound_column = f"{priority_metric}_Top_Compound"

    ranked = series.copy()
    if top_value_column in ranked.columns:
        ranked[top_value_column] = pd.to_numeric(
            ranked[top_value_column], errors="coerce"
        )
        ranked = ranked.sort_values(top_value_column, ascending=False, na_position="last")
    top = ranked.iloc[0]

    return [
        _card(
            source="sar",
            learning_type="series_lead",
            title=f"{top[group_column]} carries the strongest {priority_metric} signal",
            evidence=(
                f"Top {priority_metric} {_as_text(top.get(top_value_column), '(missing)')} "
                f"with {_as_text(top.get('Count'), '0')} compounds; lead {_as_text(top.get(top_compound_column), '(unknown)')}"
            ),
            recommendation=(
                "Use this series as the reference branch for the next analog round, but keep the analog set narrow enough to test one substituent axis at a time."
            ),
            priority="promote",
        )
    ]


def _mmp_cards(path: Path) -> list[dict[str, str]]:
    summary = _read_json(path / "summary.json")
    transforms = _read_csv(path / "transform_summary.csv")
    if transforms.empty:
        return []

    priority_metric = _as_text(summary.get("priority_metric"), "priority_metric")
    mean_col = f"{priority_metric}_mean_delta"
    if mean_col not in transforms.columns:
        return []

    ranked = transforms.copy()
    ranked[mean_col] = pd.to_numeric(ranked[mean_col], errors="coerce")
    positive = ranked.loc[ranked[mean_col] > 0].sort_values(mean_col, ascending=False)
    negative = ranked.loc[ranked[mean_col] < 0].sort_values(mean_col, ascending=True)

    cards: list[dict[str, str]] = []
    if not positive.empty:
        row = positive.iloc[0]
        cards.append(
            _card(
                source="mmp",
                learning_type="transform_gain",
                title=f"Transform {row['transform']} tends to improve {priority_metric}",
                evidence=(
                    f"Mean delta {row[mean_col]:.4g} across {int(row['pair_count'])} matched pairs"
                ),
                recommendation=(
                    "Repeat this transform in one more constant environment before treating it as a reusable SAR rule."
                ),
                priority="promote",
            )
        )
    if not negative.empty:
        row = negative.iloc[0]
        cards.append(
            _card(
                source="mmp",
                learning_type="transform_risk",
                title=f"Transform {row['transform']} tends to hurt {priority_metric}",
                evidence=(
                    f"Mean delta {row[mean_col]:.4g} across {int(row['pair_count'])} matched pairs"
                ),
                recommendation=(
                    "Avoid scaling this transform blindly; use it only if another project objective clearly outweighs the regression."
                ),
                priority="caution",
            )
        )
    return cards


def _picklists_cards(path: Path) -> tuple[list[dict[str, str]], str | None]:
    metadata = _read_json(path / "picklists_run.metadata.json")
    stage = _as_text(metadata.get("stage")) or None
    cards: list[dict[str, str]] = []

    for queue_key in ("make_next", "test_next", "de_risk_next", "backup_series", "hold"):
        queue_path = path / f"{queue_key}.csv"
        queue_df = _read_csv(queue_path)
        if queue_df.empty:
            continue

        top = queue_df.iloc[0]
        compound_id = _as_text(top.get("Compound_ID"), "(unknown)")
        queue_name = _as_text(top.get("Decision_Queue"), queue_key)
        cards.append(
            _card(
                source="picklists",
                learning_type=queue_name,
                title=f"{compound_id} leads the {queue_name} queue",
                evidence=_as_text(
                    top.get("Decision_Rationale"),
                    _as_text(
                        top.get("Decision_Evidence_Summary"),
                        "Decision evidence not available.",
                    ),
                ),
                recommendation=(
                    f"Minimal test set: {_as_text(top.get('Decision_Minimal_Test_Set'), '(not set)')} "
                    f"| Minimal analog set: {_as_text(top.get('Decision_Minimal_Analog_Set'), '(not set)')}"
                ),
                priority="action",
            )
        )

    return cards[:5], stage


def _write_markdown(
    path: Path,
    *,
    title: str,
    cards: list[dict[str, str]],
    sources: list[str],
    stage: str | None,
) -> None:
    lines = [
        f"# {title}",
        "",
        "## Overview",
        f"- learning cards: {len(cards)}",
        f"- sources: {', '.join(sources) if sources else '(none)'}",
        f"- stage context: {stage or '(none)'}",
        "",
    ]

    grouped: dict[str, list[dict[str, str]]] = {}
    for card in cards:
        grouped.setdefault(card["source"], []).append(card)

    for source in ("compare", "sar", "mmp", "picklists"):
        source_cards = grouped.get(source, [])
        if not source_cards:
            continue
        lines.extend([f"## {source.title()} learnings", ""])
        for card in source_cards:
            lines.append(f"### {card['title']}")
            lines.append("")
            lines.append(f"- type: {card['learning_type']}")
            lines.append(f"- priority: {card['priority']}")
            lines.append(f"- evidence: {card['evidence']}")
            lines.append(f"- recommendation: {card['recommendation']}")
            lines.append("")

    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def _write_html(
    path: Path,
    *,
    title: str,
    cards: list[dict[str, str]],
    stage: str | None,
) -> None:
    card_html = []
    for card in cards:
        card_html.append(
            "<section class='card'>"
            f"<div class='eyebrow'>{_html_escape(card['source'])} • {_html_escape(card['priority'])}</div>"
            f"<h2>{_html_escape(card['title'])}</h2>"
            f"<p><strong>Type:</strong> {_html_escape(card['learning_type'])}</p>"
            f"<p><strong>Evidence:</strong> {_html_escape(card['evidence'])}</p>"
            f"<p><strong>Recommendation:</strong> {_html_escape(card['recommendation'])}</p>"
            "</section>"
        )

    html = f"""<!doctype html>
<html lang='en'>
<head>
  <meta charset='utf-8'>
  <meta name='viewport' content='width=device-width, initial-scale=1'>
  <title>{_html_escape(title)}</title>
  <style>
    :root {{
      --bg: #f5f7f8;
      --surface: #ffffff;
      --text: #17212b;
      --muted: #5f6b76;
      --border: rgba(15, 23, 42, 0.10);
      --brand: #0f766e;
      --shadow: 0 20px 50px rgba(15, 23, 42, 0.08);
    }}
    body {{
      margin: 0;
      font-family: system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      background: linear-gradient(180deg, #f8fafb 0%, #f3f6f7 100%);
      color: var(--text);
    }}
    .shell {{
      max-width: 1120px;
      margin: 0 auto;
      padding: 40px 20px 56px;
    }}
    .hero {{
      background: radial-gradient(circle at top right, rgba(15, 118, 110, 0.10), transparent 38%), var(--surface);
      border: 1px solid var(--border);
      border-radius: 28px;
      padding: 28px 32px;
      box-shadow: var(--shadow);
      margin-bottom: 24px;
    }}
    .hero h1 {{
      margin: 0 0 8px;
      font-size: 36px;
      line-height: 1.1;
    }}
    .hero p {{
      margin: 6px 0 0;
      max-width: 76ch;
      color: var(--muted);
    }}
    .grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
      gap: 18px;
    }}
    .card {{
      background: var(--surface);
      border: 1px solid var(--border);
      border-radius: 22px;
      padding: 22px;
      box-shadow: var(--shadow);
    }}
    .card h2 {{
      margin: 0 0 10px;
      font-size: 22px;
    }}
    .eyebrow {{
      display: inline-block;
      margin-bottom: 10px;
      padding: 4px 10px;
      border-radius: 999px;
      background: rgba(15, 118, 110, 0.10);
      color: var(--brand);
      font-size: 12px;
      font-weight: 700;
      letter-spacing: 0.04em;
      text-transform: uppercase;
    }}
  </style>
</head>
<body>
  <main class='shell'>
    <section class='hero'>
      <div class='eyebrow'>Campaign learnings</div>
      <h1>{_html_escape(title)}</h1>
      <p>
        This bundle distills compare, SAR, MMP, and decision outputs into reusable project lessons instead of leaving the
        learning trapped inside separate artifacts.
      </p>
      <p>Stage context: {_html_escape(stage or '(none)')}</p>
    </section>
    <div class='grid'>
      {''.join(card_html)}
    </div>
  </main>
</body>
</html>
"""
    path.write_text(html, encoding="utf-8")


def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(
        prog="molscope learnings",
        description="Synthesize compare, SAR, MMP, and decision outputs into a reusable campaign-learnings bundle.",
    )
    ap.add_argument(
        "-o",
        "--output-dir",
        default=None,
        help="Output directory for the learnings bundle (default: learnings/campaign_learnings_<timestamp>)",
    )
    ap.add_argument("--title", default="MolScope Campaign Learnings", help="Bundle title")
    ap.add_argument("--compare-dir", default=None, help="Directory containing molscope compare outputs")
    ap.add_argument("--sar-dir", default=None, help="Directory containing molscope sar outputs")
    ap.add_argument("--mmp-dir", default=None, help="Directory containing molscope mmp outputs")
    ap.add_argument("--picklists-dir", default=None, help="Directory containing molscope picklists outputs")
    args = ap.parse_args(list(argv) if argv is not None else None)

    output_dir = Path(args.output_dir or _default_output_dir())
    output_dir.mkdir(parents=True, exist_ok=True)

    cards: list[dict[str, str]] = []
    sources: list[str] = []
    stage: str | None = None

    if args.compare_dir:
        compare_dir = Path(args.compare_dir)
        if not compare_dir.exists():
            raise SystemExit(f"molscope learnings: compare directory not found: {compare_dir}")
        sources.append("compare")
        cards.extend(_compare_cards(compare_dir))

    if args.sar_dir:
        sar_dir = Path(args.sar_dir)
        if not sar_dir.exists():
            raise SystemExit(f"molscope learnings: SAR directory not found: {sar_dir}")
        sources.append("sar")
        cards.extend(_sar_cards(sar_dir))

    if args.mmp_dir:
        mmp_dir = Path(args.mmp_dir)
        if not mmp_dir.exists():
            raise SystemExit(f"molscope learnings: MMP directory not found: {mmp_dir}")
        sources.append("mmp")
        cards.extend(_mmp_cards(mmp_dir))

    if args.picklists_dir:
        picklists_dir = Path(args.picklists_dir)
        if not picklists_dir.exists():
            raise SystemExit(f"molscope learnings: picklists directory not found: {picklists_dir}")
        sources.append("picklists")
        picklist_cards, detected_stage = _picklists_cards(picklists_dir)
        cards.extend(picklist_cards)
        stage = detected_stage or stage

    if not sources:
        raise SystemExit("molscope learnings: provide at least one artifact directory")

    cards_df = pd.DataFrame(cards)
    cards_path = output_dir / "learning_cards.csv"
    summary_json_path = output_dir / "summary.json"
    summary_md_path = output_dir / "summary.md"
    index_html_path = output_dir / "index.html"

    cards_df.to_csv(cards_path, index=False)

    summary = {
        "title": args.title,
        "sources": sources,
        "stage": stage,
        "learning_cards": int(len(cards_df)),
        "cards_by_source": cards_df["source"].value_counts().to_dict() if not cards_df.empty else {},
    }
    summary_json_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    _write_markdown(
        summary_md_path,
        title=args.title,
        cards=cards,
        sources=sources,
        stage=stage,
    )
    _write_html(index_html_path, title=args.title, cards=cards, stage=stage)

    write_run_metadata(
        tool="molscope-learnings",
        output_table_path=cards_path,
        parameters={
            "title": args.title,
            "sources": sources,
            "stage": stage,
        },
        artifacts={
            "summary_json": str(summary_json_path.resolve()),
            "summary_markdown": str(summary_md_path.resolve()),
            "index_html": str(index_html_path.resolve()),
            "compare_dir": args.compare_dir,
            "sar_dir": args.sar_dir,
            "mmp_dir": args.mmp_dir,
            "picklists_dir": args.picklists_dir,
        },
        notes="Campaign learnings bundle synthesized from compare, SAR, MMP, and decision outputs.",
    )

    print(f"[ok] wrote learnings bundle to {output_dir}")


if __name__ == "__main__":
    main()
