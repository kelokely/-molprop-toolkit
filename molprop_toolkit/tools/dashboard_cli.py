#!/usr/bin/env python3
"""Build a persistent campaign dashboard from MolScope learnings bundles."""

from __future__ import annotations

import argparse
import json
import re
import shutil
from datetime import datetime
from pathlib import Path
from typing import Any, Optional, Sequence

import pandas as pd

from molprop_toolkit.core.metadata import write_run_metadata


def _default_output_dir() -> str:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return str(Path("dashboard") / f"campaign_dashboard_{timestamp}")


def _read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def _read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path)


def _as_text(value: Any, default: str = "") -> str:
    if value is None:
        return default
    try:
        if pd.isna(value):
            return default
    except TypeError:
        pass
    return str(value)


def _html_escape(text: Any) -> str:
    return (
        str(text)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def _sanitize_snapshot_id(raw: str) -> str:
    cleaned = re.sub(r"[^a-zA-Z0-9]+", "_", raw).strip("_").lower()
    return cleaned or "snapshot"


def _snapshot_label(path: Path, summary: dict[str, Any]) -> str:
    title = _as_text(summary.get("title"))
    if title and title.lower() not in {"campaign learnings", "molscope campaign learnings"}:
        return title
    return path.name


def _snapshot_id(path: Path, summary: dict[str, Any]) -> str:
    preferred = _as_text(summary.get("snapshot_id")) or _snapshot_label(path, summary)
    return _sanitize_snapshot_id(preferred)


def _replace_tree(src: Path, dst: Path) -> None:
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(src, dst)


def _write_csv(df: pd.DataFrame, path: Path) -> None:
    if df.empty:
        path.write_text("", encoding="utf-8")
        return
    df.to_csv(path, index=False)


def _sort_cards(cards: pd.DataFrame) -> pd.DataFrame:
    if cards.empty:
        return cards
    priority_rank = {
        "action": 0,
        "promote": 1,
        "follow_up": 2,
        "caution": 3,
    }
    ranked = cards.copy()
    ranked["_priority_rank"] = ranked["priority"].map(priority_rank).fillna(9)
    ranked = ranked.sort_values(
        ["_priority_rank", "source", "title"],
        ascending=[True, True, True],
        na_position="last",
    )
    return ranked.drop(columns=["_priority_rank"])


def _render_dataframe(
    df: pd.DataFrame,
    *,
    columns: Sequence[str],
    empty_message: str,
    link_columns: Optional[dict[str, str]] = None,
) -> str:
    if df.empty:
        return f"<p class='muted'>{_html_escape(empty_message)}</p>"

    headers = "".join(f"<th>{_html_escape(col)}</th>" for col in columns)
    body_rows = []
    for _, row in df.iterrows():
        cells = []
        for col in columns:
            value = _as_text(row.get(col), "(missing)")
            if link_columns and col in link_columns and value not in {"", "(missing)"}:
                href = _html_escape(_as_text(row.get(link_columns[col]), "#"))
                cells.append(f"<td><a href='{href}'>{_html_escape(value)}</a></td>")
            else:
                cells.append(f"<td>{_html_escape(value)}</td>")
        body_rows.append("<tr>" + "".join(cells) + "</tr>")
    return (
        "<table class='summary-table'><thead><tr>"
        + headers
        + "</tr></thead><tbody>"
        + "".join(body_rows)
        + "</tbody></table>"
    )


def _build_recurring_lessons(cards: pd.DataFrame) -> pd.DataFrame:
    if cards.empty:
        return pd.DataFrame(
            columns=[
                "title",
                "source",
                "learning_type",
                "priority",
                "snapshot_count",
                "total_mentions",
                "latest_snapshot",
                "latest_stage",
                "latest_recommendation",
            ]
        )

    grouped = (
        cards.groupby(["title", "source", "learning_type", "priority"], dropna=False)
        .agg(
            snapshot_count=("snapshot_id", "nunique"),
            total_mentions=("snapshot_id", "size"),
            latest_snapshot=("snapshot_label", "last"),
            latest_stage=("snapshot_stage", "last"),
            latest_recommendation=("recommendation", "last"),
        )
        .reset_index()
    )
    grouped = grouped.loc[grouped["snapshot_count"] > 1]
    if grouped.empty:
        return grouped
    return grouped.sort_values(
        ["snapshot_count", "total_mentions", "priority", "title"],
        ascending=[False, False, True, True],
    )


def _write_markdown(
    output_dir: Path,
    *,
    title: str,
    snapshot_registry: pd.DataFrame,
    active_cards: pd.DataFrame,
    recurring: pd.DataFrame,
    summary: dict[str, Any],
) -> None:
    lines = [
        f"# {title}",
        "",
        f"- Snapshots tracked: {summary['snapshots']}",
        f"- Learning cards tracked: {summary['total_learning_cards']}",
        f"- Latest snapshot: {summary.get('latest_snapshot') or '(none)'}",
        f"- Latest stage: {summary.get('latest_stage') or '(none)'}",
        "",
        "## Recent snapshots",
        "",
    ]

    if snapshot_registry.empty:
        lines.append("- No snapshots tracked yet.")
    else:
        for _, row in snapshot_registry.head(8).iterrows():
            lines.append(
                "- "
                + f"{row['snapshot_label']} ({row['stage']}) — "
                + f"{int(row['learning_cards'])} learning cards "
                + f"[index]({row['snapshot_link']})"
            )

    lines.extend(["", "## Active lessons", ""])
    if active_cards.empty:
        lines.append("- No active lessons available.")
    else:
        for _, row in active_cards.iterrows():
            lines.append(
                "- "
                + f"[{row['priority']}] {row['title']} — {row['recommendation']}"
            )

    lines.extend(["", "## Recurring lessons", ""])
    if recurring.empty:
        lines.append("- No recurring lessons detected yet.")
    else:
        for _, row in recurring.head(10).iterrows():
            lines.append(
                "- "
                + f"{row['title']} ({row['source']}) repeated across {int(row['snapshot_count'])} snapshots"
            )

    (output_dir / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_html(
    output_dir: Path,
    *,
    title: str,
    summary: dict[str, Any],
    snapshot_registry: pd.DataFrame,
    active_cards: pd.DataFrame,
    recurring: pd.DataFrame,
    source_summary: pd.DataFrame,
) -> None:
    cards = [
        ("Snapshots tracked", str(summary["snapshots"])),
        ("Learning cards", str(summary["total_learning_cards"])),
        ("Latest snapshot", summary.get("latest_snapshot") or "(none)"),
        ("Latest stage", summary.get("latest_stage") or "(none)"),
    ]
    hero_cards = "".join(
        "<div class='hero-card'>"
        f"<div class='hero-card-label'>{_html_escape(label)}</div>"
        f"<div class='hero-card-value'>{_html_escape(value)}</div>"
        "</div>"
        for label, value in cards
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
      max-width: 1180px;
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
      margin: 0 0 18px;
      max-width: 80ch;
      color: var(--muted);
    }}
    .hero-grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
      gap: 14px;
    }}
    .hero-card {{
      background: rgba(15, 118, 110, 0.06);
      border: 1px solid rgba(15, 118, 110, 0.10);
      border-radius: 18px;
      padding: 16px;
    }}
    .hero-card-label {{
      color: var(--muted);
      font-size: 12px;
      font-weight: 700;
      letter-spacing: 0.04em;
      text-transform: uppercase;
      margin-bottom: 6px;
    }}
    .hero-card-value {{
      font-size: 22px;
      font-weight: 700;
    }}
    .grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(320px, 1fr));
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
      margin: 0 0 8px;
      font-size: 22px;
    }}
    .muted {{
      color: var(--muted);
    }}
    .summary-table {{
      width: 100%;
      border-collapse: collapse;
      margin-top: 14px;
      font-size: 14px;
    }}
    .summary-table th,
    .summary-table td {{
      text-align: left;
      vertical-align: top;
      padding: 8px 10px 8px 0;
      border-bottom: 1px solid var(--border);
    }}
    .summary-table th {{
      color: var(--muted);
      font-weight: 700;
    }}
    .summary-table a {{
      color: var(--brand);
      text-decoration: none;
      font-weight: 600;
    }}
    .summary-table a:hover {{
      text-decoration: underline;
    }}
  </style>
</head>
<body>
  <main class='shell'>
    <section class='hero'>
      <h1>{_html_escape(title)}</h1>
      <p>
        This dashboard accumulates campaign learnings across rounds so the team can keep a durable registry of what changed,
        what patterns keep recurring, and which lessons should still shape the next design cycle.
      </p>
      <div class='hero-grid'>
        {hero_cards}
      </div>
    </section>
    <div class='grid'>
      <section class='card'>
        <h2>Snapshot history</h2>
        <p class='muted'>Each imported learnings bundle is copied into the dashboard so round-level context stays browsable.</p>
        {_render_dataframe(snapshot_registry.head(12), columns=["snapshot_label", "stage", "learning_cards", "sources", "imported_at"], empty_message="No snapshots tracked yet.", link_columns={"snapshot_label": "snapshot_link"})}
      </section>
      <section class='card'>
        <h2>Active lessons</h2>
        <p class='muted'>The highest-signal cards from the latest snapshot, sorted by actionability and priority.</p>
        {_render_dataframe(active_cards.head(12), columns=["priority", "source", "title", "recommendation"], empty_message="No active lessons available.")}
      </section>
      <section class='card'>
        <h2>Recurring lessons</h2>
        <p class='muted'>Signals that have repeated across multiple snapshots and are worth treating as campaign memory.</p>
        {_render_dataframe(recurring.head(10), columns=["title", "source", "snapshot_count", "latest_snapshot", "latest_recommendation"], empty_message="No recurring lessons detected yet.")}
      </section>
      <section class='card'>
        <h2>Coverage by source</h2>
        <p class='muted'>Where the dashboard’s learning density is coming from across compare, SAR, MMP, and decision outputs.</p>
        {_render_dataframe(source_summary, columns=["source", "card_count", "snapshot_count"], empty_message="No source coverage available.")}
      </section>
    </div>
  </main>
</body>
</html>
"""
    (output_dir / "index.html").write_text(html, encoding="utf-8")


def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(
        prog="molscope dashboard",
        description="Build or update a persistent campaign dashboard from MolScope learnings bundles.",
    )
    ap.add_argument(
        "--learnings-dir",
        action="append",
        default=[],
        help="Directory containing a campaign learnings bundle. Repeat to add multiple snapshots.",
    )
    ap.add_argument(
        "-o",
        "--output-dir",
        default=None,
        help="Persistent dashboard directory (default: dashboard/campaign_dashboard_<timestamp>)",
    )
    ap.add_argument("--title", default="MolScope Campaign Dashboard", help="Dashboard title")
    ap.add_argument(
        "--max-active-cards",
        type=int,
        default=12,
        help="Maximum number of active learning cards to highlight from the latest snapshot",
    )
    ap.add_argument(
        "--max-recurring",
        type=int,
        default=10,
        help="Maximum number of recurring lessons to highlight",
    )
    args = ap.parse_args(list(argv) if argv is not None else None)

    output_dir = Path(args.output_dir or _default_output_dir())
    output_dir.mkdir(parents=True, exist_ok=True)
    snapshots_root = output_dir / "snapshots"
    snapshots_root.mkdir(parents=True, exist_ok=True)

    snapshot_registry_path = output_dir / "snapshot_registry.csv"
    knowledge_registry_path = output_dir / "knowledge_registry.csv"

    snapshot_registry = _read_csv(snapshot_registry_path)
    knowledge_registry = _read_csv(knowledge_registry_path)

    imported_snapshot_ids: list[str] = []
    imported_sources: list[str] = []
    imported_dirs = [Path(path) for path in args.learnings_dir]

    if not imported_dirs and snapshot_registry.empty and knowledge_registry.empty:
        raise SystemExit(
            "molscope dashboard: provide at least one --learnings-dir or point -o to an existing dashboard"
        )

    for learnings_dir in imported_dirs:
        if not learnings_dir.exists() or not learnings_dir.is_dir():
            raise SystemExit(
                f"molscope dashboard: learnings directory not found: {learnings_dir}"
            )
        cards_path = learnings_dir / "learning_cards.csv"
        if not cards_path.exists():
            raise SystemExit(
                f"molscope dashboard: missing learning_cards.csv in {learnings_dir}"
            )

        summary = _read_json(learnings_dir / "summary.json")
        cards = _read_csv(cards_path)
        snapshot_id = _snapshot_id(learnings_dir, summary)
        snapshot_label = _snapshot_label(learnings_dir, summary)
        stage = _as_text(summary.get("stage"), "(none)")
        sources = summary.get("sources") or sorted(
            cards.get("source", pd.Series(dtype=str)).dropna().astype(str).unique().tolist()
        )
        imported_at = datetime.now().isoformat(timespec="microseconds")
        snapshot_link = f"snapshots/{snapshot_id}/index.html"

        if not snapshot_registry.empty:
            snapshot_registry = snapshot_registry.loc[
                snapshot_registry["snapshot_id"] != snapshot_id
            ]
        if not knowledge_registry.empty:
            knowledge_registry = knowledge_registry.loc[
                knowledge_registry["snapshot_id"] != snapshot_id
            ]

        _replace_tree(learnings_dir, snapshots_root / snapshot_id)

        cards = cards.copy()
        if cards.empty:
            cards = pd.DataFrame(
                columns=[
                    "source",
                    "learning_type",
                    "title",
                    "evidence",
                    "recommendation",
                    "priority",
                ]
            )
        cards.insert(0, "snapshot_id", snapshot_id)
        cards.insert(1, "snapshot_label", snapshot_label)
        cards.insert(2, "snapshot_stage", stage)
        cards.insert(3, "snapshot_imported_at", imported_at)
        cards.insert(4, "snapshot_link", snapshot_link)

        snapshot_row = pd.DataFrame(
            [
                {
                    "snapshot_id": snapshot_id,
                    "snapshot_label": snapshot_label,
                    "stage": stage,
                    "learning_cards": int(summary.get("learning_cards") or len(cards)),
                    "sources": ", ".join(sources) if isinstance(sources, list) else _as_text(sources),
                    "imported_at": imported_at,
                    "snapshot_link": snapshot_link,
                    "source_dir": str(learnings_dir.resolve()),
                }
            ]
        )

        snapshot_registry = pd.concat([snapshot_registry, snapshot_row], ignore_index=True)
        knowledge_registry = pd.concat([knowledge_registry, cards], ignore_index=True)
        imported_snapshot_ids.append(snapshot_id)
        imported_sources.append(str(learnings_dir.resolve()))

    if snapshot_registry.empty and knowledge_registry.empty:
        raise SystemExit("molscope dashboard: no learnings bundles available to index")

    snapshot_registry = snapshot_registry.sort_values(
        ["imported_at", "snapshot_label"], ascending=[False, True], na_position="last"
    ).reset_index(drop=True)
    knowledge_registry = knowledge_registry.sort_values(
        ["snapshot_imported_at", "snapshot_label", "priority", "title"],
        ascending=[False, True, True, True],
        na_position="last",
    ).reset_index(drop=True)

    latest_snapshot = (
        snapshot_registry.iloc[0] if not snapshot_registry.empty else pd.Series(dtype=object)
    )
    latest_snapshot_id = _as_text(latest_snapshot.get("snapshot_id"))
    active_cards = knowledge_registry.loc[
        knowledge_registry["snapshot_id"] == latest_snapshot_id
    ].copy()
    active_cards = _sort_cards(active_cards).head(args.max_active_cards)

    recurring = _build_recurring_lessons(knowledge_registry).head(args.max_recurring)
    source_summary = (
        knowledge_registry.groupby("source", dropna=False)
        .agg(card_count=("title", "size"), snapshot_count=("snapshot_id", "nunique"))
        .reset_index()
        .sort_values(["card_count", "snapshot_count", "source"], ascending=[False, False, True])
    )

    _write_csv(snapshot_registry, snapshot_registry_path)
    _write_csv(knowledge_registry, knowledge_registry_path)
    _write_csv(active_cards, output_dir / "active_learning_cards.csv")
    _write_csv(recurring, output_dir / "recurring_lessons.csv")
    _write_csv(source_summary, output_dir / "source_summary.csv")

    summary = {
        "title": args.title,
        "snapshots": int(len(snapshot_registry)),
        "total_learning_cards": int(len(knowledge_registry)),
        "latest_snapshot": _as_text(latest_snapshot.get("snapshot_label")),
        "latest_stage": _as_text(latest_snapshot.get("stage")),
        "sources": source_summary["source"].astype(str).tolist() if not source_summary.empty else [],
        "recurring_lessons": int(len(recurring)),
        "imported_snapshot_ids": imported_snapshot_ids,
    }
    (output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    _write_markdown(
        output_dir,
        title=args.title,
        snapshot_registry=snapshot_registry,
        active_cards=active_cards,
        recurring=recurring,
        summary=summary,
    )
    _write_html(
        output_dir,
        title=args.title,
        summary=summary,
        snapshot_registry=snapshot_registry,
        active_cards=active_cards,
        recurring=recurring,
        source_summary=source_summary,
    )

    write_run_metadata(
        tool="molscope-dashboard",
        output_table_path=output_dir / "index.html",
        parameters={
            "title": args.title,
            "learnings_dirs": args.learnings_dir,
            "max_active_cards": args.max_active_cards,
            "max_recurring": args.max_recurring,
        },
        artifacts={
            "snapshot_registry": str(snapshot_registry_path),
            "knowledge_registry": str(knowledge_registry_path),
            "active_learning_cards": str(output_dir / "active_learning_cards.csv"),
            "recurring_lessons": str(output_dir / "recurring_lessons.csv"),
            "source_summary": str(output_dir / "source_summary.csv"),
            "imported_sources": imported_sources,
        },
        notes="Persistent campaign dashboard built from MolScope learnings bundles.",
    )

    print(f"Campaign dashboard written to: {output_dir}")
    print(f"- Snapshots tracked: {summary['snapshots']}")
    print("- Entry point: index.html")


if __name__ == "__main__":
    main()
