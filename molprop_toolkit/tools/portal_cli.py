#!/usr/bin/env python3
"""Assemble MolScope artifact bundles into a single campaign portal."""

from __future__ import annotations

import argparse
import csv
import json
import shutil
from datetime import datetime
from pathlib import Path
from typing import Any, Optional, Sequence

from molprop_toolkit.core.metadata import write_run_metadata


def _default_output_dir() -> str:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return str(Path("portal") / f"campaign_portal_{timestamp}")


def _copy_tree(src: Path, dst: Path) -> None:
    shutil.copytree(src, dst, dirs_exist_ok=True)


def _read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def _html_escape(text: Any) -> str:
    return (
        str(text)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def _count_csv_rows(path: Path) -> Optional[int]:
    if not path.exists():
        return None
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle)
        try:
            next(reader)
        except StopIteration:
            return 0
        return sum(1 for _ in reader)


def _picklists_summary(path: Path) -> tuple[list[tuple[str, str]], list[tuple[str, str]]]:
    metadata = _read_json(path / "picklists_run.metadata.json")
    csv_count = len(list(path.glob("*.csv")))
    html_available = (path / "index.html").exists()
    summary_rows = [
        ("stage", str(metadata.get("stage") or "(none)")),
        ("csv_files", str(csv_count)),
        ("html_index", "yes" if html_available else "no"),
    ]
    links = []
    if html_available:
        links.append(("Picklists index", "picklists/index.html"))
    return summary_rows, links


def _report_summary(path: Path) -> tuple[list[tuple[str, str]], list[tuple[str, str]]]:
    plots_dir = path / "plots"
    plot_count = len(list(plots_dir.glob("*.png"))) if plots_dir.exists() else 0
    summary_rows = [
        ("html_report", "yes" if (path / "report.html").exists() else "no"),
        ("markdown_report", "yes" if (path / "report.md").exists() else "no"),
        ("plot_count", str(plot_count)),
    ]
    links = []
    if (path / "report.html").exists():
        links.append(("Report HTML", "report/report.html"))
    if (path / "report.md").exists():
        links.append(("Report Markdown", "report/report.md"))
    return summary_rows, links


def _compare_summary(path: Path) -> tuple[list[tuple[str, str]], list[tuple[str, str]]]:
    summary = _read_json(path / "summary.json")
    round_review = _read_json(path / "round_review.json")
    summary_rows = [
        ("reference", str(summary.get("reference_name") or "(unknown)")),
        ("candidate", str(summary.get("candidate_name") or "(unknown)")),
        ("shared_ids", str(summary.get("shared_ids") or 0)),
        ("group_column", str(summary.get("group_column") or "(none)")),
        ("positive_shifts", str(len(round_review.get("positive_shifts") or []))),
        ("negative_shifts", str(len(round_review.get("negative_shifts") or []))),
    ]
    links = []
    if (path / "summary.md").exists():
        links.append(("Summary Markdown", "compare/summary.md"))
    if (path / "round_review.md").exists():
        links.append(("Round review", "compare/round_review.md"))
    if (path / "endpoint_deltas.csv").exists():
        links.append(("Endpoint deltas", "compare/endpoint_deltas.csv"))
    return summary_rows, links


def _learnings_summary(path: Path) -> tuple[list[tuple[str, str]], list[tuple[str, str]]]:
    summary = _read_json(path / "summary.json")
    summary_rows = [
        ("learning_cards", str(summary.get("learning_cards") or 0)),
        ("sources", ", ".join(summary.get("sources") or []) or "(none)"),
        ("stage", str(summary.get("stage") or "(none)")),
    ]
    links = []
    if (path / "index.html").exists():
        links.append(("Learnings index", "learnings/index.html"))
    if (path / "summary.md").exists():
        links.append(("Summary Markdown", "learnings/summary.md"))
    if (path / "learning_cards.csv").exists():
        links.append(("Learning cards", "learnings/learning_cards.csv"))
    return summary_rows, links


def _dashboard_summary(path: Path) -> tuple[list[tuple[str, str]], list[tuple[str, str]]]:
    summary = _read_json(path / "summary.json")
    summary_rows = [
        ("snapshots", str(summary.get("snapshots") or 0)),
        ("learning_cards", str(summary.get("total_learning_cards") or 0)),
        ("latest_snapshot", str(summary.get("latest_snapshot") or "(none)")),
        ("latest_stage", str(summary.get("latest_stage") or "(none)")),
    ]
    links = []
    if (path / "index.html").exists():
        links.append(("Dashboard index", "dashboard/index.html"))
    if (path / "summary.md").exists():
        links.append(("Summary Markdown", "dashboard/summary.md"))
    if (path / "snapshot_registry.csv").exists():
        links.append(("Snapshot registry", "dashboard/snapshot_registry.csv"))
    if (path / "knowledge_registry.csv").exists():
        links.append(("Knowledge registry", "dashboard/knowledge_registry.csv"))
    return summary_rows, links


def _sar_summary(path: Path) -> tuple[list[tuple[str, str]], list[tuple[str, str]]]:
    summary = _read_json(path / "summary.json")
    summary_rows = [
        ("groups", str(summary.get("groups") or 0)),
        ("rows", str(summary.get("rows") or 0)),
        ("priority_metric", str(summary.get("priority_metric") or "(none)")),
    ]
    links = []
    if (path / "summary.md").exists():
        links.append(("Summary Markdown", "sar/summary.md"))
    if (path / "series_summary.csv").exists():
        links.append(("Series summary", "sar/series_summary.csv"))
    if (path / "top_compounds.csv").exists():
        links.append(("Top compounds", "sar/top_compounds.csv"))
    return summary_rows, links


def _mmp_summary(path: Path) -> tuple[list[tuple[str, str]], list[tuple[str, str]]]:
    summary = _read_json(path / "summary.json")
    summary_rows = [
        ("pair_rows", str(summary.get("pair_rows") or 0)),
        ("transforms", str(summary.get("transforms") or 0)),
        ("priority_metric", str(summary.get("priority_metric") or "(none)")),
    ]
    links = []
    if (path / "summary.md").exists():
        links.append(("Summary Markdown", "mmp/summary.md"))
    if (path / "transform_summary.csv").exists():
        links.append(("Transform summary", "mmp/transform_summary.csv"))
    return summary_rows, links


def _retro_summary(path: Path) -> tuple[list[tuple[str, str]], list[tuple[str, str]]]:
    routes_csv = path / "routes_summary.csv"
    summary_rows = [
        ("route_rows", str(_count_csv_rows(routes_csv) or 0)),
        ("html_index", "yes" if (path / "index.html").exists() else "no"),
        ("retro_enriched_table", "yes" if (path / "retro_enriched.csv").exists() else "no"),
    ]
    links = []
    if (path / "index.html").exists():
        links.append(("Retro index", "retro/index.html"))
    if routes_csv.exists():
        links.append(("Route summary", "retro/routes_summary.csv"))
    if (path / "retro_enriched.csv").exists():
        links.append(("Retro-enriched table", "retro/retro_enriched.csv"))
    return summary_rows, links


SECTION_CONFIG: tuple[tuple[str, str, str, Any], ...] = (
    ("report", "Report bundle", "Shareable campaign report with stage-aware framing and decision products.", _report_summary),
    ("picklists", "Decision picklists", "Operational queues for make/test/de-risk/hold workflows.", _picklists_summary),
    ("compare", "Campaign compare", "Round-over-round or lead-set deltas for endpoints and groups.", _compare_summary),
    ("learnings", "Campaign learnings", "Reusable project lessons distilled from compare, SAR, MMP, and decision outputs.", _learnings_summary),
    ("dashboard", "Campaign dashboard", "Persistent cross-round dashboard and knowledge registry built from learnings bundles.", _dashboard_summary),
    ("sar", "SAR bundle", "Scaffold-centric series summary and top compounds.", _sar_summary),
    ("mmp", "Matched molecular pairs", "Transform summaries and pair-level deltas.", _mmp_summary),
    ("retro", "Retrosynthesis", "Route-aware feasibility bundle for synthesis-facing prioritization.", _retro_summary),
)


def _render_table(rows: Sequence[tuple[str, str]]) -> str:
    items = "".join(
        f"<tr><th>{_html_escape(key)}</th><td>{_html_escape(value)}</td></tr>"
        for key, value in rows
    )
    return f"<table class='summary-table'><tbody>{items}</tbody></table>"


def _render_links(links: Sequence[tuple[str, str]]) -> str:
    if not links:
        return "<p class='muted'>No linked artifacts detected.</p>"
    return "<div class='link-list'>" + "".join(
        f"<a href='{_html_escape(href)}'>{_html_escape(label)}</a>"
        for label, href in links
    ) + "</div>"


def _write_portal_index(
    output_dir: Path,
    *,
    title: str,
    cards: Sequence[dict[str, Any]],
) -> None:
    card_html = []
    for card in cards:
        card_html.append(
            "<section class='card'>"
            f"<div class='eyebrow'>{_html_escape(card['slug'])}</div>"
            f"<h2>{_html_escape(card['title'])}</h2>"
            f"<p class='muted'>{_html_escape(card['description'])}</p>"
            f"{_render_table(card['summary_rows'])}"
            f"{_render_links(card['links'])}"
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
      margin: 0;
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
      margin: 0 0 8px;
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
    .muted {{
      color: var(--muted);
    }}
    .summary-table {{
      width: 100%;
      border-collapse: collapse;
      margin: 14px 0;
      font-size: 14px;
    }}
    .summary-table th,
    .summary-table td {{
      text-align: left;
      vertical-align: top;
      padding: 8px 0;
      border-bottom: 1px solid var(--border);
    }}
    .summary-table th {{
      width: 42%;
      color: var(--muted);
      font-weight: 600;
      padding-right: 12px;
    }}
    .link-list {{
      display: flex;
      flex-wrap: wrap;
      gap: 10px;
      margin-top: 16px;
    }}
    .link-list a {{
      color: var(--brand);
      text-decoration: none;
      font-weight: 600;
    }}
    .link-list a:hover {{
      text-decoration: underline;
    }}
  </style>
</head>
<body>
  <main class='shell'>
    <section class='hero'>
      <div class='eyebrow'>Campaign portal</div>
      <h1>{_html_escape(title)}</h1>
      <p>
        This bundle assembles MolScope artifacts into one workspace so a team can move across report, picklists,
        compare, learnings, dashboard, SAR, MMP, and retrosynthesis outputs without juggling separate directories.
      </p>
    </section>
    <div class='grid'>
      {''.join(card_html)}
    </div>
  </main>
</body>
</html>
"""
    (output_dir / "index.html").write_text(html, encoding="utf-8")


def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(
        prog="molscope portal",
        description="Assemble existing MolScope artifacts into a single linked campaign portal.",
    )
    ap.add_argument(
        "-o",
        "--output-dir",
        default=None,
        help="Output directory for the assembled portal (default: portal/campaign_portal_<timestamp>)",
    )
    ap.add_argument("--title", default="MolScope Campaign Portal", help="Portal title")
    ap.add_argument("--report-dir", default=None, help="Directory containing report.html/report.md")
    ap.add_argument("--picklists-dir", default=None, help="Directory containing picklist CSV/HTML outputs")
    ap.add_argument("--compare-dir", default=None, help="Directory containing compare summary outputs")
    ap.add_argument("--learnings-dir", default=None, help="Directory containing campaign learnings outputs")
    ap.add_argument("--dashboard-dir", default=None, help="Directory containing persistent campaign dashboard outputs")
    ap.add_argument("--sar-dir", default=None, help="Directory containing SAR bundle outputs")
    ap.add_argument("--mmp-dir", default=None, help="Directory containing MMP bundle outputs")
    ap.add_argument("--retro-dir", default=None, help="Directory containing retrosynthesis outputs")
    args = ap.parse_args(list(argv) if argv is not None else None)

    output_dir = Path(args.output_dir or _default_output_dir())
    output_dir.mkdir(parents=True, exist_ok=True)

    section_paths = {
        "report": args.report_dir,
        "picklists": args.picklists_dir,
        "compare": args.compare_dir,
        "learnings": args.learnings_dir,
        "dashboard": args.dashboard_dir,
        "sar": args.sar_dir,
        "mmp": args.mmp_dir,
        "retro": args.retro_dir,
    }
    selected = {slug: Path(path) for slug, path in section_paths.items() if path}
    if not selected:
        raise SystemExit("molscope portal: provide at least one artifact directory")

    cards: list[dict[str, Any]] = []
    for slug, title, description, summarizer in SECTION_CONFIG:
        src = selected.get(slug)
        if src is None:
            continue
        if not src.exists() or not src.is_dir():
            raise SystemExit(f"molscope portal: artifact directory not found: {src}")
        dst = output_dir / slug
        _copy_tree(src, dst)
        summary_rows, links = summarizer(dst)
        cards.append(
            {
                "slug": slug,
                "title": title,
                "description": description,
                "summary_rows": summary_rows,
                "links": links,
            }
        )

    _write_portal_index(output_dir, title=args.title, cards=cards)

    write_run_metadata(
        tool="molscope-portal",
        output_table_path=output_dir / "index.html",
        parameters={
            "title": args.title,
            "sections": list(selected.keys()),
        },
        artifacts={
            "report_dir": args.report_dir,
            "picklists_dir": args.picklists_dir,
            "compare_dir": args.compare_dir,
            "learnings_dir": args.learnings_dir,
            "dashboard_dir": args.dashboard_dir,
            "sar_dir": args.sar_dir,
            "mmp_dir": args.mmp_dir,
            "retro_dir": args.retro_dir,
        },
        notes="Portal bundle assembled from existing MolScope artifact directories.",
    )

    print(f"Campaign portal written to: {output_dir}")
    print("- Entry point: index.html")


if __name__ == "__main__":
    main()
