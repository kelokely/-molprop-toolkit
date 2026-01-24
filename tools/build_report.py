#!/usr/bin/env python3
"""Build shareable report artifacts from a results CSV.

Outputs a directory containing:
- report.md
- report.html
- plots/*.png

Designed to be registry-driven via molprop_toolkit.core.registry.CATEGORY_SPECS.

Example:
  python tools/build_report.py results.csv --outdir reports/run1

If you want to limit categories:
  python tools/build_report.py results.csv --categories cns_mpo,toxicity
"""

from __future__ import annotations

import argparse
import sys
from datetime import datetime
from pathlib import Path

# Allow running directly without installing the package.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from molprop_toolkit.reporting import ReportBuilder  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser(description="MolProp Toolkit report builder (Markdown + HTML + plots)")
    ap.add_argument("input", help="Results CSV")
    ap.add_argument(
        "--outdir",
        help="Output directory (default: reports/<csvstem>_report_<timestamp>)",
        default=None,
    )
    ap.add_argument(
        "--categories",
        help="Comma-separated category keys to include (default: all)",
        default=None,
    )
    ap.add_argument("--top", type=int, default=10, help="Top rows to show per category")
    ap.add_argument("--no-plots", action="store_true", help="Disable plot generation")
    ap.add_argument("--no-md", action="store_true", help="Do not write report.md")
    ap.add_argument("--no-html", action="store_true", help="Do not write report.html")
    ap.add_argument("--title", default="MolProp Toolkit Report", help="Report title")
    args = ap.parse_args()

    csv_path = args.input
    stem = Path(csv_path).stem
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    outdir = args.outdir or str(Path("reports") / f"{stem}_report_{ts}")

    categories = None
    if args.categories:
        categories = [c.strip() for c in args.categories.split(",") if c.strip()]

    builder = ReportBuilder(csv_path=csv_path, categories=categories, top_n=args.top)
    artifacts = builder.build(
        outdir=outdir,
        write_markdown=not args.no_md,
        write_html=not args.no_html,
        make_plots=not args.no_plots,
        title=args.title,
    )

    print(f"Report written to: {artifacts.outdir}")
    if not args.no_md:
        print(f"- Markdown: {artifacts.markdown_path}")
    if not args.no_html:
        print(f"- HTML: {artifacts.html_path}")
    if not args.no_plots:
        print(f"- Plots dir: {artifacts.plots_dir}")


if __name__ == "__main__":
    main()

