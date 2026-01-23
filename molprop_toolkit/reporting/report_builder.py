"""Registry-driven report builder (Task 12).

Outputs:
- report.md (Markdown)
- report.html (self-contained HTML referencing local plot images)
- plots/*.png

The builder is intentionally generic:
- It discovers which categories are present by scanning the dataset columns and
  matching against CATEGORY_SPECS.
- It creates summaries + plots primarily using each category's `key_column`
  (when present) and the available columns.

This file is dependency-light:
- Requires pandas (already required by analyzers)
- Uses matplotlib if plots are enabled

If matplotlib is not available, plot generation can be disabled.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

import pandas as pd

from molprop_toolkit.core import detect_id_column, read_csv
from molprop_toolkit.core.registry import CATEGORY_SPECS, CategorySpec

try:  # pragma: no cover
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except Exception:  # pragma: no cover
    plt = None  # type: ignore[assignment]


@dataclass(frozen=True)
class ReportArtifacts:
    outdir: Path
    markdown_path: Path
    html_path: Path
    plots_dir: Path


def _slug(s: str) -> str:
    s2 = "".join(ch if (ch.isalnum() or ch in "-_") else "_" for ch in (s or ""))
    s2 = s2.strip("_-") or "item"
    return s2[:80]


def _is_numeric(series: pd.Series) -> bool:
    # Treat bool-like columns as categorical for reporting purposes.
    return pd.api.types.is_numeric_dtype(series) and not pd.api.types.is_bool_dtype(series)


def _safe_value(v: object) -> str:
    if v is None:
        return ""
    if isinstance(v, float):
        if pd.isna(v):
            return ""
        return f"{v:.3g}"
    return str(v)


def _markdown_table(headers: Sequence[str], rows: Sequence[Sequence[object]]) -> str:
    # Very small markdown table helper (no external deps).
    def esc(x: object) -> str:
        return _safe_value(x).replace("|", "\\|")

    head = "| " + " | ".join(esc(h) for h in headers) + " |"
    sep = "| " + " | ".join(["---"] * len(headers)) + " |"
    body = "\n".join("| " + " | ".join(esc(c) for c in r) + " |" for r in rows)
    return "\n".join([head, sep, body])


def _html_escape(s: str) -> str:
    return (
        (s or "")
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def _html_table(headers: Sequence[str], rows: Sequence[Sequence[object]]) -> str:
    ths = "".join(f"<th>{_html_escape(str(h))}</th>" for h in headers)
    trs = []
    for r in rows:
        tds = "".join(f"<td>{_html_escape(_safe_value(c))}</td>" for c in r)
        trs.append(f"<tr>{tds}</tr>")
    return f"<table><thead><tr>{ths}</tr></thead><tbody>{''.join(trs)}</tbody></table>"


class ReportBuilder:
    def __init__(
        self,
        csv_path: str,
        categories: Optional[Sequence[str]] = None,
        top_n: int = 10,
    ):
        self.csv_path = csv_path
        self.df = read_csv(csv_path)
        self.id_col = detect_id_column(self.df)
        self.top_n = int(top_n)

        if categories is None:
            self.categories = list(CATEGORY_SPECS.keys())
        else:
            unknown = [c for c in categories if c not in CATEGORY_SPECS]
            if unknown:
                raise ValueError(f"Unknown categories: {unknown}")
            self.categories = list(categories)

    def category_columns(self, spec: CategorySpec) -> List[str]:
        cols: List[str] = []
        cols.extend([c for c in spec.columns if c in self.df.columns])
        if spec.prefix:
            for c in self.df.columns:
                if c.startswith(spec.prefix) and c not in cols:
                    cols.append(c)
        return cols

    def categories_present(self) -> List[str]:
        present: List[str] = []
        for k in self.categories:
            spec = CATEGORY_SPECS[k]
            if self.category_columns(spec):
                present.append(k)
        return present

    def _key_series(self, spec: CategorySpec, available_cols: Sequence[str]) -> Tuple[Optional[str], Optional[pd.Series]]:
        key = spec.key_column
        if key and key in self.df.columns:
            return key, self.df[key]

        # Fallback: first available numeric column.
        for c in available_cols:
            if c in self.df.columns and _is_numeric(self.df[c]):
                return c, self.df[c]

        # Fallback: first available col.
        if available_cols:
            c0 = available_cols[0]
            return c0, self.df[c0]

        return None, None

    def _summarize_key(self, key: str, s: pd.Series) -> Tuple[str, List[Tuple[str, str]]]:
        # Returns (kind, summary_items)
        if _is_numeric(s):
            vals = pd.to_numeric(s, errors="coerce").dropna()
            if len(vals) == 0:
                return "numeric", [("non-null", "0")]
            summary = [
                ("non-null", str(len(vals))),
                ("mean", f"{vals.mean():.3g}"),
                ("median", f"{vals.median():.3g}"),
                ("min", f"{vals.min():.3g}"),
                ("max", f"{vals.max():.3g}"),
                ("std", f"{vals.std():.3g}"),
            ]
            return "numeric", summary

        vc = s.fillna("(missing)").astype(str).value_counts(dropna=False)
        summary = [(f"count[{k}]", str(v)) for k, v in vc.head(8).items()]
        return "categorical", summary

    def _plot_numeric_hist(self, title: str, s: pd.Series, outpath: Path) -> None:
        if plt is None:
            raise RuntimeError("matplotlib is not available; disable plots or install matplotlib")

        vals = pd.to_numeric(s, errors="coerce").dropna()
        if len(vals) == 0:
            return

        fig = plt.figure(figsize=(7.2, 4.2), dpi=160)
        ax = fig.add_subplot(1, 1, 1)
        ax.hist(vals, bins=30, color="#2E86C1", alpha=0.85, edgecolor="white")
        mean = float(vals.mean())
        med = float(vals.median())
        ax.axvline(mean, color="#C0392B", linestyle="--", linewidth=1.5, label=f"mean={mean:.3g}")
        ax.axvline(med, color="#27AE60", linestyle=":", linewidth=1.5, label=f"median={med:.3g}")
        ax.set_title(title)
        ax.set_xlabel("value")
        ax.set_ylabel("count")
        ax.legend(frameon=False)
        fig.tight_layout()
        outpath.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(outpath)
        plt.close(fig)

    def _plot_categorical_bar(self, title: str, s: pd.Series, outpath: Path) -> None:
        if plt is None:
            raise RuntimeError("matplotlib is not available; disable plots or install matplotlib")

        vc = s.fillna("(missing)").astype(str).value_counts(dropna=False).head(12)
        if len(vc) == 0:
            return

        fig = plt.figure(figsize=(7.2, 4.2), dpi=160)
        ax = fig.add_subplot(1, 1, 1)
        ax.bar(vc.index.tolist(), vc.values.tolist(), color="#8E44AD", alpha=0.9)
        ax.set_title(title)
        ax.set_xlabel("class")
        ax.set_ylabel("count")
        ax.tick_params(axis="x", rotation=35)
        fig.tight_layout()
        outpath.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(outpath)
        plt.close(fig)

    def _top_table(self, key: str, s: pd.Series, available_cols: Sequence[str]) -> Tuple[List[str], List[List[object]]]:
        base_cols = [self.id_col]
        if key not in base_cols:
            base_cols.append(key)

        # Add up to 4 extra columns for context.
        for c in available_cols:
            if c not in base_cols:
                base_cols.append(c)
            if len(base_cols) >= 6:
                break

        df2 = self.df[[c for c in base_cols if c in self.df.columns]].copy()

        if _is_numeric(s):
            df2[key] = pd.to_numeric(df2[key], errors="coerce")
            df2 = df2.sort_values(key, ascending=False)
        else:
            # For categorical, just show the most frequent class first.
            df2[key] = df2[key].fillna("(missing)").astype(str)
            top_class = df2[key].value_counts().index[0] if len(df2) else "(missing)"
            df2 = df2.sort_values(by=[key, self.id_col])
            df2 = pd.concat([df2[df2[key] == top_class], df2[df2[key] != top_class]], axis=0)

        df2 = df2.head(self.top_n)
        headers = list(df2.columns)
        rows = df2.values.tolist()
        return headers, rows

    def build(
        self,
        outdir: str | Path,
        write_markdown: bool = True,
        write_html: bool = True,
        make_plots: bool = True,
        title: str = "MolProp Toolkit Report",
    ) -> ReportArtifacts:
        outdir_p = Path(outdir)
        outdir_p.mkdir(parents=True, exist_ok=True)
        plots_dir = outdir_p / "plots"
        plots_dir.mkdir(parents=True, exist_ok=True)

        md_path = outdir_p / "report.md"
        html_path = outdir_p / "report.html"

        present = self.categories_present()
        generated = datetime.now().isoformat(timespec="seconds")

        # Collect per-category artifacts for both MD and HTML.
        cat_blocks_md: List[str] = []
        cat_blocks_html: List[str] = []

        for cat_key in present:
            spec = CATEGORY_SPECS[cat_key]
            cols = self.category_columns(spec)
            key, series = self._key_series(spec, cols)
            if key is None or series is None:
                continue

            kind, summary_items = self._summarize_key(key, series)

            plot_rel: Optional[str] = None
            if make_plots:
                plot_name = f"{_slug(cat_key)}__{_slug(key)}.png"
                plot_path = plots_dir / plot_name
                try:
                    if kind == "numeric":
                        self._plot_numeric_hist(f"{spec.name}: {key}", series, plot_path)
                    else:
                        self._plot_categorical_bar(f"{spec.name}: {key}", series, plot_path)
                    plot_rel = f"plots/{plot_name}"
                except Exception as e:
                    plot_rel = None
                    # We keep going and just note plot failure in the report.
                    summary_items.append(("plot_error", str(e)))

            top_headers, top_rows = self._top_table(key, series, cols)

            md = []
            md.append(f"## {spec.name} (`{cat_key}`)")
            md.append("")
            md.append(spec.description.strip())
            md.append("")
            md.append(f"Key metric: `{key}`")
            md.append("")
            md.append(_markdown_table(["metric", "value"], summary_items))
            md.append("")
            if plot_rel:
                md.append(f"![]({plot_rel})")
                md.append("")
            md.append("Top compounds (preview):")
            md.append("")
            md.append(_markdown_table(top_headers, top_rows))
            md.append("")
            md.append(f"Available columns in dataset for this category: {len(cols)}")
            md.append("")
            cat_blocks_md.append("\n".join(md))

            h = []
            h.append(f"<section class='category' id='{_html_escape(cat_key)}'>")
            h.append(f"<h2>{_html_escape(spec.name)} <code>{_html_escape(cat_key)}</code></h2>")
            h.append(f"<p class='desc'>{_html_escape(spec.description.strip())}</p>")
            h.append(f"<p><strong>Key metric:</strong> <code>{_html_escape(key)}</code></p>")
            h.append(_html_table(["metric", "value"], summary_items))
            if plot_rel:
                h.append(
                    f"<div class='plot'><img src='{_html_escape(plot_rel)}' alt='{_html_escape(spec.name)} plot'/></div>"
                )
            h.append("<h3>Top compounds (preview)</h3>")
            h.append(_html_table(top_headers, top_rows))
            h.append(f"<p class='meta'>Available columns in dataset for this category: {len(cols)}</p>")
            h.append("</section>")
            cat_blocks_html.append("\n".join(h))

        # Overview table
        overview_rows: List[List[object]] = []
        for cat_key in present:
            spec = CATEGORY_SPECS[cat_key]
            cols = self.category_columns(spec)
            key = spec.key_column if (spec.key_column and spec.key_column in self.df.columns) else ""
            overview_rows.append([cat_key, spec.name, key, len(cols)])

        md_overview = _markdown_table(
            ["category", "name", "key_column", "available_cols"],
            overview_rows,
        )
        html_overview = _html_table(
            ["category", "name", "key_column", "available_cols"],
            overview_rows,
        )

        md_doc = "\n".join(
            [
                f"# {title}",
                "",
                f"Dataset: `{self.csv_path}`",
                "",
                f"Generated: `{generated}`",
                "",
                f"Compounds: `{len(self.df)}`  ",
                f"Columns: `{len(self.df.columns)}`  ",
                f"Detected ID column: `{self.id_col}`",
                "",
                "## Category coverage",
                "",
                md_overview,
                "",
                "---",
                "",
                "\n\n".join(cat_blocks_md) if cat_blocks_md else "(No categories detected in dataset)",
                "",
            ]
        )

        css = """
        :root{--fg:#111;--muted:#555;--border:#e6e6e6;--bg:#fff;--mono:ui-monospace,SFMono-Regular,Menlo,Monaco,Consolas,monospace;}
        body{font-family:system-ui,-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;margin:24px;line-height:1.45;color:var(--fg);background:var(--bg);}
        h1{margin:0 0 6px 0;font-size:22px;}
        .sub{color:var(--muted);margin:0 0 18px 0;}
        code{font-family:var(--mono);font-size:0.95em;}
        table{border-collapse:collapse;width:100%;margin:10px 0 18px 0;}
        th,td{border:1px solid var(--border);padding:8px 10px;text-align:left;vertical-align:top;}
        th{background:#fafafa;}
        section.category{padding:14px 0;border-top:2px solid #f0f0f0;}
        .desc{white-space:pre-wrap;}
        .meta{color:var(--muted);}
        .plot img{max-width:980px;width:100%;height:auto;border:1px solid var(--border);border-radius:8px;background:#fff;}
        @media (max-width:600px){body{margin:14px;} th,td{font-size:13px;}}
        """

        html_doc = f"""<!doctype html>
<html>
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width,initial-scale=1"/>
<title>{_html_escape(title)}</title>
<style>{css}</style>
</head>
<body>
<h1>{_html_escape(title)}</h1>
<p class="sub">Dataset: <code>{_html_escape(self.csv_path)}</code> | Generated: <code>{_html_escape(generated)}</code> | Compounds: <code>{len(self.df)}</code> | Columns: <code>{len(self.df.columns)}</code> | ID: <code>{_html_escape(self.id_col)}</code></p>
<h2>Category coverage</h2>
{html_overview}
{''.join(cat_blocks_html) if cat_blocks_html else '<p>(No categories detected in dataset)</p>'}
</body>
</html>
"""

        if write_markdown:
            md_path.write_text(md_doc, encoding="utf-8")
        if write_html:
            html_path.write_text(html_doc, encoding="utf-8")

        return ReportArtifacts(outdir=outdir_p, markdown_path=md_path, html_path=html_path, plots_dir=plots_dir)


