"""Registry-driven report builder.

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

from molprop_toolkit.core import detect_id_column, read_table
from molprop_toolkit.core.registry import CATEGORY_SPECS, CategorySpec
from molprop_toolkit.decision_presets import get_stage_preset
from molprop_toolkit.picklists.builtins import BUILTIN_PICKLISTS
from molprop_toolkit.picklists.engine import decision_evidence_buckets, run_picklist

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
    return pd.api.types.is_numeric_dtype(series) and not pd.api.types.is_bool_dtype(
        series
    )


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


def _candidate_column(df: pd.DataFrame, *names: str) -> Optional[str]:
    for name in names:
        if name in df.columns:
            return name
    return None


DEFAULT_DECISION_QUEUE_KEYS: tuple[str, ...] = (
    "make_next",
    "test_next",
    "de_risk_next",
    "hold",
    "backup_series",
)


class ReportBuilder:
    def __init__(
        self,
        csv_path: Optional[str] = None,
        *,
        df: Optional[pd.DataFrame] = None,
        categories: Optional[Sequence[str]] = None,
        top_n: int = 10,
        stage: Optional[str] = None,
    ):
        if df is None:
            if not csv_path:
                raise ValueError("csv_path is required when df is not provided")
            self.csv_path = str(csv_path)
            self.df = read_table(self.csv_path)
        else:
            self.csv_path = str(csv_path or "(in-memory)")
            self.df = df
        self.id_col = detect_id_column(self.df)
        self.top_n = int(top_n)
        self.stage = stage

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

    def _key_series(
        self, spec: CategorySpec, available_cols: Sequence[str]
    ) -> Tuple[Optional[str], Optional[pd.Series]]:
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

    def _summarize_key(
        self, key: str, s: pd.Series
    ) -> Tuple[str, List[Tuple[str, str]]]:
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
            raise RuntimeError(
                "matplotlib is not available; disable plots or install matplotlib"
            )

        vals = pd.to_numeric(s, errors="coerce").dropna()
        if len(vals) == 0:
            return

        fig = plt.figure(figsize=(7.2, 4.2), dpi=160)
        ax = fig.add_subplot(1, 1, 1)
        ax.hist(vals, bins=30, color="#2E86C1", alpha=0.85, edgecolor="white")
        mean = float(vals.mean())
        med = float(vals.median())
        ax.axvline(
            mean,
            color="#C0392B",
            linestyle="--",
            linewidth=1.5,
            label=f"mean={mean:.3g}",
        )
        ax.axvline(
            med,
            color="#27AE60",
            linestyle=":",
            linewidth=1.5,
            label=f"median={med:.3g}",
        )
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
            raise RuntimeError(
                "matplotlib is not available; disable plots or install matplotlib"
            )

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

    def _top_table(
        self, key: str, s: pd.Series, available_cols: Sequence[str]
    ) -> Tuple[List[str], List[List[object]]]:
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
            df2 = pd.concat(
                [df2[df2[key] == top_class], df2[df2[key] != top_class]], axis=0
            )

        df2 = df2.head(self.top_n)
        headers = list(df2.columns)
        rows = df2.values.tolist()
        return headers, rows

    def _retro_section(self) -> tuple[str, str]:
        rerank_col = _candidate_column(self.df, "Retro_Rerank_Score")
        rank_col = _candidate_column(self.df, "Retro_Rerank_Rank")
        solved_col = _candidate_column(self.df, "Retro_IsSolved", "is_solved")
        steps_col = _candidate_column(self.df, "Retro_StepCount", "number_of_steps")
        stock_col = _candidate_column(
            self.df, "Retro_Stock_Coverage", "precursors_in_stock"
        )
        top_score_col = _candidate_column(self.df, "Retro_TopScore", "top_score")

        if rerank_col is None and rank_col is None and solved_col is None:
            return "", ""

        summary_items: list[tuple[str, str]] = []

        if solved_col is not None:
            solved = (
                self.df[solved_col]
                .fillna(False)
                .astype("string")
                .str.lower()
                .isin({"true", "1", "yes"})
            )
            summary_items.append(("solved_targets", str(int(solved.sum()))))
            summary_items.append(("unsolved_targets", str(int((~solved).sum()))))

        if rerank_col is not None:
            rerank_values = pd.to_numeric(self.df[rerank_col], errors="coerce").dropna()
            if not rerank_values.empty:
                summary_items.extend(
                    [
                        ("rerank_score_mean", f"{rerank_values.mean():.3g}"),
                        ("rerank_score_max", f"{rerank_values.max():.3g}"),
                    ]
                )

        if steps_col is not None:
            step_values = pd.to_numeric(self.df[steps_col], errors="coerce").dropna()
            if not step_values.empty:
                summary_items.append(("step_count_median", f"{step_values.median():.3g}"))

        if stock_col is not None and stock_col == "Retro_Stock_Coverage":
            stock_values = pd.to_numeric(self.df[stock_col], errors="coerce").dropna()
            if not stock_values.empty:
                summary_items.append(("stock_coverage_mean", f"{stock_values.mean():.3g}"))

        priority_sort_col = rank_col or rerank_col
        top_cols = [self.id_col]
        for col in (rank_col, rerank_col, steps_col, stock_col, top_score_col):
            if col and col not in top_cols:
                top_cols.append(col)

        df_top = self.df[[c for c in top_cols if c in self.df.columns]].copy()
        if priority_sort_col and priority_sort_col in df_top.columns:
            df_top[priority_sort_col] = pd.to_numeric(
                df_top[priority_sort_col], errors="coerce"
            )
            ascending = priority_sort_col == rank_col
            df_top = df_top.sort_values(
                by=priority_sort_col,
                ascending=ascending,
                na_position="last",
            )
        elif rerank_col and rerank_col in df_top.columns:
            df_top[rerank_col] = pd.to_numeric(df_top[rerank_col], errors="coerce")
            df_top = df_top.sort_values(
                by=rerank_col,
                ascending=False,
                na_position="last",
            )
        df_top = df_top.head(self.top_n)

        headers = list(df_top.columns)
        rows = df_top.values.tolist()

        md = []
        md.append("## Retrosynthesis feasibility")
        md.append("")
        md.append(
            "This section summarizes retrosynthesis-aware prioritization fields when a retro-enriched campaign table is present. "
            "Use it to weigh route support, step burden, and stock coverage alongside the rest of the report."
        )
        md.append("")
        md.append(_markdown_table(["metric", "value"], summary_items))
        md.append("")
        md.append("Top retro-prioritized compounds:")
        md.append("")
        md.append(_markdown_table(headers, rows))
        md.append("")

        h = []
        h.append("<section class='category' id='retrosynthesis'>")
        h.append("<h2>Retrosynthesis feasibility</h2>")
        h.append(
            "<p class='desc'>This section summarizes retrosynthesis-aware prioritization fields when a retro-enriched campaign table is present. "
            "Use it to weigh route support, step burden, and stock coverage alongside the rest of the report.</p>"
        )
        h.append(_html_table(["metric", "value"], summary_items))
        h.append("<h3>Top retro-prioritized compounds</h3>")
        h.append(_html_table(headers, rows))
        h.append("</section>")

        return "\n".join(md), "\n".join(h)

    def _stage_section(self) -> tuple[str, str]:
        preset = get_stage_preset(self.stage)
        if preset is None:
            return "", ""

        available_columns = [
            column for column in preset.priority_columns if column in self.df.columns
        ]
        recommended_picklists = ", ".join(
            f"`{key}`" for key in preset.recommended_picklists
        )

        preview_cols = [self.id_col]
        for column in available_columns[:5]:
            if column not in preview_cols:
                preview_cols.append(column)

        preview_df = self.df[[c for c in preview_cols if c in self.df.columns]].copy()
        sort_column = next(
            (
                column
                for column in (
                    "Measured_pIC50",
                    "Lead_LLE",
                    "Retro_Rerank_Rank",
                    "Retro_Rerank_Score",
                    "CNS_MPO",
                    "Bioavail_Score",
                    "Tox_Total_Alerts",
                    "MedChem_Alerts",
                )
                if column in preview_df.columns
            ),
            None,
        )
        if sort_column is not None:
            preview_df[sort_column] = pd.to_numeric(
                preview_df[sort_column], errors="coerce"
            )
            preview_df = preview_df.sort_values(
                by=sort_column,
                ascending=(sort_column == "Retro_Rerank_Rank"),
                na_position="last",
            )
        preview_df = preview_df.head(self.top_n)

        summary_rows = [
            ("stage", preset.label),
            ("recommended_picklists", recommended_picklists),
            (
                "available_priority_columns",
                ", ".join(f"`{col}`" for col in available_columns)
                if available_columns
                else "(none detected in this dataset)",
            ),
        ]

        md = []
        md.append(f"## Decision stage: {preset.label}")
        md.append("")
        md.append(preset.summary)
        md.append("")
        md.append(_markdown_table(["item", "value"], summary_rows))
        md.append("")
        md.append("Stage focus:")
        md.append("")
        for point in preset.focus_points:
            md.append(f"- {point}")
        md.append("")
        if not preview_df.empty:
            md.append("Priority preview:")
            md.append("")
            md.append(
                _markdown_table(list(preview_df.columns), preview_df.values.tolist())
            )
            md.append("")

        h = []
        h.append(f"<section class='category' id='stage_{_html_escape(preset.key)}'>")
        h.append(f"<h2>Decision stage: {_html_escape(preset.label)}</h2>")
        h.append(f"<p class='desc'>{_html_escape(preset.summary)}</p>")
        h.append(_html_table(["item", "value"], summary_rows))
        h.append("<h3>Stage focus</h3>")
        h.append(
            "<ul>"
            + "".join(f"<li>{_html_escape(point)}</li>" for point in preset.focus_points)
            + "</ul>"
        )
        if not preview_df.empty:
            h.append("<h3>Priority preview</h3>")
            h.append(
                _html_table(list(preview_df.columns), preview_df.values.tolist())
            )
        h.append("</section>")

        return "\n".join(md), "\n".join(h)

    def _decision_queue_keys(self) -> tuple[str, ...]:
        preset = get_stage_preset(self.stage)
        if preset is not None:
            return tuple(
                key for key in preset.recommended_picklists if key in BUILTIN_PICKLISTS
            )
        return tuple(key for key in DEFAULT_DECISION_QUEUE_KEYS if key in BUILTIN_PICKLISTS)

    def _decision_priority_columns(
        self,
        queue_keys: Sequence[str],
    ) -> list[str]:
        preset = get_stage_preset(self.stage)
        candidate_columns: list[str] = []
        if preset is not None:
            candidate_columns.extend(preset.priority_columns)
        for key in queue_keys:
            definition = BUILTIN_PICKLISTS.get(key)
            if definition is None:
                continue
            candidate_columns.extend(definition.evidence_columns)

        columns: list[str] = []
        for column in candidate_columns:
            if column not in columns and (
                column in self.df.columns or f"{column}_Provenance" in self.df.columns
            ):
                columns.append(column)
        return columns

    def _decision_evidence_mix(
        self,
        priority_columns: Sequence[str],
    ) -> list[tuple[str, str]]:
        if not priority_columns:
            return [("priority_evidence", "No configured priority evidence columns detected.")]

        measured_rows = 0
        calculated_rows = 0
        missing_rows = 0

        for _, row in self.df.iterrows():
            buckets = decision_evidence_buckets(row, priority_columns)
            if buckets["measured"]:
                measured_rows += 1
            elif buckets["calculated"] or buckets["heuristic"]:
                calculated_rows += 1
            else:
                missing_rows += 1

        return [
            ("rows_with_measured_priority_evidence", str(measured_rows)),
            ("rows_with_calculated_or_heuristic_priority_evidence", str(calculated_rows)),
            ("rows_missing_priority_evidence", str(missing_rows)),
            ("priority_columns", ", ".join(f"`{c}`" for c in priority_columns)),
        ]

    def _decision_section(self) -> tuple[str, str]:
        queue_keys = self._decision_queue_keys()
        if not queue_keys:
            return "", ""

        queue_results = [
            run_picklist(self.df, BUILTIN_PICKLISTS[key], stage=self.stage)
            for key in queue_keys
            if key in BUILTIN_PICKLISTS
        ]
        if not queue_results:
            return "", ""

        priority_columns = self._decision_priority_columns(queue_keys)
        evidence_rows = self._decision_evidence_mix(priority_columns)

        queue_rows: list[list[object]] = []
        preview_blocks_md: list[str] = []
        preview_blocks_html: list[str] = []
        experiment_rows: list[list[object]] = []

        for result in queue_results:
            definition = result.definition
            top_ids = ", ".join(
                result.df[self.id_col].astype(str).head(3).tolist()
            ) if not result.df.empty else "(none)"
            queue_rows.append(
                [
                    definition.queue or definition.key,
                    definition.name,
                    len(result.df),
                    definition.next_action or "Review with project context.",
                    top_ids,
                ]
            )

            if result.df.empty:
                continue

            first_row = result.df.iloc[0]
            experiment_rows.append(
                [
                    definition.queue or definition.key,
                    first_row.get("Decision_Hypothesis", ""),
                    first_row.get("Decision_Minimal_Test_Set", ""),
                    first_row.get("Decision_Minimal_Analog_Set", ""),
                    top_ids,
                ]
            )

            preview_cols = [self.id_col]
            for column in definition.evidence_columns[:2]:
                if column in result.df.columns and column not in preview_cols:
                    preview_cols.append(column)
            for column in ("Decision_Confidence", "Decision_Rationale"):
                if column in result.df.columns and column not in preview_cols:
                    preview_cols.append(column)

            preview_df = result.df[[c for c in preview_cols if c in result.df.columns]].head(
                min(self.top_n, 3)
            )
            if preview_df.empty:
                continue

            preview_blocks_md.append(
                "\n".join(
                    [
                        f"### Queue: {definition.name}",
                        "",
                        definition.description.strip(),
                        "",
                        _markdown_table(
                            list(preview_df.columns),
                            preview_df.values.tolist(),
                        ),
                        "",
                    ]
                )
            )
            preview_blocks_html.append(
                "\n".join(
                    [
                        f"<h3>Queue: {_html_escape(definition.name)}</h3>",
                        f"<p class='desc'>{_html_escape(definition.description.strip())}</p>",
                        _html_table(list(preview_df.columns), preview_df.values.tolist()),
                    ]
                )
            )

        md = []
        md.append("## Decision products")
        md.append("")
        md.append(
            "This section turns the current campaign table into operational queues so the output says not only what is present, "
            "but what should move next, what should be tested next, and what should stay on hold."
        )
        md.append("")
        md.append(_markdown_table(["item", "value"], evidence_rows))
        md.append("")
        md.append("Operational queues:")
        md.append("")
        md.append(
            _markdown_table(
                ["queue", "name", "count", "next_action", "top_compounds"],
                queue_rows,
            )
        )
        md.append("")
        if experiment_rows:
            md.append("Experiment prompts:")
            md.append("")
            md.append(
                _markdown_table(
                    [
                        "queue",
                        "hypothesis",
                        "minimal_test_set",
                        "minimal_analog_set",
                        "top_compounds",
                    ],
                    experiment_rows,
                )
            )
            md.append("")
        if preview_blocks_md:
            md.append("\n\n".join(preview_blocks_md))

        h = []
        h.append("<section class='category' id='decision_products'>")
        h.append("<h2>Decision products</h2>")
        h.append(
            "<p class='desc'>This section turns the current campaign table into operational queues so the output says not only what is present, "
            "but what should move next, what should be tested next, and what should stay on hold.</p>"
        )
        h.append(_html_table(["item", "value"], evidence_rows))
        h.append("<h3>Operational queues</h3>")
        h.append(
            _html_table(
                ["queue", "name", "count", "next_action", "top_compounds"],
                queue_rows,
            )
        )
        if experiment_rows:
            h.append("<h3>Experiment prompts</h3>")
            h.append(
                _html_table(
                    [
                        "queue",
                        "hypothesis",
                        "minimal_test_set",
                        "minimal_analog_set",
                        "top_compounds",
                    ],
                    experiment_rows,
                )
            )
        if preview_blocks_html:
            h.extend(preview_blocks_html)
        h.append("</section>")

        return "\n".join(md), "\n".join(h)

    def build(
        self,
        outdir: str | Path,
        write_markdown: bool = True,
        write_html: bool = True,
        make_plots: bool = True,
        title: str = "MolScope Toolkit Report",
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
                        self._plot_numeric_hist(
                            f"{spec.name}: {key}", series, plot_path
                        )
                    else:
                        self._plot_categorical_bar(
                            f"{spec.name}: {key}", series, plot_path
                        )
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
            h.append(
                f"<h2>{_html_escape(spec.name)} <code>{_html_escape(cat_key)}</code></h2>"
            )
            h.append(f"<p class='desc'>{_html_escape(spec.description.strip())}</p>")
            h.append(
                f"<p><strong>Key metric:</strong> <code>{_html_escape(key)}</code></p>"
            )
            h.append(_html_table(["metric", "value"], summary_items))
            if plot_rel:
                h.append(
                    f"<div class='plot'><img src='{_html_escape(plot_rel)}' alt='{_html_escape(spec.name)} plot'/></div>"
                )
            h.append("<h3>Top compounds (preview)</h3>")
            h.append(_html_table(top_headers, top_rows))
            h.append(
                f"<p class='meta'>Available columns in dataset for this category: {len(cols)}</p>"
            )
            h.append("</section>")
            cat_blocks_html.append("\n".join(h))

        # Optional: similarity columns (computed by molprop-report and/or molprop-picklists)
        sim_blocks_md: List[str] = []
        sim_blocks_html: List[str] = []

        sim_numeric_cols = [
            c
            for c in self.df.columns
            if str(c).startswith("Sim_") and _is_numeric(self.df[c])
        ]
        if sim_numeric_cols:
            md = []
            md.append("## Similarity to reference set")
            md.append("")
            md.append(
                "This section summarizes any similarity-to-reference columns present in the dataset (columns starting with `Sim_`). "
                "These are typically produced by `molscope report --fp ... --sim-ref-*` or `molscope picklists --fp ... --sim-ref-*`."
            )
            md.append("")

            h = []
            h.append("<section class='category' id='similarity'>")
            h.append("<h2>Similarity to reference set</h2>")
            h.append(
                "<p class='desc'>This section summarizes similarity-to-reference columns present in the dataset (columns starting with <code>Sim_</code>). "
                "These are typically produced by <code>molscope report --fp ... --sim-ref-*</code> or <code>molscope picklists --fp ... --sim-ref-*</code>."
                "</p>"
            )

            for c in sim_numeric_cols[:8]:
                s = self.df[c]
                kind, summary_items = self._summarize_key(str(c), s)

                plot_rel: Optional[str] = None
                if make_plots and kind == "numeric":
                    plot_name = f"similarity__{_slug(str(c))}.png"
                    plot_path = plots_dir / plot_name
                    try:
                        self._plot_numeric_hist(f"Similarity: {c}", s, plot_path)
                        plot_rel = f"plots/{plot_name}"
                    except Exception as e:
                        summary_items.append(("plot_error", str(e)))

                # Top-N by similarity
                df_top = self.df[
                    [self.id_col, c]
                    + ([f"{c}_BestRef"] if f"{c}_BestRef" in self.df.columns else [])
                ].copy()
                df_top[c] = pd.to_numeric(df_top[c], errors="coerce")
                df_top = df_top.sort_values(by=c, ascending=False).head(self.top_n)
                headers = list(df_top.columns)
                rows = df_top.values.tolist()

                md.append(f"### `{c}`")
                md.append("")
                md.append(_markdown_table(["metric", "value"], summary_items))
                md.append("")
                if plot_rel:
                    md.append(f"![]({plot_rel})")
                    md.append("")
                md.append("Top compounds (preview):")
                md.append("")
                md.append(_markdown_table(headers, rows))
                md.append("")

                h.append(f"<h3><code>{_html_escape(str(c))}</code></h3>")
                h.append(_html_table(["metric", "value"], summary_items))
                if plot_rel:
                    h.append(
                        f"<div class='plot'><img src='{_html_escape(plot_rel)}' alt='Similarity plot'/></div>"
                    )
                h.append("<h3>Top compounds (preview)</h3>")
                h.append(_html_table(headers, rows))

            h.append("</section>")

            sim_blocks_md.append("\n".join(md))
            sim_blocks_html.append("\n".join(h))

        retro_block_md, retro_block_html = self._retro_section()
        stage_block_md, stage_block_html = self._stage_section()
        decision_block_md, decision_block_html = self._decision_section()

        # Overview table
        overview_rows: List[List[object]] = []
        for cat_key in present:
            spec = CATEGORY_SPECS[cat_key]
            cols = self.category_columns(spec)
            key = (
                spec.key_column
                if (spec.key_column and spec.key_column in self.df.columns)
                else ""
            )
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
                "\n\n".join([stage_block_md]) if stage_block_md else "",
                "\n\n".join([decision_block_md]) if decision_block_md else "",
                "\n\n".join(sim_blocks_md) if sim_blocks_md else "",
                "\n\n".join([retro_block_md]) if retro_block_md else "",
                "---",
                "",
                (
                    "\n\n".join(cat_blocks_md)
                    if cat_blocks_md
                    else "(No categories detected in dataset)"
                ),
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
<meta charset=\"utf-8\"/>
<meta name=\"viewport\" content=\"width=device-width,initial-scale=1\"/>
<title>{_html_escape(title)}</title>
<style>{css}</style>
</head>
<body>
<h1>{_html_escape(title)}</h1>
<p class=\"sub\">Dataset: <code>{_html_escape(self.csv_path)}</code> | Generated: <code>{_html_escape(generated)}</code> | Compounds: <code>{len(self.df)}</code> | Columns: <code>{len(self.df.columns)}</code> | ID: <code>{_html_escape(self.id_col)}</code></p>
<h2>Category coverage</h2>
{html_overview}
{stage_block_html}
{decision_block_html}
{''.join(sim_blocks_html)}
{retro_block_html}
{''.join(cat_blocks_html) if cat_blocks_html else '<p>(No categories detected in dataset)</p>'}
</body>
</html>
"""

        if write_markdown:
            md_path.write_text(md_doc, encoding="utf-8")
        if write_html:
            html_path.write_text(html_doc, encoding="utf-8")

        return ReportArtifacts(
            outdir=outdir_p,
            markdown_path=md_path,
            html_path=html_path,
            plots_dir=plots_dir,
        )
