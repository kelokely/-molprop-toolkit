#!/usr/bin/env python3
"""Build a scaffold-centric SAR review bundle from a MolScope results table."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Optional, Sequence

import pandas as pd
from pandas.api.types import is_numeric_dtype

from molprop_toolkit.core import detect_id_column, read_table
from molprop_toolkit.core.metadata import sha256_file, write_run_metadata

PREFERRED_METRICS: tuple[str, ...] = (
    "Measured_pIC50",
    "Lead_LLE",
    "Lead_LE",
    "Lead_BEI",
    "Lead_SEI",
    "CNS_MPO",
    "QED",
    "Sol_ESOL_LogS",
    "MolWt",
    "LogP",
    "TPSA",
)


def _default_output_dir(input_path: str) -> str:
    p = Path(input_path)
    return str(p.with_name(f"{p.stem}_sar"))


def _parse_csv_list(value: Optional[str]) -> list[str]:
    if not value:
        return []
    return [item.strip() for item in str(value).split(",") if item.strip()]


def _normalize_group_series(series: pd.Series) -> pd.Series:
    normalized = series.astype("string").str.strip()
    normalized = normalized.replace({"": pd.NA, "nan": pd.NA, "None": pd.NA})
    return normalized


def _detect_group_column(df: pd.DataFrame, requested: Optional[str]) -> str:
    if requested:
        if requested not in df.columns:
            raise SystemExit(
                f"molscope sar: group column '{requested}' was not found in the input table"
            )
        return requested

    for candidate in ("Scaffold_ID", "Series_ID", "Cluster_ID"):
        if candidate in df.columns:
            return candidate

    raise SystemExit(
        "molscope sar: no series grouping column found. Run 'molscope series' first "
        "or pass --group-column."
    )


def _select_metrics(
    df: pd.DataFrame,
    *,
    requested: Optional[str],
    exclude: set[str],
) -> list[str]:
    if requested:
        metrics = _parse_csv_list(requested)
        missing = [metric for metric in metrics if metric not in df.columns]
        if missing:
            raise SystemExit(
                "molscope sar: requested metrics not found: " + ", ".join(missing)
            )
        non_numeric = [metric for metric in metrics if not is_numeric_dtype(df[metric])]
        if non_numeric:
            raise SystemExit(
                "molscope sar: requested metrics must be numeric: "
                + ", ".join(non_numeric)
            )
        return metrics

    metrics: list[str] = []
    for metric in PREFERRED_METRICS:
        if metric in df.columns and metric not in exclude and is_numeric_dtype(df[metric]):
            metrics.append(metric)

    if metrics:
        return metrics

    for column in df.columns:
        if column in exclude:
            continue
        if is_numeric_dtype(df[column]):
            metrics.append(column)
    return metrics


def _choose_priority_metric(
    df: pd.DataFrame,
    *,
    metrics: list[str],
    requested: Optional[str],
) -> str:
    if requested:
        if requested not in df.columns:
            raise SystemExit(
                f"molscope sar: priority metric '{requested}' was not found in the input table"
            )
        if not is_numeric_dtype(df[requested]):
            raise SystemExit(
                f"molscope sar: priority metric '{requested}' must be numeric"
            )
        return requested

    for metric in ("Measured_pIC50", "Lead_LLE", "CNS_MPO", "QED"):
        if metric in metrics:
            return metric
    if metrics:
        return metrics[0]

    raise SystemExit("molscope sar: no numeric metrics available for SAR review")


def _series_examples(series: pd.Series, *, limit: int = 5) -> str:
    values = series.dropna().astype(str).head(limit).tolist()
    return ", ".join(values)


def _safe_stat(series: pd.Series, fn: str) -> Optional[float]:
    values = pd.to_numeric(series, errors="coerce").dropna()
    if values.empty:
        return None
    if fn == "mean":
        return round(float(values.mean()), 4)
    if fn == "median":
        return round(float(values.median()), 4)
    if fn == "min":
        return round(float(values.min()), 4)
    if fn == "max":
        return round(float(values.max()), 4)
    raise ValueError(f"Unsupported stat function: {fn}")


def _build_series_summary(
    df: pd.DataFrame,
    *,
    id_column: str,
    group_column: str,
    metrics: list[str],
    priority_metric: str,
    top_n: int,
    sort_ascending: bool,
) -> pd.DataFrame:
    grouped = (
        df.assign(_group=_normalize_group_series(df[group_column]))
        .dropna(subset=["_group"])
        .groupby("_group", dropna=False)
    )

    rows: list[dict[str, object]] = []
    for group_value, sub in grouped:
        sub_sorted = sub.sort_values(priority_metric, ascending=sort_ascending)
        best_row = sub_sorted.iloc[0] if len(sub_sorted) else None

        row: dict[str, object] = {
            group_column: group_value,
            "Count": int(len(sub)),
            "Example_Compounds": _series_examples(sub[id_column], limit=top_n),
            f"{priority_metric}_Top_Compound": (
                str(best_row[id_column]) if best_row is not None else None
            ),
            f"{priority_metric}_Top_Value": (
                round(float(best_row[priority_metric]), 4)
                if best_row is not None and pd.notna(best_row[priority_metric])
                else None
            ),
        }

        if "Scaffold_SMILES" in sub.columns:
            scaffold_smiles = (
                sub["Scaffold_SMILES"].dropna().astype(str).iloc[0]
                if sub["Scaffold_SMILES"].dropna().any()
                else ""
            )
            row["Scaffold_SMILES"] = scaffold_smiles or None

        for metric in metrics:
            row[f"{metric}_mean"] = _safe_stat(sub[metric], "mean")
            row[f"{metric}_median"] = _safe_stat(sub[metric], "median")
            row[f"{metric}_min"] = _safe_stat(sub[metric], "min")
            row[f"{metric}_max"] = _safe_stat(sub[metric], "max")

        rows.append(row)

    summary = pd.DataFrame(rows)
    if summary.empty:
        return summary

    sort_columns = ["Count"]
    ascending = [False]
    top_metric_column = f"{priority_metric}_Top_Value"
    if top_metric_column in summary.columns:
        sort_columns.append(top_metric_column)
        ascending.append(sort_ascending)
    return summary.sort_values(sort_columns, ascending=ascending, na_position="last")


def _build_top_compounds(
    df: pd.DataFrame,
    *,
    id_column: str,
    group_column: str,
    priority_metric: str,
    top_n: int,
    sort_ascending: bool,
    metrics: list[str],
) -> pd.DataFrame:
    working = df.assign(_group=_normalize_group_series(df[group_column])).dropna(
        subset=["_group"]
    )
    rows: list[dict[str, object]] = []

    extra_columns = [
        column
        for column in ("SMILES", "Calc_Canonical_SMILES", "Scaffold_SMILES")
        if column in working.columns
    ]

    for group_value, sub in working.groupby("_group", dropna=False):
        ranked = sub.sort_values(priority_metric, ascending=sort_ascending).head(top_n)
        for rank, (_, row) in enumerate(ranked.iterrows(), start=1):
            out_row: dict[str, object] = {
                group_column: group_value,
                "Rank": rank,
                id_column: row[id_column],
                priority_metric: row[priority_metric],
            }
            for metric in metrics:
                if metric == priority_metric:
                    continue
                out_row[metric] = row[metric] if metric in row else None
            for column in extra_columns:
                out_row[column] = row[column]
            rows.append(out_row)

    return pd.DataFrame(rows)


def _build_metric_summary(
    series_summary: pd.DataFrame,
    *,
    group_column: str,
    metrics: list[str],
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    if series_summary.empty:
        return pd.DataFrame(rows)

    for _, row in series_summary.iterrows():
        for metric in metrics:
            rows.append(
                {
                    group_column: row[group_column],
                    "metric": metric,
                    "mean": row.get(f"{metric}_mean"),
                    "median": row.get(f"{metric}_median"),
                    "min": row.get(f"{metric}_min"),
                    "max": row.get(f"{metric}_max"),
                }
            )
    return pd.DataFrame(rows)


def _write_summary_markdown(
    path: Path,
    *,
    input_name: str,
    group_column: str,
    priority_metric: str,
    summary: dict[str, object],
    series_summary: pd.DataFrame,
) -> None:
    lines = [
        f"# SAR Summary: {input_name}",
        "",
        "## Dataset scope",
        f"- rows: {summary['rows']}",
        f"- groups: {summary['groups']}",
        f"- grouping column: {group_column}",
        f"- priority metric: {priority_metric}",
        "",
        "## Largest groups",
    ]

    top_groups = series_summary.head(5)
    if top_groups.empty:
        lines.append("- no valid groups were found in the input table")
    else:
        top_value_column = f"{priority_metric}_Top_Value"
        top_compound_column = f"{priority_metric}_Top_Compound"
        for _, row in top_groups.iterrows():
            lines.append(
                f"- {row[group_column]}: {int(row['Count'])} compounds; top "
                f"{priority_metric} {row.get(top_value_column)} ({row.get(top_compound_column)})"
            )

    lines.extend(["", "## Suggested use", "- Review the largest groups first, then inspect the top compounds table to see which analogs are carrying the current series signal.", ""])
    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(
        prog="molscope sar",
        description="Build a scaffold-centric SAR review bundle from a MolScope-compatible table.",
    )
    ap.add_argument("input", help="Results table (ideally series-enriched)")
    ap.add_argument(
        "-o",
        "--output-dir",
        default=None,
        help="Output directory (default: <input>_sar)",
    )
    ap.add_argument(
        "--id-column",
        default=None,
        help="Compound ID column (default: auto-detect)",
    )
    ap.add_argument(
        "--group-column",
        default=None,
        help="Grouping column for SAR review (default: Scaffold_ID, then Series_ID or Cluster_ID)",
    )
    ap.add_argument(
        "--metrics",
        default=None,
        help="Comma-separated numeric metrics to summarize (default: common medchem metrics present in the table)",
    )
    ap.add_argument(
        "--priority-metric",
        default=None,
        help="Metric used to rank top compounds within each group (default: auto-select)",
    )
    ap.add_argument(
        "--top-compounds-per-group",
        type=int,
        default=5,
        help="Number of top compounds to retain per group (default: 5)",
    )
    ap.add_argument(
        "--priority-ascending",
        action="store_true",
        help="Treat lower values of the priority metric as better",
    )
    args = ap.parse_args(list(argv) if argv is not None else None)

    df = read_table(args.input)
    id_column = args.id_column or detect_id_column(df)
    if not id_column:
        raise SystemExit("molscope sar: could not detect a compound ID column")

    group_column = _detect_group_column(df, args.group_column)
    metrics = _select_metrics(df, requested=args.metrics, exclude={id_column, group_column})
    priority_metric = _choose_priority_metric(
        df,
        metrics=metrics,
        requested=args.priority_metric,
    )
    output_dir = Path(args.output_dir or _default_output_dir(args.input))
    output_dir.mkdir(parents=True, exist_ok=True)

    series_summary = _build_series_summary(
        df,
        id_column=id_column,
        group_column=group_column,
        metrics=metrics,
        priority_metric=priority_metric,
        top_n=max(args.top_compounds_per_group, 1),
        sort_ascending=bool(args.priority_ascending),
    )
    top_compounds = _build_top_compounds(
        df,
        id_column=id_column,
        group_column=group_column,
        priority_metric=priority_metric,
        top_n=max(args.top_compounds_per_group, 1),
        sort_ascending=bool(args.priority_ascending),
        metrics=metrics,
    )
    metric_summary = _build_metric_summary(
        series_summary,
        group_column=group_column,
        metrics=metrics,
    )

    series_summary_path = output_dir / "series_summary.csv"
    top_compounds_path = output_dir / "top_compounds.csv"
    metric_summary_path = output_dir / "metric_summary.csv"
    summary_json_path = output_dir / "summary.json"
    summary_md_path = output_dir / "summary.md"

    series_summary.to_csv(series_summary_path, index=False)
    top_compounds.to_csv(top_compounds_path, index=False)
    metric_summary.to_csv(metric_summary_path, index=False)

    summary = {
        "input_name": Path(args.input).stem,
        "rows": int(len(df)),
        "groups": int(series_summary[group_column].nunique()) if not series_summary.empty else 0,
        "id_column": id_column,
        "group_column": group_column,
        "metrics": metrics,
        "priority_metric": priority_metric,
        "priority_ascending": bool(args.priority_ascending),
        "top_compounds_per_group": int(max(args.top_compounds_per_group, 1)),
    }
    summary_json_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    _write_summary_markdown(
        summary_md_path,
        input_name=Path(args.input).stem,
        group_column=group_column,
        priority_metric=priority_metric,
        summary=summary,
        series_summary=series_summary,
    )

    write_run_metadata(
        tool="molscope-sar",
        output_table_path=series_summary_path,
        input_path=args.input,
        parameters={
            "id_column": id_column,
            "group_column": group_column,
            "metrics": metrics,
            "priority_metric": priority_metric,
            "priority_ascending": bool(args.priority_ascending),
            "top_compounds_per_group": int(max(args.top_compounds_per_group, 1)),
        },
        artifacts={
            "top_compounds": str(top_compounds_path.resolve()),
            "metric_summary": str(metric_summary_path.resolve()),
            "summary_json": str(summary_json_path.resolve()),
            "summary_markdown": str(summary_md_path.resolve()),
            "input_sha256": sha256_file(Path(args.input)),
        },
        notes="First SAR slice: scaffold-centric series summary and top-compound review bundle.",
    )

    print(f"[ok] wrote SAR bundle to {output_dir}")


if __name__ == "__main__":
    main()
