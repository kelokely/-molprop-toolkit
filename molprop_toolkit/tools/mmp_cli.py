#!/usr/bin/env python3
"""Build a matched molecular pair bundle from a MolScope-compatible table."""

from __future__ import annotations

import argparse
import json
import re
from itertools import combinations
from pathlib import Path
from typing import Any, Optional, Sequence

import pandas as pd
from pandas.api.types import is_numeric_dtype

from molprop_toolkit.core import detect_best_smiles_column, detect_id_column, read_table
from molprop_toolkit.core.metadata import sha256_file, write_run_metadata

MMP_IMPORT_ERROR: ModuleNotFoundError | None = None

try:
    from rdkit import Chem
    from rdkit.Chem import rdMMPA
except ModuleNotFoundError as exc:
    if exc.name and exc.name.startswith("rdkit"):
        MMP_IMPORT_ERROR = exc
    else:  # pragma: no cover
        raise


PREFERRED_METRICS: tuple[str, ...] = (
    "Measured_pIC50",
    "Lead_LLE",
    "Lead_LE",
    "Lead_BEI",
    "Lead_SEI",
    "CNS_MPO",
    "QED",
    "Sol_ESOL_LogS",
    "Perm_Caco2_Class",
    "MolWt",
    "LogP",
    "TPSA",
)

ATTACHMENT_RE = re.compile(r"\[\*:\d+\]")


def _require_mmp_dependencies() -> None:
    if MMP_IMPORT_ERROR is None:
        return
    raise SystemExit(
        "molscope mmp requires RDKit. Install RDKit, preferably via conda-forge, "
        "before running this command."
    ) from MMP_IMPORT_ERROR


def _default_output_dir(input_path: str) -> str:
    p = Path(input_path)
    return str(p.with_name(f"{p.stem}_mmp"))


def _parse_csv_list(value: Optional[str]) -> list[str]:
    if not value:
        return []
    return [item.strip() for item in str(value).split(",") if item.strip()]


def _canonical_fragment_key(fragment: str) -> str:
    parts = [part for part in str(fragment).split(".") if part]
    normalized = [ATTACHMENT_RE.sub("[*]", part) for part in parts]
    return ".".join(sorted(normalized))


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
                "molscope mmp: requested metrics not found: " + ", ".join(missing)
            )
        non_numeric = [metric for metric in metrics if not is_numeric_dtype(df[metric])]
        if non_numeric:
            raise SystemExit(
                "molscope mmp: requested metrics must be numeric: "
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
                f"molscope mmp: priority metric '{requested}' was not found in the input table"
            )
        if not is_numeric_dtype(df[requested]):
            raise SystemExit(
                f"molscope mmp: priority metric '{requested}' must be numeric"
            )
        return requested

    for metric in ("Measured_pIC50", "Lead_LLE", "CNS_MPO", "QED"):
        if metric in metrics:
            return metric
    if metrics:
        return metrics[0]
    raise SystemExit("molscope mmp: no numeric metrics available for delta summaries")


def _build_fragment_records(
    df: pd.DataFrame,
    *,
    id_column: str,
    smiles_column: str,
    metrics: list[str],
    min_cuts: int,
    max_cuts: int,
    max_cut_bonds: int,
) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []

    extra_context_columns = [
        column
        for column in ("Scaffold_ID", "Scaffold_SMILES", "Cluster_ID")
        if column in df.columns
    ]

    for _, row in df.iterrows():
        smiles = row.get(smiles_column)
        if pd.isna(smiles):
            continue
        mol = Chem.MolFromSmiles(str(smiles))
        if mol is None:
            continue

        fragments = rdMMPA.FragmentMol(
            mol,
            minCuts=min_cuts,
            maxCuts=max_cuts,
            maxCutBonds=max_cut_bonds,
            resultsAsMols=False,
        )
        seen_keys: set[tuple[str, str]] = set()
        for variable_fragment, constant_environment in fragments:
            if not variable_fragment or not constant_environment:
                continue
            variable_key = _canonical_fragment_key(variable_fragment)
            constant_key = _canonical_fragment_key(constant_environment)
            if not variable_key or not constant_key:
                continue
            dedupe_key = (variable_key, constant_key)
            if dedupe_key in seen_keys:
                continue
            seen_keys.add(dedupe_key)

            record: dict[str, Any] = {
                "compound_id": row[id_column],
                "smiles": row[smiles_column],
                "variable_fragment": variable_key,
                "constant_environment": constant_key,
                "variable_fragment_raw": variable_fragment,
                "constant_environment_raw": constant_environment,
            }
            for metric in metrics:
                record[metric] = row.get(metric)
            for column in extra_context_columns:
                record[column] = row.get(column)
            records.append(record)

    return records


def _pair_rows_for_constant(
    records: list[dict[str, Any]],
    *,
    metrics: list[str],
) -> list[dict[str, Any]]:
    by_fragment: dict[str, list[dict[str, Any]]] = {}
    for record in records:
        by_fragment.setdefault(str(record["variable_fragment"]), []).append(record)

    rows: list[dict[str, Any]] = []
    for fragment_a, fragment_b in combinations(sorted(by_fragment.keys()), 2):
        left_fragment, right_fragment = sorted((fragment_a, fragment_b))
        left_records = by_fragment[left_fragment]
        right_records = by_fragment[right_fragment]

        for left in left_records:
            for right in right_records:
                pair_row: dict[str, Any] = {
                    "transform": f"{left_fragment}>>{right_fragment}",
                    "constant_environment": left["constant_environment"],
                    "left_fragment": left_fragment,
                    "right_fragment": right_fragment,
                    "left_compound_id": left["compound_id"],
                    "right_compound_id": right["compound_id"],
                    "left_smiles": left["smiles"],
                    "right_smiles": right["smiles"],
                }
                for context_column in ("Scaffold_ID", "Scaffold_SMILES", "Cluster_ID"):
                    if context_column in left:
                        pair_row[f"left_{context_column}"] = left.get(context_column)
                    if context_column in right:
                        pair_row[f"right_{context_column}"] = right.get(context_column)

                for metric in metrics:
                    left_value = pd.to_numeric([left.get(metric)], errors="coerce")[0]
                    right_value = pd.to_numeric([right.get(metric)], errors="coerce")[0]
                    pair_row[f"{metric}_left"] = left_value
                    pair_row[f"{metric}_right"] = right_value
                    pair_row[f"{metric}_delta"] = (
                        round(float(right_value - left_value), 4)
                        if pd.notna(left_value) and pd.notna(right_value)
                        else None
                    )
                rows.append(pair_row)

    return rows


def _build_pair_table(
    records: list[dict[str, Any]],
    *,
    metrics: list[str],
    max_pairs_per_constant: int,
) -> pd.DataFrame:
    by_constant: dict[str, list[dict[str, Any]]] = {}
    for record in records:
        by_constant.setdefault(str(record["constant_environment"]), []).append(record)

    rows: list[dict[str, Any]] = []
    for _constant_environment, constant_records in sorted(by_constant.items()):
        constant_rows = _pair_rows_for_constant(constant_records, metrics=metrics)
        constant_rows = sorted(
            constant_rows,
            key=lambda row: (
                str(row["transform"]),
                str(row["left_compound_id"]),
                str(row["right_compound_id"]),
            ),
        )
        rows.extend(constant_rows[:max_pairs_per_constant])

    return pd.DataFrame(rows)


def _safe_series_stat(series: pd.Series, fn: str) -> Optional[float]:
    values = pd.to_numeric(series, errors="coerce").dropna()
    if values.empty:
        return None
    if fn == "mean":
        return round(float(values.mean()), 4)
    if fn == "median":
        return round(float(values.median()), 4)
    raise ValueError(f"Unsupported stat: {fn}")


def _build_transform_summary(
    pair_table: pd.DataFrame,
    *,
    metrics: list[str],
) -> pd.DataFrame:
    if pair_table.empty:
        return pd.DataFrame()

    rows: list[dict[str, Any]] = []
    grouped = pair_table.groupby("transform", dropna=False)
    for transform, sub in grouped:
        row: dict[str, Any] = {
            "transform": transform,
            "pair_count": int(len(sub)),
            "constant_contexts": int(sub["constant_environment"].nunique()),
            "left_fragments": int(sub["left_fragment"].nunique()),
            "right_fragments": int(sub["right_fragment"].nunique()),
        }
        for metric in metrics:
            delta_column = f"{metric}_delta"
            delta_series = pd.to_numeric(sub[delta_column], errors="coerce")
            row[f"{metric}_mean_delta"] = _safe_series_stat(delta_series, "mean")
            row[f"{metric}_median_delta"] = _safe_series_stat(delta_series, "median")
            valid = delta_series.dropna()
            row[f"{metric}_positive_fraction"] = (
                round(float((valid > 0).mean()), 4) if not valid.empty else None
            )
            row[f"{metric}_negative_fraction"] = (
                round(float((valid < 0).mean()), 4) if not valid.empty else None
            )
        rows.append(row)

    summary = pd.DataFrame(rows)
    if "pair_count" in summary.columns:
        summary = summary.sort_values(
            ["pair_count", "constant_contexts", "transform"],
            ascending=[False, False, True],
            na_position="last",
        )
    return summary


def _write_summary_markdown(
    path: Path,
    *,
    input_name: str,
    summary: dict[str, Any],
    transform_summary: pd.DataFrame,
    priority_metric: str,
) -> None:
    lines = [
        f"# MMP Summary: {input_name}",
        "",
        "## Dataset scope",
        f"- rows: {summary['rows']}",
        f"- fragmented compounds: {summary['fragmented_compounds']}",
        f"- matched-pair rows: {summary['pair_rows']}",
        f"- transform summaries: {summary['transforms']}",
        f"- priority metric: {priority_metric}",
        "",
    ]

    if transform_summary.empty:
        lines.extend(["## Transform highlights", "- no matched-pair transforms were found", ""])
        path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")
        return

    mean_delta_column = f"{priority_metric}_mean_delta"
    ranked = transform_summary.dropna(subset=[mean_delta_column]).copy()

    lines.append("## Largest positive shifts")
    if ranked.empty:
        lines.append(f"- no transforms had valid deltas for {priority_metric}")
    else:
        for _, row in ranked.sort_values(mean_delta_column, ascending=False).head(5).iterrows():
            lines.append(
                f"- {row['transform']}: mean delta {row[mean_delta_column]} "
                f"across {int(row['pair_count'])} pairs"
            )

    lines.extend(["", "## Largest negative shifts"])
    if ranked.empty:
        lines.append(f"- no transforms had valid deltas for {priority_metric}")
    else:
        for _, row in ranked.sort_values(mean_delta_column, ascending=True).head(5).iterrows():
            lines.append(
                f"- {row['transform']}: mean delta {row[mean_delta_column]} "
                f"across {int(row['pair_count'])} pairs"
            )

    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def main(argv: Optional[Sequence[str]] = None) -> None:
    _require_mmp_dependencies()

    ap = argparse.ArgumentParser(
        prog="molscope mmp",
        description="Build matched molecular pair tables and transform summaries from a MolScope-compatible table.",
    )
    ap.add_argument("input", help="Results table")
    ap.add_argument(
        "-o",
        "--output-dir",
        default=None,
        help="Output directory (default: <input>_mmp)",
    )
    ap.add_argument(
        "--id-column",
        default=None,
        help="Compound ID column (default: auto-detect)",
    )
    ap.add_argument(
        "--smiles-column",
        default=None,
        help="SMILES column to fragment (default: auto-detect)",
    )
    ap.add_argument(
        "--metrics",
        default=None,
        help="Comma-separated numeric metrics to summarize as pair deltas",
    )
    ap.add_argument(
        "--priority-metric",
        default=None,
        help="Metric used for markdown transform ranking (default: auto-select)",
    )
    ap.add_argument(
        "--min-cuts",
        type=int,
        default=2,
        help="Minimum number of cuts for RDKit MMP fragmentation (default: 2)",
    )
    ap.add_argument(
        "--max-cuts",
        type=int,
        default=2,
        help="Maximum number of cuts for RDKit MMP fragmentation (default: 2)",
    )
    ap.add_argument(
        "--max-cut-bonds",
        type=int,
        default=20,
        help="Maximum candidate cut bonds per molecule (default: 20)",
    )
    ap.add_argument(
        "--max-pairs-per-constant",
        type=int,
        default=500,
        help="Cap pair rows retained per constant environment (default: 500)",
    )
    args = ap.parse_args(list(argv) if argv is not None else None)

    df = read_table(args.input)
    id_column = args.id_column or detect_id_column(df)
    smiles_column = args.smiles_column or detect_best_smiles_column(df.columns)

    if not id_column:
        raise SystemExit("molscope mmp: could not detect a compound ID column")
    if not smiles_column or smiles_column not in df.columns:
        raise SystemExit("molscope mmp: could not detect a valid SMILES column")

    metrics = _select_metrics(
        df,
        requested=args.metrics,
        exclude={id_column, smiles_column},
    )
    priority_metric = _choose_priority_metric(
        df,
        metrics=metrics,
        requested=args.priority_metric,
    )

    output_dir = Path(args.output_dir or _default_output_dir(args.input))
    output_dir.mkdir(parents=True, exist_ok=True)

    fragment_records = _build_fragment_records(
        df,
        id_column=id_column,
        smiles_column=smiles_column,
        metrics=metrics,
        min_cuts=int(args.min_cuts),
        max_cuts=int(args.max_cuts),
        max_cut_bonds=int(args.max_cut_bonds),
    )
    pair_table = _build_pair_table(
        fragment_records,
        metrics=metrics,
        max_pairs_per_constant=max(int(args.max_pairs_per_constant), 1),
    )
    transform_summary = _build_transform_summary(pair_table, metrics=metrics)

    pair_table_path = output_dir / "pair_table.csv"
    transform_summary_path = output_dir / "transform_summary.csv"
    summary_json_path = output_dir / "summary.json"
    summary_md_path = output_dir / "summary.md"

    pair_table.to_csv(pair_table_path, index=False)
    transform_summary.to_csv(transform_summary_path, index=False)

    summary = {
        "input_name": Path(args.input).stem,
        "rows": int(len(df)),
        "fragment_records": int(len(fragment_records)),
        "fragmented_compounds": int(
            len({record["compound_id"] for record in fragment_records})
        ),
        "pair_rows": int(len(pair_table)),
        "transforms": int(transform_summary["transform"].nunique())
        if not transform_summary.empty
        else 0,
        "id_column": id_column,
        "smiles_column": smiles_column,
        "metrics": metrics,
        "priority_metric": priority_metric,
        "min_cuts": int(args.min_cuts),
        "max_cuts": int(args.max_cuts),
        "max_cut_bonds": int(args.max_cut_bonds),
    }
    summary_json_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    _write_summary_markdown(
        summary_md_path,
        input_name=Path(args.input).stem,
        summary=summary,
        transform_summary=transform_summary,
        priority_metric=priority_metric,
    )

    write_run_metadata(
        tool="molscope-mmp",
        output_table_path=pair_table_path,
        input_path=args.input,
        smiles_column=smiles_column,
        parameters={
            "id_column": id_column,
            "smiles_column": smiles_column,
            "metrics": metrics,
            "priority_metric": priority_metric,
            "min_cuts": int(args.min_cuts),
            "max_cuts": int(args.max_cuts),
            "max_cut_bonds": int(args.max_cut_bonds),
            "max_pairs_per_constant": int(max(args.max_pairs_per_constant, 1)),
        },
        artifacts={
            "transform_summary": str(transform_summary_path.resolve()),
            "summary_json": str(summary_json_path.resolve()),
            "summary_markdown": str(summary_md_path.resolve()),
            "input_sha256": sha256_file(Path(args.input)),
        },
        notes="First MMP slice: RDKit fragmentation-based pair rows and aggregate transform summaries.",
    )

    print(f"[ok] wrote MMP bundle to {output_dir}")


if __name__ == "__main__":
    main()
