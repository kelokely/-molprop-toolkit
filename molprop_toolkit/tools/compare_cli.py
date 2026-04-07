#!/usr/bin/env python3
"""Compare two MolScope-compatible result tables.

This first compare slice is intentionally lightweight:
- detect the table ID columns,
- summarize compound overlap and numeric endpoint shifts,
- break out scaffold-level count changes when scaffold assignments exist,
- write a small bundle that later campaign views can build on.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Optional, Sequence

import pandas as pd
from pandas.api.types import is_bool_dtype, is_numeric_dtype

from molprop_toolkit.core import detect_id_column, read_table
from molprop_toolkit.core.metadata import sha256_file, write_run_metadata


def _default_output_dir(reference_path: str, candidate_path: str) -> str:
    ref = Path(reference_path)
    cand = Path(candidate_path)
    return str(cand.with_name(f"{cand.stem}_vs_{ref.stem}_compare"))


def _parse_csv_list(value: Optional[str]) -> list[str]:
    if not value:
        return []
    return [item.strip() for item in str(value).split(",") if item.strip()]


def _normalize_id_series(series: pd.Series) -> pd.Series:
    normalized = series.astype("string").str.strip()
    normalized = normalized.replace({"": pd.NA, "nan": pd.NA, "None": pd.NA})
    return normalized


def _detect_group_column(
    reference_df: pd.DataFrame,
    candidate_df: pd.DataFrame,
    requested: Optional[str],
) -> Optional[str]:
    if requested:
        if requested not in reference_df.columns or requested not in candidate_df.columns:
            raise SystemExit(
                f"molscope compare: group column '{requested}' must exist in both tables"
            )
        return requested

    for candidate in ("Scaffold_ID", "Series_ID", "Cluster_ID"):
        if candidate in reference_df.columns and candidate in candidate_df.columns:
            return candidate
    return None


def _select_numeric_columns(
    reference_df: pd.DataFrame,
    candidate_df: pd.DataFrame,
    requested: Optional[str],
    *,
    exclude: set[str],
) -> list[str]:
    if requested:
        numeric_columns = _parse_csv_list(requested)
        missing = [
            col
            for col in numeric_columns
            if col not in reference_df.columns or col not in candidate_df.columns
        ]
        if missing:
            raise SystemExit(
                "molscope compare: numeric columns not found in both tables: "
                + ", ".join(missing)
            )
        return numeric_columns

    columns: list[str] = []
    for column in reference_df.columns:
        if column in exclude or column not in candidate_df.columns:
            continue
        if is_bool_dtype(reference_df[column]) or is_bool_dtype(candidate_df[column]):
            continue
        if is_numeric_dtype(reference_df[column]) and is_numeric_dtype(
            candidate_df[column]
        ):
            columns.append(column)
    return columns


def _series_to_mean(series: pd.Series) -> Optional[float]:
    values = pd.to_numeric(series, errors="coerce").dropna()
    if values.empty:
        return None
    return round(float(values.mean()), 4)


def _series_to_median(series: pd.Series) -> Optional[float]:
    values = pd.to_numeric(series, errors="coerce").dropna()
    if values.empty:
        return None
    return round(float(values.median()), 4)


def _delta(a: Optional[float], b: Optional[float]) -> Optional[float]:
    if a is None or b is None:
        return None
    return round(b - a, 4)


def _build_endpoint_summary(
    reference_df: pd.DataFrame,
    candidate_df: pd.DataFrame,
    *,
    reference_id: str,
    candidate_id: str,
    numeric_columns: list[str],
) -> pd.DataFrame:
    ref_subset = reference_df[[reference_id, *numeric_columns]].copy()
    cand_subset = candidate_df[[candidate_id, *numeric_columns]].copy()
    ref_subset["_join_id"] = _normalize_id_series(ref_subset[reference_id])
    cand_subset["_join_id"] = _normalize_id_series(cand_subset[candidate_id])

    ref_subset = ref_subset.dropna(subset=["_join_id"])
    cand_subset = cand_subset.dropna(subset=["_join_id"])

    shared = ref_subset.merge(
        cand_subset,
        on="_join_id",
        how="inner",
        suffixes=("_reference", "_candidate"),
    )

    rows: list[dict[str, object]] = []
    for column in numeric_columns:
        ref_values = pd.to_numeric(reference_df[column], errors="coerce")
        cand_values = pd.to_numeric(candidate_df[column], errors="coerce")
        shared_ref = pd.to_numeric(shared.get(f"{column}_reference"), errors="coerce")
        shared_cand = pd.to_numeric(shared.get(f"{column}_candidate"), errors="coerce")
        shared_delta = shared_cand - shared_ref
        shared_mask = shared_ref.notna() & shared_cand.notna()

        reference_mean = _series_to_mean(ref_values)
        candidate_mean = _series_to_mean(cand_values)
        reference_median = _series_to_median(ref_values)
        candidate_median = _series_to_median(cand_values)
        paired_mean_delta = _series_to_mean(shared_delta.loc[shared_mask])
        paired_median_delta = _series_to_median(shared_delta.loc[shared_mask])

        rows.append(
            {
                "endpoint": column,
                "reference_non_null": int(ref_values.notna().sum()),
                "candidate_non_null": int(cand_values.notna().sum()),
                "shared_non_null": int(shared_mask.sum()),
                "reference_mean": reference_mean,
                "candidate_mean": candidate_mean,
                "mean_delta": _delta(reference_mean, candidate_mean),
                "reference_median": reference_median,
                "candidate_median": candidate_median,
                "median_delta": _delta(reference_median, candidate_median),
                "paired_mean_delta": paired_mean_delta,
                "paired_median_delta": paired_median_delta,
            }
        )

    return pd.DataFrame(rows)


def _build_group_summary(
    reference_df: pd.DataFrame,
    candidate_df: pd.DataFrame,
    *,
    reference_id: str,
    candidate_id: str,
    group_column: str,
) -> pd.DataFrame:
    ref_groups = (
        reference_df[[reference_id, group_column]]
        .assign(
            _join_id=lambda df: _normalize_id_series(df[reference_id]),
            _group=lambda df: _normalize_id_series(df[group_column]),
        )
        .dropna(subset=["_group"])
    )
    cand_groups = (
        candidate_df[[candidate_id, group_column]]
        .assign(
            _join_id=lambda df: _normalize_id_series(df[candidate_id]),
            _group=lambda df: _normalize_id_series(df[group_column]),
        )
        .dropna(subset=["_group"])
    )

    ref_counts = ref_groups["_group"].value_counts(dropna=False)
    cand_counts = cand_groups["_group"].value_counts(dropna=False)

    shared = ref_groups[["_join_id", "_group"]].merge(
        cand_groups[["_join_id", "_group"]],
        on="_join_id",
        how="inner",
        suffixes=("_reference", "_candidate"),
    )
    shared = shared[
        shared["_group_reference"].notna()
        & (shared["_group_reference"] == shared["_group_candidate"])
    ]
    shared_counts = shared["_group_reference"].value_counts(dropna=False)

    groups = sorted(set(ref_counts.index.tolist()) | set(cand_counts.index.tolist()))
    rows = []
    for group in groups:
        ref_count = int(ref_counts.get(group, 0))
        cand_count = int(cand_counts.get(group, 0))
        rows.append(
            {
                group_column: group,
                "reference_count": ref_count,
                "candidate_count": cand_count,
                "count_delta": cand_count - ref_count,
                "shared_compounds": int(shared_counts.get(group, 0)),
            }
        )

    return pd.DataFrame(rows)


def _write_summary_markdown(
    summary_path: Path,
    *,
    reference_name: str,
    candidate_name: str,
    summary: dict[str, object],
    endpoint_summary: pd.DataFrame,
    group_summary: Optional[pd.DataFrame],
    group_column: Optional[str],
) -> None:
    lines = [
        f"# Compare Summary: {candidate_name} vs {reference_name}",
        "",
        "## Dataset overlap",
        f"- reference rows: {summary['reference_rows']}",
        f"- candidate rows: {summary['candidate_rows']}",
        f"- shared IDs: {summary['shared_ids']}",
        f"- candidate-only IDs: {summary['candidate_only_ids']}",
        f"- reference-only IDs: {summary['reference_only_ids']}",
        "",
    ]

    if not endpoint_summary.empty:
        top_shifts = endpoint_summary.copy()
        top_shifts["abs_paired_mean_delta"] = (
            pd.to_numeric(top_shifts["paired_mean_delta"], errors="coerce").abs()
        )
        top_shifts = top_shifts.sort_values(
            by="abs_paired_mean_delta", ascending=False, na_position="last"
        ).head(5)
        lines.append("## Largest paired endpoint shifts")
        if top_shifts.empty:
            lines.append("- no numeric endpoints had paired values in both tables")
        else:
            for _, row in top_shifts.iterrows():
                lines.append(
                    f"- {row['endpoint']}: paired mean delta {row['paired_mean_delta']}"
                )
        lines.append("")

    if group_summary is not None and not group_summary.empty and group_column:
        top_groups = group_summary.reindex(
            group_summary["count_delta"].abs().sort_values(ascending=False).index
        ).head(5)
        lines.append(f"## Largest {group_column} count changes")
        for _, row in top_groups.iterrows():
            lines.append(
                f"- {row[group_column]}: {row['reference_count']} -> {row['candidate_count']} "
                f"(delta {row['count_delta']})"
            )
        lines.append("")

    summary_path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def _top_endpoint_rows(
    endpoint_summary: pd.DataFrame,
    *,
    positive: bool,
    limit: int = 3,
) -> list[dict[str, object]]:
    if endpoint_summary.empty:
        return []

    ranked = endpoint_summary.copy()
    ranked["paired_mean_delta"] = pd.to_numeric(
        ranked["paired_mean_delta"], errors="coerce"
    )
    ranked["shared_non_null"] = pd.to_numeric(
        ranked["shared_non_null"], errors="coerce"
    ).fillna(0)
    if positive:
        ranked = ranked[ranked["paired_mean_delta"] > 0]
        ranked = ranked.sort_values("paired_mean_delta", ascending=False)
    else:
        ranked = ranked[ranked["paired_mean_delta"] < 0]
        ranked = ranked.sort_values("paired_mean_delta", ascending=True)

    rows: list[dict[str, object]] = []
    for _, row in ranked.head(limit).iterrows():
        rows.append(
            {
                "endpoint": str(row["endpoint"]),
                "paired_mean_delta": round(float(row["paired_mean_delta"]), 4),
                "shared_non_null": int(row["shared_non_null"]),
            }
        )
    return rows


def _low_coverage_rows(
    endpoint_summary: pd.DataFrame,
    *,
    limit: int = 3,
) -> list[dict[str, object]]:
    if endpoint_summary.empty:
        return []

    ranked = endpoint_summary.copy()
    ranked["shared_non_null"] = pd.to_numeric(
        ranked["shared_non_null"], errors="coerce"
    ).fillna(0)
    ranked = ranked[ranked["shared_non_null"] <= 1]
    ranked = ranked.sort_values(["shared_non_null", "endpoint"], ascending=[True, True])

    rows: list[dict[str, object]] = []
    for _, row in ranked.head(limit).iterrows():
        rows.append(
            {
                "endpoint": str(row["endpoint"]),
                "shared_non_null": int(row["shared_non_null"]),
            }
        )
    return rows


def _top_group_rows(
    group_summary: Optional[pd.DataFrame],
    *,
    group_column: Optional[str],
    positive: bool,
    limit: int = 3,
) -> list[dict[str, object]]:
    if group_summary is None or group_summary.empty or not group_column:
        return []

    ranked = group_summary.copy()
    ranked["count_delta"] = pd.to_numeric(ranked["count_delta"], errors="coerce").fillna(0)
    if positive:
        ranked = ranked[ranked["count_delta"] > 0]
        ranked = ranked.sort_values("count_delta", ascending=False)
    else:
        ranked = ranked[ranked["count_delta"] < 0]
        ranked = ranked.sort_values("count_delta", ascending=True)

    rows: list[dict[str, object]] = []
    for _, row in ranked.head(limit).iterrows():
        rows.append(
            {
                "group": str(row[group_column]),
                "count_delta": int(row["count_delta"]),
                "reference_count": int(row["reference_count"]),
                "candidate_count": int(row["candidate_count"]),
            }
        )
    return rows


def _build_round_review(
    *,
    reference_name: str,
    candidate_name: str,
    endpoint_summary: pd.DataFrame,
    group_summary: Optional[pd.DataFrame],
    group_column: Optional[str],
) -> dict[str, object]:
    positive_shifts = _top_endpoint_rows(endpoint_summary, positive=True)
    negative_shifts = _top_endpoint_rows(endpoint_summary, positive=False)
    low_coverage = _low_coverage_rows(endpoint_summary)
    expanded_groups = _top_group_rows(
        group_summary, group_column=group_column, positive=True
    )
    contracted_groups = _top_group_rows(
        group_summary, group_column=group_column, positive=False
    )

    open_questions: list[str] = []
    for row in low_coverage:
        open_questions.append(
            f"Add more overlapping measurements for {row['endpoint']} before treating it as a stable round-over-round signal."
        )
    if group_column and expanded_groups:
        open_questions.append(
            f"Check whether the expanded {group_column.lower()} groups preserved property balance or only increased campaign breadth."
        )
    if group_column and contracted_groups:
        open_questions.append(
            f"Review whether contracted {group_column.lower()} groups were intentionally deprioritized or regressed under the current objective."
        )

    return {
        "reference_name": reference_name,
        "candidate_name": candidate_name,
        "group_column": group_column,
        "positive_shifts": positive_shifts,
        "negative_shifts": negative_shifts,
        "low_coverage_endpoints": low_coverage,
        "expanded_groups": expanded_groups,
        "contracted_groups": contracted_groups,
        "open_questions": open_questions,
    }


def _write_round_review_markdown(
    path: Path,
    *,
    review: dict[str, object],
) -> None:
    lines = [
        f"# Round Review: {review['candidate_name']} vs {review['reference_name']}",
        "",
    ]

    def _append_rows(title: str, rows: list[dict[str, object]], formatter) -> None:
        lines.append(f"## {title}")
        if not rows:
            lines.append("- none surfaced from the current comparison")
        else:
            for row in rows:
                lines.append(f"- {formatter(row)}")
        lines.append("")

    _append_rows(
        "Largest positive paired shifts",
        review["positive_shifts"],
        lambda row: f"{row['endpoint']}: paired mean delta {row['paired_mean_delta']} across {row['shared_non_null']} shared measurements",
    )
    _append_rows(
        "Largest negative paired shifts",
        review["negative_shifts"],
        lambda row: f"{row['endpoint']}: paired mean delta {row['paired_mean_delta']} across {row['shared_non_null']} shared measurements",
    )
    _append_rows(
        "Lowest paired coverage",
        review["low_coverage_endpoints"],
        lambda row: f"{row['endpoint']}: only {row['shared_non_null']} shared non-null measurement(s)",
    )

    group_column = review.get("group_column")
    if group_column:
        _append_rows(
            f"Largest {group_column} expansions",
            review["expanded_groups"],
            lambda row: f"{row['group']}: {row['reference_count']} -> {row['candidate_count']} (delta {row['count_delta']})",
        )
        _append_rows(
            f"Largest {group_column} contractions",
            review["contracted_groups"],
            lambda row: f"{row['group']}: {row['reference_count']} -> {row['candidate_count']} (delta {row['count_delta']})",
        )

    _append_rows(
        "Open questions",
        review["open_questions"],
        lambda row: str(row),
    )

    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(
        prog="molscope compare",
        description="Compare two MolScope-compatible results tables and summarize campaign shifts.",
    )
    ap.add_argument("reference", help="Reference table (for example: prior round)")
    ap.add_argument("candidate", help="Candidate table (for example: current round)")
    ap.add_argument(
        "-o",
        "--output-dir",
        default=None,
        help="Output directory (default: <candidate>_vs_<reference>_compare)",
    )
    ap.add_argument(
        "--reference-id",
        default=None,
        help="ID column in the reference table (default: auto-detect)",
    )
    ap.add_argument(
        "--candidate-id",
        default=None,
        help="ID column in the candidate table (default: auto-detect)",
    )
    ap.add_argument(
        "--numeric-columns",
        default=None,
        help="Comma-separated numeric columns to compare (default: common numeric columns)",
    )
    ap.add_argument(
        "--group-column",
        default=None,
        help="Optional grouping column for count deltas (default: auto-detect Scaffold_ID/Series_ID/Cluster_ID)",
    )
    ap.add_argument(
        "--reference-name",
        default=None,
        help="Human-readable label for the reference table (default: file stem)",
    )
    ap.add_argument(
        "--candidate-name",
        default=None,
        help="Human-readable label for the candidate table (default: file stem)",
    )
    args = ap.parse_args(list(argv) if argv is not None else None)

    reference_df = read_table(args.reference)
    candidate_df = read_table(args.candidate)

    reference_id = args.reference_id or detect_id_column(reference_df)
    candidate_id = args.candidate_id or detect_id_column(candidate_df)

    if not reference_id or not candidate_id:
        raise SystemExit(
            "molscope compare: could not detect an ID column in both input tables"
        )

    group_column = _detect_group_column(reference_df, candidate_df, args.group_column)
    numeric_columns = _select_numeric_columns(
        reference_df,
        candidate_df,
        args.numeric_columns,
        exclude={reference_id, candidate_id, *( [group_column] if group_column else [] )},
    )

    output_dir = Path(args.output_dir or _default_output_dir(args.reference, args.candidate))
    output_dir.mkdir(parents=True, exist_ok=True)

    reference_name = args.reference_name or Path(args.reference).stem
    candidate_name = args.candidate_name or Path(args.candidate).stem

    reference_ids = set(
        _normalize_id_series(reference_df[reference_id]).dropna().astype(str).tolist()
    )
    candidate_ids = set(
        _normalize_id_series(candidate_df[candidate_id]).dropna().astype(str).tolist()
    )
    shared_ids = reference_ids & candidate_ids

    endpoint_summary = _build_endpoint_summary(
        reference_df,
        candidate_df,
        reference_id=reference_id,
        candidate_id=candidate_id,
        numeric_columns=numeric_columns,
    )
    endpoint_path = output_dir / "endpoint_deltas.csv"
    endpoint_summary.to_csv(endpoint_path, index=False)

    group_summary: Optional[pd.DataFrame] = None
    group_path: Optional[Path] = None
    if group_column:
        group_summary = _build_group_summary(
            reference_df,
            candidate_df,
            reference_id=reference_id,
            candidate_id=candidate_id,
            group_column=group_column,
        )
        group_path = output_dir / f"{group_column.lower()}_deltas.csv"
        group_summary.to_csv(group_path, index=False)

    summary = {
        "reference_name": reference_name,
        "candidate_name": candidate_name,
        "reference_rows": int(len(reference_df)),
        "candidate_rows": int(len(candidate_df)),
        "reference_unique_ids": int(len(reference_ids)),
        "candidate_unique_ids": int(len(candidate_ids)),
        "shared_ids": int(len(shared_ids)),
        "candidate_only_ids": int(len(candidate_ids - reference_ids)),
        "reference_only_ids": int(len(reference_ids - candidate_ids)),
        "reference_id_column": reference_id,
        "candidate_id_column": candidate_id,
        "numeric_columns_compared": numeric_columns,
        "group_column": group_column,
    }

    summary_path = output_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    markdown_path = output_dir / "summary.md"
    _write_summary_markdown(
        markdown_path,
        reference_name=reference_name,
        candidate_name=candidate_name,
        summary=summary,
        endpoint_summary=endpoint_summary,
        group_summary=group_summary,
        group_column=group_column,
    )

    round_review = _build_round_review(
        reference_name=reference_name,
        candidate_name=candidate_name,
        endpoint_summary=endpoint_summary,
        group_summary=group_summary,
        group_column=group_column,
    )
    round_review_json_path = output_dir / "round_review.json"
    round_review_json_path.write_text(
        json.dumps(round_review, indent=2) + "\n", encoding="utf-8"
    )
    round_review_markdown_path = output_dir / "round_review.md"
    _write_round_review_markdown(
        round_review_markdown_path,
        review=round_review,
    )

    artifacts = {
        "endpoint_summary": str(endpoint_path.resolve()),
        "summary_markdown": str(markdown_path.resolve()),
        "round_review_json": str(round_review_json_path.resolve()),
        "round_review_markdown": str(round_review_markdown_path.resolve()),
        "reference_sha256": sha256_file(Path(args.reference)),
        "candidate_sha256": sha256_file(Path(args.candidate)),
    }
    if group_path is not None:
        artifacts["group_summary"] = str(group_path.resolve())

    write_run_metadata(
        tool="molscope-compare",
        output_table_path=endpoint_path,
        input_path=args.candidate,
        parameters={
            "reference": str(Path(args.reference).resolve()),
            "candidate": str(Path(args.candidate).resolve()),
            "reference_id": reference_id,
            "candidate_id": candidate_id,
            "group_column": group_column,
            "numeric_columns": numeric_columns,
        },
        artifacts=artifacts,
        notes="Compare bundle with dataset overlap, endpoint deltas, group-level count changes, and a round-review summary.",
    )

    print(f"[ok] wrote compare bundle to {output_dir}")


if __name__ == "__main__":
    main()
