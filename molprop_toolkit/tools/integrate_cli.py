#!/usr/bin/env python3
"""Merge measured assay data into a MolScope results table.

This command keeps the existing table-first workflow intact:
- read a MolScope-compatible results table,
- left-join measured assay columns,
- stamp simple provenance fields for imported values,
- optionally recompute canonical lead metrics from measured potency.

The output stays compatible with the existing report and picklist tools because
it preserves the original results table shape and only adds or updates columns.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Any, Dict, Optional, Sequence

import numpy as np
import pandas as pd

from molprop_toolkit.core import detect_id_column, read_table, write_table
from molprop_toolkit.core.metadata import sha256_file, write_run_metadata


LEAD_METRIC_COLUMNS: tuple[str, ...] = (
    "Lead_LE",
    "Lead_LLE",
    "Lead_LELP",
    "Lead_BEI",
    "Lead_SEI",
    "Lead_LE_Quality",
    "Lead_LLE_Quality",
    "Lead_LELP_Quality",
)


def _parse_csv_list(value: Optional[str]) -> list[str]:
    if not value:
        return []
    return [item.strip() for item in str(value).split(",") if item.strip()]


def _default_output_path(results_path: str) -> str:
    p = Path(results_path)
    suffix = p.suffix or ".csv"
    return str(p.with_name(f"{p.stem}_integrated{suffix}"))


def _prefixed_name(column: str, prefix: str) -> str:
    if not prefix:
        return column
    marker = f"{prefix}_"
    if column.startswith(marker):
        return column
    return f"{marker}{column}"


def _select_assay_columns(
    assay_df: pd.DataFrame,
    *,
    right_id: str,
    assay_columns_arg: Optional[str],
    potency_column: Optional[str],
) -> list[str]:
    cols = _parse_csv_list(assay_columns_arg)
    if not cols:
        cols = [c for c in assay_df.columns if c != right_id]

    if potency_column and potency_column != right_id and potency_column not in cols:
        cols.append(potency_column)

    missing = [c for c in cols if c not in assay_df.columns]
    if missing:
        raise SystemExit(
            "molscope integrate: assay columns not found: " + ", ".join(missing)
        )
    return cols


def _normalize_potency_series(series: pd.Series, scale: str) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce")
    values = values.where(values > 0)

    if scale == "pic50":
        return values
    if scale == "ic50_nm":
        return 9.0 - np.log10(values)
    if scale == "ic50_um":
        return 6.0 - np.log10(values)
    if scale == "ic50_m":
        return -np.log10(values)

    raise SystemExit(f"molscope integrate: unsupported potency scale: {scale}")


def _ensure_provenance(
    df: pd.DataFrame,
    *,
    column: str,
    measured_mask: Optional[pd.Series] = None,
    heuristic_mask: Optional[pd.Series] = None,
) -> None:
    provenance_column = f"{column}_Provenance"
    if provenance_column not in df.columns:
        df[provenance_column] = pd.Series(pd.NA, index=df.index, dtype="object")

    if heuristic_mask is not None:
        df.loc[heuristic_mask & df[column].notna(), provenance_column] = "heuristic"
    if measured_mask is not None:
        df.loc[measured_mask & df[column].notna(), provenance_column] = "measured"


def _recompute_lead_metrics(
    df: pd.DataFrame,
    *,
    pic50_column: str,
) -> int:
    if pic50_column not in df.columns:
        return 0

    pic50 = pd.to_numeric(df[pic50_column], errors="coerce")
    heavy_atoms = pd.to_numeric(df.get("HeavyAtoms"), errors="coerce")
    logp = pd.to_numeric(df.get("LogP"), errors="coerce")
    mw = pd.to_numeric(df.get("MolWt"), errors="coerce")
    tpsa = pd.to_numeric(df.get("TPSA"), errors="coerce")

    if not any(
        col in df.columns for col in ("HeavyAtoms", "LogP", "MolWt", "TPSA")
    ):
        return 0

    for col in LEAD_METRIC_COLUMNS:
        if col in df.columns and f"Heuristic_{col}" not in df.columns:
            df[f"Heuristic_{col}"] = df[col]

    le_mask = pic50.notna() & heavy_atoms.notna()
    lle_mask = pic50.notna() & logp.notna()
    bei_mask = pic50.notna() & mw.notna() & (mw > 0)
    sei_mask = pic50.notna() & tpsa.notna() & (tpsa > 0)

    if le_mask.any():
        heavy = heavy_atoms.where(heavy_atoms > 0, 1.0)
        le = (pic50 / heavy).round(3)
        df.loc[le_mask, "Lead_LE"] = le.loc[le_mask]
        df.loc[le_mask, "Lead_LE_Quality"] = np.where(
            le.loc[le_mask] >= 0.3, "Good", "Poor"
        )

    if lle_mask.any():
        lle = (pic50 - logp).round(3)
        df.loc[lle_mask, "Lead_LLE"] = lle.loc[lle_mask]
        df.loc[lle_mask, "Lead_LLE_Quality"] = np.where(
            lle.loc[lle_mask] >= 7,
            "Excellent",
            np.where(lle.loc[lle_mask] >= 5, "Good", "Poor"),
        )

    if le_mask.any() and logp is not None:
        le_vals = pd.to_numeric(df.get("Lead_LE"), errors="coerce")
        lelp_mask = le_vals.notna() & (le_vals > 0) & logp.notna() & pic50.notna()
        if lelp_mask.any():
            lelp = (logp / le_vals).round(3)
            df.loc[lelp_mask, "Lead_LELP"] = lelp.loc[lelp_mask]
            df.loc[lelp_mask, "Lead_LELP_Quality"] = np.where(
                lelp.loc[lelp_mask] < 10, "Good", "Poor"
            )
    else:
        lelp_mask = pd.Series(False, index=df.index)

    if bei_mask.any():
        bei = (pic50 / (mw / 1000.0)).round(3)
        df.loc[bei_mask, "Lead_BEI"] = bei.loc[bei_mask]

    if sei_mask.any():
        sei = (pic50 / (tpsa / 100.0)).round(3)
        df.loc[sei_mask, "Lead_SEI"] = sei.loc[sei_mask]

    heuristic_rows = pic50.isna()
    for col in (
        "Lead_LE",
        "Lead_LLE",
        "Lead_LELP",
        "Lead_BEI",
        "Lead_SEI",
        "Lead_LE_Quality",
        "Lead_LLE_Quality",
        "Lead_LELP_Quality",
    ):
        if col in df.columns:
            measured_mask = pd.Series(False, index=df.index)
            if col in {"Lead_LE", "Lead_LE_Quality"}:
                measured_mask = le_mask
            elif col in {"Lead_LLE", "Lead_LLE_Quality"}:
                measured_mask = lle_mask
            elif col in {"Lead_LELP", "Lead_LELP_Quality"}:
                measured_mask = lelp_mask
            elif col == "Lead_BEI":
                measured_mask = bei_mask
            elif col == "Lead_SEI":
                measured_mask = sei_mask

            _ensure_provenance(
                df,
                column=col,
                measured_mask=measured_mask,
                heuristic_mask=heuristic_rows,
            )

    return int(pic50.notna().sum())


def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(
        prog="molscope integrate",
        description="Merge measured assay data into a MolScope-compatible results table.",
    )
    ap.add_argument("results", help="MolScope results table (CSV/TSV/Parquet)")
    ap.add_argument("assays", help="Measured assay table (CSV/TSV/Parquet)")
    ap.add_argument(
        "-o",
        "--output",
        default=None,
        help="Output table path (default: <results>_integrated.<ext>)",
    )
    ap.add_argument(
        "--left-id",
        default=None,
        help="ID column in the results table (default: auto-detect)",
    )
    ap.add_argument(
        "--right-id",
        default=None,
        help="ID column in the assay table (default: auto-detect)",
    )
    ap.add_argument(
        "--assay-columns",
        default=None,
        help="Comma-separated assay columns to import (default: all assay columns except the join key)",
    )
    ap.add_argument(
        "--prefix",
        default="Measured",
        help="Prefix for imported assay columns (default: Measured)",
    )
    ap.add_argument(
        "--potency-column",
        default=None,
        help="Assay-table potency column to normalize and use for lead-metric recomputation",
    )
    ap.add_argument(
        "--potency-scale",
        choices=["pic50", "ic50_nm", "ic50_um", "ic50_m"],
        default="pic50",
        help="Scale of --potency-column values (default: pic50)",
    )
    ap.add_argument(
        "--normalized-potency-column",
        default="Measured_pIC50",
        help="Canonical normalized potency column to write (default: Measured_pIC50)",
    )
    args = ap.parse_args(list(argv) if argv is not None else None)

    output_path = args.output or _default_output_path(args.results)

    results_df = read_table(args.results)
    assay_df = read_table(args.assays)

    left_id = args.left_id or detect_id_column(results_df)
    right_id = args.right_id or detect_id_column(assay_df)

    if left_id not in results_df.columns:
        raise SystemExit(f"molscope integrate: results ID column not found: {left_id}")
    if right_id not in assay_df.columns:
        raise SystemExit(f"molscope integrate: assay ID column not found: {right_id}")

    duplicate_mask = assay_df[right_id].duplicated(keep=False)
    if duplicate_mask.any():
        dupes = sorted(
            str(v) for v in assay_df.loc[duplicate_mask, right_id].dropna().unique()
        )
        preview = ", ".join(dupes[:5])
        raise SystemExit(
            "molscope integrate: assay table has duplicate join keys. "
            f"Please deduplicate before merging. Example IDs: {preview}"
        )

    assay_columns = _select_assay_columns(
        assay_df,
        right_id=right_id,
        assay_columns_arg=args.assay_columns,
        potency_column=args.potency_column,
    )

    renamed_columns: dict[str, str] = {
        col: _prefixed_name(col, args.prefix) for col in assay_columns
    }
    collisions = sorted(
        new_name
        for new_name in renamed_columns.values()
        if new_name in results_df.columns and new_name != left_id
    )
    if collisions:
        raise SystemExit(
            "molscope integrate: imported assay columns would collide with existing result columns: "
            + ", ".join(collisions)
        )

    assay_payload = assay_df[[right_id, *assay_columns]].copy()
    assay_payload = assay_payload.rename(
        columns={right_id: "_AssayJoinKey", **renamed_columns}
    )

    merged = results_df.merge(
        assay_payload,
        how="left",
        left_on=left_id,
        right_on="_AssayJoinKey",
        indicator="_AssayMerge",
    )
    if "_AssayJoinKey" in merged.columns:
        merged = merged.drop(columns=["_AssayJoinKey"])

    matched_mask = merged["_AssayMerge"].eq("both")
    merged["Measured_Record_Matched"] = matched_mask
    merged = merged.drop(columns=["_AssayMerge"])

    for imported_col in renamed_columns.values():
        _ensure_provenance(
            merged,
            column=imported_col,
            measured_mask=merged[imported_col].notna(),
        )

    recomputed_rows = 0
    normalized_potency_column = args.normalized_potency_column
    if args.potency_column:
        imported_potency_column = renamed_columns[args.potency_column]
        normalized = _normalize_potency_series(
            merged[imported_potency_column], args.potency_scale
        ).round(3)

        merged[normalized_potency_column] = normalized
        _ensure_provenance(
            merged,
            column=normalized_potency_column,
            measured_mask=normalized.notna(),
        )
        recomputed_rows = _recompute_lead_metrics(
            merged, pic50_column=normalized_potency_column
        )

    write_table(merged, output_path)

    artifacts: Dict[str, Any] = {
        "assay_input": {
            "path": str(Path(args.assays).resolve()),
            "sha256": sha256_file(Path(args.assays)),
        },
        "matched_rows": int(matched_mask.sum()),
        "unmatched_rows": int((~matched_mask).sum()),
        "imported_columns": list(renamed_columns.values()),
        "lead_metrics_recomputed_rows": recomputed_rows,
    }
    write_run_metadata(
        tool="molscope-integrate",
        input_path=args.results,
        output_table_path=output_path,
        parameters={
            "left_id": left_id,
            "right_id": right_id,
            "assay_columns": assay_columns,
            "prefix": args.prefix,
            "potency_column": args.potency_column,
            "potency_scale": args.potency_scale,
            "normalized_potency_column": (
                normalized_potency_column if args.potency_column else None
            ),
        },
        artifacts=artifacts,
        notes="Measured assay values are left-joined onto the canonical MolScope results table.",
    )

    print(f"Wrote integrated table: {output_path}")
    print(f"Matched rows: {int(matched_mask.sum())} / {len(merged)}")
    if args.potency_column:
        print(
            f"Normalized potency column: {normalized_potency_column} ({args.potency_scale})"
        )
        print(f"Lead-metric rows updated from measured potency: {recomputed_rows}")
