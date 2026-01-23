#!/usr/bin/env python3
"""Validate a results CSV against a MolProp Toolkit schema (Task 11A).

This tool is designed for pragmatic validation:
- Results tables are often partial (not every calculator computes every column).
- Unknown columns are allowed by default (warn-only).
- Regex patterns in the schema provide metadata for prefixed column families.

Checks performed:
- Required columns missing
- Optional requirement of registry categories (via --require-category)
- Type coercion sanity (numeric/bool/string)
- Unknown columns (warn or fail)

Example:
  python tools/validate_csv_schema.py results.csv
  python tools/validate_csv_schema.py results.csv --fail-unknown
  python tools/validate_csv_schema.py results.csv --require-category cns_mpo --require-category toxicity
"""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

# Allow running directly without installing the package.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd

from molprop_toolkit.core import read_csv
from molprop_toolkit.core.registry import CATEGORY_SPECS
from molprop_toolkit.schema import load_schema


@dataclass(frozen=True)
class ColMeta:
    name: str
    type: str
    required: bool
    units: Optional[str]
    description: str
    produced_by: Sequence[str]
    tags: Sequence[str]
    source: str  # explicit|pattern


def _load_schema(path: str | Path | None) -> Dict[str, Any]:
    # Delegate to the shared loader so installed packages can use bundled schema.
    return load_schema(path)


def _compile_patterns(schema: Dict[str, Any]) -> List[Tuple[re.Pattern[str], Dict[str, Any]]]:
    out: List[Tuple[re.Pattern[str], Dict[str, Any]]] = []
    for pat in schema.get("patterns", []) or []:
        rx = re.compile(str(pat.get("match")))
        out.append((rx, pat))
    return out


def _meta_for(col: str, schema: Dict[str, Any], patterns: List[Tuple[re.Pattern[str], Dict[str, Any]]]) -> Optional[ColMeta]:
    cols = schema.get("columns", {}) or {}
    if col in cols:
        m = cols[col]
        return ColMeta(
            name=col,
            type=str(m.get("type", "auto")),
            required=bool(m.get("required", False)),
            units=m.get("units"),
            description=str(m.get("description", "")),
            produced_by=list(m.get("produced_by", []) or []),
            tags=list(m.get("tags", []) or []),
            source="explicit",
        )

    for rx, pat in patterns:
        if rx.search(col):
            return ColMeta(
                name=col,
                type=str(pat.get("type", "auto")),
                required=False,
                units=pat.get("units"),
                description=str(pat.get("description", "")),
                produced_by=list(pat.get("produced_by", []) or []),
                tags=list(pat.get("tags", []) or []),
                source="pattern",
            )

    return None


def _coerce_bool(s: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(s):
        return s
    # tolerate common encodings
    truthy = {"true", "t", "1", "yes", "y"}
    falsy = {"false", "f", "0", "no", "n"}

    def conv(v: Any) -> Optional[bool]:
        if v is None or (isinstance(v, float) and pd.isna(v)):
            return None
        if isinstance(v, bool):
            return v
        if isinstance(v, (int, float)):
            if v == 1:
                return True
            if v == 0:
                return False
        vv = str(v).strip().lower()
        if vv in truthy:
            return True
        if vv in falsy:
            return False
        return None

    return s.map(conv)


def _coercion_ok_ratio(meta_type: str, s: pd.Series) -> float:
    nonnull = s[~s.isna()]
    if len(nonnull) == 0:
        return 1.0

    t = meta_type
    if t in ("auto", "any"):
        return 1.0

    if t == "string":
        return 1.0

    if t == "float":
        sn = pd.to_numeric(nonnull, errors="coerce")
        return float((~sn.isna()).mean())

    if t == "int":
        sn = pd.to_numeric(nonnull, errors="coerce")
        ok = (~sn.isna()) & (sn.round(0) == sn)
        return float(ok.mean())

    if t == "bool":
        sb = _coerce_bool(nonnull)
        return float((~sb.isna()).mean())

    return 1.0


def _category_columns_present(df: pd.DataFrame, category_key: str) -> List[str]:
    if category_key not in CATEGORY_SPECS:
        raise ValueError(f"Unknown registry category: {category_key}")

    spec = CATEGORY_SPECS[category_key]
    cols: List[str] = [c for c in spec.columns if c in df.columns]
    if spec.prefix:
        for c in df.columns:
            if str(c).startswith(spec.prefix) and c not in cols:
                cols.append(c)
    return cols


def main() -> None:
    ap = argparse.ArgumentParser(description="Validate a MolProp Toolkit results CSV against docs/schema.json")
    ap.add_argument("input", help="Results CSV")
    ap.add_argument(
        "--schema",
        default=None,
        help="Schema JSON path (defaults to bundled package schema; use docs/schema.json when running from source)",
    )
    ap.add_argument("--fail-unknown", action="store_true", help="Fail if unknown columns are found")
    ap.add_argument("--warn-unknown", action="store_true", help="Warn on unknown columns (default)")
    ap.add_argument(
        "--require-category",
        action="append",
        default=[],
        help="Require at least one column from a registry category to be present (repeatable)",
    )
    ap.add_argument(
        "--min-coercion",
        type=float,
        default=0.95,
        help="Minimum fraction of non-null values that must coerce to the schema type for typed columns",
    )
    args = ap.parse_args()

    schema = _load_schema(args.schema)
    patterns = _compile_patterns(schema)

    df = read_csv(args.input)

    # 1) Required explicit columns
    required_missing: List[str] = []
    for name, meta in (schema.get("columns", {}) or {}).items():
        if meta.get("required") and name not in df.columns:
            required_missing.append(name)

    # 2) Require categories
    category_missing: List[str] = []
    for cat in args.require_category or []:
        cols = _category_columns_present(df, cat)
        if not cols:
            category_missing.append(cat)

    # 3) Known vs unknown + type coercion
    unknown: List[str] = []
    type_issues: List[str] = []
    known_explicit = 0
    known_pattern = 0

    for c in df.columns:
        meta = _meta_for(str(c), schema, patterns)
        if meta is None:
            unknown.append(str(c))
            continue

        if meta.source == "explicit":
            known_explicit += 1
        else:
            known_pattern += 1

        # For pattern columns with type=auto, we do not check coercion.
        if meta.type in ("auto", "any"):
            continue

        ratio = _coercion_ok_ratio(meta.type, df[c])
        if ratio < float(args.min_coercion):
            type_issues.append(f"{c}: expected {meta.type}, coercion_ok={ratio:.2f}")

    # Reporting
    print("=" * 80)
    print("MolProp Toolkit CSV schema validation")
    print("=" * 80)
    print(f"csv: {args.input}")
    print(f"schema: {args.schema}")
    print(f"rows: {len(df)}")
    print(f"columns: {len(df.columns)}")
    print(f"known columns: {known_explicit + known_pattern} (explicit={known_explicit}, pattern={known_pattern})")
    print(f"unknown columns: {len(unknown)}")

    ok = True

    if required_missing:
        ok = False
        print("\nMissing required columns:")
        for c in required_missing:
            print(f"- {c}")

    if category_missing:
        ok = False
        print("\nMissing required categories (no matching columns present):")
        for cat in category_missing:
            print(f"- {cat}")

    if type_issues:
        ok = False
        print("\nType coercion issues:")
        for t in type_issues[:200]:
            print(f"- {t}")
        if len(type_issues) > 200:
            print(f"... +{len(type_issues) - 200} more")

    if unknown:
        if args.fail_unknown:
            ok = False
        if args.warn_unknown or (not args.fail_unknown):
            print("\nUnknown columns:")
            for c in unknown[:200]:
                print(f"- {c}")
            if len(unknown) > 200:
                print(f"... +{len(unknown) - 200} more")

    print("\nResult:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()

