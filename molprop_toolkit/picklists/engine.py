"""Picklist engine.

A picklist is defined as:
- a list of filters (column/operator/value)
- optional sorting
- optional row limit

The engine is intentionally conservative:
- Filters are structured objects (not arbitrary pandas query strings)
- Missing columns can be treated as warnings (default) or failures
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, Sequence, Tuple

import pandas as pd

from molprop_toolkit.core import detect_id_column, read_table
from molprop_toolkit.core.registry import CATEGORY_SPECS

FilterOp = Literal[
    "==",
    "!=",
    ">",
    ">=",
    "<",
    "<=",
    "in",
    "notin",
    "contains",
    "notcontains",
    "isnull",
    "notnull",
    "between",
]


MissingPolicy = Literal["warn_skip", "fail", "false"]
LogicMode = Literal["all", "any"]


@dataclass(frozen=True)
class FilterRule:
    column: str
    op: FilterOp
    value: Any = None
    value2: Any = None  # for between
    missing: MissingPolicy = "warn_skip"


@dataclass(frozen=True)
class SortRule:
    column: str
    ascending: bool = False


@dataclass(frozen=True)
class PicklistDefinition:
    key: str
    name: str
    description: str = ""
    filters: Sequence[FilterRule] = ()
    logic: LogicMode = "all"  # "all" = AND; "any" = OR
    sort: Sequence[SortRule] = ()
    limit: int = 0
    include_categories: Sequence[str] = ()
    output_columns: Sequence[str] = ()  # explicit override


@dataclass(frozen=True)
class PicklistResult:
    definition: PicklistDefinition
    df: pd.DataFrame
    warnings: Sequence[str]


def _coerce_numeric(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce")


def _apply_filter(
    df: pd.DataFrame,
    rule: FilterRule,
    neutral_on_skip: bool,
) -> Tuple[pd.Series, Optional[str]]:
    """Apply a single filter rule.

    `neutral_on_skip` controls the behavior when a column is missing and the
    rule's missing policy is warn_skip.

    - For AND logic, neutral should be True (skipping a condition keeps rows).
    - For OR logic, neutral should be False (skipping a condition should not add rows).
    """

    if rule.column not in df.columns:
        if rule.missing == "fail":
            raise ValueError(f"Missing required column for filter: {rule.column}")
        if rule.missing == "false":
            return pd.Series([False] * len(df), index=df.index), None
        # warn_skip
        return (
            pd.Series([bool(neutral_on_skip)] * len(df), index=df.index),
            f"Missing column '{rule.column}' → skipped filter ({rule.op})",
        )

    s = df[rule.column]
    op = rule.op

    if op in ("==", "!="):
        mask = s == rule.value if op == "==" else s != rule.value
        return mask.fillna(False), None

    if op in (">", ">=", "<", "<="):
        sn = _coerce_numeric(s)
        vn = float(rule.value)
        if op == ">":
            mask = sn > vn
        elif op == ">=":
            mask = sn >= vn
        elif op == "<":
            mask = sn < vn
        else:
            mask = sn <= vn
        return mask.fillna(False), None

    if op in ("in", "notin"):
        values = rule.value
        if not isinstance(values, (list, tuple, set)):
            values = [values]
        mask = s.isin(list(values))
        if op == "notin":
            mask = ~mask
        return mask.fillna(False), None

    if op in ("contains", "notcontains"):
        patt = str(rule.value or "")
        ss = s.fillna("").astype(str)
        mask = ss.str.contains(patt, case=False, regex=False)
        if op == "notcontains":
            mask = ~mask
        return mask.fillna(False), None

    if op == "isnull":
        return s.isna(), None

    if op == "notnull":
        return ~s.isna(), None

    if op == "between":
        sn = _coerce_numeric(s)
        lo = float(rule.value)
        hi = float(rule.value2)
        mask = (sn >= lo) & (sn <= hi)
        return mask.fillna(False), None

    raise ValueError(f"Unknown filter op: {op}")


def registry_category_columns(df: pd.DataFrame, category_key: str) -> List[str]:
    """Return dataset columns that belong to a registry category.

    Uses CategorySpec.columns plus optional prefix matching.
    """

    if category_key not in CATEGORY_SPECS:
        return []
    spec = CATEGORY_SPECS[category_key]

    cols: List[str] = [c for c in spec.columns if c in df.columns]
    if spec.prefix:
        for c in df.columns:
            if str(c).startswith(spec.prefix) and c not in cols:
                cols.append(c)
    return cols


def select_output_columns(
    df: pd.DataFrame,
    id_col: str,
    definition: PicklistDefinition,
) -> List[str]:
    if definition.output_columns:
        cols = [id_col] + [c for c in definition.output_columns if c != id_col]
        return [c for c in cols if c in df.columns]

    # Minimal default: id + columns referenced by filters/sort + key columns for included categories.
    cols: List[str] = [id_col]

    for f in definition.filters:
        if f.column in df.columns and f.column not in cols:
            cols.append(f.column)

    for s in definition.sort:
        if s.column in df.columns and s.column not in cols:
            cols.append(s.column)

    for cat in definition.include_categories:
        if cat not in CATEGORY_SPECS:
            continue
        spec = CATEGORY_SPECS[cat]
        if (
            spec.key_column
            and spec.key_column in df.columns
            and spec.key_column not in cols
        ):
            cols.append(spec.key_column)

        for c in registry_category_columns(df, cat):
            if c not in cols:
                cols.append(c)
            if len(cols) >= 24:
                break

    return cols


def run_picklist(df: pd.DataFrame, definition: PicklistDefinition) -> PicklistResult:
    warnings: List[str] = []

    logic = definition.logic or "all"
    if logic not in ("all", "any"):
        raise ValueError(f"Unknown picklist logic: {logic}")

    if logic == "all":
        mask = pd.Series([True] * len(df), index=df.index)
    else:
        mask = pd.Series([False] * len(df), index=df.index)

    for rule in definition.filters:
        m, w = _apply_filter(df, rule, neutral_on_skip=(logic == "all"))
        if w:
            warnings.append(w)
        mask = (mask & m) if logic == "all" else (mask | m)

    out = df[mask].copy()

    # Sorting
    if definition.sort:
        sort_cols = [s.column for s in definition.sort if s.column in out.columns]
        ascending = [s.ascending for s in definition.sort if s.column in out.columns]
        if sort_cols:
            # Try numeric sort first for each column.
            for c in sort_cols:
                if c in out.columns and pd.api.types.is_object_dtype(out[c]):
                    # If object but numeric-looking, coerce.
                    out[c] = pd.to_numeric(out[c], errors="ignore")
            out = out.sort_values(by=sort_cols, ascending=ascending, kind="mergesort")

    # Limit
    if definition.limit and definition.limit > 0:
        out = out.head(int(definition.limit))

    return PicklistResult(definition=definition, df=out, warnings=warnings)


def build_picklists_from_df(
    df: pd.DataFrame,
    definitions: Sequence[PicklistDefinition],
    *,
    id_col: Optional[str] = None,
) -> Tuple[pd.DataFrame, str, List[PicklistResult]]:
    """Run picklists using an in-memory dataframe.

    This is useful when the caller wants to add computed columns (for example, similarity-to-reference columns) before
    applying picklist rules.
    """

    if id_col is None:
        id_col = detect_id_column(df)

    results: List[PicklistResult] = []
    for d in definitions:
        results.append(run_picklist(df, d))
    return df, id_col, results


def build_picklists(
    csv_path: str,
    definitions: Sequence[PicklistDefinition],
) -> Tuple[pd.DataFrame, str, List[PicklistResult]]:
    """Backwards-compatible helper that reads a CSV and runs picklists."""

    df = read_table(csv_path)
    return build_picklists_from_df(df, definitions)


def _load_json_or_yaml(path: str | Path) -> Dict[str, Any]:
    p = Path(path)
    txt = p.read_text(encoding="utf-8")

    if p.suffix.lower() in (".yaml", ".yml"):
        try:
            import yaml  # type: ignore
        except Exception as e:  # pragma: no cover
            raise RuntimeError(
                "YAML config requested but PyYAML is not installed. Install with: pip install pyyaml"
            ) from e
        data = yaml.safe_load(txt)
        return data or {}

    return json.loads(txt)


def load_picklist_config(path: str | Path) -> List[PicklistDefinition]:
    """Load picklists from a JSON/YAML config.

    Expected shape:
    {
      "picklists": [
        {
          "key": "my_list",
          "name": "My List",
          "description": "...",
          "filters": [{"column":"CNS_MPO","op":">=","value":4.5}],
          "sort": [{"column":"CNS_MPO","ascending": false}],
          "limit": 50,
          "include_categories": ["cns_mpo","toxicity"]
        }
      ]
    }

    or with a dict mapping keys → picklist objects:
    {"picklists": {"my_list": {...}}}
    """

    raw = _load_json_or_yaml(path)
    pls = raw.get("picklists", [])

    items: List[Dict[str, Any]]
    if isinstance(pls, dict):
        items = []
        for k, v in pls.items():
            vv = dict(v)
            vv.setdefault("key", k)
            items.append(vv)
    elif isinstance(pls, list):
        items = pls
    else:
        raise ValueError("Config 'picklists' must be a list or dict")

    defs: List[PicklistDefinition] = []
    for it in items:
        key = str(it.get("key") or "").strip()
        name = str(it.get("name") or key).strip()
        if not key:
            raise ValueError("Picklist is missing 'key'")

        filters_in = it.get("filters", []) or []
        filters: List[FilterRule] = []
        for f in filters_in:
            filters.append(
                FilterRule(
                    column=str(f.get("column")),
                    op=str(f.get("op")),  # type: ignore[arg-type]
                    value=f.get("value"),
                    value2=f.get("value2"),
                    missing=str(f.get("missing") or "warn_skip"),  # type: ignore[arg-type]
                )
            )

        sort_in = it.get("sort", []) or []
        sort: List[SortRule] = []
        for s in sort_in:
            sort.append(
                SortRule(
                    column=str(s.get("column")),
                    ascending=bool(s.get("ascending", False)),
                )
            )

        defs.append(
            PicklistDefinition(
                key=key,
                name=name,
                description=str(it.get("description") or ""),
                filters=filters,
                logic=str(it.get("logic") or "all"),  # type: ignore[arg-type]
                sort=sort,
                limit=int(it.get("limit") or 0),
                include_categories=list(it.get("include_categories") or []),
                output_columns=list(it.get("output_columns") or []),
            )
        )

    return defs
