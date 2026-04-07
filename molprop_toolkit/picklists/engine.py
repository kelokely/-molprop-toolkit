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
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, Sequence, Tuple

import pandas as pd

from molprop_toolkit.core import detect_id_column, read_table
from molprop_toolkit.core.registry import CATEGORY_SPECS
from molprop_toolkit.experiment_design import suggest_experiment_plan

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
    queue: str = ""
    next_action: str = ""
    rationale: str = ""
    evidence_columns: Sequence[str] = ()
    diversity_columns: Sequence[str] = ()
    max_per_diversity_group: int = 0


@dataclass(frozen=True)
class PicklistResult:
    definition: PicklistDefinition
    df: pd.DataFrame
    warnings: Sequence[str]


DECISION_OUTPUT_COLUMNS: tuple[str, ...] = (
    "Decision_Queue",
    "Decision_Rationale",
    "Decision_Next_Action",
    "Decision_Confidence",
    "Decision_Evidence_Summary",
    "Decision_Hypothesis",
    "Decision_Minimal_Test_Set",
    "Decision_Minimal_Analog_Set",
)


def _coerce_numeric(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce")


def _value_text(value: Any) -> str:
    if value is None or pd.isna(value):
        return "(missing)"
    if isinstance(value, float):
        return f"{value:.3g}"
    return str(value)


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
        cols.extend(c for c in DECISION_OUTPUT_COLUMNS if c not in cols)
        return [c for c in cols if c in df.columns]

    # Minimal default: id + columns referenced by filters/sort + key columns for included categories.
    cols: List[str] = [id_col]

    for f in definition.filters:
        if f.column in df.columns and f.column not in cols:
            cols.append(f.column)

    for s in definition.sort:
        if s.column in df.columns and s.column not in cols:
            cols.append(s.column)

    for c in definition.diversity_columns:
        if c in df.columns and c not in cols:
            cols.append(c)
            break

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

    for c in DECISION_OUTPUT_COLUMNS:
        if c in df.columns and c not in cols:
            cols.append(c)

    return cols


def _resolve_diversity_column(
    df: pd.DataFrame,
    definition: PicklistDefinition,
) -> str | None:
    for column in definition.diversity_columns:
        if column in df.columns:
            return column
    return None


def _apply_diversity_limit(
    df: pd.DataFrame,
    *,
    group_column: str,
    max_per_group: int,
) -> pd.DataFrame:
    if max_per_group <= 0 or group_column not in df.columns or df.empty:
        return df

    counts: dict[str, int] = defaultdict(int)
    keep_index: list[object] = []

    for index, value in df[group_column].items():
        if pd.isna(value) or str(value).strip() == "":
            keep_index.append(index)
            continue

        key = str(value).strip()
        if counts[key] >= max_per_group:
            continue
        counts[key] += 1
        keep_index.append(index)

    return df.loc[keep_index].copy()


def _row_matches_rule(row: pd.Series, rule: FilterRule) -> bool:
    if rule.column not in row.index:
        return False

    value = row.get(rule.column)
    op = rule.op

    if op == "==":
        return bool(value == rule.value)
    if op == "!=":
        return bool(value != rule.value)
    if op in (">", ">=", "<", "<="):
        numeric = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
        if pd.isna(numeric):
            return False
        threshold = float(rule.value)
        if op == ">":
            return bool(numeric > threshold)
        if op == ">=":
            return bool(numeric >= threshold)
        if op == "<":
            return bool(numeric < threshold)
        return bool(numeric <= threshold)
    if op in ("in", "notin"):
        values = rule.value if isinstance(rule.value, (list, tuple, set)) else [rule.value]
        matched = value in values
        return matched if op == "in" else not matched
    if op in ("contains", "notcontains"):
        text = "" if value is None or pd.isna(value) else str(value)
        matched = str(rule.value or "").lower() in text.lower()
        return matched if op == "contains" else not matched
    if op == "isnull":
        return bool(pd.isna(value))
    if op == "notnull":
        return bool(not pd.isna(value))
    if op == "between":
        numeric = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
        if pd.isna(numeric):
            return False
        return bool(float(rule.value) <= numeric <= float(rule.value2))
    raise ValueError(f"Unknown filter op: {op}")


def _describe_rule_match(row: pd.Series, rule: FilterRule) -> str:
    value = row.get(rule.column)
    rendered_value = _value_text(value)
    if rule.op == ">=":
        return f"{rule.column}={rendered_value} meets >= {rule.value}"
    if rule.op == ">":
        return f"{rule.column}={rendered_value} exceeds {rule.value}"
    if rule.op == "<=":
        return f"{rule.column}={rendered_value} stays <= {rule.value}"
    if rule.op == "<":
        return f"{rule.column}={rendered_value} stays < {rule.value}"
    if rule.op == "==":
        return f"{rule.column}={rendered_value} matches {rule.value}"
    if rule.op == "!=":
        return f"{rule.column}={rendered_value} differs from {rule.value}"
    if rule.op == "in":
        return f"{rule.column}={rendered_value} matches allowed set"
    if rule.op == "notin":
        return f"{rule.column}={rendered_value} avoids {', '.join(str(v) for v in (rule.value if isinstance(rule.value, (list, tuple, set)) else [rule.value]))}"
    if rule.op == "contains":
        return f"{rule.column} contains '{rule.value}'"
    if rule.op == "notcontains":
        return f"{rule.column} avoids '{rule.value}'"
    if rule.op == "isnull":
        return f"{rule.column} is still missing"
    if rule.op == "notnull":
        return f"{rule.column} is available"
    if rule.op == "between":
        return f"{rule.column}={rendered_value} stays between {rule.value} and {rule.value2}"
    return f"{rule.column} matches {rule.op}"


def _column_evidence_state(row: pd.Series, column: str) -> str:
    value = row.get(column)
    if pd.isna(value):
        return "missing"

    provenance_column = f"{column}_Provenance"
    provenance = row.get(provenance_column)
    if provenance is not None and not pd.isna(provenance):
        provenance_text = str(provenance).strip().lower()
        if provenance_text in {"measured", "derived_from_measured"}:
            return "measured"
        if provenance_text in {"heuristic"}:
            return "heuristic"
        return "calculated"

    if str(column).startswith("Measured_"):
        return "measured"
    return "calculated"


def decision_evidence_summary(
    row: pd.Series,
    evidence_columns: Sequence[str],
) -> str:
    buckets = decision_evidence_buckets(row, evidence_columns)
    parts: list[str] = []
    for key in ("measured", "calculated", "heuristic", "missing"):
        cols = buckets.get(key, [])
        if cols:
            parts.append(f"{key}: {', '.join(cols)}")
    return "; ".join(parts) if parts else "No evidence signals detected."


def decision_evidence_buckets(
    row: pd.Series,
    evidence_columns: Sequence[str],
) -> dict[str, list[str]]:
    if not evidence_columns:
        return {"measured": [], "calculated": [], "heuristic": [], "missing": []}

    buckets: dict[str, list[str]] = {
        "measured": [],
        "calculated": [],
        "heuristic": [],
        "missing": [],
    }

    for column in evidence_columns:
        if column not in row.index and f"{column}_Provenance" not in row.index:
            buckets["missing"].append(column)
            continue
        state = _column_evidence_state(row, column)
        buckets.setdefault(state, []).append(column)

    return buckets


def decision_confidence(
    row: pd.Series,
    evidence_columns: Sequence[str],
) -> str:
    if not evidence_columns:
        return "medium"

    buckets = decision_evidence_buckets(row, evidence_columns)
    measured = len(buckets["measured"])
    missing = len(buckets["missing"])
    available = measured + len(buckets["calculated"]) + len(buckets["heuristic"])

    if available == 0:
        return "low"
    if measured >= 2 or (measured >= 1 and missing == 0):
        return "high"
    if measured >= 1 or available >= 2:
        return "medium"
    return "low"


def _decision_rationale(row: pd.Series, definition: PicklistDefinition) -> str:
    matched_reasons = [
        _describe_rule_match(row, rule)
        for rule in definition.filters
        if _row_matches_rule(row, rule)
    ]
    matched_reasons = matched_reasons[:4]
    base = definition.rationale.strip() or f"Selected for {definition.name.lower()}."
    if not matched_reasons:
        return base
    return base + " Why: " + "; ".join(matched_reasons) + "."


def annotate_picklist_decisions(
    df: pd.DataFrame,
    definition: PicklistDefinition,
    *,
    stage: str | None = None,
) -> pd.DataFrame:
    out = df.copy()
    queue = definition.queue or definition.key
    out["Decision_Queue"] = queue
    out["Decision_Next_Action"] = (
        definition.next_action.strip()
        or f"Review {definition.name.lower()} as the next operational queue."
    )
    evidence_columns = tuple(definition.evidence_columns or ())
    out["Decision_Rationale"] = out.apply(
        lambda row: _decision_rationale(row, definition),
        axis=1,
    )
    out["Decision_Confidence"] = out.apply(
        lambda row: decision_confidence(row, evidence_columns),
        axis=1,
    )
    out["Decision_Evidence_Summary"] = out.apply(
        lambda row: decision_evidence_summary(row, evidence_columns),
        axis=1,
    )
    plans = out.apply(
        lambda row: suggest_experiment_plan(row, queue=queue, stage=stage),
        axis=1,
    )
    out["Decision_Hypothesis"] = plans.apply(lambda plan: plan.hypothesis)
    out["Decision_Minimal_Test_Set"] = plans.apply(
        lambda plan: plan.minimal_test_set
    )
    out["Decision_Minimal_Analog_Set"] = plans.apply(
        lambda plan: plan.minimal_analog_set
    )
    return out


def run_picklist(
    df: pd.DataFrame,
    definition: PicklistDefinition,
    *,
    stage: str | None = None,
) -> PicklistResult:
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

    diversity_column = _resolve_diversity_column(out, definition)
    if definition.max_per_diversity_group > 0:
        if diversity_column is None and definition.diversity_columns:
            warnings.append(
                "No diversity column found among "
                + ", ".join(definition.diversity_columns)
                + "; queue remains score-sorted only."
            )
        elif diversity_column is not None:
            out = _apply_diversity_limit(
                out,
                group_column=diversity_column,
                max_per_group=int(definition.max_per_diversity_group),
            )

    # Limit
    if definition.limit and definition.limit > 0:
        out = out.head(int(definition.limit))

    out = annotate_picklist_decisions(out, definition, stage=stage)

    return PicklistResult(definition=definition, df=out, warnings=warnings)


def build_picklists_from_df(
    df: pd.DataFrame,
    definitions: Sequence[PicklistDefinition],
    *,
    id_col: Optional[str] = None,
    stage: str | None = None,
) -> Tuple[pd.DataFrame, str, List[PicklistResult]]:
    """Run picklists using an in-memory dataframe.

    This is useful when the caller wants to add computed columns (for example, similarity-to-reference columns) before
    applying picklist rules.
    """

    if id_col is None:
        id_col = detect_id_column(df)

    results: List[PicklistResult] = []
    for d in definitions:
        results.append(run_picklist(df, d, stage=stage))
    return df, id_col, results


def build_picklists(
    csv_path: str,
    definitions: Sequence[PicklistDefinition],
    *,
    stage: str | None = None,
) -> Tuple[pd.DataFrame, str, List[PicklistResult]]:
    """Backwards-compatible helper that reads a CSV and runs picklists."""

    df = read_table(csv_path)
    return build_picklists_from_df(df, definitions, stage=stage)


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
                queue=str(it.get("queue") or ""),
                next_action=str(it.get("next_action") or ""),
                rationale=str(it.get("rationale") or ""),
                evidence_columns=list(it.get("evidence_columns") or []),
            )
        )

    return defs
