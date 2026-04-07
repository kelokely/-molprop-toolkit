"""Table IO helpers (CSV/TSV/Parquet).

MolProp Toolkit is "table-first": calculators and downstream tools exchange results through a tabular artifact.

These helpers provide:
- consistent defaults for large descriptor tables,
- transparent support for CSV/TSV and Parquet,
- centralized format detection so all CLIs behave the same.

Parquet support requires `pyarrow` (recommended) or another pandas parquet engine.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Literal, Optional

import pandas as pd

TableFormat = Literal["csv", "tsv", "parquet"]


def detect_table_format(path: str, fmt: Optional[str] = None) -> TableFormat:
    """Detect table format.

    If fmt is provided and not 'auto', it takes precedence.
    Otherwise, detect from file extension.
    """

    if fmt and fmt.lower() != "auto":
        f = fmt.lower()
        if f in ("csv", "tsv", "parquet"):
            return f  # type: ignore[return-value]
        raise ValueError(f"Unknown table format: {fmt}")

    ext = Path(path).suffix.lower()
    if ext in (".parquet", ".pq"):
        return "parquet"
    if ext in (".tsv", ".tab"):
        return "tsv"
    return "csv"


def read_table(path: str, *, fmt: Optional[str] = None, **kwargs: Any) -> pd.DataFrame:
    """Read a table (CSV/TSV/Parquet) with consistent defaults."""

    f = detect_table_format(path, fmt)

    if f == "parquet":
        try:
            return pd.read_parquet(path, **kwargs)
        except ImportError as e:  # pragma: no cover
            raise RuntimeError(
                "Reading Parquet requires pyarrow (recommended). Install with: pip install pyarrow\n"
                "or install molprop-toolkit with the parquet extra: pip install 'molprop-toolkit[parquet]'"
            ) from e

    # CSV/TSV defaults
    if "low_memory" not in kwargs:
        # Avoid mixed-type inference surprises on large descriptor tables.
        kwargs["low_memory"] = False

    if f == "tsv":
        kwargs.setdefault("sep", "\t")

    return pd.read_csv(path, **kwargs)


def write_table(
    df: pd.DataFrame, path: str, *, fmt: Optional[str] = None, **kwargs: Any
) -> None:
    """Write a table (CSV/TSV/Parquet), creating parent directories."""

    Path(path).parent.mkdir(parents=True, exist_ok=True)

    f = detect_table_format(path, fmt)

    if f == "parquet":
        try:
            # Parquet should preserve dtypes; never store the index.
            return df.to_parquet(path, index=False, **kwargs)
        except ImportError as e:  # pragma: no cover
            raise RuntimeError(
                "Writing Parquet requires pyarrow (recommended). Install with: pip install pyarrow\n"
                "or install molprop-toolkit with the parquet extra: pip install 'molprop-toolkit[parquet]'"
            ) from e

    # CSV/TSV
    if f == "tsv":
        kwargs.setdefault("sep", "\t")

    return df.to_csv(path, index=False, **kwargs)


# Backwards-compatible wrappers (older tools import these names)


def read_csv(path: str, **kwargs: Any) -> pd.DataFrame:
    return read_table(path, fmt="csv", **kwargs)


def write_csv(df: pd.DataFrame, path: str) -> None:
    return write_table(df, path, fmt="csv")
