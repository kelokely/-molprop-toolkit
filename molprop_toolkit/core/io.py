"""CSV IO helpers.

The repo uses pandas throughout analyzers; these helpers standardize read/write
settings and ensure output directories exist.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd


def read_csv(path: str, **kwargs: Any) -> pd.DataFrame:
    """Read a CSV with consistent defaults."""

    # low_memory=False avoids mixed-type inference surprises on large descriptor tables.
    if "low_memory" not in kwargs:
        kwargs["low_memory"] = False
    return pd.read_csv(path, **kwargs)


def write_csv(df: pd.DataFrame, path: str) -> None:
    """Write a CSV creating the parent directory if needed."""

    Path(path).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)

