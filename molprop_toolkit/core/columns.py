"""Column detection helpers.

The toolkit accepts many possible upstream CSV formats, so we try multiple common identifier and SMILES column names.
"""

from __future__ import annotations

from typing import Optional, Sequence, Union


ID_CANDIDATES: Sequence[str] = (
    "Compound_ID",
    "compound_id",
    "Compound",
    "compound",
    "Name",
    "name",
    "ID",
    "id",
)

SMILES_PRIORITY: Sequence[str] = (
    # Most traceable (calculation structure)
    "Calc_Canonical_SMILES",
    # Pre-protomer calc base (tautomer+stereo)
    "Calc_Base_SMILES",
    # Prepared (standardized)
    "Canonical_SMILES",
    # Input canonicalized
    "Input_Canonical_SMILES",
    # Raw
    "SMILES",
)


def _as_columns(df_or_columns: Union[Sequence[str], object]) -> list[str]:
    # Avoid importing pandas here; accept anything with a .columns attribute.
    cols = None
    if hasattr(df_or_columns, "columns"):
        cols = getattr(df_or_columns, "columns")
    else:
        cols = df_or_columns

    try:
        return list(cols)  # type: ignore[arg-type]
    except Exception:
        return []


def detect_id_column(df_or_columns: Union[Sequence[str], object]) -> str:
    """Detect a compound identifier column.

    Accepts a pandas DataFrame or a list of column names.
    """

    cols = _as_columns(df_or_columns)
    for c in ID_CANDIDATES:
        if c in cols:
            return c
    return cols[0] if cols else "Compound_ID"


def detect_best_smiles_column(
    columns: Sequence[str],
    priority: Sequence[str] = SMILES_PRIORITY,
) -> Optional[str]:
    """Pick the best SMILES column for rendering or downstream operations."""

    cols = list(columns)
    for c in priority:
        if c in cols:
            return c
    for c in cols:
        if "smiles" in c.lower():
            return c
    return None

