# ruff: noqa: F401

"""Decision filters and picklist generation.

This subpackage turns a wide results CSV into operational outputs: one or more filtered picklist CSVs plus optional
lightweight HTML views. If a structure-sketch directory is provided, required structure images can be copied into the picklist
output folder so the HTML is portable.

The picklist engine is designed to work without RDKit; it only depends on pandas and the results table.
"""

from __future__ import annotations

from .builtins import BUILTIN_PICKLISTS
from .engine import (
    PicklistDefinition,
    PicklistResult,
    build_picklists,
    load_picklist_config,
)
