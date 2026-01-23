# ruff: noqa: F401

"""Decision filters + pick list generation (Task 14).

This subpackage turns a wide results CSV into operational outputs:
- one or more filtered picklist CSVs
- optional lightweight HTML views (with optional structure images copied from
  Task 5 sketch assets)

The picklist engine is designed to work without RDKit; it only depends on
pandas and the results table.
"""

from __future__ import annotations

from .builtins import BUILTIN_PICKLISTS
from .engine import PicklistDefinition, PicklistResult, build_picklists, load_picklist_config

