# ruff: noqa: F401

"""Shared utilities for MolProp Toolkit (Task 10).

This subpackage is intentionally lightweight and dependency-minimal; analyzers
already depend on pandas/numpy.
"""

from __future__ import annotations

from .columns import detect_best_smiles_column, detect_id_column
from .io import read_csv, write_csv
from .printing import print_banner, print_section
from .registry import CATEGORY_SPECS, CategorySpec

