# ruff: noqa: F401

"""Shared utilities for MolProp Toolkit.

This subpackage is intentionally lightweight and dependency-minimal; analyzers already depend on pandas/numpy.
"""

from __future__ import annotations

from .columns import detect_best_smiles_column, detect_id_column
from .io import detect_table_format, read_csv, read_table, write_csv, write_table
<<<<<<< HEAD
from .metadata import metadata_sidecar_path, sha256_file, write_run_metadata
=======
>>>>>>> 5d70153 (Core: prefer Calc_Canonical_SMILES + add Parquet table IO)
from .printing import print_banner, print_section
from .registry import CATEGORY_SPECS, CategorySpec

