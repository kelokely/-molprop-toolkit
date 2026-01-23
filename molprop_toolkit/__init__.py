"""MolProp Toolkit (importable package).

The repository originally started as a collection of CLI scripts.

Task 10 introduces a small internal package namespace so shared logic (column
name detection, CSV IO, printing, and registries) can live in one place and be
used consistently by analyzers/tools.

The CLI scripts remain runnable directly (e.g. `python analyzers/analyzer_cns.py`).
"""

from __future__ import annotations

