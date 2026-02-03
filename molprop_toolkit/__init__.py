"""MolProp Toolkit (importable package).

The repository started as a collection of CLI scripts. The `molprop_toolkit` package provides a small internal namespace so
shared logic (column detection, CSV IO, printing, and registries) lives in one place and is used consistently by analyzers and
tools.

The original scripts remain runnable directly (for example, `python analyzers/analyzer_cns.py`).
"""

from __future__ import annotations
