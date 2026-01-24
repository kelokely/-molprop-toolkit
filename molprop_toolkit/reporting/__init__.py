# ruff: noqa: F401

"""Report generation utilities.

This subpackage produces shareable report artifacts (Markdown + HTML + plots) from a results CSV.

Design goal: remain registry-driven by consuming `molprop_toolkit.core.registry.CATEGORY_SPECS` rather than hard-coding
category knowledge into the report builder.
"""

from __future__ import annotations

from .report_builder import ReportArtifacts, ReportBuilder

