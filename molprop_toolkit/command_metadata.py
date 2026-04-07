"""Shared command metadata for the public MolScope workflow.

This module is the lightweight source of truth for the command names we want to
surface publicly. The repository still contains several legacy entrypoints for
compatibility, but the site, README, and CLI help should all converge on the
same top-level workflow language.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class CommandSpec:
    command: str
    summary: str


WORKFLOW_COMMANDS: tuple[CommandSpec, ...] = (
    CommandSpec("molscope calc", "Create the recommended extended results table from SMILES."),
    CommandSpec("molscope analyze", "Inspect categories, summarize liabilities, and export focused slices."),
    CommandSpec("molscope report", "Build a shareable Markdown and HTML report bundle."),
    CommandSpec("molscope picklists", "Generate operational decision lists and optional HTML outputs."),
    CommandSpec("molscope portal", "Assemble report, compare, SAR, MMP, picklist, and retro bundles into one workspace."),
)


ADDITIONAL_COMMANDS: tuple[CommandSpec, ...] = (
    CommandSpec("molscope prep", "Apply named ligand-preparation presets and write a prep audit table."),
    CommandSpec("molscope sketch", "Render 2D structures and browseable sketch sheets."),
    CommandSpec("molscope integrate", "Merge measured assay endpoints into the canonical MolScope results table."),
    CommandSpec("molscope compare", "Compare two campaign tables and summarize dataset, endpoint, and scaffold shifts."),
    CommandSpec("molscope learnings", "Synthesize compare, SAR, MMP, and decision outputs into reusable campaign lessons."),
    CommandSpec("molscope dashboard", "Accumulate learnings bundles into a persistent campaign dashboard and knowledge registry."),
    CommandSpec("molscope sar", "Build a scaffold-centric SAR bundle with series summaries and top compounds."),
    CommandSpec("molscope mmp", "Build matched-pair tables and aggregate transform summaries from a campaign table."),
    CommandSpec("molscope search", "Search campaign tables by exact SMILES, SMARTS substructure, or motif panel."),
    CommandSpec("molscope series", "Annotate scaffolds and similarity-based series clusters."),
    CommandSpec("molscope similarity", "Run similarity search, matrices, diversity picking, and clustering."),
    CommandSpec("molscope featurize", "Export ML-ready feature packs and optional deep-learning assets."),
    CommandSpec("molscope retro", "Run offline retrosynthesis via AiZynthFinder and rerank route summaries."),
    CommandSpec("molscope schema", "Validate results tables against the published column schema."),
)


CALCULATOR_PROFILES: dict[str, str] = {
    "extended": "Recommended default profile for full medchem and developability triage.",
    "baseline": "Smaller medchem-focused profile for quicker, lighter-weight tables.",
}


LEGACY_ALIASES: dict[str, str] = {
    "molscope calc": "molprop-calc-v5 (default) or molprop-calc-v4 (--profile baseline)",
    "molscope prep": "molprop-prep",
    "molscope analyze": "molprop-analyze",
    "molscope report": "molprop-report",
    "molscope picklists": "molprop-picklists",
    "molscope portal": "molprop-portal",
    "molscope sketch": "molprop-sketch",
    "molscope integrate": "molprop-integrate",
    "molscope compare": "molprop-compare",
    "molscope learnings": "molprop-learnings",
    "molscope dashboard": "molprop-dashboard",
    "molscope sar": "molprop-sar",
    "molscope mmp": "molprop-mmp",
    "molscope search": "molprop-search",
    "molscope series": "molprop-series",
    "molscope similarity": "molprop-similarity",
    "molscope featurize": "molprop-featurize",
    "molscope retro": "molprop-retro",
    "molscope schema": "python tools/validate_csv_schema.py",
}


def workflow_help_lines() -> list[str]:
    lines = ["Core workflow:"]
    for spec in WORKFLOW_COMMANDS:
        if spec.command == "molscope calc":
            lines.append("  molscope calc [--profile baseline|extended] <input> -o <output>")
        elif spec.command == "molscope picklists":
            lines.append("  molscope picklists <results.csv> --html")
        elif spec.command == "molscope portal":
            lines.append("  molscope portal --report-dir <dir> [--picklists-dir <dir> ...]")
        else:
            lines.append(f"  {spec.command} <results.csv>")
    return lines
