#!/usr/bin/env python3
"""Generate a small command summary page from the shared command metadata."""

from __future__ import annotations

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from molprop_toolkit.command_metadata import (
    ADDITIONAL_COMMANDS,
    CALCULATOR_PROFILES,
    LEGACY_ALIASES,
    WORKFLOW_COMMANDS,
)


def render_markdown() -> str:
    lines: list[str] = []
    lines.append("# MolScope command summary")
    lines.append("")
    lines.append("This page is generated from `molprop_toolkit.command_metadata` and summarizes the public command surface.")
    lines.append("")
    lines.append("## Core workflow")
    lines.append("")
    for spec in WORKFLOW_COMMANDS:
        lines.append(f"- `{spec.command}`: {spec.summary}")
    lines.append("")
    lines.append("## Additional commands")
    lines.append("")
    for spec in ADDITIONAL_COMMANDS:
        lines.append(f"- `{spec.command}`: {spec.summary}")
    lines.append("")
    lines.append("## Calculator profiles")
    lines.append("")
    for key, description in CALCULATOR_PROFILES.items():
        lines.append(f"- `{key}`: {description}")
    lines.append("")
    lines.append("## Legacy aliases")
    lines.append("")
    for canonical, legacy in LEGACY_ALIASES.items():
        lines.append(f"- `{canonical}` → `{legacy}`")
    lines.append("")
    return "\n".join(lines)


def main() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    out_path = repo_root / "docs" / "commands.md"
    out_path.write_text(render_markdown(), encoding="utf-8")
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
