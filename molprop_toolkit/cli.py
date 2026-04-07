"""Console-script entrypoints.

These wrappers keep existing workflows intact.

The original scripts remain runnable as `python calculators/...` and `python tools/...` while the package also exposes
installable console scripts that call the same underlying `main()` functions.

Design note: these entrypoints intentionally avoid re-implementing the full
argument parsing of each tool. The new ``molscope`` wrapper only normalizes the
top-level workflow and then delegates to the original CLIs.
"""

from __future__ import annotations

import sys
from typing import Callable, Sequence

from molprop_toolkit.command_metadata import ADDITIONAL_COMMANDS, workflow_help_lines


def _run_with_argv(func: Callable[[], None], argv0: str, argv: Sequence[str]) -> None:
    original_argv = sys.argv[:]
    try:
        sys.argv = [argv0, *argv]
        func()
    finally:
        sys.argv = original_argv


def _parse_calc_profile(argv: Sequence[str]) -> tuple[str, list[str]]:
    profile = "extended"
    forwarded: list[str] = []

    i = 0
    while i < len(argv):
        arg = argv[i]
        if arg == "--profile":
            if i + 1 >= len(argv):
                raise SystemExit("molscope calc: expected a value after --profile")
            profile = argv[i + 1].strip().lower()
            i += 2
            continue
        if arg.startswith("--profile="):
            profile = arg.split("=", 1)[1].strip().lower()
            i += 1
            continue
        forwarded.append(arg)
        i += 1

    if profile not in {"baseline", "extended"}:
        raise SystemExit(
            "molscope calc: --profile must be one of: baseline, extended"
        )

    return profile, forwarded


def _print_main_help() -> None:
    lines = [
        "MolScope Toolkit CLI",
        "",
        "Usage:",
        "  molscope <command> [options]",
        "",
        *workflow_help_lines(),
        "",
        "Additional commands:",
    ]
    for spec in ADDITIONAL_COMMANDS:
        lines.append(f"  {spec.command}")
    lines.extend(
        [
            "",
            "Notes:",
            "  - The default calculator profile is 'extended'.",
            "  - Legacy molprop-* console scripts remain supported for compatibility.",
        ]
    )
    print("\n".join(lines))


def calc_v4() -> None:
    from calculators.mpo_v4 import main

    main()


def calc_v5() -> None:
    from calculators.mpo_v5 import main

    main()


def analyze() -> None:
    # Delegates to the registry-driven master analyzer.
    from analyzers.analyzer_master import main

    main()


def report() -> None:
    from tools.build_report import main

    main()


def portal() -> None:
    from molprop_toolkit.tools.portal_cli import main

    main()


def learnings() -> None:
    from molprop_toolkit.tools.learnings_cli import main

    main()


def dashboard() -> None:
    from molprop_toolkit.tools.dashboard_cli import main

    main()


def prep() -> None:
    from molprop_toolkit.tools.prep_cli import main

    main()


def sketch() -> None:
    from tools.sketch_structures import main

    main()


def picklists() -> None:
    from tools.picklists import main

    main()


def integrate() -> None:
    from molprop_toolkit.tools.integrate_cli import main

    main()


def compare() -> None:
    from molprop_toolkit.tools.compare_cli import main

    main()


def sar() -> None:
    from molprop_toolkit.tools.sar_cli import main

    main()


def mmp() -> None:
    from molprop_toolkit.tools.mmp_cli import main

    main()


def search() -> None:
    from molprop_toolkit.tools.search_cli import main

    main()


def series() -> None:
    from tools.series_analytics import main

    main()


def retro() -> None:
    from tools.retrosynthesis import main

    main()


def schema() -> None:
    from tools.validate_csv_schema import main

    main()


def featurize() -> None:
    from tools.featurize import main

    main()


def similarity() -> None:
    from molprop_toolkit.tools.similarity_cli import main

    main()


def main(argv: Sequence[str] | None = None) -> None:
    args = list(sys.argv[1:] if argv is None else argv)
    if not args or args[0] in {"-h", "--help", "help"}:
        _print_main_help()
        return

    command = args[0]
    remainder = args[1:]

    if command == "calc":
        profile, forwarded = _parse_calc_profile(remainder)
        if profile == "baseline":
            _run_with_argv(calc_v4, "molscope calc", forwarded)
            return
        _run_with_argv(calc_v5, "molscope calc", forwarded)
        return

    dispatch: dict[str, Callable[[], None]] = {
        "prep": prep,
        "analyze": analyze,
        "report": report,
        "portal": portal,
        "learnings": learnings,
        "dashboard": dashboard,
        "sketch": sketch,
        "picklists": picklists,
        "integrate": integrate,
        "compare": compare,
        "sar": sar,
        "mmp": mmp,
        "search": search,
        "series": series,
        "similarity": similarity,
        "featurize": featurize,
        "retro": retro,
        "schema": schema,
    }

    func = dispatch.get(command)
    if func is None:
        raise SystemExit(
            f"Unknown molscope command: {command}\n\nRun 'molscope --help' to see the supported workflow commands."
        )

    _run_with_argv(func, f"molscope {command}", remainder)


if __name__ == "__main__":
    main()
