"""Console-script entrypoints.

These wrappers keep existing workflows intact.

The original scripts remain runnable as `python calculators/...` and `python tools/...` while the package also exposes
installable console scripts that call the same underlying `main()` functions.

Design note: these entrypoints intentionally avoid re-implementing argument parsing; they delegate to the existing script CLIs.
"""

from __future__ import annotations


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


def sketch() -> None:
    from tools.sketch_structures import main

    main()


def picklists() -> None:
    from tools.picklists import main

    main()


def series() -> None:
    from tools.series_analytics import main

    main()


def retro() -> None:
    from tools.retrosynthesis import main

    main()


def featurize() -> None:
    from tools.featurize import main

    main()

