from __future__ import annotations

import sys

import molprop_toolkit.cli as cli


def test_cli_has_molscope_entrypoint() -> None:
    assert hasattr(cli, "main")


def test_molscope_calc_baseline_dispatches_without_profile_arg(monkeypatch) -> None:
    captured: dict[str, list[str]] = {}

    def fake_calc_v4() -> None:
        captured["argv"] = sys.argv[:]

    monkeypatch.setattr(cli, "calc_v4", fake_calc_v4)
    cli.main(["calc", "--profile", "baseline", "examples/example.smi", "-o", "results.csv"])

    assert captured["argv"] == [
        "molscope calc",
        "examples/example.smi",
        "-o",
        "results.csv",
    ]


def test_molscope_picklists_dispatches_subcommand(monkeypatch) -> None:
    captured: dict[str, list[str]] = {}

    def fake_picklists() -> None:
        captured["argv"] = sys.argv[:]

    monkeypatch.setattr(cli, "picklists", fake_picklists)
    cli.main(["picklists", "results.csv", "--html"])

    assert captured["argv"] == [
        "molscope picklists",
        "results.csv",
        "--html",
    ]


def test_molscope_portal_dispatches_subcommand(monkeypatch) -> None:
    captured: dict[str, list[str]] = {}

    def fake_portal() -> None:
        captured["argv"] = sys.argv[:]

    monkeypatch.setattr(cli, "portal", fake_portal)
    cli.main(["portal", "--report-dir", "reports/demo", "--picklists-dir", "picklists/demo"])

    assert captured["argv"] == [
        "molscope portal",
        "--report-dir",
        "reports/demo",
        "--picklists-dir",
        "picklists/demo",
    ]


def test_molscope_learnings_dispatches_subcommand(monkeypatch) -> None:
    captured: dict[str, list[str]] = {}

    def fake_learnings() -> None:
        captured["argv"] = sys.argv[:]

    monkeypatch.setattr(cli, "learnings", fake_learnings)
    cli.main(
        [
            "learnings",
            "--compare-dir",
            "compare/demo",
            "--picklists-dir",
            "picklists/demo",
            "-o",
            "learnings/demo",
        ]
    )

    assert captured["argv"] == [
        "molscope learnings",
        "--compare-dir",
        "compare/demo",
        "--picklists-dir",
        "picklists/demo",
        "-o",
        "learnings/demo",
    ]


def test_molscope_dashboard_dispatches_subcommand(monkeypatch) -> None:
    captured: dict[str, list[str]] = {}

    def fake_dashboard() -> None:
        captured["argv"] = sys.argv[:]

    monkeypatch.setattr(cli, "dashboard", fake_dashboard)
    cli.main(
        [
            "dashboard",
            "--learnings-dir",
            "learnings/round5",
            "-o",
            "dashboard/demo",
        ]
    )

    assert captured["argv"] == [
        "molscope dashboard",
        "--learnings-dir",
        "learnings/round5",
        "-o",
        "dashboard/demo",
    ]


def test_molscope_integrate_dispatches_subcommand(monkeypatch) -> None:
    captured: dict[str, list[str]] = {}

    def fake_integrate() -> None:
        captured["argv"] = sys.argv[:]

    monkeypatch.setattr(cli, "integrate", fake_integrate)
    cli.main(["integrate", "results.csv", "assays.csv", "-o", "integrated.csv"])

    assert captured["argv"] == [
        "molscope integrate",
        "results.csv",
        "assays.csv",
        "-o",
        "integrated.csv",
    ]


def test_molscope_compare_dispatches_subcommand(monkeypatch) -> None:
    captured: dict[str, list[str]] = {}

    def fake_compare() -> None:
        captured["argv"] = sys.argv[:]

    monkeypatch.setattr(cli, "compare", fake_compare)
    cli.main(["compare", "round4.csv", "round5.csv", "-o", "compare-dir"])

    assert captured["argv"] == [
        "molscope compare",
        "round4.csv",
        "round5.csv",
        "-o",
        "compare-dir",
    ]


def test_molscope_sar_dispatches_subcommand(monkeypatch) -> None:
    captured: dict[str, list[str]] = {}

    def fake_sar() -> None:
        captured["argv"] = sys.argv[:]

    monkeypatch.setattr(cli, "sar", fake_sar)
    cli.main(["sar", "results_series.csv", "-o", "sar-dir"])

    assert captured["argv"] == [
        "molscope sar",
        "results_series.csv",
        "-o",
        "sar-dir",
    ]


def test_molscope_mmp_dispatches_subcommand(monkeypatch) -> None:
    captured: dict[str, list[str]] = {}

    def fake_mmp() -> None:
        captured["argv"] = sys.argv[:]

    monkeypatch.setattr(cli, "mmp", fake_mmp)
    cli.main(["mmp", "results.csv", "-o", "mmp-dir"])

    assert captured["argv"] == [
        "molscope mmp",
        "results.csv",
        "-o",
        "mmp-dir",
    ]


def test_molscope_search_dispatches_subcommand(monkeypatch) -> None:
    captured: dict[str, list[str]] = {}

    def fake_search() -> None:
        captured["argv"] = sys.argv[:]

    monkeypatch.setattr(cli, "search", fake_search)
    cli.main(["search", "results.csv", "--smarts", "c1ccccc1", "-o", "search-dir"])

    assert captured["argv"] == [
        "molscope search",
        "results.csv",
        "--smarts",
        "c1ccccc1",
        "-o",
        "search-dir",
    ]


def test_molscope_retro_dispatches_subcommand(monkeypatch) -> None:
    captured: dict[str, list[str]] = {}

    def fake_retro() -> None:
        captured["argv"] = sys.argv[:]

    monkeypatch.setattr(cli, "retro", fake_retro)
    cli.main(
        [
            "retro",
            "results.csv",
            "--config",
            "aizynth_data/config.yml",
            "--rerank-policy",
            "stock-first",
        ]
    )

    assert captured["argv"] == [
        "molscope retro",
        "results.csv",
        "--config",
        "aizynth_data/config.yml",
        "--rerank-policy",
        "stock-first",
    ]


def test_molscope_prep_dispatches_subcommand(monkeypatch) -> None:
    captured: dict[str, list[str]] = {}

    def fake_prep() -> None:
        captured["argv"] = sys.argv[:]

    monkeypatch.setattr(cli, "prep", fake_prep)
    cli.main(["prep", "examples/example.smi", "-o", "prepared.csv", "--preset", "audit"])

    assert captured["argv"] == [
        "molscope prep",
        "examples/example.smi",
        "-o",
        "prepared.csv",
        "--preset",
        "audit",
    ]


def test_molscope_schema_dispatches_subcommand(monkeypatch) -> None:
    captured: dict[str, list[str]] = {}

    def fake_schema() -> None:
        captured["argv"] = sys.argv[:]

    monkeypatch.setattr(cli, "schema", fake_schema)
    cli.main(["schema", "results.csv", "--require-category", "developability"])

    assert captured["argv"] == [
        "molscope schema",
        "results.csv",
        "--require-category",
        "developability",
    ]
