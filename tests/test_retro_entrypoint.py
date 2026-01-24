from __future__ import annotations

import molprop_toolkit.cli as cli


def test_cli_has_retro_entrypoint() -> None:
    assert hasattr(cli, "retro")

