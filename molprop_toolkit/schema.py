"""Schema loader utilities.

We keep a copy of the schema inside the installable package so that a normal `pip install molprop-toolkit` works without
requiring a separate `docs/` tree. Tools default to the bundled schema, while still allowing overrides via `--schema` when
working from a source checkout.
"""

from __future__ import annotations

import json
from importlib import resources
from pathlib import Path
from typing import Any, Dict, Optional


def load_schema(schema_path: Optional[str | Path] = None) -> Dict[str, Any]:
    """Load a schema JSON.

    Precedence:
    1) schema_path if provided and exists
    2) bundled package schema (molprop_toolkit/resources/schema.json)
    3) docs/schema.json when running from a source checkout
    """

    if schema_path is not None:
        p = Path(schema_path)
        if p.exists():
            return json.loads(p.read_text(encoding="utf-8"))

    # Bundled resource
    try:
        txt = (resources.files("molprop_toolkit.resources") / "schema.json").read_text(
            encoding="utf-8"
        )
        return json.loads(txt)
    except Exception:
        pass

    # Source checkout fallback
    p2 = Path(__file__).resolve().parents[1] / "docs" / "schema.json"
    return json.loads(p2.read_text(encoding="utf-8"))
