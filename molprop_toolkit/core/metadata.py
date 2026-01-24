"""Run metadata helpers.

Task 3: Every command that produces an output table should emit a small, machine-readable
metadata JSON artifact capturing provenance (input hash, parameters, versions).

Convention:
- for an output table path like `results.csv` or `results.parquet`, write a sidecar
  file next to it named `results.metadata.json`.

These helpers are intentionally dependency-light.
"""

from __future__ import annotations

import json
import os
import sys
from datetime import datetime
from importlib import metadata as importlib_metadata
from pathlib import Path
from typing import Any, Dict, Optional


def get_toolkit_version() -> str:
    """Return installed package version if available, else 'unknown'."""

    try:
        return importlib_metadata.version("molprop-toolkit")
    except Exception:
        return "unknown"


def get_rdkit_version() -> str:
    try:
        from rdkit import rdBase  # type: ignore

        return str(getattr(rdBase, "rdkitVersion", "unknown"))
    except Exception:
        return "not_installed"


def sha256_file(path: Path, *, max_bytes: int = 200 * 1024 * 1024) -> Optional[str]:
    """Compute SHA256 for a file, returning None if too large or unreadable."""

    try:
        st = path.stat()
        if st.st_size > max_bytes:
            return None

        import hashlib

        h = hashlib.sha256()
        with path.open("rb") as f:
            for chunk in iter(lambda: f.read(1024 * 1024), b""):
                h.update(chunk)
        return h.hexdigest()
    except Exception:
        return None


def metadata_sidecar_path(output_table_path: str | Path) -> Path:
    p = Path(output_table_path)
    return p.with_name(f"{p.stem}.metadata.json")


def write_run_metadata(
    *,
    tool: str,
    output_table_path: str | Path,
    input_path: Optional[str | Path] = None,
    smiles_column: Optional[str] = None,
    parameters: Optional[Dict[str, Any]] = None,
    artifacts: Optional[Dict[str, Any]] = None,
    notes: Optional[str] = None,
) -> Path:
    """Write a standardized run metadata JSON sidecar.

    This is meant to be called by CLIs after they have written the output artifact.
    """

    out_p = Path(output_table_path)
    in_p = Path(input_path) if input_path else None

    payload: Dict[str, Any] = {
        "tool": str(tool),
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "cwd": os.getcwd(),
        "argv": list(sys.argv),
        "versions": {
            "molprop_toolkit": get_toolkit_version(),
            "python": sys.version.split()[0],
            "rdkit": get_rdkit_version(),
        },
        "structure_basis": {
            "smiles_column": smiles_column,
            "policy": "auto (prefers Calc_Canonical_SMILES when present)",
        },
        "input": None,
        "output": {
            "path": str(out_p.resolve()),
            "name": out_p.name,
            "format": out_p.suffix.lower().lstrip(".") or "unknown",
            "size_bytes": int(out_p.stat().st_size) if out_p.exists() else None,
        },
        "parameters": parameters or {},
        "artifacts": artifacts or {},
    }

    if notes:
        payload["notes"] = str(notes)

    if in_p is not None:
        payload["input"] = {
            "path": str(in_p.resolve()),
            "name": in_p.name,
            "sha256": sha256_file(in_p),
            "size_bytes": int(in_p.stat().st_size) if in_p.exists() else None,
            "mtime": datetime.fromtimestamp(in_p.stat().st_mtime).isoformat(timespec="seconds") if in_p.exists() else None,
        }

    sidecar = metadata_sidecar_path(out_p)
    sidecar.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return sidecar

