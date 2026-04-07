#!/usr/bin/env python3
"""2D structure depiction utilities (Task 5).

This module provides lightweight helpers for converting SMILES to RDKit mols
and rendering 2D depictions to SVG/PNG.

We keep this separate from calculators/analyzers so it can be reused by a
future report generator (HTML/PDF) without duplicating RDKit drawing code.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Optional, Tuple

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdDepictor


def sanitize_filename(name: str, max_len: int = 120) -> str:
    name = (name or "").strip() or "compound"
    name = re.sub(r"[^A-Za-z0-9._-]+", "_", name)
    name = name.strip("._-") or "compound"
    return name[:max_len]


def mol_from_smiles(smiles: str) -> Optional[Chem.Mol]:
    if not smiles:
        return None
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    try:
        rdDepictor.Compute2DCoords(m)
    except Exception:
        pass
    return m


def draw_svg(mol: Chem.Mol, size: Tuple[int, int] = (280, 200), legend: str = "") -> str:
    w, h = int(size[0]), int(size[1])
    d2d = rdMolDraw2D.MolDraw2DSVG(w, h)
    opts = d2d.drawOptions()
    opts.addStereoAnnotation = True
    opts.addAtomIndices = False
    d2d.DrawMolecule(mol, legend=legend)
    d2d.FinishDrawing()
    svg = d2d.GetDrawingText()
    return svg


def save_svg(path: str | Path, mol: Chem.Mol, size: Tuple[int, int] = (280, 200), legend: str = "") -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    svg = draw_svg(mol, size=size, legend=legend)
    p.write_text(svg, encoding="utf-8")


def save_png(path: str | Path, mol: Chem.Mol, size: Tuple[int, int] = (280, 200), legend: str = "") -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    img = Draw.MolToImage(mol, size=(int(size[0]), int(size[1])), legend=legend)
    img.save(str(p))

