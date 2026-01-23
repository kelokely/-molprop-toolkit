#!/usr/bin/env python3
"""Protomer (pH-dependent protonation state) enumeration and selection.

Task 2B implements an enumeration-based option using Dimorphite-DL.

Dimorphite-DL generates plausible protonation variants over a user-specified
pH range. This module wraps that functionality and provides a simple selection
strategy so the calculators can record a representative protomer.

Key design points:
- Dimorphite-DL is an OPTIONAL dependency. If it is not installed and the user
  requests dimorphite mode, we raise a clear error.
- Enumeration yields SMILES; we canonicalize and de-duplicate before selection.
- Selection can be driven by an estimated target net charge (e.g., from a
  heuristic pKa model) or by preference for neutral species.

Sources:
- Dimorphite-DL project: https://github.com/durrantlab/dimorphite_dl
- Paper: Ropp PJ et al., J Cheminform 11:14 (2019). https://doi.org/10.1186/s13321-019-0336-9
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Tuple

from rdkit import Chem


@dataclass
class ProtomerConfig:
    mode: str = "none"  # none | dimorphite
    ph_min: float = 6.8
    ph_max: float = 8.0
    precision: float = 1.0
    max_variants: int = 64
    label_states: bool = False
    select: str = "closest-charge"  # closest-charge | neutral | first
    calc_on_protomer: bool = False


@dataclass
class ProtomerReport:
    mode: str
    count: int
    selected_smiles: str
    topk_smiles: str
    selected_charge: int
    target_charge: Optional[int]
    calc_on_protomer: bool
    notes: str = ""


def _try_import_dimorphite():
    try:
        from dimorphite_dl import protonate_smiles  # type: ignore

        return protonate_smiles
    except Exception as e:  # pragma: no cover
        raise ImportError(
            "Dimorphite-DL is not installed. Install with: pip install dimorphite_dl"
        ) from e


def _canonicalize_smiles_list(smiles_list: List[str]) -> List[str]:
    out: List[str] = []
    seen = set()
    for s in smiles_list:
        s = s.strip()
        if not s:
            continue
        # Dimorphite can optionally label states like "SMILES\tPROTONATED".
        if "\t" in s:
            s = s.split("\t")[0]
        m = Chem.MolFromSmiles(s)
        if m is None:
            continue
        can = Chem.MolToSmiles(m, canonical=True, isomericSmiles=True)
        if can not in seen:
            seen.add(can)
            out.append(can)
    return out


def enumerate_protomers_dimorphite(
    smiles: str,
    ph_min: float,
    ph_max: float,
    precision: float = 1.0,
    max_variants: int = 64,
    label_states: bool = False,
) -> List[str]:
    protonate_smiles = _try_import_dimorphite()

    # dimorphite_dl.protonate_smiles returns list[str]
    protomers = protonate_smiles(
        smiles,
        ph_min=float(ph_min),
        ph_max=float(ph_max),
        precision=float(precision),
        max_variants=int(max_variants),
        label_states=bool(label_states),
    )

    if not isinstance(protomers, list):
        protomers = list(protomers)

    return _canonicalize_smiles_list([str(x) for x in protomers])


def _formal_charge(smiles: str) -> Optional[int]:
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    try:
        return int(Chem.GetFormalCharge(m))
    except Exception:
        return None


def select_protomer(
    protomer_smiles: List[str],
    target_charge: Optional[int] = None,
    strategy: str = "closest-charge",
    top_k: int = 5,
) -> Tuple[str, int, str]:
    """Select a representative protomer.

    Returns (selected_smiles, selected_charge, topk_smiles).
    """

    if not protomer_smiles:
        return "", 0, ""

    topk = ";".join(protomer_smiles[: max(1, int(top_k))])

    strategy = (strategy or "first").lower()

    if strategy == "first" or target_charge is None:
        sel = protomer_smiles[0]
        ch = _formal_charge(sel) or 0
        return sel, ch, topk

    charges = [(s, _formal_charge(s)) for s in protomer_smiles]
    charges = [(s, c) for s, c in charges if c is not None]
    if not charges:
        sel = protomer_smiles[0]
        ch = _formal_charge(sel) or 0
        return sel, ch, topk

    if strategy == "neutral":
        neutral = [x for x in charges if x[1] == 0]
        if neutral:
            sel, ch = neutral[0]
            return sel, ch, topk
        # fall through to closest-charge

    # closest-charge
    ranked = sorted(
        charges,
        key=lambda x: (abs(x[1] - int(target_charge)), abs(x[1])),
    )
    sel, ch = ranked[0]
    return sel, int(ch), topk


def protomerize_smiles(
    smiles: str,
    cfg: ProtomerConfig,
    target_charge: Optional[int] = None,
    top_k: int = 5,
) -> Tuple[Optional[Chem.Mol], ProtomerReport]:
    """Enumerate protomers and return (selected_mol, report)."""

    if cfg.mode.lower() != "dimorphite":
        m = Chem.MolFromSmiles(smiles)
        rep = ProtomerReport(
            mode=cfg.mode,
            count=1,
            selected_smiles=smiles,
            topk_smiles=smiles,
            selected_charge=int(Chem.GetFormalCharge(m)) if m is not None else 0,
            target_charge=target_charge,
            calc_on_protomer=cfg.calc_on_protomer,
            notes="",
        )
        return m, rep

    protomers = enumerate_protomers_dimorphite(
        smiles,
        ph_min=cfg.ph_min,
        ph_max=cfg.ph_max,
        precision=cfg.precision,
        max_variants=cfg.max_variants,
        label_states=cfg.label_states,
    )

    if not protomers:
        m = Chem.MolFromSmiles(smiles)
        rep = ProtomerReport(
            mode="dimorphite",
            count=0,
            selected_smiles=smiles,
            topk_smiles=smiles,
            selected_charge=int(Chem.GetFormalCharge(m)) if m is not None else 0,
            target_charge=target_charge,
            calc_on_protomer=cfg.calc_on_protomer,
            notes="no_protomers_returned; fell back to input",
        )
        return m, rep

    sel_smi, sel_charge, topk = select_protomer(
        protomers,
        target_charge=target_charge,
        strategy=cfg.select,
        top_k=top_k,
    )

    m = Chem.MolFromSmiles(sel_smi) if sel_smi else Chem.MolFromSmiles(smiles)

    rep = ProtomerReport(
        mode="dimorphite",
        count=len(protomers),
        selected_smiles=sel_smi or smiles,
        topk_smiles=topk,
        selected_charge=int(sel_charge),
        target_charge=target_charge,
        calc_on_protomer=cfg.calc_on_protomer,
        notes="",
    )

    return m, rep

