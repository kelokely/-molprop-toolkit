#!/usr/bin/env python3
"""RDKit structure preparation utilities.

This module centralizes structure preparation steps that should happen before
property calculation.

Design goals:
- Keep the original input SMILES for traceability.
- Produce a standardized parent structure (salt/solvent stripped when possible).
- Canonicalize tautomer form (optional).
- Assign stereochemistry and report unassigned stereocenters/bonds.

Important limitations:
- This module does NOT perform true pH-dependent protonation or rigorous pKa
  prediction. That is handled separately (or via optional external tools).
- RDKit standardization tools aim for a reasonable "parent" form, but may not
  match the preferred state for a specific project.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

from rdkit import Chem

try:
    # RDKit standardization utilities
    from rdkit.Chem.MolStandardize import rdMolStandardize
except Exception as e:  # pragma: no cover
    rdMolStandardize = None  # type: ignore


@dataclass
class PrepConfig:
    """Configuration for structure preparation."""

    strip_fragments: bool = True
    cleanup: bool = True
    uncharge: bool = True
    reionize: bool = True
    canonicalize_tautomer: bool = True


@dataclass
class PrepReport:
    """What happened during preparation."""

    changed: bool
    salts_stripped: bool
    uncharged: bool
    reionized: bool
    tautomer_changed: bool
    stereo_status: str
    unassigned_chiral_centers: int
    unassigned_stereo_bonds: int
    notes: str


def _safe_mol_to_smiles(mol: Chem.Mol) -> str:
    return Chem.MolToSmiles(mol, canonical=True)


def stereo_summary(mol: Chem.Mol) -> Tuple[str, int, int]:
    """Return (status, n_unassigned_centers, n_unassigned_double_bonds)."""

    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)

    centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    n_unassigned_centers = sum(1 for _, tag in centers if tag == "?")

    n_unassigned_bonds = 0
    try:
        # FindPotentialStereo reports possible stereo elements and whether specified.
        stereos = Chem.FindPotentialStereo(mol)
        for si in stereos:
            # Depending on RDKit version, attributes are slightly different.
            # We defensively check for 'specified' flag.
            specified = getattr(si, "specified", True)
            stype = getattr(si, "type", None)
            if stype is None:
                continue
            if str(stype).lower().endswith("bond_double") and not specified:
                n_unassigned_bonds += 1
    except Exception:
        n_unassigned_bonds = 0

    if len(centers) == 0 and n_unassigned_bonds == 0:
        status = "None"
    elif n_unassigned_centers == 0 and n_unassigned_bonds == 0:
        status = "Defined"
    else:
        status = "Undefined"

    return status, n_unassigned_centers, n_unassigned_bonds


def prepare_mol(mol: Chem.Mol, cfg: Optional[PrepConfig] = None) -> Tuple[Chem.Mol, PrepReport]:
    """Prepare an RDKit molecule and return (prepared_mol, report)."""

    if cfg is None:
        cfg = PrepConfig()

    in_smiles = _safe_mol_to_smiles(mol)

    salts_stripped = False
    uncharged = False
    reionized = False
    tautomer_changed = False
    notes = []

    m = Chem.Mol(mol)

    if rdMolStandardize is None:
        # RDKit built without MolStandardize
        status, uc, ub = stereo_summary(m)
        rep = PrepReport(
            changed=False,
            salts_stripped=False,
            uncharged=False,
            reionized=False,
            tautomer_changed=False,
            stereo_status=status,
            unassigned_chiral_centers=uc,
            unassigned_stereo_bonds=ub,
            notes="MolStandardize unavailable; stereo assigned only",
        )
        return m, rep

    if cfg.cleanup:
        try:
            m2 = rdMolStandardize.Cleanup(m)
            m = m2
        except Exception:
            notes.append("cleanup_failed")

    if cfg.strip_fragments:
        try:
            m2 = rdMolStandardize.FragmentParent(m)
            if _safe_mol_to_smiles(m2) != _safe_mol_to_smiles(m):
                salts_stripped = True
            m = m2
        except Exception:
            notes.append("fragment_parent_failed")

    if cfg.uncharge:
        try:
            uncharger = rdMolStandardize.Uncharger()
            m2 = uncharger.uncharge(m)
            if _safe_mol_to_smiles(m2) != _safe_mol_to_smiles(m):
                uncharged = True
            m = m2
        except Exception:
            notes.append("uncharge_failed")

    if cfg.reionize:
        try:
            reionizer = rdMolStandardize.Reionizer()
            m2 = reionizer.reionize(m)
            if _safe_mol_to_smiles(m2) != _safe_mol_to_smiles(m):
                reionized = True
            m = m2
        except Exception:
            notes.append("reionize_failed")

    if cfg.canonicalize_tautomer:
        try:
            te = rdMolStandardize.TautomerEnumerator()
            m2 = te.Canonicalize(m)
            if _safe_mol_to_smiles(m2) != _safe_mol_to_smiles(m):
                tautomer_changed = True
            m = m2
        except Exception:
            notes.append("tautomer_failed")

    status, uc, ub = stereo_summary(m)

    out_smiles = _safe_mol_to_smiles(m)
    changed = out_smiles != in_smiles

    rep = PrepReport(
        changed=changed,
        salts_stripped=salts_stripped,
        uncharged=uncharged,
        reionized=reionized,
        tautomer_changed=tautomer_changed,
        stereo_status=status,
        unassigned_chiral_centers=uc,
        unassigned_stereo_bonds=ub,
        notes=",".join(notes) if notes else "",
    )

    return m, rep

