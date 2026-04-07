"""Preparation-first CLI for MolScope.

This command promotes the documented ligand-preparation modes into a real
workflow step. It writes a preparation-focused table plus a summary JSON so
teams can lock a structure-handling policy before descriptor calculation.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Optional

import pandas as pd
from molprop_toolkit.core.io import write_table
from molprop_toolkit.core.metadata import write_run_metadata

PREP_IMPORT_ERROR: ModuleNotFoundError | None = None

try:
    from rdkit import Chem
    from rdkit.Chem import Crippen

    from calculators.ionization import estimate_net_charge, ionization_features
    from calculators.mpo_v5 import estimate_pka, parse_smiles_file
    from calculators.prep import PrepConfig, PrepReport, prepare_mol, stereo_summary
    from calculators.protomers import ProtomerConfig, protomerize_smiles
    from calculators.stereo import StereoConfig, apply_stereo
    from calculators.tautomers import TautomerConfig, select_representative_tautomer
except ModuleNotFoundError as exc:
    if exc.name and exc.name.startswith("rdkit"):
        PREP_IMPORT_ERROR = exc
    else:  # pragma: no cover - unrelated import errors should still surface
        raise


PREP_PRESETS: dict[str, dict[str, Any]] = {
    "triage": {
        "use_prep": True,
        "uncharge": True,
        "tautomer_mode": "prep-canonical",
        "tautomer_max": 64,
        "tautomer_topk": 5,
        "stereo_mode": "keep",
        "stereo_max": 32,
        "stereo_topk": 5,
        "stereo_select": "canonical",
        "ionization": "heuristic",
        "ph": 7.4,
        "ph_min": None,
        "ph_max": None,
        "dimorphite_precision": 1.0,
        "dimorphite_max_variants": 64,
        "protomer_select": "closest-charge",
        "calc_on_protomer": False,
    },
    "audit": {
        "use_prep": True,
        "uncharge": True,
        "tautomer_mode": "enumerate",
        "tautomer_max": 64,
        "tautomer_topk": 5,
        "stereo_mode": "enumerate",
        "stereo_max": 32,
        "stereo_topk": 5,
        "stereo_select": "canonical",
        "ionization": "heuristic",
        "ph": 7.4,
        "ph_min": None,
        "ph_max": None,
        "dimorphite_precision": 1.0,
        "dimorphite_max_variants": 64,
        "protomer_select": "closest-charge",
        "calc_on_protomer": False,
    },
    "3d": {
        "use_prep": True,
        "uncharge": True,
        "tautomer_mode": "prep-canonical",
        "tautomer_max": 64,
        "tautomer_topk": 5,
        "stereo_mode": "keep",
        "stereo_max": 32,
        "stereo_topk": 5,
        "stereo_select": "canonical",
        "ionization": "dimorphite",
        "ph": 7.4,
        "ph_min": None,
        "ph_max": None,
        "dimorphite_precision": 1.0,
        "dimorphite_max_variants": 64,
        "protomer_select": "closest-charge",
        "calc_on_protomer": True,
    },
    "external": {
        "use_prep": False,
        "uncharge": True,
        "tautomer_mode": "none",
        "tautomer_max": 64,
        "tautomer_topk": 5,
        "stereo_mode": "keep",
        "stereo_max": 32,
        "stereo_topk": 5,
        "stereo_select": "canonical",
        "ionization": "none",
        "ph": 7.4,
        "ph_min": None,
        "ph_max": None,
        "dimorphite_precision": 1.0,
        "dimorphite_max_variants": 64,
        "protomer_select": "closest-charge",
        "calc_on_protomer": False,
    },
}


def _require_prep_dependencies() -> None:
    if PREP_IMPORT_ERROR is None:
        return
    raise SystemExit(
        "molscope prep requires RDKit. Install RDKit, preferably via conda-forge, "
        "before running this command."
    ) from PREP_IMPORT_ERROR


def _summary_path_for_output(output_path: str) -> Path:
    out = Path(output_path)
    return out.with_name(f"{out.stem}.prep-summary.json")


def _safe_num(value: Optional[float]) -> Optional[float]:
    if value is None:
        return None
    return round(float(value), 3)


def _fallback_prep_report(mol: Chem.Mol, *, note: str) -> PrepReport:
    status, unassigned_chiral, unassigned_bonds = stereo_summary(Chem.Mol(mol))
    return PrepReport(
        changed=False,
        salts_stripped=False,
        uncharged=False,
        reionized=False,
        tautomer_changed=False,
        stereo_status=status,
        unassigned_chiral_centers=unassigned_chiral,
        unassigned_stereo_bonds=unassigned_bonds,
        notes=note,
    )


def _build_effective_settings(args: argparse.Namespace) -> dict[str, Any]:
    settings = dict(PREP_PRESETS[args.preset])

    if args.ph is not None:
        settings["ph"] = float(args.ph)
    if args.no_prep:
        settings["use_prep"] = False
    if args.keep_charges:
        settings["uncharge"] = False
    if args.tautomer_mode is not None:
        settings["tautomer_mode"] = str(args.tautomer_mode)
    if args.tautomer_max is not None:
        settings["tautomer_max"] = int(args.tautomer_max)
    if args.tautomer_topk is not None:
        settings["tautomer_topk"] = int(args.tautomer_topk)
    if args.stereo_mode is not None:
        settings["stereo_mode"] = str(args.stereo_mode)
    if args.stereo_max is not None:
        settings["stereo_max"] = int(args.stereo_max)
    if args.stereo_topk is not None:
        settings["stereo_topk"] = int(args.stereo_topk)
    if args.stereo_select is not None:
        settings["stereo_select"] = str(args.stereo_select)
    if args.ionization is not None:
        settings["ionization"] = str(args.ionization)
    if args.ph_min is not None:
        settings["ph_min"] = float(args.ph_min)
    if args.ph_max is not None:
        settings["ph_max"] = float(args.ph_max)
    if args.dimorphite_precision is not None:
        settings["dimorphite_precision"] = float(args.dimorphite_precision)
    if args.dimorphite_max_variants is not None:
        settings["dimorphite_max_variants"] = int(args.dimorphite_max_variants)
    if args.protomer_select is not None:
        settings["protomer_select"] = str(args.protomer_select)
    if args.calc_on_protomer:
        settings["calc_on_protomer"] = True

    if not settings["use_prep"] and settings["tautomer_mode"] == "prep-canonical":
        settings["tautomer_mode"] = "none"

    if settings["ionization"] != "dimorphite":
        settings["calc_on_protomer"] = False

    if settings["ionization"] == "dimorphite":
        if settings["ph_min"] is None:
            settings["ph_min"] = settings["ph"] - 0.5
        if settings["ph_max"] is None:
            settings["ph_max"] = settings["ph"] + 0.5

    return settings


def _build_configs(
    settings: dict[str, Any],
) -> tuple[Optional[PrepConfig], TautomerConfig, StereoConfig, ProtomerConfig]:
    prep_cfg: Optional[PrepConfig]
    if settings["use_prep"]:
        prep_cfg = PrepConfig(
            canonicalize_tautomer=(settings["tautomer_mode"] == "prep-canonical"),
            uncharge=bool(settings["uncharge"]),
        )
    else:
        prep_cfg = None

    taut_cfg = TautomerConfig(
        mode=str(settings["tautomer_mode"]),
        max_variants=int(settings["tautomer_max"]),
        top_k=int(settings["tautomer_topk"]),
    )
    stereo_cfg = StereoConfig(
        mode=str(settings["stereo_mode"]),
        max_isomers=int(settings["stereo_max"]),
        top_k=int(settings["stereo_topk"]),
        select=str(settings["stereo_select"]),
    )

    if settings["ionization"] == "dimorphite":
        prot_cfg = ProtomerConfig(
            mode="dimorphite",
            ph_min=float(settings["ph_min"]),
            ph_max=float(settings["ph_max"]),
            precision=float(settings["dimorphite_precision"]),
            max_variants=int(settings["dimorphite_max_variants"]),
            select=str(settings["protomer_select"]),
            calc_on_protomer=bool(settings["calc_on_protomer"]),
        )
    else:
        prot_cfg = ProtomerConfig(mode="none", calc_on_protomer=False)

    return prep_cfg, taut_cfg, stereo_cfg, prot_cfg


def _prepare_record(
    compound_id: str,
    smiles: str,
    *,
    preset: str,
    ph: float,
    ionization_mode: str,
    prep_cfg: Optional[PrepConfig],
    tautomer_cfg: TautomerConfig,
    stereo_cfg: StereoConfig,
    protomer_cfg: ProtomerConfig,
) -> Optional[dict[str, object]]:
    mol_in = Chem.MolFromSmiles(smiles)
    if mol_in is None:
        return None

    input_can = Chem.MolToSmiles(mol_in, canonical=True, isomericSmiles=True)

    if prep_cfg is None:
        mol = Chem.Mol(mol_in)
        prep_report = _fallback_prep_report(mol, note="prep_disabled")
    else:
        mol, prep_report = prepare_mol(mol_in, prep_cfg)

    prep_can = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)

    calc_mol = mol
    tautomer_report = None
    if tautomer_cfg.mode == "enumerate":
        calc_mol, tautomer_report = select_representative_tautomer(
            mol,
            max_variants=tautomer_cfg.max_variants,
            top_k=tautomer_cfg.top_k,
        )

    calc_mol, stereo_report = apply_stereo(calc_mol, stereo_cfg)
    calc_base_smiles = Chem.MolToSmiles(
        calc_mol,
        canonical=True,
        isomericSmiles=True,
    )

    basic_pka, acidic_pka, pka_type = estimate_pka(calc_mol)

    protomer_report = None
    if protomer_cfg.mode.lower() == "dimorphite":
        target_charge = int(round(estimate_net_charge(basic_pka, acidic_pka, ph)))
        prot_mol, protomer_report = protomerize_smiles(
            calc_base_smiles,
            protomer_cfg,
            target_charge=target_charge,
            top_k=5,
        )
        if protomer_cfg.calc_on_protomer and prot_mol is not None:
            calc_mol = prot_mol

    calc_can = Chem.MolToSmiles(calc_mol, canonical=True, isomericSmiles=True)
    logp = float(Crippen.MolLogP(calc_mol))

    out: dict[str, object] = {
        "Compound_ID": compound_id,
        "SMILES": smiles,
        "Input_Canonical_SMILES": input_can,
        "Canonical_SMILES": prep_can,
        "Calc_Base_SMILES": calc_base_smiles,
        "Calc_Canonical_SMILES": calc_can,
        "Prep_Changed": prep_report.changed,
        "Prep_SaltsStripped": prep_report.salts_stripped,
        "Prep_Uncharged": prep_report.uncharged,
        "Prep_Reionized": prep_report.reionized,
        "Prep_TautomerChanged": prep_report.tautomer_changed,
        "Prep_StereoStatus": prep_report.stereo_status,
        "Prep_UnassignedChiral": prep_report.unassigned_chiral_centers,
        "Prep_UnassignedStereoBonds": prep_report.unassigned_stereo_bonds,
        "Prep_Notes": prep_report.notes,
        "Tautomer_Mode": tautomer_cfg.mode,
        "Stereo_Mode": stereo_cfg.mode,
        "pKa_Basic": _safe_num(basic_pka),
        "pKa_Acidic": _safe_num(acidic_pka),
        "pKa_Type": pka_type,
        "Ion_Mode": ionization_mode,
    }

    if tautomer_cfg.mode == "enumerate" and tautomer_report is not None:
        out["Tautomer_Count"] = tautomer_report.count
        out["Tautomer_Selected_SMILES"] = tautomer_report.selected_smiles
        out["Tautomer_TopK_SMILES"] = tautomer_report.topk_smiles
        out["Tautomer_Changed"] = tautomer_report.changed
    else:
        out["Tautomer_Count"] = 1
        out["Tautomer_Selected_SMILES"] = prep_can
        out["Tautomer_TopK_SMILES"] = prep_can
        out["Tautomer_Changed"] = False

    out["Stereo_IsomerCount"] = stereo_report.count
    out["Stereo_Selected_SMILES"] = stereo_report.selected_smiles
    out["Stereo_TopK_SMILES"] = stereo_report.topk_smiles
    out["Stereo_Changed"] = stereo_report.changed
    out["Stereo_Notes"] = stereo_report.notes

    out["Protomer_Mode"] = protomer_cfg.mode
    out["Protomer_pH_Min"] = (
        float(protomer_cfg.ph_min)
        if protomer_cfg.mode == "dimorphite"
        else None
    )
    out["Protomer_pH_Max"] = (
        float(protomer_cfg.ph_max)
        if protomer_cfg.mode == "dimorphite"
        else None
    )
    out["Protomer_Precision"] = (
        float(protomer_cfg.precision)
        if protomer_cfg.mode == "dimorphite"
        else None
    )
    out["Protomer_MaxVariants"] = (
        int(protomer_cfg.max_variants)
        if protomer_cfg.mode == "dimorphite"
        else None
    )
    out["Protomer_CalcOnProtomer"] = bool(protomer_cfg.calc_on_protomer)

    if protomer_report is not None:
        out["Protomer_Count"] = protomer_report.count
        out["Protomer_Selected_SMILES"] = protomer_report.selected_smiles
        out["Protomer_Selected_Charge"] = protomer_report.selected_charge
        out["Protomer_Target_Charge"] = protomer_report.target_charge
        out["Protomer_TopK_SMILES"] = protomer_report.topk_smiles
        out["Protomer_Notes"] = protomer_report.notes
    else:
        out["Protomer_Count"] = 1
        out["Protomer_Selected_SMILES"] = calc_base_smiles
        out["Protomer_Selected_Charge"] = None
        out["Protomer_Target_Charge"] = None
        out["Protomer_TopK_SMILES"] = calc_base_smiles
        out["Protomer_Notes"] = ""

    if ionization_mode in {"heuristic", "dimorphite"}:
        out.update(ionization_features(logp, basic_pka, acidic_pka, ph=ph))

    if preset:
        out["Prep_Notes"] = ",".join(
            part for part in [str(out["Prep_Notes"]).strip(","), f"preset={preset}"] if part
        )

    return out


def _build_summary(
    rows: list[dict[str, object]],
    *,
    preset: str,
    input_path: str,
    output_path: str,
    summary_path: str,
    settings: dict[str, Any],
) -> dict[str, object]:
    stereo_status_counts = Counter(str(row.get("Prep_StereoStatus")) for row in rows)
    ion_state_counts = Counter(
        str(row.get("Ion_State"))
        for row in rows
        if row.get("Ion_State") not in (None, "")
    )
    pka_type_counts = Counter(
        str(row.get("pKa_Type"))
        for row in rows
        if row.get("pKa_Type") not in (None, "")
    )

    return {
        "tool": "molscope prep",
        "preset": preset,
        "input_path": input_path,
        "output_table_path": output_path,
        "summary_path": summary_path,
        "prepared_count": len(rows),
        "counts": {
            "prep_changed": sum(bool(row.get("Prep_Changed")) for row in rows),
            "salts_stripped": sum(bool(row.get("Prep_SaltsStripped")) for row in rows),
            "uncharged": sum(bool(row.get("Prep_Uncharged")) for row in rows),
            "reionized": sum(bool(row.get("Prep_Reionized")) for row in rows),
            "tautomer_changed": sum(bool(row.get("Tautomer_Changed")) for row in rows),
            "stereo_changed": sum(bool(row.get("Stereo_Changed")) for row in rows),
            "protomerized": sum(
                1 for row in rows if str(row.get("Protomer_Mode")) == "dimorphite"
            ),
            "calc_on_protomer": sum(
                bool(row.get("Protomer_CalcOnProtomer")) for row in rows
            ),
            "undefined_stereo": sum(
                1
                for row in rows
                if str(row.get("Prep_StereoStatus")) == "Undefined"
            ),
        },
        "totals": {
            "unassigned_chiral_centers": sum(
                int(row.get("Prep_UnassignedChiral") or 0) for row in rows
            ),
            "unassigned_stereo_bonds": sum(
                int(row.get("Prep_UnassignedStereoBonds") or 0) for row in rows
            ),
        },
        "stereo_status_counts": dict(sorted(stereo_status_counts.items())),
        "ion_state_counts": dict(sorted(ion_state_counts.items())),
        "pka_type_counts": dict(sorted(pka_type_counts.items())),
        "parameters": {
            "preset": preset,
            "ph": settings["ph"],
            "ionization": settings["ionization"],
            "tautomer_mode": settings["tautomer_mode"],
            "stereo_mode": settings["stereo_mode"],
            "use_prep": settings["use_prep"],
            "calc_on_protomer": settings["calc_on_protomer"],
        },
    }


def build_arg_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description="MolScope preparation workflow: standardize structures, resolve tautomer/stereo choices, and write a prep-focused table.",
    )
    ap.add_argument("input", help="Input SMILES file")
    ap.add_argument(
        "-o",
        "--output",
        default=None,
        help="Output prep table (CSV/TSV/Parquet inferred from extension)",
    )
    ap.add_argument(
        "--preset",
        choices=sorted(PREP_PRESETS),
        default="triage",
        help="Preparation preset: triage, audit, 3d, or external",
    )
    ap.add_argument(
        "--summary",
        default=None,
        help="Optional summary JSON path (default: alongside output table)",
    )
    ap.add_argument("--ph", type=float, default=None, help="Target pH for heuristic ionization features or dimorphite defaults")
    ap.add_argument("--no-prep", action="store_true", default=False, help="Disable RDKit preparation even if the preset enables it")
    ap.add_argument("--keep-charges", action="store_true", default=False, help="Skip uncharging during preparation")
    ap.add_argument(
        "--tautomer-mode",
        choices=["prep-canonical", "none", "enumerate"],
        default=None,
        help="Override tautomer handling mode",
    )
    ap.add_argument("--tautomer-max", type=int, default=None, help="Max tautomers when enumerating")
    ap.add_argument("--tautomer-topk", type=int, default=None, help="Top-k tautomer SMILES to record")
    ap.add_argument(
        "--stereo-mode",
        choices=["keep", "strip", "enumerate"],
        default=None,
        help="Override stereochemistry handling mode",
    )
    ap.add_argument("--stereo-max", type=int, default=None, help="Max stereoisomers when enumerating")
    ap.add_argument("--stereo-topk", type=int, default=None, help="Top-k stereoisomer SMILES to record")
    ap.add_argument(
        "--stereo-select",
        choices=["canonical", "first"],
        default=None,
        help="How to select a representative stereoisomer",
    )
    ap.add_argument(
        "--ionization",
        choices=["none", "heuristic", "dimorphite"],
        default=None,
        help="Override ionization handling mode",
    )
    ap.add_argument("--ph-min", type=float, default=None, help="Dimorphite minimum pH")
    ap.add_argument("--ph-max", type=float, default=None, help="Dimorphite maximum pH")
    ap.add_argument("--dimorphite-precision", type=float, default=None, help="Dimorphite precision")
    ap.add_argument("--dimorphite-max-variants", type=int, default=None, help="Max Dimorphite variants")
    ap.add_argument(
        "--protomer-select",
        choices=["closest-charge", "neutral", "first"],
        default=None,
        help="How to select a representative protomer",
    )
    ap.add_argument(
        "--calc-on-protomer",
        action="store_true",
        default=False,
        help="When using dimorphite, set Calc_Canonical_SMILES to the selected protomer",
    )
    return ap


def main(argv: Optional[list[str]] = None) -> None:
    ap = build_arg_parser()
    args = ap.parse_args(argv)
    _require_prep_dependencies()

    output_path = args.output or f"{Path(args.input).stem}_prepared.csv"
    summary_path = args.summary or str(_summary_path_for_output(output_path))

    items = parse_smiles_file(args.input)
    if not items:
        raise SystemExit("No valid SMILES found")

    settings = _build_effective_settings(args)
    prep_cfg, tautomer_cfg, stereo_cfg, protomer_cfg = _build_configs(settings)

    rows: list[dict[str, object]] = []
    try:
        for compound_id, smiles in items:
            row = _prepare_record(
                compound_id,
                smiles,
                preset=args.preset,
                ph=float(settings["ph"]),
                ionization_mode=str(settings["ionization"]),
                prep_cfg=prep_cfg,
                tautomer_cfg=tautomer_cfg,
                stereo_cfg=stereo_cfg,
                protomer_cfg=protomer_cfg,
            )
            if row is not None:
                rows.append(row)
    except ImportError as exc:
        raise SystemExit(str(exc)) from exc

    if not rows:
        raise SystemExit("No valid molecules could be prepared")

    df = pd.DataFrame(rows)
    write_table(df, output_path)

    summary_payload = _build_summary(
        rows,
        preset=args.preset,
        input_path=args.input,
        output_path=output_path,
        summary_path=summary_path,
        settings=settings,
    )
    Path(summary_path).write_text(
        json.dumps(summary_payload, indent=2) + "\n",
        encoding="utf-8",
    )

    try:
        write_run_metadata(
            tool="molscope-prep",
            output_table_path=output_path,
            input_path=args.input,
            smiles_column="Calc_Canonical_SMILES",
            parameters={
                "preset": args.preset,
                "ph": float(settings["ph"]),
                "ionization": str(settings["ionization"]),
                "tautomer_mode": str(settings["tautomer_mode"]),
                "stereo_mode": str(settings["stereo_mode"]),
                "use_prep": bool(settings["use_prep"]),
                "calc_on_protomer": bool(settings["calc_on_protomer"]),
            },
            artifacts={"prep_summary_json": str(Path(summary_path).resolve())},
        )
    except Exception as exc:
        print(
            f"[molscope prep] warning: failed to write metadata sidecar: {exc}",
            file=sys.stderr,
        )

    print(f"Saved: {output_path} (compounds={len(rows)}, columns={len(df.columns)})")
    print(f"Summary: {summary_path}")


if __name__ == "__main__":
    main()
