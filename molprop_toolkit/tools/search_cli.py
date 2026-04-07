#!/usr/bin/env python3
"""Search MolScope-compatible tables by exact SMILES, SMARTS, or motif panel."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any, Optional, Sequence

import pandas as pd

from molprop_toolkit.core import (
    detect_best_smiles_column,
    detect_id_column,
    read_table,
)
from molprop_toolkit.core.metadata import sha256_file, write_run_metadata

SEARCH_IMPORT_ERROR: ModuleNotFoundError | None = None

try:
    from rdkit import Chem
except ModuleNotFoundError as exc:
    if exc.name and exc.name.startswith("rdkit"):
        SEARCH_IMPORT_ERROR = exc
    else:  # pragma: no cover
        raise


def _require_search_dependencies() -> None:
    if SEARCH_IMPORT_ERROR is None:
        return
    raise SystemExit(
        "molscope search requires RDKit. Install RDKit, preferably via conda-forge, "
        "before running this command."
    ) from SEARCH_IMPORT_ERROR


def _default_output_dir(input_path: str) -> str:
    p = Path(input_path)
    return str(p.with_name(f"{p.stem}_search"))


def _canonicalize_smiles(smiles: str) -> Optional[str]:
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return None
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


def _iter_molecule_rows(
    df: pd.DataFrame,
    *,
    smiles_column: str,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for idx, row in df.iterrows():
        smiles = row.get(smiles_column)
        mol = None
        canonical_smiles = None
        if pd.notna(smiles):
            mol = Chem.MolFromSmiles(str(smiles))
            if mol is not None:
                canonical_smiles = Chem.MolToSmiles(
                    mol, canonical=True, isomericSmiles=True
                )
        rows.append(
            {
                "index": idx,
                "smiles": smiles,
                "mol": mol,
                "canonical_smiles": canonical_smiles,
            }
        )
    return rows


def _parse_motif_file(path: str) -> list[dict[str, str]]:
    motif_path = Path(path)
    suffix = motif_path.suffix.lower()

    if suffix in {".csv", ".tsv"}:
        sep = "\t" if suffix == ".tsv" else ","
        df = pd.read_csv(motif_path, sep=sep)
        if {"name", "smarts"}.issubset(df.columns):
            return [
                {"name": str(row["name"]), "smarts": str(row["smarts"])}
                for _, row in df.iterrows()
                if pd.notna(row["name"]) and pd.notna(row["smarts"])
            ]
        if df.shape[1] >= 2:
            cols = list(df.columns[:2])
            return [
                {"name": str(row[cols[0]]), "smarts": str(row[cols[1]])}
                for _, row in df.iterrows()
                if pd.notna(row[cols[0]]) and pd.notna(row[cols[1]])
            ]
        raise SystemExit(
            "molscope search: motif file must contain name/smarts columns or at least two columns"
        )

    motifs: list[dict[str, str]] = []
    with motif_path.open("r", encoding="utf-8") as handle:
        for _i, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            if "\t" in line:
                name, smarts = line.split("\t", 1)
            elif "," in line:
                name, smarts = next(csv.reader([line]))
            else:
                name = f"motif_{len(motifs) + 1:03d}"
                smarts = line
            motifs.append({"name": name.strip(), "smarts": smarts.strip()})
    return motifs


def _run_exact_search(
    df: pd.DataFrame,
    *,
    query_smiles: str,
    query_name: str,
    id_column: str,
    molecule_rows: list[dict[str, Any]],
) -> tuple[pd.DataFrame, dict[str, Any]]:
    canonical_query = _canonicalize_smiles(query_smiles)
    if canonical_query is None:
        raise SystemExit("molscope search: invalid query SMILES for --exact-smiles")

    matched_indices = [
        row["index"]
        for row in molecule_rows
        if row["canonical_smiles"] == canonical_query
    ]
    results = df.loc[matched_indices].copy()
    results["Search_Query_Name"] = query_name
    results["Search_Query"] = query_smiles
    results["Search_Query_Type"] = "exact"
    results["Search_Match_Count"] = 1

    summary = {
        "query_type": "exact",
        "query_name": query_name,
        "query": query_smiles,
        "canonical_query": canonical_query,
        "matched_rows": int(len(results)),
        "matched_ids": results[id_column].astype(str).tolist() if not results.empty else [],
    }
    return results, summary


def _run_smarts_search(
    df: pd.DataFrame,
    *,
    smarts: str,
    query_name: str,
    molecule_rows: list[dict[str, Any]],
) -> tuple[pd.DataFrame, dict[str, Any]]:
    query_mol = Chem.MolFromSmarts(smarts)
    if query_mol is None:
        raise SystemExit("molscope search: invalid SMARTS query for --smarts")

    matched_indices: list[int] = []
    match_counts: list[int] = []
    for row in molecule_rows:
        mol = row["mol"]
        if mol is None:
            continue
        matches = mol.GetSubstructMatches(query_mol)
        if matches:
            matched_indices.append(int(row["index"]))
            match_counts.append(int(len(matches)))

    results = df.loc[matched_indices].copy()
    results["Search_Query_Name"] = query_name
    results["Search_Query"] = smarts
    results["Search_Query_Type"] = "smarts"
    results["Search_Match_Count"] = match_counts

    summary = {
        "query_type": "smarts",
        "query_name": query_name,
        "query": smarts,
        "matched_rows": int(len(results)),
        "total_substructure_hits": int(sum(match_counts)),
    }
    return results, summary


def _run_motif_search(
    df: pd.DataFrame,
    *,
    motifs: list[dict[str, str]],
    id_column: str,
    smiles_column: str,
    molecule_rows: list[dict[str, Any]],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    compiled_motifs: list[dict[str, Any]] = []
    for motif in motifs:
        pattern = Chem.MolFromSmarts(motif["smarts"])
        if pattern is None:
            raise SystemExit(
                f"molscope search: invalid SMARTS in motif '{motif['name']}': {motif['smarts']}"
            )
        compiled_motifs.append(
            {
                "name": motif["name"],
                "smarts": motif["smarts"],
                "pattern": pattern,
            }
        )

    hit_rows: list[dict[str, Any]] = []
    motif_counts: dict[str, int] = {motif["name"]: 0 for motif in compiled_motifs}
    matched_motif_names: dict[int, list[str]] = {}

    for mol_row in molecule_rows:
        mol = mol_row["mol"]
        if mol is None:
            continue
        idx = int(mol_row["index"])
        row = df.loc[idx]
        for motif in compiled_motifs:
            matches = mol.GetSubstructMatches(motif["pattern"])
            if not matches:
                continue
            motif_counts[motif["name"]] += 1
            matched_motif_names.setdefault(idx, []).append(motif["name"])
            hit_rows.append(
                {
                    "motif_name": motif["name"],
                    "motif_smarts": motif["smarts"],
                    id_column: row[id_column],
                    smiles_column: row[smiles_column],
                    "match_count": int(len(matches)),
                }
            )

    motif_hits = pd.DataFrame(hit_rows)
    matched_indices = sorted(matched_motif_names.keys())
    results = df.loc[matched_indices].copy() if matched_indices else df.iloc[0:0].copy()
    if not results.empty:
        results["Search_Query_Name"] = "motif_panel"
        results["Search_Query"] = "motif_panel"
        results["Search_Query_Type"] = "motif_panel"
        results["Matched_Motif_Count"] = [
            len(matched_motif_names[int(idx)]) for idx in results.index
        ]
        results["Matched_Motifs"] = [
            ", ".join(sorted(matched_motif_names[int(idx)])) for idx in results.index
        ]
    else:
        results["Search_Query_Name"] = pd.Series(dtype="object")
        results["Search_Query"] = pd.Series(dtype="object")
        results["Search_Query_Type"] = pd.Series(dtype="object")
        results["Matched_Motif_Count"] = pd.Series(dtype="int64")
        results["Matched_Motifs"] = pd.Series(dtype="object")

    motif_summary = pd.DataFrame(
        [
            {
                "motif_name": motif["name"],
                "motif_smarts": motif["smarts"],
                "matched_compounds": int(motif_counts[motif["name"]]),
            }
            for motif in compiled_motifs
        ]
    )
    motif_summary = motif_summary.sort_values(
        ["matched_compounds", "motif_name"], ascending=[False, True]
    )

    summary = {
        "query_type": "motif_panel",
        "motif_count": int(len(compiled_motifs)),
        "matched_rows": int(len(results)),
        "motif_hits": int(len(motif_hits)),
    }
    return results, motif_hits, motif_summary, summary


def _write_summary_markdown(
    path: Path,
    *,
    input_name: str,
    summary: dict[str, Any],
    motif_summary: Optional[pd.DataFrame] = None,
) -> None:
    lines = [
        f"# Search Summary: {input_name}",
        "",
        f"- query type: {summary['query_type']}",
        f"- matched rows: {summary['matched_rows']}",
    ]
    if "query_name" in summary:
        lines.append(f"- query name: {summary['query_name']}")
    if "query" in summary:
        lines.append(f"- query: {summary['query']}")
    if "motif_count" in summary:
        lines.append(f"- motifs: {summary['motif_count']}")
        lines.append(f"- motif hits: {summary['motif_hits']}")
    lines.append("")

    if motif_summary is not None and not motif_summary.empty:
        lines.append("## Top motifs")
        for _, row in motif_summary.head(5).iterrows():
            lines.append(
                f"- {row['motif_name']}: {int(row['matched_compounds'])} compounds"
            )
        lines.append("")

    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def main(argv: Optional[Sequence[str]] = None) -> None:
    _require_search_dependencies()

    ap = argparse.ArgumentParser(
        prog="molscope search",
        description="Search MolScope-compatible tables by exact SMILES, SMARTS, or motif panel.",
    )
    ap.add_argument("input", help="Results table or compatible molecular table")
    ap.add_argument(
        "-o",
        "--output-dir",
        default=None,
        help="Output directory (default: <input>_search)",
    )
    ap.add_argument(
        "--id-column",
        default=None,
        help="Compound ID column (default: auto-detect)",
    )
    ap.add_argument(
        "--smiles-column",
        default=None,
        help="SMILES column to search (default: auto-detect)",
    )
    query_group = ap.add_mutually_exclusive_group(required=True)
    query_group.add_argument(
        "--exact-smiles",
        default=None,
        help="Exact SMILES query (canonicalized before matching)",
    )
    query_group.add_argument(
        "--smarts",
        default=None,
        help="SMARTS substructure query",
    )
    query_group.add_argument(
        "--motif-file",
        default=None,
        help="Motif panel file (TSV/CSV or text lines of name<TAB>SMARTS)",
    )
    ap.add_argument(
        "--query-name",
        default=None,
        help="Optional human-readable query name for single-query modes",
    )
    args = ap.parse_args(list(argv) if argv is not None else None)

    df = read_table(args.input)
    id_column = args.id_column or detect_id_column(df)
    smiles_column = args.smiles_column or detect_best_smiles_column(df.columns)
    if not id_column:
        raise SystemExit("molscope search: could not detect a compound ID column")
    if not smiles_column or smiles_column not in df.columns:
        raise SystemExit("molscope search: could not detect a valid SMILES column")

    output_dir = Path(args.output_dir or _default_output_dir(args.input))
    output_dir.mkdir(parents=True, exist_ok=True)

    molecule_rows = _iter_molecule_rows(df, smiles_column=smiles_column)

    search_results_path = output_dir / "search_results.csv"
    motif_hits_path = output_dir / "motif_hits.csv"
    motif_summary_path = output_dir / "motif_summary.csv"
    summary_json_path = output_dir / "summary.json"
    summary_md_path = output_dir / "summary.md"

    artifacts: dict[str, Any] = {
        "summary_json": str(summary_json_path.resolve()),
        "summary_markdown": str(summary_md_path.resolve()),
        "input_sha256": sha256_file(Path(args.input)),
    }

    motif_summary: Optional[pd.DataFrame] = None
    if args.exact_smiles:
        query_name = args.query_name or "exact_smiles"
        results, summary = _run_exact_search(
            df,
            query_smiles=args.exact_smiles,
            query_name=query_name,
            id_column=id_column,
            molecule_rows=molecule_rows,
        )
        results.to_csv(search_results_path, index=False)
        artifacts["search_results"] = str(search_results_path.resolve())
    elif args.smarts:
        query_name = args.query_name or "smarts_query"
        results, summary = _run_smarts_search(
            df,
            smarts=args.smarts,
            query_name=query_name,
            molecule_rows=molecule_rows,
        )
        results.to_csv(search_results_path, index=False)
        artifacts["search_results"] = str(search_results_path.resolve())
    else:
        motifs = _parse_motif_file(str(args.motif_file))
        results, motif_hits, motif_summary, summary = _run_motif_search(
            df,
            motifs=motifs,
            id_column=id_column,
            smiles_column=smiles_column,
            molecule_rows=molecule_rows,
        )
        results.to_csv(search_results_path, index=False)
        motif_hits.to_csv(motif_hits_path, index=False)
        motif_summary.to_csv(motif_summary_path, index=False)
        artifacts["search_results"] = str(search_results_path.resolve())
        artifacts["motif_hits"] = str(motif_hits_path.resolve())
        artifacts["motif_summary"] = str(motif_summary_path.resolve())

    summary.update(
        {
            "input_name": Path(args.input).stem,
            "rows": int(len(df)),
            "id_column": id_column,
            "smiles_column": smiles_column,
        }
    )
    summary_json_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    _write_summary_markdown(
        summary_md_path,
        input_name=Path(args.input).stem,
        summary=summary,
        motif_summary=motif_summary,
    )

    write_run_metadata(
        tool="molscope-search",
        output_table_path=search_results_path,
        input_path=args.input,
        smiles_column=smiles_column,
        parameters={
            "id_column": id_column,
            "smiles_column": smiles_column,
            "query_type": summary["query_type"],
            "query_name": summary.get("query_name"),
            "query": summary.get("query"),
            "motif_count": summary.get("motif_count"),
        },
        artifacts=artifacts,
        notes="First search slice: exact SMILES, SMARTS, and motif-panel workflows against MolScope-compatible tables.",
    )

    print(f"[ok] wrote search bundle to {output_dir}")


if __name__ == "__main__":
    main()
