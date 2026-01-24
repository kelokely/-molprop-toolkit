#!/usr/bin/env python3
"""Offline retrosynthesis wrapper using AiZynthFinder (optional feature).

This tool is intentionally thin: it delegates the actual planning to AiZynthFinder
(`aizynthcli`) and then post-processes the output into:
- a compact summary CSV/JSON for downstream tooling
- a lightweight HTML index for browsing per-target routes
- optional PNG images for the top-ranked routes

Why we wrap the CLI instead of re-implementing the Python API:
AiZynthFinder's CLI output format is documented and stable, and it supports
multiprocessing via `--nproc`.

Prerequisites (end-user):
- Install AiZynthFinder (and its dependencies) in the active environment
- Download model + stock data and generate a config.yml

References:
- AiZynthFinder CLI docs: https://molecularai.github.io/aizynthfinder/cli.html
- AiZynthFinder config docs: https://molecularai.github.io/aizynthfinder/configuration.html

Example:
  download_public_data aizynth_data
  molprop-retro results.csv --config aizynth_data/config.yml --top-routes 5 --nproc 4
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import List, Optional, Tuple

# Allow running directly without installing the package.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd

from molprop_toolkit.core import detect_id_column, read_csv


PREFERRED_SMILES_COLUMNS = [
    "Calc_Canonical_SMILES",
    "Canonical_SMILES",
    "SMILES",
    "Input_Canonical_SMILES",
]


def _detect_smiles_column(df: pd.DataFrame, smiles_col: Optional[str]) -> str:
    if smiles_col:
        if smiles_col not in df.columns:
            raise SystemExit(f"SMILES column not found: {smiles_col}. Available: {list(df.columns)[:30]}...")
        return smiles_col

    for c in PREFERRED_SMILES_COLUMNS:
        if c in df.columns:
            return c

    # Fallback: first column containing 'smiles'
    for c in df.columns:
        if "smiles" in c.lower():
            return c

    raise SystemExit(
        "Could not detect a SMILES column. Provide --smiles-col (e.g., SMILES or Calc_Canonical_SMILES)."
    )


def _read_targets(path: str, smiles_col: Optional[str], id_col: Optional[str], limit: Optional[int]) -> List[Tuple[str, str]]:
    p = Path(path)

    # If it's a CSV, use our normal column detection.
    if p.suffix.lower() in (".csv", ".tsv"):
        df = read_csv(path)
        smi_col = _detect_smiles_column(df, smiles_col)

        cid_col = id_col
        if cid_col:
            if cid_col not in df.columns:
                raise SystemExit(f"ID column not found: {cid_col}")
        else:
            cid_col = detect_id_column(df)

        out: List[Tuple[str, str]] = []
        for _, r in df.iterrows():
            cid = str(r.get(cid_col, ""))
            smi = str(r.get(smi_col, "")).strip()
            if not smi or smi.lower() in ("nan", "none"):
                continue
            out.append((cid, smi))
            if limit and len(out) >= limit:
                break
        return out

    # Otherwise treat as a SMILES list file (name<tab>smiles or smiles<tab>name).
    out2: List[Tuple[str, str]] = []
    with open(path, "r", encoding="utf-8") as f:
        for i, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) == 1:
                cid, smi = f"row{i}", parts[0]
            else:
                # Heuristic: SMILES has typical characters
                a, b = parts[0], parts[-1]
                smi = a if any(ch in a for ch in "[]=#()") else b
                cid = b if smi == a else a
            out2.append((cid, smi))
            if limit and len(out2) >= limit:
                break
    return out2


def _write_smiles_file(targets: List[Tuple[str, str]], outpath: Path) -> None:
    outpath.write_text("\n".join([smi for _, smi in targets]) + "\n", encoding="utf-8")


def _run_aizynthcli(config: str, smiles_file: Path, output_file: Path, nproc: int, extra_args: List[str]) -> None:
    exe = shutil.which("aizynthcli")
    if not exe:
        raise SystemExit(
            "AiZynthFinder is not available (missing 'aizynthcli').\n\n"
            "Install it in the active environment, for example:\n"
            "  pip install aizynthfinder[all]\n\n"
            "Then download model/stock data and create a config.yml (see AiZynthFinder docs)."
        )

    cmd = [
        exe,
        "--config",
        config,
        "--smiles",
        str(smiles_file),
        "--output",
        str(output_file),
    ]
    if nproc and nproc > 1:
        cmd += ["--nproc", str(nproc)]
    cmd += extra_args

    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)


def _safe_filename(s: str) -> str:
    keep = []
    for ch in s:
        if ch.isalnum() or ch in ("-", "_", "."):
            keep.append(ch)
        else:
            keep.append("_")
    out = "".join(keep).strip("_")
    return out[:80] or "compound"


def _write_html_index(outdir: Path, rows: List[dict]) -> None:
    items = []
    for r in rows:
        href = r.get("html_page", "")
        items.append(
            f"<tr><td><a href='{href}'>{r.get('Compound_ID','')}</a></td>"
            f"<td><code>{r.get('target','')}</code></td>"
            f"<td>{r.get('is_solved','')}</td>"
            f"<td>{r.get('top_score','')}</td>"
            f"<td>{r.get('number_of_routes','')}</td>"
            f"<td>{r.get('search_time','')}</td></tr>"
        )

    html = f"""<!doctype html>
<html lang='en'>
<head>
  <meta charset='utf-8'>
  <meta name='viewport' content='width=device-width, initial-scale=1'>
  <title>MolProp Retro • Index</title>
  <style>
    body{{font-family:system-ui, -apple-system, Segoe UI, Roboto, Helvetica, Arial; margin:28px;}}
    h1{{margin:0 0 8px;}}
    p{{color:#444; max-width: 80ch;}}
    table{{border-collapse:collapse; width:100%;}}
    th,td{{border:1px solid #ddd; padding:8px; font-size:14px; vertical-align:top;}}
    th{{background:#f7f7f7; text-align:left;}}
    code{{font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, Liberation Mono, Courier New, monospace; font-size:12px;}}
  </style>
</head>
<body>
  <h1>MolProp retrosynthesis results</h1>
  <p>This report was generated by <code>molprop-retro</code> using AiZynthFinder. Each row links to a per-target route page.</p>
  <table>
    <thead>
      <tr>
        <th>Compound_ID</th>
        <th>Target</th>
        <th>Solved</th>
        <th>Top score</th>
        <th>Routes</th>
        <th>Search time (s)</th>
      </tr>
    </thead>
    <tbody>
      {''.join(items)}
    </tbody>
  </table>
</body>
</html>"""

    (outdir / "index.html").write_text(html, encoding="utf-8")


def _write_target_page(outdir: Path, cid: str, target: str, trees: list, images: List[str]) -> str:
    rows = []
    for i, t in enumerate(trees):
        img = images[i] if i < len(images) else ""
        rows.append(
            f"<div style='border:1px solid #ddd; border-radius:10px; padding:12px; margin:12px 0'>"
            f"<div style='font-weight:700'>Route {i+1}</div>"
            + (f"<div style='margin-top:8px'><img src='{img}' style='max-width:100%; height:auto'></div>" if img else "")
            + f"<details style='margin-top:8px'><summary>Route JSON</summary><pre style='white-space:pre-wrap'>{json.dumps(t, indent=2)[:200000]}</pre></details>"
            + "</div>"
        )

    page = f"""<!doctype html>
<html lang='en'>
<head>
  <meta charset='utf-8'>
  <meta name='viewport' content='width=device-width, initial-scale=1'>
  <title>MolProp Retro • {cid}</title>
  <style>
    body{{font-family:system-ui, -apple-system, Segoe UI, Roboto, Helvetica, Arial; margin:28px;}}
    a{{color:#0b66d6; text-decoration:none;}}
    a:hover{{text-decoration:underline;}}
    code{{font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, Liberation Mono, Courier New, monospace; font-size:12px;}}
  </style>
</head>
<body>
  <div style='margin-bottom:14px'><a href='../index.html'>&larr; Back to index</a></div>
  <h1 style='margin:0 0 8px'>{cid}</h1>
  <div style='color:#444; margin-bottom:14px'>Target: <code>{target}</code></div>
  {''.join(rows) if rows else '<p>No extracted routes.</p>'}
</body>
</html>"""

    subdir = outdir / "targets"
    subdir.mkdir(parents=True, exist_ok=True)
    fn = f"{_safe_filename(cid)}.html"
    (subdir / fn).write_text(page, encoding="utf-8")
    return f"targets/{fn}"


def _render_route_images(outdir: Path, cid: str, trees: list, max_routes: int) -> List[str]:
    # Lazy import: only needed when generating images.
    try:
        from aizynthfinder.reactiontree import ReactionTree  # type: ignore
    except Exception:
        return []

    img_dir = outdir / "images" / _safe_filename(cid)
    img_dir.mkdir(parents=True, exist_ok=True)

    images: List[str] = []
    for itree, tree in enumerate(trees[:max_routes]):
        try:
            imagefile = img_dir / f"route{itree:03d}.png"
            ReactionTree.from_dict(tree).to_image().save(imagefile)  # type: ignore
            images.append(str(imagefile.relative_to(outdir)))
        except Exception:
            # Keep robust even if an individual image fails.
            continue

    return images


def main() -> None:
    ap = argparse.ArgumentParser(description="MolProp Toolkit retrosynthesis wrapper (AiZynthFinder)")
    ap.add_argument("input", help="Input results CSV OR a SMILES list file")
    ap.add_argument("--config", required=True, help="AiZynthFinder config.yml")
    ap.add_argument("--outdir", default=None, help="Output directory (default: retrosynthesis/<stem>_retro_<timestamp>)")
    ap.add_argument("--smiles-col", default=None, help="SMILES column for CSV inputs")
    ap.add_argument("--id-col", default=None, help="Compound ID column for CSV inputs")
    ap.add_argument("--limit", type=int, default=None, help="Limit number of targets")
    ap.add_argument("--top-routes", type=int, default=5, help="Number of top routes to extract/render")
    ap.add_argument("--nproc", type=int, default=1, help="Number of processes for aizynthcli")
    ap.add_argument(
        "--no-images",
        action="store_true",
        help="Disable route PNG rendering (still writes JSON/HTML summaries)",
    )
    ap.add_argument(
        "--aizynth-args",
        default=None,
        help="Extra args passed to aizynthcli (quoted string), e.g. --aizynth-args \"--cluster\"",
    )
    args = ap.parse_args()

    targets = _read_targets(args.input, args.smiles_col, args.id_col, args.limit)
    if not targets:
        raise SystemExit("No targets found")

    stem = Path(args.input).stem
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    outdir = Path(args.outdir or (Path("retrosynthesis") / f"{stem}_retro_{ts}"))
    outdir.mkdir(parents=True, exist_ok=True)

    smiles_file = outdir / "targets.smi"
    _write_smiles_file(targets, smiles_file)

    output_file = outdir / "aizynth_output.json.gz"
    extra_args = []
    if args.aizynth_args:
        # Very small parser: split on whitespace.
        extra_args = [x.strip() for x in args.aizynth_args.split() if x.strip()]

    _run_aizynthcli(args.config, smiles_file, output_file, args.nproc, extra_args)

    # Load output file (documented format: orient="table")
    df = pd.read_json(output_file, orient="table")

    # Build summary rows + HTML pages
    summary_rows: List[dict] = []
    for (cid, _smi), (_, r) in zip(targets, df.iterrows()):
        trees = r.get("trees") or []
        trees = list(trees)[: args.top_routes]

        images: List[str] = []
        if not args.no_images and trees:
            images = _render_route_images(outdir, cid, trees, args.top_routes)

        page_href = _write_target_page(outdir, cid, str(r.get("target", "")), trees, images)

        row = {
            "Compound_ID": cid,
            "target": r.get("target"),
            "is_solved": r.get("is_solved"),
            "top_score": r.get("top_score"),
            "number_of_routes": r.get("number_of_routes"),
            "number_of_solved_routes": r.get("number_of_solved_routes"),
            "search_time": r.get("search_time"),
            "max_transforms": r.get("max_transforms"),
            "number_of_steps": r.get("number_of_steps"),
            "number_of_precursors": r.get("number_of_precursors"),
            "precursors_in_stock": r.get("precursors_in_stock"),
            "precursors_not_in_stock": r.get("precursors_not_in_stock"),
            "html_page": page_href,
        }
        summary_rows.append(row)

    # Write outputs
    df_summary = pd.DataFrame(summary_rows)
    df_summary.to_csv(outdir / "routes_summary.csv", index=False)
    (outdir / "routes_summary.json").write_text(json.dumps(summary_rows, indent=2), encoding="utf-8")

    _write_html_index(outdir, summary_rows)

    print(f"Retrosynthesis outputs written to: {outdir}")
    print(f"- AiZynthFinder raw output: {output_file.name}")
    print("- Summary CSV: routes_summary.csv")
    print("- Summary JSON: routes_summary.json")
    print("- HTML index: index.html")


if __name__ == "__main__":
    main()

