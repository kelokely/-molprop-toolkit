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
import os
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

# Allow running directly without installing the package.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd

from molprop_toolkit.core import detect_id_column, read_table, write_table
from molprop_toolkit.core.metadata import write_run_metadata

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
        df = read_table(path)
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


def _copy_site_assets(outdir: Path) -> None:
    """Copy GitHub Pages site assets into a retrosynthesis output folder.

    This makes retrosynthesis HTML outputs match the look-and-feel of the Pages
    site (docs/index.html). We copy:
    - docs/assets/style.css
    - docs/assets/logo.svg

    If the repo docs assets are not available (e.g., pip install without docs),
    we fall back to a minimal stylesheet so the report remains readable.
    """

    assets_out = outdir / "assets"
    assets_out.mkdir(parents=True, exist_ok=True)

    repo_root = Path(__file__).resolve().parents[1]
    assets_src = repo_root / "docs" / "assets"

    for fn in ("style.css", "logo.svg"):
        src = assets_src / fn
        dst = assets_out / fn
        if src.exists() and not dst.exists():
            shutil.copy2(src, dst)

    # Minimal fallback if docs assets aren't present
    if not (assets_out / "style.css").exists():
        (assets_out / "style.css").write_text(
            ":root{--bg:#0b1020;--bg2:#0d1630;--card:#111a33;--text:#e8eefc;--muted:#b8c3e1;--border:rgba(255,255,255,.10);--brand:#6ee7ff;}\n"
            "body{margin:0;font-family:system-ui,-apple-system,Segoe UI,Roboto,Helvetica,Arial;background:linear-gradient(180deg,var(--bg),var(--bg2));color:var(--text);line-height:1.55}\n"
            ".container{max-width:1100px;margin:0 auto;padding:0 22px}\n"
            ".site-header{position:sticky;top:0;z-index:20;backdrop-filter:blur(14px);background:rgba(11,16,32,.55);border-bottom:1px solid var(--border)}\n"
            ".header-inner{display:flex;align-items:center;justify-content:space-between;gap:14px;padding:14px 0}\n"
            ".brand{display:flex;align-items:center;gap:12px;text-decoration:none;color:var(--text)}\n"
            ".brand-logo{width:34px;height:34px}\n"
            ".brand-text{display:flex;flex-direction:column}\n"
            ".brand-name{font-weight:700}\n"
            ".brand-tag{font-size:12px;color:var(--muted)}\n"
            ".nav{display:flex;gap:16px;flex-wrap:wrap}\n"
            ".nav a{color:var(--muted);text-decoration:none;font-size:14px}\n"
            ".nav a:hover{color:var(--text)}\n"
            ".hero{padding:34px 0 10px}\n"
            ".section{padding:22px 0}\n"
            ".card{background:rgba(17,26,51,.65);border:1px solid var(--border);border-radius:16px;padding:18px}\n"
            ".muted{color:var(--muted)}\n"
            "code,pre{font-family:ui-monospace,SFMono-Regular,Menlo,Monaco,Consolas,Liberation Mono,Courier New,monospace;font-size:12px}\n",
            encoding="utf-8",
        )

    if not (assets_out / "logo.svg").exists():
        (assets_out / "logo.svg").write_text(
            "<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 64 64'><rect width='64' height='64' rx='14' fill='#0b1020'/><path d='M18 42 L28 22 L40 26 L46 42 L34 50 Z' fill='none' stroke='#6ee7ff' stroke-width='4' stroke-linejoin='round'/></svg>",
            encoding="utf-8",
        )


def _coerce_numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def _normalize_series(series: pd.Series, *, higher_is_better: bool) -> pd.Series:
    numeric = _coerce_numeric(series)
    values = numeric.dropna()
    if values.empty:
        return pd.Series(0.0, index=series.index, dtype="float64")

    minimum = float(values.min())
    maximum = float(values.max())
    if maximum == minimum:
        normalized = pd.Series(0.0, index=series.index, dtype="float64")
        normalized.loc[numeric.notna()] = 1.0
    else:
        normalized = (numeric - minimum) / (maximum - minimum)
        normalized = normalized.fillna(0.0)

    if not higher_is_better:
        normalized = 1.0 - normalized

    return normalized.clip(lower=0.0, upper=1.0)


def _apply_retro_reranking(df_summary: pd.DataFrame, policy: str) -> pd.DataFrame:
    if df_summary.empty:
        return df_summary.copy()

    ranked = df_summary.copy()

    solved = ranked.get("is_solved", pd.Series(False, index=ranked.index)).fillna(False)
    solved = solved.astype("string").str.lower().isin({"true", "1", "yes"})

    top_score_norm = _normalize_series(ranked.get("top_score", pd.Series(dtype="float64")), higher_is_better=True)
    route_count_score = _normalize_series(
        ranked.get("number_of_solved_routes", ranked.get("number_of_routes", pd.Series(dtype="float64"))),
        higher_is_better=True,
    )
    step_score = _normalize_series(
        ranked.get("number_of_steps", pd.Series(dtype="float64")),
        higher_is_better=False,
    )
    not_in_stock_score = _normalize_series(
        ranked.get("precursors_not_in_stock", pd.Series(dtype="float64")),
        higher_is_better=False,
    )

    precursor_count = _coerce_numeric(
        ranked.get("number_of_precursors", pd.Series(dtype="float64"))
    )
    in_stock = _coerce_numeric(
        ranked.get("precursors_in_stock", pd.Series(dtype="float64"))
    ).fillna(0.0)
    stock_coverage = pd.Series(0.0, index=ranked.index, dtype="float64")
    valid_precursors = precursor_count > 0
    stock_coverage.loc[valid_precursors] = (
        in_stock.loc[valid_precursors] / precursor_count.loc[valid_precursors]
    )
    stock_coverage = stock_coverage.clip(lower=0.0, upper=1.0)

    policy_weights: dict[str, dict[str, float]] = {
        "balanced": {
            "solved": 0.30,
            "top_score": 0.25,
            "stock_coverage": 0.20,
            "steps": 0.15,
            "route_count": 0.10,
        },
        "score": {
            "solved": 0.20,
            "top_score": 0.55,
            "stock_coverage": 0.10,
            "steps": 0.05,
            "route_count": 0.10,
        },
        "stock-first": {
            "solved": 0.25,
            "top_score": 0.10,
            "stock_coverage": 0.40,
            "steps": 0.15,
            "route_count": 0.10,
        },
        "shortest": {
            "solved": 0.20,
            "top_score": 0.10,
            "stock_coverage": 0.20,
            "steps": 0.40,
            "route_count": 0.10,
        },
    }
    weights = policy_weights[policy]

    rerank_score = (
        solved.astype("float64") * weights["solved"]
        + top_score_norm * weights["top_score"]
        + stock_coverage * weights["stock_coverage"]
        + step_score * weights["steps"]
        + route_count_score * weights["route_count"]
    )

    ranked["Retro_Stock_Coverage"] = stock_coverage.round(4)
    ranked["Retro_NotInStock_Score"] = not_in_stock_score.round(4)
    ranked["Retro_Rerank_Policy"] = policy
    ranked["Retro_Rerank_Score"] = rerank_score.round(4)

    sort_frame = pd.DataFrame(
        {
            "Retro_Rerank_Score": ranked["Retro_Rerank_Score"],
            "top_score": _coerce_numeric(ranked.get("top_score", pd.Series(dtype="float64"))).fillna(-1.0),
            "number_of_steps": _coerce_numeric(
                ranked.get("number_of_steps", pd.Series(dtype="float64"))
            ).fillna(float("inf")),
            "Retro_Stock_Coverage": ranked["Retro_Stock_Coverage"].fillna(0.0),
            "precursors_not_in_stock": _coerce_numeric(
                ranked.get("precursors_not_in_stock", pd.Series(dtype="float64"))
            ).fillna(float("inf")),
        },
        index=ranked.index,
    )
    ranked = ranked.loc[
        sort_frame.sort_values(
            by=[
                "Retro_Rerank_Score",
                "top_score",
                "Retro_Stock_Coverage",
                "number_of_steps",
                "precursors_not_in_stock",
            ],
            ascending=[False, False, False, True, True],
        ).index
    ].reset_index(drop=True)
    ranked["Retro_Rerank_Rank"] = range(1, len(ranked) + 1)

    return ranked


def _write_summary_markdown(
    path: Path,
    *,
    summary_df: pd.DataFrame,
    policy: str,
    retro_table_name: Optional[str],
) -> None:
    solved_count = int(
        summary_df.get("is_solved", pd.Series(False, index=summary_df.index))
        .fillna(False)
        .astype("string")
        .str.lower()
        .isin({"true", "1", "yes"})
        .sum()
    )
    lines = [
        "# Retrosynthesis Summary",
        "",
        f"- rerank policy: `{policy}`",
        f"- targets processed: {len(summary_df)}",
        f"- solved targets: {solved_count}",
        "",
        "## Top reranked targets",
    ]

    if summary_df.empty:
        lines.append("- No retrosynthesis rows were generated.")
    else:
        for _, row in summary_df.head(5).iterrows():
            lines.append(
                "- "
                f"{row.get('Compound_ID', 'unknown')} "
                f"(rank {row.get('Retro_Rerank_Rank', '')}, "
                f"score {row.get('Retro_Rerank_Score', '')}, "
                f"steps {row.get('number_of_steps', '')}, "
                f"stock coverage {row.get('Retro_Stock_Coverage', '')})"
            )

    if retro_table_name:
        lines.extend(
            [
                "",
                "## Downstream table",
                f"- retro-enriched campaign table: `{retro_table_name}`",
            ]
        )

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _build_retro_enriched_table(
    input_path: str,
    summary_df: pd.DataFrame,
    *,
    requested_id_column: Optional[str],
) -> tuple[pd.DataFrame, str]:
    input_df = read_table(input_path)
    id_column = requested_id_column or detect_id_column(input_df)
    if id_column not in input_df.columns:
        raise SystemExit(
            f"molscope retro: could not find input ID column '{id_column}' for retro table merge"
        )

    retro_columns = [
        "Compound_ID",
        "target",
        "is_solved",
        "top_score",
        "number_of_routes",
        "number_of_solved_routes",
        "search_time",
        "max_transforms",
        "number_of_steps",
        "number_of_precursors",
        "precursors_in_stock",
        "precursors_not_in_stock",
        "Retro_Stock_Coverage",
        "Retro_Rerank_Policy",
        "Retro_Rerank_Score",
        "Retro_Rerank_Rank",
        "html_page",
    ]
    retro_summary = summary_df[retro_columns].rename(
        columns={
            "Compound_ID": "__retro_join_id",
            "target": "Retro_Target",
            "is_solved": "Retro_IsSolved",
            "top_score": "Retro_TopScore",
            "number_of_routes": "Retro_RouteCount",
            "number_of_solved_routes": "Retro_SolvedRouteCount",
            "search_time": "Retro_SearchTime",
            "max_transforms": "Retro_MaxTransforms",
            "number_of_steps": "Retro_StepCount",
            "number_of_precursors": "Retro_PrecursorCount",
            "precursors_in_stock": "Retro_InStockCount",
            "precursors_not_in_stock": "Retro_NotInStockCount",
            "html_page": "Retro_HTMLPage",
        }
    )

    merged = input_df.copy()
    merged["__retro_join_id"] = merged[id_column].astype("string").str.strip()
    retro_summary["__retro_join_id"] = retro_summary["__retro_join_id"].astype("string").str.strip()
    merged = merged.merge(retro_summary, on="__retro_join_id", how="left")
    merged.drop(columns="__retro_join_id", inplace=True)
    merged["Retro_RecordMatched"] = merged["Retro_Rerank_Rank"].notna()
    return merged, id_column


def _write_html_index(outdir: Path, rows: List[dict], site_url: str | None) -> None:
    items = []
    for r in rows:
        href = r.get("html_page", "")
        items.append(
            "<tr>"
            f"<td><a href='{href}'>{r.get('Compound_ID','')}</a></td>"
            f"<td><code>{r.get('target','')}</code></td>"
            f"<td>{r.get('is_solved','')}</td>"
            f"<td>{r.get('top_score','')}</td>"
            f"<td>{r.get('number_of_steps','')}</td>"
            f"<td>{r.get('Retro_Stock_Coverage','')}</td>"
            f"<td>{r.get('Retro_Rerank_Score','')}</td>"
            f"<td>{r.get('Retro_Rerank_Rank','')}</td>"
            "</tr>"
        )

    site_links_html = ""
    if site_url:
        site = site_url.rstrip("/") + "/"
        site_links_html = (
            f"<a href='{site}'>Site home</a>"
            f"<a href='{site}quickstart.html'>Quickstart</a>"
            f"<a href='{site}tools/'>Tools</a>"
        )

    html = f"""<!doctype html>
<html lang='en'>
<head>
  <meta charset='utf-8'>
  <meta name='viewport' content='width=device-width, initial-scale=1'>
  <title>MolScope Retro • Index</title>
  <link rel='stylesheet' href='assets/style.css'>
  <link rel='icon' href='assets/logo.svg' type='image/svg+xml'>
  <style>
    .table{{width:100%; border-collapse:collapse;}}
    .table th,.table td{{border:1px solid rgba(255,255,255,.12); padding:10px; font-size:14px; vertical-align:top;}}
    .table th{{background: rgba(0,0,0,.18); color: var(--text); text-align:left;}}
    .table td{{color: var(--muted);}}
    .table a{{color: var(--brand); text-decoration:none;}}
    .table a:hover{{text-decoration:underline;}}
    .table code{{color: var(--text);}}
  </style>
</head>
<body>
  <header class='site-header'>
    <div class='container header-inner'>
      <a class='brand' href='./'>
          <img class='brand-logo' src='assets/logo.svg' alt='MolScope Toolkit logo'>
          <span class='brand-text'>
            <span class='brand-name'>MolScope Toolkit</span>
            <span class='brand-tag'>Retrosynthesis results (AiZynthFinder)</span>
          </span>
        </a>
      <nav class='nav'>
        <a href='routes_summary.csv'>routes_summary.csv</a>
        <a href='routes_summary.json'>routes_summary.json</a>
        <a href='summary.md'>summary.md</a>
        {site_links_html}
      </nav>
    </div>
  </header>

  <section class='hero'>
    <div class='container'>
      <h1 style='margin:0 0 10px; font-size:34px; line-height:1.15'>Retrosynthesis route index</h1>
      <p class='muted' style='max-width: 90ch'>
        Generated by <code>molprop-retro</code>. Each row links to a per-target page with extracted routes.
        Scores, route counts, and rerank fields come from AiZynthFinder’s batch output plus MolScope post-processing.
      </p>
    </div>
  </section>

  <section class='section'>
    <div class='container'>
      <div class='card'>
        <table class='table'>
          <thead>
            <tr>
              <th>Compound_ID</th>
              <th>Target</th>
              <th>Solved</th>
              <th>Top score</th>
              <th>Steps</th>
              <th>Stock coverage</th>
              <th>Rerank score</th>
              <th>Rerank rank</th>
            </tr>
          </thead>
          <tbody>
            {''.join(items)}
          </tbody>
        </table>
      </div>
    </div>
  </section>
</body>
</html>"""

    (outdir / "index.html").write_text(html, encoding="utf-8")


def _write_target_page(outdir: Path, cid: str, target: str, trees: list, images: List[str], site_url: str | None) -> str:
    site_links_html = ""
    if site_url:
        site = site_url.rstrip("/") + "/"
        site_links_html = f"<a href='{site}'>Site home</a><a href='{site}tools/'>Tools</a>"

    blocks = []
    for i, t in enumerate(trees):
        img = images[i] if i < len(images) else ""
        img_src = f"../{img}" if img else ""

        blocks.append(
            "<div class='card' style='margin: 14px 0'>"
            f"<div style='font-weight:800; margin-bottom:8px'>Route {i+1}</div>"
            + (
                f"<div style='margin:10px 0'><img src='{img_src}' style='max-width:100%; height:auto; border-radius:12px; border:1px solid rgba(255,255,255,.10)'></div>"
                if img_src
                else ""
            )
            + f"<details><summary class='muted' style='cursor:pointer'>Route JSON</summary><pre style='white-space:pre-wrap; overflow:auto; margin-top:10px; padding:14px; border-radius:12px; background: rgba(0,0,0,.28); border: 1px solid rgba(255,255,255,.08);'>{json.dumps(t, indent=2)[:200000]}</pre></details>"
            + "</div>"
        )

    page = f"""<!doctype html>
<html lang='en'>
<head>
  <meta charset='utf-8'>
  <meta name='viewport' content='width=device-width, initial-scale=1'>
  <title>MolScope Retro • {cid}</title>
  <link rel='stylesheet' href='../assets/style.css'>
  <link rel='icon' href='../assets/logo.svg' type='image/svg+xml'>
</head>
<body>
  <header class='site-header'>
    <div class='container header-inner'>
      <a class='brand' href='../index.html'>
        <img class='brand-logo' src='../assets/logo.svg' alt='MolScope Toolkit logo'>
        <span class='brand-text'>
          <span class='brand-name'>MolScope Toolkit</span>
          <span class='brand-tag'>Retrosynthesis: {cid}</span>
        </span>
      </a>
      <nav class='nav'>
        <a href='../index.html'>&larr; Index</a>
        <a href='../routes_summary.csv'>Summary CSV</a>
        {site_links_html}
      </nav>
    </div>
  </header>

  <section class='hero'>
    <div class='container'>
      <h1 style='margin:0 0 10px; font-size:32px; line-height:1.15'>{cid}</h1>
      <p class='muted' style='max-width: 90ch'>Target: <code>{target}</code></p>
    </div>
  </section>

  <section class='section'>
    <div class='container'>
      {''.join(blocks) if blocks else "<div class='card'><p class='muted' style='margin:0'>No extracted routes.</p></div>"}
    </div>
  </section>
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


def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(description="MolScope Toolkit retrosynthesis wrapper (AiZynthFinder)")
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
    ap.add_argument(
        "--rerank-policy",
        default="balanced",
        choices=["balanced", "score", "stock-first", "shortest"],
        help="Policy used to rerank route summaries for downstream prioritization.",
    )
    ap.add_argument(
        "--write-retro-table",
        action="store_true",
        help="Join retro summary fields back onto the input results table and write retro_enriched.csv.",
    )
    ap.add_argument(
        "--site-url",
        default=None,
        help=(
            "Optional base URL of the MolScope Toolkit GitHub Pages site to link back to "
            "(e.g., https://kelokely.github.io/-molprop-toolkit/). You can also set "
            "MOLSCOPE_SITE_URL or MOLPROP_SITE_URL."
        ),
    )
    args = ap.parse_args(argv)

    targets = _read_targets(args.input, args.smiles_col, args.id_col, args.limit)
    if not targets:
        raise SystemExit("No targets found")

    # Site URL used for nav links in the generated HTML
    site_url = (
        args.site_url
        or os.environ.get("MOLSCOPE_SITE_URL")
        or os.environ.get("MOLPROP_SITE_URL")
        or "https://kelokely.github.io/-molprop-toolkit/"
    )

    stem = Path(args.input).stem
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    outdir = Path(args.outdir or (Path("retrosynthesis") / f"{stem}_retro_{ts}"))
    outdir.mkdir(parents=True, exist_ok=True)

    # Copy the GitHub Pages site styling into this output folder
    _copy_site_assets(outdir)

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

        page_href = _write_target_page(outdir, cid, str(r.get("target", "")), trees, images, site_url)

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
    df_summary = _apply_retro_reranking(pd.DataFrame(summary_rows), args.rerank_policy)
    df_summary.to_csv(outdir / "routes_summary.csv", index=False)
    (outdir / "routes_summary.json").write_text(
        json.dumps(df_summary.to_dict(orient="records"), indent=2), encoding="utf-8"
    )

    retro_table_path: Optional[Path] = None
    retro_table_id_column: Optional[str] = None
    if args.write_retro_table:
        if Path(args.input).suffix.lower() not in (".csv", ".tsv", ".tab"):
            raise SystemExit(
                "molscope retro: --write-retro-table requires a CSV/TSV results table input."
            )
        retro_table_path = outdir / "retro_enriched.csv"
        retro_df, retro_table_id_column = _build_retro_enriched_table(
            args.input,
            df_summary,
            requested_id_column=args.id_col,
        )
        write_table(retro_df, str(retro_table_path))

    summary_md_path = outdir / "summary.md"
    _write_summary_markdown(
        summary_md_path,
        summary_df=df_summary,
        policy=args.rerank_policy,
        retro_table_name=retro_table_path.name if retro_table_path else None,
    )

    _write_html_index(outdir, df_summary.to_dict(orient="records"), site_url)

    routes_summary_path = outdir / "routes_summary.csv"
    write_run_metadata(
        tool="molscope-retro",
        input_path=args.input,
        output_table_path=routes_summary_path,
        smiles_column=args.smiles_col,
        parameters={
            "config": args.config,
            "id_column": args.id_col,
            "limit": args.limit,
            "top_routes": args.top_routes,
            "nproc": args.nproc,
            "rerank_policy": args.rerank_policy,
            "write_retro_table": bool(args.write_retro_table),
            "no_images": bool(args.no_images),
        },
        artifacts={
            "aizynth_output": output_file.name,
            "summary_json": "routes_summary.json",
            "summary_markdown": summary_md_path.name,
            "html_index": "index.html",
            "retro_table": retro_table_path.name if retro_table_path else None,
        },
        notes="Route summaries are post-processed into rerankable prioritization signals.",
    )
    if retro_table_path is not None:
        write_run_metadata(
            tool="molscope-retro",
            input_path=args.input,
            output_table_path=retro_table_path,
            smiles_column=args.smiles_col,
            parameters={
                "config": args.config,
                "id_column": retro_table_id_column,
                "rerank_policy": args.rerank_policy,
            },
            artifacts={
                "routes_summary": routes_summary_path.name,
                "summary_markdown": summary_md_path.name,
                "html_index": "index.html",
            },
            notes="Input table enriched with retro reranking fields for downstream reports and picklists.",
        )

    print(f"Retrosynthesis outputs written to: {outdir}")
    print(f"- AiZynthFinder raw output: {output_file.name}")
    print("- Summary CSV: routes_summary.csv")
    print("- Summary JSON: routes_summary.json")
    print("- Summary Markdown: summary.md")
    if retro_table_path is not None:
        print(f"- Retro-enriched table: {retro_table_path.name}")
    print("- HTML index: index.html")


if __name__ == "__main__":
    main()
