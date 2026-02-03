"""HTML rendering for picklists.

The HTML is intentionally lightweight and self-contained. If a sketches directory is provided, we can copy the required
structure images into the picklists output directory and embed them in tables.
"""

from __future__ import annotations

import html as _html
import re
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import pandas as pd


def sanitize_filename(name: str, max_len: int = 120) -> str:
    name = (name or "").strip() or "compound"
    name = re.sub(r"[^A-Za-z0-9._-]+", "_", name)
    name = name.strip("._-") or "compound"
    return name[:max_len]


@dataclass(frozen=True)
class StructureAssetSpec:
    sketch_dir: Path
    views: Sequence[str] = ("calc",)
    formats: Sequence[str] = ("svg",)


def _candidate_structure_files(
    cid: str, views: Sequence[str], formats: Sequence[str]
) -> List[str]:
    base = sanitize_filename(cid)
    cands: List[str] = []

    # Single-view
    for fmt in formats:
        cands.append(f"structures/{base}.{fmt}")

    # Compare-view
    for v in views:
        for fmt in formats:
            cands.append(f"structures/{base}__{v}.{fmt}")

    # Keep stable order while de-duping
    seen = set()
    out: List[str] = []
    for x in cands:
        if x not in seen:
            seen.add(x)
            out.append(x)
    return out


def copy_structure_assets(
    *,
    cids: Iterable[str],
    asset: StructureAssetSpec,
    outdir: Path,
) -> Dict[str, List[str]]:
    """Copy structure files needed for a set of compound IDs.

    Returns mapping cid -> list of relative paths copied (relative to outdir).
    """

    out_struct = outdir / "structures"
    out_struct.mkdir(parents=True, exist_ok=True)

    mapping: Dict[str, List[str]] = {}

    for cid in cids:
        rels_copied: List[str] = []
        for rel in _candidate_structure_files(cid, asset.views, asset.formats):
            src = asset.sketch_dir / rel
            if not src.exists():
                continue
            dst = outdir / rel
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dst)
            rels_copied.append(rel)
        mapping[cid] = rels_copied

    return mapping


def _css() -> str:
    return """
    :root{--fg:#111;--muted:#555;--border:#e6e6e6;--bg:#fff;--mono:ui-monospace,SFMono-Regular,Menlo,Monaco,Consolas,monospace;}
    body{font-family:system-ui,-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;margin:24px;line-height:1.45;color:var(--fg);background:var(--bg);}
    h1{margin:0 0 6px 0;font-size:22px;}
    h2{margin:18px 0 10px 0;font-size:18px;}
    .sub{color:var(--muted);margin:0 0 18px 0;}
    code{font-family:var(--mono);font-size:0.95em;}
    table{border-collapse:collapse;width:100%;margin:10px 0 18px 0;}
    th,td{border:1px solid var(--border);padding:8px 10px;text-align:left;vertical-align:top;}
    th{background:#fafafa;position:sticky;top:0;}
    .imgcell img{max-width:260px;width:100%;height:auto;border:1px solid var(--border);border-radius:8px;background:#fff;}
    .warn{background:#fff7e6;border:1px solid #ffe0a3;padding:10px;border-radius:8px;margin:10px 0;color:#7a4a00;}
    .pill{display:inline-block;padding:2px 8px;border-radius:999px;background:#f2f2f2;color:#333;font-size:12px;margin-left:6px;}
    a{color:#0b5bd3;text-decoration:none;} a:hover{text-decoration:underline;}
    @media (max-width:600px){body{margin:14px;} th,td{font-size:13px;} .imgcell img{max-width:180px;}}
    """


def write_picklist_html(
    *,
    outpath: Path,
    title: str,
    subtitle: str,
    df: pd.DataFrame,
    id_col: str,
    warnings: Sequence[str] = (),
    structure_map: Optional[Dict[str, List[str]]] = None,
) -> None:
    structure_map = structure_map or {}

    cols = list(df.columns)

    # Add a structure column if we have any assets.
    add_struct = (
        any(structure_map.get(str(cid)) for cid in df[id_col].astype(str).tolist())
        if id_col in df.columns
        else False
    )

    headers = ["Structure"] + cols if add_struct else cols
    ths = "".join(f"<th>{_html.escape(h)}</th>" for h in headers)

    trs: List[str] = []
    for _, r in df.iterrows():
        tds: List[str] = []

        if add_struct:
            cid = str(r.get(id_col, ""))
            rels = structure_map.get(cid, [])
            imgs = []
            for rel in rels[:3]:
                imgs.append(
                    f"<div><img src='{_html.escape(rel)}' alt='{_html.escape(cid)}'/></div>"
                )
            tds.append(f"<td class='imgcell'>{''.join(imgs) if imgs else ''}</td>")

        for c in cols:
            v = r.get(c)
            if pd.isna(v):
                s = ""
            elif isinstance(v, float):
                s = f"{v:.3g}"
            else:
                s = str(v)
            tds.append(f"<td>{_html.escape(s)}</td>")
        trs.append(f"<tr>{''.join(tds)}</tr>")

    warn_block = "".join(f"<div class='warn'>{_html.escape(w)}</div>" for w in warnings)

    doc = f"""<!doctype html>
<html>
<head>
<meta charset='utf-8'/>
<meta name='viewport' content='width=device-width,initial-scale=1'/>
<title>{_html.escape(title)}</title>
<style>{_css()}</style>
</head>
<body>
<h1>{_html.escape(title)}</h1>
<p class='sub'>{_html.escape(subtitle)}</p>
{warn_block}
<table>
<thead><tr>{ths}</tr></thead>
<tbody>{''.join(trs)}</tbody>
</table>
</body>
</html>
"""

    outpath.write_text(doc, encoding="utf-8")


def write_index_html(
    *,
    outpath: Path,
    title: str,
    subtitle: str,
    entries: Sequence[Tuple[str, str, int, str]],
) -> None:
    # entries: (key, name, count, href)
    rows = []
    for key, name, count, href in entries:
        rows.append(
            "<tr>"
            f"<td><code>{_html.escape(key)}</code></td>"
            f"<td><a href='{_html.escape(href)}'>{_html.escape(name)}</a></td>"
            f"<td>{count}</td>"
            "</tr>"
        )

    doc = f"""<!doctype html>
<html>
<head>
<meta charset='utf-8'/>
<meta name='viewport' content='width=device-width,initial-scale=1'/>
<title>{_html.escape(title)}</title>
<style>{_css()}</style>
</head>
<body>
<h1>{_html.escape(title)}</h1>
<p class='sub'>{_html.escape(subtitle)}</p>
<h2>Picklists</h2>
<table>
<thead><tr><th>Key</th><th>Name</th><th>Count</th></tr></thead>
<tbody>
{''.join(rows)}
</tbody>
</table>
</body>
</html>
"""

    outpath.write_text(doc, encoding="utf-8")
