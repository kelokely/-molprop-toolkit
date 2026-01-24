#!/usr/bin/env python3
"""Build shareable report artifacts from a results CSV.

Outputs a directory containing:
- report.md
- report.html
- plots/*.png

Designed to be registry-driven via molprop_toolkit.core.registry.CATEGORY_SPECS.

Example:
  python tools/build_report.py results.csv --outdir reports/run1

If you want to limit categories:
  python tools/build_report.py results.csv --categories cns_mpo,toxicity
"""

from __future__ import annotations

import argparse
import sys
from datetime import datetime
from pathlib import Path

# Allow running directly without installing the package.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from molprop_toolkit.reporting import ReportBuilder  # noqa: E402
from molprop_toolkit.core import detect_id_column
from molprop_toolkit.core.columns import detect_best_smiles_column
from molprop_toolkit.fingerprints import FingerprintSpec, compute_similarity_columns, load_reference_smiles


def main() -> None:
    ap = argparse.ArgumentParser(description="MolProp Toolkit report builder (Markdown + HTML + plots)")
    ap.add_argument("input", help="Results table (CSV/TSV/Parquet)")
    ap.add_argument(
        "--outdir",
        help="Output directory (default: reports/<csvstem>_report_<timestamp>)",
        default=None,
    )
    ap.add_argument(
        "--categories",
        help="Comma-separated category keys to include (default: all)",
        default=None,
    )
    ap.add_argument("--top", type=int, default=10, help="Top rows to show per category")
    ap.add_argument("--no-plots", action="store_true", help="Disable plot generation")
    ap.add_argument("--no-md", action="store_true", help="Do not write report.md")
    ap.add_argument("--no-html", action="store_true", help="Do not write report.html")
    ap.add_argument("--title", default="MolProp Toolkit Report", help="Report title")

    # Optional similarity-to-reference columns (fingerprints)
    ap.add_argument(
        "--fp",
        default=None,
        help=(
            "Optional fingerprint kind for similarity columns: morgan, rdkit, maccs, atom_pair, torsion, usr, usrcat. "
            "If set with --sim-ref-* options, the report will compute Sim_<kind> columns and include a similarity section."
        ),
    )
    ap.add_argument("--fp-radius", type=int, default=2, help="Morgan radius (only for --fp morgan)")
    ap.add_argument("--fp-nbits", type=int, default=2048, help="Fingerprint bit length (where applicable)")
    ap.add_argument("--fp-chirality", action="store_true", help="Use chirality in Morgan fingerprints")
    ap.add_argument("--fp-features", action="store_true", help="Use feature-invariants for Morgan fingerprints")
    ap.add_argument("--fp-3d-num-confs", type=int, default=10, help="Conformers for 3D fingerprints (usr/usrcat)")
    ap.add_argument(
        "--fp-3d-minimize",
        choices=["none", "uff", "mmff"],
        default="none",
        help="Optional geometry minimization for 3D fingerprints",
    )
    ap.add_argument(
        "--sim-ref-id",
        action="append",
        default=[],
        help="Reference Compound_ID (or detected ID) in the input CSV. Can be provided multiple times.",
    )
    ap.add_argument(
        "--sim-ref-smiles",
        action="append",
        default=[],
        help="Reference SMILES string. Can be provided multiple times.",
    )
    ap.add_argument(
        "--sim-ref-file",
        default=None,
        help="Reference SMILES file (.smi). Lines can be 'name smiles' or 'smiles name'.",
    )
    ap.add_argument(
        "--sim-prefix",
        default="Sim",
        help="Column prefix for generated similarity columns (default: Sim)",
    )
    args = ap.parse_args()

    csv_path = args.input
    stem = Path(csv_path).stem
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    outdir = args.outdir or str(Path("reports") / f"{stem}_report_{ts}")

    categories = None
    if args.categories:
        categories = [c.strip() for c in args.categories.split(",") if c.strip()]

    # Optional: compute similarity columns to a reference set
    df = None
    if args.fp and (args.sim_ref_id or args.sim_ref_smiles or args.sim_ref_file):
        from molprop_toolkit.core import read_table

        df = read_table(csv_path)
        id_col = detect_id_column(df)
        smiles_col = detect_best_smiles_column(df)

        refs = load_reference_smiles(
            df=df,
            id_col=id_col,
            smiles_col=smiles_col,
            ref_ids=args.sim_ref_id,
            ref_smiles=args.sim_ref_smiles,
            ref_file=args.sim_ref_file,
        )

        spec = FingerprintSpec(
            kind=args.fp,
            radius=int(args.fp_radius),
            nbits=int(args.fp_nbits),
            use_chirality=bool(args.fp_chirality),
            use_features=bool(args.fp_features),
            num_confs=int(args.fp_3d_num_confs),
            minimize=args.fp_3d_minimize,
        )

        df = compute_similarity_columns(df, refs=refs, spec=spec, smiles_col=smiles_col, id_col=id_col, column_prefix=args.sim_prefix)

    builder = ReportBuilder(csv_path=csv_path, df=df, categories=categories, top_n=args.top)
    artifacts = builder.build(
        outdir=outdir,
        write_markdown=not args.no_md,
        write_html=not args.no_html,
        make_plots=not args.no_plots,
        title=args.title,
    )

    print(f"Report written to: {artifacts.outdir}")
    if not args.no_md:
        print(f"- Markdown: {artifacts.markdown_path}")
    if not args.no_html:
        print(f"- HTML: {artifacts.html_path}")
    if not args.no_plots:
        print(f"- Plots dir: {artifacts.plots_dir}")


if __name__ == "__main__":
    main()

