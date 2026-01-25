#!/usr/bin/env python
"""
molprop-similarity: Find similar compounds in a molecular library.

Command-line interface for fingerprint-based similarity searching.

Examples
--------
# Basic similarity search
molprop-similarity "CCO" library.smi -o results.csv

# Search with specific fingerprint and threshold
molprop-similarity "c1ccccc1" library.csv --fp morgan --threshold 0.7 --top 20

# Search from a query file (multiple queries)
molprop-similarity queries.smi library.csv --batch -o similarity_results/

# Pairwise similarity matrix
molprop-similarity --pairwise library.smi -o similarity_matrix.csv

# Diversity picking
molprop-similarity --diversity library.smi --pick 100 -o diverse_subset.csv

# Clustering by similarity
molprop-similarity --cluster library.smi --threshold 0.7 -o clusters.csv

# List available fingerprints and metrics
molprop-similarity --list-fps
molprop-similarity --list-metrics
"""

import argparse
import sys
import os
from pathlib import Path
from typing import Optional
import warnings

import pandas as pd
import numpy as np


def main():
    parser = argparse.ArgumentParser(
        description="Molecular similarity search and analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    
    # Positional arguments (optional based on mode)
    parser.add_argument(
        "query",
        nargs="?",
        help="Query SMILES string, or path to SMILES/CSV file for batch mode",
    )
    parser.add_argument(
        "library",
        nargs="?",
        help="Library to search (SMILES file, CSV, or Parquet)",
    )
    
    # Output options
    parser.add_argument(
        "-o", "--output",
        help="Output file or directory (default: stdout for single query)",
    )
    parser.add_argument(
        "--html",
        action="store_true",
        help="Generate HTML report in addition to CSV",
    )
    
    # Fingerprint options
    parser.add_argument(
        "--fp", "--fingerprint",
        dest="fp_type",
        default="morgan",
        help="Fingerprint type (default: morgan)",
    )
    parser.add_argument(
        "--radius",
        type=int,
        default=2,
        help="Morgan fingerprint radius (default: 2)",
    )
    parser.add_argument(
        "--nbits",
        type=int,
        default=2048,
        help="Fingerprint bit length (default: 2048)",
    )
    
    # Similarity options
    parser.add_argument(
        "--metric",
        default="tanimoto",
        help="Similarity metric (default: tanimoto)",
    )
    parser.add_argument(
        "-t", "--threshold",
        type=float,
        default=0.0,
        help="Minimum similarity threshold (default: 0.0)",
    )
    parser.add_argument(
        "-n", "--top",
        type=int,
        default=None,
        help="Return only top N results",
    )
    
    # Column specification
    parser.add_argument(
        "--smiles-col",
        default="SMILES",
        help="SMILES column name in library (default: SMILES)",
    )
    parser.add_argument(
        "--id-col",
        default="Compound_ID",
        help="ID column name in library (default: Compound_ID)",
    )
    
    # Mode flags
    parser.add_argument(
        "--batch",
        action="store_true",
        help="Batch mode: query is a file with multiple SMILES",
    )
    parser.add_argument(
        "--pairwise",
        action="store_true",
        help="Calculate pairwise similarity matrix",
    )
    parser.add_argument(
        "--diversity",
        action="store_true",
        help="Diversity picking mode",
    )
    parser.add_argument(
        "--cluster",
        action="store_true",
        help="Cluster molecules by similarity",
    )
    parser.add_argument(
        "--knn",
        type=int,
        default=None,
        help="Find K nearest neighbors",
    )
    
    # Diversity picking options
    parser.add_argument(
        "--pick",
        type=int,
        help="Number of compounds to pick (for --diversity)",
    )
    parser.add_argument(
        "--seed",
        help="Seed compound SMILES or index (for --diversity)",
    )
    
    # Performance options
    parser.add_argument(
        "-j", "--jobs",
        type=int,
        default=1,
        help="Number of parallel jobs (default: 1, use -1 for all CPUs)",
    )
    parser.add_argument(
        "--progress",
        action="store_true",
        help="Show progress bar",
    )
    
    # Info options
    parser.add_argument(
        "--list-fps",
        action="store_true",
        help="List available fingerprint types",
    )
    parser.add_argument(
        "--list-metrics",
        action="store_true",
        help="List available similarity metrics",
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Verbose output",
    )
    
    args = parser.parse_args()
    
    # Handle info commands
    if args.list_fps:
        _list_fingerprints()
        return 0
    
    if args.list_metrics:
        _list_metrics()
        return 0
    
    # Import here to delay RDKit load
    try:
        from molprop_toolkit.similarity import (
            similarity_search,
            pairwise_similarity,
            diversity_pick,
            cluster_by_similarity,
            find_nearest_neighbors,
            FINGERPRINT_TYPES,
            SIMILARITY_METRICS,
        )
    except ImportError as e:
        print(f"Error importing similarity module: {e}", file=sys.stderr)
        print("Make sure RDKit is installed: conda install -c conda-forge rdkit", 
              file=sys.stderr)
        return 1
    
    # Validate arguments based on mode
    if args.pairwise:
        if not args.query:
            parser.error("--pairwise requires a library file as the first argument")
        library = args.query
        return _run_pairwise(library, args)
    
    if args.diversity:
        if not args.query:
            parser.error("--diversity requires a library file as the first argument")
        if not args.pick:
            parser.error("--diversity requires --pick N to specify number of compounds")
        library = args.query
        return _run_diversity(library, args)
    
    if args.cluster:
        if not args.query:
            parser.error("--cluster requires a library file as the first argument")
        library = args.query
        return _run_cluster(library, args)
    
    # Standard similarity search
    if not args.query:
        parser.error("Query SMILES or file is required")
    if not args.library:
        parser.error("Library file is required for similarity search")
    
    if args.batch:
        return _run_batch_search(args)
    else:
        return _run_single_search(args)


def _list_fingerprints():
    """Print available fingerprint types."""
    from molprop_toolkit.similarity.fingerprints import FINGERPRINT_TYPES
    
    print("\nAvailable fingerprint types:\n")
    print(f"{'Type':<15} {'Description'}")
    print("-" * 60)
    for fp_type, info in FINGERPRINT_TYPES.items():
        print(f"{fp_type:<15} {info['description']}")
    print()


def _list_metrics():
    """Print available similarity metrics."""
    from molprop_toolkit.similarity.metrics import SIMILARITY_METRICS
    
    print("\nAvailable similarity metrics:\n")
    print(f"{'Metric':<15} {'Description'}")
    print("-" * 70)
    for metric, desc in SIMILARITY_METRICS.items():
        print(f"{metric:<15} {desc}")
    print()


def _get_fp_params(args) -> dict:
    """Build fingerprint parameters from arguments."""
    params = {}
    if args.radius:
        params["radius"] = args.radius
    if args.nbits:
        params["nBits"] = args.nbits
    return params


def _run_single_search(args):
    """Run similarity search for a single query."""
    from molprop_toolkit.similarity import similarity_search
    
    if args.verbose:
        print(f"Query: {args.query}", file=sys.stderr)
        print(f"Library: {args.library}", file=sys.stderr)
        print(f"Fingerprint: {args.fp_type}", file=sys.stderr)
        print(f"Metric: {args.metric}", file=sys.stderr)
    
    try:
        results = similarity_search(
            query=args.query,
            library=args.library,
            fp_type=args.fp_type,
            metric=args.metric,
            threshold=args.threshold,
            top_n=args.top or args.knn,
            smiles_col=args.smiles_col,
            id_col=args.id_col,
            fp_params=_get_fp_params(args),
            n_jobs=args.jobs,
            show_progress=args.progress,
        )
    except Exception as e:
        print(f"Error during search: {e}", file=sys.stderr)
        return 1
    
    # Output results
    df = results.to_dataframe()
    
    if args.verbose:
        print(f"\nSearched {results.n_searched} compounds ({results.n_valid} valid)",
              file=sys.stderr)
        print(f"Found {len(results)} results above threshold {args.threshold}",
              file=sys.stderr)
    
    if args.output:
        output_path = Path(args.output)
        if output_path.suffix == ".parquet":
            df.to_parquet(output_path, index=False)
        else:
            df.to_csv(output_path, index=False)
        print(f"Results saved to {output_path}", file=sys.stderr)
        
        if args.html:
            html_path = output_path.with_suffix(".html")
            _write_html_report(df, args.query, html_path, args)
            print(f"HTML report saved to {html_path}", file=sys.stderr)
    else:
        # Print to stdout
        print(df.to_string(index=False))
    
    return 0


def _run_batch_search(args):
    """Run similarity search for multiple queries."""
    from molprop_toolkit.similarity import similarity_search
    from molprop_toolkit.similarity.search import _load_smiles_file
    
    # Load queries
    query_path = Path(args.query)
    if query_path.suffix == ".csv":
        queries_df = pd.read_csv(query_path)
        if args.smiles_col in queries_df.columns:
            queries = queries_df[args.smiles_col].tolist()
            query_ids = queries_df.get(args.id_col, 
                                        [f"Query_{i+1}" for i in range(len(queries))])
        else:
            queries = queries_df.iloc[:, 0].tolist()
            query_ids = [f"Query_{i+1}" for i in range(len(queries))]
    else:
        queries_df = _load_smiles_file(query_path)
        queries = queries_df["SMILES"].tolist()
        query_ids = queries_df.get("Compound_ID", 
                                    [f"Query_{i+1}" for i in range(len(queries))])
    
    if args.verbose:
        print(f"Loaded {len(queries)} queries", file=sys.stderr)
    
    # Setup output directory
    if args.output:
        output_dir = Path(args.output)
        output_dir.mkdir(parents=True, exist_ok=True)
    else:
        output_dir = Path("similarity_results")
        output_dir.mkdir(exist_ok=True)
    
    # Run searches
    all_results = []
    for i, (query, query_id) in enumerate(zip(queries, query_ids)):
        if args.verbose:
            print(f"Searching query {i+1}/{len(queries)}: {query_id}", file=sys.stderr)
        
        try:
            results = similarity_search(
                query=query,
                library=args.library,
                fp_type=args.fp_type,
                metric=args.metric,
                threshold=args.threshold,
                top_n=args.top,
                smiles_col=args.smiles_col,
                id_col=args.id_col,
                query_id=str(query_id),
                fp_params=_get_fp_params(args),
                n_jobs=args.jobs,
            )
            
            df = results.to_dataframe()
            df.insert(0, "Query_ID", query_id)
            df.insert(1, "Query_SMILES", query)
            all_results.append(df)
            
        except Exception as e:
            print(f"Error searching {query_id}: {e}", file=sys.stderr)
            continue
    
    # Combine and save results
    if all_results:
        combined = pd.concat(all_results, ignore_index=True)
        output_path = output_dir / "batch_similarity_results.csv"
        combined.to_csv(output_path, index=False)
        print(f"Combined results saved to {output_path}", file=sys.stderr)
        
        if args.html:
            html_path = output_dir / "batch_similarity_results.html"
            _write_batch_html_report(combined, html_path, args)
            print(f"HTML report saved to {html_path}", file=sys.stderr)
    
    return 0


def _run_pairwise(library, args):
    """Calculate pairwise similarity matrix."""
    from molprop_toolkit.similarity import pairwise_similarity
    
    if args.verbose:
        print(f"Calculating pairwise similarity for {library}", file=sys.stderr)
    
    try:
        matrix, smiles_list = pairwise_similarity(
            molecules=library,
            fp_type=args.fp_type,
            metric=args.metric,
            smiles_col=args.smiles_col,
            fp_params=_get_fp_params(args),
            n_jobs=args.jobs,
            show_progress=args.progress,
        )
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    
    # Create DataFrame with SMILES as index
    df = pd.DataFrame(matrix, index=smiles_list, columns=smiles_list)
    
    if args.output:
        output_path = Path(args.output)
        df.to_csv(output_path)
        print(f"Similarity matrix saved to {output_path}", file=sys.stderr)
        print(f"Matrix shape: {matrix.shape}", file=sys.stderr)
    else:
        print(df.to_string())
    
    return 0


def _run_diversity(library, args):
    """Run diversity picking."""
    from molprop_toolkit.similarity import diversity_pick
    
    if args.verbose:
        print(f"Selecting {args.pick} diverse compounds from {library}", file=sys.stderr)
    
    seed = args.seed
    if seed and seed.isdigit():
        seed = int(seed)
    
    try:
        df = diversity_pick(
            molecules=library,
            n_pick=args.pick,
            seed=seed,
            fp_type=args.fp_type,
            metric=args.metric,
            smiles_col=args.smiles_col,
            id_col=args.id_col,
            fp_params=_get_fp_params(args),
            n_jobs=args.jobs,
            show_progress=args.progress,
        )
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    
    if args.output:
        output_path = Path(args.output)
        if output_path.suffix == ".parquet":
            df.to_parquet(output_path, index=False)
        else:
            df.to_csv(output_path, index=False)
        print(f"Diverse subset saved to {output_path}", file=sys.stderr)
    else:
        print(df.to_string(index=False))
    
    return 0


def _run_cluster(library, args):
    """Run similarity-based clustering."""
    from molprop_toolkit.similarity import cluster_by_similarity
    
    if args.verbose:
        print(f"Clustering {library} at threshold {args.threshold}", file=sys.stderr)
    
    try:
        df = cluster_by_similarity(
            molecules=library,
            threshold=args.threshold,
            fp_type=args.fp_type,
            metric=args.metric,
            smiles_col=args.smiles_col,
            id_col=args.id_col,
            fp_params=_get_fp_params(args),
            n_jobs=args.jobs,
        )
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    
    n_clusters = df["Cluster_ID"].max() + 1
    
    if args.output:
        output_path = Path(args.output)
        if output_path.suffix == ".parquet":
            df.to_parquet(output_path, index=False)
        else:
            df.to_csv(output_path, index=False)
        print(f"Clustering results saved to {output_path}", file=sys.stderr)
        print(f"Found {n_clusters} clusters", file=sys.stderr)
    else:
        print(df.to_string(index=False))
        print(f"\nFound {n_clusters} clusters", file=sys.stderr)
    
    return 0


def _write_html_report(df: pd.DataFrame, query: str, output_path: Path, args):
    """Generate HTML report for single query search."""
    
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Similarity Search Results</title>
    <style>
        :root {{
            --primary: #2563eb;
            --bg: #f8fafc;
            --card-bg: #ffffff;
            --text: #1e293b;
            --text-muted: #64748b;
            --border: #e2e8f0;
            --success: #22c55e;
        }}
        * {{ box-sizing: border-box; margin: 0; padding: 0; }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            background: var(--bg);
            color: var(--text);
            line-height: 1.6;
            padding: 2rem;
        }}
        .container {{ max-width: 1200px; margin: 0 auto; }}
        h1 {{ 
            font-size: 1.75rem; 
            font-weight: 600; 
            margin-bottom: 0.5rem;
            color: var(--text);
        }}
        .subtitle {{
            color: var(--text-muted);
            margin-bottom: 2rem;
        }}
        .query-card {{
            background: var(--card-bg);
            border: 1px solid var(--border);
            border-radius: 8px;
            padding: 1.5rem;
            margin-bottom: 2rem;
        }}
        .query-card h2 {{
            font-size: 0.875rem;
            text-transform: uppercase;
            letter-spacing: 0.05em;
            color: var(--text-muted);
            margin-bottom: 0.5rem;
        }}
        .query-smiles {{
            font-family: 'Courier New', monospace;
            font-size: 1rem;
            background: var(--bg);
            padding: 0.75rem 1rem;
            border-radius: 4px;
            word-break: break-all;
        }}
        .stats {{
            display: flex;
            gap: 2rem;
            margin-top: 1rem;
            font-size: 0.875rem;
            color: var(--text-muted);
        }}
        .stats span {{ font-weight: 600; color: var(--text); }}
        table {{
            width: 100%;
            background: var(--card-bg);
            border: 1px solid var(--border);
            border-radius: 8px;
            border-collapse: collapse;
            overflow: hidden;
        }}
        th, td {{
            padding: 0.75rem 1rem;
            text-align: left;
            border-bottom: 1px solid var(--border);
        }}
        th {{
            background: var(--bg);
            font-weight: 600;
            font-size: 0.75rem;
            text-transform: uppercase;
            letter-spacing: 0.05em;
            color: var(--text-muted);
        }}
        tr:last-child td {{ border-bottom: none; }}
        tr:hover td {{ background: var(--bg); }}
        .smiles {{
            font-family: 'Courier New', monospace;
            font-size: 0.875rem;
            max-width: 400px;
            overflow: hidden;
            text-overflow: ellipsis;
            white-space: nowrap;
        }}
        .similarity {{
            font-weight: 600;
            color: var(--primary);
        }}
        .high {{ color: var(--success); }}
        .footer {{
            margin-top: 2rem;
            text-align: center;
            font-size: 0.75rem;
            color: var(--text-muted);
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>Similarity Search Results</h1>
        <p class="subtitle">MolProp Toolkit</p>
        
        <div class="query-card">
            <h2>Query Compound</h2>
            <div class="query-smiles">{query}</div>
            <div class="stats">
                <div>Fingerprint: <span>{args.fp_type}</span></div>
                <div>Metric: <span>{args.metric}</span></div>
                <div>Threshold: <span>{args.threshold}</span></div>
                <div>Results: <span>{len(df)}</span></div>
            </div>
        </div>
        
        <table>
            <thead>
                <tr>
                    <th>Rank</th>
                    <th>Compound ID</th>
                    <th>SMILES</th>
                    <th>Similarity</th>
                </tr>
            </thead>
            <tbody>
"""
    
    for _, row in df.iterrows():
        sim = row["Similarity"]
        sim_class = "high" if sim >= 0.7 else ""
        html += f"""                <tr>
                    <td>{row['Rank']}</td>
                    <td>{row['Compound_ID']}</td>
                    <td class="smiles" title="{row['SMILES']}">{row['SMILES']}</td>
                    <td class="similarity {sim_class}">{sim:.4f}</td>
                </tr>
"""
    
    html += """            </tbody>
        </table>
        
        <p class="footer">Generated by MolProp Toolkit</p>
    </div>
</body>
</html>
"""
    
    with open(output_path, "w") as f:
        f.write(html)


def _write_batch_html_report(df: pd.DataFrame, output_path: Path, args):
    """Generate HTML report for batch search results."""
    
    queries = df.groupby(["Query_ID", "Query_SMILES"]).size().reset_index(name="Count")
    
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Batch Similarity Search Results</title>
    <style>
        :root {{
            --primary: #2563eb;
            --bg: #f8fafc;
            --card-bg: #ffffff;
            --text: #1e293b;
            --text-muted: #64748b;
            --border: #e2e8f0;
        }}
        * {{ box-sizing: border-box; margin: 0; padding: 0; }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            background: var(--bg);
            color: var(--text);
            line-height: 1.6;
            padding: 2rem;
        }}
        .container {{ max-width: 1400px; margin: 0 auto; }}
        h1 {{ font-size: 1.75rem; font-weight: 600; margin-bottom: 0.5rem; }}
        .subtitle {{ color: var(--text-muted); margin-bottom: 2rem; }}
        .summary {{
            display: flex;
            gap: 1rem;
            margin-bottom: 2rem;
        }}
        .stat-card {{
            background: var(--card-bg);
            border: 1px solid var(--border);
            border-radius: 8px;
            padding: 1rem 1.5rem;
            flex: 1;
        }}
        .stat-card .value {{
            font-size: 2rem;
            font-weight: 700;
            color: var(--primary);
        }}
        .stat-card .label {{
            font-size: 0.875rem;
            color: var(--text-muted);
        }}
        table {{
            width: 100%;
            background: var(--card-bg);
            border: 1px solid var(--border);
            border-radius: 8px;
            border-collapse: collapse;
            font-size: 0.875rem;
        }}
        th, td {{
            padding: 0.5rem 0.75rem;
            text-align: left;
            border-bottom: 1px solid var(--border);
        }}
        th {{
            background: var(--bg);
            font-weight: 600;
            font-size: 0.75rem;
            text-transform: uppercase;
            letter-spacing: 0.05em;
            color: var(--text-muted);
            position: sticky;
            top: 0;
        }}
        .smiles {{
            font-family: monospace;
            max-width: 300px;
            overflow: hidden;
            text-overflow: ellipsis;
            white-space: nowrap;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>Batch Similarity Search Results</h1>
        <p class="subtitle">MolProp Toolkit</p>
        
        <div class="summary">
            <div class="stat-card">
                <div class="value">{len(queries)}</div>
                <div class="label">Queries</div>
            </div>
            <div class="stat-card">
                <div class="value">{len(df)}</div>
                <div class="label">Total Hits</div>
            </div>
            <div class="stat-card">
                <div class="value">{args.fp_type}</div>
                <div class="label">Fingerprint</div>
            </div>
            <div class="stat-card">
                <div class="value">{args.threshold}</div>
                <div class="label">Threshold</div>
            </div>
        </div>
        
        <table>
            <thead>
                <tr>
                    <th>Query ID</th>
                    <th>Query SMILES</th>
                    <th>Rank</th>
                    <th>Hit ID</th>
                    <th>Hit SMILES</th>
                    <th>Similarity</th>
                </tr>
            </thead>
            <tbody>
"""
    
    for _, row in df.head(1000).iterrows():  # Limit for large results
        html += f"""                <tr>
                    <td>{row['Query_ID']}</td>
                    <td class="smiles" title="{row['Query_SMILES']}">{row['Query_SMILES']}</td>
                    <td>{row['Rank']}</td>
                    <td>{row['Compound_ID']}</td>
                    <td class="smiles" title="{row['SMILES']}">{row['SMILES']}</td>
                    <td>{row['Similarity']:.4f}</td>
                </tr>
"""
    
    html += """            </tbody>
        </table>
    </div>
</body>
</html>
"""
    
    with open(output_path, "w") as f:
        f.write(html)


if __name__ == "__main__":
    sys.exit(main())
