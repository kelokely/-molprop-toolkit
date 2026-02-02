# MolProp Toolkit — Similarity Search

`molprop-similarity` is now integrated into the MolProp Toolkit package (no manual copying of scripts required).

## What’s included

- Python module: `molprop_toolkit/similarity/`
  - fingerprint generation (`fingerprints.py`)
  - similarity metrics (`metrics.py`)
  - similarity search / pairwise / diversity pick / clustering (`search.py`)

- CLI entry point: `molprop-similarity`
  - implemented in `molprop_toolkit/tools/similarity_cli.py`
  - registered in `pyproject.toml` under `[project.scripts]`

- Docs: `docs/tools/molprop-similarity.html` (short tool page)

## Install (recommended)

```bash
conda env create -f environment.yml
conda activate molprop-toolkit
pip install -e ".[dev,parallel]"
```

RDKit is required for similarity search:

```bash
conda install -c conda-forge rdkit
```

## Quick usage

```bash
molprop-similarity --help
molprop-similarity --list-fps
molprop-similarity --list-metrics

# Basic search
molprop-similarity "CCO" library.smi --top 25 -o hits.csv

# Diversity picking
molprop-similarity --diversity library.smi --pick 100 -o diverse.csv

# Clustering
molprop-similarity --cluster library.smi --threshold 0.7 -o clusters.csv
```

## Notes

Older experimental scripts that lived under the top-level `tools/` folder have been removed to avoid drift. The supported CLI is
`molprop_toolkit/tools/similarity_cli.py` via the `molprop-similarity` console script.

