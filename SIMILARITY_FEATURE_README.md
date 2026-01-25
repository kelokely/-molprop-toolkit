# MolProp Toolkit - Similarity Search Feature

This directory contains the new similarity search module for MolProp Toolkit.

## Files to Add to Your Repository

```
molprop_toolkit/
├── similarity/
│   ├── __init__.py          # Module exports
│   ├── fingerprints.py      # Fingerprint generation
│   ├── metrics.py           # Similarity metrics
│   └── search.py            # Search functions
└── tools/
    └── similarity_cli.py    # CLI entry point

tools/
└── molprop_similarity.py    # Full CLI implementation

tests/
└── test_similarity.py       # Comprehensive tests

docs/tools/
└── similarity.md            # Documentation

examples/
└── similarity_example.smi   # Example compounds
```

## Integration Steps

### 1. Copy the similarity module

```bash
# From this directory, copy to your repo
cp -r molprop_toolkit/similarity /path/to/your/repo/molprop_toolkit/
cp molprop_toolkit/tools/similarity_cli.py /path/to/your/repo/molprop_toolkit/tools/
cp tools/molprop_similarity.py /path/to/your/repo/tools/
```

### 2. Update pyproject.toml

Add this entry to `[project.scripts]`:

```toml
molprop-similarity = "molprop_toolkit.tools.similarity_cli:main"
```

Add similarity to packages in `[tool.setuptools]`:

```toml
packages = ["molprop_toolkit", "molprop_toolkit.similarity", ...]
```

Add optional dependencies:

```toml
[project.optional-dependencies]
parallel = [
    "joblib>=1.2.0",
    "tqdm>=4.64.0",
]
```

### 3. Reinstall the package

```bash
pip install -e '.[dev,parallel]'
```

### 4. Run tests

```bash
pytest tests/test_similarity.py -v
```

### 5. Verify CLI

```bash
molprop-similarity --help
molprop-similarity --list-fps
molprop-similarity "CCO" examples/similarity_example.smi --top 5
```

## Quick Usage Examples

```bash
# Basic similarity search
molprop-similarity "c1ccccc1O" library.smi -o phenol_analogs.csv

# Search with threshold
molprop-similarity "c1ccccc1O" library.smi --threshold 0.6 --top 20

# Different fingerprint
molprop-similarity "c1ccccc1O" library.smi --fp maccs

# HTML report
molprop-similarity "c1ccccc1O" library.smi --top 50 -o results.csv --html

# Batch search
molprop-similarity queries.smi library.smi --batch -o results/

# Diversity picking
molprop-similarity --diversity library.smi --pick 100 -o diverse.csv

# Clustering
molprop-similarity --cluster library.smi --threshold 0.7 -o clusters.csv
```

## Python API

```python
from molprop_toolkit.similarity import (
    similarity_search,
    diversity_pick,
    get_fingerprint,
    tanimoto_similarity,
)

# Search
results = similarity_search("CCO", "library.smi", threshold=0.7)
df = results.to_dataframe()

# Fingerprints
fp = get_fingerprint("CCO", fp_type="morgan")
```

## Features

- **7 fingerprint types**: Morgan, MACCS, RDKit, AtomPair, Torsion, Pattern, Morgan+Features
- **9 similarity metrics**: Tanimoto, Dice, Cosine, Sokal, Russel, Kulczynski, McConnaughey, Asymmetric, Braun-Blanquet
- **Multiple modes**: Single query, batch, pairwise matrix, diversity picking, clustering
- **Flexible input**: SMILES files, CSV, Parquet, DataFrames, lists
- **Parallel processing**: Optional joblib support for large libraries
- **HTML reports**: Browseable output for results

## Version

This feature is designed for MolProp Toolkit v0.6.0.
