# Similarity Search Tool

`molprop-similarity` provides fingerprint-based similarity searching for molecular libraries. It supports multiple fingerprint types, similarity metrics, and operational modes including batch search, pairwise matrices, diversity picking, and clustering.

## Quick Start

```bash
# Basic similarity search
molprop-similarity "CCO" library.smi -o similar_to_ethanol.csv

# Find top 10 compounds similar to benzene with Tanimoto ≥ 0.7
molprop-similarity "c1ccccc1" library.csv --fp morgan --threshold 0.7 --top 10

# List available fingerprint types
molprop-similarity --list-fps

# List available similarity metrics
molprop-similarity --list-metrics
```

## Features

### Similarity Search

Find compounds similar to a query molecule:

```bash
# Basic search
molprop-similarity "CC(=O)Oc1ccccc1C(=O)O" results.csv --top 20

# With threshold filter
molprop-similarity "CC(=O)Oc1ccccc1C(=O)O" results.csv --threshold 0.6

# Using different fingerprint
molprop-similarity "CC(=O)Oc1ccccc1C(=O)O" results.csv --fp maccs

# Using different metric
molprop-similarity "CC(=O)Oc1ccccc1C(=O)O" results.csv --metric dice

# Generate HTML report
molprop-similarity "CC(=O)Oc1ccccc1C(=O)O" results.csv --top 50 -o results.csv --html
```

### Batch Search

Search multiple queries against a library:

```bash
# queries.smi contains one SMILES per line
molprop-similarity queries.smi library.csv --batch -o results/

# With all options
molprop-similarity queries.smi library.csv --batch --fp morgan --threshold 0.5 --top 10 -o results/ --html
```

### Pairwise Similarity Matrix

Calculate similarity between all pairs of compounds:

```bash
molprop-similarity --pairwise library.smi -o similarity_matrix.csv

# With specific fingerprint
molprop-similarity --pairwise library.smi --fp maccs -o matrix_maccs.csv
```

### Diversity Picking

Select a diverse subset of compounds using the MaxMin algorithm:

```bash
# Pick 100 diverse compounds
molprop-similarity --diversity library.smi --pick 100 -o diverse_100.csv

# With specific seed compound
molprop-similarity --diversity library.smi --pick 50 --seed "c1ccccc1" -o diverse_50.csv
```

### Clustering

Cluster compounds by similarity using leader-follower algorithm:

```bash
# Cluster at 70% similarity
molprop-similarity --cluster library.smi --threshold 0.7 -o clusters.csv

# More stringent clustering
molprop-similarity --cluster library.smi --threshold 0.85 -o clusters_strict.csv
```

### K-Nearest Neighbors

Find the K most similar compounds:

```bash
molprop-similarity "CCO" library.smi --knn 5
```

## Fingerprint Types

| Type | Description | Default Bits |
|------|-------------|--------------|
| `morgan` | Morgan circular fingerprint (ECFP-like) | 2048 |
| `morgan_feat` | Morgan with pharmacophoric features (FCFP-like) | 2048 |
| `maccs` | MACCS 166-bit structural keys | 167 |
| `rdkit` | RDKit topological fingerprint | 2048 |
| `atompair` | Atom pair fingerprint | 2048 |
| `torsion` | Topological torsion fingerprint | 2048 |
| `pattern` | Pattern fingerprint (substructure keys) | 2048 |

### Fingerprint Parameters

```bash
# Morgan with radius 3 and 1024 bits
molprop-similarity "CCO" library.smi --fp morgan --radius 3 --nbits 1024

# FCFP4 equivalent
molprop-similarity "CCO" library.smi --fp morgan_feat --radius 2
```

## Similarity Metrics

| Metric | Description |
|--------|-------------|
| `tanimoto` | Tanimoto coefficient (Jaccard index) - most common |
| `dice` | Dice coefficient - weights matches more heavily |
| `cosine` | Cosine similarity - treats fingerprints as vectors |
| `sokal` | Sokal-Michener coefficient - includes shared absences |
| `russel` | Russel-Rao coefficient - considers only shared bits |
| `kulczynski` | Kulczynski coefficient - asymmetric measure |
| `mcconnaughey` | McConnaughey coefficient - balanced measure |
| `asymmetric` | Tversky index (a=1, b=0) - substructure similarity |
| `braun_blanquet` | Braun-Blanquet coefficient |

## Input Formats

The tool accepts multiple input formats:

### SMILES Files

```
# name<tab>smiles format
Aspirin	CC(=O)Oc1ccccc1C(=O)O
Ibuprofen	CC(C)Cc1ccc(cc1)C(C)C(=O)O

# smiles<tab>name format  
CC(=O)Oc1ccccc1C(=O)O	Aspirin
CC(C)Cc1ccc(cc1)C(C)C(=O)O	Ibuprofen

# SMILES only
CC(=O)Oc1ccccc1C(=O)O
CC(C)Cc1ccc(cc1)C(C)C(=O)O
```

### CSV Files

```csv
Compound_ID,SMILES,Other_Column
Aspirin,CC(=O)Oc1ccccc1C(=O)O,value1
Ibuprofen,CC(C)Cc1ccc(cc1)C(C)C(=O)O,value2
```

### MolProp Results Tables

Works directly with output from `molprop-calc-v5`:

```bash
molprop-calc-v5 library.smi -o results.csv
molprop-similarity "CCO" results.csv --smiles-col Canonical_SMILES
```

## Output Format

### Search Results CSV

```csv
Rank,Compound_ID,SMILES,Similarity
1,Cpd_042,CCCO,0.8571
2,Cpd_017,CCCCO,0.7500
3,Cpd_089,CCO,1.0000
```

### Diversity Pick CSV

```csv
Compound_ID,SMILES,Pick_Order,Min_Distance
Cpd_001,c1ccccc1,1,0.0000
Cpd_089,CCCCCCCC,2,0.9231
Cpd_042,CC(=O)O,3,0.8750
```

### Cluster CSV

```csv
Compound_ID,SMILES,Cluster_ID,Is_Centroid
Cpd_001,c1ccccc1,0,True
Cpd_002,c1ccccc1O,0,False
Cpd_003,CCCCCC,1,True
```

## Performance Options

```bash
# Use 4 parallel workers (requires joblib)
molprop-similarity "CCO" library.smi -j 4

# Use all CPU cores
molprop-similarity "CCO" library.smi -j -1

# Show progress bar (requires tqdm)
molprop-similarity "CCO" library.smi --progress
```

## Python API

The similarity module can also be used programmatically:

```python
from molprop_toolkit.similarity import (
    similarity_search,
    pairwise_similarity,
    diversity_pick,
    cluster_by_similarity,
    get_fingerprint,
    tanimoto_similarity,
)

# Generate fingerprints
fp1 = get_fingerprint("CCO", fp_type="morgan")
fp2 = get_fingerprint("CCCO", fp_type="morgan")

# Calculate similarity
sim = tanimoto_similarity(fp1, fp2)
print(f"Similarity: {sim:.4f}")

# Search a library
results = similarity_search(
    query="CCO",
    library="library.smi",
    fp_type="morgan",
    threshold=0.7,
    top_n=10,
)

# Convert to DataFrame
df = results.to_dataframe()
print(df)

# Diversity picking
diverse = diversity_pick(
    molecules="library.smi",
    n_pick=50,
    fp_type="morgan",
)

# Clustering
clusters = cluster_by_similarity(
    molecules="library.smi",
    threshold=0.7,
    fp_type="morgan",
)
```

## Use Cases

### Lead Hopping

Find compounds similar to a lead but from different scaffolds:

```bash
# Search with moderate threshold
molprop-similarity "$LEAD_SMILES" library.csv --threshold 0.5 --top 100 -o lead_analogs.csv

# Then filter by scaffold in downstream analysis
```

### Library Deduplication

Remove near-duplicates from a compound library:

```bash
# Cluster at high similarity
molprop-similarity --cluster library.smi --threshold 0.95 -o clusters.csv

# Keep only centroids
```

### Diversity Analysis

Analyze chemical diversity of a library:

```bash
# Pairwise matrix for small libraries
molprop-similarity --pairwise library.smi -o matrix.csv

# Diversity picking for coverage analysis
molprop-similarity --diversity library.smi --pick 100 -o diverse_subset.csv
```

### SAR Exploration

Find all close analogs of a compound for SAR:

```bash
# High threshold for close analogs
molprop-similarity "$COMPOUND" library.csv --threshold 0.8 -o close_analogs.csv --html
```

## Tips

1. **Start with Morgan fingerprints** - they work well for most applications
2. **Use threshold 0.7-0.8** for finding similar compounds with same scaffold
3. **Use threshold 0.4-0.6** for finding related compounds with different scaffolds
4. **MACCS keys are faster** but less discriminative for large libraries
5. **Dice similarity gives higher values** than Tanimoto for same pairs
6. **Use parallel processing** (`-j -1`) for libraries > 10,000 compounds
