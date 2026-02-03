# MolProp Toolkit

[![CI](https://github.com/kelokely/-molprop-toolkit/actions/workflows/ci.yml/badge.svg)](https://github.com/kelokely/-molprop-toolkit/actions/workflows/ci.yml)
[![Releases](https://img.shields.io/github/v/release/kelokely/-molprop-toolkit)](https://github.com/kelokely/-molprop-toolkit/releases)
[![Downloads](https://img.shields.io/badge/download-latest-blue)](https://github.com/kelokely/-molprop-toolkit/releases/latest/download/molprop-toolkit_latest.zip)

MolProp Toolkit is a small, practical repository for turning SMILES lists into analysis-ready medicinal chemistry property tables and then generating category-focused summaries, reports, and operational picklists from those tables. It is designed for the early decision-making loop where you want the “medchem-first” columns surfaced up front, you want identity preserved in every output, and you want a repeatable way to move from a raw SMILES list to a shareable triage package.

Created by Elokely Lab at the University of Wyoming School of Pharmacy: https://www.uwyo.edu/pharmacy/faculty-and-staff-directory/faculty/khaled-elokely.html

The project is organized around two ideas. First, compute a wide descriptor table where the first columns are the ones medicinal chemists typically reach for when scanning a dataset, including CNS MPO, oral/RO5-style filters, structural alert flags, toxicity and hERG risk indicators, metabolism soft-spot summaries, solubility and permeability heuristics, and simple PK-oriented signals. Second, treat that wide CSV as a dataset and provide lightweight analyzers, report builders, and picklist generators that are consistent about columns, consistent about identifiers, and easy to integrate into everyday workflows.

## Quickstart

If you only want to evaluate the repository quickly, run the example file shipped in `examples/example.smi`. The sequence is always the same: calculate a results table, inspect categories, then generate outputs that you can share.

```bash
# 1) Create an environment (recommended for RDKit reliability)
conda env create -f environment.yml
conda activate molprop-toolkit

# 2) Install console scripts + Parquet support
pip install -e ".[dev,parquet]"

# 3) Calculate a results table
molprop-calc-v5 examples/example.smi -o results.csv

# 4) Inspect available analysis categories
molprop-analyze results.csv --list

# 5) Generate a report bundle (Markdown + HTML + plots)
molprop-report results.csv

# 6) Generate picklists (CSV and optional HTML)
molprop-picklists results.csv --html
```

If you already have a compatible results CSV from another pipeline, you can skip RDKit and run analyzers and report tooling directly.

## Repository layout

The calculators live in `calculators/` and produce a wide results CSV from an input SMILES file. The analyzers live in `analyzers/` and read a results CSV to print category summaries and optionally export smaller category CSVs. The “packaging” tools that generate shareable artifacts live in `tools/` (report bundles, picklists, structure sketches, series analytics, and retrosynthesis wrappers). The `docs/` folder contains schema and narrative documentation, and it is also used to publish the GitHub Pages site at https://kelokely.github.io/-molprop-toolkit/.

## Installation

The core analysis layer depends only on `pandas` and `numpy`. The property calculators require RDKit (and a few optional scientific/cheminformatics libraries). RDKit is intentionally not forced into `requirements.txt` because it is environment dependent and is most reliably installed via conda.

If you want the command-line tools (`molprop-calc-v5`, `molprop-analyze`, `molprop-report`, `molprop-sketch`, `molprop-series`, `molprop-picklists`, `molprop-retro`), install the repository in editable mode.

```bash
pip install -e .
```

For development and tests:

```bash
pip install -e '.[dev]'
```

If you want calculators, install RDKit using conda-forge.

```bash
conda install -c conda-forge rdkit
```

## Input and output conventions

MolProp Toolkit assumes a simple workflow: an input SMILES file produces a single wide CSV, and all downstream tools treat that CSV as the canonical dataset. The analyzers and generators are intentionally tolerant of missing columns so that you can use them with partial tables, but the “best experience” is achieved when the results include a stable compound identifier column and a canonical SMILES column.

SMILES input is flexible. Files like `name<tab>smiles`, `smiles<tab>name`, `name smiles`, or `smiles name` are supported, and the calculators try to detect which token is the SMILES by looking for typical SMILES characters. The results CSV should contain a compound identifier column, ideally `Compound_ID`. The analyzers attempt to detect an ID column automatically and will keep that identity visible in terminal output and exported tables.

## Standard workflow (recommended sequence)

Start with a calculator, then use analyzers and generators to produce outputs for humans and downstream systems.

```bash
# Calculate a wide table
molprop-calc-v5 my_library.smi -o results.csv

# Targeted category summaries
molprop-analyze results.csv --category cns_mpo
molprop-analyze results.csv --category toxicity

# Shareable report bundle
molprop-report results.csv

# Quick 2D depictions + browseable HTML index
molprop-sketch results.csv --outdir sketches --format svg --html

# Operational picklists (filter/sort outputs)
molprop-picklists results.csv --html

# Series analytics (scaffolds + similarity clusters)
molprop-series results.csv
```

## Structure handling options (standardization, tautomers, stereo, ionization)

The calculators perform an RDKit-based preparation step to make downstream tables more consistent when input lists contain salts, mixtures, charge variants, or tautomeric forms. Preparation writes traceability metadata into columns such as `Input_Canonical_SMILES`, `Canonical_SMILES`, and `Prep_*` columns so that you can debug what changed and why. You can disable preparation or parts of it if you need a minimal table or want to preserve input forms.

```bash
# Disable preparation entirely
molprop-calc-v5 examples/example.smi -o results.csv --no-prep

# Disable tautomer canonicalization
molprop-calc-v5 examples/example.smi -o results.csv --tautomer-mode none

# Keep charges (skip uncharging)
molprop-calc-v5 examples/example.smi -o results.csv --keep-charges
```

Tautomers and stereochemistry can be treated explicitly when you need deterministic behavior or auditability. You can keep, strip, or enumerate stereoisomers for unassigned stereo, and you can use a single canonical tautomer (fast) or enumerate tautomers (more expensive) while recording summary columns.

```bash
# Enumerate tautomers, record the set size, and select a representative
molprop-calc-v5 examples/example.smi -o results.csv --tautomer-mode enumerate --tautomer-max 64 --tautomer-topk 5

# Enumerate stereoisomers for unassigned stereo and select a representative
molprop-calc-v5 examples/example.smi -o results.csv --stereo-mode enumerate --stereo-max 32 --stereo-topk 5 --stereo-select canonical
```

Ionization handling is available in two modes. The default is a lightweight, pH-aware feature set that adds `Ion_*` columns without changing the structure used for descriptor calculation. If you need explicit protomer enumeration, Dimorphite-DL can be enabled as an optional dependency.

```bash
# Add pH-dependent ionization features (no structure change)
molprop-calc-v5 examples/example.smi -o results.csv --ph 7.4 --ionization heuristic

# Disable ionization features entirely
molprop-calc-v5 examples/example.smi -o results.csv --ionization none
```

For Dimorphite-DL enumeration:

```bash
pip install dimorphite-dl
molprop-calc-v5 examples/example.smi -o results.csv --ionization dimorphite --ph 7.4
```

## Calculators

The calculators generate a wide results table and keep the most decision-relevant columns near the front.

`molprop-calc-v4` is the baseline table focused on 2D descriptors, drug-likeness style filters, alerts, and interpretability. `molprop-calc-v5` extends v4 with additional families such as solubility and permeability heuristics, PK-oriented features, and optional 3D descriptors.

```bash
molprop-calc-v4 examples/example.smi -o results.csv
molprop-calc-v5 examples/example.smi -o results.csv
```

3D features are optional because they require conformer generation and can significantly increase runtime.

```bash
molprop-calc-v5 examples/example.smi -o results.csv --3d --3d-num-confs 10 --3d-minimize mmff
```

## Reports, picklists, and series analytics

`molprop-report` builds a timestamped report bundle (Markdown + HTML + plots) under `reports/`. `molprop-picklists` produces CSV picklists and can also write a small HTML site for browsing picklist outputs. `molprop-series` runs series-level enrichment (scaffolds + similarity clusters) and writes joinable tables plus an enriched CSV.

```bash
molprop-report results.csv
molprop-picklists results.csv --html
molprop-series results.csv
```

## Optional: offline retrosynthesis planning

`molprop-retro` is an offline wrapper around AiZynthFinder. It delegates planning to `aizynthcli` and post-processes outputs into a summary CSV/JSON plus an HTML route browser that reuses the MolProp Toolkit GitHub Pages styling. This feature is intentionally optional because AiZynthFinder has heavier dependencies and currently supports Python 3.10–3.12.

```bash
# Download AiZynthFinder public model/stock bundle (create the folder first)
mkdir -p aizynth_data
download_public_data aizynth_data

# Run retrosynthesis and write a browseable HTML output folder
molprop-retro results.csv --config aizynth_data/config.yml --outdir retrosynthesis/demo
```

If you want to host a demo report on GitHub Pages, you can set the output directory under `docs/` and commit it, but be aware that route images and large route sets can make the repository grow quickly.

## CSV schema and documentation

MolProp Toolkit ships a lightweight schema at `docs/schema.json`. The schema is used to document columns and can be used to validate that a results CSV contains the expected fields. The GitHub Pages site mirrors the docs in a browsable format at https://kelokely.github.io/-molprop-toolkit/.

```bash
# Validate a CSV against the schema
python tools/validate_csv_schema.py results.csv

# Generate the Markdown column reference from the schema
python tools/generate_columns_md_from_schema.py
```

The primary docs entry points are `docs/columns_schema.md` (generated column reference), `docs/feature_families.md` (prefix-based feature grouping), and `docs/categories.md` (narrative interpretation notes).

## Notes on missing values (`--fill-missing`)

Many descriptor families are only defined when a corresponding structural feature exists. In a raw table, those cells should be treated as “not applicable” rather than numeric zeros. Some modeling pipelines require dense matrices, so the calculators provide `--fill-missing` to replace missing values with a constant. If you choose to use `--fill-missing 0`, treat it as a modeling decision and document it explicitly, because it converts “not applicable” into a numeric value that some models may interpret as meaningful signal.

## Reproducibility and disclaimers

The outputs in this repository are heuristic and intended for triage and prioritization. They are not a substitute for experimental ADME/Tox, validated predictive models, or expert review. Alert- and rule-based features can overflag or underflag depending on chemistry, dose, exposure, and context. Use the toolkit as a screening and communication layer, then confirm key decisions with appropriate experimental or validated computational follow-up.

## References

Wager, T. T. et al., “Defining Desirable Central Nervous System Drug Space through the Alignment of Molecular Properties, In Vitro ADME, and Safety Attributes,” ACS Chemical Neuroscience (2010). https://pubs.acs.org/doi/10.1021/cn100007x

Bickerton, G. R. et al., “Quantifying the chemical beauty of drugs,” Nature Chemistry (2012). https://www.nature.com/articles/nchem.1243

