# MolScope Toolkit

[![CI](https://github.com/kelokely/-molprop-toolkit/actions/workflows/ci.yml/badge.svg)](https://github.com/kelokely/-molprop-toolkit/actions/workflows/ci.yml)
[![Releases](https://img.shields.io/github/v/release/kelokely/-molprop-toolkit)](https://github.com/kelokely/-molprop-toolkit/releases)
[![Downloads](https://img.shields.io/badge/download-latest-blue)](https://github.com/kelokely/-molprop-toolkit/releases/latest/download/molprop-toolkit_latest.zip)

MolScope Toolkit is a small, practical repository for turning SMILES lists into analysis-ready medicinal chemistry property tables and then generating category-focused summaries, reports, and operational picklists from those tables. It is designed for the early decision-making loop where you want the “medchem-first” columns surfaced up front, you want identity preserved in every output, and you want a repeatable way to move from a raw SMILES list to a shareable triage package.

Created by Elokely Lab at the University of Wyoming School of Pharmacy: https://www.uwyo.edu/pharmacy/faculty-and-staff-directory/faculty/khaled-elokely.html

The project is organized around two ideas. First, compute a wide descriptor table where the first columns are the ones medicinal chemists typically reach for when scanning a dataset, including CNS MPO, oral/RO5-style filters, structural alert flags, toxicity and hERG risk indicators, metabolism soft-spot summaries, solubility and permeability heuristics, and simple PK-oriented signals. Second, treat that wide CSV as a dataset and provide lightweight analyzers, report builders, and picklist generators that are consistent about columns, consistent about identifiers, and easy to integrate into everyday workflows.

## Quickstart

If you only want to evaluate the repository quickly, run the example file shipped in `examples/example.smi`. The recommended public workflow now starts from a single top-level command: calculate a results table, inspect categories, then generate outputs that you can share.

```bash
# 1) Create an environment (recommended for RDKit reliability)
conda env create -f environment.yml
conda activate molprop-toolkit

# 2) Install console scripts + Parquet support
pip install -e ".[dev,parquet]"

# 3) Optional: lock a preparation policy first
molscope prep examples/example.smi -o prepared.csv --preset triage

# 4) Calculate a results table
molscope calc examples/example.smi -o results.csv

# 5) Inspect available analysis categories
molscope analyze results.csv --list

# 6) Generate a report bundle (Markdown + HTML + plots)
molscope report results.csv

# 7) Generate picklists (CSV and optional HTML)
molscope picklists results.csv --html
```

If you already have a compatible results CSV from another pipeline, you can skip RDKit and run analyzers and report tooling directly.

Compatibility aliases remain available, but `molscope` is now the simplest entry point to the canonical workflow.

A compact command summary is also available in [docs/commands.md](docs/commands.md).
Scientific grounding and release-gate pages are published in [docs/docs/methods.html](docs/docs/methods.html), [docs/docs/validation.html](docs/docs/validation.html), [docs/docs/benchmarks.html](docs/docs/benchmarks.html), and [docs/docs/uncertainty.html](docs/docs/uncertainty.html).

If you want to validate the second-level campaign workflows on a richer example than the single-file quickstart, the synthetic Campaign Alpha pack lives under `examples/campaign_alpha/` and can be built end to end with:

```bash
python examples/build_campaign_alpha_demo.py
```

## Repository layout

The calculators live in `calculators/` and produce a wide results CSV from an input SMILES file. The analyzers live in `analyzers/` and read a results CSV to print category summaries and optionally export smaller category CSVs. The “packaging” tools that generate shareable artifacts live in `tools/` (report bundles, picklists, structure sketches, series analytics, and retrosynthesis wrappers). The `docs/` folder contains schema and narrative documentation, and it is also used to publish the GitHub Pages site at https://kelokely.github.io/-molprop-toolkit/.

## Installation

The core analysis layer depends only on `pandas` and `numpy`. The property calculators require RDKit (and a few optional scientific/cheminformatics libraries). RDKit is intentionally not forced into `requirements.txt` because it is environment dependent and is most reliably installed via conda.

If you want the command-line tools (`molscope` plus compatibility aliases for existing scripts), install the repository in editable mode.

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

MolScope Toolkit assumes a simple workflow: an input SMILES file produces a single wide CSV, and all downstream tools treat that CSV as the canonical dataset. The analyzers and generators are intentionally tolerant of missing columns so that you can use them with partial tables, but the “best experience” is achieved when the results include a stable compound identifier column and a canonical SMILES column.

SMILES input is flexible. Files like `name<tab>smiles`, `smiles<tab>name`, `name smiles`, or `smiles name` are supported, and the calculators try to detect which token is the SMILES by looking for typical SMILES characters. The results CSV should contain a compound identifier column, ideally `Compound_ID`. The analyzers attempt to detect an ID column automatically and will keep that identity visible in terminal output and exported tables.

## Standard workflow (recommended sequence)

Start with a calculator, then use analyzers and generators to produce outputs for humans and downstream systems.

```bash
# Calculate a wide table
molscope calc my_library.smi -o results.csv

# Targeted category summaries
molscope analyze results.csv --category cns_mpo
molscope analyze results.csv --category toxicity

# Shareable report bundle
molscope report results.csv

# Quick 2D depictions + browseable HTML index
molscope sketch results.csv --outdir sketches --format svg --html

# Operational picklists (filter/sort outputs)
molscope picklists results.csv --html

# Series analytics (scaffolds + similarity clusters)
molscope series results.csv
```

For the calculator profile, `molscope calc` defaults to the recommended extended table. If you want the smaller baseline table, use:

```bash
molscope calc --profile baseline my_library.smi -o results.csv
```

## Structure handling options (standardization, tautomers, stereo, ionization)

The calculators perform an RDKit-based preparation step to make downstream tables more consistent when input lists contain salts, mixtures, charge variants, or tautomeric forms. Preparation writes traceability metadata into columns such as `Input_Canonical_SMILES`, `Canonical_SMILES`, and `Prep_*` columns so that you can debug what changed and why. You can disable preparation or parts of it if you need a minimal table or want to preserve input forms.

If you want preparation to be a first-class workflow step before descriptor calculation, use `molscope prep`. It writes a preparation-focused table plus a summary JSON sidecar so teams can standardize the structure-of-record policy before they move into the full calculator workflow.

```bash
# standard 2D triage
molscope prep examples/example.smi -o prepared.csv --preset triage

# ambiguity audit
molscope prep examples/example.smi -o prepared_audit.csv --preset audit

# 3D-ready protomer-aware preparation
molscope prep examples/example.smi -o prepared_3d.csv --preset 3d

# minimal / externally prepared structures
molscope prep examples/example.smi -o prepared_external.csv --preset external
```

`molprop-prep` is also installed as a direct console script, but `molscope prep` is the recommended public entry point.

If you need fine-grained control inside the descriptor calculator itself, the lower-level `molscope calc` flags still work:

```bash
# Disable preparation entirely
molscope calc examples/example.smi -o results.csv --no-prep

# Disable tautomer canonicalization
molscope calc examples/example.smi -o results.csv --tautomer-mode none

# Keep charges (skip uncharging)
molscope calc examples/example.smi -o results.csv --keep-charges
```

Tautomers and stereochemistry can be treated explicitly when you need deterministic behavior or auditability. You can keep, strip, or enumerate stereoisomers for unassigned stereo, and you can use a single canonical tautomer (fast) or enumerate tautomers (more expensive) while recording summary columns.

```bash
# Enumerate tautomers, record the set size, and select a representative
molscope calc examples/example.smi -o results.csv --tautomer-mode enumerate --tautomer-max 64 --tautomer-topk 5

# Enumerate stereoisomers for unassigned stereo and select a representative
molscope calc examples/example.smi -o results.csv --stereo-mode enumerate --stereo-max 32 --stereo-topk 5 --stereo-select canonical
```

Ionization handling is available in two modes. The default is a lightweight, pH-aware feature set that adds `Ion_*` columns without changing the structure used for descriptor calculation. If you need explicit protomer enumeration, Dimorphite-DL can be enabled as an optional dependency.

```bash
# Add pH-dependent ionization features (no structure change)
molscope calc examples/example.smi -o results.csv --ph 7.4 --ionization heuristic

# Disable ionization features entirely
molscope calc examples/example.smi -o results.csv --ionization none
```

For Dimorphite-DL enumeration:

```bash
pip install dimorphite-dl
molscope calc examples/example.smi -o results.csv --ionization dimorphite --ph 7.4
```

## Calculators

The calculators generate a wide results table and keep the most decision-relevant columns near the front.

`molscope calc --profile baseline` is the baseline table focused on 2D descriptors, drug-likeness style filters, alerts, and interpretability. `molscope calc` defaults to the extended profile, which adds solubility and permeability heuristics, PK-oriented features, and optional 3D descriptors.

```bash
molscope calc --profile baseline examples/example.smi -o results.csv
molscope calc examples/example.smi -o results.csv
```

3D features are optional because they require conformer generation and can significantly increase runtime.

```bash
molscope calc examples/example.smi -o results.csv --3d --3d-num-confs 10 --3d-minimize mmff
```

## Reports, picklists, and series analytics

`molscope report` builds a timestamped report bundle (Markdown + HTML + plots) under `reports/`. `molscope picklists` produces CSV picklists and can also write a small HTML site for browsing picklist outputs. `molscope series` runs series-level enrichment (scaffolds + similarity clusters) and writes joinable tables plus an enriched CSV.

```bash
molscope report results.csv
molscope picklists results.csv --html
molscope series results.csv
```

When you want a single handoff bundle for a project review, point `molscope portal` at any combination of previously generated artifact directories:

```bash
molscope report results.csv --outdir reports/demo
molscope picklists results.csv --outdir picklists/demo --html
molscope learnings --picklists-dir picklists/demo -o learnings/demo
molscope dashboard --learnings-dir learnings/demo -o dashboard/campaign_alpha
molscope portal --report-dir reports/demo --picklists-dir picklists/demo --learnings-dir learnings/demo --dashboard-dir dashboard/campaign_alpha -o portal/demo
```

`molscope learnings` turns compare, SAR, MMP, and decision outputs into reusable project lessons (`learning_cards.csv`, `summary.md`, `index.html`) so teams can capture what changed, what worked, and what should be tested next. `molscope dashboard` then accumulates those learnings bundles across rounds into a persistent registry (`snapshot_registry.csv`, `knowledge_registry.csv`, `recurring_lessons.csv`, `index.html`). `molscope portal` can include both the current learnings bundle and the cross-round dashboard so teams get one entry point instead of a pile of folders.

You can also frame those outputs around a campaign stage:

```bash
molscope report results.csv --stage lead-optimization
molscope picklists results.csv --stage hit-expansion --html
molscope picklists retrosynthesis/demo/retro_enriched.csv --stage backup-series --html
```

The current stage presets are `hit-expansion`, `lead-optimization`, `backup-series`, and `liability-de-risking`. They do two things:
- add stage-aware framing and a priority preview to reports
- choose a smaller, stage-appropriate default subset of built-in picklists when you do not explicitly pass `--picklists`

The newer decision-oriented picklists are designed to read more like operational queues than filtered exports: `make_next`, `test_next`, `de_risk_next`, `hold`, and `backup_series`. Their CSV and HTML outputs now carry `Decision_Queue`, `Decision_Rationale`, `Decision_Next_Action`, `Decision_Confidence`, `Decision_Evidence_Summary`, `Decision_Hypothesis`, `Decision_Minimal_Test_Set`, and `Decision_Minimal_Analog_Set` fields so the outputs are easier to use directly in project reviews and follow-up planning.

`molscope report` now mirrors that decision layer. In addition to the queue summary, the report includes experiment prompts that state the current hypothesis, the smallest next test set that should resolve the uncertainty, and the narrowest analog set worth making before the next round.

## Preclinical profiling

MolScope Toolkit can also be read as a lightweight preclinical profiling layer built from structure-derived heuristics. The fastest way to use it is to separate the outputs into in vitro-relevant modules and in vivo-relevant modules, then use the analyzers, reports, and picklists to turn those signals into concrete experiment priorities.

The in vitro-relevant modules are the ones most useful for assay planning and early liability review: solubility, permeability, toxicity alerts, med-chem flags, hERG risk, CYP inhibition burden, and metabolism soft spots. These are ranking and communication tools rather than replacements for experiments, but they help answer which compounds should move first into solubility assays, microsomal stability, CYP inhibition panels, hERG counterscreens, or orthogonal artifact checks.

The in vivo-relevant modules are the PK-oriented summaries and developability views: protein binding, distribution class, clearance risk, oral exposure class, and composite developability burden. In practice these are most useful when combined with `molscope report`, `molscope picklists`, and `molscope series`, because those tools make it easier to compare compounds or scaffolds before committing to more expensive PK and efficacy studies.

## Optional: offline retrosynthesis planning

`molscope retro` is an offline wrapper around AiZynthFinder. It delegates planning to `aizynthcli` and post-processes outputs into a reranked route summary bundle (`routes_summary.csv`, `routes_summary.json`, `summary.md`, `index.html`) and, when requested, a `retro_enriched.csv` table that can feed reports and picklists. This feature is intentionally optional because AiZynthFinder has heavier dependencies and currently supports Python 3.10–3.12.

```bash
# Download AiZynthFinder public model/stock bundle (create the folder first)
mkdir -p aizynth_data
download_public_data aizynth_data

# Run retrosynthesis and write a browseable HTML output folder
molscope retro results.csv --config aizynth_data/config.yml --outdir retrosynthesis/demo

# Prefer stock-supported routes and write a downstream-friendly retro table
molscope retro results.csv --config aizynth_data/config.yml --rerank-policy stock-first --write-retro-table

# Use the retro-enriched campaign table in downstream decision products
molscope report retrosynthesis/demo/retro_enriched.csv
molscope picklists retrosynthesis/demo/retro_enriched.csv --picklists make_next_synthesis_ready,route_risk_deprioritize --html
```

If you want to host a demo report on GitHub Pages, you can set the output directory under `docs/` and commit it, but be aware that route images and large route sets can make the repository grow quickly.

## CSV schema and documentation

MolScope Toolkit ships a lightweight schema at `docs/schema.json`. The schema is used to document columns and can be used to validate that a results CSV contains the expected fields. The GitHub Pages site mirrors the docs in a browsable format at https://kelokely.github.io/-molprop-toolkit/.

```bash
# Validate a CSV against the schema
molscope schema results.csv

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
