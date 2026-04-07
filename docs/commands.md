# MolScope command summary

This page is generated from `molprop_toolkit.command_metadata` and summarizes the public command surface.

## Core workflow

- `molscope calc`: Create the recommended extended results table from SMILES.
- `molscope analyze`: Inspect categories, summarize liabilities, and export focused slices.
- `molscope report`: Build a shareable Markdown and HTML report bundle.
- `molscope picklists`: Generate operational decision lists and optional HTML outputs.
- `molscope portal`: Assemble report, compare, SAR, MMP, picklist, and retro bundles into one workspace.

## Additional commands

- `molscope prep`: Apply named ligand-preparation presets and write a prep audit table.
- `molscope sketch`: Render 2D structures and browseable sketch sheets.
- `molscope integrate`: Merge measured assay endpoints into the canonical MolScope results table.
- `molscope compare`: Compare two campaign tables and summarize dataset, endpoint, and scaffold shifts.
- `molscope learnings`: Synthesize compare, SAR, MMP, and decision outputs into reusable campaign lessons.
- `molscope dashboard`: Accumulate learnings bundles into a persistent campaign dashboard and knowledge registry.
- `molscope sar`: Build a scaffold-centric SAR bundle with series summaries and top compounds.
- `molscope mmp`: Build matched-pair tables and aggregate transform summaries from a campaign table.
- `molscope search`: Search campaign tables by exact SMILES, SMARTS substructure, or motif panel.
- `molscope series`: Annotate scaffolds and similarity-based series clusters.
- `molscope similarity`: Run similarity search, matrices, diversity picking, and clustering.
- `molscope featurize`: Export ML-ready feature packs and optional deep-learning assets.
- `molscope retro`: Run offline retrosynthesis via AiZynthFinder and rerank route summaries.
- `molscope schema`: Validate results tables against the published column schema.

## Calculator profiles

- `extended`: Recommended default profile for full medchem and developability triage.
- `baseline`: Smaller medchem-focused profile for quicker, lighter-weight tables.

## Legacy aliases

- `molscope calc` → `molprop-calc-v5 (default) or molprop-calc-v4 (--profile baseline)`
- `molscope prep` → `molprop-prep`
- `molscope analyze` → `molprop-analyze`
- `molscope report` → `molprop-report`
- `molscope picklists` → `molprop-picklists`
- `molscope portal` → `molprop-portal`
- `molscope sketch` → `molprop-sketch`
- `molscope integrate` → `molprop-integrate`
- `molscope compare` → `molprop-compare`
- `molscope learnings` → `molprop-learnings`
- `molscope dashboard` → `molprop-dashboard`
- `molscope sar` → `molprop-sar`
- `molscope mmp` → `molprop-mmp`
- `molscope search` → `molprop-search`
- `molscope series` → `molprop-series`
- `molscope similarity` → `molprop-similarity`
- `molscope featurize` → `molprop-featurize`
- `molscope retro` → `molprop-retro`
- `molscope schema` → `python tools/validate_csv_schema.py`
