# Column conventions and grouping

MolProp Toolkit uses simple naming conventions so that analysis scripts can discover groups of columns even if you add your own custom properties. The analyzers are tolerant to missing columns, so you can run them on partial outputs. The master analyzer prefers explicit column lists for the most important metrics, but several analyzers will also include all columns that share a prefix.

## Feature-to-columns reference

A compact “column families” cheat sheet that groups common `*_` prefixes and explains how to interpret them is available here: `docs/feature_families.md`.

## Identifier columns

The analyzers attempt to display a compound identifier column in every table they print. The preferred column name is `Compound_ID`. If it is not present, most analyzers fall back to the first column in the CSV, which is a common pattern when the first column contains names. If your dataset uses a different name, such as `Name`, you can rename it to `Compound_ID` in the CSV for best compatibility.

## Preparation metadata columns

When you run the calculators with structure preparation enabled (default), additional traceability columns may be present.

`Input_Canonical_SMILES` is the canonical SMILES for the molecule as RDKit parsed it from the input. `Canonical_SMILES` is the canonical SMILES after preparation. The `Prep_*` fields describe which operations changed the structure and whether stereochemistry is fully specified.

Common fields include `Prep_SaltsStripped`, `Prep_TautomerChanged`, `Prep_StereoStatus`, `Prep_UnassignedChiral`, and `Prep_Notes`.

## Ionization metadata columns

When heuristic ionization is enabled in the calculators (default), additional columns prefixed with `Ion_` may be present. These are computed from estimated pKa values using Henderson–Hasselbalch at the chosen pH.

Typical fields include `Ion_Mode` (none/heuristic/dimorphite), `Ion_pH`, `Ion_State`, `Ion_NetCharge_Est`, and `Ion_LogD_pH_Est`. These are intended as ranking signals and for selecting among enumerated protomers, not as a replacement for rigorous pKa/protomer modeling.

## Stereochemistry metadata columns (Task 4)

The calculators may add `Stereo_*` columns when stereochemistry handling is enabled (default mode is `keep`). These columns describe whether stereochemistry was preserved, stripped, or enumerated.

In `strip` mode, stereochemistry is removed before calculations so the structure is treated as unspecified/racemic.

In `enumerate` mode, the toolkit enumerates stereoisomers for unassigned stereocenters/bonds (bounded by a maximum) and selects a representative stereoisomer for calculations, while recording the total count and a short list of example stereoisomer SMILES.

Note that preparation-level stereo diagnostics are also recorded in the `Prep_*` columns (e.g., `Prep_StereoStatus`). The `Stereo_*` columns reflect what the calculator actually did with stereochemistry for the *calculation structure*.

## Protomer (protonation) metadata columns (Task 2B)

When Dimorphite-DL protomer enumeration is enabled (`--ionization dimorphite`), the calculators will add `Protomer_*` columns describing the enumeration and the representative protomer selected. `Protomer_Selected_SMILES` stores the selected protomer. `Protomer_Count` stores how many unique protomers were generated (bounded by the configured maximum). `Protomer_Selected_Charge` stores the RDKit formal charge of the selected protomer.

By default, the calculators only record the representative protomer. If you pass `--calc-on-protomer`, the calculators will compute descriptors on the selected protomer and `Calc_Canonical_SMILES` will reflect that protomer.

## Tautomer metadata columns (Task 3)

When tautomer handling is enabled, the calculators may add `Tautomer_*` columns that describe how tautomer canonicalization/enumeration was handled. In `prep-canonical` mode, RDKit canonicalizes a single tautomer during preparation. In `enumerate` mode, RDKit enumerates a set of possible tautomers (up to a maximum) and the calculator selects a representative tautomer (currently the RDKit canonical tautomer) while also reporting the size of the tautomer set.

`Tautomer_Count` is the number of enumerated tautomers (bounded by the max). `Tautomer_Selected_SMILES` is the selected representative tautomer. `Tautomer_TopK_SMILES` is a semicolon-separated list of up to K tautomer SMILES for quick inspection.

To avoid ambiguity, `Canonical_SMILES` refers to the prepared molecule after standardization, while `Calc_Canonical_SMILES` refers to the exact structure used for descriptor calculations (which may be a selected tautomer or, in future, a selected protomer).

## Prefix based groups

Several property families are produced as a set of columns with shared prefixes.

`Tox_` columns represent toxicity structural alert counts or flags. `Tox_Risk` and `Tox_Total_Alerts` are summary columns. The toxicity analyzer will include every column that starts with `Tox_` so that if you add new alerts you do not need to change the analyzer.

`Met_` columns represent metabolic soft spot counts, conjugation flags, or summary metrics. `Met_Stability` and `Met_Total_Sites` are summary columns.

`hERG_` columns represent hERG risk features and risk classification.

`CYP` columns represent cytochrome P450 inhibition risk flags and summary metrics.

`Sol_` columns represent solubility predictions and classes.

`Perm_` columns represent permeability predictions, classes, and transporter likelihood heuristics.

`Lead_` columns represent lead optimization and efficiency metrics.

`SA_` columns represent synthetic accessibility summary outputs. `SA_Score` is a 1–10 score (lower is easier). `SA_Class` is a coarse label.

`Complexity_` columns represent a small bundle of interpretable complexity metrics (graph indices, stereocenter counts, spiro/bridgeheads, aromaticity fraction).

`Scaffold_` columns represent Murcko scaffold assignments used for series analysis. `Scaffold_ID` is a stable identifier assigned in descending frequency order within the analyzed dataset; `Scaffold_Size` is the number of compounds sharing that scaffold.

`Cluster_` columns represent similarity clustering assignments (Butina clustering using Morgan fingerprints). `Cluster_ID` is the cluster assignment and `Cluster_Size` is the number of compounds in the cluster.

`3D_` columns represent optional 3D conformer-based shape descriptors computed by `calculators/mpo_v5.py` when run with `--3d`. They are based on RDKit ETKDG conformer generation and optional MMFF/UFF minimization, so treat them as qualitative ranking features rather than ground truth physical properties.

`Dev_` columns represent developability-oriented indices and rule-style filters (Task 8), including the GSK 4/400 and Pfizer 3/75 checks plus an aggregated risk burden and score for triage.

## Short vs long property names

The calculators usually provide a short set of convenience columns such as `MolWt`, `LogP`, `TPSA`, `HBD`, and `HBA`, and also provide longer, namespaced versions such as `RDKit_MolWt` for the full RDKit descriptor table. This is helpful because the short names are easy to read in the early columns, while the namespaced versions avoid collisions and are convenient for machine learning.

## Missing values

Some descriptor values are undefined when the corresponding structural feature does not exist in a molecule. In raw outputs, these cells are empty. In the calculators, `--fill-missing 0` can be used to replace empties with a constant, which is often necessary for model training. When you do this, it is good practice to track which descriptors are truly “not applicable” rather than “measured as zero,” because those are different concepts.


## Schema metadata (Task 11)

If you need a machine-readable manifest of column types, units, and short guidance, see `docs/schema.json`. You can validate a results table with `tools/validate_csv_schema.py` and you can generate a schema-derived column reference doc (`docs/columns_schema.md`) with `tools/generate_columns_md_from_schema.py`.

