# Column reference (generated)

This file is generated from `docs/schema.json`. Do not edit by hand; update the schema and re-run the generator.

Generated: `2026-01-24T08:46:09`

## 3d

| column | type | units | required | description |
| --- | --- | --- | --- | --- |
| `3D_Asphericity` | `float` | `` | no | RDKit 3D Asphericity descriptor. |
| `3D_BestEnergy` | `float` | `kcal/mol` | no | Force-field energy of the selected conformer (RDKit FF units; qualitative). |
| `3D_Converged` | `bool` | `` | no | Whether the selected conformer minimization converged (best-effort). |
| `3D_Eccentricity` | `float` | `` | no | RDKit 3D Eccentricity descriptor. |
| `3D_InertialShapeFactor` | `float` | `` | no | RDKit 3D InertialShapeFactor descriptor. |
| `3D_Minimizer` | `string` | `` | no | Minimizer used for conformers (mmff/uff/none). |
| `3D_Mode` | `string` | `` | no | 3D descriptor mode/status (off/etkdg/fail). |
| `3D_NPR1` | `float` | `` | no | RDKit 3D NPR1 (normalized principal moment ratio 1). |
| `3D_NPR2` | `float` | `` | no | RDKit 3D NPR2 (normalized principal moment ratio 2). |
| `3D_NumConfs` | `int` | `count` | no | Number of conformers embedded for 3D descriptors. |
| `3D_PBF` | `float` | `` | no | RDKit 3D Plane of Best Fit (PBF) descriptor (lower is more planar). |
| `3D_PMI1` | `float` | `` | no | RDKit 3D PMI1 (principal moment of inertia 1). |
| `3D_PMI2` | `float` | `` | no | RDKit 3D PMI2 (principal moment of inertia 2). |
| `3D_PMI3` | `float` | `` | no | RDKit 3D PMI3 (principal moment of inertia 3). |
| `3D_RadiusOfGyration` | `float` | `angstrom` | no | RDKit 3D RadiusOfGyration descriptor. |
| `3D_SelectedConfId` | `int` | `` | no | Conformer ID selected for descriptor calculation (lowest energy when minimized). |
| `3D_SpherocityIndex` | `float` | `` | no | RDKit 3D SpherocityIndex descriptor. |

## calc

| column | type | units | required | description |
| --- | --- | --- | --- | --- |
| `Calc_Base_SMILES` | `string` | `` | no | Base structure SMILES used as the starting point for calculation structure selection (e.g., after prep + tautomer + stereo handling, before optional protomer selection). |
| `Calc_Canonical_SMILES` | `string` | `` | no | Canonical SMILES for the exact structure used to compute descriptors. If calc-on-protomer is enabled, this may reflect a selected protomer. |

## cns_mpo

| column | type | units | required | description |
| --- | --- | --- | --- | --- |
| `BBB_Likely` | `bool` | `` | no | Heuristic BBB penetration likelihood flag derived from MPO and property bounds. |
| `CNS_MPO` | `float` | `score` | no | CNS MPO total score (0–6). |
| `CNS_MPO_Class` | `string` | `` | no | CNS MPO qualitative class (Excellent/Good/Moderate/Poor). |
| `MPO_HBD` | `float` | `desirability` | no | CNS MPO desirability subscore for HBD. |
| `MPO_LogD` | `float` | `desirability` | no | CNS MPO desirability subscore for logD. |
| `MPO_LogP` | `float` | `desirability` | no | CNS MPO desirability subscore for LogP. |
| `MPO_MW` | `float` | `desirability` | no | CNS MPO desirability subscore for MW. |
| `MPO_TPSA` | `float` | `desirability` | no | CNS MPO desirability subscore for TPSA. |
| `MPO_pKa` | `float` | `desirability` | no | CNS MPO desirability subscore for basic pKa. |

## core

| column | type | units | required | description |
| --- | --- | --- | --- | --- |
| `AromaticRings` | `int` | `count` | no | Aromatic ring count. |
| `FractionCSP3` | `float` | `fraction` | no | Fraction sp3 carbon (0–1). |
| `HBA` | `int` | `count` | no | H-bond acceptor count. |
| `HBD` | `int` | `count` | no | H-bond donor count. |
| `HeavyAtoms` | `int` | `count` | no | Heavy atom count. |
| `LogP` | `float` | `log10` | no | Crippen logP estimate (unitless log10). |
| `MolWt` | `float` | `Da` | no | Molecular weight (Daltons) of the calculation structure. |
| `Rings` | `int` | `count` | no | Ring count. |
| `RotBonds` | `int` | `count` | no | Rotatable bond count. |
| `TPSA` | `float` | `A^2` | no | Topological polar surface area (Å^2). |

## cyp

| column | type | units | required | description |
| --- | --- | --- | --- | --- |
| `CYP2C9_Inhibitor` | `bool` | `` | no | CYP2C9 inhibition risk flag (heuristic). |
| `CYP2D6_Inhibitor` | `bool` | `` | no | CYP2D6 inhibition risk flag (heuristic). |
| `CYP3A4_Inhibitor` | `bool` | `` | no | CYP3A4 inhibition risk flag (heuristic). |
| `CYP_MBI_Risk` | `bool` | `` | no | Mechanism-based inactivation motif risk flag. |
| `CYP_Total_Inhibition` | `int` | `count` | no | Total CYP inhibition burden (sum of flags). |

## developability

| column | type | units | required | description |
| --- | --- | --- | --- | --- |
| `Dev_BasePoints` | `int` | `points` | no | Base points used to build Dev_Score (counts property-space passes and other positive signals). |
| `Dev_Class` | `string` | `` | no | Developability class derived from Dev_Score and Dev_Risk_Class (Low/Medium/High). |
| `Dev_Flags` | `string` | `` | no | Comma-separated short flags explaining key contributors to risk burden. |
| `Dev_GSK_4_400_Pass` | `bool` | `` | no | GSK 4/400 rule: MW <= 400 and LogP <= 4. |
| `Dev_GoldenTriangle_Pass` | `bool` | `` | no | Coarse golden triangle check: 200 <= MW <= 500 and 1 <= cLogD_7.4 <= 3. |
| `Dev_Pfizer_3_75_Pass` | `bool` | `` | no | Pfizer 3/75 rule: LogP <= 3 and TPSA <= 75. |
| `Dev_Risk_Burden` | `int` | `points` | no | Aggregated risk burden points from available red flags (tox/hERG/etc.). |
| `Dev_Risk_Class` | `string` | `` | no | Risk class derived from Dev_Risk_Burden (Low/Medium/High). |
| `Dev_Score` | `int` | `points` | no | Composite developability score (heuristic): max(0, Dev_BasePoints - Dev_Risk_Burden). |
| `Dev_Solubility_Risk` | `string` | `` | no | Solubility risk class derived from Sol_ESOL_LogS when available (Low/Medium/High). |

## herg

| column | type | units | required | description |
| --- | --- | --- | --- | --- |
| `hERG_Alerts` | `int` | `count` | no | hERG alert feature count. |
| `hERG_Risk` | `string` | `` | no | hERG risk class (Low/Medium/High) from heuristic features. |

## identifier

| column | type | units | required | description |
| --- | --- | --- | --- | --- |
| `Compound_ID` | `string` | `` | yes | Compound identifier used for display and joining. Preferred ID column name for all analyzers/tools. |
| `Name` | `string` | `` | no | Alternative compound identifier sometimes used by upstream pipelines. |

## ionization

| column | type | units | required | description |
| --- | --- | --- | --- | --- |
| `Ion_FracDeprotonatedAcid` | `float` | `fraction` | no | Fraction deprotonated for an acid site at Ion_pH. |
| `Ion_FracProtonatedBase` | `float` | `fraction` | no | Fraction protonated for a base site at Ion_pH. |
| `Ion_LogD_pH_Est` | `float` | `log10` | no | Estimated logD at Ion_pH. |
| `Ion_Mode` | `string` | `` | no | Ionization feature mode (none/heuristic/dimorphite). |
| `Ion_NetCharge_Est` | `float` | `e` | no | Estimated net charge at Ion_pH (single-site approximation). |
| `Ion_NetCharge_Round` | `int` | `e` | no | Rounded net charge estimate at Ion_pH. |
| `Ion_State` | `string` | `` | no | Dominant ionization class at Ion_pH (neutral/cationic/anionic/zwitterionic). |
| `Ion_pH` | `float` | `pH` | no | pH used for heuristic ionization features. |
| `cLogD_7.4` | `float` | `log10` | no | Approximate logD at pH 7.4 derived from LogP and estimated pKa. |

## lead

| column | type | units | required | description |
| --- | --- | --- | --- | --- |
| `Lead_BEI` | `float` | `pIC50/(MW/1000)` | no | Binding efficiency index estimate. |
| `Lead_LE` | `float` | `pIC50/HA` | no | Ligand efficiency estimate (requires meaningful pIC50). |
| `Lead_LELP` | `float` | `logP/LE` | no | LELP estimate. |
| `Lead_LELP_Quality` | `string` | `` | no | LELP quality label. |
| `Lead_LE_Quality` | `string` | `` | no | LE quality label. |
| `Lead_LLE` | `float` | `pIC50-logP` | no | Lipophilic ligand efficiency estimate. |
| `Lead_LLE_Quality` | `string` | `` | no | LLE quality label. |
| `Lead_SEI` | `float` | `pIC50/(TPSA/100)` | no | Surface efficiency index estimate. |

## medchem_flags

| column | type | units | required | description |
| --- | --- | --- | --- | --- |
| `BRENK` | `bool` | `` | no | BRENK alert flag (RDKit FilterCatalog). |
| `MedChem_Alerts` | `int` | `count` | no | Total medchem alert count (PAINS/BRENK/NIH/ZINC). |
| `NIH` | `bool` | `` | no | NIH alert flag (RDKit FilterCatalog). |
| `PAINS` | `bool` | `` | no | PAINS alert flag (RDKit FilterCatalog). |
| `PAINS_Type` | `string` | `` | no | PAINS alert description if flagged. |
| `ZINC` | `bool` | `` | no | ZINC alert flag (RDKit FilterCatalog). |

## metabolism

| column | type | units | required | description |
| --- | --- | --- | --- | --- |
| `Met_Acetylation` | `bool` | `` | no | Phase II acetylation handle heuristic. |
| `Met_Glucuronidation` | `bool` | `` | no | Phase II glucuronidation handle heuristic. |
| `Met_Methylation` | `bool` | `` | no | Phase II methylation handle heuristic. |
| `Met_Stability` | `string` | `` | no | Metabolic stability class (High/Medium/Low) from soft-spot burden heuristic. |
| `Met_Sulfation` | `bool` | `` | no | Phase II sulfation handle heuristic. |
| `Met_Total_Sites` | `int` | `count` | no | Total metabolic soft-spot count. |

## oral_bioavailability

| column | type | units | required | description |
| --- | --- | --- | --- | --- |
| `Bioavail_Score` | `int` | `count` | no | Count of passed oral filters (0–4). |
| `Egan_Pass` | `bool` | `` | no | Egan oral filter pass flag. |
| `Ghose_Pass` | `bool` | `` | no | Ghose filter pass flag. |
| `Muegge_Pass` | `bool` | `` | no | Muegge filter pass flag. |
| `Veber_Pass` | `bool` | `` | no | Veber oral filter pass flag. |

## permeability

| column | type | units | required | description |
| --- | --- | --- | --- | --- |
| `Perm_Caco2_Class` | `string` | `` | no | Caco-2 permeability class label. |
| `Perm_Caco2_LogPapp` | `float` | `log10(Papp)` | no | Heuristic Caco-2 log(Papp) estimate. |
| `Perm_GI_Absorption` | `string` | `` | no | GI absorption class label. |
| `Perm_PAMPA_Class` | `string` | `` | no | PAMPA permeability class label. |
| `Perm_PAMPA_LogPe` | `float` | `log10(Pe)` | no | Heuristic PAMPA log(Pe) estimate. |
| `Perm_Pgp_Likelihood` | `string` | `` | no | P-gp substrate likelihood class (Low/Medium/High). |
| `Perm_Skin_LogKp` | `float` | `log10(cm/s)` | no | Heuristic skin permeability logKp estimate. |

## pk

| column | type | units | required | description |
| --- | --- | --- | --- | --- |
| `CL_Risk` | `string` | `` | no | Clearance risk class label. |
| `CL_SoftSpots` | `int` | `count` | no | Clearance soft-spot burden (copied from Met_Total_Sites). |
| `F_Class` | `string` | `` | no | Oral bioavailability class label. |
| `F_Est_%` | `string` | `%` | no | Oral F estimate band (text, e.g., '>70'). |
| `F_Score` | `int` | `count` | no | Oral F scoring heuristic (0–4). |
| `PPB_Class` | `string` | `` | no | Protein binding class label. |
| `PPB_Est_%` | `float` | `%` | no | Estimated percent protein bound (heuristic). |
| `PPB_FreeEst_%` | `float` | `%` | no | Estimated percent unbound (100 - PPB_Est_%). |
| `T_half_Est` | `string` | `` | no | Half-life class estimate. |
| `Vd_Class` | `string` | `` | no | Vd class label. |
| `Vd_Est_L_kg` | `float` | `L/kg` | no | Heuristic volume of distribution estimate. |

## pka

| column | type | units | required | description |
| --- | --- | --- | --- | --- |
| `pKa_Acidic` | `float` | `pKa` | no | Estimated acidic pKa (heuristic/fragment-based). |
| `pKa_Basic` | `float` | `pKa` | no | Estimated basic pKa (heuristic/fragment-based). |
| `pKa_Type` | `string` | `` | no | Dominant pKa type label (basic/acidic/neutral) used for simplified modeling. |

## prep

| column | type | units | required | description |
| --- | --- | --- | --- | --- |
| `Canonical_SMILES` | `string` | `` | no | RDKit-canonicalized SMILES after preparation (standardization/parent selection/optional uncharge+reionize/tautomer handling). |
| `Input_Canonical_SMILES` | `string` | `` | no | RDKit-canonicalized SMILES as parsed from the input before preparation. |
| `Prep_Changed` | `bool` | `` | no | Whether preparation modified the structure. |
| `Prep_Notes` | `string` | `` | no | Free-text notes from the preparation step. |
| `Prep_Reionized` | `bool` | `` | no | Whether the structure was reionized after neutralization. |
| `Prep_SaltsStripped` | `bool` | `` | no | Whether salts/fragments were stripped during parent selection. |
| `Prep_StereoStatus` | `string` | `` | no | Stereo completeness diagnostic from preparation (e.g., specified vs unassigned). |
| `Prep_TautomerChanged` | `bool` | `` | no | Whether preparation changed the tautomer form. |
| `Prep_UnassignedChiral` | `int` | `count` | no | Number of unassigned chiral centers detected during preparation. |
| `Prep_UnassignedStereoBonds` | `int` | `count` | no | Number of unassigned stereogenic double bonds detected during preparation. |
| `Prep_Uncharged` | `bool` | `` | no | Whether the structure was neutralized during preparation. |

## protomer

| column | type | units | required | description |
| --- | --- | --- | --- | --- |
| `Protomer_CalcOnProtomer` | `bool` | `` | no | Whether descriptors were computed on the selected protomer. |
| `Protomer_Count` | `int` | `count` | no | Number of protomers enumerated (bounded by config). |
| `Protomer_MaxVariants` | `int` | `count` | no | Maximum protomer variants retained. |
| `Protomer_Mode` | `string` | `` | no | Protomer handling mode (none/heuristic/dimorphite). |
| `Protomer_Notes` | `string` | `` | no | Free-text notes from protomer enumeration/selection. |
| `Protomer_Precision` | `float` | `pH` | no | Dimorphite precision parameter controlling microstate generation. |
| `Protomer_Selected_Charge` | `int` | `e` | no | Formal charge of the selected protomer (RDKit). |
| `Protomer_Selected_SMILES` | `string` | `` | no | Representative protomer SMILES selected. |
| `Protomer_Target_Charge` | `int` | `e` | no | Target charge used by selection logic (if applicable). |
| `Protomer_TopK_SMILES` | `string` | `` | no | Semicolon-separated list of up to K protomer SMILES. |
| `Protomer_pH_Max` | `float` | `pH` | no | Maximum pH used for protomer enumeration window. |
| `Protomer_pH_Min` | `float` | `pH` | no | Minimum pH used for protomer enumeration window. |

## qed

| column | type | units | required | description |
| --- | --- | --- | --- | --- |
| `QED` | `float` | `score` | no | Quantitative Estimate of Drug-likeness (0–1). |

## rule_of_3

| column | type | units | required | description |
| --- | --- | --- | --- | --- |
| `RO3_HBA` | `bool` | `` | no | RO3 HBA criterion (HBA<=3). |
| `RO3_HBD` | `bool` | `` | no | RO3 HBD criterion (HBD<=3). |
| `RO3_LogP` | `bool` | `` | no | RO3 LogP criterion (LogP<=3). |
| `RO3_MW` | `bool` | `` | no | RO3 MW criterion (MW<=300). |
| `RO3_Pass` | `bool` | `` | no | Rule-of-3 pass flag (0 violations). |
| `RO3_RotBonds` | `bool` | `` | no | RO3 rotatable bonds criterion (RotBonds<=3). |
| `RO3_TPSA` | `bool` | `` | no | RO3 TPSA criterion (TPSA<=60). |
| `RO3_Violations` | `int` | `count` | no | Number of RO3 violations. |

## rule_of_5

| column | type | units | required | description |
| --- | --- | --- | --- | --- |
| `RO5_HBA` | `bool` | `` | no | RO5 HBA criterion (HBA<=10). |
| `RO5_HBD` | `bool` | `` | no | RO5 HBD criterion (HBD<=5). |
| `RO5_LogP` | `bool` | `` | no | RO5 LogP criterion (LogP<=5). |
| `RO5_MW` | `bool` | `` | no | RO5 MW criterion (MW<=500). |
| `RO5_Pass` | `bool` | `` | no | Lipinski Rule-of-5 pass flag (≤1 violation). |
| `RO5_Violations` | `int` | `count` | no | Number of RO5 violations. |

## sa_complexity

| column | type | units | required | description |
| --- | --- | --- | --- | --- |
| `Complexity_AromaticAtomFrac` | `float` | `fraction` | no | Fraction of heavy atoms that are aromatic (0–1). |
| `Complexity_BalabanJ` | `float` | `index` | no | Balaban J index. |
| `Complexity_BertzCT` | `float` | `index` | no | Bertz complexity index (graph complexity). |
| `Complexity_HallKierAlpha` | `float` | `index` | no | Hall–Kier alpha value. |
| `Complexity_NumBridgeheadAtoms` | `int` | `count` | no | Number of bridgehead atoms. |
| `Complexity_NumHeteroAtoms` | `int` | `count` | no | Number of heteroatoms. |
| `Complexity_NumSpiroAtoms` | `int` | `count` | no | Number of spiro atoms. |
| `Complexity_NumStereoCenters` | `int` | `count` | no | Number of atom stereocenters (specified). |
| `Complexity_NumUnspecifiedStereoCenters` | `int` | `count` | no | Number of unspecified atom stereocenters. |
| `SA_Class` | `string` | `` | no | Coarse SA score class (Easy/Moderate/Hard). |
| `SA_Score` | `float` | `score` | no | Synthetic accessibility score (Ertl–Schuffenhauer), scaled 1–10 where lower is easier. |

## series

| column | type | units | required | description |
| --- | --- | --- | --- | --- |
| `Cluster_ID` | `int` | `` | no | Similarity cluster assignment (Butina clustering on Morgan fingerprints). |
| `Cluster_MedoidMeanSim` | `float` | `tanimoto` | no | Mean Tanimoto similarity of the medoid to other cluster members. |
| `Cluster_Representative_ID` | `string` | `` | no | Compound_ID of the cluster medoid (representative). |
| `Cluster_Size` | `int` | `count` | no | Number of compounds in the similarity cluster. |
| `Scaffold_ID` | `string` | `` | no | Murcko scaffold identifier assigned in descending frequency order within the analyzed dataset. |
| `Scaffold_SMILES` | `string` | `` | no | Canonical Murcko scaffold SMILES. |
| `Scaffold_Size` | `int` | `count` | no | Number of compounds sharing the scaffold. |

## solubility

| column | type | units | required | description |
| --- | --- | --- | --- | --- |
| `Sol_ESOL_Class` | `string` | `` | no | ESOL solubility class label. |
| `Sol_ESOL_LogS` | `float` | `log10(mol/L)` | no | ESOL predicted logS. |
| `Sol_ESOL_mg_L` | `float` | `mg/L` | no | ESOL predicted solubility converted to mg/L. |

## stereo

| column | type | units | required | description |
| --- | --- | --- | --- | --- |
| `Stereo_Changed` | `bool` | `` | no | Whether stereochemistry handling changed the calculation structure. |
| `Stereo_IsomerCount` | `int` | `count` | no | Number of stereoisomers considered in enumerate mode (bounded by config). |
| `Stereo_Mode` | `string` | `` | no | Stereochemistry handling mode (keep/strip/enumerate). |
| `Stereo_Notes` | `string` | `` | no | Free-text notes from stereochemistry handling. |
| `Stereo_Selected_SMILES` | `string` | `` | no | Representative stereoisomer SMILES selected for calculations. |
| `Stereo_TopK_SMILES` | `string` | `` | no | Semicolon-separated list of up to K stereoisomer SMILES for quick inspection. |

## structure

| column | type | units | required | description |
| --- | --- | --- | --- | --- |
| `SMILES` | `string` | `` | no | Raw input SMILES if provided in the upstream table. |

## tautomer

| column | type | units | required | description |
| --- | --- | --- | --- | --- |
| `Tautomer_Changed` | `bool` | `` | no | Whether tautomer handling changed the calculation structure. |
| `Tautomer_Count` | `int` | `count` | no | Number of tautomers enumerated (bounded by config). |
| `Tautomer_Mode` | `string` | `` | no | Tautomer handling mode (none/prep-canonical/enumerate). |
| `Tautomer_Selected_SMILES` | `string` | `` | no | Representative tautomer SMILES selected. |
| `Tautomer_TopK_SMILES` | `string` | `` | no | Semicolon-separated list of up to K tautomer SMILES. |

## toxicity

| column | type | units | required | description |
| --- | --- | --- | --- | --- |
| `Tox_Risk` | `string` | `` | no | Toxicity risk class (Low/Medium/High) derived from alert families. |
| `Tox_Total_Alerts` | `int` | `count` | no | Total toxicity alert family count. |

## Pattern-based column families

Patterns define metadata for whole families of columns, typically produced with shared prefixes.

| regex | type | tags | description |
| --- | --- | --- | --- |
| `^Tox_` | `auto` | toxicity | Toxicity structural alert family columns. These may be boolean flags or alert counts depending on implementation. |
| `^Met_` | `auto` | metabolism | Metabolism/soft-spot columns (counts, boolean handles, or summary fields). |
| `^hERG_` | `auto` | herg | hERG risk feature columns. |
| `^CYP` | `auto` | cyp | CYP inhibition risk columns. |
| `^Sol_` | `auto` | solubility | Solubility prediction columns. |
| `^Perm_` | `auto` | permeability | Permeability prediction columns. |
| `^Lead_` | `auto` | lead | Lead optimization metric columns. |
| `^SA_` | `auto` | sa_complexity | Synthetic accessibility family columns. |
| `^Complexity_` | `auto` | sa_complexity | Complexity descriptor family columns. |
| `^Scaffold_` | `auto` | series | Series analytics scaffold columns (Murcko scaffolds). |
| `^Cluster_` | `auto` | series | Series analytics similarity clustering columns. |
| `^3D_` | `auto` | 3d | 3D conformer-based descriptor columns (v5 calculator optional). |
| `^Dev_` | `auto` | developability | Developability indices and aggregated risk/score columns. |
| `^RDKit_` | `float` | descriptors | Namespaced RDKit descriptor columns (long-form descriptor table). |
| `^Mordred_` | `float` | descriptors | Namespaced Mordred descriptor columns (optional family, if enabled). |
