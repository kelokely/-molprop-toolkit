# Feature-to-columns (column families)

This reference groups related output fields by their shared column prefixes (or common naming patterns) and gives a one-line interpretation for how each family is typically used during triage.

| Column family (pattern) | Example columns | One-line interpretation (how to use it) |
|---|---|---|
| `Prep_*` | `Prep_Changed`, `Prep_SaltsStripped`, `Prep_Uncharged`, `Prep_Notes` | Audit trail of structure standardization; use to explain unexpected property shifts vs the input SMILES. |
| `Tautomer_*` | `Tautomer_Mode`, `Tautomer_Count`, `Tautomer_Selected_SMILES`, `Tautomer_Changed` | Records tautomer handling so you can judge whether tautomer choice might be affecting ranking/SAR interpretation. |
| `Stereo_*` | `Stereo_Mode`, `Stereo_IsomerCount`, `Stereo_Selected_SMILES`, `Stereo_Notes` | Captures stereochemistry handling; use to flag ambiguity/enumeration and interpret variability driven by stereo selection. |
| `Protomer_*` | `Protomer_Mode`, `Protomer_Count`, `Protomer_Selected_SMILES`, `Protomer_Selected_Charge` | Summarizes protomer enumeration/selection; use to understand charge-state sensitivity and pH-driven ADME risks. |
| `Ion_*` | `Ion_State`, `Ion_NetCharge_Est`, `Ion_LogD_pH_Est`, `Ion_FracProtonatedBase` | Heuristic ionization at a chosen pH; use as qualitative guidance for permeability/BBB/PPB risk and LogD context. |
| `MPO_*` | `MPO_MW`, `MPO_LogP`, `MPO_LogD`, `MPO_TPSA`, `MPO_HBD`, `MPO_pKa` | Component-level breakdown of CNS MPO; use to diagnose why CNS_MPO is high/low and what to optimize. |
| `RO5_*` | `RO5_Pass`, `RO5_Violations`, `RO5_MW`, `RO5_LogP` | Lipinski Rule-of-5 compliance; use violations as risk indicators (not absolute gates) for oral developability. |
| `RO3_*` | `RO3_Pass`, `RO3_Violations`, `RO3_TPSA`, `RO3_RotBonds` | Fragment-likeness heuristic; use when triaging fragment libraries or checking lead-like vs fragment-like space. |
| `Tox_*` | `Tox_Risk`, `Tox_Total_Alerts`, `Tox_michael_acceptor`, `Tox_nitro_aromatic` | Structural alert counts and a coarse risk class; use to prioritize follow-up safety/reactivity assessment early. |
| `Met_*` | `Met_Total_Sites`, `Met_Stability`, `Met_Glucuronidation`, `Met_Sulfation` | Heuristic metabolic soft-spot and conjugation signals; use for clearance-risk ranking and liability flagging. |
| `hERG_*` | `hERG_Risk`, `hERG_Alerts`, `hERG_basic_nitrogen`, `hERG_tertiary_amine` | Coarse hERG liability heuristics; use to flag compounds for earlier electrophysiology counterscreens. |
| `CYP_*` | `CYP3A4_Inhibitor`, `CYP2D6_Inhibitor`, `CYP_MBI_Risk`, `CYP_Total_Inhibition` | DDI risk heuristics; use total burden plus MBI motifs to prioritize cleaner profiles and avoid late surprises. |
| `Sol_*` | `Sol_ESOL_LogS`, `Sol_ESOL_mg_L`, `Sol_ESOL_Class` | Predicted solubility signals; best used for within-series ranking and anticipating formulation risk. |
| `Perm_*` | `Perm_Caco2_Class`, `Perm_PAMPA_Class`, `Perm_GI_Absorption`, `Perm_Pgp_Likelihood` | Qualitative permeability/absorption/efflux triage; use to spot polarity/efflux-driven exposure and BBB issues. |
| `PPB_*` | `PPB_Est_%`, `PPB_FreeEst_%`, `PPB_Class` | Heuristic protein binding; use as a quick exposure headwind signal (interpret qualitatively). |
| `Vd_*` | `Vd_Est_L_kg`, `Vd_Class` | Rough distribution tendency (plasma vs tissue); use as qualitative PK context, not a predictor. |
| `CL_*` | `CL_Risk`, `CL_SoftSpots` | Clearance-risk heuristic driven by lipophilicity and soft spots; use to prioritize metabolic stability work. |
| `F_*` | `F_Class`, `F_Est_%`, `F_Score` | Coarse oral bioavailability class; use as a triage hint alongside solubility/permeability flags. |
| `Lead_*` | `Lead_LE`, `Lead_LLE`, `Lead_LELP`, `Lead_LLE_Quality` | Potency-efficiency metrics (meaningful with real potency); use to prevent size/lipophilicity creep during lead optimization. |
| `SA_*` | `SA_Score`, `SA_Class` | Synthetic accessibility estimate; use to avoid drifting into hard-to-make chemical space unintentionally. |
| `Complexity_*` | `Complexity_BertzCT`, `Complexity_NumStereoCenters`, `Complexity_NumSpiroAtoms` | Structural complexity descriptors; use for risk-of-synthesis awareness and cross-series comparability. |
| `3D_*` | `3D_Mode`, `3D_NPR1`, `3D_NPR2`, `3D_RadiusOfGyration`, `3D_PBF` | Conformer-derived 3D shape descriptors (v5 `--3d`); use for qualitative shape diversity/planarity/3D-ness ranking. |
| `Dev_*` | `Dev_GSK_4_400_Pass`, `Dev_Pfizer_3_75_Pass`, `Dev_Risk_Burden`, `Dev_Score`, `Dev_Class` | Developability-oriented rule checks plus a heuristic score/burden; use for operational triage and picklist filtering. |
| `Scaffold_*` | `Scaffold_ID`, `Scaffold_SMILES`, `Scaffold_Size` | Murcko scaffold grouping; use for series-level summaries and to spot scaffold-specific outliers/opportunities. |
| `Cluster_*` | `Cluster_ID`, `Cluster_Size`, `Cluster_Representative_ID`, `Cluster_MedoidMeanSim` | Similarity clustering; use to understand diversity, redundancy, and analog-series structure in the dataset. |

