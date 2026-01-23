"""Category registry (Task 10C).

The registry is the canonical source of category metadata used by:
- analyzer_master (interactive CLI)
- future report builders (Task 12)
- future schema/manifest generation (Task 11)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence


@dataclass(frozen=True)
class CategorySpec:
    key: str
    name: str
    description: str
    columns: Sequence[str]
    key_column: Optional[str] = None
    prefix: Optional[str] = None


# NOTE: This was originally embedded in analyzers/analyzer_master.py.
# It is moved here so the rest of the toolkit can share a single registry.
CATEGORY_SPECS: Dict[str, CategorySpec] = {
    "cns_mpo": CategorySpec(
        key="cns_mpo",
        name="CNS MPO Analysis",
        description="""CNS Multiparameter Optimization (MPO) Score

The CNS MPO score was developed to capture a region of property space where
central nervous system (CNS) drugs tend to succeed. The score combines six
components, each mapped to a desirability score, and summed to a total score.

The components are molecular weight, LogP, LogD at pH 7.4, TPSA, hydrogen-bond
Donors, and the most basic pKa. The total score ranges from 0 to 6. In many
workflows, values around 4 and above are treated as a practical threshold for
acceptable CNS balance, and values around 5 and above are considered strong.

The utility of CNS MPO is not that it replaces experimental BBB assays, but
that it helps teams see why a compound is failing. If the total is low, the
component breakdown often tells you whether the issue is excessive polarity,
excessive basicity, excessive lipophilicity, or size.
""",
        columns=(
            "CNS_MPO",
            "CNS_MPO_Class",
            "BBB_Likely",
            "MPO_MW",
            "MPO_LogP",
            "MPO_LogD",
            "MPO_TPSA",
            "MPO_HBD",
            "MPO_pKa",
            "cLogD_7.4",
            "pKa_Basic",
            "pKa_Acidic",
            "MolWt",
            "LogP",
            "TPSA",
            "HBD",
        ),
        key_column="CNS_MPO",
    ),
    "oral_bioavailability": CategorySpec(
        key="oral_bioavailability",
        name="Oral Bioavailability Analysis",
        description="""Oral bioavailability filters

Oral bioavailability depends on solubility, permeability, stability, and
first-pass metabolism. Because these factors are complex, multiple filters are
commonly used together. Veber emphasizes TPSA and rotatable bonds, Egan uses a
TPSA–LogP boundary, and Ghose and Muegge are broader property filters derived
from distributions observed in known drugs.

The calculators often include a simple Bioavailability Score computed as the
count of passed filters. This is primarily a ranking tool rather than an
absolute predictor.
""",
        columns=(
            "Bioavail_Score",
            "Veber_Pass",
            "Egan_Pass",
            "Ghose_Pass",
            "Muegge_Pass",
            "RO5_Pass",
            "TPSA",
            "RotBonds",
            "MolWt",
            "LogP",
        ),
        key_column="Bioavail_Score",
    ),
    "rule_of_5": CategorySpec(
        key="rule_of_5",
        name="Rule of 5 (Lipinski) Analysis",
        description="""Lipinski Rule of 5

The Rule of 5 is a classic heuristic for oral small-molecule drug likeness. It
considers molecular weight, LogP, hydrogen-bond donors, and hydrogen-bond
acceptors. The typical interpretation is that 0 violations is excellent, 1
violation is often acceptable, and 2 or more violations increase the risk of
poor absorption and/or solubility.
""",
        columns=(
            "RO5_Pass",
            "RO5_Violations",
            "RO5_MW",
            "RO5_LogP",
            "RO5_HBD",
            "RO5_HBA",
            "MolWt",
            "LogP",
            "HBD",
            "HBA",
        ),
        key_column="RO5_Pass",
    ),
    "rule_of_3": CategorySpec(
        key="rule_of_3",
        name="Rule of 3 (Fragment) Analysis",
        description="""Rule of 3

The Rule of 3 is associated with fragment-based drug discovery. It helps
identify small, simple compounds that often serve as efficient starting points
for optimization. Rule of 3 compliance is not required for a drug, but it is
useful for triaging fragment libraries or fragment-like hits.
""",
        columns=(
            "RO3_Pass",
            "RO3_Violations",
            "RO3_MW",
            "RO3_LogP",
            "RO3_HBD",
            "RO3_HBA",
            "RO3_RotBonds",
            "RO3_TPSA",
        ),
        key_column="RO3_Pass",
    ),
    "medchem_flags": CategorySpec(
        key="medchem_flags",
        name="MedChem Flags (PAINS/BRENK/NIH/ZINC)",
        description="""Med chem structural flags

These flags are substructure filters intended to highlight motifs that are
commonly problematic, either because they are frequent assay interferers, are
reactive, or historically correlate with instability or toxicity. A flag is a
risk indicator, not an automatic rejection, but it should prompt orthogonal
validation and careful interpretation.
""",
        columns=("PAINS", "PAINS_Type", "BRENK", "NIH", "ZINC", "MedChem_Alerts"),
        key_column="MedChem_Alerts",
    ),
    "toxicity": CategorySpec(
        key="toxicity",
        name="Toxicity Structural Alerts",
        description="""Toxicity alert summary

This toolkit primarily uses structural alerting rather than trained predictive
models. Alerts are intended to prompt investigation and prioritization of
experiments, not to claim clinical toxicity. The summary columns include a risk
class and a total alert count, while the per-alert columns are prefixed with
Tox_.
""",
        columns=("Tox_Risk", "Tox_Total_Alerts"),
        prefix="Tox_",
        key_column="Tox_Risk",
    ),
    "metabolism": CategorySpec(
        key="metabolism",
        name="Metabolism and Soft Spots",
        description="""Metabolism and soft spots

Metabolic soft spots are a practical way to anticipate clearance and reactive
metabolite risk. The columns prefixed with Met_ represent soft spot counts,
phase II conjugation flags, and summary outputs such as total sites and
stability class.
""",
        columns=(
            "Met_Stability",
            "Met_Total_Sites",
            "Met_Glucuronidation",
            "Met_Sulfation",
            "Met_Acetylation",
            "Met_Methylation",
        ),
        prefix="Met_",
        key_column="Met_Stability",
    ),
    "herg": CategorySpec(
        key="herg",
        name="hERG Risk",
        description="""hERG risk

hERG liability is a key cardiac safety risk. The toolkit provides heuristic
signals based on common structural and property correlations. Experimental
testing is required for decisions.
""",
        columns=(
            "hERG_Risk",
            "hERG_Alerts",
            "hERG_basic_nitrogen",
            "hERG_aromatic_ring",
            "hERG_hydrophobic_center",
            "hERG_tertiary_amine",
        ),
        key_column="hERG_Risk",
    ),
    "cyp": CategorySpec(
        key="cyp",
        name="CYP Inhibition Risk",
        description="""CYP inhibition risk

Cytochrome P450 inhibition can create drug–drug interaction risk. The toolkit
reports simple inhibitor flags for major isoforms and an additional flag for
mechanism-based inactivation motifs.
""",
        columns=(
            "CYP3A4_Inhibitor",
            "CYP2D6_Inhibitor",
            "CYP2C9_Inhibitor",
            "CYP_MBI_Risk",
            "CYP_Total_Inhibition",
        ),
        key_column="CYP_Total_Inhibition",
    ),
    "solubility": CategorySpec(
        key="solubility",
        name="Solubility",
        description="""Solubility

Predicted solubility metrics are best treated as ranking signals within a
series. The toolkit includes logS estimators such as ESOL and a class label.
""",
        columns=(
            "Sol_ESOL_LogS",
            "Sol_ESOL_mg_L",
            "Sol_ESOL_Class",
            "Sol_Ali_LogS",
            "Sol_SilicosIT_LogS",
            "Sol_Consensus_LogS",
        ),
        key_column="Sol_Consensus_LogS",
    ),
    "permeability": CategorySpec(
        key="permeability",
        name="Permeability",
        description="""Permeability

Permeability estimates and qualitative absorption flags help triage whether a
series is likely to be absorption limited. P-gp likelihood is included as a
separate indicator because it can cause cell efflux and reduced brain
penetration.
""",
        columns=(
            "Perm_GI_Absorption",
            "Perm_Caco2_LogPapp",
            "Perm_Caco2_Class",
            "Perm_PAMPA_LogPe",
            "Perm_PAMPA_Class",
            "Perm_Pgp_Likelihood",
            "Perm_Skin_LogKp",
        ),
        key_column="Perm_GI_Absorption",
    ),
    "pk": CategorySpec(
        key="pk",
        name="PK heuristics",
        description="""PK heuristics

The PK columns provide coarse summaries of protein binding, distribution,
clearance risk, half-life class, and an oral bioavailability class. These are
intended as triage signals.
""",
        columns=(
            "PPB_Est_%",
            "PPB_Class",
            "Vd_Est_L_kg",
            "Vd_Class",
            "CL_Risk",
            "CL_SoftSpots",
            "T_half_Est",
            "F_Class",
            "F_Est_%",
        ),
        key_column="F_Class",
    ),
    "lead": CategorySpec(
        key="lead",
        name="Lead optimization metrics",
        description="""Lead optimization metrics

Lead metrics such as LE, LLE, and LELP are most meaningful when experimental
potency is available. In this toolkit they may be computed using default
assumptions if potency is not provided. These metrics help avoid making
molecules larger and more lipophilic without gaining potency in an efficient
manner.
""",
        columns=(
            "Lead_LE",
            "Lead_LE_Quality",
            "Lead_LLE",
            "Lead_LLE_Quality",
            "Lead_LELP",
            "Lead_LELP_Quality",
            "Lead_BEI",
            "Lead_SEI",
            "Lead_FitQuality",
        ),
        key_column="Lead_LLE",
    ),
    "qed": CategorySpec(
        key="qed",
        name="QED",
        description="""QED

QED is a continuous drug-likeness score (0 to 1) derived from desirability
functions fit to distributions of known drugs. It is often used for ranking and
triage. Because it is continuous, it can be more informative than pass/fail
filters, but it should still be used alongside other project constraints.
""",
        columns=(
            "QED",
            "QED_MW",
            "QED_ALOGP",
            "QED_HBA",
            "QED_HBD",
            "QED_PSA",
            "QED_ROTB",
            "QED_AROM",
            "QED_ALERTS",
        ),
        key_column="QED",
    ),
    "sa_complexity": CategorySpec(
        key="sa_complexity",
        name="Synthetic accessibility + complexity",
        description="""Synthetic accessibility and complexity

This category summarizes the Ertl–Schuffenhauer SA score (1–10; lower is easier)
along with a small bundle of interpretable RDKit complexity descriptors.

SA_Score is best used as a ranking signal. The Complexity_* columns help
interpret why a structure may be challenging: spiro/bridgehead content,
stereocenters (including unspecified stereocenters), aromaticity fraction, and
graph complexity indices (e.g., BertzCT).
""",
        columns=(
            "SA_Score",
            "SA_Class",
            "Complexity_BertzCT",
            "Complexity_BalabanJ",
            "Complexity_HallKierAlpha",
            "Complexity_NumSpiroAtoms",
            "Complexity_NumBridgeheadAtoms",
            "Complexity_NumStereoCenters",
            "Complexity_NumUnspecifiedStereoCenters",
            "Complexity_NumHeteroAtoms",
            "Complexity_AromaticAtomFrac",
        ),
        prefix="Complexity_",
        key_column="SA_Score",
    ),
    "scaffolds": CategorySpec(
        key="scaffolds",
        name="Murcko scaffolds",
        description="""Murcko scaffolds

Murcko scaffolds collapse compounds to a core ring/linker framework so that
series can be summarized at the scaffold level. In practice this helps answer
questions like: "Are my top MPO compounds all on the same chemotype?" and
"Which scaffold families have the best balance of potency proxies vs risk?"

Run `python analyzers/analyzer_scaffolds.py results.csv` or the combined runner
`python tools/series_analytics.py results.csv` to generate these columns.
""",
        columns=(
            "Scaffold_ID",
            "Scaffold_SMILES",
            "Scaffold_Size",
        ),
        prefix="Scaffold_",
        key_column="Scaffold_Size",
    ),
    "clustering": CategorySpec(
        key="clustering",
        name="Similarity clusters",
        description="""Similarity clustering

Similarity clustering groups compounds using Morgan fingerprints and Tanimoto
similarity (Butina clustering). This helps identify closely related analog
series and provides quick diversity context.

Run `python analyzers/analyzer_clustering.py results.csv` or the combined runner
`python tools/series_analytics.py results.csv` to generate these columns.
""",
        columns=(
            "Cluster_ID",
            "Cluster_Size",
            "Cluster_Representative_ID",
            "Cluster_MedoidMeanSim",
        ),
        prefix="Cluster_",
        key_column="Cluster_Size",
    ),
    "developability": CategorySpec(
        key="developability",
        name="Developability indices",
        description="""Developability indices

This category bundles a few widely used, interpretable property-space rules and
an aggregated heuristic risk/score intended for early developability triage.

Implemented rules:
- GSK 4/400: MW <= 400 and LogP <= 4
- Pfizer 3/75: TPSA <= 75 and LogP <= 3
- Golden Triangle (coarse): 200 <= MW <= 500 and 1 <= cLogD_7.4 <= 3

The composite Dev_Risk_Burden and Dev_Score incorporate available red flags
(toxicity class, hERG class, CYP burden, clearance/solubility risk when present,
Ro5 violations, and med-chem alert counts). These are intended for ranking and
triage rather than as hard decision gates.
""",
        columns=(
            "Dev_Score",
            "Dev_Class",
            "Dev_Risk_Class",
            "Dev_Risk_Burden",
            "Dev_BasePoints",
            "Dev_GSK_4_400_Pass",
            "Dev_Pfizer_3_75_Pass",
            "Dev_GoldenTriangle_Pass",
            "Dev_Solubility_Risk",
            "Dev_Flags",
        ),
        prefix="Dev_",
        key_column="Dev_Score",
    ),
    "shape3d": CategorySpec(
        key="shape3d",
        name="3D shape descriptors",
        description="""3D shape descriptors

These columns are optional outputs from `calculators/mpo_v5.py` when 3D mode is
enabled (`--3d`). The calculator embeds conformers with RDKit ETKDG, optionally
minimizes them (MMFF94s or UFF), selects the lowest-energy conformer, and then
computes a small set of RDKit 3D shape descriptors (NPR1/NPR2, PMI1-3, radius of
gyration, etc.).

These descriptors are useful for rough shape comparisons and to help flag very
flat vs very globular molecules, but they are sensitive to conformer generation
settings and should be treated as qualitative ranking features.
""",
        columns=(
            "3D_Mode",
            "3D_NumConfs",
            "3D_SelectedConfId",
            "3D_Minimizer",
            "3D_Converged",
            "3D_BestEnergy",
            "3D_Asphericity",
            "3D_Eccentricity",
            "3D_InertialShapeFactor",
            "3D_NPR1",
            "3D_NPR2",
            "3D_PMI1",
            "3D_PMI2",
            "3D_PMI3",
            "3D_RadiusOfGyration",
            "3D_SpherocityIndex",
            "3D_PBF",
        ),
        prefix="3D_",
        key_column="3D_RadiusOfGyration",
    ),
}


def list_categories() -> List[str]:
    return list(CATEGORY_SPECS.keys())

