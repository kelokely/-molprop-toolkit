# Category guide

This document explains what each analysis category in MolProp Toolkit is trying to capture, why the category exists in early drug discovery, how to interpret the common thresholds, and how to use the output to decide what to do next. The analyzers are intentionally pragmatic. They are meant to help you find compounds that are clean and promising, compounds that are likely to create trouble later, and compounds that are worth optimizing because the issues are addressable.

## CNS MPO

CNS MPO is a multi-parameter score created to describe a region of property space where central nervous system drugs tend to succeed. It combines molecular weight, lipophilicity, distribution coefficient at physiological pH, polar surface area, hydrogen-bond donors, and basic pKa. The appeal of this score is not that it is perfect, but that it encodes a set of trade-offs that medicinal chemists repeatedly see in practice: as lipophilicity and basicity increase, passive permeability can improve but off-target binding, hERG risk, and metabolic clearance often worsen; as polarity increases, solubility can improve but brain penetration and cellular permeability often worsen. In the analyzer output, a high score is a signal that the compound is balanced for CNS delivery. A low score is a signal that one or more properties are far from typical CNS drug space, and the component breakdown helps you decide whether the issue is a heavy molecule, an overly polar molecule, an overly basic molecule, or a lipophilicity problem.

## Oral bioavailability filters

Oral bioavailability is governed by dissolution, permeability, stability in the gastrointestinal tract, and first-pass metabolism in gut and liver. No single rule predicts all of this, so the toolkit includes several widely used filters. Veber emphasizes polar surface area and rotatable bonds, which approximate permeability and conformational flexibility. Egan provides a slightly different boundary in TPSA–LogP space that correlates with absorption. Ghose and Muegge add additional constraints motivated by distributions seen in known drugs. The oral bioavailability score reported by the calculator is a simple count of how many of these filters pass, which is convenient for ranking and triage. A low score does not guarantee failure, but it should prompt you to ask whether the project can accept more formulation effort, higher doses, or alternative routes.

## Rule of 5 and Rule of 3

The Rule of 5 is a classic heuristic for oral small molecules. It is useful because it is fast, interpretable, and often correlates with permeability and solubility constraints. It is not a law, and there are many important drug classes that violate it, but it is still a strong early filter for many discovery programs. The Rule of 3 is oriented toward fragment-based drug discovery, where small, simple compounds are used as starting points and are grown or linked during optimization. Rule of 3 compliance is not required for a compound to be a good drug, but it is a useful way to identify molecules that are small enough to be efficient starting points, especially when paired with ligand efficiency metrics.

## Med-chem flags

Med-chem flags are structural filters that highlight motifs associated with assay interference, instability, or recurring development problems. PAINS patterns were popularized as a way to reduce false positives in screening campaigns, particularly when readouts are sensitive to redox cycling, aggregation, metal chelation, or covalent reactivity. Brenk and related filters represent collections of fragments that medicinal chemists often avoid because they are reactive, unstable, or historically problematic. In this toolkit, these flags are treated as risk indicators rather than absolute rejection criteria. The most important best practice is that a flagged compound should not be trusted until it has passed orthogonal assays that eliminate common artifacts.

## Toxicity alerts

The toxicity analysis in this repository is alert based rather than model based. It looks for functional groups that are electrophilic, groups that are known mutagenicity red flags, and motifs that can be bioactivated to reactive metabolites. These alerts are not predictions of clinical toxicity; they are prompts to investigate. An epoxide alert does not mean a compound will be toxic at therapeutic exposure, but it does mean the chemist should consider whether the epoxide is essential, whether it can be replaced, and whether glutathione trapping experiments should be prioritized. Similarly, nitro-aromatic motifs are historically associated with mutagenicity risk and should trigger Ames testing early.

## Metabolism and soft spots

Metabolic soft spots identify likely sites of CYP oxidation, dealkylation, or hydrolysis, as well as obvious phase II conjugation handles. The purpose is not to produce a precise clearance number, but to anticipate liabilities and prioritize experiments such as microsomal stability, hepatocyte stability, metabolite identification, and reactive metabolite trapping. A compound with many benzylic or allylic sites, multiple N-dealkylation handles, and high lipophilicity is often a candidate for rapid clearance. When a project needs longer exposure, these features suggest concrete medicinal chemistry moves, such as blocking benzylic positions, replacing methoxy groups, reducing LogP, or changing ester linkers.

## hERG

hERG liability remains one of the most important safety risks for small molecules. The toolkit’s hERG analysis is again heuristic, using the presence of a basic nitrogen, aromatic or hydrophobic features, and property thresholds that often correlate with higher hERG binding probability. The output is best used as a screening signal to prioritize early patch clamp or binding assays. When a series shows persistent medium or high risk, common mitigation levers include reducing basicity, reducing lipophilicity, increasing three-dimensionality, and avoiding motifs that strongly bind into the channel pharmacophore.

## CYP inhibition and mechanism-based inactivation

Inhibition of major CYP isoforms can create clinically meaningful drug–drug interactions. The toolkit provides binary alerts for inhibition risk and a separate signal for mechanism-based inactivation patterns, which are especially important because they can lead to time-dependent inhibition. If a compound is flagged, the most appropriate next steps are experimental CYP inhibition panels and, for MBI-like motifs, time-dependent inhibition experiments. From a design standpoint, reducing basicity and lipophilicity can help for CYP2D6-like inhibitors, while removal of certain heteroaromatic motifs can reduce CYP3A4 inhibition risk.

## Solubility

Solubility is frequently the limiting factor for oral exposure, especially for lipophilic, aromatic, and high-molecular-weight compounds. The toolkit includes simple logS estimators that are most useful for ranking within a series rather than for making absolute claims. A low predicted logS should trigger thinking about salt formation, pKa engineering, reduction of aromatic surface, increasing sp3 character, and formulation options such as amorphous dispersions, lipid formulations, or particle size reduction.

## Permeability

Permeability predictions in the toolkit are coarse estimates intended for prioritization. The most useful signals are the GI absorption classification and the P-gp substrate likelihood. A compound with good solubility but poor permeability may require polarity reduction or prodrug strategies. A compound with high P-gp likelihood may show good permeability in PAMPA but poor cellular permeability, and may fail to reach the brain even if it has otherwise acceptable CNS MPO scores.

## Pharmacokinetics heuristics

The PK section includes simple estimates and categories for protein binding, distribution, clearance risk, half-life class, and an oral bioavailability class. These are not substitutes for measured PK, but they are useful for communicating the expected direction of risk and for tracking progress across optimization cycles. When the PK summary flags multiple issues, such as low solubility, high clearance, and high P-gp likelihood, the program likely needs a more substantial redesign rather than incremental optimization.

## Lead optimization metrics

Lead metrics such as ligand efficiency and lipophilic ligand efficiency are most meaningful when you have experimental potency values, but they can also be used as conceptual targets to avoid making molecules bigger and more lipophilic without real gains. If the same potency can be achieved at lower LogP, the resulting compound is often easier to develop. If potency gains are accompanied by large increases in heavy atom count and decreases in ligand efficiency, the series is often drifting into an unattractive region of chemical space.



## Synthetic accessibility and complexity

Synthetic accessibility is an operational concern: even when a compound has attractive ADME or potency signals, the project may stall if the route is long, the chemistry is fragile, or key building blocks are expensive or unavailable. MolProp Toolkit includes the Ertl–Schuffenhauer SA score (scaled 1–10, where lower is easier) as a fast, ranking-friendly indicator. It combines fragment contributions (frequent medicinal chemistry fragments tend to be easier) with penalties for structural complexity such as spiro/bridgehead content, stereocenters, and macrocycles. The score is not a substitute for retrosynthetic analysis, but it is useful for highlighting molecules that are likely to require more chemistry effort, and for tracking whether a series is becoming harder to make over time.

The accompanying Complexity_* descriptors are intended to explain the direction of the SA score by pointing to specific complexity drivers: graph complexity indices (e.g., BertzCT), spiro and bridgehead atom counts, stereocenter burden (including unspecified stereocenters), and aromatic atom fraction as a rough proxy for flatness. When a compound looks "too hard" by these heuristics, common mitigation strategies include reducing ring fusion complexity, simplifying stereochemistry, avoiding spiro/bridged scaffolds when not essential, and controlling aromaticity by introducing saturated ring replacements where appropriate.



## Developability indices

Developability is a pragmatic umbrella concept that covers the many reasons a compound may be difficult to turn into a robust, manufacturable, and clinically usable drug candidate. Some of the most common early drivers are property-space issues like excessive size, lipophilicity, or polarity that create solubility and permeability challenges, as well as red flags like structural alerts, high hERG liability risk, high CYP inhibition burden, or rapid metabolic clearance. Because these factors trade off against each other and vary by modality, the toolkit implements a small set of interpretable “rule-style” checks and a lightweight composite score intended for ranking and triage rather than as a hard gate.

The rule-style checks include the classic GSK “4/400” filter (MW ≤ 400 and LogP ≤ 4), the Pfizer “3/75” filter (LogP ≤ 3 and TPSA ≤ 75), and a coarse “golden triangle” check using MW and cLogD at pH 7.4. The composite developability outputs (`Dev_*` columns) then combine these property-space checks with whichever risk flags are available in the dataset (toxicity class, hERG class, CYP burden, Ro5 violations, med‑chem alerts, and, when present, solubility and clearance risk) to create an interpretable risk burden and a simple score.

## 3D shape descriptors

The optional 3D descriptor mode in the v5 calculator adds a small set of conformer-derived shape descriptors that are useful for rough geometry comparisons across a compound set. This mode uses RDKit ETKDG to generate multiple conformers per compound, optionally minimizes them with MMFF94s or UFF, selects the lowest-energy conformer, and then computes principal-moment and normalized principal moment (PMI/NPR) metrics along with other compact shape descriptors such as radius of gyration and asphericity. These descriptors can be helpful for flagging overly flat molecules, prioritizing shape diversity, or sanity-checking whether two 2D-similar compounds occupy similar 3D shape space. Because the results depend on conformer generation settings and force-field assumptions, they should be interpreted as qualitative ranking signals rather than accurate physical measurements.

## Series analytics: scaffolds and clustering

Series analytics is where the toolkit shifts from per-compound scoring to questions that matter in real medicinal chemistry cycles: whether a hit list represents one chemotype or several, whether improvements are happening within a consistent series, and whether different scaffolds exhibit systematically different liability profiles. The scaffold and clustering tooling is designed to be a lightweight way to add this context to a results table without requiring a full cheminformatics platform.

Murcko scaffolds reduce each compound to a core ring–linker framework, which makes it easier to see which chemotypes dominate a dataset and to compute per-scaffold property summaries. A common workflow is to examine the largest scaffolds first, then compare the mean or median CNS MPO, solubility, and risk flags across scaffolds to decide which chemotypes deserve the most investment. Outlier inspection within a scaffold can reveal whether the series has "headroom" to improve, for example when one analog has much higher MPO or lower risk than its close neighbors.

Similarity clustering groups compounds based on Morgan fingerprints and Tanimoto similarity. In practice this is a way to find close analog series even when the Murcko scaffold definition splits or merges chemotypes in unintuitive ways. Cluster size distributions can highlight whether a dataset is diverse or dominated by a few analog series, and cluster representatives can be used as a compact visual index for chemistry review.

