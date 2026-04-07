# MolScope Toolkit Roadmap

This roadmap turns the current product direction into issue-ready work. The focus is to make MolScope more useful for medicinal-chemistry decisions rather than adding another descriptor-heavy calculator profile.

## Milestone 0: Cleanup the public command surface

Title: `Finish molscope rename cleanup across detail docs and examples`

Problem:
Public docs still leak legacy command names and uneven validator examples, which makes the product feel half-migrated even though `molscope` is already the canonical workflow wrapper.

Scope:
- Make `molscope` the first command shown on README and detail pages.
- Keep compatibility aliases working, but stop leading with them in public docs.
- Expose schema validation as `molscope schema`.
- Remove accidental local-path references and remaining MolProp branding from public-facing outputs.

Acceptance criteria:
- README, quickstart, and detail pages lead with `molscope`.
- Schema validation docs use `molscope schema`.
- Generated report and similarity HTML use `MolScope Toolkit` branding.
- Compatibility aliases still function.

## Milestone 1: Add scaffold-centric SAR

Title: `Add molscope sar for scaffold-centric R-group SAR tables`

Problem:
`molscope series` currently stops at scaffolds and clusters. Med-chem users still need a way to inspect substitution patterns and property deltas around a common core.

Scope:
- Add `molscope sar` as a new workflow command.
- Use RDKit scaffold decomposition / R-group decomposition to build per-series tables.
- Export CSV or Parquet outputs and an HTML SAR view with aligned depictions.
- Surface per-position deltas for key endpoints such as MPO, solubility, permeability, hERG risk, and developability.

Acceptance criteria:
- Users can provide a scaffold or let the tool infer series cores.
- Output contains R-group columns and per-series summary metrics.
- HTML output is readable for medicinal chemists, not just raw descriptor tables.
- SAR outputs can be linked from reports and picklists.

## Milestone 2: Add matched molecular pair analysis

Title: `Add molscope mmp for matched molecular pair analysis`

Problem:
The toolkit can score compounds, but it does not yet turn those scores into “what change improved the property?” design rules.

Scope:
- Add `molscope mmp` using RDKit MMP fragmentation.
- Build pair tables with structural transform notation and endpoint deltas.
- Add aggregate transform summaries for repeated changes across a library.
- Prioritize med-chem endpoints such as solubility, clearance heuristics, hERG risk, CYP risk, MPO, and lipophilicity.

Acceptance criteria:
- Users can generate matched-pair tables from a results dataset.
- Output includes transform identifiers, left/right compounds, and delta columns.
- Aggregate transform summaries rank recurring edits by consistency and effect size.
- Reports can reference top beneficial and harmful transforms.

## Milestone 3: Add exact and substructure search

Title: `Add molscope search for exact and substructure search`

Problem:
Fingerprint similarity is useful, but medicinal chemists also need exact SMARTS and motif-driven retrieval.

Scope:
- Add `molscope search` for exact SMILES, SMARTS substructure, and motif-panel workflows.
- Support searching MolScope result tables and plain SMILES libraries.
- Reuse RDKit substructure and, where appropriate, scalable library search primitives.
- Allow custom motif panels for project-specific alerts or chemotype checks.

Acceptance criteria:
- Users can run exact and SMARTS-based search from the CLI.
- Outputs include match counts and matching compounds.
- Search supports both simple one-off queries and reusable motif files.
- Search results can feed into picklists or reports.

## Milestone 4: Integrate measured assay data

Title: `Add molscope integrate for assay joins and real LE/LLE recalculation`

Problem:
Potency-dependent metrics are most useful when measured potency and assay data are joined into the same table, and that workflow is currently missing.

Scope:
- Add `molscope integrate` to merge measured assay endpoints into MolScope outputs.
- Support joins by compound ID and configurable alternate identifiers.
- Recalculate potency-dependent metrics such as LE, LLE, LELP, BEI, and SEI from measured values when available.
- Mark whether a field is measured, inferred, or heuristic.

Acceptance criteria:
- Users can merge a measured assay table into a MolScope table with one command.
- Recomputed potency-dependent metrics are emitted in a clear, auditable way.
- Provenance columns distinguish measured values from heuristic outputs.
- Integrated outputs remain compatible with report and picklist generation.

## Milestone 5: Compare design rounds

Title: `Add molscope compare for round-over-round campaign comparison`

Problem:
Teams need to compare round N versus round N+1 or series A versus series B, not just inspect one static table at a time.

Scope:
- Add `molscope compare` for dataset-versus-dataset or round-versus-round comparison.
- Summarize changes by scaffold, cluster, and priority endpoints.
- Highlight where potency improved while liabilities worsened, or vice versa.
- Reuse integrated assay-aware outputs when available.

Acceptance criteria:
- Users can compare two result tables directly.
- Output includes both summary statistics and scaffold- or cluster-level breakdowns.
- Reports clearly surface positive and negative tradeoffs between rounds.
- The command works with both heuristic-only tables and assay-integrated tables.

## Milestone 6: Make retrosynthesis actionable

Title: `Add retro reranking and report/picklist integration`

Problem:
`molscope retro` already produces route summaries, but synthesis feasibility is not yet threaded back into prioritization.

Scope:
- Add route reranking by score, stock availability, precursor burden, and step count.
- Write route burden metrics back into the enriched results table.
- Surface synthesis-aware sorting and flags in reports and picklists.
- Keep the existing AiZynthFinder-based workflow intact.

Acceptance criteria:
- Users can rerank routes with a clear policy.
- Route summaries can be merged back into report and picklist workflows.
- Picklists can sort or filter by synthesis feasibility alongside ADME and tox criteria.
- Outputs remain reproducible from route summary files.

## Milestone 7: Add focused analog design

Title: `Add molscope design for focused analog enumeration`

Problem:
After SAR and MMP workflows exist, the next logical step is suggesting what to make next.

Scope:
- Add `molscope design` for focused analog enumeration.
- Start with bounded workflows such as reaction-based enumeration, stereoisomer generation, and heterocycle hopping.
- Feed generated designs back into `molscope calc`, `molscope report`, and `molscope picklists`.
- Keep the output reviewable and chemistry-aware rather than unconstrained enumeration.

Acceptance criteria:
- Users can generate a focused analog set from a seed or scaffold.
- Generated designs can move through the existing scoring and reporting workflow.
- Enumeration rules are explicit and reproducible.
- The command defaults to conservative, reviewable proposal counts.

## Recommended order

1. Cleanup command surface and finish the rename.
2. Add `molscope integrate`.
3. Add `molscope compare`.
4. Add `molscope sar`.
5. Add `molscope mmp`.
6. Add `molscope search`.
7. Add retro reranking.
8. Add `molscope design`.
