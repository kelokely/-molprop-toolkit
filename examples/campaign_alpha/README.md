# Campaign Alpha Validation Pack

This synthetic example pack is meant to exercise the second-level MolScope workflows on something closer to a real campaign loop than `examples/example.smi`.

Included inputs:

- `round4_results.csv` and `round5_results.csv`: MolScope-compatible campaign tables with scaffold, cluster, route, and liability context.
- `round4_assays.csv` and `round5_assays.csv`: measured potency and small ADME/safety follow-up panels that can be joined back into those tables with `molscope integrate`.

The intended validation story is:

1. Join measured assay data onto both rounds with `molscope integrate`.
2. Compare the rounds with `molscope compare`.
3. Build a report and decision picklists for the current round.
4. Generate SAR and matched-pair bundles.
5. Distill those artifacts into `molscope learnings`, persist them with `molscope dashboard`, and assemble everything with `molscope portal`.

The benchmark pack is designed so the expected high-level outcomes are easy to inspect:

- shared compounds improve on measured potency and potency-normalized balance in round 5
- solubility gets worse for some shared compounds, so compare should surface at least one regression
- `Scaffold_ID = SCF_C` expands from round 4 to round 5
- low coverage remains on at least one measured endpoint, so round-review should still emit an open question

To build the full demo bundle:

```bash
python examples/build_campaign_alpha_demo.py
```
