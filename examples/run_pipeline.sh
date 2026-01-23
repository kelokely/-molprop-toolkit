#!/usr/bin/env bash
set -euo pipefail

# Example end-to-end run.
# This script assumes you have RDKit installed for the calculators.

INPUT=${1:-examples/example.smi}
OUTDIR=${2:-outputs}

mkdir -p "$OUTDIR"

python calculators/mpo_v4.py "$INPUT" -o "$OUTDIR/results.csv" --fill-missing 0

python analyzers/analyzer_master.py "$OUTDIR/results.csv" --report -o "$OUTDIR/report.txt"
python analyzers/analyzer_toxicity.py "$OUTDIR/results.csv" -o "$OUTDIR/toxicity.csv"
python analyzers/analyzer_medchem.py "$OUTDIR/results.csv" -o "$OUTDIR/medchem.csv"
python analyzers/analyzer_cns.py "$OUTDIR/results.csv" -o "$OUTDIR/cns.csv"

