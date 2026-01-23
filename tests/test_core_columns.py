from __future__ import annotations

import pandas as pd

from molprop_toolkit.core.columns import detect_best_smiles_column, detect_id_column


def test_detect_id_column_prefers_compound_id() -> None:
    df = pd.DataFrame({"Compound_ID": ["A"], "SMILES": ["CCO"]})
    assert detect_id_column(df) == "Compound_ID"


def test_detect_best_smiles_column_priority() -> None:
    cols = ["SMILES", "Canonical_SMILES", "Calc_Canonical_SMILES"]
    assert detect_best_smiles_column(cols) == "Calc_Canonical_SMILES"

