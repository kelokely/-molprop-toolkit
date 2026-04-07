from __future__ import annotations

import pandas as pd

from molprop_toolkit.picklists.builtins import BUILTIN_PICKLISTS
from molprop_toolkit.picklists.engine import run_picklist


def test_retro_picklist_builtins_use_retro_columns() -> None:
    df = pd.DataFrame(
        {
            "Compound_ID": ["CMPD_001", "CMPD_002", "CMPD_003"],
            "Retro_RecordMatched": [True, True, False],
            "Retro_IsSolved": [True, False, False],
            "Retro_Rerank_Rank": [3, 24, 80],
            "Retro_Rerank_Score": [0.88, 0.22, 0.05],
            "Retro_Stock_Coverage": [1.0, 0.2, 0.0],
            "Retro_StepCount": [3, 9, 10],
            "Retro_NotInStockCount": [0, 4, 5],
            "Retro_HTMLPage": ["targets/a.html", "targets/b.html", "targets/c.html"],
            "Lead_LLE": [4.8, 3.1, 2.2],
            "CNS_MPO": [4.2, 3.5, 2.1],
        }
    )

    make_next = run_picklist(df, BUILTIN_PICKLISTS["make_next_synthesis_ready"])
    route_risk = run_picklist(df, BUILTIN_PICKLISTS["route_risk_deprioritize"])

    assert list(make_next.df["Compound_ID"]) == ["CMPD_001"]
    assert set(route_risk.df["Compound_ID"]) == {"CMPD_002", "CMPD_003"}
    assert make_next.df.loc[make_next.df.index[0], "Decision_Queue"] == "make_next"
    assert "Decision_Rationale" in make_next.df.columns
    assert "Retro_Rerank_Rank" in make_next.df.loc[make_next.df.index[0], "Decision_Evidence_Summary"]
