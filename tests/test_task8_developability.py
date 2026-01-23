from __future__ import annotations

from molprop_toolkit.core.registry import CATEGORY_SPECS
from molprop_toolkit.schema import load_schema


def test_registry_has_developability_category() -> None:
    assert "developability" in CATEGORY_SPECS


def test_schema_has_dev_columns() -> None:
    schema = load_schema(None)
    cols = schema.get("columns", {})
    for c in [
        "Dev_Score",
        "Dev_Class",
        "Dev_Risk_Burden",
        "Dev_GSK_4_400_Pass",
        "Dev_Pfizer_3_75_Pass",
    ]:
        assert c in cols

