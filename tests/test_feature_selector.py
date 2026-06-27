import numpy as np
import pandas as pd

from agvarpred_core.feature_selector import FeatureSelector


def test_feature_selector_adds_missing_columns():
    selector = FeatureSelector(["a", "b", "c"])
    df = pd.DataFrame({"a": [1.0, 2.0], "b": [3.0, np.nan]})
    out = selector.select(df)
    assert list(out.columns) == ["a", "b", "c"]
    assert np.isnan(out["c"].iloc[0])


def test_feature_selector_drops_extra_columns():
    selector = FeatureSelector(["a", "b"])
    df = pd.DataFrame({"a": [1], "b": [2], "d": [3]})
    out = selector.select(df)
    assert list(out.columns) == ["a", "b"]
    assert "d" not in out.columns
