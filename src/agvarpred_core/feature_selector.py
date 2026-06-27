"""FeatureSelector: reduce the full engineered matrix to a model's subset."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from .utils import clean_name, load_selected_features


class FeatureSelector:
    """Select the model-specific feature subset from a full feature matrix.

    Parameters
    ----------
    selected_features:
        Either a list of feature names or a path to a newline-delimited
        feature list file.
    """

    def __init__(self, selected_features: Iterable[str] | str | Path):
        if isinstance(selected_features, (str, Path)):
            self.selected_features = load_selected_features(selected_features)
        else:
            self.selected_features = [clean_name(f) for f in selected_features]

    def select(self, df: pd.DataFrame) -> pd.DataFrame:
        """Return ``df`` with only the selected feature columns, in order.

        Missing selected columns are filled with NaN so that downstream
        imputers can handle them consistently.
        """
        df = df.copy()
        for col in self.selected_features:
            if col not in df.columns:
                df[col] = np.nan
        return df[self.selected_features]

    @property
    def n_features(self) -> int:
        return len(self.selected_features)
