"""Unit tests for parity-comparison helpers."""

import numpy as np
import pandas as pd

from webserver.validate_feature_parity import compare_selected_features


def test_feature_comparison_accepts_equal_values_and_matching_nans():
    reference = pd.DataFrame(
        {"track_a": [0.1, np.nan], "vep_SIFT_score": [0.8, np.nan]},
        index=["v1", "v2"],
    )
    generated = reference.copy()

    result = compare_selected_features(
        reference, generated, ["track_a", "vep_SIFT_score"], tolerance=1e-8
    )

    assert result["feature_mismatches"] == 0
    assert result["feature_match_rate"] == 1.0


def test_feature_comparison_detects_feature_drift():
    reference = pd.DataFrame({"track_a": [0.1], "gnomAD_AF": [0.01]}, index=["v1"])
    generated = pd.DataFrame({"track_a": [0.2], "gnomAD_AF": [0.01]}, index=["v1"])

    result = compare_selected_features(
        reference, generated, ["track_a", "gnomAD_AF"], tolerance=1e-6
    )

    assert result["feature_mismatches"] == 1
    assert result["mismatch_examples"][0]["feature"] == "track_a"
