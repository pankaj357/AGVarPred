"""Tests for live AGVarPred feature-schema and annotation guards."""

import numpy as np
import pandas as pd

from webserver.prediction_service import (
    alpha_feature_names,
    apply_training_missing_annotation_defaults,
    missing_model_feature_columns,
)


def test_alpha_feature_names_excludes_vep_and_frequency_features():
    selected = [
        "gnomAD_AF",
        "log10_gnomAD_AF",
        "vep_IMPACT_score",
        "ATAC_track_1",
        "SPLICE_track_2",
    ]

    assert alpha_feature_names(selected) == ["ATAC_track_1", "SPLICE_track_2"]


def test_missing_alpha_feature_column_is_reported():
    features = pd.DataFrame({"ATAC_track_1": [0.2]})
    selected = ["ATAC_track_1", "SPLICE_track_2"]

    assert missing_model_feature_columns(features, selected) == ["SPLICE_track_2"]


def test_training_vEP_defaults_are_applied_when_entire_batch_lacks_vEP():
    selected = [
        "vep_IMPACT_score",
        "vep_LoF",
        "vep_PolyPhen_score",
        "vep_SIFT_score",
        "vep_is_synonymous",
        "vep_Consequence_synonymous_variant",
    ]
    features = pd.DataFrame({"ATAC_track_1": [0.2]}, index=["chr1_10_A_G"])

    result = apply_training_missing_annotation_defaults(features, selected)
    row = result.loc["chr1_10_A_G"]

    assert row["vep_IMPACT_score"] == 0
    assert row["vep_LoF"] == -1
    assert np.isnan(row["vep_PolyPhen_score"])
    assert np.isnan(row["vep_SIFT_score"])
    assert row["vep_is_synonymous"] == 0
    assert row["vep_Consequence_synonymous_variant"] == 0
    assert missing_model_feature_columns(result, selected) == []


def test_missing_vep_row_has_same_defaults_in_mixed_batch_and_alone():
    selected = [
        "vep_IMPACT_score",
        "vep_LoF",
        "vep_PolyPhen_score",
        "vep_SIFT_score",
        "vep_is_synonymous",
        "vep_Consequence_synonymous_variant",
    ]
    mixed = pd.DataFrame(
        {
            "vep_IMPACT_score": [4, np.nan],
            "vep_LoF": [1, np.nan],
            "vep_PolyPhen_score": [0.01, np.nan],
            "vep_SIFT_score": [0.02, np.nan],
            "vep_is_synonymous": [0, np.nan],
            "vep_Consequence_synonymous_variant": [1, np.nan],
        },
        index=["annotated", "missing"],
    )

    result = apply_training_missing_annotation_defaults(mixed, selected)
    row = result.loc["missing"]

    assert row["vep_IMPACT_score"] == 0
    assert row["vep_LoF"] == -1
    assert np.isnan(row["vep_PolyPhen_score"])
    assert np.isnan(row["vep_SIFT_score"])
    assert row["vep_is_synonymous"] == 0
    assert row["vep_Consequence_synonymous_variant"] == 0


def test_unknown_vep_feature_missing_rule_fails_closed():
    features = pd.DataFrame(index=["variant"])
    selected = ["vep_new_unverified_feature"]

    try:
        apply_training_missing_annotation_defaults(features, selected)
    except RuntimeError as exc:
        assert "No training-compatible missing-value rule" in str(exc)
    else:
        raise AssertionError("An unknown VEP feature must not receive a guessed default")
