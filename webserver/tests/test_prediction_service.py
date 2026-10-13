"""Tests for live AGVarPred feature-schema and annotation guards."""

import pandas as pd

from webserver.prediction_service import (
    alpha_feature_names,
    has_vep_consequence_annotation,
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


def test_absent_consequence_one_hot_column_can_be_zero_filled():
    features = pd.DataFrame({"vep_Consequence": ["missense_variant"]})
    selected = ["vep_Consequence_intron_variant"]

    assert missing_model_feature_columns(features, selected) == []


def test_vep_guard_requires_a_nonempty_raw_consequence():
    features = pd.DataFrame(
        {
            "vep_Consequence": ["missense_variant", None, ""],
        },
        index=["annotated", "missing", "empty"],
    )

    assert has_vep_consequence_annotation(features, "annotated")
    assert not has_vep_consequence_annotation(features, "missing")
    assert not has_vep_consequence_annotation(features, "empty")
    assert not has_vep_consequence_annotation(features, "not_in_matrix")
