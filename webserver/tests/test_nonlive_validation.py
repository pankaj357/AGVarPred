"""Non-live model tests using a synthetic in-memory feature matrix.

These tests check schema compatibility and inference mechanics only.
They do not validate live AlphaGenome parity or biological accuracy.
"""

import re

import numpy as np
import pandas as pd
import pytest

from AGVarPred import AGVarPredPredictor
from agvarpred_core.utils import parse_variant_id


@pytest.fixture
def predictor():
    return AGVarPredPredictor(
        model_name="model_full",
        af_source_name="webserver-test",
    )


@pytest.fixture
def reference(predictor):
    variant_ids = ["chr1_12345_A_G", "chr13_200_C_T"]
    features = pd.DataFrame(
        np.zeros((len(variant_ids), predictor.n_features), dtype=float),
        columns=predictor.selected_features,
    )
    features.insert(0, "variant_id", variant_ids)
    return features


def test_synthetic_matrix_has_unique_variant_ids(reference):
    assert "variant_id" in reference.columns
    assert not reference["variant_id"].duplicated().any()


def test_matrix_contains_all_selected_model_features(predictor, reference):
    missing = [
        name for name in predictor.selected_features
        if name not in reference.columns
    ]
    assert not missing, f"Missing {len(missing)} selected features: {missing[:8]}"


def test_variant_ids_are_canonical(reference):
    pattern = re.compile(r"^chr[0-9XYM]+_\d+_[ACGTN]+_[ACGTN]+$")
    for variant_id in reference["variant_id"]:
        assert pattern.fullmatch(variant_id), variant_id
        chrom, pos, ref, alt = parse_variant_id(variant_id)
        assert chrom.startswith("chr")
        assert pos > 0
        assert set(ref) <= set("ACGTN")
        assert set(alt) <= set("ACGTN")
        assert ref != alt


def test_predictor_returns_finite_scores(predictor, reference):
    indexed = reference.set_index("variant_id", drop=False)
    predictions = predictor.predict(indexed)

    assert len(predictions) == len(reference)
    assert np.isfinite(predictions["probability"].to_numpy()).all()
    assert predictions["predicted_class"].isin([0, 1]).all()
    assert list(predictions["variant_id"]) == list(reference["variant_id"])


def test_prediction_metadata_columns(predictor, reference):
    indexed = reference.set_index("variant_id", drop=False)
    predictions = predictor.predict(indexed)

    expected = {
        "variant_id", "probability", "predicted_class", "model_version",
        "model_type", "af_source", "feature_version", "alphagenome_version",
        "gnomAD_version", "VEP_version", "prediction_timestamp",
    }
    assert expected.issubset(predictions.columns)
    assert predictions["model_version"].iloc[0] == "model_full"
    assert predictions["model_type"].iloc[0] == "full"


def test_manifest_matches_predictor(predictor):
    assert predictor.model_version == "model_full"
    assert predictor.model_type == "full"
    assert predictor.n_features == 120
    assert predictor.threshold == pytest.approx(0.42)
