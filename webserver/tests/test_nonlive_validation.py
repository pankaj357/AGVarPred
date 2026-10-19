"""Non-live validation checks that do not need AlphaGenome API access.

These tests use the bundled frozen model and the local reference feature matrix
(`parity-reference-test-part0.parquet`) to verify that:
- the reference matrix matches the model's selected-feature schema,
- variant IDs follow the canonical format,
- the frozen model returns finite, calibrated scores on real reference rows,
- the predictor output contains the expected metadata columns.

They do not validate live AlphaGenome feature generation or claim end-to-end
live parity.
"""

from __future__ import annotations

import re

import numpy as np
import pandas as pd
import pytest

from AGVarPred import AGVarPredPredictor
from agvarpred_core.utils import parse_variant_id

REFERENCE_PATH = "parity-reference-test-part0.parquet"


@pytest.fixture
def predictor():
    return AGVarPredPredictor(model_name="model_full", af_source_name="webserver")


@pytest.fixture
def reference():
    df = pd.read_parquet(REFERENCE_PATH)
    df["variant_id"] = df["variant_id"].astype(str)
    return df


def test_reference_matrix_exists_and_has_variant_id():
    df = pd.read_parquet(REFERENCE_PATH)
    assert "variant_id" in df.columns
    assert not df["variant_id"].duplicated().any()


def test_reference_contains_all_selected_model_features(predictor, reference):
    missing = [name for name in predictor.selected_features if name not in reference.columns]
    assert not missing, f"Reference is missing {len(missing)} selected features"


def test_reference_variant_ids_are_canonical(reference):
    pattern = re.compile(r"^chr[0-9XYM]+_\d+_[ACGTN]+_[ACGTN]+$")
    for variant_id in reference["variant_id"]:
        assert pattern.match(variant_id), f"Unexpected variant_id format: {variant_id}"
        chrom, pos, ref, alt = parse_variant_id(variant_id)
        assert chrom.startswith("chr")
        assert pos > 0
        assert set(ref) <= set("ACGTN")
        assert set(alt) <= set("ACGTN")
        assert ref != alt


def test_predictor_on_reference_features_returns_finite_scores(predictor, reference):
    """Score reference rows directly; the frozen model must return finite probabilities."""
    indexed = reference.set_index("variant_id", drop=False)
    predictions = predictor.predict(indexed)

    assert len(predictions) == len(reference)
    assert np.isfinite(predictions["probability"].to_numpy()).all()
    assert predictions["predicted_class"].isin([0, 1]).all()
    assert list(predictions["variant_id"]) == list(reference["variant_id"])


def test_predictor_output_includes_expected_metadata_columns(predictor, reference):
    indexed = reference.set_index("variant_id", drop=False)
    predictions = predictor.predict(indexed)

    expected_columns = {
        "variant_id",
        "probability",
        "predicted_class",
        "model_version",
        "model_type",
        "af_source",
        "feature_version",
        "alphagenome_version",
        "gnomAD_version",
        "VEP_version",
        "prediction_timestamp",
    }
    assert expected_columns.issubset(set(predictions.columns))
    assert predictions["model_version"].iloc[0] == "model_full"
    assert predictions["model_type"].iloc[0] == "full"


def test_manifest_matches_predictor_properties(predictor):
    assert predictor.model_version == "model_full"
    assert predictor.model_type == "full"
    assert predictor.n_features == 120
    assert predictor.threshold == pytest.approx(0.42)
