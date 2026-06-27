"""Shared test fixtures."""

import json
import shutil
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import pytest
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import RobustScaler


def _make_synthetic_pipeline(selected_features, threshold=0.5):
    rng = np.random.RandomState(42)
    n = 200
    X = pd.DataFrame(
        rng.randn(n, len(selected_features)),
        columns=selected_features,
    )
    # Make first feature predictive
    y = (X.iloc[:, 0] + X.iloc[:, 1] * 0.5 > 0).astype(int)

    imputer = X.median()
    scaler = RobustScaler()
    X_s = pd.DataFrame(scaler.fit_transform(X.fillna(imputer)), columns=selected_features)

    model = LogisticRegression(max_iter=1000)
    model.fit(X_s, y)

    calibrator = IsotonicRegression(out_of_bounds="clip")
    cal_probs = model.predict_proba(X_s)[:, 1]
    calibrator.fit(cal_probs, y)

    return {
        "model": model,
        "calibrator": calibrator,
        "imputer": imputer,
        "scaler": scaler,
        "features": selected_features,
        "threshold": threshold,
    }


def _write_manifest(model_dir, model_version, model_type, selected_features, threshold, sha256=""):
    import yaml

    manifest = {
        "model_version": model_version,
        "model_name": f"Synthetic {model_type} model",
        "model_type": model_type,
        "training_dataset": "synthetic",
        "feature_list_file": "selected_features.txt",
        "selected_feature_count": len(selected_features),
        "threshold": threshold,
        "calibration": "isotonic_regression",
        "scaler": "robust_scaler",
        "imputation": "median",
        "alphagenome_version": "test",
        "gnomad_version": "test",
        "vep_version": "test",
        "training_date": "2026-06-26",
        "git_commit": "abc123",
        "package_version": "test",
        "pipeline_file": "final_pipeline.pkl",
        "pipeline_sha256": sha256,
    }
    with open(model_dir / "manifest.yaml", "w") as fh:
        yaml.dump(manifest, fh)


@pytest.fixture
def synthetic_model_dir(tmp_path: Path):
    """Create a temporary model root with synthetic full and no-AF models."""
    model_root = tmp_path / "model"
    model_root.mkdir()

    (model_root / "active_model.json").write_text(
        json.dumps({"default": "model_full", "no_af": "model_no_af"})
    )

    # Full model
    full_dir = model_root / "model_full"
    full_dir.mkdir()
    full_features = ["feat_a", "feat_b", "gnomAD_AF", "vep_IMPACT_score"]
    full_pipeline = _make_synthetic_pipeline(full_features, threshold=0.5)
    joblib.dump(full_pipeline, full_dir / "final_pipeline.pkl")
    with open(full_dir / "selected_features.txt", "w") as fh:
        fh.write("\n".join(full_features))
    _write_manifest(full_dir, "model_full", "full", full_features, 0.5)

    # No-AF model
    noaf_dir = model_root / "model_no_af"
    noaf_dir.mkdir()
    noaf_features = ["feat_a", "feat_b", "vep_IMPACT_score"]
    noaf_pipeline = _make_synthetic_pipeline(noaf_features, threshold=0.45)
    joblib.dump(noaf_pipeline, noaf_dir / "final_pipeline.pkl")
    with open(noaf_dir / "selected_features.txt", "w") as fh:
        fh.write("\n".join(noaf_features))
    _write_manifest(noaf_dir, "model_no_af", "no_AF", noaf_features, 0.45)

    yield model_root

    shutil.rmtree(model_root, ignore_errors=True)
