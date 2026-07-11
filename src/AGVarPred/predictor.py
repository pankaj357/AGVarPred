"""AGVarPred predictor: load a frozen model and produce predictions."""

from __future__ import annotations

import warnings
from datetime import datetime, timezone
from importlib.resources import as_file
from pathlib import Path
from typing import Any

import joblib
import numpy as np
import pandas as pd

from agvarpred_core.feature_selector import FeatureSelector

from .pipeline import get_model_root, validate_manifest


class AGVarPredPredictor:
    """Load a frozen AGVarPred model and run predictions.

    Parameters
    ----------
    model_name:
        Model version directory name (e.g. ``model_full``). Required.
    model_dir:
        Root directory containing model versions. If None, resolved via
        ``AGVARPRED_MODEL_DIR`` env var or the editable-install repo root.
    af_source_name:
        Name of the AF/VEP source used for this prediction, stored in output
        metadata (e.g. ``local_gnomad``, ``online``, ``none``).
    """

    def __init__(
        self,
        model_name: str,
        model_dir: str | Path | None = None,
        af_source_name: str = "unknown",
    ):
        if not model_name:
            raise ValueError("model_name is required")
        self.model_root = get_model_root(model_dir)
        self.model_name = model_name
        self.af_source_name = af_source_name
        self.manifest = validate_manifest(self.model_root, self.model_name)

        pipeline_path = self.model_root / self.model_name / self.manifest["pipeline_file"]
        feature_path = self.model_root / self.model_name / self.manifest["feature_list_file"]
        with (
            as_file(pipeline_path) as pipeline_file,
            as_file(feature_path) as feature_file,
        ):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=UserWarning)
                self.pipeline = joblib.load(pipeline_file)

            self.feature_selector = FeatureSelector(feature_file)
        self.threshold = float(self.manifest["threshold"])
        self.selected_features = list(self.pipeline.get("features", self.feature_selector.selected_features))
        self.model_type = self.manifest.get("model_type", "unknown")

        # Sanity check: feature list in pipeline vs manifest should match
        if set(self.selected_features) != set(self.feature_selector.selected_features):
            warnings.warn(
                "Feature list in pipeline differs from manifest feature list. "
                "Using the pipeline's feature list."
            )
            self.feature_selector = FeatureSelector(self.selected_features)

    def predict(self, features_df: pd.DataFrame) -> pd.DataFrame:
        """Predict pathogenicity from a full or selected feature DataFrame.

        Parameters
        ----------
        features_df:
            DataFrame indexed by ``variant_id`` containing the engineered
            features. It may be the full matrix or already restricted to the
            selected features.

        Returns
        -------
        DataFrame with columns:
            - variant_id
            - probability (calibrated)
            - predicted_class (0/1)
            - model_version, model_type, af_source, feature_version,
              alphagenome_version, gnomAD_version, VEP_version,
              prediction_timestamp
        """
        if features_df.empty:
            return pd.DataFrame(
                columns=[
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
                ]
            )

        # Select model features
        X = self.feature_selector.select(features_df)
        variant_ids = X.index.tolist()

        # Preprocessing: median imputation + robust scaling
        imputer = self.pipeline["imputer"]
        scaler = self.pipeline["scaler"]
        model = self.pipeline["model"]
        calibrator = self.pipeline["calibrator"]

        X_imp = X.fillna(imputer)
        X_s = pd.DataFrame(scaler.transform(X_imp), index=X_imp.index, columns=X_imp.columns)

        # Raw model probability
        raw_prob = model.predict_proba(X_s)[:, 1]

        # Calibration and thresholding
        cal_prob = calibrator.transform(raw_prob)
        predicted_class = (cal_prob > self.threshold).astype(int)

        timestamp = datetime.now(timezone.utc).isoformat()
        result = pd.DataFrame(
            {
                "variant_id": variant_ids,
                "probability": cal_prob,
                "predicted_class": predicted_class,
                "model_version": self.manifest["model_version"],
                "model_type": self.model_type,
                "af_source": self.af_source_name,
                "feature_version": self.manifest.get("package_version", ""),
                "alphagenome_version": self.manifest.get("alphagenome_version", ""),
                "gnomAD_version": self.manifest.get("gnomad_version", ""),
                "VEP_version": self.manifest.get("vep_version", ""),
                "prediction_timestamp": timestamp,
            }
        )
        return result

    @property
    def model_version(self) -> str:
        return str(self.manifest["model_version"])

    @property
    def n_features(self) -> int:
        return len(self.selected_features)


class AGVarPredAutoPredictor:
    """Automatically choose the best available model and run predictions.

    Resolution order:
        1. If a local gnomAD VCF is available, use the full model.
        2. If an online AF/VEP provider is enabled, use the full model.
        3. Otherwise, fall back to the bundled no-AF model with a warning.

    Parameters
    ----------
    model_dir:
        Root directory containing model versions.
    requested_model:
        ``full``, ``no_af``, or a concrete model directory name. If None,
        automatic selection is used.
    gnomad_vcf:
        Explicit path to a local gnomAD VCF. Falls back to ``GNOMAD_VCF``.
    online_config:
        Configuration for a future online provider.
    """

    def __init__(
        self,
        model_dir: str | Path | None = None,
        requested_model: str | None = None,
        gnomad_vcf: str | Path | None = None,
        online_config: dict[str, Any] | None = None,
    ):
        from agvarpred_core.af_source import resolve_af_source

        self.af_source = resolve_af_source(gnomad_vcf, online_config)
        self.af_source_name = self.af_source.name

        model_root = get_model_root(model_dir)
        from .pipeline import resolve_model_name

        self.model_name = resolve_model_name(
            model_root, requested_model, self.af_source_name
        )

        if requested_model is None and self.af_source_name == "none":
            warnings.warn(
                "Population allele-frequency annotations were unavailable. "
                f"Using the bundled {self.model_name} model instead. "
                "Prediction accuracy may differ slightly from the primary model.",
                UserWarning,
                stacklevel=2,
            )

        self.predictor = AGVarPredPredictor(
            model_name=self.model_name,
            model_dir=model_dir,
            af_source_name=self.af_source_name,
        )

    def predict(self, features_df: pd.DataFrame) -> pd.DataFrame:
        """Run prediction using the automatically selected model."""
        return self.predictor.predict(features_df)

    @property
    def model_version(self) -> str:
        return self.predictor.model_version

    @property
    def model_type(self) -> str:
        return self.predictor.model_type

    @property
    def n_features(self) -> int:
        return self.predictor.n_features
