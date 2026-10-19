"""Live inference adapter that reuses the frozen AGVarPred package.

Live calls are gated by deployment configuration and should only be enabled
after the AlphaGenome provider terms have been reviewed for the intended
service. No provider credentials are stored in source control.
"""

from __future__ import annotations

import importlib.metadata
import importlib.util
import os
from functools import lru_cache
from pathlib import Path
from typing import Any

import pandas as pd

from AGVarPred.predictor import AGVarPredPredictor
from agvarpred_core.af_source import resolve_af_source
from agvarpred_core.feature_generator import FeatureGenerator
from agvarpred_core.utils import make_variant_id


ANNOTATION_ONLY_FEATURES = {
    "gnomAD_AF",
    "log10_gnomAD_AF",
    "AF_missing",
    "is_ultra_rare",
}


def env_enabled(name: str) -> bool:
    """Return true only for an explicit truthy deployment setting."""
    return os.getenv(name, "").strip().lower() in {"1", "true", "yes", "on"}


def alpha_feature_names(selected_features: list[str]) -> list[str]:
    """Return frozen-model AlphaGenome features, excluding VEP/AF features."""
    return [
        name
        for name in selected_features
        if not name.startswith("vep_") and name not in ANNOTATION_ONLY_FEATURES
    ]


def apply_training_missing_annotation_defaults(
    features: pd.DataFrame, selected_features: list[str]
) -> pd.DataFrame:
    """Apply the same VEP-missing defaults used when preparing training data.

    Missing gnomAD VEP annotations are a supported input pattern in the training
    pipeline. Keep those values batch-independent instead of letting an entirely
    unannotated batch turn selected VEP columns into median-imputed values.
    Numeric SIFT/PolyPhen scores remain NaN and are imputed by the frozen model.
    """
    features = features.copy()
    ordinal_missing = {"vep_LoF", "vep_SIFT_pred", "vep_PolyPhen_pred", "vep_IMPACT"}
    numeric_missing = {"vep_SIFT_score", "vep_PolyPhen_score", "vep_Protein_position"}

    for name in selected_features:
        if not name.startswith("vep_"):
            continue

        if name.startswith("vep_Consequence_"):
            if name not in features.columns:
                features[name] = 0
            else:
                features[name] = features[name].fillna(0)
        elif name == "vep_IMPACT_score" or name.startswith("vep_is_") or name.startswith("vep_has_"):
            if name not in features.columns:
                features[name] = 0
            else:
                features[name] = features[name].fillna(0)
        elif name in ordinal_missing:
            if name not in features.columns:
                features[name] = -1
            else:
                features[name] = pd.to_numeric(features[name], errors="coerce").fillna(-1)
        elif name in numeric_missing:
            if name not in features.columns:
                features[name] = float("nan")
            else:
                features[name] = pd.to_numeric(features[name], errors="coerce")
        else:
            raise RuntimeError(
                "No training-compatible missing-value rule is defined for selected "
                f"VEP feature {name!r}; refusing to score."
            )

    return features


def missing_model_feature_columns(
    features: pd.DataFrame, selected_features: list[str]
) -> list[str]:
    """List required selected columns that are absent after normalization."""
    return [name for name in selected_features if name not in features.columns]

def live_provider_status(
    expected_alphagenome_version: str | None = None,
) -> tuple[bool, str]:
    """Check non-secret prerequisites for live API-based feature extraction.

    Refuse to serve with a different AlphaGenome SDK version than the one
    recorded for the frozen model until feature parity is re-validated.
    """
    api_key = os.getenv("ALPHAGENOME_API_KEY")
    if not api_key:
        return False, "AlphaGenome API key is not configured."

    gnomad_tsv = os.getenv("GNOMAD_TSV")
    gnomad_vcf = os.getenv("GNOMAD_VCF")
    if gnomad_tsv and Path(gnomad_tsv).is_file():
        if not Path(f"{gnomad_tsv}.tbi").is_file():
            return False, "The configured GNOMAD_TSV index is missing."
    elif gnomad_vcf and Path(gnomad_vcf).is_file():
        index_candidates = (Path(f"{gnomad_vcf}.tbi"), Path(f"{gnomad_vcf}.csi"))
        if not any(path.is_file() for path in index_candidates):
            return False, "The configured GNOMAD_VCF index is missing."
    else:
        return False, "Neither GNOMAD_TSV nor GNOMAD_VCF is configured."

    try:
        sdk_available = importlib.util.find_spec("alphagenome") is not None
    except (ImportError, ValueError):
        sdk_available = False
    if not sdk_available:
        return False, "The AlphaGenome Python SDK is not installed."

    try:
        installed_version = importlib.metadata.version("alphagenome")
    except importlib.metadata.PackageNotFoundError:
        return False, "AlphaGenome SDK package metadata is unavailable."

    if expected_alphagenome_version and installed_version != str(expected_alphagenome_version):
        return (
            False,
            "Installed AlphaGenome SDK version "
            f"{installed_version} differs from the model's recorded feature version "
            f"{expected_alphagenome_version}; validate feature parity before serving.",
        )

    return True, "ready"


class LivePredictionService:
    """Generate live AlphaGenome + gnomAD/VEP features and score with frozen model."""

    def __init__(self, predictor: AGVarPredPredictor):
        ready, message = live_provider_status(
            predictor.manifest.get("alphagenome_version")
        )
        if not ready:
            raise RuntimeError(message)

        api_key = os.environ["ALPHAGENOME_API_KEY"]

        # Prefer the reduced tabix-indexed TSV (GNOMAD_TSV) if available, but
        # fall back to the full local VCF (GNOMAD_VCF) for backward compatibility.
        af_source = resolve_af_source()
        if not af_source.is_available():
            raise RuntimeError("Neither GNOMAD_TSV nor GNOMAD_VCF is configured.")

        self.predictor = predictor
        if predictor.model_type != "full" or predictor.n_features != 120:
            raise RuntimeError("The frozen full AGVarPred model is not the expected version.")

        self.feature_generator = FeatureGenerator(
            af_source=af_source,
            alpha_mode="sdk",
            alpha_api_key=api_key,
        )
        self.gnomad_release = self.predictor.manifest.get("gnomad_version", "")
        self.alphagenome_release = self.predictor.manifest.get("alphagenome_version", "")

    def predict(self, variants: list[Any]) -> dict[str, Any]:
        """Score input variants; return explicit per-variant failures if API scoring fails."""
        rows = []
        original_by_id: dict[str, Any] = {}
        for variant in variants:
            variant_id = make_variant_id(
                variant.chrom, variant.pos, variant.ref, variant.alt
            )
            original_by_id[variant_id] = variant
            rows.append(
                {
                    "chrom": variant.chrom,
                    "pos": int(variant.pos),
                    "ref": variant.ref.upper(),
                    "alt": variant.alt.upper(),
                    "gene": variant.gene or "",
                    "variant_id": variant_id,
                }
            )

        variants_df = pd.DataFrame(rows)
        features = self.feature_generator.from_variants(variants_df)

        selected = list(self.predictor.selected_features)
        alpha_features = alpha_feature_names(selected)
        if not alpha_features:
            raise RuntimeError("The frozen model's AlphaGenome feature schema is unavailable.")

        splice_prefixes = (
            "SPLICE_SITES_",
            "SPLICE_SITE_USAGE_",
            "SPLICE_JUNCTIONS_",
        )
        
        missing_alpha = [
            name for name in alpha_features
            if name not in features.columns
        ]
        
        blocking_missing = [
            name for name in missing_alpha
            if not name.startswith(splice_prefixes)
        ]
        
        # If a splice modality has returned some columns but selected tracks
        # are missing, treat this as a schema mismatch rather than hiding it.
        for prefix in splice_prefixes:
            missing_for_modality = [
                name for name in missing_alpha
                if name.startswith(prefix)
            ]
        
            modality_present = any(
                str(col).startswith(prefix)
                for col in features.columns
            )
        
            if modality_present:
                blocking_missing.extend(missing_for_modality)
        
        if blocking_missing:
            examples = ", ".join(blocking_missing[:8])
            more = " ..." if len(blocking_missing) > 8 else ""
        
            raise RuntimeError(
                "AlphaGenome feature-schema mismatch: "
                f"{len(blocking_missing)} selected tracks are absent "
                f"(e.g. {examples}{more}). No predictions were generated."
            )
        
        # If an entire splice modality is absent, preserve its missingness.
        # The frozen predictor's FeatureSelector fills missing selected columns
        # with NaN, and the trained imputer handles them during inference.
        for name in missing_alpha:
            features[name] = float("nan")

        # Training permits absent gnomAD/VEP annotations. Apply its exact missing
        # value policy to selected VEP inputs so the result does not depend on
        # whether this request contains another variant with a VEP record.
        features = apply_training_missing_annotation_defaults(features, selected)
        missing_required = missing_model_feature_columns(features, selected)
        if missing_required:
            examples = ", ".join(missing_required[:8])
            more = " ..." if len(missing_required) > 8 else ""
            raise RuntimeError(
                "Frozen-model feature-schema mismatch: "
                f"{len(missing_required)} required columns are absent "
                f"(e.g. {examples}{more}). No predictions were generated."
            )

        successful_ids: list[str] = []
        row_errors: dict[str, str] = {}
        for variant_id in original_by_id:
            if variant_id not in features.index:
                row_errors[variant_id] = (
                    "No usable AlphaGenome scores were returned for this variant."
                )
                continue

            row = features.loc[variant_id, alpha_features]
            if not row.notna().any():
                row_errors[variant_id] = (
                    "No usable AlphaGenome scores were returned for this variant."
                )
                continue
            successful_ids.append(variant_id)

        prediction_map: dict[str, dict[str, Any]] = {}
        if successful_ids:
            model_input = features.loc[successful_ids].copy()
            predictions = self.predictor.predict(model_input)
            for record in predictions.to_dict(orient="records"):
                variant_id = str(record["variant_id"])
                feature_row = features.loc[variant_id]
                variant = original_by_id[variant_id]
                predicted_class = int(record["predicted_class"])
                af_value = feature_row.get("gnomAD_AF")
                consequence = feature_row.get("vep_Consequence")
                prediction_map[variant_id] = {
                    "variant_id": variant_id,
                    "chrom": variant.chrom,
                    "pos": int(variant.pos),
                    "ref": variant.ref,
                    "alt": variant.alt,
                    "gene": variant.gene,
                    "status": "ok",
                    "score": float(record["probability"]),
                    "predicted_class": (
                        "Pathogenic" if predicted_class == 1 else "Benign"
                    ),
                    "threshold": float(self.predictor.threshold),
                    "gnomAD_AF": (
                        float(af_value)
                        if pd.notna(af_value)
                        else None
                    ),
                    "consequence": (
                        str(consequence)
                        if isinstance(consequence, str) and consequence
                        else None
                    ),
                    "model_version": self.predictor.model_version,
                    "model_type": self.predictor.model_type,
                    "gnomAD_release": self.gnomad_release,
                    "AlphaGenome_version": self.alphagenome_release,
                }

        results = []
        for variant_id, variant in original_by_id.items():
            if variant_id in prediction_map:
                results.append(prediction_map[variant_id])
            else:
                results.append(
                    {
                        "variant_id": variant_id,
                        "chrom": variant.chrom,
                        "pos": int(variant.pos),
                        "ref": variant.ref,
                        "alt": variant.alt,
                        "gene": variant.gene,
                        "status": "error",
                        "error": row_errors.get(
                            variant_id, "No prediction was generated for this variant."
                        ),
                    }
                )

        return {
            "model_version": self.predictor.model_version,
            "model_type": self.predictor.model_type,
            "threshold": float(self.predictor.threshold),
            "requested_variants": len(variants),
            "scored_variants": len(prediction_map),
            "results": results,
            "notice": (
                "Research use only. Scores are calibrated to the ClinVar "
                "development distribution and are not clinical diagnoses."
            ),
        }


@lru_cache(maxsize=1)
def get_live_service(predictor: AGVarPredPredictor) -> LivePredictionService:
    """Build the API-backed feature service once per worker process."""
    return LivePredictionService(predictor)
