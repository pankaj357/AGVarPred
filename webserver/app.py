"""FastAPI entry point for the AGVarPred webserver prototype.

Model inference is intentionally gated until an authorized and scientifically
equivalent live feature-provider adapter has been implemented and validated.
"""

from __future__ import annotations

import os
from functools import lru_cache
from typing import Annotated, Literal

from fastapi import FastAPI, HTTPException, status
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, ConfigDict, Field, field_validator

from AGVarPred import AGVarPredPredictor

APP_NAME = "AGVarPred Web Server"
MAX_VARIANTS_PER_REQUEST = 25

app = FastAPI(
    title=APP_NAME,
    version="0.1.0",
    description=(
        "Prototype API for germline variant pathogenicity research. "
        "Not for clinical decision-making. Live prediction remains disabled "
        "until the live feature provider is implemented and validated."
    ),
)


class VariantInput(BaseModel):
    """One small variant on the GRCh38 reference assembly."""

    model_config = ConfigDict(extra="forbid")

    chrom: str = Field(
        description="Chromosome, e.g. chr1 or 1. GRCh38 coordinates are required."
    )
    pos: int = Field(gt=0, description="1-based reference position.")
    ref: str = Field(min_length=1, max_length=50)
    alt: str = Field(min_length=1, max_length=50)
    gene: str | None = Field(default=None, max_length=64)

    @field_validator("chrom")
    @classmethod
    def validate_chrom(cls, value: str) -> str:
        value = value.strip()
        if value.lower().startswith("chr"):
            value = value[3:]
        if value.upper() not in {
            *(str(i) for i in range(1, 23)),
            "X", "Y", "M", "MT",
        }:
            raise ValueError("chrom must be 1-22, X, Y, M, or MT")
        # Keep canonical casing; downstream AGVarPred code is chromosome-aware.
        return "chr" + value.upper().replace("MT", "M")

    @field_validator("ref", "alt")
    @classmethod
    def validate_allele(cls, value: str) -> str:
        value = value.strip().upper()
        if not value or any(base not in "ACGTN" for base in value):
            raise ValueError("ref and alt must contain only A, C, G, T, or N")
        if value == "." or value.startswith("<") or "[" in value or "]" in value:
            raise ValueError("symbolic alleles and breakends are not supported")
        return value

    @field_validator("gene")
    @classmethod
    def clean_gene(cls, value: str | None) -> str | None:
        if value is None:
            return None
        value = value.strip()
        return value or None


class PredictRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    genome_build: Literal["GRCh38"] = "GRCh38"
    variants: list[VariantInput] = Field(min_length=1, max_length=MAX_VARIANTS_PER_REQUEST)


class HealthResponse(BaseModel):
    status: str
    service: str
    version: str
    model_version: str
    model_type: str
    feature_count: int
    threshold: float
    live_prediction_enabled: bool
    live_provider_ready: bool


@lru_cache(maxsize=1)
def get_frozen_predictor() -> AGVarPredPredictor:
    """Load the exact bundled full model; no automatic no-AF fallback."""
    return AGVarPredPredictor(model_name="model_full", af_source_name="webserver")


@app.get("/api/health", response_model=HealthResponse, tags=["service"])
def health() -> HealthResponse:
    """Return service and frozen-model metadata without exposing secrets."""
    predictor = get_frozen_predictor()
    live_enabled = os.getenv("AGVARPRED_ENABLE_LIVE_PREDICTION", "false").lower() == "true"
    # This remains false until a real provider adapter passes integration tests.
    provider_ready = False
    return HealthResponse(
        status="ok",
        service=APP_NAME,
        version=app.version,
        model_version=predictor.model_version,
        model_type=predictor.model_type,
        feature_count=predictor.n_features,
        threshold=predictor.threshold,
        live_prediction_enabled=live_enabled and provider_ready,
        live_provider_ready=provider_ready,
    )


@app.post("/api/predict", tags=["prediction"])
def predict(request: PredictRequest) -> dict:
    """Validate a request and fail closed until live feature generation is ready.

    This route must never fabricate a score or substitute an incompatible
    annotation pipeline. The next development phase will implement the live
    AlphaGenome + VEP + gnomAD feature adapter, guarded by parity tests against
    the frozen research workflow.
    """
    # Force early model integrity/checksum verification before accepting work.
    get_frozen_predictor()

    live_enabled = os.getenv("AGVARPRED_ENABLE_LIVE_PREDICTION", "false").lower() == "true"
    if not live_enabled:
        raise HTTPException(
            status_code=status.HTTP_503_SERVICE_UNAVAILABLE,
            detail=(
                "Prediction is not enabled on this prototype. Live AlphaGenome, "
                "VEP, and gnomAD feature extraction must be configured and validated "
                "before scores can be returned."
            ),
        )

    raise HTTPException(
        status_code=status.HTTP_503_SERVICE_UNAVAILABLE,
        detail=(
            "Live prediction provider is not implemented yet. No predictions were "
            "generated; this endpoint will be enabled only after feature parity tests."
        ),
    )
