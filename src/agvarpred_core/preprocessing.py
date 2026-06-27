"""Preprocessing helpers: name cleaning and VEP feature encoding."""

from __future__ import annotations

import re
from typing import Iterable

import numpy as np
import pandas as pd

from .utils import clean_name


# Fixed ordinal mappings used in the training pipeline. Unknown categories map to -1.
# The integer order is least/most benign -> most damaging where applicable.
VEP_ORDINAL_MAPS = {
    "vep_IMPACT": {"MODIFIER": 0, "LOW": 1, "MODERATE": 2, "HIGH": 3},
    "vep_SIFT_pred": {
        "tolerated": 0,
        "tolerated_low_confidence": 1,
        "deleterious_low_confidence": 2,
        "deleterious": 3,
    },
    "vep_PolyPhen_pred": {
        "benign": 0,
        "possibly_damaging": 1,
        "probably_damaging": 2,
        "unknown": 3,
    },
    "vep_LoF": {"LC": 0, "HC": 1},
}


def clean_dataframe_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Return a copy of ``df`` with column names normalized via ``clean_name``."""
    df = df.copy()
    df.columns = [clean_name(c) for c in df.columns]
    return df


def encode_vep_features(
    df: pd.DataFrame,
    consequence_categories: Iterable[str] | None = None,
) -> pd.DataFrame:
    """Add derived VEP columns to a DataFrame containing raw VEP annotations.

    Input columns (created by :mod:`agvarpred_core.vep`) are expected:
        - vep_SIFT_score, vep_PolyPhen_score, vep_Protein_position (numeric)
        - vep_IMPACT_score (numeric)
        - vep_is_missense, vep_is_synonymous, vep_is_stop_gained,
          vep_is_frameshift, vep_is_splice, vep_is_LoF_HC,
          vep_has_SIFT, vep_has_PolyPhen (binary int)
        - vep_IMPACT, vep_SIFT_pred, vep_PolyPhen_pred, vep_LoF (categorical)
        - vep_Consequence (raw consequence string)

    The function adds / ensures correct dtypes and one-hot encodes the
    consequence categories. If ``consequence_categories`` is None, one-hot
    columns are generated for every unique non-null consequence observed in
    ``df``.
    """
    df = df.copy()

    # Numeric features
    for col in ["vep_SIFT_score", "vep_PolyPhen_score", "vep_Protein_position"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").astype("float32")

    if "vep_IMPACT_score" in df.columns:
        df["vep_IMPACT_score"] = (
            pd.to_numeric(df["vep_IMPACT_score"], errors="coerce").fillna(0).astype("int8")
        )

    # Binary flags
    binary_flags = [
        "vep_is_missense",
        "vep_is_synonymous",
        "vep_is_stop_gained",
        "vep_is_frameshift",
        "vep_is_splice",
        "vep_is_LoF_HC",
        "vep_has_SIFT",
        "vep_has_PolyPhen",
    ]
    for flag_col in binary_flags:
        if flag_col not in df.columns:
            if flag_col == "vep_has_SIFT" and "vep_SIFT_score" in df.columns:
                df[flag_col] = df["vep_SIFT_score"].notna().astype("int8")
            elif flag_col == "vep_has_PolyPhen" and "vep_PolyPhen_score" in df.columns:
                df[flag_col] = df["vep_PolyPhen_score"].notna().astype("int8")
            else:
                df[flag_col] = 0
        df[flag_col] = pd.to_numeric(df[flag_col], errors="coerce").fillna(0).astype("int8")

    # Ordinal categorical columns
    for col, mapping in VEP_ORDINAL_MAPS.items():
        if col in df.columns:
            df[col] = df[col].map(mapping).fillna(-1).astype("int16")

    # One-hot encode consequences
    if "vep_Consequence" in df.columns:
        if consequence_categories is None:
            consequence_categories = sorted(
                df["vep_Consequence"].dropna().unique()
            )
        for cat_val in consequence_categories:
            col_name = f"vep_Consequence_{cat_val}"
            df[col_name] = (df["vep_Consequence"] == cat_val).astype("int8")

    return df


def add_af_features(df: pd.DataFrame) -> pd.DataFrame:
    """Add gnomAD AF-derived features used by the model.

    Expects ``gnomAD_AF`` to be present (0.0 for missing).
    Adds:
        - log10_gnomAD_AF
        - AF_missing
        - is_ultra_rare
    """
    df = df.copy()
    if "gnomAD_AF" not in df.columns:
        df["gnomAD_AF"] = 0.0
    df["gnomAD_AF"] = pd.to_numeric(df["gnomAD_AF"], errors="coerce").fillna(0.0).astype("float32")
    df["log10_gnomAD_AF"] = np.log10(df["gnomAD_AF"] + 1e-8).astype("float32")
    df["AF_missing"] = (df["gnomAD_AF"] == 0).astype("int8")
    df["is_ultra_rare"] = (df["gnomAD_AF"] < 0.0001).astype("int8")
    return df


def align_columns(df: pd.DataFrame, required_columns: Iterable[str]) -> pd.DataFrame:
    """Ensure ``df`` contains exactly ``required_columns`` in the right order.

    Missing columns are filled with NaN. Extra columns are dropped.
    """
    required_columns = list(required_columns)
    for col in required_columns:
        if col not in df.columns:
            df[col] = np.nan
    return df[required_columns]
