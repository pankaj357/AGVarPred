#!/usr/bin/env python3
"""Preprocess PanelApp benign VEP parquet to match training format."""

import pandas as pd
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
IN_PATH = ROOT / "external_validation" / "panelapp_benign_vep.parquet"
OUTDIR = ROOT / "external_validation" / "vep_preprocessed"
OUTDIR.mkdir(exist_ok=True, parents=True)
OUT_PATH = OUTDIR / "panelapp_benign_vep.parquet"
VEP_LINK = ROOT / "external_validation" / "vep" / "panelapp_benign_vep.parquet"


def preprocess_vep(df):
    df = df.copy()
    if "variant_id" in df.columns:
        df = df.set_index("variant_id")

    # Numeric features
    for col in ["vep_SIFT_score", "vep_PolyPhen_score", "vep_Protein_position"]:
        if col in df.columns:
            df[col] = df[col].astype("float32")

    # vep_IMPACT_score
    if "vep_IMPACT_score" in df.columns:
        df["vep_IMPACT_score"] = df["vep_IMPACT_score"].fillna(0).astype("int8")

    # Binary flags
    for flag_col in ["vep_is_missense", "vep_is_synonymous", "vep_is_stop_gained",
                     "vep_is_frameshift", "vep_is_splice", "vep_is_LoF_HC",
                     "vep_has_SIFT", "vep_has_PolyPhen"]:
        if flag_col not in df.columns:
            if flag_col == "vep_has_SIFT" and "vep_SIFT_score" in df.columns:
                df[flag_col] = df["vep_SIFT_score"].notna().astype("int8")
            elif flag_col == "vep_has_PolyPhen" and "vep_PolyPhen_score" in df.columns:
                df[flag_col] = df["vep_PolyPhen_score"].notna().astype("int8")
            else:
                df[flag_col] = 0
        df[flag_col] = df[flag_col].fillna(0).astype("int8")

    # vep_LoF ordinal
    if "vep_LoF" in df.columns:
        lof_map = {"LC": 0, "HC": 1}
        df["vep_LoF"] = df["vep_LoF"].map(lof_map).fillna(-1).astype("int16")

    # vep_IMPACT ordinal
    if "vep_IMPACT" in df.columns:
        impact_map = {"MODIFIER": 0, "LOW": 1, "MODERATE": 2, "HIGH": 3}
        df["vep_IMPACT"] = df["vep_IMPACT"].map(impact_map).fillna(-1).astype("int16")

    # vep_SIFT_pred ordinal
    if "vep_SIFT_pred" in df.columns:
        sift_map = {"tolerated": 0, "tolerated_low_confidence": 1, "deleterious_low_confidence": 2, "deleterious": 3}
        df["vep_SIFT_pred"] = df["vep_SIFT_pred"].map(sift_map).fillna(-1).astype("int16")

    # vep_PolyPhen_pred ordinal
    if "vep_PolyPhen_pred" in df.columns:
        pp_map = {"benign": 0, "possibly_damaging": 1, "probably_damaging": 2, "unknown": 3}
        df["vep_PolyPhen_pred"] = df["vep_PolyPhen_pred"].map(pp_map).fillna(-1).astype("int16")

    # One-hot consequences
    if "vep_Consequence" in df.columns:
        for cat_val in df["vep_Consequence"].dropna().unique():
            col_name = f"vep_Consequence_{cat_val}"
            df[col_name] = (df["vep_Consequence"] == cat_val).astype("int8")

    return df


def main():
    if not IN_PATH.exists():
        print(f"⚠ Input VEP not found: {IN_PATH}")
        return

    print(f"🔄 Preprocessing {IN_PATH.name}...")
    df = pd.read_parquet(IN_PATH)
    print(f"  Raw columns: {len(df.columns)}")

    df = preprocess_vep(df)
    print(f"  After preprocessing: {len(df.columns)}")

    df.to_parquet(OUT_PATH)
    print(f"  Saved: {OUT_PATH}")

    # Create symlink in external_validation/vep/
    VEP_LINK.parent.mkdir(exist_ok=True, parents=True)
    if VEP_LINK.exists() or VEP_LINK.is_symlink():
        VEP_LINK.unlink()
    VEP_LINK.symlink_to(Path("..") / "vep_preprocessed" / "panelapp_benign_vep.parquet")
    print(f"  Symlinked: {VEP_LINK} -> {OUT_PATH}")


if __name__ == "__main__":
    main()
