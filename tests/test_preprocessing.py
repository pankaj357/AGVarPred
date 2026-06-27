import numpy as np
import pandas as pd

from agvarpred_core.preprocessing import (
    add_af_features,
    clean_dataframe_columns,
    encode_vep_features,
)
from agvarpred_core.utils import clean_name


def test_clean_name_matches_training():
    # The original training pipeline normalizes names as follows:
    # replace colons/hyphens/spaces with _, remove other non-alnum/underscore,
    # then collapse repeated underscores.
    assert clean_name("ATAC__A549__EFO:0001086 ATAC-seq") == "ATAC_A549_EFO_0001086_ATAC_seq"
    assert clean_name("SPLICE_SITE_USAGE_testis_usage_UBERON:0000473_total_RNA-seq") == (
        "SPLICE_SITE_USAGE_testis_usage_UBERON_0000473_total_RNA_seq"
    )


def test_clean_dataframe_columns():
    df = pd.DataFrame({"A: B-C": [1], "x  y": [2]})
    df = clean_dataframe_columns(df)
    assert list(df.columns) == ["A_B_C", "x_y"]


def test_encode_vep_features():
    df = pd.DataFrame(
        {
            "vep_SIFT_score": [0.1, None],
            "vep_PolyPhen_score": [0.99, None],
            "vep_Protein_position": [123, None],
            "vep_IMPACT_score": [3, 0],
            "vep_IMPACT": ["MODERATE", "MODIFIER"],
            "vep_SIFT_pred": ["deleterious", None],
            "vep_PolyPhen_pred": ["probably_damaging", None],
            "vep_LoF": ["HC", None],
            "vep_Consequence": ["missense_variant", "synonymous_variant"],
            "vep_is_missense": [1, 0],
            "vep_is_synonymous": [0, 1],
        }
    )
    out = encode_vep_features(df)
    assert out["vep_IMPACT_score"].dtype == np.int8
    assert out["vep_LoF"].iloc[0] == 1
    assert out["vep_LoF"].iloc[1] == -1
    assert "vep_Consequence_missense_variant" in out.columns
    assert "vep_Consequence_synonymous_variant" in out.columns


def test_add_af_features():
    df = pd.DataFrame({"gnomAD_AF": [0.0, 0.00001, 0.5]})
    out = add_af_features(df)
    assert out["AF_missing"].tolist() == [1, 0, 0]
    assert out["is_ultra_rare"].tolist() == [1, 1, 0]
    assert np.allclose(out["log10_gnomAD_AF"].iloc[0], np.log10(1e-8))
