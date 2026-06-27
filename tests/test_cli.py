import os
from pathlib import Path

import pandas as pd
import pytest

from AGVarPred.cli import main


def test_cli_list_models():
    assert main(["list-models"]) == 0


def _make_minimal_vcf_alpha(tmp_path, gene="GENE1"):
    import pysam

    vcf_path = tmp_path / "input.vcf"
    header = pysam.VariantHeader()
    header.add_meta("fileformat", value="VCFv4.2")
    header.add_meta(
        "INFO",
        items=[("ID", "GENE"), ("Number", "."), ("Type", "String"), ("Description", "Gene")],
    )
    header.contigs.add("chr1")
    with pysam.VariantFile(str(vcf_path), "w", header=header) as vcf:
        rec = vcf.new_record()
        rec.chrom = "chr1"
        rec.pos = 100
        rec.ref = "A"
        rec.alts = ["G"]
        rec.info["GENE"] = gene
        vcf.write(rec)

    alpha_dir = tmp_path / "alpha"
    gene_dir = alpha_dir / f"{gene}_VEP"
    gene_dir.mkdir(parents=True)
    feats = pd.DataFrame(
        {
            "feat_a": [0.5],
            "feat_b": [-0.2],
            "gnomAD_AF": [0.0],
            "vep_IMPACT_score": [3],
        },
        index=pd.Index(["chr1_100_A_G"], name="variant_id"),
    )
    feats.to_parquet(gene_dir / f"{gene}_ALL_VEP_RAW_SCORE_MATRIX.parquet")
    return vcf_path, alpha_dir


def test_cli_predict_auto_no_af_fallback(tmp_path, synthetic_model_dir):
    vcf_path, alpha_dir = _make_minimal_vcf_alpha(tmp_path)
    out_path = tmp_path / "pred.csv"

    rc = main(
        [
            "predict",
            str(vcf_path),
            "--output",
            str(out_path),
            "--model-dir",
            str(synthetic_model_dir),
            "--alpha-mode",
            "precomputed",
            "--alpha-dir",
            str(alpha_dir),
        ]
    )
    assert rc == 0
    preds = pd.read_csv(out_path)
    assert len(preds) == 1
    assert preds["model_type"].iloc[0] == "no_AF"
    assert preds["af_source"].iloc[0] == "none"


def test_cli_predict_force_full_without_af(tmp_path, synthetic_model_dir):
    """Advanced users can force the full model even without AF; metadata records af_source=none."""
    vcf_path, alpha_dir = _make_minimal_vcf_alpha(tmp_path)
    out_path = tmp_path / "pred.csv"

    rc = main(
        [
            "predict",
            str(vcf_path),
            "--output",
            str(out_path),
            "--model-dir",
            str(synthetic_model_dir),
            "--model",
            "full",
            "--alpha-mode",
            "precomputed",
            "--alpha-dir",
            str(alpha_dir),
        ]
    )
    assert rc == 0
    preds = pd.read_csv(out_path)
    assert preds["model_type"].iloc[0] == "full"
    assert preds["af_source"].iloc[0] == "none"
