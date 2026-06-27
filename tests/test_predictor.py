import numpy as np
import pandas as pd

from AGVarPred import AGVarPredAutoPredictor, AGVarPredPredictor
from agvarpred_core.af_source import NoAFSource


def test_predictor_with_synthetic_model(synthetic_model_dir):
    predictor = AGVarPredPredictor(model_name="model_full", model_dir=synthetic_model_dir)
    assert predictor.model_version == "model_full"
    assert predictor.model_type == "full"
    assert predictor.n_features == 4

    features = pd.DataFrame(
        {
            "feat_a": [0.5, -0.5, 1.0],
            "feat_b": [0.1, 0.2, -0.3],
            "gnomAD_AF": [0.0, 0.01, 0.5],
            "vep_IMPACT_score": [3, 2, 4],
        },
        index=pd.Index(["v1", "v2", "v3"], name="variant_id"),
    )
    preds = predictor.predict(features)
    assert len(preds) == 3
    assert set(preds.columns) == {
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
    assert preds["model_version"].iloc[0] == "model_full"
    assert preds["model_type"].iloc[0] == "full"
    assert set(preds["predicted_class"]).issubset({0, 1})


def test_predictor_empty_input(synthetic_model_dir):
    predictor = AGVarPredPredictor(model_name="model_full", model_dir=synthetic_model_dir)
    empty = pd.DataFrame(index=pd.Index([], name="variant_id"))
    preds = predictor.predict(empty)
    assert preds.empty


def test_auto_predictor_selects_no_af_without_gnomad(synthetic_model_dir):
    with np.errstate(invalid="ignore"):
        auto = AGVarPredAutoPredictor(
            model_dir=synthetic_model_dir,
            gnomad_vcf=None,
        )
    assert auto.model_name == "model_no_af"
    assert auto.model_type == "no_AF"
    assert auto.af_source_name == "none"

    features = pd.DataFrame(
        {
            "feat_a": [0.5],
            "feat_b": [0.1],
            "gnomAD_AF": [0.0],  # ignored by no-AF model
            "vep_IMPACT_score": [3],
        },
        index=pd.Index(["v1"], name="variant_id"),
    )
    preds = auto.predict(features)
    assert preds["model_type"].iloc[0] == "no_AF"
    assert preds["af_source"].iloc[0] == "none"


def test_auto_predictor_selects_full_with_local_gnomad(synthetic_model_dir, tmp_path):
    # Create a tiny fake gnomAD VCF so the local source is available
    import pysam

    vcf_path = tmp_path / "fake_gnomad.vcf.gz"
    header = pysam.VariantHeader()
    header.add_meta("fileformat", value="VCFv4.2")
    header.add_meta("INFO", items=[("ID", "AF"), ("Number", "A"), ("Type", "Float"), ("Description", "AF")])
    header.add_meta("INFO", items=[("ID", "vep"), ("Number", "."), ("Type", "String"), ("Description", "VEP")])
    header.contigs.add("chr1")
    with pysam.VariantFile(str(vcf_path), "wz", header=header) as vcf:
        rec = vcf.new_record()
        rec.chrom = "chr1"
        rec.pos = 100
        rec.ref = "A"
        rec.alts = ["G"]
        rec.info["AF"] = (0.01,)
        vcf.write(rec)

    auto = AGVarPredAutoPredictor(
        model_dir=synthetic_model_dir,
        gnomad_vcf=str(vcf_path),
    )
    assert auto.model_name == "model_full"
    assert auto.model_type == "full"
    assert auto.af_source_name == "local_gnomad"
