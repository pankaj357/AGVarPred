from pathlib import Path

from AGVarPred.pipeline import list_models, load_active_model, load_active_model_map, validate_manifest
from agvarpred_core.utils import sha256_file


def _model_root() -> Path:
    """Return the bundled model directory path."""
    import AGVarPred
    return Path(AGVarPred.__file__).resolve().parent / "model"


def test_real_model_manifest():
    model_root = _model_root()
    assert model_root.exists()

    mapping = load_active_model_map(model_root)
    assert mapping["default"] == "model_full"
    assert mapping["no_af"] == "model_no_af"

    manifest = validate_manifest(model_root, "model_full")
    required_keys = {
        "model_version",
        "model_type",
        "training_dataset",
        "feature_list_file",
        "selected_feature_count",
        "threshold",
        "pipeline_file",
        "pipeline_sha256",
    }
    assert required_keys.issubset(manifest.keys())
    assert manifest["selected_feature_count"] == 120
    assert manifest["model_type"] == "full"

    pipeline_path = model_root / "model_full" / manifest["pipeline_file"]
    actual_sha = sha256_file(pipeline_path)
    assert actual_sha == manifest["pipeline_sha256"]


def test_no_af_model_manifest():
    model_root = _model_root()
    manifest = validate_manifest(model_root, "model_no_af")
    assert manifest["model_type"] == "no_AF"
    assert manifest["selected_feature_count"] == 119


def test_list_models():
    model_root = _model_root()
    models = list_models(model_root)
    assert "model_full" in models
    assert "model_no_af" in models
