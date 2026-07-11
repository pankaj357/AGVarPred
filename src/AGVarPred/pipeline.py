"""Model-discovery and manifest utilities for AGVarPred."""

from __future__ import annotations

import hashlib
import os
from importlib import resources
from pathlib import Path
from typing import Any

import yaml


def _as_resource(model_root: str | Path | resources.abc.Traversable) -> Path | resources.abc.Traversable:
    """Return a path-like or package resource for ``model_root``."""
    if isinstance(model_root, (str, Path)):
        return Path(model_root)
    return model_root


def get_model_root(model_dir: str | Path | None = None) -> Path | resources.abc.Traversable:
    """Resolve the root directory containing model versions.

    Resolution order:
        1. ``model_dir`` argument.
        2. ``AGVARPRED_MODEL_DIR`` environment variable.
        3. Bundled ``model/`` directory shipped as package data.
    """
    if model_dir is not None:
        return Path(model_dir).resolve()

    env_path = os.environ.get("AGVARPRED_MODEL_DIR")
    if env_path:
        return Path(env_path).resolve()

    # Installed-package fallback: models are shipped as package data.
    return resources.files("AGVarPred") / "model"


def load_active_model_map(model_root: str | Path | resources.abc.Traversable) -> dict[str, str]:
    """Load ``active_model.json`` and return the full mapping."""
    model_root = _as_resource(model_root)
    active_path = model_root / "active_model.json"
    if active_path.exists():
        with active_path.open("r") as fh:
            data = yaml.safe_load(fh)
        return {
            "default": data.get("default"),
            "full": data.get("full") or data.get("default"),
            "no_af": data.get("no_af"),
        }

    # Fallback: infer from directories
    versions = [
        d.name
        for d in model_root.iterdir()
        if d.is_dir() and (d / "manifest.yaml").exists()
    ]
    return {
        "default": versions[0] if versions else None,
        "full": versions[0] if versions else None,
        "no_af": None,
    }


def load_active_model(model_root: str | Path | resources.abc.Traversable) -> str:
    """Return the default model version from ``active_model.json``."""
    mapping = load_active_model_map(model_root)
    default = mapping.get("default")
    if not default:
        raise FileNotFoundError(
            f"No active_model.json found in {model_root} and could not infer a default model."
        )
    return default


def resolve_model_name(
    model_root: str | Path | resources.abc.Traversable,
    requested: str | None,
    af_source_name: str,
) -> str:
    """Resolve which model version to load.

    Parameters
    ----------
    requested:
        Explicit model request: ``full``, ``no_af``, or a directory name.
    af_source_name:
        Name of the resolved AF source (``local_gnomad``, ``online``, ``none``).
    """
    model_root = _as_resource(model_root)
    mapping = load_active_model_map(model_root)

    if requested:
        requested = requested.lower()
        if requested == "full":
            name = mapping.get("full") or mapping.get("default")
            if not name:
                raise ValueError("No full model configured in active_model.json")
            return name
        if requested in ("no_af", "no-af", "noaf"):
            name = mapping.get("no_af")
            if not name:
                raise ValueError("No no-AF model configured in active_model.json")
            return name
        # Assume it is a concrete directory name
        return requested

    # Automatic selection
    if af_source_name in ("local_gnomad", "online"):
        return mapping.get("full") or mapping.get("default")

    return mapping.get("no_af") or mapping.get("default")


def load_manifest(
    model_root: str | Path | resources.abc.Traversable,
    model_name: str,
) -> dict[str, Any]:
    """Load ``manifest.yaml`` for a specific model version."""
    model_root = _as_resource(model_root)
    manifest_path = model_root / model_name / "manifest.yaml"
    if not manifest_path.exists():
        raise FileNotFoundError(f"Model manifest not found: {manifest_path}")
    with manifest_path.open("r") as fh:
        return yaml.safe_load(fh)


def _sha256_of_resource(resource: resources.abc.Traversable) -> str:
    """Compute the SHA256 hex digest of a package resource or file."""
    h = hashlib.sha256()
    with resource.open("rb") as fh:
        for chunk in iter(lambda: fh.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def validate_manifest(
    model_root: str | Path | resources.abc.Traversable,
    model_name: str,
) -> dict[str, Any]:
    """Load manifest and verify the pipeline file checksum."""
    model_root = _as_resource(model_root)
    manifest = load_manifest(model_root, model_name)
    pipeline_file = model_root / model_name / manifest["pipeline_file"]
    expected = manifest.get("pipeline_sha256")
    if expected:
        actual = _sha256_of_resource(pipeline_file)
        if actual != expected:
            raise ValueError(
                f"SHA256 mismatch for {pipeline_file}: expected {expected}, got {actual}"
            )
    return manifest


def list_models(model_root: str | Path | resources.abc.Traversable | None = None) -> list[str]:
    """Return the names of available model versions."""
    model_root = get_model_root(model_root) if model_root is None else _as_resource(model_root)
    return sorted(
        d.name
        for d in model_root.iterdir()
        if d.is_dir() and (d / "manifest.yaml").exists()
    )
