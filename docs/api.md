# API Reference

## `agvarpred_core.af_source`

Provider interface for AF/VEP annotations.

- `AFSource` — abstract base class.
- `LocalGnomADSource(vcf_path)` — query a local gnomAD VCF.
- `OnlineAFSource(config)` — stub for a future online provider.
- `NoAFSource()` — no annotation source (fallback).
- `resolve_af_source(gnomad_vcf=None, online_config=None)` — choose the best available provider.

## `agvarpred_core.feature_generator.FeatureGenerator`

```python
FeatureGenerator(
    af_source=None,
    alpha_mode="auto",
    alpha_api_key=None,
    alpha_dir=None,
)
```

- `from_vcf(path) -> pd.DataFrame`: generate the full engineered feature matrix from a VCF.
- `from_variants(df) -> pd.DataFrame`: generate features from a DataFrame with columns `chrom`, `pos`, `ref`, `alt`, `gene`.

## `agvarpred_core.feature_selector.FeatureSelector`

```python
FeatureSelector(selected_features)
```

- `select(df) -> pd.DataFrame`: return only the selected feature columns.

## `AGVarPred.AGVarPredPredictor`

Load an explicit model version and predict.

```python
AGVarPredPredictor(
    model_name="model_full",
    model_dir=None,
    af_source_name="unknown",
)
```

- `model_name` is required and must be a valid model directory name.
- `predict(features_df) -> pd.DataFrame`: run prediction with metadata.

## `AGVarPred.AGVarPredAutoPredictor`

Automatically select and load the best available model.

```python
AGVarPredAutoPredictor(
    model_dir=None,
    requested_model=None,  # "auto", "full", "no_af", or a directory name
    gnomad_vcf=None,
    online_config=None,
)
```

- `predict(features_df) -> pd.DataFrame`: run prediction using the selected model.
- `af_source`: resolved `AFSource` instance.
- `model_name`, `model_type`, `n_features`: properties of the selected model.

## `AGVarPred.cli.main`

Entry point for the `AGVarPred` console script.
