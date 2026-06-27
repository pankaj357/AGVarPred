# AGVarPred Zenodo Archive v1.0.0

This archive contains the frozen reproducibility artifacts for the AGVarPred
manuscript. It is designed to be cited via its Zenodo DOI.

## Contents

```
AGVarPred-zenodo/
├── model/
│   ├── active_model.json
│   ├── model_full/                       # Primary production model (with AF)
│   │   ├── final_pipeline.pkl
│   │   ├── selected_features.txt         # 120 selected features
│   │   ├── selected_features_with_importance.csv
│   │   └── manifest.yaml
│   └── model_no_af/                      # Fallback model (without AF)
│       ├── final_pipeline.pkl
│       ├── selected_features.txt         # 119 selected features
│       └── manifest.yaml
├── splits/
│   ├── train.csv                         # Train/cal/test variant definitions
│   ├── cal.csv
│   ├── test.csv
│   └── train_gene_set.json
├── benchmarks/                           # External validation benchmark definitions
├── metrics/                              # Final model performance metrics
├── data_manifest.json                    # External data download instructions
├── checksums.sha256                      # SHA256 checksums of archive files
├── code_reference.txt                    # Link to the GitHub release
├── CITATION.cff
└── LICENSE
```

## Models

- **model_full** uses 120 features including gnomAD allele frequency. It is the
  primary model used when a local gnomAD VCF or online AF provider is available.
- **model_no_af** uses 119 features (gnomAD_AF excluded) and is the automatic
  fallback when no AF/VEP source is available.

See `model/active_model.json` for the bundled model mapping.

## Citation

Please cite this Zenodo archive using the DOI provided on the Zenodo record.
Also cite the AGVarPred GitHub repository (see `code_reference.txt`).

The machine-readable citation metadata is in `CITATION.cff`.

## License

Apache-2.0. See `LICENSE`.

## Notes

- The large gnomAD VCF is not included; download instructions are in
  `data_manifest.json`.
- Source code is not duplicated here; use the GitHub release linked in
  `code_reference.txt`.
