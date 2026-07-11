# Reproducibility

The complete research pipeline is in `AGVarPred-training/`.

## Steps

1. **Prepare data**
   - Obtain ClinVar `variant_summary.txt` from
     <https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz>
     and place it at the repository root.
   - Download the gnomAD exomes r2.1.1 liftover to GRCh38 VCF (see
     [`installation.md`](installation.md) for the exact `wget` commands).

2. **Install training dependencies**

   ```bash
   cd AGVarPred
   pip install -e ".[test]"
   pip install -r AGVarPred-training/requirements.txt
   ```

3. **Run the pipeline**

   ```bash
   cd AGVarPred-training
   make features   # slow; runs AlphaGenome feature extraction
   make train
   make evaluate
   make benchmark
   ```

   > `make evaluate` and `make benchmark` are thin wrappers that print the
   > location of the relevant scripts and outputs. See
   > `AGVarPred-training/scripts/model_training/` and
   > `AGVarPred-training/scripts/external_validation/` for the full set of
   > evaluation and benchmark scripts.

## Outputs

- `final_dataset_parts_*/` — assembled train/cal/test feature matrices.
- `feature_selection_output_nested/` — selected features and selection curves.
- `final_model_output_regularized/` — final calibrated model and metrics.
- `external_validation/results/` — benchmark scores and comparisons.

## Zenodo archive

Frozen model artifacts, splits, benchmark definitions, and metrics are archived
in `AGVarPred-zenodo/` and uploaded to Zenodo with a DOI. See
`AGVarPred-zenodo/code_reference.txt` for the GitHub release link and
`RELEASE_CHECKLIST.md` for the remaining Zenodo DOI step.
