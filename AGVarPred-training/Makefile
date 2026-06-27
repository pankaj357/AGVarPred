# Makefile for reproducing the AGVarPred research pipeline.
# All script commands are run from the repository root so that hardcoded
# relative paths (e.g., external_data/, final_dataset_parts_*/) resolve correctly.

SHELL := /bin/bash
PYTHON := python
PROJECT_ROOT := $(abspath ..)

# Allow external scripts to find the repository root.
export PROJECT_ROOT := $(PROJECT_ROOT)

.PHONY: help features train evaluate benchmark clean

help:
	@echo "AGVarPred-training make targets:"
	@echo "  make features   - prepare chunks and extract AlphaGenome features"
	@echo "  make train      - run feature selection and train all models"
	@echo "  make evaluate   - report the location of the final model metrics"
	@echo "  make benchmark  - run external validation benchmarks (best effort)"
	@echo "  make clean      - remove intermediate outputs"

features:
	@echo "==> Building clean ClinVar gold set and train/cal/test splits..."
	cd $(PROJECT_ROOT) && $(PYTHON) AGVarPred-training/scripts/data_preparation/build_clean_clinvar_gold.py
	cd $(PROJECT_ROOT) && $(PYTHON) AGVarPred-training/scripts/data_preparation/split_clinvar_clean.py
	cd $(PROJECT_ROOT) && $(PYTHON) AGVarPred-training/scripts/data_preparation/build_alphagenome_input_split.py
	@echo "==> Creating chunks for parallel AlphaGenome scoring..."
	cd $(PROJECT_ROOT) && $(PYTHON) AGVarPred-training/scripts/feature_extraction/create_chunks_for_train.py
	cd $(PROJECT_ROOT) && $(PYTHON) AGVarPred-training/scripts/feature_extraction/create_chunks_cal.py
	cd $(PROJECT_ROOT) && $(PYTHON) AGVarPred-training/scripts/feature_extraction/create_chunks_test.py
	@echo "==> Running AlphaGenome feature extraction (requires SDK + API keys)..."
	@echo "Run the following in parallel (e.g., with GNU parallel or a job array):"
	@echo "  for i in {0..19}; do RUN_ID=\$$i python AGVarPred-training/scripts/feature_extraction/code.py; done"
	@echo "==> Assembling final dataset parts..."
	cd $(PROJECT_ROOT) && $(PYTHON) AGVarPred-training/scripts/feature_extraction/feature_labbling.py

train:
	@echo "==> Feature selection..."
	cd $(PROJECT_ROOT) && $(PYTHON) AGVarPred-training/scripts/feature_selection/feature_selection.py
	@echo "==> Training baselines..."
	cd $(PROJECT_ROOT) && $(PYTHON) AGVarPred-training/scripts/model_training/baseline_models_final.py
	@echo "==> Training XGBoost..."
	cd $(PROJECT_ROOT) && $(PYTHON) AGVarPred-training/scripts/model_training/train_xgboost_model.py
	@echo "==> Training final models..."
	cd $(PROJECT_ROOT) && $(PYTHON) AGVarPred-training/scripts/model_training/final_model_training.py
	cd $(PROJECT_ROOT) && $(PYTHON) AGVarPred-training/scripts/model_training/regularized_final_model.py

evaluate:
	@echo "==> Evaluating final regularized model..."
	@if [ -f $(PROJECT_ROOT)/final_model_output_regularized/final_metrics.csv ]; then \
		echo "Metrics available at: $(PROJECT_ROOT)/final_model_output_regularized/final_metrics.csv"; \
		cat $(PROJECT_ROOT)/final_model_output_regularized/final_metrics.csv; \
	else \
		echo "Metrics not found. Run 'make train' first."; \
	fi

benchmark:
	@echo "==> Running external validation benchmarks (best effort; individual scripts may be skipped if input data are missing)..."
	cd $(PROJECT_ROOT) && $(PYTHON) AGVarPred-training/scripts/external_validation/score_humsavar.py || true
	cd $(PROJECT_ROOT) && $(PYTHON) AGVarPred-training/scripts/external_validation/score_clinvar_3star.py || true
	cd $(PROJECT_ROOT) && $(PYTHON) AGVarPred-training/scripts/external_validation/score_mave_independent.py || true
	cd $(PROJECT_ROOT) && $(PYTHON) AGVarPred-training/scripts/external_validation/score_vip.py || true
	cd $(PROJECT_ROOT) && $(PYTHON) AGVarPred-training/scripts/external_validation/score_grimm2015.py || true
	cd $(PROJECT_ROOT) && $(PYTHON) AGVarPred-training/scripts/external_validation/score_gnomad_benign.py || true
	cd $(PROJECT_ROOT) && $(PYTHON) AGVarPred-training/scripts/external_validation/score_panelapp_benign.py || true
	cd $(PROJECT_ROOT) && $(PYTHON) AGVarPred-training/scripts/external_validation/score_dvd.py || true
	@echo "See $(PROJECT_ROOT)/external_validation/results/ for outputs and"
	@echo "$(PROJECT_ROOT)/AGVarPred-training/scripts/external_validation/ for the full list of benchmark scripts."

clean:
	@echo "==> Removing intermediate outputs..."
	cd $(PROJECT_ROOT) && rm -rf feature_extraction_train_dataset feature_extraction_cal_dataset feature_extraction_test_dataset
	cd $(PROJECT_ROOT) && rm -rf chunks_input_train chunks_input_cal chunks_input_test
	cd $(PROJECT_ROOT) && rm -rf alphagenome_input
	cd $(PROJECT_ROOT) && rm -rf final_dataset_parts_train final_dataset_parts_cal final_dataset_parts_test
	cd $(PROJECT_ROOT) && rm -rf logs_train logs_cal logs_test
	find $(PROJECT_ROOT) -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
