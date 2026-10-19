# AGVarPred Web Server (prototype)

This directory contains an initial FastAPI service for the AGVarPred webserver project.

## Safety and scientific validity status

- This is a prototype, not a production or clinical service.
- Live prediction remains disabled by default. Runtime configuration must explicitly set `AGVARPRED_ENABLE_LIVE_PREDICTION=true` and `AGVARPRED_ALPHA_TERMS_REVIEWED=true`, and provide the AlphaGenome SDK/API key plus an indexed GRCh38 gnomAD exomes r2.1.1 VCF. The confirmation flag is an operator-controlled deployment gate; it is not a substitute for following the research team's approval and applicable provider requirements.
- No API keys belong in source control. Supply `ALPHAGENOME_API_KEY` and any other secrets through deployment environment variables/secrets.
- The service is designed to reuse the existing frozen AGVarPred package and bundled `model_full` model. It must not retrain, mutate, or silently fall back to the no-AF model.
- Requests accept normalized GRCh38 small variants as structured JSON (currently up to 25 per request). Feature generation uses the package’s AlphaGenome SDK source and local gnomAD VCF source, preserving the training-era VEP/AF source where configured. VCF upload, background jobs and richer visualizations remain follow-up milestones.

## Local development

From the repository root:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e ".[test]"
pip install -r webserver/requirements.txt
uvicorn webserver.app:app --reload
```

Open http://127.0.0.1:8000/docs to inspect the API.

For a first local run, keep live prediction disabled. The health endpoint and input validation are available without provider credentials. Live prediction requires the frozen feature schema to be present and the annotation defaults to follow the training pipeline; keep it disabled until the parity check below passes on real, release-matched data.

In a second terminal, start the React interface:

```bash
cd webserver/frontend
npm install
cp .env.example .env.local
npm run dev
```

Open http://127.0.0.1:5173. The form validates and submits inputs to the API. The adapter is implemented, but real-provider feature/prediction parity has not yet been demonstrated, so leave live prediction disabled until validation passes.

## Initial API

- `GET /api/health`: service status, frozen-model metadata and provider readiness (never returns secrets).
- `POST /api/predict`: accepts structured GRCh38 variants, performs live AlphaGenome + local gnomAD/VEP feature extraction when fully configured, and returns calibrated scores plus per-variant failures.
- If the provider is unavailable, if the model's usable AlphaGenome features are not returned, or if the operator has not explicitly confirmed review of provider terms, prediction is blocked instead of fabricating a score.

## Live feature and prediction parity check

The adapter calls the existing package's AlphaGenome SDK source and local gnomAD/VEP source, then uses the bundled frozen model. A parity checker is included, but it must be run against a small reference matrix from the established research feature pipeline using the same GRCh38 variant IDs. The checker makes real AlphaGenome API calls.

1. Configure `ALPHAGENOME_API_KEY` and `GNOMAD_VCF` in your local environment, and install the SDK version recorded in `src/AGVarPred/model/model_full/manifest.yaml`.
2. Choose one small parquet part file from the existing feature-extraction/training output. It must contain `variant_id` plus the selected model feature columns. Do not point it at the full training dataset.
3. From the repository root, run:

```bash
python -m webserver.validate_feature_parity \
  --reference /path/to/final_dataset_parts_test/part_0.parquet \
  --max-variants 5 \
  --output parity-report.json
```

The command compares the selected feature values, calibrated probabilities and predicted classes for the sampled variants and returns a non-zero exit code if it detects drift or missing AlphaGenome features. Inspect `parity-report.json` and resolve any discrepancies before enabling public live predictions. Running this check uses the AlphaGenome API and consumes provider quota.

## Current validation status

- Unit/API tests, existing package regression tests on Python 3.11–3.13, and the frontend build passed in CI for the preceding commits; check the latest commit's workflow status before merging.
- Real AlphaGenome API + local gnomAD/VEP end-to-end parity has **not** yet been run from a configured research environment.
- No production deployment has been performed. Keep live prediction disabled until the real-provider parity report has been reviewed and passes.
