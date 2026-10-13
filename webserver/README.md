# AGVarPred Web Server (prototype)

This directory contains an initial FastAPI service for the AGVarPred webserver project.

## Safety and scientific validity status

- This is a prototype, not a production or clinical service.
- Live prediction is disabled by default. It requires `AGVARPRED_ENABLE_LIVE_PREDICTION=true`, `AGVARPRED_ALPHA_TERMS_REVIEWED=true`, the AlphaGenome SDK and API key, and an indexed GRCh38 gnomAD exomes r2.1.1 VCF. Do not enable it until the AlphaGenome API/model terms and deployment permissions have been reviewed.
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

For a first local run, keep live prediction disabled. The health endpoint and input validation are available without provider credentials, but prediction requests deliberately return `503` until the feature provider is ready and explicitly enabled.

In a second terminal, start the React interface:

```bash
cd webserver/frontend
npm install
cp .env.example .env.local
npm run dev
```

Open http://127.0.0.1:5173. The form currently validates and submits inputs to the API; the backend intentionally returns a clear `503` response until the live feature-provider adapter is implemented and verified.

## Initial API

- `GET /api/health`: service status, frozen-model metadata and provider readiness (never returns secrets).
- `POST /api/predict`: accepts structured GRCh38 variants, performs live AlphaGenome + local gnomAD/VEP feature extraction when fully configured, and returns calibrated scores plus per-variant failures.
- If the provider is unavailable, if the model's usable AlphaGenome features are not returned, or if the operator has not explicitly confirmed review of provider terms, prediction is blocked instead of fabricating a score.

The feature-extraction adapter is new and still requires successful integration/parity testing with real authorized provider access before any public deployment.
