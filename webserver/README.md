# AGVarPred Web Server (prototype)

This directory contains an initial FastAPI service for the AGVarPred webserver project.

## Safety and scientific validity status

- This is a prototype, not a production or clinical service.
- Live prediction is disabled by default. Set `AGVARPRED_ENABLE_LIVE_PREDICTION=true` only in an environment where the AlphaGenome API/model terms and deployment permissions have been reviewed.
- No API keys belong in source control. Supply `ALPHAGENOME_API_KEY` and any other secrets through deployment environment variables/secrets.
- The service is designed to reuse the existing frozen AGVarPred package and bundled `model_full` model. It must not retrain, mutate, or silently fall back to the no-AF model.
- Initial requests accept normalized GRCh38 variants as structured JSON. VCF upload, background jobs, live external annotation, and richer visualizations are follow-up milestones.

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

In a second terminal, start the React interface:

```bash
cd webserver/frontend
npm install
cp .env.example .env.local
npm run dev
```

Open http://127.0.0.1:5173. The form currently validates and submits inputs to the API; the backend intentionally returns a clear `503` response until the live feature-provider adapter is implemented and verified.

## Initial API

- `GET /api/health`: service status and model artifact metadata (no secret values).
- `POST /api/predict`: validates a small batch, then deliberately returns `503` unless live prediction has been explicitly enabled and a verified live feature-provider adapter is configured.

The initial implementation does not pretend to return predictions without the required AlphaGenome, VEP and gnomAD-derived features. That avoids scientific misrepresentation while we integrate and test a real provider.
