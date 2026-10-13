from fastapi.testclient import TestClient

from webserver.app import app

client = TestClient(app)


def test_health_reports_frozen_full_model():
    response = client.get("/api/health")
    assert response.status_code == 200
    body = response.json()
    assert body["status"] == "ok"
    assert body["model_version"] == "model_full"
    assert body["model_type"] == "full"
    assert body["feature_count"] == 120
    assert body["threshold"] == 0.42
    assert body["live_provider_ready"] is False


def test_predict_is_disabled_by_default(monkeypatch):
    monkeypatch.delenv("AGVARPRED_ENABLE_LIVE_PREDICTION", raising=False)
    response = client.post(
        "/api/predict",
        json={
            "genome_build": "GRCh38",
            "variants": [{"chrom": "1", "pos": 12345, "ref": "A", "alt": "G", "gene": "TEST1"}],
        },
    )
    assert response.status_code == 503
    assert "not enabled" in response.json()["detail"].lower()


def test_predict_rejects_wrong_build():
    response = client.post(
        "/api/predict",
        json={
            "genome_build": "GRCh37",
            "variants": [{"chrom": "1", "pos": 12345, "ref": "A", "alt": "G"}],
        },
    )
    assert response.status_code == 422


def test_predict_rejects_too_many_variants():
    variants = [{"chrom": "1", "pos": i + 1, "ref": "A", "alt": "G"} for i in range(26)]
    response = client.post(
        "/api/predict",
        json={"genome_build": "GRCh38", "variants": variants},
    )
    assert response.status_code == 422


def test_predict_rejects_symbolic_alleles():
    response = client.post(
        "/api/predict",
        json={
            "genome_build": "GRCh38",
            "variants": [{"chrom": "1", "pos": 12345, "ref": "A", "alt": "<DEL>"}],
        },
    )
    assert response.status_code == 422
