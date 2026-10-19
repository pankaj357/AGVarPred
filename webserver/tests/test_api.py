from fastapi.testclient import TestClient

from webserver.app import app

client = TestClient(app)


def test_health_reports_frozen_full_model(monkeypatch):
    # Isolate the test from any local gnomAD configuration.
    monkeypatch.delenv("GNOMAD_TSV", raising=False)
    monkeypatch.delenv("GNOMAD_VCF", raising=False)

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

def test_predict_rejects_identical_ref_alt():
    response = client.post(
        "/api/predict",
        json={
            "genome_build": "GRCh38",
            "variants": [{"chrom": "chr1", "pos": 12345, "ref": "A", "alt": "A"}],
        },
    )
    assert response.status_code == 422


def test_predict_rejects_invalid_chromosome():
    response = client.post(
        "/api/predict",
        json={
            "genome_build": "GRCh38",
            "variants": [{"chrom": "chr23", "pos": 12345, "ref": "A", "alt": "G"}],
        },
    )
    assert response.status_code == 422

def test_health_never_enables_live_mode_without_terms_confirmation(monkeypatch):
    monkeypatch.setenv("AGVARPRED_ENABLE_LIVE_PREDICTION", "true")
    monkeypatch.delenv("AGVARPRED_ALPHA_TERMS_REVIEWED", raising=False)
    monkeypatch.delenv("ALPHAGENOME_API_KEY", raising=False)
    monkeypatch.delenv("GNOMAD_TSV", raising=False)
    monkeypatch.delenv("GNOMAD_VCF", raising=False)

    response = client.get("/api/health")
    assert response.status_code == 200
    body = response.json()
    assert body["live_prediction_enabled"] is False
    assert body["live_provider_ready"] is False


def test_predict_requires_explicit_terms_confirmation(monkeypatch):
    monkeypatch.setenv("AGVARPRED_ENABLE_LIVE_PREDICTION", "true")
    monkeypatch.delenv("AGVARPRED_ALPHA_TERMS_REVIEWED", raising=False)
    response = client.post(
        "/api/predict",
        json={
            "genome_build": "GRCh38",
            "variants": [{"chrom": "1", "pos": 12345, "ref": "A", "alt": "G"}],
        },
    )
    assert response.status_code == 503
    assert "terms" in response.json()["detail"].lower()


def test_predict_rejects_duplicate_variants():
    variant = {"chrom": "chr1", "pos": 12345, "ref": "A", "alt": "G"}
    response = client.post(
        "/api/predict",
        json={"genome_build": "GRCh38", "variants": [variant, variant]},
    )
    assert response.status_code == 422


def test_cors_preflight_from_allowed_origin():
    response = client.options(
        "/api/predict",
        headers={
            "Origin": "http://127.0.0.1:5173",
            "Access-Control-Request-Method": "POST",
            "Access-Control-Request-Headers": "Content-Type",
        },
    )
    assert response.status_code == 200
    assert response.headers.get("access-control-allow-origin") == "http://127.0.0.1:5173"
    assert "POST" in (response.headers.get("access-control-allow-methods") or "")


def test_cors_headers_on_cross_origin_health_request():
    response = client.get("/api/health", headers={"Origin": "http://localhost:5173"})
    assert response.status_code == 200
    assert response.headers.get("access-control-allow-origin") == "http://localhost:5173"


def test_cors_helper_parses_origin_env(monkeypatch):
    """_cors_origins reads the comma-separated deployment variable."""
    from webserver.app import _cors_origins

    monkeypatch.setenv(
        "AGVARPRED_CORS_ORIGINS",
        "https://trusted.example.com, https://app.example.com",
    )
    assert _cors_origins() == [
        "https://trusted.example.com",
        "https://app.example.com",
    ]


def test_health_rate_limit_enforced():
    """Repeated health calls from the same address receive HTTP 429."""
    from webserver.app import app

    app.state.limiter.reset()
    # Exhaust the 10/minute limit for the test client address.
    for _ in range(10):
        response = client.get("/api/health")
        assert response.status_code == 200

    response = client.get("/api/health")
    assert response.status_code == 429
    assert "rate limit" in response.json()["detail"].lower()
