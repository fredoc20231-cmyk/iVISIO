"""End-to-end API test: demo load through the core analysis workflow."""

import pytest

pytest.importorskip("fastapi")
from fastapi.testclient import TestClient  # noqa: E402

from app.main import app  # noqa: E402

client = TestClient(app)


def _new_session() -> str:
    r = client.post("/api/session")
    assert r.status_code == 200
    return r.json()["session_id"]


def test_health():
    assert client.get("/api/health").json()["status"] == "ok"


def test_full_workflow_on_demo():
    sid = _new_session()

    assert client.post(f"/api/{sid}/load-demo").status_code == 200

    st = client.get(f"/api/{sid}/status").json()
    assert st["loaded"] and st["n_spots"] > 0 and st["has_image"]

    assert client.post(f"/api/{sid}/qc", json={"mito_prefix": "MT-"}).status_code == 200
    assert client.post(f"/api/{sid}/normalize",
                       json={"method": "log", "n_top_genes": 100, "n_pcs": 10}).status_code == 200
    assert client.post(f"/api/{sid}/cluster",
                       json={"n_neighbors": 15, "n_pcs_low": 1, "n_pcs_high": 10,
                             "resolution": 0.5, "metric": "euclidean"}).status_code == 200

    st = client.get(f"/api/{sid}/status").json()
    assert st["has_clusters"] and st["has_umap"]

    mk = client.post(f"/api/{sid}/markers", json={"group_col": "clusters"})
    assert mk.status_code == 200 and mk.json()["total"] > 0

    # Spatial overlay + downloads.
    assert client.get(f"/api/{sid}/image").status_code == 200
    assert client.get(f"/api/{sid}/download/adata").status_code == 200
    assert client.get(f"/api/{sid}/download/markers").status_code == 200
    assert client.get(f"/api/{sid}/download/metadata").status_code == 200


def test_stie_endpoint_on_demo(tmp_path):
    sid = _new_session()
    client.post(f"/api/{sid}/load-demo")
    client.post(f"/api/{sid}/normalize", json={"method": "log", "n_top_genes": 100, "n_pcs": 10})

    # Pull the generated demo cells + signature CSVs and feed them to STIE.
    cells_csv = client.get(f"/api/{sid}/download/demo_cells").content
    sig_csv = client.get(f"/api/{sid}/download/demo_signature").content
    files = {
        "cells_csv": ("cells.csv", cells_csv, "text/csv"),
        "signature_csv": ("sig.csv", sig_csv, "text/csv"),
    }
    data = {"features": "area,perimeter,circularity", "gamma": "2.5", "lam": "0", "steps": "8"}
    r = client.post(f"/api/{sid}/stie/deconvolve", files=files, data=data)
    assert r.status_code == 200, r.text
    body = r.json()
    assert body["n_cells"] > 0
    assert client.get(f"/api/{sid}/stie/cells").status_code == 200
    assert client.get(f"/api/{sid}/download/stie_cells").status_code == 200
