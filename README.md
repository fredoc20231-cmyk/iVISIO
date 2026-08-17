# iVisio — Spatial Omics Platform

iVisio is a full-stack web application for end-to-end analysis of **10x Genomics
Visium** spatial transcriptomics data. It takes Space Ranger output and guides
you through the complete workflow — load, QC, normalization, clustering, spatial
mapping, marker/DE discovery, cell-type deconvolution, AI-assisted annotation,
pathway enrichment, and cell-cell communication — in an interactive, tabbed UI,
with every result table and object available for download.

> **Note:** iVisio was originally prototyped as an R Shiny app (preserved for
> reference at [`legacy/app.R`](legacy/app.R)). It has been re-platformed to a
> **Python (FastAPI + scanpy/squidpy) backend** and a **React + TypeScript
> (Vite) frontend**.

## Architecture

```
┌──────────────────────────┐        REST / JSON        ┌───────────────────────────┐
│  React + TypeScript (Vite)│  ───────────────────────▶ │  FastAPI + scanpy/squidpy  │
│  Plotly visualizations    │  ◀─────────────────────── │  AnnData session store     │
│  12-step tabbed workflow  │        CSV / .h5ad         │  scanpy / squidpy / gseapy │
└──────────────────────────┘                            └───────────────────────────┘
        frontend/                                                backend/
```

- **backend/** — FastAPI service wrapping a scanpy pipeline. See
  [`backend/README.md`](backend/README.md) for the full R→Python feature-parity
  table and API reference.
- **frontend/** — React + TypeScript SPA (Vite) with Plotly charts and an
  interactive tissue overlay.

## Workflow (frontend tabs)

| # | Tab | Backend capability |
|---|-----|--------------------|
| 1 | Load | Assemble AnnData from the 4 Space Ranger files + optional metadata |
| 2 | QC | Counts / genes / mito %, distributions, and filtering |
| 3 | Normalize | LogNormalize or Pearson residuals (SCT analog) + PCA |
| 4 | Clusters | Leiden clustering + UMAP |
| 5 | Spatial map | Expression / label overlay on the tissue image |
| 6 | Explore | Feature-on-UMAP + custom signature scoring |
| 7 | Markers / DE | `rank_genes_groups` markers and pairwise DE |
| 8 | Spatially variable | Moran's I via squidpy |
| 9 | Deconvolution | NNLS reference deconvolution |
| 10 | AI & Enrichment | Local Jaccard cell-type annotation + gseapy enrichment |
| 11 | Cell-Cell Comm. | squidpy `ligrec` (CellChat analog) |
| 12 | Export | Download every table (CSV), the AnnData (.h5ad), and an HTML report |

## Quick start (Docker)

```bash
docker compose up --build
```

- Frontend: http://localhost:8080
- Backend API docs: http://localhost:8000/docs

## Local development

**Backend**
```bash
cd backend
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
uvicorn app.main:app --reload --port 8000
```

**Frontend**
```bash
cd frontend
npm install
npm run dev      # http://localhost:5173 (proxies /api to :8000)
```

## Downloads

Everything the platform computes is downloadable from the **Export** tab (or
directly at `GET /api/{session}/download/{artifact}`): marker/DE/SVG/AI/
enrichment/ligand-receptor/deconvolution tables and reference signature (CSV),
the spot metadata table (CSV), the full **AnnData `.h5ad`** object, and a
self-contained **HTML report**.

## Deployment notes

- Visium uploads are large. The bundled nginx config sets `client_max_body_size 0`;
  set a finite value if your platform requires one, and raise equivalent limits on
  any upstream load balancer / ingress.
- Optional analysis packages (`squidpy`, `gseapy`, `harmonypy`) degrade gracefully:
  if one is missing, only its endpoint returns a clear message — the core workflow
  keeps working.
