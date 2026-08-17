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

The UI opens on an **Overview dashboard** (KPI tiles, cluster sizes, PCA scree,
AI-predicted composition) and includes a selectable **dark mode** whose chart
palette is validated for colorblind-safety in both light and dark surfaces.

| # | Tab | Backend capability / visualization |
|---|-----|--------------------|
| — | Overview | KPI tiles + cluster-size bar + PCA scree + composition donut |
| 1 | Load | Assemble AnnData from the 4 Space Ranger files + optional metadata |
| 2 | QC | Counts / genes / mito %, distributions, scatter, and filtering |
| 3 | Normalize | LogNormalize or Pearson residuals (SCT analog) + PCA elbow |
| 4 | Clusters | Leiden clustering + interactive UMAP |
| 5 | Spatial map | Expression / label overlay on the tissue image |
| 6 | Explore | Feature-on-UMAP + custom signature scoring |
| 7 | Markers / DE | `rank_genes_groups` markers, pairwise DE, **volcano plot** |
| 8 | Expression matrix | **Dot plot**, **marker heatmap**, **stacked violin** |
| 9 | Spatially variable | Moran's I via squidpy |
| 10 | Spatial statistics | **Neighborhood enrichment**, **co-occurrence**, cluster correlation |
| 11 | Trajectory | **PAGA graph** + **diffusion pseudotime** |
| 12 | Deconvolution | NNLS reference deconvolution + composition |
| 13 | **STIE (single-cell)** | **Native EM implementation of Zhu et al. 2024** — single-cell deconvolution / convolution / clustering integrating spot expression with nuclear morphology, with gap-area cell recovery |
| 14 | AI & Enrichment | Local Jaccard cell-type annotation + gseapy enrichment |
| 15 | Cell-Cell Comm. | squidpy `ligrec` (CellChat analog) |
| 16 | Integration | Multi-slice ZIP upload + harmonypy batch correction (Harmony analog) |
| 17 | Export | Download every table (CSV), the AnnData (.h5ad), and an HTML report |

### STIE — single-cell level (Nature Communications 2024)

Tab 13 is a faithful native-Python implementation of the **STIE** Expectation-
Maximization algorithm from Zhu, Kubota, Wang, Wang, Xiao & Hoshida,
*"STIE: Single-cell level deconvolution, convolution, and clustering in in situ
capturing-based spatial transcriptomics,"* Nat Commun **15**, 7559 (2024). It
jointly models spot-level gene expression and matched histology-image nuclear
morphology to reach the single-cell *level* (not just resolution):

- **Deconvolution / convolution** given a gene × cell-type signature, and
  **signature-free clustering** when none is supplied (with per-iteration
  signature re-estimation, yielding CAGE).
- **Bona-fide spot area** hyperparameter **γ** (default 2.5× the reported spot
  diameter) and morphology shrinkage penalty **λ** (Eq. 6–8).
- Per-spot penalized NNLS for the cell-type coefficients β, soft-EM refresh of
  the Gaussian morphology parameters μ_k / σ_k, and combined
  morphology + expression cell-type assignment (Eq. 9).
- **Gap-area recovery** of cells missed by all spots via nearest-spot
  neighborhood information.

Input is a nucleus-segmentation table (`cell_id`, pixel coordinates, and
morphology columns such as area / perimeter / circularity) — produced by any
segmentation tool (e.g. the paper's DeepImageJ Multi-Organ model, StarDist, or
Cellpose). The tab visualizes single cells on the tissue image, EM convergence,
composition, and the learned per-type nuclear-morphology profiles; all outputs
(per-cell assignments, spot proportions, morphology profiles, CAGE) are
downloadable.

### Advanced analytical concepts included

Dot plots · marker heatmaps · stacked violins · volcano plots · Moran's I
spatially-variable genes · **neighborhood enrichment** · **co-occurrence vs
distance** · **PAGA connectivity** · **diffusion pseudotime** · cluster
correlation dendrograms · ligand-receptor communication · NNLS deconvolution ·
local AI cell-type annotation · gseapy pathway enrichment · Harmony multi-slice
integration. Every visualization uses a colorblind-safe palette validated in
light and dark modes, always paired with a legend and a downloadable table.

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
