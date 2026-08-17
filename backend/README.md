# iVisio backend (FastAPI + scanpy)

Python re-implementation of the iVisio Visium spatial-omics workflow, exposed
as a REST API for the React frontend.

## Feature parity with the original R app

| Workflow step | R (Seurat) | Python analog used here |
|---------------|-----------|--------------------------|
| Load Space Ranger output | `Load10X_Spatial` / manual | `scanpy.read_10x_h5` + custom overlay assembly (`io_visium.py`) |
| QC | `PercentageFeatureSet` | `sc.pp.calculate_qc_metrics` |
| Normalize | LogNormalize / SCTransform | `normalize_total`+`log1p` / analytic Pearson residuals |
| PCA / clustering | `RunPCA`/`FindClusters` (Louvain) | `sc.tl.pca` / `sc.tl.leiden` |
| UMAP | `RunUMAP` | `sc.tl.umap` |
| Markers / DE | `FindAllMarkers` / `FindMarkers` | `sc.tl.rank_genes_groups` |
| Signature scoring | `AddModuleScore` | `sc.tl.score_genes` |
| Spatially variable | markvariogram | **Moran's I** via `squidpy.gr.spatial_autocorr` |
| Reference deconvolution | projected-gradient NNLS | same algorithm, vectorized in NumPy |
| AI cluster annotation | local Jaccard dictionary | ported verbatim (`knowledge.py`) |
| Pathway enrichment | clusterProfiler | `gseapy.enrichr` |
| Cell-cell communication | CellChat | `squidpy.gr.ligrec` |
| Multi-slice integration | Harmony | `harmonypy` via `sc.external.pp.harmony_integrate` |
| **STIE single-cell EM** | STIE R package | **native NumPy/SciPy EM** (`stie.py`) — no R dependency |

### STIE (`app/stie.py`)

A faithful re-implementation of the STIE EM algorithm (Zhu et al., *Nat Commun*
2024, doi:10.1038/s41467-024-51728-5). It jointly models spot expression (Λ, X, β)
and Gaussian nuclear morphology (μ_k, σ_k), using the bona-fide spot area γ and
morphology shrinkage penalty λ. The M-step gene-expression update (Eq. 6–8) is a
per-spot penalized NNLS surrogate of the paper's quadprog formulation; morphology
parameters are refreshed from soft responsibilities (Eq. 4); cells are typed by
the combined morphology + expression score (Eq. 9); the signature is re-estimated
by NNLS in clustering mode (Eq. 10); and gap-area cells are recovered via their
nearest spot. Endpoints: `POST /stie/deconvolve`, `POST /stie/cluster`,
`GET /stie/cells`, `GET /stie/morphology`; downloads: `stie_cells`,
`stie_morphology`, `stie_signature` (CAGE), `stie_spot_proportions`.

> Faithfulness note: the geometry uses `spot_diameter_fullres` from the
> scale-factor JSON and the full-resolution tissue-position pixel coordinates, so
> the cells-on-image CSV must be in the same full-resolution pixel space as the
> histology image. The penalized-NNLS M-step is a convex surrogate of the exact
> quadprog objective; results are directionally equivalent but not bit-identical
> to the reference R package.

Optional packages (`squidpy`, `gseapy`, `harmonypy`) are imported lazily; if a
package is missing, the relevant endpoint returns a clear 400 rather than
crashing, exactly like the R app's `requireNamespace()` guards.

## Run locally

```bash
cd backend
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
uvicorn app.main:app --reload --port 8000
```

Interactive API docs are then at http://localhost:8000/docs.

## Downloads

Every result is downloadable via `GET /api/{sid}/download/{artifact}`:
`markers`, `de`, `svg`, `ai`, `enrichment`, `ligrec`, `deconv`, `signature`
(CSV), `metadata` (spot table CSV), `adata` (full `.h5ad`), and `report`
(self-contained HTML).
