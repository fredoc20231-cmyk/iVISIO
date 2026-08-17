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
