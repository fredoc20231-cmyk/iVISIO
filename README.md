# iVisio — 10x Visium Spatial-Omics Platform

iVisio is a self-contained [R Shiny](https://shiny.posit.co/) application for
end-to-end analysis of 10x Genomics **Visium** spatial transcriptomics data. It
takes Space Ranger output and walks you through a complete workflow — from
loading and QC through clustering, spatial mapping, marker discovery, cell-type
deconvolution, and report export — in an interactive, tabbed interface.

## Workflow

The app is organized as a numbered set of tabs that mirror a standard Visium
analysis pipeline:

| # | Tab | What it does |
|---|-----|--------------|
| 1 | **Load** | Upload the filtered `.h5` matrix, tissue image, scale factors, and tissue positions (plus optional metadata CSV) |
| 2 | **QC** | Compute per-spot counts/features/mitochondrial % and apply filters |
| 3 | **Normalize** | LogNormalize or SCTransform, variable feature selection, and PCA |
| 4 | **Clusters** | Nearest-neighbor graph, UMAP, and Louvain clustering |
| 5 | **Spatial map** | Overlay gene expression or cluster labels on the tissue image |
| 6 | **Explore** | Plot individual features and score custom gene signatures |
| 7 | **Markers / DE** | Cluster markers and pairwise differential expression (fast vectorized or Seurat Wilcoxon) |
| 8 | **STIE** | Integrate spot expression with nuclear morphology for cell-level typing (requires the STIE package) |
| 9 | **Reference / Deconvolution** | Estimate cell-type proportions per spot from a reference signature (NNLS) |
| 10 | **Spatially variable** | Identify spatially variable features |
| 11 | **Export** | Download the Seurat object, tables, spot metadata, and an HTML report |

## Input

The **Load** tab accepts each Space Ranger file separately (this avoids large
multi-file requests and works with GSM-prefixed filenames):

1. `filtered_feature_bc_matrix.h5` — filtered expression matrix
2. Tissue image — `tissue_hires_image.png` (or `.jpg`/`.jpeg`)
3. `scalefactors_json.json` — scale factors
4. `tissue_positions_list.csv` / `tissue_positions.csv` — spot positions

An optional metadata CSV (with a `barcode` column) can be joined onto the spots.

## Requirements

Dependencies are **not** installed by the app — install them once in a separate
R session before launching:

```r
install.packages(c(
  "shiny", "bslib", "shinycssloaders", "DT", "plotly", "ggplot2",
  "dplyr", "readr", "tidyr", "zip", "rmarkdown", "png", "jsonlite"
))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("Seurat", "SeuratObject", "sp", "sctransform", "SummarizedExperiment"))
```

The optional **STIE** integration can be installed from within the app (tab 8)
or manually:

```r
install.packages(c("quadprog", "magick", "remotes"))
remotes::install_github("zhushijia/STIE")
```

## Running

```r
shiny::runApp("app.R")
```

### Deployment notes

- The app disables Shiny's internal request-size cap (`shiny.maxRequestSize`),
  but reverse proxies and hosting platforms enforce their own limits. For nginx,
  set `client_max_body_size 0;` (or a large finite value). Raise the equivalent
  limit for Apache (`LimitRequestBody`), Posit Connect, shinyapps.io, Kubernetes
  ingress, or your cloud load balancer.
- The app intentionally does not install packages at startup or execute uploaded
  files.
