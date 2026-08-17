"""Analysis pipeline for iVisio, built on scanpy + optional squidpy/gseapy.

Every function is a pure operation on an AnnData object (or its derived
tables) so it can be called from the API layer without reaching into Shiny-
style reactive state. Advanced features that depend on optional packages
(squidpy, gseapy, harmonypy) import lazily and raise a clear error if the
package is absent, so the core workflow never breaks because of a missing
optional dependency.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from scipy.sparse import issparse

from .knowledge import AI_CELL_TYPE_DICT


# --------------------------------------------------------------------------- #
# QC
# --------------------------------------------------------------------------- #
def compute_qc(adata: AnnData, mito_prefix: str = "MT-") -> dict:
    """Compute per-spot counts, detected genes, and mitochondrial percentage."""
    adata.var["mt"] = adata.var_names.str.upper().str.startswith(mito_prefix.upper().lstrip("^"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False)
    # Friendly aliases matching the UI vocabulary.
    adata.obs["n_counts"] = adata.obs["total_counts"]
    adata.obs["n_genes"] = adata.obs["n_genes_by_counts"]
    adata.obs["percent_mt"] = adata.obs.get("pct_counts_mt", 0.0)
    return {
        "n_spots": int(adata.n_obs),
        "median_counts": float(np.median(adata.obs["n_counts"])),
        "median_genes": float(np.median(adata.obs["n_genes"])),
        "median_percent_mt": float(np.median(adata.obs["percent_mt"])),
        "mito_genes_found": int(adata.var["mt"].sum()),
    }


def filter_spots(
    adata: AnnData,
    min_counts: float,
    max_counts: float,
    min_genes: float,
    max_percent_mt: float,
) -> AnnData:
    if "n_counts" not in adata.obs:
        raise ValueError("Run QC metrics before applying filters.")
    keep = (
        (adata.obs["n_counts"] >= min_counts)
        & (adata.obs["n_counts"] <= max_counts)
        & (adata.obs["n_genes"] >= min_genes)
    )
    if "percent_mt" in adata.obs:
        keep &= adata.obs["percent_mt"] <= max_percent_mt
    return adata[keep.values].copy()


# --------------------------------------------------------------------------- #
# Normalization + PCA
# --------------------------------------------------------------------------- #
def normalize(adata: AnnData, method: str, n_top_genes: int, n_pcs: int) -> dict:
    """LogNormalize (scanpy) or an SCTransform-style analytic Pearson residual.

    Keeps a raw counts copy in ``adata.layers['counts']`` so downstream steps
    (marker tests, deconvolution) can access untransformed values.
    """
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()

    n_top_genes = int(min(n_top_genes, adata.n_vars))
    max_pcs = int(max(2, min(n_pcs, adata.n_obs - 1, adata.n_vars - 1)))

    if method == "sct":
        # Analytic Pearson residuals are scanpy's closest analog to SCTransform.
        sc.experimental.pp.normalize_pearson_residuals(adata)
        sc.experimental.pp.highly_variable_genes(
            adata, flavor="pearson_residuals", n_top_genes=n_top_genes
        )
        used = "pearson_residuals"
    else:
        adata.X = adata.layers["counts"].copy()
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        adata.raw = adata
        sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=n_top_genes)
        used = "lognormalize"

    sc.pp.scale(adata, max_value=10, zero_center=True)
    sc.tl.pca(adata, n_comps=max_pcs, svd_solver="arpack")
    variance_ratio = adata.uns["pca"]["variance_ratio"].tolist()
    return {"method": used, "n_pcs": max_pcs, "n_hvg": int(n_top_genes),
            "variance_ratio": variance_ratio}


# --------------------------------------------------------------------------- #
# Neighbors + clustering + UMAP
# --------------------------------------------------------------------------- #
def cluster(
    adata: AnnData,
    n_neighbors: int,
    n_pcs_low: int,
    n_pcs_high: int,
    resolution: float,
    metric: str = "cosine",
) -> dict:
    if "X_pca" not in adata.obsm:
        raise ValueError("Run normalization and PCA before clustering.")
    total_pcs = adata.obsm["X_pca"].shape[1]
    low = max(1, min(int(n_pcs_low), total_pcs))
    high = max(low, min(int(n_pcs_high), total_pcs))
    n_use = high  # scanpy uses leading PCs; honor the upper bound of the range
    k = int(min(max(5, n_neighbors), adata.n_obs - 1))

    sc.pp.neighbors(adata, n_neighbors=k, n_pcs=n_use, metric=metric, random_state=1234)
    sc.tl.leiden(adata, resolution=float(resolution), key_added="clusters", random_state=1234)
    sc.tl.umap(adata, min_dist=0.3, random_state=1234)

    counts = adata.obs["clusters"].value_counts().sort_index()
    return {
        "n_clusters": int(counts.shape[0]),
        "pcs_used": [low, high],
        "cluster_sizes": {str(k): int(v) for k, v in counts.items()},
    }


# --------------------------------------------------------------------------- #
# Markers + differential expression
# --------------------------------------------------------------------------- #
def _expr_layer(adata: AnnData) -> AnnData:
    """Return an AnnData view whose ``X`` is log-normalized for DE tests."""
    if adata.raw is not None:
        return adata.raw.to_adata()
    return adata


def find_markers(adata: AnnData, group_col: str, method: str = "wilcoxon",
                 min_pct: float = 0.1, logfc: float = 0.25, top_n: int = 0) -> pd.DataFrame:
    if group_col not in adata.obs:
        raise ValueError(f"Group column '{group_col}' not found. Run clustering first.")
    src = _expr_layer(adata)
    src.obs[group_col] = adata.obs[group_col].values
    if src.obs[group_col].nunique() < 2:
        raise ValueError(f"Group column '{group_col}' must contain at least two groups.")
    sc.tl.rank_genes_groups(src, groupby=group_col, method=method, pts=True)
    df = sc.get.rank_genes_groups_df(src, group=None)
    df = df.rename(columns={
        "names": "feature", "group": "cluster", "logfoldchanges": "avg_log2FC",
        "pvals": "p_val", "pvals_adj": "p_val_adj", "pct_nz_group": "pct.1",
        "pct_nz_reference": "pct.2",
    })
    if "pct.1" in df:
        df = df[(df["pct.1"] >= min_pct) | (df.get("pct.2", 0) >= min_pct)]
    df = df[df["avg_log2FC"].abs() >= logfc]
    df = df.sort_values(["cluster", "p_val_adj", "avg_log2FC"],
                        ascending=[True, True, False])
    if top_n and top_n > 0:
        df = df.groupby("cluster", group_keys=False).head(top_n)
    if df.empty:
        raise ValueError("No marker genes passed the thresholds. Lower min_pct or logFC.")
    return df.reset_index(drop=True)


def pairwise_de(adata: AnnData, group_col: str, group1: str, group2: str,
                method: str = "wilcoxon", min_pct: float = 0.1, logfc: float = 0.25) -> pd.DataFrame:
    if group_col not in adata.obs:
        raise ValueError(f"Group column '{group_col}' not found. Run clustering first.")
    if group1 == group2:
        raise ValueError("Group 1 and Group 2 must be different.")
    src = _expr_layer(adata)
    src.obs[group_col] = adata.obs[group_col].astype(str).values
    mask = src.obs[group_col].isin([str(group1), str(group2)])
    sub = src[mask.values].copy()
    sc.tl.rank_genes_groups(sub, groupby=group_col, groups=[str(group1)],
                            reference=str(group2), method=method, pts=True)
    df = sc.get.rank_genes_groups_df(sub, group=str(group1))
    df = df.rename(columns={
        "names": "feature", "logfoldchanges": "avg_log2FC", "pvals": "p_val",
        "pvals_adj": "p_val_adj", "pct_nz_group": "pct.1", "pct_nz_reference": "pct.2",
    })
    if "pct.1" in df:
        df = df[(df["pct.1"] >= min_pct) | (df.get("pct.2", 0) >= min_pct)]
    df = df[df["avg_log2FC"].abs() >= logfc]
    df = df.sort_values(["p_val_adj", "avg_log2FC"], ascending=[True, False])
    if df.empty:
        raise ValueError("No differential genes passed the thresholds.")
    return df.reset_index(drop=True)


# --------------------------------------------------------------------------- #
# Signature scoring
# --------------------------------------------------------------------------- #
def score_signature(adata: AnnData, name: str, genes: list[str]) -> dict:
    present = [g for g in genes if g in adata.var_names]
    if len(present) < 2:
        raise ValueError("Fewer than two signature genes were found in the dataset.")
    score_key = f"Signature_{name}"
    sc.tl.score_genes(adata, gene_list=present, score_name=score_key, use_raw=adata.raw is not None)
    return {"name": name, "column": score_key, "genes": present}


# --------------------------------------------------------------------------- #
# Spatially variable features (Moran's I via squidpy; markvariogram analog)
# --------------------------------------------------------------------------- #
def spatially_variable(adata: AnnData, n_features: int, n_neighbors: int = 6) -> pd.DataFrame:
    try:
        import squidpy as sq  # noqa: WPS433
    except Exception as exc:  # pragma: no cover - optional dependency
        raise ValueError(
            "Spatially-variable analysis requires the 'squidpy' package. "
            f"Install it to enable Moran's I. ({exc})"
        )
    if "spatial" not in adata.obsm:
        raise ValueError("Spatial coordinates are unavailable; reload the Visium files.")
    hvg = adata.var_names[adata.var["highly_variable"]] if "highly_variable" in adata.var else adata.var_names
    genes = list(hvg[: int(n_features)])
    sq.gr.spatial_neighbors(adata, n_neighs=int(n_neighbors), coord_type="generic")
    sq.gr.spatial_autocorr(adata, mode="moran", genes=genes, n_perms=100, seed=1234)
    res = adata.uns["moranI"].copy()
    res = res.reset_index().rename(columns={"index": "feature", "I": "morans_i",
                                            "pval_norm": "p_value", "pval_sim": "p_value_sim"})
    keep = [c for c in ["feature", "morans_i", "p_value", "p_value_sim", "pval_norm_fdr_bh"] if c in res.columns]
    return res[keep].sort_values("morans_i", ascending=False).reset_index(drop=True)


# --------------------------------------------------------------------------- #
# Reference deconvolution (projected-gradient NNLS, ported from the R app)
# --------------------------------------------------------------------------- #
def build_reference_signature(expr: pd.DataFrame, groups: list[str]) -> pd.DataFrame:
    """Mean expression per cell type -> genes x cell-type signature matrix."""
    groups = np.asarray([str(g) for g in groups])
    keep = (groups != "nan") & (groups != "")
    expr = expr.loc[:, keep]
    groups = groups[keep]
    levels = sorted(set(groups))
    sig = pd.DataFrame(
        {lev: expr.loc[:, groups == lev].mean(axis=1) for lev in levels}
    )
    return sig


def nnls_deconvolution(st_expr: np.ndarray, signature: np.ndarray,
                       max_iter: int = 300, tol: float = 1e-5) -> np.ndarray:
    """Simplex-constrained projected-gradient NNLS, one row of proportions per spot."""
    a = np.clip(signature.astype(float), 0, None)
    y = np.clip(st_expr.astype(float), 0, None)
    scale = np.sqrt((a ** 2).sum(axis=0))
    scale[scale == 0] = 1.0
    a = a / scale
    ata = a.T @ a
    step = 1.0 / max(1.0, ata.sum(axis=0).max())
    p = np.full((y.shape[0], a.shape[1]), 1.0 / a.shape[1])
    ya = y @ a
    for _ in range(max_iter):
        old = p
        grad = (p @ ata - ya) * step
        p = np.clip(p - grad, 0, None)
        rs = p.sum(axis=1, keepdims=True)
        rs[rs == 0] = 1.0
        p = p / rs
        if np.abs(p - old).max() < tol:
            break
    return p


def run_reference_deconvolution(adata: AnnData, signature: pd.DataFrame) -> pd.DataFrame:
    src = _expr_layer(adata)
    genes = [g for g in signature.index if g in src.var_names]
    if len(genes) < 10:
        raise ValueError("Fewer than 10 genes overlap the Visium matrix and reference signature.")
    sub = src[:, genes]
    x = sub.X.toarray() if issparse(sub.X) else np.asarray(sub.X)
    props = nnls_deconvolution(x, signature.loc[genes].to_numpy())
    return pd.DataFrame(props, index=adata.obs_names, columns=list(signature.columns))


# --------------------------------------------------------------------------- #
# AI heuristic cluster annotation (local Jaccard matching)
# --------------------------------------------------------------------------- #
def ai_annotate(markers: pd.DataFrame, top_n: int = 10,
                dictionary: dict[str, list[str]] | None = None) -> pd.DataFrame:
    dictionary = dictionary or AI_CELL_TYPE_DICT
    if markers is None or markers.empty:
        raise ValueError("Run marker analysis before AI annotation.")
    if "cluster" not in markers.columns:
        raise ValueError("Marker table has no 'cluster' column. Run 'find markers for all groups'.")
    gene_col = "feature" if "feature" in markers.columns else "gene"
    rows = []
    for cl, sub in markers.groupby("cluster"):
        sub = sub.sort_values("avg_log2FC", ascending=False)
        top = [str(g).upper() for g in sub[gene_col].head(top_n)]
        best_label, best_score = "Unknown / Novel", 0.0
        for label, dict_genes in dictionary.items():
            dg = {g.upper() for g in dict_genes}
            inter = len(set(top) & dg)
            union = len(set(top) | dg)
            jac = inter / union if union else 0.0
            if jac > best_score:
                best_score, best_label = jac, label
        rows.append({
            "cluster": str(cl),
            "predicted_cell_type": best_label,
            "confidence_pct": round(best_score * 100, 1),
            "top_markers": ", ".join(top[:5]),
        })
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# Pathway enrichment (gseapy Enrichr; clusterProfiler analog)
# --------------------------------------------------------------------------- #
def pathway_enrichment(genes: list[str], gene_set: str = "GO_Biological_Process_2021",
                       organism: str = "human", top_n: int = 50) -> pd.DataFrame:
    try:
        import gseapy as gp  # noqa: WPS433
    except Exception as exc:  # pragma: no cover - optional dependency
        raise ValueError(f"Pathway enrichment requires the 'gseapy' package. ({exc})")
    genes = sorted({str(g).upper().strip() for g in genes if str(g).strip()})
    if len(genes) < 3:
        raise ValueError("Fewer than three genes were supplied for enrichment.")
    enr = gp.enrichr(gene_list=genes, gene_sets=[gene_set], organism=organism, outdir=None)
    df = enr.results.sort_values("Adjusted P-value").head(top_n)
    keep = [c for c in ["Term", "Overlap", "P-value", "Adjusted P-value",
                        "Odds Ratio", "Combined Score", "Genes"] if c in df.columns]
    return df[keep].reset_index(drop=True)


# --------------------------------------------------------------------------- #
# Cell-cell communication (squidpy ligrec; CellChat analog)
# --------------------------------------------------------------------------- #
def ligand_receptor(adata: AnnData, group_col: str, n_perms: int = 100) -> pd.DataFrame:
    try:
        import squidpy as sq  # noqa: WPS433
    except Exception as exc:  # pragma: no cover - optional dependency
        raise ValueError(f"Cell-cell communication requires the 'squidpy' package. ({exc})")
    if group_col not in adata.obs:
        raise ValueError(f"Group column '{group_col}' not found. Run clustering first.")
    src = _expr_layer(adata)
    src.obs[group_col] = adata.obs[group_col].astype("category").values
    res = sq.gr.ligrec(
        src, cluster_key=group_col, n_perms=int(n_perms), seed=1234,
        use_raw=False, copy=True, threshold=0.01,
    )
    means = res["means"].copy()
    pvals = res["pvalues"].copy()
    means.columns = [f"{a}|{b}" for a, b in means.columns]
    pvals.columns = [f"{a}|{b}" for a, b in pvals.columns]
    long = means.stack().reset_index()
    long.columns = ["ligand", "receptor", "pair", "mean_expr"] if long.shape[1] == 4 else long.columns
    # Flatten ligand/receptor multiindex robustly:
    m = means.reset_index()
    id_cols = [c for c in m.columns if c not in means.columns]
    tidy = m.melt(id_vars=id_cols, var_name="source_target", value_name="mean_expr")
    tidy = tidy.rename(columns={id_cols[0]: "ligand", id_cols[1]: "receptor"}) if len(id_cols) >= 2 else tidy
    tidy = tidy.dropna(subset=["mean_expr"]).sort_values("mean_expr", ascending=False)
    return tidy.head(500).reset_index(drop=True)


# --------------------------------------------------------------------------- #
# Multi-slice integration (harmonypy; Harmony analog)
# --------------------------------------------------------------------------- #
def integrate_slices(adatas: list[AnnData], n_pcs: int = 20, use_harmony: bool = True) -> AnnData:
    if len(adatas) < 2:
        raise ValueError("At least two slices are required for integration.")
    for i, a in enumerate(adatas):
        a.obs["ivisio_slice"] = f"slice{i + 1}"
    merged = adatas[0].concatenate(*adatas[1:], batch_key="ivisio_slice_batch",
                                   index_unique="-")
    merged.layers["counts"] = merged.X.copy()
    sc.pp.normalize_total(merged, target_sum=1e4)
    sc.pp.log1p(merged)
    merged.raw = merged
    sc.pp.highly_variable_genes(merged, n_top_genes=2000)
    sc.pp.scale(merged, max_value=10)
    max_pcs = int(max(2, min(n_pcs, merged.n_obs - 1)))
    sc.tl.pca(merged, n_comps=max_pcs)
    if use_harmony:
        try:
            sc.external.pp.harmony_integrate(merged, "ivisio_slice")
            merged.uns["ivisio_integration"] = "harmony"
        except Exception:  # pragma: no cover - optional dependency
            merged.uns["ivisio_integration"] = "uncorrected (harmonypy unavailable)"
    else:
        merged.uns["ivisio_integration"] = "uncorrected merge"
    return merged
