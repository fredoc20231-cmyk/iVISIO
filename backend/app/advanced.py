"""Advanced analytical concepts for iVisio.

These functions power the platform's advanced visualization tabs: expression
matrix views (dot plot, marker heatmap, stacked violin), spatial statistics
(neighborhood enrichment, co-occurrence), trajectory inference (PAGA +
diffusion pseudotime), cluster relationships (dendrogram / correlation), and
multi-slice batch integration. Each returns plain JSON-serializable structures
so the React frontend can render them with Plotly.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from scipy.sparse import issparse


def _log_expr(adata: AnnData) -> AnnData:
    """Log-normalized expression view (raw if present)."""
    return adata.raw.to_adata() if adata.raw is not None else adata


def _dense(mat) -> np.ndarray:
    return mat.toarray() if issparse(mat) else np.asarray(mat)


def _resolve_genes(adata: AnnData, genes: list[str]) -> list[str]:
    present = [g for g in genes if g in adata.var_names]
    if not present:
        raise ValueError("None of the requested genes were found in the dataset.")
    return present


# --------------------------------------------------------------------------- #
# Expression matrix views
# --------------------------------------------------------------------------- #
def dot_plot(adata: AnnData, genes: list[str], group_col: str) -> dict:
    """Mean expression (z-scored per gene) and fraction expressing, per group."""
    if group_col not in adata.obs:
        raise ValueError(f"Group column '{group_col}' not found. Run clustering first.")
    src = _log_expr(adata)
    genes = _resolve_genes(src, genes)
    groups = sorted(adata.obs[group_col].astype(str).unique(), key=_natural_key)
    labels = adata.obs[group_col].astype(str).values

    mean = np.zeros((len(groups), len(genes)))
    frac = np.zeros((len(groups), len(genes)))
    x = _dense(src[:, genes].X)
    for i, g in enumerate(groups):
        rows = labels == g
        block = x[rows]
        mean[i] = block.mean(axis=0)
        frac[i] = (block > 0).mean(axis=0)
    # z-score each gene (column) across groups for comparable color.
    mu = mean.mean(axis=0, keepdims=True)
    sd = mean.std(axis=0, keepdims=True)
    sd[sd == 0] = 1.0
    zmean = (mean - mu) / sd
    return {
        "genes": genes, "groups": groups,
        "mean_z": zmean.round(4).tolist(),
        "fraction": frac.round(4).tolist(),
    }


def marker_heatmap(adata: AnnData, markers: pd.DataFrame, group_col: str,
                   top_n: int = 5) -> dict:
    """Top-N markers per cluster -> genes x groups mean expression (z per gene)."""
    if markers is None or markers.empty:
        raise ValueError("Run marker analysis before building the heatmap.")
    gene_col = "feature" if "feature" in markers.columns else "gene"
    genes: list[str] = []
    for _, sub in markers.groupby("cluster"):
        sub = sub.sort_values("avg_log2FC", ascending=False)
        for g in sub[gene_col].head(top_n):
            if g not in genes and g in adata.var_names:
                genes.append(str(g))
    if not genes:
        raise ValueError("No marker genes available for the heatmap.")
    return _matrix_zscore(adata, genes, group_col)


def _matrix_zscore(adata: AnnData, genes: list[str], group_col: str) -> dict:
    src = _log_expr(adata)
    groups = sorted(adata.obs[group_col].astype(str).unique(), key=_natural_key)
    labels = adata.obs[group_col].astype(str).values
    x = _dense(src[:, genes].X)
    mat = np.zeros((len(genes), len(groups)))
    for j, g in enumerate(groups):
        mat[:, j] = x[labels == g].mean(axis=0)
    mu = mat.mean(axis=1, keepdims=True)
    sd = mat.std(axis=1, keepdims=True)
    sd[sd == 0] = 1.0
    z = (mat - mu) / sd
    return {"genes": genes, "groups": groups, "z": z.round(4).tolist()}


def violin(adata: AnnData, genes: list[str], group_col: str, max_per_group: int = 400) -> dict:
    """Per-group expression samples for stacked violin / box plots."""
    if group_col not in adata.obs:
        raise ValueError(f"Group column '{group_col}' not found. Run clustering first.")
    src = _log_expr(adata)
    genes = _resolve_genes(src, genes)
    groups = sorted(adata.obs[group_col].astype(str).unique(), key=_natural_key)
    labels = adata.obs[group_col].astype(str).values
    x = _dense(src[:, genes].X)
    rng = np.random.default_rng(1234)
    out: dict[str, dict[str, list[float]]] = {}
    for gi, gene in enumerate(genes):
        out[gene] = {}
        for grp in groups:
            vals = x[labels == grp, gi]
            if vals.shape[0] > max_per_group:
                vals = vals[rng.choice(vals.shape[0], max_per_group, replace=False)]
            out[gene][grp] = vals.round(4).tolist()
    return {"genes": genes, "groups": groups, "values": out}


# --------------------------------------------------------------------------- #
# Spatial statistics (squidpy)
# --------------------------------------------------------------------------- #
def neighborhood_enrichment(adata: AnnData, group_col: str, n_neighs: int = 6) -> dict:
    sq = _require_squidpy()
    if group_col not in adata.obs:
        raise ValueError(f"Group column '{group_col}' not found. Run clustering first.")
    adata.obs[group_col] = adata.obs[group_col].astype("category")
    sq.gr.spatial_neighbors(adata, n_neighs=int(n_neighs), coord_type="generic")
    sq.gr.nhood_enrichment(adata, cluster_key=group_col, seed=1234)
    z = adata.uns[f"{group_col}_nhood_enrichment"]["zscore"]
    cats = list(adata.obs[group_col].cat.categories.astype(str))
    return {"groups": cats, "zscore": np.asarray(z).round(3).tolist()}


def co_occurrence(adata: AnnData, group_col: str) -> dict:
    sq = _require_squidpy()
    if group_col not in adata.obs:
        raise ValueError(f"Group column '{group_col}' not found. Run clustering first.")
    adata.obs[group_col] = adata.obs[group_col].astype("category")
    sq.gr.co_occurrence(adata, cluster_key=group_col)
    res = adata.uns[f"{group_col}_co_occurrence"]
    occ = np.asarray(res["occ"])           # (n_cat, n_cat, n_intervals)
    interval = np.asarray(res["interval"])
    cats = list(adata.obs[group_col].cat.categories.astype(str))
    # Midpoints of distance bins for plotting.
    dist = ((interval[:-1] + interval[1:]) / 2).round(2).tolist()
    return {"groups": cats, "distance": dist, "occ": occ.round(4).tolist()}


# --------------------------------------------------------------------------- #
# Trajectory inference: PAGA + diffusion pseudotime
# --------------------------------------------------------------------------- #
def paga_trajectory(adata: AnnData, group_col: str) -> dict:
    if group_col not in adata.obs:
        raise ValueError(f"Group column '{group_col}' not found. Run clustering first.")
    if "neighbors" not in adata.uns:
        sc.pp.neighbors(adata, random_state=1234)
    adata.obs[group_col] = adata.obs[group_col].astype("category")
    sc.tl.paga(adata, groups=group_col)
    conn = _dense(adata.uns["paga"]["connectivities"])
    cats = list(adata.obs[group_col].cat.categories.astype(str))
    # Node positions = mean UMAP coordinate per group (falls back to PCA).
    key = "X_umap" if "X_umap" in adata.obsm else "X_pca"
    coords = adata.obsm[key][:, :2]
    labels = adata.obs[group_col].astype(str).values
    pos = {c: coords[labels == c].mean(axis=0).round(4).tolist() for c in cats}
    edges = []
    for i in range(len(cats)):
        for j in range(i + 1, len(cats)):
            w = float(conn[i, j])
            if w > 0.01:
                edges.append({"source": cats[i], "target": cats[j], "weight": round(w, 4)})
    return {"groups": cats, "positions": pos, "edges": edges}


def diffusion_pseudotime(adata: AnnData, root_group: str, group_col: str) -> dict:
    if group_col not in adata.obs:
        raise ValueError(f"Group column '{group_col}' not found. Run clustering first.")
    if "neighbors" not in adata.uns:
        sc.pp.neighbors(adata, random_state=1234)
    labels = adata.obs[group_col].astype(str).values
    root_cells = np.where(labels == str(root_group))[0]
    if root_cells.size == 0:
        raise ValueError(f"Root group '{root_group}' not found in '{group_col}'.")
    adata.uns["iroot"] = int(root_cells[0])
    sc.tl.diffmap(adata)
    sc.tl.dpt(adata)
    pt = adata.obs["dpt_pseudotime"].to_numpy()
    pt = np.where(np.isfinite(pt), pt, np.nan)
    payload = {"pseudotime": np.nan_to_num(pt, nan=0.0).round(4).tolist(),
               "barcodes": list(adata.obs_names)}
    if "X_umap" in adata.obsm:
        payload["umap"] = adata.obsm["X_umap"][:, :2].round(4).tolist()
    return payload


# --------------------------------------------------------------------------- #
# Cluster relationships
# --------------------------------------------------------------------------- #
def cluster_correlation(adata: AnnData, group_col: str) -> dict:
    """Pairwise correlation of per-cluster mean HVG profiles (dendrogram input)."""
    if group_col not in adata.obs:
        raise ValueError(f"Group column '{group_col}' not found. Run clustering first.")
    src = _log_expr(adata)
    if "highly_variable" in adata.var.columns:
        hvg = adata.var_names[adata.var["highly_variable"]]
        genes = [g for g in hvg if g in src.var_names][:2000]
    else:
        genes = list(src.var_names[:2000])
    x = _dense(src[:, genes].X)
    labels = adata.obs[group_col].astype(str).values
    groups = sorted(set(labels), key=_natural_key)
    profiles = np.vstack([x[labels == g].mean(axis=0) for g in groups])
    corr = np.corrcoef(profiles)
    return {"groups": groups, "correlation": np.nan_to_num(corr).round(4).tolist()}


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def _require_squidpy():
    try:
        import squidpy as sq  # noqa: WPS433
        return sq
    except Exception as exc:  # pragma: no cover - optional dependency
        raise ValueError(f"This analysis requires the 'squidpy' package. ({exc})")


def _natural_key(value: str):
    """Sort '0','1','2','10' numerically when possible, else lexically."""
    try:
        return (0, float(value))
    except (TypeError, ValueError):
        return (1, str(value))
