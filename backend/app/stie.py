"""STIE — Spatial Transcriptomics + Image-based nuclear morphology EM.

A faithful Python re-implementation of the algorithm in

    Zhu, Kubota, Wang, Wang, Xiao & Hoshida.
    "STIE: Single-cell level deconvolution, convolution, and clustering in in
    situ capturing-based spatial transcriptomics." Nat Commun 15, 7559 (2024).
    https://doi.org/10.1038/s41467-024-51728-5

The method bridges the gap between spot resolution and single-cell *level* by
jointly modeling (a) spot-level gene expression and (b) matched histology-image
nuclear morphology through an Expectation-Maximization algorithm. Given a
cell-type transcriptomic signature it performs single-cell **deconvolution**
(splitting low-resolution spots into cells) and **convolution** (assembling
high-resolution sub-cell pieces); given no signature it performs signature-free
single-cell **clustering**. It also recovers cells in the ~70% gap area between
spots via spot-neighborhood information.

Correspondence to the paper's Methods
-------------------------------------
* Λ  (`signature`)          genes × K cell-type expression signature.
* X  (`expr`)               spots × genes spatial expression.
* β  (`spot_beta`)          per-spot non-negative cell-type coefficients
                            (∝ cell counts); proportions Pᵉˣᵖʳ = β / Σβ.
* μ_k, σ_k (`mu`,`sigma`)   Gaussian nuclear-morphology parameters per type.
* γ  (`gamma`)              bona-fide spot area = γ × reported spot diameter
                            (paper default 2.5×).
* λ  (`lam`)                morphology shrinkage penalty balancing the
                            transcriptome vs. morphology contribution (Eq. 6–8).

The M-step gene-expression update (Eq. 6–8) is a per-spot penalized NNLS:
minimize ‖X_s − Λβ_s‖² + λ‖β_s − nₛ·Pᵐᵒʳᵖʰ_s‖² over β_s ≥ 0, solved as an
augmented non-negative least squares problem — a convex surrogate of the paper's
quadprog formulation that couples the expression fit to the morphology-derived
expected cell counts. Morphology μ_k, σ_k are refreshed from soft responsibilities
(Eq. 4); cell types are assigned by the combined morphology + expression score
(Eq. 9); the signature is re-estimated by NNLS in clustering mode (Eq. 10).
Deviations from the reference R package are documented inline and in the README.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy.optimize import nnls
from scipy.spatial import cKDTree
from sklearn.cluster import KMeans

_EPS = 1e-9
_VAR_FLOOR = 1e-3


@dataclass
class STIEResult:
    mode: str
    cell_types: list[str]
    cells: pd.DataFrame                      # cell_id, spot, x, y, cell_type, prob, recovered
    spot_beta: pd.DataFrame                  # spots × K (estimated cell counts)
    spot_prop: pd.DataFrame                  # spots × K (proportions)
    morphology: pd.DataFrame                 # cell_type × feature (mean ± sd)
    signature: pd.DataFrame | None           # genes × K (CAGE, clustering mode)
    rmse_trace: list[float] = field(default_factory=list)
    n_recovered: int = 0
    features: list[str] = field(default_factory=list)
    gamma: float = 2.5
    lam: float = 0.0


# --------------------------------------------------------------------------- #
# Geometry: assign cells to the bona-fide area of each spot
# --------------------------------------------------------------------------- #
def _coverage(cell_xy: np.ndarray, spot_xy: np.ndarray, radius: float):
    """Return (cells_per_spot, spot_of_cell, nearest_spot) using a KD-tree.

    A cell is covered by a spot when the distance between their centroids is
    below ``radius`` = γ × (spot_diameter / 2). Cells covered by no spot are
    later recovered via their nearest spot (paper's gap-filling step).
    """
    tree = cKDTree(spot_xy)
    covered = tree.query_ball_point(cell_xy, r=radius)   # list per cell
    _, nearest = tree.query(cell_xy, k=1)
    cells_per_spot: dict[int, list[int]] = {s: [] for s in range(len(spot_xy))}
    spot_of_cell = np.full(len(cell_xy), -1, dtype=int)
    for ci, spots in enumerate(covered):
        if spots:
            # Assign the cell to its closest covering spot for the E/M coupling.
            d = np.linalg.norm(spot_xy[spots] - cell_xy[ci], axis=1)
            s = spots[int(np.argmin(d))]
            spot_of_cell[ci] = s
            cells_per_spot[s].append(ci)
    return cells_per_spot, spot_of_cell, nearest


# --------------------------------------------------------------------------- #
# Morphology likelihood (diagonal Gaussian)
# --------------------------------------------------------------------------- #
def _morph_loglik(feats: np.ndarray, mu: np.ndarray, sigma: np.ndarray) -> np.ndarray:
    """C×K matrix of log N(feature_c | μ_k, σ_k) for diagonal covariance."""
    var = np.clip(sigma ** 2, _VAR_FLOOR, None)             # K×M
    # (C,1,M) - (1,K,M) -> (C,K,M)
    diff2 = (feats[:, None, :] - mu[None, :, :]) ** 2
    ll = -0.5 * (np.log(2 * np.pi * var)[None] + diff2 / var[None])
    return ll.sum(axis=2)                                    # C×K


# --------------------------------------------------------------------------- #
# M-step: per-spot penalized NNLS for β (Eq. 6–8 surrogate)
# --------------------------------------------------------------------------- #
def _solve_beta(x_s: np.ndarray, lam_mat: np.ndarray, pmorph_s: np.ndarray,
                n_s: float, lam: float) -> np.ndarray:
    """min ‖x_s − Λβ‖² + λ‖β − nₛ·Pᵐᵒʳᵖʰ‖² s.t. β ≥ 0 via augmented NNLS."""
    K = lam_mat.shape[1]
    if n_s <= 0:
        return np.zeros(K)
    if lam > 0:
        aug_a = np.vstack([lam_mat, np.sqrt(lam) * np.eye(K)])
        aug_b = np.concatenate([x_s, np.sqrt(lam) * pmorph_s * n_s])
    else:
        aug_a, aug_b = lam_mat, x_s
    beta, _ = nnls(aug_a, aug_b)
    # Rescale so the coefficients sum to the observed cell count in the spot,
    # keeping β interpretable as per-type cell counts.
    tot = beta.sum()
    if tot > _EPS:
        beta = beta / tot * n_s
    return beta


# --------------------------------------------------------------------------- #
# Core EM shared by deconvolution and clustering
# --------------------------------------------------------------------------- #
def _run_em(expr: pd.DataFrame, signature: pd.DataFrame, cells: pd.DataFrame,
            spot_index: list[str], spot_xy: np.ndarray, feature_cols: list[str],
            radius: float, lam: float, steps: int, update_signature: bool,
            cell_type_names: list[str]):
    genes = [g for g in signature.index if g in expr.columns]
    if len(genes) < 5:
        raise ValueError("Fewer than 5 genes overlap the spot matrix and signature.")
    lam_mat = signature.loc[genes].to_numpy(dtype=float)          # G×K
    x = expr.loc[spot_index, genes].to_numpy(dtype=float)         # S×G
    K = lam_mat.shape[1]

    cell_xy = cells[["_x", "_y"]].to_numpy(dtype=float)
    raw_feats = cells[feature_cols].to_numpy(dtype=float)
    fmean = raw_feats.mean(axis=0)
    fstd = raw_feats.std(axis=0)
    fstd[fstd == 0] = 1.0
    feats = (raw_feats - fmean) / fstd                            # standardized

    cells_per_spot, spot_of_cell, nearest = _coverage(cell_xy, spot_xy, radius)
    spot_pos = {s: i for i, s in enumerate(spot_index)}

    # Initialize morphology params from a K-means partition of the nuclei.
    init = KMeans(n_clusters=K, n_init=10, random_state=1234).fit_predict(feats)
    mu = np.vstack([feats[init == k].mean(axis=0) if (init == k).any() else feats.mean(axis=0)
                    for k in range(K)])
    sigma = np.vstack([feats[init == k].std(axis=0) if (init == k).sum() > 1 else np.ones(feats.shape[1])
                       for k in range(K)])
    resp = np.eye(K)[init]                                         # C×K soft responsibilities
    n_s = np.array([len(cells_per_spot[s]) for s in range(len(spot_xy))], dtype=float)

    rmse_trace: list[float] = []
    beta = np.tile(np.where(n_s[:, None] > 0, n_s[:, None] / K, 0.0), (1, K))

    for _ in range(max(1, steps)):
        # --- Morphology-derived proportions per spot (Pᵐᵒʳᵖʰ) ---
        pmorph = np.zeros((len(spot_xy), K))
        for s, members in cells_per_spot.items():
            if members:
                pmorph[s] = resp[members].mean(axis=0)

        # --- M-step (expression): per-spot penalized NNLS for β ---
        beta = np.zeros((len(spot_xy), K))
        for s in range(len(spot_xy)):
            if n_s[s] > 0:
                beta[s] = _solve_beta(x[s], lam_mat, pmorph[s], n_s[s], lam)
        row = beta.sum(axis=1, keepdims=True)
        pexpr = np.divide(beta, np.where(row > _EPS, row, 1.0))

        # --- E-step: combine morphology likelihood with expression prior ---
        loglik = _morph_loglik(feats, mu, sigma)                  # C×K
        log_prior = np.zeros_like(loglik)
        for ci in range(len(cell_xy)):
            s = spot_of_cell[ci]
            if s >= 0:
                log_prior[ci] = np.log(pexpr[s] + _EPS)
        combined = loglik + log_prior
        combined -= combined.max(axis=1, keepdims=True)
        resp = np.exp(combined)
        resp /= resp.sum(axis=1, keepdims=True) + _EPS

        # --- M-step (morphology): refresh μ_k, σ_k from responsibilities ---
        wsum = resp.sum(axis=0) + _EPS
        mu = (resp.T @ feats) / wsum[:, None]
        var = (resp.T @ (feats ** 2)) / wsum[:, None] - mu ** 2
        sigma = np.sqrt(np.clip(var, _VAR_FLOOR, None))

        # --- M-step (signature): re-estimate Λ by NNLS in clustering mode ---
        if update_signature:
            new_sig = np.zeros_like(lam_mat)
            for gi in range(len(genes)):
                coef, _ = nnls(pexpr, x[:, gi])
                new_sig[gi] = coef
            lam_mat = new_sig

        # --- Track reconstruction RMSE ‖X − PᵉˣᵖʳΛᵀ‖ ---
        recon = pexpr @ lam_mat.T
        rmse_trace.append(float(np.sqrt(np.mean((x - recon) ** 2))))

    # --- Final single-cell assignment (Eq. 9), incl. gap-area recovery ---
    loglik = _morph_loglik(feats, mu, sigma)
    row = beta.sum(axis=1, keepdims=True)
    pexpr = np.divide(beta, np.where(row > _EPS, row, 1.0))
    assign = np.zeros(len(cell_xy), dtype=int)
    prob = np.zeros(len(cell_xy))
    recovered = np.zeros(len(cell_xy), dtype=bool)
    for ci in range(len(cell_xy)):
        s = spot_of_cell[ci]
        if s < 0:                       # gap-area cell: recover via nearest spot
            s = int(nearest[ci])
            recovered[ci] = True
        score = loglik[ci] + np.log(pexpr[s] + _EPS)
        assign[ci] = int(np.argmax(score))
        p = np.exp(score - score.max())
        prob[ci] = float((p / p.sum())[assign[ci]])

    names = cell_type_names
    cell_df = pd.DataFrame({
        "cell_id": cells["cell_id"].values,
        "spot": [spot_index[spot_of_cell[ci]] if spot_of_cell[ci] >= 0
                 else spot_index[int(nearest[ci])] for ci in range(len(cell_xy))],
        "x": cells["_x"].values, "y": cells["_y"].values,
        "cell_type": [names[a] for a in assign],
        "probability": prob.round(4), "recovered": recovered,
    })
    beta_df = pd.DataFrame(beta, index=spot_index, columns=names)
    prop_df = pd.DataFrame(pexpr, index=spot_index, columns=names)

    # Morphology params back in original feature units for interpretability.
    morph_rows = []
    for k, name in enumerate(names):
        for mi, feat in enumerate(feature_cols):
            morph_rows.append({
                "cell_type": name, "feature": feat,
                "mean": float(mu[k, mi] * fstd[mi] + fmean[mi]),
                "sd": float(sigma[k, mi] * fstd[mi]),
            })
    morph_df = pd.DataFrame(morph_rows)
    sig_df = pd.DataFrame(lam_mat, index=genes, columns=names) if update_signature else None
    return cell_df, beta_df, prop_df, morph_df, sig_df, rmse_trace, int(recovered.sum())


# --------------------------------------------------------------------------- #
# Public entry points
# --------------------------------------------------------------------------- #
def deconvolve(expr: pd.DataFrame, signature: pd.DataFrame, cells: pd.DataFrame,
               spot_index: list[str], spot_xy: np.ndarray, feature_cols: list[str],
               spot_diameter_fullres: float, gamma: float = 2.5, lam: float = 0.0,
               steps: int = 20) -> STIEResult:
    """Signature-based single-cell deconvolution / convolution."""
    radius = gamma * spot_diameter_fullres / 2.0
    names = [str(c) for c in signature.columns]
    out = _run_em(expr, signature, cells, spot_index, spot_xy, feature_cols,
                  radius, lam, steps, update_signature=False, cell_type_names=names)
    cell_df, beta_df, prop_df, morph_df, _sig, trace, n_rec = out
    return STIEResult("deconvolution", names, cell_df, beta_df, prop_df, morph_df,
                      None, trace, n_rec, feature_cols, gamma, lam)


def cluster(expr: pd.DataFrame, cells: pd.DataFrame, spot_index: list[str],
            spot_xy: np.ndarray, feature_cols: list[str], k: int,
            spot_diameter_fullres: float, gamma: float = 2.5, lam: float = 0.0,
            steps: int = 20, n_hvg: int = 2000) -> STIEResult:
    """Signature-free single-cell clustering (EM with signature re-estimation)."""
    radius = gamma * spot_diameter_fullres / 2.0
    # Use the most variable genes to bound the signature re-estimation cost.
    variances = expr.var(axis=0).sort_values(ascending=False)
    genes = list(variances.index[:n_hvg])
    sub = expr[genes]
    # Initialize the signature from a spot-level K-means partition (paper's init).
    init_labels = KMeans(n_clusters=k, n_init=10, random_state=1234).fit_predict(sub.to_numpy())
    names = [f"cluster{c + 1}" for c in range(k)]
    init_sig = np.vstack([sub.to_numpy()[init_labels == c].mean(axis=0) if (init_labels == c).any()
                          else sub.to_numpy().mean(axis=0) for c in range(k)]).T   # G×K
    signature = pd.DataFrame(init_sig, index=genes, columns=names)
    out = _run_em(sub, signature, cells, spot_index, spot_xy, feature_cols,
                  radius, lam, steps, update_signature=True, cell_type_names=names)
    cell_df, beta_df, prop_df, morph_df, sig_df, trace, n_rec = out
    return STIEResult("clustering", names, cell_df, beta_df, prop_df, morph_df,
                      sig_df, trace, n_rec, feature_cols, gamma, lam)
