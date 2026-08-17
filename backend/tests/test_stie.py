"""Tests for the STIE EM implementation against synthetic ground truth."""

import numpy as np
import pandas as pd

from app import demo, stie


def _inputs(seed=0):
    d = demo.generate(seed=seed)
    adata = d["adata"]
    # STIE expects log-normalized-ish expression; use log1p of counts here.
    expr = pd.DataFrame(np.log1p(adata.X), index=list(adata.obs_names),
                        columns=list(adata.var_names))
    spot_index = list(adata.obs_names)
    spot_xy = d["spots_fullres"].loc[spot_index, ["x", "y"]].to_numpy()
    cells = d["cells"].copy()
    cells["_x"] = cells["pixel_x"]
    cells["_y"] = cells["pixel_y"]
    return d, expr, cells, spot_index, spot_xy


def test_deconvolve_recovers_cell_types():
    d, expr, cells, spot_index, spot_xy = _inputs()
    res = stie.deconvolve(
        expr, d["signature"], cells, spot_index, spot_xy,
        feature_cols=["area", "perimeter", "circularity"],
        spot_diameter_fullres=d["spot_diameter_fullres"], gamma=2.5, lam=0.0, steps=15,
    )
    assert res.mode == "deconvolution"
    assert set(res.cells.columns) >= {"cell_id", "spot", "x", "y", "cell_type", "probability", "recovered"}
    assert res.cells.shape[0] == cells.shape[0]
    assert list(res.spot_prop.columns) == list(d["signature"].columns)
    # Proportions per spot are a simplex.
    row_sums = res.spot_prop.to_numpy().sum(axis=1)
    assert np.allclose(row_sums[row_sums > 0], 1.0, atol=1e-6)
    # EM should not increase reconstruction error overall.
    assert res.rmse_trace[-1] <= res.rmse_trace[0] + 1e-6
    # Cell typing should beat chance against the known truth.
    merged = res.cells.merge(cells[["cell_id", "true_type"]], on="cell_id")
    accuracy = (merged["cell_type"] == merged["true_type"]).mean()
    assert accuracy > 0.5


def test_cluster_produces_signature_and_labels():
    d, expr, cells, spot_index, spot_xy = _inputs(seed=1)
    res = stie.cluster(
        expr, cells, spot_index, spot_xy,
        feature_cols=["area", "perimeter", "circularity"], k=4,
        spot_diameter_fullres=d["spot_diameter_fullres"], gamma=2.5, lam=0.0,
        steps=10, n_hvg=200,
    )
    assert res.mode == "clustering"
    assert res.signature is not None                     # CAGE estimated
    assert res.signature.shape[1] == 4
    assert res.cells["cell_type"].nunique() <= 4
    assert res.cells.shape[0] == cells.shape[0]


def test_gap_area_recovery_flags_uncovered_cells():
    d, expr, cells, spot_index, spot_xy = _inputs(seed=2)
    # Shrink the bona-fide area so some cells fall outside every spot and must
    # be recovered via their nearest spot.
    res = stie.deconvolve(
        expr, d["signature"], cells, spot_index, spot_xy,
        feature_cols=["area", "perimeter"],
        spot_diameter_fullres=d["spot_diameter_fullres"], gamma=0.5, lam=0.0, steps=5,
    )
    assert res.n_recovered >= 0
    assert res.cells["recovered"].sum() == res.n_recovered


def test_lambda_penalty_runs():
    d, expr, cells, spot_index, spot_xy = _inputs(seed=3)
    res = stie.deconvolve(
        expr, d["signature"], cells, spot_index, spot_xy,
        feature_cols=["area", "circularity"],
        spot_diameter_fullres=d["spot_diameter_fullres"], gamma=2.5, lam=1e3, steps=8,
    )
    assert res.lam == 1e3
    assert np.isfinite(res.rmse_trace).all()
