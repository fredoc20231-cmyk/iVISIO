"""Unit tests for deconvolution math and AI annotation."""

import numpy as np
import pandas as pd

from app import pipeline
from app.knowledge import AI_CELL_TYPE_DICT


def test_nnls_deconvolution_simplex():
    rng = np.random.default_rng(0)
    n_genes, n_types, n_spots = 50, 3, 40
    signature = np.abs(rng.normal(1, 0.5, (n_genes, n_types)))
    true_prop = rng.dirichlet(np.ones(n_types), size=n_spots)
    st_expr = true_prop @ signature.T
    props = pipeline.nnls_deconvolution(st_expr, signature, max_iter=500)
    assert props.shape == (n_spots, n_types)
    assert np.allclose(props.sum(axis=1), 1.0, atol=1e-4)
    # Recovered proportions should correlate with the truth.
    corr = np.corrcoef(props.ravel(), true_prop.ravel())[0, 1]
    assert corr > 0.7


def test_ai_annotate_matches_known_signature():
    # A cluster whose top markers are canonical fibroblast genes should map to CAFs.
    caf = AI_CELL_TYPE_DICT["Fibroblasts / CAFs"]
    markers = pd.DataFrame({
        "feature": caf + ["RANDOM1", "RANDOM2"],
        "cluster": ["0"] * (len(caf) + 2),
        "avg_log2FC": [3.0] * len(caf) + [0.1, 0.1],
    })
    out = pipeline.ai_annotate(markers, top_n=len(caf))
    assert out.loc[0, "predicted_cell_type"] == "Fibroblasts / CAFs"
    assert out.loc[0, "confidence_pct"] > 0


def test_build_reference_signature_means():
    expr = pd.DataFrame(
        {"c1": [1.0, 3.0], "c2": [1.0, 3.0], "c3": [5.0, 7.0]},
        index=["GENE1", "GENE2"],
    )
    sig = pipeline.build_reference_signature(expr, ["A", "A", "B"])
    assert list(sig.columns) == ["A", "B"]
    assert np.isclose(sig.loc["GENE1", "A"], 1.0)
    assert np.isclose(sig.loc["GENE2", "B"], 7.0)
