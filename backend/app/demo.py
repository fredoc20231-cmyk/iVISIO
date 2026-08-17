"""Synthetic Visium dataset generator for demos and tests.

Produces a self-consistent toy 10x Visium library — spot grid, expression
matrix with a known spatial cell-type pattern, a matched H&E-like image, a
nucleus-segmentation table with morphology signal, and a gene x cell-type
signature — so the whole platform (including STIE) can be exercised end-to-end
without uploading real data. Deterministic given a seed, so it also backs the
test suite.
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd


def generate(seed: int = 0, grid: int = 24, n_types: int = 4,
             genes_per_type: int = 30, noise_genes: int = 20):
    """Build a synthetic Visium dataset.

    Returns a dict with: ``adata`` (spots x genes counts), ``image`` (HxWx3
    uint8), ``overlay`` (scaled spot coords), ``spots_fullres`` (full-res spot
    centres), ``spot_diameter_fullres``, ``scale_factor``, ``cells`` (nucleus
    table with morphology), and ``signature`` (genes x cell types).
    """
    rng = np.random.default_rng(seed)

    # --- Spot grid in full-resolution pixel space ---
    spacing = 100.0
    spot_diameter = 60.0
    margin = 80.0
    xs = np.arange(grid) * spacing + margin
    ys = np.arange(grid) * spacing + margin
    gx, gy = np.meshgrid(xs, ys)
    spot_x = gx.ravel()
    spot_y = gy.ravel()
    n_spots = spot_x.size
    barcodes = [f"S{ i }-1" for i in range(n_spots)]

    # --- Spatial domains: quadrant-based dominant cell type ---
    cx, cy = margin + grid * spacing / 2, margin + grid * spacing / 2
    domain = ((spot_x > cx).astype(int) + 2 * (spot_y > cy).astype(int))
    domain = domain % n_types                                    # 0..n_types-1

    # --- Gene signature: block-structured marker expression ---
    n_marker = genes_per_type * n_types
    n_genes = n_marker + noise_genes
    gene_names = [f"GENE{i:04d}" for i in range(n_genes)]
    signature = np.full((n_genes, n_types), 0.5)
    for k in range(n_types):
        block = slice(k * genes_per_type, (k + 1) * genes_per_type)
        signature[block, k] = 6.0
    signature += rng.normal(0, 0.1, signature.shape).clip(-0.4, 0.4)
    signature = signature.clip(0.01, None)

    # --- Spot proportions: dominant domain type + contamination ---
    prop = np.full((n_spots, n_types), 0.08)
    prop[np.arange(n_spots), domain] = 1.0
    prop += rng.normal(0, 0.03, prop.shape).clip(0, None)
    prop /= prop.sum(axis=1, keepdims=True)

    # --- Counts ~ Poisson(depth * proportion @ signature^T) ---
    depth = 200.0
    mean_expr = (prop @ signature.T) * depth / signature.sum()
    counts = rng.poisson(mean_expr).astype(np.float32)

    adata = ad.AnnData(X=counts)
    adata.obs_names = barcodes
    adata.var_names = gene_names
    adata.obsm["spatial"] = np.column_stack([spot_x, spot_y])
    adata.obs["barcode"] = barcodes

    # --- Nuclei: cells scattered within each spot, morphology by type ---
    # Distinct nuclear-size means per type give STIE a morphology signal.
    type_area = np.linspace(40, 120, n_types)
    cell_rows = []
    radius = spot_diameter / 2
    for si in range(n_spots):
        n_cells = int(rng.poisson(4) + 1)
        types_here = rng.choice(n_types, size=n_cells, p=prop[si])
        for j, t in enumerate(types_here):
            ang = rng.uniform(0, 2 * np.pi)
            rr = radius * np.sqrt(rng.uniform(0, 1))
            area = max(10.0, rng.normal(type_area[t], 8))
            perim = 2 * np.sqrt(np.pi * area) * rng.normal(1.0, 0.05)
            circ = np.clip(4 * np.pi * area / (perim ** 2), 0.1, 1.0)
            cell_rows.append({
                "cell_id": f"C{si}_{j}",
                "spot": barcodes[si],
                "pixel_x": spot_x[si] + rr * np.cos(ang),
                "pixel_y": spot_y[si] + rr * np.sin(ang),
                "area": area, "perimeter": perim, "circularity": circ,
                "true_type": f"type{t + 1}",
            })
    cells = pd.DataFrame(cell_rows)

    # --- Synthetic H&E-like image with nuclei dots ---
    height = int(ys.max() + margin)
    width = int(xs.max() + margin)
    image = np.full((height, width, 3), (238, 224, 234), dtype=np.uint8)   # pale eosin
    for _, c in cells.iterrows():
        x, y = int(c["pixel_x"]), int(c["pixel_y"])
        r = max(2, int(np.sqrt(c["area"]) / 2))
        y0, y1 = max(0, y - r), min(height, y + r)
        x0, x1 = max(0, x - r), min(width, x + r)
        image[y0:y1, x0:x1] = (110, 70, 140)                               # haematoxylin nuclei

    scale_factor = 1.0
    overlay = pd.DataFrame({
        "barcode": barcodes,
        "x": spot_x * scale_factor,
        "y": height - spot_y * scale_factor,
    }).set_index("barcode")
    spots_fullres = pd.DataFrame({"barcode": barcodes, "x": spot_x, "y": spot_y}).set_index("barcode")

    adata.uns["ivisio"] = {
        "project": "Demo_Visium",
        "image_shape": [int(height), int(width)],
        "scale_factor": scale_factor,
        "spot_diameter_fullres": spot_diameter,
    }
    signature_df = pd.DataFrame(signature, index=gene_names,
                                columns=[f"type{k + 1}" for k in range(n_types)])
    return {
        "adata": adata, "image": image, "overlay": overlay,
        "spots_fullres": spots_fullres, "spot_diameter_fullres": spot_diameter,
        "scale_factor": scale_factor, "cells": cells, "signature": signature_df,
    }
