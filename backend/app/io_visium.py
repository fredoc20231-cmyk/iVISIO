"""Loading 10x Visium Space Ranger output into an AnnData object.

Mirrors the flexible loader from the original R app: it accepts the four
Space Ranger files uploaded separately (filtered `.h5` matrix, tissue image,
scale-factor JSON, tissue-positions CSV) rather than requiring a rigid
directory layout, and it tolerates the common barcode-suffix mismatches
(`-1`) seen across public datasets.
"""

from __future__ import annotations

import json
import re

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from PIL import Image


def _normalize_barcode(values: pd.Series | list[str]) -> list[str]:
    """Strip a trailing ``-<n>`` suffix so barcodes match across pipelines."""
    return [re.sub(r"-[0-9]+$", "", str(v).strip()) for v in values]


def read_positions(path: str) -> pd.DataFrame:
    """Read tissue_positions(_list).csv, header-optional, into a tidy frame."""
    raw = pd.read_csv(path, header=None)
    # Space Ranger v2 files carry a header row starting with "barcode".
    if str(raw.iloc[0, 0]).strip().lower() == "barcode":
        raw = pd.read_csv(path, header=0)
    if raw.shape[1] < 6:
        raise ValueError("Tissue-position CSV must contain at least six columns.")
    pos = raw.iloc[:, :6].copy()
    pos.columns = ["barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"]
    pos["barcode"] = pos["barcode"].astype(str)
    pos["in_tissue"] = pd.to_numeric(pos["in_tissue"], errors="coerce").fillna(0).astype(int)
    for col in ("pxl_row", "pxl_col"):
        pos[col] = pd.to_numeric(pos[col], errors="coerce")
    return pos


def read_scalefactors(path: str) -> dict:
    with open(path) as fh:
        return json.load(fh)


def load_visium(
    matrix_h5: str,
    image_path: str,
    scalefactors_path: str,
    positions_path: str,
    project: str = "Visium_Project",
) -> tuple[ad.AnnData, np.ndarray, float, pd.DataFrame]:
    """Assemble an AnnData plus spatial overlay assets from four files.

    Returns ``(adata, image_array, scale_factor, positions)`` where
    ``positions`` pixel coordinates are already multiplied by the appropriate
    scale factor and flipped into image space, and ``adata.obsm['spatial']``
    holds the scaled (x, y) coordinates for every retained spot.
    """
    adata = sc.read_10x_h5(matrix_h5)
    adata.var_names_make_unique()
    adata.obs_names = [str(b) for b in adata.obs_names]

    scale = read_scalefactors(scalefactors_path)
    # Choose the scale factor matching the supplied image resolution.
    lower = image_path.lower()
    if "hires" in lower:
        scale_factor = float(scale.get("tissue_hires_scalef", 1.0))
    elif "lowres" in lower:
        scale_factor = float(scale.get("tissue_lowres_scalef", 1.0))
    else:
        scale_factor = float(scale.get("tissue_hires_scalef") or scale.get("tissue_lowres_scalef") or 1.0)

    img = np.asarray(Image.open(image_path).convert("RGB"))
    height, width = img.shape[0], img.shape[1]

    pos = read_positions(positions_path)
    pos = pos[pos["in_tissue"] == 1].copy()
    pos = pos[np.isfinite(pos["pxl_row"]) & np.isfinite(pos["pxl_col"])]

    # Match positions to the expression matrix, tolerating a barcode-suffix mismatch.
    obs_index = list(adata.obs_names)
    common = set(pos["barcode"]) & set(obs_index)
    if not common:
        norm_obs = _normalize_barcode(obs_index)
        norm_pos = _normalize_barcode(pos["barcode"])
        obs_map = {n: o for n, o in zip(norm_obs, obs_index)}
        pos = pos.assign(_norm=norm_pos)
        pos = pos[pos["_norm"].isin(obs_map)]
        pos["barcode"] = pos["_norm"].map(obs_map)
        pos = pos.drop(columns="_norm")
        common = set(pos["barcode"]) & set(obs_index)
    if not common:
        raise ValueError(
            "No tissue-position barcodes matched the expression matrix. "
            "Check that the H5 matrix and tissue-position CSV come from the same Visium library."
        )

    pos = pos.drop_duplicates("barcode").set_index("barcode")
    keep = [b for b in obs_index if b in pos.index]
    adata = adata[keep].copy()
    pos = pos.loc[keep]

    # Scaled pixel coordinates in image space (y flipped so origin is top-left).
    x = pos["pxl_col"].to_numpy() * scale_factor
    y = pos["pxl_row"].to_numpy() * scale_factor
    adata.obsm["spatial"] = np.column_stack([x, y])

    overlay = pd.DataFrame({
        "barcode": keep,
        "x": x,
        "y": height - y,          # flip for plotting over the raster image
    }).set_index("barcode")

    adata.uns["ivisio"] = {
        "project": project,
        "image_shape": [int(height), int(width)],
        "scale_factor": scale_factor,
    }
    adata.obs["barcode"] = adata.obs_names

    return adata, img, scale_factor, overlay


def attach_metadata(adata: ad.AnnData, metadata_path: str) -> ad.AnnData:
    """Join an optional user metadata CSV (must contain a ``barcode`` column)."""
    md = pd.read_csv(metadata_path)
    if "barcode" not in md.columns:
        raise ValueError("Metadata CSV must contain a 'barcode' column.")
    md["barcode"] = md["barcode"].astype(str)
    md = md.drop_duplicates("barcode").set_index("barcode")
    common = [b for b in adata.obs_names if b in md.index]
    if not common:
        raise ValueError("No metadata barcodes matched the Visium object.")
    add = md.reindex(adata.obs_names)
    for col in add.columns:
        adata.obs[col] = add[col].values
    return adata
