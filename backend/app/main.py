"""iVisio FastAPI application.

Exposes the full Visium spatial-omics workflow as a JSON/REST API consumed by
the React frontend: load -> QC -> normalize -> cluster -> spatial maps ->
markers/DE -> signatures -> spatially variable -> reference deconvolution ->
AI annotation -> pathway enrichment -> cell-cell communication -> exports.

Design notes
------------
* Heavy analysis steps run in a threadpool (``run_in_threadpool``) so a single
  worker stays responsive while scanpy crunches.
* All result tables are downloadable as CSV; the AnnData object and a summary
  HTML report are downloadable too, satisfying the "download of relevant data
  and files" requirement.
* Optional features degrade with a clear 400 message instead of a 500 crash.
"""

from __future__ import annotations

import io
import tempfile

import numpy as np
import pandas as pd
import scanpy as sc
from fastapi import Body, FastAPI, File, Form, HTTPException, UploadFile
from fastapi.concurrency import run_in_threadpool
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse, Response, StreamingResponse
from PIL import Image
from scipy.sparse import issparse

from . import pipeline
from .io_visium import attach_metadata, load_visium
from .knowledge import AI_CELL_TYPE_DICT, REFERENCE_CATALOG
from .schemas import (
    AIParams, ClusterParams, DEParams, EnrichmentParams, FilterParams,
    LigRecParams, MarkerParams, NormalizeParams, QCParams, SVGParams,
    SignatureParams,
)
from .state import Session, store

app = FastAPI(title="iVisio Spatial Omics Platform", version="2.0.0")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def _session(sid: str) -> Session:
    try:
        return store.get(sid)
    except KeyError as exc:
        raise HTTPException(status_code=404, detail=str(exc))


def _require_adata(sess: Session):
    if sess.adata is None:
        raise HTTPException(status_code=400, detail="Load a Visium dataset first.")
    return sess.adata


def _guard(fn, sess: Session, label: str):
    """Run a compute closure, logging success/failure and mapping errors to 400."""
    try:
        result = fn()
        return result
    except HTTPException:
        raise
    except Exception as exc:  # surface analysis errors as clean 400s
        sess.note(f"{label} ERROR: {exc}")
        raise HTTPException(status_code=400, detail=str(exc))


def _save_upload(upload: UploadFile, suffix: str) -> str:
    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=suffix)
    tmp.write(upload.file.read())
    tmp.close()
    return tmp.name


def _df_response(df: pd.DataFrame, limit: int = 2000) -> dict:
    out = df.head(limit).replace({np.nan: None})
    return {"columns": list(out.columns), "rows": out.to_dict(orient="records"),
            "total": int(df.shape[0])}


# --------------------------------------------------------------------------- #
# Session lifecycle
# --------------------------------------------------------------------------- #
@app.get("/api/health")
def health() -> dict:
    return {"status": "ok", "version": app.version}


@app.post("/api/session")
def create_session() -> dict:
    sess = store.create()
    return {"session_id": sess.id}


@app.get("/api/reference-catalog")
def reference_catalog() -> dict:
    return {"catalog": REFERENCE_CATALOG}


@app.get("/api/ai-dictionary")
def ai_dictionary() -> dict:
    return {"dictionary": AI_CELL_TYPE_DICT}


@app.get("/api/{sid}/status")
def status(sid: str) -> dict:
    sess = _session(sid)
    a = sess.adata
    return {
        "loaded": a is not None,
        "project": sess.project,
        "n_spots": int(a.n_obs) if a is not None else 0,
        "n_features": int(a.n_vars) if a is not None else 0,
        "has_pca": a is not None and "X_pca" in a.obsm,
        "has_clusters": a is not None and "clusters" in a.obs,
        "has_umap": a is not None and "X_umap" in a.obsm,
        "has_image": sess.image is not None,
        "obs_columns": list(a.obs.columns) if a is not None else [],
        "results": {
            "markers": sess.markers is not None,
            "de": sess.de is not None,
            "svg": sess.svg is not None,
            "ai": sess.ai_annotations is not None,
            "enrichment": sess.enrichment is not None,
            "ligrec": sess.ligrec is not None,
            "deconv": sess.reference_deconv is not None,
        },
        "log": sess.log[-200:],
    }


# --------------------------------------------------------------------------- #
# 1. Load
# --------------------------------------------------------------------------- #
@app.post("/api/{sid}/load")
async def load(
    sid: str,
    matrix_h5: UploadFile = File(...),
    tissue_image: UploadFile = File(...),
    scale_factors: UploadFile = File(...),
    tissue_positions: UploadFile = File(...),
    metadata: UploadFile | None = File(None),
    project: str = Form("Visium_Project"),
) -> dict:
    sess = _session(sid)

    def _do():
        h5 = _save_upload(matrix_h5, ".h5")
        img = _save_upload(tissue_image, "." + tissue_image.filename.rsplit(".", 1)[-1])
        # Ensure the loader picks the hires scale factor for an arbitrary image name.
        if "hires" not in img.lower() and "lowres" not in img.lower():
            import os
            hires = img + "_hires.png"
            os.rename(img, hires)
            img_path = hires
        else:
            img_path = img
        scale = _save_upload(scale_factors, ".json")
        pos = _save_upload(tissue_positions, ".csv")
        adata, image, sf, overlay = load_visium(h5, img_path, scale, pos, project=project)
        if metadata is not None:
            md = _save_upload(metadata, ".csv")
            adata = attach_metadata(adata, md)
        sess.adata = adata
        sess.image = image
        sess.scale_factor = sf
        sess.positions = overlay
        sess.project = project
        sess.note(f"Loaded {adata.n_obs} spots and {adata.n_vars} features.")
        return {"n_spots": int(adata.n_obs), "n_features": int(adata.n_vars)}

    return _guard(_do, sess, "LOAD")


@app.get("/api/{sid}/image")
def tissue_image(sid: str):
    sess = _session(sid)
    if sess.image is None:
        raise HTTPException(status_code=404, detail="No tissue image loaded.")
    buf = io.BytesIO()
    Image.fromarray(sess.image).save(buf, format="PNG")
    buf.seek(0)
    return Response(content=buf.read(), media_type="image/png")


@app.get("/api/{sid}/features")
def features(sid: str, q: str = "", limit: int = 50) -> dict:
    a = _require_adata(_session(sid))
    names = a.var_names
    if q:
        ql = q.upper()
        names = [n for n in names if ql in n.upper()]
    return {"features": list(names[:limit])}


# --------------------------------------------------------------------------- #
# 2. QC
# --------------------------------------------------------------------------- #
@app.post("/api/{sid}/qc")
def qc(sid: str, params: QCParams = Body(default=QCParams())) -> dict:
    sess = _session(sid)
    a = _require_adata(sess)

    def _do():
        summary = pipeline.compute_qc(a, params.mito_prefix)
        sess.note("QC metrics calculated.")
        obs = a.obs[["n_counts", "n_genes", "percent_mt"]].copy()
        return {"summary": summary,
                "distributions": {
                    "n_counts": obs["n_counts"].tolist(),
                    "n_genes": obs["n_genes"].tolist(),
                    "percent_mt": obs["percent_mt"].tolist(),
                }}

    return _guard(_do, sess, "QC")


@app.post("/api/{sid}/qc/filter")
def qc_filter(sid: str, params: FilterParams) -> dict:
    sess = _session(sid)
    a = _require_adata(sess)

    def _do():
        sess.adata = pipeline.filter_spots(a, params.min_counts, params.max_counts,
                                            params.min_genes, params.max_percent_mt)
        sess.note(f"QC filters applied; {sess.adata.n_obs} spots remain.")
        return {"n_spots": int(sess.adata.n_obs)}

    return _guard(_do, sess, "FILTER")


# --------------------------------------------------------------------------- #
# 3. Normalize
# --------------------------------------------------------------------------- #
@app.post("/api/{sid}/normalize")
async def normalize(sid: str, params: NormalizeParams) -> dict:
    sess = _session(sid)
    a = _require_adata(sess)

    def _do():
        info = pipeline.normalize(a, params.method, params.n_top_genes, params.n_pcs)
        sess.note(f"{info['method']} normalization + PCA done ({info['n_pcs']} PCs).")
        return info

    return await run_in_threadpool(_guard, _do, sess, "NORMALIZE")


# --------------------------------------------------------------------------- #
# 4. Cluster
# --------------------------------------------------------------------------- #
@app.post("/api/{sid}/cluster")
async def cluster(sid: str, params: ClusterParams) -> dict:
    sess = _session(sid)
    a = _require_adata(sess)

    def _do():
        info = pipeline.cluster(a, params.n_neighbors, params.n_pcs_low,
                                params.n_pcs_high, params.resolution, params.metric)
        sess.note(f"Clustering done: {info['n_clusters']} clusters (Leiden).")
        return info

    return await run_in_threadpool(_guard, _do, sess, "CLUSTER")


@app.get("/api/{sid}/umap")
def umap(sid: str, color: str = "clusters") -> dict:
    sess = _session(sid)
    a = _require_adata(sess)
    if "X_umap" not in a.obsm:
        raise HTTPException(status_code=400, detail="Run clustering first.")
    coords = a.obsm["X_umap"]
    if color not in a.obs:
        color = "clusters" if "clusters" in a.obs else a.obs.columns[0]
    return {
        "x": coords[:, 0].tolist(),
        "y": coords[:, 1].tolist(),
        "labels": a.obs[color].astype(str).tolist(),
        "color_by": color,
    }


# --------------------------------------------------------------------------- #
# 5. Spatial map
# --------------------------------------------------------------------------- #
@app.get("/api/{sid}/spatial")
def spatial(sid: str, mode: str = "expression", feature: str = "", group: str = "clusters") -> dict:
    sess = _session(sid)
    a = _require_adata(sess)
    if sess.positions is None:
        raise HTTPException(status_code=400, detail="Spatial coordinates unavailable; reload the Visium files.")
    pos = sess.positions
    common = [b for b in a.obs_names if b in pos.index]
    p = pos.loc[common]
    payload = {"x": p["x"].tolist(), "y": p["y"].tolist(), "barcodes": common,
               "image_shape": a.uns["ivisio"]["image_shape"], "mode": mode}
    if mode == "labels":
        col = group if group in a.obs else "clusters"
        if col not in a.obs:
            raise HTTPException(status_code=400, detail="No cluster/label column; run clustering.")
        payload["labels"] = a.obs.loc[common, col].astype(str).tolist()
        payload["legend"] = col
    else:
        if feature not in a.var_names:
            raise HTTPException(status_code=400, detail=f"Feature '{feature}' not found.")
        src = a.raw.to_adata() if a.raw is not None else a
        vals = src[common, feature].X
        vals = vals.toarray().ravel() if issparse(vals) else np.asarray(vals).ravel()
        payload["values"] = vals.tolist()
        payload["legend"] = feature
    return payload


# --------------------------------------------------------------------------- #
# 6. Explore (feature on UMAP) + signature scoring
# --------------------------------------------------------------------------- #
@app.get("/api/{sid}/feature")
def feature_on_umap(sid: str, feature: str, reduction: str = "umap") -> dict:
    a = _require_adata(_session(sid))
    key = "X_umap" if reduction == "umap" else "X_pca"
    if key not in a.obsm:
        raise HTTPException(status_code=400, detail=f"Run the step producing {reduction} first.")
    if feature not in a.var_names:
        raise HTTPException(status_code=400, detail=f"Feature '{feature}' not found.")
    coords = a.obsm[key]
    src = a.raw.to_adata() if a.raw is not None else a
    vals = src[:, feature].X
    vals = vals.toarray().ravel() if issparse(vals) else np.asarray(vals).ravel()
    return {"x": coords[:, 0].tolist(), "y": coords[:, 1].tolist(),
            "values": vals.tolist(), "feature": feature}


@app.post("/api/{sid}/signature")
def signature(sid: str, params: SignatureParams) -> dict:
    sess = _session(sid)
    a = _require_adata(sess)

    def _do():
        info = pipeline.score_signature(a, params.name, params.genes)
        sess.signatures[params.name] = info
        sess.note(f"Signature scored: {params.name} ({len(info['genes'])} genes).")
        return info

    return _guard(_do, sess, "SIGNATURE")


# --------------------------------------------------------------------------- #
# 7. Markers / DE
# --------------------------------------------------------------------------- #
@app.post("/api/{sid}/markers")
async def markers(sid: str, params: MarkerParams) -> dict:
    sess = _session(sid)
    a = _require_adata(sess)

    def _do():
        df = pipeline.find_markers(a, params.group_col, params.method,
                                   params.min_pct, params.logfc, params.top_n)
        sess.markers = df
        sess.note(f"Marker analysis done: {df.shape[0]} rows.")
        return _df_response(df)

    return await run_in_threadpool(_guard, _do, sess, "MARKERS")


@app.post("/api/{sid}/de")
async def de(sid: str, params: DEParams) -> dict:
    sess = _session(sid)
    a = _require_adata(sess)

    def _do():
        df = pipeline.pairwise_de(a, params.group_col, params.group1, params.group2,
                                  params.method, params.min_pct, params.logfc)
        sess.de = df
        sess.note(f"Pairwise DE done: {df.shape[0]} rows.")
        return _df_response(df)

    return await run_in_threadpool(_guard, _do, sess, "DE")


@app.get("/api/{sid}/groups")
def groups(sid: str, group_col: str = "clusters") -> dict:
    a = _require_adata(_session(sid))
    if group_col not in a.obs:
        return {"groups": []}
    vals = sorted(set(a.obs[group_col].astype(str)))
    return {"groups": vals}


# --------------------------------------------------------------------------- #
# 8. Spatially variable
# --------------------------------------------------------------------------- #
@app.post("/api/{sid}/svg")
async def svg(sid: str, params: SVGParams) -> dict:
    sess = _session(sid)
    a = _require_adata(sess)

    def _do():
        df = pipeline.spatially_variable(a, params.n_features, params.n_neighbors)
        sess.svg = df
        sess.note(f"Moran's I computed for {df.shape[0]} features.")
        return _df_response(df)

    return await run_in_threadpool(_guard, _do, sess, "SVG")


# --------------------------------------------------------------------------- #
# 9. Reference deconvolution
# --------------------------------------------------------------------------- #
@app.post("/api/{sid}/reference/load")
def load_reference(sid: str, signature_csv: UploadFile = File(...)) -> dict:
    sess = _session(sid)

    def _do():
        path = _save_upload(signature_csv, ".csv")
        df = pd.read_csv(path, index_col=0)
        df = df.apply(pd.to_numeric, errors="coerce").dropna(how="all")
        if df.isna().any().any():
            raise ValueError("Signature CSV must be numeric after its first (gene) column.")
        sess.reference_signature = df
        sess.note(f"Reference signature loaded: {df.shape[0]} genes x {df.shape[1]} cell types.")
        return {"n_genes": int(df.shape[0]), "cell_types": list(df.columns)}

    return _guard(_do, sess, "REFERENCE")


@app.post("/api/{sid}/reference/deconv")
async def reference_deconv(sid: str) -> dict:
    sess = _session(sid)
    a = _require_adata(sess)
    if sess.reference_signature is None:
        raise HTTPException(status_code=400, detail="Load a reference signature first.")

    def _do():
        df = pipeline.run_reference_deconvolution(a, sess.reference_signature)
        sess.reference_deconv = df
        # Attach proportions to obs so they can be viewed spatially.
        for col in df.columns:
            a.obs[f"Reference_{col}"] = df[col].values
        sess.note(f"Deconvolution done: {df.shape[0]} spots x {df.shape[1]} cell types.")
        comp = df.mean(axis=0).sort_values(ascending=False)
        return {"cell_types": list(df.columns),
                "mean_composition": {k: float(v) for k, v in comp.items()},
                "preview": _df_response(df.reset_index().rename(columns={"index": "barcode"}), 500)}

    return await run_in_threadpool(_guard, _do, sess, "DECONV")


# --------------------------------------------------------------------------- #
# 10. AI annotation + 11. enrichment + 12. ligand-receptor
# --------------------------------------------------------------------------- #
@app.post("/api/{sid}/ai-annotate")
def ai_annotate(sid: str, params: AIParams = Body(default=AIParams())) -> dict:
    sess = _session(sid)

    def _do():
        df = pipeline.ai_annotate(sess.markers, params.top_n)
        sess.ai_annotations = df
        sess.note(f"AI annotation done for {df.shape[0]} clusters.")
        return _df_response(df)

    return _guard(_do, sess, "AI")


@app.post("/api/{sid}/enrichment")
async def enrichment(sid: str, params: EnrichmentParams) -> dict:
    sess = _session(sid)
    if sess.markers is None:
        raise HTTPException(status_code=400, detail="Run marker analysis first.")

    def _do():
        m = sess.markers
        gene_col = "feature" if "feature" in m.columns else "gene"
        if params.source == "cluster" and params.cluster:
            genes = m.loc[m["cluster"].astype(str) == str(params.cluster), gene_col].tolist()
        else:
            genes = m[gene_col].tolist()
        df = pipeline.pathway_enrichment(genes, params.gene_set, params.organism)
        sess.enrichment = df
        sess.note(f"Pathway enrichment done: {df.shape[0]} terms.")
        return _df_response(df)

    return await run_in_threadpool(_guard, _do, sess, "ENRICHMENT")


@app.post("/api/{sid}/ligrec")
async def ligrec(sid: str, params: LigRecParams) -> dict:
    sess = _session(sid)
    a = _require_adata(sess)

    def _do():
        df = pipeline.ligand_receptor(a, params.group_col, params.n_perms)
        sess.ligrec = df
        sess.note(f"Ligand-receptor analysis done: {df.shape[0]} pairs.")
        return _df_response(df)

    return await run_in_threadpool(_guard, _do, sess, "LIGREC")


# --------------------------------------------------------------------------- #
# 13. Exports / downloads
# --------------------------------------------------------------------------- #
def _csv_stream(df: pd.DataFrame, filename: str) -> StreamingResponse:
    buf = io.StringIO()
    df.to_csv(buf)
    buf.seek(0)
    return StreamingResponse(iter([buf.getvalue()]), media_type="text/csv",
                             headers={"Content-Disposition": f"attachment; filename={filename}"})


@app.get("/api/{sid}/download/{artifact}")
def download(sid: str, artifact: str, prefix: str = "visium_analysis"):
    sess = _session(sid)
    tables = {
        "markers": sess.markers, "de": sess.de, "svg": sess.svg,
        "ai": sess.ai_annotations, "enrichment": sess.enrichment,
        "ligrec": sess.ligrec, "deconv": sess.reference_deconv,
        "signature": sess.reference_signature,
    }
    if artifact in tables:
        df = tables[artifact]
        if df is None:
            raise HTTPException(status_code=400, detail=f"No '{artifact}' result available yet.")
        return _csv_stream(df, f"{prefix}.{artifact}.csv")

    if artifact == "metadata":
        a = _require_adata(sess)
        return _csv_stream(a.obs, f"{prefix}.spot_metadata.csv")

    if artifact == "adata":
        a = _require_adata(sess)
        tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".h5ad")
        a.write_h5ad(tmp.name)
        with open(tmp.name, "rb") as fh:
            data = fh.read()
        return Response(content=data, media_type="application/octet-stream",
                        headers={"Content-Disposition": f"attachment; filename={prefix}.h5ad"})

    if artifact == "report":
        return Response(content=_build_report(sess), media_type="text/html",
                        headers={"Content-Disposition": f"attachment; filename={prefix}_report.html"})

    raise HTTPException(status_code=404, detail=f"Unknown artifact '{artifact}'.")


def _build_report(sess: Session) -> str:
    a = sess.adata
    rows = "".join(f"<li>{line}</li>" for line in sess.log[-100:])
    summary = "No dataset loaded."
    if a is not None:
        summary = (f"Project: {sess.project}<br>Spots: {a.n_obs}<br>"
                   f"Features: {a.n_vars}<br>Clusters: "
                   f"{a.obs['clusters'].nunique() if 'clusters' in a.obs else 'n/a'}")
    return f"""<!doctype html><html><head><meta charset="utf-8">
<title>iVisio Spatial Omics Report</title>
<style>body{{font-family:system-ui,sans-serif;margin:2rem;color:#1a2733}}
h1{{color:#2166AC}} .card{{border:1px solid #e2e8f0;border-radius:10px;padding:1rem 1.25rem;margin:1rem 0}}</style>
</head><body><h1>iVisio Spatial Omics Report</h1>
<div class="card"><h2>Summary</h2><p>{summary}</p></div>
<div class="card"><h2>Run log</h2><ol>{rows}</ol></div>
</body></html>"""
