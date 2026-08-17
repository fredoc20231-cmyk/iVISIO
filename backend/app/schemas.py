"""Pydantic request models for the iVisio API."""

from __future__ import annotations

from pydantic import BaseModel, Field


class QCParams(BaseModel):
    mito_prefix: str = "MT-"


class FilterParams(BaseModel):
    min_counts: float = 500
    max_counts: float = 100000
    min_genes: float = 100
    max_percent_mt: float = 25


class NormalizeParams(BaseModel):
    method: str = Field("log", pattern="^(log|sct)$")
    n_top_genes: int = 2000
    n_pcs: int = 20


class ClusterParams(BaseModel):
    n_neighbors: int = 30
    n_pcs_low: int = 1
    n_pcs_high: int = 30
    resolution: float = 0.5
    metric: str = "cosine"


class MarkerParams(BaseModel):
    group_col: str = "clusters"
    method: str = "wilcoxon"
    min_pct: float = 0.1
    logfc: float = 0.25
    top_n: int = 0


class DEParams(BaseModel):
    group_col: str = "clusters"
    group1: str
    group2: str
    method: str = "wilcoxon"
    min_pct: float = 0.1
    logfc: float = 0.25


class SignatureParams(BaseModel):
    name: str = "My_signature"
    genes: list[str]


class SVGParams(BaseModel):
    n_features: int = 100
    n_neighbors: int = 6


class AIParams(BaseModel):
    top_n: int = 10


class EnrichmentParams(BaseModel):
    source: str = "all"           # "all" or "cluster"
    cluster: str | None = None
    gene_set: str = "GO_Biological_Process_2021"
    organism: str = "human"


class LigRecParams(BaseModel):
    group_col: str = "clusters"
    n_perms: int = 100
