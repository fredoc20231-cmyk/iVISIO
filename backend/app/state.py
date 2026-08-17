"""In-memory session state for the iVisio backend.

Each browser session owns one :class:`Session` holding its AnnData object,
the tissue image + scale factor for spatial overlays, and any computed result
tables (markers, DE, spatially variable features, AI annotations, enrichment,
ligand-receptor results). Sessions are kept in a process-local dict; for a
single-instance deployment this is sufficient. For horizontal scaling, swap
this for a shared store (Redis + on-disk AnnData) behind the same interface.
"""

from __future__ import annotations

import threading
import time
import uuid
from dataclasses import dataclass, field
from typing import Any

import numpy as np
import pandas as pd


@dataclass
class Session:
    """All state for a single user's analysis session."""

    id: str
    created: float
    log: list[str] = field(default_factory=list)

    # Core data
    adata: Any | None = None                 # anndata.AnnData
    project: str = "Visium_Project"

    # Spatial overlay assets
    image: np.ndarray | None = None          # H x W x C uint8 tissue image
    scale_factor: float = 1.0                # hires/lowres scalefactor applied to pixel coords
    positions: pd.DataFrame | None = None    # barcode -> pixel coords (already scaled)

    # Result tables (pandas DataFrames) and objects
    markers: pd.DataFrame | None = None
    de: pd.DataFrame | None = None
    svg: pd.DataFrame | None = None
    reference_signature: pd.DataFrame | None = None
    reference_deconv: pd.DataFrame | None = None
    ai_annotations: pd.DataFrame | None = None
    enrichment: pd.DataFrame | None = None
    ligrec: pd.DataFrame | None = None
    signatures: dict[str, dict] = field(default_factory=dict)

    def note(self, text: str) -> None:
        stamp = time.strftime("%Y-%m-%d %H:%M:%S")
        self.log.append(f"[{stamp}] {text}")


class SessionStore:
    """Thread-safe registry of active sessions with lazy TTL eviction."""

    def __init__(self, ttl_seconds: int = 6 * 3600) -> None:
        self._sessions: dict[str, Session] = {}
        self._lock = threading.Lock()
        self._ttl = ttl_seconds

    def create(self) -> Session:
        self._evict_expired()
        sid = uuid.uuid4().hex[:16]
        sess = Session(id=sid, created=time.time())
        sess.note("Session created.")
        with self._lock:
            self._sessions[sid] = sess
        return sess

    def get(self, sid: str) -> Session:
        with self._lock:
            sess = self._sessions.get(sid)
        if sess is None:
            raise KeyError(f"Unknown or expired session '{sid}'.")
        return sess

    def drop(self, sid: str) -> None:
        with self._lock:
            self._sessions.pop(sid, None)

    def _evict_expired(self) -> None:
        cutoff = time.time() - self._ttl
        with self._lock:
            stale = [sid for sid, s in self._sessions.items() if s.created < cutoff]
            for sid in stale:
                self._sessions.pop(sid, None)


# Module-level singleton store used by the API layer.
store = SessionStore()
