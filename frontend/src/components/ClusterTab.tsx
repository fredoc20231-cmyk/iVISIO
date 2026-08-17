import { useEffect, useState } from "react";
import { api } from "../api";
import { useStore } from "../store";
import Plot, { colorFor } from "./Plot";
import { Button, Card, Field } from "./ui";
import type { UmapResponse } from "../types";

export default function ClusterTab() {
  const { sid, status, run, busy } = useStore();
  const [p, setP] = useState({ n_neighbors: 30, n_pcs_low: 1, n_pcs_high: 30, resolution: 0.5, metric: "cosine" });
  const [umap, setUmap] = useState<UmapResponse | null>(null);

  const runCluster = async () => {
    if (!sid) return;
    const res = await run(() => api.cluster(sid, p), "Clustering complete.");
    if (res) loadUmap();
  };

  const loadUmap = async () => {
    if (!sid) return;
    try {
      setUmap(await api.umap(sid, "clusters"));
    } catch {
      /* not clustered yet */
    }
  };

  useEffect(() => {
    if (status?.has_umap) loadUmap();
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [status?.has_umap]);

  const traces = umap ? buildTraces(umap) : [];

  return (
    <div className="layout">
      <div className="sidebar">
        <Card title="Clustering (Leiden)">
          <Field label="First PC"><input type="number" value={p.n_pcs_low} onChange={(e) => setP({ ...p, n_pcs_low: Number(e.target.value) })} /></Field>
          <Field label="Last PC"><input type="number" value={p.n_pcs_high} onChange={(e) => setP({ ...p, n_pcs_high: Number(e.target.value) })} /></Field>
          <Field label="Resolution"><input type="number" step="0.05" value={p.resolution} onChange={(e) => setP({ ...p, resolution: Number(e.target.value) })} /></Field>
          <Field label="Neighbors"><input type="number" value={p.n_neighbors} onChange={(e) => setP({ ...p, n_neighbors: Number(e.target.value) })} /></Field>
          <Field label="UMAP metric">
            <select value={p.metric} onChange={(e) => setP({ ...p, metric: e.target.value })}>
              <option>cosine</option><option>euclidean</option><option>manhattan</option>
            </select>
          </Field>
          <Button className="btn-success btn-block" disabled={busy} onClick={runCluster}>
            Neighbors, UMAP & clusters
          </Button>
        </Card>
      </div>

      <div className="content">
        <Card title="UMAP">
          {umap ? (
            <Plot
              data={traces}
              layout={{ height: 560, margin: { t: 20, r: 10, b: 40, l: 40 }, legend: { title: { text: "cluster" } }, xaxis: { title: "UMAP1" }, yaxis: { title: "UMAP2" } }}
              style={{ width: "100%" }}
              config={{ displaylogo: false, responsive: true }}
            />
          ) : (
            <p className="muted">Run clustering to render the UMAP.</p>
          )}
        </Card>
      </div>
    </div>
  );
}

function buildTraces(u: UmapResponse) {
  const groups = Array.from(new Set(u.labels)).sort((a, b) => Number(a) - Number(b) || a.localeCompare(b));
  return groups.map((g, i) => {
    const idx = u.labels.map((l, j) => (l === g ? j : -1)).filter((j) => j >= 0);
    return {
      x: idx.map((j) => u.x[j]),
      y: idx.map((j) => u.y[j]),
      mode: "markers",
      type: "scattergl",
      name: g,
      marker: { size: 4, color: colorFor(i) },
    };
  });
}
