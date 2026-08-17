import { useEffect, useState } from "react";
import { api } from "../api";
import { useStore } from "../store";
import Plot from "./Plot";
import { Button, Card, Field } from "./ui";
import type { FeatureResponse } from "../types";

export default function ExploreTab() {
  const { sid, status, run, busy } = useStore();
  const [feature, setFeature] = useState("");
  const [options, setOptions] = useState<string[]>([]);
  const [reduction, setReduction] = useState("umap");
  const [data, setData] = useState<FeatureResponse | null>(null);

  const [sigName, setSigName] = useState("My_signature");
  const [sigGenes, setSigGenes] = useState("COL1A1\nCOL1A2\nDCN");

  useEffect(() => {
    if (sid && status?.loaded) api.features(sid, "").then(setOptions).catch(() => {});
  }, [sid, status?.loaded]);

  const plot = async () => {
    if (!sid) return;
    const res = await run(() => api.feature(sid, feature || options[0], reduction));
    if (res) setData(res as FeatureResponse);
  };

  const score = async () => {
    if (!sid) return;
    const genes = sigGenes.split(/[\s,;]+/).map((g) => g.trim()).filter(Boolean);
    await run(() => api.signature(sid, { name: sigName, genes }), "Signature scored.");
  };

  return (
    <div className="layout">
      <div className="sidebar">
        <Card title="Feature expression">
          <Field label="Feature">
            <input list="explore-feats" value={feature} onChange={(e) => setFeature(e.target.value)} placeholder="Type a gene…" />
            <datalist id="explore-feats">{options.map((f) => <option key={f} value={f} />)}</datalist>
          </Field>
          <Field label="Reduction">
            <select value={reduction} onChange={(e) => setReduction(e.target.value)}>
              <option value="umap">UMAP</option><option value="pca">PCA</option>
            </select>
          </Field>
          <Button className="btn-block" disabled={busy} onClick={plot}>Plot feature</Button>
        </Card>
        <Card title="Signature scoring">
          <Field label="Signature name"><input value={sigName} onChange={(e) => setSigName(e.target.value)} /></Field>
          <Field label="Genes (one per line / comma separated)">
            <textarea rows={5} value={sigGenes} onChange={(e) => setSigGenes(e.target.value)} />
          </Field>
          <Button className="btn-success btn-block" disabled={busy} onClick={score}>Score signature</Button>
        </Card>
      </div>

      <div className="content">
        <Card title={data ? `Expression: ${data.feature}` : "Feature expression"}>
          {data ? (
            <Plot
              data={[{
                x: data.x, y: data.y, mode: "markers", type: "scattergl",
                marker: { size: 5, color: data.values, colorscale: "Viridis", showscale: true },
              }]}
              layout={{ height: 620, margin: { t: 20, r: 10, b: 40, l: 40 }, xaxis: { title: `${reduction}1` }, yaxis: { title: `${reduction}2` } }}
              style={{ width: "100%" }}
              config={{ displaylogo: false, responsive: true }}
            />
          ) : (
            <p className="muted">Plot a feature to view its expression.</p>
          )}
        </Card>
      </div>
    </div>
  );
}
