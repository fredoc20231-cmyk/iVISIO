import { useState } from "react";
import { api } from "../api";
import { useStore } from "../store";
import Plot from "./Plot";
import { Button, Card, Field } from "./ui";

export default function NormalizeTab() {
  const { sid, run, busy } = useStore();
  const [method, setMethod] = useState("log");
  const [nTop, setNTop] = useState(2000);
  const [nPcs, setNPcs] = useState(20);
  const [info, setInfo] = useState<any>(null);

  const runNorm = async () => {
    if (!sid) return;
    const res = await run(
      () => api.normalize(sid, { method, n_top_genes: nTop, n_pcs: nPcs }),
      "Normalization + PCA complete."
    );
    if (res) setInfo(res);
  };

  return (
    <div className="layout">
      <div className="sidebar">
        <Card title="Normalization">
          <Field label="Method">
            <select value={method} onChange={(e) => setMethod(e.target.value)}>
              <option value="log">LogNormalize</option>
              <option value="sct">Pearson residuals (SCT analog)</option>
            </select>
          </Field>
          <Field label="Highly variable genes">
            <input type="number" value={nTop} onChange={(e) => setNTop(Number(e.target.value))} />
          </Field>
          <Field label="Principal components">
            <input type="number" value={nPcs} onChange={(e) => setNPcs(Number(e.target.value))} />
          </Field>
          <p className="help">
            Pearson residuals (scanpy's analytic normalization) are the closest analog to
            Seurat's SCTransform.
          </p>
          <Button className="btn-success btn-block" disabled={busy} onClick={runNorm}>
            Normalize and run PCA
          </Button>
        </Card>
      </div>

      <div className="content">
        <Card title="PCA variance (elbow)">
          {info?.variance_ratio ? (
            <Plot
              data={[{
                x: info.variance_ratio.map((_: number, i: number) => i + 1),
                y: info.variance_ratio,
                type: "scatter", mode: "lines+markers", marker: { color: "#2166ac" },
              }]}
              layout={{ height: 420, margin: { t: 20, r: 10, b: 45, l: 55 }, xaxis: { title: "PC" }, yaxis: { title: "Variance ratio" } }}
              style={{ width: "100%" }}
              config={{ displaylogo: false, responsive: true }}
            />
          ) : (
            <p className="muted">Run normalization to view the PCA elbow.</p>
          )}
        </Card>
        {info && (
          <Card title="Result">
            <p className="muted">
              Method: <b>{info.method}</b> · PCs: <b>{info.n_pcs}</b> · HVGs: <b>{info.n_hvg}</b>
            </p>
          </Card>
        )}
      </div>
    </div>
  );
}
