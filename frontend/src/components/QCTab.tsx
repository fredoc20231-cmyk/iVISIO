import { useState } from "react";
import { api } from "../api";
import { useStore } from "../store";
import Plot from "./Plot";
import { Button, Card, Field } from "./ui";

export default function QCTab() {
  const { sid, run, busy } = useStore();
  const [mito, setMito] = useState("MT-");
  const [dist, setDist] = useState<any>(null);
  const [filter, setFilter] = useState({ min_counts: 500, max_counts: 100000, min_genes: 100, max_percent_mt: 25 });

  const compute = async () => {
    if (!sid) return;
    const res = await run(() => api.qc(sid, { mito_prefix: mito }), "QC metrics computed.");
    if (res) setDist(res.distributions);
  };

  const applyFilter = async () => {
    if (!sid) return;
    await run(() => api.qcFilter(sid, filter), "QC filters applied.");
    await compute();
  };

  return (
    <div className="layout">
      <div className="sidebar">
        <Card title="QC parameters">
          <Field label="Mitochondrial gene prefix">
            <input value={mito} onChange={(e) => setMito(e.target.value)} />
          </Field>
          <Button className="btn-success btn-block" disabled={busy} onClick={compute}>
            Calculate QC metrics
          </Button>
          <hr />
          {(["min_counts", "max_counts", "min_genes", "max_percent_mt"] as const).map((k) => (
            <Field key={k} label={k.replace(/_/g, " ")}>
              <input
                type="number"
                value={(filter as any)[k]}
                onChange={(e) => setFilter({ ...filter, [k]: Number(e.target.value) })}
              />
            </Field>
          ))}
          <Button className="btn-warn btn-block" disabled={busy} onClick={applyFilter}>
            Apply QC filters
          </Button>
        </Card>
      </div>

      <div className="content">
        <Card title="QC distributions">
          {dist ? (
            <Plot
              data={[
                { x: dist.n_counts, type: "histogram", name: "counts", marker: { color: "#2166ac" }, xaxis: "x1", yaxis: "y1" },
                { x: dist.n_genes, type: "histogram", name: "genes", marker: { color: "#2f9e44" }, xaxis: "x2", yaxis: "y2" },
                { x: dist.percent_mt, type: "histogram", name: "% mito", marker: { color: "#e8590c" }, xaxis: "x3", yaxis: "y3" },
              ]}
              layout={{
                grid: { rows: 1, columns: 3, pattern: "independent" },
                height: 340,
                margin: { t: 20, r: 10, b: 40, l: 40 },
                showlegend: false,
              }}
              style={{ width: "100%" }}
              config={{ displaylogo: false, responsive: true }}
            />
          ) : (
            <p className="muted">Run QC metrics to see distributions.</p>
          )}
        </Card>

        {dist && (
          <Card title="Counts vs genes per spot">
            <Plot
              data={[{
                x: dist.n_counts, y: dist.n_genes, mode: "markers", type: "scattergl",
                marker: { size: 4, color: dist.percent_mt, colorscale: "Viridis", showscale: true, colorbar: { title: "% mito" } },
              }]}
              layout={{ height: 420, margin: { t: 20, r: 10, b: 45, l: 55 }, xaxis: { title: "Total counts" }, yaxis: { title: "Detected genes" } }}
              style={{ width: "100%" }}
              config={{ displaylogo: false, responsive: true }}
            />
          </Card>
        )}
      </div>
    </div>
  );
}
