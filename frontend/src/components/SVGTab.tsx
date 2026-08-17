import { useState } from "react";
import { api } from "../api";
import { useStore } from "../store";
import { Button, Card, DataTable, Field } from "./ui";
import type { TableResponse } from "../types";

export default function SVGTab() {
  const { sid, status, run, busy } = useStore();
  const [n, setN] = useState(100);
  const [k, setK] = useState(6);
  const [table, setTable] = useState<TableResponse | null>(null);

  const runSvg = async () => {
    if (!sid) return;
    const res = await run(() => api.svg(sid, { n_features: n, n_neighbors: k }), "Moran's I computed.");
    if (res) setTable(res as TableResponse);
  };

  return (
    <div className="layout">
      <div className="sidebar">
        <Card title="Spatially variable features">
          <p className="help">Moran's I spatial autocorrelation via squidpy — the Python analog of Seurat's markvariogram.</p>
          <Field label="Number of features"><input type="number" value={n} onChange={(e) => setN(Number(e.target.value))} /></Field>
          <Field label="Spatial neighbors (k)"><input type="number" value={k} onChange={(e) => setK(Number(e.target.value))} /></Field>
          <Button className="btn-success btn-block" disabled={busy} onClick={runSvg}>Find spatially variable features</Button>
        </Card>
      </div>
      <div className="content">
        <Card
          title="Moran's I results"
          actions={sid && table ? <a className="pill" href={api.downloadUrl(sid, "svg", status?.project ?? "visium")}>Download CSV</a> : undefined}
        >
          <DataTable data={table} maxRows={40} />
        </Card>
      </div>
    </div>
  );
}
