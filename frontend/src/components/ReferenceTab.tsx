import { useEffect, useRef, useState } from "react";
import { api } from "../api";
import { useStore } from "../store";
import Plot, { colorFor } from "./Plot";
import { Button, Card, DataTable, Field } from "./ui";
import type { TableResponse } from "../types";

export default function ReferenceTab() {
  const { sid, status, run, busy } = useStore();
  const sig = useRef<HTMLInputElement>(null);
  const [catalog, setCatalog] = useState<any[]>([]);
  const [comp, setComp] = useState<Record<string, number> | null>(null);
  const [preview, setPreview] = useState<TableResponse | null>(null);

  useEffect(() => { api.referenceCatalog().then(setCatalog).catch(() => {}); }, []);

  const loadSig = async () => {
    if (!sid || !sig.current?.files?.[0]) return alert("Choose a signature CSV first.");
    const f = new FormData();
    f.append("signature_csv", sig.current.files[0]);
    await run(() => api.loadReference(sid, f), "Reference signature loaded.");
  };

  const deconv = async () => {
    if (!sid) return;
    const res = await run(() => api.deconv(sid), "Deconvolution complete.");
    if (res) { setComp(res.mean_composition); setPreview(res.preview); }
  };

  const compData = comp
    ? Object.entries(comp).sort((a, b) => b[1] - a[1])
    : [];

  return (
    <div className="layout">
      <div className="sidebar">
        <Card title="Reference / deconvolution">
          <p className="help">
            NNLS deconvolution using a gene × cell-type signature CSV (genes in the first
            column, cell types across the remaining columns). References are never
            downloaded silently — grab them from the cited source and upload.
          </p>
          <Field label="Signature CSV"><input ref={sig} type="file" accept=".csv" /></Field>
          <Button className="btn-block" disabled={busy} onClick={loadSig}>Load signature</Button>
          <Button className="btn-success btn-block" disabled={busy} onClick={deconv}>Run reference deconvolution</Button>
        </Card>
      </div>

      <div className="content">
        {compData.length > 0 && (
          <Card title="Average composition across spots">
            <Plot
              data={[{
                x: compData.map((d) => d[0]), y: compData.map((d) => d[1] * 100),
                type: "bar", marker: { color: compData.map((_, i) => colorFor(i)) },
              }]}
              layout={{ height: 380, margin: { t: 20, r: 10, b: 90, l: 50 }, yaxis: { title: "Mean %" }, xaxis: { tickangle: -40 } }}
              style={{ width: "100%" }}
              config={{ displaylogo: false, responsive: true }}
            />
          </Card>
        )}
        <Card
          title="Per-spot proportions"
          actions={sid && preview ? <a className="pill" href={api.downloadUrl(sid, "deconv", status?.project ?? "visium")}>Download CSV</a> : undefined}
        >
          <DataTable data={preview} maxRows={25} />
        </Card>
        <Card title="Reference catalog (public sources)">
          <DataTable data={catalog.length ? { columns: Object.keys(catalog[0]), rows: catalog, total: catalog.length } : null} maxRows={20} />
        </Card>
      </div>
    </div>
  );
}
