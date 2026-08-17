import { useRef, useState } from "react";
import { api } from "../api";
import { useStore } from "../store";
import Plot from "./Plot";
import { baseLayout, colorFor, plotConfig } from "../theme";
import { Button, Card, Field } from "./ui";

export default function IntegrationTab() {
  const { sid, run, busy, dark } = useStore();
  const zip = useRef<HTMLInputElement>(null);
  const [harmony, setHarmony] = useState(true);
  const [result, setResult] = useState<any>(null);
  const [umap, setUmap] = useState<any>(null);

  const submit = async () => {
    if (!sid || !zip.current?.files?.[0]) return alert("Choose a ZIP of Space Ranger folders.");
    const f = new FormData();
    f.append("slices_zip", zip.current.files[0]);
    f.append("use_harmony", String(harmony));
    f.append("n_pcs", "20");
    const r = await run(() => api.integrate(sid, f), "Slices integrated.");
    if (r) setResult(r);
  };

  const showUmap = async () => {
    if (!sid) return;
    const r = await run(() => api.umap(sid, "ivisio_slice"));
    if (r) setUmap(r);
  };

  const traces = umap ? buildTraces(umap) : [];

  return (
    <div className="layout">
      <div className="sidebar">
        <Card title="Multi-slice integration">
          <p className="help">
            Upload one ZIP containing two or more Space Ranger output folders (each with a
            <code> filtered_feature_bc_matrix.h5</code>). Batch effects are corrected with
            harmonypy (the Harmony analog); without it, slices are still merged.
          </p>
          <Field label="ZIP of Space Ranger folders"><input ref={zip} type="file" accept=".zip" /></Field>
          <Field label="Correct batch with Harmony">
            <input type="checkbox" checked={harmony} onChange={(e) => setHarmony(e.target.checked)} />
          </Field>
          <Button className="btn-block" disabled={busy} onClick={submit}>Load & integrate slices</Button>
          <Button className="btn-ghost btn-block" disabled={busy} onClick={showUmap}>Show integrated UMAP</Button>
        </Card>
      </div>

      <div className="content">
        {result && (
          <Card title="Integration result">
            <p className="muted">
              Slices: <b>{result.n_slices}</b> · Spots: <b>{result.n_spots?.toLocaleString()}</b> ·
              Mode: <b>{result.integration}</b>
            </p>
            <p className="help">Now run normalization/clustering on the merged object; colour the UMAP by <code>ivisio_slice</code> to inspect mixing.</p>
          </Card>
        )}
        <Card title="Integrated UMAP by slice">
          {umap ? (
            <Plot
              data={traces}
              layout={baseLayout({ height: 560, legend: { title: { text: "slice" } }, xaxis: { title: "UMAP1" }, yaxis: { title: "UMAP2" } })}
              config={plotConfig}
              style={{ width: "100%" }}
              key={dark ? "id" : "il"}
            />
          ) : <p className="muted">After clustering the merged object, show the UMAP to assess batch mixing.</p>}
        </Card>
      </div>
    </div>
  );
}

function buildTraces(u: any) {
  const groups = Array.from(new Set<string>(u.labels)).sort();
  return groups.map((g, i) => {
    const idx = u.labels.map((l: string, j: number) => (l === g ? j : -1)).filter((j: number) => j >= 0);
    return {
      x: idx.map((j: number) => u.x[j]), y: idx.map((j: number) => u.y[j]),
      mode: "markers", type: "scattergl", name: g, marker: { size: 4, color: colorFor(i) },
    };
  });
}
