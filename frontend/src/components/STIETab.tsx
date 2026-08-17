import { useRef, useState } from "react";
import { api } from "../api";
import { useStore } from "../store";
import Plot from "./Plot";
import { baseLayout, colorFor, plotConfig } from "../theme";
import { Button, Card, DataTable, Field } from "./ui";
import type { TableResponse } from "../types";

// STIE — single-cell deconvolution/convolution & clustering (Zhu et al.,
// Nat Commun 2024). Integrates spot expression with matched histology
// nuclear morphology via an EM algorithm.
export default function STIETab() {
  const { sid, run, busy, dark } = useStore();
  const cellsRef = useRef<HTMLInputElement>(null);
  const sigRef = useRef<HTMLInputElement>(null);
  const [mode, setMode] = useState<"deconvolve" | "cluster">("deconvolve");
  const [features, setFeatures] = useState("area,perimeter,circularity");
  const [gamma, setGamma] = useState(2.5);
  const [lam, setLam] = useState(0);
  const [steps, setSteps] = useState(20);
  const [k, setK] = useState(5);

  const [summary, setSummary] = useState<any>(null);
  const [cells, setCells] = useState<any>(null);
  const [morph, setMorph] = useState<TableResponse | null>(null);

  const runStie = async () => {
    if (!sid || !cellsRef.current?.files?.[0]) return alert("Upload a cells-on-image CSV (cell_id, pixel_x, pixel_y, morphology features).");
    const f = new FormData();
    f.append("cells_csv", cellsRef.current.files[0]);
    f.append("features", features);
    f.append("gamma", String(gamma));
    f.append("lam", String(lam));
    f.append("steps", String(steps));
    let res;
    if (mode === "deconvolve") {
      if (!sigRef.current?.files?.[0]) return alert("Deconvolution needs a gene × cell-type signature CSV.");
      f.append("signature_csv", sigRef.current.files[0]);
      res = await run(() => api.stieDeconvolve(sid, f), "STIE deconvolution complete.");
    } else {
      f.append("k", String(k));
      res = await run(() => api.stieCluster(sid, f), "STIE clustering complete.");
    }
    if (res) {
      setSummary(res);
      setCells(await api.stieCells(sid));
      setMorph(await api.stieMorphology(sid));
    }
  };

  return (
    <div className="layout">
      <div className="sidebar">
        <Card title="STIE — single-cell level">
          <p className="help">
            Integrates spot expression with histology nuclear morphology via an EM algorithm
            to reach the single-cell <em>level</em> (Zhu et al., Nat Commun 2024). Upload a
            nucleus-segmentation table (cell_id, pixel_x, pixel_y + morphology columns).
          </p>
          <Field label="Mode">
            <select value={mode} onChange={(e) => setMode(e.target.value as any)}>
              <option value="deconvolve">Deconvolution / convolution (signature)</option>
              <option value="cluster">Signature-free clustering</option>
            </select>
          </Field>
          <Field label="Cells-on-image CSV"><input ref={cellsRef} type="file" accept=".csv" /></Field>
          {mode === "deconvolve" ? (
            <Field label="Signature CSV (genes × cell types)"><input ref={sigRef} type="file" accept=".csv" /></Field>
          ) : (
            <Field label="Number of clusters (k)"><input type="number" value={k} onChange={(e) => setK(Number(e.target.value))} /></Field>
          )}
          <Field label="Morphology features (comma separated)">
            <input value={features} onChange={(e) => setFeatures(e.target.value)} />
          </Field>
          <Field label="γ · bona-fide spot area (× diameter)">
            <input type="number" step="0.5" value={gamma} onChange={(e) => setGamma(Number(e.target.value))} />
          </Field>
          <Field label="λ · morphology shrinkage penalty">
            <input type="number" step="100" value={lam} onChange={(e) => setLam(Number(e.target.value))} />
          </Field>
          <Field label="EM iterations"><input type="number" value={steps} onChange={(e) => setSteps(Number(e.target.value))} /></Field>
          <Button className="btn-danger btn-block" disabled={busy} onClick={runStie}>
            {busy ? "Running EM…" : "Run STIE"}
          </Button>
          <p className="help">
            γ default 2.5× (paper's bona-fide spot size). λ=0 still uses morphology in cell typing;
            larger λ weights morphology more in the proportion fit.
          </p>
        </Card>
      </div>

      <div className="content">
        {summary && (
          <div className="kpi-row" style={{ gridTemplateColumns: "repeat(4,1fr)" }}>
            <div className="kpi"><div className="kpi-value">{summary.n_cells?.toLocaleString()}</div><div className="kpi-label">Single cells</div></div>
            <div className="kpi"><div className="kpi-value">{summary.n_recovered?.toLocaleString()}</div><div className="kpi-label">Recovered (gap area)</div></div>
            <div className="kpi"><div className="kpi-value">{summary.cell_types?.length}</div><div className="kpi-label">{summary.mode === "clustering" ? "Clusters" : "Cell types"}</div></div>
            <div className="kpi"><div className="kpi-value">{summary.rmse_trace?.length}</div><div className="kpi-label">EM iterations</div></div>
          </div>
        )}

        <Card
          title="Single cells on tissue (STIE)"
          actions={sid && summary ? <a className="pill" href={api.downloadUrl(sid, "stie_cells", "stie")}>Download cells CSV</a> : undefined}
        >
          {cells ? <CellOverlay sid={sid!} cells={cells} keyProp={dark ? "d" : "l"} /> : <p className="muted">Run STIE to map single cells onto the histology image.</p>}
        </Card>

        <div className="grid-2">
          <Card title="EM convergence (reconstruction RMSE)">
            {summary?.rmse_trace ? (
              <Plot
                data={[{
                  x: summary.rmse_trace.map((_: number, i: number) => i + 1), y: summary.rmse_trace,
                  type: "scatter", mode: "lines+markers", line: { color: colorFor(0), width: 2 },
                }]}
                layout={baseLayout({ height: 320, xaxis: { title: "EM iteration" }, yaxis: { title: "RMSE" } })}
                config={plotConfig}
                style={{ width: "100%" }}
                key={dark ? "rd" : "rl"}
              />
            ) : <p className="muted">Convergence trace appears after a run.</p>}
          </Card>
          <Card title="Cell-type composition">
            {summary?.composition ? <Composition comp={summary.composition} keyProp={dark ? "cd" : "cl"} /> : <p className="muted">—</p>}
          </Card>
        </div>

        <Card
          title="Learned nuclear-morphology profiles"
          actions={sid && morph ? <a className="pill" href={api.downloadUrl(sid, "stie_morphology", "stie")}>Download CSV</a> : undefined}
        >
          {morph ? <MorphologyChart data={morph} keyProp={dark ? "md" : "ml"} /> : <p className="muted">STIE learns a Gaussian morphology profile per cell type (Fig. 3j/m of the paper).</p>}
        </Card>

        {summary?.mode === "clustering" && sid && (
          <Card
            title="Cluster average gene expression (CAGE)"
            actions={<a className="pill" href={api.downloadUrl(sid, "stie_signature", "stie")}>Download CAGE CSV</a>}
          >
            <p className="muted">The single-cell-level clustering signature is downloadable as CAGE — compare it against scRNA-seq signatures as in the paper.</p>
          </Card>
        )}
      </div>
    </div>
  );
}

function CellOverlay({ sid, cells, keyProp }: { sid: string; cells: any; keyProp: string }) {
  const [h, w] = cells.image_shape;
  const groups = Array.from(new Set<string>(cells.labels)).sort();
  const traces = groups.map((g, i) => {
    const idx = cells.labels.map((l: string, j: number) => (l === g ? j : -1)).filter((j: number) => j >= 0);
    return {
      x: idx.map((j: number) => cells.x[j]), y: idx.map((j: number) => cells.y[j]),
      text: idx.map((j: number) => `${cells.cell_ids[j]}${cells.recovered[j] ? " (recovered)" : ""}`),
      mode: "markers", type: "scattergl", name: g,
      marker: { size: 5, color: colorFor(i), line: { width: 0 } },
      hovertemplate: "%{text}<br>" + g + "<extra></extra>",
    };
  });
  return (
    <Plot
      data={traces}
      layout={baseLayout({
        height: 700,
        xaxis: { visible: false, range: [0, w] },
        yaxis: { visible: false, range: [0, h], scaleanchor: "x", scaleratio: 1 },
        images: [{ source: api.imageUrl(sid), xref: "x", yref: "y", x: 0, y: h, sizex: w, sizey: h, sizing: "stretch", layer: "below", opacity: 0.85 }],
        legend: { orientation: "v" },
        margin: { t: 10, r: 10, b: 10, l: 10 },
      })}
      config={plotConfig}
      style={{ width: "100%" }}
      key={keyProp}
    />
  );
}

function Composition({ comp, keyProp }: { comp: Record<string, number>; keyProp: string }) {
  const entries = Object.entries(comp).sort((a, b) => b[1] - a[1]);
  return (
    <Plot
      data={[{ x: entries.map((e) => e[0]), y: entries.map((e) => e[1]), type: "bar", marker: { color: entries.map((_, i) => colorFor(i)) } }]}
      layout={baseLayout({ height: 320, xaxis: { tickangle: -35 }, yaxis: { title: "Cells" } })}
      config={plotConfig}
      style={{ width: "100%" }}
      key={keyProp}
    />
  );
}

function MorphologyChart({ data, keyProp }: { data: TableResponse; keyProp: string }) {
  // Grouped bars of mean ± sd per (cell_type, feature).
  const rows = data.rows as any[];
  const features = Array.from(new Set(rows.map((r) => r.feature)));
  const types = Array.from(new Set(rows.map((r) => r.cell_type)));
  const traces = types.map((t, i) => {
    const sub = features.map((f) => rows.find((r) => r.cell_type === t && r.feature === f));
    return {
      x: features, y: sub.map((r) => (r ? r.mean : 0)),
      error_y: { type: "data", array: sub.map((r) => (r ? r.sd : 0)), visible: true },
      type: "bar", name: String(t), marker: { color: colorFor(i) },
    };
  });
  return (
    <Plot
      data={traces}
      layout={baseLayout({ height: 380, barmode: "group", xaxis: { title: "Nuclear feature" }, yaxis: { title: "Mean ± sd" } })}
      config={plotConfig}
      style={{ width: "100%" }}
      key={keyProp}
    />
  );
}
