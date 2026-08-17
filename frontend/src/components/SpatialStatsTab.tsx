import { useState } from "react";
import { api } from "../api";
import { useStore } from "../store";
import Plot from "./Plot";
import { baseLayout, colorFor, divColorscale, plotConfig } from "../theme";
import { Button, Card, Field } from "./ui";

export default function SpatialStatsTab() {
  const { sid, status, run, busy, dark } = useStore();
  const [groupCol, setGroupCol] = useState("clusters");
  const [k, setK] = useState(6);
  const [nhood, setNhood] = useState<any>(null);
  const [cooc, setCooc] = useState<any>(null);
  const [corr, setCorr] = useState<any>(null);

  const cols = status?.obs_columns?.length ? status.obs_columns : ["clusters"];

  const runNhood = async () => {
    if (!sid) return;
    const r = await run(() => api.nhood(sid, { group_col: groupCol, n_neighbors: k }), "Neighborhood enrichment done.");
    if (r) setNhood(r);
  };
  const runCooc = async () => {
    if (!sid) return;
    const r = await run(() => api.cooccurrence(sid, { group_col: groupCol, n_neighbors: k }), "Co-occurrence done.");
    if (r) setCooc(r);
  };
  const runCorr = async () => {
    if (!sid) return;
    const r = await run(() => api.dendrogram(sid, { group_col: groupCol, top_n: 5 }));
    if (r) setCorr(r);
  };

  return (
    <div className="layout">
      <div className="sidebar">
        <Card title="Spatial statistics">
          <p className="help">Squidpy graph statistics over the tissue neighbor graph — the spatial analytical core.</p>
          <Field label="Group column">
            <select value={groupCol} onChange={(e) => setGroupCol(e.target.value)}>
              {cols.map((c) => <option key={c}>{c}</option>)}
            </select>
          </Field>
          <Field label="Spatial neighbors (k)"><input type="number" value={k} onChange={(e) => setK(Number(e.target.value))} /></Field>
          <Button className="btn-block" disabled={busy} onClick={runNhood}>Neighborhood enrichment</Button>
          <Button className="btn-block btn-ghost" disabled={busy} onClick={runCooc}>Co-occurrence vs distance</Button>
          <Button className="btn-success btn-block" disabled={busy} onClick={runCorr}>Cluster correlation</Button>
        </Card>
      </div>

      <div className="content">
        <Card title="Neighborhood enrichment (z-score)">
          {nhood ? (
            <Plot
              data={[{
                z: nhood.zscore, x: nhood.groups, y: nhood.groups, type: "heatmap",
                colorscale: divColorscale(), zmid: 0, colorbar: { title: "z" },
                hovertemplate: "%{y} ↔ %{x}<br>z=%{z:.1f}<extra></extra>",
              }]}
              layout={baseLayout({ height: 460, xaxis: { title: groupCol }, yaxis: { title: groupCol, automargin: true } })}
              config={plotConfig}
              style={{ width: "100%" }}
              key={dark ? "nd" : "nl"}
            />
          ) : <p className="muted">Enriched (red) / depleted (blue) spatial adjacency between groups.</p>}
        </Card>

        <Card title="Co-occurrence probability vs distance">
          {cooc ? (
            <Plot
              data={cooc.groups.map((g: string, i: number) => ({
                x: cooc.distance,
                // occ[i][i] = self co-occurrence ratio across distance bins.
                y: cooc.occ[i][i],
                type: "scatter", mode: "lines", name: g,
                line: { color: colorFor(i), width: 2 },
              }))}
              layout={baseLayout({ height: 460, xaxis: { title: "Distance" }, yaxis: { title: "Co-occurrence ratio" } })}
              config={plotConfig}
              style={{ width: "100%" }}
              key={dark ? "cd" : "cl"}
            />
          ) : <p className="muted">Self co-occurrence ratio of each group as a function of distance.</p>}
        </Card>

        <Card title="Cluster–cluster correlation">
          {corr ? (
            <Plot
              data={[{
                z: corr.correlation, x: corr.groups, y: corr.groups, type: "heatmap",
                colorscale: divColorscale(), zmid: 0, colorbar: { title: "r" },
                hovertemplate: "%{y} ↔ %{x}<br>r=%{z:.2f}<extra></extra>",
              }]}
              layout={baseLayout({ height: 440, xaxis: { title: groupCol }, yaxis: { title: groupCol, automargin: true } })}
              config={plotConfig}
              style={{ width: "100%" }}
              key={dark ? "rd" : "rl"}
            />
          ) : <p className="muted">Transcriptomic similarity between clusters (dendrogram input).</p>}
        </Card>
      </div>
    </div>
  );
}
