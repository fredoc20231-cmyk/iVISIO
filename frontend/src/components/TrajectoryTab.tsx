import { useEffect, useState } from "react";
import { api } from "../api";
import { useStore } from "../store";
import Plot from "./Plot";
import { baseLayout, colorFor, plotConfig, seqColorscale } from "../theme";
import { Button, Card, Field } from "./ui";

export default function TrajectoryTab() {
  const { sid, status, run, busy, dark } = useStore();
  const [groupCol, setGroupCol] = useState("clusters");
  const [root, setRoot] = useState("");
  const [groups, setGroups] = useState<string[]>([]);
  const [paga, setPaga] = useState<any>(null);
  const [pt, setPt] = useState<any>(null);

  const cols = status?.obs_columns?.length ? status.obs_columns : ["clusters"];

  useEffect(() => {
    if (sid && status?.has_clusters) api.groups(sid, groupCol).then((g) => {
      setGroups(g);
      if (g.length && !root) setRoot(g[0]);
    }).catch(() => {});
  }, [sid, status?.has_clusters, groupCol]);

  const runPaga = async () => {
    if (!sid) return;
    const r = await run(() => api.paga(sid, { group_col: groupCol, top_n: 5 }), "PAGA graph computed.");
    if (r) setPaga(r);
  };
  const runPt = async () => {
    if (!sid) return;
    const r = await run(() => api.pseudotime(sid, { group_col: groupCol, root_group: root }), "Pseudotime computed.");
    if (r) setPt(r);
  };

  return (
    <div className="layout">
      <div className="sidebar">
        <Card title="Trajectory inference">
          <p className="help">PAGA connectivity + diffusion pseudotime — the scanpy analog of a lineage/pseudotime workflow.</p>
          <Field label="Group column">
            <select value={groupCol} onChange={(e) => setGroupCol(e.target.value)}>
              {cols.map((c) => <option key={c}>{c}</option>)}
            </select>
          </Field>
          <Button className="btn-block" disabled={busy} onClick={runPaga}>Compute PAGA graph</Button>
          <hr />
          <Field label="Root group (pseudotime origin)">
            <select value={root} onChange={(e) => setRoot(e.target.value)}>
              {groups.map((g) => <option key={g}>{g}</option>)}
            </select>
          </Field>
          <Button className="btn-success btn-block" disabled={busy} onClick={runPt}>Compute pseudotime</Button>
        </Card>
      </div>

      <div className="content">
        <Card title="PAGA connectivity graph">
          {paga ? <PagaGraph data={paga} keyProp={dark ? "pd" : "pl"} /> : <p className="muted">Partition-based graph abstraction: nodes are clusters, edge weight is connectivity.</p>}
        </Card>
        <Card title="Diffusion pseudotime (UMAP)">
          {pt?.umap ? (
            <Plot
              data={[{
                x: pt.umap.map((c: number[]) => c[0]), y: pt.umap.map((c: number[]) => c[1]),
                mode: "markers", type: "scattergl",
                marker: { size: 5, color: pt.pseudotime, colorscale: seqColorscale(), showscale: true, colorbar: { title: "pseudotime" } },
                hovertemplate: "dpt=%{marker.color:.3f}<extra></extra>",
              }]}
              layout={baseLayout({ height: 560, xaxis: { title: "UMAP1" }, yaxis: { title: "UMAP2" } })}
              config={plotConfig}
              style={{ width: "100%" }}
              key={dark ? "ptd" : "ptl"}
            />
          ) : <p className="muted">Choose a root group and compute pseudotime.</p>}
        </Card>
      </div>
    </div>
  );
}

function PagaGraph({ data, keyProp }: { data: any; keyProp: string }) {
  const maxW = Math.max(...data.edges.map((e: any) => e.weight), 1);
  const edgeTraces = data.edges.map((e: any) => {
    const s = data.positions[e.source];
    const t = data.positions[e.target];
    return {
      x: [s[0], t[0]], y: [s[1], t[1]], mode: "lines", type: "scatter",
      line: { width: 1 + (e.weight / maxW) * 8, color: "rgba(120,140,160,0.45)" },
      hoverinfo: "text", text: `${e.source}–${e.target}: ${e.weight.toFixed(3)}`, showlegend: false,
    };
  });
  const nodeTrace = {
    x: data.groups.map((g: string) => data.positions[g][0]),
    y: data.groups.map((g: string) => data.positions[g][1]),
    mode: "markers+text", type: "scatter",
    marker: { size: 26, color: data.groups.map((_: string, i: number) => colorFor(i)), line: { width: 2, color: "#fff" } },
    text: data.groups, textposition: "middle center", textfont: { color: "#fff", size: 11 },
    hovertemplate: "Cluster %{text}<extra></extra>", showlegend: false,
  };
  return (
    <Plot
      data={[...edgeTraces, nodeTrace]}
      layout={baseLayout({ height: 560, xaxis: { visible: false }, yaxis: { visible: false } })}
      config={plotConfig}
      style={{ width: "100%" }}
      key={keyProp}
    />
  );
}
