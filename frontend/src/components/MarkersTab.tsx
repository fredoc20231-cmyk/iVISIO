import { useEffect, useState } from "react";
import { api } from "../api";
import { useStore } from "../store";
import Plot from "./Plot";
import { baseLayout, plotConfig } from "../theme";
import { Button, Card, DataTable, Field } from "./ui";
import type { TableResponse } from "../types";

export default function MarkersTab() {
  const { sid, status, run, busy, dark } = useStore();
  const [groupCol, setGroupCol] = useState("clusters");
  const [minPct, setMinPct] = useState(0.1);
  const [logfc, setLogfc] = useState(0.25);
  const [method, setMethod] = useState("wilcoxon");
  const [groups, setGroups] = useState<string[]>([]);
  const [g1, setG1] = useState("");
  const [g2, setG2] = useState("");
  const [table, setTable] = useState<TableResponse | null>(null);

  useEffect(() => {
    if (sid && status?.has_clusters) {
      api.groups(sid, groupCol).then((gs) => {
        setGroups(gs);
        if (gs.length) { setG1(gs[0]); setG2(gs[1] ?? gs[0]); }
      });
    }
  }, [sid, status?.has_clusters, groupCol]);

  const findMarkers = async () => {
    if (!sid) return;
    const res = await run(() => api.markers(sid, { group_col: groupCol, method, min_pct: minPct, logfc }), "Markers computed.");
    if (res) setTable(res as TableResponse);
  };

  const runDe = async () => {
    if (!sid) return;
    const res = await run(() => api.de(sid, { group_col: groupCol, group1: g1, group2: g2, method, min_pct: minPct, logfc }), "Pairwise DE complete.");
    if (res) setTable(res as TableResponse);
  };

  return (
    <div className="layout">
      <div className="sidebar">
        <Card title="Markers / DE">
          <Field label="Group column">
            <select value={groupCol} onChange={(e) => setGroupCol(e.target.value)}>
              {(status?.obs_columns ?? ["clusters"]).map((c) => <option key={c}>{c}</option>)}
            </select>
          </Field>
          <Field label="Test method">
            <select value={method} onChange={(e) => setMethod(e.target.value)}>
              <option value="wilcoxon">Wilcoxon</option>
              <option value="t-test">t-test</option>
              <option value="logreg">Logistic regression</option>
            </select>
          </Field>
          <Field label="Min fraction expressing"><input type="number" step="0.05" value={minPct} onChange={(e) => setMinPct(Number(e.target.value))} /></Field>
          <Field label="Log2FC threshold"><input type="number" step="0.05" value={logfc} onChange={(e) => setLogfc(Number(e.target.value))} /></Field>
          <Button className="btn-success btn-block" disabled={busy} onClick={findMarkers}>Find markers for all groups</Button>
          <hr />
          <Field label="Group 1">
            <select value={g1} onChange={(e) => setG1(e.target.value)}>{groups.map((g) => <option key={g}>{g}</option>)}</select>
          </Field>
          <Field label="Group 2">
            <select value={g2} onChange={(e) => setG2(e.target.value)}>{groups.map((g) => <option key={g}>{g}</option>)}</select>
          </Field>
          <Button className="btn-warn btn-block" disabled={busy} onClick={runDe}>Run pairwise DE</Button>
        </Card>
      </div>
      <div className="content">
        {table && <Card title="Volcano plot"><Volcano table={table} dark={dark} /></Card>}
        <Card
          title="Differential results"
          actions={sid && table ? <a className="pill" href={api.downloadUrl(sid, "markers", status?.project ?? "visium")}>Download CSV</a> : undefined}
        >
          <DataTable data={table} maxRows={40} />
        </Card>
      </div>
    </div>
  );
}

function Volcano({ table, dark }: { table: TableResponse; dark: boolean }) {
  const rows = table.rows;
  const lfc = rows.map((r) => Number(r["avg_log2FC"] ?? 0));
  const padj = rows.map((r) => Number(r["p_val_adj"] ?? 1));
  const y = padj.map((p) => -Math.log10(Math.max(p, 1e-300)));
  const names = rows.map((r) => String(r["feature"] ?? r["gene"] ?? ""));
  // Significance/direction categories drive color AND are labeled in the legend,
  // so identity never rests on hue alone.
  const cat = lfc.map((f, i) =>
    padj[i] < 0.05 && Math.abs(f) >= 1 ? (f > 0 ? "Up" : "Down") : "n.s."
  );
  const groups: { name: string; color: string }[] = [
    { name: "Up", color: "#e34948" },
    { name: "Down", color: "#2a78d6" },
    { name: "n.s.", color: dark ? "#6b7686" : "#b8beca" },
  ];
  const traces = groups.map((g) => {
    const idx = cat.map((c, i) => (c === g.name ? i : -1)).filter((i) => i >= 0);
    return {
      x: idx.map((i) => lfc[i]), y: idx.map((i) => y[i]),
      text: idx.map((i) => names[i]), mode: "markers", type: "scattergl", name: g.name,
      marker: { size: 6, color: g.color, opacity: 0.8 },
      hovertemplate: "%{text}<br>log2FC=%{x:.2f}<br>-log10 p.adj=%{y:.1f}<extra></extra>",
    };
  });
  return (
    <Plot
      data={traces}
      layout={baseLayout({
        height: 460,
        xaxis: { title: "log2 fold change" },
        yaxis: { title: "-log10 adjusted p-value" },
        shapes: [
          { type: "line", x0: 1, x1: 1, yref: "paper", y0: 0, y1: 1, line: { dash: "dot", width: 1, color: "#98a2b3" } },
          { type: "line", x0: -1, x1: -1, yref: "paper", y0: 0, y1: 1, line: { dash: "dot", width: 1, color: "#98a2b3" } },
        ],
      })}
      config={plotConfig}
      style={{ width: "100%" }}
      key={dark ? "vd" : "vl"}
    />
  );
}
