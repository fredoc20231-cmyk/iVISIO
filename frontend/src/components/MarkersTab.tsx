import { useEffect, useState } from "react";
import { api } from "../api";
import { useStore } from "../store";
import { Button, Card, DataTable, Field } from "./ui";
import type { TableResponse } from "../types";

export default function MarkersTab() {
  const { sid, status, run, busy } = useStore();
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
