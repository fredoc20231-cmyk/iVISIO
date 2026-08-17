import { useState } from "react";
import { api } from "../api";
import { useStore } from "../store";
import { Button, Card, DataTable, Field } from "./ui";
import type { TableResponse } from "../types";

export default function LigRecTab() {
  const { sid, status, run, busy } = useStore();
  const [groupCol, setGroupCol] = useState("clusters");
  const [nPerms, setNPerms] = useState(100);
  const [table, setTable] = useState<TableResponse | null>(null);

  const runLigrec = async () => {
    if (!sid) return;
    const res = await run(() => api.ligrec(sid, { group_col: groupCol, n_perms: nPerms }), "Ligand-receptor analysis complete.");
    if (res) setTable(res as TableResponse);
  };

  return (
    <div className="layout">
      <div className="sidebar">
        <Card title="Cell-cell communication">
          <p className="help">
            Ligand-receptor interaction analysis via squidpy's <code>ligrec</code> — the Python
            analog of CellChat. Computes interaction strength between spatial domains/clusters.
          </p>
          <Field label="Group by">
            <select value={groupCol} onChange={(e) => setGroupCol(e.target.value)}>
              {(status?.obs_columns ?? ["clusters"]).map((c) => <option key={c}>{c}</option>)}
            </select>
          </Field>
          <Field label="Permutations"><input type="number" value={nPerms} onChange={(e) => setNPerms(Number(e.target.value))} /></Field>
          <Button className="btn-block" disabled={busy} onClick={runLigrec}>Run ligand-receptor analysis</Button>
        </Card>
      </div>
      <div className="content">
        <Card
          title="Ligand-receptor pairs"
          actions={sid && table ? <a className="pill" href={api.downloadUrl(sid, "ligrec", status?.project ?? "visium")}>Download CSV</a> : undefined}
        >
          <DataTable data={table} maxRows={40} />
        </Card>
      </div>
    </div>
  );
}
