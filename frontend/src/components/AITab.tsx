import { useEffect, useState } from "react";
import { api } from "../api";
import { useStore } from "../store";
import { Button, Card, DataTable, Field } from "./ui";
import type { TableResponse } from "../types";

export default function AITab() {
  const { sid, status, run, busy } = useStore();
  const [topN, setTopN] = useState(10);
  const [ai, setAi] = useState<TableResponse | null>(null);

  const [source, setSource] = useState("all");
  const [cluster, setCluster] = useState("");
  const [geneSet, setGeneSet] = useState("GO_Biological_Process_2021");
  const [organism, setOrganism] = useState("human");
  const [clusters, setClusters] = useState<string[]>([]);
  const [enrich, setEnrich] = useState<TableResponse | null>(null);

  useEffect(() => {
    if (sid && status?.results?.markers) api.groups(sid, "clusters").then(setClusters).catch(() => {});
  }, [sid, status?.results?.markers]);

  const runAi = async () => {
    if (!sid) return;
    const res = await run(() => api.aiAnnotate(sid, { top_n: topN }), "AI annotation complete.");
    if (res) setAi(res as TableResponse);
  };

  const runEnrich = async () => {
    if (!sid) return;
    const res = await run(
      () => api.enrichment(sid, { source, cluster, gene_set: geneSet, organism }),
      "Enrichment complete."
    );
    if (res) setEnrich(res as TableResponse);
  };

  return (
    <div className="layout">
      <div className="sidebar">
        <Card title="AI heuristic annotation">
          <p className="help">
            Maps each cluster's top markers onto known cell-type signatures by Jaccard
            similarity — fully local, no external API calls. Run markers first.
          </p>
          <Field label="Top markers per cluster"><input type="number" value={topN} onChange={(e) => setTopN(Number(e.target.value))} /></Field>
          <Button className="btn-danger btn-block" disabled={busy} onClick={runAi}>Run AI annotation engine</Button>
        </Card>
        <Card title="Pathway enrichment">
          <p className="help">GO/pathway enrichment via gseapy Enrichr (needs network access at run time).</p>
          <Field label="Gene list source">
            <select value={source} onChange={(e) => setSource(e.target.value)}>
              <option value="all">All current markers</option>
              <option value="cluster">Selected cluster only</option>
            </select>
          </Field>
          {source === "cluster" && (
            <Field label="Cluster">
              <select value={cluster} onChange={(e) => setCluster(e.target.value)}>{clusters.map((c) => <option key={c}>{c}</option>)}</select>
            </Field>
          )}
          <Field label="Gene set">
            <select value={geneSet} onChange={(e) => setGeneSet(e.target.value)}>
              <option>GO_Biological_Process_2021</option>
              <option>GO_Molecular_Function_2021</option>
              <option>GO_Cellular_Component_2021</option>
              <option>KEGG_2021_Human</option>
            </select>
          </Field>
          <Field label="Organism">
            <select value={organism} onChange={(e) => setOrganism(e.target.value)}>
              <option value="human">human</option><option value="mouse">mouse</option>
            </select>
          </Field>
          <Button className="btn-warn btn-block" disabled={busy} onClick={runEnrich}>Run pathway enrichment</Button>
        </Card>
      </div>

      <div className="content">
        <Card
          title="AI predicted cell types"
          actions={sid && ai ? <a className="pill" href={api.downloadUrl(sid, "ai", status?.project ?? "visium")}>Download CSV</a> : undefined}
        >
          <DataTable data={ai} maxRows={30} />
        </Card>
        <Card
          title="Top enriched pathways"
          actions={sid && enrich ? <a className="pill" href={api.downloadUrl(sid, "enrichment", status?.project ?? "visium")}>Download CSV</a> : undefined}
        >
          <DataTable data={enrich} maxRows={30} />
        </Card>
      </div>
    </div>
  );
}
