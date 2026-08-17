import { useEffect, useState } from "react";
import { api } from "../api";
import { useStore } from "../store";
import Plot from "./Plot";
import { baseLayout, colorFor, plotConfig } from "../theme";
import { Card, DataTable } from "./ui";

export default function Dashboard() {
  const { sid, status, dark } = useStore();
  const [data, setData] = useState<any>(null);

  useEffect(() => {
    if (sid && status?.loaded) api.dashboard(sid).then(setData).catch(() => setData(null));
  }, [sid, status?.loaded, status?.has_clusters, status?.results?.ai]);

  if (!status?.loaded) {
    return <div className="layout"><div className="content"><Card title="Overview"><p className="muted">Load a dataset to see the overview dashboard.</p></Card></div></div>;
  }

  const tiles = [
    { label: "Spots", value: data?.n_spots?.toLocaleString() ?? "—" },
    { label: "Features", value: data?.n_features?.toLocaleString() ?? "—" },
    { label: "Clusters", value: data?.n_clusters ?? "—" },
    { label: "Median genes/spot", value: data?.median_genes ? Math.round(data.median_genes) : "—" },
    { label: "Median counts/spot", value: data?.median_counts ? Math.round(data.median_counts) : "—" },
    { label: "Median % mito", value: data?.median_percent_mt != null ? data.median_percent_mt.toFixed(1) : "—" },
  ];

  const sizes = data?.cluster_sizes ? Object.entries<number>(data.cluster_sizes) : [];
  const variance = data?.variance_ratio as number[] | undefined;

  return (
    <div style={{ padding: "1rem", display: "flex", flexDirection: "column", gap: "1rem" }}>
      <div className="kpi-row">
        {tiles.map((t) => (
          <div className="kpi" key={t.label}>
            <div className="kpi-value">{t.value}</div>
            <div className="kpi-label">{t.label}</div>
          </div>
        ))}
      </div>

      <div className="grid-2">
        <Card title="Cluster sizes">
          {sizes.length ? (
            <Plot
              data={[{
                x: sizes.map((s) => s[0]), y: sizes.map((s) => s[1]),
                type: "bar",
                marker: { color: sizes.map((_, i) => colorFor(i)) },
                hovertemplate: "Cluster %{x}<br>%{y} spots<extra></extra>",
              }]}
              layout={baseLayout({ height: 360, xaxis: { title: "Cluster" }, yaxis: { title: "Spots" } })}
              config={plotConfig}
              style={{ width: "100%" }}
            />
          ) : <p className="muted">Run clustering to see cluster sizes.</p>}
        </Card>

        <Card title="PCA variance explained">
          {variance ? (
            <Plot
              data={[{
                x: variance.map((_, i) => i + 1), y: variance,
                type: "scatter", mode: "lines+markers", fill: "tozeroy",
                line: { color: colorFor(0), width: 2 }, marker: { size: 6, color: colorFor(0) },
                hovertemplate: "PC%{x}<br>%{y:.3f}<extra></extra>",
              }]}
              layout={baseLayout({ height: 360, xaxis: { title: "Principal component" }, yaxis: { title: "Variance ratio" } })}
              config={plotConfig}
              style={{ width: "100%" }}
            />
          ) : <p className="muted">Run normalization to see the PCA scree plot.</p>}
        </Card>
      </div>

      {data?.ai_annotations && (
        <Card title="AI-predicted cell-type composition">
          <div className="grid-2">
            <Plot
              data={[{
                labels: data.ai_annotations.map((r: any) => `${r.cluster}: ${r.predicted_cell_type}`),
                values: sizes.map((s) => s[1]),
                type: "pie", hole: 0.55,
                marker: { colors: sizes.map((_, i) => colorFor(i)) },
                textinfo: "label", textposition: "outside",
              }]}
              layout={baseLayout({ height: 380, showlegend: false, margin: { t: 20, r: 10, b: 10, l: 10 } })}
              config={plotConfig}
              style={{ width: "100%" }}
              key={dark ? "d" : "l"}
            />
            <DataTable data={{ columns: Object.keys(data.ai_annotations[0]), rows: data.ai_annotations, total: data.ai_annotations.length }} maxRows={20} />
          </div>
        </Card>
      )}
    </div>
  );
}
