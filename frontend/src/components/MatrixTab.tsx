import { useState } from "react";
import { api } from "../api";
import { useStore } from "../store";
import Plot from "./Plot";
import { baseLayout, colorFor, divColorscale, plotConfig, seqColorscale } from "../theme";
import { Button, Card, Field } from "./ui";

function parseGenes(text: string): string[] {
  return text.split(/[\s,;]+/).map((g) => g.trim()).filter(Boolean);
}

export default function MatrixTab() {
  const { sid, status, run, busy, dark } = useStore();
  const [genesText, setGenesText] = useState("COL1A1, DCN, EPCAM, PECAM1, CD68, GFAP, MBP, GAD1");
  const [groupCol, setGroupCol] = useState("clusters");
  const [topN, setTopN] = useState(5);
  const [dot, setDot] = useState<any>(null);
  const [heat, setHeat] = useState<any>(null);
  const [vio, setVio] = useState<any>(null);

  const cols = status?.obs_columns?.length ? status.obs_columns : ["clusters"];

  const runDot = async () => {
    if (!sid) return;
    const r = await run(() => api.dotplot(sid, { genes: parseGenes(genesText), group_col: groupCol }));
    if (r) setDot(r);
  };
  const runHeat = async () => {
    if (!sid) return;
    const r = await run(() => api.heatmap(sid, { group_col: groupCol, top_n: topN }));
    if (r) setHeat(r);
  };
  const runVio = async () => {
    if (!sid) return;
    const r = await run(() => api.violin(sid, { genes: parseGenes(genesText).slice(0, 4), group_col: groupCol }));
    if (r) setVio(r);
  };

  return (
    <div className="layout">
      <div className="sidebar">
        <Card title="Expression matrix views">
          <Field label="Group column">
            <select value={groupCol} onChange={(e) => setGroupCol(e.target.value)}>
              {cols.map((c) => <option key={c}>{c}</option>)}
            </select>
          </Field>
          <Field label="Genes (comma / space separated)">
            <textarea rows={4} value={genesText} onChange={(e) => setGenesText(e.target.value)} />
          </Field>
          <Button className="btn-block" disabled={busy} onClick={runDot}>Dot plot</Button>
          <Button className="btn-block btn-ghost" disabled={busy} onClick={runVio}>Stacked violin (first 4 genes)</Button>
          <hr />
          <Field label="Top markers / cluster (heatmap)">
            <input type="number" value={topN} onChange={(e) => setTopN(Number(e.target.value))} />
          </Field>
          <Button className="btn-success btn-block" disabled={busy} onClick={runHeat}>Marker heatmap</Button>
          <p className="help">The marker heatmap uses the results from the Markers/DE tab; run those first.</p>
        </Card>
      </div>

      <div className="content">
        <Card title="Dot plot — mean expression (z) & fraction expressing">
          {dot ? <DotPlot data={dot} keyProp={dark ? "d" : "l"} /> : <p className="muted">Enter genes and click Dot plot.</p>}
        </Card>
        <Card title="Marker heatmap (z-scored)">
          {heat ? (
            <Plot
              data={[{
                z: heat.z, x: heat.groups, y: heat.genes, type: "heatmap",
                colorscale: divColorscale(), zmid: 0, colorbar: { title: "z" },
                hovertemplate: "%{y} · %{x}<br>z=%{z:.2f}<extra></extra>",
              }]}
              layout={baseLayout({ height: Math.max(360, heat.genes.length * 16 + 120), xaxis: { title: groupCol }, yaxis: { automargin: true } })}
              config={plotConfig}
              style={{ width: "100%" }}
              key={dark ? "hd" : "hl"}
            />
          ) : <p className="muted">Run marker analysis, then build the heatmap.</p>}
        </Card>
        <Card title="Stacked violin">
          {vio ? <StackedViolin data={vio} keyProp={dark ? "vd" : "vl"} /> : <p className="muted">Click Stacked violin to plot the first four genes.</p>}
        </Card>
      </div>
    </div>
  );
}

function DotPlot({ data, keyProp }: { data: any; keyProp: string }) {
  const x: string[] = [];
  const y: string[] = [];
  const size: number[] = [];
  const color: number[] = [];
  data.groups.forEach((g: string, gi: number) => {
    data.genes.forEach((gene: string, ni: number) => {
      x.push(g);
      y.push(gene);
      size.push(6 + data.fraction[gi][ni] * 26);
      color.push(data.mean_z[gi][ni]);
    });
  });
  return (
    <Plot
      data={[{
        x, y, mode: "markers", type: "scatter",
        marker: { size, color, colorscale: divColorscale(), cmid: 0, showscale: true, colorbar: { title: "mean z" }, line: { width: 0 } },
        hovertemplate: "%{y} · %{x}<br>z=%{marker.color:.2f}<extra></extra>",
      }]}
      layout={baseLayout({ height: Math.max(340, data.genes.length * 26 + 120), xaxis: { title: "group" }, yaxis: { automargin: true } })}
      config={plotConfig}
      style={{ width: "100%" }}
      key={keyProp}
    />
  );
}

function StackedViolin({ data, keyProp }: { data: any; keyProp: string }) {
  const traces: any[] = [];
  data.genes.forEach((gene: string, gi: number) => {
    data.groups.forEach((grp: string, i: number) => {
      const vals: number[] = data.values[gene][grp];
      traces.push({
        type: "violin", y: vals, name: grp, legendgroup: grp, showlegend: gi === 0,
        x: vals.map(() => grp), xaxis: `x${gi + 1}`, yaxis: `y${gi + 1}`,
        box: { visible: true }, meanline: { visible: true },
        line: { color: colorFor(i), width: 1 }, fillcolor: colorFor(i), opacity: 0.65,
        points: false,
      });
    });
  });
  return (
    <Plot
      data={traces}
      layout={baseLayout({
        height: 220 * data.genes.length,
        grid: { rows: data.genes.length, columns: 1, pattern: "independent" },
        annotations: data.genes.map((gene: string, i: number) => ({
          text: gene, showarrow: false, x: 0, xref: "paper",
          y: 1 - i / data.genes.length, yref: "paper", xanchor: "left", font: { size: 12 },
        })),
        margin: { t: 30, r: 16, b: 50, l: 60 },
      })}
      config={plotConfig}
      style={{ width: "100%" }}
      key={keyProp}
    />
  );
}
