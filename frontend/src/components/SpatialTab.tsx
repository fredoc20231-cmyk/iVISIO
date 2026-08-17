import { useEffect, useState } from "react";
import { api } from "../api";
import { useStore } from "../store";
import Plot, { colorFor } from "./Plot";
import { Button, Card, Field } from "./ui";
import type { SpatialResponse } from "../types";

export default function SpatialTab() {
  const { sid, status, run, busy } = useStore();
  const [mode, setMode] = useState("expression");
  const [feature, setFeature] = useState("");
  const [group, setGroup] = useState("clusters");
  const [featOptions, setFeatOptions] = useState<string[]>([]);
  const [data, setData] = useState<SpatialResponse | null>(null);
  const [size, setSize] = useState(6);

  useEffect(() => {
    if (sid && status?.loaded) api.features(sid, "").then(setFeatOptions).catch(() => {});
  }, [sid, status?.loaded]);

  const refresh = async () => {
    if (!sid) return;
    const res = await run(() =>
      api.spatial(sid, { mode, feature: feature || featOptions[0], group })
    );
    if (res) setData(res as SpatialResponse);
  };

  const searchFeatures = async (q: string) => {
    setFeature(q);
    if (sid && q.length >= 2) setFeatOptions(await api.features(sid, q));
  };

  return (
    <div className="layout">
      <div className="sidebar">
        <Card title="Spatial overlay">
          <Field label="Overlay mode">
            <select value={mode} onChange={(e) => setMode(e.target.value)}>
              <option value="expression">Expression feature</option>
              <option value="labels">Metadata clusters/labels</option>
            </select>
          </Field>
          {mode === "expression" ? (
            <Field label="Feature (gene)">
              <input list="feat-list" value={feature} onChange={(e) => searchFeatures(e.target.value)} placeholder="Type a gene…" />
              <datalist id="feat-list">{featOptions.map((f) => <option key={f} value={f} />)}</datalist>
            </Field>
          ) : (
            <Field label="Label column">
              <select value={group} onChange={(e) => setGroup(e.target.value)}>
                {(status?.obs_columns ?? ["clusters"]).map((c) => <option key={c}>{c}</option>)}
              </select>
            </Field>
          )}
          <Field label="Spot size">
            <input type="range" min={2} max={16} value={size} onChange={(e) => setSize(Number(e.target.value))} />
          </Field>
          <Button className="btn-block" disabled={busy} onClick={refresh}>Refresh spatial map</Button>
        </Card>
      </div>

      <div className="content">
        <Card title="Spatial feature map">
          {data ? <SpatialFigure sid={sid!} data={data} size={size} /> : <p className="muted">Choose a feature or label column and refresh.</p>}
        </Card>
      </div>
    </div>
  );
}

function SpatialFigure({ sid, data, size }: { sid: string; data: SpatialResponse; size: number }) {
  const [h, w] = data.image_shape;
  let traces: any[];
  if (data.mode === "labels" && data.labels) {
    const groups = Array.from(new Set(data.labels)).sort();
    traces = groups.map((g, i) => {
      const idx = data.labels!.map((l, j) => (l === g ? j : -1)).filter((j) => j >= 0);
      return {
        x: idx.map((j) => data.x[j]), y: idx.map((j) => data.y[j]),
        mode: "markers", type: "scattergl", name: g,
        marker: { size, color: colorFor(i) },
      };
    });
  } else {
    traces = [{
      x: data.x, y: data.y, mode: "markers", type: "scattergl",
      marker: { size, color: data.values, colorscale: "Magma", showscale: true, colorbar: { title: data.legend } },
      text: data.barcodes,
    }];
  }
  return (
    <Plot
      data={traces}
      layout={{
        height: 720,
        margin: { t: 10, r: 10, b: 10, l: 10 },
        xaxis: { visible: false, range: [0, w], constrain: "domain" },
        yaxis: { visible: false, range: [0, h], scaleanchor: "x", scaleratio: 1 },
        images: [{
          source: api.imageUrl(sid),
          xref: "x", yref: "y", x: 0, y: h, sizex: w, sizey: h,
          sizing: "stretch", layer: "below", opacity: 0.9,
        }],
        legend: { orientation: "v" },
      }}
      style={{ width: "100%" }}
      config={{ displaylogo: false, responsive: true }}
    />
  );
}
