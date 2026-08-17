import { useState } from "react";
import { api } from "../api";
import { useStore } from "../store";
import { Card, Field } from "./ui";

const ARTIFACTS: { key: string; label: string; needs?: string }[] = [
  { key: "adata", label: "AnnData object (.h5ad)" },
  { key: "metadata", label: "Spot metadata (CSV)" },
  { key: "markers", label: "Marker table (CSV)", needs: "markers" },
  { key: "de", label: "Pairwise DE (CSV)", needs: "de" },
  { key: "svg", label: "Spatially variable (CSV)", needs: "svg" },
  { key: "deconv", label: "Deconvolution (CSV)", needs: "deconv" },
  { key: "signature", label: "Reference signature (CSV)" },
  { key: "ai", label: "AI annotations (CSV)", needs: "ai" },
  { key: "enrichment", label: "Pathway enrichment (CSV)", needs: "enrichment" },
  { key: "ligrec", label: "Ligand-receptor pairs (CSV)", needs: "ligrec" },
  { key: "report", label: "HTML report" },
];

export default function ExportTab() {
  const { sid, status } = useStore();
  const [prefix, setPrefix] = useState("visium_analysis");

  return (
    <div className="layout">
      <div className="sidebar">
        <Card title="Export options">
          <Field label="Output filename prefix">
            <input value={prefix} onChange={(e) => setPrefix(e.target.value)} />
          </Field>
          <p className="help">
            Every result computed in this session is downloadable here. Greyed-out items
            become available once their analysis has been run.
          </p>
        </Card>
      </div>
      <div className="content">
        <Card title="Downloads">
          <div className="download-grid">
            {ARTIFACTS.map((a) => {
              const ready = !a.needs || status?.results?.[a.needs] || status?.loaded;
              const enabled = sid && status?.loaded && ready;
              return enabled ? (
                <a key={a.key} className="btn btn-ghost" href={api.downloadUrl(sid!, a.key, prefix)}>
                  ⬇ {a.label}
                </a>
              ) : (
                <button key={a.key} className="btn btn-ghost" disabled>
                  {a.label}
                </button>
              );
            })}
          </div>
        </Card>
        <Card title="Export status">
          <pre className="log">
            {status?.loaded
              ? [
                  `Prefix: ${prefix}`,
                  `Spots: ${status.n_spots}`,
                  `Markers: ${status.results.markers}`,
                  `DE: ${status.results.de}`,
                  `SVG: ${status.results.svg}`,
                  `AI: ${status.results.ai}`,
                  `Enrichment: ${status.results.enrichment}`,
                  `Ligand-receptor: ${status.results.ligrec}`,
                  `Deconvolution: ${status.results.deconv}`,
                ].join("\n")
              : "Load a dataset to enable exports."}
          </pre>
        </Card>
      </div>
    </div>
  );
}
