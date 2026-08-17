import { useRef, useState } from "react";
import { api } from "../api";
import { useStore } from "../store";
import { Button, Card, Field } from "./ui";

export default function LoadTab() {
  const { sid, status, run, busy } = useStore();
  const [project, setProject] = useState("Visium_Project");
  const h5 = useRef<HTMLInputElement>(null);
  const img = useRef<HTMLInputElement>(null);
  const scale = useRef<HTMLInputElement>(null);
  const pos = useRef<HTMLInputElement>(null);
  const meta = useRef<HTMLInputElement>(null);

  const submit = async () => {
    if (!sid) return;
    const f = new FormData();
    const need = { matrix_h5: h5, tissue_image: img, scale_factors: scale, tissue_positions: pos };
    for (const [name, ref] of Object.entries(need)) {
      const file = ref.current?.files?.[0];
      if (!file) return alert(`Please choose the ${name.replace("_", " ")} file.`);
      f.append(name, file);
    }
    if (meta.current?.files?.[0]) f.append("metadata", meta.current.files[0]);
    f.append("project", project);
    await run(() => api.load(sid, f), "Visium data loaded.");
  };

  return (
    <div className="layout">
      <div className="sidebar">
        <Card title="Space Ranger input">
          <Field label="1. Filtered expression matrix (.h5)">
            <input ref={h5} type="file" accept=".h5" />
          </Field>
          <Field label="2. Tissue image (.png/.jpg)">
            <input ref={img} type="file" accept=".png,.jpg,.jpeg" />
          </Field>
          <Field label="3. Scale factors (.json)">
            <input ref={scale} type="file" accept=".json" />
          </Field>
          <Field label="4. Tissue positions (.csv)">
            <input ref={pos} type="file" accept=".csv,.tsv" />
          </Field>
          <p className="help">
            Upload each Visium file separately. Barcode suffixes (e.g. <code>-1</code>) are
            reconciled automatically, and GSM-prefixed filenames are accepted.
          </p>
          <Field label="Optional metadata CSV (barcode column)">
            <input ref={meta} type="file" accept=".csv" />
          </Field>
          <Field label="Project name">
            <input value={project} onChange={(e) => setProject(e.target.value)} />
          </Field>
          <Button className="btn-block" disabled={busy} onClick={submit}>
            {busy ? "Loading…" : "Load Visium data"}
          </Button>
        </Card>
      </div>

      <div className="content">
        <Card title="Dataset status">
          {status?.loaded ? (
            <div className="status-line">
              <span>Project: <b>{status.project}</b></span>
              <span>Spots: <b>{status.n_spots.toLocaleString()}</b></span>
              <span>Features: <b>{status.n_features.toLocaleString()}</b></span>
              <span>Image: <b>{status.has_image ? "yes" : "no"}</b></span>
              <span>PCA: <b>{status.has_pca ? "yes" : "no"}</b></span>
              <span>Clusters: <b>{status.has_clusters ? "yes" : "no"}</b></span>
            </div>
          ) : (
            <p className="muted">No dataset loaded yet.</p>
          )}
        </Card>

        {status?.has_image && (
          <Card title="Tissue image">
            <img
              src={api.imageUrl(sid!)}
              alt="tissue"
              style={{ maxWidth: "100%", borderRadius: 8, border: "1px solid var(--border)" }}
            />
          </Card>
        )}

        <Card title="Run log">
          <div className="log">{(status?.log ?? []).join("\n") || "No activity yet."}</div>
        </Card>
      </div>
    </div>
  );
}
