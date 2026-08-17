import { useState } from "react";
import { useStore } from "./store";
import LoadTab from "./components/LoadTab";
import QCTab from "./components/QCTab";
import NormalizeTab from "./components/NormalizeTab";
import ClusterTab from "./components/ClusterTab";
import SpatialTab from "./components/SpatialTab";
import ExploreTab from "./components/ExploreTab";
import MarkersTab from "./components/MarkersTab";
import SVGTab from "./components/SVGTab";
import ReferenceTab from "./components/ReferenceTab";
import AITab from "./components/AITab";
import LigRecTab from "./components/LigRecTab";
import ExportTab from "./components/ExportTab";

interface TabDef {
  id: string;
  label: string;
  render: () => JSX.Element;
  requires?: "loaded" | "pca" | "clusters";
}

const TABS: TabDef[] = [
  { id: "load", label: "1 · Load", render: () => <LoadTab /> },
  { id: "qc", label: "2 · QC", render: () => <QCTab />, requires: "loaded" },
  { id: "normalize", label: "3 · Normalize", render: () => <NormalizeTab />, requires: "loaded" },
  { id: "cluster", label: "4 · Clusters", render: () => <ClusterTab />, requires: "pca" },
  { id: "spatial", label: "5 · Spatial map", render: () => <SpatialTab />, requires: "loaded" },
  { id: "explore", label: "6 · Explore", render: () => <ExploreTab />, requires: "loaded" },
  { id: "markers", label: "7 · Markers / DE", render: () => <MarkersTab />, requires: "clusters" },
  { id: "svg", label: "8 · Spatially variable", render: () => <SVGTab />, requires: "loaded" },
  { id: "reference", label: "9 · Deconvolution", render: () => <ReferenceTab />, requires: "loaded" },
  { id: "ai", label: "10 · AI & Enrichment", render: () => <AITab />, requires: "clusters" },
  { id: "ligrec", label: "11 · Cell-Cell Comm.", render: () => <LigRecTab />, requires: "clusters" },
  { id: "export", label: "12 · Export", render: () => <ExportTab />, requires: "loaded" },
];

export default function App() {
  const { status, busy, toast } = useStore();
  const [active, setActive] = useState("load");

  const isEnabled = (t: TabDef): boolean => {
    if (!t.requires) return true;
    if (!status) return false;
    if (t.requires === "loaded") return status.loaded;
    if (t.requires === "pca") return status.has_pca;
    if (t.requires === "clusters") return status.has_clusters;
    return true;
  };

  const current = TABS.find((t) => t.id === active) ?? TABS[0];

  return (
    <div className="app">
      <div className="topbar">
        <div>
          <h1>iVisio · Spatial Omics Platform</h1>
          <div className="sub">FastAPI + scanpy backend · React + TypeScript frontend</div>
        </div>
        <div className="sub">
          {status?.loaded
            ? `${status.n_spots.toLocaleString()} spots · ${status.n_features.toLocaleString()} features`
            : "No dataset loaded"}
        </div>
      </div>

      {busy && <div className="busybar" />}

      <div className="tabs">
        {TABS.map((t) => (
          <button
            key={t.id}
            className={`tab ${active === t.id ? "active" : ""}`}
            disabled={!isEnabled(t)}
            onClick={() => setActive(t.id)}
          >
            {t.label}
          </button>
        ))}
      </div>

      {current.render()}

      {toast && <div className="toast">{toast}</div>}
    </div>
  );
}
