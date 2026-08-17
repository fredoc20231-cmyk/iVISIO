import { useState } from "react";
import { useStore } from "./store";
import Dashboard from "./components/Dashboard";
import LoadTab from "./components/LoadTab";
import QCTab from "./components/QCTab";
import NormalizeTab from "./components/NormalizeTab";
import ClusterTab from "./components/ClusterTab";
import SpatialTab from "./components/SpatialTab";
import ExploreTab from "./components/ExploreTab";
import MarkersTab from "./components/MarkersTab";
import MatrixTab from "./components/MatrixTab";
import SVGTab from "./components/SVGTab";
import SpatialStatsTab from "./components/SpatialStatsTab";
import TrajectoryTab from "./components/TrajectoryTab";
import ReferenceTab from "./components/ReferenceTab";
import STIETab from "./components/STIETab";
import AITab from "./components/AITab";
import LigRecTab from "./components/LigRecTab";
import IntegrationTab from "./components/IntegrationTab";
import ExportTab from "./components/ExportTab";

interface TabDef {
  id: string;
  label: string;
  render: () => JSX.Element;
  requires?: "loaded" | "pca" | "clusters";
}

const TABS: TabDef[] = [
  { id: "dashboard", label: "Overview", render: () => <Dashboard />, requires: "loaded" },
  { id: "load", label: "1 · Load", render: () => <LoadTab /> },
  { id: "qc", label: "2 · QC", render: () => <QCTab />, requires: "loaded" },
  { id: "normalize", label: "3 · Normalize", render: () => <NormalizeTab />, requires: "loaded" },
  { id: "cluster", label: "4 · Clusters", render: () => <ClusterTab />, requires: "pca" },
  { id: "spatial", label: "5 · Spatial map", render: () => <SpatialTab />, requires: "loaded" },
  { id: "explore", label: "6 · Explore", render: () => <ExploreTab />, requires: "loaded" },
  { id: "markers", label: "7 · Markers / DE", render: () => <MarkersTab />, requires: "clusters" },
  { id: "matrix", label: "8 · Expression matrix", render: () => <MatrixTab />, requires: "clusters" },
  { id: "svg", label: "9 · Spatially variable", render: () => <SVGTab />, requires: "loaded" },
  { id: "spstats", label: "10 · Spatial statistics", render: () => <SpatialStatsTab />, requires: "clusters" },
  { id: "trajectory", label: "11 · Trajectory", render: () => <TrajectoryTab />, requires: "clusters" },
  { id: "reference", label: "12 · Deconvolution", render: () => <ReferenceTab />, requires: "loaded" },
  { id: "stie", label: "13 · STIE (single-cell)", render: () => <STIETab />, requires: "loaded" },
  { id: "ai", label: "14 · AI & Enrichment", render: () => <AITab />, requires: "clusters" },
  { id: "ligrec", label: "15 · Cell-Cell Comm.", render: () => <LigRecTab />, requires: "clusters" },
  { id: "integration", label: "16 · Integration", render: () => <IntegrationTab />, requires: "loaded" },
  { id: "export", label: "17 · Export", render: () => <ExportTab />, requires: "loaded" },
];

export default function App() {
  const { status, busy, toast, dark, toggleDark } = useStore();
  const [active, setActive] = useState("load");

  const isEnabled = (t: TabDef): boolean => {
    if (!t.requires) return true;
    if (!status) return false;
    if (t.requires === "loaded") return status.loaded;
    if (t.requires === "pca") return status.has_pca;
    if (t.requires === "clusters") return status.has_clusters;
    return true;
  };

  const current = TABS.find((t) => t.id === active) ?? TABS[1];

  return (
    <div className="app">
      <div className="topbar">
        <div>
          <h1>iVisio · Spatial Omics Platform</h1>
          <div className="sub">FastAPI + scanpy/squidpy · React + TypeScript · advanced spatial analytics</div>
        </div>
        <div style={{ display: "flex", alignItems: "center", gap: "1rem" }}>
          <span className="sub">
            {status?.loaded
              ? `${status.n_spots.toLocaleString()} spots · ${status.n_features.toLocaleString()} features`
              : "No dataset loaded"}
          </span>
          <button className="theme-toggle" onClick={toggleDark} title="Toggle dark mode">
            {dark ? "☀ Light" : "☾ Dark"}
          </button>
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
