// Central visualization theme for iVisio.
//
// The categorical palette is the dataviz reference theme, validated with the
// skill's validate_palette.js in BOTH modes (light: all checks pass with the
// relief rule — we always ship legends + table views; dark: all checks pass).
// Assign hues in fixed order, never cycled: a 9th series folds to "Other".

const CATEGORICAL_LIGHT = [
  "#2a78d6", "#eb6834", "#1baf7a", "#eda100",
  "#e87ba4", "#008300", "#4a3aa7", "#e34948",
];
const CATEGORICAL_DARK = [
  "#3987e5", "#d95926", "#199e70", "#c98500",
  "#d55181", "#008300", "#9085e9", "#e66767",
];

// Sequential blue ramp (magnitude — heatmaps, expression) light→dark.
const SEQ_STOPS = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"];

// Diverging blue↔red with a neutral gray midpoint (z-scores, enrichment).
const DIV_LIGHT: [number, string][] = [
  [0, "#184f95"], [0.25, "#6da7ec"], [0.5, "#f0efec"], [0.75, "#eb6834"], [1, "#c92a2a"],
];
const DIV_DARK: [number, string][] = [
  [0, "#3987e5"], [0.25, "#256abf"], [0.5, "#383835"], [0.75, "#d95926"], [1, "#e66767"],
];

const INK = {
  light: { text: "#0b0b0b", muted: "#898781", grid: "#e1e0d9", surface: "#ffffff" },
  dark: { text: "#f4f6f8", muted: "#a7b0bb", grid: "#2c2c2a", surface: "#141a22" },
};

let _dark = false;

export function setDarkMode(d: boolean): void {
  _dark = d;
}
export function isDark(): boolean {
  return _dark;
}
export function palette(): string[] {
  return _dark ? CATEGORICAL_DARK : CATEGORICAL_LIGHT;
}
export function colorFor(i: number): string {
  const p = palette();
  return p[i % p.length];
}

export function seqColorscale(): [number, string][] {
  return SEQ_STOPS.map((c, i) => [i / (SEQ_STOPS.length - 1), c]);
}
export function divColorscale(): [number, string][] {
  return _dark ? DIV_DARK : DIV_LIGHT;
}

// Shared Plotly layout so every chart reads as one system in both modes.
export function baseLayout(overrides: Record<string, any> = {}): Record<string, any> {
  const c = _dark ? INK.dark : INK.light;
  const axis = {
    gridcolor: c.grid,
    zerolinecolor: c.grid,
    linecolor: c.muted,
    tickfont: { color: c.muted, size: 11 },
    titlefont: { color: c.text, size: 12 },
  };
  return {
    paper_bgcolor: "rgba(0,0,0,0)",
    plot_bgcolor: "rgba(0,0,0,0)",
    font: { family: 'system-ui, -apple-system, "Segoe UI", sans-serif', color: c.text, size: 12 },
    margin: { t: 24, r: 16, b: 44, l: 52 },
    xaxis: { ...axis, ...(overrides.xaxis || {}) },
    yaxis: { ...axis, ...(overrides.yaxis || {}) },
    legend: { font: { color: c.text, size: 11 }, ...(overrides.legend || {}) },
    hoverlabel: { font: { family: 'system-ui, sans-serif' } },
    ...stripAxisKeys(overrides),
  };
}

function stripAxisKeys(o: Record<string, any>): Record<string, any> {
  const { xaxis, yaxis, legend, ...rest } = o;
  return rest;
}

export const plotConfig = { displaylogo: false, responsive: true, displayModeBar: true };
