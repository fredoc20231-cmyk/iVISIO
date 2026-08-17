import Plotly from "plotly.js-dist-min";
import createPlotlyComponent from "react-plotly.js/factory";

// Bind react-plotly.js to the lightweight dist-min bundle so we don't pull in
// the full plotly.js build. Used across every visualization tab.
const Plot = createPlotlyComponent(Plotly);

export default Plot;

// A colorblind-friendly categorical palette matching the original R app's
// ivisio_discrete_palette(), so clusters look consistent across ports.
export const IVISIO_PALETTE = [
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFD92F",
  "#A65628", "#F781BF", "#66C2A5", "#8DA0CB", "#E78AC3", "#A6D854",
  "#1B9E77", "#D95F02", "#7570B3", "#66A61E", "#E6AB02", "#A6761D",
];

export function colorFor(index: number): string {
  return IVISIO_PALETTE[index % IVISIO_PALETTE.length];
}
