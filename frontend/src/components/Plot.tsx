import Plotly from "plotly.js-dist-min";
import createPlotlyComponent from "react-plotly.js/factory";

// Bind react-plotly.js to the lightweight dist-min bundle so we don't pull in
// the full plotly.js build. Used across every visualization tab.
const Plot = createPlotlyComponent(Plotly);

export default Plot;

// Re-export the shared theme helpers so chart components have one import site.
export { colorFor, palette, seqColorscale, divColorscale, baseLayout, plotConfig } from "../theme";
