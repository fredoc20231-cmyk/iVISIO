// react-plotly.js ships its factory without bundled types for the
// dist-min build; declare the minimal surface we use.
declare module "react-plotly.js/factory" {
  import type { ComponentType } from "react";
  const createPlotlyComponent: (plotly: unknown) => ComponentType<any>;
  export default createPlotlyComponent;
}

declare module "plotly.js-dist-min";
