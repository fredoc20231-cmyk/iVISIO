export interface StatusResponse {
  loaded: boolean;
  project: string;
  n_spots: number;
  n_features: number;
  has_pca: boolean;
  has_clusters: boolean;
  has_umap: boolean;
  has_image: boolean;
  obs_columns: string[];
  results: Record<string, boolean>;
  log: string[];
}

export interface TableResponse {
  columns: string[];
  rows: Record<string, unknown>[];
  total: number;
}

export interface SpatialResponse {
  x: number[];
  y: number[];
  barcodes: string[];
  image_shape: [number, number];
  mode: string;
  legend: string;
  values?: number[];
  labels?: string[];
}

export interface UmapResponse {
  x: number[];
  y: number[];
  labels: string[];
  color_by: string;
}

export interface FeatureResponse {
  x: number[];
  y: number[];
  values: number[];
  feature: string;
}
