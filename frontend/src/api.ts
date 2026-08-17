import axios from "axios";
import type {
  FeatureResponse,
  SpatialResponse,
  StatusResponse,
  TableResponse,
  UmapResponse,
} from "./types";

const http = axios.create({ baseURL: "/api" });

// Extract a clean, human-readable message from a FastAPI error response.
export function apiError(err: unknown): string {
  if (axios.isAxiosError(err)) {
    const detail = err.response?.data?.detail;
    if (typeof detail === "string") return detail;
    if (detail) return JSON.stringify(detail);
    return err.message;
  }
  return String(err);
}

export const api = {
  createSession: async (): Promise<string> => {
    const { data } = await http.post("/session");
    return data.session_id;
  },

  status: async (sid: string): Promise<StatusResponse> =>
    (await http.get(`/${sid}/status`)).data,

  referenceCatalog: async () => (await http.get("/reference-catalog")).data.catalog,
  aiDictionary: async () => (await http.get("/ai-dictionary")).data.dictionary,

  load: async (sid: string, form: FormData) =>
    (await http.post(`/${sid}/load`, form)).data,

  features: async (sid: string, q = ""): Promise<string[]> =>
    (await http.get(`/${sid}/features`, { params: { q, limit: 50 } })).data.features,

  qc: async (sid: string, body: object) => (await http.post(`/${sid}/qc`, body)).data,
  qcFilter: async (sid: string, body: object) =>
    (await http.post(`/${sid}/qc/filter`, body)).data,

  normalize: async (sid: string, body: object) =>
    (await http.post(`/${sid}/normalize`, body)).data,
  cluster: async (sid: string, body: object) =>
    (await http.post(`/${sid}/cluster`, body)).data,

  umap: async (sid: string, color = "clusters"): Promise<UmapResponse> =>
    (await http.get(`/${sid}/umap`, { params: { color } })).data,
  spatial: async (
    sid: string,
    params: { mode: string; feature?: string; group?: string }
  ): Promise<SpatialResponse> =>
    (await http.get(`/${sid}/spatial`, { params })).data,
  feature: async (sid: string, feature: string, reduction = "umap"): Promise<FeatureResponse> =>
    (await http.get(`/${sid}/feature`, { params: { feature, reduction } })).data,

  signature: async (sid: string, body: object) =>
    (await http.post(`/${sid}/signature`, body)).data,

  markers: async (sid: string, body: object): Promise<TableResponse> =>
    (await http.post(`/${sid}/markers`, body)).data,
  de: async (sid: string, body: object): Promise<TableResponse> =>
    (await http.post(`/${sid}/de`, body)).data,
  groups: async (sid: string, group_col = "clusters"): Promise<string[]> =>
    (await http.get(`/${sid}/groups`, { params: { group_col } })).data.groups,

  svg: async (sid: string, body: object): Promise<TableResponse> =>
    (await http.post(`/${sid}/svg`, body)).data,

  loadReference: async (sid: string, form: FormData) =>
    (await http.post(`/${sid}/reference/load`, form)).data,
  deconv: async (sid: string) => (await http.post(`/${sid}/reference/deconv`)).data,

  aiAnnotate: async (sid: string, body: object): Promise<TableResponse> =>
    (await http.post(`/${sid}/ai-annotate`, body)).data,
  enrichment: async (sid: string, body: object): Promise<TableResponse> =>
    (await http.post(`/${sid}/enrichment`, body)).data,
  ligrec: async (sid: string, body: object): Promise<TableResponse> =>
    (await http.post(`/${sid}/ligrec`, body)).data,

  imageUrl: (sid: string) => `/api/${sid}/image`,
  downloadUrl: (sid: string, artifact: string, prefix: string) =>
    `/api/${sid}/download/${artifact}?prefix=${encodeURIComponent(prefix)}`,
};
