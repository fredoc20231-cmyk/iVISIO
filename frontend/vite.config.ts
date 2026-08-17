import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";

// During development the API base is proxied to the FastAPI server so the
// browser talks to a single origin (no CORS complications). In production the
// built static assets are served behind the same reverse proxy as /api.
export default defineConfig({
  plugins: [react()],
  server: {
    port: 5173,
    proxy: {
      "/api": {
        target: process.env.VITE_API_TARGET || "http://localhost:8000",
        changeOrigin: true,
      },
    },
  },
});
