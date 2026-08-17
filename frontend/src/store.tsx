import { createContext, useContext, useEffect, useState, type ReactNode } from "react";
import { api, apiError } from "./api";
import { setDarkMode } from "./theme";
import type { StatusResponse } from "./types";

interface Store {
  sid: string | null;
  status: StatusResponse | null;
  busy: boolean;
  toast: string | null;
  dark: boolean;
  toggleDark: () => void;
  refresh: () => Promise<void>;
  setBusy: (b: boolean) => void;
  notify: (msg: string) => void;
  run: <T>(fn: () => Promise<T>, okMsg?: string) => Promise<T | undefined>;
}

const Ctx = createContext<Store | null>(null);

export function StoreProvider({ children }: { children: ReactNode }) {
  const [sid, setSid] = useState<string | null>(null);
  const [status, setStatus] = useState<StatusResponse | null>(null);
  const [busy, setBusy] = useState(false);
  const [toast, setToast] = useState<string | null>(null);
  const [dark, setDark] = useState<boolean>(
    () => window.matchMedia?.("(prefers-color-scheme: dark)").matches ?? false
  );

  useEffect(() => {
    api.createSession().then(setSid).catch((e) => setToast(apiError(e)));
  }, []);

  // Keep the theme module + document attribute in sync so Plotly charts and
  // CSS both respond to the toggle. Dark mode is explicitly selected, not an
  // automatic invert.
  useEffect(() => {
    setDarkMode(dark);
    document.body.dataset.theme = dark ? "dark" : "light";
  }, [dark]);

  const toggleDark = () => setDark((d) => !d);

  const refresh = async () => {
    if (!sid) return;
    try {
      setStatus(await api.status(sid));
    } catch (e) {
      setToast(apiError(e));
    }
  };

  useEffect(() => {
    if (sid) refresh();
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [sid]);

  const notify = (msg: string) => {
    setToast(msg);
    window.clearTimeout((notify as any)._t);
    (notify as any)._t = window.setTimeout(() => setToast(null), 4000);
  };

  // Wraps an async API call with busy state, error surfacing, and a status
  // refresh so every tab shares consistent feedback.
  const run = async <T,>(fn: () => Promise<T>, okMsg?: string): Promise<T | undefined> => {
    setBusy(true);
    try {
      const result = await fn();
      await refresh();
      if (okMsg) notify(okMsg);
      return result;
    } catch (e) {
      notify(apiError(e));
      return undefined;
    } finally {
      setBusy(false);
    }
  };

  return (
    <Ctx.Provider value={{ sid, status, busy, toast, dark, toggleDark, refresh, setBusy, notify, run }}>
      {children}
    </Ctx.Provider>
  );
}

export function useStore(): Store {
  const ctx = useContext(Ctx);
  if (!ctx) throw new Error("useStore must be used within StoreProvider");
  return ctx;
}
