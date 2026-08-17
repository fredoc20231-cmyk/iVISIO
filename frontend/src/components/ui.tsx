import type { ButtonHTMLAttributes, ReactNode } from "react";
import type { TableResponse } from "../types";

export function Card({ title, children, actions }: { title: string; children: ReactNode; actions?: ReactNode }) {
  return (
    <div className="card">
      <div className="card-header">
        <span>{title}</span>
        {actions}
      </div>
      <div className="card-body">{children}</div>
    </div>
  );
}

export function Field({ label, children }: { label: string; children: ReactNode }) {
  return (
    <label className="field">
      <span>{label}</span>
      {children}
    </label>
  );
}

export function Button(props: ButtonHTMLAttributes<HTMLButtonElement> & { variant?: string }) {
  const { variant = "primary", className = "", ...rest } = props;
  return <button className={`btn btn-${variant} ${className}`} {...rest} />;
}

export function DataTable({ data, maxRows = 25 }: { data: TableResponse | null; maxRows?: number }) {
  if (!data) return <p className="muted">No results yet.</p>;
  if (!data.rows.length) return <p className="muted">No rows returned.</p>;
  const rows = data.rows.slice(0, maxRows);
  return (
    <div className="table-wrap">
      <table>
        <thead>
          <tr>{data.columns.map((c) => <th key={c}>{c}</th>)}</tr>
        </thead>
        <tbody>
          {rows.map((r, i) => (
            <tr key={i}>
              {data.columns.map((c) => (
                <td key={c}>{formatCell(r[c])}</td>
              ))}
            </tr>
          ))}
        </tbody>
      </table>
      <p className="muted">
        Showing {rows.length} of {data.total} rows.
      </p>
    </div>
  );
}

function formatCell(v: unknown): string {
  if (v === null || v === undefined) return "";
  if (typeof v === "number") return Number.isInteger(v) ? String(v) : v.toPrecision(4);
  return String(v);
}
