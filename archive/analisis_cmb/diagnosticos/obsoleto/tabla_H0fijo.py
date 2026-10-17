#!/usr/bin/env python3
"""Markdown table for one sigma_X from resultados/sigma_<s>/resumen.json -> tabla.md (same folder)."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import planck_datos as P

HERE = Path(__file__).resolve().parent
OBS = {"TT": "TT 2≤l<30", "pp": "C_L^φφ (L 8–400)", "thetastar": "θ_*"}
FILE = {"TT": "COM_PowerSpect_CMB-TT-full_R3.01.txt",
        "pp": "…_consext8_bandpowers.dat",
        "thetastar": "Planck18 VI tabla 2"}


def fmt(x, nd=3):
    return "—" if x is None else f"{x:.{nd}g}"


def main(sigma_dir: str):
    d = Path(sigma_dir)
    s = json.loads((d / "resumen.json").read_text())
    refs = s["reference_observables"]
    lines = [f"# σ_X = {s['sigma_X']:.4e} — {s['n_real']} realizaciones, semilla {s['seed']}", "",
             "Δχ² = χ²(curva) − χ²(ΛCDM) frente a datos Planck. Fuente: `resumen.json` y `dchi2_*.txt` de esta carpeta.", "",
             "## Curvas de referencia (Δχ² frente a ΛCDM; θ_* como 100θ_*)", "",
             "| observable | fichero Planck | ΛCDM | CPL repo (diag.) | CPL joint | DESI DR2+CMB (params Planck fijos) |",
             "|---|---|---|---|---|---|"]
    for o, name in OBS.items():
        cells = []
        for k in ("LCDM", "CPL_repo", "CPL_joint", "DESI_DR2_CMB"):
            r = refs[k]
            dc = r[f"dchi2_{o}_vs_LCDM"]
            cells.append(f"{r['thetastar100']:.5f} (Δχ² {dc:.3g})" if o == "thetastar" else f"Δχ² {dc:.3g}")
        lines.append(f"| {name} | {FILE[o]} | " + " | ".join(cells) + " |")

    lines += ["", "## OU por (θ, Δ)", "",
              "| observable | θ | Δ | ρ≤0 | CAMB falla | Δχ² p50 | Δχ² p95 | frac Δχ²>4 | hueco banda68–ΛCDM [σ] | borde banda68 [σ] | p95 vs CPL repo | p95 vs CPL joint | nota |",
              "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for o, name in OBS.items():
        for tag, c in s["cases"].items():
            th = tag.split("_")[0].replace("theta", "")
            de = tag.split("Delta")[1]
            e = c["obs"][o]
            note = "dominado por el suavizado (1/θ = Δ)" if th == "20" else ""
            if th == "0.001":
                note = "θ preferido por DESI a esta σ_X (borde de malla)"
            if "error" in e:
                lines.append(f"| {name} | {th} | {de} | {c['n_bad_rho_le_0']} | {c.get('n_camb_fail', 0)} | — | — | — | — | — | — | — | {e['error']} |")
                continue
            extra = ""
            if o == "thetastar":
                extra = f" 100θ_* media {e['mean']:.5f} [{e['p16']:.5f}, {e['p84']:.5f}]"
            lines.append(
                f"| {name} | {th} | {de} | {c['n_bad_rho_le_0']} | {c.get('n_camb_fail', 0)} | {fmt(e['dchi2_p50'])} | "
                f"{fmt(e['dchi2_p95'])} | {e['frac_dchi2_gt_4']:.3f} | {fmt(e['band68_gap_vs_LCDM_sigma'])} | "
                f"{fmt(e['band68_edge_vs_LCDM_sigma'])} | {fmt(e['dchi2_vs_CPL_repo_p95'])} | "
                f"{fmt(e['dchi2_vs_CPL_joint_p95'])} | {note}{extra} |")

    lines += ["", "## Sensibilidad a Δ: p95(Δ=0.05) − p95(Δ=0.10)", "", "| observable | θ | Δp95 |", "|---|---|---|"]
    for o, name in OBS.items():
        for th in ("0.001", "1", "20"):
            a = s["cases"].get(f"theta{th}_Delta0.05", {}).get("obs", {}).get(o, {})
            b = s["cases"].get(f"theta{th}_Delta0.10", {}).get("obs", {}).get(o, {})
            if "dchi2_p95" in a and "dchi2_p95" in b:
                lines.append(f"| {name} | {th} | {a['dchi2_p95'] - b['dchi2_p95']:+.3g} |")
            else:
                lines.append(f"| {name} | {th} | — (sin datos válidos) |")
    lines.append("")
    lines.append(f"θ_*: σ_Planck = {P.THETASTAR100_ERR} en 100θ_*.")
    (d / "tabla.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else str(HERE / "resultados" / "sigma_2.506e-02"))
