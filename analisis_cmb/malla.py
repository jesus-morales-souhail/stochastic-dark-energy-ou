#!/usr/bin/env python3
"""
sigma_X grid summary: resultados/takahashi/malla/sigma_*/resumen.json -> malla.json, malla.md.

Per point: p50/p95 of Delta chi^2 (TT, phiphi) and p16/p50/p84 of H0.
Scaling fits (least squares in log10-log10 over the grid):
    p95(Delta chi^2) ∝ sigma_X^n,   sigma(H0) ≡ (p84 - p16)/2 ∝ sigma_X^m
Residual = rms of log10 residuals. Points with p95 <= 0 cannot enter a power law: counted, excluded.
Bound: largest sigma_X with p95 < 4 (checked over the whole grid).
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
GRID = HERE / "resultados" / "takahashi" / "malla"
TAGS = {"0.001": "theta0.001_Delta0.05", "1": "theta1_Delta0.05"}


def fit_power(x, y):
    x, y = np.asarray(x), np.asarray(y)
    ok = y > 0
    if ok.sum() < 2:
        return {"exponent": None, "n_points": int(ok.sum()), "n_excluded_nonpositive": int((~ok).sum())}
    lx, ly = np.log10(x[ok]), np.log10(y[ok])
    slope, icpt = np.polyfit(lx, ly, 1)
    res = ly - (slope * lx + icpt)
    return {"exponent": float(slope), "log10_prefactor": float(icpt),
            "rms_log10_residual": float(np.sqrt(np.mean(res**2))),
            "max_abs_log10_residual": float(np.max(np.abs(res))),
            "n_points": int(ok.sum()), "n_excluded_nonpositive": int((~ok).sum())}


def main():
    rows = {th: [] for th in TAGS}
    for d in sorted(GRID.glob("sigma_*"), key=lambda p: float(p.name[6:])):
        s = json.loads((d / "resumen.json").read_text())
        sig = s["meta"]["sigma_X"]
        for th, tag in TAGS.items():
            c = s["cases"][tag]
            h = c["H0"]["p2.5_16_50_84_97.5"]
            rows[th].append({
                "sigma_X": sig, "dir": d.name,
                "TT_p50": c["obs"]["TT"]["dchi2_p50"], "TT_p95": c["obs"]["TT"]["dchi2_p95"],
                "pp_p50": c["obs"]["pp"]["dchi2_p50"], "pp_p95": c["obs"]["pp"]["dchi2_p95"],
                "H0_p16": h[1], "H0_p50": h[2], "H0_p84": h[3], "sigma_H0": 0.5 * (h[3] - h[1]),
                "fails": {o: c["obs"][o]["n_fail_camb"] + c["obs"][o]["n_fail_silent_nan"] for o in ("TT", "pp")}
                         | {"H0": c["H0"]["n_fail"], "rho_le_0": c["n_rho_le_0"]},
            })
    out = {"grid_dir": str(GRID.relative_to(HERE)), "rows": rows, "fits": {}, "bound": {}}
    for th, r in rows.items():
        x = [q["sigma_X"] for q in r]
        out["fits"][th] = {
            "TT_p95": fit_power(x, [q["TT_p95"] for q in r]),
            "pp_p95": fit_power(x, [q["pp_p95"] for q in r]),
            "sigma_H0": fit_power(x, [q["sigma_H0"] for q in r]),
        }
        for o in ("TT", "pp"):
            p95 = np.array([q[f"{o}_p95"] for q in r])
            out["bound"][f"{th}_{o}"] = {"all_p95_lt_4": bool(np.all(p95 < 4)), "max_p95": float(p95.max()),
                                         "sigma_X_at_max": float(x[int(np.argmax(p95))])}
    (GRID / "malla.json").write_text(json.dumps(out, indent=2), encoding="utf-8")

    L = ["# Malla σ_X — takahashi, Δ = 0.05, 200 realizaciones, semilla 20261003", "",
         f"Fuente: `{out['grid_dir']}/sigma_*/resumen.json`. σ(H0) ≡ (p84 − p16)/2.", ""]
    for th, r in rows.items():
        L += [f"## θ = {th}", "",
              "| σ_X | TT p50 | TT p95 | φφ p50 | φφ p95 | H0 p16 | H0 p50 | H0 p84 | σ(H0) | fallos |",
              "|---|---|---|---|---|---|---|---|---|---|"]
        for q in r:
            nf = sum(q["fails"].values())
            L.append(f"| {q['sigma_X']:.4e} | {q['TT_p50']:.3g} | {q['TT_p95']:.3g} | {q['pp_p50']:.3g} | "
                     f"{q['pp_p95']:.3g} | {q['H0_p16']:.3f} | {q['H0_p50']:.3f} | {q['H0_p84']:.3f} | "
                     f"{q['sigma_H0']:.4f} | {nf} |")
        L.append("")
        for k, f in out["fits"][th].items():
            if f["exponent"] is None:
                L.append(f"- {k}: sin ajuste ({f['n_points']} puntos > 0)")
            else:
                L.append(f"- {k} ∝ σ_X^{f['exponent']:.3f} (rms residuo log10 {f['rms_log10_residual']:.3g}, "
                         f"máx {f['max_abs_log10_residual']:.3g}; {f['n_points']} puntos, "
                         f"{f['n_excluded_nonpositive']} excluidos por p95 ≤ 0)")
        L.append("")
    L.append("## Cota")
    for k, b in out["bound"].items():
        L.append(f"- {k}: p95 < 4 en toda la malla: {b['all_p95_lt_4']}; máx p95 = {b['max_p95']:.3g} "
                 f"(σ_X = {b['sigma_X_at_max']:.4e})")
    (GRID / "malla.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print("\n".join(L))


if __name__ == "__main__":
    main()
