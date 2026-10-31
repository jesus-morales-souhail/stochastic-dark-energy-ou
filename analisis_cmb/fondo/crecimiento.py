#!/usr/bin/env python3
"""
Exploration only (CRECIMIENTO.md): relative change of f sigma_8(z) when the residual is a constant shift
delta Omega_L = +-2.633e-2 (the theta -> 0 mode that sets the C2 bound) on Planck flat LCDM, H0 fixed, so
Omega_m -> Omega_m -+ delta.

  f(z)       = Omega_m(z)^0.55 (as requested); the exact f = d ln D / d ln a from the ODE is also given as a check.
  sigma_8(z) = sigma_8 D(z), with D from the linear growth equation in ln a,
               D'' + (2 + d ln H / d ln a) D' - 3/2 Omega_m(a) D = 0, D = a deep in matter domination (a = 1e-3),
               flat LCDM without radiation.
  Baseline normalisation: D -> a in matter domination, i.e. the early universe and the primordial amplitude are
  unchanged (the residual is a late-time dark-energy shift). Variant: sigma_8(z=0) held fixed.

Errors: DESI DR1 full-shape, DESI 2024 V (arXiv:2411.12021), Sec. 7.1, Table 9, ShapeFit-only column, f sigma_s8 /
(f sigma_s8)_fid per bin (MAP values; asymmetric errors symmetrised as the mean). arXiv:2411.12022 (DESI 2024 VII)
gives only cosmological parameters, no per-bin f sigma_8.

Writes analisis_cmb/fondo/crecimiento.{json,txt}; nothing in results/ or figures/.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent
OM_P = 0.3153          # Planck 2018 TT,TE,EE+lowE+lensing; Omega_L = 0.6847 as in the S_z table
DELTA = 2.633e-2       # C2 bound, theta -> 0 mode = constant delta Omega_L
Z_OUT = [0.3, 0.5, 0.7, 1.0]
# DESI 2024 V, Table 9 (ShapeFit only): tracer, z_eff, ratio, +err, -err
DESI_FS = [("BGS", 0.295, 0.80, 0.20, 0.20), ("LRG1", 0.510, 1.09, 0.12, 0.14), ("LRG2", 0.706, 1.05, 0.12, 0.12),
           ("LRG3", 0.919, 0.96, 0.11, 0.10), ("ELG2", 1.317, 0.95, 0.11, 0.08), ("QSO", 1.491, 1.16, 0.12, 0.12)]
NEAREST = {0.3: "BGS", 0.5: "LRG1", 0.7: "LRG2", 1.0: "LRG3"}


def growth(om: float):
    """Return functions D(z) (D = a at a = 1e-3) and f_exact(z)."""
    E2 = lambda a: om * a**-3 + 1 - om
    def rhs(lna, y):
        a = np.exp(lna)
        om_a = om * a**-3 / E2(a)
        dlnh = -1.5 * om_a
        return [y[1], -(2 + dlnh) * y[1] + 1.5 * om_a * y[0]]
    a0 = 1e-3
    sol = solve_ivp(rhs, [np.log(a0), 0.0], [a0, a0], rtol=1e-11, atol=1e-14, dense_output=True)
    D = lambda z: sol.sol(-np.log1p(z))[0]
    f = lambda z: sol.sol(-np.log1p(z))[1] / sol.sol(-np.log1p(z))[0]
    return D, f


def om_z(om: float, z: float) -> float:
    return om * (1 + z) ** 3 / (om * (1 + z) ** 3 + 1 - om)


def fs8(om: float, z: float, which_f: str, norm: str, D0_fid: float) -> float:
    D, fex = growth(om)
    f = om_z(om, z) ** 0.55 if which_f == "approx" else fex(z)
    d = D(z) if norm == "early" else D(z) / D(0.0) * D0_fid
    return f * d


def main() -> None:
    Dfid, ffid = growth(OM_P)
    D0 = Dfid(0.0)
    err = {t: 0.5 * (ep + em) for t, _, _, ep, em in DESI_FS}
    rows = []
    for z in Z_OUT:
        r = {"z": z}
        for sgn, lab in ((+1, "dOL_plus"), (-1, "dOL_minus")):
            om = OM_P - sgn * DELTA
            for wf in ("approx", "exact"):
                for norm in ("early", "sigma8_0_fixed"):
                    r[f"{lab}|{wf}|{norm}"] = fs8(om, z, wf, norm, D0) / fs8(OM_P, z, wf, norm, D0) - 1
        r["parts_dOL_plus"] = {"f_approx": (om_z(OM_P - DELTA, z) / om_z(OM_P, z)) ** 0.55 - 1,
                               "D_early": growth(OM_P - DELTA)[0](z) / Dfid(z) - 1}
        t = NEAREST[z]
        r["desi_bin"], r["desi_frac_err"] = t, err[t]
        base = max(abs(r["dOL_plus|approx|early"]), abs(r["dOL_minus|approx|early"]))
        r["error_over_shift"] = err[t] / base
        rows.append(r)
    out = {"Omega_m_Planck": OM_P, "delta_Omega_L": DELTA, "f_exact_vs_approx_fid":
           {str(z): ffid(z) / om_z(OM_P, z) ** 0.55 - 1 for z in Z_OUT},
           "desi_dr1_fs_table9_shapefit": [dict(zip(("tracer", "z_eff", "ratio", "err_plus", "err_minus"), d)) for d in DESI_FS],
           "rows": rows}
    (HERE / "crecimiento.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    lines = [f"Planck Omega_m = {OM_P}, delta Omega_L = +-{DELTA}", "",
             "z    | dfs8/fs8 (+dOL; -dOL), f=Om^0.55, D early | f exacta | sigma8(0) fijo | partes (+dOL): f, D | DESI tramo, error | error/cambio"]
    for r in rows:
        lines.append(f"{r['z']:.1f}  | {100*r['dOL_plus|approx|early']:+.2f} %; {100*r['dOL_minus|approx|early']:+.2f} % | "
                     f"{100*r['dOL_plus|exact|early']:+.2f} % | {100*r['dOL_plus|approx|sigma8_0_fixed']:+.2f} % | "
                     f"{100*r['parts_dOL_plus']['f_approx']:+.2f} %, {100*r['parts_dOL_plus']['D_early']:+.2f} % | "
                     f"{r['desi_bin']} {100*r['desi_frac_err']:.1f} % | {r['error_over_shift']:.1f}")
    lines.append("f exacta / Om(z)^0.55 - 1 (fiducial): " + ", ".join(f"z={k}: {v:+.4f}" for k, v in out["f_exact_vs_approx_fid"].items()))
    (HERE / "crecimiento.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
