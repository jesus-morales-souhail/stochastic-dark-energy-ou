#!/usr/bin/env python3
"""
Exploration only (RELOJ.md): the 95 % profile bound on sigma_X with four "clocks" of stationarity for the OU
kernel. Same data, same official covariance, same sigma_X grid, same theta search and the same criterion
(Delta lnL = -1.92 from the profile maximum, root by brentq) as scripts/profile_sigma_x_desi.py, whose functions
are reused: only the coordinate x in add_ou_kernel (exp(-theta |x_i - x_j|)) changes.

  (a) x = ln(1+z) = -ln a        (current; must reproduce results/profile_sigma_x/profile_sigma_x.json)
  (b) cosmic time t(z), units 1/H0
  (c) conformal time eta(z), units 1/H0
  (d) x = z
Clocks (b), (c) from LCDM-Planck: H0 = 67.36, Omega_m = 0.3153, Omega_r = omega_r / h^2 with omega_r of
analisis_cmb/everpresent/resultados/rd_planck.json (photons + 3.046 massless neutrinos), flat.
theta is searched in the same fixed range (1e-3, 20) for every clock, as in the original script.

Writes nothing outside analisis_cmb/reloj/: perfil_reloj.{json,txt}.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import profile_sigma_x_desi as PS  # noqa: E402  (module import loads the data; main() is not called)

H0, OM = 67.36, 0.3153
OR = json.loads((ROOT / "analisis_cmb" / "everpresent" / "resultados" / "rd_planck.json").read_text())["omega_r"] / (H0 / 100) ** 2
OL = 1.0 - OM - OR
THR = -1.92


def E(a: float) -> float:
    return np.sqrt(OR / a**4 + OM / a**3 + OL)


def cosmic_time(z: float) -> float:
    return quad(lambda a: 1.0 / (a * E(a)), 0.0, 1.0 / (1.0 + z), epsabs=0, epsrel=1e-12, limit=200)[0]


def conformal_time(z: float) -> float:
    return quad(lambda a: 1.0 / (a * a * E(a)), 0.0, 1.0 / (1.0 + z), epsabs=0, epsrel=1e-12, limit=200)[0]


def clocks(z: np.ndarray) -> dict[str, np.ndarray]:
    return {"a_ln_a": np.log1p(z),
            "b_cosmic_time": np.array([cosmic_time(v) for v in z]),
            "c_conformal_time": np.array([conformal_time(v) for v in z]),
            "d_redshift": np.asarray(z, dtype=float).copy()}


def profile(x: np.ndarray) -> dict:
    """Same steps as profile_sigma_x_desi.main(), without writing anything."""
    PS.x = x
    ll_lcdm = PS.logL_gaussian(PS.residuals, PS.C_MEAS)
    sig_grid = np.geomspace(1e-6, 5e-2, 48)
    prof = [PS.max_logL_over_theta(float(s)) for s in sig_grid]
    dlogL = np.array([p[0] for p in prof]) - ll_lcdm
    dlogL_rel = dlogL - np.max(dlogL)
    above = np.where(dlogL_rel >= THR)[0]
    s_95_grid = float(sig_grid[above[-1]]) if len(above) else float(sig_grid[0])
    ll_max = float(np.max(dlogL)) + ll_lcdm
    j = above[-1] if len(above) else 0
    s_95 = s_95_grid
    if j + 1 < len(sig_grid):
        g = lambda s: PS.max_logL_over_theta(float(s))[0] - ll_max - THR
        s_95 = float(brentq(g, sig_grid[j], sig_grid[j + 1], xtol=1e-10 * sig_grid[j]))
    th_95 = PS.max_logL_over_theta(s_95)[1]
    return {"sigma_X_95": s_95, "sigma_X_95_grid": s_95_grid, "theta_at_bound": th_95,
            "max_dlogL_vs_lcdm": float(np.max(dlogL)), "sigma_X_at_max": float(sig_grid[int(np.argmax(dlogL))]),
            "theta_best_profile": [float(p[1]) for p in prof], "dlogL_profile": [float(v) for v in dlogL],
            "theta_dx_max_at_bound": float(th_95 * (x.max() - x.min()))}


def main() -> None:
    z = np.asarray(PS.z_eff, dtype=float)
    ref = json.loads((ROOT / "results" / "profile_sigma_x" / "profile_sigma_x.json").read_text())
    X = clocks(z)
    out = {"z_eff": z.tolist(), "cosmology": {"H0": H0, "Omega_m": OM, "Omega_r": OR, "Omega_L": OL},
           "theta_search": "(1e-3, 20), same for every clock", "reference_sigma_X_95": ref["sigma_X_95CL_upper_profile"],
           "clocks": {}}
    for name, x in X.items():
        r = profile(x)
        d = np.abs(x[:, None] - x[None, :])
        r.update(x=x.tolist(), dx_adjacent=np.abs(np.diff(x)).tolist(), dx_matrix=d.tolist(),  # z_eff is sorted
                 span=float(x.max() - x.min()), rel_to_reference=r["sigma_X_95"] / ref["sigma_X_95CL_upper_profile"] - 1)
        out["clocks"][name] = r
        print(f"{name:18s} sigma_X < {r['sigma_X_95']:.6e} ({r['rel_to_reference']:+.3%})  theta = {r['theta_at_bound']:.4g}"
              f"  theta*span = {r['theta_dx_max_at_bound']:.3g}  max dlnL = {r['max_dlogL_vs_lcdm']:+.4f}", flush=True)
    # theta -> 0: the kernel tends to the rank-1 S_i S_j sigma_X^2 for any clock; bound with theta fixed at 1e-12
    def g0(s):
        return PS.logL(1e-12, float(s)) - PS.logL_gaussian(PS.residuals, PS.C_MEAS) - THR
    sg = np.geomspace(1e-6, 5e-2, 48)
    k = np.nonzero([g0(v) >= 0 for v in sg])[0][-1]
    out["fully_correlated_limit_theta_1e-12"] = float(brentq(g0, sg[k], sg[k + 1], xtol=1e-10 * sg[k]))
    out["theta_at_lower_edge_near_bound"] = {n_: bool(abs(r_["theta_at_bound"] - 1e-3) < 1e-6) for n_, r_ in out["clocks"].items()}
    out["check_a_reproduces_reference"] = {"ours": out["clocks"]["a_ln_a"]["sigma_X_95"],
                                           "reference": ref["sigma_X_95CL_upper_profile"],
                                           "identical": out["clocks"]["a_ln_a"]["sigma_X_95"] == ref["sigma_X_95CL_upper_profile"]}
    (HERE / "perfil_reloj.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    lines = [f"z_eff: {', '.join(f'{v:.3f}' for v in z)}", "",
             f"Limite totalmente correlacionado (theta = 1e-12, igual para todo reloj): sigma_X < {out['fully_correlated_limit_theta_1e-12']:.5e}",
             f"theta en la cota en el borde inferior de la busqueda (1e-3): {json.dumps(out['theta_at_lower_edge_near_bound'])}", "",
             f"{'reloj':18s} {'sigma_X 95%':>12s} {'vs 2.633e-2':>11s} {'theta cota':>10s} {'theta*span':>10s} {'max dlnL':>9s}"]
    for name, r in out["clocks"].items():
        lines.append(f"{name:18s} {r['sigma_X_95']:12.5e} {r['rel_to_reference']:+11.3%} {r['theta_at_bound']:10.4g}"
                     f" {r['theta_dx_max_at_bound']:10.3g} {r['max_dlogL_vs_lcdm']:+9.2e}")
    lines += ["", "Separaciones entre bins consecutivos en orden de z (dx) y extension total:"]
    for name, r in out["clocks"].items():
        lines.append(f"{name:18s} " + " ".join(f"{v:.4f}" for v in r["dx_adjacent"]) + f"   total {r['span']:.4f}")
    lines += ["", f"(a) reproduce la referencia: {json.dumps(out['check_a_reproduces_reference'])}"]
    (HERE / "perfil_reloj.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
