#!/usr/bin/env python3
"""
Exploration only (EXACTO.md): exact (non-linear) bound on the OU tremor in convention B. Nothing in results/ or papers changes.

Convention B: rho_DE(z)/rho_crit,0 = Omega_L + X(z) - X(0) (the constant part of the residual is Lambda; omega_m is the
matter density and does not respond to X). H^2 = 100^2 [omega_m (1+z)^3 + omega_r (1+z)^4 + (h^2 - omega_m - omega_r)
+ h^2 (X(z) - X(0))], as escalon.py (B) with the step replaced by the OU path. X: stationary OU in u = ln(1+z) on the grid
of cmb_fondo.py (uniform, du = 2e-3, 0 <= z <= z_*) with the DESI z inserted as exact nodes (AR(1) with the step-dependent
decay). Distances by trapezoid on the grid (exact nodes at the data z), the 13 DESI numbers and the calibrated theta_*
(r_s(omega_m), z_*(omega_m), CAL of escalon.py). For sigma_X in {0.05, 0.1, 0.2, 0.3, 0.5, 0.7} x theta in {0.1, 0.3, 1, 3, 10}:
N_REAL = 2000 realisations (seed per cell); a realisation is discarded if H^2 <= 0 somewhere at the LCDM best fit of the 15
rows (h, omega_m fixed there); during the fit, (h, omega_m) with H^2 <= 0 are penalised. For every surviving realisation the
background (h, omega_m) is fitted to BAO + theta_* + omega_m (15 rows) and to BAO only (13 rows); chi^2_min stored.
Bound (i): the sigma_X at which P(chi^2_min <= chi^2_obs) = 5 %, chi^2_obs = 17.24 (15 rows) and 10.28 (BAO only),
interpolated linearly in (ln sigma_X, P). Bound (ii): the linear convention-B bounds of cmb_fondo.json (marginalised and
profiled, -1.92) in the same cells.

    ../venv/bin/python exacto.py     # writes exacto.{json,txt} and exacto_chi2.npz
"""

from __future__ import annotations

import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "escalon"))
import escalon as E  # noqa: E402

C_KMS = E.C_KMS
OBS13, COV13, Z13, KINDS13 = E.OBS, E.COV, E.Z, E.CU.KINDS
TH_P, SIG_TH, WM_P, SIG_WM, OMR, ZSTAR = E.TH_P, E.SIG_TH, E.WM_P, E.SIG_WM, E.OMR, E.ZSTAR
SIGMAS = (0.05, 0.1, 0.2, 0.3, 0.5, 0.7)
THETAS = (0.1, 0.3, 1.0, 3.0, 10.0)
N_REAL = 2000
SEED = 20261112
DU = 2e-3
CHI2_OBS = {"BAO+CMB": None, "BAO": None}  # filled in main from the X = 0 fits (must be 17.24 and 10.28)
BOUNDS_P = ([0.4, 0.05], [1.0, 0.30])

# grid: uniform in u plus the exact data nodes
_u = np.linspace(0.0, np.log1p(ZSTAR), int(np.ceil(np.log1p(ZSTAR) / DU)) + 1)
U = np.unique(np.concatenate([_u, np.log1p(np.unique(Z13))]))
ZG = np.expm1(U)
N = len(U)
DUK = np.diff(U)
IDX = [int(np.argmin(np.abs(U - np.log1p(z)))) for z in Z13]
assert max(abs(U[i] - np.log1p(z)) for i, z in zip(IDX, Z13)) < 1e-12
_1PZ3, _1PZ4, _1PZ = (1 + ZG) ** 3, (1 + ZG) ** 4, 1 + ZG
IS_DM, IS_DH = np.array([k == "DM" for k in KINDS13]), np.array([k == "DH" for k in KINDS13])
Y15 = np.concatenate([OBS13, [TH_P, WM_P]])
L13 = np.linalg.cholesky(COV13)


def ou_paths(sigma: float, theta: float, n: int, rng) -> np.ndarray:
    """n stationary OU paths on U (AR(1) with step-dependent decay), shape (n, N)."""
    X = np.empty((n, N))
    X[:, 0] = sigma * rng.standard_normal(n)
    rho = np.exp(-theta * DUK)
    s = sigma * np.sqrt(1.0 - rho**2)
    xi = rng.standard_normal((n, N - 1))
    for k in range(N - 1):
        X[:, k + 1] = rho[k] * X[:, k] + s[k] * xi[:, k]
    return X


def model(h: float, wm: float, dX: np.ndarray):
    """dX = X(z) - X(0) on the grid. Returns (15-vector, min E^2 h^2) or (None, min) if H^2 <= 0."""
    h2 = h * h
    e2 = wm * _1PZ3 + OMR * _1PZ4 + (h2 - wm - OMR) + h2 * dX  # H^2 / 100^2
    m = float(e2.min())
    if m <= 0.0:
        return None, m
    inv = 1.0 / np.sqrt(e2)  # 100 / H
    g = _1PZ * inv  # dz / (H/100) per du (dz = e^u du)
    dm_u = np.concatenate([[0.0], np.cumsum(0.5 * (g[:-1] + g[1:]) * DUK)])  # cumulative trapezoid, exact at the nodes
    dm = np.array([dm_u[i] for i in IDX])
    dh = np.array([inv[i] for i in IDX])
    d = np.where(IS_DM, dm, np.where(IS_DH, dh, np.cbrt(Z13 * dm * dm * dh)))
    bao = C_KMS / (100.0 * float(E.RD_SPL(wm))) * d
    dm_star = C_KMS / 100.0 * (dm_u[-1] + inv[-1] * (E.zstar_of(wm) - ZSTAR))
    th = 100.0 * float(E.RS_SPL(wm)) / dm_star / E.CAL
    return np.concatenate([bao, [th, wm]]), m


def resid(p, dX, rows):
    v, m = model(p[0], p[1], dX)
    if v is None:
        return np.full(len(rows), 1e3 * (1.0 - m))
    r = Y15 - v
    out = np.linalg.solve(L13, r[:13])
    if len(rows) == 15:
        out = np.concatenate([out, [r[13] / SIG_TH, r[14] / SIG_WM]])
    return out


def fit(dX: np.ndarray, rows: int, x0=(0.682, 0.1407)) -> tuple[float, float, float]:
    r = least_squares(resid, np.asarray(x0), bounds=BOUNDS_P, args=(dX, range(rows)), method="trf",
                      x_scale=[0.01, 0.001], ftol=1e-10, xtol=1e-10, max_nfev=300)
    return float(2.0 * r.cost), float(r.x[0]), float(r.x[1])


P0 = None  # LCDM best fit (h, omega_m) of the 15 rows, set in main


def _cell(args):
    sigma, theta, i_cell = args
    rng = np.random.default_rng([SEED, i_cell])
    X = ou_paths(sigma, theta, N_REAL, rng)
    dX = X - X[:, [0]]
    h0, wm0 = P0
    e2min = (wm0 * _1PZ3 + OMR * _1PZ4 + (h0**2 - wm0 - OMR) + h0**2 * dX).min(axis=1)
    ok = e2min > 0.0
    chi_cmb, chi_bao = np.full(N_REAL, np.nan), np.full(N_REAL, np.nan)
    for i in np.nonzero(ok)[0]:
        chi_cmb[i] = fit(dX[i], 15)[0]
        chi_bao[i] = fit(dX[i], 13)[0]
    return sigma, theta, ok, chi_cmb, chi_bao


def crossing(sig: np.ndarray, p: np.ndarray, target: float = 0.05):
    """sigma at which P crosses target (linear in ln sigma); None if never below target; '<min' if below at the first."""
    if p[0] <= target:
        return "<" + f"{sig[0]:g}"
    for k in range(len(sig) - 1):
        if p[k] > target >= p[k + 1]:
            t = (p[k] - target) / (p[k] - p[k + 1])
            return float(np.exp(np.log(sig[k]) + t * (np.log(sig[k + 1]) - np.log(sig[k]))))
    return None


def main() -> None:
    global P0
    t0 = time.time()
    zero = np.zeros(N)
    c15, h0, wm0 = fit(zero, 15)
    c13, hb, wmb = fit(zero, 13)
    P0 = (h0, wm0)
    CHI2_OBS["BAO+CMB"], CHI2_OBS["BAO"] = c15, c13
    checks = {"lcdm_15": {"chi2": c15, "H0": 100 * h0, "omega_m": wm0, "escalon_B": {"chi2": 17.238, "H0": 68.19, "omega_m": 0.14072}},
              "lcdm_13": {"chi2": c13, "H0": 100 * hb, "omega_m": wmb, "cmb_fondo": {"chi2": 10.28, "H0": 68.67, "omega_m": 0.14014}},
              "grid": {"nodes": int(N), "du": DU}}
    v15, _ = model(h0, wm0, zero)
    checks["theta_at_lcdm15"] = float(v15[13])
    print(json.dumps(checks, indent=1), flush=True)
    cells = [(s, t, i) for i, (s, t) in enumerate((s, t) for t in THETAS for s in SIGMAS)]
    with Pool() as pool:
        res = pool.map(_cell, cells, chunksize=1)
    print(f"celdas: {time.time() - t0:.0f} s", flush=True)
    lin = json.load(open(HERE / "cmb_fondo.json"))["fits"]
    out = {"sigmas": SIGMAS, "thetas": THETAS, "n_real": N_REAL, "seed": SEED, "chi2_obs": CHI2_OBS, "checks": checks, "cells": {}, "bounds": {}}
    raw = {}
    for sigma, theta, ok, cc, cb in res:
        key = f"theta_{theta:g}|sigma_{sigma:g}"
        raw[key.replace("|", "__")] = np.vstack([ok, cc, cb])
        fin = ok
        out["cells"][key] = {"frac_H2_le_0": float(1 - ok.mean()), "n_ok": int(ok.sum()),
                             "P_cmb_le_obs": float(np.mean(cc[fin] <= c15)) if ok.any() else None,
                             "P_bao_le_obs": float(np.mean(cb[fin] <= c13)) if ok.any() else None,
                             "chi2_cmb_p5_p50_p95": [float(np.percentile(cc[fin], q)) for q in (5, 50, 95)] if ok.any() else None,
                             "chi2_bao_p5_p50_p95": [float(np.percentile(cb[fin], q)) for q in (5, 50, 95)] if ok.any() else None,
                             "P_cmb_le_obs_discarded_as_fail": float(np.mean(np.where(ok, cc, np.inf) <= c15)),
                             "P_bao_le_obs_discarded_as_fail": float(np.mean(np.where(ok, cb, np.inf) <= c13))}
    np.savez_compressed(HERE / "exacto_chi2.npz", **raw, columns=np.array(["ok", "chi2_cmb", "chi2_bao"]))
    sig = np.array(SIGMAS)
    for theta in THETAS:
        k = f"theta_{theta:g}"
        rec = {}
        for lab, pk in (("BAO+CMB", "P_cmb_le_obs"), ("BAO", "P_bao_le_obs")):
            p = np.array([out["cells"][f"{k}|sigma_{s:g}"][pk] if out["cells"][f"{k}|sigma_{s:g}"][pk] is not None else 0.0 for s in SIGMAS])
            pf = np.array([out["cells"][f"{k}|sigma_{s:g}"][pk + "_discarded_as_fail"] for s in SIGMAS])
            rec[lab] = {"P_vs_sigma": p.tolist(), "P_vs_sigma_discarded_as_fail": pf.tolist(),
                        "exact_5pct": crossing(sig, p), "exact_5pct_discarded_as_fail": crossing(sig, pf)}
            lkey = "BAO+CMB|Lambda_absorbs_X0" if lab == "BAO+CMB" else "BAO|Lambda_absorbs_X0"
            rec[lab]["linear_marginalised_192"] = lin[lkey]["marginalised"][f"theta_{float(theta)}"]["192"]
            rec[lab]["linear_marginalised_135"] = lin[lkey]["marginalised"][f"theta_{float(theta)}"]["135"]
            if "profiled" in lin[lkey]:
                rec[lab]["linear_profiled_192"] = lin[lkey]["profiled"][f"theta_{float(theta)}"]["192"]
        rec["frac_H2_le_0_vs_sigma"] = [out["cells"][f"{k}|sigma_{s:g}"]["frac_H2_le_0"] for s in SIGMAS]
        out["bounds"][k] = rec
    out["runtime_s"] = time.time() - t0
    (HERE / "exacto.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    f = lambda v: ("sin cota (<= 0.7)" if v is None else (v if isinstance(v, str) else f"{v:.3g}"))
    g = lambda v: "sin cota" if v is None else f"{v:.3g}"
    L = [f"chi2_obs: BAO+CMB = {c15:.3f} (H0 = {100 * h0:.2f}, omega_m = {wm0:.5f}); BAO = {c13:.3f}; malla {N} nodos; {N_REAL} realizaciones por celda", ""]
    for lab in ("BAO+CMB", "BAO"):
        L.append(f"{lab}: {'theta':>6s} {'exacta (5 %)':>14s} {'exacta, descartes=fallo':>24s} {'lineal marg -1.92':>18s} {'lineal marg -1.35':>18s} {'P(chi2<=obs) por sigma ' + str(SIGMAS):>60s}  frac H2<=0")
        for theta in THETAS:
            r = out["bounds"][f"theta_{theta:g}"]
            L.append(f"{'':9s}{theta:6g} {f(r[lab]['exact_5pct']):>14s} {f(r[lab]['exact_5pct_discarded_as_fail']):>24s} {g(r[lab]['linear_marginalised_192']):>18s} "
                     f"{g(r[lab]['linear_marginalised_135']):>18s} {str([round(x, 3) for x in r[lab]['P_vs_sigma']]):>60s}  {[round(x, 3) for x in r['frac_H2_le_0_vs_sigma']]}")
        L.append("")
    (HERE / "exacto.txt").write_text("\n".join(L) + "\n", encoding="utf-8")
    print("\n".join(L))


if __name__ == "__main__":
    main()
