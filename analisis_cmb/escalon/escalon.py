#!/usr/bin/env python3
"""
Exploration only (ESCALON.md): how late can dark energy switch on? Sweep of a deterministic step.

Model (flat):  rho_DE(z) / rho_crit,0 = Omega_L [1 - f s(z)],   s(z) = 1/2 [1 + tanh((ln(1+z) - ln(1+z_t)) / w)],
s -> 1 in the past: f = 1 is dark energy absent before z_t, f = 2 a sign change, f = 0 LCDM. Flatness fixes Omega_L today,
Omega_L [1 - f s(0)] = 1 - Omega_m (- Omega_r in B), so the shape actually used is g(z) = [1 - f s(z)] / [1 - f s(0)]
(s(0) < 4e-5 for w = 0.05; up to 0.15 for w = 0.3, z_t = 0.3). Distances are exact integrals (Simpson on fixed grids, checked
against quad), H^2 > 0 is checked on every grid point (a violation gives chi^2 = inf and is reported).
Grid: z_t in {0.3, 0.5, 0.7, 1.0, 1.3, 1.6, 2.0, 2.4, 3.0} x w in {0.05, 0.3}; profile of chi^2 in f in [-2, 3] (step 0.05)
with the 13 DESI DR2 numbers and the official covariance, in two treatments:
  (A) BAO only: Omega_m and h r_d free (h r_d analytic, Omega_m by bounded Brent), no radiation, as cuesta.py.
  (B) BAO + compressed CMB: parameters (h, omega_m, f). r_d(omega_m) and r_s(z_*; omega_m) from the no-dark-energy integral of
      everpresent/rdcalc.py (omega_b = 0.02237, z_d scaled with Eisenstein-Hu, normalised to CAMB at Planck as rd_planck.json),
      z_* scaled with omega_m by Hu-Sugiyama (theta_star.py); D_M(z_*) with the modified history (radiation included);
      100 theta_* = 100 r_s(z_*) / D_M(z_*), calibrated so that the Planck LCDM history gives CAMB's value (rd_planck.json,
      the "calibrated" variant of theta_star.py). chi^2 adds [(100 theta_* - 1.04110) / 0.00031]^2 (theta_star.py, Planck 2018
      Table 2) and [(omega_m - 0.14301) / 0.0011]^2 (omega_m incl. nu of rd_planck.json; the error is not in the repo files:
      Planck 2018 VI, arXiv:1807.06209, Table 2, TT,TE,EE+lowE+lensing, Omega_m h^2 = 0.1430 +- 0.0011, checked on arXiv).
      (h, omega_m) by least squares (trf, bounds) on the whitened residuals.
Per cell: best f, Delta chi^2 = chi^2(f_best) - chi^2(f = 0) (<= 0), the 95 % interval of f (Delta chi^2 = 3.84 from the
minimum, connected interval around the best fit, edges / H^2 > 0 wall flagged), Delta chi^2 at f = 1, and in (B) H0 and
Omega_m of the best fit. For comparison, Delta chi^2 of flat w0waCDM with the same data: (A) cuesta.fit_both (same priors);
(B) (h, omega_m, w0, wa) with the same CMB terms, multi-start least squares with the DESI priors w0 in [-3, 1], wa in [-3, 2],
w0 + wa < 0.
Look-elsewhere: N_UNI = 2000 noise-only LCDM universes, the noise vectors of cuesta.py (seed 20261105, official covariance);
truth = the LCDM best fit of each treatment ((A): Omega_m, h r_d; (B): h, omega_m, which also sets the "observed" theta_* and
omega_m, perturbed with their Planck errors, seed SEED_CMB). Same sweep per universe (coarse f grid, step 0.25, plus bounded
refinement); distribution of the best Delta chi^2 over (z_t, w, f); global p-value = fraction of universes with a best
Delta chi^2 at least as negative as the observed one; local p-values per cell as well.

    ../venv/bin/python escalon.py        # writes escalon.{json,txt,png} and escalon_universos.npz; nothing in results/ or papers
"""

from __future__ import annotations

import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.interpolate import CubicSpline
from scipy.linalg import solve_triangular
from scipy.optimize import brentq, least_squares, minimize_scalar

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "cuesta"))
sys.path.insert(0, str(HERE.parent / "everpresent"))
import cuesta as CU  # noqa: E402  (data, Simpson grid, LCDM / w0wa fitters, noise seed)
import rdcalc as RC  # noqa: E402  (c_s, rd_standard, zd, NORM, OMB; reads rd_planck.json)

C_KMS = CU.C_KMS
RDJ = RC.RDJ
ZT_GRID = (0.3, 0.5, 0.7, 1.0, 1.3, 1.6, 2.0, 2.4, 3.0)
W_GRID = (0.05, 0.3)
F_FINE = np.round(np.arange(-2.0, 3.0 + 1e-9, 0.05), 3)
F_COARSE = np.round(np.arange(-2.0, 3.0 + 1e-9, 0.25), 3)
F_LO, F_HI = -2.0, 3.0
DCHI2_95 = 3.84
N_UNI = CU.N_UNI
SEED_CMB = 20261111
# Planck 2018, TT,TE,EE+lowE+lensing (arXiv:1807.06209 Table 2)
TH_P, SIG_TH = 1.04110, 0.00031  # as everpresent/cmb_real/theta_star.py
WM_P = RDJ["omega_m_incl_nu"]  # 0.14301 (omega_b + omega_c + omega_nu of ou_cmb.PLANCK18)
SIG_WM = 0.0011  # Omega_m h^2 = 0.1430 +- 0.0011; not in rd_planck.json / theta_star.py, taken from Table 2 (checked)
H0_P = RDJ["H0"]
OMR = RDJ["omega_r"]
ZSTAR = RDJ["zstar"]
TH_CAMB = RDJ["thetastar100"]
BOUNDS_B = ([0.4, 0.05], [1.0, 0.30])  # h, omega_m
BOUNDS_W = ([0.4, 0.05, -3.0, -3.0], [1.0, 0.30, 1.0, 2.0])  # h, omega_m, w0, wa (DESI DR2 priors on w0, wa)
W_STARTS = ((-1.0, 0.0), (-0.5, -1.6), (-1.2, 0.8), (-0.75, -1.0), (-0.3, -2.5), (-0.9, -0.3))
BIG = 1e8

# ----------------------------------------------------------------------------------------------------------------------
# data and grids (low z from cuesta.py: 7 Simpson segments [z_{k-1}, z_k], 128 subintervals; high z: Simpson in u = ln(1+z))
# ----------------------------------------------------------------------------------------------------------------------
OBS, COV, CINV, L_COV = CU.OBS, CU.COV, CU.CINV, CU.L_COV
Z, ZU, ROW_Z, SEG, ZG, WG = CU.Z, CU.ZU, CU.ROW_Z, CU.SEG, CU.ZG, CU.WG
IS_DM, IS_DH = CU.IS_DM, CU.IS_DH
_1PZ, _1PZ3, _1PZ4 = 1.0 + ZG, (1.0 + ZG) ** 3, (1.0 + ZG) ** 4
_1PZU, _1PZU3, _1PZU4 = 1.0 + ZU, (1.0 + ZU) ** 3, (1.0 + ZU) ** 4
N_HI = 4001
U_HI = np.linspace(np.log1p(ZU[-1]), np.log1p(ZSTAR), N_HI)
_W_HI = np.ones(N_HI)
_W_HI[1:-1:2], _W_HI[2:-1:2] = 4.0, 2.0
_W_HI *= (U_HI[-1] - U_HI[0]) / (N_HI - 1) / 3.0
_EU, _E3U, _E4U = np.exp(U_HI), np.exp(3 * U_HI), np.exp(4 * U_HI)
Z_HI = _EU - 1.0


def step_s(z, zt: float, w: float):
    return 0.5 * (1.0 + np.tanh((np.log1p(z) - np.log1p(zt)) / w))


class Cell:
    """s(z) of one (z_t, w) on the three grids and at z = 0."""

    def __init__(self, zt: float, w: float):
        self.zt, self.w = zt, w
        self.s_lo, self.s_u, self.s_hi, self.s0 = step_s(ZG, zt, w), step_s(ZU, zt, w), step_s(Z_HI, zt, w), float(step_s(0.0, zt, w))
        self.name = f"zt_{zt:g}|w_{w:g}"

    def shapes(self, f: float):
        n = 1.0 - f * self.s0
        return (1.0 - f * self.s_lo) / n, (1.0 - f * self.s_u) / n, (1.0 - f * self.s_hi) / n


CELLS = [Cell(zt, w) for w in W_GRID for zt in ZT_GRID]


def cpl_shapes(w0: float, wa: float):
    p = 3.0 * (1.0 + w0 + wa)
    g_lo = _1PZ**p * np.exp(-3.0 * wa * ZG / _1PZ)
    g_u = _1PZU**p * np.exp(-3.0 * wa * ZU / _1PZU)
    g_hi = np.exp(p * U_HI) * np.exp(-3.0 * wa * (1.0 - 1.0 / _EU))
    return g_lo, g_u, g_hi


def dist_lo(om: float, orad: float, g_lo, g_u):
    """Dimensionless D_type H0/c at the 13 rows and the cumulative D_M at the 7 z; None if H^2 <= 0 somewhere."""
    ode = 1.0 - om - orad
    e2 = om * _1PZ3 + orad * _1PZ4 + ode * g_lo
    e2u = om * _1PZU3 + orad * _1PZU4 + ode * g_u
    if e2.min() <= 0.0 or e2u.min() <= 0.0:
        return None, None, min(e2.min(), e2u.min())
    dm_u = np.cumsum(np.bincount(SEG, weights=WG / np.sqrt(e2), minlength=len(ZU)))
    dh_u = 1.0 / np.sqrt(e2u)
    dm, dh = dm_u[ROW_Z], dh_u[ROW_Z]
    return np.where(IS_DM, dm, np.where(IS_DH, dh, np.cbrt(Z * dm * dm * dh))), dm_u, 0.0


def dist_hi(om: float, orad: float, g_hi):
    """Dimensionless int_{z_7}^{z_*} dz / E(z) and 1/E(z_*); None if H^2 <= 0."""
    e2 = om * _E3U + orad * _E4U + (1.0 - om - orad) * g_hi
    if e2.min() <= 0.0:
        return None, None, e2.min()
    e = np.sqrt(e2)
    return float(np.sum(_W_HI * _EU / e)), float(1.0 / e[-1]), 0.0


# ----------------------------------------------------------------------------------------------------------------------
# sound horizons and z_* as functions of omega_m (splines on quad results)
# ----------------------------------------------------------------------------------------------------------------------
def zstar_hs(wm: float, wb: float = RC.OMB) -> float:
    """Hu & Sugiyama (1996) fitting formula for z_* (copied from theta_star.py), used only as a ratio to scale CAMB's z_*."""
    g1 = 0.0783 * wb**-0.238 / (1 + 39.5 * wb**0.763)
    g2 = 0.560 / (1 + 21.1 * wb**1.81)
    return 1048 * (1 + 0.00124 * wb**-0.738) * (1 + g1 * wm**g2)


def zstar_of(wm: float) -> float:
    return ZSTAR * zstar_hs(wm) / zstar_hs(WM_P)


_WM_GRID = np.linspace(0.08, 0.26, 361)
RD_SPL = CubicSpline(_WM_GRID, [RC.NORM * RC.rd_standard(w, float(RC.zd(w))) for w in _WM_GRID])
RS_SPL = CubicSpline(_WM_GRID, [RC.NORM * RC.rd_standard(w, zstar_of(w)) for w in _WM_GRID])


def theta_raw(h: float, wm: float, dm7: float, hi: float, inv_e_star: float) -> float:
    """100 r_s(z_*) / D_M(z_*), uncalibrated. dm7, hi: dimensionless integrals to z_7 and from z_7 to CAMB's z_*;
    first-order correction for z_*(omega_m) - z_*^CAMB."""
    dm_star = C_KMS / (100.0 * h) * (dm7 + hi + inv_e_star * (zstar_of(wm) - ZSTAR))
    return 100.0 * float(RS_SPL(wm)) / dm_star


def _cal() -> float:
    om, orad = WM_P / (H0_P / 100) ** 2, OMR / (H0_P / 100) ** 2
    _, dm_u, _ = dist_lo(om, orad, 1.0, 1.0)
    hi, ie, _ = dist_hi(om, orad, 1.0)
    return theta_raw(H0_P / 100, WM_P, dm_u[-1], hi, ie) / TH_CAMB


CAL = _cal()  # theta_* of the Planck LCDM history through this code over CAMB's


# ----------------------------------------------------------------------------------------------------------------------
# treatment (A): BAO only
# ----------------------------------------------------------------------------------------------------------------------
def chi2_A_om(om: float, g_lo, g_u, y) -> tuple[float, float]:
    d, _, bad = dist_lo(om, 0.0, g_lo, g_u)
    if d is None:
        return BIG * (1.0 - bad), np.nan
    cd = CINV @ d
    s = (cd @ y) / (cd @ d)
    r = y - s * d
    return float(r @ CINV @ r), float(C_KMS / (100.0 * s))


_OM_GRID = np.linspace(0.05, 0.95, 19)


def fit_A(f: float, cell: Cell, y) -> tuple[float, float, float]:
    """chi^2 minimised over Omega_m (grid + bounded Brent) and h r_d (analytic) at fixed f. Returns (chi2, Omega_m, h r_d)."""
    g_lo, g_u, _ = cell.shapes(f)
    fun = lambda om: chi2_A_om(om, g_lo, g_u, y)[0]
    vals = [fun(o) for o in _OM_GRID]
    j = int(np.argmin(vals))
    if vals[j] >= BIG:
        return np.inf, np.nan, np.nan
    lo, hi = _OM_GRID[max(j - 1, 0)], _OM_GRID[min(j + 1, len(_OM_GRID) - 1)]
    r = minimize_scalar(fun, bounds=(lo, hi), method="bounded", options={"xatol": 1e-7})
    c, hrd = chi2_A_om(float(r.x), g_lo, g_u, y)
    if c >= BIG:
        return np.inf, np.nan, np.nan
    return c, float(r.x), hrd


# ----------------------------------------------------------------------------------------------------------------------
# treatment (B): BAO + compressed CMB
# ----------------------------------------------------------------------------------------------------------------------
def resid_B(p, shapes, y, th_obs: float, wm_obs: float) -> np.ndarray:
    h, wm = p[0], p[1]
    g_lo, g_u, g_hi = shapes(p)
    om, orad = wm / h**2, OMR / h**2
    d, dm_u, bad1 = dist_lo(om, orad, g_lo, g_u)
    hi, ie, bad2 = dist_hi(om, orad, g_hi)
    if d is None or hi is None:
        return np.full(15, 1e3 * (1.0 + abs(min(bad1, bad2))))
    pred = C_KMS / (100.0 * h * float(RD_SPL(wm))) * d
    r = solve_triangular(L_COV, y - pred, lower=True)
    th = theta_raw(h, wm, dm_u[-1], hi, ie) / CAL
    return np.concatenate([r, [(th - th_obs) / SIG_TH, (wm - wm_obs) / SIG_WM]])


def _ls(fun, x0, bounds, args):
    return least_squares(fun, np.asarray(x0, float), bounds=bounds, args=args, method="trf",
                         x_scale=np.array([0.01, 0.001] + [0.1] * (len(x0) - 2)), ftol=1e-10, xtol=1e-10, gtol=1e-10, max_nfev=400)


def fit_B(f: float, cell: Cell, y, th_obs: float, wm_obs: float, x0=(0.68, 0.143)) -> tuple[float, float, float]:
    """chi^2 minimised over (h, omega_m) at fixed f. Returns (chi2, h, omega_m); chi2 = inf if H^2 <= 0 at the solution."""
    sh = lambda p, f=f: cell.shapes(f)
    r = _ls(resid_B, x0, BOUNDS_B, (sh, y, th_obs, wm_obs))
    if r.fun[0] >= 1e3:
        return np.inf, np.nan, np.nan
    return float(2.0 * r.cost), float(r.x[0]), float(r.x[1])


def theta_of_B(h: float, wm: float, f: float, cell: Cell | None) -> float:
    g_lo, g_u, g_hi = cell.shapes(f) if cell is not None else (1.0, 1.0, 1.0)
    om, orad = wm / h**2, OMR / h**2
    _, dm_u, _ = dist_lo(om, orad, g_lo, g_u)
    hi, ie, _ = dist_hi(om, orad, g_hi)
    return theta_raw(h, wm, dm_u[-1], hi, ie) / CAL


def fit_w0wa_B(y, th_obs: float, wm_obs: float, x0_lcdm) -> dict:
    def res(p, y, th, wm):
        base = resid_B(p, lambda q: cpl_shapes(q[2], q[3]), y, th, wm)
        pen = 1e3 * max(0.0, p[2] + p[3])  # DESI prior w0 + wa < 0
        return np.concatenate([base, [pen]])

    best = None
    for w0, wa in W_STARTS:
        r = _ls(res, [x0_lcdm[0], x0_lcdm[1], w0, wa], BOUNDS_W, (y, th_obs, wm_obs))
        if best is None or r.cost < best.cost:
            best = r
    h, wm, w0, wa = (float(v) for v in best.x)
    return {"chi2": float(2.0 * best.cost), "H0": 100 * h, "omega_m": wm, "Omega_m": wm / h**2, "w0": w0, "wa": wa}


# ----------------------------------------------------------------------------------------------------------------------
# profiles, intervals, sweeps
# ----------------------------------------------------------------------------------------------------------------------
def profile_fun(treat: str, cell: Cell, y, th_obs=None, wm_obs=None):
    """Returns chi2(f) (float, inf if H^2 <= 0) with warm start in (B)."""
    if treat == "A":
        return lambda f: fit_A(f, cell, y)[0]
    state = {"x0": (0.68, 0.143)}

    def fun(f):
        c, h, wm = fit_B(f, cell, y, th_obs, wm_obs, state["x0"])
        if np.isfinite(c):
            state["x0"] = (h, wm)
        return c
    return fun


def refine(fun, fgrid, vals) -> tuple[float, float]:
    """Bounded minimisation between the neighbours of the grid minimum. Returns (f_best, chi2_best)."""
    fin = np.isfinite(vals)
    if not fin.any():
        return np.nan, np.inf
    j = int(np.nanargmin(np.where(fin, vals, np.inf)))
    lo = fgrid[j - 1] if j > 0 and fin[j - 1] else fgrid[j]
    hi = fgrid[j + 1] if j < len(fgrid) - 1 and fin[j + 1] else fgrid[j]
    if hi - lo < 1e-9:
        return float(fgrid[j]), float(vals[j])
    r = minimize_scalar(fun, bounds=(lo, hi), method="bounded", options={"xatol": 1e-4})
    if np.isfinite(r.fun) and r.fun < vals[j]:
        return float(r.x), float(r.fun)
    return float(fgrid[j]), float(vals[j])


def interval95(fun, fgrid, vals, f_best: float, c_best: float) -> dict:
    """Connected 95 % interval around f_best: Delta chi^2 = chi2 - c_best <= 3.84; brentq on the re-fitted profile."""
    out = {}
    g = lambda f: fun(f) - c_best - DCHI2_95
    for side, order in (("lo", range(len(fgrid) - 1, -1, -1)), ("hi", range(len(fgrid)))):
        sel = [k for k in order if (fgrid[k] < f_best if side == "lo" else fgrid[k] > f_best)]
        last_ok, first_bad, wall = f_best, None, False
        for k in (sel[::-1] if side == "lo" else sel):  # walking away from f_best
            v = vals[k]
            if not np.isfinite(v):
                first_bad, wall = fgrid[k], True
                break
            if v - c_best > DCHI2_95:
                first_bad = fgrid[k]
                break
            last_ok = fgrid[k]
        if first_bad is None:
            out[side] = {"f": F_LO if side == "lo" else F_HI, "flag": "borde"}
        elif wall:
            # H^2 <= 0 beyond last_ok: bisect to the wall, then the 3.84 crossing if it is before the wall
            a, b = (first_bad, last_ok) if side == "lo" else (last_ok, first_bad)
            for _ in range(30):
                m = 0.5 * (a + b)
                ok = np.isfinite(fun(m))
                if side == "lo":
                    a, b = (a, m) if ok else (m, b)
                else:
                    a, b = (m, b) if ok else (a, m)
            fw = b if side == "lo" else a
            if g(fw) > 0:
                out[side] = {"f": float(brentq(g, *(sorted((last_ok, fw))), xtol=1e-5)), "flag": ""}
            else:
                out[side] = {"f": float(fw), "flag": "muro H2"}
        else:
            out[side] = {"f": float(brentq(g, *(sorted((last_ok, first_bad))), xtol=1e-5)), "flag": ""}
    dv = vals - c_best
    inside = (fgrid >= out["lo"]["f"] - 1e-9) & (fgrid <= out["hi"]["f"] + 1e-9)
    out["n_fine_points_allowed_outside"] = int(np.sum(np.isfinite(dv) & (dv <= DCHI2_95) & ~inside))
    return out


def sweep(treat: str, y, th_obs=None, wm_obs=None, fgrid=F_COARSE, full: bool = False) -> dict:
    """All cells. full=True: fine grid, intervals, f = 1, H0 (real data); else best f and Delta chi^2 only (universes)."""
    if treat == "A":
        c0 = CU.fit_lcdm(y)[0]
    else:
        c0, h0, wm0 = fit_B(0.0, CELLS[0], y, th_obs, wm_obs)
    res = {"chi2_0": float(c0), "cells": {}}
    if treat == "B":
        res["lcdm"] = {"H0": 100 * h0, "omega_m": wm0, "Omega_m": wm0 / h0**2}
    for cell in CELLS:
        fun = profile_fun(treat, cell, y, th_obs, wm_obs)
        vals = np.array([fun(f) for f in fgrid])
        fb, cb = refine(fun, fgrid, vals)
        d = {"f_best": fb, "dchi2": cb - c0}
        if full:
            d["profile"] = [None if not np.isfinite(v) else float(v - c0) for v in vals]
            d["interval95"] = interval95(fun, fgrid, vals, fb, cb)
            d["dchi2_f1"] = float(fun(1.0) - c0) if np.isfinite(fun(1.0)) else None
            d["f1_excluded95"] = (d["dchi2_f1"] is None) or (d["dchi2_f1"] - (cb - c0) > DCHI2_95)
            d["chi2_at_f0_check"] = float(fun(0.0) - c0)
            if treat == "A":
                _, om, hrd = fit_A(fb, cell, y)
                d.update(Omega_m=om, h_rd=hrd)
            else:
                _, h, wm = fit_B(fb, cell, y, th_obs, wm_obs)
                d.update(H0=100 * h, omega_m=wm, Omega_m=wm / h**2, theta100=theta_of_B(h, wm, fb, cell))
        res["cells"][cell.name] = d
    return res


# ----------------------------------------------------------------------------------------------------------------------
# checks
# ----------------------------------------------------------------------------------------------------------------------
def check_distances() -> dict:
    """Grid distances against quad for a few (cell, f, Omega_m), low z (13 rows) and D_M(z_*)."""
    worst_lo = worst_hi = 0.0
    for cell, f, om in ((CELLS[0], 1.0, 0.3), (CELLS[9], 2.0, 0.3), (CELLS[4], -1.5, 0.25), (CELLS[16], 0.8, 0.32)):
        g_lo, g_u, g_hi = cell.shapes(f)
        orad = OMR / 0.68**2
        d, dm_u, _ = dist_lo(om, orad, g_lo, g_u)
        hi, _, _ = dist_hi(om, orad, g_hi)
        n = 1.0 - f * cell.s0
        E = lambda x: np.sqrt(om * (1 + x) ** 3 + orad * (1 + x) ** 4 + (1 - om - orad) * (1 - f * step_s(x, cell.zt, cell.w)) / n)
        for i, (k, z) in enumerate(zip(CU.KINDS, Z)):
            dm = quad(lambda x: 1 / E(x), 0, z, epsrel=1e-12, points=[cell.zt], limit=200)[0]
            ref = {"DM": dm, "DH": 1 / E(z), "DV": (z * dm * dm / E(z)) ** (1 / 3)}[k]
            worst_lo = max(worst_lo, abs(d[i] / ref - 1))
        ref_hi = quad(lambda x: 1 / E(x), ZU[-1], ZSTAR, epsrel=1e-12, points=[cell.zt] if cell.zt > ZU[-1] else None, limit=400)[0]
        worst_hi = max(worst_hi, abs((dm_u[-1] + hi) / (quad(lambda x: 1 / E(x), 0, ZU[-1], epsrel=1e-12, limit=200)[0] + ref_hi) - 1))
    return {"grid_vs_quad_max_rel_13rows": worst_lo, "grid_vs_quad_max_rel_DM_zstar": worst_hi}


def check_sound_horizons() -> dict:
    rd_aub = lambda wm: 55.154 * np.exp(-72.3 * (RDJ["rd_aubourg"]["omega_nu"] + 0.0006) ** 2) / ((wm - RDJ["rd_aubourg"]["omega_nu"]) ** 0.25351 * RC.OMB**0.12807)
    return {"rd_at_planck": float(RD_SPL(WM_P)), "rd_camb": RDJ["rdrag_Mpc"], "rs_zstar_at_planck": float(RS_SPL(WM_P)),
            "rs_zstar_camb": RDJ["rs_zstar_Mpc"], "rel_rs": float(RS_SPL(WM_P) / RDJ["rs_zstar_Mpc"] - 1),
            "theta_calibration_CAL": CAL, "theta_raw_minus_camb_in_sigma": (CAL - 1) * TH_CAMB / SIG_TH,
            "dlnrd_dlnwm_here": float(RD_SPL.derivative()(WM_P) * WM_P / RD_SPL(WM_P)),
            "dlnrd_dlnwm_aubourg": float((np.log(rd_aub(WM_P * 1.01)) - np.log(rd_aub(WM_P * 0.99))) / (np.log(1.01) - np.log(0.99))),
            "zstar_planck": ZSTAR, "zstar_at_wm_plus_1pct": zstar_of(WM_P * 1.01)}


# ----------------------------------------------------------------------------------------------------------------------
# universes
# ----------------------------------------------------------------------------------------------------------------------
_U = {}


def _worker(i: int):
    yA = _U["predA"] + _U["noise"][i]
    yB = _U["predB"] + _U["noise"][i]
    th = _U["thB"] + SIG_TH * _U["cmb"][i, 0]
    wm = _U["wmB"] + SIG_WM * _U["cmb"][i, 1]
    rA = sweep("A", yA)
    rB = sweep("B", yB, th, wm)
    out = []
    for r in (rA, rB):
        out.append([r["cells"][c.name]["dchi2"] for c in CELLS])
        out.append([r["cells"][c.name]["f_best"] for c in CELLS])
    return np.array(out)  # (4, 18): dchi2_A, f_A, dchi2_B, f_B


def main() -> None:
    t0 = time.time()
    checks = {**check_distances(), **check_sound_horizons()}
    print("comprobaciones:", json.dumps(checks, indent=1), flush=True)

    # real data
    realA = sweep("A", OBS, fgrid=F_FINE, full=True)
    cA = CU.fit_both(OBS)
    realA["w0wa"] = {"dchi2": -cA[0], "chi2": cA[4], "Omega_m": cA[5], "h_rd": cA[6], "w0": cA[7], "wa": cA[8]}
    print(f"(A) chi2_0 = {realA['chi2_0']:.3f}; w0wa dchi2 = {-cA[0]:.3f}  [{time.time() - t0:.0f} s]", flush=True)
    realB = sweep("B", OBS, TH_P, WM_P, fgrid=F_FINE, full=True)
    lB = realB["lcdm"]
    wB = fit_w0wa_B(OBS, TH_P, WM_P, (lB["H0"] / 100, lB["omega_m"]))
    realB["w0wa"] = {**wB, "dchi2": wB["chi2"] - realB["chi2_0"]}
    realB["lcdm"]["theta100"] = theta_of_B(lB["H0"] / 100, lB["omega_m"], 0.0, None)
    realB["lcdm"]["chi2_parts"] = {"theta": ((realB["lcdm"]["theta100"] - TH_P) / SIG_TH) ** 2, "omega_m": ((lB["omega_m"] - WM_P) / SIG_WM) ** 2}
    print(f"(B) LCDM: H0 = {lB['H0']:.2f}, omega_m = {lB['omega_m']:.5f}, chi2_0 = {realB['chi2_0']:.3f}; "
          f"w0wa: dchi2 = {realB['w0wa']['dchi2']:.3f}, H0 = {wB['H0']:.2f}, w0 = {wB['w0']:.3f}, wa = {wB['wa']:.3f}  [{time.time() - t0:.0f} s]", flush=True)
    for name, r in (("A", realA), ("B", realB)):
        for c in CELLS:
            d = r["cells"][c.name]
            iv = d["interval95"]
            extra = f"  H0 = {d['H0']:.2f}" if name == "B" else ""
            print(f"({name}) {c.name:14s} f = {d['f_best']:+.3f}  dchi2 = {d['dchi2']:+.3f}  95%: [{iv['lo']['f']:+.3f}{iv['lo']['flag'][:1]}, "
                  f"{iv['hi']['f']:+.3f}{iv['hi']['flag'][:1]}]  dchi2(f=1) = {d['dchi2_f1']}  f=1 excl: {d['f1_excluded95']}{extra}", flush=True)
    bestA = min(realA["cells"].values(), key=lambda d: d["dchi2"])["dchi2"]
    bestB = min(realB["cells"].values(), key=lambda d: d["dchi2"])["dchi2"]
    # the same observed sweep with the coarse grid of the universes (for the p-value on equal footing)
    obsA_c = min(d["dchi2"] for d in sweep("A", OBS)["cells"].values())
    obsB_c = min(d["dchi2"] for d in sweep("B", OBS, TH_P, WM_P)["cells"].values())

    # universes
    rng_n = np.random.default_rng(CU.SEED_NOISE)
    noise = rng_n.standard_normal((N_UNI, 13)) @ L_COV.T
    cmb = np.random.default_rng(SEED_CMB).standard_normal((N_UNI, 2))
    _, omA, hrdA = CU.fit_lcdm(OBS)
    predA = C_KMS / (100.0 * hrdA) * CU.dhat13(omA, -1.0, 0.0)
    hB, wmB = lB["H0"] / 100, lB["omega_m"]
    dB, _, _ = dist_lo(wmB / hB**2, OMR / hB**2, 1.0, 1.0)
    predB = C_KMS / (100.0 * hB * float(RD_SPL(wmB))) * dB
    _U.update(noise=noise, cmb=cmb, predA=predA, predB=predB, thB=realB["lcdm"]["theta100"], wmB=wmB)
    t1 = time.time()
    with Pool() as pool:
        U = np.array(pool.map(_worker, range(N_UNI), chunksize=10))  # (N_UNI, 4, 18)
    print(f"universos: {time.time() - t1:.0f} s", flush=True)
    np.savez_compressed(HERE / "escalon_universos.npz", cells=np.array([c.name for c in CELLS]),
                        columns=np.array(["dchi2_A", "f_A", "dchi2_B", "f_B"]), data=U)
    uni = {}
    for name, k, obs_fine, obs_c in (("A", 0, bestA, obsA_c), ("B", 2, bestB, obsB_c)):
        best = U[:, k, :].min(axis=1)
        jbest = U[:, k, :].argmin(axis=1)
        uni[name] = {"observed_best_dchi2_fine": obs_fine, "observed_best_dchi2_coarse": obs_c,
                     "p_global": float(np.mean(best <= obs_c)), "p_global_vs_fine": float(np.mean(best <= obs_fine)),
                     "p_err": float(np.sqrt(np.mean(best <= obs_c) * (1 - np.mean(best <= obs_c)) / N_UNI)),
                     "best_dchi2_percentiles_5_50_95": [float(np.percentile(best, q)) for q in (5, 50, 95)],
                     "best_dchi2_mean": float(best.mean()), "frac_best_le_m6": float(np.mean(best <= -6)),
                     "frac_best_le_m3.84": float(np.mean(best <= -3.84)),
                     "best_cell_counts": {c.name: int(np.sum(jbest == j)) for j, c in enumerate(CELLS)},
                     "local": {c.name: {"p_local": float(np.mean(U[:, k, j] <= (realA if name == "A" else realB)["cells"][c.name]["dchi2"])),
                                        "median_dchi2": float(np.median(U[:, k, j])), "p5_dchi2": float(np.percentile(U[:, k, j], 5))}
                               for j, c in enumerate(CELLS)}}
        print(f"({name}) observado {obs_c:+.3f} (fino {obs_fine:+.3f}); p global = {uni[name]['p_global']:.4f}; "
              f"percentiles 5/50/95 del mejor dchi2: {uni[name]['best_dchi2_percentiles_5_50_95']}", flush=True)

    out = {"model": "rho_DE/rho_crit0 = Omega_L [1 - f s(z)], s = 1/2 [1 + tanh((ln(1+z) - ln(1+z_t))/w)], flat, Omega_L from flatness today",
           "zt_grid": ZT_GRID, "w_grid": W_GRID, "f_fine_step": 0.05, "f_coarse_step": 0.25, "f_range": [F_LO, F_HI],
           "planck": {"theta100": TH_P, "sigma_theta100": SIG_TH, "omega_m": WM_P, "sigma_omega_m": SIG_WM,
                      "source": "arXiv:1807.06209 Table 2, TT,TE,EE+lowE+lensing (theta: theta_star.py; omega_m error: Table 2, not in the repo files)",
                      "omega_b": RC.OMB, "omega_r": OMR, "zstar_camb": ZSTAR, "theta100_camb": TH_CAMB},
           "n_universes": N_UNI, "seeds": {"noise": CU.SEED_NOISE, "cmb": SEED_CMB}, "checks": checks,
           "real": {"A": realA, "B": realB}, "universes": uni, "runtime_s": time.time() - t0}
    (HERE / "escalon.json").write_text(json.dumps(out, indent=1), encoding="utf-8")

    lines = [f"escalon: f in [{F_LO}, {F_HI}], paso fino {0.05}; dchi2 = chi2(f_best) - chi2(f = 0); intervalo 95 %: dchi2 = 3.84 desde el minimo",
             f"(A) BAO solo: chi2(LCDM) = {realA['chi2_0']:.3f}; w0wa: dchi2 = {realA['w0wa']['dchi2']:+.3f}",
             f"(B) BAO + CMB: chi2(LCDM) = {realB['chi2_0']:.3f}, H0 = {lB['H0']:.2f}, omega_m = {lB['omega_m']:.5f}, Omega_m = {lB['Omega_m']:.4f}; "
             f"w0wa: dchi2 = {realB['w0wa']['dchi2']:+.3f}, H0 = {wB['H0']:.2f}, w0 = {wB['w0']:.2f}, wa = {wB['wa']:.2f}", ""]
    for name, r in (("A", realA), ("B", realB)):
        lines.append(f"({name}) {'celda':14s} {'f_best':>7s} {'dchi2':>7s} {'f_lo':>8s} {'f_hi':>8s} {'dchi2(f=1)':>10s} {'f=1 excl':>8s}" + ("  H0     Om" if name == "B" else "  Om     h rd"))
        for c in CELLS:
            d = r["cells"][c.name]
            iv = d["interval95"]
            d1 = f"{d['dchi2_f1']:10.2f}" if d["dchi2_f1"] is not None else f"{'H2<0':>10s}"
            tail = f"  {d['H0']:.2f} {d['Omega_m']:.4f}" if name == "B" else f"  {d['Omega_m']:.4f} {d['h_rd']:.2f}"
            lines.append(f"    {c.name:14s} {d['f_best']:+7.3f} {d['dchi2']:+7.2f} {iv['lo']['f']:+7.3f}{(iv['lo']['flag'] or ' ')[0]} "
                         f"{iv['hi']['f']:+7.3f}{(iv['hi']['flag'] or ' ')[0]} {d1} {str(d['f1_excluded95']):>8s}{tail}")
        u = uni[name]
        lines.append(f"    universos: observado {u['observed_best_dchi2_coarse']:+.3f}, p global = {u['p_global']:.4f} +- {u['p_err']:.4f}; "
                     f"mejor dchi2 p5/p50/p95 = {u['best_dchi2_percentiles_5_50_95'][0]:.2f}/{u['best_dchi2_percentiles_5_50_95'][1]:.2f}/{u['best_dchi2_percentiles_5_50_95'][2]:.2f}")
        lines.append("")
    lines.append("comprobaciones: " + json.dumps(checks))
    (HERE / "escalon.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    figure(out)


def figure(out: dict) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    col = {"A": "#2a78d6", "B": "#eb6834"}
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.2), sharex=True)
    for j, w in enumerate(W_GRID):
        ax, axd = axes[0, j], axes[1, j]
        for name in ("A", "B"):
            r = out["real"][name]["cells"]
            lo = np.array([r[f"zt_{zt:g}|w_{w:g}"]["interval95"]["lo"]["f"] for zt in ZT_GRID])
            hi = np.array([r[f"zt_{zt:g}|w_{w:g}"]["interval95"]["hi"]["f"] for zt in ZT_GRID])
            fb = np.array([r[f"zt_{zt:g}|w_{w:g}"]["f_best"] for zt in ZT_GRID])
            dc = np.array([r[f"zt_{zt:g}|w_{w:g}"]["dchi2"] for zt in ZT_GRID])
            lab = {"A": "(A) BAO solo", "B": "(B) BAO + CMB comprimido"}[name]
            ax.fill_between(ZT_GRID, lo, hi, color=col[name], alpha=0.18, lw=0)
            ax.plot(ZT_GRID, lo, color=col[name], lw=1.2)
            ax.plot(ZT_GRID, hi, color=col[name], lw=1.2)
            ax.plot(ZT_GRID, fb, color=col[name], lw=2, marker="o", ms=5, label=lab)
            axd.plot(ZT_GRID, dc, color=col[name], lw=2, marker="o", ms=5, label=lab)
            dw = out["real"][name]["w0wa"]["dchi2"]
            axd.axhline(dw, color=col[name], lw=1, ls="--")
            axd.text(ZT_GRID[0], dw, f" w₀wₐ {name}: {dw:+.1f}", color=col[name], va="bottom", fontsize=8)
        ax.axhline(1.0, color="#52514e", lw=1, ls=":")
        ax.axhline(0.0, color="#52514e", lw=1, ls="-", alpha=0.5)
        ax.text(ZT_GRID[0], 1.0, " f = 1 (sin Λ antes de z_t)", color="#52514e", va="bottom", ha="left", fontsize=8)
        ax.set_title(f"w = {w:g}: f de mejor ajuste e intervalo al 95 %")
        ax.set_ylabel("f")
        ax.set_ylim(F_LO - 0.1, F_HI + 0.1)
        axd.set_ylabel("Δχ² = χ²(f̂) − χ²(f = 0)")
        axd.set_xlabel("z_t")
        axd.set_title(f"w = {w:g}: mejora frente a ΛCDM")
        for a in (ax, axd):
            a.grid(alpha=0.25)
            a.set_xticks(ZT_GRID)
        ax.legend(loc="lower left", fontsize=8, frameon=False)
    fig.suptitle("Escalón tardío de la energía oscura, ρ_DE = Ω_Λ[1 − f s(z)]: DESI DR2 (13 números) solo y con θ*, ω_m de Planck")
    fig.tight_layout()
    fig.savefig(HERE / "escalon.png", dpi=150)


if __name__ == "__main__":
    main()
