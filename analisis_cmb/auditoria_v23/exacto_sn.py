#!/usr/bin/env python3
"""
Exploration only (EXACTO_SN.md): how much do the ceilings of EXACTO.md drop with the Pantheon+ supernovae? Nothing in
results/ or papers changes.

Set-up of exacto.py unchanged (convention B, rho_DE = Omega_L + X(z) - X(0); same grid in u = ln(1+z) with the DESI z as exact
nodes; same seed per cell; discard if H^2 <= 0 at the LCDM best fit of the 15 rows; criterion P(chi^2_min <= chi^2_obs) = 5 %,
interpolated in ln sigma_X), plus:
- SN data as analisis_cmb/sn/perfil_sn.py: Pantheon+ (datos/Pantheon+SH0ES.dat, STAT+SYS covariance), m_b_corr with
  z_HD > 0.01 and IS_CALIBRATOR = 0 (1580 SN); M marginalised analytically (flat prior):
  chi^2_SN = r'C^-1 r - (1'C^-1 r)^2 / (1'C^-1 1) (model-independent constants dropped; they cancel in P).
- SN model: m = 5 log10[(1 + z_HEL) D_M(z_HD)] + M, D_M from the same exact history (the H^2 of exacto.py), interpolated
  linearly in u on the grid (the SN z are not nodes; the interpolation error is < 1e-6 relative).
- Cells: those of exacto.py (sigma_X in {0.05, 0.1, 0.2, 0.3, 0.5, 0.7}) plus sigma_X = 0.01, 0.02, 0.03; theta in {0.1, 0.3, 1, 3, 10}.
  New cells get seeds [SEED, 100 + index]; the old ones keep exacto.py's.
- Combinations, each with its own chi^2_obs (fit with X = 0 to the same data): SN alone (h degenerate with M: Omega_m fitted,
  h fixed at 0.68 for Omega_r only); BAO + SN (h, omega_m); BAO + theta_* + omega_m + SN (h, omega_m). In the least-squares
  fits M is a third parameter (equivalent to the analytic marginalisation at the minimum). Variant: the last combination with
  z_HD > 0.023 (1371 SN).
Checks: (1) with X = 0, SN alone must give Omega_m = 0.3326 (SN.md); (2) grid: the theta = 10 cells repeated with step 1e-3 in u
below z = 0.3 (same seeds, different nodes: P must agree within the Monte Carlo error); (3) linear marginalised bounds
(-1.92 and -1.35) with the full-response kernel of convention B, as cmb_fondo.py, for the three combinations (M as a
linear nuisance template on the SN rows), theta in {0.1, 0.3, 1, 3, 10}.
Also: fraction of survivors with rho_DE <= 0 somewhere on the grid (rho_DE = Omega_L + X - X(0) at the LCDM best fit) at the
new ceilings of BAO + theta_* + omega_m + SN (paths regenerated at the ceiling with the seed of the nearest cell).

    ../venv/bin/python exacto_sn.py     # writes exacto_sn.{json,txt} and exacto_sn_chi2.npz
"""

from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")  # one BLAS thread per pool process (1580 x 1580 SN products)

import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy.linalg import solve_triangular
from scipy.optimize import least_squares, minimize_scalar
from scipy.signal import lfilter

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import exacto as X  # noqa: E402
import cmb_fondo as CF  # noqa: E402

E = X.E
C_KMS, OMR, ZSTAR = X.C_KMS, X.OMR, X.ZSTAR
SIGMAS = (0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7)
THETAS = (0.1, 0.3, 1.0, 3.0, 10.0)
N_REAL = X.N_REAL
SEED = X.SEED
MAG = 5.0 / np.log(10.0)
H_SN = 0.68
SN_DIR = HERE.parent / "sn" / "datos"
BOUNDS_P = X.BOUNDS_P
COMBOS = ("SN", "BAO+SN", "BAO+CMB+SN", "BAO+CMB+SN_z023")


# ----------------------------------------------------------------------------------------------------------------------
# grids (exacto's, and the fine variant below z = 0.3)
# ----------------------------------------------------------------------------------------------------------------------
class Grid:
    def __init__(self, du_low: float | None = None):
        u = X._u.copy()
        if du_low is not None:
            u = np.unique(np.concatenate([u[u > np.log1p(0.3)], np.arange(0.0, np.log1p(0.3) + 1e-12, du_low)]))
        self.U = np.unique(np.concatenate([u, np.log1p(np.unique(X.Z13))]))
        self.Z = np.expm1(self.U)
        self.N = len(self.U)
        self.DUK = np.diff(self.U)
        self.IDX = [int(np.argmin(np.abs(self.U - np.log1p(z)))) for z in X.Z13]
        assert max(abs(self.U[i] - np.log1p(z)) for i, z in zip(self.IDX, X.Z13)) < 1e-12
        self.p3, self.p4, self.p1 = (1 + self.Z) ** 3, (1 + self.Z) ** 4, 1 + self.Z


G0 = Grid()
GF = Grid(du_low=1e-3)


def ou_paths(g: Grid, sigma: float, theta: float, n: int, rng) -> np.ndarray:
    Xp = np.empty((n, g.N))
    Xp[:, 0] = sigma * rng.standard_normal(n)
    rho = np.exp(-theta * g.DUK)
    s = sigma * np.sqrt(1.0 - rho**2)
    xi = rng.standard_normal((n, g.N - 1))
    for k in range(g.N - 1):
        Xp[:, k + 1] = rho[k] * Xp[:, k] + s[k] * xi[:, k]
    return Xp


# ----------------------------------------------------------------------------------------------------------------------
# SN data
# ----------------------------------------------------------------------------------------------------------------------
def load_sn():
    dat, cov = SN_DIR / "Pantheon+SH0ES.dat", SN_DIR / "Pantheon+SH0ES_STAT+SYS.cov"
    with open(dat) as f:
        head = f.readline().split()
    raw = np.genfromtxt(dat, skip_header=1, dtype=None, encoding="ascii", names=head)
    n = int(open(cov).readline())
    C = np.loadtxt(cov, skiprows=1).reshape(n, n)
    assert n == len(raw)
    out = {}
    for lab, zmin in (("all", 0.01), ("z023", 0.023)):
        m = (raw["zHD"] > zmin) & (raw["IS_CALIBRATOR"] == 0)
        idx = np.nonzero(m)[0]
        Cs = C[np.ix_(idx, idx)]
        L = np.linalg.cholesky(Cs)
        Ci = np.linalg.inv(Cs)
        one = np.ones(len(idx))
        out[lab] = {"zhd": raw["zHD"][idx].astype(float), "zhel": raw["zHEL"][idx].astype(float), "mb": raw["m_b_corr"][idx].astype(float),
                    "L": L, "Ci": Ci, "Ci1": Ci @ one, "e": float(one @ Ci @ one), "n": len(idx), "uhd": np.log1p(raw["zHD"][idx].astype(float))}
    return out


SN = load_sn()


def chi2_sn(sn: dict, mu: np.ndarray) -> float:
    r = sn["mb"] - mu
    cr = sn["Ci"] @ r
    b = sn["Ci1"] @ r
    return float(r @ cr - b * b / sn["e"])


# ----------------------------------------------------------------------------------------------------------------------
# history and models on a grid
# ----------------------------------------------------------------------------------------------------------------------
def history(g: Grid, om: float, orad: float, dX: np.ndarray):
    """E^2 = om (1+z)^3 + orad (1+z)^4 + (1 - om - orad) + dX; returns (dm_u, inv) dimensionless (D_M H0/c, H0/H) or None."""
    e2 = om * g.p3 + orad * g.p4 + (1.0 - om - orad) + dX
    if e2.min() <= 0.0:
        return None, float(e2.min())
    inv = 1.0 / np.sqrt(e2)
    f = g.p1 * inv
    dm_u = np.concatenate([[0.0], np.cumsum(0.5 * (f[:-1] + f[1:]) * g.DUK)])
    return (dm_u, inv), 0.0


def mu_sn(g: Grid, sn: dict, dm_u: np.ndarray) -> np.ndarray:
    return 5.0 * np.log10((1.0 + sn["zhel"]) * np.interp(sn["uhd"], g.U, dm_u))


def vec15(g: Grid, h: float, wm: float, dm_u: np.ndarray, inv: np.ndarray) -> np.ndarray:
    dm = np.array([dm_u[i] for i in g.IDX])
    dh = np.array([inv[i] for i in g.IDX])
    d = np.where(X.IS_DM, dm, np.where(X.IS_DH, dh, np.cbrt(X.Z13 * dm * dm * dh)))
    bao = C_KMS / (100.0 * h * float(E.RD_SPL(wm))) * d
    dm_star = C_KMS / (100.0 * h) * (dm_u[-1] + inv[-1] * (E.zstar_of(wm) - ZSTAR))
    th = 100.0 * float(E.RS_SPL(wm)) / dm_star / E.CAL
    return np.concatenate([bao, [th, wm]])


# ----------------------------------------------------------------------------------------------------------------------
# fits
# ----------------------------------------------------------------------------------------------------------------------
def fit_sn_only(g: Grid, dX: np.ndarray, sn: dict) -> tuple[float, float]:
    orad = OMR / H_SN**2

    def f(om):
        hist, m = history(g, om, orad, dX)
        if hist is None:
            return 1e8 * (1.0 - m)
        return chi2_sn(sn, mu_sn(g, sn, hist[0]))

    r = minimize_scalar(f, bounds=(0.05, 0.95), method="bounded", options={"xatol": 1e-6})
    return float(r.fun), float(r.x)


def resid_combo(p, g: Grid, dX: np.ndarray, sn: dict, with_cmb: bool):
    h, wm, M = p
    hist, m = history(g, wm / h**2, OMR / h**2, dX)
    nr = 13 + (2 if with_cmb else 0) + sn["n"]
    if hist is None:
        return np.full(nr, 1e3 * (1.0 - m))
    dm_u, inv = hist
    v = vec15(g, h, wm, dm_u, inv)
    r = X.Y15 - v
    out = [solve_triangular(X.L13, r[:13], lower=True, check_finite=False)]
    if with_cmb:
        out.append([r[13] / X.SIG_TH, r[14] / X.SIG_WM])
    rs = sn["mb"] - mu_sn(g, sn, dm_u) - M
    out.append(solve_triangular(sn["L"], rs, lower=True, check_finite=False))
    return np.concatenate(out)


def fit_combo(g: Grid, dX: np.ndarray, sn: dict, with_cmb: bool, x0=(0.682, 0.1407)) -> tuple[float, float, float, float]:
    # M start: analytic optimum at x0
    hist, _ = history(g, x0[1] / x0[0] ** 2, OMR / x0[0] ** 2, dX)
    M0 = 0.0
    if hist is not None:
        rr = sn["mb"] - mu_sn(g, sn, hist[0])
        M0 = float((sn["Ci1"] @ rr) / sn["e"])
    r = least_squares(resid_combo, np.array([x0[0], x0[1], M0]), bounds=([0.4, 0.05, -50.0], [1.0, 0.30, 50.0]), args=(g, dX, sn, with_cmb),
                      method="trf", x_scale=[0.01, 0.001, 0.01], ftol=1e-10, xtol=1e-10, max_nfev=200)
    return float(2.0 * r.cost), float(r.x[0]), float(r.x[1]), float(r.x[2])


def chi2_all(g: Grid, dX: np.ndarray) -> np.ndarray:
    """chi^2_min of the four combinations for one path."""
    return np.array([fit_sn_only(g, dX, SN["all"])[0], fit_combo(g, dX, SN["all"], False)[0], fit_combo(g, dX, SN["all"], True)[0],
                     fit_combo(g, dX, SN["z023"], True)[0]])


P0 = (0.681929893729943, 0.14071555057242932)  # LCDM best fit of the 15 rows (exacto.py)


def cell_seed(sigma: float, theta: float) -> list[int]:
    if sigma in X.SIGMAS:
        i = [(s, t) for t in X.THETAS for s in X.SIGMAS].index((sigma, theta))
        return [SEED, i]
    i = [(s, t) for t in THETAS for s in (0.01, 0.02, 0.03)].index((sigma, theta))
    return [SEED, 100 + i]


def _cell(args):
    sigma, theta, fine = args
    g = GF if fine else G0
    rng = np.random.default_rng(cell_seed(sigma, theta))
    Xp = ou_paths(g, sigma, theta, N_REAL, rng)
    dX = Xp - Xp[:, [0]]
    h0, wm0 = P0
    om0, or0 = wm0 / h0**2, OMR / h0**2
    e2min = (om0 * g.p3 + or0 * g.p4 + (1 - om0 - or0) + dX).min(axis=1)
    ok = e2min > 0.0
    rho_neg = ((1 - om0 - or0) + dX).min(axis=1) <= 0.0
    chi = np.full((N_REAL, 4), np.nan)
    for i in np.nonzero(ok)[0]:
        chi[i] = chi2_all(g, dX[i])
    return sigma, theta, fine, ok, rho_neg, chi


def _ceiling(args):
    theta, combo, b, chi_obs = args
    s_near = min(SIGMAS, key=lambda s: abs(np.log(s / b)))
    rng = np.random.default_rng(cell_seed(s_near, theta))
    Xp = ou_paths(G0, b, theta, N_REAL, rng)
    dX = Xp - Xp[:, [0]]
    h0, wm0 = P0
    om0, or0 = wm0 / h0**2, OMR / h0**2
    ok = (om0 * G0.p3 + or0 * G0.p4 + (1 - om0 - or0) + dX).min(axis=1) > 0
    rho_neg = ((1 - om0 - or0) + dX).min(axis=1) <= 0
    sn = SN["z023"] if combo.endswith("z023") else SN["all"]
    chi = np.full(N_REAL, np.inf)
    for i in np.nonzero(ok)[0]:
        chi[i] = fit_combo(G0, dX[i], sn, True)[0]
    surv = chi <= chi_obs
    return {"sigma_ceiling": b, "seed_cell_sigma": s_near, "n_survivors": int(surv.sum()), "frac_survive": float(surv.mean()),
            "frac_rho_le_0_among_survivors": float(rho_neg[surv].mean()) if surv.any() else None, "frac_H2_le_0": float(1 - ok.mean())}


def crossing(sig, p, target=0.05):
    sig, p = np.asarray(sig), np.asarray(p)
    if p[0] <= target:
        return "<" + f"{sig[0]:g}"
    for k in range(len(sig) - 1):
        if p[k] > target >= p[k + 1]:
            t = (p[k] - target) / (p[k] - p[k + 1])
            return float(np.exp(np.log(sig[k]) + t * (np.log(sig[k + 1]) - np.log(sig[k]))))
    return None


# ----------------------------------------------------------------------------------------------------------------------
# linear marginalised bounds, convention B, full response (cmb_fondo grid), SN rows added
# ----------------------------------------------------------------------------------------------------------------------
def linear_bounds() -> dict:
    R = CF.Response(lambda_absorbs_x0=True)  # 15 rows on the uniform grid
    u, z = R.u, R.z
    om, orad = CF.OM_FID, OMR / CF.H_FID**2
    E2 = om * (1 + z) ** 3 + orad * (1 + z) ** 4 + (1 - om - orad)
    Es = np.sqrt(E2)
    base = R.du * (1 + z)
    N = len(u)

    def wM_at(k):
        cw = np.zeros(N)
        cw[: k + 1] = base[: k + 1]
        cw[0] *= 0.5
        cw[k] *= 0.5
        dm = float(np.sum(cw / Es))
        return -cw / (2 * Es**3) / dm

    out = {}
    sn = SN["all"]
    # SN rows: linear interpolation of the row between bracketing nodes
    Wsn = np.zeros((sn["n"], N))
    for i, uu in enumerate(sn["uhd"]):
        k = int(np.searchsorted(u, uu)) - 1
        t = (uu - u[k]) / (u[k + 1] - u[k])
        row = (1 - t) * wM_at(k) + t * wM_at(k + 1)
        row[0] -= row.sum()  # Lambda absorbs X(0)
        Wsn[i] = MAG * row
    W = np.vstack([R.W[:15], Wsn])  # BAO (13, in delta ln D), theta_* (delta ln theta), omega_m (0), SN (delta m)
    # scaling to raw units: BAO rows x prediction, theta row x theta; SN rows already in mag
    fit15 = CF.Fit(CF.SETS["BAO+CMB"], R)
    scale = np.concatenate([fit15.scale[:15], np.ones(sn["n"])])
    Ws = W * scale[:, None]
    # data residual r0 and templates around the LCDM best fit of each combination
    y = np.concatenate([X.Y15, sn["mb"]])
    Cm = np.zeros((15 + sn["n"], 15 + sn["n"]))
    Cm[:15, :15] = X.CMEAS15 if hasattr(X, "CMEAS15") else CF.CMEAS15
    Cm[15:, 15:] = np.linalg.inv(sn["Ci"])
    rows_of = {"SN": list(range(15, 15 + sn["n"])), "BAO+SN": list(range(13)) + list(range(15, 15 + sn["n"])),
               "BAO+CMB+SN": list(range(15 + sn["n"]))}
    dX0 = np.zeros(G0.N)

    def model_full(h, wm, M):
        hist, _ = history(G0, wm / h**2, OMR / h**2, dX0)
        return np.concatenate([vec15(G0, h, wm, *hist), mu_sn(G0, sn, hist[0]) + M])

    for combo, rows in rows_of.items():
        rows = np.array(rows)
        if combo == "SN":
            c0, om_b = fit_sn_only(G0, dX0, sn)
            h_b, wm_b = H_SN, om_b * H_SN**2
            hist, _ = history(G0, om_b, OMR / H_SN**2, dX0)
            rr = sn["mb"] - mu_sn(G0, sn, hist[0])
            M_b = float((sn["Ci1"] @ rr) / sn["e"])
            pars = ["Omega_m", "M"]
        else:
            c0, h_b, wm_b, M_b = fit_combo(G0, dX0, sn, combo == "BAO+CMB+SN")
            pars = ["h", "omega_m", "M"]
        m0 = model_full(h_b, wm_b, M_b)
        T = []
        for pname in pars:
            eps = {"h": 1e-4, "omega_m": 1e-5, "Omega_m": 1e-4, "M": 1e-3}[pname]
            def mp(sign):
                if pname == "h":
                    return model_full(h_b + sign * eps, wm_b, M_b)
                if pname == "omega_m":
                    return model_full(h_b, wm_b + sign * eps, M_b)
                if pname == "Omega_m":
                    return model_full(h_b, (om_b + sign * eps) * H_SN**2, M_b)
                return model_full(h_b, wm_b, M_b + sign * eps)
            T.append((mp(+1) - mp(-1)) / (2 * eps))
        T = np.array(T).T[rows]
        r0 = (y - m0)[rows]
        C = Cm[np.ix_(rows, rows)]
        Wr = Ws[rows]
        rec = {"lcdm_chi2": c0, "params_best": {"h": h_b, "omega_m": wm_b, "M": M_b}}
        for theta in THETAS:
            Kw = R.Kw(theta, Wr)
            Mth = Wr @ Kw.T

            def ll_c(s, Mth=Mth):
                L = np.linalg.cholesky(C + s * s * Mth)
                solve = lambda b: np.linalg.solve(L.T, np.linalg.solve(L, b))
                cr, cT = solve(r0), solve(T)
                F = T.T @ cT
                b = T.T @ cr
                return float(-0.5 * (r0 @ cr - b @ np.linalg.solve(F, b) + 2 * np.sum(np.log(np.diag(L))) + np.linalg.slogdet(F)[1]))

            rec[f"theta_{theta:g}"] = {name: CF.V.bound(ll_c, thr)["sigma_95"] for name, thr in CF.THRS.items()}
        out[combo] = rec
        print("linear", combo, json.dumps(rec), flush=True)
    return out


def main() -> None:
    t0 = time.time()
    zero = np.zeros(G0.N)
    chk = {}
    c_sn, om_sn = fit_sn_only(G0, zero, SN["all"])
    chk["sn_only_lcdm"] = {"Omega_m": om_sn, "chi2_rel": c_sn, "expected_SN_md": 0.3326, "n_sn": SN["all"]["n"], "n_sn_z023": SN["z023"]["n"]}
    obs = chi2_all(G0, zero)
    chk["chi2_obs"] = dict(zip(COMBOS, obs.tolist()))
    c15, h15, wm15, _ = fit_combo(G0, zero, SN["all"], True)
    chk["bao_cmb_sn_lcdm"] = {"chi2": c15, "H0": 100 * h15, "omega_m": wm15}
    cb, hb, wmb, _ = fit_combo(G0, zero, SN["all"], False)
    chk["bao_sn_lcdm"] = {"chi2": cb, "H0": 100 * hb, "omega_m": wmb}
    chk["grid"] = {"nodes": G0.N, "nodes_fine": GF.N}
    print(json.dumps(chk, indent=1), flush=True)
    cells = [(s, t, False) for t in THETAS for s in SIGMAS] + [(s, 10.0, True) for s in SIGMAS]
    with Pool() as pool:
        res = pool.map(_cell, cells, chunksize=1)
    print(f"celdas: {time.time() - t0:.0f} s", flush=True)
    out = {"sigmas": SIGMAS, "thetas": THETAS, "n_real": N_REAL, "seed": SEED, "combos": COMBOS, "checks": chk, "cells": {}, "bounds": {}}
    raw = {}
    for sigma, theta, fine, ok, rho_neg, chi in res:
        key = f"theta_{theta:g}|sigma_{sigma:g}" + ("|fine" if fine else "")
        raw[key.replace("|", "__")] = np.column_stack([ok, rho_neg, chi])
        rec = {"frac_H2_le_0": float(1 - ok.mean()), "frac_rho_le_0_all": float(rho_neg.mean())}
        for j, combo in enumerate(COMBOS):
            c = chi[ok, j]
            rec[f"P_{combo}"] = float(np.mean(c <= obs[j]))
            rec[f"P_{combo}_err"] = float(np.sqrt(rec[f"P_{combo}"] * (1 - rec[f"P_{combo}"]) / max(ok.sum(), 1)))
            rec[f"chi2_{combo}_p5_p50"] = [float(np.percentile(c, 5)), float(np.percentile(c, 50))]
        out["cells"][key] = rec
    np.savez_compressed(HERE / "exacto_sn_chi2.npz", **raw, columns=np.array(["ok", "rho_neg"] + [f"chi2_{c}" for c in COMBOS]))
    for theta in THETAS:
        k = f"theta_{theta:g}"
        rec = {}
        for combo in COMBOS:
            p = [out["cells"][f"{k}|sigma_{s:g}"][f"P_{combo}"] for s in SIGMAS]
            rec[combo] = {"P_vs_sigma": p, "exact_5pct": crossing(SIGMAS, p)}
            if theta == 10.0:
                pf = [out["cells"][f"{k}|sigma_{s:g}|fine"][f"P_{combo}"] for s in SIGMAS]
                rec[combo]["P_vs_sigma_fine"] = pf
                rec[combo]["exact_5pct_fine"] = crossing(SIGMAS, pf)
                err = [out["cells"][f"{k}|sigma_{s:g}"][f"P_{combo}_err"] for s in SIGMAS]
                rec[combo]["max_abs_dP_fine_over_err"] = float(np.max(np.abs(np.array(p) - np.array(pf)) / np.maximum(np.sqrt(2) * np.array(err), 1e-9)))
        rec["frac_H2_le_0_vs_sigma"] = [out["cells"][f"{k}|sigma_{s:g}"]["frac_H2_le_0"] for s in SIGMAS]
        out["bounds"][k] = rec
    # rho_DE <= 0 among survivors at the new ceilings (BAO+CMB+SN and z023), in parallel
    jobs = []
    for theta in THETAS:
        for combo in ("BAO+CMB+SN", "BAO+CMB+SN_z023"):
            b = out["bounds"][f"theta_{theta:g}"][combo]["exact_5pct"]
            if isinstance(b, float):
                jobs.append((theta, combo, b, float(obs[COMBOS.index(combo)])))
    with Pool() as pool:
        cr = pool.map(_ceiling, jobs, chunksize=1)
    out["ceiling_rho"] = {f"theta_{j[0]:g}|{j[1]}": r for j, r in zip(jobs, cr)}
    print("ceiling rho:", json.dumps(out["ceiling_rho"], indent=1), flush=True)
    out["linear"] = linear_bounds()
    out["runtime_s"] = time.time() - t0
    (HERE / "exacto_sn.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    summary(out)


def summary(out: dict) -> None:
    f = lambda v: "sin cota (<= 0.7)" if v is None else (v if isinstance(v, str) else f"{v:.3g}")
    g = lambda v: "sin cota" if v is None else f"{v:.3g}"
    L = [f"chi2_obs: {json.dumps(out['checks']['chi2_obs'])}; SN solas LCDM Omega_m = {out['checks']['sn_only_lcdm']['Omega_m']:.4f} (SN.md 0.3326)", ""]
    for combo in COMBOS:
        L.append(f"{combo}: theta  exacta(5%)  [fina]  lineal marg -1.92 / -1.35   P por sigma {SIGMAS}")
        for theta in THETAS:
            r = out["bounds"][f"theta_{theta:g}"][combo]
            lin = out["linear"].get(combo.replace("_z023", ""), {}).get(f"theta_{theta:g}", {})
            fine = f"[{f(r['exact_5pct_fine'])}]" if "exact_5pct_fine" in r else ""
            L.append(f"   {theta:5g}  {f(r['exact_5pct']):>16s} {fine:>10s}  {g(lin.get('192')):>8s} / {g(lin.get('135')):>8s}   {[round(x, 3) for x in r['P_vs_sigma']]}")
        L.append("")
    L.append("rho_DE <= 0 entre supervivientes en los techos nuevos: " + json.dumps(out["ceiling_rho"]))
    (HERE / "exacto_sn.txt").write_text("\n".join(L) + "\n", encoding="utf-8")
    print("\n".join(L))


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--txt-only":
        summary(json.loads((HERE / "exacto_sn.json").read_text()))
    else:
        main()
