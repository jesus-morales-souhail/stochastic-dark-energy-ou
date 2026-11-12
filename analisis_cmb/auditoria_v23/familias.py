#!/usr/bin/env python3
"""
Exploration only (FAMILIAS.md): is the allowed change of the dark energy between today and z ~ 0.5 the same for different
families of fluctuations, even if their sigma_X bounds differ? Nothing in results/ or papers changes.

Set-up of exacto_sn.py unchanged (convention B, exact history, BAO + theta_* + omega_m + Pantheon+, same data, covariances
and treatment of h, omega_m and M; same grid in u = ln(1+z) with the DESI z as exact nodes; discard if H^2 <= 0 at the LCDM
best fit of the 15 rows; P(chi^2_min <= chi^2_obs) = 5 % interpolated in ln sigma_X). Only the process X(u), stationary,
zero mean, variance sigma_X^2, changes:
- OU:        K = exp(-|du|/l) (theta = 1/l); exacto_sn.ou_paths with exacto_sn's seeds where the cell exists there.
- Matern 3/2: K = (1 + sqrt3 r/l) exp(-sqrt3 r/l); exact two-state Markov process (X, X') with the transition exp(F du),
             F = [[0, 1], [-lam^2, -2 lam]], lam = sqrt3/l, innovation Q = P_inf - Phi P_inf Phi', P_inf = diag(1, lam^2) sigma^2
             (exact on the non-uniform grid: Phi(d1) Phi(d2) = Phi(d1 + d2)).
- Gaussian:  K = exp(-r^2/(2 l^2)); Cholesky of the correlation matrix with the smallest jitter in {1e-12, ..., 1e-6} that
             succeeds (reported; the run must stop if > 1e-6).
- Steps:     independent N(0, sigma^2) values on intervals of width l with a uniform random offset per path; K = (1 - |du|/l)_+;
             the path is a function of u (offset and values drawn first), so the grid only enters through the integration.
Cells: l in {0.1, 0.3, 1, 3, 10} x sigma_X in {0.02, ..., 1.0}; 2000 paths per cell, seed per cell and family.
Per family and l: (a) 5 % bound with bootstrap, and the drop of 1.92 / 1.35 of the analytic linear marginal (convention B,
full response A of marginal_guiada.build_A, kernel of the family) from the GLOBAL maximum over the extended sigma grid of
maximos_globales.py, with the class of the maximum; (b) at each bound, Delta rho(z) = X(z) - X(0) at z in Z_EVAL: prior rms
(analytic, sigma sqrt(2 (1 - K(u_z))), checked against Monte Carlo) and p95 of |Delta rho| among survivors (paths regenerated
at the bound with the seed of the nearest cell); (c) z_* where the spread of the prior rms across families and l is minimal;
(d) fraction of survivors with rho_DE <= 0; (e) guided exact marginal (marginal_guiada.py, N = 5000) at l = 1 and sigma_X at
0.6-1.5 x the linear bound, four families (for Steps, the Gaussian process with the triangular kernel: same covariance as the
linear marginal); (f) Steps on the half-step grid at l = 0.3 and 1 (same offsets and values); extra: the Gaussianised Steps
(Cholesky of the triangular kernel) through the 5 % rule at l = 0.3 and 1, to see whether the non-Gaussianity matters.
Checks: OU cells shared with exacto_sn.py must reproduce its chi^2 exactly; sample variance and lag correlations of each
sampler against K; P at sigma_X = 0.02 in 0.40-0.50 for every family; effective scales (e^-1 separation and integral scale).

    ../venv/bin/python familias.py            # writes familias.{json,txt} and familias_chi2.npz
    ../venv/bin/python familias.py --txt-only
"""

from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.linalg import cho_factor, cho_solve, solve_triangular
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import exacto as X  # noqa: E402
import exacto_sn as S  # noqa: E402
import marginal_guiada as MG  # noqa: E402
import cmb_fondo as CF  # noqa: E402
import precision_sn as PS  # noqa: E402

ELLS = (0.1, 0.3, 1.0, 3.0, 10.0)
SIGMAS = (0.02, 0.03, 0.05, 0.07, 0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 1.0)
FAMILIES = ("OU", "Matern32", "Gauss", "Steps")
EXTRA = {"StepsGauss": (0.3, 1.0)}  # Gaussianised steps (triangular kernel, Cholesky), 5 % rule only
HALF_ELLS = (0.3, 1.0)  # Steps on the half-step grid
N_REAL = S.N_REAL
N_BOOT = 1000
SEED_F = 20261116
Z_EVAL = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.33)
U_EVAL = np.log1p(np.array(Z_EVAL))
RULES = ("5pct", "192", "135")
GUIDED_FACTORS = (0.6, 0.8, 1.0, 1.25, 1.5)
JITTERS = (1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6)
SN = S.SN["all"]
P0 = S.P0
G0 = S.G0
OMR = S.OMR
GRID_LIN = np.concatenate([CF.V.SIG, np.geomspace(CF.V.SIG[-1], 1.0, 40)[1:]])
THRS = {"192": 1.92, "135": 1.35}


class GridU(S.Grid):
    """exacto_sn.Grid from an explicit base u array (plus the DESI nodes)."""

    def __init__(self, u_base):
        self.U = np.unique(np.concatenate([u_base, np.log1p(np.unique(X.Z13))]))
        self.Z = np.expm1(self.U)
        self.N = len(self.U)
        self.DUK = np.diff(self.U)
        self.IDX = [int(np.argmin(np.abs(self.U - np.log1p(z)))) for z in X.Z13]
        assert max(abs(self.U[i] - np.log1p(z)) for i, z in zip(self.IDX, X.Z13)) < 1e-12
        self.p3, self.p4, self.p1 = (1 + self.Z) ** 3, (1 + self.Z) ** 4, 1 + self.Z


GH = GridU(np.linspace(0.0, X._u[-1], 2 * (len(X._u) - 1) + 1))  # half step everywhere
assert abs(np.median(np.diff(GH.U)) / (np.median(np.diff(X._u)) / 2) - 1) < 1e-6


# ----------------------------------------------------------------------------------------------------------------------
# kernels (unit variance) and samplers
# ----------------------------------------------------------------------------------------------------------------------
def kfun(fam: str, r, l: float):
    r = np.abs(np.asarray(r, dtype=float))
    if fam == "OU":
        return np.exp(-r / l)
    if fam == "Matern32":
        a = np.sqrt(3.0) * r / l
        return (1.0 + a) * np.exp(-a)
    if fam == "Gauss":
        return np.exp(-(r * r) / (2.0 * l * l))
    if fam in ("Steps", "StepsGauss"):
        return np.maximum(1.0 - r / l, 0.0)
    raise ValueError(fam)


def kernel_matrix(fam: str, g, l: float) -> np.ndarray:
    return kfun(fam, g.U[:, None] - g.U[None, :], l)


def chol_corr(fam: str, g, l: float):
    """Cholesky of the correlation matrix with the smallest jitter that works; returns (L, jitter)."""
    K = kernel_matrix(fam, g, l)
    for jit in JITTERS:
        try:
            L = np.linalg.cholesky(K + jit * np.eye(g.N))
            return L, jit
        except np.linalg.LinAlgError:
            continue
    raise RuntimeError(f"Cholesky of {fam} l={l} failed even with jitter 1e-6")


def paths_matern(g, sigma: float, l: float, n: int, rng) -> np.ndarray:
    lam = np.sqrt(3.0) / l
    d = g.DUK
    e = np.exp(-lam * d)
    # Phi_k (2x2 per step) and Q_k = P - Phi P Phi', P = diag(1, lam^2); 2x2 Cholesky of Q
    P = np.diag([1.0, lam * lam])
    Phi = np.empty((len(d), 2, 2))
    Phi[:, 0, 0], Phi[:, 0, 1] = e * (1 + lam * d), e * d
    Phi[:, 1, 0], Phi[:, 1, 1] = -e * lam * lam * d, e * (1 - lam * d)
    Q = P[None] - np.einsum("kij,jl,kml->kim", Phi, P, Phi)
    a, b, c = np.maximum(Q[:, 0, 0], 0.0), Q[:, 0, 1], np.maximum(Q[:, 1, 1], 0.0)
    Lq = np.zeros_like(Q)
    Lq[:, 0, 0] = np.sqrt(a)
    Lq[:, 1, 0] = np.where(a > 0, b / np.sqrt(np.where(a > 0, a, 1.0)), 0.0)
    Lq[:, 1, 1] = np.sqrt(np.maximum(c - Lq[:, 1, 0] ** 2, 0.0))
    Y = np.empty((n, 2))
    Y[:, 0] = rng.standard_normal(n)
    Y[:, 1] = lam * rng.standard_normal(n)
    out = np.empty((n, g.N))
    out[:, 0] = Y[:, 0]
    xi = rng.standard_normal((n, g.N - 1, 2))
    for k in range(g.N - 1):
        Y = Y @ Phi[k].T + xi[:, k, :] @ Lq[k].T
        out[:, k + 1] = Y[:, 0]
    return sigma * out


def paths_chol(fam: str, g, sigma: float, l: float, n: int, rng):
    L, jit = chol_corr(fam, g, l)
    return sigma * (rng.standard_normal((n, g.N)) @ L.T), jit


def steps_draw(l: float, n: int, rng, umax: float):
    """Offsets and step values first (grid-independent); returns (phi, V)."""
    phi = rng.uniform(0.0, l, n)
    ncell = int(np.ceil(umax / l)) + 2
    V = rng.standard_normal((n, ncell))
    return phi, V


def steps_eval(g, sigma: float, l: float, phi, V) -> np.ndarray:
    idx = np.floor((g.U[None, :] + phi[:, None]) / l).astype(int)
    return sigma * V[np.arange(len(phi))[:, None], idx]


def sample(fam: str, g, sigma: float, l: float, n: int, rng):
    """Paths (n, g.N) and a dict of extras (jitter)."""
    if fam == "OU":
        return S.ou_paths(g, sigma, 1.0 / l, n, rng), {}
    if fam == "Matern32":
        return paths_matern(g, sigma, l, n, rng), {}
    if fam in ("Gauss", "StepsGauss"):
        Xp, jit = paths_chol(fam, g, sigma, l, n, rng)
        return Xp, {"jitter_rel": jit}
    if fam == "Steps":
        phi, V = steps_draw(l, n, rng, X._u[-1])
        return steps_eval(g, sigma, l, phi, V), {}
    raise ValueError(fam)


def cell_seed(fam: str, sigma: float, l: float):
    """exacto_sn's seed for the OU cells that exist there; otherwise [SEED_F, family, l, sigma] indices."""
    if fam == "OU":
        theta = 1.0 / l
        for t in S.THETAS:
            if abs(theta - t) < 1e-9 and (sigma in X.SIGMAS or sigma in (0.01, 0.02, 0.03)):
                return S.cell_seed(sigma, t)
    fi = (FAMILIES + tuple(EXTRA)).index(fam)
    return [SEED_F, fi, ELLS.index(l), SIGMAS.index(sigma)]


# ----------------------------------------------------------------------------------------------------------------------
# cells
# ----------------------------------------------------------------------------------------------------------------------
def drho_at(g, dX: np.ndarray) -> np.ndarray:
    """Delta rho = X(z) - X(0) at Z_EVAL by linear interpolation in u (n, len(Z_EVAL))."""
    k = np.searchsorted(g.U, U_EVAL) - 1
    t = (U_EVAL - g.U[k]) / (g.U[k + 1] - g.U[k])
    return (1 - t)[None, :] * dX[:, k] + t[None, :] * dX[:, k + 1]


def run_paths(g, dX: np.ndarray):
    h0, wm0 = P0
    om0, or0 = wm0 / h0**2, OMR / h0**2
    ok = (om0 * g.p3 + or0 * g.p4 + (1 - om0 - or0) + dX).min(axis=1) > 0.0
    rho_neg = ((1 - om0 - or0) + dX).min(axis=1) <= 0.0
    chi = np.full(len(dX), np.nan)
    for i in np.nonzero(ok)[0]:
        chi[i] = S.fit_combo(g, dX[i], SN, True)[0]
    return ok, rho_neg, chi


def _cell(args):
    fam, l, sigma, half = args
    g = GH if half else G0
    t0 = time.time()
    rng = np.random.default_rng(cell_seed(fam, sigma, l))
    Xp, extra = sample(fam, g, sigma, l, N_REAL, rng)
    dX = Xp - Xp[:, [0]]
    ok, rho_neg, chi = run_paths(g, dX)
    return {"fam": fam, "l": l, "sigma": sigma, "half": half, "ok": ok, "rho_neg": rho_neg, "chi": chi,
            "drho": drho_at(g, dX).astype(np.float32), "extra": extra, "seconds": time.time() - t0}


def _regen(args):
    """Paths regenerated at a bound sigma_b (seed of the nearest cell): survivors' statistics."""
    fam, l, rule, sigma_b, chi_obs = args
    s_near = min(SIGMAS, key=lambda s: abs(np.log(s / sigma_b)))
    rng = np.random.default_rng(cell_seed(fam, s_near, l))
    Xp, _ = sample(fam, G0, sigma_b, l, N_REAL, rng)
    dX = Xp - Xp[:, [0]]
    dr = drho_at(G0, dX)
    ok, rho_neg, chi = run_paths(G0, dX)
    surv = ok & (chi <= chi_obs)
    rec = {"fam": fam, "l": l, "rule": rule, "sigma_b": sigma_b, "seed_cell_sigma": s_near, "n_ok": int(ok.sum()), "n_surv": int(surv.sum()),
           "frac_survive_among_ok": float(surv.sum() / max(ok.sum(), 1)), "frac_H2_le_0": float(1 - ok.mean()),
           "frac_rho_le_0_among_survivors": float(rho_neg[surv].mean()) if surv.any() else None,
           "frac_rho_le_0_all": float(rho_neg.mean()),
           "rms_drho_prior_mc": np.sqrt(np.mean(dr**2, axis=0)).tolist(),
           "rms_drho_prior_analytic": (sigma_b * np.sqrt(2.0 * (1.0 - kfun(fam, U_EVAL, l)))).tolist(),
           "p95_abs_drho_survivors": np.percentile(np.abs(dr[surv]), 95, axis=0).tolist() if surv.sum() >= 20 else None,
           "rms_drho_survivors": np.sqrt(np.mean(dr[surv] ** 2, axis=0)).tolist() if surv.any() else None}
    return rec


def crossing_boot(sig, P_boot):
    """P_boot (B, n_sigma) -> bounds (B,) (nan = none below 5 %, -1 = below at the first node)."""
    out = np.full(P_boot.shape[0], np.nan)
    ls = np.log(np.asarray(sig))
    for b in range(P_boot.shape[0]):
        c = S.crossing(sig, P_boot[b])
        out[b] = -1.0 if isinstance(c, str) else (np.nan if c is None else c)
    return out


# ----------------------------------------------------------------------------------------------------------------------
# linear marginal (analytic, convention B, full response), kernel of the family
# ----------------------------------------------------------------------------------------------------------------------
def _linear(args):
    fam, l = args
    t0 = time.time()
    A, Bm = LIN["A"], LIN["Bm"]
    K = kernel_matrix(fam, G0, l)
    AK = A @ K
    Msig = AK @ A.T
    Msig = 0.5 * (Msig + Msig.T)
    r = MG.D - MG.B0

    def lnz(s):
        Sf = cho_factor(MG.C + s * s * Msig, lower=True)
        return MG.lnZ_lin(Sf, Bm, r)

    z0 = lnz(0.0)
    prof = np.array([lnz(float(s)) - z0 for s in GRID_LIN])
    vals = np.concatenate([[0.0], prof])
    nodes = np.concatenate([[0.0], GRID_LIN])
    k0 = int(np.argmax(vals))
    rec = {"fam": fam, "l": l, "grid": GRID_LIN.tolist(), "lnZ_rel": prof.tolist(), "sigma_at_max": float(nodes[k0]), "lnZ_max_rel": float(vals[k0]),
           "max_class": "sigma->0" if k0 == 0 else ("upper edge (1.0)" if k0 == len(nodes) - 1 else "interior")}
    for name, drop in THRS.items():
        thr = vals[k0] - drop
        below = np.nonzero((vals < thr) & (np.arange(len(vals)) > k0))[0]
        if len(below) == 0:
            rec[name] = None
        else:
            k = int(below[0])
            lo = nodes[k - 1]
            if lo == 0.0:
                rec[name] = -1.0
            else:
                rec[name] = float(brentq(lambda s: lnz(s) - z0 - thr, lo, nodes[k], xtol=1e-10 * lo))
    rec["seconds"] = time.time() - t0
    return rec


LIN = {}


def init_linear():
    LIN["A"], LIN["Bm"] = MG.build_A(), MG.build_B()


# ----------------------------------------------------------------------------------------------------------------------
# guided exact marginal (marginal_guiada.cell with the family's prior sampler)
# ----------------------------------------------------------------------------------------------------------------------
def guided_cell(fam, l, sigma, K, AK, Msig, A, Bm, rng, Lcorr=None):
    N_SAMP, NROWS = MG.N_SAMP, MG.NROWS
    r = MG.D - MG.B0
    Sf = cho_factor(MG.C + sigma**2 * Msig, lower=True)
    lz = MG.lnZ_lin(Sf, Bm, r)
    SB = cho_solve(Sf, Bm)
    Sp = np.linalg.inv(Bm.T @ SB)
    dp_hat = Sp @ (Bm.T @ cho_solve(Sf, r))
    dp = dp_hat[None, :] + rng.standard_normal((N_SAMP, 3)) @ np.linalg.cholesky(Sp).T
    if sigma > 0:
        if fam in ("Gauss", "StepsGauss"):
            x0 = sigma * (rng.standard_normal((N_SAMP, G0.N)) @ Lcorr.T)
        else:
            x0, _ = sample(fam, G0, sigma, l, N_SAMP, rng)
    else:
        x0 = np.zeros((N_SAMP, G0.N))
    n0 = rng.standard_normal((N_SAMP, NROWS)) @ MG.LC.T
    R = r[None, :] - dp @ Bm.T - x0 @ A.T - n0
    x = x0 + sigma**2 * cho_solve(Sf, R.T).T @ AK if sigma > 0 else x0
    m_lin = MG.B0[None, :] + dp @ Bm.T + x @ A.T
    chi_lin = MG.chi2_of(m_lin)
    p = MG.P0[None, :] + dp
    m_ex, bad = MG.exact_model(p[:, 0], p[:, 1], p[:, 2], x - x[:, [0]])
    chi_ex = np.where(bad, np.inf, MG.chi2_of(np.where(bad[:, None], 0.0, m_ex)))
    lw = -0.5 * (chi_ex - chi_lin)
    w = np.where(bad, 0.0, np.exp(np.minimum(lw, 700)))
    idx = rng.integers(0, N_SAMP, size=(MG.N_BOOT, N_SAMP))
    lnw_b = np.log(np.maximum(w[idx].mean(axis=1), 1e-300))
    return {"sigma": sigma, "lnZ_lin": lz, "ln_mean_w": float(np.log(max(w.mean(), 1e-300))), "ln_mean_w_boot": lnw_b.tolist(),
            "ess": float(w.sum() ** 2 / np.sum(w**2)) if w.sum() > 0 else 0.0, "w_max_frac": float(w.max() / w.sum()) if w.sum() > 0 else 1.0,
            "frac_H2_le_0": float(bad.mean()), "chi_lin_mean": float(chi_lin.mean()), "chi_ex_mean_finite": float(np.mean(chi_ex[np.isfinite(chi_ex)])),
            "lw_p16_p50_p84": [float(v) for v in np.percentile(lw[np.isfinite(lw)], [16, 50, 84])]}


def _guided(args):
    fam, l, b_lin, s_max = args
    t0 = time.time()
    A, Bm = LIN["A"], LIN["Bm"]
    K = kernel_matrix(fam, G0, l)
    AK = A @ K
    Msig = AK @ A.T
    Msig = 0.5 * (Msig + Msig.T)
    Lcorr = chol_corr(fam, G0, l)[0] if fam in ("Gauss", "StepsGauss") else None
    rng = np.random.default_rng([SEED_F, 999, FAMILIES.index(fam)])
    # nodes: around the linear maximum (so the drop is measured from the maximum) and 0.6-1.5 x the linear bound
    sigs = sorted({0.0} | {f * s_max for f in (0.5, 1.0, 1.5) if 0 < f * s_max < 0.6 * b_lin} | {f * b_lin for f in GUIDED_FACTORS})
    cells = [guided_cell(fam, l, s, K, AK, Msig, A, Bm, rng, Lcorr) for s in sigs]
    lz_lin = np.array([c["lnZ_lin"] for c in cells])
    lnw = np.array([c["ln_mean_w"] for c in cells])
    lz_lin_rel = lz_lin - lz_lin[0]
    lz_ex_rel = (lz_lin + lnw) - (lz_lin[0] + lnw[0])
    sig_pos = np.array(sigs[1:])
    rec = {"fam": fam, "l": l, "sigma_linear_192": b_lin, "sigmas": sigs, "cells": cells, "lnZ_lin_rel": lz_lin_rel.tolist(), "lnZ_exact_rel": lz_ex_rel.tolist(),
           "correction": (lnw - lnw[0]).tolist()}
    for name, drop in THRS.items():
        rec[f"lin_{name}_from_nodes"] = PS.lnl_bound(sig_pos, lz_lin_rel[1:], drop)
        rec[f"exact_{name}"] = PS.lnl_bound(sig_pos, lz_ex_rel[1:], drop)
        bb = []
        for b in range(MG.N_BOOT):
            lnw_b = np.array([c["ln_mean_w_boot"][b] for c in cells])
            bb.append(PS.lnl_bound(sig_pos, ((lz_lin + lnw_b) - (lz_lin[0] + lnw_b[0]))[1:], drop))
        bb = np.array(bb)
        fin = np.isfinite(bb) & (bb > 0)
        rec[f"exact_{name}_boot_p16_p84"] = [float(np.percentile(bb[fin], 16)), float(np.percentile(bb[fin], 84))] if fin.sum() > 10 else None
        rec[f"exact_{name}_boot_frac_no_bound"] = float(np.mean(~np.isfinite(bb)))
    ess = np.array([c["ess"] for c in cells[1:]])
    trusted = ess >= 0.01 * MG.N_SAMP
    rec["ess_vs_sigma"] = ess.tolist()
    rec["nodes_trusted_ess_ge_1pct"] = sig_pos[trusted].tolist()
    for name, drop in THRS.items():
        rec[f"exact_{name}_trusted"] = PS.lnl_bound(sig_pos[trusted], lz_ex_rel[1:][trusted], drop) if trusted.sum() >= 2 else None
    kb = int(np.argmin(np.abs(np.array(sigs) - b_lin)))
    rec["at_linear_bound"] = {"ess": cells[kb]["ess"], "ess_frac_N": cells[kb]["ess"] / MG.N_SAMP, "correction": float(lnw[kb] - lnw[0]),
                              "w_max_frac": cells[kb]["w_max_frac"], "frac_H2_le_0": cells[kb]["frac_H2_le_0"]}
    rec["ess_min"] = float(min(c["ess"] for c in cells[1:]))
    rec["seconds"] = time.time() - t0
    return rec


# ----------------------------------------------------------------------------------------------------------------------
# checks: samplers against K, effective scales
# ----------------------------------------------------------------------------------------------------------------------
def sampler_check(fam: str, l: float, n: int = 4000) -> dict:
    rng = np.random.default_rng([SEED_F, 7, FAMILIES.index(fam) if fam in FAMILIES else 9, ELLS.index(l)])
    Xp, extra = sample(fam, G0, 1.0, l, n, rng)
    var = float(np.mean(Xp**2))
    lags = (0.05, 0.1, 0.3, 1.0, 3.0)
    res = {}
    for r in lags:
        if r >= G0.U[-1]:
            continue
        k = int(np.searchsorted(G0.U, r))
        # correlation at lag r averaged over starting nodes (every 50th)
        i0 = np.arange(0, G0.N - k, 50)
        c = float(np.mean(Xp[:, i0] * Xp[:, i0 + k]) / var)  # correlation (normalised by the sample variance)
        res[f"lag_{r:g}"] = {"sample": c, "K": float(kfun(fam, G0.U[i0[0] + k] - G0.U[i0[0]], l)), "mc_err": float(1.0 / np.sqrt(n))}
    worst = max(abs(v["sample"] - v["K"]) / v["mc_err"] for v in res.values())
    return {"var_sample": var, "var_err": float(np.sqrt(2.0 / n)), "lags": res, "worst_lag_dev_over_err": worst, **extra}


def effective_scales(fam: str, l: float) -> dict:
    f = lambda r: float(kfun(fam, r, l))
    r_e = brentq(lambda r: f(r) - np.exp(-1.0), 1e-9 * l, 50 * l)
    upper = l if fam in ("Steps", "StepsGauss") else 60 * l
    integ = quad(f, 0.0, upper, limit=200)[0]
    return {"r_e_over_l": r_e / l, "integral_over_l": integ / l}


EXPECTED_SCALES = {"OU": (1.0, 1.0), "Matern32": (1.24, 1.15), "Gauss": (1.41, 1.25), "Steps": (0.63, 0.5)}


# ----------------------------------------------------------------------------------------------------------------------
def main() -> None:
    t0 = time.time()
    out = {"ells": ELLS, "sigmas": SIGMAS, "families": FAMILIES, "extra": EXTRA, "n_real": N_REAL, "seed": SEED_F, "z_eval": Z_EVAL, "checks": {}}
    zero = np.zeros(G0.N)
    chi_obs = float(S.fit_combo(G0, zero, SN, True)[0])
    chi_obs_half = float(S.fit_combo(GH, np.zeros(GH.N), SN, True)[0])
    out["chi2_obs"] = chi_obs
    out["chi2_obs_half_grid"] = chi_obs_half
    out["checks"]["grid"] = {"nodes": G0.N, "nodes_half": GH.N, "chi2_obs_exacto_sn": 1406.3271474297196}
    # effective scales and sampler checks
    out["checks"]["scales"] = {fam: {f"l_{l:g}": effective_scales(fam, l) for l in ELLS} for fam in FAMILIES}
    out["checks"]["scales_expected"] = EXPECTED_SCALES
    out["checks"]["samplers"] = {fam: {f"l_{l:g}": sampler_check(fam, l) for l in ELLS} for fam in FAMILIES + tuple(EXTRA)}
    print(json.dumps({"chi2_obs": chi_obs, "chi2_obs_half": chi_obs_half, "scales": out["checks"]["scales"]}, indent=1), flush=True)
    for fam, d in out["checks"]["samplers"].items():
        for lk, v in d.items():
            print(f"sampler {fam:10s} {lk:6s} var={v['var_sample']:.3f} worst lag dev/err={v['worst_lag_dev_over_err']:.1f} {v.get('jitter_rel', '')}", flush=True)
    # cells (longest first)
    cells = [("Gauss", l, s, False) for l in ELLS for s in SIGMAS] + [("Steps", l, s, True) for l in HALF_ELLS for s in SIGMAS]
    cells += [(fam, l, s, False) for fam in ("OU", "Matern32", "Steps") for l in ELLS for s in SIGMAS]
    cells += [(fam, l, s, False) for fam, ls in EXTRA.items() for l in ls for s in SIGMAS]
    res = []
    with Pool() as pool:
        for i, r in enumerate(pool.imap_unordered(_cell, cells, chunksize=1)):
            res.append(r)
            if i % 20 == 0:
                print(f"  celda {i + 1}/{len(cells)} ({time.time() - t0:.0f} s): {r['fam']} l={r['l']} s={r['sigma']} half={r['half']} {r['seconds']:.0f} s", flush=True)
    print(f"celdas: {time.time() - t0:.0f} s", flush=True)
    raw, out["cells"] = {}, {}
    for r in res:
        key = f"{r['fam']}|l_{r['l']:g}|sigma_{r['sigma']:g}" + ("|half" if r["half"] else "")
        raw[key.replace("|", "__")] = np.column_stack([r["ok"], r["rho_neg"], r["chi"], r["drho"]])
        ok, chi = r["ok"], r["chi"]
        co = chi_obs_half if r["half"] else chi_obs
        surv = ok & (chi <= co)
        out["cells"][key] = {"P": float(np.mean(chi[ok] <= co)), "P_err": float(np.sqrt(np.mean(chi[ok] <= co) * (1 - np.mean(chi[ok] <= co)) / max(ok.sum(), 1))),
                             "n_ok": int(ok.sum()), "frac_H2_le_0": float(1 - ok.mean()), "frac_rho_le_0_all": float(r["rho_neg"].mean()),
                             "frac_rho_le_0_among_survivors": float(r["rho_neg"][surv].mean()) if surv.any() else None,
                             "n_chi_equal_obs": int(np.sum(np.abs(chi[ok] - co) < 1e-6)),
                             "chi2_p5_p50": [float(np.percentile(chi[ok], 5)), float(np.percentile(chi[ok], 50))], **r["extra"], "seconds": r["seconds"]}
    np.savez_compressed(HERE / "familias_chi2.npz", **raw, columns=np.array(["ok", "rho_neg", "chi2_BAO+CMB+SN"] + [f"drho_z{z:g}" for z in Z_EVAL]))
    # reproduction of exacto_sn (OU cells with the same seeds)
    ref = np.load(HERE / "exacto_sn_chi2.npz")
    rep = {}
    for l in (0.1, 1.0, 10.0):
        theta = 1.0 / l
        for s in SIGMAS:
            kk = f"theta_{theta:g}__sigma_{s:g}"
            if kk in ref.files:
                a = raw[f"OU__l_{l:g}__sigma_{s:g}"][:, 2]
                b = ref[kk][:, 4]  # chi2_BAO+CMB+SN column
                m = np.isfinite(a) & np.isfinite(b)
                rep[f"l_{l:g}|sigma_{s:g}"] = {"max_abs_dchi2": float(np.max(np.abs(a[m] - b[m]))), "same_ok": bool(np.array_equal(np.isfinite(a), np.isfinite(b)))}
    out["checks"]["reproduce_exacto_sn"] = rep
    print("reproduce exacto_sn: max |dchi2| =", max(v["max_abs_dchi2"] for v in rep.values()), "all same ok:", all(v["same_ok"] for v in rep.values()), flush=True)
    # 5 % bounds with bootstrap
    rng_b = np.random.default_rng([SEED_F, 1])
    out["bounds"] = {}
    for fam in FAMILIES + tuple(EXTRA):
        for l in (ELLS if fam in FAMILIES else EXTRA[fam]):
            for half in ((False, True) if fam == "Steps" and l in HALF_ELLS else (False,)):
                keys = [f"{fam}|l_{l:g}|sigma_{s:g}" + ("|half" if half else "") for s in SIGMAS]
                P = np.array([out["cells"][k]["P"] for k in keys])
                b5 = S.crossing(SIGMAS, P)
                Pb = np.empty((N_BOOT, len(SIGMAS)))
                for j, k in enumerate(keys):
                    arr = raw[k.replace("|", "__")]
                    ok = arr[:, 0] > 0
                    ind = (arr[ok, 2] <= (chi_obs_half if half else chi_obs)).astype(float)
                    idx = rng_b.integers(0, len(ind), size=(N_BOOT, len(ind)))
                    Pb[:, j] = ind[idx].mean(axis=1)
                bb = crossing_boot(SIGMAS, Pb)
                fin = np.isfinite(bb) & (bb > 0)
                rec = {"P_vs_sigma": P.tolist(), "exact_5pct": b5, "boot_p16_p84": [float(np.percentile(bb[fin], 16)), float(np.percentile(bb[fin], 84))] if fin.sum() > 10 else None,
                       "boot_frac_no_bound": float(np.mean(np.isnan(bb))), "frac_H2_le_0_vs_sigma": [out["cells"][k]["frac_H2_le_0"] for k in keys],
                       "n_chi_equal_obs_vs_sigma": [out["cells"][k]["n_chi_equal_obs"] for k in keys]}
                out["bounds"][f"{fam}|l_{l:g}" + ("|half" if half else "")] = rec
    # small-sigma check
    out["checks"]["P_at_sigma_0.02"] = {f"{fam}|l_{l:g}": out["cells"][f"{fam}|l_{l:g}|sigma_0.02"]["P"] for fam in FAMILIES for l in ELLS}
    # linear marginal bounds (global maximum), in parallel over (fam, l)
    init_linear()
    pairs = [(fam, l) for fam in FAMILIES for l in ELLS]
    with Pool(min(16, len(pairs))) as pool:
        lin = pool.map(_linear, pairs, chunksize=1)
    out["linear"] = {f"{r['fam']}|l_{r['l']:g}": r for r in lin}
    print(f"lineal: {time.time() - t0:.0f} s", flush=True)
    for r in lin:
        print(f"  lineal {r['fam']:9s} l={r['l']:<4g} max@{r['sigma_at_max']:g} ({r['max_class']}) 192={PS.fmt(r['192'])} 135={PS.fmt(r['135'])}", flush=True)
    # regeneration at the bounds (survivors' Delta rho, rho <= 0)
    jobs = []
    for fam in FAMILIES:
        for l in ELLS:
            b5 = out["bounds"][f"{fam}|l_{l:g}"]["exact_5pct"]
            for rule in RULES:
                b = b5 if rule == "5pct" else out["linear"][f"{fam}|l_{l:g}"][rule]
                if isinstance(b, float) and b > 0:
                    jobs.append((fam, l, rule, b, chi_obs))
    with Pool() as pool:
        rg = pool.map(_regen, jobs, chunksize=1)
    out["at_bound"] = {f"{r['fam']}|l_{r['l']:g}|{r['rule']}": r for r in rg}
    print(f"regeneracion en las cotas: {time.time() - t0:.0f} s", flush=True)
    # z_* and the universality table
    out["universality"] = universality(out)
    # guided exact at l = 1
    run_guided(out)
    out["runtime_s"] = time.time() - t0
    (HERE / "familias.json").write_text(json.dumps(out, indent=1, default=_json_default), encoding="utf-8")
    summary(out)


def run_guided(out: dict) -> None:
    gj = [(fam, 1.0, out["linear"][f"{fam}|l_1"]["192"], out["linear"][f"{fam}|l_1"]["sigma_at_max"]) for fam in FAMILIES
          if isinstance(out["linear"][f"{fam}|l_1"]["192"], float) and out["linear"][f"{fam}|l_1"]["192"] > 0]
    with Pool(len(gj)) as pool:
        gr = pool.map(_guided, gj, chunksize=1)
    out["guided_l1"] = {r["fam"]: r for r in gr}
    for r in gr:
        print(f"  guiada {r['fam']:9s} lin192={PS.fmt(r['sigma_linear_192'])} exact192={PS.fmt(r['exact_192'])} {r['exact_192_boot_p16_p84']} corr@lin={r['at_linear_bound']['correction']:+.3f} ESS={r['at_linear_bound']['ess']:.0f} ESSmin={r['ess_min']:.0f}", flush=True)


def _json_default(o):
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(str(type(o)))


def rms_table(out):
    """rms of the prior Delta rho at each z, at each bound: {fam|l|rule: array over Z_EVAL or None}."""
    T = {}
    for fam in FAMILIES:
        for l in ELLS:
            for rule in RULES:
                b = out["bounds"][f"{fam}|l_{l:g}"]["exact_5pct"] if rule == "5pct" else out["linear"][f"{fam}|l_{l:g}"][rule]
                T[f"{fam}|l_{l:g}|{rule}"] = (b * np.sqrt(2.0 * (1.0 - kfun(fam, U_EVAL, l)))) if isinstance(b, float) and b > 0 else None
    return T


def universality(out) -> dict:
    T = rms_table(out)
    rec = {"rms_prior_at_bound": {k: (v.tolist() if v is not None else None) for k, v in T.items()}}
    for rule in RULES:
        sets = {"l_0.3_1_3": [(f, l) for f in FAMILIES for l in (0.3, 1.0, 3.0)], "all_l": [(f, l) for f in FAMILIES for l in ELLS]}
        R = {}
        for sname, members in sets.items():
            M = [T[f"{f}|l_{l:g}|{rule}"] for f, l in members]
            missing = [f"{f}|l_{l:g}" for (f, l), m in zip(members, M) if m is None]
            M = np.array([m for m in M if m is not None])
            if len(M) < 2:
                R[sname] = {"n": int(len(M)), "missing": missing, "ratio_vs_z": None}
                continue
            ratio = (M.max(axis=0) / np.maximum(M.min(axis=0), 1e-300)).tolist()
            k = int(np.argmin(ratio))
            R[sname] = {"n": int(len(M)), "missing": missing, "ratio_vs_z": ratio, "z_star": Z_EVAL[k], "ratio_at_z_star": ratio[k],
                        "rms_range_at_z_star": [float(M[:, k].min()), float(M[:, k].max())], "median_at_z_star": float(np.median(M[:, k]))}
        # ratio of the all-l set at the z_* of the 0.3-1-3 set
        if R["l_0.3_1_3"].get("ratio_vs_z") and R["all_l"].get("ratio_vs_z"):
            k = Z_EVAL.index(R["l_0.3_1_3"]["z_star"])
            R["all_l"]["ratio_at_z_star_of_l_0.3_1_3"] = R["all_l"]["ratio_vs_z"][k]
        rec[rule] = R
    return rec


def summary(out: dict) -> None:
    f5 = lambda v: "sin cota (> 1.0)" if v is None else (f"< {v[1:]}" if isinstance(v, str) else f"{v:.3f}")
    L = [f"chi2_obs = {out['chi2_obs']:.3f} (malla media {out['chi2_obs_half_grid']:.3f}); {out['n_real']} caminos por celda; sigma_X en {out['sigmas']}", ""]
    L.append("Escalas efectivas (en unidades de l): familia  r(e^-1)/l  integral/l  [esperado]")
    for fam in FAMILIES:
        s = out["checks"]["scales"][fam]["l_1"]
        L.append(f"   {fam:9s} {s['r_e_over_l']:.3f} {s['integral_over_l']:.3f}  {out['checks']['scales_expected'][fam]}")
    L.append("")
    L.append("Cotas de sigma_X: familia  l  | 5 % (boot 16-84) | lineal marg -1.92 | -1.35 | max lineal | P por sigma | H2<=0 por sigma")
    for fam in FAMILIES + tuple(EXTRA):
        for l in (ELLS if fam in FAMILIES else EXTRA[fam]):
            b = out["bounds"][f"{fam}|l_{l:g}"]
            lin = out["linear"].get(f"{'Steps' if fam == 'StepsGauss' else fam}|l_{l:g}", {})  # StepsGauss: same kernel as Steps
            bp = b["boot_p16_p84"]
            bps = f"({bp[0]:.3f}-{bp[1]:.3f})" if bp else "(-)"
            L.append(f"   {fam:10s} {l:<4g} | {f5(b['exact_5pct']):>16s} {bps:>14s} | {PS.fmt(lin.get('192')):>8s} | {PS.fmt(lin.get('135')):>8s} | {lin.get('sigma_at_max', float('nan')):g} {lin.get('max_class', '')} | "
                     f"{[round(p, 3) for p in b['P_vs_sigma']]} | {[round(p, 2) for p in b['frac_H2_le_0_vs_sigma']]}")
            if fam == "Steps" and l in HALF_ELLS:
                bh = out["bounds"][f"{fam}|l_{l:g}|half"]
                L.append(f"   {'Steps/2':10s} {l:<4g} | {f5(bh['exact_5pct']):>16s} {'':>14s} |          |          | malla media | {[round(p, 3) for p in bh['P_vs_sigma']]}")
    L.append("")
    U = out["universality"]
    for rule in RULES:
        R = U[rule]
        for sname in ("l_0.3_1_3", "all_l"):
            r = R[sname]
            if r.get("ratio_vs_z"):
                L.append(f"regla {rule}, conjunto {sname} (n={r['n']}, faltan {r['missing']}): z_* = {r['z_star']}, factor max/min = {r['ratio_at_z_star']:.2f}, "
                         f"rms en z_* {r['rms_range_at_z_star'][0]:.3f}-{r['rms_range_at_z_star'][1]:.3f} (mediana {r['median_at_z_star']:.3f}); factor por z: {[round(x, 2) for x in r['ratio_vs_z']]}"
                         + (f"; en el z_* del conjunto 0.3-1-3: {r['ratio_at_z_star_of_l_0.3_1_3']:.2f}" if "ratio_at_z_star_of_l_0.3_1_3" in r else ""))
            else:
                L.append(f"regla {rule}, conjunto {sname}: n={r['n']} < 2 (faltan {r['missing']})")
    L.append("")
    zs = U["5pct"]["l_0.3_1_3"].get("z_star")
    if zs is not None:
        k = Z_EVAL.index(zs)
        L.append(f"rms del prior de Delta rho en z_* = {zs} (regla 5 % / -1.92 / -1.35) y p95 |Delta rho| de los supervivientes en la cota 5 %, rho<=0 entre supervivientes:")
        for fam in FAMILIES:
            for l in ELLS:
                vals = []
                for rule in RULES:
                    v = U["rms_prior_at_bound"][f"{fam}|l_{l:g}|{rule}"]
                    vals.append("  -  " if v is None else f"{v[k]:.3f}")
                ab = out["at_bound"].get(f"{fam}|l_{l:g}|5pct")
                p95 = f"{ab['p95_abs_drho_survivors'][k]:.3f}" if ab and ab["p95_abs_drho_survivors"] else "  -  "
                rn = f"{ab['frac_rho_le_0_among_survivors']:.3f}" if ab and ab["frac_rho_le_0_among_survivors"] is not None else "  -  "
                L.append(f"   {fam:9s} l={l:<4g} rms: {' / '.join(vals)}   p95 surv: {p95}   rho<=0 surv: {rn}")
    L.append("")
    L.append("Exacta guiada en l = 1 (N = 5000): familia | lineal -1.92 (analitica) | lineal en los mismos nodos | exacta -1.92 (boot 16-84) | exacta solo nodos con ESS >= 1 % N | -1.35 | correccion ln(Zex/Zlin) en la cota lineal | ESS en la cota (min) | nodos sigma")
    for fam, r in out.get("guided_l1", {}).items():
        a = r["at_linear_bound"]
        bp = r["exact_192_boot_p16_p84"]
        L.append(f"   {fam:9s} | {PS.fmt(r['sigma_linear_192'])} | {PS.fmt(r['lin_192_from_nodes'])} | {PS.fmt(r['exact_192'])} {[round(v, 3) for v in bp] if bp else None} | {PS.fmt(r.get('exact_192_trusted'))} | {PS.fmt(r['exact_135'])} | {a['correction']:+.3f} | {a['ess']:.0f} ({r['ess_min']:.0f}) | {[round(x, 3) for x in r['sigmas']]} ESS {[round(e) for e in r['ess_vs_sigma']]}")
    L.append("")
    rep = out["checks"]["reproduce_exacto_sn"]
    L.append(f"reproduce exacto_sn (OU, {len(rep)} celdas): max |dchi2| = {max(v['max_abs_dchi2'] for v in rep.values()):.2e}")
    L.append("P en sigma = 0.02: " + json.dumps({k: round(v, 3) for k, v in out["checks"]["P_at_sigma_0.02"].items()}))
    L.append("jitter gaussiano: " + json.dumps({k: v.get("jitter_rel") for k, v in out["checks"]["samplers"]["Gauss"].items()}))
    L.append(f"tiempo: {out['runtime_s']:.0f} s")
    (HERE / "familias.txt").write_text("\n".join(L) + "\n", encoding="utf-8")
    print("\n".join(L))


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--txt-only":
        summary(json.loads((HERE / "familias.json").read_text()))
    elif len(sys.argv) > 1 and sys.argv[1] == "--guided-only":
        out = json.loads((HERE / "familias.json").read_text())
        init_linear()
        run_guided(out)
        (HERE / "familias.json").write_text(json.dumps(out, indent=1, default=_json_default), encoding="utf-8")
        summary(out)
    else:
        main()
