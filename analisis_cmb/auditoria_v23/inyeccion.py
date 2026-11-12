#!/usr/bin/env python3
"""
Exploration only (INYECCION.md): fake universes where the truth is known. Which rule gives correct bounds and which one
confuses tension with jitter? Nothing in results/ or papers changes.

Data vector: BAO (13) + theta_* + omega_m + Pantheon+ (1580 rows), generated from the LCDM best fit of the real data (muestras_sn
make_linear: P0, B0, Bm, A, C) with Gaussian noise from the real covariances; 200 universes per scenario, the SAME noise vector per
universe in every scenario.
Scenarios:
  A   pure LCDM.
  B   an OU theta = 1 path (one per universe, convention B: X - X(0)) injected through the EXACT model, sigma_true in {0.1, 0.2, 0.3}.
  C   tension of the DES-SN5YR kind, no jitter: the SN block generated with Omega_m = Omega_m(best fit) + delta, delta chosen so that
      the mean tension Delta chi^2 (joint - 15 rows alone - SN alone) is 8.7; BAO and CMB unchanged.
  D   the SN of A compressed into 22 redshift bins (edges at the geometric midpoints of the Union3 nodes) with the GLS compression
      T = (E' W E)^-1 E' W, covariance (E' W E)^-1 (no tension, less information).
Per universe, OU theta = 1, linear version (chi^2_min of each path as a quadratic form with (h, omega_m, M) profiled, which equals
the Gaussian marginalisation up to a constant; the SAME 2000 unit paths scaled by sigma for all universes):
  - 5 % rule: chi^2_i(sigma) - chi^2_LCDM = -2 sigma b_i + sigma^2 c_i <= 0 iff sigma <= s_i = 2 b_i / c_i, so
    P(sigma) = frac(s_i >= sigma) and the 5 % bound is the 95th percentile of max(s_i, 0) (closed form, no grid).
  - marginal -1.92 (primary): analytic linear marginal of marginal_guiada.lnZ_lin (C + sigma^2 A K A', flat priors on dp) on a
    dense sigma grid, cubic spline, drop 1.92 from the global maximum; also the Monte Carlo marginal over the 2000 paths
    (logmeanexp of sigma b_i - sigma^2 c_i / 2) with its ESS, as a check.
  - Delta ln L max = max over sigma of ln L(sigma) - ln L(0); "sigma_X = 0 excluded" iff Delta ln L max > 1.92.
  - tension Delta chi^2 (linear) for every universe.
Checks: (1) real data: the linear 5 % and marginal reproduce the exact bounds (0.169, 0.287) within 10 % (stop otherwise);
(2) 5 universes per scenario: the exact 5 % rule (nonlinear fit per path, familias.py style, discarding H^2 <= 0) reproduces the
linear one within 10 %.

    ../venv/bin/python -I inyeccion.py --check     # real-data check only
    ../venv/bin/python -I inyeccion.py             # everything: writes inyeccion.{json,txt}
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
from scipy.interpolate import CubicSpline
from scipy.linalg import cho_factor, cho_solve, solve_triangular
from scipy.optimize import brentq, least_squares
from scipy.special import logsumexp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import exacto as X  # noqa: E402
import exacto_sn as S  # noqa: E402
import marginal_guiada as MG  # noqa: E402
import muestras_sn as MS  # noqa: E402

G0 = S.G0
OMR = S.OMR
SEED_I = 20261118
N_UNI, N_PATH, N_EXACT = 200, 2000, 5
SIG_TRUE = (0.1, 0.2, 0.3)
TENSION_TARGET = 8.7
GRID = np.concatenate([[0.0], np.geomspace(0.02, 1.5, 90)])
REAL = {"5pct_exact": 0.169, "marg_exact": 0.287}
DROP = 1.92


# ----------------------------------------------------------------------------------------------------------------------
# settings (full Pantheon+ and binned)
# ----------------------------------------------------------------------------------------------------------------------
def quad_perp(CF, Bm, R):
    """r' W_perp r for every column of R, W_perp = W - W Bm (Bm' W Bm)^-1 Bm' W, W = C^-1 (profile over the columns of Bm)."""
    WR = cho_solve(CF, R)
    WB = cho_solve(CF, Bm)
    Fm = Bm.T @ WB
    BWR = Bm.T @ WR
    return np.einsum("ij,ij->j", R, WR) - np.einsum("ij,ij->j", BWR, np.linalg.solve(Fm, BWR))


def wperp(CF, Bm, R):
    WR = cho_solve(CF, R)
    WB = cho_solve(CF, Bm)
    return WR - WB @ np.linalg.solve(Bm.T @ WB, Bm.T @ WR)


def make_full():
    sn = S.SN["all"]
    L = MS.make_linear(sn)
    return {"name": "full", "sn": sn, "D": L["D"], "C": L["C"], "CF": L["CF"], "P0": L["P0"], "B0": L["B0"], "Bm": L["Bm"], "A": L["A"],
            "nsn": sn["n"], "T": None, "Lc": np.linalg.cholesky(L["C"])}


def make_binned(full):
    sn = full["sn"]
    znodes = MS.load_union3()["zhd"]
    edges = np.concatenate([[0.01], np.sqrt(znodes[:-1] * znodes[1:]), [max(sn["zhd"].max(), znodes[-1]) + 1e-6]])
    idx = np.digitize(sn["zhd"], edges) - 1
    counts = np.bincount(idx, minlength=len(edges) - 1)
    keep = np.nonzero(counts > 0)[0]
    remap = {b: i for i, b in enumerate(keep)}
    nb = len(keep)
    E = np.zeros((sn["n"], nb))
    E[np.arange(sn["n"]), [remap[b] for b in idx]] = 1.0
    W = sn["Ci"]
    Cb = np.linalg.inv(E.T @ W @ E)
    Cb = 0.5 * (Cb + Cb.T)
    T = Cb @ E.T @ W
    n15 = 15
    C = np.zeros((n15 + nb, n15 + nb))
    C[:n15, :n15] = full["C"][:n15, :n15]
    C[n15:, n15:] = Cb
    B0 = np.concatenate([full["B0"][:n15], T @ full["B0"][n15:]])
    Bm = np.vstack([full["Bm"][:n15], T @ full["Bm"][n15:]])
    A = np.vstack([full["A"][:n15], T @ full["A"][n15:]])
    return {"name": "binned", "sn": sn, "C": C, "CF": cho_factor(C, lower=True), "P0": full["P0"], "B0": B0, "Bm": Bm, "A": A, "nsn": nb,
            "T": T, "Lb": np.linalg.cholesky(Cb), "edges": edges.tolist(), "counts": counts.tolist(), "nb": nb}


def compress(setting, Dfull):
    """Columns of Dfull (nrows_full x k) -> binned data vector(s)."""
    T = setting["T"]
    if T is None:
        return Dfull
    return np.vstack([Dfull[:15], T @ Dfull[15:]])


# ----------------------------------------------------------------------------------------------------------------------
# linear machinery for one setting
# ----------------------------------------------------------------------------------------------------------------------
class Linear:
    def __init__(self, setting, Xp):
        self.s = setting
        A, Bm, CF, C = setting["A"], setting["Bm"], setting["CF"], setting["C"]
        K = MG.kernel(1.0)
        Msig = (A @ K) @ A.T
        self.Msig = 0.5 * (Msig + Msig.T)
        self.Y = A @ Xp.T  # (nrows, N_PATH): response of each unit path
        WY = wperp(CF, Bm, self.Y)
        self.c = np.einsum("ij,ij->j", self.Y, WY)
        # analytic marginal: per grid sigma, factorisations shared by all universes
        self.fac = []
        for sg in GRID:
            Sf = cho_factor(C + sg * sg * self.Msig, lower=True)
            SB = cho_solve(Sf, Bm)
            Fm = Bm.T @ SB
            self.fac.append((Sf, SB, Fm, 2 * np.sum(np.log(np.diag(Sf[0]))) + np.linalg.slogdet(Fm)[1]))
        # tension blocks
        self.CF15 = cho_factor(C[:15, :15], lower=True)
        self.CFsn = cho_factor(C[15:, 15:], lower=True)

    def lnz_profile(self, R):
        """ln Z(sigma) - ln Z(0) for the columns of R (residuals to B0): (len(GRID), k)."""
        Bm = self.s["Bm"]
        out = np.zeros((len(GRID), R.shape[1]))
        for g, (Sf, SB, Fm, ld) in enumerate(self.fac):
            SR = cho_solve(Sf, R)
            b = Bm.T @ SR
            out[g] = -0.5 * (np.einsum("ij,ij->j", R, SR) - np.einsum("ij,ij->j", b, np.linalg.solve(Fm, b)) + ld)
        return out - out[0]

    def lnz_single(self, r, sg):
        Sf = cho_factor(self.s["C"] + sg * sg * self.Msig, lower=True)
        return MG.lnZ_lin(Sf, self.s["Bm"], r)

    def analyse(self, R):
        """All linear statistics for the universes in the columns of R (residuals to B0)."""
        s = self.s
        Wr = wperp(s["CF"], s["Bm"], R)  # (nrows, k)
        a = np.einsum("ij,ij->j", R, Wr)
        b = self.Y.T @ Wr  # (N_PATH, k)
        sstar = np.maximum(2.0 * b / self.c[:, None], 0.0)
        b5 = np.quantile(sstar, 0.95, axis=0)
        prof = self.lnz_profile(R)
        k = R.shape[1]
        res = {"chi2_lcdm": a, "b5": b5, "marg": np.full(k, np.nan), "marg_class": [""] * k, "dlnL_max": prof.max(axis=0),
               "sigma_max": GRID[np.argmax(prof, axis=0)], "excl0": prof.max(axis=0) > DROP,
               "mc_marg": np.full(k, np.nan), "mc_dlnL_max": np.zeros(k), "mc_ess_at_bound": np.zeros(k), "mc_excl0": np.zeros(k, bool)}
        for j in range(k):
            res["marg"][j], res["marg_class"][j] = drop_bound(prof[:, j])
            # Monte Carlo marginal over the same paths
            lw = GRID[:, None] * b[:, j][None, :] - 0.5 * GRID[:, None] ** 2 * self.c[None, :]  # (grid, paths)
            lz = logsumexp(lw, axis=1) - np.log(N_PATH)
            lz -= lz[0]
            res["mc_dlnL_max"][j] = lz.max()
            res["mc_excl0"][j] = lz.max() > DROP
            res["mc_marg"][j], _ = drop_bound(lz)
            if np.isfinite(res["mc_marg"][j]):
                g = int(np.argmin(np.abs(GRID - res["mc_marg"][j])))
                w = np.exp(lw[g] - lw[g].max())
                res["mc_ess_at_bound"][j] = w.sum() ** 2 / np.sum(w**2)
        # tension: joint (3 params) - 15 rows alone (h, omega_m) - SN alone (omega_m, M at fixed h)
        c_joint = a
        c_15 = quad_perp(self.CF15, s["Bm"][:15, :2], R[:15])
        c_sn = quad_perp(self.CFsn, s["Bm"][15:, 1:], R[15:])
        res["tension"] = c_joint - c_15 - c_sn
        res["chi2_15"] = c_15
        res["chi2_sn"] = c_sn
        return res


def drop_bound(prof):
    """Upper bound where ln L drops by DROP from its global maximum (cubic spline on GRID). Returns (sigma or nan, class)."""
    k0 = int(np.argmax(prof))
    cls = "sigma->0" if k0 == 0 else ("upper edge" if k0 == len(GRID) - 1 else "interior")
    thr = prof[k0] - DROP
    below = np.nonzero((prof < thr) & (np.arange(len(prof)) > k0))[0]
    if len(below) == 0:
        return np.nan, cls + " (no crossing)"
    kk = int(below[0])
    cs = CubicSpline(GRID, prof)
    return float(brentq(lambda s_: cs(s_) - thr, GRID[kk - 1], GRID[kk], xtol=1e-8)), cls


# ----------------------------------------------------------------------------------------------------------------------
# exact (nonlinear) fits for a fake universe
# ----------------------------------------------------------------------------------------------------------------------
EX = {}  # filled in workers by fork: paths, settings


def resid_exact(p, dX, data, setting):
    h, wm, M = p
    sn = setting["sn"]
    hist, m = S.history(G0, wm / h**2, OMR / h**2, dX)
    nr = 15 + setting["nsn"]
    if hist is None:
        return np.full(nr, 1e3 * (1.0 - m))
    dm_u, inv = hist
    v = S.vec15(G0, h, wm, dm_u, inv)
    r = data[:15] - v
    out = [solve_triangular(X.L13, r[:13], lower=True, check_finite=False), [r[13] / X.SIG_TH, r[14] / X.SIG_WM]]
    mu = S.mu_sn(G0, sn, dm_u) + M
    if setting["T"] is None:
        out.append(solve_triangular(setting["sn"]["L"], data[15:] - mu, lower=True, check_finite=False))
    else:
        out.append(solve_triangular(setting["Lb"], data[15:] - setting["T"] @ mu, lower=True, check_finite=False))
    return np.concatenate(out)


def fit_exact(dX, data, setting, x0):
    r = least_squares(resid_exact, np.array(x0), bounds=([0.4, 0.05, -50.0], [1.0, 0.30, 50.0]), args=(dX, data, setting),
                      method="trf", x_scale=[0.01, 0.001, 0.01], ftol=1e-10, xtol=1e-10, max_nfev=200)
    return float(2.0 * r.cost), r.x


def _exact_task(args):
    scen, u, sg = args
    setting = EX["settings"][EX["scen_setting"][scen]]
    data = EX["data"][scen][:, u]
    Xp = EX["Xp"]
    t0 = time.time()
    c0, p0 = fit_exact(np.zeros(G0.N), data, setting, setting["P0"])
    dX = sg * (Xp - Xp[:, [0]])
    h0, wm0 = p0[0], p0[1]
    om0, or0 = wm0 / h0**2, OMR / h0**2
    ok = (om0 * G0.p3 + or0 * G0.p4 + (1 - om0 - or0) + dX).min(axis=1) > 0.0
    chi = np.full(N_PATH, np.nan)
    for i in np.nonzero(ok)[0]:
        chi[i] = fit_exact(dX[i], data, setting, p0)[0]
    P = float(np.mean(chi[ok] <= c0))
    return {"scen": scen, "u": u, "sigma": sg, "chi2_lcdm": c0, "P": P, "n_ok": int(ok.sum()), "seconds": time.time() - t0}


def cross_5pct(sigs, P):
    """5 % crossing interpolated in ln sigma (as exacto_sn); None if P never crosses 0.05 inside the range."""
    sigs, P = np.asarray(sigs), np.asarray(P)
    for k in range(1, len(sigs)):
        if P[k - 1] >= 0.05 > P[k]:
            x0, x1 = np.log(sigs[k - 1]), np.log(sigs[k])
            return float(np.exp(x0 + (0.05 - P[k - 1]) * (x1 - x0) / (P[k] - P[k - 1])))
    return None


# ----------------------------------------------------------------------------------------------------------------------
# scenarios
# ----------------------------------------------------------------------------------------------------------------------
def tension_noise_free(full, lin, delta):
    sn, P0 = full["sn"], full["P0"]
    h0, wm0, M0 = P0
    om0 = wm0 / h0**2
    wm_c = (om0 + delta) * h0**2
    Bc = MS.exact_model(sn, h0, wm_c, M0, np.zeros((1, G0.N)))[0][0]
    Bc[:15] = full["B0"][:15]
    r = (Bc - full["B0"])[:, None]
    return float(quad_perp(full["CF"], full["Bm"], r)[0] - quad_perp(lin.CF15, full["Bm"][:15, :2], r[:15])[0]
                 - quad_perp(lin.CFsn, full["Bm"][15:, 1:], r[15:])[0]), Bc


def summarise(name, res, sig_true=None):
    def q(v):
        v = np.asarray(v, float)
        fin = v[np.isfinite(v)]
        return {"median": float(np.median(fin)) if len(fin) else None, "p16": float(np.quantile(fin, 0.16)) if len(fin) else None,
                "p84": float(np.quantile(fin, 0.84)) if len(fin) else None, "n_finite": int(len(fin)), "n": int(len(v))}
    out = {"scenario": name, "sigma_true": sig_true, "n_universes": int(len(res["b5"])),
           "5pct": q(res["b5"]), "marg": q(res["marg"]), "mc_marg": q(res["mc_marg"]),
           "frac_5pct_above_1": float(np.mean(res["b5"] > 1.0)), "marg_no_crossing": int(np.sum(~np.isfinite(res["marg"]))),
           "dlnL_max": q(res["dlnL_max"]), "frac_excl0": float(np.mean(res["excl0"])), "frac_excl0_mc": float(np.mean(res["mc_excl0"])),
           "tension": {"mean": float(np.mean(res["tension"])), "median": float(np.median(res["tension"]))},
           "chi2_lcdm_mean": float(np.mean(res["chi2_lcdm"])), "mc_ess_at_bound_median": float(np.nanmedian(res["mc_ess_at_bound"]))}
    if sig_true is not None:
        out["coverage_5pct"] = float(np.mean(res["b5"] >= sig_true))
        out["coverage_marg"] = float(np.mean(~np.isfinite(res["marg"]) | (res["marg"] >= sig_true)))
        out["coverage_mc_marg"] = float(np.mean(~np.isfinite(res["mc_marg"]) | (res["mc_marg"] >= sig_true)))
        out["5pct_over_true"] = q(res["b5"] / sig_true)
        out["marg_over_true"] = q(res["marg"] / sig_true)
    return out


def main(check_only: bool) -> None:
    t_start = time.time()
    full = make_full()
    rng = np.random.default_rng([SEED_I, 1])
    Xp = S.ou_paths(G0, 1.0, 1.0, N_PATH, rng)  # unit paths, convention B applied by A (rows sum to zero) / explicitly in exact
    lin = Linear(full, Xp)
    # ---- real-data check
    r_real = (full["D"] - full["B0"])[:, None]
    real = lin.analyse(r_real)
    real_rec = {"5pct_linear": float(real["b5"][0]), "marg_linear": float(real["marg"][0]), "marg_linear_class": real["marg_class"][0],
                "mc_marg_linear": float(real["mc_marg"][0]), "mc_ess_at_bound": float(real["mc_ess_at_bound"][0]),
                "dlnL_max_linear": float(real["dlnL_max"][0]), "sigma_max": float(real["sigma_max"][0]), "tension_linear": float(real["tension"][0]),
                "chi2_lcdm_linear": float(real["chi2_lcdm"][0]), "chi2_lcdm_exact": float(S.fit_combo(G0, np.zeros(G0.N), full["sn"], True)[0]),
                "ref": REAL, "muestras_sn_linear_192": 0.277}
    real_rec["5pct_ratio"] = real_rec["5pct_linear"] / REAL["5pct_exact"]
    real_rec["marg_ratio"] = real_rec["marg_linear"] / REAL["marg_exact"]
    # brentq on the exact linear profile (as muestras_sn) to validate the spline
    thr = real_rec["dlnL_max_linear"] - DROP
    z0 = lin.lnz_single(r_real[:, 0], 0.0)
    prof = lin.lnz_profile(r_real)[:, 0]
    k0 = int(np.argmax(prof))
    kk = int(np.nonzero((prof < thr) & (np.arange(len(prof)) > k0))[0][0])
    real_rec["marg_linear_brentq"] = float(brentq(lambda s_: lin.lnz_single(r_real[:, 0], s_) - z0 - thr, GRID[kk - 1], GRID[kk], xtol=1e-9))
    real_rec["ok_within_10pct"] = bool(abs(real_rec["5pct_ratio"] - 1) < 0.10 and abs(real_rec["marg_ratio"] - 1) < 0.10)
    print("REAL DATA CHECK:", json.dumps(real_rec, indent=1), flush=True)
    if check_only:
        (HERE / "inyeccion_check.json").write_text(json.dumps(real_rec, indent=1), encoding="utf-8")
        return
    if not real_rec["ok_within_10pct"]:
        print("STOP: linear version does not reproduce the real bounds within 10 %", flush=True)
        (HERE / "inyeccion_check.json").write_text(json.dumps(real_rec, indent=1), encoding="utf-8")
        return
    # ---- binned setting
    binned = make_binned(full)
    lin_b = Linear(binned, Xp)
    # ---- noise and injected paths
    rng_n = np.random.default_rng([SEED_I, 2])
    Z = rng_n.standard_normal((len(full["D"]), N_UNI))
    noise = full["Lc"] @ Z  # same noise per universe in every scenario
    rng_p = np.random.default_rng([SEED_I, 3])
    Xinj = S.ou_paths(G0, 1.0, 1.0, N_UNI, rng_p)  # unit injected paths, one per universe
    dXinj = Xinj - Xinj[:, [0]]
    # ---- scenario C: delta
    f = lambda d: tension_noise_free(full, lin, d)[0] - (TENSION_TARGET - 1.0)
    delta = brentq(f, 1e-4, 0.2, xtol=1e-8)
    lam, Bc = tension_noise_free(full, lin, delta)
    om0 = full["P0"][1] / full["P0"][0] ** 2
    print(f"delta = {delta:.5f} (Omega_m {om0:.4f} -> {om0 + delta:.4f}); noise-free tension {lam:.3f}, expected mean {lam + 1:.3f}", flush=True)
    # ---- data vectors per scenario (full rows), then residuals to B0 of the setting
    data = {"A": full["B0"][:, None] + noise}
    inj_bad = {}
    for st in SIG_TRUE:
        sig = np.zeros_like(noise)
        nbad = 0
        for u in range(N_UNI):
            out, bad = MS.exact_model(full["sn"], full["P0"][0], full["P0"][1], full["P0"][2], (st * dXinj[u])[None, :])
            if bad[0]:
                nbad += 1
                out = full["B0"][None, :]
            sig[:, u] = out[0] - full["B0"]
        inj_bad[f"B_{st:g}"] = nbad
        data[f"B_{st:g}"] = full["B0"][:, None] + sig + noise
    data["C"] = Bc[:, None] + noise
    data["D"] = compress(binned, data["A"])
    scen_setting = {k: ("binned" if k == "D" else "full") for k in data}
    settings = {"full": full, "binned": binned}
    lins = {"full": lin, "binned": lin_b}
    # ---- linear analysis
    out = {"seed": SEED_I, "n_universes": N_UNI, "n_paths": N_PATH, "grid": GRID.tolist(), "real_check": real_rec, "delta": delta,
           "omega_m_fit": om0, "omega_m_C": om0 + delta, "tension_noise_free_C": lam, "tension_target": TENSION_TARGET,
           "injection_H2_le_0": inj_bad, "bins": {"edges": binned["edges"], "counts": binned["counts"], "nb": binned["nb"]},
           "scenarios": {}, "per_universe": {}}
    for scen, Dm in data.items():
        st = settings[scen_setting[scen]]
        R = Dm - st["B0"][:, None]
        t0 = time.time()
        res = lins[scen_setting[scen]].analyse(R)
        sig_true = float(scen.split("_")[1]) if scen.startswith("B_") else None
        out["scenarios"][scen] = summarise(scen, res, sig_true)
        out["scenarios"][scen]["seconds"] = time.time() - t0
        out["per_universe"][scen] = {k: (np.asarray(v).tolist() if not isinstance(v, list) else v) for k, v in res.items()}
        print(f"{scen}: {json.dumps({k: v for k, v in out['scenarios'][scen].items() if k in ('5pct', 'marg', 'frac_excl0', 'tension', 'coverage_5pct', 'coverage_marg')})}", flush=True)
    # ratios vs A
    A5, Am = out["scenarios"]["A"]["5pct"]["median"], out["scenarios"]["A"]["marg"]["median"]
    for scen in out["scenarios"]:
        out["scenarios"][scen]["ratio_vs_A"] = {"5pct": out["scenarios"][scen]["5pct"]["median"] / A5,
                                                "marg": (out["scenarios"][scen]["marg"]["median"] / Am) if out["scenarios"][scen]["marg"]["median"] else None}
    # ---- exact check on N_EXACT universes per scenario
    EX["settings"], EX["scen_setting"], EX["data"], EX["Xp"] = settings, scen_setting, data, Xp
    tasks = []
    for scen in data:
        for u in range(N_EXACT):
            b5 = out["per_universe"][scen]["b5"][u]
            for sg in np.geomspace(0.6 * b5, 1.6 * b5, 6):
                tasks.append((scen, u, float(sg)))
    t0 = time.time()
    with Pool(min(16, os.cpu_count() or 1)) as pool:
        recs = pool.map(_exact_task, tasks, chunksize=1)
    exact = {}
    for scen in data:
        exact[scen] = []
        for u in range(N_EXACT):
            rr = sorted([r for r in recs if r["scen"] == scen and r["u"] == u], key=lambda r: r["sigma"])
            sigs, P = [r["sigma"] for r in rr], [r["P"] for r in rr]
            b_ex = cross_5pct(sigs, P)
            b_lin = out["per_universe"][scen]["b5"][u]
            exact[scen].append({"u": u, "sigmas": sigs, "P": P, "n_ok_min": min(r["n_ok"] for r in rr), "chi2_lcdm_exact": rr[0]["chi2_lcdm"],
                                "chi2_lcdm_linear": out["per_universe"][scen]["chi2_lcdm"][u], "5pct_exact": b_ex, "5pct_linear": b_lin,
                                "ratio": (b_ex / b_lin) if b_ex else None})
    ratios = [e["ratio"] for s in exact.values() for e in s if e["ratio"]]
    out["exact_check"] = {"per_scenario": exact, "seconds": time.time() - t0, "n_with_crossing": len(ratios), "n_total": N_EXACT * len(data),
                          "ratio_min": min(ratios) if ratios else None, "ratio_max": max(ratios) if ratios else None,
                          "all_within_10pct": bool(ratios and max(abs(r - 1) for r in ratios) < 0.10)}
    out["runtime_s"] = time.time() - t_start
    (HERE / "inyeccion.json").write_text(json.dumps(out, indent=1, default=_json_default), encoding="utf-8")
    write_txt(out)
    print(f"done in {out['runtime_s']:.0f} s", flush=True)


def _json_default(o):
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(type(o))


def fmt(v, d=3):
    return "—" if v is None or (isinstance(v, float) and not np.isfinite(v)) else f"{v:.{d}f}"


def write_txt(out):
    L = [f"INYECCION: {out['n_universes']} universos por escenario, {out['n_paths']} caminos OU theta = 1, semilla {out['seed']}",
         f"real: 5 % lineal {out['real_check']['5pct_linear']:.4f} (exacta {REAL['5pct_exact']}), marginal lineal {out['real_check']['marg_linear']:.4f} "
         f"(exacta {REAL['marg_exact']}; brentq {out['real_check']['marg_linear_brentq']:.4f}; MC {out['real_check']['mc_marg_linear']:.4f}), "
         f"dlnL max {out['real_check']['dlnL_max_linear']:+.3f}, tension lineal {out['real_check']['tension_linear']:.2f}",
         f"delta (C) = {out['delta']:.5f}; Omega_m {out['omega_m_fit']:.4f} -> {out['omega_m_C']:.4f}; tension sin ruido {out['tension_noise_free_C']:.3f}",
         f"tramos (D): {out['bins']['nb']} con SN; cuentas {out['bins']['counts']}",
         "", "escenario  sigma_true | 5 % mediana (16-84)        | marginal mediana (16-84)    | MC marg.  | cobertura 5 % / marg | frac sigma=0 excl. | dlnL max med | tension media | /A 5 % | /A marg"]
    for scen, s in out["scenarios"].items():
        cov = f"{s['coverage_5pct']:.3f} / {s['coverage_marg']:.3f}" if s.get("coverage_5pct") is not None else "—"
        L.append(f"{scen:9s} {fmt(s['sigma_true'], 1):>5s} | {fmt(s['5pct']['median'])} ({fmt(s['5pct']['p16'])}-{fmt(s['5pct']['p84'])}) | "
                 f"{fmt(s['marg']['median'])} ({fmt(s['marg']['p16'])}-{fmt(s['marg']['p84'])}) | {fmt(s['mc_marg']['median'])} | {cov:>15s} | "
                 f"{s['frac_excl0']:.3f} (MC {s['frac_excl0_mc']:.3f}) | {s['dlnL_max']['median']:+.2f} | {s['tension']['mean']:.2f} | "
                 f"{s['ratio_vs_A']['5pct']:.3f} | {fmt(s['ratio_vs_A']['marg'])}")
    ex = out["exact_check"]
    L += ["", f"comprobacion exacta (5 universos/escenario, 5 % exacta / lineal): {ex['n_with_crossing']}/{ex['n_total']} con cruce, "
          f"cociente {fmt(ex['ratio_min'])}-{fmt(ex['ratio_max'])}, todos a <10 %: {ex['all_within_10pct']}; {ex['seconds']:.0f} s"]
    for scen, lst in ex["per_scenario"].items():
        L.append(f"  {scen:7s} " + "  ".join(f"u{e['u']}: {fmt(e['5pct_exact'])}/{e['5pct_linear']:.3f} ({fmt(e['ratio'], 2)})" for e in lst))
    L.append(f"tiempo total {out['runtime_s']:.0f} s")
    (HERE / "inyeccion.txt").write_text("\n".join(L) + "\n", encoding="utf-8")
    print("\n".join(L))


if __name__ == "__main__":
    main("--check" in sys.argv[1:])
