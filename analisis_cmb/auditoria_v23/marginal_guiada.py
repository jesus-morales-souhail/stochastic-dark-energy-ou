#!/usr/bin/env python3
"""
Exploration only (MARGINAL_GUIADA.md): does the exact marginal likelihood agree with the analytic linear one? Importance
sampling guided by the linear solution (not prior paths). Nothing in results/ or papers changes.

Combination BAO + theta_* + omega_m + SN (1595 rows), theta in {0.3, 1, 3, 10}, the sigma_X cells of convergencia_sn.py plus
sigma_X = 0 and 0.01 (checks), N = 5000 samples per cell.
1. Linear model of convention B around the LCDM best fit p0 = (h, omega_m, M) of this data set: d = b0 + B dp + A x + n, with x
   the path on the exacto.py grid (3506 nodes, DESI z exact), A the full linear response in raw units (delta y = y delta ln D;
   SN rows (5/ln 10) delta ln D_M interpolated between bracketing nodes; rows sum to zero: Lambda absorbs X(0)), B the finite-
   difference derivatives of the exact model, prior x ~ N(0, sigma_X^2 K_theta) with K the OU kernel on the grid, flat priors
   on dp, n ~ N(0, C), C = blockdiag(official 13x13, sigma_theta^2, sigma_omega^2, C_SN STAT+SYS).
2. (dp, x) sampled from the linear posterior: dp from its Gaussian marginal N(dp_hat, (B' S^-1 B)^-1); x by Matheron's rule,
   x = x0 + sigma^2 K A' S^-1 (d - b0 - B dp - A x0 - n0), x0 from the prior (exact AR(1)), n0 from the noise, S = C + sigma^2 A K A'.
3. Weight w = exp(-(chi2_exact - chi2_linear)/2): chi2_exact is the non-linear model of exacto_sn.py at (p0 + dp, x), no refit,
   M explicit; chi2_linear the linear model at the same point, same data, same C, same treatment of (h, omega_m, M).
   H^2 <= 0 -> w = 0.
4. ln Z_exact(sigma) = ln Z_linear(sigma) + ln <w>; ln Z_linear from the analytic marginal (flat priors on dp). Bounds: drop of
   1.92 and 1.35 from the maximum over the grid including sigma = 0, interpolated in ln sigma. ESS = (sum w)^2 / sum w^2.
   Bootstrap error (1000 resamples of the weights).
Checks: at sigma = 0.01 the weights must be ~ 1; ln Z_linear must reproduce the analytic linear bound of EXACTO_SN.md; the
samples must fit the data (mean chi2_linear against chi2_linear(E[m]) + tr(C^-1 Cov_post(m))); with X = 0 and dp = 0 the two
chi2 coincide (difference reported).

    ../venv/bin/python marginal_guiada.py     # writes marginal_guiada.{json,txt}
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.linalg import cho_factor, cho_solve, solve_triangular

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import exacto_sn as S  # noqa: E402
import precision_sn as PS  # noqa: E402

E, X = S.E, S.X
G = S.G0
SN = S.SN["all"]
NSN = SN["n"]
MAG = S.MAG
OMR, ZSTAR, C_KMS = S.OMR, S.ZSTAR, S.C_KMS
THETAS = (0.3, 1.0, 3.0, 10.0)
SIG_MAIN = (0.03, 0.05, 0.07, 0.1, 0.15, 0.2, 0.3, 0.4)
SIG_CHECK = (0.0, 0.01)
N_SAMP, N_BOOT, SEED = 5000, 1000, 20261115
NROWS = 15 + NSN

# ----------------------------------------------------------------------------------------------------------------------
# data, covariance, best fit, exact model (vectorised over samples)
# ----------------------------------------------------------------------------------------------------------------------
D = np.concatenate([X.Y15, SN["mb"]])
C = np.zeros((NROWS, NROWS))
C[:13, :13] = X.COV13
C[13, 13], C[14, 14] = X.SIG_TH**2, X.SIG_WM**2
C[15:, 15:] = np.linalg.inv(SN["Ci"])
LC = np.linalg.cholesky(C)
CF = cho_factor(C, lower=True)
# SN interpolation indices on the grid
_k = np.searchsorted(G.U, SN["uhd"]) - 1
_t = (SN["uhd"] - G.U[_k]) / (G.U[_k + 1] - G.U[_k])


def exact_model(h, wm, M, dX):
    """Exact model vector(s) for parameter arrays h, wm, M (n,) and paths dX (n, N) [dX = X(z) - X(0)]; rows with H^2 <= 0 -> NaN."""
    h, wm, M = np.atleast_1d(h), np.atleast_1d(wm), np.atleast_1d(M)
    om, orad = (wm / h**2)[:, None], (OMR / h**2)[:, None]
    e2 = om * G.p3 + orad * G.p4 + (1.0 - om - orad) + dX
    bad = e2.min(axis=1) <= 0.0
    e2 = np.where(e2 > 0, e2, 1.0)
    inv = 1.0 / np.sqrt(e2)
    f = G.p1 * inv
    dm_u = np.concatenate([np.zeros((len(h), 1)), np.cumsum(0.5 * (f[:, :-1] + f[:, 1:]) * G.DUK, axis=1)], axis=1)
    dm, dh = dm_u[:, G.IDX], inv[:, G.IDX]
    d = np.where(X.IS_DM, dm, np.where(X.IS_DH, dh, np.cbrt(X.Z13 * dm * dm * dh)))
    bao = (C_KMS / (100.0 * h * E.RD_SPL(wm)))[:, None] * d
    dm_star = C_KMS / (100.0 * h) * (dm_u[:, -1] + inv[:, -1] * (np.array([E.zstar_of(w) for w in wm]) - ZSTAR))
    th = 100.0 * E.RS_SPL(wm) / dm_star / E.CAL
    mu = 5.0 * np.log10((1.0 + SN["zhel"])[None, :] * ((1 - _t)[None, :] * dm_u[:, _k] + _t[None, :] * dm_u[:, _k + 1])) + M[:, None]
    out = np.concatenate([bao, th[:, None], wm[:, None], mu], axis=1)
    out[bad] = np.nan
    return out, bad


def chi2_of(models):
    r = D[None, :] - models
    y = solve_triangular(LC, r.T, lower=True, check_finite=False)
    return np.sum(y * y, axis=0)


# best fit (h, omega_m, M) of this data set with X = 0: from exacto_sn.fit_combo (M explicit)
_c0, _h0, _wm0, _M0 = S.fit_combo(G, np.zeros(G.N), SN, True)
P0 = np.array([_h0, _wm0, _M0])
B0 = exact_model(P0[0], P0[1], P0[2], np.zeros((1, G.N)))[0][0]


def build_B():
    B = np.zeros((NROWS, 3))
    for j, eps in enumerate((1e-4, 1e-5, 1e-3)):
        pp, pm = P0.copy(), P0.copy()
        pp[j] += eps
        pm[j] -= eps
        B[:, j] = (exact_model(pp[0], pp[1], pp[2], np.zeros((1, G.N)))[0][0] - exact_model(pm[0], pm[1], pm[2], np.zeros((1, G.N)))[0][0]) / (2 * eps)
    return B


def build_A():
    """Full linear response (convention B, Lambda absorbs X(0)) on the exacto grid at p0, raw units."""
    h, wm = P0[0], P0[1]
    om, orad = wm / h**2, OMR / h**2
    E2 = om * G.p3 + orad * G.p4 + (1 - om - orad)
    Es = np.sqrt(E2)
    N = G.N
    wint = np.zeros(N)
    wint[1:] += 0.5 * G.DUK
    wint[:-1] += 0.5 * G.DUK  # interior trapezoid weights (full)
    g = G.p1 / (2 * Es**3)
    f = G.p1 / Es
    dm = np.concatenate([[0.0], np.cumsum(0.5 * (f[:-1] + f[1:]) * G.DUK)])
    # Wdm[k, j]: weight of X(z_j) in delta ln D_M(z_k): -tw(k)_j g_j / dm_k, tw(k)_j = wint_j for j < k, 0.5 DUK_{k-1} for j = k
    T = np.tril(np.tile(wint, (N, 1)), -1)
    T[np.arange(1, N), np.arange(1, N)] = 0.5 * G.DUK
    T[0, 0] = 0.0
    Wdm = -(T * g[None, :]) / np.where(dm > 0, dm, 1.0)[:, None]
    Wdm[:, 0] -= Wdm.sum(axis=1)  # Lambda absorbs X(0): rows sum to zero
    WdH = np.zeros((len(G.IDX), N))
    for i, k in enumerate(G.IDX):
        WdH[i, k] = -1.0 / (2 * E2[k])
        WdH[i, 0] += 1.0 / (2 * E2[k])
    rows = []
    for i, (k, kind) in enumerate(zip(G.IDX, X.KINDS13)):
        if kind == "DM":
            rows.append(Wdm[k])
        elif kind == "DH":
            rows.append(WdH[i])
        else:
            rows.append((2 * Wdm[k] + WdH[i]) / 3)
    A = np.zeros((NROWS, N))
    A[:13] = np.array(rows) * B0[:13, None]
    A[13] = -Wdm[-1] * B0[13]
    A[15:] = MAG * ((1 - _t)[:, None] * Wdm[_k] + _t[:, None] * Wdm[_k + 1])
    return A


def kernel(theta):
    return np.exp(-theta * np.abs(G.U[:, None] - G.U[None, :]))


def lnZ_lin(Sf, Bm, r):
    """Analytic linear marginal with flat priors on dp (constants dropped). Sf: cho_factor of S."""
    Sr, SB = cho_solve(Sf, r), cho_solve(Sf, Bm)
    F = Bm.T @ SB
    b = Bm.T @ Sr
    logdetS = 2 * np.sum(np.log(np.diag(Sf[0])))
    return float(-0.5 * (r @ Sr - b @ np.linalg.solve(F, b) + logdetS + np.linalg.slogdet(F)[1]))


def cell(theta, sigma, K, AK, Msig, A, Bm, rng):
    """One (theta, sigma) cell: samples, weights, statistics."""
    r = D - B0
    Sm = C + sigma**2 * Msig
    Sf = cho_factor(Sm, lower=True)
    lz = lnZ_lin(Sf, Bm, r)
    SB = cho_solve(Sf, Bm)
    Sp = np.linalg.inv(Bm.T @ SB)
    dp_hat = Sp @ (Bm.T @ cho_solve(Sf, r))
    dp = dp_hat[None, :] + rng.standard_normal((N_SAMP, 3)) @ np.linalg.cholesky(Sp).T
    if sigma > 0:
        x0 = S.ou_paths(G, sigma, theta, N_SAMP, rng)
    else:
        x0 = np.zeros((N_SAMP, G.N))
    n0 = rng.standard_normal((N_SAMP, NROWS)) @ LC.T
    R = r[None, :] - dp @ Bm.T - x0 @ A.T - n0
    if sigma > 0:
        x = x0 + sigma**2 * cho_solve(Sf, R.T).T @ AK  # sigma^2 K A' S^-1 R  (AK = A K, symmetric K)
    else:
        x = x0
    # linear chi2 at (dp, x)
    m_lin = B0[None, :] + dp @ Bm.T + x @ A.T
    chi_lin = chi2_of(m_lin)
    # exact chi2 at (p0 + dp, x): dX = x - x(0)
    p = P0[None, :] + dp
    m_ex, bad = exact_model(p[:, 0], p[:, 1], p[:, 2], x - x[:, [0]])
    chi_ex = np.where(bad, np.inf, chi2_of(np.where(bad[:, None], 0.0, m_ex)))
    lw = -0.5 * (chi_ex - chi_lin)
    w = np.where(bad, 0.0, np.exp(np.minimum(lw, 700)))
    # expected linear chi2: chi2(E[m]) + tr(C^-1 Cov_post(m))
    P_ = sigma**2 * cho_solve(Sf, Msig).T  # M S^-1  (Msig = A K A')
    Em = Bm @ dp_hat + P_ @ (r - Bm @ dp_hat)
    IP = np.eye(NROWS) - P_
    cov_m = IP @ (Bm @ Sp @ Bm.T) @ IP.T + sigma**2 * Msig - P_ @ (sigma**2 * Msig)
    y = solve_triangular(LC, (r - Em), lower=True, check_finite=False)
    n_eff = float(np.trace(cho_solve(CF, cov_m)))
    exp_chi = float(y @ y + n_eff)
    # bootstrap of ln <w>
    idx = rng.integers(0, N_SAMP, size=(N_BOOT, N_SAMP))
    lnw_b = np.log(np.maximum(w[idx].mean(axis=1), 1e-300))
    return {"sigma": sigma, "lnZ_lin": lz, "ln_mean_w": float(np.log(max(w.mean(), 1e-300))), "ln_mean_w_boot": lnw_b.tolist(),
            "ess": float(w.sum() ** 2 / np.sum(w**2)) if w.sum() > 0 else 0.0, "w_max_frac": float(w.max() / w.sum()) if w.sum() > 0 else 1.0,
            "frac_H2_le_0": float(bad.mean()), "chi_lin_mean": float(chi_lin.mean()), "chi_lin_expected": exp_chi, "n_eff": n_eff,
            "chi_ex_mean_finite": float(np.mean(chi_ex[np.isfinite(chi_ex)])), "lw_p16_p50_p84": [float(v) for v in np.percentile(lw[np.isfinite(lw)], [16, 50, 84])],
            "dp_hat": dp_hat.tolist(), "dp_std": np.sqrt(np.diag(Sp)).tolist()}


def bounds_from(sig, lz_rel, drop):
    return PS.lnl_bound(np.array(sig), np.array(lz_rel), drop)


def main():
    t0 = time.time()
    Bm, A = build_B(), build_A()
    rowsum = float(np.max(np.abs(A.sum(axis=1))))
    chi_ex0 = float(chi2_of(B0[None, :])[0])
    chi_lin0 = float(chi2_of(B0[None, :])[0])  # same point: identical by construction
    # independent check of exact vs linear at (dp = 0, X = 0): the exact model recomputed
    m0, _ = exact_model(P0[0], P0[1], P0[2], np.zeros((1, G.N)))
    checks = {"p0": {"h": P0[0], "omega_m": P0[1], "M": P0[2], "chi2_exact_fit": _c0}, "chi2_exact_at_p0_X0": float(chi2_of(m0)[0]),
              "chi2_linear_at_p0_X0": chi_lin0, "difference": float(chi2_of(m0)[0] - chi_lin0), "A_rows_sum_max_abs": rowsum,
              "n_rows": NROWS}
    print(json.dumps(checks, indent=1), flush=True)
    lin_ref = json.loads((HERE / "exacto_sn.json").read_text())["linear"]["BAO+CMB+SN"]
    p5_ref = json.loads((HERE / "convergencia_sn.json").read_text())["thetas"]
    out = {"n_samples": N_SAMP, "n_boot": N_BOOT, "seed": SEED, "checks": checks, "thetas": {}}
    rng = np.random.default_rng(SEED)
    for theta in THETAS:
        t1 = time.time()
        K = kernel(theta)
        AK = A @ K
        Msig = AK @ A.T
        Msig = 0.5 * (Msig + Msig.T)
        sigs = SIG_CHECK + SIG_MAIN + ((0.5, 0.7) if theta == 0.3 else ())
        cells = [cell(theta, s, K, AK, Msig, A, Bm, rng) for s in sigs]
        sig = [c["sigma"] for c in cells]
        lz_lin = np.array([c["lnZ_lin"] for c in cells])
        lnw = np.array([c["ln_mean_w"] for c in cells])
        lz_lin_rel = lz_lin - lz_lin[0]
        lz_ex_rel = (lz_lin + lnw) - (lz_lin[0] + lnw[0])
        sig_pos = np.array(sig[1:])
        rec = {"sigmas": sig, "cells": cells, "lnZ_lin_rel": lz_lin_rel.tolist(), "lnZ_exact_rel": lz_ex_rel.tolist(),
               "correction_ln_Zex_over_Zlin": (lnw - lnw[0]).tolist()}
        for name, drop in (("192", 1.92), ("135", 1.35)):
            rec[f"lin_{name}"] = bounds_from(sig_pos, lz_lin_rel[1:], drop)
            rec[f"exact_{name}"] = bounds_from(sig_pos, lz_ex_rel[1:], drop)
            bb = []
            for b in range(N_BOOT):
                lnw_b = np.array([c["ln_mean_w_boot"][b] for c in cells])
                bb.append(bounds_from(sig_pos, ((lz_lin + lnw_b) - (lz_lin[0] + lnw_b[0]))[1:], drop))
            bb = np.array(bb)
            fin = np.isfinite(bb) & (bb > 0)
            rec[f"exact_{name}_boot_p16_p84"] = [float(np.percentile(bb[fin], 16)), float(np.percentile(bb[fin], 84))] if fin.sum() > 10 else None
            rec[f"exact_{name}_boot_frac_no_bound"] = float(np.mean(~np.isfinite(bb)))
        rec["lnZ_exact_max_at_sigma"] = float(np.array(sig)[int(np.argmax(lz_ex_rel))])
        rec["lnZ_exact_max"] = float(lz_ex_rel.max())
        rec["lnZ_lin_max_at_sigma"] = float(np.array(sig)[int(np.argmax(lz_lin_rel))])
        for key in ("lin_192", "exact_192"):
            b = rec[key]
            if np.isfinite(b) and b > 0:
                k = int(np.argmin(np.abs(np.log(sig_pos / b)))) + 1
                rec[f"at_{key}"] = {"sigma_cell": sig[k], "ess": cells[k]["ess"], "ess_frac_N": cells[k]["ess"] / N_SAMP, "correction": float(lnw[k] - lnw[0]),
                                    "w_max_frac": cells[k]["w_max_frac"], "frac_H2_le_0": cells[k]["frac_H2_le_0"]}
        rec["reference"] = {"linear_exacto_sn_192": lin_ref[f"theta_{theta:g}"]["192"], "linear_exacto_sn_135": lin_ref[f"theta_{theta:g}"]["135"],
                            "p5_N20000": p5_ref[f"theta_{theta:g}"]["by_n"]["20000"]["p5"],
                            "mc_direct_192_N20000": p5_ref[f"theta_{theta:g}"]["by_n"]["20000"]["marg_192"],
                            "mc_direct_ess": p5_ref[f"theta_{theta:g}"]["by_n"]["20000"].get("ess_at_marg_192", {}).get("ess")}
        rec["seconds"] = time.time() - t1
        out["thetas"][f"theta_{theta:g}"] = rec
        print(f"theta {theta:g} [{rec['seconds']:.0f} s]: lin192={PS.fmt(rec['lin_192'])} (exacto_sn {lin_ref[f'theta_{theta:g}']['192']:.3f}) exact192={PS.fmt(rec['exact_192'])} "
              f"{rec['exact_192_boot_p16_p84']} exact135={PS.fmt(rec['exact_135'])} max@{rec['lnZ_exact_max_at_sigma']:g} ({rec['lnZ_exact_max']:+.2f}) "
              f"at_exact={rec.get('at_exact_192')}", flush=True)
        for c in cells:
            print(f"   sigma={c['sigma']:<5g} lnZlin={c['lnZ_lin'] - lz_lin[0]:+8.3f} ln<w>={c['ln_mean_w']:+7.3f} ess={c['ess']:7.0f} wmax={c['w_max_frac']:.3f} H2<=0={c['frac_H2_le_0']:.4f} "
                  f"chi_lin mean={c['chi_lin_mean']:.1f} exp={c['chi_lin_expected']:.1f} (n_eff {c['n_eff']:.1f}) lw p16/50/84={[round(v, 2) for v in c['lw_p16_p50_p84']]}", flush=True)
    out["runtime_s"] = time.time() - t0
    (HERE / "marginal_guiada.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    L = [f"BAO+CMB+SN, N = {N_SAMP} muestras por celda; checks: {json.dumps(checks)}", "",
         "theta | 5 % (N=20000) | lineal analitica (exacto_sn) | lineal (este A, ln Z_lin) | exacta guiada -1.92 (16-84) | -1.35 | ESS en la cota | correccion ln(Zex/Zlin) en la cota | MC directo -1.92 (N=20000, ESS)"]
    for theta in THETAS:
        r = out["thetas"][f"theta_{theta:g}"]; ref = r["reference"]; a = r.get("at_exact_192", {})
        L.append(f"{theta:5g} | {PS.fmt(ref['p5_N20000'])} | {PS.fmt(ref['linear_exacto_sn_192'])} | {PS.fmt(r['lin_192'])} | {PS.fmt(r['exact_192'])} {r['exact_192_boot_p16_p84']} | "
                 f"{PS.fmt(r['exact_135'])} | {a.get('ess', float('nan')):.0f} ({100 * a.get('ess_frac_N', float('nan')):.0f} %) | {a.get('correction', float('nan')):+.2f} | {PS.fmt(ref['mc_direct_192_N20000'])} ({ref['mc_direct_ess']})")
    (HERE / "marginal_guiada.txt").write_text("\n".join(L) + "\n", encoding="utf-8")
    print("\n".join(L))


if __name__ == "__main__":
    main()
