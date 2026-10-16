#!/usr/bin/env python3
"""
r_d computed from each realisation's own expansion history (rdcalc.py), CMB proxy on a window, and
the final criterion. Writes resultados/rd_real.json.

Model 2 (Zwane, alpha = 0.8824, mu = 0.9804): the 20000 realisations of bao_sim.py (same seed, same
401 points) are extended to lambda in [-ln 1e5, -ln 1000] (step 0.1) plus the 10 window points
z in [1000, 5000] (log spaced). The new points are sampled conditioned on the existing lambda_rec =
-ln 1101 (the only correlated one: every other existing point is > 5.6 e-folds away, correlation
< 2e-7): mean k Omega_hat(lambda_rec), covariance K - k k^T, new normals from seed 20261017.
H^2 = 100^2 [w_m (1+z)^3 + w_r (1+z)^4] / (1 - Omega_de(z)); r_d = int_{z_d}^{1e5} c_s/H dz (trapezoid on
4000 log-z points, Omega_hat interpolated linearly in lambda) + standard tail z > 1e5 (Omega_de = 0
beyond 1e5, 2.65 Mpc of ~147), times the CAMB normalisation of rdcalc.

Model 1 (Das, n = 1e5): r_d accumulated inside each run for 25 omega_m values (bao_sim.py m1), for
z < 9999 (start a = 1e-4); for z > 9999 the standard integrand (Lambda = 0, consistent with S_0 = 0
at the start), which is 24.4 Mpc of 146.9 (16.6 %). Proxy on the stored Omega_de at the 25 log points
of z in [1000, 1e4] restricted to [1000, 5000].

Fractions (A: chi2 <= 29.84; B: chi2 - 29.84 < 4; A': chi2 <= 10.28):
  (a) r_d per realisation, omega_m fixed (Planck)
  (b) r_d per realisation, omega_m free in [0.10, 0.20]
  (c) r_d per realisation, window proxy, omega_m within +-1.6 % (2 sigma Planck)   <- final criterion
Low-z distances scale as 1/sqrt(omega_m) (radiation not rescaled; error <= 0.13 in chi2 for (b) in
Model 2, robustez.json; negligible within +-1.6 %).
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline

import bao_sim as B
import rdcalc as RC
import robustez as RB

HERE = Path(__file__).resolve().parent
OUT = HERE / "resultados"
OM0 = RC.OM0
REF_AP = json.loads((OUT / "robustez.json").read_text())["lcdm_desi"]["chi2_min"]
WM_RD = np.unique(np.append(np.geomspace(0.10, 0.20, 24), OM0))
Z_WIN10 = np.geomspace(1000.0, 5000.0, 10)
WM_FREE = np.linspace(0.10, 0.20, 801)
WM_C = OM0 * np.linspace(0.984, 1.016, 161)


def fractions(chi: np.ndarray, n_total: int, mask=None) -> dict:
    ok = np.isfinite(chi) if mask is None else (np.isfinite(chi) & mask)
    A = ok & (chi <= B.CHI2_LCDM)
    Bc = ok & (chi - B.CHI2_LCDM < 4)
    Ap = ok & (chi <= REF_AP)
    out = {"n_total": n_total, "nA": int(A.sum()), "fA": A.sum() / n_total, "nB": int(Bc.sum()),
           "fB": Bc.sum() / n_total, "nA_prime": int(Ap.sum()), "fA_prime": Ap.sum() / n_total}
    if ok.any():
        out["chi2_min"] = float(np.min(chi[ok]))
    return out


def minimise(t0: np.ndarray, rd: np.ndarray, wm_grid: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """t0: (13, n) D/r_d at omega_m0 with r_d = rd(omega_m0); rd: (n, 25) on WM_RD (NaN allowed).
    Returns chi2 min over wm_grid and the argmin omega_m."""
    m = B.BAO["mean"][:, None]
    a = np.einsum("in,ij,jn->n", t0, B.CINV, t0)
    b = np.einsum("in,ij,jn->n", t0, B.CINV, np.broadcast_to(m, t0.shape))
    c = float(B.BAO["mean"] @ B.CINV @ B.BAO["mean"])
    lr = np.log(rd)
    good = np.all(np.isfinite(lr), axis=1)
    rd_g = np.full((rd.shape[0], len(wm_grid)), np.nan)
    if good.any():
        cs = CubicSpline(np.log(WM_RD), lr[good], axis=1)
        rd_g[good] = np.exp(cs(np.log(wm_grid)))
    # rows with NaN (H^2 < 0 for some omega_m): linear interpolation over the finite nodes only
    for j in np.nonzero(~good)[0]:
        f = np.isfinite(lr[j])
        if f.sum() >= 2:
            inside = (wm_grid >= WM_RD[f].min()) & (wm_grid <= WM_RD[f].max())
            rd_g[j, inside] = np.exp(np.interp(np.log(wm_grid[inside]), np.log(WM_RD[f]), lr[j, f]))
    i0 = int(np.argmin(np.abs(WM_RD - OM0)))
    s = np.sqrt(OM0 / wm_grid)[None, :] * rd[:, [i0]] / rd_g
    chi = s * s * a[:, None] - 2 * s * b[:, None] + c
    chi = np.where(np.isfinite(chi), chi, np.inf)
    k = np.argmin(chi, axis=1)
    return chi[np.arange(len(k)), k], wm_grid[k]


def int_from_zd(zq: np.ndarray, g: np.ndarray, zdk: float) -> np.ndarray:
    """Trapezoid of g (len(zq), ...) from zdk to zq[-1]; the first partial interval is interpolated."""
    i = int(np.searchsorted(zq, zdk))  # zq[i-1] < zdk <= zq[i]
    w = (zdk - zq[i - 1]) / (zq[i] - zq[i - 1])
    g_zd = (1 - w) * g[i - 1] + w * g[i]
    zz, gg = zq[i:], g[i:]
    body = np.sum(0.5 * (gg[1:] + gg[:-1]) * np.diff(zz).reshape((-1,) + (1,) * (g.ndim - 1)), axis=0)
    return body + 0.5 * (g_zd + gg[0]) * (zz[0] - zdk)


def model2() -> dict:
    mu, alpha = 0.9804, 0.8824
    L = B.chol_m2(mu)
    rng = np.random.default_rng(B.SEED)
    Z = rng.standard_normal((401, 20000))
    Ohat = alpha * L @ Z
    lam_rec = B.LAM_ALL[400]
    lam_new = np.unique(np.concatenate([np.arange(-np.log(1e5), -np.log(1000.0) + 1e-9, 0.1), -np.log1p(Z_WIN10)]))
    k = np.exp(-mu * (lam_new - lam_rec) ** 2 / 2)
    Knn = np.exp(-mu * (lam_new[:, None] - lam_new[None, :]) ** 2 / 2) - np.outer(k, k)
    w, V = np.linalg.eigh(Knn)
    Ls = V * np.sqrt(np.clip(w, 0, None))[None, :]
    Zn = np.random.default_rng(B.SEED + 1).standard_normal((len(lam_new), 20000))
    Onew = alpha * (k[:, None] * (Ohat[400][None, :] / alpha) + Ls @ Zn)
    out = {"lambda_new": {"n": int(len(lam_new)), "min": float(lam_new.min()), "max": float(lam_new.max())},
           "cond_check_max_corr_other_points": float(np.exp(-mu * (lam_new.max() - B.LAM[0]) ** 2 / 2))}

    # low-z distances at omega_m0 (same code as robustez.py)
    Ode = np.tanh(Ohat)
    zg = np.exp(-B.LAM) - 1.0
    H = 100.0 * np.sqrt((OM0 * (1 + zg[:, None]) ** 3 + B.ORAD * (1 + zg[:, None]) ** 4) / (1.0 - Ode[:400]))
    f = B.C * np.exp(-B.LAM)[:, None] / H
    cum = np.concatenate([np.zeros((1, 20000)), np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(B.LAM)[:, None], axis=0)])
    DMg = cum[-1][None, :] - cum
    zu = np.unique(B.ZR)
    lu = -np.log1p(zu)
    kk = np.clip(np.searchsorted(B.LAM, lu) - 1, 0, len(B.LAM) - 2)
    ww = ((lu - B.LAM[kk]) / (B.LAM[kk + 1] - B.LAM[kk]))[:, None]
    DMu = (1 - ww) * DMg[kk] + ww * DMg[kk + 1]
    Hu = (1 - ww) * H[kk] + ww * H[kk + 1]
    H0 = H[-1]
    t_rd0 = RB.theory13((DMu * H0[None, :] / B.C).T, (H0[None, :] / Hu).T, H0, zu)  # D / RD_CAMB

    # r_d per realisation for WM_RD
    zq = np.geomspace(RC.zd(WM_RD).min() * 0.999, 1e5, 4000)
    lq = -np.log1p(zq)
    order = np.argsort(lam_new)
    Oq = np.empty((len(zq), 20000))
    for j0 in range(0, 20000, 2000):
        sl = slice(j0, j0 + 2000)
        Oq[:, sl] = np.array([np.interp(lq, lam_new[order], Onew[order, j]) for j in range(j0, min(j0 + 2000, 20000))]).T
    Odq = np.tanh(Oq)
    rd = np.empty((20000, len(WM_RD)))
    for i, wm in enumerate(WM_RD):
        Hq = 100.0 * np.sqrt((wm * (1 + zq[:, None]) ** 3 + B.ORAD * (1 + zq[:, None]) ** 4) / (1.0 - Odq))
        g = RC.cs(zq)[:, None] / Hq
        rd[:, i] = RC.NORM * (int_from_zd(zq, g, float(RC.zd(wm))) + RC.rd_standard(wm, 1e5))
    # validation: Omega_de = 0 through the same quadrature
    g0 = RC.cs(zq) / (100.0 * np.sqrt(OM0 * (1 + zq) ** 3 + B.ORAD * (1 + zq) ** 4))
    v0 = RC.NORM * (int_from_zd(zq, g0, float(RC.zd(OM0))) + RC.rd_standard(OM0, 1e5))
    out["check_rd_no_de"] = {"rd_quadrature": float(v0), "rd_camb": RC.RD_CAMB, "rel": float(v0 / RC.RD_CAMB - 1),
                             "rd_unnormalised_vs_aubourg": float(v0 / RC.NORM / 147.06303008488774 - 1)}

    i0 = int(np.argmin(np.abs(WM_RD - OM0)))
    t0 = t_rd0 * (RC.RD_CAMB / rd[:, i0])[None, :]
    chi_a = np.einsum("in,ij,jn->n", t0 - B.BAO["mean"][:, None], B.CINV, t0 - B.BAO["mean"][:, None])
    chi_b, wm_b = minimise(t0, rd, WM_FREE)
    chi_c, wm_c = minimise(t0, rd, WM_C)
    lw = -np.log1p(Z_WIN10)
    iw = np.array([int(np.argmin(np.abs(lam_new - l))) for l in lw])
    proxy_win = np.all(np.abs(np.tanh(Onew[iw])) < 0.25, axis=0)
    proxy_rec = np.abs(np.tanh(Ohat[400])) < 0.25
    out["rd_stats"] = {"p16_p50_p84_at_wm0": [float(v) for v in np.percentile(rd[:, i0], [16, 50, 84])]}
    out["proxy"] = {"rec_only": float(proxy_rec.mean()), "window_1000_5000": float(proxy_win.mean()),
                    "n_window": int(proxy_win.sum())}
    out["a_rd_real_wm_fixed"] = fractions(chi_a, 20000)
    out["a_rd_real_wm_fixed_window"] = fractions(chi_a, 20000, proxy_win)
    out["b_rd_real_wm_free"] = fractions(chi_b, 20000)
    out["b_rd_real_wm_free_window"] = fractions(chi_b, 20000, proxy_win)
    out["c_rd_real_window_wm_2sigma"] = fractions(chi_c, 20000, proxy_win)
    out["c_without_window"] = fractions(chi_c, 20000)
    jb = int(np.argmin(np.where(proxy_win, chi_c, np.inf)))
    out["c_best"] = {"index": jb, "chi2": float(chi_c[jb]), "omega_m": float(wm_c[jb]), "rd": float(rd[jb, i0]),
                     "Ode0": float(Ode[399, jb])}
    jb = int(np.argmin(chi_b))
    out["b_best"] = {"index": jb, "chi2": float(chi_b[jb]), "omega_m": float(wm_b[jb]), "window": bool(proxy_win[jb])}
    np.savez_compressed(OUT / "rd_real_m2.npz", rd=rd, wm_rd=WM_RD, chi_a=chi_a, chi_b=chi_b, wm_b=wm_b,
                        chi_c=chi_c, wm_c=wm_c, proxy_win=proxy_win)
    return out


def model1(alpha: float) -> dict:
    R = np.load(OUT / f"bao_m1_alpha{alpha}.npz")
    if "rd_part" not in R:
        return {"note": "npz without rd_part"}
    assert np.allclose(R["wm_rd"], WM_RD)
    tail = np.array([RC.rd_standard(w, 1.0 / 1e-4 - 1.0) for w in WM_RD])
    rd = RC.NORM * (R["rd_part"] + tail[None, :])
    t_rd0 = RB.theory13(R["dmh"], R["eh"], R["H0"], R["z_unique"])
    i0 = int(np.argmin(np.abs(WM_RD - OM0)))
    t0 = t_rd0 * (RC.RD_CAMB / rd[:, i0])[None, :]
    chi_a = np.einsum("in,ij,jn->n", t0 - B.BAO["mean"][:, None], B.CINV, t0 - B.BAO["mean"][:, None])
    chi_b = np.empty(len(chi_a)); wm_b = np.empty(len(chi_a))
    chi_c = np.empty(len(chi_a)); wm_c = np.empty(len(chi_a))
    for j0 in range(0, len(chi_a), 20000):
        sl = slice(j0, j0 + 20000)
        chi_b[sl], wm_b[sl] = minimise(t0[:, sl], rd[sl], WM_FREE)
        chi_c[sl], wm_c[sl] = minimise(t0[:, sl], rd[sl], WM_C)
    zw = R["z_win"]
    inwin = (zw >= 1000.0) & (zw <= 5000.0)
    proxy_win = np.all(np.abs(R["win_ode"][:, inwin]) < 0.25, axis=1)
    n_seeds = json.loads((OUT / f"bao_m1_alpha{alpha}.json").read_text())["n_seeds"]
    out = {"n_seeds": n_seeds, "n_complete": int(len(chi_a)), "n_window_points": int(inwin.sum()),
           "tail_z_gt_9999_Mpc_at_wm0": float(RC.NORM * tail[i0]),
           "rd_p16_p50_p84_at_wm0": [float(v) for v in np.nanpercentile(rd[:, i0], [16, 50, 84])],
           "frac_rd_nan_any_wm": float(np.mean(~np.all(np.isfinite(rd), axis=1))),
           "proxy_window_frac_of_complete": float(proxy_win.mean()), "n_window": int(proxy_win.sum())}
    out["a_rd_real_wm_fixed"] = fractions(chi_a, n_seeds)
    out["a_rd_real_wm_fixed_window"] = fractions(chi_a, n_seeds, proxy_win)
    out["b_rd_real_wm_free"] = fractions(chi_b, n_seeds)
    out["b_rd_real_wm_free_window"] = fractions(chi_b, n_seeds, proxy_win)
    out["c_rd_real_window_wm_2sigma"] = fractions(chi_c, n_seeds, proxy_win)
    out["c_without_window"] = fractions(chi_c, n_seeds)
    jb = int(np.argmin(chi_b))
    out["b_best"] = {"chi2": float(chi_b[jb]), "omega_m": float(wm_b[jb]), "window": bool(proxy_win[jb]),
                     "H0": float(R["H0"][jb] * np.sqrt(wm_b[jb] / OM0))}
    jc = int(np.argmin(np.where(proxy_win, chi_c, np.inf)))
    out["c_best"] = {"chi2": float(chi_c[jc]), "omega_m": float(wm_c[jc]),
                     "H0": float(R["H0"][jc] * np.sqrt(wm_c[jc] / OM0)), "OmL0": float(R["OmL0"][jc])}
    np.savez_compressed(OUT / f"rd_real_m1_alpha{alpha}.npz", rd=rd, chi_a=chi_a, chi_b=chi_b, wm_b=wm_b,
                        chi_c=chi_c, wm_c=wm_c, proxy_win=proxy_win)
    return out


def main() -> None:
    import sys

    p = OUT / "rd_real.json"
    res = json.loads(p.read_text()) if p.exists() else {}
    res["criteria"] = {"A": B.CHI2_LCDM, "B": B.CHI2_LCDM + 4, "A_prime": REF_AP,
                       "wm_free": [0.10, 0.20], "wm_c": [float(WM_C[0]), float(WM_C[-1])],
                       "rd_norm_camb": RC.NORM, "zd_camb": RC.ZD_CAMB}
    which = sys.argv[1:] or ["m2", "m1"]
    if "m2" in which:
        res["model2"] = model2()
        print(json.dumps(res["model2"], indent=1))
    if "m1" in which:
        for a in (0.01, 0.005):
            res[f"model1_alpha{a}"] = model1(a)
            print(a, json.dumps(res[f"model1_alpha{a}"], indent=1))
    p.write_text(json.dumps(res, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
