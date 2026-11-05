#!/usr/bin/env python3
"""
Exploration only (CUESTA.md): false positives of w0-wa produced by an OU residual, BAO only.

Base universe: flat LCDM at the DESI DR2 best fit (Omega_m,f = 0.29740, h r_d,f = 101.538 Mpc, FONDO.md), predicting the
13 DESI DR2 BAO numbers (D_V/r_d, D_M/r_d, D_H/r_d; desi_dr2_data.load_gaussian_bao_vector). For each (sigma_X, theta) of
the grid theta in {0.1, 0.3, 1, 3} x sigma_X in {0, 95% bound with the background fixed at that theta, 95% bound with the
background marginalised at that theta} (bounds recomputed here with perfil_13's machinery on the same 13 numbers; sigma_X = 0
is the noise-only control, independent of theta), N_UNI universes are generated:
  obs_i = pred_i(ref) * (1 + S_i X(z_i)) + n,  X ~ N(0, sigma_X^2 e^{-theta |x_i - x_j|}), x = ln(1 + z),  n ~ N(0, C_13x13),
with the same kernel as the nucleus (perfil_fondo.py) and S by type (D_V, D_M, D_H) as in perfil_13.py. Rows at the same z share
X. The noise vectors are the same N_UNI draws in every cell (common random numbers), so "excess over the control" is a paired
difference; the OU draws use an independent stream per cell.
Each universe is fitted with flat LCDM (Omega_m, h r_d) and flat w0waCDM (Omega_m, h r_d, w0, wa) by chi^2 minimisation with
the 13x13 covariance: h r_d analytically (the prediction is linear in 1/(h r_d)), Omega_m by bounded Brent (LCDM) and
(Omega_m, w0, wa) by multi-start Nelder-Mead (w0waCDM) with the DESI DR2 priors w0 in [-3, 1], wa in [-3, 2], w0 + wa < 0,
Omega_m in [0.01, 0.99]. Same ingredients as scripts/eos_efectiva.py (CPL E(z), chi^2 with the official covariance,
multi-start Nelder-Mead), but with the 13 numbers and Omega_m, h r_d free instead of the 7 alpha with Omega_m, H0 fixed.
Saved per universe: Delta chi^2 = chi^2(LCDM) - chi^2(w0wa), w0, wa, Omega_m of both fits.
DESI DR2 reference (arXiv:2503.14738, TeX, tab:dchisq and eq:w0wa_DESI): Delta chi^2_MAP = -4.7 (1.7 sigma) for DESI alone;
w0 = -0.48 (+0.35, -0.17), wa < -1.34 (68%), marginalised means, cut by the priors. False positive = universe with
Delta chi^2 >= 4.7 and w0 > -1 and wa < 0.
Checks: the fitter on the real 13 numbers must give the DESI LCDM best fit (0.2975, 101.54 Mpc) and a Delta chi^2 close to 4.7;
the grid distances agree with quad to <= 1e-7.

Writes analisis_cmb/cuesta/cuesta.{json,txt} and cuesta_universos.npz; nothing in results/, figures/ or papers.
"""

from __future__ import annotations

import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize, minimize_scalar

HERE = Path(__file__).resolve().parent
FONDO = HERE.parent / "fondo"
sys.path.insert(0, str(FONDO))
import perfil_13 as P13  # noqa: E402  (imports perfil_fondo: reference background, SIG grid, bound())
import perfil_fondo as PF  # noqa: E402
from desi_dr2_data import load_gaussian_bao_vector  # noqa: E402

OM_F, HRD_F = PF.OM_F, PF.HRD_F
C_KMS = 299792.458
THETAS = (0.1, 0.3, 1.0, 3.0)
N_UNI = 2000
SEED_NOISE, SEED_OU = 20261105, 7_000_000
DCHI2_DESI, SIG_DESI = 4.7, 1.7  # tab:dchisq, DESI alone (minus sign: chi^2(w0wa) - chi^2(LCDM) = -4.7)
W0_DESI, WA_DESI_68 = -0.48, -1.34  # eq:w0wa_DESI
PRIOR = {"om": (0.01, 0.99), "w0": (-3.0, 1.0), "wa": (-3.0, 2.0)}
CHI2_LCDM_OK = 19.675  # 95% point of chi^2 with 11 dof (13 numbers - 2 parameters): "the LCDM fit is still acceptable"
NM_STARTS = ((0.35, -0.5, -1.6), (0.25, -1.2, 0.8), (0.30, -0.75, -1.0), (0.40, -0.3, -2.5))

# ----------------------------------------------------------------------------------------------------------------------
# data
# ----------------------------------------------------------------------------------------------------------------------
G = load_gaussian_bao_vector()
KINDS = [lab[:2] for lab in G["labels"]]
Z = np.asarray(G["z_rows"], float)
OBS, COV = np.asarray(G["mean"], float), np.asarray(G["cov"], float)
CINV = np.linalg.inv(COV)
L_COV = np.linalg.cholesky(COV)
ZU = np.unique(Z)  # 7 effective redshifts
ROW_Z = np.array([int(np.nonzero(ZU == z)[0][0]) for z in Z])  # row -> unique z index
X_U = np.log1p(ZU)
S13 = np.array([P13.s_type(k, z) for k, z in zip(KINDS, Z)])

# Simpson grid: 7 segments [z_{k-1}, z_k] (z_0 = 0), 128 subintervals each
_NSUB = 128
_seg_z, _seg_w = [], []
_edges = np.concatenate([[0.0], ZU])
for a, b in zip(_edges[:-1], _edges[1:]):
    zz = np.linspace(a, b, _NSUB + 1)
    w = np.ones(_NSUB + 1)
    w[1:-1:2], w[2:-1:2] = 4.0, 2.0
    _seg_z.append(zz)
    _seg_w.append(w * (b - a) / _NSUB / 3.0)
ZG = np.concatenate(_seg_z)
WG = np.concatenate(_seg_w)
SEG = np.repeat(np.arange(len(ZU)), _NSUB + 1)
_1PZ = 1.0 + ZG
_1PZ3 = _1PZ**3
_EXPARG = -3.0 * ZG / _1PZ
_1PZU, _1PZU3, _EXPARGU = 1.0 + ZU, (1.0 + ZU) ** 3, -3.0 * ZU / (1.0 + ZU)
IS_DV = np.array([k == "DV" for k in KINDS])
IS_DM = np.array([k == "DM" for k in KINDS])
IS_DH = np.array([k == "DH" for k in KINDS])


def dhat13(om: float, w0: float, wa: float) -> np.ndarray:
    """Dimensionless D_type H0/c at the 13 rows, flat CPL without radiation."""
    p = 3.0 * (1.0 + w0 + wa)
    e_g = np.sqrt(om * _1PZ3 + (1.0 - om) * _1PZ**p * np.exp(wa * _EXPARG))
    dm_u = np.cumsum(np.bincount(SEG, weights=WG / e_g, minlength=len(ZU)))
    dh_u = 1.0 / np.sqrt(om * _1PZU3 + (1.0 - om) * _1PZU**p * np.exp(wa * _EXPARGU))
    dm, dh = dm_u[ROW_Z], dh_u[ROW_Z]
    out = np.where(IS_DM, dm, np.where(IS_DH, dh, np.cbrt(Z * dm * dm * dh)))
    return out


def chi2_profiled(om: float, w0: float, wa: float, y: np.ndarray) -> tuple[float, float]:
    """chi^2 minimised analytically over the scale s = c/(100 h r_d); returns (chi^2, h r_d)."""
    d = dhat13(om, w0, wa)
    cd = CINV @ d
    s = (cd @ y) / (cd @ d)
    r = y - s * d
    return float(r @ CINV @ r), float(C_KMS / (100.0 * s))


def fit_lcdm(y: np.ndarray) -> tuple[float, float, float]:
    f = lambda om: chi2_profiled(om, -1.0, 0.0, y)[0]
    grid = np.linspace(0.05, 0.95, 37)
    j = int(np.argmin([f(o) for o in grid]))
    lo, hi = max(PRIOR["om"][0], grid[max(j - 1, 0)]), min(PRIOR["om"][1], grid[min(j + 1, len(grid) - 1)])
    r = minimize_scalar(f, bounds=(lo, hi), method="bounded", options={"xatol": 1e-7})
    chi2, hrd = chi2_profiled(float(r.x), -1.0, 0.0, y)
    return chi2, float(r.x), hrd


def _obj(p: np.ndarray, y: np.ndarray) -> float:
    om, w0, wa = p
    if not (PRIOR["om"][0] <= om <= PRIOR["om"][1] and PRIOR["w0"][0] <= w0 <= PRIOR["w0"][1]
            and PRIOR["wa"][0] <= wa <= PRIOR["wa"][1]):
        return 1e8
    pen = 0.0 if w0 + wa < 0 else 1e6 * (1.0 + w0 + wa) ** 2  # DESI prior w0 + wa < 0
    return chi2_profiled(om, w0, wa, y)[0] + pen


def fit_w0wa(y: np.ndarray, om_lcdm: float) -> tuple[float, float, float, float, float]:
    best = None
    for p0 in ((om_lcdm, -1.0, 0.0),) + NM_STARTS:
        r = minimize(_obj, np.array(p0), args=(y,), method="Nelder-Mead",
                     options={"xatol": 1e-5, "fatol": 1e-6, "maxiter": 4000, "initial_simplex": None})
        if best is None or r.fun < best.fun:
            best = r
    # polish from the best
    r = minimize(_obj, best.x, args=(y,), method="Nelder-Mead", options={"xatol": 1e-6, "fatol": 1e-8, "maxiter": 4000})
    if r.fun < best.fun:
        best = r
    om, w0, wa = (float(v) for v in best.x)
    chi2, hrd = chi2_profiled(om, w0, wa, y)
    return chi2, om, hrd, w0, wa


def fit_both(y: np.ndarray) -> tuple:
    c_l, om_l, hrd_l = fit_lcdm(y)
    c_w, om_w, hrd_w, w0, wa = fit_w0wa(y, om_l)
    return c_l - c_w, c_l, om_l, hrd_l, c_w, om_w, hrd_w, w0, wa


# ----------------------------------------------------------------------------------------------------------------------
# bounds at the four theta (perfil_13 machinery, 13 numbers), checks
# ----------------------------------------------------------------------------------------------------------------------
def bounds_13() -> tuple[dict, P13.Data]:
    pred = np.array([C_KMS / (100 * HRD_F) * P13.dhat(k, OM_F, z) for k, z in zip(KINDS, Z)])
    D = P13.Data(Z, KINDS, OBS / pred, COV / np.outer(pred, pred), S13)
    out = {}
    for th in THETAS:
        out[f"theta_{th:g}"] = {"a_fixed": PF.bound(lambda s, t=th: D.ll_a(s, t))["sigma_95"],
                                "c_marginalised": PF.bound(lambda s, t=th: D.ll_c(s, t))["sigma_95"]}
    return out, D


def check_distances() -> float:
    worst = 0.0
    for om, w0, wa in ((OM_F, -1.0, 0.0), (0.35, -0.5, -1.6), (0.25, -1.2, 0.8)):
        d = dhat13(om, w0, wa)
        for i, (k, z) in enumerate(zip(KINDS, Z)):
            E = lambda x: np.sqrt(om * (1 + x) ** 3 + (1 - om) * (1 + x) ** (3 * (1 + w0 + wa)) * np.exp(-3 * wa * x / (1 + x)))
            dm = quad(lambda x: 1 / E(x), 0, z, epsrel=1e-12)[0]
            ref = {"DM": dm, "DH": 1 / E(z), "DV": (z * dm * dm / E(z)) ** (1 / 3)}[k]
            worst = max(worst, abs(d[i] / ref - 1))
    return worst


# ----------------------------------------------------------------------------------------------------------------------
# universes
# ----------------------------------------------------------------------------------------------------------------------
PRED_REF = np.array([C_KMS / (100 * HRD_F) * P13.dhat(k, OM_F, z) for k, z in zip(KINDS, Z)])


def universes(sigma: float, theta: float, noise: np.ndarray, seed_ou: int) -> np.ndarray:
    """N_UNI x 13 observed vectors. noise: N_UNI x 13 draws with the official covariance (shared across cells)."""
    if sigma == 0.0:
        return PRED_REF[None, :] + noise
    K = sigma**2 * np.exp(-theta * np.abs(X_U[:, None] - X_U[None, :]))
    L = np.linalg.cholesky(K + 1e-14 * sigma**2 * np.eye(len(ZU)))
    rng = np.random.default_rng(seed_ou)
    X_u = rng.standard_normal((noise.shape[0], len(ZU))) @ L.T
    X = X_u[:, ROW_Z]
    return PRED_REF[None, :] * (1.0 + S13[None, :] * X) + noise


def run_cell(ys: np.ndarray, pool: Pool) -> np.ndarray:
    return np.array(pool.map(fit_both, list(ys), chunksize=25))


def summarise(res: np.ndarray, ctrl: np.ndarray | None) -> dict:
    dchi2, chi2_l, w0, wa = res[:, 0], res[:, 1], res[:, 7], res[:, 8]
    q = (w0 > -1) & (wa < 0)
    fp = (dchi2 >= DCHI2_DESI) & q
    n = len(dchi2)
    d = {"n": int(n), "frac_dchi2_ge_desi": float(np.mean(dchi2 >= DCHI2_DESI)),
         "frac_quadrant": float(np.mean(q)),
         "false_positive": float(fp.mean()), "fp_err": float(np.sqrt(fp.mean() * (1 - fp.mean()) / n)),
         "p_dchi2_ge_desi_given_quadrant": float(np.mean(dchi2[q] >= DCHI2_DESI)),
         "p_dchi2_ge_desi_given_not_quadrant": float(np.mean(dchi2[~q] >= DCHI2_DESI)),
         "frac_wa_at_prior_edge": float(np.mean(wa < PRIOR["wa"][0] + 1e-3)),
         "mean_chi2_lcdm": float(chi2_l.mean()), "median_chi2_lcdm": float(np.median(chi2_l)),
         "frac_chi2_lcdm_ok": float(np.mean(chi2_l <= CHI2_LCDM_OK)),
         "false_positive_chi2_ok": float(np.mean(fp & (chi2_l <= CHI2_LCDM_OK))),
         "median_dchi2": float(np.median(dchi2)), "mean_dchi2": float(np.mean(dchi2)),
         "median_w0_fp": float(np.median(w0[fp])) if fp.any() else None,
         "median_wa_fp": float(np.median(wa[fp])) if fp.any() else None,
         "frac_dchi2_ge_desi_in_quadrant": float(np.mean(fp[(w0 > -1) & (wa < 0)])) if ((w0 > -1) & (wa < 0)).any() else None}
    if ctrl is not None:
        fpc = (ctrl[:, 0] >= DCHI2_DESI) & (ctrl[:, 7] > -1) & (ctrl[:, 8] < 0)
        diff = fp.astype(float) - fpc.astype(float)
        d["excess_over_control"] = float(diff.mean())
        d["excess_err_paired"] = float(diff.std(ddof=1) / np.sqrt(n))
    return d


def main() -> None:
    t0 = time.time()
    dist_err = check_distances()
    real = fit_both(OBS)
    print(f"distancias de malla frente a quad: max dif rel = {dist_err:.1e}", flush=True)
    print(f"datos reales: LCDM Om = {real[2]:.4f}, h r_d = {real[3]:.2f}; w0wa Om = {real[5]:.4f}, h r_d = {real[6]:.2f}, "
          f"w0 = {real[7]:.3f}, wa = {real[8]:.3f}; Delta chi2 = {real[0]:.3f} (DESI: {DCHI2_DESI})", flush=True)
    bnd, _ = bounds_13()
    print("cotas (13 numeros):", json.dumps(bnd), flush=True)

    rng_n = np.random.default_rng(SEED_NOISE)
    noise = rng_n.standard_normal((N_UNI, 13)) @ L_COV.T
    cells = [("control", 0.0, 0.0)]
    for th in THETAS:
        cells.append((f"theta_{th:g}|a_fixed", bnd[f"theta_{th:g}"]["a_fixed"], th))
        cells.append((f"theta_{th:g}|c_marginalised", bnd[f"theta_{th:g}"]["c_marginalised"], th))
    results, raw = {}, {}
    with Pool() as pool:
        for i, (name, sigma, th) in enumerate(cells):
            t1 = time.time()
            ys = universes(sigma, th, noise, SEED_OU + i)
            res = run_cell(ys, pool)
            raw[name] = res
            results[name] = {"sigma_X": sigma, "theta": th, **summarise(res, raw.get("control") if name != "control" else None)}
            print(f"{name:28s} sigma_X = {sigma:.4g}: FP = {results[name]['false_positive']:.4f}  "
                  f"P(dchi2>=4.7) = {results[name]['frac_dchi2_ge_desi']:.4f}  cuadrante = {results[name]['frac_quadrant']:.3f}  "
                  f"[{time.time() - t1:.0f} s]", flush=True)
    np.savez_compressed(HERE / "cuesta_universos.npz", cells=np.array([c[0] for c in cells]),
                        columns=np.array(["dchi2", "chi2_lcdm", "om_lcdm", "hrd_lcdm", "chi2_w0wa", "om_w0wa", "hrd_w0wa", "w0", "wa"]),
                        **{c[0].replace("|", "__"): raw[c[0]] for c in cells})
    out = {"reference_background": {"Omega_m": OM_F, "h_rd_Mpc": HRD_F}, "n_universes": N_UNI, "thetas": THETAS,
           "seeds": {"noise": SEED_NOISE, "ou_base": SEED_OU}, "priors": PRIOR,
           "desi_dr2": {"dchi2_MAP_desi_alone": -DCHI2_DESI, "significance_sigma": SIG_DESI, "dDIC": -0.8,
                        "w0_marg": "-0.48 (+0.35, -0.17)", "wa_marg_68": "< -1.34",
                        "source": "arXiv:2503.14738 TeX: tab:dchisq, eq:w0wa_DESI (DESI BAO alone)"},
           "false_positive_definition": "dchi2 = chi2(LCDM) - chi2(w0wa) >= 4.7 and w0 > -1 and wa < 0",
           "chi2_lcdm_ok_threshold": CHI2_LCDM_OK,
           "checks": {"grid_vs_quad_max_rel": dist_err,
                      "real_data": {"lcdm": {"Omega_m": real[2], "h_rd": real[3], "chi2": real[1]},
                                    "w0wa": {"Omega_m": real[5], "h_rd": real[6], "w0": real[7], "wa": real[8], "chi2": real[4]},
                                    "dchi2": real[0]}},
           "bounds_13": bnd, "cells": results, "runtime_s": time.time() - t0}
    (HERE / "cuesta.json").write_text(json.dumps(out, indent=1), encoding="utf-8")

    ctrl = results["control"]
    lines = [f"universo base: Omega_m = {OM_F:.5f}, h r_d = {HRD_F:.3f} Mpc; {N_UNI} universos por celda; "
             f"falso positivo: dchi2 >= {DCHI2_DESI} y w0 > -1 y wa < 0",
             f"datos reales: dchi2 = {real[0]:.2f} (DESI -{DCHI2_DESI}), w0 = {real[7]:.2f}, wa = {real[8]:.2f}",
             f"control (sigma_X = 0): FP = {ctrl['false_positive']:.4f} +- {ctrl['fp_err']:.4f}, "
             f"P(dchi2 >= 4.7) = {ctrl['frac_dchi2_ge_desi']:.4f}, cuadrante = {ctrl['frac_quadrant']:.3f}, "
             f"<chi2 LCDM> = {ctrl['mean_chi2_lcdm']:.2f}, FP con chi2 LCDM <= {CHI2_LCDM_OK} = {ctrl['false_positive_chi2_ok']:.4f}", "",
             f"{'theta':>6s} {'tratamiento':16s} {'sigma_X':>8s} {'FP':>7s} {'control':>8s} {'exceso':>8s} {'+-':>7s} "
             f"{'P(dchi2>=4.7)':>14s} {'cuadrante':>9s} {'med dchi2':>9s} {'<chi2 LCDM>':>11s} {'FP&chi2ok':>9s}"]
    for th in THETAS:
        for tr in ("a_fixed", "c_marginalised"):
            r = results[f"theta_{th:g}|{tr}"]
            lines.append(f"{th:6g} {tr:16s} {r['sigma_X']:8.4f} {r['false_positive']:7.4f} {ctrl['false_positive']:8.4f} "
                         f"{r['excess_over_control']:+8.4f} {r['excess_err_paired']:7.4f} {r['frac_dchi2_ge_desi']:14.4f} "
                         f"{r['frac_quadrant']:9.3f} {r['median_dchi2']:9.2f} {r['mean_chi2_lcdm']:11.2f} "
                         f"{r['false_positive_chi2_ok']:9.4f}")
    (HERE / "cuesta.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
