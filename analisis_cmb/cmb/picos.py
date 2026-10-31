#!/usr/bin/env python3
"""
Exploration only (PICOS.md): does the OU residual at the DESI ceiling leave a trace in the high-l CMB spectra?

Same set-up as realizaciones.py (section 3.5 of the preprint): path on the fine grid, Gaussian smoothing Delta = 0.05,
2000-point w(a) table in DarkEnergyPPF, Takahashi halofit, H0 solved so that 100 theta_* equals the LCDM value, the
standard normals of seed 20261003. Section 3.5 used 200 realisations; here the first 60 (indices 0-59, the same
normals):
  sigma_X = 2.633e-2 (C2), theta = 1 and theta = 0.001   -> H0 from resultados/takahashi/sigma_2.633e-02/ (checked)
  sigma_X = 4.611e-2 (marginalised ceiling at theta = 1, analisis_cmb/fondo/perfil_fondo.json), theta = 1 -> H0 solved
Reference: LCDM-Planck through the same code (PPF with w = -1, H0 solved: the LCDM reference curve of realizaciones.py).

Statistic (a theory distance, not a fit to data): lensed C_l^TT, TE, EE for l in [30, 2500];
  Var(C_TT) = 2 (C_TT + N_T)^2 / nu,  Var(C_EE) = 2 (C_EE + N_P)^2 / nu,  Var(C_TE) = [C_TE^2 + (C_TT + N_T)(C_EE + N_P)] / nu,
  nu = (2l + 1) f_sky, f_sky = 0.6; per-spectrum Delta chi^2 = sum (Delta C / sigma)^2, and the joint Delta chi^2 with the
  full 3x3 covariance per l, Cov(C_XY, C_WZ) = [C_XW C_YZ + C_XZ C_YW] / nu (Gaussian, cosmic variance + noise).
  Noise: the Planck 2018 I N_l are not in the repository; white noise 33 muK arcmin (temperature), sqrt(2) x 33 for
  polarisation, Gaussian beam FWHM 7 arcmin.
Also the same distance with unlensed spectra, to separate lensing from the primary anisotropies, and the l of the
maximum |Delta C_l / sigma_l| per spectrum.

    OMP_NUM_THREADS=2 ../venv/bin/python picos.py     # writes cmb/picos.json
"""

from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "2")

import json
import multiprocessing as mp
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

N_USE = 60
DELTA = 0.05
LMIN, LMAX_USE = 30, 2500
FSKY = 0.6
NOISE_T = 33.0  # muK arcmin
BEAM = 7.0  # arcmin FWHM
SIG_C2 = 0.0263282756994459
SIG_MARG = json.loads((HERE.parent / "fondo" / "perfil_fondo.json").read_text())["treatments"]["c_marginalised"]["theta_1"]["sigma_95"]
RES_C2 = HERE.parent / "resultados" / "takahashi" / "sigma_2.633e-02"

_TARGET = None


def _init(target):
    global _TARGET
    _TARGET = target


def spectra(task):
    """task = ('cpl', H0) for the LCDM reference or ('table', a, w, H0 or None). Returns lensed/unlensed C_l and H0."""
    import ou_cmb as M
    import camb

    if task[0] == "cpl":
        make = lambda h: M.make_params(w0=-1.0, wa=0.0, H0=h, halofit="takahashi")
        h0 = task[1]
    else:
        make = lambda h: M.make_params(a=task[1], w=task[2], H0=h, halofit="takahashi")
        h0 = task[3]
        if h0 is None:
            h0 = M.solve_H0_for_thetastar(_TARGET, make)[0]
    pars = make(h0)
    res = camb.get_results(pars)
    cl = res.get_cmb_power_spectra(pars, CMB_unit="muK", raw_cl=True)
    L = cl["total"].shape[0]
    return {"H0": float(h0), "lensed": cl["total"][: LMAX_USE + 1, :4], "unlensed": cl["unlensed_total"][: LMAX_USE + 1, :4],
            "thetastar": float(res.get_derived_params()["thetastar"]), "lmax_camb": int(L - 1)}


def noise():
    ell = np.arange(LMAX_USE + 1)
    arc = np.pi / (180 * 60)
    bl2 = np.exp(ell * (ell + 1) * (BEAM * arc) ** 2 / (8 * np.log(2)))
    nt = (NOISE_T * arc) ** 2 * bl2
    return nt, 2 * nt


def distance(ref: np.ndarray, mod: np.ndarray) -> dict:
    """Gaussian Delta chi^2 (TT, TE, EE, joint) of mod against ref (columns TT, EE, BB, TE), l in [LMIN, LMAX_USE]."""
    nt, npol = noise()
    ell = np.arange(LMAX_USE + 1)
    sl = slice(LMIN, LMAX_USE + 1)
    nu = (2 * ell[sl] + 1) * FSKY
    TT, EE, TE = ref[sl, 0] + nt[sl], ref[sl, 1] + npol[sl], ref[sl, 3]
    d = (mod - ref)[sl]
    dTT, dEE, dTE = d[:, 0], d[:, 1], d[:, 3]
    vTT, vEE, vTE = 2 * TT**2 / nu, 2 * EE**2 / nu, (TE**2 + TT * EE) / nu
    zTT, zEE, zTE = dTT / np.sqrt(vTT), dEE / np.sqrt(vEE), dTE / np.sqrt(vTE)
    # joint 3x3 covariance of (TT, EE, TE) per l
    cov = np.empty((len(nu), 3, 3))
    cov[:, 0, 0], cov[:, 1, 1], cov[:, 2, 2] = vTT, vEE, vTE
    cov[:, 0, 1] = cov[:, 1, 0] = 2 * TE**2 / nu
    cov[:, 0, 2] = cov[:, 2, 0] = 2 * TT * TE / nu
    cov[:, 1, 2] = cov[:, 2, 1] = 2 * EE * TE / nu
    v = np.stack([dTT, dEE, dTE], axis=1)
    joint_l = np.einsum("li,li->l", v, np.linalg.solve(cov, v[:, :, None])[:, :, 0])
    ls = ell[sl]
    hi = ls > 1000
    return {"TT": float(np.sum(zTT**2)), "TE": float(np.sum(zTE**2)), "EE": float(np.sum(zEE**2)),
            "joint": float(np.sum(joint_l)), "joint_frac_l_gt_1000": float(np.sum(joint_l[hi]) / np.sum(joint_l)),
            "lmax_TT": int(ls[np.argmax(np.abs(zTT))]), "lmax_TE": int(ls[np.argmax(np.abs(zTE))]),
            "lmax_EE": int(ls[np.argmax(np.abs(zEE))]), "lmax_joint": int(ls[np.argmax(joint_l)])}


def tables(theta: float, sigma: float):
    import ou_cmb as M

    lna_f, tab_idx, main = M.production_grids()
    h = lna_f[1] - lna_f[0]
    z = np.random.default_rng(M.SEED).standard_normal((200, len(lna_f)))[:N_USE]  # the first 60 of section 3.5
    xs = M.smooth_uniform(sigma * M.ou_uniform(z, h, theta, 1.0), h, DELTA)
    _, w, bad = M.rho_w_table(lna_f[tab_idx], main, xs[:, tab_idx], M.omega_l0())
    return np.exp(lna_f[tab_idx][main]), w, bad


def summarise(rows: list[dict]) -> dict:
    out = {}
    for k in ("TT", "TE", "EE", "joint"):
        v = np.array([r[k] for r in rows])
        out[k] = {"p50": float(np.median(v)), "p95": float(np.percentile(v, 95)), "max": float(v.max())}
    for k in ("lmax_TT", "lmax_TE", "lmax_EE", "lmax_joint"):
        v = np.array([r[k] for r in rows])
        out[k] = {"p16_p50_p84": [float(x) for x in np.percentile(v, [16, 50, 84])], "frac_gt_1000": float(np.mean(v > 1000))}
    out["joint_frac_l_gt_1000_median"] = float(np.median([r["joint_frac_l_gt_1000"] for r in rows]))
    return out


def main() -> None:
    import ou_cmb as M

    t0 = time.time()
    target = M.thetastar_background(M.make_params(halofit="takahashi"))
    meta = json.loads((RES_C2 / "meta.json").read_text())
    h0_lcdm = meta["references"]["LCDM"]["H0"]
    cases = [("sigma2.633e-02_theta1", 1.0, SIG_C2, RES_C2 / "observables_theta1_Delta0.05.npz"),
             ("sigma2.633e-02_theta0.001", 0.001, SIG_C2, RES_C2 / "observables_theta0.001_Delta0.05.npz"),
             (f"sigma{SIG_MARG:.3e}_theta1", 1.0, SIG_MARG, None)]
    out = {"n_realisations": N_USE, "realisations": "indices 0-59 of the 200 of section 3.5 (seed 20261003)",
           "delta": DELTA, "l_range": [LMIN, LMAX_USE], "f_sky": FSKY,
           "noise": f"white {NOISE_T} muK arcmin (T), sqrt(2) x for P, Gaussian beam {BEAM} arcmin FWHM; Planck 2018 I N_l not in the repository",
           "sigma_marg": SIG_MARG, "H0_lcdm": h0_lcdm, "cases": {}}
    with mp.get_context("spawn").Pool(8, initializer=_init, initargs=(target,)) as pool:
        ref = pool.map(spectra, [("cpl", h0_lcdm)])[0]
        out["reference"] = {"H0": ref["H0"], "thetastar": ref["thetastar"], "lmax_camb": ref["lmax_camb"]}
        for name, theta, sigma, npz in cases:
            a, w, bad = tables(theta, sigma)
            assert not bad.any()
            h0s = np.load(npz)["H0"][:N_USE] if npz is not None else [None] * N_USE
            res = pool.map(spectra, [("table", a, w[i], h0s[i]) for i in range(N_USE)], chunksize=2)
            rows = [distance(ref["lensed"], r["lensed"]) for r in res]
            rows_un = [distance(ref["unlensed"], r["unlensed"]) for r in res]
            c = {"theta": theta, "sigma_X": sigma, "summary_lensed": summarise(rows),
                 "summary_unlensed": summarise(rows_un), "per_realisation": rows,
                 "H0": [r["H0"] for r in res], "thetastar_max_abs_diff": float(np.max(np.abs(np.array([r["thetastar"] for r in res]) - target)))}
            if npz is not None:
                # check: re-solve H0 for realisation 0 and compare with the stored value
                import ou_cmb as M2
                make = lambda h: M2.make_params(a=a, w=w[0], H0=h, halofit="takahashi")
                c["check_H0_resolved_r0"] = {"stored": float(h0s[0]), "resolved": float(M2.solve_H0_for_thetastar(target, make)[0])}
            out["cases"][name] = c
            print(name, json.dumps({k: v for k, v in c.items() if k not in ("per_realisation", "H0")}), f"{time.time() - t0:.0f} s", flush=True)
            (HERE / "picos.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    # scaling between the two sigma at theta = 1 (same normals: per-realisation ratio too)
    a1 = out["cases"]["sigma2.633e-02_theta1"]
    a2 = out["cases"][f"sigma{SIG_MARG:.3e}_theta1"]
    sc = {"sigma_ratio": SIG_MARG / SIG_C2}
    for k in ("TT", "TE", "EE", "joint"):
        r1 = np.array([r[k] for r in a1["per_realisation"]])
        r2 = np.array([r[k] for r in a2["per_realisation"]])
        sc[k] = {"p95_ratio": a2["summary_lensed"][k]["p95"] / a1["summary_lensed"][k]["p95"],
                 "median_ratio_per_realisation": float(np.median(r2 / r1)),
                 "exponent_from_median_ratio": float(np.log(np.median(r2 / r1)) / np.log(SIG_MARG / SIG_C2))}
    out["scaling_theta1"] = sc
    out["seconds"] = time.time() - t0
    (HERE / "picos.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    print(json.dumps(sc, indent=1))


if __name__ == "__main__":
    main()
