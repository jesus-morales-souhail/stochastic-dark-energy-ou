#!/usr/bin/env python3
"""
200 OU realisations per (theta, Delta, sigma_X) through CAMB, with H0 re-solved per curve.

Scheme (NOTAS_ANALISIS_CMB.md): path on a fine grid (step h_tab/16) padded to ln a = 0.4,
Gaussian smoothing Delta, 2000-point w(a) table for CAMB (DarkEnergyPPF.set_w_a_table).
Common random numbers: the same standard normals (seed 20261003) for every theta, Delta and
sigma_X, so X is exactly proportional to sigma_X realisation by realisation.

Per curve (every OU realisation and every reference curve):
  1. 100 theta_* at the Planck 2018 H0 (background only) — kept for information.
  2. H0 by bisection on get_background so that 100 theta_* equals the LCDM reference value
     (omega_b, omega_c, mnu fixed; flat).
  3. TT (2 <= l < 30) and binned C_L^phiphi at that H0. Non-linear lensing: --halofit.
Failures (H0 not bracketed, CAMB error, silent NaN spectra) are flagged per observable and
counted; nothing is dropped silently. Statistics: estadisticos.py.

The previous version (H0 fixed, HMCode mead2020) is kept in diagnosticos/obsoleto/.

Usage:
    ./venv/bin/python realizaciones.py --halofit takahashi                 # sigma_X = 2.506e-2
    ./venv/bin/python realizaciones.py --halofit takahashi --sigma 1e-3 5e-3
    ./venv/bin/python realizaciones.py --halofit takahashi --subdir malla --thetas 0.001 1 --deltas 0.05 \
        --sigma $(python3 -c "import numpy as n;print(*n.geomspace(0.025063155506929828,1e-4,10))")
Outputs: resultados/<halofit>/sigma_<sigma>/
"""

from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "2")

import argparse
import json
import multiprocessing as mp
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
N_REAL = int(os.environ.get("N_REAL", 200))  # override only for smoke tests
THETAS = (0.001, 1.0, 20.0)
DELTAS = (0.05, 0.10)
N_WORKERS = 8

# Reference curves (w0, wa)
REFS = {
    "LCDM": (-1.0, 0.0),
    # results/eos_cpl_desi_dr2/eos_cpl_summary.json (pure_CPL; DIAGONAL errors)
    "CPL_repo": (-0.9898956968059271, -0.015902984131323503),
    # results/joint_w0wa_sigma/joint_w0wa_sigma.json (model 'cpl'; cov_mode contradictory)
    "CPL_joint": (-1.0016540106177185, 0.0172470528926633),
    # arXiv:2503.14738 eq:w0wa_DESI_CMB (DESI+CMB, no SNe)
    "DESI_DR2_CMB": (-0.42, -1.75),
}

_LENS = None
_HALOFIT = None
_TARGET = None


def _init_worker(halofit, target):
    global _LENS, _HALOFIT, _TARGET
    import planck_datos as P

    _LENS, _HALOFIT, _TARGET = P.load_lensing(), halofit, target


def _curve(task):
    """task = ('cpl', w0, wa) or ('table', a, w). Returns dict with H0, theta_*, TT, pp and failure messages."""
    import ou_cmb as M
    import planck_datos as P

    if task[0] == "cpl":
        def make(h):
            return M.make_params(w0=task[1], wa=task[2], H0=h, halofit=_HALOFIT)
    else:
        def make(h):
            return M.make_params(a=task[1], w=task[2], H0=h, halofit=_HALOFIT)

    out = {"H0": np.nan, "thetastar_H0planck": np.nan, "thetastar_solved": np.nan,
           "TT": np.full(28, np.nan), "pp": np.full(9, np.nan), "fail": {}}
    try:
        out["thetastar_H0planck"] = M.thetastar_background(make(M.PLANCK18["H0"]))
        out["H0"], out["thetastar_solved"], _ = M.solve_H0_for_thetastar(_TARGET, make)
    except Exception as exc:  # noqa: BLE001
        msg = f"H0: {type(exc).__name__}: {exc}".replace("\n", " ")
        out["fail"] = {"TT": msg, "pp": msg, "H0": msg}
        return out
    try:
        r = M.run(make(out["H0"]), keep_results=False)
    except Exception as exc:  # noqa: BLE001
        msg = f"CAMB: {type(exc).__name__}: {exc}".replace("\n", " ")
        out["fail"] = {"TT": msg, "pp": msg}
        return out
    tt, pp = r["DlTT"][2:30], P.bin_pp(r["Clpp"], _LENS)
    for name, v in (("TT", tt), ("pp", pp)):
        if np.all(np.isfinite(v)):
            out[name] = v
        else:
            out["fail"][name] = f"silent NaN ({int(np.sum(~np.isfinite(v)))} bins)"
    return out


def main():
    import ou_cmb as M
    import estadisticos as E

    ap = argparse.ArgumentParser()
    ap.add_argument("--sigma", type=float, nargs="+", default=[M.SIGMA_X])
    ap.add_argument("--halofit", default="takahashi", choices=["takahashi", "lineal", "mead2020"])
    ap.add_argument("--thetas", type=float, nargs="+", default=list(THETAS))
    ap.add_argument("--deltas", type=float, nargs="+", default=list(DELTAS))
    ap.add_argument("--subdir", default="", help="extra folder under resultados/<halofit>/ (e.g. malla)")
    args = ap.parse_args()
    halofit = None if args.halofit == "mead2020" else args.halofit

    om_l = M.omega_l0()
    target = M.thetastar_background(M.make_params(halofit=halofit))  # LCDM reference 100 theta_*
    lna_f, tab_idx, main_mask = M.production_grids()
    h = lna_f[1] - lna_f[0]
    lna_tab = lna_f[tab_idx]
    a_tab = np.exp(lna_tab[main_mask])
    z = np.random.default_rng(M.SEED).standard_normal((N_REAL, len(lna_f)))
    base = HERE / "resultados" / args.halofit / args.subdir

    ctx = mp.get_context("spawn")
    with ctx.Pool(N_WORKERS, initializer=_init_worker, initargs=(halofit, target)) as pool:
        refs = dict(zip(REFS, pool.map(_curve, [("cpl", *v) for v in REFS.values()])))
        for sigma in args.sigma:
            out = base / f"sigma_{sigma:.3e}"
            out.mkdir(parents=True, exist_ok=True)
            meta = {"sigma_X": sigma, "n_real": N_REAL, "seed": M.SEED, "Omega_L0_planck": om_l,
                    "n_tab": M.N_TAB, "k_path": M.K_PATH, "pad_lna": M.PAD, "halofit": args.halofit,
                    "thetastar100_target_LCDM": target, "references_w0_wa": REFS,
                    "references": {k: {kk: (vv.tolist() if isinstance(vv, np.ndarray) else vv)
                                       for kk, vv in v.items()} for k, v in refs.items()}}
            (out / "meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
            for theta in args.thetas:
                x_unit = M.ou_uniform(z, h, theta, 1.0)  # X = sigma * x_unit (common random numbers)
                for delta in args.deltas:
                    xs = M.smooth_uniform(sigma * x_unit, h, delta)
                    _, w, rho_bad = M.rho_w_table(lna_tab, main_mask, xs[:, tab_idx], om_l)
                    good = np.flatnonzero(~rho_bad)
                    res = pool.map(_curve, [("table", a_tab, w[i]) for i in good], chunksize=4)
                    arr = {"H0": np.full(N_REAL, np.nan), "thetastar_H0planck": np.full(N_REAL, np.nan),
                           "thetastar_solved": np.full(N_REAL, np.nan),
                           "TT": np.full((N_REAL, 28), np.nan), "pp": np.full((N_REAL, 9), np.nan)}
                    fails = {}
                    for i, r in zip(good, res):
                        for k in arr:
                            arr[k][i] = r[k]
                        if r["fail"]:
                            fails[int(i)] = r["fail"]
                    tag = f"theta{theta:g}_Delta{delta:.2f}"
                    np.savez(out / f"observables_{tag}.npz", rho_bad=rho_bad,
                             w_min=np.nanmin(w, axis=1), w_max=np.nanmax(w, axis=1), **arr)
                    (out / f"fallos_{tag}.json").write_text(
                        json.dumps({"rho_le_0": np.flatnonzero(rho_bad).tolist(), "per_realisation": fails},
                                   indent=1), encoding="utf-8")
                    nf = {o: sum(o in f for f in fails.values()) for o in ("TT", "pp", "H0")}
                    print(f"[{args.halofit}] sigma={sigma:.3e} {tag}: rho<=0={rho_bad.sum()} "
                          f"fallos TT={nf['TT']} pp={nf['pp']} H0={nf['H0']} "
                          f"H0 p50={np.nanmedian(arr['H0']):.3f}", flush=True)
            E.resumen(out)


if __name__ == "__main__":
    main()
