#!/usr/bin/env python3
"""
200 OU realisations per (theta, Delta, sigma_X) through CAMB; Delta chi^2 against LCDM on Planck bins.

Scheme (NOTAS_ANALISIS_CMB.md): path on a fine grid (step h_tab/16) padded to ln a = 0.4,
Gaussian smoothing Delta, table of 2000 points for CAMB (DarkEnergyPPF.set_w_a_table).
Common random numbers: the same standard normals (seed 20261003) for every theta, Delta
and sigma_X, so X is exactly proportional to sigma_X realisation by realisation.

Per realisation r and observable o:
    Delta chi^2_r = chi^2(OU_r) - chi^2(LCDM), sigma = Planck errors
    (TT: upper error if model > data, lower error otherwise).

Usage:
    ./venv/bin/python realizaciones.py                 # sigma_X = 2.506e-2 only
    ./venv/bin/python realizaciones.py --sigma 1e-3 5e-3
Outputs: resultados/sigma_<sigma>/
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
    # arXiv:2503.14738 eq:w0wa_DESI_CMB (DESI+CMB, no SNe), evaluated at fixed Planck 2018 params
    "DESI_DR2_CMB": (-0.42, -1.75),
}

_LENS = None


def _init_worker():
    global _LENS
    import planck_datos as P

    _LENS = P.load_lensing()


def _observables(args):
    """(a, w) table or ('cpl', w0, wa) -> (Dl_TT[2..29], pp_binned[9], 100 theta_*)."""
    import ou_cmb as M
    import planck_datos as P

    if isinstance(args[0], str):
        pars = M.make_params(w0=args[1], wa=args[2])
    else:
        pars = M.make_params(a=args[0], w=args[1])
    try:
        r = M.run(pars, keep_results=False)
    except Exception as exc:  # CAMB failure: counted and reported by the caller, never dropped silently
        return f"{type(exc).__name__}: {exc}".replace("\n", " ")
    return r["DlTT"][2:30], P.bin_pp(r["Clpp"], _LENS), r["thetastar100"]


def chi2_terms(model: np.ndarray, data: np.ndarray, err_lo: np.ndarray, err_hi: np.ndarray) -> np.ndarray:
    """Sum over the last axis of ((model - data)/sigma_side)^2."""
    d = model - data
    s = np.where(d > 0, err_hi, err_lo)
    return np.sum((d / s) ** 2, axis=-1)


def planck_targets():
    import planck_datos as P

    tt = P.load_tt()
    lens = P.load_lensing()
    return {
        "TT": (tt["Dl"], tt["err_lo"], tt["err_hi"]),
        "pp": (lens["pp"], lens["err"], lens["err"]),
        "thetastar": (np.array([P.THETASTAR100]), np.array([P.THETASTAR100_ERR]), np.array([P.THETASTAR100_ERR])),
    }


def band_gap_and_edge(vals: np.ndarray, ref: np.ndarray, err_lo: np.ndarray, err_hi: np.ndarray):
    """
    68 % band [p16, p84] across realisations vs a reference curve, per bin, in sigma_Planck.
    gap: distance from ref to the band (0 if ref inside); edge: farthest band edge from ref.
    Error side: upper error where the band lies above ref, lower where below.
    Returns (max gap, max edge) over bins.
    """
    lo, hi = np.nanpercentile(vals, [16, 84], axis=0)
    gap = np.where(lo > ref, (lo - ref) / err_hi, np.where(hi < ref, (ref - hi) / err_lo, 0.0))
    edge = np.maximum((hi - ref) / err_hi, (ref - lo) / err_lo)
    return float(np.max(gap)), float(np.max(edge))


def main():
    import ou_cmb as M

    ap = argparse.ArgumentParser()
    ap.add_argument("--sigma", type=float, nargs="+", default=[M.SIGMA_X])
    args = ap.parse_args()

    om_l = M.omega_l0()
    lna_f, tab_idx, main_mask = M.production_grids()
    h = lna_f[1] - lna_f[0]
    lna_tab = lna_f[tab_idx]
    a_tab = np.exp(lna_tab[main_mask])
    z = np.random.default_rng(M.SEED).standard_normal((N_REAL, len(lna_f)))
    targets = planck_targets()

    ctx = mp.get_context("spawn")
    with ctx.Pool(N_WORKERS, initializer=_init_worker) as pool:
        ref_obs = dict(zip(REFS, pool.map(_observables, [("cpl", *v) for v in REFS.values()])))

        for sigma in args.sigma:
            out = HERE / "resultados" / f"sigma_{sigma:.3e}"
            out.mkdir(parents=True, exist_ok=True)
            summary = {
                "sigma_X": sigma, "n_real": N_REAL, "seed": M.SEED, "Omega_L0": om_l,
                "n_tab": M.N_TAB, "k_path": M.K_PATH, "pad_lna": M.PAD,
                "references_wa_w0": REFS, "cases": {},
            }
            for theta in THETAS:
                x_unit = M.ou_uniform(z, h, theta, 1.0)  # X = sigma * x_unit (common random numbers)
                for delta in DELTAS:
                    xs = M.smooth_uniform(sigma * x_unit, h, delta)
                    rho, w, bad = M.rho_w_table(lna_tab, main_mask, xs[:, tab_idx], om_l)
                    good = np.flatnonzero(~bad)
                    res = pool.map(_observables, [(a_tab, w[i]) for i in good], chunksize=4)
                    obs = {
                        "TT": np.full((N_REAL, 28), np.nan),
                        "pp": np.full((N_REAL, 9), np.nan),
                        "thetastar": np.full((N_REAL, 1), np.nan),
                    }
                    camb_fail = np.zeros(N_REAL, bool)
                    fail_msgs = {}
                    for i, rr in zip(good, res):
                        if isinstance(rr, str):
                            camb_fail[i] = True
                            fail_msgs[int(i)] = rr
                            continue
                        obs["TT"][i], obs["pp"][i], obs["thetastar"][i] = rr
                    tag = f"theta{theta:g}_Delta{delta:.2f}"
                    np.savez(out / f"observables_{tag}.npz", w_min=np.nanmin(w, axis=1),
                             w_max=np.nanmax(w, axis=1), bad=bad, camb_fail=camb_fail, **obs)
                    rho_bad = bad.copy()
                    bad = bad | camb_fail

                    case = {"n_bad_rho_le_0": int(rho_bad.sum()), "rho_bad_indices": np.flatnonzero(rho_bad).tolist(),
                            "n_camb_fail": int(camb_fail.sum()), "camb_fail_messages": fail_msgs,
                            "w_range_all": [float(np.nanmin(w)), float(np.nanmax(w))], "obs": {}}
                    for o, (data, elo, ehi) in targets.items():
                        chi_ref = {k: float(chi2_terms(v[{"TT": 0, "pp": 1, "thetastar": 2}[o]]
                                                      if o != "thetastar" else np.array([v[2]]), data, elo, ehi))
                                   for k, v in ref_obs.items()}
                        chi_r = chi2_terms(obs[o], data, elo, ehi)
                        chi_r[bad] = np.nan
                        dchi = chi_r - chi_ref["LCDM"]
                        np.savetxt(out / f"dchi2_{o}_{tag}.txt", dchi,
                                   header=f"Delta chi2_r = chi2(OU_r) - chi2(LCDM), {o}, {tag}, sigma_X={sigma:.6e}, "
                                          f"NaN = rho_DE<=0 (reading (a) broken) or CAMB failure (see resumen.json)")
                        ok = dchi[~np.isnan(dchi)]
                        if len(ok) == 0:
                            case["obs"][o] = {"error": "no valid realisations"}
                            continue
                        ref_lcdm = ref_obs["LCDM"][{"TT": 0, "pp": 1, "thetastar": 2}[o]]
                        gap, edge = band_gap_and_edge(obs[o][~bad], np.atleast_1d(ref_lcdm), elo, ehi)
                        entry = {
                            "dchi2_p50": float(np.percentile(ok, 50)),
                            "dchi2_p95": float(np.percentile(ok, 95)),
                            "frac_dchi2_gt_4": float(np.mean(ok > 4)),
                            "band68_gap_vs_LCDM_sigma": gap,
                            "band68_edge_vs_LCDM_sigma": edge,
                            "dchi2_vs_CPL_repo_p95": float(np.percentile(chi_r[~bad] - chi_ref["CPL_repo"], 95)),
                            "dchi2_vs_CPL_joint_p95": float(np.percentile(chi_r[~bad] - chi_ref["CPL_joint"], 95)),
                            "chi2_refs": chi_ref,
                        }
                        if o == "thetastar":
                            v = obs[o][~bad, 0]
                            entry.update(mean=float(v.mean()), p16=float(np.percentile(v, 16)),
                                         p84=float(np.percentile(v, 84)))
                        case["obs"][o] = entry
                    summary["cases"][tag] = case
                    print(f"sigma={sigma:.3e} {tag}: rho<=0={rho_bad.sum()} camb_fail={camb_fail.sum()} "
                          + " ".join(f"{o}:p95={case['obs'][o].get('dchi2_p95', float('nan')):.3g}" for o in targets), flush=True)

            summary["reference_observables"] = {
                k: {"thetastar100": float(v[2]),
                    **{f"dchi2_{o}_vs_LCDM": float(
                        chi2_terms(v[i] if o != "thetastar" else np.array([v[2]]), *targets[o])
                        - chi2_terms(ref_obs["LCDM"][i] if o != "thetastar" else np.array([ref_obs["LCDM"][2]]),
                                     *targets[o]))
                       for i, o in enumerate(targets)}}
                for k, v in ref_obs.items()
            }
            (out / "resumen.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
