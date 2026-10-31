#!/usr/bin/env python3
"""
Exploration only (EXPLORACION.md, point 4): CAMB cost per realisation with the realizaciones.py setup
(OMP_NUM_THREADS = 2 per worker, 8 workers; Takahashi halofit; lmax 2500 with lens_potential_accuracy = 1;
H0 solved by bisection on get_background so that 100 theta_* equals the LCDM value).

History: survivor 5492 with the gap draw of prueba_camb.py (same seed, same pick). Also times the
LCDM reference with the same calls, to compare with the OU runs of realizaciones.py.

    OMP_NUM_THREADS=2 ../../venv/bin/python tiempos.py

Writes cmb_real/tiempos.json.
"""

from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "2")

import json
import time

import numpy as np

import prueba_camb as P


def main() -> None:
    sv = {r["index"]: r for r in json.loads((P.HERE / "supervivientes.json").read_text())["model2"]}
    lam_all, O = P.stored_points()
    lam_gap = np.arange(-np.log(1000.0) + 0.02, -np.log(4.0) - 0.01, 0.02)
    rng = np.random.default_rng(P.GAP_SEED)
    j = 5492
    wm = sv[j]["passes"]["3d_wm_free"]["omega_m"]
    mean, Ls = P.gap_sampler(lam_all, O[:, j], lam_gap)
    draws = mean[:, None] + Ls @ rng.standard_normal((len(lam_gap), P.N_GAP))
    pick = int(np.nonzero(np.all(draws > 0, axis=0))[0][0])
    a, ode, rho, w = P.history(lam_all, O[:, j], lam_gap, draws[:, pick], wm)
    H0r = 100.0 * np.sqrt((wm + P.B.ORAD) / (1.0 - ode[-1]))

    def make(h):
        pars = P.M.make_params(H0=h, halofit="takahashi")
        pars.set_cosmology(H0=h, ombh2=P.OMB, omch2=wm - P.OMB - P.OMNU, mnu=0.06, omk=0.0, tau=P.M.PLANCK18["tau"])
        de = P.DarkEnergyPPF()
        de.set_w_a_table(a, w)
        pars.DarkEnergy = de
        return pars

    out = {"omp_threads": os.environ["OMP_NUM_THREADS"], "index": j, "omega_m": wm, "H0_realisation": H0r}
    target = P.M.thetastar_background(P.M.make_params(halofit="takahashi"))
    out["target_thetastar100"] = target
    scan = []
    for h in (40, 45, 50, 55, 59.82, 65, 70, 80, 90, 100, 120, 150):
        t0 = time.time()
        try:
            th = P.M.thetastar_background(make(h))
            scan.append({"H0": h, "thetastar100": th, "s": time.time() - t0})
        except Exception as exc:  # noqa: BLE001
            scan.append({"H0": h, "error": f"{type(exc).__name__}: {exc}"[:160], "s": time.time() - t0})
        print(scan[-1], flush=True)
    out["thetastar_scan"] = scan
    ok = [s for s in scan if "thetastar100" in s]
    out["t_background_mean_s"] = float(np.mean([s["s"] for s in ok]))
    lo = min(s["H0"] for s in ok)
    try:
        t0 = time.time()
        h, th, it = P.M.solve_H0_for_thetastar(target, make, lo=lo, hi=150.0)
        out["solve"] = {"H0": h, "thetastar100": th, "iterations": it, "s": time.time() - t0}
    except Exception as exc:  # noqa: BLE001
        out["solve"] = {"error": f"{type(exc).__name__}: {exc}"[:300], "bracket": [lo, 150.0]}
    print(out["solve"], flush=True)
    t0 = time.time()
    r = P.M.run(make(H0r), keep_results=False)
    out["t_full_run_realisation_s"] = time.time() - t0
    out["full_run_finite"] = bool(np.all(np.isfinite(r["DlTT"][2:30])) and np.all(np.isfinite(r["Clpp"][2:2001])))
    # LCDM reference with the same calls
    t0 = time.time()
    P.M.thetastar_background(P.M.make_params(halofit="takahashi"))
    out["t_background_lcdm_s"] = time.time() - t0
    t0 = time.time()
    P.M.solve_H0_for_thetastar(target, lambda hh: P.M.make_params(H0=hh, halofit="takahashi"))
    out["t_solve_lcdm_s"] = time.time() - t0
    t0 = time.time()
    P.M.run(P.M.make_params(halofit="takahashi"), keep_results=False)
    out["t_full_run_lcdm_s"] = time.time() - t0
    print(json.dumps(out, indent=1))
    (P.HERE / "tiempos.json").write_text(json.dumps(out, indent=1), encoding="utf-8")


if __name__ == "__main__":
    main()
