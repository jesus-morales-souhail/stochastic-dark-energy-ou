#!/usr/bin/env python3
"""
Crossings of w = -1 (sign changes of d rho_DE / d ln a on 0 < z < z*) in the 40 posterior histories per ceiling of cmb_final
(theta = 1, 3, 10; same seeds, no CAMB) as a function of the Gaussian smoothing scale Delta in ln a. The raw posterior paths are
regenerated with the sampling code of cmb_final.prep (identical seeds); the saved smoothed histories of cmb_final_theta_*.npz are
reproduced bit for bit as a check (Delta of cmb_final). Writes cruces_w.json.

    OMP_NUM_THREADS=2 ../venv/bin/python cruces_w.py
"""
from __future__ import annotations
import os
os.environ.setdefault("OMP_NUM_THREADS", "2")
import json, sys
from pathlib import Path
import numpy as np
from scipy.linalg import cho_factor, cho_solve
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent)); sys.path.insert(0, str(HERE.parent / "cmb"))
import cmb_final as CF  # noqa: E402
import marginal_guiada as MG  # noqa: E402
import exacto as X  # noqa: E402

DELTAS = (0.0125, 0.025, 0.05, 0.1, 0.2)
out = {"deltas": list(DELTAS), "n_run": CF.N_RUN, "seed": CF.SEED, "definition": "crossing = sign change of d rho_DE/d ln a between consecutive grid nodes on 0 <= u <= u* (0 < z < z*), rho_DE = Omega_Lambda + Xs - Xs(0)", "thetas": {}}
Bm, A = MG.build_B(), MG.build_A()
for ic, (theta, sigma, delta_cf) in enumerate(CF.CELLS):
    key = f"theta_{theta:g}"
    K = MG.kernel(theta); AK = A @ K; Msig = AK @ A.T
    rng = np.random.default_rng([CF.SEED, ic])
    r = MG.D - MG.B0
    Sf = cho_factor(MG.C + sigma**2 * Msig, lower=True)
    SB = cho_solve(Sf, Bm); Sp = np.linalg.inv(Bm.T @ SB); dp_hat = Sp @ (Bm.T @ cho_solve(Sf, r))
    dp = dp_hat[None, :] + rng.standard_normal((CF.N_SAMP, 3)) @ np.linalg.cholesky(Sp).T
    x0 = MG.S.ou_paths(CF.G, sigma, theta, CF.N_SAMP, rng)
    n0 = rng.standard_normal((CF.N_SAMP, MG.NROWS)) @ MG.LC.T
    R = r[None, :] - dp @ Bm.T - x0 @ A.T - n0
    x = x0 + sigma**2 * cho_solve(Sf, R.T).T @ AK
    p = MG.P0[None, :] + dp
    h, wm = p[:, 0], p[:, 1]
    ol = (h * h - wm - X.OMR) / (h * h)
    saved = np.load(HERE / f"cmb_final_{key}.npz")
    idx_run = saved["idx_run"]
    # reproduce prep's candidate set (first 4*N_RUN with rho > 0 and H^2 > 0) so that the backward extension uses the same normals
    rho_raw = ol[:, None] + x - x[:, [0]]
    neg_raw = (rho_raw <= 0).any(axis=1)
    e2 = (wm / h**2)[:, None] * CF.G.p3 + (X.OMR / h**2)[:, None] * CF.G.p4 + (1 - wm / h**2 - X.OMR / h**2)[:, None] + x - x[:, [0]]
    bad = e2.min(axis=1) <= 0
    idx_keep = np.nonzero(~bad & ~neg_raw)[0][: 4 * CF.N_RUN]
    pos_in_keep = np.array([int(np.nonzero(idx_keep == i)[0][0]) for i in idx_run])
    rec = {"sigma_X": sigma, "delta_cmb_final": delta_cf, "n": int(len(idx_run)), "per_delta": {}}
    for delta in DELTAS:
        rng2 = np.random.default_rng([CF.SEED, ic, 1])
        Xs_all = CF.smooth_ext(x[idx_keep], theta, sigma, delta, rng2)
        Xs = Xs_all[pos_in_keep]
        if abs(delta - delta_cf) < 1e-12:
            rec["check_saved_max_abs_diff"] = float(np.max(np.abs(Xs - saved["Xs"])))
        m = CF.MAIN
        rho = ol[idx_run][:, None] + Xs - Xs[:, [CF.I0]]
        d = np.diff(rho[:, m], axis=1)  # d rho / du up to sign; crossings of w = -1 are sign changes of d rho/d ln a = -d rho/du
        s = np.sign(d); s[s == 0] = 1
        ncross = np.sum(s[:, 1:] * s[:, :-1] < 0, axis=1)
        rec["per_delta"][f"{delta:g}"] = {"frac_at_least_one": float(np.mean(ncross >= 1)), "n_cross_median": float(np.median(ncross)),
                                          "n_cross_p16_p84": [float(np.percentile(ncross, 16)), float(np.percentile(ncross, 84))], "n_cross_min_max": [int(ncross.min()), int(ncross.max())],
                                          "abs_wp1_max_median": float(np.median(np.nanmax(np.abs(-1.0 + np.gradient(np.log(np.where(rho > 0, rho, np.nan)), CF.U_EXT, axis=1)[:, m] / 3.0 + 1.0), axis=1)))}
    out["thetas"][key] = rec
    print(key, "check", rec.get("check_saved_max_abs_diff"), {k: (v["frac_at_least_one"], v["n_cross_median"]) for k, v in rec["per_delta"].items()}, flush=True)
(HERE / "cruces_w.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
