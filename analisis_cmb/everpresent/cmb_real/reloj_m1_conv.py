#!/usr/bin/env python3
"""
Exploration only (RELOJ_M1.md, "Convergencia en ln a"): survival of Model 1 (alpha = 0.01, the 2000 seeds of
reloj_m1.py) with steps uniform in ln a for n = 3e5 and 1e6, and two estimates of the continuous-monitoring limit:

 (i)  weighted fit p(Delta) = p_inf + c sqrt(Delta) to n = 1e4, 1e5 (from reloj_m1.json), 3e5, 1e6, with
      Delta = ln(1e4) / n_steps and binomial errors. The four runs share the seeds, so their errors are correlated;
      the fit treats them as independent (stated), and the error is also given scaled by sqrt(chi2/dof) when > 1.
 (ii) Broadie-Glasserman-Kou continuity correction for discretely monitored barriers: the discrete run with the barrier
      moved towards the process by beta sigma_step, beta = -zeta(1/2)/sqrt(2 pi) = 0.5826, approximates continuous
      monitoring. Here the monitored process is u = K alpha S / V with S a Wiener process in V, and sigma_step =
      K alpha sqrt(dV) / V; the barrier -(w_m a^-3 + w_r a^-4) moves with time, u is not a pure Brownian motion and
      Lambda is held fixed between steps, so the correction is approximate. Applied to n = 1e5 (and, as a check,
      to n = 1e4).

    python3 reloj_m1_conv.py    # writes cmb_real/reloj_m1_conv.json
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
from scipy.special import zeta

import reloj_m1 as R

HERE = Path(__file__).resolve().parent
BETA = float(-zeta(0.5) / np.sqrt(2 * np.pi))


def main() -> None:
    prev = json.loads((HERE / "reloj_m1.json").read_text())["clocks"]
    out = {"beta_bgk": BETA, "runs": {}}
    for n in (300000, 1000000):
        t0 = time.time()
        s = R.summary(R.run_grid(R.grid("lna", n), R.Streams(R.N_REAL), R.N_REAL))
        s["seconds"] = time.time() - t0
        out["runs"][f"lna_n{n}"] = s
        print(n, json.dumps(s), flush=True)
        (HERE / "reloj_m1_conv.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    for n in (10000, 100000):
        t0 = time.time()
        s = R.summary(R.run_grid(R.grid("lna", n), R.Streams(R.N_REAL), R.N_REAL, bgk=BETA))
        s["seconds"] = time.time() - t0
        out["runs"][f"lna_n{n}_bgk"] = s
        print("bgk", n, json.dumps(s), flush=True)
    # (i) fit
    ns = [10000, 100000, 300000, 1000000]
    p = np.array([prev["lna_n10000"]["frac_alive"], prev["lna_n100000"]["frac_alive"],
                  out["runs"]["lna_n300000"]["frac_alive"], out["runs"]["lna_n1000000"]["frac_alive"]])
    nsteps = np.array([len(R.grid("lna", n)) - 1 for n in ns])
    d = np.log(1.0 / R.A_START) / nsteps
    sd = np.sqrt(p * (1 - p) / R.N_REAL)
    X = np.column_stack([np.ones_like(d), np.sqrt(d)])
    W = np.diag(1 / sd**2)
    cov = np.linalg.inv(X.T @ W @ X)
    beta = cov @ X.T @ W @ p
    res = p - X @ beta
    chi2 = float(res @ W @ res)
    scale = max(1.0, np.sqrt(chi2 / (len(p) - 2)))
    out["fit_sqrt"] = {"n": ns, "delta_lna": d.tolist(), "p": p.tolist(), "sd": sd.tolist(),
                       "p_inf": float(beta[0]), "p_inf_err": float(np.sqrt(cov[0, 0])),
                       "p_inf_err_scaled": float(np.sqrt(cov[0, 0]) * scale), "c": float(beta[1]),
                       "chi2": chi2, "dof": len(p) - 2, "residuals": res.tolist(),
                       "note": "binomial errors treated as independent although the runs share seeds"}
    print(json.dumps(out["fit_sqrt"], indent=1))
    (HERE / "reloj_m1_conv.json").write_text(json.dumps(out, indent=1), encoding="utf-8")


if __name__ == "__main__":
    main()
