#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Route 2 — Late-time horizon-exit analogue for an OU residual X(x).

Hypothesis: mean-reversion rate theta (or gamma) becomes mode-dependent and
drops when a mode is "super-horizon", freezing the amplitude — analogous to
inflationary freeze-out, but applied to a residual on x = ln a during late
acceleration.

We simulate the SDE (Euler–Maruyama, many paths in parallel):

    dX = -theta(x) X dx + sigma(x) dW

with theta(x) = theta_sub * S(x) + theta_super * (1-S(x)),
S = smoothstep for "horizon exit" of a fiducial comoving scale k.

WARNING: This is a *toy* analogue. It does NOT derive that DE residuals
freeze like inflaton modes. It quantifies how much frozen amplitude you get
from a given seed once theta→0 outside the horizon.

Run:
  python scripts/amplification/route2_late_horizon_exit.py
  python scripts/amplification/route2_late_horizon_exit.py --heavy
"""

from __future__ import annotations

import argparse
import csv
import os
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, "results", "amplification_routes")
os.makedirs(OUT, exist_ok=True)

OM, OL = 0.315, 0.685
H0 = 67.4  # km/s/Mpc only for ratios


def E_z(z: float) -> float:
    return float(np.sqrt(OM * (1 + z) ** 3 + OL))


def a_of_x(x: float) -> float:
    return float(np.exp(x))


def z_of_x(x: float) -> float:
    return float(np.exp(-x) - 1.0)


def horizon_exit_switch(x: float, x_exit: float, width: float = 0.05) -> float:
    """
    S→1 deep inside horizon (restoring), S→0 outside (freeze).
    x increases toward future (a grows). Exit when x > x_exit.
    """
    # logistic: S = 1/(1+exp((x-x_exit)/width))  → 1 before exit, 0 after
    return float(1.0 / (1.0 + np.exp((x - x_exit) / width)))


def simulate_path(
    x_grid: np.ndarray,
    theta_sub: float,
    theta_super: float,
    sigma_amp: float,
    x_exit: float,
    seed: int,
) -> dict:
    rng = np.random.default_rng(seed)
    X = 0.0
    xs = []
    Xs = []
    for i in range(len(x_grid) - 1):
        x = x_grid[i]
        dx = x_grid[i + 1] - x
        S = horizon_exit_switch(x, x_exit)
        th = theta_sub * S + theta_super * (1.0 - S)
        # diffusion may also freeze outside horizon
        sig = sigma_amp * S  # no new kicks outside
        dW = rng.normal() * np.sqrt(max(dx, 0.0))
        X = X + (-th * X) * dx + sig * dW
        xs.append(x)
        Xs.append(X)
    xs.append(x_grid[-1])
    Xs.append(X)
    arr = np.asarray(Xs)
    return {
        "seed": seed,
        "X_final": float(arr[-1]),
        "X_max_abs": float(np.max(np.abs(arr))),
        "X_rms": float(np.sqrt(np.mean(arr**2))),
        "theta_sub": theta_sub,
        "theta_super": theta_super,
        "sigma_amp": sigma_amp,
        "x_exit": x_exit,
    }


def _job(payload):
    return simulate_path(**payload)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--heavy", action="store_true", help="more paths / denser grid")
    ap.add_argument("--paths", type=int, default=None)
    args = ap.parse_args()

    n_paths = args.paths or (4000 if args.heavy else 800)
    n_x = 4000 if args.heavy else 1500

    # Path covering DESI-like span and late acceleration: z~1100 → 0
    # x = ln a = -ln(1+z)
    z_hi, z_lo = 1100.0, 0.0
    x_start, x_end = -np.log(1 + z_hi), -np.log(1 + z_lo)
    x_grid = np.linspace(x_start, x_end, n_x)

    # Exit near z~0.5 (late acceleration) — illustrative
    x_exit = -np.log(1.5)

    print("=" * 72)
    print("ROUTE 2: Late horizon-exit analogue (OU with theta(x) → freeze)")
    print("=" * 72)
    print(f"x: {x_start:.3f} → {x_end:.3f}  (z {z_hi} → {z_lo})")
    print(f"x_exit ≈ {x_exit:.3f}  (z≈0.5),  paths={n_paths},  grid={n_x}")
    print()

    # Scans: seed-level diffusion vs freeze strength
    # sigma_amp chosen so integrated noise without freeze is tiny or moderate
    scan = []
    for theta_sub in [0.5, 1.2, 3.0]:
        for sigma_amp in [1e-6, 1e-5, 1e-4, 1e-3]:
            for theta_super in [0.0, 1e-4]:
                scan.append((theta_sub, theta_super, sigma_amp))

    # Parallel ensembles
    jobs = []
    job_id = 0
    for theta_sub, theta_super, sigma_amp in scan:
        for p in range(n_paths):
            jobs.append(
                {
                    "x_grid": x_grid,
                    "theta_sub": theta_sub,
                    "theta_super": theta_super,
                    "sigma_amp": sigma_amp,
                    "x_exit": x_exit,
                    "seed": 1000 * job_id + p,
                    "_meta": (theta_sub, theta_super, sigma_amp),
                }
            )
        job_id += 1

    # Don't pass _meta to worker — strip
    print(f"Total path simulations: {len(jobs)} on {os.cpu_count()} cores")
    # Group by meta for efficiency: run in batches per parameter set
    summary_rows = []
    for theta_sub, theta_super, sigma_amp in scan:
        payloads = [
            {
                "x_grid": x_grid,
                "theta_sub": theta_sub,
                "theta_super": theta_super,
                "sigma_amp": sigma_amp,
                "x_exit": x_exit,
                "seed": 10_000_000
                + int(theta_sub * 1000)
                + int(sigma_amp * 1e9) % 100000
                + p,
            }
            for p in range(n_paths)
        ]
        finals = []
        maxabs = []
        with ProcessPoolExecutor(max_workers=os.cpu_count() or 8) as ex:
            for res in ex.map(_job, payloads, chunksize=max(1, n_paths // 64)):
                finals.append(res["X_final"])
                maxabs.append(res["X_max_abs"])
        finals = np.asarray(finals)
        maxabs = np.asarray(maxabs)
        row = {
            "theta_sub": theta_sub,
            "theta_super": theta_super,
            "sigma_amp": sigma_amp,
            "x_exit": x_exit,
            "n_paths": n_paths,
            "mean_X_final": float(np.mean(finals)),
            "std_X_final": float(np.std(finals)),
            "rms_X_final": float(np.sqrt(np.mean(finals**2))),
            "mean_X_max_abs": float(np.mean(maxabs)),
            "p95_abs_X_final": float(np.percentile(np.abs(finals), 95)),
        }
        summary_rows.append(row)
        print(
            f"  θ_sub={theta_sub:g} θ_sup={theta_super:g} σ={sigma_amp:.0e}  "
            f"→ rms|X_f|={row['rms_X_final']:.3e}  p95={row['p95_abs_X_final']:.3e}"
        )

    path = os.path.join(OUT, "route2_horizon_exit_scan.csv")
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()))
        w.writeheader()
        w.writerows(summary_rows)
    print()
    print(f"Wrote {path}")

    # Analytic freeze insight
    print()
    print("ANALYTIC INSIGHT:")
    print("  If sigma_amp is the bare kick scale and freeze sets theta→0 after exit,")
    print("  frozen rms is set by integrated noise *before* exit — order sigma_amp")
    print("  times sqrt(effective x-duration inside horizon), NOT automatically 1e-5")
    print("  from a 1e-61 seed unless sigma_amp itself is already raised (Route 1)")
    print("  or a large number of coherent e-folds acts on a different seed (H/MPl).")
    print("  Late acceleration has Δx = ln(a0/a_exit) ≲ O(1), not 60 e-folds.")
    dx_late = float(x_end - x_exit)
    print(f"  Δx after z=0.5 exit in this setup: {dx_late:.3f}  (≪ 60)")
    print("=" * 72)


if __name__ == "__main__":
    main()
