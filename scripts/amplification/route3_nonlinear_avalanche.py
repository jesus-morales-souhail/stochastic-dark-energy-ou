#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Route 3 — Nonlinear "avalanche" / threshold trigger (toy Langevin).

Hypothesis: tiny seed noise does not need linear gain 1e56; it only needs to
kick the system over an unstable threshold. The large residual is paid by
the background potential (phase-transition / bifurcation), not by multiplying
white noise.

Toy model (overdamped Langevin in x = ln a):

    dX = -V'(X) dx + sigma dW

with double-well / inverted barrier potential:

    V(X) = (1/4) a X^4 - (1/2) b X^2     (a>0, b>0)
    V'(X) = a X^3 - b X

Unstable at X=0 when linearized: V''(0) = -b < 0 → exponential push away
from 0 once |X| exceeds noise-driven escape.

We measure:
  - escape probability from a neighborhood of 0
  - mean |X| after fixed path length
  - fraction of paths with |X| > X_obs (e.g. 1e-5, 1e-4)

This does NOT prove a microphysical DE potential; it maps which (b, sigma)
produce rare large residuals vs smooth BAO-like smallness.

Run:
  python scripts/amplification/route3_nonlinear_avalanche.py
  python scripts/amplification/route3_nonlinear_avalanche.py --heavy
"""

from __future__ import annotations

import argparse
import csv
import os
from concurrent.futures import ProcessPoolExecutor

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, "results", "amplification_routes")
os.makedirs(OUT, exist_ok=True)


def V_prime(X: float, a: float, b: float) -> float:
    return a * X**3 - b * X


def simulate_ensemble(
    a: float,
    b: float,
    sigma: float,
    n_paths: int,
    n_steps: int,
    x_span: float,
    seed0: int,
    thresholds: tuple[float, ...],
) -> dict:
    rng = np.random.default_rng(seed0)
    dx = x_span / n_steps
    X = np.zeros(n_paths)
    # small initial seed ~ sigma scale
    X += rng.normal(0.0, sigma * np.sqrt(dx), size=n_paths)

    for _ in range(n_steps):
        dW = rng.normal(0.0, np.sqrt(dx), size=n_paths)
        X = X + (-V_prime(X, a, b)) * dx + sigma * dW

    absX = np.abs(X)
    row = {
        "a": a,
        "b": b,
        "sigma": sigma,
        "n_paths": n_paths,
        "n_steps": n_steps,
        "x_span": x_span,
        "mean_abs_X": float(np.mean(absX)),
        "rms_X": float(np.sqrt(np.mean(X**2))),
        "p95_abs_X": float(np.percentile(absX, 95)),
        "p99_abs_X": float(np.percentile(absX, 99)),
        "median_abs_X": float(np.median(absX)),
    }
    for thr in thresholds:
        row[f"frac_abs_gt_{thr:.0e}"] = float(np.mean(absX > thr))
    return row


def _job(kwargs):
    return simulate_ensemble(**kwargs)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--heavy", action="store_true")
    args = ap.parse_args()

    n_paths = 20000 if args.heavy else 5000
    n_steps = 8000 if args.heavy else 2500
    x_span = 1.0  # Δx = ln a span O(1) late universe
    thresholds = (1e-6, 1e-5, 1e-4, 1e-3)

    # Grid: weak barrier vs strong; tiny vs moderate noise
    a_vals = [1.0, 10.0]
    b_vals = np.logspace(-6, -1, 12)  # instability strength
    sigma_vals = np.logspace(-8, -3, 12)

    jobs = []
    sid = 0
    for a in a_vals:
        for b in b_vals:
            for sig in sigma_vals:
                jobs.append(
                    {
                        "a": float(a),
                        "b": float(b),
                        "sigma": float(sig),
                        "n_paths": n_paths,
                        "n_steps": n_steps,
                        "x_span": x_span,
                        "seed0": 10_000 + sid,
                        "thresholds": thresholds,
                    }
                )
                sid += 1

    print("=" * 72)
    print("ROUTE 3: Nonlinear avalanche / double-well Langevin (toy)")
    print("=" * 72)
    print(f"Jobs: {len(jobs)}, paths/job={n_paths}, steps={n_steps}, cores={os.cpu_count()}")
    print(f"Potential: V' = a X^3 - b X  (unstable at 0 if b>0)")
    print()

    rows = []
    # ProcessPool — each job is already multi-path vectorized; limit workers
    # to avoid oversubscription: use fewer workers, fat jobs
    workers = max(1, (os.cpu_count() or 8) // 2)
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for i, row in enumerate(ex.map(_job, jobs, chunksize=1)):
            rows.append(row)
            if (i + 1) % 20 == 0 or i == 0:
                print(
                    f"  [{i+1}/{len(jobs)}] a={row['a']} b={row['b']:.2e} σ={row['sigma']:.2e} "
                    f"rms={row['rms_X']:.3e} frac>1e-5={row['frac_abs_gt_1e-05']:.3f}"
                )

    path = os.path.join(OUT, "route3_avalanche_scan.csv")
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print()
    print(f"Wrote {path}")

    # Highlight regimes: mostly small but rare large (avalanche-like)
    # vs always large (excluded by BAO smoothness)
    rare = [
        r
        for r in rows
        if r["frac_abs_gt_1e-05"] > 0.01
        and r["frac_abs_gt_1e-05"] < 0.5
        and r["median_abs_X"] < 1e-4
    ]
    always = [r for r in rows if r["frac_abs_gt_1e-04"] > 0.8]
    never = [r for r in rows if r["frac_abs_gt_1e-05"] < 1e-3]

    print()
    print(f"Regimes: rare large (>1% & <50% paths |X|>1e-5, median small): {len(rare)}")
    print(f"         almost always |X|>1e-4: {len(always)}  (BAO-tension candidates)")
    print(f"         almost never |X|>1e-5: {len(never)}")
    if rare:
        r = rare[len(rare) // 2]
        print(
            f"  Example rare: a={r['a']} b={r['b']:.2e} σ={r['sigma']:.2e} "
            f"median={r['median_abs_X']:.2e} frac>1e-5={r['frac_abs_gt_1e-05']:.3f}"
        )
    print()
    print("INTERPRETATION:")
    print("  Avalanche = threshold dynamics can produce LARGE X without linear")
    print("  gain on white noise — but only in regions of (b,σ) that must still")
    print("  obey BAO smoothness (not always-large). This is a *hypothesis class*,")
    print("  not a derived DE potential. Next theory step: microphysical V(X).")
    print("=" * 72)


if __name__ == "__main__":
    main()
