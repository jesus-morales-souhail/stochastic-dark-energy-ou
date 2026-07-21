#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Route 1 — Scale-invariant / local causal-set seed (NOT a magic amplifier).

Hypothesis: the Sorkin count N ~ 10^122 uses the *full* observable horizon.
If the DE residual couples only to a *local correlation patch* of scale L
(or an effective N_eff << N_BH), then

    sigma_0_eff = 1 / sqrt(N_eff)

can be much larger than 10^{-61} without inventing a linear gain factor.

This script:
  - Scans N_eff and L/R_H for area (p=2), volume (p=3), 4-volume (p=4) counting
  - Reports N_eff needed for sigma_0_eff in {1e-5, 1.5e-4, 1e-4}
  - Does NOT claim the correct causal-set counting is known

Run (uses all CPU cores for grid):
  python scripts/amplification/route1_local_causal_set_seed.py
"""

from __future__ import annotations

import csv
import os
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, "results", "amplification_routes")
os.makedirs(OUT, exist_ok=True)

N_BH = 1.0e122
SIGMA_0_GLOBAL = 1.0 / np.sqrt(N_BH)  # ~1e-61
R_H_M = 1.4e26  # ~ Hubble radius, metres (order of magnitude)
L_P_M = 1.616255e-35

TARGETS = {
    "euclid_1e-5": 1.0e-5,
    "desi_1p5e-4": 1.5e-4,
    "window_1e-4": 1.0e-4,
}


def sigma_from_N(N_eff: float) -> float:
    return 1.0 / np.sqrt(max(N_eff, 1.0))


def N_from_sigma(sig: float) -> float:
    return 1.0 / (sig * sig)


def N_eff_from_fraction(frac: float, power: int) -> float:
    """N_eff = N_BH * (L/R_H)^power  (same counting dimension as global N)."""
    return N_BH * (frac ** power)


def worker_row(args):
    power, frac = args
    Ne = N_eff_from_fraction(frac, power)
    return {
        "power": power,
        "L_over_R_H": frac,
        "N_eff": Ne,
        "sigma_0_eff": sigma_from_N(Ne),
        "L_metres": frac * R_H_M,
        "L_over_l_P": (frac * R_H_M) / L_P_M,
    }


def main() -> None:
    print("=" * 72)
    print("ROUTE 1: Local / meso-scale causal-set seed (redefine N_eff)")
    print("=" * 72)
    print(f"Global N_BH = {N_BH:.3e}  →  sigma_0_global = {SIGMA_0_GLOBAL:.3e}")
    print()

    print("--- N_eff required for target seeds (pure 1/sqrt(N) counting) ---")
    need_rows = []
    for name, sig in TARGETS.items():
        Ne = N_from_sigma(sig)
        print(f"  {name:16s}  sigma={sig:.1e}  ⇒  N_eff = {Ne:.3e}")
        need_rows.append({"target": name, "sigma": sig, "N_eff_required": Ne})

    # Scan L/R_H log-spaced deep enough that N_eff = N_BH f^p can hit 1e10
    fracs = np.logspace(-60, 0, 600)
    powers = [2, 3, 4]
    jobs = [(p, float(f)) for p in powers for f in fracs]

    print()
    print(f"Scanning {len(jobs)} (power, L/R_H) points on {os.cpu_count()} cores...")
    rows = []
    with ProcessPoolExecutor(max_workers=os.cpu_count() or 8) as ex:
        futs = [ex.submit(worker_row, j) for j in jobs]
        for fut in as_completed(futs):
            rows.append(fut.result())

    rows.sort(key=lambda r: (r["power"], r["L_over_R_H"]))

    # For each power, find L/R_H where sigma_0_eff crosses 1e-5 from below
    # (smaller L → smaller N_eff → larger sigma)
    print()
    print("--- L/R_H so that sigma_0_eff ≥ 1e-5 under N_eff = N_BH (L/R_H)^p ---")
    print("    (smaller patch ⇒ larger seed; these L are typically ≪ any meso scale)")
    for p in powers:
        subset = [r for r in rows if r["power"] == p and r["sigma_0_eff"] >= 1e-5]
        if not subset:
            print(f"  p={p}: not reached in scan")
            continue
        # largest L (most "macro") that still has sigma >= 1e-5
        best = max(subset, key=lambda r: r["L_over_R_H"])
        print(
            f"  p={p}: L/R_H ≲ {best['L_over_R_H']:.3e}  "
            f"(L ≈ {best['L_metres']:.3e} m ≪ R_H, N_eff ≈ {best['N_eff']:.3e}, "
            f"σ_eff ≈ {best['sigma_0_eff']:.3e})"
        )
        print(
            f"       → Under global-N rescaling, a 'local' seed of 1e-5 forces "
            f"sub-physical L; Route 1 needs a *different counting rule*, not a fraction of N_BH."
        )
    # Also: pure Planck counting without anchoring to N_BH
    print()
    print("--- Pure Planck counting N=(L/l_P)^p for σ_eff=1e-5 ---")
    Ne_need = N_from_sigma(1e-5)
    for p in powers:
        L_over_lp = Ne_need ** (1.0 / p)
        L = L_over_lp * L_P_M
        print(f"  p={p}: L/l_P = {L_over_lp:.3e}  ⇒  L = {L:.3e} m  "
              f"({L/R_H_M:.3e} R_H)")

    # CSV dumps
    path_grid = os.path.join(OUT, "route1_N_eff_grid.csv")
    with open(path_grid, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    path_need = os.path.join(OUT, "route1_N_eff_required.csv")
    with open(path_need, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(need_rows[0].keys()))
        w.writeheader()
        w.writerows(need_rows)

    print()
    print(f"Wrote {path_grid}")
    print(f"Wrote {path_need}")
    print()
    print("INTERPRETATION:")
    print("  Route 1 does NOT multiply sigma_0 by 1e56.")
    print("  It *redefines* the counting volume so the bare seed is larger.")
    print("  Price: justify local causal-set / meso-scale N_eff from first principles.")
    print("  Under N_eff = N_BH (L/R_H)^p with p=2 (area-like), L must be ~ R_H to")
    print("  recover global counting; smaller patches raise sigma_0_eff.")
    print("=" * 72)


if __name__ == "__main__":
    main()
