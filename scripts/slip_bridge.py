#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Stochastic-DE anisotropic stress → gravitational-slip bridge (Option 0).

SDiff / unimodular gravity cancels isotropic vacuum stress of the form
T_μν = V(x) g_μν. An anisotropic (shear) piece is not of that form, so the
trace-free projection need not remove it. If a fraction eps of the
BAO-bounded residual σ_X appears as shear, it sources gravitational slip
γ = Φ/Ψ ≠ 1.

This script puts order-of-magnitude numbers on that map. It is NOT a
Boltzmann solver (no hi_class / MGCAMB). Verdict: amplitude-starved.

Physics (Newtonian gauge, flat, sub-horizon k ≫ aH):
  Poisson   :  k² Ψ = −4π G a² ρ_m δ_m            (GR, μ=1)
  Anisotropy:  k² (Φ − Ψ) = 8π G a² π_T           (traceless ij)
  ⇒  (Φ − Ψ)/Ψ = −2 π_T / (ρ_m δ_m)

  π_T = eps · σ_X · ρ_X
  |γ − 1| = 2 · eps · σ_X · (ρ_X/ρ_m) / |δ_m|

Independent of (ρ+p) and of w: sidesteps the (1+w)→0 issue for a Λ-like
fluid because π_T is the raw shear amplitude.

Run:
  python scripts/slip_bridge.py

Related:
  papers/anisotropic-slip-option0.md
  papers/amplification-gap.md
  papers/resume.txt
  papers/data-pack-option0-internet.md
"""

from __future__ import annotations

import math

# ---- Cosmology (Planck / DESI-ish) ----
OM_M0 = 0.315
OM_L0 = 0.685

# ---- Repo bound (DESI DR2 BAO, Morales Souhail 2026) ----
SIGMA_X_95 = 1.5e-4  # 95% CL upper limit on fractional DE residual

# ---- Observational targets ----
# Maus et al. 2025 (arXiv:2505.20656): γ ≡ Φ/Ψ, GR = 1
GAMMA_OBS = 1.17
GAMMA_OBS_ERR = 0.11
DEV_OBS = abs(GAMMA_OBS - 1.0)  # 0.17

# Sakr et al. 2025 (arXiv:2501.07477): indicative floors on |η−1|
FORECAST_FLOOR_CONST = 0.05  # constant η
FORECAST_FLOOR_ZK = 0.30  # z,k-dependent η

# Matter contrast on measurement scales (sub-horizon, k ~ 0.1 h/Mpc)
DELTA_M = 1.0


def rho_ratio_x_over_m(z: float) -> float:
    """ρ_DE / ρ_m at redshift z (DE ~ cosmological constant)."""
    return OM_L0 / (OM_M0 * (1.0 + z) ** 3)


def slip_deviation(
    eps: float,
    sigma_x: float,
    z: float,
    delta_m: float = DELTA_M,
) -> float:
    """|γ − 1| from anisotropic fraction eps of stochastic σ_X."""
    return 2.0 * eps * sigma_x * rho_ratio_x_over_m(z) / abs(delta_m)


def eps_needed(target_dev: float, z: float, sigma_x: float = SIGMA_X_95) -> float:
    """eps required to reach target |γ−1| at fixed σ_X (can exceed 1)."""
    denom = 2.0 * sigma_x * rho_ratio_x_over_m(z)
    return target_dev / denom if denom > 0 else math.inf


def main() -> None:
    print("=" * 72)
    print("STOCHASTIC-DE → GRAVITATIONAL SLIP  (Option 0 bridge)")
    print("=" * 72)
    print(f"σ_X 95% CL bound (repo, DESI DR2 BAO) : {SIGMA_X_95:.2e}")
    print(
        f"Observed slip  γ = {GAMMA_OBS} ± {GAMMA_OBS_ERR}  "
        f"(Maus+ 2025, arXiv:2505.20656)  → |γ−1| = {DEV_OBS:.3f}"
    )
    print(
        f"Forecast floor on |η−1|: {FORECAST_FLOOR_CONST:.2f} (const), "
        f"{FORECAST_FLOOR_ZK:.2f} (z,k)  (Sakr+ 2025, arXiv:2501.07477)"
    )
    print(f"δ_m (fiducial linear rms)             : {DELTA_M:.1f}")
    print()

    for z in (0.5, 1.0, 1.5):
        r = rho_ratio_x_over_m(z)
        print(f"--- z = {z}   (ρ_X/ρ_m = {r:.3f}) ---")
        for eps in (1.0, 0.1, 0.01):
            dev = slip_deviation(eps, SIGMA_X_95, z)
            gap_obs = DEV_OBS / dev if dev > 0 else math.inf
            gap_fore = FORECAST_FLOOR_CONST / dev if dev > 0 else math.inf
            print(
                f"  eps={eps:<5} : |γ−1| = {dev:.3e}"
                f"   gap vs observed({DEV_OBS:.2f}): {gap_obs:8.0f}×  "
                f"vs forecast({FORECAST_FLOOR_CONST:.2f}): {gap_fore:8.0f}×"
            )
        print()

    print("=" * 72)
    print("BEST CASE: eps = 1 (all bounded residual is uncanceled shear)")
    print("=" * 72)
    for z in (0.5, 1.0, 1.5):
        dev_max = slip_deviation(1.0, SIGMA_X_95, z)
        eps_obs = eps_needed(DEV_OBS, z)
        eps_fore = eps_needed(FORECAST_FLOOR_CONST, z)
        print(f"  z={z}: max |γ−1| = {dev_max:.3e}")
        print(
            f"         to reach observed {DEV_OBS:.2f} need eps = {eps_obs:.0f}  "
            f"({'impossible' if eps_obs > 1 else 'allowed'})"
        )
        print(
            f"         to reach forecast {FORECAST_FLOOR_CONST:.2f} need eps = {eps_fore:.0f}  "
            f"({'impossible' if eps_fore > 1 else 'allowed'})"
        )
    print()

    # Robustness: δ_m bracket
    print("=" * 72)
    print("ROBUSTNESS: δ_m bracket at z=0.5, eps=1")
    print("=" * 72)
    for dm in (0.5, 1.0, 2.0):
        dev = slip_deviation(1.0, SIGMA_X_95, 0.5, delta_m=dm)
        print(f"  δ_m={dm:<4} → |γ−1| = {dev:.3e}  (gap vs 0.05: {FORECAST_FLOOR_CONST/dev:.0f}×)")
    print()

    print("VERDICT:")
    print("  The SDiff 'grieta' (anisotropic stress is not canceled) is real")
    print("  in principle, but σ_X < 1.5×10⁻⁴ caps the mechanism. Even eps=1")
    print("  gives |γ−1| ~ 10⁻⁴ — O(10²–10³)× below Maus/Sakr targets.")
    print("  Observable slip requires shear DECOUPLED from the BAO-bounded")
    print("  isotropic residual (or amplification of A0). That is no longer")
    print("  this stochastic mechanism alone — it is generic anisotropic")
    print("  stress / MG (Clifton+Ferreira+Padilla+Skordis degeneracy).")
    print()
    print("  Act IV inherits Act III: amplification gap on amplitude (no free gain).")


if __name__ == "__main__":
    main()
