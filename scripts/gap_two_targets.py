#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Two amplification gaps — do not mix labels.

  G_Euclid = 1e-5 / sigma_0     ~ 10^56
  G_DESI   = 2.5e-2 / sigma_0   ~ 10^57

sigma_0 = L_P / L_H (holographic Sorkin, H0=67.4).
See papers/amplification-gap.md §2.1
"""

from __future__ import annotations

import math

C_M_S = 299_792_458.0
MPC_M = 3.085_677_581e22
L_P_M = 1.616_255e-35
H0 = 67.4

SIGMA_EUCLID = 1e-5
SIGMA_DESI = 2.5e-2


def sorkin_sigma0(H0_km_s_mpc: float = H0) -> float:
    H0_si = (H0_km_s_mpc * 1e3) / MPC_M
    L_H = C_M_S / H0_si
    return L_P_M / L_H


def gap(target: float, sigma0: float) -> float:
    return target / sigma0


def r_soft(target: float, sigma0: float) -> float:
    return 0.5 * math.log(target / sigma0)


def main() -> None:
    s0 = sorkin_sigma0()
    g_e = gap(SIGMA_EUCLID, s0)
    g_d = gap(SIGMA_DESI, s0)
    print("Two amplification gaps (exact under holographic Sorkin)")
    print(f"  sigma_0 = L_P/L_H = {s0:.6e}  (H0={H0})")
    print()
    print(f"  G_Euclid = {SIGMA_EUCLID:.0e}/sigma_0 = {g_e:.3e}  ~ 10^{{{math.log10(g_e):.0f}}}")
    print(f"    soft r needed = {r_soft(SIGMA_EUCLID, s0):.2f}")
    print(f"  G_DESI   = {SIGMA_DESI:.1e}/sigma_0 = {g_d:.3e}  ~ 10^{{{math.log10(g_d):.0f}}}")
    print(f"    soft r needed = {r_soft(SIGMA_DESI, s0):.2f}")
    print()
    print(f"  Ratio G_DESI/G_Euclid = {g_d/g_e:.2f}  (exactly sigma_DESI/sigma_Euclid = {SIGMA_DESI/SIGMA_EUCLID:.1f})")
    print()
    print("Label rule:")
    print("  '~10^56' ONLY for Euclid-scale target ~1e-5")
    print("  '~10^59' for measured DESI profile ceiling 2.5e-2")
    print("Illegal: compute density with 2.5e-2 and call the DESI gap 10^56")
    print("Not claimed: a new value of rho_Lambda (imported cosmology, not derived here)")


if __name__ == "__main__":
    main()
