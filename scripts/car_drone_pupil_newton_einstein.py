#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Car + drone @ 120 km/h + diffraction in a 1 mm pupil:
Newton/Galileo vs Einstein (SR) + local Hubble expansion.

No tesseract. Numbers only.

Run:
    python scripts/car_drone_pupil_newton_einstein.py
"""

from __future__ import annotations

import numpy as np

c = 299_792_458.0  # m/s


def main() -> None:
    v = 120_000 / 3600  # m/s
    print("=" * 72)
    print("EXPERIMENT: car + drone @ 120 km/h, diffraction in a pupil")
    print("=" * 72)
    print(f"v_car = v_drone = {v:.6f} m/s = {v/c:.3e} c")
    print("v_rel (car-drone) = 0  (same velocity vector)")
    print()

    H0 = 70e3 / (3.085677581e22)  # 70 km/s/Mpc -> 1/s
    print("--- A) Expansion of space vs EM rigidity ---")
    print(f"H0 ≈ {H0:.3e} s^-1")
    for L in [1e-10, 1e-9, 1e-6, 1.0, 1e6]:
        print(f"  Hubble recession over L={L:.0e} m: v_H = H0*L = {H0*L:.3e} m/s")
    print("  Bound systems (atoms, car, eye) do NOT participate in Hubble flow.")
    print()

    print("--- B) Diffraction in the pupil (car rest frame) ---")
    lam = 550e-9
    f_eye = 0.017
    for D in [1e-3, 4e-3]:
        theta = 1.22 * lam / D
        r_airy = 1.22 * lam * f_eye / D
        print(f"  D = {D*1e3:.1f} mm, λ = {lam*1e9:.0f} nm:")
        print(f"    θ_Airy ≈ {theta:.3e} rad = {np.degrees(theta)*3600:.3f} arcsec")
        print(f"    Airy radius (f={f_eye*1e3:.0f} mm) ≈ {r_airy*1e6:.3f} μm")
    print()

    print("--- C) Newton / Galileo (inertial, constant v) ---")
    print("  Passenger: inertial ⇒ same Airy pattern as parked car.")
    print("  Drone co-moving: car at rest relative to drone ⇒ identical pattern.")
    print("  No 4D optical 'tesseract' projection appears.")
    print()

    print("--- D) Einstein special relativity ---")
    beta = v / c
    gamma = 1 / np.sqrt(1 - beta**2)
    D = 1e-3
    print(f"  β = {beta:.3e},  γ-1 = {gamma-1:.3e}")
    print(f"  Light-crossing time of 1 mm: D/c = {D/c:.3e} s")
    print("  Correct setup (matched velocities): relative β=0")
    print("  ⇒ no contraction/Doppler between car and drone; same diffraction.")
    print("  Wrong setup (drone on road):")
    print(f"    |Δf/f|_transverse ≈ {abs(1/gamma - 1):.3e}")
    print(f"    aberration ~ β ≈ {beta:.3e} rad = {np.degrees(beta)*3600:.4f} arcsec")
    theta_A = 1.22 * lam / 1e-3
    print(f"    Airy angle ≈ {theta_A:.3e} rad; ratio β/θ_Airy ≈ {beta/theta_A:.3e}")
    print()

    print("--- E) When SR would matter for the optical pattern ---")
    v_need = theta_A * c
    print(f"  Aberration ~ Airy (1 mm): v ~ {v_need/1000:.1f} km/s "
          f"({v_need/v:.0f}× 120 km/h)")
    v_01 = c * np.sqrt(1 - 1 / 1.01**2)
    print(f"  For γ-1 ~ 1%: v ≈ {v_01/c:.3f} c ≈ {v_01/1000:.0f} km/s")
    print()

    print("--- F) Hubble stretch during 60 s ---")
    t_exp = 60.0
    for L in [1.0, 1e-3, 1e-10]:
        print(f"  ΔL(L={L:.0e} m) = H0 L t = {H0*L*t_exp:.3e} m")
    print()

    print("--- G) 'Mechanism blinds the question' ---")
    print("  Local EM + local inertial physics shield the pupil experiment from")
    print("  cosmology. Expansion is measured with free photons over Gpc (redshift),")
    print("  not with a 1 mm pupil at 120 km/h — wrong scale / wrong operator.")
    print("=" * 72)


if __name__ == "__main__":
    main()
