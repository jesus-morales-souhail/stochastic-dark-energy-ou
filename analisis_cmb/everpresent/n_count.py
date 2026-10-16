#!/usr/bin/env python3
"""
Everpresent Lambda: Planck-unit count N for today, and the naive amplitude 1/sqrt(N)
mapped to the repo's X = delta Omega_L (scripts/desi_dr2_data.py::add_ou_kernel).

Background: flat LCDM, Planck 2018 (TT,TE,EE+lowE+lensing) H0 = 67.36, ombh2 = 0.02237,
omch2 = 0.1200; radiation = photons (T0 = 2.7255 K) + 3.046 massless neutrinos.
The 0.06 eV neutrino mass is ignored (counted as matter would change N by << 1 %).

(1) 4-volume of the causal past of (here, now), Ahmed et al. eq. (3) with t = t0:
        V4 = (4 pi / 3) c * int_0^t0 dt a(t)^3 r_c(t)^3,
        r_c(t) = c int_t^t0 dt' / a(t')        (comoving radius of the past light cone at t)
    In terms of a, with E = H/H0 and R_H = c/H0:
        V4 = (4 pi / 3) R_H^4 int_0^1 da a^2 rhat(a)^3 / E(a),   rhat(a) = int_a^1 da'/(a'^2 E(a'))
    N_4 = V4 / l_P^4.
    Variant printed for reference only: Zwane et al. eq. (2.6) writes the inner integral as
    int_0^t' (particle horizon at t' instead of light-cone radius); not used in the table.
(2) N_area = A_H / (4 l_P^2), A_H = 4 pi R_H^2.
(3) delta Lambda / Lambda ~ 1/sqrt(N);  X = Omega_L * delta Lambda / Lambda;  ratio to 2.506e-2.
(4) Papers' own normalisation (Ahmed eq. (1)-(2), (5); Zwane eq. (2.3)): Delta rho_L ~ 1/sqrt(V)
    in units 8 pi G = hbar = c = 1, i.e. an absolute fluctuation. Relative to the observed
    rho_L0 = 3 Omega_L H0^2 (same units) this is (c/H0)^2 / (3 Omega_L sqrt(V4)), independent
    of l_P. Printed with unit coefficient, in reduced units and in G = c = hbar = 1 units
    (where rho_crit = 3H^2/8pi, so the same 1/sqrt(V) is 8 pi larger relative to rho_L0).
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import quad

C = 299_792_458.0  # m/s
MPC = 3.0856775814913673e22  # m
HBAR = 1.054571817e-34
G = 6.67430e-11
L_P = np.sqrt(HBAR * G / C**3)  # 1.616e-35 m
L_P_RED = np.sqrt(8 * np.pi * HBAR * G / C**3)

H0_KMS = 67.36
H = H0_KMS / 100.0
OMBH2, OMCH2 = 0.02237, 0.1200
T_CMB, N_EFF = 2.7255, 3.046
SIGMA_X = 0.025063155506929828  # analisis_cmb/ou_cmb.py

OMG_H2 = 2.4728e-5 * (T_CMB / 2.7255) ** 4  # photons
OMR = OMG_H2 * (1.0 + 0.22710731766 * N_EFF) / H**2
OMM = (OMBH2 + OMCH2) / H**2
OML = 1.0 - OMM - OMR

H0_SI = H0_KMS * 1e3 / MPC
R_H = C / H0_SI  # m


def E(a: float) -> float:
    return np.sqrt(OMR / a**4 + OMM / a**3 + OML)


def eta(a: float) -> float:
    """Conformal time int_0^a da'/(a'^2 E) in units of 1/H0 (c = 1)."""
    return quad(lambda u: 1.0 / (u * u * E(u)), 0.0, a, epsabs=0, epsrel=1e-11, limit=400)[0]


ETA0 = eta(1.0)


def v4_hat(inner) -> float:
    """(V4 / R_H^4) for a given comoving-radius function rhat(a)."""
    f = lambda la: np.exp(3 * la) * inner(np.exp(la)) ** 3 / E(np.exp(la))  # a^2 da = a^3 dln a
    return 4 * np.pi / 3 * quad(f, np.log(1e-12), 0.0, epsabs=0, epsrel=1e-9, limit=400)[0]


def main() -> None:
    v_cone = v4_hat(lambda a: ETA0 - eta(a))  # Ahmed eq. (3): past light cone
    v_zwane = v4_hat(eta)  # Zwane eq. (2.6) as printed

    V4 = v_cone * R_H**4
    N4 = V4 / L_P**4
    A_H = 4 * np.pi * R_H**2
    NA = A_H / (4 * L_P**2)

    print(f"Planck 2018: H0 = {H0_KMS}, Om = {OMM:.5f}, Or = {OMR:.3e}, OL = {OML:.5f}")
    print(f"R_H = c/H0 = {R_H:.4e} m,  l_P = {L_P:.6e} m")
    print(f"particle horizon today eta0 = {ETA0:.4f} c/H0 = {ETA0 * R_H / MPC:.1f} Mpc")
    print(f"V4 (past light cone, Ahmed eq.3) = {v_cone:.5f} (c/H0)^4 = {V4:.4e} m^4")
    print(f"   [Zwane eq.2.6 as printed: {v_zwane:.5f} (c/H0)^4 -- reference only]")
    print()

    rows = [("4-volumen V4/l_P^4", N4), ("area A_H/(4 l_P^2)", NA)]
    print(f"{'conteo':<22} {'N':>12} {'dL/L':>12} {'X equiv':>12} {'X/2.506e-2':>12}")
    for name, n in rows:
        d = 1 / np.sqrt(n)
        x = OML * d
        print(f"{name:<22} {n:12.3e} {d:12.3e} {x:12.3e} {x / SIGMA_X:12.3e}")
    print()

    # Papers' normalisation: absolute Delta rho_L ~ 1/sqrt(V) in 8piG = c = hbar = 1 units.
    rel = 1.0 / (3 * OML * np.sqrt(v_cone))  # Delta rho_L / rho_L0, l_P-independent
    print("Normalizacion de los papers (Delta rho_L ~ 1/sqrt(V) absoluto, coeficiente 1):")
    for lab, f in [("unidades reducidas (8piG=c=hbar=1)", 1.0), ("unidades sin reducir (G=c=hbar=1)", 8 * np.pi)]:
        r = rel * f
        x = OML * r
        print(f"  {lab:<34} dL/L0 = {r:.3f}  X = {x:.3f}  X/2.506e-2 = {x / SIGMA_X:.1f}")


if __name__ == "__main__":
    main()
