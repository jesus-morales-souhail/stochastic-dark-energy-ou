#!/usr/bin/env python3
"""
Sound horizon at the drag epoch for a given expansion history:

    r_d = int_{z_d}^{inf} c_s(z) / H(z) dz,   c_s = c / sqrt(3 (1 + R)),  R = (3 omega_b / 4 omega_gamma) / (1 + z)

- omega_b = 0.02237, omega_gamma for T0 = 2.7255 K.
- z_d: CAMB value at Planck 2018 (rd_planck.json, 1059.93), scaled with omega_m by the Eisenstein & Hu
  (1998) fitting formula, z_d(omega_m) = z_d^CAMB * EH(omega_m) / EH(omega_m0). z_d is not recomputed for
  the dark-energy history of each realisation (the expansion rate near z_d changes it slightly).
- The integral is normalised to CAMB: the plain integral with Omega_de = 0 gives r_d^int(omega_m0); every
  r_d is multiplied by r_d^CAMB / r_d^int(omega_m0) (one constant, recorded in rd_check). The validation
  of the unnormalised integral against CAMB and Aubourg is in main().
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad

HERE = Path(__file__).resolve().parent
RDJ = json.loads((HERE / "resultados" / "rd_planck.json").read_text())
C = 299792.458
OMB = 0.02237
OMG = 2.4728e-5  # photons, T0 = 2.7255 K
OMR = RDJ["omega_r"]
OM0 = RDJ["omega_m_incl_nu"]
ZD_CAMB = RDJ["zdrag"]
RD_CAMB = RDJ["rdrag_Mpc"]
RB = 3.0 * OMB / (4.0 * OMG)


def cs(z):
    return C / np.sqrt(3.0 * (1.0 + RB / (1.0 + np.asarray(z))))


def zd_eh(wm, wb=OMB):
    b1 = 0.313 * wm**-0.419 * (1 + 0.607 * wm**0.674)
    b2 = 0.238 * wm**0.223
    return 1291 * wm**0.251 / (1 + 0.659 * wm**0.828) * (1 + b1 * wb**b2)


def zd(wm):
    return ZD_CAMB * zd_eh(np.asarray(wm, dtype=float)) / zd_eh(OM0)


def rd_standard(wm: float, z_lo: float, z_hi: float = np.inf) -> float:
    """int_{z_lo}^{z_hi} c_s/H dz without dark energy (H = 100 sqrt(wm (1+z)^3 + wr (1+z)^4)), Mpc."""
    f = lambda lz: np.exp(lz) * cs(np.exp(lz) - 1) / (100 * np.sqrt(wm * np.exp(3 * lz) + OMR * np.exp(4 * lz)))
    hi = np.log1p(z_hi) if np.isfinite(z_hi) else np.log(1e12)
    return quad(f, np.log1p(z_lo), hi, epsabs=0, epsrel=1e-11, limit=400)[0]


NORM = RD_CAMB / rd_standard(OM0, ZD_CAMB)


if __name__ == "__main__":
    raw = rd_standard(OM0, ZD_CAMB)
    aub = 55.154 * np.exp(-72.3 * (RDJ["rd_aubourg"]["omega_nu"] + 0.0006) ** 2) / (
        RDJ["rd_aubourg"]["omega_cb"] ** 0.25351 * OMB**0.12807)
    print(f"r_d integral (Omega_de = 0, z_d CAMB) = {raw:.3f} Mpc; CAMB {RD_CAMB:.3f} ({raw / RD_CAMB - 1:+.4%}); "
          f"Aubourg {aub:.3f} ({raw / aub - 1:+.4%}); norm = {NORM:.6f}")
    print(f"z_d EH(om0) = {zd_eh(OM0):.2f}  vs CAMB {ZD_CAMB:.2f}")
