#!/usr/bin/env python3
"""
r_d and LCDM chi^2_BAO for Planck 2018 (ou_cmb.make_params, as in oscilador/), with CAMB, plus a
check of the analytic distance formula used in bao_sim.py. Run with analisis_cmb/venv/bin/python.
Writes resultados/rd_planck.json.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import quad

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import camb  # noqa: E402
import ou_cmb as M  # noqa: E402
from desi_dr2_data import load_gaussian_bao_vector  # noqa: E402  (local copy)

C = 299792.458


def rd_aubourg(omega_b: float, omega_cb, omega_nu: float):
    """Aubourg et al. 2015 (arXiv:1411.1074), Sec. II, eq. labelled eqn:rd; quoted accuracy 0.021 %."""
    return 55.154 * np.exp(-72.3 * (omega_nu + 0.0006) ** 2) / (np.asarray(omega_cb) ** 0.25351 * omega_b**0.12807)


def vec(DM, DH, bao, rd):
    out = []
    for z, lab in zip(bao["z_rows"], bao["labels"]):
        dm, dh = DM(z), DH(z)
        out.append({"DM": dm, "DH": dh, "DV": (z * dm * dm * dh) ** (1 / 3)}[lab[:2]] / rd)
    return np.array(out)


def main() -> None:
    p = M.make_params(halofit="takahashi")
    r = camb.get_background(p)
    rd = r.get_derived_params()["rdrag"]
    bao = load_gaussian_bao_vector()
    ci = np.linalg.inv(bao["cov"])
    v1 = vec(r.comoving_radial_distance, lambda z: C / r.hubble_parameter(z), bao, rd)
    om = p.ombh2 + p.omch2 + p.omnuh2
    orad = 2.4728e-5 * (1 + 0.22710731766 * 3.046)
    h = p.H0 / 100
    H = lambda z: 100 * np.sqrt(om * (1 + z) ** 3 + orad * (1 + z) ** 4 + (h * h - om - orad))
    v2 = vec(lambda z: quad(lambda x: C / H(x), 0, z, epsrel=1e-10)[0], lambda z: C / H(z), bao, rd)
    d1, d2 = v1 - bao["mean"], v2 - bao["mean"]
    out = {"rdrag_Mpc": rd, "zdrag": r.get_derived_params()["zdrag"], "H0": p.H0, "omega_m_incl_nu": om, "omega_r": orad,
           "chi2_lcdm_camb": float(d1 @ ci @ d1), "chi2_lcdm_analytic": float(d2 @ ci @ d2),
           "max_rel_diff_distances": float(np.max(np.abs(v2 / v1 - 1)))}
    rf = float(rd_aubourg(p.ombh2, p.ombh2 + p.omch2, p.omnuh2))
    out["rd_aubourg"] = {"omega_b": p.ombh2, "omega_cb": p.ombh2 + p.omch2, "omega_nu": p.omnuh2,
                         "rd_formula_Mpc": rf, "rel_diff_vs_camb": rf / rd - 1}
    (HERE / "resultados" / "rd_planck.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
