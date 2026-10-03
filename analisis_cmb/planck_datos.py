#!/usr/bin/env python3
"""Planck 2018 data used as error bars (see README for download)."""

from __future__ import annotations

from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
TT_FILE = HERE / "datos_planck" / "COM_PowerSpect_CMB-TT-full_R3.01.txt"
LENS_DIR = HERE / "datos_planck" / "cobaya_packages" / "data" / "planck_supp_data_and_covmats" / "lensing" / "2018"
LENS_BASE = "smicadx12_Dec5_ftl_mv2_ndclpp_p_teb_consext8"

# Planck 2018 VI (arXiv:1807.06209) Table 2, column TT,TE,EE+lowE+lensing
THETASTAR100 = 1.04110
THETASTAR100_ERR = 0.00031

FILES = {
    "TT": "COM_PowerSpect_CMB-TT-full_R3.01.txt",
    "pp": f"{LENS_BASE}_bandpowers.dat (planck_2018_lensing.native, conservative L=8-400)",
    "thetastar": "Planck 2018 VI Table 2 (TT,TE,EE+lowE+lensing): 100theta_* = 1.04110 +- 0.00031",
}


def load_tt(lmin: int = 2, lmax_excl: int = 30) -> dict:
    """D_l (muK^2) with asymmetric errors for lmin <= l < lmax_excl."""
    t = np.loadtxt(TT_FILE)
    ell = t[:, 0].astype(int)
    s = (ell >= lmin) & (ell < lmax_excl)
    return {"ell": ell[s], "Dl": t[s, 1], "err_lo": t[s, 2], "err_hi": t[s, 3]}


def load_lensing() -> dict:
    """Conservative MV bandpowers ([L(L+1)]^2 C_L/2pi) and their bin windows."""
    bp = np.loadtxt(LENS_DIR / f"{LENS_BASE}_bandpowers.dat")
    wins = []
    for b in range(1, len(bp) + 1):
        w = np.loadtxt(LENS_DIR / f"{LENS_BASE}_window" / f"window{b}.dat")
        wins.append((w[:, 0].astype(int), w[:, 1]))
    return {"L_av": bp[:, 3], "pp": bp[:, 4], "err": bp[:, 5], "windows": wins}


def bin_pp(clpp: np.ndarray, lens: dict) -> np.ndarray:
    """Bin [L(L+1)]^2 C_L/2pi theory with the Planck bin windows (no linear correction)."""
    return np.array([np.sum(w * clpp[L]) for L, w in lens["windows"]])
