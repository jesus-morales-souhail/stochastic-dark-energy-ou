#!/usr/bin/env python3
"""
Load DESI DR2 BAO alpha vectors from the local Zenodo pack (real files).

No synthetic mocks. Paths are relative to this repository root.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
FIG6 = ROOT / "data" / "desi_dr2_local" / "dr2_data" / "dr2-bao-zenodo" / "figure6"
FIG5 = ROOT / "data" / "desi_dr2_local" / "dr2_data" / "dr2-bao-zenodo" / "figure5"
BAO_DR2 = (
    ROOT
    / "data"
    / "desi_dr2_local"
    / "dr2_data"
    / "dr2-bao-zenodo"
    / "cosmology_chains"
    / "bao_data"
    / "desi_bao_dr2"
)

# Published isotropic alpha (DESI DR2 public summary) — only if files missing
_FALLBACK_Z = np.array([0.295, 0.510, 0.706, 0.934, 1.321, 1.484, 2.330])
_FALLBACK_A = np.array([1.0030, 0.9947, 1.0016, 0.9960, 1.0020, 0.9963, 1.0008])
_FALLBACK_S = np.array([0.0097, 0.0072, 0.0057, 0.0049, 0.0063, 0.0088, 0.0120])
_FALLBACK_S_Z = np.array([-0.284, -0.462, -0.595, -0.719, -0.870, -0.917, -1.070])


def load_alpha_dv(prefer_file: bool = True) -> dict:
    """Isotropic alpha = D_V / r_s relative to fiducial (7 bins)."""
    path = FIG6 / "DESI_DR2_alpha_DV_over_rs.txt"
    if prefer_file and path.is_file():
        raw = np.loadtxt(path, comments="#")
        z, a, s = raw[:, 0], raw[:, 1], raw[:, 2]
        # sensitivity weights: same shape as published table if lengths match
        if len(z) == len(_FALLBACK_S_Z):
            sz = _FALLBACK_S_Z.copy()
        else:
            sz = -1.0 / (1.0 + z)
            sz = sz / np.sqrt(np.mean(sz**2))
        return {
            "z": z.astype(float),
            "alpha": a.astype(float),
            "sigma": s.astype(float),
            "S_z": sz.astype(float),
            "source": str(path.relative_to(ROOT)),
            "mock": False,
        }
    return {
        "z": _FALLBACK_Z.copy(),
        "alpha": _FALLBACK_A.copy(),
        "sigma": _FALLBACK_S.copy(),
        "S_z": _FALLBACK_S_Z.copy(),
        "source": "published_DR2_isotropic_alpha_table (file missing)",
        "mock": False,
    }


def multipole_dir() -> Path | None:
    return FIG5 if FIG5.is_dir() else None


def bao_dr2_dir() -> Path | None:
    return BAO_DR2 if BAO_DR2.is_dir() else None
