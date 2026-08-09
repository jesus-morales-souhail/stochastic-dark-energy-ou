#!/usr/bin/env python3
"""
Load DESI DR2 BAO data from the local Zenodo pack (real files).

Public products (CobayaSampler/bao_data + DESI DR2 Results II, arXiv:2503.14738):
  - Gaussian BAO mean vector (13 distances: DV, DM, DH per tracer bin)
  - Full 13×13 covariance (block-diagonal across tracers; DM–DH correlated within bin)
  - Isotropic α from figure6 and its projected 7×7 covariance

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
# Small tracked copies for pipelines that must not depend on the huge local pack path
COV_PUBLIC = ROOT / "results" / "desi_cov"
COV_PUBLIC.mkdir(parents=True, exist_ok=True)

# Published isotropic alpha fallback (DESI DR2 public summary)
_FALLBACK_Z = np.array([0.295, 0.510, 0.706, 0.934, 1.321, 1.484, 2.330])
_FALLBACK_A = np.array([1.0030, 0.9947, 1.0016, 0.9960, 1.0020, 0.9963, 1.0008])
_FALLBACK_S = np.array([0.0097, 0.0072, 0.0057, 0.0049, 0.0063, 0.0088, 0.0120])
_FALLBACK_S_Z = np.array([-0.284, -0.462, -0.595, -0.719, -0.870, -0.917, -1.070])

# Labels for the official 13-vector (ALL_GCcomb mean file order)
_GAUSS_LABELS = [
    "DV_BGS_z0.295",
    "DM_LRG_z0.510",
    "DH_LRG_z0.510",
    "DM_LRG_z0.706",
    "DH_LRG_z0.706",
    "DM_LRG+ELG_z0.934",
    "DH_LRG+ELG_z0.934",
    "DM_ELG_z1.321",
    "DH_ELG_z1.321",
    "DM_QSO_z1.484",
    "DH_QSO_z1.484",
    "DH_Lya_z2.33",
    "DM_Lya_z2.33",
]


def multipole_dir() -> Path | None:
    return FIG5 if FIG5.is_dir() else None


def bao_dr2_dir() -> Path | None:
    return BAO_DR2 if BAO_DR2.is_dir() else None


def load_gaussian_bao_vector() -> dict:
    """
    Official DESI DR2 Gaussian BAO mean + full 13×13 covariance.

    Source: desi_gaussian_bao_ALL_GCcomb_{mean,cov}.txt in the Zenodo / Cobaya bao_data pack.
    Structure: block-diagonal across tracer bins; within each anisotropic bin, DM and DH
    are anti-correlated (typical |ρ| ~ 0.35–0.50). Cross-bin blocks are zero in this public
    likelihood (DESI DR2 Results II compressed Gaussian BAO).
    """
    mean_path = BAO_DR2 / "desi_gaussian_bao_ALL_GCcomb_mean.txt"
    cov_path = BAO_DR2 / "desi_gaussian_bao_ALL_GCcomb_cov.txt"
    if not (mean_path.is_file() and cov_path.is_file()):
        raise FileNotFoundError(
            "Missing DESI Gaussian BAO mean/cov under data/.../desi_bao_dr2/. "
            "Unpack the Zenodo dr2-bao-zenodo pack or clone CobayaSampler/bao_data."
        )
    # columns: z, value, quantity_name (string) — load only numeric cols
    mean_table = np.loadtxt(mean_path, comments="#", usecols=(0, 1))
    z_rows = mean_table[:, 0]
    mean = mean_table[:, 1].astype(float)
    cov = np.loadtxt(cov_path).astype(float)
    if cov.shape != (13, 13) or mean.shape != (13,):
        raise ValueError(f"Unexpected shapes mean={mean.shape} cov={cov.shape}")
    cov = 0.5 * (cov + cov.T)
    # persist a small public copy
    np.savetxt(
        COV_PUBLIC / "desi_gaussian_bao_ALL_GCcomb_cov.txt",
        cov,
        header="DESI DR2 Gaussian BAO 13x13 cov (public; block-diagonal across bins)",
    )
    return {
        "mean": mean,
        "cov": cov,
        "z_rows": z_rows.astype(float),
        "labels": list(_GAUSS_LABELS),
        "source": str(mean_path.relative_to(ROOT)),
        "cov_source": str(cov_path.relative_to(ROOT)),
        "mock": False,
        "note": (
            "Full public Gaussian BAO covariance. Block-diagonal across redshift tracers; "
            "DM–DH correlated within anisotropic bins. See arXiv:2503.14738 and CobayaSampler/bao_data."
        ),
    }


def project_alpha_covariance(alphas: np.ndarray, mean: np.ndarray, cov13: np.ndarray) -> np.ndarray:
    """
    Project the 13-distance Gaussian BAO covariance onto isotropic α via Jacobian.

    α_iso ≈ DV/DV_fid for BGS; for anisotropic bins
    α_iso ∝ DM^{2/3} DH^{1/3}, so
      dα/α = (2/3) dDM/DM + (1/3) dDH/DH.

    Because the public 13×13 is block-diagonal across bins, the projected 7×7 α
    covariance is diagonal across bins (off-block zeros). Within-bin DM–DH
    correlation is folded into each α variance.
    """
    alphas = np.asarray(alphas, dtype=float)
    assert alphas.shape == (7,)
    assert mean.shape == (13,) and cov13.shape == (13, 13)

    J = np.zeros((7, 13))
    # BGS DV
    J[0, 0] = alphas[0] / mean[0]
    # anisotropic pairs (DM, DH) indices
    pairs = [(1, 1, 2), (2, 3, 4), (3, 5, 6), (4, 7, 8), (5, 9, 10)]
    for i_alpha, i_dm, i_dh in pairs:
        J[i_alpha, i_dm] = (2.0 / 3.0) * alphas[i_alpha] / mean[i_dm]
        J[i_alpha, i_dh] = (1.0 / 3.0) * alphas[i_alpha] / mean[i_dh]
    # Lya: mean file order is DH then DM
    J[6, 11] = (1.0 / 3.0) * alphas[6] / mean[11]
    J[6, 12] = (2.0 / 3.0) * alphas[6] / mean[12]

    C = J @ cov13 @ J.T
    return 0.5 * (C + C.T)


def load_alpha_dv(prefer_file: bool = True, use_full_cov: bool = True) -> dict:
    """
    Isotropic alpha = D_V / r_s relative to fiducial (7 bins).

    If use_full_cov and the Gaussian BAO pack is present, attach
    cov (7×7) from the official 13×13 via Jacobian projection.
    """
    path = FIG6 / "DESI_DR2_alpha_DV_over_rs.txt"
    if prefer_file and path.is_file():
        raw = np.loadtxt(path, comments="#")
        z, a, s = raw[:, 0], raw[:, 1], raw[:, 2]
        if len(z) == len(_FALLBACK_S_Z):
            sz = _FALLBACK_S_Z.copy()
        else:
            sz = -1.0 / (1.0 + z)
            sz = sz / np.sqrt(np.mean(sz**2))
        out = {
            "z": z.astype(float),
            "alpha": a.astype(float),
            "sigma": s.astype(float),
            "S_z": sz.astype(float),
            "source": str(path.relative_to(ROOT)),
            "mock": False,
            "cov_mode": "diagonal_figure6_errors",
        }
    else:
        out = {
            "z": _FALLBACK_Z.copy(),
            "alpha": _FALLBACK_A.copy(),
            "sigma": _FALLBACK_S.copy(),
            "S_z": _FALLBACK_S_Z.copy(),
            "source": "published_DR2_isotropic_alpha_table (file missing)",
            "mock": False,
            "cov_mode": "diagonal_fallback",
        }

    # Full measurement covariance
    C = np.diag(out["sigma"] ** 2)
    cov_mode = out["cov_mode"]
    if use_full_cov and BAO_DR2.is_dir():
        try:
            g = load_gaussian_bao_vector()
            C_proj = project_alpha_covariance(out["alpha"], g["mean"], g["cov"])
            # Prefer projected cov; keep figure6 sigma as cross-check
            # If projection fails PD, fall back
            eig = np.linalg.eigvalsh(C_proj)
            if np.all(eig > 0):
                C = C_proj
                out["sigma"] = np.sqrt(np.diag(C))
                cov_mode = "full_gaussian_bao_projected_13x13"
                out["gaussian_bao_source"] = g["source"]
                out["cov_note"] = g["note"]
                # save public 7x7
                np.savetxt(
                    COV_PUBLIC / "desi_cov_alpha_iso_7x7.txt",
                    C,
                    header=(
                        "DESI DR2 isotropic alpha 7x7 cov = J C_13 J^T "
                        "from desi_gaussian_bao_ALL_GCcomb (block-diagonal across bins)"
                    ),
                )
                np.save(BAO_DR2 / "desi_cov_alpha_iso.npy", C)
        except Exception as exc:  # noqa: BLE001
            out["cov_fallback_reason"] = str(exc)

    out["cov"] = C
    out["cov_mode"] = cov_mode
    return out


def measurement_cov(data: dict | None = None) -> np.ndarray:
    """Return 7×7 measurement covariance for residual pipelines."""
    if data is None:
        data = load_alpha_dv(prefer_file=True, use_full_cov=True)
    if "cov" in data and data["cov"] is not None:
        return np.asarray(data["cov"], dtype=float)
    return np.diag(np.asarray(data["sigma"], dtype=float) ** 2)


def add_ou_kernel(C_meas: np.ndarray, S_z: np.ndarray, x: np.ndarray, theta: float, sigma_X: float) -> np.ndarray:
    """C_total = C_meas + OU residual signal."""
    C = np.array(C_meas, dtype=float, copy=True)
    s2 = float(sigma_X) ** 2
    n = len(S_z)
    for i in range(n):
        for j in range(n):
            dx = abs(float(x[i]) - float(x[j]))
            C[i, j] += S_z[i] * S_z[j] * s2 * np.exp(-theta * dx)
    return C


def logL_gaussian(residuals: np.ndarray, cov: np.ndarray) -> float:
    """Multivariate Gaussian log-likelihood."""
    from scipy.linalg import cho_factor, cho_solve

    r = np.asarray(residuals, dtype=float)
    C = np.asarray(cov, dtype=float)
    n = len(r)
    if not np.all(np.isfinite(C)):
        return -1e30
    try:
        c, lower = cho_factor(C, lower=True, check_finite=False)
        y = cho_solve((c, lower), r, check_finite=False)
        logdet = 2.0 * np.sum(np.log(np.diag(c)))
        return float(-0.5 * (r @ y + logdet + n * np.log(2.0 * np.pi)))
    except Exception:
        return -1e30
