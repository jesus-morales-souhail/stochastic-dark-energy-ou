#!/usr/bin/env python3
"""
CPL background and optional nested effective-EoS extension on DESI DR2 BAO.

Fits free (w0, wa) to the public DESI DR2 BAO alpha vector (diagonal errors),
then a nested extension with mean-field amplitude σ and rate θ. Uses multi-start
Nelder–Mead (derivative-free). Outputs go to ``results/eos_cpl_desi_dr2/``.

Data: same BAO summary statistics as ``ou_bao_likelihood.py`` (arXiv:2503.14738).

Author: Jesús Morales Souhail
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize

# ---------------------------------------------------------------------------
# Real DESI DR2 BAO data (repo standard)
# ---------------------------------------------------------------------------
z_eff = np.array([0.295, 0.510, 0.706, 0.934, 1.321, 1.484, 2.330])
alpha_obs = np.array([1.0030, 0.9947, 1.0016, 0.9960, 1.0020, 0.9963, 1.0008])
sigma_obs = np.array([0.0097, 0.0072, 0.0057, 0.0049, 0.0063, 0.0088, 0.0120])

H0 = 67.4
Om = 0.315
c_kms = 299792.458

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "results" / "eos_cpl_desi_dr2"
OUT_DIR.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Pure CPL cosmology (standard textbook expressions)
# ---------------------------------------------------------------------------
def E2_cpl(z: float, w0: float, wa: float, Om_: float = Om) -> float:
    """E(z)^2 for flat CPL dark energy."""
    rho_de = (1.0 + z) ** (3.0 * (1.0 + w0 + wa)) * np.exp(-3.0 * wa * z / (1.0 + z))
    return Om_ * (1.0 + z) ** 3 + (1.0 - Om_) * rho_de


def H_cpl(z: float, w0: float, wa: float) -> float:
    e2 = E2_cpl(z, w0, wa)
    if e2 <= 0.0 or not np.isfinite(e2):
        return np.nan
    return H0 * np.sqrt(e2)


def DM_cpl(z: float, w0: float, wa: float) -> float:
    def integrand(zp: float) -> float:
        h = H_cpl(zp, w0, wa)
        if not np.isfinite(h) or h <= 0.0:
            raise ValueError("non-positive H")
        return 1.0 / h

    return c_kms * quad(integrand, 0.0, z, limit=120, epsabs=1e-6)[0]


def DH_cpl(z: float, w0: float, wa: float) -> float:
    h = H_cpl(z, w0, wa)
    if not np.isfinite(h) or h <= 0.0:
        raise ValueError("non-positive H")
    return c_kms / h


def DV_cpl(z: float, w0: float, wa: float) -> float:
    dm = DM_cpl(z, w0, wa)
    dh = DH_cpl(z, w0, wa)
    return (dm**2 * dh * z) ** (1.0 / 3.0)


# Fiducial = flat ΛCDM (w0=-1, wa=0), matching alpha≈1 definition in the repo
W0_FID, WA_FID = -1.0, 0.0
DV_fid = DV_cpl(float(z_eff[0]), W0_FID, WA_FID)
DM_fid = np.array([DM_cpl(float(z), W0_FID, WA_FID) for z in z_eff[1:]])
DH_fid = np.array([DH_cpl(float(z), W0_FID, WA_FID) for z in z_eff[1:]])


def alpha_pred_cpl(w0: float, wa: float) -> np.ndarray:
    """Model alpha relative to flat-ΛCDM fiducial (same construction as DR2 tests)."""
    out = np.zeros(len(z_eff))
    out[0] = DV_cpl(float(z_eff[0]), w0, wa) / DV_fid
    for i in range(1, len(z_eff)):
        z = float(z_eff[i])
        dm = DM_cpl(z, w0, wa)
        dh = DH_cpl(z, w0, wa)
        out[i] = (dm / DM_fid[i - 1]) ** (2.0 / 3.0) * (dh / DH_fid[i - 1]) ** (1.0 / 3.0)
    return out


def chi2_cpl(w0: float, wa: float) -> float:
    pred = alpha_pred_cpl(w0, wa)
    if not np.all(np.isfinite(pred)):
        return np.inf
    return float(np.sum(((alpha_obs - pred) / sigma_obs) ** 2))


def logL_cpl(w0: float, wa: float) -> float:
    c2 = chi2_cpl(w0, wa)
    if not np.isfinite(c2):
        return -np.inf
    return -0.5 * c2


# ---------------------------------------------------------------------------
# Optional OU-style additive correction on w (nested extension)
# w_eff(z) = w_CPL(z) + δw(z),  δw ~ σ²/(2θ) (1 - e^{-2θ t}) / Ω_DE  (toy)
# No hard phantom clamp — that created a flat ΛCDM plateau.
# ---------------------------------------------------------------------------
def conformal_time_to_z(z: float) -> float:
    def integrand(zp: float) -> float:
        # background for time coordinate: flat ΛCDM
        h = H_cpl(zp, W0_FID, WA_FID)
        return 1.0 / ((1.0 + zp) * h)

    return quad(integrand, 0.0, z, limit=80, epsabs=1e-6)[0]


def Omega_DE_eff(z: float, w0: float, wa: float, sigma: float, theta: float) -> float:
    """Integrate continuity equation for DE with optional mean-field OU correction."""

    def integrand(zp: float) -> float:
        w_cpl = w0 + wa * zp / (1.0 + zp)
        if sigma <= 0.0:
            w_eff = w_cpl
        else:
            tz = conformal_time_to_z(zp)
            # background DE fraction (ΛCDM) only as a scale in the toy correction
            e2_bg = E2_cpl(zp, W0_FID, WA_FID)
            omega_bg = (1.0 - Om) / e2_bg
            if theta < 1e-8:
                corr = (sigma**2) * tz / max(omega_bg, 1e-12)
            else:
                corr = (sigma**2 / (2.0 * theta * max(omega_bg, 1e-12))) * (
                    1.0 - np.exp(-2.0 * theta * tz)
                )
            w_eff = w_cpl + corr
        return (1.0 + w_eff) / (1.0 + zp)

    exponent = 3.0 * quad(integrand, 0.0, z, limit=100, epsabs=1e-6)[0]
    return (1.0 - Om) * np.exp(exponent)


def E2_eff(z: float, w0: float, wa: float, sigma: float, theta: float) -> float:
    return Om * (1.0 + z) ** 3 + Omega_DE_eff(z, w0, wa, sigma, theta)


def H_eff(z: float, w0: float, wa: float, sigma: float, theta: float) -> float:
    e2 = E2_eff(z, w0, wa, sigma, theta)
    if e2 <= 0.0 or not np.isfinite(e2):
        return np.nan
    return H0 * np.sqrt(e2)


def alpha_pred_eff(w0: float, wa: float, sigma: float, theta: float) -> np.ndarray:
    def DM(z: float) -> float:
        def integrand(zp: float) -> float:
            h = H_eff(zp, w0, wa, sigma, theta)
            if not np.isfinite(h) or h <= 0.0:
                raise ValueError("non-positive H_eff")
            return 1.0 / h

        return c_kms * quad(integrand, 0.0, z, limit=100, epsabs=1e-6)[0]

    def DH(z: float) -> float:
        h = H_eff(z, w0, wa, sigma, theta)
        if not np.isfinite(h) or h <= 0.0:
            raise ValueError("non-positive H_eff")
        return c_kms / h

    out = np.zeros(len(z_eff))
    z0 = float(z_eff[0])
    dm0 = DM(z0)
    dh0 = DH(z0)
    out[0] = (dm0**2 * dh0 * z0) ** (1.0 / 3.0) / DV_fid
    for i in range(1, len(z_eff)):
        z = float(z_eff[i])
        out[i] = (DM(z) / DM_fid[i - 1]) ** (2.0 / 3.0) * (DH(z) / DH_fid[i - 1]) ** (
            1.0 / 3.0
        )
    return out


def chi2_eff(w0: float, wa: float, sigma: float, theta: float) -> float:
    try:
        pred = alpha_pred_eff(w0, wa, sigma, theta)
    except (ValueError, FloatingPointError, ZeroDivisionError, OverflowError):
        return np.inf
    if not np.all(np.isfinite(pred)):
        return np.inf
    return float(np.sum(((alpha_obs - pred) / sigma_obs) ** 2))


def logL_eff(w0: float, wa: float, sigma: float, theta: float) -> float:
    c2 = chi2_eff(w0, wa, sigma, theta)
    if not np.isfinite(c2):
        return -np.inf
    return -0.5 * c2


# ---------------------------------------------------------------------------
# Optimizers (derivative-free; multi-start)
# ---------------------------------------------------------------------------
def fit_cpl_pure(starts=None):
    if starts is None:
        starts = [
            [-1.0, 0.0],
            [-0.9, -0.3],
            [-0.87, -0.41],
            [-0.7, -0.5],
            [-1.2, 0.2],
            [-0.5, -1.0],
        ]

    def neg(p):
        w0, wa = p
        if not (-3.0 < w0 < 0.5 and -5.0 < wa < 5.0):
            return 1e12
        ll = logL_cpl(w0, wa)
        return -ll if np.isfinite(ll) else 1e12

    best = None
    for x0 in starts:
        res = minimize(
            neg,
            x0,
            method="Nelder-Mead",
            options={"xatol": 1e-7, "fatol": 1e-7, "maxiter": 4000},
        )
        if best is None or res.fun < best.fun:
            best = res
    w0, wa = best.x
    return {
        "w0": float(w0),
        "wa": float(wa),
        "chi2": float(chi2_cpl(w0, wa)),
        "logL": float(logL_cpl(w0, wa)),
        "success": bool(best.success),
        "nfev": int(best.nfev),
    }


def fit_eff_nested(starts=None):
    """Joint {w0, wa, sigma, theta}; sigma>=0, theta>0."""
    if starts is None:
        starts = [
            [-1.0, 0.0, 1e-5, 0.1],
            [-0.99, -0.02, 1e-4, 0.5],
            [-0.87, -0.41, 1e-3, 1.0],
            [-0.9, -0.2, 5e-3, 0.2],
            [-1.1, 0.1, 1e-6, 1.0],
        ]

    def neg(p):
        w0, wa, sig, th = p
        if not (-3.0 < w0 < 0.5 and -5.0 < wa < 5.0):
            return 1e12
        if sig < 0.0 or th < 1e-6 or th > 50.0 or sig > 0.05:
            return 1e12
        ll = logL_eff(w0, wa, sig, th)
        return -ll if np.isfinite(ll) else 1e12

    best = None
    for x0 in starts:
        res = minimize(
            neg,
            x0,
            method="Nelder-Mead",
            options={"xatol": 1e-7, "fatol": 1e-7, "maxiter": 6000},
        )
        if best is None or res.fun < best.fun:
            best = res
    w0, wa, sig, th = best.x
    return {
        "w0": float(w0),
        "wa": float(wa),
        "sigma": float(sig),
        "theta": float(th),
        "chi2": float(chi2_eff(w0, wa, sig, th)),
        "logL": float(logL_eff(w0, wa, sig, th)),
        "success": bool(best.success),
        "nfev": int(best.nfev),
    }


def aic(logL: float, k: int) -> float:
    return 2 * k - 2 * logL


def bic(logL: float, k: int, n: int) -> float:
    return k * np.log(n) - 2 * logL


def main():
    print("=" * 70)
    print("DESI DR2 BAO — CPL background and nested effective-EoS extension")
    print("=" * 70)
    print(f"N bins = {len(z_eff)}  |  fiducial = flat ΛCDM (w0=-1, wa=0)")
    print()

    ll_lcdm = logL_cpl(-1.0, 0.0)
    print("Reference ΛCDM (w0=-1, wa=0):")
    print(f"  chi2 = {chi2_cpl(-1.0, 0.0):.4f}   logL = {ll_lcdm:.4f}")
    print()

    cpl = fit_cpl_pure()
    print("Pure CPL MLE (multi-start Nelder–Mead):")
    print(f"  w0   = {cpl['w0']:+.4f}")
    print(f"  wa   = {cpl['wa']:+.4f}")
    print(f"  chi2 = {cpl['chi2']:.4f}")
    print(f"  logL = {cpl['logL']:.4f}")
    print(f"  AIC  = {aic(cpl['logL'], 2):.3f}   BIC = {bic(cpl['logL'], 2, 7):.3f}")
    ll_check = logL_cpl(cpl["w0"], cpl["wa"])
    if abs(ll_check - cpl["logL"]) > 1e-6:
        print(f"  WARNING: re-evaluated logL differs ({ll_check:.4f})")
    print()

    eff = fit_eff_nested()
    print("Nested extension (CPL + σ, θ):")
    print(f"  w0    = {eff['w0']:+.4f}")
    print(f"  wa    = {eff['wa']:+.4f}")
    print(f"  sigma = {eff['sigma']:.3e}")
    print(f"  theta = {eff['theta']:.4f}")
    print(f"  chi2  = {eff['chi2']:.4f}")
    print(f"  logL  = {eff['logL']:.4f}")
    print(f"  AIC   = {aic(eff['logL'], 4):.3f}   BIC = {bic(eff['logL'], 4, 7):.3f}")
    print()

    dlogL = eff["logL"] - cpl["logL"]
    daic = aic(eff["logL"], 4) - aic(cpl["logL"], 2)
    print("Nested comparison:")
    print(f"  ΔlogL (extension − CPL) = {dlogL:+.4f}")
    print(f"  ΔAIC  (extension − CPL) = {daic:+.3f}")
    if daic > 0:
        print("  → Pure CPL preferred (no improvement from σ, θ).")
    print()

    print("Profile in σ at best pure-CPL (w0, wa), θ = 1:")
    w0b, wab = cpl["w0"], cpl["wa"]
    ll0 = logL_eff(w0b, wab, 0.0, 1.0)
    for sig in [0.0, 1e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2]:
        ll = logL_eff(w0b, wab, sig, 1.0)
        print(f"  σ = {sig:.1e}:  logL = {ll:.4f}  (Δ = {ll - ll0:+.4f})")
    print()

    summary = {
        "data": "DESI DR2 BAO summary statistics (7 bins, diagonal errors)",
        "fiducial": "flat LCDM, Om=0.315, H0=67.4, w0=-1, wa=0",
        "pure_CPL": {
            "w0": cpl["w0"],
            "wa": cpl["wa"],
            "chi2": cpl["chi2"],
            "logL": cpl["logL"],
            "AIC": aic(cpl["logL"], 2),
            "BIC": bic(cpl["logL"], 2, 7),
        },
        "lcdm_reference": {
            "w0": -1.0,
            "wa": 0.0,
            "chi2": float(chi2_cpl(-1.0, 0.0)),
            "logL": float(ll_lcdm),
        },
        "nested_effective_EoS": {
            "w0": eff["w0"],
            "wa": eff["wa"],
            "sigma": eff["sigma"],
            "theta": eff["theta"],
            "chi2": eff["chi2"],
            "logL": eff["logL"],
            "AIC": aic(eff["logL"], 4),
            "BIC": bic(eff["logL"], 4, 7),
        },
        "delta_logL_extension_minus_cpl": float(dlogL),
        "delta_AIC_extension_minus_cpl": float(daic),
        "notes": [
            "BAO-only fit with diagonal measurement errors.",
            "Nested stochastic amplitude is consistent with zero; pure CPL preferred by AIC.",
        ],
    }
    out_json = OUT_DIR / "eos_cpl_summary.json"
    out_txt = OUT_DIR / "eos_cpl_summary.txt"
    out_json.write_text(json.dumps(summary, indent=2) + "\n")
    lines = [
        "DESI DR2 BAO — CPL background and nested effective-EoS extension",
        "=" * 64,
        "Data: public DESI DR2 BAO summary statistics (7 bins; diagonal errors).",
        "Fiducial distances: flat ΛCDM (Ωm = 0.315, H0 = 67.4, w0 = −1, wa = 0).",
        "",
        "Pure CPL (free w0, wa)",
        f"  w0   = {cpl['w0']:+.4f}",
        f"  wa   = {cpl['wa']:+.4f}",
        f"  χ²   = {cpl['chi2']:.4f}",
        f"  logL = {cpl['logL']:.4f}",
        f"  AIC  = {aic(cpl['logL'], 2):.3f}",
        f"  BIC  = {bic(cpl['logL'], 2, 7):.3f}",
        "",
        "Reference ΛCDM (w0 = −1, wa = 0)",
        f"  χ²   = {chi2_cpl(-1.0, 0.0):.4f}",
        f"  logL = {ll_lcdm:.4f}",
        "",
        "Nested extension (CPL + σ, θ)",
        f"  w0    = {eff['w0']:+.4f}",
        f"  wa    = {eff['wa']:+.4f}",
        f"  σ     = {eff['sigma']:.3e}",
        f"  θ     = {eff['theta']:.4f}",
        f"  logL  = {eff['logL']:.4f}",
        f"  ΔlogL (extension − CPL) = {dlogL:+.4f}",
        f"  ΔAIC  (extension − CPL) = {daic:+.3f}",
        "",
        "Interpretation",
        "  BAO-only fits prefer a background close to flat ΛCDM.",
        "  The nested stochastic correction is not favoured (σ consistent with zero;",
        "  positive AIC penalty). Full DESI analyses use the survey covariance matrix",
        "  and multi-probe combinations; these numbers are a transparent BAO-only",
        "  baseline under the diagonal approximation used in this repository.",
        "",
    ]
    out_txt.write_text("\n".join(lines))
    print(f"Wrote {out_json}")
    print(f"Wrote {out_txt}")
    print("=" * 70)


if __name__ == "__main__":
    main()
