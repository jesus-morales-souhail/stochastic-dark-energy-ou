#!/usr/bin/env python3
"""
Profile likelihood for sigma_X on DESI DR2 BAO (full measurement cov when available).

For each fixed sigma_X, maximize logL over theta (OU kernel) with free
optional (w0, wa) in a second pass (background fixed to LCDM first).

This fills the gap noted in the paper: a formal 95% CL profile over sigma_X
rather than only the optimizer floor estimate.

Public-repo compatible. Author: Jesús Morales Souhail
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import minimize_scalar

ROOT = Path(__file__).resolve().parents[1]
import sys
sys.path.insert(0, str(ROOT / "scripts"))
from desi_dr2_data import load_alpha_dv, measurement_cov, add_ou_kernel, logL_gaussian
OUT = ROOT / "results" / "profile_sigma_x"
OUT.mkdir(parents=True, exist_ok=True)

_d = load_alpha_dv(prefer_file=True, use_full_cov=True)
z_eff = _d["z"]
alpha = _d["alpha"]
sigma_obs = _d["sigma"]
S_z = _d["S_z"]
C_MEAS = measurement_cov(_d)
print("profile_sigma_x data:", _d["source"], "cov_mode:", _d.get("cov_mode"))

residuals = alpha - 1.0
x = np.log(1.0 + z_eff)
n = len(z_eff)


def build_C(theta: float, sigma_X: float) -> np.ndarray:
    return add_ou_kernel(C_MEAS, S_z, x, theta, sigma_X)


def logL(theta: float, sigma_X: float) -> float:
    if theta <= 0 or sigma_X < 0:
        return -1e30
    return logL_gaussian(residuals, build_C(theta, sigma_X))


def max_logL_over_theta(sigma_X: float) -> tuple[float, float]:
    """Maximize logL over theta in (1e-3, 20) for fixed sigma_X."""
    def neg(th):
        return -logL(th, sigma_X)

    # multi-bracket
    best_ll, best_th = -1e30, 1.0
    for bracket in [(1e-3, 0.5), (0.1, 3.0), (1.0, 15.0)]:
        try:
            res = minimize_scalar(neg, bounds=bracket, method="bounded", options={"xatol": 1e-5})
            ll = -res.fun
            if ll > best_ll:
                best_ll, best_th = ll, float(res.x)
        except Exception:
            pass
    # also evaluate grid
    for th in np.geomspace(1e-3, 20, 40):
        ll = logL(th, sigma_X)
        if ll > best_ll:
            best_ll, best_th = ll, float(th)
    return best_ll, best_th


def main():
    # LCDM baseline (measurement cov only)
    ll_lcdm = logL_gaussian(residuals, C_MEAS)

    # profile over sigma_X
    sig_grid = np.geomspace(1e-6, 5e-2, 48)
    profile = []
    for s in sig_grid:
        ll, th = max_logL_over_theta(float(s))
        profile.append({"sigma_X": float(s), "theta_best": th, "logL": ll, "dlogL": ll - ll_lcdm})
        print(f"  sigma_X={s:.3e}  theta={th:.4f}  dlogL={ll-ll_lcdm:+.4f}")

    dlogL = np.array([p["dlogL"] for p in profile])
    # 95% CL for 1 parameter: ΔlogL = -1.92 (≈ half of 3.84 chi2)
    thr = -1.92
    # find largest sigma_X with dlogL >= thr (from left, null is at small s)
    # For upper limit: where profile drops below thr relative to max
    dlogL_rel = dlogL - np.max(dlogL)
    above = np.where(dlogL_rel >= thr)[0]
    if len(above):
        s_95 = float(sig_grid[above[-1]])
    else:
        s_95 = float(sig_grid[0])

    # also report where dlogL vs LCDM crosses -1.92 (if max is at LCDM)
    summary = {
        "data": "DESI DR2 BAO 7 bins, full measurement cov when available",
        "cov_mode": _d.get("cov_mode"),
        "ll_lcdm": ll_lcdm,
        "max_dlogL_vs_lcdm": float(np.max(dlogL)),
        "sigma_X_at_max": float(sig_grid[int(np.argmax(dlogL))]),
        "sigma_X_95CL_upper_profile": s_95,
        "criterion": "95% CL: Delta logL >= -1.92 from profile max (1 dof)",
        "paper_working_limit": 1.5e-4,
        "profile": profile,
        "notes": [
            "Uses official Gaussian BAO 13x13 projected to alpha (block-diagonal across bins).",
            "Background fixed to alpha=1 (LCDM fiducial). Free {w0,wa} is a separate scan.",
        ],
    }

    out_json = OUT / "profile_sigma_x.json"
    out_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    # figure
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.semilogx(sig_grid, dlogL_rel, "b-o", ms=3, lw=1.5)
    ax.axhline(0, color="k", lw=0.8)
    ax.axhline(thr, color="r", ls="--", label=r"$95\%$ CL ($\Delta\ln\mathcal{L}=-1.92$)")
    ax.axvline(1.5e-4, color="orange", ls=":", label=r"paper working limit $1.5\times10^{-4}$")
    ax.axvline(s_95, color="green", ls="-.", label=rf"profile $95\%$ $\sigma_X\leq{s_95:.2e}$")
    ax.set_xlabel(r"$\sigma_X$")
    ax.set_ylabel(r"$\Delta\ln\mathcal{L}$ (profile over $\theta$)")
    ax.set_title(r"DESI DR2 BAO: profile likelihood for $\sigma_X$ (OU, full meas. cov)")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "profile_sigma_x.png", dpi=150)
    fig.savefig(ROOT / "figures" / "profile_sigma_x.png", dpi=150)
    plt.close(fig)

    lines = [
        "Profile likelihood σ_X — DESI DR2 BAO (full meas. cov)",
        "=" * 50,
        f"ll_LCDM = {ll_lcdm:.4f}",
        f"max ΔlnL vs LCDM = {np.max(dlogL):+.4f} at σ_X = {sig_grid[int(np.argmax(dlogL))]:.3e}",
        f"95% CL upper (profile): σ_X ≤ {s_95:.3e}",
        f"Paper working limit:     σ_X < 1.5e-4",
        f"Wrote results/profile_sigma_x/profile_sigma_x.json",
        f"Wrote results/profile_sigma_x/profile_sigma_x.png",
    ]
    text = "\n".join(lines) + "\n"
    (OUT / "profile_sigma_x.txt").write_text(text, encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
