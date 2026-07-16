#!/usr/bin/env python3
"""
Joint MLE on DESI DR2 BAO: {w0, wa, theta, sigma_X} vs nested baselines.

Uses the same alpha construction as eos_efectiva (CPL distances relative to
flat LCDM fiducial) plus the OU residual kernel of ou_bao_likelihood.

Author: Jesús Morales Souhail
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "joint_w0wa_sigma"
OUT.mkdir(parents=True, exist_ok=True)
FIG = ROOT / "figures"

# DESI DR2 BAO (repo standard)
z_eff = np.array([0.295, 0.510, 0.706, 0.934, 1.321, 1.484, 2.330])
alpha_obs = np.array([1.0030, 0.9947, 1.0016, 0.9960, 1.0020, 0.9963, 1.0008])
sigma_obs = np.array([0.0097, 0.0072, 0.0057, 0.0049, 0.0063, 0.0088, 0.0120])
S_z = np.array([-0.284, -0.462, -0.595, -0.719, -0.870, -0.917, -1.070])

H0, Om, c_kms = 67.4, 0.315, 299792.458
n = len(z_eff)
xbins = np.log(1.0 + z_eff)


def E2_cpl(z, w0, wa):
    rho = (1 + z) ** (3 * (1 + w0 + wa)) * np.exp(-3 * wa * z / (1 + z))
    return Om * (1 + z) ** 3 + (1 - Om) * rho


def H_cpl(z, w0, wa):
    e2 = E2_cpl(z, w0, wa)
    if e2 <= 0 or not np.isfinite(e2):
        return np.nan
    return H0 * np.sqrt(e2)


def DM(z, w0, wa):
    def integrand(zp):
        h = H_cpl(zp, w0, wa)
        if not np.isfinite(h) or h <= 0:
            raise ValueError("bad H")
        return 1.0 / h

    return c_kms * quad(integrand, 0, z, limit=100, epsabs=1e-5)[0]


def DH(z, w0, wa):
    h = H_cpl(z, w0, wa)
    if not np.isfinite(h) or h <= 0:
        raise ValueError("bad H")
    return c_kms / h


def DV(z, w0, wa):
    return (DM(z, w0, wa) ** 2 * DH(z, w0, wa) * z) ** (1.0 / 3.0)


# Fiducial LCDM distances (cache)
_DV_fid = None
_DM_fid = None
_DH_fid = None


def _ensure_fid():
    global _DV_fid, _DM_fid, _DH_fid
    if _DV_fid is None:
        _DV_fid = DV(float(z_eff[0]), -1.0, 0.0)
        _DM_fid = np.array([DM(float(z), -1.0, 0.0) for z in z_eff[1:]])
        _DH_fid = np.array([DH(float(z), -1.0, 0.0) for z in z_eff[1:]])


def alpha_pred_cpl(w0, wa):
    _ensure_fid()
    out = np.zeros(n)
    out[0] = DV(float(z_eff[0]), w0, wa) / _DV_fid
    for i in range(1, n):
        z = float(z_eff[i])
        out[i] = (DM(z, w0, wa) / _DM_fid[i - 1]) ** (2 / 3) * (
            DH(z, w0, wa) / _DH_fid[i - 1]
        ) ** (1 / 3)
    return out


def logL_total(w0, wa, theta, sigma_X):
    """Gaussian logL with OU kernel on residuals alpha_obs - alpha_pred(w0,wa)."""
    if theta < 1e-4 or sigma_X < 0:
        return -1e30
    if not (-2.5 < w0 < 0.2 and -3.0 < wa < 2.0):
        return -1e30
    try:
        pred = alpha_pred_cpl(w0, wa)
    except Exception:
        return -1e30
    if not np.all(np.isfinite(pred)):
        return -1e30
    r = alpha_obs - pred
    C = np.diag(sigma_obs**2)
    s2 = sigma_X**2
    for i in range(n):
        for j in range(n):
            dx = abs(xbins[i] - xbins[j])
            C[i, j] += S_z[i] * S_z[j] * s2 * np.exp(-theta * dx)
    try:
        c, lower = cho_factor(C, lower=True, check_finite=False)
        y = cho_solve((c, lower), r, check_finite=False)
        logdet = 2.0 * np.sum(np.log(np.diag(c)))
        return float(-0.5 * (r @ y + logdet + n * np.log(2 * np.pi)))
    except Exception:
        return -1e30


def fit(model: str):
    """
    model:
      lcdm: w0=-1, wa=0, sigma_X=0
      cpl: free w0,wa; sigma_X=0
      ou_fixed_bg: w0=-1,wa=0; free theta, sigma_X
      joint: free w0, wa, theta, sigma_X
    """
    starts = []
    if model == "lcdm":
        ll = logL_total(-1.0, 0.0, 1.0, 0.0)
        return {
            "model": model,
            "w0": -1.0,
            "wa": 0.0,
            "theta": None,
            "sigma_X": 0.0,
            "logL": ll,
            "k": 0,
        }
    if model == "cpl":
        starts = [
            [-1.0, 0.0],
            [-0.9, -0.3],
            [-0.87, -0.41],
            [-0.99, -0.02],
            [-0.7, -0.5],
        ]

        def neg(p):
            return -logL_total(p[0], p[1], 1.0, 0.0)

        best = None
        for x0 in starts:
            res = minimize(neg, x0, method="Nelder-Mead", options={"maxiter": 4000, "xatol": 1e-6, "fatol": 1e-6})
            if best is None or res.fun < best.fun:
                best = res
        w0, wa = best.x
        ll = -best.fun
        return {"model": model, "w0": float(w0), "wa": float(wa), "theta": None, "sigma_X": 0.0, "logL": ll, "k": 2}

    if model == "ou_fixed_bg":
        starts = [
            [1.0, 1e-5],
            [0.5, 1e-4],
            [2.0, 1e-3],
            [0.1, 1e-3],
            [5.0, 1e-2],
        ]

        def neg(p):
            th, s = p
            if th < 1e-3 or s < 0 or s > 0.1:
                return 1e30
            return -logL_total(-1.0, 0.0, th, s)

        best = None
        for x0 in starts:
            res = minimize(neg, x0, method="Nelder-Mead", options={"maxiter": 5000, "xatol": 1e-6, "fatol": 1e-6})
            if best is None or res.fun < best.fun:
                best = res
        th, s = best.x
        ll = -best.fun
        return {
            "model": model,
            "w0": -1.0,
            "wa": 0.0,
            "theta": float(th),
            "sigma_X": float(max(s, 0.0)),
            "logL": ll,
            "k": 2,
        }

    # joint
    starts = [
        [-1.0, 0.0, 1.0, 1e-5],
        [-0.99, -0.02, 0.5, 1e-4],
        [-0.87, -0.41, 1.0, 1e-3],
        [-0.9, -0.2, 2.0, 1e-3],
        [-1.1, 0.1, 0.2, 5e-3],
        [-0.95, -0.1, 5.0, 1e-2],
    ]

    def neg(p):
        w0, wa, th, s = p
        if th < 1e-3 or s < 0 or s > 0.1:
            return 1e30
        return -logL_total(w0, wa, th, s)

    best = None
    for x0 in starts:
        res = minimize(neg, x0, method="Nelder-Mead", options={"maxiter": 8000, "xatol": 1e-6, "fatol": 1e-6})
        if best is None or res.fun < best.fun:
            best = res
    w0, wa, th, s = best.x
    ll = -best.fun
    return {
        "model": model,
        "w0": float(w0),
        "wa": float(wa),
        "theta": float(th),
        "sigma_X": float(max(s, 0.0)),
        "logL": ll,
        "k": 4,
    }


def aic(ll, k):
    return 2 * k - 2 * ll


def bic(ll, k, npts=n):
    return k * np.log(npts) - 2 * ll


def main():
    print("Joint DESI DR2 BAO fits...")
    results = {}
    for m in ["lcdm", "cpl", "ou_fixed_bg", "joint"]:
        print(f"  fitting {m}...")
        results[m] = fit(m)
        r = results[m]
        print(
            f"    w0={r['w0']:+.4f} wa={r.get('wa') if r['wa'] is not None else 0:+.4f} "
            f"th={r['theta']} sX={r['sigma_X']:.3e} logL={r['logL']:.4f}"
        )

    ll0 = results["lcdm"]["logL"]
    table = []
    for m, r in results.items():
        k = r["k"]
        table.append(
            {
                **r,
                "dlogL_vs_lcdm": r["logL"] - ll0,
                "AIC": aic(r["logL"], k) if k > 0 else None,
                "BIC": bic(r["logL"], k) if k > 0 else None,
                "dAIC_vs_cpl": (aic(r["logL"], k) - aic(results["cpl"]["logL"], 2))
                if k > 0
                else None,
            }
        )

    summary = {
        "data": "DESI DR2 BAO 7 bins, diagonal cov, S(z) fixed",
        "models": table,
        "reading": [
            "If joint drives sigma_X→0 and (w0,wa) near LCDM, null stochastic is robust.",
            "If joint prefers large sigma_X only with exotic (w0,wa), degeneracy warning.",
            "Compare AIC: joint vs cpl (k=4 vs k=2).",
        ],
    }
    (OUT / "joint_w0wa_sigma.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    lines = ["Joint MLE DESI DR2 BAO {w0,wa,theta,sigma_X}", "=" * 55]
    for r in table:
        lines.append(
            f"{r['model']:12s}  w0={r['w0']:+.4f}  wa={r['wa'] if r['wa'] is not None else 0:+.4f}  "
            f"θ={str(r['theta'])[:8]:>8}  σ_X={r['sigma_X']:.3e}  "
            f"logL={r['logL']:.4f}  ΔlnL={r['dlogL_vs_lcdm']:+.4f}  "
            f"AIC={r['AIC']}"
        )
    j, c = results["joint"], results["cpl"]
    lines.append("")
    lines.append(f"ΔAIC(joint − cpl) = {aic(j['logL'],4) - aic(c['logL'],2):+.3f}")
    lines.append(f"σ_X joint = {j['sigma_X']:.3e}  (→0 means stochastic not needed)")
    text = "\n".join(lines) + "\n"
    (OUT / "joint_w0wa_sigma.txt").write_text(text, encoding="utf-8")
    print(text)

    # bar figure of dlogL
    import matplotlib.pyplot as plt

    names = [r["model"] for r in table]
    dls = [r["dlogL_vs_lcdm"] for r in table]
    fig, ax = plt.subplots(figsize=(7, 4))
    colors = ["#455a64", "#1565c0", "#6a1b9a", "#c62828"]
    ax.bar(names, dls, color=colors[: len(names)])
    ax.axhline(0, color="k", lw=0.8)
    ax.set_ylabel(r"$\Delta\ln\mathcal{L}$ vs $\Lambda$CDM")
    ax.set_title("DESI DR2 BAO joint model comparison (diagonal cov)")
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "joint_w0wa_sigma.png", dpi=150)
    fig.savefig(FIG / "joint_w0wa_sigma.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
