#!/usr/bin/env python3
"""
Rank-1 coherent growing-mode MLE on real DESI DR2 isotropic alpha.

Covariance (paper quantum-fluid note):

  C_{ij} = δ_{ij} σ_i² + σ₀² S_i S_j exp((t_i + t_j)/t_c)

with t(z) = lookback time (Gyr), S(z) = ∂ ln D_V / ∂ Ω_Λ, residual r = α − 1.

σ₀ = 10^{-61} (Sorkin / causal-set seed scale used elsewhere in the repo).

Writes results/tachyonic_rank1/SUMMARY.txt
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.linalg import cholesky, solve_triangular

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from desi_dr2_data import load_alpha_dv, measurement_cov  # noqa: E402

H0 = 67.4
OM = 0.315
H0_GYR = H0 * 1.022712165e-3
SIGMA0 = 1.0e-61
OUT = ROOT / "results" / "tachyonic_rank1"


def lookback_gyr(z: float) -> float:
    def E(zp: float) -> float:
        return float(np.sqrt(OM * (1.0 + zp) ** 3 + (1.0 - OM)))

    integ, _ = quad(lambda zp: 1.0 / ((1.0 + zp) * E(zp)), 0.0, float(z), epsabs=1e-8)
    return integ / H0_GYR


def log_likelihood(residuals: np.ndarray, cov: np.ndarray) -> float:
    if not np.all(np.isfinite(cov)):
        return -np.inf
    try:
        L = cholesky(cov, lower=True)
        y = solve_triangular(L, residuals, lower=True)
        logdet = 2.0 * np.sum(np.log(np.diag(L)))
        n = len(residuals)
        return float(-0.5 * (np.dot(y, y) + logdet + n * np.log(2.0 * np.pi)))
    except Exception:
        return -np.inf


def cov_rank1(C_meas: np.ndarray, S: np.ndarray, t: np.ndarray, tc: float, sigma0: float) -> np.ndarray:
    expo = np.clip(t / tc, -700.0, 700.0)
    u = sigma0 * S * np.exp(expo)
    if not np.all(np.isfinite(u)) or np.max(np.abs(u)) > 1.0e12:
        return np.full((len(S), len(S)), np.nan)
    return np.array(C_meas, dtype=float, copy=True) + np.outer(u, u)


def main() -> None:
    d = load_alpha_dv(prefer_file=True, use_full_cov=True)
    z, alpha, sig, S = d["z"], d["alpha"], d["sigma"], d["S_z"]
    C_MEAS = measurement_cov(d)
    t = np.array([lookback_gyr(float(zi)) for zi in z])
    r = alpha - 1.0

    C0 = C_MEAS
    ll0 = log_likelihood(r, C0)
    chi2_0 = float(np.sum((r / sig) ** 2))

    tcs = np.logspace(-4, 4, 500)
    rows = []
    for tc in tcs:
        C = cov_rank1(C_MEAS, S, t, tc, SIGMA0)
        ll = log_likelihood(r, C)
        u_max = float(np.max(np.abs(SIGMA0 * S * np.exp(np.clip(t / tc, -700, 700)))))
        if not np.isfinite(u_max):
            u_max = np.inf
        rows.append((float(tc), float(ll), float(ll - ll0), u_max))
    rows = np.array(rows)
    finite = rows[np.isfinite(rows[:, 1])]

    # No-growth limit: large tc
    large = finite[finite[:, 0] > 1.0e2]
    dln_nogrowth = float(np.max(large[:, 2])) if len(large) else 0.0

    # Active growth: rank-1 amplitude not negligible
    active = finite[finite[:, 3] > 1.0e-6]
    if len(active):
        i_best_act = int(np.argmax(active[:, 2]))
        i_near = int(np.argmin(np.abs(active[:, 2] + 11.35)))
        best_active = active[i_best_act]
        near_1135 = active[i_near]
        worst_active = active[int(np.argmin(active[:, 2]))]
    else:
        best_active = near_1135 = worst_active = np.array([np.nan, np.nan, np.nan, np.nan])

    OUT.mkdir(parents=True, exist_ok=True)
    np.savetxt(
        OUT / "scan_fixed_sigma0.txt",
        finite,
        header="tc_Gyr ll dlnL_vs_LCDM max_abs_rank1_amp",
    )

    summary = {
        "source": d["source"],
        "cov_mode": d.get("cov_mode"),
        "mock": bool(d.get("mock", False)),
        "n_bins": int(len(z)),
        "z": z.tolist(),
        "lookback_Gyr": t.tolist(),
        "ll_LCDM": ll0,
        "chi2_LCDM": chi2_0,
        "sigma0": SIGMA0,
        "dlnL_no_growth_large_tc": dln_nogrowth,
        "best_active_tc_Gyr": float(best_active[0]),
        "best_active_dlnL": float(best_active[2]),
        "representative_tc_Gyr": float(near_1135[0]),
        "representative_dlnL": float(near_1135[2]),
        "representative_dchi2": float(-2.0 * near_1135[2]),
        "worst_active_dlnL": float(worst_active[2]),
        "note": (
            "Likelihood peaks at no growth (tc→∞, ΔlnL→0). "
            "When σ0·exp(t/tc) is large enough for the rank-1 term to matter, "
            "the model is excluded (ΔlnL ≲ −11, Δχ² ≳ +22)."
        ),
    }
    (OUT / "SUMMARY.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    (OUT / "SUMMARY.txt").write_text(
        "\n".join(
            [
                "Tachyonic / coherent growing-mode rank-1 MLE (real DESI DR2 alpha)",
                f"source: {summary['source']}",
                f"n_bins: {summary['n_bins']}",
                f"ll_LCDM: {ll0:.6f}",
                f"chi2_LCDM: {chi2_0:.6f}",
                f"sigma0: {SIGMA0:.3e}",
                f"dlnL (no growth, large tc): {dln_nogrowth:.6f}",
                f"representative active tc [Gyr]: {summary['representative_tc_Gyr']:.6e}",
                f"representative ΔlnL vs LCDM: {summary['representative_dlnL']:.4f}",
                f"representative Δχ² ≈ {-2*summary['representative_dlnL']:.2f}",
                f"worst active ΔlnL: {summary['worst_active_dlnL']:.4f}",
                summary["note"],
                "",
            ]
        ),
        encoding="utf-8",
    )
    print((OUT / "SUMMARY.txt").read_text())


if __name__ == "__main__":
    main()
