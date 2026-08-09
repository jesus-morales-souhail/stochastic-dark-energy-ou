#!/usr/bin/env python3
"""
DESI DR2 BAO residual test on public files on disk (not a synthetic mock).

Uses isotropic alpha from the Zenodo DR2 pack (figure6) and the official
Gaussian BAO 13×13 covariance projected to alpha (full measurement cov).

Runs OU residual MLE + profile-style scan of sigma_X.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar

ROOT = Path(__file__).resolve().parents[1]
import sys

sys.path.insert(0, str(ROOT / "scripts"))
from desi_dr2_data import (  # noqa: E402
    add_ou_kernel,
    load_alpha_dv,
    logL_gaussian,
    measurement_cov,
)

OUT = ROOT / "results" / "desi_dr2_real_bao"
OUT.mkdir(parents=True, exist_ok=True)


def max_over_theta(res, x, S_z, C_meas, sigma_X: float) -> tuple[float, float]:
    best_ll, best_th = -1e30, 1.0

    def neg(th):
        if th <= 0:
            return 1e30
        C = add_ou_kernel(C_meas, S_z, x, float(th), sigma_X)
        return -logL_gaussian(res, C)

    for bracket in [(1e-3, 0.5), (0.1, 3.0), (1.0, 15.0)]:
        try:
            r = minimize_scalar(neg, bounds=bracket, method="bounded", options={"xatol": 1e-5})
            ll = -r.fun
            if ll > best_ll:
                best_ll, best_th = ll, float(r.x)
        except Exception:
            pass
    for th in np.geomspace(1e-3, 20, 40):
        ll = -neg(th)
        if ll > best_ll:
            best_ll, best_th = ll, float(th)
    return best_ll, best_th


def main() -> int:
    print("=== DESI DR2 REAL BAO residual test (full meas. cov) ===")
    data = load_alpha_dv(prefer_file=True, use_full_cov=True)
    data["n_bins"] = int(len(data["z"]))
    z = data["z"]
    alpha = data["alpha"]
    sigma = data["sigma"]
    S_z = data["S_z"]
    C_meas = measurement_cov(data)
    print("source:", data["source"])
    print("cov_mode:", data.get("cov_mode"))
    print("n_bins:", data["n_bins"])

    res = alpha - 1.0
    x = np.log(1.0 + z)

    ll_meas = logL_gaussian(res, C_meas)

    grid = np.geomspace(1e-6, 5e-2, 40)
    profile = []
    best = {"sigma_X": 0.0, "ll": ll_meas, "theta": None}
    for sx in grid:
        ll, th = max_over_theta(res, x, S_z, C_meas, float(sx))
        profile.append({"sigma_X": float(sx), "ll": ll, "theta": th, "dlnL": ll - ll_meas})
        if ll > best["ll"]:
            best = {"sigma_X": float(sx), "ll": ll, "theta": th}

    thr = ll_meas - 1.92
    allowed = [p for p in profile if p["ll"] >= thr]
    sx_up = max(p["sigma_X"] for p in allowed) if allowed else float(grid[0])

    out = {
        "data_source": data["source"],
        "cov_mode": data.get("cov_mode"),
        "n_bins": data["n_bins"],
        "z": z.tolist(),
        "alpha": alpha.tolist(),
        "sigma": sigma.tolist(),
        "ll_measurement_cov_LCDM": ll_meas,
        "best_ou": best,
        "delta_lnL_best_vs_meas": best["ll"] - ll_meas,
        "sigma_X_95_profile_upper": sx_up,
        "profile": profile,
        "mock": False,
        "note": (
            "Real DESI DR2 alpha; measurement cov = official Gaussian BAO 13x13 "
            "projected to isotropic alpha (block-diagonal across bins)."
        ),
    }
    (OUT / "desi_dr2_real_bao.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    txt = "\n".join(
        [
            "DESI DR2 REAL BAO residual test",
            f"source: {data['source']}",
            f"cov_mode: {data.get('cov_mode')}",
            f"n_bins: {data['n_bins']}",
            f"ll_meas: {ll_meas:.4f}",
            f"best OU: sigma_X={best['sigma_X']:.3e} theta={best['theta']} ll={best['ll']:.4f}",
            f"dlnL best vs meas: {best['ll'] - ll_meas:.4f}",
            f"sigma_X 95% profile upper (approx): {sx_up:.3e}",
            "mock: false",
            "",
        ]
    )
    (OUT / "desi_dr2_real_bao.txt").write_text(txt, encoding="utf-8")
    print(txt)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
