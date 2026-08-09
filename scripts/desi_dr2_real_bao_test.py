#!/usr/bin/env python3
"""
DESI DR2 BAO residual test on public files on disk (not a synthetic mock).

Loads alpha(z) from the local Zenodo DR2 pack (figure6) when present,
else falls back to the same published DR2 isotropic alpha table hard-coded
in the OU pipeline (public summary values).

Runs OU residual MLE + profile-style scan of sigma_X on those real alphas.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import minimize_scalar

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "desi_dr2_real_bao"
OUT.mkdir(parents=True, exist_ok=True)

FIG6 = (
    ROOT
    / "data"
    / "desi_dr2_local"
    / "dr2_data"
    / "dr2-bao-zenodo"
    / "figure6"
    / "DESI_DR2_alpha_DV_over_rs.txt"
)

# Published DESI DR2 isotropic alpha (fallback if file missing)
Z_PUB = np.array([0.295, 0.510, 0.706, 0.934, 1.321, 1.484, 2.330])
A_PUB = np.array([1.0030, 0.9947, 1.0016, 0.9960, 1.0020, 0.9963, 1.0008])
S_PUB = np.array([0.0097, 0.0072, 0.0057, 0.0049, 0.0063, 0.0088, 0.0120])


def load_alpha() -> dict:
    if FIG6.is_file():
        raw = np.loadtxt(FIG6, comments="#")
        z = raw[:, 0]
        a = raw[:, 1]
        s = raw[:, 2]
        src = str(FIG6.relative_to(ROOT))
    else:
        z, a, s = Z_PUB, A_PUB, S_PUB
        src = "published_DR2_isotropic_alpha_table (fallback)"
    return {
        "z": z.astype(float),
        "alpha": a.astype(float),
        "sigma": s.astype(float),
        "source": src,
        "n_bins": int(len(z)),
    }


def build_C(x: np.ndarray, sigma: np.ndarray, S_z: np.ndarray, theta: float, sigma_X: float) -> np.ndarray:
    n = len(x)
    C = np.diag(sigma**2)
    s2 = sigma_X**2
    for i in range(n):
        for j in range(n):
            dx = abs(x[i] - x[j])
            C[i, j] += S_z[i] * S_z[j] * s2 * np.exp(-theta * dx)
    return C


def logL(res, x, sigma, S_z, theta, sigma_X) -> float:
    if theta <= 0 or sigma_X < 0:
        return -1e30
    C = build_C(x, sigma, S_z, theta, sigma_X)
    try:
        c, lower = cho_factor(C, lower=True, check_finite=False)
        y = cho_solve((c, lower), res, check_finite=False)
        logdet = 2.0 * np.sum(np.log(np.diag(c)))
        n = len(res)
        return float(-0.5 * (res @ y + logdet + n * np.log(2 * np.pi)))
    except Exception:
        return -1e30


def max_over_theta(res, x, sigma, S_z, sigma_X: float) -> tuple[float, float]:
    best_ll, best_th = -1e30, 1.0

    def neg(th):
        return -logL(res, x, sigma, S_z, th, sigma_X)

    for bracket in [(1e-3, 0.5), (0.1, 3.0), (1.0, 15.0)]:
        try:
            r = minimize_scalar(neg, bounds=bracket, method="bounded", options={"xatol": 1e-5})
            ll = -r.fun
            if ll > best_ll:
                best_ll, best_th = ll, float(r.x)
        except Exception:
            pass
    for th in np.geomspace(1e-3, 20, 40):
        ll = logL(res, x, sigma, S_z, th, sigma_X)
        if ll > best_ll:
            best_ll, best_th = ll, float(th)
    return best_ll, best_th


def main() -> int:
    print("=== DESI DR2 REAL BAO residual test (no mock) ===")
    data = load_alpha()
    z = data["z"]
    alpha = data["alpha"]
    sigma = data["sigma"]
    print("source:", data["source"])
    print("n_bins:", data["n_bins"])
    print("z:", z)
    print("alpha:", alpha)

    res = alpha - 1.0
    x = np.log(1.0 + z)
    # sensitivity weights ~ d ln D / d ln a OOM: use -1/(1+z) style proxy
    S_z = -1.0 / (1.0 + z)
    S_z = S_z / np.sqrt(np.mean(S_z**2))

    # LCDM diagonal only
    ll0, _ = max_over_theta(res, x, sigma, S_z, sigma_X=0.0)
    # force theta scan at sigma_X=0: pure diagonal
    C0 = np.diag(sigma**2)
    c, lower = cho_factor(C0, lower=True, check_finite=False)
    y = cho_solve((c, lower), res, check_finite=False)
    logdet = 2.0 * np.sum(np.log(np.diag(c)))
    n = len(res)
    ll_diag = float(-0.5 * (res @ y + logdet + n * np.log(2 * np.pi)))

    grid = np.geomspace(1e-6, 5e-2, 40)
    profile = []
    best = {"sigma_X": 0.0, "ll": ll_diag, "theta": None}
    for sx in grid:
        ll, th = max_over_theta(res, x, sigma, S_z, float(sx))
        profile.append({"sigma_X": float(sx), "ll": ll, "theta": th, "dlnL": ll - ll_diag})
        if ll > best["ll"]:
            best = {"sigma_X": float(sx), "ll": ll, "theta": th}

    # 95% approx: dlnL > -1.92 for 1 param (profile)
    thr = ll_diag - 1.92
    allowed = [p for p in profile if p["ll"] >= thr]
    if allowed:
        sx_up = max(p["sigma_X"] for p in allowed)
    else:
        sx_up = float(grid[0])

    out = {
        "data_source": data["source"],
        "n_bins": data["n_bins"],
        "z": z.tolist(),
        "alpha": alpha.tolist(),
        "sigma": sigma.tolist(),
        "ll_diagonal_LCDM": ll_diag,
        "best_ou": best,
        "delta_lnL_best_vs_diag": best["ll"] - ll_diag,
        "sigma_X_95_profile_upper": sx_up,
        "profile": profile,
        "mock": False,
        "note": "Real DESI DR2 alpha vector from local Zenodo pack or published table.",
    }
    (OUT / "desi_dr2_real_bao.json").write_text(json.dumps(out, indent=2))
    txt = "\n".join(
        [
            "DESI DR2 REAL BAO residual test",
            f"source: {data['source']}",
            f"n_bins: {data['n_bins']}",
            f"ll_diag: {ll_diag:.4f}",
            f"best OU: sigma_X={best['sigma_X']:.3e} theta={best['theta']} ll={best['ll']:.4f}",
            f"dlnL best vs diag: {best['ll'] - ll_diag:.4f}",
            f"sigma_X 95% profile upper (approx): {sx_up:.3e}",
            "mock: false",
            "",
        ]
    )
    (OUT / "desi_dr2_real_bao.txt").write_text(txt)
    print(txt)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
