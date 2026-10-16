#!/usr/bin/env python3
"""
DESI DR2 BAO profile likelihoods for Everpresent Lambda kernels.

Copy of scripts/profile_sigma_x_desi.py (original untouched), same data, same residuals
(alpha_iso - 1), same full measurement covariance (cov_mode full_gaussian_bao_projected_13x13),
same 95 % criterion (Delta lnL = -1.92 from the profile maximum, 1 dof). Uses the local
copy everpresent/desi_dr2_data.py (writes nothing outside everpresent/).

Runs:
  T0  OU, theta profiled in (1e-3, 20), sigma grid geomspace(1e-6, 5e-2, 48): must give the
      current 2.506e-2 (test that the copy is not broken).
  T1  OU, theta = 1 fixed (information).
  T2  OU, theta fixed to theta_fit of the Ahmed process (resultados/ahmed_sim.json).
  Z1  Zwane eq. (2.9), mu = 0.9804 fixed, profile in alpha.
  Z2  Zwane eq. (2.9), mu profiled in [0.1, 10], profile in alpha.
Limits are given on the 48-point grid (as the current profile) and refined by root finding
between grid neighbours.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, minimize_scalar

from desi_dr2_data import add_ou_kernel, add_zwane_kernel, load_alpha_dv, logL_gaussian, measurement_cov

HERE = Path(__file__).resolve().parent
OUT = HERE / "resultados"
THR = -1.92
GRID = np.geomspace(1e-6, 5e-2, 48)  # same as scripts/profile_sigma_x_desi.py
GRID_ZW = np.geomspace(1e-6, 2.0, 80)
ALPHA_ZW, MU_ZW = 0.8824, 0.9804  # Zwane et al. Table 1

_d = load_alpha_dv(prefer_file=True, use_full_cov=True)
Z, S_Z = _d["z"], _d["S_z"]
C_MEAS = measurement_cov(_d)
RES = _d["alpha"] - 1.0
X = np.log(1.0 + Z)
LL0 = logL_gaussian(RES, C_MEAS)


def ll_ou(theta: float, s: float) -> float:
    if theta <= 0 or s < 0:
        return -1e30
    return logL_gaussian(RES, add_ou_kernel(C_MEAS, S_Z, X, theta, s))


def ll_zw(mu: float, a: float) -> float:
    return logL_gaussian(RES, add_zwane_kernel(C_MEAS, S_Z, Z, a, mu))


def prof_ou_theta(s: float) -> tuple[float, float]:
    """Same maximisation over theta as scripts/profile_sigma_x_desi.py."""
    best_ll, best_th = -1e30, 1.0
    for br in [(1e-3, 0.5), (0.1, 3.0), (1.0, 15.0)]:
        r = minimize_scalar(lambda th: -ll_ou(th, s), bounds=br, method="bounded", options={"xatol": 1e-5})
        if -r.fun > best_ll:
            best_ll, best_th = -r.fun, float(r.x)
    for th in np.geomspace(1e-3, 20, 40):
        ll = ll_ou(th, s)
        if ll > best_ll:
            best_ll, best_th = ll, float(th)
    return best_ll, best_th


def tanh_kernel(a: float, mu: float, n: int = 60) -> np.ndarray:
    """Cov[tanh(O_i), tanh(O_j)] for the Gaussian O of eq. (2.9) (Gauss-Hermite, 2D)."""
    t, w = np.polynomial.hermite_e.hermegauss(n)
    w = w / w.sum()
    lam = -np.log1p(Z)
    rho = np.exp(-mu * (lam[:, None] - lam[None, :]) ** 2 / 2)
    K = np.empty_like(rho)
    for i in range(len(Z)):
        for j in range(len(Z)):
            r = min(rho[i, j], 1 - 1e-12)
            x1 = a * t[:, None]
            x2 = a * (r * t[:, None] + np.sqrt(1 - r * r) * t[None, :])
            K[i, j] = np.sum(w[:, None] * w[None, :] * np.tanh(x1) * np.tanh(x2))
    return K


def prof_zw_mu(a: float) -> tuple[float, float]:
    """Maximise over ln mu in [ln 0.1, ln 10]."""
    lg = np.linspace(np.log(0.1), np.log(10.0), 41)
    lls = [ll_zw(np.exp(v), a) for v in lg]
    k = int(np.argmax(lls))
    lo, hi = lg[max(k - 1, 0)], lg[min(k + 1, len(lg) - 1)]
    r = minimize_scalar(lambda v: -ll_zw(np.exp(v), a), bounds=(lo, hi), method="bounded", options={"xatol": 1e-6})
    if -r.fun >= lls[k]:
        return -r.fun, float(np.exp(r.x))
    return lls[k], float(np.exp(lg[k]))


def profile(f, grid) -> dict:
    """f(amp) -> (lnL, nuisance). Returns grid profile, grid and refined 95 % upper limits."""
    rows = [(float(g),) + tuple(f(float(g))) for g in grid]
    ll = np.array([r[1] for r in rows])
    imax = int(np.argmax(ll))
    llmax = ll[imax]
    rel = ll - llmax
    above = np.where(rel >= THR)[0]
    up_grid = float(grid[above[-1]])
    up_ref = None
    j = above[-1]
    if j + 1 < len(grid):
        g = lambda s: f(s)[0] - llmax - THR
        up_ref = float(brentq(g, grid[j], grid[j + 1], xtol=1e-10 * grid[j]))
    return {
        "best_amp": float(grid[imax]),
        "best_nuisance": rows[imax][2],
        "max_dlnL_vs_null": float(llmax - LL0),
        "upper95_grid": up_grid,
        "upper95_refined": up_ref,
        "nuisance_at_upper95": float(f(up_ref)[1]) if up_ref else None,
        "profile": [{"amp": r[0], "lnL": r[1], "nuisance": r[2], "dlnL_rel": float(r[1] - llmax)} for r in rows],
    }


def main() -> None:
    OUT.mkdir(exist_ok=True)
    print("data:", _d["source"], "| cov_mode:", _d.get("cov_mode"), f"| lnL(null) = {LL0:.6f}")
    theta_fit = json.loads((OUT / "ahmed_sim.json").read_text())["process"]["R_alpha0.02"]["theta_fit"]

    res = {"cov_mode": _d.get("cov_mode"), "lnL_null": LL0, "threshold": THR}
    res["T0_ou_theta_profiled"] = profile(prof_ou_theta, GRID)
    res["T1_ou_theta1"] = profile(lambda s: (ll_ou(1.0, s), 1.0), GRID)
    res["T2_ou_theta_fit"] = profile(lambda s: (ll_ou(theta_fit, s), theta_fit), GRID)
    res["T2_ou_theta_fit"]["theta"] = theta_fit
    res["Z1_zwane_mu_fixed"] = profile(lambda a: (ll_zw(MU_ZW, a), MU_ZW), GRID_ZW)
    res["Z2_zwane_mu_free"] = profile(prof_zw_mu, GRID_ZW)

    ll_a_fix = ll_zw(MU_ZW, ALPHA_ZW)
    ll_a_free, mu_a_free = prof_zw_mu(ALPHA_ZW)
    for key, ll, mu in (("Z1_zwane_mu_fixed", ll_a_fix, MU_ZW), ("Z2_zwane_mu_free", ll_a_free, mu_a_free)):
        d = ll - LL0
        res[key]["alpha_zwane"] = {"alpha": ALPHA_ZW, "mu": mu, "dlnL_vs_alpha0": float(d), "dchi2": float(-2 * d),
                                   "nsigma_sqrt_dchi2": float(np.sqrt(max(-2 * d, 0.0)))}
    Kt = tanh_kernel(ALPHA_ZW, MU_ZW)
    d_t = logL_gaussian(RES, C_MEAS + np.outer(S_Z, S_Z) * Kt) - LL0
    res["Z1_zwane_mu_fixed"]["alpha_zwane_with_tanh"] = {
        "note": "Gaussian likelihood with the exact covariance of tanh(O_hat), eq. (2.10)",
        "std_tanh": float(np.sqrt(Kt[0, 0])), "dlnL_vs_alpha0": float(d_t), "dchi2": float(-2 * d_t),
        "nsigma_sqrt_dchi2": float(np.sqrt(max(-2 * d_t, 0.0)))}
    # where mu sits along the free-mu profile
    res["Z2_zwane_mu_free"]["mu_along_profile"] = [
        {"alpha": p["amp"], "mu": p["nuisance"]} for p in res["Z2_zwane_mu_free"]["profile"] if p["amp"] >= 1e-3
    ]

    (OUT / "profile_desi_everpresent.json").write_text(json.dumps(res, indent=2), encoding="utf-8")

    vig = json.loads((HERE.parents[1] / "results" / "profile_sigma_x" / "profile_sigma_x.json").read_text())
    print(f"T0 OU theta perfilado: 95% grid {res['T0_ou_theta_profiled']['upper95_grid']:.6e}"
          f"  (vigente {vig['sigma_X_95CL_upper_profile']:.6e})  refinado {res['T0_ou_theta_profiled']['upper95_refined']:.4e}")
    for k in ("T1_ou_theta1", "T2_ou_theta_fit"):
        r = res[k]
        print(f"{k}: max dlnL {r['max_dlnL_vs_null']:+.2e} en {r['best_amp']:.2e}; 95% grid {r['upper95_grid']:.4e}"
              f"  refinado {r['upper95_refined']:.4e}" + (f"  (theta = {theta_fit:.4f})" if k.startswith("T2") else ""))
    for k in ("Z1_zwane_mu_fixed", "Z2_zwane_mu_free"):
        r = res[k]
        az = r["alpha_zwane"]
        print(f"{k}: mejor alpha {r['best_amp']:.2e} (mu {r['best_nuisance']:.3g}), max dlnL {r['max_dlnL_vs_null']:+.2e};"
              f" 95% grid {r['upper95_grid']:.4e} refinado {r['upper95_refined']:.4e} (mu ahi {r['nuisance_at_upper95']:.4g});"
              f" alpha=0.8824: mu {az['mu']:.4g} dlnL {az['dlnL_vs_alpha0']:+.2f} dchi2 {az['dchi2']:.2f}"
              f" -> {az['nsigma_sqrt_dchi2']:.2f} sigma")
    t = res["Z1_zwane_mu_fixed"]["alpha_zwane_with_tanh"]
    print(f"Z1 con tanh (cov exacta de tanh, alpha=0.8824, mu=0.9804): std {t['std_tanh']:.4f}"
          f" dlnL {t['dlnL_vs_alpha0']:+.2f} dchi2 {t['dchi2']:.2f} -> {t['nsigma_sqrt_dchi2']:.2f} sigma")
    print("Z2 mu a lo largo del perfil (alpha >= 1e-3):")
    for p in res["Z2_zwane_mu_free"]["mu_along_profile"][::6]:
        print(f"   alpha {p['alpha']:.3e}  mu {p['mu']:.4g}")
    print("Escrito resultados/profile_desi_everpresent.json")


if __name__ == "__main__":
    main()
