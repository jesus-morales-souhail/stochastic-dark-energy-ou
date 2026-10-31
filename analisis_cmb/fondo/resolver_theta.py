#!/usr/bin/env python3
"""
Exploration only (RESOLVER_THETA.md): how large must sigma_X be for the data to tell two memories of the residual
apart? Expected separation between H1 = N(0, C0 + K(sigma_X, theta1)) and H2 = N(0, C0 + K(sigma_X, theta2)),
K_ij = S_i S_j sigma_X^2 exp(-theta |x_i - x_j|), x = ln(1+z), as the Kullback-Leibler divergence
    KL(1||2) = 1/2 [tr(S2^-1 S1) - n + ln det S2 - ln det S1]   (= expected Delta lnL if H1 is true)
in both directions, their mean ("symmetric": the expected Delta lnL averaged over which hypothesis is true, compared
with 1.92) and their sum (Jeffreys J). Pairs (theta1, theta2) = (0.001, 1), (0.001, 10), (1, 10); sigma_X on
geomspace(1e-3, 1e-1, 241); the sigma_X where each divergence reaches 1.92 by brentq on log sigma.

Surveys:
  DESI DR2  the 7 alpha_DV of scripts/profile_sigma_x_desi.py: z_eff, S_z table, official covariance projected (C0).
  Euclid-like  13 bins equispaced in ln(1+z) between z = 0.1 and 2.0, 0.5 % uncorrelated errors on alpha,
            S(z) = -d ln D_V / d Omega_L (flat, H0 fixed, Omega_L = 0.6847): the definition that reproduces the DESI S_z
            table (checked here at the DESI redshifts).
Background: fixed (as in the C2 profile), and marginalised (restricted likelihood: the data projected onto the
complement of the templates T = [S, 1], i.e. d alpha / d Omega_m and the h r_d scale, exact for the scale and to
first order for Omega_m; FONDO.md shows S and d ln alpha / d Omega_m have cosine 0.99996).

Writes analisis_cmb/fondo/resolver_theta.{json,txt}.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import profile_sigma_x_desi as PS  # noqa: E402  (loads the data; main() is not called)

THR = 1.92
PAIRS = [(0.001, 1.0), (0.001, 10.0), (1.0, 10.0)]
SIG = np.geomspace(1e-3, 1e-1, 241)
OL_FID = 0.6847


def s_of_z(z: np.ndarray) -> np.ndarray:
    e = 1e-4

    def lndv(ol, zz):
        E = lambda x: np.sqrt((1 - ol) * (1 + x) ** 3 + ol)
        dm = quad(lambda x: 1 / E(x), 0, zz, epsrel=1e-11)[0]
        return np.log((zz * dm * dm / E(zz)) ** (1 / 3))

    return -np.array([(lndv(OL_FID + e, zz) - lndv(OL_FID - e, zz)) / (2 * e) for zz in z])


def survey_desi():
    z = np.asarray(PS.z_eff, float)
    return {"z": z, "S": np.asarray(PS.S_z, float), "C0": np.asarray(PS.C_MEAS, float)}


def survey_euclid():
    x = np.linspace(np.log(1.1), np.log(3.0), 13)
    z = np.expm1(x)
    return {"z": z, "S": s_of_z(z), "C0": np.diag(np.full(13, 0.005**2))}


def projector(S: np.ndarray) -> np.ndarray:
    """Orthonormal basis of the complement of span{S, 1} (restricted likelihood)."""
    T = np.column_stack([S, np.ones_like(S)])
    q, _ = np.linalg.qr(np.column_stack([T, np.eye(len(S))]))
    return q[:, 2:]


def cov(sv, sigma, theta, P=None):
    x = np.log1p(sv["z"])
    C = sv["C0"] + np.outer(sv["S"], sv["S"]) * sigma**2 * np.exp(-theta * np.abs(x[:, None] - x[None, :]))
    return C if P is None else P.T @ C @ P


def kl(S1, S2) -> float:
    n = S1.shape[0]
    L2 = np.linalg.cholesky(S2)
    L1 = np.linalg.cholesky(S1)
    tr = np.trace(np.linalg.solve(S2, S1))
    return 0.5 * (tr - n + 2 * np.sum(np.log(np.diag(L2))) - 2 * np.sum(np.log(np.diag(L1))))


def divergences(sv, sigma, t1, t2, P=None) -> dict:
    A, B = cov(sv, sigma, t1, P), cov(sv, sigma, t2, P)
    k12, k21 = kl(A, B), kl(B, A)
    return {"KL_12": k12, "KL_21": k21, "sym_mean": 0.5 * (k12 + k21), "jeffreys": k12 + k21}


def crossing(sv, t1, t2, key, P=None):
    f = lambda ls: divergences(sv, float(np.exp(ls)), t1, t2, P)[key] - THR
    vals = np.array([f(np.log(s)) for s in SIG])
    k = np.nonzero(vals >= 0)[0]
    if not len(k):
        return None
    k = int(k[0])
    if k == 0:
        return float(SIG[0])
    return float(np.exp(brentq(f, np.log(SIG[k - 1]), np.log(SIG[k]), xtol=1e-12)))


def main() -> None:
    surveys = {"DESI_DR2": survey_desi(), "Euclid_like": survey_euclid()}
    d = surveys["DESI_DR2"]
    out = {"threshold": THR, "sigma_grid": [float(SIG[0]), float(SIG[-1]), len(SIG)],
           "check_S_definition_vs_DESI_table_max_abs": float(np.max(np.abs(s_of_z(d["z"]) - d["S"]))),
           "euclid_like": {"z": surveys["Euclid_like"]["z"].tolist(), "S": surveys["Euclid_like"]["S"].tolist(),
                           "alpha_error": 0.005},
           "results": {}}
    for name, sv in surveys.items():
        for bg, P in (("background_fixed", None), ("background_marginalised", projector(sv["S"]))):
            for t1, t2 in PAIRS:
                key = f"{name}|{bg}|{t1:g}_vs_{t2:g}"
                r = {k: crossing(sv, t1, t2, k, P) for k in ("KL_12", "KL_21", "sym_mean", "jeffreys")}
                r["at_sigma_1.5e-2"] = divergences(sv, 1.5e-2, t1, t2, P)
                r["at_sigma_2.633e-2"] = divergences(sv, 2.633e-2, t1, t2, P)
                r["curve_sym_mean"] = [divergences(sv, float(s), t1, t2, P)["sym_mean"] for s in SIG[::10]]
                out["results"][key] = r
                print(key, {k: r[k] for k in ("KL_12", "KL_21", "sym_mean", "jeffreys")},
                      "sym@1.5e-2 =", round(r["at_sigma_1.5e-2"]["sym_mean"], 4), flush=True)
    (HERE / "resolver_theta.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    fmt = lambda v: f"{v:.3g}" if v is not None else "> 0.1"
    lines = [f"S de la definicion frente a la tabla de DESI: max |dif| = {out['check_S_definition_vs_DESI_table_max_abs']:.2e}",
             "sigma_X en que la divergencia alcanza 1.92 (media simetrica | KL(1||2) | KL(2||1) | Jeffreys); "
             "y la media simetrica en sigma_X = 1.5e-2 y 2.633e-2", ""]
    for key, r in out["results"].items():
        lines.append(f"{key:52s} {fmt(r['sym_mean']):>8s} | {fmt(r['KL_12']):>8s} | {fmt(r['KL_21']):>8s} | {fmt(r['jeffreys']):>8s}"
                     f"   sym(1.5e-2) = {r['at_sigma_1.5e-2']['sym_mean']:.3g}, sym(2.633e-2) = {r['at_sigma_2.633e-2']['sym_mean']:.3g}")
    (HERE / "resolver_theta.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
