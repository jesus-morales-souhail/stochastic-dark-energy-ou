#!/usr/bin/env python3
"""
Exploration only (TRECE.md): the sigma_X profile of perfil_fondo.py with the 13 DESI DR2 BAO numbers (D_V/r_d, D_M/r_d,
D_H/r_d with their 13x13 covariance, desi_dr2_data.load_gaussian_bao_vector) instead of the 7 alpha_DV.

Data in alpha units: alpha_i = obs_i / pred_i, pred from the DESI DR2 best-fit flat LCDM recovered in perfil_fondo.py
(Omega_m,f = 0.29740, h r_d,f = 101.538 Mpc); covariance diag(1/pred) C diag(1/pred).
Model: alpha_th,i(Omega_m, h r_d) = [(h r_d)_f / (h r_d)] Dhat_type(z_i; Omega_m) / Dhat_type(z_i; Omega_m,f), with
Dhat_M = int dz/E, Dhat_H = 1/E, Dhat_V = (z Dhat_M^2 Dhat_H)^(1/3) (flat LCDM, no radiation).
Residual: delta alpha_i = S_i X(z_i), S_type = -d ln D_type / d Omega_L (flat, H0 fixed, Omega_L = 0.6847: the sign
convention of the S_z table, which S_DV reproduces); kernel K_ij = S_i S_j sigma_X^2 exp(-theta |x_i - x_j|),
x = ln(1+z) of each row (rows at the same z: Delta x = 0).
Treatments as perfil_fondo.py: (a) background fixed, (b) Omega_m and h r_d profiled, (c) marginalised (linearised, flat
prior); bounds with theta free and at theta = 0.001, 0.1, 1, 10; same sigma_X grid, theta search, Delta lnL = -1.92
with brentq, extension to sigma_X = 1 if needed.
Check first: the same generic code fed with the 7 alpha_DV, the S_z table and the projected covariance reproduces the
table of FONDO.md (perfil_fondo.json).

Writes analisis_cmb/fondo/perfil_13.{json,txt}.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq, minimize_scalar

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import perfil_fondo as PF  # noqa: E402  (imports the DESI profile module and the reference background)
from desi_dr2_data import load_gaussian_bao_vector  # noqa: E402

THR, SIG = PF.THR, PF.SIG
OM_F, HRD_F = PF.OM_F, PF.HRD_F
OL_FID = 0.6847
C_KMS = 299792.458
OM_GRID = np.linspace(0.15, 0.50, 701)


def dhat(kind: str, om: float, z: float) -> float:
    E = lambda x: np.sqrt(om * (1 + x) ** 3 + 1 - om)
    if kind == "DH":
        return 1 / E(z)
    dm = quad(lambda x: 1 / E(x), 0, z, epsrel=1e-11)[0]
    return dm if kind == "DM" else (z * dm * dm / E(z)) ** (1 / 3)


def s_type(kind: str, z: float) -> float:
    """-d ln D_type / d Omega_L, flat, H0 fixed, at Omega_L = 0.6847 (Omega_m = 1 - Omega_L)."""
    e = 1e-4
    # d ln D / d Omega_L = [ln D(Omega_m = 1 - OL - e) - ln D(Omega_m = 1 - OL + e)] / 2e
    return -(np.log(dhat(kind, 1 - OL_FID - e, z)) - np.log(dhat(kind, 1 - OL_FID + e, z))) / (2 * e)


class Data:
    def __init__(self, z, kinds, alpha, cov, S):
        self.z, self.kinds, self.alpha, self.cov, self.S = np.asarray(z, float), list(kinds), np.asarray(alpha, float), \
            np.asarray(cov, float), np.asarray(S, float)
        self.x = np.log1p(self.z)
        self.dx = np.abs(self.x[:, None] - self.x[None, :])
        f = np.array([[dhat(k, om, zz) for k, zz in zip(self.kinds, self.z)] for om in OM_GRID])
        f0 = np.array([dhat(k, OM_F, zz) for k, zz in zip(self.kinds, self.z)])
        self.fspl = CubicSpline(OM_GRID, f / f0[None, :], axis=0)
        # linear templates at the sigma_X = 0 best fit (h r_d free) for (c)
        L = np.linalg.cholesky(self.cov)
        r = minimize_scalar(lambda om: -self.ll_om(L, om)[0], bounds=(0.15, 0.5), method="bounded", options={"xatol": 1e-8})
        self.om_b = float(r.x)
        self.s_b = self.ll_om(L, self.om_b)[1]
        fb, dfb = self.fspl(self.om_b), self.fspl(self.om_b, 1)
        self.T = np.column_stack([self.s_b * dfb, fb])
        self.r0 = self.alpha - self.s_b * fb
        self.cov_p = np.linalg.inv(self.T.T @ np.linalg.solve(self.cov, self.T))

    def kern(self, sigma, theta):
        return self.cov + np.outer(self.S, self.S) * sigma**2 * np.exp(-theta * self.dx)

    @staticmethod
    def solve(L, b):
        return np.linalg.solve(L.T, np.linalg.solve(L, b))

    def ll_om(self, L, om):
        f = self.fspl(om)
        cf, ca = self.solve(L, f), self.solve(L, self.alpha)
        s = (f @ ca) / (f @ cf)
        r = self.alpha - s * f
        return float(-0.5 * (r @ self.solve(L, r) + 2 * np.sum(np.log(np.diag(L))) + len(r) * np.log(2 * np.pi))), float(s)

    def ll_a(self, sigma, theta):
        L = np.linalg.cholesky(self.kern(sigma, theta))
        r = self.alpha - 1.0
        return float(-0.5 * (r @ self.solve(L, r) + 2 * np.sum(np.log(np.diag(L))) + len(r) * np.log(2 * np.pi)))

    def ll_b(self, sigma, theta):
        L = np.linalg.cholesky(self.kern(sigma, theta))
        res = minimize_scalar(lambda om: -self.ll_om(L, om)[0], bounds=(0.15, 0.5), method="bounded", options={"xatol": 1e-6})
        return -float(res.fun)

    def ll_c(self, sigma, theta):
        L = np.linalg.cholesky(self.kern(sigma, theta))
        cr, cT = self.solve(L, self.r0), self.solve(L, self.T)
        F = self.T.T @ cT
        b = self.T.T @ cr
        return float(-0.5 * (self.r0 @ cr - b @ np.linalg.solve(F, b) + 2 * np.sum(np.log(np.diag(L)))
                             + np.linalg.slogdet(F)[1] + len(self.r0) * np.log(2 * np.pi)))


def max_theta(fun, sigma):
    best_ll, best_th = -1e30, 1.0
    for br in [(1e-3, 0.5), (0.1, 3.0), (1.0, 15.0)]:
        r = minimize_scalar(lambda th: -fun(sigma, th), bounds=br, method="bounded", options={"xatol": 1e-5})
        if -r.fun > best_ll:
            best_ll, best_th = -float(r.fun), float(r.x)
    for th in np.geomspace(1e-3, 20, 40):
        v = fun(sigma, th)
        if v > best_ll:
            best_ll, best_th = v, float(th)
    return best_ll, best_th


def table(D: Data) -> dict:
    out = {}
    for name, fun in (("a_fixed", D.ll_a), ("b_profiled", D.ll_b), ("c_marginalised", D.ll_c)):
        d = {"theta_free": PF.bound(lambda s, f=fun: max_theta(f, s)[0])}
        if d["theta_free"]["sigma_95"]:
            d["theta_free"]["theta_at_bound"] = max_theta(fun, d["theta_free"]["sigma_95"])[1]
        for th in (0.001, 0.1, 1.0, 10.0):
            d[f"theta_{th:g}"] = PF.bound(lambda s, f=fun, t=th: f(s, t))
        out[name] = d
    ci = np.linalg.solve(D.cov, D.S)
    sv = float(1 / np.sqrt(D.S @ ci))
    out["s_mode"] = {"sigma_v": sv, "X_hat_over_sigma_v": float(ci @ (D.alpha - 1) * sv),
                     "occam_bound": float(np.sqrt(np.exp(3.84) - 1) * sv),
                     "sigma_Omega_m_hrd_free": float(np.sqrt(D.cov_p[0, 0])),
                     "sigma_Omega_m_hrd_fixed": float(1 / np.sqrt(D.T[:, 0] @ np.linalg.solve(D.cov, D.T[:, 0]))),
                     "Omega_m_best_sigma0": D.om_b, "h_rd_best_sigma0": HRD_F / D.s_b}
    return out


def main() -> None:
    keys = ("theta_free", "theta_0.001", "theta_0.1", "theta_1", "theta_10")
    # check: 7 alpha_DV through the generic code against perfil_fondo.json
    D7 = Data(PF.Z, ["DV"] * 7, PF.ALPHA, PF.CM, PF.S)
    t7 = table(D7)
    ref = json.loads((HERE / "perfil_fondo.json").read_text())["treatments"]
    chk = {}
    for tr in ("a_fixed", "b_profiled", "c_marginalised"):
        for k in keys:
            a, b = t7[tr][k]["sigma_95"], ref[tr][k]["sigma_95"]
            chk[f"{tr}|{k}"] = None if (a is None and b is None) else (abs(a / b - 1) if (a and b) else "mismatch")
    print("check 7 vs FONDO:", chk, flush=True)
    # the 13 numbers
    g = load_gaussian_bao_vector()
    kinds = [lab[:2] for lab in g["labels"]]
    z = np.asarray(g["z_rows"], float)
    pred = np.array([C_KMS / (100 * HRD_F) * dhat(k, OM_F, zz) for k, zz in zip(kinds, z)])
    alpha13 = g["mean"] / pred
    cov13 = g["cov"] / np.outer(pred, pred)
    S13 = np.array([s_type(k, zz) for k, zz in zip(kinds, z)])
    s_dv_check = float(np.max(np.abs(np.array([s_type("DV", zz) for zz in PF.Z]) - PF.S)))
    D13 = Data(z, kinds, alpha13, cov13, S13)
    t13 = table(D13)
    out = {"reference": {"Omega_m": OM_F, "h_rd_Mpc": HRD_F}, "check_S_DV_vs_table_max_abs": s_dv_check,
           "check_7_vs_FONDO_rel_diff": chk, "alpha13": alpha13.tolist(), "S13": S13.tolist(), "labels": g["labels"],
           "seven": t7, "thirteen": t13}
    (HERE / "perfil_13.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    fmt = lambda v: f"{v:.4g}" if v is not None else "sin cota"
    lines = [f"S_DV frente a la tabla S_z: max |dif| = {s_dv_check:.2e}",
             f"7 alpha_DV por este codigo frente a perfil_fondo.json: max dif relativa = "
             f"{max(v for v in chk.values() if isinstance(v, float)):.2e}; casos sin cota coinciden: "
             f"{all(v is None or isinstance(v, float) for v in chk.values())}", ""]
    for lab, t in (("7 alpha_DV", t7), ("13 numeros", t13)):
        lines.append(f"{lab}: sigma_v = {t['s_mode']['sigma_v']:.5f}, X_hat/sigma_v = {t['s_mode']['X_hat_over_sigma_v']:+.3f}, "
                     f"6.75 sigma_v = {t['s_mode']['occam_bound']:.5f}, sigma(Omega_m) h r_d libre = {t['s_mode']['sigma_Omega_m_hrd_free']:.4f}, "
                     f"h r_d fijo = {t['s_mode']['sigma_Omega_m_hrd_fixed']:.4f}, Omega_m(sigma=0) = {t['s_mode']['Omega_m_best_sigma0']:.4f}")
        for tr in ("a_fixed", "b_profiled", "c_marginalised"):
            lines.append(f"  {tr:16s} " + " ".join(f"{fmt(t[tr][k]['sigma_95']):>10s}" for k in keys))
    (HERE / "perfil_13.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
