#!/usr/bin/env python3
"""
Exploration only (FONDO.md): the DESI DR2 sigma_X profile of scripts/profile_sigma_x_desi.py (same 7 alpha_DV, same
projected official covariance, same sigma_X grid, same theta search, Delta lnL = -1.92 with brentq) with the mean
background free.

The alpha_DV of figure6 are (D_V/r_d) / (D_V/r_d)^best, relative to the DESI DR2 best-fit flat LCDM (plot_figure6.ipynb
of the data release). That reference is recovered from the data themselves: fitting obs D_V/r_d (from the 13-vector)
/ alpha gives Omega_m,f = 0.29740 and h r_d,f = 101.538 Mpc, reproducing the 7 alpha to <= 3e-4. With it,
  alpha_th(Omega_m, h r_d) = [(h r_d)_f / (h r_d)] * Dhat_V(z; Omega_m) / Dhat_V(z; Omega_m,f),  Dhat_V = D_V H0 / c,
flat LCDM without radiation, so alpha_th = 1 at the reference.

  (a) background fixed (alpha_th = 1): scripts/profile_sigma_x_desi.py's own functions; must give 2.633e-2.
  (b) Omega_m and h r_d profiled at every (sigma_X, theta): h r_d analytically (alpha_th is linear in 1/(h r_d)),
      Omega_m by bounded Brent in [0.15, 0.5].
  (c) Omega_m and h r_d marginalised (flat prior), linearised around the best fit at sigma_X = 0, as in
      analisis_cmb/sn/perfil_sn.py: -2 lnL = r'C^-1 r - b'F^-1 b + ln det C + ln det F, T = [d alpha/d Omega_m, alpha_th].
For each: bound with theta free and at theta = 0.001, 0.1, 1, 10. If a profile does not cross -1.92 inside the
grid (<= 5e-2), it is extended to sigma_X = 1 and this is said.
Also: cosine between S (the S_z table of desi_dr2_data.py) and d ln alpha / d Omega_m at the reference, the best-fit
Omega_m with h r_d free and sigma_X = 0, and sigma(Omega_m) from the Fisher matrix of the 7 alpha.

Writes analisis_cmb/fondo/perfil_fondo.{json,txt}; nothing in results/ or figures/.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq, least_squares, minimize_scalar

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import profile_sigma_x_desi as PS  # noqa: E402  (module import loads the data; main() is not called)
from desi_dr2_data import load_gaussian_bao_vector  # noqa: E402

THR = -1.92
SIG = np.geomspace(1e-6, 5e-2, 48)
C_KMS = 299792.458
Z = np.asarray(PS.z_eff, dtype=float)
ALPHA = np.asarray(PS.alpha, dtype=float)
S = np.asarray(PS.S_z, dtype=float)
X = np.asarray(PS.x, dtype=float)
CM = np.asarray(PS.C_MEAS, dtype=float)


def dvhat(om: float, z: float) -> float:
    E = lambda x: np.sqrt(om * (1 + x) ** 3 + 1 - om)
    dm = quad(lambda x: 1 / E(x), 0, z, epsrel=1e-11)[0]
    return (z * dm * dm / E(z)) ** (1 / 3)


def reference() -> tuple[float, float, float]:
    g = load_gaussian_bao_vector()
    m = g["mean"]
    obs = np.array([m[0]] + [(Z[i] * m[2 * i - 1] ** 2 * m[2 * i]) ** (1 / 3) for i in range(1, 6)]
                   + [(Z[6] * m[12] ** 2 * m[11]) ** (1 / 3)])
    r = least_squares(lambda p: obs / (C_KMS / (100 * p[1]) * np.array([dvhat(p[0], z) for z in Z])) - ALPHA,
                      [0.3, 101.5], xtol=1e-14, ftol=1e-14)
    return float(r.x[0]), float(r.x[1]), float(np.max(np.abs(r.fun)))


OM_F, HRD_F, REF_RES = reference()
OM_GRID = np.linspace(0.15, 0.50, 701)
_DV = np.array([[dvhat(om, z) for z in Z] for om in OM_GRID])
_DVF = np.array([dvhat(OM_F, z) for z in Z])
F_SPL = CubicSpline(OM_GRID, _DV / _DVF[None, :], axis=0)  # f_i(Omega_m), f(OM_F) = 1


def kern(sigma: float, theta: float) -> np.ndarray:
    return CM + np.outer(S, S) * sigma**2 * np.exp(-theta * np.abs(X[:, None] - X[None, :]))


def chol(C):
    try:
        return np.linalg.cholesky(C)
    except np.linalg.LinAlgError:
        return None


def solve(L, b):
    return np.linalg.solve(L.T, np.linalg.solve(L, b))


def ll_profiled_om(L, om: float) -> tuple[float, float]:
    """lnL with h r_d profiled analytically at fixed Omega_m (constant terms as in logL_gaussian)."""
    f = F_SPL(om)
    cf, ca = solve(L, f), solve(L, ALPHA)
    s = (f @ ca) / (f @ cf)
    r = ALPHA - s * f
    return float(-0.5 * (r @ solve(L, r) + 2 * np.sum(np.log(np.diag(L))) + len(r) * np.log(2 * np.pi))), float(s)


def ll_b(sigma: float, theta: float) -> float:
    L = chol(kern(sigma, theta))
    if L is None:
        return -1e30
    res = minimize_scalar(lambda om: -ll_profiled_om(L, om)[0], bounds=(0.15, 0.5), method="bounded",
                          options={"xatol": 1e-6})
    return -float(res.fun)


def lin_templates():
    """Best fit at sigma_X = 0 (h r_d free) and the linear templates for (c)."""
    L = chol(CM)
    res = minimize_scalar(lambda om: -ll_profiled_om(L, om)[0], bounds=(0.15, 0.5), method="bounded",
                          options={"xatol": 1e-8})
    om, s = float(res.x), ll_profiled_om(L, float(res.x))[1]
    f, df = F_SPL(om), F_SPL(om, 1)
    T = np.column_stack([s * df, f])  # d alpha / d Omega_m, d alpha / d s  (s = (h r_d)_f / h r_d)
    r0 = ALPHA - s * f
    Fis = T.T @ np.linalg.solve(CM, T)
    cov = np.linalg.inv(Fis)
    return om, s, T, r0, cov


OM_B, S_B, T_C, R_C, COV_P = lin_templates()


def ll_c(sigma: float, theta: float) -> float:
    L = chol(kern(sigma, theta))
    if L is None:
        return -1e30
    cr, cT = solve(L, R_C), solve(L, T_C)
    Fm = T_C.T @ cT
    b = T_C.T @ cr
    return float(-0.5 * (R_C @ cr - b @ np.linalg.solve(Fm, b) + 2 * np.sum(np.log(np.diag(L)))
                         + np.linalg.slogdet(Fm)[1] + len(R_C) * np.log(2 * np.pi)))


def ll_a(sigma: float, theta: float) -> float:
    return PS.logL(theta, sigma)


def max_theta(fun, sigma: float) -> tuple[float, float]:
    """Same search as PS.max_logL_over_theta: three bounded searches and a 40-point log grid in (1e-3, 20)."""
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


def bound(g) -> dict:
    """g(sigma) -> lnL (already maximised over what is free). Profile on SIG, root; extension to 1 if needed."""
    ll = np.array([g(float(s)) for s in SIG])
    llmax = float(ll.max())
    rel = ll - llmax
    out = {"sigma_at_max": float(SIG[int(np.argmax(ll))]), "dlnL_max_minus_at_grid_start": float(llmax - ll[0])}
    j = int(np.nonzero(rel >= THR)[0][-1])
    if j + 1 < len(SIG):
        out.update(sigma_95=float(brentq(lambda s: g(float(s)) - llmax - THR, SIG[j], SIG[j + 1], xtol=1e-10 * SIG[j])),
                   crossing="inside grid")
        return out
    ext = np.geomspace(SIG[-1], 1.0, 25)[1:]
    lle = np.array([g(float(s)) for s in ext])
    k = np.nonzero(lle - llmax < THR)[0]
    out["dlnL_rel_at_sigma_1"] = float(lle[-1] - llmax)
    if len(k):
        k = int(k[0])
        lo = SIG[-1] if k == 0 else ext[k - 1]
        out.update(sigma_95=float(brentq(lambda s: g(float(s)) - llmax - THR, lo, ext[k], xtol=1e-10 * lo)),
                   crossing="beyond the DESI grid (extended to 1)")
    else:
        out.update(sigma_95=None, crossing="no crossing up to sigma_X = 1")
    return out


def main() -> None:
    ref = json.loads((ROOT / "results" / "profile_sigma_x" / "profile_sigma_x.json").read_text())
    dlna = F_SPL(OM_F, 1)  # d ln alpha / d Omega_m at the reference (h r_d fixed; f = 1 there)
    out = {"reference_background": {"Omega_m": OM_F, "h_rd_Mpc": HRD_F, "max_abs_alpha_residual": REF_RES,
                                    "note": "alpha_DV of figure6 are relative to the DESI DR2 best-fit LCDM (plot_figure6.ipynb)"},
           "cos_S_dlnalpha_dOm": float(S @ dlna / np.linalg.norm(S) / np.linalg.norm(dlna)),
           "fit_sigma0_hrd_free": {"Omega_m": OM_B, "h_rd_Mpc": HRD_F / S_B, "sigma_Omega_m_fisher": float(np.sqrt(COV_P[0, 0])),
                                   "desi_dr2_paper": "Omega_m = 0.2975 +- 0.0086, h r_d = 101.54 +- 0.73 Mpc (DESI only, "
                                                     "13 BAO points; arXiv:2503.14738 Sec. 6, eq. constraints_lcdm_desi)"},
           "treatments": {}}
    # the theta -> 0 mode: conditional error of its amplitude (everything else fixed), its estimate from the alpha,
    # and the bound a pure ln-det (Occam) penalty would give, sqrt(e^3.84 - 1) sigma_v (if the estimate were 0)
    ci_s = np.linalg.solve(CM, S)
    sv = float(1 / np.sqrt(S @ ci_s))
    out["s_mode"] = {"sigma_v_conditional": sv, "X_hat": float(ci_s @ (ALPHA - 1) * sv**2),
                     "X_hat_over_sigma_v": float(ci_s @ (ALPHA - 1) * sv), "occam_bound": float(np.sqrt(np.exp(3.84) - 1) * sv),
                     "sigma_Omega_m_marginalised_over_hrd": float(np.sqrt(COV_P[0, 0])),
                     "sigma_Omega_m_conditional_hrd_fixed": float(1 / np.sqrt(T_C[:, 0] @ np.linalg.solve(CM, T_C[:, 0])))}
    print(json.dumps({k: out[k] for k in ("reference_background", "cos_S_dlnalpha_dOm", "fit_sigma0_hrd_free")}, indent=1), flush=True)
    for name, fun in (("a_fixed", ll_a), ("b_profiled", ll_b), ("c_marginalised", ll_c)):
        d = {"theta_free": bound(lambda s, f=fun: max_theta(f, s)[0])}
        if d["theta_free"]["sigma_95"]:
            d["theta_free"]["theta_at_bound"] = max_theta(fun, d["theta_free"]["sigma_95"])[1]
        for th in (0.001, 0.1, 1.0, 10.0):
            d[f"theta_{th:g}"] = bound(lambda s, f=fun, t=th: f(s, t))
        out["treatments"][name] = d
        print(name, json.dumps(d), flush=True)
    a = out["treatments"]["a_fixed"]["theta_free"]["sigma_95"]
    out["check_a_reproduces_reference"] = {"ours": a, "reference": ref["sigma_X_95CL_upper_profile"],
                                           "identical": a == ref["sigma_X_95CL_upper_profile"]}
    (HERE / "perfil_fondo.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    fmt = lambda v: f"{v:.4g}" if v is not None else "sin cota"
    lines = [f"fondo de referencia de los alpha (mejor ajuste DESI DR2): Omega_m = {OM_F:.5f}, h r_d = {HRD_F:.3f} Mpc",
             f"cos(S, d ln alpha / d Omega_m) = {out['cos_S_dlnalpha_dOm']:.6f}",
             f"sigma_X = 0, h r_d libre: Omega_m = {OM_B:.4f} +- {np.sqrt(COV_P[0, 0]):.4f} (Fisher), h r_d = {HRD_F / S_B:.2f} Mpc",
             f"(a) reproduce 2.633e-2: {out['check_a_reproduces_reference']}",
             f"modo S (theta -> 0): sigma_v condicional = {out['s_mode']['sigma_v_conditional']:.5f}, X_hat/sigma_v = "
             f"{out['s_mode']['X_hat_over_sigma_v']:+.3f}, cota de Occam sqrt(e^3.84 - 1) sigma_v = {out['s_mode']['occam_bound']:.5f}; "
             f"sigma(Omega_m) con h r_d libre = {out['s_mode']['sigma_Omega_m_marginalised_over_hrd']:.4f}, con h r_d fijo = "
             f"{out['s_mode']['sigma_Omega_m_conditional_hrd_fixed']:.4f}", "",
             f"{'tratamiento':16s} {'theta libre':>12s} {'theta=0.001':>12s} {'theta=0.1':>12s} {'theta=1':>12s} {'theta=10':>12s}"]
    for name, d in out["treatments"].items():
        lines.append(f"{name:16s} " + " ".join(f"{fmt(d[k]['sigma_95']):>12s}" for k in
                                                ("theta_free", "theta_0.001", "theta_0.1", "theta_1", "theta_10")))
    (HERE / "perfil_fondo.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
