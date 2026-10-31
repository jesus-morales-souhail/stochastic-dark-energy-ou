#!/usr/bin/env python3
"""
Exploration only (SN.md): profile bound on sigma_X with Pantheon+ supernovae, with the same residual model as the DESI
profile (scripts/profile_sigma_x_desi.py):

  X(x), x = ln(1+z): delta Omega_L at redshift z; delta ln D_L(z_i) = S_i X(z_i), S_i = -d ln D_L / d Omega_L
  (flat, H0 fixed, Omega_L = 0.6847: the definition that reproduces the S_z table of scripts/desi_dr2_data.py to 0.1 %;
  the sign does not enter S_i S_j). In magnitudes the kernel is (5/ln 10)^2 S_i S_j sigma_X^2 exp(-theta |x_i - x_j|),
  added to the STAT+SYS covariance.

Data: datos/Pantheon+SH0ES.dat and datos/Pantheon+SH0ES_STAT+SYS.cov (PantheonPlusSH0ES/DataRelease), m_b_corr with
z_HD > 0.01 and IS_CALIBRATOR = 0. Model: flat LCDM, m = 5 log10[(1 + z_HEL) int_0^{z_HD} dz/E(z)] + M, M marginalised
analytically (flat prior): -2 lnL = r'C^-1 r - (1'C^-1 r)^2 / (1'C^-1 1) + ln det C + ln(1'C^-1 1); Omega_m free
(profiled in [0.1, 0.6] for every (sigma_X, theta)).
Profile as in the DESI script: sigma_X grid geomspace(1e-6, 5e-2, 48); theta maximised in (1e-3, 20) with three bounded
searches (1e-3, 0.5), (0.1, 3), (1, 15) plus a 40-point log grid; 95 %: Delta lnL = -1.92 from the profile maximum,
root by brentq between the bracketing grid nodes. Also the bound at fixed theta = 0.001, 0.1, 1, 10 (if the profile
does not cross -1.92 inside the grid, the grid is extended to sigma_X = 1 and this is said).

Run from analisis_cmb/sn/ with python3 -I (the data files are downloaded, untrusted input; they are only parsed as
numbers). Writes perfil_sn.{json,txt}; then `python3 -I perfil_sn.py --marg` appends the Omega_m-marginalised
diagnostic to both.
"""

from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")

import json
import multiprocessing as mp
import time
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import brentq, minimize_scalar

HERE = Path(__file__).resolve().parent
DAT = HERE / "datos" / "Pantheon+SH0ES.dat"
COV = HERE / "datos" / "Pantheon+SH0ES_STAT+SYS.cov"
OL_FID = 0.6847
THR = -1.92
SIG_GRID = np.geomspace(1e-6, 5e-2, 48)
MAG = 5.0 / np.log(10.0)


def load():
    with open(DAT) as f:
        head = f.readline().split()
    raw = np.genfromtxt(DAT, skip_header=1, dtype=None, encoding="ascii", names=head)
    n = int(open(COV).readline())
    cov = np.loadtxt(COV, skiprows=1).reshape(n, n)
    assert n == len(raw)
    m = (raw["zHD"] > 0.01) & (raw["IS_CALIBRATOR"] == 0)
    idx = np.nonzero(m)[0]
    return {"zhd": raw["zHD"][idx].astype(float), "zhel": raw["zHEL"][idx].astype(float),
            "mb": raw["m_b_corr"][idx].astype(float), "cov": cov[np.ix_(idx, idx)], "n_all": n, "n": len(idx)}


D = load()
ZG = np.concatenate([[0.0], np.geomspace(1e-4, D["zhd"].max() * 1.001, 4000)])


def dl_hat(om: float, z: np.ndarray, zhel: np.ndarray) -> np.ndarray:
    """(1 + z_hel) int_0^z dz/E, flat LCDM, cumulative trapezoid on a fine grid + linear interpolation."""
    f = 1.0 / np.sqrt(om * (1 + ZG) ** 3 + 1 - om)
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(ZG))])
    return (1 + zhel) * np.interp(z, ZG, cum)


def s_vector(z: np.ndarray) -> np.ndarray:
    """-d ln D_L / d Omega_L, flat, H0 fixed, at Omega_L = 0.6847 (D_L's (1+z) factor cancels)."""
    e = 1e-4

    def lndm(ol, zz):
        return np.log(quad(lambda x: 1 / np.sqrt((1 - ol) * (1 + x) ** 3 + ol), 0, zz, epsrel=1e-10)[0])

    return -np.array([(lndm(OL_FID + e, zz) - lndm(OL_FID - e, zz)) / (2 * e) for zz in z])


S = s_vector(D["zhd"])
X = np.log1p(D["zhd"])
DX = np.abs(X[:, None] - X[None, :])
SS = MAG**2 * np.outer(S, S)
ONE = np.ones(D["n"])


def loglike_fixed_C(cf, om: float) -> float:
    r = D["mb"] - 5 * np.log10(dl_hat(om, D["zhd"], D["zhel"]))
    ci_r = cho_solve(cf, r)
    ci_1 = cho_solve(cf, ONE)
    e = ONE @ ci_1
    b = ONE @ ci_r
    logdet = 2 * np.sum(np.log(np.diag(cf[0])))
    return -0.5 * (r @ ci_r - b * b / e + logdet + np.log(e))


def max_over_om(sigma: float, theta: float) -> tuple[float, float]:
    C = D["cov"] + SS * sigma**2 * np.exp(-theta * DX)
    try:
        cf = cho_factor(C, lower=True, check_finite=False)
    except np.linalg.LinAlgError:
        return -1e30, np.nan
    res = minimize_scalar(lambda om: -loglike_fixed_C(cf, om), bounds=(0.1, 0.6), method="bounded",
                          options={"xatol": 1e-5})
    return -float(res.fun), float(res.x)


def max_over_theta(sigma: float) -> tuple[float, float, float]:
    """As profile_sigma_x_desi.max_logL_over_theta, with Omega_m profiled inside."""
    best = (-1e30, 1.0, np.nan)

    def neg(th):
        return -max_over_om(sigma, th)[0]

    for br in [(1e-3, 0.5), (0.1, 3.0), (1.0, 15.0)]:
        r = minimize_scalar(neg, bounds=br, method="bounded", options={"xatol": 1e-5})
        if -r.fun > best[0]:
            best = (-float(r.fun), float(r.x), max_over_om(sigma, float(r.x))[1])
    for th in np.geomspace(1e-3, 20, 40):
        ll, om = max_over_om(sigma, th)
        if ll > best[0]:
            best = (ll, float(th), om)
    return best


def profile(fun, grid: np.ndarray, extend: bool) -> dict:
    with mp.get_context("fork").Pool(8) as pool:
        vals = pool.map(fun, list(grid))
    ll = np.array([v[0] for v in vals])
    llmax = float(ll.max())
    rel = ll - llmax
    above = np.nonzero(rel >= THR)[0]
    out = {"grid_max_sigma": float(grid[-1]), "ll_max": llmax, "sigma_at_max": float(grid[int(np.argmax(ll))]),
           "dlogL_rel": rel.tolist(), "theta_best": [float(v[1]) for v in vals], "om_best": [float(v[2]) for v in vals]}
    j = int(above[-1])
    if j + 1 < len(grid):
        g = lambda s: fun(float(s))[0] - llmax - THR
        s95 = float(brentq(g, grid[j], grid[j + 1], xtol=1e-10 * grid[j]))
        out.update(sigma_95=s95, crossing="inside grid", at_bound=list(fun(s95)))
    elif extend:
        # extension beyond the DESI grid, up to sigma_X = 1 (stated in SN.md)
        ext = np.geomspace(grid[-1], 1.0, 25)[1:]
        with mp.get_context("fork").Pool(8) as pool:
            lle = np.array([v[0] for v in pool.map(fun, list(ext))])
        out["extension_grid"] = ext.tolist()
        out["extension_dlogL_rel"] = (lle - llmax).tolist()
        k = np.nonzero(lle - llmax < THR)[0]
        if len(k):
            k = int(k[0])
            lo = grid[-1] if k == 0 else ext[k - 1]
            s95 = float(brentq(lambda s: fun(float(s))[0] - llmax - THR, lo, ext[k], xtol=1e-10 * lo))
            out.update(sigma_95=s95, crossing="beyond the DESI grid (extended to 1)", at_bound=list(fun(s95)))
        else:
            out.update(sigma_95=None, crossing="no crossing up to sigma_X = 1")
    else:
        out.update(sigma_95=None, crossing="no crossing")
    return out


def fixed_om_job(args):
    """Diagnostic: theta fixed and Omega_m fixed at the LCDM best fit (no Omega_m freedom)."""
    s, theta, om = args
    C = D["cov"] + SS * s**2 * np.exp(-theta * DX)
    cf = cho_factor(C, lower=True, check_finite=False)
    return loglike_fixed_C(cf, om), theta, om


def _fixed(theta):
    return lambda s: (*max_over_om(s, theta), theta)


def fixed_theta_job(args):
    s, theta = args
    ll, om = max_over_om(s, theta)
    return ll, theta, om


def main() -> None:
    t0 = time.time()
    out = {"data": {"file": DAT.name, "cov": COV.name, "n_all": D["n_all"], "n_used": D["n"],
                    "cut": "z_HD > 0.01 and IS_CALIBRATOR = 0", "z_range": [float(D["zhd"].min()), float(D["zhd"].max())]},
           "S_definition": "-dlnD_L/dOmega_L, flat, H0 fixed, Omega_L = 0.6847 (as the DESI S_z table)",
           "S_range": [float(S.min()), float(S.max())]}
    # LCDM baseline
    ll0, om0 = max_over_om(0.0, 1.0)
    out["lcdm"] = {"ll": ll0, "Omega_m": om0}
    print("LCDM:", out["lcdm"], flush=True)
    out["free_theta"] = profile(max_over_theta, SIG_GRID, extend=True)
    print("theta libre:", {k: v for k, v in out["free_theta"].items() if k not in ("dlogL_rel", "theta_best", "om_best")}, flush=True)
    out["fixed_theta"] = {}
    for th in (0.001, 0.1, 1.0, 10.0):
        fun = _fixed(th)
        with mp.get_context("fork").Pool(8) as pool:
            vals = pool.map(fixed_theta_job, [(float(s), th) for s in SIG_GRID])
        ll = np.array([v[0] for v in vals])
        llmax = float(ll.max())
        rel = ll - llmax
        above = np.nonzero(rel >= THR)[0]
        j = int(above[-1])
        d = {"ll_max": llmax, "sigma_at_max": float(SIG_GRID[int(np.argmax(ll))]), "dlogL_rel": rel.tolist()}
        if j + 1 < len(SIG_GRID):
            d["sigma_95"] = float(brentq(lambda s: fun(float(s))[0] - llmax - THR, SIG_GRID[j], SIG_GRID[j + 1],
                                         xtol=1e-10 * SIG_GRID[j]))
            d["crossing"] = "inside grid"
        else:
            # extension beyond the DESI grid (stated in SN.md)
            ext = np.geomspace(SIG_GRID[-1], 1.0, 25)[1:]
            lle = np.array([fun(float(s))[0] for s in ext])
            k = np.nonzero(lle - llmax < THR)[0]
            d["extension_grid"] = ext.tolist()
            d["extension_dlogL_rel"] = (lle - llmax).tolist()
            if len(k):
                k = int(k[0])
                lo = SIG_GRID[-1] if k == 0 else ext[k - 1]
                d["sigma_95"] = float(brentq(lambda s: fun(float(s))[0] - llmax - THR, lo, ext[k], xtol=1e-10 * lo))
                d["crossing"] = "beyond the DESI grid (extended to 1)"
            else:
                d["sigma_95"] = None
                d["crossing"] = "no crossing up to sigma_X = 1"
        out["fixed_theta"][str(th)] = d
        print("theta", th, {k: v for k, v in d.items() if "dlogL" not in k and "grid" not in k}, flush=True)
    # diagnostic: the theta -> 0 mode S_i X is, to first order, the Omega_m direction (S = d ln D_L / d Omega_m, flat,
    # H0 fixed), so it is degenerate with free Omega_m; repeat theta = 0.001 with Omega_m fixed at the LCDM best fit
    grid_d = np.geomspace(1e-6, 1.0, 72)
    with mp.get_context("fork").Pool(8) as pool:
        vals = pool.map(fixed_om_job, [(float(s_), 0.001, om0) for s_ in grid_d])
    ll = np.array([v[0] for v in vals])
    rel = ll - ll.max()
    j = int(np.nonzero(rel >= THR)[0][-1])
    f = lambda s_: fixed_om_job((float(s_), 0.001, om0))[0] - ll.max() - THR
    out["diagnostic_theta0.001_om_fixed"] = {"Omega_m": om0, "sigma_95": float(brentq(f, grid_d[j], grid_d[j + 1])) if j + 1 < len(grid_d) else None,
                                             "sigma_at_max": float(grid_d[int(np.argmax(ll))])}
    # alignment of S with the Omega_m direction (numerical d ln D_L / d Omega_m at the LCDM best fit)
    e = 1e-4
    dlm = (np.log(dl_hat(om0 + e, D["zhd"], D["zhel"])) - np.log(dl_hat(om0 - e, D["zhd"], D["zhel"]))) / (2 * e)
    out["diagnostic_S_vs_dlnDL_dOm"] = {"cosine": float(S @ dlm / np.linalg.norm(S) / np.linalg.norm(dlm)),
                                        "ratio_median": float(np.median(S / dlm))}
    print("diagnostico:", out["diagnostic_theta0.001_om_fixed"], out["diagnostic_S_vs_dlnDL_dOm"], flush=True)
    out["seconds"] = time.time() - t0
    (HERE / "perfil_sn.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    ft = out["free_theta"]
    lines = [f"Pantheon+ (m_b_corr, STAT+SYS), {D['n']} SNe (z_HD > 0.01, sin calibradores), M_B marginalizado, Omega_m libre",
             f"LCDM: Omega_m = {om0:.4f}",
             f"theta libre en (1e-3, 20): sigma_X < {ft['sigma_95']}  ({ft['crossing']}); en la cota (lnL, theta, Omega_m) = {ft.get('at_bound')}",
             f"  maximo del perfil en sigma_X = {ft['sigma_at_max']:.3g}"]
    for th, d in out["fixed_theta"].items():
        lines.append(f"theta = {th}: sigma_X < {d['sigma_95']}  ({d['crossing']}); maximo en sigma_X = {d['sigma_at_max']:.3g}")
    dg = out["diagnostic_theta0.001_om_fixed"]
    lines.append(f"diagnostico theta = 0.001 con Omega_m fijo en {dg['Omega_m']:.4f}: sigma_X < {dg['sigma_95']}")
    lines.append(f"coseno entre S y d ln D_L / d Omega_m: {out['diagnostic_S_vs_dlnDL_dOm']['cosine']:.6f}")
    (HERE / "perfil_sn.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))




def marg_diagnostic() -> dict:
    """Diagnostic appended to perfil_sn.json: Omega_m marginalised (flat prior) instead of profiled, linearised around
    the LCDM best fit: templates T = [1, (5/ln 10) d ln D_L / d Omega_m] marginalised analytically together with M,
    -2 lnL = r'C^-1 r - b'F^-1 b + ln det C + ln det F, b = T'C^-1 r, F = T'C^-1 T. Bound at fixed theta."""
    res = json.loads((HERE / "perfil_sn.json").read_text())
    om0 = res["lcdm"]["Omega_m"]
    e = 1e-4
    dl0 = dl_hat(om0, D["zhd"], D["zhel"])
    v = MAG * (np.log(dl_hat(om0 + e, D["zhd"], D["zhel"])) - np.log(dl_hat(om0 - e, D["zhd"], D["zhel"]))) / (2 * e)
    T = np.column_stack([ONE, v])
    r = D["mb"] - 5 * np.log10(dl0)

    def ll(s, theta):
        C = D["cov"] + SS * s**2 * np.exp(-theta * DX)
        cf = cho_factor(C, lower=True, check_finite=False)
        cr, cT = cho_solve(cf, r), cho_solve(cf, T)
        F = T.T @ cT
        b = T.T @ cr
        return -0.5 * (r @ cr - b @ np.linalg.solve(F, b) + 2 * np.sum(np.log(np.diag(cf[0]))) + np.linalg.slogdet(F)[1])

    grid = np.geomspace(1e-6, 1.0, 72)
    out = {"Omega_m_lin": om0}
    for th in (0.001, 1.0, 10.0):
        vals = np.array([ll(float(s), th) for s in grid])
        rel = vals - vals.max()
        below = np.nonzero(rel < THR)[0]
        d = {"sigma_at_max": float(grid[int(np.argmax(vals))]), "dlogL_rel_at_sigma_1": float(rel[-1])}
        if len(below):
            k = int(below[0])
            d["sigma_95"] = float(brentq(lambda s: ll(s, th) - vals.max() - THR, grid[k - 1], grid[k]))
        else:
            d["sigma_95"] = None
        out[str(th)] = d
        print("marginal Omega_m, theta", th, d, flush=True)
    res["diagnostic_om_marginalised"] = out
    (HERE / "perfil_sn.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    with open(HERE / "perfil_sn.txt", "a", encoding="utf-8") as f:
        for th in ("0.001", "1.0", "10.0"):
            f.write(f"diagnostico Omega_m marginalizado (lineal), theta = {th}: sigma_X < {out[th]['sigma_95']}"
                    f" (Delta lnL en sigma_X = 1: {out[th]['dlogL_rel_at_sigma_1']:.3f})\n")
    return out


if __name__ == "__main__":
    import sys

    if sys.argv[1:] == ["--marg"]:
        marg_diagnostic()
    else:
        main()
