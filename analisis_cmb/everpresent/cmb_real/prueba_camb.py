#!/usr/bin/env python3
"""
Exploration only (EXPLORACION.md, point 3): can the history of a Model 2 survivor be passed to CAMB with
DarkEnergyPPF.set_w_a_table, as ou_cmb.py does? Run with analisis_cmb/venv (CAMB 2.0.4):

    ../../venv/bin/python prueba_camb.py

Two survivors (cmb_real/supervivientes.json):
  5492   Omega_de > 0 on every stored point (z in [0, 3] and [1000, 1e5]); passes A in 3d (omega_m free).
  11555  Omega_de changes sign at z < 3 (min -0.23, Omega_de(0) = +0.56); passes A in 3e(b), and is the
         only survivor that passes the window proxy.

The stored history has a gap: bao_sim.py sampled lambda = ln a in [-ln 4, 0] (+ lambda_rec) and rd_real.py
added [-ln 1e5, -ln 1000]; z in (3, 1000) was never sampled, because BAO and r_d do not need it. The CMB
does (distance to last scattering, ISW). Here the gap is filled by a conditional Gaussian draw given all
458 stored points (new seed GAP_SEED): it is a new random choice, not part of the stored realisation.
The probability that the gap stays positive for 5492 is estimated from N_GAP draws.

Background: H^2 = 100^2 [w_m a^-3 + w_r a^-4] / (1 - Omega_de) (bao_sim.py), so
rho_DE / (3 (100 km/s/Mpc)^2 / 8 pi G) = Omega_de / (1 - Omega_de) (w_m a^-3 + w_r a^-4), and
w = -1 - (1/3) d ln rho_DE / d ln a, defined only while rho_DE keeps one sign.
CAMB: omega_b, omega_c (+ m_nu = 0.06) so that omega_m equals the survivor's omega_m, H0 = H(z = 0) of the
realisation, flat, DarkEnergyPPF with the w table on a in [1e-5, 1].

Writes cmb_real/prueba_camb.json and prints a log.
"""

from __future__ import annotations

import json
import sys
import time
import warnings
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline

HERE = Path(__file__).resolve().parent
EV = HERE.parent
sys.path.insert(0, str(EV))
sys.path.insert(0, str(EV.parent))

import bao_sim as B  # noqa: E402
import camb  # noqa: E402
import ou_cmb as M  # noqa: E402
from camb.dark_energy import DarkEnergyPPF  # noqa: E402

MU, ALPHA = 0.9804, 0.8824
GAP_SEED = 20261031
N_GAP = 500
OMB, OMNU = 0.02237, 0.0006448993754203254  # rd_planck.json
LNA_MIN = np.log(1e-5)
N_TAB = 4000


def stored_points() -> tuple[np.ndarray, np.ndarray]:
    """lambda (458,) and Omega_hat (458, 20000): the points of bao_sim.py and rd_real.py, same seeds."""
    L = B.chol_m2(MU)
    Ohat = ALPHA * L @ np.random.default_rng(B.SEED).standard_normal((401, 20000))
    lam_rec = B.LAM_ALL[400]
    zw = np.geomspace(1000.0, 5000.0, 10)
    lam_new = np.unique(np.concatenate([np.arange(-np.log(1e5), -np.log(1000.0) + 1e-9, 0.1), -np.log1p(zw)]))
    k = np.exp(-MU * (lam_new - lam_rec) ** 2 / 2)
    Knn = np.exp(-MU * (lam_new[:, None] - lam_new[None, :]) ** 2 / 2) - np.outer(k, k)
    w, V = np.linalg.eigh(Knn)
    Ls = V * np.sqrt(np.clip(w, 0, None))[None, :]
    Zn = np.random.default_rng(B.SEED + 1).standard_normal((len(lam_new), 20000))
    Onew = ALPHA * (k[:, None] * (Ohat[400][None, :] / ALPHA) + Ls @ Zn)
    return np.concatenate([B.LAM_ALL, lam_new]), np.vstack([Ohat, Onew])


def gap_sampler(lam_obs: np.ndarray, y_obs: np.ndarray, lam_gap: np.ndarray, nugget: float = 1e-8):
    """Conditional mean and a square root of the conditional covariance of Omega_hat on lam_gap."""
    kern = lambda x, y: ALPHA**2 * np.exp(-MU * (x[:, None] - y[None, :]) ** 2 / 2)
    Koo = kern(lam_obs, lam_obs) + nugget * ALPHA**2 * np.eye(len(lam_obs))
    Kgo = kern(lam_gap, lam_obs)
    cf = np.linalg.cholesky(Koo)
    sol = np.linalg.solve(cf.T, np.linalg.solve(cf, np.column_stack([y_obs, Kgo.T])))
    mean = Kgo @ sol[:, 0]
    cov = kern(lam_gap, lam_gap) - Kgo @ sol[:, 1:]
    w, V = np.linalg.eigh(0.5 * (cov + cov.T))
    return mean, V * np.sqrt(np.clip(w, 0, None))[None, :]


def history(lam_obs, y_obs, lam_gap, y_gap, wm):
    """Omega_de, rho_DE and w on a uniform ln a grid [ln 1e-5, 0]; cubic spline in lambda for Omega_hat."""
    lam = np.concatenate([lam_obs, lam_gap])
    y = np.concatenate([y_obs, y_gap])
    o = np.argsort(lam)
    lam, y = lam[o], y[o]
    keep = np.concatenate([[True], np.diff(lam) > 1e-9])
    lna = np.linspace(LNA_MIN, 0.0, N_TAB)
    ode = np.tanh(CubicSpline(lam[keep], y[keep])(lna))
    a = np.exp(lna)
    rho = ode / (1.0 - ode) * (wm / a**3 + B.ORAD / a**4)  # units 3 (100 km/s/Mpc)^2 / 8 pi G
    with np.errstate(divide="ignore", invalid="ignore"):
        w = -1.0 - np.gradient(np.log(np.abs(rho)), lna) / 3.0
    return a, ode, rho, w


def camb_try(label: str, a, w, wm, H0, ode, extra_runs: bool) -> dict:
    """Background + (optionally) TT l<30 / C_L^phiphi at the realisation's H0 and the theta_* solve."""
    out = {"label": label, "H0_input": H0, "omega_m": wm, "w_min": float(np.nanmin(w)), "w_max": float(np.nanmax(w)),
           "w_nonfinite": int(np.sum(~np.isfinite(w)))}

    def make(h):
        pars = M.make_params(H0=h, halofit="takahashi")
        pars.set_cosmology(H0=h, ombh2=OMB, omch2=wm - OMB - OMNU, mnu=0.06, omk=0.0, tau=M.PLANCK18["tau"])
        de = DarkEnergyPPF()
        de.set_w_a_table(a, w)
        pars.DarkEnergy = de
        return pars

    with warnings.catch_warnings(record=True) as wr:
        warnings.simplefilter("always")
        try:
            t0 = time.time()
            pars = make(H0)
            bg = camb.get_background(pars)
            out["t_background_s"] = time.time() - t0
            zc = np.array([0.0, 0.5, 1.0, 2.0, 3.0, 10.0, 100.0, 1000.0, 1100.0, 3000.0, 1e4])
            Hc = bg.hubble_parameter(zc)
            ai = 1.0 / (1.0 + zc)
            odi = np.interp(np.log(ai), np.log(a), ode)
            Hw = 100.0 * np.sqrt((wm / ai**3 + B.ORAD / ai**4) / (1.0 - odi))
            out["H_camb_over_H_realisation"] = {f"{z:g}": float(x) for z, x in zip(zc, Hc / Hw)}
            out["omega_de_camb_z0"] = float(bg.omega_de)
            out["Ode0_realisation"] = float(ode[-1])
            out["thetastar100_at_H0"] = float(bg.get_derived_params()["thetastar"])
            out["background"] = "ok" if np.all(np.isfinite(Hc)) else "NaN in H(z)"
        except Exception as exc:  # noqa: BLE001
            out["background"] = f"{type(exc).__name__}: {exc}".replace("\n", " ")
            out["warnings"] = [str(x.message) for x in wr]
            return out
        if extra_runs:
            try:
                t0 = time.time()
                r = M.run(make(H0), keep_results=False)
                out["t_full_run_s"] = time.time() - t0
                out["TT_l2_29_finite"] = bool(np.all(np.isfinite(r["DlTT"][2:30])))
                out["Clpp_finite"] = bool(np.all(np.isfinite(r["Clpp"][2:2001])))
                out["DlTT_l2_10_20_29"] = [float(r["DlTT"][l]) for l in (2, 10, 20, 29)]
            except Exception as exc:  # noqa: BLE001
                out["full_run"] = f"{type(exc).__name__}: {exc}".replace("\n", " ")
            try:
                target = M.thetastar_background(M.make_params(halofit="takahashi"))
                t0 = time.time()
                h, th, it = M.solve_H0_for_thetastar(target, make, lo=30.0, hi=120.0)
                out["thetastar_solve"] = {"H0": h, "thetastar100": th, "iterations": it, "seconds": time.time() - t0,
                                          "target": target}
            except Exception as exc:  # noqa: BLE001
                out["thetastar_solve"] = f"{type(exc).__name__}: {exc}".replace("\n", " ")
        out["warnings"] = sorted({str(x.message)[:300] for x in wr})
    return out


def main() -> None:
    sv = {r["index"]: r for r in json.loads((HERE / "supervivientes.json").read_text())["model2"]}
    lam_all, O = stored_points()
    lam_gap = np.arange(-np.log(1000.0) + 0.02, -np.log(4.0) - 0.01, 0.02)
    res = {"gap": {"lambda_range": [float(lam_gap[0]), float(lam_gap[-1])], "n": int(len(lam_gap)),
                   "seed": GAP_SEED, "z_range": [float(np.exp(-lam_gap[-1]) - 1), float(np.exp(-lam_gap[0]) - 1)]},
           "cases": []}
    rng = np.random.default_rng(GAP_SEED)
    for j, crit, sign_note in ((5492, "3d_wm_free", "positive on stored points"),
                               (11555, "3e_b_rdreal_wm_free", "sign change at z < 3")):
        wm = sv[j]["passes"][crit]["omega_m"]
        mean, Ls = gap_sampler(lam_all, O[:, j], lam_gap)
        draws = mean[:, None] + Ls @ rng.standard_normal((len(lam_gap), N_GAP))
        p_pos = float(np.mean(np.all(draws > 0, axis=0)))
        # continuity check: conditional mean next to the stored edges
        edge = {"mean_at_z~3": float(mean[-1]), "stored_at_z=3": float(O[0, j]),
                "mean_at_z~1000": float(mean[0]), "stored_at_z=1000": float(O[np.argmin(np.abs(lam_all + np.log(1000.0))), j])}
        k_pos = np.nonzero(np.all(draws > 0, axis=0))[0]
        pick = int(k_pos[0]) if len(k_pos) else 0
        a, ode, rho, w = history(lam_all, O[:, j], lam_gap, draws[:, pick], wm)
        H0 = 100.0 * np.sqrt((wm + B.ORAD) / (1.0 - ode[-1]))
        case = {"index": j, "criterion": crit, "note": sign_note, "omega_m": wm, "H0": H0,
                "gap_draws_all_positive_frac": p_pos, "gap_pick": pick, "gap_edge_check": edge,
                "Ode_min": float(ode.min()), "Ode_max": float(ode.max()),
                "rho_sign_changes": int(np.sum(np.diff(np.sign(rho)) != 0))}
        print(f"== {j} ({sign_note}): omega_m {wm:.4f}, H0 {H0:.2f}, gap all-positive in {p_pos:.1%} of {N_GAP} draws; "
              f"Ode in [{ode.min():+.2f}, {ode.max():+.2f}], sign changes of rho_DE: {case['rho_sign_changes']}", flush=True)
        if case["rho_sign_changes"] == 0:
            case["camb"] = [camb_try("w from ln rho", a, w, wm, H0, ode, extra_runs=True)]
        else:
            cross = np.zeros(len(rho), bool)
            ic = np.nonzero(np.diff(np.sign(rho)) != 0)[0]
            cross[ic] = cross[ic + 1] = True  # the two grid points around each zero of rho_DE
            case["camb"] = [
                camb_try("naive: w from ln|rho| (large but finite next to the crossings)", a, w, wm, H0, ode,
                         extra_runs=True),
                camb_try("w = +inf next to the crossings", a, np.where(cross, np.inf, w), wm, H0, ode, False),
                camb_try("w = NaN next to the crossings", a, np.where(cross, np.nan, w), wm, H0, ode, False),
            ]
        for c in case["camb"]:
            print(json.dumps(c, indent=1), flush=True)
        res["cases"].append(case)
    (HERE / "prueba_camb.json").write_text(json.dumps(res, indent=1), encoding="utf-8")


if __name__ == "__main__":
    main()
