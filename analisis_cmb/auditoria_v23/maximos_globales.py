#!/usr/bin/env python3
"""
Recomputation of the linear marginalised bounds of CMB_FONDO.md (all data sets and treatments, both conventions) and of the
linear bounds of EXACTO_SN.md (three combinations) with the drop of 1.92 / 1.35 measured from the GLOBAL maximum of the profile
over an extended sigma_X grid (perfil_fondo.SIG up to 0.05, then geomspace(0.05, 1, 40)), instead of the maximum over
perfil_fondo.SIG (sigma_X <= 0.05) that perfil_fondo.bound uses. Reports, per entry: the historical bound, sigma at the global
maximum, ln L(max) - ln L(sigma -> 0), the new bound and the change. The historical value is reproduced first with the clipped
maximum as a check. Writes maximos_globales.json. Nothing of the historical files is modified.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cmb_fondo as CF  # noqa: E402
import exacto_sn as S  # noqa: E402

SIG = CF.V.SIG
EXT = np.geomspace(SIG[-1], 1.0, 40)[1:]
GRID = np.concatenate([SIG, EXT])
THRS = {"192": -1.92, "135": -1.35}


def profile(g):
    return np.array([g(float(s)) for s in GRID])


def bound_from(ll, g, thr, clip_max: bool):
    """Bound from the profile ll on GRID: drop thr from the max over SIG (clip_max) or over GRID; crossing from the max outward."""
    ref_region = ll[: len(SIG)] if clip_max else ll
    llmax = float(ref_region.max())
    k0 = int(np.argmax(ref_region))
    below = np.nonzero((ll - llmax < thr) & (np.arange(len(ll)) > k0))[0]
    if len(below) == 0:
        return None
    k = int(below[0])
    return float(brentq(lambda s: g(float(s)) - llmax - thr, GRID[k - 1], GRID[k], xtol=1e-10 * GRID[k - 1]))


def entry(g, old_192, old_135):
    ll = profile(g)
    k = int(np.argmax(ll))
    rec = {"sigma_at_max_grid": float(GRID[k]), "lnL_max_minus_sigma0": float(ll[k] - ll[0]),
           "max_class": "sigma->0" if k == 0 else ("interior" if k < len(GRID) - 1 else "upper edge of extended grid (1.0)"),
           "clipped_max_at_grid_edge_0.05": bool(int(np.argmax(ll[: len(SIG)])) == len(SIG) - 1)}
    for name, thr in THRS.items():
        rec[f"{name}_clipped"] = bound_from(ll, g, thr, True)
        rec[f"{name}_global"] = bound_from(ll, g, thr, False)
    rec["192_historical"], rec["135_historical"] = old_192, old_135
    for name in ("192", "135"):
        h, c, gl = rec[f"{name}_historical"], rec[f"{name}_clipped"], rec[f"{name}_global"]
        rec[f"{name}_reproduced"] = (h is None and c is None) or (h is not None and c is not None and abs(c / h - 1) < 2e-3)
        rec[f"{name}_change_pct"] = (100 * (gl / h - 1)) if (h and gl) else None
    return rec


def main():
    t0 = time.time()
    old_cf = json.load(open(HERE / "cmb_fondo.json"))["fits"]
    out = {"grid": {"sig_max_historical": float(SIG[-1]), "extended_to": 1.0}, "cmb_fondo": {}, "exacto_sn": {}}
    R = CF.Response()
    Rl = CF.Response(lambda_absorbs_x0=True)
    sets = {k: (v, R) for k, v in CF.SETS.items()}
    sets["BAO+CMB|Lambda_absorbs_X0"] = (CF.SETS["BAO+CMB"], Rl)
    sets["BAO|Lambda_absorbs_X0"] = (CF.SETS["BAO"], Rl)
    for name, (rows, Rx) in sets.items():
        fit = CF.Fit(rows, Rx)
        out["cmb_fondo"][name] = {}
        for treat in ("profiled", "marginalised"):
            if treat not in old_cf[name]:
                continue
            fun = fit.ll_b if treat == "profiled" else fit.ll_c
            fit.kern = (lambda s, t, D=fit: D.cov + s**2 * D.R.M(t, D.scale)[np.ix_(D.rows, D.rows)]) if False else None
            out["cmb_fondo"][name][treat] = {}
            for th in list(CF.THETAS) + ["free"]:
                key = f"theta_{th}" if th != "free" else "theta_free"
                old = old_cf[name][treat][key]
                g = (lambda s, f=fun: CF.V.P13.max_theta(f, s)[0]) if th == "free" else (lambda s, f=fun, t=th: f(s, t))
                rec = entry(g, old["192"], old["135"])
                out["cmb_fondo"][name][treat][key] = rec
                print(f"cmb_fondo {name:26s} {treat:12s} {key:11s} max@{rec['sigma_at_max_grid']:.3g} ({rec['lnL_max_minus_sigma0']:+.2f}) "
                      f"hist={old['192']} clip={rec['192_clipped']} glob={rec['192_global']} repro={rec['192_reproduced']} dpct={rec['192_change_pct']}  [{time.time() - t0:.0f} s]", flush=True)
    # EXACTO_SN linear bounds: replicate exacto_sn.linear_bounds with the extended-max bound
    old_sn = json.load(open(HERE / "exacto_sn.json"))["linear"]
    R = CF.Response(lambda_absorbs_x0=True)
    u, z = R.u, R.z
    om, orad = CF.OM_FID, S.OMR / CF.H_FID**2
    E2 = om * (1 + z) ** 3 + orad * (1 + z) ** 4 + (1 - om - orad)
    Es = np.sqrt(E2)
    base = R.du * (1 + z)
    N = len(u)
    sn = S.SN["all"]

    def wM_at(k):
        cw = np.zeros(N)
        cw[: k + 1] = base[: k + 1]
        cw[0] *= 0.5
        cw[k] *= 0.5
        dm = float(np.sum(cw / Es))
        return -cw / (2 * Es**3) / dm

    Wsn = np.zeros((sn["n"], N))
    for i, uu in enumerate(sn["uhd"]):
        k = int(np.searchsorted(u, uu)) - 1
        t = (uu - u[k]) / (u[k + 1] - u[k])
        row = (1 - t) * wM_at(k) + t * wM_at(k + 1)
        row[0] -= row.sum()
        Wsn[i] = S.MAG * row
    W = np.vstack([R.W[:15], Wsn])
    fit15 = CF.Fit(CF.SETS["BAO+CMB"], R)
    scale = np.concatenate([fit15.scale[:15], np.ones(sn["n"])])
    Ws = W * scale[:, None]
    y = np.concatenate([S.X.Y15, sn["mb"]])
    Cm = np.zeros((15 + sn["n"], 15 + sn["n"]))
    Cm[:15, :15] = CF.CMEAS15
    Cm[15:, 15:] = np.linalg.inv(sn["Ci"])
    rows_of = {"SN": list(range(15, 15 + sn["n"])), "BAO+SN": list(range(13)) + list(range(15, 15 + sn["n"])), "BAO+CMB+SN": list(range(15 + sn["n"]))}
    G0 = S.G0
    dX0 = np.zeros(G0.N)

    def model_full(h, wm, M):
        hist, _ = S.history(G0, wm / h**2, S.OMR / h**2, dX0)
        return np.concatenate([S.vec15(G0, h, wm, *hist), S.mu_sn(G0, sn, hist[0]) + M])

    for combo, rows in rows_of.items():
        rows = np.array(rows)
        if combo == "SN":
            c0, om_b = S.fit_sn_only(G0, dX0, sn)
            h_b, wm_b = S.H_SN, om_b * S.H_SN**2
            hist, _ = S.history(G0, om_b, S.OMR / S.H_SN**2, dX0)
            rr = sn["mb"] - S.mu_sn(G0, sn, hist[0])
            M_b = float((sn["Ci1"] @ rr) / sn["e"])
            pars = ["Omega_m", "M"]
        else:
            c0, h_b, wm_b, M_b = S.fit_combo(G0, dX0, sn, combo == "BAO+CMB+SN")
            pars = ["h", "omega_m", "M"]
        m0 = model_full(h_b, wm_b, M_b)
        T = []
        for pname in pars:
            eps = {"h": 1e-4, "omega_m": 1e-5, "Omega_m": 1e-4, "M": 1e-3}[pname]

            def mp(sign):
                if pname == "h":
                    return model_full(h_b + sign * eps, wm_b, M_b)
                if pname == "omega_m":
                    return model_full(h_b, wm_b + sign * eps, M_b)
                if pname == "Omega_m":
                    return model_full(h_b, (om_b + sign * eps) * S.H_SN**2, M_b)
                return model_full(h_b, wm_b, M_b + sign * eps)
            T.append((mp(+1) - mp(-1)) / (2 * eps))
        T = np.array(T).T[rows]
        r0 = (y - m0)[rows]
        C = Cm[np.ix_(rows, rows)]
        Wr = Ws[rows]
        out["exacto_sn"][combo] = {}
        for theta in S.THETAS:
            Kw = R.Kw(theta, Wr)
            Mth = Wr @ Kw.T

            def ll_c(s, Mth=Mth):
                L = np.linalg.cholesky(C + s * s * Mth)
                solve = lambda b: np.linalg.solve(L.T, np.linalg.solve(L, b))
                cr, cT = solve(r0), solve(T)
                F = T.T @ cT
                b = T.T @ cr
                return float(-0.5 * (r0 @ cr - b @ np.linalg.solve(F, b) + 2 * np.sum(np.log(np.diag(L))) + np.linalg.slogdet(F)[1]))

            old = old_sn[combo][f"theta_{theta:g}"]
            rec = entry(ll_c, old["192"], old["135"])
            out["exacto_sn"][combo][f"theta_{theta:g}"] = rec
            print(f"exacto_sn {combo:11s} theta={theta:<4g} max@{rec['sigma_at_max_grid']:.3g} ({rec['lnL_max_minus_sigma0']:+.2f}) hist={old['192']} clip={rec['192_clipped']} "
                  f"glob={rec['192_global']} repro={rec['192_reproduced']} dpct={rec['192_change_pct']}  [{time.time() - t0:.0f} s]", flush=True)
    out["runtime_s"] = time.time() - t0
    (HERE / "maximos_globales.json").write_text(json.dumps(out, indent=1), encoding="utf-8")


if __name__ == "__main__":
    main()
