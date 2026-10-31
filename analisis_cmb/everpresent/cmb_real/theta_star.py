#!/usr/bin/env python3
"""
100 theta_* = 100 r_s(z_*) / D_M(z_*) for every survivor of EVERPRESENT_vs_DESI.md 3a/3d/3e
(supervivientes.json), computed with each realisation's own expansion history (option d1 of
EXPLORACION.md). Run with analisis_cmb/venv (imports prueba_camb, which imports camb):

    ../../venv/bin/python theta_star.py

History:
  Model 2  Omega_hat on the 458 stored points (bao_sim.py, rd_real.py; same seeds) plus N_DRAW draws of the
           unsampled stretch z in (3, 1000), conditioned on the stored points (prueba_camb.gap_sampler,
           275 points, step 0.02 in lambda); seed np.random.default_rng([20261031, index]). Omega_hat is
           interpolated by a cubic spline in lambda = ln a; Omega_de = tanh(Omega_hat); Omega_de = 0 for z > 1e5.
           H^2 = 100^2 [w_m (1+z)^3 + w_r (1+z)^4] / (1 - Omega_de).
  Model 1  full u(a) of historias_m1/ (resim_m1.py, identical to bao_sim.py), Lambda held fixed within each
           a-step as in run_das; omega_m rescaling of Das I Sec. 4.4: H^2 = 100^2 [w_m a^-3 + w_r a^-4 + u w_m / w_m0]
           (radiation not rescaled, as in rd_real.py); Lambda = 0 for z > 9999 (S_0 = 0 at the start).
omega_m: the one with which the realisation passed (optimum of 3d / 3e(b) / 3e(c)*, Planck for 3a / 3e(a));
H0 = 100 sqrt((w_m + w_r)/(1 - Omega_de(0))) (Model 2) or from the run (Model 1) -- theta_* does not depend
on it beyond H(z).
r_s(z_*) = NORM [int_{z_*}^{1e5} c_s/H dz + standard tail z > 1e5] (Model 1: the run up to z = 9999 + standard
tail), with c_s, NORM and the tail of rdcalc.py (the normalisation of rd_real.py).
D_M(z_*) = int_0^{z_*} c/H dz. z_* = CAMB value for Planck 2018 (rd_planck.json -> zstar), the same for every
realisation; the sensitivity to scaling it with omega_m (Hu & Sugiyama 1996) is recorded.
Planck 2018 (arXiv:1807.06209, Table 2, TT,TE,EE+lowE+lensing): 100 theta_* = 1.04110 +- 0.00031.

Validations (written to the json):
  1. LCDM-Planck history through the same Model 2 code: 100 theta_* and r_s(z_*) against CAMB (rd_planck.json);
     the Model 1 code with u = 0 against the plain integral; Model 1 D_M at the DESI redshifts against the run's
     conformal time (bao_m1 npz, dmh).
  2. Survivor 5492 with the gap draw of prueba_camb.py (100theta_* = 0.6845 from CAMB): our value, and CAMB's
     z_*, r_s and D_M for that cosmology.

Writes cmb_real/theta_star.{json,txt}.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline

import prueba_camb as P

HERE = Path(__file__).resolve().parent
EV = HERE.parent
sys.path.insert(0, str(EV))

import bao_sim as B  # noqa: E402
import rdcalc as RC  # noqa: E402

RDJ = json.loads((EV / "resultados" / "rd_planck.json").read_text())
ZSTAR = RDJ["zstar"]
TH_PLANCK, SIG_PLANCK = 1.04110, 0.00031  # arXiv:1807.06209 Table 2, TT,TE,EE+lowE+lensing
N_DRAW = 100
GAP_SEED = 20261031
OM0, ORAD, C = B.OM, B.ORAD, B.C
LAM_GAP = np.arange(-np.log(1000.0) + 0.02, -np.log(4.0) - 0.01, 0.02)  # as prueba_camb.py
N_RS, N_DM = 6000, 8000
CAL = 1.0
LAM_EXT = np.arange(-np.log(1e8), -np.log(1e5) - 0.05, 0.1)  # sensitivity: Model 2 process continued to z = 1e8
# 100 theta_* of the LCDM-Planck history through this code over CAMB's (validation 1): the "calibrated" variant
# divides by it


def zstar_hs(wm: float, wb: float = RC.OMB) -> float:
    """Hu & Sugiyama (1996) fitting formula for z_*, used only as a ratio to scale the CAMB value."""
    g1 = 0.0783 * wb**-0.238 / (1 + 39.5 * wb**0.763)
    g2 = 0.560 / (1 + 21.1 * wb**1.81)
    return 1048 * (1 + 0.00124 * wb**-0.738) * (1 + g1 * wm**g2)


def grids(zs: float, zhi: float = 1e5):
    lam_s = -np.log1p(zs)
    n_rs = int(N_RS * np.log1p(zhi / zs) / np.log1p(1e5 / zs))
    return np.linspace(-np.log1p(zhi), lam_s, n_rs), np.linspace(lam_s, 0.0, N_DM)


def theta_from_ode(ode_rs: np.ndarray, ode_dm: np.ndarray, wm: float, zs: float, zhi: float = 1e5) -> dict:
    """ode_*: Omega_de on grids(zs, zhi) (shape (N,) or (N, n)). Returns 100 theta_*, r_s, D_M (arrays over draws).
    Standard (no dark energy) tail of r_s beyond zhi."""
    l_rs, l_dm = grids(zs, zhi)
    if ode_rs.ndim == 2:
        l_rs, l_dm = l_rs[:, None], l_dm[:, None]
    H_rs = 100.0 * np.sqrt((wm * np.exp(-3 * l_rs) + ORAD * np.exp(-4 * l_rs)) / (1.0 - ode_rs))
    H_dm = 100.0 * np.sqrt((wm * np.exp(-3 * l_dm) + ORAD * np.exp(-4 * l_dm)) / (1.0 - ode_dm))
    rs = RC.NORM * (np.trapezoid(RC.cs(np.exp(-l_rs) - 1.0) * np.exp(-l_rs) / H_rs, l_rs, axis=0)
                    + RC.rd_standard(wm, zhi))
    dm = np.trapezoid(C * np.exp(-l_dm) / H_dm, l_dm, axis=0)
    return {"th": 100.0 * rs / dm, "rs": rs, "dm": dm}


def m2_history(lam_obs, y_obs, y_gap: np.ndarray, zs: float, y_ext: np.ndarray | None = None, zhi: float = 1e5):
    """Omega_de on grids(zs, zhi) from the stored points and gap draws y_gap (275, n), cubic spline in lambda.
    y_ext: draws on LAM_EXT (z in (1e5, 1e8]) for the sensitivity variant."""
    lam = np.concatenate([lam_obs, LAM_GAP] + ([LAM_EXT] if y_ext is not None else []))
    o = np.argsort(lam)
    lam = lam[o]
    keep = np.concatenate([[True], np.diff(lam) > 1e-9])
    Y = np.vstack([np.repeat(y_obs[:, None], y_gap.shape[1], axis=1), y_gap] + ([y_ext] if y_ext is not None else []))
    Y = Y[o][keep]
    cs = CubicSpline(lam[keep], Y, axis=0)
    l_rs, l_dm = grids(zs, zhi)
    return np.tanh(cs(l_rs)), np.tanh(cs(l_dm))


def m1_theta(a: np.ndarray, u: np.ndarray, wm: float, zs: float) -> dict:
    """Model 1: Simpson per a-step with the step's Lambda (u_step[k] = u[k], u[0] := 0), split at a_*."""
    ust = np.where(np.isfinite(u), u, 0.0)[:-1] * wm / OM0
    a0, a1 = a[:-1], a[1:]
    a_s = 1.0 / (1.0 + zs)

    def simpson(lo, hi, g):
        tot = np.zeros_like(lo)
        for aa, wgt in ((lo, 1.0), (0.5 * (lo + hi), 4.0), (hi, 1.0)):
            H2 = wm / aa**3 + ORAD / aa**4 + ust
            H = 100.0 * np.sqrt(np.where(H2 > 0, H2, np.nan))
            tot += wgt * g(aa) / (aa * aa * H)
        return np.nansum(np.where(hi > lo, tot * (hi - lo) / 6.0, 0.0)), bool(np.any((H2 <= 0) & (hi > lo)))

    rs_part, bad1 = simpson(a0, np.minimum(a1, a_s), lambda aa: RC.cs(1.0 / aa - 1.0))
    dm, bad2 = simpson(np.maximum(a0, a_s), a1, lambda aa: C + 0.0 * aa)
    rs = RC.NORM * (rs_part + RC.rd_standard(wm, 1.0 / a[0] - 1.0))
    th = np.nan if (bad1 or bad2) else 100.0 * rs / dm
    return {"th": np.array([th]), "rs": np.array([rs]), "dm": np.array([dm]), "H2_negative": bad1 or bad2}


def wm_of(passes: dict) -> dict:
    out = {}
    for k, v in passes.items():
        out[k] = v.get("omega_m", OM0) if k not in ("3a_fixed", "3c_fixed", "3e_a_rdreal_wm_fixed") else OM0
    return out


def summarise(th: np.ndarray) -> dict:
    d = np.abs(th - TH_PLANCK) / SIG_PLANCK
    f = np.isfinite(th)
    return {"n": int(len(th)), "n_finite": int(f.sum()),
            "th_min": float(np.nanmin(th)) if f.any() else None, "th_med": float(np.nanmedian(th)) if f.any() else None,
            "th_max": float(np.nanmax(th)) if f.any() else None,
            "nsigma_min": float(np.nanmin(d)) if f.any() else None,
            "frac_5sigma": float(np.mean(f & (d < 5))), "frac_10sigma": float(np.mean(f & (d < 10))),
            "draws_10sigma": [int(i) for i in np.nonzero(f & (d < 10))[0]]}


def validations(lam_all, O) -> dict:
    v = {}
    # 1a. LCDM-Planck history through the Model 2 code
    h2 = RDJ["H0"] ** 2 / 1e4
    l_rs, l_dm = grids(ZSTAR)

    def ode_l(l):
        z = np.exp(-l) - 1.0
        de = h2 - OM0 - ORAD
        return de / (OM0 * (1 + z) ** 3 + ORAD * (1 + z) ** 4 + de)

    r = theta_from_ode(ode_l(l_rs), ode_l(l_dm), OM0, ZSTAR)
    v["lcdm_model2_code"] = {"th100": float(r["th"]), "th100_camb": RDJ["thetastar100"],
                             "rel_th_vs_camb": float(r["th"] / RDJ["thetastar100"] - 1),
                             "th_minus_camb_in_planck_sigma": float((r["th"] - RDJ["thetastar100"]) / SIG_PLANCK),
                             "rs": float(r["rs"]), "rs_camb": RDJ["rs_zstar_Mpc"], "rel_rs": float(r["rs"] / RDJ["rs_zstar_Mpc"] - 1),
                             "dm": float(r["dm"]), "dm_camb": RDJ["DMstar_Mpc"], "rel_dm": float(r["dm"] / RDJ["DMstar_Mpc"] - 1),
                             "note": "omega_m incl. massive nu treated as matter, omega_r = photons + 3.046 massless nu; "
                                     "dark energy of LCDM included at z > z_* in r_s up to 1e5"}
    # 1b. Model 1 code with u = 0 against the plain integrals
    a = np.arange(10, 100001) / 1e5
    r1 = m1_theta(a, np.zeros_like(a), OM0, ZSTAR)
    rs_plain = RC.NORM * RC.rd_standard(OM0, ZSTAR)
    from scipy.integrate import quad
    dm_plain = quad(lambda z: C / (100 * np.sqrt(OM0 * (1 + z) ** 3 + ORAD * (1 + z) ** 4)), 0, ZSTAR, epsrel=1e-11, limit=200)[0]
    v["model1_code_no_de"] = {"rel_rs": float(r1["rs"][0] / rs_plain - 1), "rel_dm": float(r1["dm"][0] / dm_plain - 1)}
    # 1c. Model 1 D_M at the DESI redshifts against the run's own conformal time (bao_m1 npz, dmh = D_M H0 / c)
    R = np.load(EV / "resultados" / "bao_m1_alpha0.01.npz")
    h = np.load(HERE / "historias_m1" / "alpha0.01_row5367.npz")
    zu = R["z_unique"]
    errs = []
    for zz, dmh in zip(zu, R["dmh"][5367]):
        rr = m1_theta(h["a"], h["u"], OM0, float(zz))
        errs.append(float(rr["dm"][0] * R["H0"][5367] / C / dmh - 1))
    v["model1_dm_vs_run_row5367"] = {"z": [float(x) for x in zu], "rel_diff": errs}
    # 2. survivor 5492 with the gap draw of prueba_camb.py, against CAMB
    j = 5492
    sv = {r_["index"]: r_ for r_ in json.loads((HERE / "supervivientes.json").read_text())["model2"]}
    wm = sv[j]["passes"]["3d_wm_free"]["omega_m"]
    mean, Ls = P.gap_sampler(lam_all, O[:, j], LAM_GAP)
    draws = mean[:, None] + Ls @ np.random.default_rng(GAP_SEED).standard_normal((len(LAM_GAP), P.N_GAP))
    pick = int(np.nonzero(np.all(draws > 0, axis=0))[0][0])
    ours = theta_from_ode(*m2_history(lam_all, O[:, j], draws[:, [pick]], ZSTAR), wm, ZSTAR)
    a_t, ode_t, _, w_t = P.history(lam_all, O[:, j], LAM_GAP, draws[:, pick], wm)
    H0 = 100.0 * np.sqrt((wm + ORAD) / (1.0 - ode_t[-1]))
    pars = P.M.make_params(H0=H0, halofit="takahashi")
    pars.set_cosmology(H0=H0, ombh2=P.OMB, omch2=wm - P.OMB - P.OMNU, mnu=0.06, omk=0.0, tau=P.M.PLANCK18["tau"])
    de = P.DarkEnergyPPF()
    de.set_w_a_table(a_t, w_t)
    pars.DarkEnergy = de
    der = P.camb.get_background(pars).get_derived_params()
    at_camb_z = theta_from_ode(*m2_history(lam_all, O[:, j], draws[:, [pick]], der["zstar"]), wm, der["zstar"])
    v["survivor5492_vs_camb"] = {
        "gap_pick": pick, "omega_m": wm, "H0": float(H0),
        "ours_at_zstar_planck": {"th100": float(ours["th"][0]), "rs": float(ours["rs"][0]), "dm": float(ours["dm"][0])},
        "camb": {"th100": der["thetastar"], "zstar": der["zstar"], "rs": der["rstar"], "dm": der["DAstar"] * 1000.0},
        "ours_at_zstar_camb": {"th100": float(at_camb_z["th"][0]), "rs": float(at_camb_z["rs"][0]), "dm": float(at_camb_z["dm"][0])},
        "note": "CAMB extrapolates the w table beyond z = 1e5 (rho_DE continues); here Omega_de = 0 for z > 1e5"}
    return v


def main() -> None:
    S = json.loads((HERE / "supervivientes.json").read_text())
    lam_all, O = P.stored_points()
    res = {"planck": {"th100": TH_PLANCK, "sigma": SIG_PLANCK, "source": "arXiv:1807.06209 Table 2, TT,TE,EE+lowE+lensing"},
           "zstar": ZSTAR, "n_draw_model2": N_DRAW, "gap_seed": f"default_rng([{GAP_SEED}, index])",
           "validations": validations(lam_all, O), "model2": [], "model1": []}
    print(json.dumps(res["validations"], indent=1), flush=True)
    global CAL
    CAL = res["validations"]["lcdm_model2_code"]["th100"] / res["validations"]["lcdm_model2_code"]["th100_camb"]
    dz_max = 0.0
    for s in S["model2"]:
        j = s["index"]
        mean, Ls = P.gap_sampler(lam_all, O[:, j], LAM_GAP)
        Yg = mean[:, None] + Ls @ np.random.default_rng([GAP_SEED, j]).standard_normal((len(LAM_GAP), N_DRAW))
        me, Le = P.gap_sampler(lam_all, O[:, j], LAM_EXT)
        Ye = me[:, None] + Le @ np.random.default_rng([GAP_SEED, j, 1]).standard_normal((len(LAM_EXT), N_DRAW))
        row = {"index": j, "by_criterion": {}}
        for crit, wm in wm_of(s["passes"]).items():
            r = theta_from_ode(*m2_history(lam_all, O[:, j], Yg, ZSTAR), wm, ZSTAR)
            zs2 = ZSTAR * zstar_hs(wm) / zstar_hs(OM0)
            r2 = theta_from_ode(*m2_history(lam_all, O[:, j], Yg, zs2), wm, zs2)
            r3 = theta_from_ode(*m2_history(lam_all, O[:, j], Yg, ZSTAR, Ye, 1e8), wm, ZSTAR, 1e8)
            dz_max = max(dz_max, float(np.nanmax(np.abs(r2["th"] - r["th"]))))
            row["by_criterion"][crit] = {"omega_m": wm, **summarise(r["th"]),
                                         "variants": {"zstar_hs": summarise(r2["th"]), "calibrated": summarise(r["th"] / CAL),
                                                      "process_to_1e8": summarise(r3["th"]),
                                                      "max_abs_dth_process_to_1e8": float(np.max(np.abs(r3["th"] - r["th"])))}}
        res["model2"].append(row)
    for alpha in ("0.01", "0.005"):
        for s in S[f"model1_alpha{alpha}"]:
            h = np.load(HERE / "historias_m1" / f"alpha{alpha}_row{s['npz_row']}.npz")
            row = {"alpha": float(alpha), "npz_row": s["npz_row"], "chunk": s["chunk"], "idx_in_chunk": s["idx_in_chunk"],
                   "by_criterion": {}}
            for crit, wm in wm_of(s["passes"]).items():
                r = m1_theta(h["a"], h["u"], wm, ZSTAR)
                zs2 = ZSTAR * zstar_hs(wm) / zstar_hs(OM0)
                r2 = m1_theta(h["a"], h["u"], wm, zs2)
                if np.isfinite(r["th"][0]):
                    dz_max = max(dz_max, float(abs(r2["th"][0] - r["th"][0])))
                row["by_criterion"][crit] = {"omega_m": wm, "H2_negative": r["H2_negative"], **summarise(r["th"]),
                                             "variants": {"zstar_hs": summarise(r2["th"]), "calibrated": summarise(r["th"] / CAL)}}
            res["model1"].append(row)
    res["zstar_hs_scaling_max_abs_dth100"] = dz_max

    # per survivor: the criterion that brings it closest to Planck
    def best(row):
        return min(row["by_criterion"].items(), key=lambda kv: kv[1]["nsigma_min"] if kv[1]["nsigma_min"] is not None else np.inf)

    allrows = [("M2", r) for r in res["model2"]] + [("M1", r) for r in res["model1"]]
    tot = {"n_survivors": len(allrows), "n_any_draw_5sigma": 0, "n_any_draw_10sigma": 0, "passing_10sigma": []}
    for m, r in allrows:
        for crit, d in r["by_criterion"].items():
            if d["draws_10sigma"]:
                tot["passing_10sigma"].append({"model": m, "id": r.get("index", r.get("npz_row")), "criterion": crit,
                                               "draws": d["draws_10sigma"]})
        k, d = best(r)
        tot["n_any_draw_5sigma"] += int(any(x["frac_5sigma"] > 0 for x in r["by_criterion"].values()))
        tot["n_any_draw_10sigma"] += int(any(x["frac_10sigma"] > 0 for x in r["by_criterion"].values()))
    allth = np.concatenate([[x for x in [d["th_min"], d["th_max"]] if x is not None] for _, r in allrows
                            for d in r["by_criterion"].values()])
    tot["th100_range_all"] = [float(allth.min()), float(allth.max())]
    nsig = [d["nsigma_min"] for _, r in allrows for d in r["by_criterion"].values() if d["nsigma_min"] is not None]
    tot["closest_nsigma"] = float(min(nsig))
    for var in ("zstar_hs", "calibrated", "process_to_1e8"):
        n5 = n10 = n10c = 0
        for m, r in allrows:
            vs = [d["variants"][var] for d in r["by_criterion"].values() if var in d["variants"]]
            if not vs:
                continue
            n5 += int(any(x["frac_5sigma"] > 0 for x in vs))
            n10 += int(any(x["frac_10sigma"] > 0 for x in vs))
        for m, r in allrows:
            d = r["by_criterion"].get("3e_c_rdreal_wm_2sigma_nowindow")
            if d and var in d["variants"]:
                n10c += int(d["variants"][var]["frac_10sigma"] > 0)
        tot[f"variant_{var}"] = {"n_any_draw_5sigma": n5, "n_any_draw_10sigma": n10, "n_3e_c_any_draw_10sigma": n10c}
    tot["n_3e_c_any_draw_5sigma"] = sum(int(r["by_criterion"].get("3e_c_rdreal_wm_2sigma_nowindow", {"frac_5sigma": 0})["frac_5sigma"] > 0)
                                        for _, r in allrows)
    tot["n_3e_c_any_draw_10sigma"] = sum(int(r["by_criterion"].get("3e_c_rdreal_wm_2sigma_nowindow", {"frac_10sigma": 0})["frac_10sigma"] > 0)
                                         for _, r in allrows)
    # criterion (c) with theta_* in place of the window proxy: expected fraction of the 20000 Model 2 realisations
    # (each realisation's unsampled stretch is random, so a realisation passes with probability = its fraction of draws)
    for k in ("5", "10"):
        tot[f"model2_3e_c_expected_fraction_{k}sigma"] = sum(
            r["by_criterion"].get("3e_c_rdreal_wm_2sigma_nowindow", {f"frac_{k}sigma": 0.0})[f"frac_{k}sigma"]
            for r in res["model2"]) / 20000
    sv = {x["index"]: x for x in S["model2"]}
    near = sorted({r["index"] for r in res["model2"] if any(d["frac_10sigma"] > 0 for d in r["by_criterion"].values())})
    tot["model2_10sigma_with_window_proxy"] = [i for i in near if sv[i]["window_proxy"]]
    tot["n_model1_any_draw_10sigma"] = sum(int(any(d["frac_10sigma"] > 0 for d in r["by_criterion"].values()))
                                           for m, r in allrows if m == "M1")
    tot["max_abs_dth_process_to_1e8"] = max(d["variants"]["max_abs_dth_process_to_1e8"] for m, r in allrows if m == "M2"
                                            for d in r["by_criterion"].values())
    res["totals"] = tot
    (HERE / "theta_star.json").write_text(json.dumps(res, indent=1), encoding="utf-8")

    lines = [f"100 theta_* por superviviente (Planck 2018: {TH_PLANCK} +- {SIG_PLANCK}); z_* = {ZSTAR:.3f} (CAMB, Planck 2018)",
             f"Modelo 2: {N_DRAW} sorteos del hueco z in (3, 1000) por superviviente; Modelo 1: historia completa (1 valor).",
             "Una fila por (superviviente, criterio con el que pasa). n_sigma_min: distancia minima a Planck en sigmas.",
             "", f"{'modelo':7s} {'id':>7s} {'criterio':32s} {'omega_m':>8s} {'min':>8s} {'mediana':>8s} {'max':>8s}"
             f" {'n_sig_min':>9s} {'f<5s':>5s} {'f<10s':>5s}"]
    for m, r in allrows:
        idn = str(r["index"]) if m == "M2" else f"{r['alpha']:g}/{r['npz_row']}"
        for crit, d in r["by_criterion"].items():
            fmt = lambda x: f"{x:8.4f}" if x is not None else f"{'NaN':>8s}"
            ns = f"{d['nsigma_min']:9.0f}" if d["nsigma_min"] is not None else f"{'H2<0':>9s}"
            lines.append(f"{m:7s} {idn:>7s} {crit:32s} {d['omega_m']:8.4f} {fmt(d['th_min'])} {fmt(d['th_med'])} {fmt(d['th_max'])}"
                         f" {ns} {d['frac_5sigma']:5.2f} {d['frac_10sigma']:5.2f}")
    lines += ["", "Totales: " + json.dumps(tot), f"Escalar z_* con omega_m (Hu-Sugiyama): max |d 100theta_*| = {dz_max:.2e}",
              "", "Validaciones: " + json.dumps(res["validations"], indent=1)]
    (HERE / "theta_star.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines[:8]))
    print(json.dumps(tot, indent=1), f"\nzstar HS max dth {dz_max:.2e}")


if __name__ == "__main__":
    main()
