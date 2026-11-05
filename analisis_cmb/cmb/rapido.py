#!/usr/bin/env python3
"""
Exploration only (RAPIDO.md): the fast tremor (theta = 3, 10) in Planck 2018 and in the acoustic peaks.

Set-up of section 3.5 / realizaciones.py and of picos.py: OU path on the fine grid, Gaussian smoothing Delta, 2000-point w(a)
table in DarkEnergyPPF, Takahashi halofit, H0 solved so that 100 theta_* equals the LCDM value, the first 60 standard-normal
realisations of seed 20261003. Cases: theta in {3, 10} x sigma_X in {1.290e-2, 2.204e-2} (bounds at theta = 10 with the
background fixed and marginalised, FONDO.md / perfil_fondo.json).
Smoothing check: a Gaussian smoothing of std Delta keeps a fraction e^{(theta Delta)^2} erfc(theta Delta) of the variance of a
stationary OU with rate theta. Delta = 0.05 keeps 85% at theta = 3 and 62% at theta = 10, so theta = 10 is run with
Delta = 0.0125 (87%, the same theta Delta as theta = 3 with 0.05); the numerical fraction on the grid and CAMB's rho_DE against
ours are reported. Sensitivity: theta = 10 also at Delta = 0.05 and 0.025, theta = 3 also at Delta = 0.025. Every case with Takahashi halofit
(the choice of section 3.5) and with linear lensing, because Takahashi's formula depends explicitly on w(a) and is not valid for
the fast tremor (see CASES).
Per curve: Delta chi^2 against Planck 2018 TT (2 <= l < 30) and phiphi bandpowers (L = 8-400) relative to LCDM (estadisticos.chi2,
asymmetric errors), the l / bin of the largest |Delta D_l / sigma_l| (Delta = curve - LCDM, sigma = Planck error), and the
theory distance of picos.py in TT/TE/EE l in [30, 2500] (lensed and unlensed) with its l of maximum. Also the amplitude of
w_eff(a) = -1 - X'/(3 (Omega_L + X)): p95 of |w + 1| over the table grid per realisation.

    OMP_NUM_THREADS=2 ../venv/bin/python rapido.py     # writes cmb/rapido.{json,txt}
"""

from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "2")

import json
import multiprocessing as mp
import sys
import time
from pathlib import Path

import numpy as np
from scipy.special import erfc

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
import picos as PI  # noqa: E402  (distance(), N_USE, LMAX_USE, RES_C2)

N_USE = PI.N_USE
_PF = json.loads((HERE.parent / "fondo" / "perfil_fondo.json").read_text())["treatments"]
SIG_FIX, SIG_MARG = _PF["a_fixed"]["theta_10"]["sigma_95"], _PF["c_marginalised"]["theta_10"]["sigma_95"]
DELTA_MAIN = {3.0: 0.05, 10.0: 0.0125}
# (theta, Delta, halofit, main). Takahashi is the halofit of section 3.5 / picos.py; its fitting formula depends explicitly on
# w(a) and, with the |w + 1| ~ 0.4 swings of theta = 10 at Delta = 0.0125, inflates C_L^phiphi by ~1.5% while the Weyl potential
# and P(k, z) are unchanged (< 5e-4 at every z): a fit-formula artifact (HMcode fails outright). Every case is therefore run
# also with linear lensing ("lineal"), where the lensing response is converged in Delta.
CASES = [(3.0, 0.05, "takahashi", True), (10.0, 0.0125, "takahashi", True), (3.0, 0.05, "lineal", True), (10.0, 0.0125, "lineal", True),
         (10.0, 0.05, "takahashi", False), (10.0, 0.025, "takahashi", False), (3.0, 0.025, "takahashi", False),
         (10.0, 0.05, "lineal", False), (10.0, 0.025, "lineal", False), (3.0, 0.025, "lineal", False)]

_TARGET = None


def _init(target):
    global _TARGET
    _TARGET = target


def curve(task):
    """task = ('cpl', H0) or ('table', a, w, None). Returns H0, theta_*, TT D_l (2..29), binned pp, lensed/unlensed raw C_l."""
    import camb
    import ou_cmb as M
    import planck_datos as P

    hf = task[-1]
    if task[0] == "cpl":
        make = lambda h: M.make_params(w0=-1.0, wa=0.0, H0=h, halofit=hf)
        h0 = task[1]
    else:
        make = lambda h: M.make_params(a=task[1], w=task[2], H0=h, halofit=hf)
        h0 = M.solve_H0_for_thetastar(_TARGET, make)[0]
    pars = make(h0)
    res = camb.get_results(pars)
    cl = res.get_cmb_power_spectra(pars, CMB_unit="muK", raw_cl=True)
    dl = res.get_cmb_power_spectra(pars, CMB_unit="muK", spectra=["total"])["total"][:, 0]
    pp = res.get_lens_potential_cls(lmax=M.LMAX)[:, 0]
    return {"H0": float(h0), "thetastar": float(res.get_derived_params()["thetastar"]), "TT": dl[2:30],
            "pp": P.bin_pp(pp, P.load_lensing()), "lensed": cl["total"][: PI.LMAX_USE + 1, :4],
            "unlensed": cl["unlensed_total"][: PI.LMAX_USE + 1, :4]}


def tables(theta, sigma, delta):
    import ou_cmb as M

    lna_f, tab_idx, main = M.production_grids()
    h = lna_f[1] - lna_f[0]
    z = np.random.default_rng(M.SEED).standard_normal((200, len(lna_f)))[:N_USE]
    x = sigma * M.ou_uniform(z, h, theta, 1.0)
    xs = M.smooth_uniform(x, h, delta)
    rho, w, bad = M.rho_w_table(lna_f[tab_idx], main, xs[:, tab_idx], M.omega_l0())
    var_kept = float(np.mean(np.var(xs[:, tab_idx][:, main], axis=1) / np.var(x[:, tab_idx][:, main], axis=1)))
    return np.exp(lna_f[tab_idx][main]), w, bad, var_kept, rho


def rho_check(a, w, rho):
    import camb
    import ou_cmb as M

    rc = camb.get_background(M.make_params(a=a, w=w, halofit="takahashi")).get_dark_energy_rho_w(a)[0]
    return float(np.max(np.abs(rc / rc[-1] / rho - 1)))


def pct(v):
    v = np.asarray(v, float)
    return {"p50": float(np.median(v)), "p95": float(np.percentile(v, 95)), "max": float(v.max()), "min": float(v.min())}


def main() -> None:
    import estadisticos as E
    import ou_cmb as M

    t0 = time.time()
    target = M.thetastar_background(M.make_params(halofit="takahashi"))
    h0_lcdm = json.loads((PI.RES_C2 / "meta.json").read_text())["references"]["LCDM"]["H0"]
    tg = E.targets()
    out = {"n_realisations": N_USE, "realisations": "indices 0-59 of the 200 of section 3.5 (seed 20261003)",
           "sigmas": {"fixed_theta10": SIG_FIX, "marginalised_theta10": SIG_MARG}, "delta_main": DELTA_MAIN,
           "smoothing_variance_kept_analytic": {f"theta_{th:g}": {f"Delta_{D:g}": float(np.exp((th * D) ** 2) * erfc(th * D))
                                                                  for D in (0.05, 0.025, 0.0125)} for th in (1.0, 3.0, 10.0)},
           "cases": {}}
    with mp.get_context("spawn").Pool(8, initializer=_init, initargs=(target,)) as pool:
        refs, chi_refs = {}, {}
        for hf in ("takahashi", "lineal"):
            refs[hf] = pool.map(curve, [("cpl", h0_lcdm, hf)])[0]
            chi_refs[hf] = {"TT": float(E.chi2(refs[hf]["TT"], *tg["TT"])), "pp": float(E.chi2(refs[hf]["pp"], *tg["pp"]))}
        out["reference"] = {hf: {"H0": refs[hf]["H0"], "thetastar": refs[hf]["thetastar"], "chi2_TT": chi_refs[hf]["TT"],
                                 "chi2_pp": chi_refs[hf]["pp"]} for hf in refs}
        sig_tt = 0.5 * (tg["TT"][1] + tg["TT"][2])
        sig_pp = tg["pp"][1]
        for theta, delta, hf, is_main in CASES:
            ref, chi_ref = refs[hf], chi_refs[hf]
            for sigma in (SIG_FIX, SIG_MARG):
                a, w, bad, var_kept, rho = tables(theta, sigma, delta)
                assert not bad.any()
                res = pool.map(curve, [("table", a, w[i], hf) for i in range(N_USE)], chunksize=2)
                d_tt = np.array([E.chi2(r["TT"], *tg["TT"]) for r in res]) - chi_ref["TT"]
                d_pp = np.array([E.chi2(r["pp"], *tg["pp"]) for r in res]) - chi_ref["pp"]
                z_tt = np.array([(r["TT"] - ref["TT"]) / sig_tt for r in res])  # 60 x 28, l = 2..29
                z_pp = np.array([(r["pp"] - ref["pp"]) / sig_pp for r in res])  # 60 x 9
                rows = [PI.distance(ref["lensed"], r["lensed"]) for r in res]
                rows_un = [PI.distance(ref["unlensed"], r["unlensed"]) for r in res]
                wamp = np.percentile(np.abs(w + 1), 95, axis=1)
                name = f"theta{theta:g}_Delta{delta:g}_{hf}_sigma{sigma:.3e}"
                c = {"theta": theta, "Delta": delta, "halofit": hf, "sigma_X": sigma, "main": is_main,
                     "variance_kept_numerical": var_kept, "rho_camb_vs_ours_max_rel_r0": rho_check(a, w[0], rho[0]),
                     "w_p95_abs_wp1": pct(wamp), "w_std_wp1_median": float(np.median(np.std(w + 1, axis=1))),
                     "H0": pct([r["H0"] for r in res]), "H0_std": float(np.std([r["H0"] for r in res])),
                     "thetastar_max_abs_diff": float(np.max(np.abs([r["thetastar"] - target for r in res]))),
                     "dchi2_TT": pct(d_tt), "dchi2_pp": pct(d_pp), "dchi2_TT_plus_pp": pct(d_tt + d_pp),
                     "frac_TT_plus_pp_gt_0.5": float(np.mean(d_tt + d_pp > 0.5)),
                     "max_abs_z_TT": pct(np.max(np.abs(z_tt), axis=1)), "max_abs_z_pp": pct(np.max(np.abs(z_pp), axis=1)),
                     "l_of_max_z_TT_p16_p50_p84": [float(x) for x in np.percentile(2 + np.argmax(np.abs(z_tt), axis=1), [16, 50, 84])],
                     "Lav_of_max_z_pp_p16_p50_p84": [float(x) for x in np.percentile(np.array([P_LAV[i] for i in np.argmax(np.abs(z_pp), axis=1)]), [16, 50, 84])],
                     "frac_max_z_in_pp": float(np.mean(np.max(np.abs(z_pp), axis=1) > np.max(np.abs(z_tt), axis=1))),
                     "frac_dchi2_pp_gt_TT": float(np.mean(np.abs(d_pp) > np.abs(d_tt))),
                     "picos_lensed": PI.summarise(rows), "picos_unlensed": PI.summarise(rows_un),
                     "per_realisation": {"dchi2_TT": d_tt.tolist(), "dchi2_pp": d_pp.tolist(), "picos_joint": [r["joint"] for r in rows],
                                         "w_p95": wamp.tolist(), "H0": [r["H0"] for r in res]}}
                out["cases"][name] = c
                print(name, f"var_kept={var_kept:.3f} w95={c['w_p95_abs_wp1']['p50']:.3f} TT p95={c['dchi2_TT']['p95']:.4f} "
                      f"pp p95={c['dchi2_pp']['p95']:.4f} TT+pp p95={c['dchi2_TT_plus_pp']['p95']:.4f} picos p95={c['picos_lensed']['joint']['p95']:.2e} "
                      f"H0 std={c['H0_std']:.3f} [{time.time() - t0:.0f} s]", flush=True)
                (HERE / "rapido.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    out["seconds"] = time.time() - t0
    (HERE / "rapido.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    lines = [f"60 realizaciones; sigma_X = {SIG_FIX:.4e} (fijo, theta=10) y {SIG_MARG:.4e} (marginalizado, theta=10); "
             f"chi2 ref (takahashi) TT={chi_refs['takahashi']['TT']:.3f} pp={chi_refs['takahashi']['pp']:.3f}; "
             f"(lineal) TT={chi_refs['lineal']['TT']:.3f} pp={chi_refs['lineal']['pp']:.3f}",
             f"{'caso':42s} {'varK':>5s} {'w95':>6s} {'TT p50':>8s} {'TT p95':>8s} {'pp p50':>8s} {'pp p95':>8s} {'TT+pp p95':>9s} "
             f"{'picos p50':>9s} {'picos p95':>9s} {'sinlente p95':>12s} {'l max TT':>8s} {'frac pp>TT':>10s} {'H0 std':>6s}"]
    for name, c in out["cases"].items():
        lines.append(f"{name:42s} {c['variance_kept_numerical']:5.3f} {c['w_p95_abs_wp1']['p50']:6.3f} {c['dchi2_TT']['p50']:8.4f} "
                     f"{c['dchi2_TT']['p95']:8.4f} {c['dchi2_pp']['p50']:8.4f} {c['dchi2_pp']['p95']:8.4f} {c['dchi2_TT_plus_pp']['p95']:9.4f} "
                     f"{c['picos_lensed']['joint']['p50']:9.2e} {c['picos_lensed']['joint']['p95']:9.2e} {c['picos_unlensed']['joint']['p95']:12.2e} "
                     f"{c['l_of_max_z_TT_p16_p50_p84'][1]:8.0f} {c['frac_dchi2_pp_gt_TT']:10.2f} {c['H0_std']:6.3f}")
    (HERE / "rapido.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


import planck_datos as _P  # noqa: E402

P_LAV = _P.load_lensing()["L_av"]

if __name__ == "__main__":
    main()
