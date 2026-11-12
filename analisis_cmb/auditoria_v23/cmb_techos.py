#!/usr/bin/env python3
"""
Exploration only (CMB_TECHOS.md): do the paths that survive BAO + theta_* + omega_m at the ceilings of EXACTO.md get
discarded by the full CMB (late ISW and lensing)? Nothing in results/ or papers changes.

Cells (theta, sigma_X) = (1, 0.43), (3, 0.29), (10, 0.19), the exact 5 % ceilings of EXACTO.md (BAO + theta_* + omega_m).
1. The 2000 paths of each cell are regenerated with exacto.py (same grid; seed of the exacto.py cell with the same theta and the
   nearest sigma_X, the paths being linear in sigma_X at fixed normals) and the survivors chi^2_min <= chi^2_obs
   (17.238, 15 rows) are kept.
2. Smoothing as rapido.py: Gaussian of std Delta in ln a (Delta = 0.05 for theta = 1, 3; 0.0125 for theta = 10), kernel
   renormalised on the grid. The exacto grid starts at a = 1, so each path is continued to ln a in (0, 0.4] (u < 0) with the
   same AR(1) recursion (a stationary OU is reversible; fresh normals, seed [SEED, cell, 1]) so that a = 1 is interior to the
   kernel, as the padding of ou_cmb.py. chi^2_min is recomputed with the smoothed path (exacto.fit) and the passers counted.
3. rho_DE/rho_crit,0 = (h^2 - omega_m - omega_r)/h^2 + X_s(z) - X_s(0) with the (h, omega_m) fitted to the smoothed path;
   realisations with rho_DE <= 0 anywhere on the CAMB table (1e-3 <= a <= 1) cannot be run (DarkEnergyPPF integrates rho from
   w) and are counted. This biases the sample towards the tamer paths.
4. The first 40 remaining (by index): CAMB, DarkEnergyPPF with the w(a) table (w = -1 - (1/3) d ln rho / d ln a on the padded
   grid, table on 1e-3 <= a <= 1), linear lensing, H0 = 100 h and omega_m of that realisation, omega_b = 0.02237 (rdcalc.OMB,
   the Planck value of exacto.py), omega_c = omega_m - omega_b - omega_nu (mnu = 0.06 eV), tau, n_s, A_s of Planck 2018.
5. Reference per realisation: LCDM with the same omega_b, omega_c, mnu, and H0 solved so that CAMB's 100 theta_* equals the
   realisation's (ou_cmb.solve_H0_for_thetastar).
6. Statistics, realisation minus reference: Delta chi^2 against Planck 2018 TT (2 <= l < 30) and phiphi bandpowers
   (estadisticos.chi2, planck_datos), and the peak distance TT/TE/EE l in [30, 2500] of picos.distance (joint). Total = sum.
   Per cell: median, p95, fractions with total > 1 and > 4, and where the maximum |Delta/sigma| falls (TT l < 30 bin, phiphi
   bin, or peaks l).
Check by two routes, every run realisation: CAMB's rho_DE(a) (get_dark_energy_rho_w) against the path's rho on the table;
faithful if max_a |rho_DE(CAMB) - rho_DE(path)| / rho_total(a) < 1e-3 (the relative error is also recorded). Faithful and
unfaithful realisations are reported separately; only the faithful ones enter the conclusion.

    OMP_NUM_THREADS=2 ../venv/bin/python cmb_techos.py prep   # steps 1-3 + one trial CAMB run, writes cmb_techos_prep.json
    OMP_NUM_THREADS=2 ../venv/bin/python cmb_techos.py run    # steps 4-6, writes cmb_techos.json
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

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parent / "cmb"))
import exacto as X  # noqa: E402

CELLS = [(1.0, 0.43, 0.05), (3.0, 0.29, 0.05), (10.0, 0.19, 0.0125)]  # theta, sigma_X, Delta
# the ceilings are interpolated crossings, not cells of exacto.py; the paths are linear in sigma_X at fixed normals, so each
# ceiling is generated with the seed of the exacto.py cell with the same theta and the nearest sigma_X (0.5, 0.3, 0.2)
_IDX = {(s, t): i for i, (s, t) in enumerate((s, t) for t in X.THETAS for s in X.SIGMAS)}
CELL_INDEX = {(0.43, 1.0): _IDX[(0.5, 1.0)], (0.29, 3.0): _IDX[(0.3, 3.0)], (0.19, 10.0): _IDX[(0.2, 10.0)]}
N_RUN = 40
A_MIN = 1e-3
PAD = 0.4
OMB = X.E.RC.OMB  # 0.02237
OMNU = X.E.RDJ["rd_aubourg"]["omega_nu"]  # 0.000645 (mnu = 0.06 eV), as rd_planck.py
CHI2_OBS = 17.238  # LCDM, 15 rows (EXACTO.md / ESCALON.md)
LMAX = 2500

# padded grid: u in [-PAD, 0) with the uniform step, then the exacto grid U (u >= 0)
_npad = int(np.ceil(PAD / X.DU))
U_PAD = np.concatenate([-X.DU * np.arange(_npad, 0, -1), X.U])
N_PAD = _npad
MAIN = np.arange(len(U_PAD)) >= N_PAD
TABLE = (U_PAD >= 0) & (U_PAD <= -np.log(A_MIN) + 1e-12)  # 1e-3 <= a <= 1


def smoothing_matrix(u: np.ndarray, delta: float) -> np.ndarray:
    d = u[:, None] - u[None, :]
    w = np.exp(-0.5 * (d / delta) ** 2)
    w[np.abs(d) > 6 * delta] = 0.0
    return w / w.sum(axis=1, keepdims=True)


def extend_backwards(Xm: np.ndarray, theta: float, sigma: float, rng) -> np.ndarray:
    """Continue the paths (rows) from u = 0 to u = -PAD with the AR(1) recursion; returns (n, N_PAD + N)."""
    n = Xm.shape[0]
    rho = np.exp(-theta * X.DU)
    s = sigma * np.sqrt(1 - rho**2)
    back = np.empty((n, N_PAD))
    prev = Xm[:, 0]
    xi = rng.standard_normal((n, N_PAD))
    for k in range(N_PAD):
        prev = rho * prev + s * xi[:, k]
        back[:, N_PAD - 1 - k] = prev
    return np.concatenate([back, Xm], axis=1)


def rho_w(Xs_pad: np.ndarray, h: float, wm: float):
    """rho_DE/rho_crit,0 and w on the padded grid for one smoothed path; (rho, w)."""
    ol = (h * h - wm - X.OMR) / (h * h)
    i0 = N_PAD
    rho = ol + Xs_pad - Xs_pad[i0]
    with np.errstate(invalid="ignore", divide="ignore"):
        w = -1.0 + np.gradient(np.log(np.where(rho > 0, rho, np.nan)), U_PAD) / 3.0  # d ln a = -du
    return rho, w


def prep() -> dict:
    out = {"cells": {}, "chi2_obs": CHI2_OBS, "pad": PAD, "n_pad": int(N_PAD), "A_min": A_MIN}
    X.P0 = (0.681929893729943, 0.14071555057242932)
    for theta, sigma, delta in CELLS:
        t0 = time.time()
        i_cell = CELL_INDEX[(sigma, theta)]
        rng = np.random.default_rng([X.SEED, i_cell])
        P = X.ou_paths(sigma, theta, X.N_REAL, rng)
        dX = P - P[:, [0]]
        e2min = (X.P0[1] * X._1PZ3 + X.OMR * X._1PZ4 + (X.P0[0] ** 2 - X.P0[1] - X.OMR) + X.P0[0] ** 2 * dX).min(axis=1)
        ok = e2min > 0
        chi = np.full(X.N_REAL, np.inf)
        for i in np.nonzero(ok)[0]:
            chi[i] = X.fit(dX[i], 15)[0]
        surv = np.nonzero(chi <= CHI2_OBS)[0]
        rec = {"theta": theta, "sigma_X": sigma, "Delta": delta, "n_ok_H2": int(ok.sum()), "n_survivors": int(len(surv)),
               "P_survive": float(len(surv) / X.N_REAL)}
        # smoothing with backward extension
        rng2 = np.random.default_rng([X.SEED, i_cell, 1])
        Pp = extend_backwards(P[surv], theta, sigma, rng2)
        S = smoothing_matrix(U_PAD, delta)
        Xs = Pp @ S.T
        var_kept = float(np.mean(np.var(Xs[:, MAIN], axis=1) / np.var(Pp[:, MAIN], axis=1)))
        rec["variance_kept_numerical"] = var_kept
        rec["variance_kept_analytic"] = float(np.exp((theta * delta) ** 2) * __import__("scipy.special", fromlist=["erfc"]).erfc(theta * delta))
        dXs = Xs[:, MAIN] - Xs[:, [N_PAD]]
        chi_s, hs, wms = np.empty(len(surv)), np.empty(len(surv)), np.empty(len(surv))
        for j in range(len(surv)):
            chi_s[j], hs[j], wms[j] = X.fit(dXs[j], 15)
        pass_s = chi_s <= CHI2_OBS
        rec["n_pass_smoothed"] = int(pass_s.sum())
        rec["chi2_smoothed_minus_raw_p16_p50_p84"] = [float(v) for v in np.percentile(chi_s - chi[surv], [16, 50, 84])]
        # rho_DE > 0 on the CAMB table
        neg = np.zeros(len(surv), bool)
        wmax = np.full(len(surv), np.nan)
        for j in range(len(surv)):
            rho, w = rho_w(Xs[j], hs[j], wms[j])
            neg[j] = bool(np.any(rho[TABLE] <= 0))
            if not neg[j]:
                wmax[j] = float(np.nanmax(np.abs(w[TABLE] + 1)))
        runnable = pass_s & ~neg
        rec["n_rho_le_0_among_passers"] = int(np.sum(pass_s & neg))
        rec["frac_rho_le_0_among_passers"] = float(np.mean(neg[pass_s])) if pass_s.any() else None
        rec["n_runnable"] = int(runnable.sum())
        rec["abs_wp1_max_p50_p95_among_runnable"] = [float(v) for v in np.percentile(wmax[runnable], [50, 95])] if runnable.any() else None
        idx_run = surv[runnable][:N_RUN]
        rec["indices_run"] = [int(i) for i in idx_run]
        np.savez_compressed(HERE / f"cmb_techos_cell{i_cell}.npz", surv=surv, chi_raw=chi[surv], chi_s=chi_s, h=hs, wm=wms, pass_s=pass_s,
                            neg=neg, Xs=Xs[runnable][:N_RUN], idx_run=idx_run, U_PAD=U_PAD)
        rec["seconds"] = time.time() - t0
        out["cells"][f"theta_{theta:g}"] = rec
        print(json.dumps(rec), flush=True)
    (HERE / "cmb_techos_prep.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    return out


# ----------------------------------------------------------------------------------------------------------------------
# CAMB
# ----------------------------------------------------------------------------------------------------------------------
def make_params(H0: float, wm: float, a=None, w=None):
    import camb
    from camb.dark_energy import DarkEnergyPPF
    import ou_cmb as M

    p = M.PLANCK18
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=H0, ombh2=OMB, omch2=wm - OMB - OMNU, tau=p["tau"], mnu=p["mnu"], omk=0.0)
    pars.InitPower.set_params(As=p["As"], ns=p["ns"])
    pars.set_for_lmax(LMAX, lens_potential_accuracy=1)
    de = DarkEnergyPPF()
    if a is not None:
        de.set_w_a_table(a, w)
    else:
        de.set_params(w=-1.0, wa=0.0)
    pars.DarkEnergy = de
    pars.NonLinear = camb.model.NonLinear_none
    return pars


def spectra(pars):
    import camb
    import planck_datos as P
    import picos as PI

    res = camb.get_results(pars)
    cl = res.get_cmb_power_spectra(pars, CMB_unit="muK", raw_cl=True)
    dl = res.get_cmb_power_spectra(pars, CMB_unit="muK", spectra=["total"])["total"][:, 0]
    pp = res.get_lens_potential_cls(lmax=LMAX)[:, 0]
    return {"TT": dl[2:30], "pp": P.bin_pp(pp, P.load_lensing()), "lensed": cl["total"][: PI.LMAX_USE + 1, :4],
            "thetastar": float(res.get_derived_params()["thetastar"])}


def joint_per_l(ref: np.ndarray, mod: np.ndarray) -> np.ndarray:
    """Per-l joint (TT, EE, TE) Delta chi^2 of picos.distance (same covariance), l in [LMIN, LMAX_USE]."""
    import picos as PI

    nt, npol = PI.noise()
    ell = np.arange(PI.LMAX_USE + 1)
    sl = slice(PI.LMIN, PI.LMAX_USE + 1)
    nu = (2 * ell[sl] + 1) * PI.FSKY
    TT, EE, TE = ref[sl, 0] + nt[sl], ref[sl, 1] + npol[sl], ref[sl, 3]
    d = (mod - ref)[sl]
    cov = np.empty((len(nu), 3, 3))
    cov[:, 0, 0], cov[:, 1, 1], cov[:, 2, 2] = 2 * TT**2 / nu, 2 * EE**2 / nu, (TE**2 + TT * EE) / nu
    cov[:, 0, 1] = cov[:, 1, 0] = 2 * TE**2 / nu
    cov[:, 0, 2] = cov[:, 2, 0] = 2 * TT * TE / nu
    cov[:, 1, 2] = cov[:, 2, 1] = 2 * EE * TE / nu
    v = np.stack([d[:, 0], d[:, 1], d[:, 3]], axis=1)
    return np.einsum("li,li->l", v, np.linalg.solve(cov, v[:, :, None])[:, :, 0])


def one(task):
    """task = (cell key, j, Xs_pad, h, wm). Realisation + its reference; rho_DE check in H^2 units; statistics."""
    import camb
    import ou_cmb as M
    import estadisticos as ES
    import picos as PI

    key, j, Xs, h, wm = task[:5]
    out = {"cell": key, "j": int(j), "fail": None}
    try:
        rho, w = rho_w(Xs, h, wm)
        a = np.exp(-U_PAD[TABLE])[::-1]
        wt = w[TABLE][::-1]
        rp = rho[TABLE][::-1]  # path rho_DE / rho_crit,0, increasing a
        pars = make_params(100 * h, wm, a, wt)
        rc = camb.get_background(pars).get_dark_energy_rho_w(a)[0]
        rc = rc / rc[-1] * rp[-1]  # CAMB rho_DE in units of rho_crit,0 (normalised today)
        om, orad = wm / h**2, X.OMR / h**2
        rtot = om / a**3 + orad / a**4 + rp
        out["rho_err_H2"] = float(np.max(np.abs(rc - rp) / rtot))
        out["rho_err_rel"] = float(np.max(np.abs(rc / rp - 1)))
        out["faithful"] = bool(out["rho_err_H2"] < 1e-3)
        out["abs_wp1_max"] = float(np.nanmax(np.abs(wt + 1)))
        out["rho_min_over_rho0"] = float(rp.min() / rp[-1])
        dxs = Xs[TABLE] - Xs[N_PAD]
        out["path_amp_std"] = float(np.std(dxs))
        out["path_amp_max"] = float(np.max(np.abs(dxs)))
        r = spectra(pars)
        out["thetastar"] = r["thetastar"]
        out["H0"] = 100 * h
        make = lambda H: make_params(H, wm)
        H0_ref, th_ref, _ = M.solve_H0_for_thetastar(r["thetastar"], make)
        ref = spectra(make(H0_ref))
        out["H0_ref"], out["thetastar_ref"] = float(H0_ref), ref["thetastar"]
        tg = ES.targets()
        d_tt = float(ES.chi2(r["TT"], *tg["TT"]) - ES.chi2(ref["TT"], *tg["TT"]))
        d_pp = float(ES.chi2(r["pp"], *tg["pp"]) - ES.chi2(ref["pp"], *tg["pp"]))
        sig_tt = 0.5 * (tg["TT"][1] + tg["TT"][2])
        z_tt = (r["TT"] - ref["TT"]) / sig_tt
        z_pp = (r["pp"] - ref["pp"]) / tg["pp"][1]
        pk = PI.distance(ref["lensed"], r["lensed"])
        jl = joint_per_l(ref["lensed"], r["lensed"])
        out.update(dchi2_TT=d_tt, dchi2_pp=d_pp, picos_joint=pk["joint"], picos_TT=pk["TT"], picos_lmax_joint=pk["lmax_joint"],
                   total=d_tt + d_pp + pk["joint"], max_z_TT=float(np.max(np.abs(z_tt))), l_max_z_TT=int(2 + np.argmax(np.abs(z_tt))),
                   max_z_pp=float(np.max(np.abs(z_pp))), max_z_picos=float(np.sqrt(jl.max())), l_max_z_picos=int(PI.LMIN + np.argmax(jl)))
    except Exception as exc:  # noqa: BLE001
        out["fail"] = f"{type(exc).__name__}: {exc}".replace("\n", " ")[:300]
    return out


def _stats(rs: list[dict]) -> dict:
    if not rs:
        return {"n": 0}
    pct = lambda v: {"p50": float(np.median(v)), "p95": float(np.percentile(v, 95)), "max": float(np.max(v)), "min": float(np.min(v))}
    tot = np.array([r["total"] for r in rs])
    c = {"n": len(rs), "total": pct(tot), "dchi2_TT": pct([r["dchi2_TT"] for r in rs]), "dchi2_pp": pct([r["dchi2_pp"] for r in rs]),
         "picos_joint": pct([r["picos_joint"] for r in rs]), "frac_total_gt_1": float(np.mean(tot > 1)), "frac_total_gt_4": float(np.mean(tot > 4)),
         "abs_wp1_max": pct([r["abs_wp1_max"] for r in rs]), "rho_err_H2": pct([r["rho_err_H2"] for r in rs]),
         "rho_err_rel": pct([r["rho_err_rel"] for r in rs]), "H0": pct([r["H0"] for r in rs]),
         "H0_ref_minus_H0": pct([r["H0_ref"] - r["H0"] for r in rs]), "path_amp_std": pct([r["path_amp_std"] for r in rs])}
    where = []
    for r in rs:
        zs = {"TT_l<30": r["max_z_TT"], "phiphi": r["max_z_pp"], "picos": r["max_z_picos"]}
        where.append(max(zs, key=zs.get))
    c["where_max_abs_z"] = {k: float(np.mean([w == k for w in where])) for k in ("TT_l<30", "phiphi", "picos")}
    c["l_of_max_z_TT_p16_p50_p84"] = [float(v) for v in np.percentile([r["l_max_z_TT"] for r in rs], [16, 50, 84])]
    c["frac_l_max_TT_lt_10"] = float(np.mean([r["l_max_z_TT"] < 10 for r in rs]))
    c["max_z_TT_p50_p95"] = [float(v) for v in np.percentile([r["max_z_TT"] for r in rs], [50, 95])]
    c["max_z_pp_p50_p95"] = [float(v) for v in np.percentile([r["max_z_pp"] for r in rs], [50, 95])]
    c["max_z_picos_p50_p95"] = [float(v) for v in np.percentile([r["max_z_picos"] for r in rs], [50, 95])]
    if len(rs) >= 5:
        from scipy.stats import pearsonr, spearmanr

        w1 = np.array([r["abs_wp1_max"] for r in rs])
        amp = np.array([r["path_amp_std"] for r in rs])
        ampm = np.array([r["path_amp_max"] for r in rs])
        c["corr_total_vs_abs_wp1_max"] = {"spearman": float(spearmanr(tot, w1)[0]), "pearson": float(pearsonr(tot, w1)[0])}
        c["corr_total_vs_path_amp_std"] = {"spearman": float(spearmanr(tot, amp)[0]), "pearson": float(pearsonr(tot, amp)[0])}
        c["corr_total_vs_path_amp_max"] = {"spearman": float(spearmanr(tot, ampm)[0]), "pearson": float(pearsonr(tot, ampm)[0])}
        c["corr_abs_wp1_max_vs_path_amp_std"] = {"spearman": float(spearmanr(w1, amp)[0])}
    return c


def _child(task, q):
    q.put(one(task))


def run_isolated(tasks, n_par: int = 8, timeout: float = 900.0) -> list[dict]:
    """One spawned process per task (a Fortran STOP inside CAMB kills the worker; a Pool would hang), n_par at a time,
    with a time limit; progress to cmb_techos_progress.log."""
    ctx = mp.get_context("spawn")
    log = open(HERE / "cmb_techos_progress.log", "w")
    res = [None] * len(tasks)
    pending = list(range(len(tasks)))
    active = {}
    t0 = time.time()
    while pending or active:
        while pending and len(active) < n_par:
            i = pending.pop(0)
            q = ctx.Queue()
            pr = ctx.Process(target=_child, args=(tasks[i], q))
            pr.start()
            active[i] = (pr, q, time.time())
        time.sleep(0.5)
        for i in list(active):
            pr, q, ts = active[i]
            if not q.empty():
                res[i] = q.get()
                pr.join(5)
            elif not pr.is_alive():
                res[i] = {"cell": tasks[i][0], "j": int(tasks[i][1]), "fail": f"worker died (exit code {pr.exitcode})"}
            elif time.time() - ts > timeout:
                pr.terminate()
                res[i] = {"cell": tasks[i][0], "j": int(tasks[i][1]), "fail": f"timeout {timeout:.0f} s"}
            else:
                continue
            del active[i]
            r = res[i]
            log.write(f"{time.time() - t0:7.0f} s  {r['cell']} j={r['j']}  " + (f"FAIL {r['fail']}" if r.get("fail") else
                      f"faithful={r['faithful']} errH2={r['rho_err_H2']:.2e} total={r['total']:+.3f}") + "\n")
            log.flush()
    log.close()
    return res


def run() -> None:
    import picos as PI  # noqa: F401
    t0 = time.time()
    prep_j = json.loads((HERE / "cmb_techos_prep.json").read_text())
    tasks = []
    for theta, sigma, delta in CELLS:
        i_cell = CELL_INDEX[(sigma, theta)]
        d = np.load(HERE / f"cmb_techos_cell{i_cell}.npz")
        key = f"theta_{theta:g}"
        for j in range(len(d["idx_run"])):
            k = np.nonzero((d["pass_s"] & ~d["neg"]))[0][j]
            tasks.append((key, int(d["idx_run"][j]), d["Xs"][j], float(d["h"][k]), float(d["wm"][k])))
    res = run_isolated(tasks)
    out = {"prep": prep_j, "faithful_criterion": "max_a |rho_DE(CAMB) - rho_DE(path)| / rho_total(a) < 1e-3", "runs": res, "cells": {}}
    for theta, sigma, delta in CELLS:
        key = f"theta_{theta:g}"
        pr = prep_j["cells"][key]
        rs = [r for r in res if r["cell"] == key]
        okr = [r for r in rs if r["fail"] is None]
        fa = [r for r in okr if r["faithful"]]
        nf = [r for r in okr if not r["faithful"]]
        c = {"n_run": len(rs), "n_fail": len(rs) - len(okr), "fails": [r["fail"] for r in rs if r["fail"]][:5],
             "n_faithful": len(fa), "n_unfaithful": len(nf), "faithful": _stats(fa), "unfaithful": _stats(nf), "all_run": _stats(okr)}
        # biases as fractions of the smoothed passers and of the raw survivors
        n_pass, n_surv, n_neg = pr["n_pass_smoothed"], pr["n_survivors"], pr["n_rho_le_0_among_passers"]
        f_unf = len(nf) / len(okr) if okr else float("nan")
        n_unf_est = f_unf * pr["n_runnable"]  # unfaithful fraction of the run sample applied to all runnable
        d = np.load(HERE / f"cmb_techos_cell{CELL_INDEX[(sigma, theta)]}.npz")
        neg_idx = set(int(v) for v in d["surv"][d["pass_s"] & d["neg"]])
        unf_idx = set(r["j"] for r in nf)
        c["bias_sets"] = {"n_rho_le_0": len(neg_idx), "n_unfaithful_run": len(unf_idx), "intersection": len(neg_idx & unf_idx),
                          "union": len(neg_idx | unf_idx), "union_frac_of_survivors": len(neg_idx | unf_idx) / n_surv,
                          "union_frac_of_passers": len(neg_idx | unf_idx) / n_pass}
        c["bias"] = {"frac_rho_le_0_of_passers": n_neg / n_pass, "frac_unfaithful_of_run": f_unf,
                     "est_unfaithful_of_runnable": n_unf_est,
                     "frac_excluded_of_passers": (n_neg + n_unf_est) / n_pass,
                     "frac_excluded_of_survivors": (n_surv - n_pass + n_neg + n_unf_est) / n_surv,
                     "n_faithful_of_run": len(fa), "n_runnable": pr["n_runnable"], "n_pass_smoothed": n_pass, "n_survivors": n_surv}
        out["cells"][key] = c
        print(key, json.dumps({k: c[k] for k in ("n_run", "n_fail", "n_faithful", "n_unfaithful", "bias")}), flush=True)
        print("   faithful:", json.dumps(c["faithful"]), flush=True)
    out["seconds"] = time.time() - t0
    (HERE / "cmb_techos.json").write_text(json.dumps(out, indent=1), encoding="utf-8")


def trial() -> None:
    """One CAMB run of the first runnable realisation of each cell (w range, rho check, timing)."""
    for theta, sigma, delta in CELLS:
        i_cell = CELL_INDEX[(sigma, theta)]
        d = np.load(HERE / f"cmb_techos_cell{i_cell}.npz")
        if len(d["idx_run"]) == 0:
            print(f"theta={theta:g}: nothing runnable")
            continue
        k = np.nonzero((d["pass_s"] & ~d["neg"]))[0][0]
        t0 = time.time()
        r = one((f"theta_{theta:g}", int(d["idx_run"][0]), d["Xs"][0], float(d["h"][k]), float(d["wm"][k])))
        r.pop("cell")
        print(f"theta={theta:g}: {json.dumps(r)}  [{time.time() - t0:.0f} s]", flush=True)



def check_edge() -> None:
    """Extend the fidelity check below the table: CAMB extrapolates w(a < 1e-3) = w(1e-3); the path's rho_DE at z > 999 is
    Omega_L + X(z) - X(0) (X of the exacto grid up to z_*, held at its last value beyond). Error in H^2 units on
    a in [1e-6, 1], and w at the table edge. Re-classifies faithful_ext = faithful (table) and error below the table < 1e-3."""
    import camb

    d = json.loads((HERE / "cmb_techos.json").read_text())
    runs = d["runs"]
    cache = {}
    for theta, sigma, delta in CELLS:
        i_cell = CELL_INDEX[(sigma, theta)]
        dd = np.load(HERE / f"cmb_techos_cell{i_cell}.npz")
        cache[f"theta_{theta:g}"] = dd
    for r in runs:
        if r.get("fail"):
            continue
        dd = cache[r["cell"]]
        jj = int(np.nonzero(dd["idx_run"] == r["j"])[0][0])
        k = np.nonzero((dd["pass_s"] & ~dd["neg"]))[0][jj]
        Xs, h, wm = dd["Xs"][jj], float(dd["h"][k]), float(dd["wm"][k])
        rho, w = rho_w(Xs, h, wm)
        a_tab = np.exp(-U_PAD[TABLE])[::-1]
        wt = w[TABLE][::-1]
        pars = make_params(100 * h, wm, a_tab, wt)
        a_lo = np.geomspace(1e-6, 1e-3, 300)
        rc = camb.get_background(pars).get_dark_energy_rho_w(np.concatenate([a_lo, [1.0]]))[0]
        rc = rc[:-1] / rc[-1] * rho[N_PAD]  # units of rho_crit,0
        # path rho_DE below the table: X at u = -ln a on the exacto grid (u <= u_*), last value beyond
        u_lo = -np.log(a_lo)
        xs = np.interp(np.minimum(u_lo, U_PAD[-1]), U_PAD, Xs)
        rp = (h * h - wm - X.OMR) / (h * h) + xs - Xs[N_PAD]
        om, orad = wm / h**2, X.OMR / h**2
        rtot = om / a_lo**3 + orad / a_lo**4 + rp
        r["w_at_table_edge"] = float(wt[0])
        r["rho_err_H2_below_table"] = float(np.max(np.abs(rc - rp) / rtot))
        r["rho_DE_over_rho_tot_camb_at_1e-5"] = float((rc / rtot)[np.argmin(np.abs(a_lo - 1e-5))])
        r["faithful_table_only"] = r["faithful"]
        r["faithful"] = bool(r["faithful_table_only"] and r["rho_err_H2_below_table"] < 1e-3)
    d["faithful_criterion"] = ("max_a |rho_DE(CAMB) - rho_DE(path)| / rho_total(a) < 1e-3 on the table (1e-3 <= a <= 1) AND below it "
                               "(1e-6 <= a < 1e-3, where CAMB extrapolates w = w(1e-3))")
    # recompute the per-cell summaries
    for theta, sigma, delta in CELLS:
        key = f"theta_{theta:g}"
        pr = d["prep"]["cells"][key]
        rs = [r for r in runs if r["cell"] == key]
        okr = [r for r in rs if not r.get("fail")]
        fa = [r for r in okr if r["faithful"]]
        nf = [r for r in okr if not r["faithful"]]
        c = d["cells"][key]
        c.update(n_faithful=len(fa), n_unfaithful=len(nf), faithful=_stats(fa), unfaithful=_stats(nf), all_run=_stats(okr),
                 n_unfaithful_table_only=int(sum(not r["faithful_table_only"] for r in okr)),
                 unfaithful_w_at_edge=[round(r["w_at_table_edge"], 2) for r in nf],
                 unfaithful_rho_frac_1e5=[r["rho_DE_over_rho_tot_camb_at_1e-5"] for r in nf],
                 unfaithful_total=[round(r["total"], 2) for r in nf],
                 w_at_table_edge_faithful_p50_max=[float(v) for v in np.percentile([r["w_at_table_edge"] for r in fa], [50, 100])] if fa else None)
        n_pass, n_surv, n_neg = pr["n_pass_smoothed"], pr["n_survivors"], pr["n_rho_le_0_among_passers"]
        f_unf = len(nf) / len(okr) if okr else float("nan")
        n_unf_est = f_unf * pr["n_runnable"]
        dd = cache[key]
        neg_idx = set(int(v) for v in dd["surv"][dd["pass_s"] & dd["neg"]])
        unf_idx = set(r["j"] for r in nf)
        c["bias_sets"] = {"n_rho_le_0": len(neg_idx), "n_unfaithful_run": len(unf_idx), "intersection": len(neg_idx & unf_idx),
                          "union": len(neg_idx | unf_idx), "union_frac_of_survivors": len(neg_idx | unf_idx) / n_surv,
                          "union_frac_of_passers": len(neg_idx | unf_idx) / n_pass}
        c["bias"] = {"frac_rho_le_0_of_passers": n_neg / n_pass, "frac_unfaithful_of_run": f_unf, "est_unfaithful_of_runnable": n_unf_est,
                     "frac_excluded_of_passers": (n_neg + n_unf_est) / n_pass,
                     "frac_excluded_of_survivors": (n_surv - n_pass + n_neg + n_unf_est) / n_surv,
                     "n_faithful_of_run": len(fa), "n_runnable": pr["n_runnable"], "n_pass_smoothed": n_pass, "n_survivors": n_surv}
        print(key, "faithful", len(fa), "unfaithful", len(nf), "table-only unfaithful", c["n_unfaithful_table_only"],
              "| unfaithful: w_edge", c["unfaithful_w_at_edge"], "total", c["unfaithful_total"], "rho frac 1e-5", [f"{v:.1e}" for v in c["unfaithful_rho_frac_1e5"]])
        print("   faithful total p50/p95/max", [round(c["faithful"]["total"][q], 3) for q in ("p50", "p95", "max")], ">1", c["faithful"]["frac_total_gt_1"], ">4", c["faithful"]["frac_total_gt_4"],
              "picos p95", round(c["faithful"]["picos_joint"]["p95"], 4), "where", c["faithful"]["where_max_abs_z"], "w_edge faithful p50/max", c["w_at_table_edge_faithful_p50_max"])
        for key2 in ("corr_total_vs_abs_wp1_max", "corr_total_vs_path_amp_std", "corr_total_vs_path_amp_max"):
            print("     ", key2, {a: round(b, 2) for a, b in c["faithful"][key2].items()})
    (HERE / "cmb_techos.json").write_text(json.dumps(d, indent=1), encoding="utf-8")


if __name__ == "__main__":
    {"prep": prep, "run": run, "trial": trial, "check_edge": check_edge}[sys.argv[1]]()
