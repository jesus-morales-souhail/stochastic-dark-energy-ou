#!/usr/bin/env python3
"""
Exploration only (CMB_FINAL.md): do the histories that BAO + theta_* + omega_m + Pantheon+ allow at the ceilings of the exact
marginal -1.92 (OU theta = 1: 0.287; 3: 0.177; 10: 0.116; marginal_guiada.json) get discarded by the full CMB we have
(Planck TT l < 30, phiphi bandpowers, peak distance)? Nothing in results/ or papers changes.

1. Histories: N_SAMP draws of (dp, x) from the LINEAR POSTERIOR of marginal_guiada.py (Matheron's rule) with sigma_X at the
   ceiling, i.e. histories compatible with the data (not prior paths); the exact/linear weight w of marginal_guiada is recorded
   as information. rho_DE/rho_crit,0 = Omega_Lambda(h, omega_m) + x(z) - x(0) with (h, omega_m) = p0 + dp of the draw; draws
   with rho_DE <= 0 anywhere on the grid (0 <= z <= z*, constant beyond) are counted and do not enter CAMB; H^2 <= 0 likewise.
2. Smoothing as cmb_techos.py (Gaussian of std Delta in ln a, Delta = 0.05 for theta = 1, 3 and 0.0125 for theta = 10),
   with the backward AR(1) extension to u in [-PAD, 0) (prior conditional of a Markov OU, fresh normals) and, NEW, a forward
   extension to u in (u*, u* + PAD] with x = x(u*) constant.
3. Border fix (NEW, this test): above z* the path is held constant, X(z) = X(z*), so rho_DE does not change backwards and
   w = -1 there; the w(a) table given to CAMB covers the whole smoothed grid (a from 1 down to exp(-(u* + PAD))) and is extended
   to a = 1e-6 with w = -1. cmb_techos.py stopped the table at a = 1e-3 and CAMB extrapolated w(1e-3) backwards, which for
   paths with w(1e-3) > 1/3 made rho_DE grow into the past (the "border artefact" of CMB_TECHOS.md).
4. CAMB as cmb_techos.py: DarkEnergyPPF with the w table, linear lensing, H0 = 100 h and omega_m of the draw, omega_b, tau,
   n_s, A_s of Planck 2018; reference per history: LCDM with the same omega_b, omega_m and H0 solved for its 100 theta_*.
5. Faithfulness by two routes: CAMB's rho_DE(a) against the path's rho on 1e-6 <= a <= 1; faithful if
   max_a |rho_DE(CAMB) - rho_DE(path)| / rho_total(a) < 1e-3. Statistics (history minus reference): Delta chi^2 TT (2 <= l < 30),
   phiphi bandpowers, peak distance (picos.distance joint, l in [30, 2500]); total = sum; where the max |Delta/sigma| falls.

    OMP_NUM_THREADS=2 ../venv/bin/python cmb_final.py prep   # draws, counts, smoothing, one trial CAMB run -> cmb_final_prep.json
    OMP_NUM_THREADS=2 ../venv/bin/python cmb_final.py run    # CAMB on the first N_RUN runnable draws per ceiling -> cmb_final.json
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
from scipy.linalg import cho_factor, cho_solve

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parent / "cmb"))
import cmb_techos as CT  # noqa: E402
import exacto as X  # noqa: E402
import marginal_guiada as MG  # noqa: E402

CEIL = json.loads((HERE / "marginal_guiada.json").read_text(encoding="utf-8"))["thetas"]
CELLS = [(1.0, CEIL["theta_1"]["exact_192"], 0.05), (3.0, CEIL["theta_3"]["exact_192"], 0.05), (10.0, CEIL["theta_10"]["exact_192"], 0.0125)]
N_SAMP, N_RUN, SEED = 5000, 40, 20261119
PAD = CT.PAD  # 0.4 e-folds, backward and forward
A_FLOOR = 1e-6
G = MG.G
U = G.U  # exacto grid, 0 <= u <= u* = ln(1 + z*)
DU = X.DU
_npad = int(np.ceil(PAD / DU))
U_EXT = np.concatenate([-DU * np.arange(_npad, 0, -1), U, U[-1] + DU * np.arange(1, _npad + 1)])
N_PAD = _npad
I0 = N_PAD  # index of u = 0
MAIN = (np.arange(len(U_EXT)) >= N_PAD) & (np.arange(len(U_EXT)) < N_PAD + len(U))
TABLE = U_EXT >= 0  # the whole smoothed grid from a = 1 down to exp(-(u* + PAD))
A_EXTRA = np.geomspace(A_FLOOR, float(np.exp(-U_EXT[-1])), 60)[:-1]  # w = -1 nodes below the smoothed grid


def smooth_ext(Xr: np.ndarray, theta: float, sigma: float, delta: float, rng) -> np.ndarray:
    """Rows of x on the exacto grid -> smoothed rows on U_EXT (backward AR(1) extension, forward constant extension)."""
    back = CT.extend_backwards(Xr, theta, sigma, rng)[:, :N_PAD]
    fwd = np.repeat(Xr[:, [-1]], _npad, axis=1)
    Pp = np.concatenate([back, Xr, fwd], axis=1)
    S = CT.smoothing_matrix(U_EXT, delta)
    return Pp @ S.T


def rho_w(Xs: np.ndarray, h: float, wm: float):
    """rho_DE/rho_crit,0 and w on U_EXT for one smoothed path (rows of U_EXT)."""
    ol = (h * h - wm - X.OMR) / (h * h)
    rho = ol + Xs - Xs[I0]
    with np.errstate(invalid="ignore", divide="ignore"):
        w = -1.0 + np.gradient(np.log(np.where(rho > 0, rho, np.nan)), U_EXT) / 3.0  # d ln a = -du
    return rho, w


def prep() -> dict:
    Bm, A = MG.build_B(), MG.build_A()
    out = {"cells": {}, "n_samp": N_SAMP, "n_run": N_RUN, "seed": SEED, "pad": PAD, "a_floor": A_FLOOR, "grid": {"u_star": float(U[-1]), "z_star": float(np.exp(U[-1]) - 1), "n_ext": int(len(U_EXT))}}
    for ic, (theta, sigma, delta) in enumerate(CELLS):
        t0 = time.time()
        key = f"theta_{theta:g}"
        K = MG.kernel(theta)
        AK = A @ K
        Msig = AK @ A.T
        rng = np.random.default_rng([SEED, ic])
        r = MG.D - MG.B0
        Sf = cho_factor(MG.C + sigma**2 * Msig, lower=True)
        SB = cho_solve(Sf, Bm)
        Sp = np.linalg.inv(Bm.T @ SB)
        dp_hat = Sp @ (Bm.T @ cho_solve(Sf, r))
        dp = dp_hat[None, :] + rng.standard_normal((N_SAMP, 3)) @ np.linalg.cholesky(Sp).T
        x0 = MG.S.ou_paths(G, sigma, theta, N_SAMP, rng)
        n0 = rng.standard_normal((N_SAMP, MG.NROWS)) @ MG.LC.T
        R = r[None, :] - dp @ Bm.T - x0 @ A.T - n0
        x = x0 + sigma**2 * cho_solve(Sf, R.T).T @ AK
        # linear and exact chi2, weight (information)
        chi_lin = MG.chi2_of(MG.B0[None, :] + dp @ Bm.T + x @ A.T)
        p = MG.P0[None, :] + dp
        m_ex, bad = MG.exact_model(p[:, 0], p[:, 1], p[:, 2], x - x[:, [0]])
        chi_ex = np.where(bad, np.inf, MG.chi2_of(np.where(bad[:, None], 0.0, m_ex)))
        lw = -0.5 * (chi_ex - chi_lin)
        w = np.where(bad, 0.0, np.exp(np.minimum(lw, 700)))
        h, wm = p[:, 0], p[:, 1]
        ol = (h * h - wm - X.OMR) / (h * h)
        rho_raw = ol[:, None] + x - x[:, [0]]
        neg_raw = (rho_raw <= 0).any(axis=1)
        rec = {"theta": theta, "sigma_X": sigma, "Delta": delta, "n_samp": N_SAMP, "dp_hat": dp_hat.tolist(), "dp_std": np.sqrt(np.diag(Sp)).tolist(),
               "n_H2_le_0": int(bad.sum()), "n_rho_le_0_raw": int(neg_raw.sum()), "frac_rho_le_0_raw": float(neg_raw.mean()),
               "rho_min_over_ol_p05_p50": [float(v) for v in np.percentile((rho_raw.min(axis=1) / ol)[~bad], [5, 50])],
               "ess_all": float(w.sum() ** 2 / np.sum(w**2)) if w.sum() > 0 else 0.0, "chi2_lin_mean": float(chi_lin.mean()),
               "chi2_ex_mean_finite": float(np.mean(chi_ex[np.isfinite(chi_ex)])), "chi2_lcdm": float(MG._c0)}
        # smoothing + extension; rho on the smoothed grid
        keep = ~bad & ~neg_raw
        idx_keep = np.nonzero(keep)[0][: 4 * N_RUN]  # enough candidates
        rng2 = np.random.default_rng([SEED, ic, 1])
        Xs = smooth_ext(x[idx_keep], theta, sigma, delta, rng2)
        neg_s = np.zeros(len(idx_keep), bool)
        wmax = np.full(len(idx_keep), np.nan)
        for j in range(len(idx_keep)):
            rho, ww = rho_w(Xs[j], h[idx_keep[j]], wm[idx_keep[j]])
            neg_s[j] = bool(np.any(rho[TABLE] <= 0))
            if not neg_s[j]:
                wmax[j] = float(np.nanmax(np.abs(ww[TABLE] + 1)))
        rec["n_candidates_smoothed"] = int(len(idx_keep))
        rec["n_rho_le_0_smoothed_among_candidates"] = int(neg_s.sum())
        run_j = np.nonzero(~neg_s)[0][:N_RUN]
        idx_run = idx_keep[run_j]
        rec["n_run"] = int(len(idx_run))
        rec["abs_wp1_max_p50_p95_run"] = [float(v) for v in np.percentile(wmax[run_j], [50, 95])]
        rec["w_exact_over_linear_run_p16_p50_p84"] = [float(v) for v in np.percentile(w[idx_run], [16, 50, 84])]
        rec["chi2_ex_minus_lin_run_p16_p50_p84"] = [float(v) for v in np.percentile((chi_ex - chi_lin)[idx_run], [16, 50, 84])]
        rec["variance_kept_numerical"] = float(np.mean(np.var(Xs[run_j][:, MAIN], axis=1) / np.var(x[idx_run], axis=1)))
        np.savez_compressed(HERE / f"cmb_final_{key}.npz", idx_run=idx_run, Xs=Xs[run_j], h=h[idx_run], wm=wm[idx_run], w=w[idx_run],
                            chi_lin=chi_lin[idx_run], chi_ex=chi_ex[idx_run], U_EXT=U_EXT, neg_raw=neg_raw, bad=bad)
        rec["seconds"] = time.time() - t0
        out["cells"][key] = rec
        print(json.dumps(rec), flush=True)
    (HERE / "cmb_final_prep.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    return out


def one(task):
    """task = (cell key, idx, Xs on U_EXT, h, wm). History + reference, faithfulness on 1e-6 <= a <= 1, statistics."""
    import camb
    import ou_cmb as M
    import estadisticos as ES
    import picos as PI

    key, j, Xs, h, wm = task[:5]
    out = {"cell": key, "j": int(j), "fail": None}
    try:
        rho, w = rho_w(Xs, h, wm)
        a_grid = np.exp(-U_EXT[TABLE])[::-1]  # increasing a, from exp(-(u*+PAD)) to 1
        w_grid = w[TABLE][::-1]
        rp_grid = rho[TABLE][::-1]
        a = np.concatenate([A_EXTRA, a_grid])
        wt = np.concatenate([-np.ones(len(A_EXTRA)), w_grid])
        rp = np.concatenate([np.full(len(A_EXTRA), rp_grid[0]), rp_grid])  # constant rho below the grid
        pars = CT.make_params(100 * h, wm, a, wt)
        rc = camb.get_background(pars).get_dark_energy_rho_w(a)[0]
        rc = rc / rc[-1] * rp[-1]
        om, orad = wm / h**2, X.OMR / h**2
        rtot = om / a**3 + orad / a**4 + rp
        err = np.abs(rc - rp) / rtot
        out["rho_err_H2"] = float(np.max(err))
        out["rho_err_H2_table"] = float(np.max(err[len(A_EXTRA):]))
        out["rho_err_rel"] = float(np.max(np.abs(rc / rp - 1)))
        out["faithful"] = bool(out["rho_err_H2"] < 1e-3)
        out["abs_wp1_max"] = float(np.nanmax(np.abs(w_grid + 1)))
        out["w_at_grid_floor"] = float(w_grid[0])
        out["rho_min_over_rho0"] = float(rp.min() / rp[-1])
        dxs = Xs[MAIN] - Xs[I0]
        out["path_amp_std"] = float(np.std(dxs))
        out["path_amp_max"] = float(np.max(np.abs(dxs)))
        r = CT.spectra(pars)
        out["thetastar"] = r["thetastar"]
        out["H0"] = 100 * h
        make = lambda H: CT.make_params(H, wm)
        H0_ref, th_ref, _ = M.solve_H0_for_thetastar(r["thetastar"], make)
        ref = CT.spectra(make(H0_ref))
        out["H0_ref"], out["thetastar_ref"] = float(H0_ref), ref["thetastar"]
        tg = ES.targets()
        d_tt = float(ES.chi2(r["TT"], *tg["TT"]) - ES.chi2(ref["TT"], *tg["TT"]))
        d_pp = float(ES.chi2(r["pp"], *tg["pp"]) - ES.chi2(ref["pp"], *tg["pp"]))
        sig_tt = 0.5 * (tg["TT"][1] + tg["TT"][2])
        z_tt = (r["TT"] - ref["TT"]) / sig_tt
        z_pp = (r["pp"] - ref["pp"]) / tg["pp"][1]
        pk = PI.distance(ref["lensed"], r["lensed"])
        jl = CT.joint_per_l(ref["lensed"], r["lensed"])
        out.update(dchi2_TT=d_tt, dchi2_pp=d_pp, picos_joint=pk["joint"], picos_TT=pk["TT"], picos_lmax_joint=pk["lmax_joint"],
                   total=d_tt + d_pp + pk["joint"], max_z_TT=float(np.max(np.abs(z_tt))), l_max_z_TT=int(2 + np.argmax(np.abs(z_tt))),
                   max_z_pp=float(np.max(np.abs(z_pp))), max_z_picos=float(np.sqrt(jl.max())), l_max_z_picos=int(PI.LMIN + np.argmax(jl)))
    except Exception as exc:  # noqa: BLE001
        out["fail"] = f"{type(exc).__name__}: {exc}".replace("\n", " ")[:300]
    return out


def _child(task, q):
    q.put(one(task))


def run_isolated(tasks, n_par: int = 8, timeout: float = 900.0) -> list[dict]:
    """As cmb_techos.run_isolated (one spawned process per task), with this test's `one` and log."""
    ctx = mp.get_context("spawn")
    log = open(HERE / "cmb_final_progress.log", "w")
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
            rr = res[i]
            log.write(f"{time.time() - t0:7.0f} s  {rr['cell']} j={rr['j']}  " + (f"FAIL {rr['fail']}" if rr.get("fail") else
                      f"faithful={rr['faithful']} errH2={rr['rho_err_H2']:.2e} total={rr['total']:+.3f}") + "\n")
            log.flush()
    log.close()
    return res


def tasks_from_prep():
    tasks = []
    for theta, sigma, delta in CELLS:
        key = f"theta_{theta:g}"
        d = np.load(HERE / f"cmb_final_{key}.npz")
        for j in range(len(d["idx_run"])):
            tasks.append((key, int(d["idx_run"][j]), d["Xs"][j], float(d["h"][j]), float(d["wm"][j])))
    return tasks


def run() -> None:
    t0 = time.time()
    prep_j = json.loads((HERE / "cmb_final_prep.json").read_text())
    res = run_isolated(tasks_from_prep())
    out = {"prep": prep_j, "faithful_criterion": "max_a |rho_DE(CAMB) - rho_DE(path)| / rho_total(a) < 1e-3 on 1e-6 <= a <= 1",
           "border_fix": "X(z) = X(z*) above z*; w table to a = 1e-6 with w = -1 below the smoothed grid; forward constant padding of PAD e-folds before smoothing",
           "runs": res, "cells": {}}
    for theta, sigma, delta in CELLS:
        key = f"theta_{theta:g}"
        rs = [r for r in res if r["cell"] == key]
        okr = [r for r in rs if r["fail"] is None]
        fa = [r for r in okr if r["faithful"]]
        nf = [r for r in okr if not r["faithful"]]
        c = {"n_run": len(rs), "n_fail": len(rs) - len(okr), "fails": [r["fail"] for r in rs if r["fail"]][:5], "n_faithful": len(fa),
             "n_unfaithful": len(nf), "frac_faithful_of_run": (len(fa) / len(okr)) if okr else None, "faithful": CT._stats(fa),
             "unfaithful": CT._stats(nf), "all_run": CT._stats(okr)}
        if okr:
            c["rho_err_H2_table_p50_max"] = [float(np.median([r["rho_err_H2_table"] for r in okr])), float(max(r["rho_err_H2_table"] for r in okr))]
            c["w_at_grid_floor_p50_min_max"] = [float(np.median([r["w_at_grid_floor"] for r in okr])), float(min(r["w_at_grid_floor"] for r in okr)), float(max(r["w_at_grid_floor"] for r in okr))]
        out["cells"][key] = c
        print(key, json.dumps({k: c[k] for k in ("n_run", "n_fail", "n_faithful", "n_unfaithful")}), flush=True)
        print("   faithful:", json.dumps(c["faithful"]), flush=True)
    out["seconds"] = time.time() - t0
    (HERE / "cmb_final.json").write_text(json.dumps(out, indent=1), encoding="utf-8")


def trial() -> None:
    for task in tasks_from_prep()[::N_RUN]:
        t0 = time.time()
        r = one(task)
        r.pop("cell")
        print(f"{task[0]}: {json.dumps(r)}  [{time.time() - t0:.0f} s]", flush=True)


if __name__ == "__main__":
    {"prep": prep, "run": run, "trial": trial}[sys.argv[1]]()
