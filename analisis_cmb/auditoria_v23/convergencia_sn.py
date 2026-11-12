#!/usr/bin/env python3
"""
Exploration only (CONVERGENCIA_SN.md): does the marginal-likelihood bound converge with the number of paths? Nothing in
results/ or papers changes; the bounds of EXACTO_SN.md and PRECISION_SN.md are not touched.

Set-up of precision_sn.py unchanged: combination BAO + theta_* + omega_m + SN only, the five theta, N = 20000 paths per cell
(new seeds [SEED, 300 + index], different from those of N = 2000), cells sigma_X in {0.03, 0.05, 0.07, 0.1, 0.15, 0.2, 0.3, 0.4}
(plus 0.5 and 0.7 for theta = 0.1 and 0.3). Paths generated in chunks of 2000 inside each worker (memory).
1. 5 % bound, marginal -1.92 and -1.35, maximum of ln L and ESS with the first 2000, 5000, 10000 and all 20000 paths.
2. Bootstrap (1000) of the three bounds with N = 20000.
3. Convergence: change of the marginal bound between 10000 and 20000 (< 5 % = converged).

    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 ../venv/bin/python convergencia_sn.py    # writes convergencia_sn.{json,txt}, convergencia_sn_chi2.npz
"""

from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import exacto_sn as S  # noqa: E402
import precision_sn as PS  # noqa: E402

N_TOT, CHUNK = 20000, 2000
NS = (2000, 5000, 10000, 20000)
SIG_MAIN = (0.03, 0.05, 0.07, 0.1, 0.15, 0.2, 0.3, 0.4)
SIG_EXT = SIG_MAIN + (0.5, 0.7)
CELLS = [(s, t) for t in S.THETAS for s in (SIG_EXT if t in (0.1, 0.3) else SIG_MAIN)]
N_BOOT, SEED_BOOT = 1000, 20261114
COMBO = "BAO+CMB+SN"


def _cell(args):
    sigma, theta = args
    i = CELLS.index((sigma, theta))
    rng = np.random.default_rng([S.SEED, 300 + i])
    h0, wm0 = S.P0
    om0, or0 = wm0 / h0**2, S.OMR / h0**2
    ok = np.zeros(N_TOT, bool)
    chi = np.full(N_TOT, np.nan)
    sn = S.SN["all"]
    for c0 in range(0, N_TOT, CHUNK):
        Xp = S.ou_paths(S.G0, sigma, theta, CHUNK, rng)
        dX = Xp - Xp[:, [0]]
        okc = (om0 * S.G0.p3 + or0 * S.G0.p4 + (1 - om0 - or0) + dX).min(axis=1) > 0.0
        ok[c0:c0 + CHUNK] = okc
        for j in np.nonzero(okc)[0]:
            chi[c0 + j] = S.fit_combo(S.G0, dX[j], sn, True)[0]
    return sigma, theta, ok, chi


def stats(sig, oks, chis, obs, n):
    P, lnL, ess, wfrac = [], [], [], []
    for ok, chi in zip(oks, chis):
        ok, chi = ok[:n], chi[:n]
        w = np.where(ok, np.exp(-0.5 * (np.where(ok, chi, 0.0) - obs)), 0.0)
        P.append(float(np.mean(chi[ok] <= obs)) if ok.any() else 0.0)
        lnL.append(float(np.log(max(w.mean(), 1e-300))))
        ess.append(float(w.sum() ** 2 / np.sum(w**2)) if w.sum() > 0 else 0.0)
        wfrac.append(float(w.max() / w.sum()) if w.sum() > 0 else 1.0)
    sig = np.array(sig)
    rec = {"n": n, "P": P, "lnL": lnL, "ess": ess, "w_max_frac": wfrac, "p5": PS.crossing(sig, P, 0.05),
           "marg_192": PS.lnl_bound(sig, np.array(lnL), 1.92), "marg_135": PS.lnl_bound(sig, np.array(lnL), 1.35)}
    allv = np.concatenate([[0.0], lnL])
    rec["lnL_max_at_sigma"] = float(np.concatenate([[0.0], sig])[int(np.argmax(allv))])
    rec["lnL_max"] = float(allv.max())
    for key in ("marg_192", "p5"):
        b = rec[key]
        if np.isfinite(b) and b > 0:
            k = int(np.argmin(np.abs(np.log(sig / b))))
            rec[f"ess_at_{key}"] = {"sigma_cell": float(sig[k]), "ess": ess[k], "w_max_frac": wfrac[k], "n_hits": int(round(P[k] * oks[k][:n].sum()))}
    return rec


def main() -> None:
    t0 = time.time()
    obs = json.loads((HERE / "exacto_sn.json").read_text())["checks"]["chi2_obs"][COMBO]
    if "--from-npz" in sys.argv:  # analysis only, from a saved run
        z = np.load(HERE / "convergencia_sn_chi2.npz")
        cells = {(s, t): (z[f"theta_{t:g}__sigma_{s:g}"][:, 0].astype(bool), z[f"theta_{t:g}__sigma_{s:g}"][:, 1]) for s, t in CELLS}
    else:
        with Pool() as pool:
            res = pool.map(_cell, CELLS, chunksize=1)
        print(f"celdas: {time.time() - t0:.0f} s", flush=True)
        cells = {(s, t): (ok, chi) for s, t, ok, chi in res}
        np.savez_compressed(HERE / "convergencia_sn_chi2.npz", columns=np.array(["ok", "chi2"]),
                            **{f"theta_{t:g}__sigma_{s:g}": np.column_stack([v[0], v[1]]) for (s, t), v in cells.items()})
    prev = json.loads((HERE / "precision_sn.json").read_text())["table"]
    rng = np.random.default_rng(SEED_BOOT)
    out = {"combo": COMBO, "chi2_obs": obs, "n_tot": N_TOT, "ns": NS, "cells": [list(c) for c in CELLS], "n_boot": N_BOOT, "thetas": {}}
    for t in S.THETAS:
        sig = [s for s, tt in CELLS if tt == t]
        oks = [cells[(s, t)][0] for s in sig]
        chis = [cells[(s, t)][1] for s in sig]
        rec = {"sigmas": sig, "frac_H2_le_0": [float(1 - ok.mean()) for ok in oks], "by_n": {}}
        for n in NS:
            rec["by_n"][str(n)] = stats(sig, oks, chis, obs, n)
        # bootstrap at N = 20000
        sigA = np.array(sig)
        boots = {"p5": [], "marg_192": [], "marg_135": []}
        hits = [np.where(ok, chi <= obs, False).astype(float) for ok, chi in zip(oks, chis)]
        ws = [np.where(ok, np.exp(-0.5 * (np.where(ok, chi, 0.0) - obs)), 0.0) for ok, chi in zip(oks, chis)]
        for b in range(N_BOOT):
            P, lnL = [], []
            for ok, h, w in zip(oks, hits, ws):
                idx = rng.integers(0, N_TOT, size=N_TOT)
                P.append(h[idx].sum() / max(ok[idx].sum(), 1))
                lnL.append(np.log(max(w[idx].mean(), 1e-300)))
            boots["p5"].append(PS.crossing(sigA, P, 0.05))
            boots["marg_192"].append(PS.lnl_bound(sigA, np.array(lnL), 1.92))
            boots["marg_135"].append(PS.lnl_bound(sigA, np.array(lnL), 1.35))
        rec["boot_20000"] = {}
        for k, v in boots.items():
            v = np.array(v)
            fin = np.isfinite(v) & (v > 0)
            rec["boot_20000"][k] = {"p16_p84": [float(np.percentile(v[fin], 16)), float(np.percentile(v[fin], 84))] if fin.sum() > 10 else None,
                                    "frac_no_bound": float(np.mean(~np.isfinite(v)))}
        r10, r20 = rec["by_n"]["10000"], rec["by_n"]["20000"]
        rec["converged_marg_192"] = bool(np.isfinite(r10["marg_192"]) and np.isfinite(r20["marg_192"]) and r10["marg_192"] > 0 and r20["marg_192"] > 0
                                     and abs(r20["marg_192"] / r10["marg_192"] - 1) < 0.05)
        rec["change_10000_to_20000"] = (float(r20["marg_192"] / r10["marg_192"] - 1) if np.isfinite(r10["marg_192"]) and np.isfinite(r20["marg_192"]) and r10["marg_192"] > 0 else None)
        p = prev[f"theta_{t:g}|{COMBO}"]
        rec["precision_sn_N2000"] = {"p5_all_cells": p["p5_all_cells"], "marg_192": p["marg_192"], "marg_135": p["marg_135"], "lnL_max_at_sigma": p["lnL_max_at_sigma"],
                                     "lnL_max": p["lnL_max"], "ess_at_marg_bound": p.get("ess_at_marg_bound")}
        out["thetas"][f"theta_{t:g}"] = rec
        print(f"theta {t:g}: " + " | ".join(f"N={n}: p5={PS.fmt(rec['by_n'][str(n)]['p5'])} m192={PS.fmt(rec['by_n'][str(n)]['marg_192'])} m135={PS.fmt(rec['by_n'][str(n)]['marg_135'])} "
                                            f"max@{rec['by_n'][str(n)]['lnL_max_at_sigma']:g}({rec['by_n'][str(n)]['lnL_max']:+.2f}) ess={rec['by_n'][str(n)].get('ess_at_marg_192', {}).get('ess', float('nan')):.0f}" for n in NS)
              + f" | boot20000 {rec['boot_20000']} | conv {rec['converged_marg_192']} ({rec['change_10000_to_20000']})", flush=True)
    out["runtime_s"] = time.time() - t0
    (HERE / "convergencia_sn.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    L = [f"{COMBO}, chi2_obs = {obs:.3f}; N = {N_TOT} caminos por celda, {len(CELLS)} celdas", ""]
    for t in S.THETAS:
        rec = out["thetas"][f"theta_{t:g}"]
        L.append(f"theta = {t:g}: N | p5 | marg -1.92 | marg -1.35 | max lnL (sigma) | ESS en la cota marg | ESS en la cota p5")
        for n in NS:
            r = rec["by_n"][str(n)]
            L.append(f"   {n:6d} | {PS.fmt(r['p5'])} | {PS.fmt(r['marg_192'])} | {PS.fmt(r['marg_135'])} | {r['lnL_max']:+.2f} ({r['lnL_max_at_sigma']:g}) | "
                     f"{r.get('ess_at_marg_192', {}).get('ess', float('nan')):.0f} | {r.get('ess_at_p5', {}).get('ess', float('nan')):.0f}")
        L.append(f"   bootstrap N=20000: {json.dumps(rec['boot_20000'])}; cambio 10000->20000 de marg -1.92: {rec['change_10000_to_20000']}; convergida: {rec['converged_marg_192']}")
        L.append("")
    (HERE / "convergencia_sn.txt").write_text("\n".join(L) + "\n", encoding="utf-8")
    print("\n".join(L))


if __name__ == "__main__":
    main()
