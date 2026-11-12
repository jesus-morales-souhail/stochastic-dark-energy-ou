#!/usr/bin/env python3
"""
Exploration only (PRECISION_SN.md): Monte Carlo error of the 5 % ceilings of EXACTO_SN.md and their dependence on the 5 % rule.
Nothing in results/ or papers changes; the 5 % rule and the bounds of EXACTO_SN.md are not touched (this goes alongside).

Set-up of exacto_sn.py unchanged, same seeds. New cells sigma_X in {0.07, 0.15, 0.4} for every theta (seeds [SEED, 200 + index]);
for them the BAO + theta_* + omega_m (no SN) chi^2 is also computed (exacto.fit), but the no-SN column of the table uses the
cells of exacto_chi2.npz only, as asked. chi^2_min of every realisation, cell (old and new) and combination -> precision_sn_chi2.npz.
1. 5 % bound on the extended sigma grid (12 values) and its bootstrap error over the realisations (1000 resamples per cell,
   independent per cell; percentiles 16 and 84 of the bound).
2. Exact marginal likelihood by Monte Carlo: L(sigma, theta) = mean over the 2000 paths of exp(-(chi^2_min - chi^2_obs)/2),
   discards (H^2 <= 0) with weight 0; ln L on the grid including sigma = 0 (L = 1); bound where ln L drops 1.92 (and 1.35)
   below its maximum over the grid, interpolated in ln sigma (if the crossing is before the first node, "< 0.01"); error with
   the same bootstrap; effective sample size (sum w)^2 / sum w^2 at the cell nearest the bound; location of the maximum.
3. Table per theta and combination (SN; BAO + SN; BAO + theta_* + omega_m + SN; BAO + theta_* + omega_m without SN), with the
   linear marginalised bounds of exacto_sn.json and cmb_fondo.json (Lambda variant).

    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 ../venv/bin/python precision_sn.py    # writes precision_sn.{json,txt}, precision_sn_chi2.npz
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
import exacto as X  # noqa: E402

NEW_SIGMAS = (0.07, 0.15, 0.4)
SIG_ALL = tuple(sorted(set(S.SIGMAS) | set(NEW_SIGMAS)))
THETAS = S.THETAS
COMBOS = ("SN", "BAO+SN", "BAO+CMB+SN", "BAO+CMB")
N_BOOT, SEED_BOOT = 1000, 20261113
THR = {"192": 1.92, "135": 1.35}


def _cell_new(args):
    sigma, theta = args
    i = [(s, t) for t in THETAS for s in NEW_SIGMAS].index((sigma, theta))
    rng = np.random.default_rng([S.SEED, 200 + i])
    Xp = S.ou_paths(S.G0, sigma, theta, S.N_REAL, rng)
    dX = Xp - Xp[:, [0]]
    h0, wm0 = S.P0
    om0, or0 = wm0 / h0**2, S.OMR / h0**2
    ok = (om0 * S.G0.p3 + or0 * S.G0.p4 + (1 - om0 - or0) + dX).min(axis=1) > 0.0
    chi = np.full((S.N_REAL, 4), np.nan)
    X.P0 = S.P0
    for j in np.nonzero(ok)[0]:
        c = S.chi2_all(S.G0, dX[j])[:3]
        chi[j, :3] = c
        chi[j, 3] = X.fit(dX[j], 15)[0]
    return sigma, theta, ok, chi


def crossing(sig, p, target):
    sig, p = np.asarray(sig, float), np.asarray(p, float)
    if p[0] <= target:
        return -1.0  # below the first node
    for k in range(len(sig) - 1):
        if p[k] > target >= p[k + 1]:
            t = (p[k] - target) / (p[k] - p[k + 1])
            return float(np.exp(np.log(sig[k]) + t * (np.log(sig[k + 1]) - np.log(sig[k]))))
    return np.inf  # never


def lnl_bound(sig, lnL, drop):
    """lnL on nodes sig (sigma > 0); the sigma = 0 node has lnL = 0. Bound where lnL falls `drop` below the max (max over
    all nodes incl. sigma = 0), from the max outwards, interpolated in ln sigma."""
    vals = np.concatenate([[0.0], lnL])
    nodes = np.concatenate([[0.0], sig])
    k0 = int(np.argmax(vals))
    thr = vals[k0] - drop
    for k in range(k0, len(nodes) - 1):
        if vals[k] >= thr > vals[k + 1]:
            if nodes[k] == 0.0:
                return -1.0  # crossing before the first node
            t = (vals[k] - thr) / (vals[k] - vals[k + 1])
            return float(np.exp(np.log(nodes[k]) + t * (np.log(nodes[k + 1]) - np.log(nodes[k]))))
    return np.inf


def fmt(v):
    if v is None or (isinstance(v, float) and np.isinf(v)):
        return "sin cota"
    if v == -1.0:
        return "< 0.01"
    return f"{v:.3f}"


def main() -> None:
    t0 = time.time()
    # old cells
    old_sn = np.load(HERE / "exacto_sn_chi2.npz")
    old_15 = np.load(HERE / "exacto_chi2.npz")
    obs_sn = json.loads((HERE / "exacto_sn.json").read_text())["checks"]["chi2_obs"]
    obs = {"SN": obs_sn["SN"], "BAO+SN": obs_sn["BAO+SN"], "BAO+CMB+SN": obs_sn["BAO+CMB+SN"], "BAO+CMB": json.loads((HERE / "exacto.json").read_text())["chi2_obs"]["BAO+CMB"]}
    cells = {}  # (theta, sigma, combo) -> (ok, chi)
    for t in THETAS:
        for s in S.SIGMAS:
            a = old_sn[f"theta_{t:g}__sigma_{s:g}"]
            ok = a[:, 0].astype(bool)
            for j, c in enumerate(("SN", "BAO+SN", "BAO+CMB+SN")):
                cells[(t, s, c)] = (ok, a[:, 2 + j])
            if s in X.SIGMAS:
                b = old_15[f"theta_{t:g}__sigma_{s:g}"]
                cells[(t, s, "BAO+CMB")] = (b[0].astype(bool), b[1])
    with Pool() as pool:
        new = pool.map(_cell_new, [(s, t) for t in THETAS for s in NEW_SIGMAS], chunksize=1)
    for s, t, ok, chi in new:
        for j, c in enumerate(COMBOS):
            cells[(t, s, c)] = (ok, chi[:, j])
    print(f"celdas nuevas: {time.time() - t0:.0f} s", flush=True)
    np.savez_compressed(HERE / "precision_sn_chi2.npz", columns=np.array(["ok", "chi2"]),
                        **{f"theta_{t:g}__sigma_{s:g}__{c}": np.column_stack([v[0], v[1]]) for (t, s, c), v in cells.items()})
    lin_sn = json.loads((HERE / "exacto_sn.json").read_text())["linear"]
    lin_cf = json.loads((HERE / "cmb_fondo.json").read_text())["fits"]["BAO+CMB|Lambda_absorbs_X0"]["marginalised"]
    rng = np.random.default_rng(SEED_BOOT)
    out = {"sigmas_all": SIG_ALL, "new_sigmas": NEW_SIGMAS, "thetas": THETAS, "combos": COMBOS, "chi2_obs": obs, "n_boot": N_BOOT, "table": {}}
    for t in THETAS:
        for c in COMBOS:
            sig = [s for s in SIG_ALL if (t, s, c) in cells]
            sig_old = [s for s in sig if s in (X.SIGMAS if c == "BAO+CMB" else S.SIGMAS)]
            P, L, ess, wmax_frac = [], [], [], []
            Pb, Lb = np.empty((N_BOOT, len(sig))), np.empty((N_BOOT, len(sig)))
            for k, s in enumerate(sig):
                ok, chi = cells[(t, s, c)]
                w = np.where(ok, np.exp(-0.5 * (np.where(ok, chi, 0.0) - obs[c])), 0.0)
                hit = np.where(ok, chi <= obs[c], False).astype(float)
                P.append(float(hit[ok].mean()) if ok.any() else 0.0)
                L.append(float(w.mean()))
                ess.append(float(w.sum() ** 2 / np.sum(w**2)) if w.sum() > 0 else 0.0)
                wmax_frac.append(float(w.max() / w.sum()) if w.sum() > 0 else 1.0)
                idx = rng.integers(0, len(chi), size=(N_BOOT, len(chi)))
                okb = ok[idx]
                Pb[:, k] = (hit[idx]).sum(axis=1) / np.maximum(okb.sum(axis=1), 1)
                Lb[:, k] = w[idx].mean(axis=1)
            lnL = np.log(np.maximum(L, 1e-300))
            rec = {"sigmas": sig, "P": P, "lnL": lnL.tolist(), "ess": ess, "w_max_frac": wmax_frac,
                   "lnL_max_at_sigma": float(np.concatenate([[0.0], sig])[int(np.argmax(np.concatenate([[0.0], lnL])))]),
                   "lnL_max": float(max(0.0, lnL.max()))}
            # 5 % bound: old cells (as EXACTO_SN) and all cells
            io = [sig.index(s) for s in sig_old]
            rec["p5_old_cells"] = crossing(sig_old, [P[i] for i in io], 0.05)
            rec["p5_all_cells"] = crossing(sig, P, 0.05)
            bb = np.array([crossing(sig, Pb[b], 0.05) for b in range(N_BOOT)])
            fin = np.isfinite(bb) & (bb > 0)
            rec["p5_boot_p16_p84"] = [float(np.percentile(bb[fin], 16)), float(np.percentile(bb[fin], 84))] if fin.sum() > 10 else None
            rec["p5_boot_frac_no_bound"] = float(np.mean(~np.isfinite(bb)))
            for name, drop in THR.items():
                rec[f"marg_{name}"] = lnl_bound(np.array(sig), lnL, drop)
                lb = np.array([lnl_bound(np.array(sig), np.log(np.maximum(Lb[b], 1e-300)), drop) for b in range(N_BOOT)])
                fin = np.isfinite(lb) & (lb > 0)
                rec[f"marg_{name}_boot_p16_p84"] = [float(np.percentile(lb[fin], 16)), float(np.percentile(lb[fin], 84))] if fin.sum() > 10 else None
                rec[f"marg_{name}_boot_frac_no_bound"] = float(np.mean(~np.isfinite(lb)))
            b192 = rec["marg_192"]
            if np.isfinite(b192) and b192 > 0:
                kk = int(np.argmin(np.abs(np.log(np.array(sig) / b192))))
                rec["ess_at_marg_bound"] = {"sigma_cell": sig[kk], "ess": ess[kk], "w_max_frac": wmax_frac[kk]}
            b5 = rec["p5_all_cells"]
            if np.isfinite(b5) and b5 > 0:
                kk = int(np.argmin(np.abs(np.log(np.array(sig) / b5))))
                rec["ess_at_p5_bound"] = {"sigma_cell": sig[kk], "ess": ess[kk], "n_hits": int(round(P[kk] * cells[(t, sig[kk], c)][0].sum()))}
            key = f"theta_{t:g}"
            rec["linear_192"] = (lin_cf[f"theta_{float(t)}"]["192"] if c == "BAO+CMB" else lin_sn[c][key]["192"])
            rec["linear_135"] = (lin_cf[f"theta_{float(t)}"]["135"] if c == "BAO+CMB" else lin_sn[c][key]["135"])
            out["table"][f"{key}|{c}"] = rec
            print(f"{key:9s} {c:11s} p5 old={fmt(rec['p5_old_cells']):>8s} all={fmt(rec['p5_all_cells']):>8s} boot={rec['p5_boot_p16_p84']}  "
                  f"marg192={fmt(rec['marg_192']):>8s} {rec['marg_192_boot_p16_p84']} marg135={fmt(rec['marg_135']):>8s}  Lmax@{rec['lnL_max_at_sigma']:g} ({rec['lnL_max']:+.2f})  "
                  f"ess={rec.get('ess_at_marg_bound')}", flush=True)
    out["runtime_s"] = time.time() - t0
    (HERE / "precision_sn.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    lines = []
    for c in COMBOS:
        lines.append(f"{c}: theta | p5 viejas | p5 todas (16-84) | marg -1.92 (16-84) | marg -1.35 | lineal -1.92 / -1.35 | max lnL en sigma | ESS en la cota marg")
        for t in THETAS:
            r = out["table"][f"theta_{t:g}|{c}"]
            e = r.get("ess_at_marg_bound")
            lines.append(f"  {t:5g} | {fmt(r['p5_old_cells'])} | {fmt(r['p5_all_cells'])} ({r['p5_boot_p16_p84']}) | {fmt(r['marg_192'])} ({r['marg_192_boot_p16_p84']}) | {fmt(r['marg_135'])} | "
                         f"{fmt(r['linear_192'])} / {fmt(r['linear_135'])} | {r['lnL_max_at_sigma']:g} ({r['lnL_max']:+.2f}) | {e['ess'] if e else None:}")
        lines.append("")
    (HERE / "precision_sn.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
