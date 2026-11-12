#!/usr/bin/env python3
"""
Exploration only (POSITIVIDAD.md, version 2): how much do the exact marginal -1.92 bounds tighten if rho_DE > 0 is required at
every z (a physical prior)? Nothing in results/ or papers changes.

BAO + theta_* + omega_m + Pantheon+ (1595 rows), OU theta in {0.3, 1, 3, 10}, the sigma_X cells of marginal_guiada.py plus
sigma_X = 0.6 for theta = 0.3. Flat priors on (h, omega_m, M): the truncation defines a path prior for each Omega_Lambda with its own
normalisation P(sigma, Omega_Lambda) = P_prior(rho_DE > 0 | sigma, Omega_Lambda).
  Z_pos(sigma) = Z_lin(sigma) * E_q[ w * 1(rho_DE > 0) / P(sigma, Omega_Lambda) ]
with q the linear posterior (Matheron's rule, marginal_guiada.cell), w = exp(-(chi2_exact - chi2_linear)/2), positivity and P
evaluated with each sample's own Omega_Lambda = (h^2 - omega_m - omega_r)/h^2 (h, omega_m = p0 + dp). rho_DE/rho_crit,0 =
Omega_Lambda + x(z) - x(0) on the exacto grid (0 <= z <= z*; constant beyond): positive iff min_grid(x - x(0)) > -Omega_Lambda.
Version 2 (design fixed before any bound was seen):
- Pass 1 samples every cell; the Omega_Lambda grid is ONE grid for all cells: from the minimum to the maximum Omega_Lambda of all
  posterior samples of all cells, plus 5 % of that range on each side, 15 nodes. Check: no sample outside (stop otherwise).
- P(sigma, Omega_Lambda): N_PRIOR prior OU paths per (theta, sigma); counts of positive paths and Monte Carlo error per node;
  linear interpolation between nodes. In the two cells bracketing each bound the interpolation is checked against the direct
  count with the same paths at 5 Omega_Lambda values between nodes (max difference reported).
- Bound: drop of 1.92 from the global maximum of ln Z_pos over the grid including sigma = 0 (precision_sn.lnl_bound, interpolated in
  ln sigma). Per cell: ESS of the truncated weights v = w 1(pos)/P and the number of positive samples. If the ESS of the truncated
  part is below 1 % of N in a cell bracketing the bound, that theta is declared "not evaluable with this method" (no bound); the
  other thetas continue.
- Sensitivity to the error of P: the N_PRIOR minima are resampled N_BOOT times, P and the bound recomputed; 16-84 % range.
- Prediction 3: P_prior(rho_DE > 0 | sigma = new bound of theta = 1, Omega_Lambda of the LCDM best fit), N_PRIOR fresh paths.

    ../venv/bin/python positividad.py        # writes positividad.json (and positividad_stop.json if a stop condition fires)
"""

from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "4")

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.linalg import cho_factor, cho_solve

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import marginal_guiada as MG  # noqa: E402
import precision_sn as PS  # noqa: E402

G, S, OMR = MG.G, MG.S, MG.OMR
MGJ = json.loads((HERE / "marginal_guiada.json").read_text(encoding="utf-8"))
THETAS = (0.3, 1.0, 3.0, 10.0)
SIGMAS = {t: sorted(set(MGJ["thetas"][f"theta_{t:g}"]["sigmas"]) | ({0.6} if t == 0.3 else set())) for t in THETAS}
N_SAMP, N_PRIOR, N_BOOT, N_NODES, MARGIN, SEED, DROP, ESS_MIN = 5000, 20000, 200, 15, 0.05, 20261120, 1.92, 0.01


def omega_l(h, wm):
    return (h * h - wm - OMR) / (h * h)


OL0 = float(omega_l(MG.P0[0], MG.P0[1]))


def sample_cell(theta, sigma, K, AK, Msig, A, Bm, rng):
    """Linear-posterior samples (Matheron) and exact/linear weights, as marginal_guiada.cell."""
    r = MG.D - MG.B0
    Sf = cho_factor(MG.C + sigma**2 * Msig, lower=True)
    lz = MG.lnZ_lin(Sf, Bm, r)
    SB = cho_solve(Sf, Bm)
    Sp = np.linalg.inv(Bm.T @ SB)
    dp_hat = Sp @ (Bm.T @ cho_solve(Sf, r))
    dp = dp_hat[None, :] + rng.standard_normal((N_SAMP, 3)) @ np.linalg.cholesky(Sp).T
    x0 = S.ou_paths(G, sigma, theta, N_SAMP, rng) if sigma > 0 else np.zeros((N_SAMP, G.N))
    n0 = rng.standard_normal((N_SAMP, MG.NROWS)) @ MG.LC.T
    R = r[None, :] - dp @ Bm.T - x0 @ A.T - n0
    x = x0 + sigma**2 * cho_solve(Sf, R.T).T @ AK if sigma > 0 else x0
    chi_lin = MG.chi2_of(MG.B0[None, :] + dp @ Bm.T + x @ A.T)
    p = MG.P0[None, :] + dp
    dX = x - x[:, [0]]
    m_ex, bad = MG.exact_model(p[:, 0], p[:, 1], p[:, 2], dX)
    chi_ex = np.where(bad, np.inf, MG.chi2_of(np.where(bad[:, None], 0.0, m_ex)))
    w = np.where(bad, 0.0, np.exp(np.minimum(-0.5 * (chi_ex - chi_lin), 700)))
    return {"lnZ_lin": lz, "dp_hat": dp_hat, "Sp": Sp, "h": p[:, 0], "wm": p[:, 1], "min_dX": dX.min(axis=1), "w": w, "bad": bad,
            "chi_lin": chi_lin, "chi_ex": chi_ex}


def P_of(mp: np.ndarray, ol):
    """Prior positivity for Omega_Lambda value(s) ol from the prior minima mp of x - x(0)."""
    ol = np.atleast_1d(ol)
    return np.array([np.mean(mp > -o) for o in ol])


def main() -> None:
    t_start = time.time()
    Bm, A = MG.build_B(), MG.build_A()
    out = {"version": 2, "n_samp": N_SAMP, "n_prior": N_PRIOR, "n_boot": N_BOOT, "n_nodes": N_NODES, "margin_frac_of_range": MARGIN, "seed": SEED,
           "drop": DROP, "ess_min_frac": ESS_MIN, "omega_l_bestfit_lcdm": OL0, "p0": MG.P0.tolist(), "sigmas": {f"theta_{t:g}": SIGMAS[t] for t in THETAS}, "thetas": {}}
    # ---- pass 1: samples and prior minima of every cell
    cells = {}
    for it, theta in enumerate(THETAS):
        key = f"theta_{theta:g}"
        K = MG.kernel(theta)
        AK = A @ K
        Msig = AK @ A.T
        for js, sigma in enumerate(SIGMAS[theta]):
            t0 = time.time()
            rng = np.random.default_rng([SEED, it, js, int(round(1000 * sigma))])
            smp = sample_cell(theta, sigma, K, AK, Msig, A, Bm, rng)
            rng_p = np.random.default_rng([SEED, it, js, int(round(1000 * sigma)), 7])
            if sigma > 0:
                Pp = S.ou_paths(G, sigma, theta, N_PRIOR, rng_p)
                mp = (Pp - Pp[:, [0]]).min(axis=1)
                del Pp
            else:
                mp = np.zeros(N_PRIOR)
            ok = ~smp["bad"]
            ol = omega_l(smp["h"], smp["wm"])
            cells[(theta, sigma)] = {"lnZ_lin": smp["lnZ_lin"], "ol": ol, "ok": ok, "pos": ok & (ol + smp["min_dX"] > 0.0), "w": smp["w"], "mp": mp,
                                     "chi_diff": (smp["chi_ex"] - smp["chi_lin"]), "dp_hat": smp["dp_hat"], "ol_c": float(omega_l(MG.P0[0] + smp["dp_hat"][0], MG.P0[1] + smp["dp_hat"][1])),
                                     "seconds": time.time() - t0}
            print(f"sampled {key} sigma={sigma:g}: ok={int(ok.sum())} pos={int(cells[(theta, sigma)]['pos'].sum())} [{time.time() - t0:.0f} s]", flush=True)
    # ---- the single Omega_Lambda grid
    all_ol = np.concatenate([c["ol"][c["ok"]] for c in cells.values()])
    lo, hi = float(all_ol.min()), float(all_ol.max())
    rng_ = hi - lo
    nodes = np.linspace(lo - MARGIN * rng_, hi + MARGIN * rng_, N_NODES)
    n_out = int(np.sum((all_ol < nodes[0]) | (all_ol > nodes[-1])))
    out["omega_l_grid"] = {"min_samples": lo, "max_samples": hi, "nodes": nodes.tolist(), "n_samples_total": int(len(all_ol)), "n_outside": n_out}
    stop = []
    if n_out > 0:
        stop.append(f"{n_out} posterior samples outside the Omega_Lambda grid")
    # ---- pass 2: P per node, truncated weights, profiles, bounds
    for it, theta in enumerate(THETAS):
        key = f"theta_{theta:g}"
        recs = []
        for js, sigma in enumerate(SIGMAS[theta]):
            c = cells[(theta, sigma)]
            mp, ol, ok, pos, w = c["mp"], c["ol"], c["ok"], c["pos"], c["w"]
            counts = np.array([int(np.sum(mp > -nd)) for nd in nodes])
            P_nodes = counts / N_PRIOR
            P_err = np.sqrt(np.maximum(P_nodes * (1 - P_nodes), 1.0 / N_PRIOR) / N_PRIOR)
            P_i = np.interp(ol, nodes, P_nodes)
            v = np.where(pos & (P_i > 0), w * pos / np.where(P_i > 0, P_i, 1.0), 0.0)
            ess_w = float(w.sum() ** 2 / np.sum(w**2)) if w.sum() > 0 else 0.0
            ess_v = float(v.sum() ** 2 / np.sum(v**2)) if v.sum() > 0 else 0.0
            rng_b = np.random.default_rng([SEED, it, js, int(round(1000 * sigma)), 11])
            ln_v_boot = np.empty(N_BOOT)
            for b in range(N_BOOT):
                mb = mp[rng_b.integers(0, N_PRIOR, N_PRIOR)]
                Pb = np.array([np.mean(mb > -nd) for nd in nodes])
                Pib = np.interp(ol, nodes, Pb)
                vb = np.where(pos & (Pib > 0), w * pos / np.where(Pib > 0, Pib, 1.0), 0.0)
                ln_v_boot[b] = np.log(max(vb.mean(), 1e-300))
            rec = {"sigma": sigma, "lnZ_lin": c["lnZ_lin"], "ln_mean_w": float(np.log(max(w.mean(), 1e-300))), "ln_mean_v": float(np.log(max(v.mean(), 1e-300))),
                   "ln_mean_v_boot": ln_v_boot.tolist(), "ess_w": ess_w, "ess_v": ess_v, "ess_v_frac_N": ess_v / N_SAMP, "n_ok": int(ok.sum()), "n_pos": int(pos.sum()),
                   "frac_pos_of_ok": float(pos.sum() / max(ok.sum(), 1)), "omega_l_center": c["ol_c"], "omega_l_samples_min_max": [float(ol[ok].min()), float(ol[ok].max())],
                   "omega_l_sd_samples": float(np.std(ol[ok])), "P_nodes": P_nodes.tolist(), "P_counts": counts.tolist(), "P_err_mc": P_err.tolist(),
                   "P_at_center": float(np.interp(c["ol_c"], nodes, P_nodes)), "P_err_at_center": float(np.interp(c["ol_c"], nodes, P_err)),
                   "P_at_samples_p16_p50_p84": [float(q) for q in np.percentile(P_i[ok], [16, 50, 84])],
                   "chi2_ex_minus_lin_p16_p50_p84": [float(q) for q in np.percentile(c["chi_diff"][ok], [16, 50, 84])], "seconds_sampling": c["seconds"]}
            recs.append(rec)
            print(key, json.dumps({k: (round(rec[k], 4) if isinstance(rec[k], float) else rec[k]) for k in ("sigma", "ln_mean_w", "ln_mean_v", "ess_w", "ess_v", "n_pos", "P_at_center")}), flush=True)
        sig = np.array([r["sigma"] for r in recs])
        lz_lin = np.array([r["lnZ_lin"] for r in recs])
        lz_ex = lz_lin + np.array([r["ln_mean_w"] for r in recs])
        lz_pos = lz_lin + np.array([r["ln_mean_v"] for r in recs])
        m = sig > 0
        trec = {"sigmas": SIGMAS[theta], "cells": recs, "reference_exact_192_marginal_guiada": MGJ["thetas"][key]["exact_192"],
                "lnZ_exact_rel": (lz_ex - lz_ex[0]).tolist(), "lnZ_pos_rel": (lz_pos - lz_pos[0]).tolist(),
                "bound_exact_192": PS.lnl_bound(sig[m], (lz_ex - lz_ex[0])[m], DROP), "bound_exact_135": PS.lnl_bound(sig[m], (lz_ex - lz_ex[0])[m], 1.35),
                "bound_pos_192": PS.lnl_bound(sig[m], (lz_pos - lz_pos[0])[m], DROP), "bound_pos_135": PS.lnl_bound(sig[m], (lz_pos - lz_pos[0])[m], 1.35),
                "max_exact_at_sigma": float(np.concatenate([[0.0], sig[m]])[int(np.argmax(np.concatenate([[0.0], (lz_ex - lz_ex[0])[m]])))]),
                "max_pos_at_sigma": float(np.concatenate([[0.0], sig[m]])[int(np.argmax(np.concatenate([[0.0], (lz_pos - lz_pos[0])[m]])))])}
        b = trec["bound_pos_192"]
        if np.isfinite(b) and b > 0:
            k = int(np.searchsorted(sig, b)) - 1
            br = [recs[k], recs[min(k + 1, len(recs) - 1)]]
            trec["bound_cells"] = []
            for r_ in br:
                cc = cells[(theta, r_["sigma"])]
                # interpolation check: 5 Omega_Lambda values between nodes around the cell centre, direct vs interpolated
                kk = int(np.clip(np.searchsorted(nodes, r_["omega_l_center"]) - 1, 0, N_NODES - 2))
                test = np.linspace(nodes[max(kk - 2, 0)], nodes[min(kk + 3, N_NODES - 1)], 7)[1:-1]
                direct = P_of(cc["mp"], test)
                interp = np.interp(test, nodes, np.array(r_["P_nodes"]))
                trec["bound_cells"].append({"sigma": r_["sigma"], "ess_v": r_["ess_v"], "ess_v_frac_N": r_["ess_v_frac_N"], "n_pos": r_["n_pos"],
                                            "P_at_center": r_["P_at_center"], "P_err_at_center": r_["P_err_at_center"],
                                            "interp_check": {"omega_l": test.tolist(), "direct": direct.tolist(), "interp": interp.tolist(), "max_abs_diff": float(np.max(np.abs(direct - interp)))}})
            trec["interp_check_max_abs_diff"] = max(bc["interp_check"]["max_abs_diff"] for bc in trec["bound_cells"])
            trec["evaluable"] = bool(min(bc["ess_v_frac_N"] for bc in trec["bound_cells"]) >= ESS_MIN)
            if not trec["evaluable"]:
                trec["not_evaluable_reason"] = f"ESS of the truncated part below {100 * ESS_MIN:g} % of N in a cell bracketing the bound: " + str([(bc["sigma"], round(bc["ess_v"], 1)) for bc in trec["bound_cells"]])
        else:
            trec["evaluable"] = False
            trec["not_evaluable_reason"] = "no crossing of the -1.92 drop inside the sigma grid"
        if trec["evaluable"]:
            boots = np.empty(N_BOOT)
            LB = np.array([r["ln_mean_v_boot"] for r in recs])
            for bb in range(N_BOOT):
                lzb = lz_lin + LB[:, bb]
                boots[bb] = PS.lnl_bound(sig[m], (lzb - lzb[0])[m], DROP)
            fin = boots[np.isfinite(boots) & (boots > 0)]
            trec["bound_pos_192_boot_p16_p84"] = [float(np.percentile(fin, 16)), float(np.percentile(fin, 84))] if len(fin) else None
            trec["bound_pos_192_boot_frac_no_bound"] = float(np.mean(~(np.isfinite(boots) & (boots > 0))))
            trec["ratio_pos_over_exact"] = trec["bound_pos_192"] / trec["bound_exact_192"] if np.isfinite(trec["bound_exact_192"]) and trec["bound_exact_192"] > 0 else None
            trec["ratio_pos_over_reference"] = trec["bound_pos_192"] / trec["reference_exact_192_marginal_guiada"]
        out["thetas"][key] = trec
        print(key, "exact", trec["bound_exact_192"], "pos", trec["bound_pos_192"], "evaluable", trec["evaluable"], trec.get("bound_pos_192_boot_p16_p84"), flush=True)
    # ---- prediction 3
    t1 = out["thetas"]["theta_1"]
    if t1["evaluable"]:
        b1 = t1["bound_pos_192"]
        rng_p = np.random.default_rng([SEED, 99])
        Pp = S.ou_paths(G, b1, 1.0, N_PRIOR, rng_p)
        mp = (Pp - Pp[:, [0]]).min(axis=1)
        P3 = float(np.mean(mp > -OL0))
        out["prediction_3"] = {"sigma": b1, "omega_l": OL0, "P_prior_pos": P3, "err_mc": float(np.sqrt(P3 * (1 - P3) / N_PRIOR)), "n_pos": int(np.sum(mp > -OL0))}
    out["stop"] = stop
    out["runtime_s"] = time.time() - t_start
    (HERE / "positividad.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")
    if stop:
        (HERE / "positividad_stop.json").write_text(json.dumps(stop, indent=1), encoding="utf-8")
        print("STOP CONDITIONS:", json.dumps(stop, indent=1), flush=True)
    elif (HERE / "positividad_stop.json").exists():
        (HERE / "positividad_stop.json").unlink()
    print(f"done in {out['runtime_s']:.0f} s", flush=True)


if __name__ == "__main__":
    main()
