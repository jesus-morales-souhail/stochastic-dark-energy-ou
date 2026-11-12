#!/usr/bin/env python3
"""
Exploration only (POSITIVIDAD.md, version 3, linear): how much does the LINEAR marginal -1.92 bound tighten if rho_DE > 0 is
required at every z? Nothing in results/ or papers changes.

Z_pos,lin(sigma) = Z_lin(sigma) * E_q[ 1(rho_DE > 0) / P(sigma, Omega_Lambda) ], q the linear posterior (Matheron, no exact/linear
weights), positivity with each sample's own Omega_Lambda, P(sigma, Omega_Lambda) the prior positivity at that Omega_Lambda (grid of
version 2: min-max of all posterior samples of all cells + 5 % margin, 15 nodes, interpolated). Cells of version 2 (positividad.py,
same seeds, so the same samples). Bound: drop of 1.92 from the global maximum of ln Z_pos,lin (precision_sn.lnl_bound). Bootstrap
(N_BOOT) over the posterior samples and over P (joint, and each alone). Linearity check among the positive samples in the cells
bracketing the bound of theta = 1 and 3: median and 16-84 % of the exact/linear weight. Stop: fewer than 100 positive posterior
samples in a cell bracketing a bound.

    ../venv/bin/python positividad_lin.py        # writes positividad_lin.json
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import positividad as PV  # noqa: E402  (sample_cell, grids, constants; its main() is not run)
import marginal_guiada as MG  # noqa: E402
import precision_sn as PS  # noqa: E402

S, G = PV.S, PV.G
THETAS, SIGMAS, N_SAMP, N_PRIOR, N_BOOT, N_NODES, MARGIN, SEED, DROP = PV.THETAS, PV.SIGMAS, PV.N_SAMP, PV.N_PRIOR, PV.N_BOOT, PV.N_NODES, PV.MARGIN, PV.SEED, PV.DROP
MIN_POS = 100
CMB_FINAL_FRAC_NEG = {"theta_1": 0.4492, "theta_3": 0.0828, "theta_10": 0.0}


def main() -> None:
    t_start = time.time()
    Bm, A = MG.build_B(), MG.build_A()
    cells = {}
    for it, theta in enumerate(THETAS):
        key = f"theta_{theta:g}"
        K = MG.kernel(theta)
        AK = A @ K
        Msig = AK @ A.T
        for js, sigma in enumerate(SIGMAS[theta]):
            t0 = time.time()
            rng = np.random.default_rng([SEED, it, js, int(round(1000 * sigma))])
            smp = PV.sample_cell(theta, sigma, K, AK, Msig, A, Bm, rng)
            rng_p = np.random.default_rng([SEED, it, js, int(round(1000 * sigma)), 7])
            if sigma > 0:
                Pp = S.ou_paths(G, sigma, theta, N_PRIOR, rng_p)
                mp = (Pp - Pp[:, [0]]).min(axis=1)
                del Pp
            else:
                mp = np.zeros(N_PRIOR)
            ok = ~smp["bad"]
            ol = PV.omega_l(smp["h"], smp["wm"])
            cells[(theta, sigma)] = {"lnZ_lin": smp["lnZ_lin"], "ol": ol, "ok": ok, "pos": ok & (ol + smp["min_dX"] > 0.0), "w": smp["w"], "mp": mp,
                                     "ol_c": float(PV.omega_l(MG.P0[0] + smp["dp_hat"][0], MG.P0[1] + smp["dp_hat"][1])), "seconds": time.time() - t0}
            print(f"sampled {key} sigma={sigma:g}: ok={int(ok.sum())} pos={int(cells[(theta, sigma)]['pos'].sum())} [{time.time() - t0:.0f} s]", flush=True)
    all_ol = np.concatenate([c["ol"][c["ok"]] for c in cells.values()])
    lo, hi = float(all_ol.min()), float(all_ol.max())
    nodes = np.linspace(lo - MARGIN * (hi - lo), hi + MARGIN * (hi - lo), N_NODES)
    n_out = int(np.sum((all_ol < nodes[0]) | (all_ol > nodes[-1])))
    out = {"version": 3, "method": "Z_pos,lin = Z_lin * E_q[1(rho_DE>0)/P(sigma,Omega_L)], no exact/linear weights", "n_samp": N_SAMP, "n_prior": N_PRIOR, "n_boot": N_BOOT,
           "seed": SEED, "drop": DROP, "min_pos": MIN_POS, "omega_l_grid": {"min_samples": lo, "max_samples": hi, "nodes": nodes.tolist(), "n_samples_total": int(len(all_ol)), "n_outside": n_out},
           "omega_l_bestfit_lcdm": PV.OL0, "thetas": {}, "stop": []}
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
            v = np.where(pos & (P_i > 0), 1.0 / np.where(P_i > 0, P_i, 1.0), 0.0)
            ess_v = float(v.sum() ** 2 / np.sum(v**2)) if v.sum() > 0 else 0.0
            # bootstraps of ln E[v]: joint (samples and P), samples only, P only
            rng_b = np.random.default_rng([SEED, it, js, int(round(1000 * sigma)), 13])
            lb = {"joint": np.empty(N_BOOT), "samples": np.empty(N_BOOT), "P": np.empty(N_BOOT)}
            for b in range(N_BOOT):
                idx = rng_b.integers(0, N_SAMP, N_SAMP)
                mb = mp[rng_b.integers(0, N_PRIOR, N_PRIOR)]
                Pb = np.array([np.mean(mb > -nd) for nd in nodes])
                for name, (ii, PP) in {"joint": (idx, Pb), "samples": (idx, P_nodes), "P": (np.arange(N_SAMP), Pb)}.items():
                    Pib = np.interp(ol[ii], nodes, PP)
                    vb = np.where(pos[ii] & (Pib > 0), 1.0 / np.where(Pib > 0, Pib, 1.0), 0.0)
                    lb[name][b] = np.log(max(vb.mean(), 1e-300))
            wpos = w[pos]
            rec = {"sigma": sigma, "lnZ_lin": c["lnZ_lin"], "ln_mean_v": float(np.log(max(v.mean(), 1e-300))), "ln_mean_v_boot": {k: val.tolist() for k, val in lb.items()},
                   "ess_v": ess_v, "n_ok": int(ok.sum()), "n_pos": int(pos.sum()), "frac_pos_of_ok": float(pos.sum() / max(ok.sum(), 1)), "frac_neg_of_ok": float(1 - pos.sum() / max(ok.sum(), 1)),
                   "omega_l_center": c["ol_c"], "P_nodes": P_nodes.tolist(), "P_counts": counts.tolist(), "P_err_mc": P_err.tolist(),
                   "P_at_center": float(np.interp(c["ol_c"], nodes, P_nodes)), "P_err_at_center": float(np.interp(c["ol_c"], nodes, P_err)),
                   "P_at_samples_p16_p50_p84": [float(q) for q in np.percentile(P_i[ok], [16, 50, 84])],
                   "w_exact_over_linear_positive_p16_p50_p84": [float(q) for q in np.percentile(wpos, [16, 50, 84])] if len(wpos) else None,
                   "w_exact_over_linear_all_ok_p16_p50_p84": [float(q) for q in np.percentile(w[ok], [16, 50, 84])],
                   "ess_w_positive": float(wpos.sum() ** 2 / np.sum(wpos**2)) if len(wpos) and wpos.sum() > 0 else 0.0}
            recs.append(rec)
            print(key, json.dumps({k: (round(rec[k], 4) if isinstance(rec[k], float) else rec[k]) for k in ("sigma", "ln_mean_v", "ess_v", "n_pos", "P_at_center")}), flush=True)
        sig = np.array([r["sigma"] for r in recs])
        lz_lin = np.array([r["lnZ_lin"] for r in recs])
        lz_pos = lz_lin + np.array([r["ln_mean_v"] for r in recs])
        m = sig > 0
        trec = {"sigmas": SIGMAS[theta], "cells": recs, "lnZ_lin_rel": (lz_lin - lz_lin[0]).tolist(), "lnZ_pos_lin_rel": (lz_pos - lz_pos[0]).tolist(),
                "bound_lin_192": PS.lnl_bound(sig[m], (lz_lin - lz_lin[0])[m], DROP), "bound_lin_135": PS.lnl_bound(sig[m], (lz_lin - lz_lin[0])[m], 1.35),
                "bound_pos_lin_192": PS.lnl_bound(sig[m], (lz_pos - lz_pos[0])[m], DROP), "bound_pos_lin_135": PS.lnl_bound(sig[m], (lz_pos - lz_pos[0])[m], 1.35),
                "max_lin_at_sigma": float(np.concatenate([[0.0], sig[m]])[int(np.argmax(np.concatenate([[0.0], (lz_lin - lz_lin[0])[m]])))]),
                "max_pos_at_sigma": float(np.concatenate([[0.0], sig[m]])[int(np.argmax(np.concatenate([[0.0], (lz_pos - lz_pos[0])[m]])))]),
                "reference_linear_192_marginal_guiada": MGJ_LIN.get(key), "reference_exact_192_marginal_guiada": MGJ["thetas"][key]["exact_192"]}

        def interp_ln(sig_, vals, s):
            return float(np.interp(np.log(s), np.log(sig_[m]), np.asarray(vals)[m]))

        b_lin, b_pos = trec["bound_lin_192"], trec["bound_pos_lin_192"]
        if np.isfinite(b_lin) and b_lin > 0:
            trec["frac_neg_at_lin_ceiling"] = interp_ln(sig, [r["frac_neg_of_ok"] for r in recs], b_lin)
        ex = trec["reference_exact_192_marginal_guiada"]
        trec["frac_neg_at_exact_ceiling"] = interp_ln(sig, [r["frac_neg_of_ok"] for r in recs], ex)
        trec["frac_neg_cmb_final_at_exact_ceiling"] = CMB_FINAL_FRAC_NEG.get(key)
        if np.isfinite(b_pos) and b_pos > 0:
            k = int(np.searchsorted(sig, b_pos)) - 1
            br = [recs[k], recs[min(k + 1, len(recs) - 1)]]
            f = (np.log(b_pos) - np.log(br[0]["sigma"])) / (np.log(br[1]["sigma"]) - np.log(br[0]["sigma"]))
            trec["bound_cells"] = [{"sigma": r["sigma"], "n_pos": r["n_pos"], "frac_pos_of_ok": r["frac_pos_of_ok"], "ess_v": r["ess_v"], "P_at_center": r["P_at_center"],
                                    "P_err_at_center": r["P_err_at_center"], "w_positive_p16_p50_p84": r["w_exact_over_linear_positive_p16_p50_p84"], "ess_w_positive": r["ess_w_positive"]} for r in br]
            trec["frac_pos_at_bound"] = (1 - f) * br[0]["frac_pos_of_ok"] + f * br[1]["frac_pos_of_ok"]
            trec["P_at_bound"] = (1 - f) * br[0]["P_at_center"] + f * br[1]["P_at_center"]
            trec["P_err_at_bound"] = max(br[0]["P_err_at_center"], br[1]["P_err_at_center"])
            if min(r["n_pos"] for r in br) < MIN_POS:
                out["stop"].append(f"{key}: fewer than {MIN_POS} positive posterior samples in a cell bracketing the bound: {[(r['sigma'], r['n_pos']) for r in br]}")
            trec["ratio_pos_over_lin"] = b_pos / b_lin if np.isfinite(b_lin) and b_lin > 0 else None
            boots = {}
            for name in ("joint", "samples", "P"):
                LB = np.array([r["ln_mean_v_boot"][name] for r in recs])
                bb = np.empty(N_BOOT)
                for j in range(N_BOOT):
                    lzb = lz_lin + LB[:, j]
                    bb[j] = PS.lnl_bound(sig[m], (lzb - lzb[0])[m], DROP)
                fin = bb[np.isfinite(bb) & (bb > 0)]
                boots[name] = {"p16_p84": [float(np.percentile(fin, 16)), float(np.percentile(fin, 84))] if len(fin) else None, "frac_no_bound": float(np.mean(~(np.isfinite(bb) & (bb > 0))))}
            trec["bound_pos_lin_192_boot"] = boots
        out["thetas"][key] = trec
        print(key, "lin", trec["bound_lin_192"], "pos", trec["bound_pos_lin_192"], "ratio", trec.get("ratio_pos_over_lin"), "boot joint", trec.get("bound_pos_lin_192_boot", {}).get("joint"), flush=True)
    t1 = out["thetas"]["theta_1"]
    b1 = t1.get("bound_pos_lin_192", np.nan)
    if np.isfinite(b1) and b1 > 0 and not out["stop"]:
        rng_p = np.random.default_rng([SEED, 99])
        Pp = S.ou_paths(G, b1, 1.0, N_PRIOR, rng_p)
        mp = (Pp - Pp[:, [0]]).min(axis=1)
        P3 = float(np.mean(mp > -PV.OL0))
        out["prediction_3"] = {"sigma": b1, "omega_l": PV.OL0, "P_prior_pos": P3, "err_mc": float(np.sqrt(P3 * (1 - P3) / N_PRIOR)), "n_pos": int(np.sum(mp > -PV.OL0))}
    out["runtime_s"] = time.time() - t_start
    (HERE / "positividad_lin.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")
    if out["stop"]:
        print("STOP CONDITIONS:", json.dumps(out["stop"], indent=1), flush=True)
    print(f"done in {out['runtime_s']:.0f} s", flush=True)


MGJ = json.loads((HERE / "marginal_guiada.json").read_text(encoding="utf-8"))
MGJ_LIN = {}
for _k, _v in MGJ["thetas"].items():
    for _name in ("linear_192", "lin_192", "bound_linear_192"):
        if _name in _v:
            MGJ_LIN[_k] = _v[_name]

if __name__ == "__main__":
    main()
