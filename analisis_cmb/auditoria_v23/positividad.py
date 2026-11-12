#!/usr/bin/env python3
"""
Exploration only (POSITIVIDAD.md): how much do the exact marginal -1.92 bounds tighten if rho_DE > 0 is required at every z
(a physical prior)? Nothing in results/ or papers changes.

BAO + theta_* + omega_m + Pantheon+ (1595 rows), OU theta in {0.3, 1, 3, 10}, the sigma_X cells of marginal_guiada.py.
Flat priors on (h, omega_m, M): the truncation defines a path prior for each Omega_Lambda with its own normalisation
P(sigma, Omega_Lambda) = P_prior(rho_DE > 0 | sigma, Omega_Lambda).
  Z_pos(sigma) = Z_lin(sigma) * E_q[ w * 1(rho_DE > 0) / P(sigma, Omega_Lambda) ]
with q the linear posterior (Matheron's rule, marginal_guiada.cell), w = exp(-(chi2_exact - chi2_linear)/2), positivity and P
evaluated with each sample's own Omega_Lambda = (h^2 - omega_m - omega_r)/h^2 (h, omega_m = p0 + dp). rho_DE/rho_crit,0 =
Omega_Lambda + x(z) - x(0) on the exacto grid (0 <= z <= z*; constant beyond): positive iff min_grid(x - x(0)) > -Omega_Lambda.
P(sigma, Omega_Lambda): N_PRIOR prior OU paths per (theta, sigma); nodes on a 7-point Omega_Lambda grid (best fit of the cell,
Omega_Lambda(p0 + dp_hat), +- NSIG times its posterior error by the delta method on the dp covariance), counts of positive paths
and Monte Carlo error per node; linear interpolation between nodes; a sample outside the grid is a stop condition (no extrapolation).
Bound: drop of 1.92 from the global maximum of ln Z_pos over the grid including sigma = 0 (precision_sn.lnl_bound, interpolated in
ln sigma). Per cell: ESS of the truncated weights v = w 1(pos)/P and the number of positive samples. Sensitivity to the error of P:
the N_PRIOR minima are resampled N_BOOT times, P and the bound recomputed; 16-84 % range of the bound. Stop conditions: ESS of
the truncated part below 1 % of N in the cell of the bound, or any posterior sample outside the Omega_Lambda grid.
Prediction 3: P_prior(rho_DE > 0 | sigma = new bound of theta = 1, Omega_Lambda of the LCDM best fit) with N_PRIOR fresh paths.

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
SIGMAS = {t: [s for s in MGJ["thetas"][f"theta_{t:g}"]["sigmas"]] for t in THETAS}
N_SAMP, N_PRIOR, N_BOOT, NSIG, N_NODES, SEED, DROP = 5000, 20000, 200, 3.0, 7, 20261120, 1.92


def omega_l(h, wm):
    return (h * h - wm - OMR) / (h * h)


OL0 = float(omega_l(MG.P0[0], MG.P0[1]))


def sample_cell(theta, sigma, K, AK, Msig, A, Bm, rng):
    """Linear-posterior samples (Matheron) and exact/linear weights, as marginal_guiada.cell; returns a dict of arrays."""
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


def main() -> None:
    t_start = time.time()
    Bm, A = MG.build_B(), MG.build_A()
    out = {"n_samp": N_SAMP, "n_prior": N_PRIOR, "n_boot": N_BOOT, "nsig_grid": NSIG, "n_nodes": N_NODES, "seed": SEED, "drop": DROP,
           "omega_l_bestfit_lcdm": OL0, "p0": MG.P0.tolist(), "thetas": {}}
    stop = []
    for it, theta in enumerate(THETAS):
        key = f"theta_{theta:g}"
        K = MG.kernel(theta)
        AK = A @ K
        Msig = AK @ A.T
        cells = []
        for js, sigma in enumerate(SIGMAS[theta]):
            t0 = time.time()
            rng = np.random.default_rng([SEED, it, js])
            smp = sample_cell(theta, sigma, K, AK, Msig, A, Bm, rng)
            ok = ~smp["bad"]
            ol = omega_l(smp["h"], smp["wm"])
            pos = ok & (ol + smp["min_dX"] > 0.0)
            # Omega_Lambda grid of the cell: best fit +- NSIG * posterior error (delta method on Sp)
            h0, wm0 = MG.P0[0] + smp["dp_hat"][0], MG.P0[1] + smp["dp_hat"][1]
            ol_c = float(omega_l(h0, wm0))
            grad = np.array([2 * (wm0 + OMR) / h0**3, -1.0 / h0**2])
            ol_sd = float(np.sqrt(grad @ smp["Sp"][:2, :2] @ grad))
            nodes = ol_c + NSIG * ol_sd * np.linspace(-1.0, 1.0, N_NODES)
            z = (ol[ok] - ol_c) / ol_sd
            n_out = int(np.sum((ol[ok] < nodes[0]) | (ol[ok] > nodes[-1])))
            # prior positivity on the same grid, same X(0), same border (constant beyond z*): min(x - x0) > -Omega_Lambda
            rng_p = np.random.default_rng([SEED, it, js, 7])
            if sigma > 0:
                Pp = S.ou_paths(G, sigma, theta, N_PRIOR, rng_p)
                mp = (Pp - Pp[:, [0]]).min(axis=1)
                del Pp
            else:
                mp = np.zeros(N_PRIOR)
            counts = np.array([int(np.sum(mp > -nd)) for nd in nodes])
            P_nodes = counts / N_PRIOR
            P_err = np.sqrt(np.maximum(P_nodes * (1 - P_nodes), 1.0 / N_PRIOR) / N_PRIOR)
            P_i = np.interp(ol, nodes, P_nodes)
            P_i_direct = np.array([np.mean(mp > -o) for o in ol[ok][:200]])  # check of the interpolation on the first 200 samples
            v = np.where(pos & (P_i > 0), smp["w"] * pos / np.where(P_i > 0, P_i, 1.0), 0.0)
            ess_w = float(smp["w"].sum() ** 2 / np.sum(smp["w"] ** 2)) if smp["w"].sum() > 0 else 0.0
            ess_v = float(v.sum() ** 2 / np.sum(v**2)) if v.sum() > 0 else 0.0
            ln_mean_w = float(np.log(max(smp["w"].mean(), 1e-300)))
            ln_mean_v = float(np.log(max(v.mean(), 1e-300)))
            # bootstrap of P (resample the prior minima) -> ln E[v]
            rng_b = np.random.default_rng([SEED, it, js, 11])
            ln_v_boot = np.empty(N_BOOT)
            for b in range(N_BOOT):
                mb = mp[rng_b.integers(0, N_PRIOR, N_PRIOR)]
                Pb = np.array([np.mean(mb > -nd) for nd in nodes])
                Pib = np.interp(ol, nodes, Pb)
                vb = np.where(pos & (Pib > 0), smp["w"] * pos / np.where(Pib > 0, Pib, 1.0), 0.0)
                ln_v_boot[b] = np.log(max(vb.mean(), 1e-300))
            cell = {"sigma": sigma, "lnZ_lin": smp["lnZ_lin"], "ln_mean_w": ln_mean_w, "ln_mean_v": ln_mean_v, "ln_mean_v_boot": ln_v_boot.tolist(),
                    "ess_w": ess_w, "ess_v": ess_v, "ess_v_frac_N": ess_v / N_SAMP, "n_H2_le_0": int(smp["bad"].sum()), "n_pos": int(pos.sum()),
                    "frac_pos_of_ok": float(pos.sum() / max(ok.sum(), 1)), "omega_l": {"center": ol_c, "sd_delta": ol_sd, "sd_samples": float(np.std(ol[ok])),
                    "nodes": nodes.tolist(), "max_abs_z": float(np.max(np.abs(z))), "n_outside_grid": n_out, "frac_outside": n_out / int(ok.sum())},
                    "P_nodes": P_nodes.tolist(), "P_counts": counts.tolist(), "P_err_mc": P_err.tolist(),
                    "P_interp_vs_direct_max_abs_diff": float(np.max(np.abs(P_i[ok][:200] - P_i_direct))) if sigma > 0 else 0.0,
                    "P_at_samples_p16_p50_p84": [float(q) for q in np.percentile(P_i[ok], [16, 50, 84])],
                    "chi2_ex_minus_lin_p16_p50_p84": [float(q) for q in np.percentile((smp["chi_ex"] - smp["chi_lin"])[ok], [16, 50, 84])],
                    "seconds": time.time() - t0}
            if n_out > 0:
                stop.append(f"{key} sigma={sigma:g}: {n_out} posterior samples outside the Omega_Lambda grid (max |z| = {cell['omega_l']['max_abs_z']:.2f})")
            cells.append(cell)
            print(key, json.dumps({k: cell[k] for k in ("sigma", "ln_mean_w", "ln_mean_v", "ess_w", "ess_v", "n_pos", "n_H2_le_0")}),
                  "Omega_L:", json.dumps({k: cell["omega_l"][k] for k in ("center", "sd_delta", "max_abs_z", "n_outside_grid")}),
                  "P nodes:", [round(q, 4) for q in P_nodes], f"[{cell['seconds']:.0f} s]", flush=True)
        rec = {"sigmas": SIGMAS[theta], "cells": cells, "reference_exact_192_marginal_guiada": MGJ["thetas"][key]["exact_192"]}
        if any(s.startswith(key) for s in stop):
            rec["bounds_skipped"] = "stop condition (posterior samples outside the Omega_Lambda grid): no bound computed"
            out["thetas"][key] = rec
            print(key, "STOP: bounds skipped", flush=True)
            continue
        # profiles and bounds
        sig = np.array([c["sigma"] for c in cells])
        lz_lin = np.array([c["lnZ_lin"] for c in cells])
        lz_ex = lz_lin + np.array([c["ln_mean_w"] for c in cells])
        lz_pos = lz_lin + np.array([c["ln_mean_v"] for c in cells])
        m = sig > 0
        rec["lnZ_exact_rel"] = (lz_ex - lz_ex[0]).tolist()
        rec["lnZ_pos_rel"] = (lz_pos - lz_pos[0]).tolist()
        rec["bound_exact_192"] = PS.lnl_bound(sig[m], (lz_ex - lz_ex[0])[m], DROP)
        rec["bound_pos_192"] = PS.lnl_bound(sig[m], (lz_pos - lz_pos[0])[m], DROP)
        rec["bound_exact_135"] = PS.lnl_bound(sig[m], (lz_ex - lz_ex[0])[m], 1.35)
        rec["bound_pos_135"] = PS.lnl_bound(sig[m], (lz_pos - lz_pos[0])[m], 1.35)
        boots = np.empty(N_BOOT)
        LB = np.array([c["ln_mean_v_boot"] for c in cells])  # (cells, boot)
        for b in range(N_BOOT):
            lzb = lz_lin + LB[:, b]
            boots[b] = PS.lnl_bound(sig[m], (lzb - lzb[0])[m], DROP)
        fin = boots[np.isfinite(boots) & (boots > 0)]
        rec["bound_pos_192_boot_p16_p84"] = [float(np.percentile(fin, 16)), float(np.percentile(fin, 84))] if len(fin) else None
        rec["bound_pos_192_boot_frac_no_bound"] = float(np.mean(~(np.isfinite(boots) & (boots > 0))))
        rec["max_pos_at_sigma"] = float(np.concatenate([[0.0], sig[m]])[int(np.argmax(np.concatenate([[0.0], (lz_pos - lz_pos[0])[m]])))])
        # cell of the bound: the two nodes bracketing the crossing
        b = rec["bound_pos_192"]
        if np.isfinite(b) and b > 0:
            k = int(np.searchsorted(sig, b)) - 1
            br = [cells[k], cells[min(k + 1, len(cells) - 1)]]
            rec["bound_cells"] = [{"sigma": c["sigma"], "ess_v": c["ess_v"], "ess_v_frac_N": c["ess_v_frac_N"], "n_pos": c["n_pos"],
                                   "P_at_center": float(np.interp(c["omega_l"]["center"], c["omega_l"]["nodes"], c["P_nodes"])),
                                   "P_err_at_center": float(np.interp(c["omega_l"]["center"], c["omega_l"]["nodes"], c["P_err_mc"]))} for c in br]
            if min(c["ess_v_frac_N"] for c in br) < 0.01:
                stop.append(f"{key}: ESS of the truncated part below 1 % of N in the cell of the bound ({[c['ess_v'] for c in br]})")
        out["thetas"][key] = rec
        print(key, "bound exact", rec["bound_exact_192"], "bound pos", rec["bound_pos_192"], "boot", rec["bound_pos_192_boot_p16_p84"], flush=True)
    # prediction 3: prior positivity at the new bound of theta = 1 with the LCDM best-fit Omega_Lambda
    b1 = out["thetas"]["theta_1"].get("bound_pos_192", np.nan)
    if np.isfinite(b1) and b1 > 0:
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
    print(f"done in {out['runtime_s']:.0f} s", flush=True)


if __name__ == "__main__":
    main()
