#!/usr/bin/env python3
"""
One OU realisation (theta = 1, sigma_X = 2.506e-2, seed 20261003) end to end through CAMB,
plus the grid-convergence test (500 vs 2000 points) for smoothing Delta = 0.05 and 0.1.

The process is padded beyond a = 1 up to ln a = 0.4 (4 x Delta_max). The same continuous
path is used for both grids: X is drawn exactly on the union of the two padded grids, then
restricted to each.

Outputs: resultados/realizacion_unica/{realizacion_theta1.png, convergencia.json, convergencia.txt}
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import ou_cmb as M
import planck_datos as P

HERE = Path(__file__).resolve().parent
OUT = HERE / "resultados" / "realizacion_unica"
OUT.mkdir(parents=True, exist_ok=True)

THETA = 1.0
NS = (500, 2000)
DELTAS = (0.05, 0.1)


def check_ar1_equals_cholesky() -> float:
    x = -M.lna_grid(500)
    z = np.random.default_rng(1).standard_normal(len(x))
    xs = np.sort(x)
    C = M.SIGMA_X**2 * np.exp(-THETA * np.abs(xs[:, None] - xs[None, :]))
    chol = np.linalg.cholesky(C) @ z
    ar = M.ou_paths(xs, THETA, M.SIGMA_X, z[None, :])[0]
    return float(np.max(np.abs(chol - ar)) / M.SIGMA_X)


def main():
    ar_vs_chol = check_ar1_equals_cholesky()
    om_l = M.omega_l0()
    tt = P.load_tt()
    err_tt = np.minimum(tt["err_lo"], tt["err_hi"])  # conservative for a grid difference
    lens = P.load_lensing()

    grids = {n: M.padded_grid(n) for n in NS}
    union = np.unique(np.concatenate([g[0] for g in grids.values()]))
    z = np.random.default_rng(M.SEED).standard_normal((1, len(union)))
    x_union = M.ou_paths(union, THETA, M.SIGMA_X, z)[0]
    lookup = dict(zip(union.tolist(), x_union))

    runs = {}
    for delta in DELTAS:
        for n in NS:
            lna, mask = grids[n]
            x = np.array([lookup[v] for v in lna.tolist()])
            xs = M.smoothing_matrix(lna, delta) @ x
            rho, w = M.rho_and_w(lna, mask, xs[None, :], om_l)
            rho, w = rho[0], w[0]
            a = np.exp(lna[mask])
            r = M.run(M.make_params(a=a, w=w))
            rho_camb, _ = r["results"].get_dark_energy_rho_w(a)
            r.update(lna=lna[mask], x=x[mask], xs=xs[mask], rho=rho, w=w,
                     pp_binned=P.bin_pp(r["Clpp"], lens),
                     rho_check=float(np.max(np.abs(np.asarray(rho_camb) / rho - 1.0))))
            runs[(delta, n)] = r
            print(f"Delta={delta} N={n}: Xs(1)={xs[mask][-1]:+.5f} 100theta*={r['thetastar100']:.6f} "
                  f"w in [{w.min():.4f},{w.max():.4f}] rho CAMB/ours-1 max={r['rho_check']:.2e}")

    conv = {}
    for delta in DELTAS:
        r1, r2 = runs[(delta, 500)], runs[(delta, 2000)]
        d_tt = np.abs(r1["DlTT"][tt["ell"]] - r2["DlTT"][tt["ell"]]) / err_tt
        d_pp = np.abs(r1["pp_binned"] - r2["pp_binned"]) / lens["err"]
        d_th = abs(r1["thetastar100"] - r2["thetastar100"]) / P.THETASTAR100_ERR
        conv[str(delta)] = {
            "TT_max_diff_over_planck_err": float(d_tt.max()),
            "TT_l_at_max": int(tt["ell"][np.argmax(d_tt)]),
            "pp_max_diff_over_planck_err": float(d_pp.max()),
            "pp_bin_at_max": int(np.argmax(d_pp)) + 1,
            "thetastar_diff_over_planck_err": float(d_th),
            "Xs_at_a1_500_vs_2000": [float(r1["xs"][-1]), float(r2["xs"][-1])],
            "pass_all_lt_1pct": bool(max(d_tt.max(), d_pp.max(), d_th) < 0.01),
        }

    lcdm = M.run(M.make_params())
    ref = runs[(0.05, 2000)]
    summary = {
        "theta": THETA, "sigma_X": M.SIGMA_X, "seed": M.SEED, "pad_lna": M.PAD,
        "sigma_X_source": "results/profile_sigma_x/profile_sigma_x.json",
        "Omega_L0": om_l,
        "ar1_vs_cholesky_max_rel": ar_vs_chol,
        "planck_files": P.FILES,
        "tt_error_convention": "min(-dDl, +dDl) per l (grid difference has no sign)",
        "convergence_500_vs_2000": conv,
        "rho_camb_vs_ours_max_rel": {f"{d}_{n}": runs[(d, n)]["rho_check"] for d in DELTAS for n in NS},
        "thetastar100_by_case": {f"{d}_{n}": runs[(d, n)]["thetastar100"] for d in DELTAS for n in NS},
        "thetastar100_Delta0.1_minus_Delta0.05_over_err_N2000":
            (runs[(0.1, 2000)]["thetastar100"] - runs[(0.05, 2000)]["thetastar100"]) / P.THETASTAR100_ERR,
        "lcdm": {"thetastar100": lcdm["thetastar100"], "cosmomc_theta100": lcdm["cosmomc_theta100"]},
    }
    (OUT / "convergencia.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    lines = [f"Realizacion unica theta={THETA} sigma_X={M.SIGMA_X:.4e} seed={M.SEED} pad={M.PAD}",
             f"AR(1) vs Cholesky max rel diff: {ar_vs_chol:.2e}"]
    for d, c in conv.items():
        lines.append(f"Delta={d}: TT {c['TT_max_diff_over_planck_err']:.3e} (l={c['TT_l_at_max']}), "
                     f"pp {c['pp_max_diff_over_planck_err']:.3e} (bin {c['pp_bin_at_max']}), "
                     f"theta* {c['thetastar_diff_over_planck_err']:.3e}  -> pass={c['pass_all_lt_1pct']}")
    lines.append(f"theta*(Delta=0.1) - theta*(Delta=0.05), N=2000: "
                 f"{summary['thetastar100_Delta0.1_minus_Delta0.05_over_err_N2000']:+.3f} sigma_Planck")
    (OUT / "convergencia.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))

    a = np.exp(ref["lna"])
    fig, ax = plt.subplots(3, 1, figsize=(7.5, 9), sharex=True)
    ax[0].plot(a, ref["x"], color="0.75", lw=0.6, label="X crudo (malla 2000)")
    ax[0].plot(a, ref["xs"], color="C0", lw=1.5, label=r"X suavizado ($\Delta=0.05$)")
    ax[0].set_ylabel(r"$X(a)=\delta\Omega_\Lambda$")
    ax[0].legend(fontsize=8)
    ax[1].plot(a, ref["rho"], color="C1")
    ax[1].axhline(1, color="k", lw=0.6)
    ax[1].set_ylabel(r"$\rho_{\rm DE}/\rho_{\Lambda 0}$")
    ax[2].plot(a, ref["w"], color="C2")
    ax[2].axhline(-1, color="k", lw=0.6)
    ax[2].set_ylabel(r"$w(a)$")
    ax[2].set_xlabel("a")
    ax[2].set_xscale("log")
    for axi in ax:
        axi.axvspan(1 / (1 + 2.330), 1 / (1 + 0.295), color="C3", alpha=0.08)
        axi.grid(alpha=0.3)
    ax[0].set_title(rf"OU $\theta={THETA}$, $\sigma_X={M.SIGMA_X:.3e}$, semilla {M.SEED}, "
                    "relleno a>1 (franja: rango z DESI)", fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / "realizacion_theta1.png", dpi=150)


if __name__ == "__main__":
    main()
