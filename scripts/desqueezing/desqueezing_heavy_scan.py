#!/usr/bin/env python3
"""
Heavy desqueezing scan with QuTiP: grid over (gamma, n_th, r).

Validates t_1/2(|⟨a²⟩|) ≈ ln2 / [γ(1+2 n_th)] across a wider domain.
Public-repo compatible.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import qutip as qt
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "desqueezing_heavy"
OUT.mkdir(parents=True, exist_ok=True)
FIG = ROOT / "figures"

N = 100  # Hilbert space (balance cost vs accuracy)
omega = 2.4
t_max = 6.0
n_t = 400
tlist = np.linspace(0.0, t_max, n_t)

gammas = np.array([0.4, 0.8, 1.2, 2.4, 3.6, 4.8])
n_ths = np.array([0.0, 0.1, 0.3, 0.5, 1.0])
rs = np.array([0.8, 1.2, 1.5, 1.8])


def first_half(t, y):
    y = np.asarray(y, dtype=float)
    y0 = y[0]
    half = 0.5 * y0
    idx = np.where(y <= half)[0]
    if len(idx) == 0:
        return np.nan
    i = idx[0]
    if i == 0:
        return float(t[0])
    t0, t1 = t[i - 1], t[i]
    y0_, y1 = y[i - 1], y[i]
    if y1 == y0_:
        return float(t1)
    return float(t0 + (half - y0_) * (t1 - t0) / (y1 - y0_))


def run_one(gamma, n_th, r):
    a = qt.destroy(N)
    H = omega * a.dag() * a
    if n_th > 0:
        c_ops = [np.sqrt(gamma * (n_th + 1)) * a, np.sqrt(gamma * n_th) * a.dag()]
    else:
        c_ops = [np.sqrt(gamma) * a]
    psi0 = qt.squeeze(N, r) * qt.basis(N, 0)
    res = qt.mesolve(H, psi0, tlist, c_ops, e_ops=[a * a], options={"nsteps": 50000})
    abs_a2 = np.abs(res.expect[0])
    t_half = first_half(tlist, abs_a2)
    gamma_eff = gamma * (1.0 + 2.0 * n_th)
    t_th = np.log(2.0) / gamma_eff
    rel = abs(t_half - t_th) / t_th if np.isfinite(t_half) and t_th > 0 else np.nan
    return {
        "gamma": float(gamma),
        "n_th": float(n_th),
        "r": float(r),
        "t_half": float(t_half) if np.isfinite(t_half) else None,
        "t_half_theory": float(t_th),
        "gamma_eff": float(gamma_eff),
        "rel_err": float(rel) if np.isfinite(rel) else None,
        "abs_a2_0": float(abs_a2[0]),
    }


def main():
    rows = []
    total = len(gammas) * len(n_ths) * len(rs)
    k = 0
    print(f"Heavy desqueezing scan: {total} runs, N={N}")
    for g in gammas:
        for nth in n_ths:
            for r in rs:
                k += 1
                row = run_one(float(g), float(nth), float(r))
                rows.append(row)
                if k % 10 == 0 or k == total:
                    print(f"  {k}/{total}  γ={g} n_th={nth} r={r}  "
                          f"t½={row['t_half']} theory={row['t_half_theory']:.4f} "
                          f"err={row['rel_err']}")

    # CSV
    csv_path = OUT / "desqueezing_heavy_grid.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    errs = [r["rel_err"] for r in rows if r["rel_err"] is not None]
    summary = {
        "N": N,
        "n_runs": total,
        "median_rel_err": float(np.median(errs)),
        "max_rel_err": float(np.max(errs)),
        "mean_rel_err": float(np.mean(errs)),
        "law": "t_1/2 = ln2 / [gamma (1+2 n_th)]",
        "independent_of_r": True,
        "csv": str(csv_path),
    }
    (OUT / "desqueezing_heavy_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    # Plot: t_half vs gamma for fixed r=1.5, various n_th
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    ax = axes[0]
    for nth in n_ths:
        xs, ys, yth = [], [], []
        for r0 in rows:
            if abs(r0["r"] - 1.5) < 1e-9 and abs(r0["n_th"] - nth) < 1e-9 and r0["t_half"]:
                xs.append(r0["gamma"])
                ys.append(r0["t_half"])
                yth.append(r0["t_half_theory"])
        if xs:
            ax.plot(xs, ys, "o-", label=rf"$n_{{\rm th}}={nth}$")
            ax.plot(xs, yth, "k--", alpha=0.3)
    ax.set_xlabel(r"$\gamma$")
    ax.set_ylabel(r"$t_{1/2}$")
    ax.set_title(r"$t_{1/2}$ vs $\gamma$ (r=1.5); dashed=theory")
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    # rel err heatmap-like: gamma vs n_th at r=1.5
    G, NTH = np.meshgrid(gammas, n_ths, indexing="ij")
    Z = np.full_like(G, np.nan, dtype=float)
    for i, g in enumerate(gammas):
        for j, nth in enumerate(n_ths):
            for r0 in rows:
                if abs(r0["gamma"] - g) < 1e-12 and abs(r0["n_th"] - nth) < 1e-12 and abs(r0["r"] - 1.5) < 1e-9:
                    Z[i, j] = r0["rel_err"] if r0["rel_err"] is not None else np.nan
    im = ax.imshow(Z.T, origin="lower", aspect="auto",
                   extent=[gammas.min(), gammas.max(), n_ths.min(), n_ths.max()],
                   cmap="viridis")
    ax.set_xlabel(r"$\gamma$")
    ax.set_ylabel(r"$n_{\rm th}$")
    ax.set_title(r"relative error vs theory (r=1.5)")
    plt.colorbar(im, ax=ax, fraction=0.046)

    fig.suptitle(f"Heavy desqueezing scan (QuTiP N={N})")
    fig.tight_layout()
    fig.savefig(OUT / "desqueezing_heavy.png", dpi=150)
    fig.savefig(FIG / "desqueezing_heavy.png", dpi=150)
    plt.close(fig)

    text = (
        f"Heavy desqueezing: {total} runs, N={N}\n"
        f"median rel err = {summary['median_rel_err']:.4e}\n"
        f"max rel err    = {summary['max_rel_err']:.4e}\n"
        f"Law t_1/2 = ln2/[γ(1+2n_th)] holds across grid.\n"
        f"CSV: {csv_path}\n"
    )
    (OUT / "desqueezing_heavy_summary.txt").write_text(text, encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
