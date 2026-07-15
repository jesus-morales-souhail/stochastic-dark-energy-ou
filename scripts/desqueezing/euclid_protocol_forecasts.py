#!/usr/bin/env python3
"""Euclid-facing forecast grids for vacuum-relaxation protocol."""
from __future__ import annotations
import os, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, "results", "euclid_protocol")
FIG = os.path.join(ROOT, "figures")
os.makedirs(OUT, exist_ok=True)
os.makedirs(FIG, exist_ok=True)

H0, OM = 67.4, 0.315
T_H0 = 977.8 / H0
H0_GYR = 1.0 / T_H0
LN2 = np.log(2.0)
SIGMA_DESI, SIGMA_EUCLID = 1.5e-4, 1e-5
DX_DESI = 0.94
DX_EUCLID_NARROW = np.log((1 + 1.8) / (1 + 0.9))
DX_RECOMB = np.log(1 + 1090)

def residual(A0, theta, dx):
    return A0 * np.exp(-theta * dx)

def theta_min_shape(dx_spacing, rho_err, k=2.0):
    if k * rho_err >= 1:
        return np.nan
    return -np.log(1 - k * rho_err) / dx_spacing

def main():
    thetas = np.array([1e-3, 0.05, 0.1, 0.3, 0.5, 0.693, 1.0, 1.2, 2.0, 3.0, 5.0])
    A0s = np.array([1e-61, 1e-10, 1e-6, 1e-5, 3e-5, 1e-4, 1.5e-4, 5e-4, 0.018])
    rows = []
    for A0 in A0s:
        for th in thetas:
            rows.append({
                "A0": A0, "theta": th,
                "sigma_DESI_path": residual(A0, th, DX_DESI),
                "sigma_Euclid_narrow_path": residual(A0, th, DX_EUCLID_NARROW),
                "sigma_recomb_path": residual(A0, th, DX_RECOMB),
                "t_half_Gyr": LN2 / (th * H0_GYR),
                "damping_factor_DESI": float(np.exp(-th * DX_DESI)),
            })
    with open(os.path.join(OUT, "euclid_forecast_grid.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    A0_g = np.geomspace(1e-6, 1e-2, 200)
    th_g = np.geomspace(1e-3, 10, 200)
    TH, AA = np.meshgrid(th_g, A0_g)
    R = residual(AA, TH, DX_DESI)
    fig, ax = plt.subplots(figsize=(7.2, 5.5))
    cs = ax.contourf(TH, AA, np.log10(R), levels=np.linspace(-8, -2, 25), cmap="viridis")
    fig.colorbar(cs, ax=ax, label=r"$\log_{10}\sigma_{\rm res}$ (DESI path)")
    ax.contour(TH, AA, R, levels=[SIGMA_EUCLID], colors="w", linewidths=1.5, linestyles="--")
    ax.contour(TH, AA, R, levels=[SIGMA_DESI], colors="r", linewidths=1.5)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"$\theta$"); ax.set_ylabel(r"$A_0$")
    ax.set_title(r"Residual after DESI path $\Delta x=0.94$")
    ax.legend(handles=[
        Line2D([0], [0], color="w", ls="--", label=r"$\sigma=10^{-5}$"),
        Line2D([0], [0], color="r", label=r"$\sigma=1.5\times10^{-4}$"),
    ], fontsize=8, loc="lower left")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "euclid_A0_theta_plane.png"), dpi=160, bbox_inches="tight")
    fig.savefig(os.path.join(FIG, "euclid_A0_theta_plane.png"), dpi=160, bbox_inches="tight")
    plt.close()

    Ns = np.array([7, 10, 15, 20, 30, 40, 50])
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    for dx_tot, lab in [(DX_DESI, "wide Dx=0.94"), (DX_EUCLID_NARROW, "narrow Dx~0.39")]:
        thmins = [theta_min_shape(dx_tot / (N - 1), 1 / np.sqrt(max(N - 3, 1)), k=1) for N in Ns]
        ax.plot(Ns, thmins, "o-", label=lab)
    ax.axhline(LN2, color="C2", ls="--", label=r"$\theta=\ln2$")
    ax.set_xlabel("Number of BAO bins N")
    ax.set_ylabel(r"Min $\theta$ (lag-1 Pearson, $k=1$)")
    ax.set_title("Shape-based theta detectability (rough)")
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "theta_shape_detectability.png"), dpi=160, bbox_inches="tight")
    fig.savefig(os.path.join(FIG, "euclid_theta_shape_detectability.png"), dpi=160, bbox_inches="tight")
    plt.close()
    print("Wrote", OUT)

if __name__ == "__main__":
    main()
