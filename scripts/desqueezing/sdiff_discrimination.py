#!/usr/bin/env python3
"""Path-integrated residual forecasts for fundamental vs emergent vacuum smoothness."""
from __future__ import annotations
import os, csv
import numpy as np
from scipy.integrate import quad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, "results", "sdiff_discrimination")
FIG = os.path.join(ROOT, "figures")
os.makedirs(OUT, exist_ok=True)
os.makedirs(FIG, exist_ok=True)

H0, OM, OL = 67.4, 0.315, 0.685
T_H0 = 977.8 / H0
H0_GYR = 1.0 / T_H0
LN2 = np.log(2.0)
SIGMA_X_DESI, SIGMA_X_EUCLID, SIGMA_0 = 2.5e-2, 1.0e-5, 1.0e-61
DX_DESI, Z_RECOMB = 0.94, 1090.0
DX_RECOMB = np.log(1.0 + Z_RECOMB)

def H_gyr(z):
    return H0_GYR * np.sqrt(OM * (1 + z) ** 3 + OL)

def age_to_z(z):
    return quad(lambda zp: 1.0 / ((1.0 + zp) * H_gyr(zp)), z, 1e5, epsabs=1e-10)[0]

def residual(A0, theta, dx):
    return A0 * np.exp(-theta * dx)

def main():
    t0 = age_to_z(0.0)
    thetas = np.array([1e-4, 1e-3, 0.01, 0.05, 0.1, 0.3, 0.693, 1.0, 1.2, 2.0, 5.0, 10.0])
    A0_scenarios = {
        "Sorkin_1e-61": 1e-61,
        "amp_1e-10": 1e-10,
        "amp_1e-5": 1e-5,
        "amp_1e-4": 1e-4,
        "old_calib_0.018": 0.018,
    }
    rows = []
    for name, A0 in A0_scenarios.items():
        for th in thetas:
            g = th * H0_GYR
            rows.append({
                "A0_name": name, "A0": A0, "theta": th,
                "gamma_per_Gyr": g,
                "t_half_Gyr": LN2 / g if g > 0 else np.nan,
                "sigma_res_DESI_path": residual(A0, th, DX_DESI),
                "sigma_res_since_recomb": residual(A0, th, DX_RECOMB),
            })
    path = os.path.join(OUT, "sdiff_discrimination_grid.csv")
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    theta_plot = np.geomspace(1e-3, 20, 200)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    ax = axes[0]
    for name, A0 in A0_scenarios.items():
        if A0 < 1e-20:
            continue
        ax.loglog(theta_plot, residual(A0, theta_plot, DX_DESI), label=name)
    ax.axhline(SIGMA_X_DESI, color="k", ls="--", lw=1, label="DESI 2.5e-2")
    ax.axhline(SIGMA_X_EUCLID, color="k", ls=":", lw=1, label="Euclid ~1e-5")
    ax.set_xlabel(r"$\theta$"); ax.set_ylabel(r"$\sigma_{\rm res}$")
    ax.set_title(rf"DESI path $\Delta x={DX_DESI}$"); ax.legend(fontsize=7)
    ax.grid(True, which="both", alpha=0.3); ax.set_ylim(1e-12, 1)
    ax = axes[1]
    for name, A0 in [("Sorkin", 1e-61), ("1e-10", 1e-10), ("1e-5", 1e-5)]:
        ax.loglog(theta_plot, residual(A0, theta_plot, DX_RECOMB), label=name)
    ax.axhline(SIGMA_X_EUCLID, color="k", ls=":"); ax.axhline(SIGMA_X_DESI, color="k", ls="--")
    ax.set_xlabel(r"$\theta$"); ax.set_ylabel(r"$\sigma_{\rm res}$")
    ax.set_title(rf"Post-recomb. $\Delta x\approx{DX_RECOMB:.2f}$"); ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "residual_vs_theta.png"), dpi=160, bbox_inches="tight")
    fig.savefig(os.path.join(FIG, "sdiff_residual_vs_theta.png"), dpi=160, bbox_inches="tight")
    plt.close()

    def A0_min(target, th, dx):
        return target * np.exp(th * dx)
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.loglog(theta_plot, A0_min(SIGMA_X_EUCLID, theta_plot, DX_DESI), label="Euclid on DESI path")
    ax.loglog(theta_plot, A0_min(SIGMA_X_DESI, theta_plot, DX_DESI), label="DESI limit on DESI path")
    ax.loglog(theta_plot, A0_min(SIGMA_X_EUCLID, theta_plot, DX_RECOMB), label="Euclid after recomb. path")
    ax.axhline(SIGMA_0, color="C3", ls="--", label=r"Sorkin $10^{-61}$")
    ax.set_xlabel(r"$\theta$"); ax.set_ylabel(r"Minimum $A_0$")
    ax.set_title("Amplification required for observability"); ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "A0_min_for_detection.png"), dpi=160, bbox_inches="tight")
    fig.savefig(os.path.join(FIG, "sdiff_A0_min_for_detection.png"), dpi=160, bbox_inches="tight")
    plt.close()

    th2 = np.geomspace(1e-3, 10, 100)
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.loglog(th2, LN2 / (th2 * H0_GYR))
    ax.axhline(t0, color="k", ls="--", label=f"Age ~ {t0:.1f} Gyr")
    ax.axhline(T_H0, color="C1", ls=":", label=f"Hubble ~ {T_H0:.1f} Gyr")
    ax.axvline(LN2, color="C2", ls="-.", label=r"$\theta=\ln 2$")
    ax.set_xlabel(r"$\theta$"); ax.set_ylabel(r"$t_{1/2}$ [Gyr]")
    ax.set_title("Physical half-life (Mapping A)"); ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "t_half_vs_theta.png"), dpi=160, bbox_inches="tight")
    fig.savefig(os.path.join(FIG, "sdiff_t_half_vs_theta.png"), dpi=160, bbox_inches="tight")
    plt.close()
    print("Wrote", OUT)

if __name__ == "__main__":
    main()
