#!/usr/bin/env python3
"""
Magnetar-scale extreme conditions as a laboratory for vacuum response,
compared with the open-system desqueezing law of this repository.

Physics used (standard astrophysics / QED; no free matching to DESI numbers):
  - Magnetar B, radius, mass, spin, giant-flare timescales from the literature
  - Schwinger critical field B_QED
  - Order-of-magnitude Euler–Heisenberg vacuum birefringence ~ alpha (B/B_QED)^2
  - Gravitational redshift for a distant observer (Schwarzschild exterior)

Open-system sector (same law as scripts/desqueezing/desqueezing_relax_time.py):
  |⟨a²⟩|(t) ∝ exp(-γ_eff t),  γ_eff = γ (1 + 2 n_th)
  t_1/2 = ln(2) / γ_eff

What this is NOT:
  - Not an identification of magnetar flares with dark-energy σ_X
  - Not a fit of DESI BAO to magnetar frequencies
  - Not a claim that QED birefringence is SDiff or unimodular gravity

What this IS:
  - A controlled comparison of timescales and dimensionless vacuum-response
    parameters across three regimes: QED/magnetar, laboratory Lindblad
    desqueezing, and the cosmological OU mean-reversion dictionary of the repo.

Author: Jesús Morales Souhail
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "magnetar_vacuum_response"
OUT.mkdir(parents=True, exist_ok=True)
FIG = ROOT / "figures"
FIG.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Constants (SI)
# ---------------------------------------------------------------------------
C = 2.99792458e8
G = 6.67430e-11
HBAR = 1.054571817e-34
E_CHARGE = 1.602176634e-19
M_E = 9.1093837015e-31
ALPHA = 7.2973525693e-3
M_SUN = 1.98847e30
PC = 3.085677581e16

# Schwinger critical field in SI: B_Q = m_e^2 c^2 / (e ħ) ≈ 4.4e9 T ≈ 4.4e13 G
B_QED_T = M_E**2 * C**2 / (E_CHARGE * HBAR)  # Tesla
B_QED_G = B_QED_T * 1e4  # Gauss


def schwarzschild_redshift(M_msun: float, R_km: float) -> dict:
    """Gravitational time-dilation factor for static observer at surface."""
    Rs = 2.0 * G * (M_msun * M_SUN) / C**2
    R = R_km * 1e3
    x = Rs / R
    if x >= 1.0:
        raise ValueError("R inside Schwarzschild radius")
    alpha = np.sqrt(1.0 - x)  # dt_proper / dt_infty at surface
    return {
        "R_s_km": Rs / 1e3,
        "compactness_Rs_over_R": x,
        "alpha_surface": alpha,  # dτ/dt_∞
        "time_dilation_infty_over_local": 1.0 / alpha,
        "redshift_z": 1.0 / alpha - 1.0,
    }


def birefringence_delta_n(B_G: float) -> float:
    """
    Order-of-magnitude vacuum birefringence (Euler–Heisenberg), low-energy limit:
      Δn ~ (α / 4π) * (B / B_QED)^2   (schematic coefficient O(α))
    Used only as a dimensionless measure of vacuum modification, not a full
    radiative-transfer calculation.
    """
    chi = B_G / B_QED_G
    return (ALPHA / (4.0 * np.pi)) * chi**2


def desqueezing_half_life(gamma: float, n_th: float = 0.0) -> float:
    """Analytic half-life of |⟨a²⟩| under thermal Lindblad damping (repo law)."""
    gamma_eff = gamma * (1.0 + 2.0 * n_th)
    if gamma_eff <= 0:
        return np.inf
    return np.log(2.0) / gamma_eff


def desqueezing_curve(t: np.ndarray, gamma: float, n_th: float = 0.0, r: float = 1.5) -> dict:
    """
    Analytic desqueezing of the anomalous correlator.
    For a squeezed vacuum under Markovian thermal damping,
      |⟨a²⟩|(t) = |⟨a²⟩|_0  exp(-γ(1+2 n_th) t)
    with |⟨a²⟩|_0 = sinh(r) cosh(r)  (vacuum squeeze; phase dropped).
    """
    gamma_eff = gamma * (1.0 + 2.0 * n_th)
    a2_0 = np.sinh(r) * np.cosh(r)
    abs_a2 = a2_0 * np.exp(-gamma_eff * t)
    return {
        "abs_a2": abs_a2,
        "abs_a2_0": a2_0,
        "gamma_eff": gamma_eff,
        "t_half": desqueezing_half_life(gamma, n_th),
    }


def magnetar_anchors() -> dict:
    """
    Literature anchors (see literature/magnetar_literature_anchors.md).
    Not free parameters of the DESI BAO analysis.
    """
    catalog = ROOT / "literature" / "magnetar_literature_anchors.json"
    if catalog.is_file():
        import json as _json
        data = _json.loads(catalog.read_text(encoding="utf-8"))
        sgr = data["objects"]["SGR_1806-20"]
        return {
            "B_typical_G": 1e14,
            "B_strong_G": float(sgr.get("B_dipole_G_order", 1e15)),
            "B_QED_G": float(data["qed"]["B_QED_G"]),
            "R_km": float(data["class_properties"]["radius_km"]["typical"]),
            "M_msun": float(data["class_properties"]["mass_Msun"]["typical"]),
            "P_spin_s": float(sgr["P_s_flare_tail"]),
            "flare_spike_s": float(sgr["giant_flare"]["spike_s"]),
            "flare_tail_s": float(sgr["giant_flare"]["tail_s_range"][0]),
            "E_iso_erg": float(sgr["giant_flare"]["E_iso_erg"]),
            "distance_SGR1806_kpc": float(sgr["distance_kpc_range"][1]),
            "source": "literature/magnetar_literature_anchors.json",
            "notes": [
                "Values from published magnetar / QED literature (McGill, Palmer 2005).",
                "Not fitted to DESI BAO.",
            ],
        }
    return {
        "B_typical_G": 1e14,
        "B_strong_G": 1e15,
        "B_QED_G": B_QED_G,
        "R_km": 12.0,
        "M_msun": 1.4,
        "P_spin_s": 7.56,
        "flare_spike_s": 0.2,
        "flare_tail_s": 400.0,
        "E_iso_erg": 2e46,
        "distance_SGR1806_kpc": 15.0,
        "source": "hardcoded fallback",
        "notes": ["Fallback if literature JSON missing."],
    }


def main() -> None:
    anchors = magnetar_anchors()
    grav = schwarzschild_redshift(anchors["M_msun"], anchors["R_km"])

    # --- Magnetar vacuum response table ---
    B_grid = np.array([1e12, 1e13, 1e14, 1e15, 2e15], dtype=float)
    rows = []
    for B in B_grid:
        chi = B / B_QED_G
        dn = birefringence_delta_n(B)
        rows.append(
            {
                "B_G": B,
                "chi_B_over_B_QED": chi,
                "delta_n_birefringence": dn,
                "QED_regime": "subcritical" if chi < 1 else "supercritical_order",
            }
        )

    # --- Desqueezing curves at magnetar-inspired timescales ---
    # Choose γ so that t_1/2 matches characteristic magnetar times (analogy only).
    timescales = {
        "spike_0.2s": anchors["flare_spike_s"],
        "spin_7.56s": anchors["P_spin_s"],
        "tail_400s": anchors["flare_tail_s"],
    }
    gamma_for_timescale = {
        name: np.log(2.0) / T for name, T in timescales.items()
    }

    t = np.linspace(0.0, 50.0, 2000)
    curves = {}
    for name, g in gamma_for_timescale.items():
        curves[name] = desqueezing_curve(t, g, n_th=0.0, r=1.5)

    # --- Cosmological dictionary from the repository ---
    # H0 ≈ 67.4 km/s/Mpc → 1/s
    H0_si = 67.4 * 1e3 / (1e6 * PC)  # s^-1
    # DESI-style OU floor θ ~ 10^{-3} and a representative θ ~ 1
    theta_floor = 1e-3
    theta_order1 = 1.0
    gamma_cosmo = {
        "theta_floor_1e-3": theta_floor * H0_si,
        "theta_O1": theta_order1 * H0_si,
    }
    t_half_cosmo = {k: np.log(2.0) / v for k, v in gamma_cosmo.items()}

    # --- Gravitational time map for a spike ---
    # Local proper duration τ → distant observer Δt_∞ = τ / α
    tau_spike = anchors["flare_spike_s"]
    t_infty_spike = tau_spike / grav["alpha_surface"]

    summary = {
        "anchors": anchors,
        "B_QED_G": B_QED_G,
        "B_QED_T": B_QED_T,
        "gravity": grav,
        "vacuum_response_vs_B": rows,
        "desqueezing_gamma_for_magnetar_timescales_s": gamma_for_timescale,
        "desqueezing_t_half_s": {k: desqueezing_half_life(g) for k, g in gamma_for_timescale.items()},
        "cosmological_gamma_si": gamma_cosmo,
        "cosmological_t_half_s": t_half_cosmo,
        "cosmological_t_half_Gyr": {k: v / (3600 * 24 * 365.25 * 1e9) for k, v in t_half_cosmo.items()},
        "flare_spike_local_s": tau_spike,
        "flare_spike_distant_observer_s": t_infty_spike,
        "interpretation": {
            "shared_idea": (
                "Both magnetar QED vacuum and open-system desqueezing exhibit a "
                "non-trivial vacuum response to extreme conditions, with a "
                "characteristic recovery timescale."
            ),
            "not_identified": (
                "Magnetar γ (if interpreted as desqueezing) is ~1/s; cosmological "
                "γ ~ θ H0 is ~1e-18/s. They are not the same physical process."
            ),
            "use_of_repo": (
                "The analytic desqueezing law t_1/2 = ln2/γ_eff is the same "
                "relation validated by QuTiP in this repository; it is applied "
                "here as a portable open-system formula, not as a DESI fit."
            ),
        },
    }

    # Write JSON
    out_json = OUT / "magnetar_vacuum_summary.json"
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, default=float)

    # Write human-readable summary
    lines = []
    lines.append("Magnetar extreme vacuum response vs repository desqueezing")
    lines.append("=" * 64)
    lines.append("")
    lines.append("1. Magnetar anchors (literature scales)")
    lines.append(f"  B_typical       = {anchors['B_typical_G']:.1e} G")
    lines.append(f"  B_strong        = {anchors['B_strong_G']:.1e} G")
    lines.append(f"  B_QED           = {B_QED_G:.3e} G")
    lines.append(f"  R, M            = {anchors['R_km']} km, {anchors['M_msun']} M_sun")
    lines.append(f"  P_spin          = {anchors['P_spin_s']} s")
    lines.append(f"  flare spike/tail= {anchors['flare_spike_s']} s / {anchors['flare_tail_s']} s")
    lines.append("")
    lines.append("2. Gravity (Schwarzschild exterior)")
    lines.append(f"  R_s / R         = {grav['compactness_Rs_over_R']:.4f}")
    lines.append(f"  α = dτ/dt_∞     = {grav['alpha_surface']:.4f}")
    lines.append(f"  1+z             = {1.0/grav['alpha_surface']:.4f}")
    lines.append(f"  spike 0.2 s local → {t_infty_spike:.3f} s at infinity")
    lines.append("")
    lines.append("3. Vacuum modification (order-of-magnitude QED)")
    lines.append(f"  {'B [G]':>12}  {'B/B_QED':>12}  {'Δn ~ α(B/B_Q)^2/4π':>20}")
    for r in rows:
        lines.append(
            f"  {r['B_G']:12.3e}  {r['chi_B_over_B_QED']:12.3e}  {r['delta_n_birefringence']:20.3e}"
        )
    lines.append("")
    lines.append("4. Open-system desqueezing (repository law)")
    lines.append("  |⟨a²⟩|(t) ∝ exp(-γ_eff t),  t_1/2 = ln2 / γ_eff")
    lines.append("  γ chosen so that t_1/2 equals a magnetar timescale (analogy only):")
    for name, g in gamma_for_timescale.items():
        lines.append(f"  {name:12s}:  γ = {g:.4e} s^-1,  t_1/2 = {desqueezing_half_life(g):.3f} s")
    lines.append("")
    lines.append("5. Cosmological OU dictionary (this repository)")
    lines.append(f"  H0             = {H0_si:.4e} s^-1")
    for k, g in gamma_cosmo.items():
        th = t_half_cosmo[k]
        lines.append(
            f"  {k:18s}: γ = θ H0 = {g:.4e} s^-1,  t_1/2 = {th:.3e} s "
            f"({th/(3600*24*365.25*1e9):.3e} Gyr)"
        )
    lines.append("")
    lines.append("6. Scale separation (the honest conclusion)")
    lines.append("  Magnetar / QED recovery, if cast as desqueezing, has γ ~ 10^{-3}–10 s^-1.")
    lines.append("  Cosmological mean-reversion γ ~ θ H0 is ~ 10^{-21}–10^{-18} s^-1.")
    lines.append("  Ratio ~ 10^{18}–10^{22}: same mathematical structure, different physics.")
    lines.append("  DESI σ_X limits are not translated into magnetar B or flare energy.")
    lines.append("")
    lines.append("Files: results/magnetar_vacuum_response/")
    text = "\n".join(lines) + "\n"
    out_txt = OUT / "magnetar_vacuum_summary.txt"
    out_txt.write_text(text, encoding="utf-8")
    print(text)

    # ---- Figures ----
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))

    # Panel 1: Δn vs B
    ax = axes[0]
    Bs = np.array([r["B_G"] for r in rows])
    dns = np.array([r["delta_n_birefringence"] for r in rows])
    ax.loglog(Bs, dns, "o-", color="steelblue", lw=2)
    ax.axvline(B_QED_G, color="gray", ls="--", label=r"$B_{\mathrm{QED}}$")
    ax.axvline(1e14, color="orange", ls=":", alpha=0.8, label=r"$10^{14}\,\mathrm{G}$")
    ax.axvline(1e15, color="crimson", ls=":", alpha=0.8, label=r"$10^{15}\,\mathrm{G}$")
    ax.set_xlabel(r"$B$ [G]")
    ax.set_ylabel(r"$\Delta n \sim \frac{\alpha}{4\pi}(B/B_{\mathrm{QED}})^2$")
    ax.set_title("Vacuum birefringence (order of magnitude)")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.3)

    # Panel 2: desqueezing curves
    ax = axes[1]
    colors = {"spike_0.2s": "crimson", "spin_7.56s": "darkorange", "tail_400s": "teal"}
    for name, cur in curves.items():
        # plot vs t/t_half so curves overlay if pure exponential
        th = cur["t_half"]
        tt = np.linspace(0, 5 * th, 500)
        yy = cur["abs_a2_0"] * np.exp(-cur["gamma_eff"] * tt)
        ax.semilogy(tt, yy / cur["abs_a2_0"], color=colors[name], lw=2, label=name)
    ax.axhline(0.5, color="k", ls="--", alpha=0.5, label="half-life")
    ax.set_xlabel(r"$t$ [s]")
    ax.set_ylabel(r"$|\langle a^2\rangle|(t)\,/\,|\langle a^2\rangle|_0$")
    ax.set_title(r"Desqueezing with $\gamma=\ln 2/t_{\mathrm{char}}$")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.3)
    ax.set_xlim(0, 40)

    # Panel 3: timescale hierarchy
    ax = axes[2]
    labels = [
        "QED Compton time",
        "flare spike 0.2 s",
        "spin 7.56 s",
        "flare tail 400 s",
        "cosmo t1/2 (theta=1)",
        "cosmo t1/2 (theta=1e-3)",
    ]
    t_qed = HBAR / (M_E * C**2)
    vals = [
        t_qed,
        anchors["flare_spike_s"],
        anchors["P_spin_s"],
        anchors["flare_tail_s"],
        t_half_cosmo["theta_O1"],
        t_half_cosmo["theta_floor_1e-3"],
    ]
    ypos = np.arange(len(vals))
    ax.barh(ypos, vals, color=["#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3", "#937860"])
    ax.set_yticks(ypos)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xscale("log")
    ax.set_xlabel("timescale [s]")
    ax.set_title("Scale separation (log s)")
    ax.grid(True, axis="x", alpha=0.3)

    fig.suptitle(
        "Magnetar vacuum response vs open-system desqueezing",
        fontsize=11,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig_path = FIG / "magnetar_vacuum_desqueezing.png"
    fig.savefig(fig_path, dpi=160, bbox_inches="tight")
    fig.savefig(OUT / "magnetar_vacuum_desqueezing.png", dpi=160, bbox_inches="tight")
    plt.close(fig)

    print(f"Wrote {out_json}")
    print(f"Wrote {out_txt}")
    print(f"Wrote {fig_path}")


if __name__ == "__main__":
    main()
