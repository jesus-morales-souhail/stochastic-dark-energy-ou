"""
Cuantificación del des-exprimido (decoherence del squeezing).

Pregunta:
  ¿Cuánto tarda el vacío squeezed en perder su compresión por el baño?
  ¿Cómo depende t_relax de γ y n_th?

Métrica principal (robusta frente a rotación libre de cuadraturas):
  |⟨a²⟩|(t) — correlador anómalo; decae sin oscilar por la rotación de q/p.
  t_half: primer tiempo con |⟨a²⟩|(t) = 0.5 |⟨a²⟩|(0).

También se reporta:
  t_half_var: primer cruce de Var(q) por el punto medio entre V0 y V_ss
              (puede ser temprano por rotación libre; se documenta).

Conexión cosmológica (repo stochastic-dark-energy-ou):
  θ_OU  ↔  γ_eff ≈ γ(1 + 2 n_th)   (tasa de mean-reversion / des-squeezing)
  σ_X   ↔  amplitud residual de fluctuaciones (aquí: profundidad de squeeze)
  SDiff emergente  →  t_relax finito tras perturbación del vacío
  SDiff fundamental →  σ_X = 0 exacto (sin “ruido de baño” efectivo)

Parámetros base (estrategia validada):
  N=150, r=1.5, omega=2.4
"""

from __future__ import annotations

import os
import csv
import numpy as np
import matplotlib.pyplot as plt
import qutip as qt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results", "desqueezing_relax_time")
os.makedirs(RESULTS, exist_ok=True)

# -------------------- parámetros base --------------------
N = 150
r = 1.5
omega = 2.4
t_max = 8.0
n_t = 800
tlist = np.linspace(0.0, t_max, n_t)

# Barridos
gammas = np.array([0.6, 1.2, 2.4, 3.6, 4.8])
n_ths = np.array([0.0, 0.05, 0.1, 0.3, 0.5, 1.0])
gamma_ref = 2.4  # para barrido en n_th
n_th_ref = 0.1   # para barrido en gamma


def var_q_from_expect(x, x2):
    return np.real(x2) / 2.0 - (np.real(x) ** 2) / 2.0


def first_half_time(t, y, y_half, decreasing=True):
    """Primer cruce de y por y_half."""
    y = np.asarray(y, dtype=float)
    if decreasing:
        idx = np.where(y <= y_half)[0]
    else:
        idx = np.where(y >= y_half)[0]
    if len(idx) == 0:
        return np.nan
    i = idx[0]
    if i == 0:
        return float(t[0])
    # interpolación lineal
    t0, t1 = t[i - 1], t[i]
    y0, y1 = y[i - 1], y[i]
    if y1 == y0:
        return float(t1)
    return float(t0 + (y_half - y0) * (t1 - t0) / (y1 - y0))


def run_squeezed(gamma: float, n_th: float):
    a = qt.destroy(N)
    X = a + a.dag()
    H0 = omega * a.dag() * a
    c_ops = [
        np.sqrt(gamma * (n_th + 1.0)) * a,
        np.sqrt(gamma * n_th) * a.dag() if n_th > 0 else 0 * a,
    ]
    if n_th == 0:
        c_ops = [np.sqrt(gamma) * a]

    e_ops = [a.dag() * a, X, X**2, a * a]
    psi0 = qt.squeeze(N, r) * qt.basis(N, 0)
    res = qt.mesolve(H0, psi0, tlist, c_ops, e_ops=e_ops)

    n_t_series = np.real(res.expect[0])
    var_q = var_q_from_expect(res.expect[1], res.expect[2])
    a2 = res.expect[3]
    abs_a2 = np.abs(a2)

    V0 = float(var_q[0])
    Vss = float(n_th + 0.5)
    Vmid = V0 + 0.5 * (Vss - V0)
    abs_a2_0 = float(abs_a2[0])
    half_a2 = 0.5 * abs_a2_0

    t_half_a2 = first_half_time(tlist, abs_a2, half_a2, decreasing=True)
    t_half_var = first_half_time(tlist, var_q, Vmid, decreasing=False)

    # tasa efectiva teórica (desfase / damping del correlador anómalo)
    # Para baño térmico Markoviano, |⟨a²⟩| ~ e^{-γ(1+2 n_th) t} * e^{-i 2 ω t}
    gamma_eff = gamma * (1.0 + 2.0 * n_th)
    t_half_theory = np.log(2.0) / gamma_eff if gamma_eff > 0 else np.nan

    # fit log |⟨a²⟩| en ventana temprana (antes de ruido numérico)
    mask = (tlist > 0) & (abs_a2 > 1e-6 * max(abs_a2_0, 1e-30)) & (tlist < min(3.0, t_max))
    if np.count_nonzero(mask) > 10:
        coeff = np.polyfit(tlist[mask], np.log(abs_a2[mask]), 1)
        rate_fit = -float(coeff[0])
        t_half_fit = np.log(2.0) / rate_fit if rate_fit > 0 else np.nan
    else:
        rate_fit = np.nan
        t_half_fit = np.nan

    return {
        "gamma": gamma,
        "n_th": n_th,
        "gamma_eff": gamma_eff,
        "V0": V0,
        "Vss": Vss,
        "abs_a2_0": abs_a2_0,
        "t_half_a2": t_half_a2,
        "t_half_var": t_half_var,
        "t_half_theory": t_half_theory,
        "t_half_fit": t_half_fit,
        "rate_fit": rate_fit,
        "n": n_t_series,
        "var_q": var_q,
        "abs_a2": abs_a2,
    }


def main():
    print("=" * 60)
    print("DES-SQUEEZING: t_relax(γ, n_th)")
    print(f"N={N}, r={r}, omega={omega}, t_max={t_max}")
    print("=" * 60)

    rows = []
    series_store = {}

    # --- barrido en gamma (n_th fijo) ---
    print("\n--- Barrido γ (n_th = %.3f) ---" % n_th_ref)
    for g in gammas:
        print(f"  mesolve γ={g} ...", flush=True)
        out = run_squeezed(g, n_th_ref)
        rows.append(out)
        series_store[f"g{g}_nth{n_th_ref}"] = out
        print(
            f"    t_½(|a²|)={out['t_half_a2']:.4f}  "
            f"t_½ theory={out['t_half_theory']:.4f}  "
            f"t_½ fit={out['t_half_fit']:.4f}  "
            f"γ_eff={out['gamma_eff']:.3f}"
        )

    # --- barrido en n_th (gamma fijo) ---
    print("\n--- Barrido n_th (γ = %.3f) ---" % gamma_ref)
    for nth in n_ths:
        key = f"g{gamma_ref}_nth{nth}"
        if key in series_store:
            out = series_store[key]
            print(f"  (reutilizado) n_th={nth}")
        else:
            print(f"  mesolve n_th={nth} ...", flush=True)
            out = run_squeezed(gamma_ref, nth)
            rows.append(out)
            series_store[key] = out
        print(
            f"    t_½(|a²|)={out['t_half_a2']:.4f}  "
            f"t_½ theory={out['t_half_theory']:.4f}  "
            f"γ_eff={out['gamma_eff']:.3f}"
        )

    # -------------------- CSV / NPZ --------------------
    csv_path = os.path.join(RESULTS, "t_relax_tabla.csv")
    fieldnames = [
        "gamma",
        "n_th",
        "gamma_eff",
        "V0",
        "Vss",
        "abs_a2_0",
        "t_half_a2",
        "t_half_var",
        "t_half_theory",
        "t_half_fit",
        "rate_fit",
    ]
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        for r0 in rows:
            w.writerow({k: r0[k] for k in fieldnames})

    npz_path = os.path.join(RESULTS, "t_relax_series.npz")
    np.savez(
        npz_path,
        t=tlist,
        N=N,
        r=r,
        omega=omega,
        gammas=gammas,
        n_ths=n_ths,
        gamma_ref=gamma_ref,
        n_th_ref=n_th_ref,
        **{
            f"{k}_abs_a2": v["abs_a2"]
            for k, v in series_store.items()
        },
        **{
            f"{k}_var_q": v["var_q"]
            for k, v in series_store.items()
        },
    )

    # -------------------- FIGURAS --------------------
    # 1) |⟨a²⟩|(t) para varios γ
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    ax = axes[0, 0]
    for g in gammas:
        out = series_store[f"g{g}_nth{n_th_ref}"]
        ax.semilogy(tlist, out["abs_a2"], label=f"γ={g}")
    ax.set_xlabel("t")
    ax.set_ylabel("|⟨a²⟩|(t)")
    ax.set_title(f"Des-squeezing vs γ  (n_th={n_th_ref}, r={r})")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.3)

    # 2) |⟨a²⟩|(t) para varios n_th
    ax = axes[0, 1]
    for nth in n_ths:
        out = series_store[f"g{gamma_ref}_nth{nth}"]
        ax.semilogy(tlist, out["abs_a2"], label=f"n_th={nth}")
    ax.set_xlabel("t")
    ax.set_ylabel("|⟨a²⟩|(t)")
    ax.set_title(f"Des-squeezing vs n_th  (γ={gamma_ref}, r={r})")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.3)

    # 3) t_half vs 1/γ_eff  (todos los puntos)
    ax = axes[1, 0]
    geff = np.array([r0["gamma_eff"] for r0 in rows])
    th = np.array([r0["t_half_a2"] for r0 in rows])
    th_th = np.array([r0["t_half_theory"] for r0 in rows])
    inv = 1.0 / geff
    ax.scatter(inv, th, s=50, label="simulación t_½(|a²|)", zorder=3)
    # línea teórica
    inv_line = np.linspace(inv.min() * 0.8, inv.max() * 1.1, 50)
    ax.plot(inv_line, np.log(2.0) * inv_line, "k--", label="teoría ln2 / γ_eff")
    ax.set_xlabel("1 / γ_eff   [γ_eff = γ(1+2n_th)]")
    ax.set_ylabel("t_½")
    ax.set_title(r"$t_{\rm relax} \propto 1/\gamma_{\rm eff}$")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # fit lineal t = a / γ_eff
    coeff = np.polyfit(inv, th, 1)
    slope, intercept = coeff[0], coeff[1]
    ax.plot(inv_line, slope * inv_line + intercept, "r:", alpha=0.8,
            label=f"fit: {slope:.3f}/γ_eff + {intercept:.3f}")
    ax.legend(fontsize=8)

    # 4) Var(q) ejemplo + mapa conceptual
    ax = axes[1, 1]
    out_ref = series_store[f"g{gamma_ref}_nth{n_th_ref}"]
    ax.plot(tlist, out_ref["var_q"], label="Var(q)")
    ax.axhline(out_ref["Vss"], color="gray", ls=":", label=f"V_ss={out_ref['Vss']:.2f}")
    ax.axvline(out_ref["t_half_a2"], color="C1", ls="--",
               label=f"t_½(|a²|)={out_ref['t_half_a2']:.3f}")
    ax.set_xlabel("t")
    ax.set_ylabel("Var(q)")
    ax.set_title(f"Var(q) ref (γ={gamma_ref}, n_th={n_th_ref})")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    fig_path = os.path.join(RESULTS, "desqueezing_relax_time.png")
    fig.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close(fig)

    # -------------------- RESUMEN TEXTO --------------------
    # ratio medio sim/teoría
    ratios = th / th_th
    ratio_mean = float(np.nanmean(ratios))
    ratio_std = float(np.nanstd(ratios))

    lines = [
        "=== Des-squeezing / tiempo de relajación ===",
        f"N={N}, r={r}, omega={omega}",
        f"Métrica: t_½ de |⟨a²⟩| (correlador anómalo)",
        f"Teoría Markoviana: |⟨a²⟩|(t) ∝ exp[-γ(1+2n_th) t]",
        f"  ⇒ t_½ = ln(2) / [γ(1+2n_th)] = ln2 / γ_eff",
        "",
        f"Fit lineal t_½ ≈ {slope:.4f} / γ_eff + {intercept:.4f}",
        f"  (teoría: pendiente = ln2 ≈ {np.log(2):.4f})",
        f"Ratio sim/teoría (media ± std): {ratio_mean:.4f} ± {ratio_std:.4f}",
        "",
        "--- Conexión con stochastic-dark-energy-ou ---",
        "OU: dX = -θ X dx + σ dW  →  tiempo de correlación τ_OU = 1/θ",
        "Aquí: γ_eff = γ(1+2n_th) juega el rol de tasa de restauración de lisura",
        "  tras una perturbación del vacío (squeezing = 'arranque' fuera del vacío liso).",
        "SDiff emergente ⇒ τ_relax finito (como este t_½).",
        "SDiff fundamental ⇒ σ_X=0 exacto (equivalente a no acoplar ruido local a la geometría).",
        f"Límite cosmológico memorizado: σ_X < 1.5e-4 (95% CL, DESI DR2).",
        "",
        "Tabla (gamma, n_th, gamma_eff, t_half_a2, t_half_theory):",
    ]
    for r0 in rows:
        lines.append(
            f"  γ={r0['gamma']:5.2f}  n_th={r0['n_th']:4.2f}  "
            f"γ_eff={r0['gamma_eff']:6.3f}  "
            f"t½={r0['t_half_a2']:.5f}  t½_th={r0['t_half_theory']:.5f}"
        )
    lines += [
        "",
        f"figura: {fig_path}",
        f"csv:    {csv_path}",
        f"npz:    {npz_path}",
        f"memoria repo: docs/MEMORIA_REPO_stochastic-dark-energy-ou.md",
    ]
    summary = "\n".join(lines)
    summary_path = os.path.join(RESULTS, "desqueezing_relax_time_resumen.txt")
    with open(summary_path, "w", encoding="utf-8") as f:
        f.write(summary + "\n")

    print("\n" + summary)
    print(f"\nGuardado en {RESULTS}")


if __name__ == "__main__":
    main()
