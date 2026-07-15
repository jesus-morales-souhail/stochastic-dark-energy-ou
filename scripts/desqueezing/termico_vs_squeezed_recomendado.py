"""
Térmico vs squeezed vacuum — estrategia recomendada.

Barrido: r ∈ {1.5, 1.8} con N = 150.
Muestra cómo Var(q) del squeezed sube por disipación hacia el baño.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import qutip as qt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results", "termico_vs_squeezed_recomendado")
os.makedirs(RESULTS, exist_ok=True)

# Parámetros fijos
N = 150
omega = 2.4
gamma = 2.4
n_th = 0.1
t = np.linspace(0, 10, 500)
r_values = [1.5, 1.8]

print(f"N={N}, omega={omega}, gamma={gamma}, n_th={n_th}")
print(f"r values: {r_values}")
print(f"t: 0 → {t[-1]}, {len(t)} puntos\n")

a = qt.destroy(N)
X = a + a.dag()
H0 = omega * a.dag() * a
c_ops = [
    np.sqrt(gamma * (n_th + 1)) * a,
    np.sqrt(gamma * n_th) * a.dag(),
]
e_ops = [a.dag() * a, X, X**2]

# Estado térmico inicial (igual en ambos casos)
rho_thermal_init = qt.thermal_dm(N, 2.25)

print("Resolviendo mesolve (térmico)...")
res_thermal = qt.mesolve(H0, rho_thermal_init, t, c_ops, e_ops=e_ops)
n_thermal = np.real(res_thermal.expect[0])
var_q_thermal = (
    np.real(res_thermal.expect[2]) / 2.0
    - (np.real(res_thermal.expect[1]) ** 2) / 2.0
)

all_data = {
    "t": t,
    "n_thermal": n_thermal,
    "var_q_thermal": var_q_thermal,
    "N": N,
    "omega": omega,
    "gamma": gamma,
    "n_th": n_th,
}

summary_lines = [
    "=== Térmico vs Squeezed (estrategia recomendada) ===",
    f"N={N}, omega={omega}, gamma={gamma}, n_th={n_th}",
    f"Var(q) final térmico: {var_q_thermal[-1]:.6f}",
    f"<n> final térmico:    {n_thermal[-1]:.6f}",
    "",
]

# Colores / estilos
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
axes[0].plot(t, n_thermal, "k-", lw=1.8, label="Térmico (n̄₀=2.25)")
axes[1].plot(t, var_q_thermal, "k-", lw=1.8, label="Térmico")

for r in r_values:
    n_theory = np.sinh(r) ** 2
    var_theory = 0.5 * np.exp(-2 * r)
    print(f"--- r={r} | sinh²(r)={n_theory:.4f} | Var_th(0)={var_theory:.6f} ---")
    print(f"Resolviendo mesolve (squeezed r={r})...")

    psi_sq = qt.squeeze(N, r) * qt.basis(N, 0)
    res_sq = qt.mesolve(H0, psi_sq, t, c_ops, e_ops=e_ops)

    n_sq = np.real(res_sq.expect[0])
    var_sq = (
        np.real(res_sq.expect[2]) / 2.0
        - (np.real(res_sq.expect[1]) ** 2) / 2.0
    )
    err_rel = abs(var_sq[0] - var_theory) / var_theory

    key = f"r{r}".replace(".", "p")
    all_data[f"n_squeezed_{key}"] = n_sq
    all_data[f"var_q_squeezed_{key}"] = var_sq
    all_data[f"var_theory_{key}"] = var_theory
    all_data[f"r_{key}"] = r

    block = [
        f"--- r = {r} ---",
        f"  sinh²(r) teórico:           {n_theory:.6f}",
        f"  <n> inicial simulado:       {n_sq[0]:.6f}",
        f"  Var(q) teórica t=0:         {var_theory:.6f}",
        f"  Var(q) simulada t=0:        {var_sq[0]:.6f}",
        f"  Error relativo Var(0):      {err_rel:.4%}",
        f"  Var(q) final squeezed:      {var_sq[-1]:.6f}",
        f"  <n> final squeezed:         {n_sq[-1]:.6f}",
    ]
    summary_lines.extend(block)
    print("\n".join(block))
    print()

    axes[0].plot(t, n_sq, lw=1.6, label=f"Squeezed r={r}")
    axes[1].plot(t, var_sq, lw=1.6, label=f"Squeezed r={r}")
    axes[1].axhline(
        var_theory,
        ls="--",
        lw=1,
        alpha=0.7,
        label=f"Var teórica t=0 (r={r})",
    )

# Línea del estacionario del baño: Var(q) = n_th + 1/2
var_ss = n_th + 0.5
axes[1].axhline(var_ss, color="gray", ls=":", lw=1.2, label=f"Estacionario ({var_ss})")

axes[0].set_title("Número medio de excitaciones ⟨n⟩(t)")
axes[0].set_xlabel("Tiempo")
axes[0].set_ylabel("⟨n⟩")
axes[0].legend()
axes[0].grid(True)

axes[1].set_title("Varianza de posición Var(q)(t)")
axes[1].set_xlabel("Tiempo")
axes[1].set_ylabel("Var(q)")
axes[1].legend(fontsize=8)
axes[1].grid(True)

plt.tight_layout()
fig_path = os.path.join(RESULTS, "termico_vs_squeezed_recomendado.png")
fig.savefig(fig_path, dpi=150, bbox_inches="tight")
plt.close(fig)

data_path = os.path.join(RESULTS, "termico_vs_squeezed_recomendado_datos.npz")
np.savez(data_path, **all_data)

summary_path = os.path.join(RESULTS, "termico_vs_squeezed_recomendado_resumen.txt")
summary_lines += ["", f"figura: {fig_path}", f"datos:  {data_path}"]
summary_text = "\n".join(summary_lines)
with open(summary_path, "w", encoding="utf-8") as f:
    f.write(summary_text + "\n")

print(summary_text)
print(f"\nFigura: {fig_path}")
print(f"Datos:  {data_path}")
print(f"Resumen:{summary_path}")
