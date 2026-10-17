#!/usr/bin/env python3
"""
Recalcula los números marcados (verificar) en analisis_cmb/AUDITORIA.md con la cota vigente
(results/profile_sigma_x/profile_sigma_x.json -> sigma_X_95CL_upper_profile = 2.633e-2, raíz del perfil desde
2026-10-17; theta_best_at_upper = 0.001, borde inferior de la malla de theta). Sin escrituras: sólo imprime. La tabla está en
verificacion_numeros.md.

Constantes: las del script del repo del que sale cada número (slip_bridge.py, gap_two_targets.py,
amplification/route1_local_causal_set_seed.py) o, si el repo no las fija, las indicadas en cada bloque.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PROF = json.loads((ROOT / "results" / "profile_sigma_x" / "profile_sigma_x.json").read_text())
SX = PROF["sigma_X_95CL_upper_profile"]  # 2.633e-2 (root); grid node: sigma_X_95CL_upper_profile_grid
TH = PROF["theta_best_at_upper"]
SX_OLD = 1.5e-4
OL0 = 0.685  # AUDITORIA.md l.19; slip_bridge.py l.42
OM0 = 0.315  # slip_bridge.py l.41


def show(label, val, extra=""):
    print(f"{label:58s} {val}" + (f"   [{extra}]" if extra else ""))


print(f"sigma_X = {SX:.6e}, theta_best(profile, sigma_X) = {TH}")

# ---------------------------------------------------------------- 1. slip (slip_bridge.py l.61-73)
print("\n# 1. slip |gamma-1| = 2 eps sigma_X rho_X/rho_m/|delta_m|, eps = delta_m = 1 (slip_bridge.py)")
MAUS, SAKR = 0.17, 0.05  # slip_bridge.py l.49-54
for z in (0.5, 1.0, 1.5):
    rr = OL0 / (OM0 * (1 + z) ** 3)
    d = 2 * SX * rr
    show(f"z={z}: |gamma-1|, Maus/|gamma-1|, Sakr/|gamma-1|", f"{d:.3e}  {MAUS / d:.2f}x  {SAKR / d:.2f}x",
         f"con 1.5e-4: {2 * SX_OLD * rr:.3e}  {MAUS / (2 * SX_OLD * rr):.0f}x  {SAKR / (2 * SX_OLD * rr):.0f}x")
for eps in (0.1, 0.01):
    show(f"z=0.5, eps={eps}", f"{2 * eps * SX * OL0 / (OM0 * 1.5**3):.3e}")
# data-pack-option0-internet.md l.83-94: |eta-1| ~ f (Omega_DE/Omega_m) sigma_X, Omega_DE/Omega_m ~ 2, f = 1
e = 1 * 2 * SX
show("data-pack: |eta-1| ~ 2 sigma_X; 0.11/.. ; 0.05/..", f"{e:.3e}  {0.11 / e:.2f}x  {0.05 / e:.2f}x",
     "Maus sigma(gamma)=0.11 como en data-pack l.93")

# ---------------------------------------------------------------- 2. gap and soft r (gap_two_targets.py)
print("\n# 2. G = sigma_tgt/sigma0, r = 0.5 ln(sigma_tgt/sigma0), sigma0 = L_P/L_H (gap_two_targets.py)")
C, MPC, LP, H0 = 299_792_458.0, 3.085_677_581e22, 1.616_255e-35, 67.4
s0 = LP / (C / (H0 * 1e3 / MPC))
show("sigma0", f"{s0:.4e}")
for tgt, name in ((1e-5, "Euclid"), (SX, "DESI"), (SX_OLD, "1.5e-4 (old)")):
    show(f"G_{name}, r_{name}", f"{tgt / s0:.3e}  r = {0.5 * math.log(tgt / s0):.2f}")
show("densidad: sigma_X/(Omega_L0 sigma0)", f"{SX / (OL0 * s0):.3e}",
     f"AUDITORIA 1.85e59 = 1.27e59 x 1.46; con 1.5e-4: {SX_OLD / (OL0 * s0):.3e}")

# ---------------------------------------------------------------- 3. N_eff and lengths (route1 script)
print("\n# 3. N_eff = 1/sigma^2; L = R_H (N_eff/N_BH)^(1/p) (route1_local_causal_set_seed.py l.35-57)")
N_BH, R_H = 1.0e122, 1.4e26
for sig, name in ((SX, "DESI"), (SX_OLD, "1.5e-4 (old)"), (1e-5, "Euclid")):
    Ne = 1 / sig**2
    Ls = [R_H * (Ne / N_BH) ** (1 / p) for p in (2, 3, 4)]
    show(f"{name}: N_eff; log10(N_BH/N_eff); L(p=2,3,4) [m]",
         f"{Ne:.3e}  {math.log10(N_BH / Ne):.1f}  " + "  ".join(f"{L:.2e}" for L in Ls))

# ---------------------------------------------------------------- 4. theta_min (euclid-protocol l.137)
print("\n# 4. theta_min = -ln(1 - s/A0)/dx, dx = 0.94 (euclid-protocol-vacuum-relaxation.md l.137-141)")
for s, A0 in ((1e-5, 3e-5), (1e-5, 1e-4), (3e-5, 1e-4), (5e-5, 2.5e-2), (5e-5, SX), (5e-5, SX_OLD)):  # 2.5e-2 = value printed at l.148
    show(f"s={s:g}, A0={A0:g}", f"{-math.log1p(-s / A0) / 0.94:.3e}")

# ---------------------------------------------------------------- 5. H.5, lambda, zeta (stochastic App. H)
print("\n# 5. H.5: |zeta V/(M V')| <~ 2.7e-6/sigma_X; lambda = 1/that; M_cut = M_Pl/|zeta|")
b = 2.7e-6 / SX
show("bound H.5", f"{b:.3e}", f"con 1.5e-4: {2.7e-6 / SX_OLD:.3e}")
show("lambda >~ (Case 1, zeta ~ 1)", f"{1 / b:.3e}", f"con 1.5e-4: {SX_OLD / 2.7e-6:.1f}")
show("M_cut/M_Pl >~ (Case 2)", f"{1 / b:.3e}")
show("D3 (sigma_X/Omega_L0 en H.2): bound, lambda", f"{2.7e-6 * OL0 / SX:.3e}  {SX / (2.7e-6 * OL0):.3e}")
show("H.4 entrada: 2 sqrt(1.3^2+0.4^2); 1.3 + eso", f"{2 * math.hypot(1.3, 0.4):.2f}e-6  {1.3 + 2 * math.hypot(1.3, 0.4):.2f}e-6")

# ---------------------------------------------------------------- 6. m_phi (stochastic l.621)
print("\n# 6. m_phi: sigma_X ~ m^2 M_Pl H0/(2 pi rho_L) -> m = sqrt(2 pi rho_L sigma_X/(M_Pl H0))")
HBAR_EVS = 6.582_119_569e-16
H0_eV = H0 * 1e3 / MPC * HBAR_EVS
MRED = 2.435_32e27  # eV, reduced Planck mass
MPL = 1.220_89e28  # eV, Planck mass
rho_L = OL0 * 3 * H0_eV**2 * MRED**2
show("H0 [eV], rho_L [eV^4]", f"{H0_eV:.4e}  {rho_L:.4e}")
for M, name in ((MRED, "M_Pl reducida"), (MPL, "M_Pl = 1.22e19 GeV")):
    m = lambda s: math.sqrt(2 * math.pi * rho_L * s / (M * H0_eV))
    show(f"m_phi [eV], {name}", f"{m(SX):.3e}", f"con 1.5e-4: {m(SX_OLD):.3e}")
show("escalado sqrt desde 9.45e-5 (si se hizo con 1.5e-4)", f"{9.45e-5 * math.sqrt(SX / SX_OLD):.3e}")

# ---------------------------------------------------------------- 7. orders of magnitude (stochastic l.372)
print("\n# 7. 'suppresses local fluctuations by N orders relative to the background density'")
show("log10(1/sigma_X) (X = dOmega_L, D1)", f"{math.log10(1 / SX):.2f}", f"con 1.5e-4: {math.log10(1 / SX_OLD):.2f}")
show("log10(Omega_L0/sigma_X) (rho_Lambda de fondo, D3)", f"{math.log10(OL0 / SX):.2f}")

# ---------------------------------------------------------------- 8. noise normalisation (stochastic l.164-176)
print("\n# 8. d(drho)/dt = -theta H drho + xi, <xi xi> = K delta(t-t'), H constante; Var(drho/rho_ref) estacionaria")
print("    Euler-Maruyama en t; se compara con sigma^2/(2 theta) (l.182) para K = sigma^2 H rho^2 y K = sigma^2 H^2 rho^2/(2 theta)")
rng = np.random.default_rng(20261017)
Hc, sig, rho = 0.7, SX * math.sqrt(2.0), 1.0
for th in (TH, 1.0, 20.0):
    dt = min(0.01 / (th * Hc), 0.1)
    nstep = int(8.0 / (th * Hc * dt))
    out = {}
    for name, K in (("sigma^2 H rho^2", sig**2 * Hc * rho**2),
                    ("sigma^2 H^2 rho^2/(2theta) (publicada)", sig**2 * Hc**2 * rho**2 / (2 * th))):
        x = np.zeros(4000)
        acc = []
        for i in range(nstep):
            x += -th * Hc * x * dt + math.sqrt(K * dt) * rng.standard_normal(x.size)
            if i > nstep // 2 and i % max(1, nstep // 200) == 0:
                acc.append(np.var(x))
        out[name] = np.mean(acc) / rho**2
    show(f"theta={th}: Var/[sigma^2/2theta] para cada K",
         "  ".join(f"{k}: {v / (sig**2 / (2 * th)):.3f}" for k, v in out.items()), f"H={Hc}, {nstep} pasos")
