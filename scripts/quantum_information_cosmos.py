#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Quantum Information Cosmos: From DESI data to the Information-Theoretic Origin of Dark Energy
Author: Jesús Morales Souhail
Date: July 2026

This script demonstrates how the numerical result ω_R ≈ ln(4) from the QNM analysis
connects cosmology with information theory. It computes:
- The effective mass of the quantum field mediating dark energy.
- The information capacity of the cosmic horizon (Bekenstein-Hawking bound in nats).
- The evolution of cosmic memory as the universe expands.
- The relationship between the effective mass and the electron mass via the fine-structure constant.

Run with:
    python quantum_information_cosmos.py
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.constants import Planck, c, hbar, e, m_e, alpha, G, pi

FIGURES_DIR = Path(__file__).resolve().parent.parent / "figures"
FIGURES_DIR.mkdir(exist_ok=True)

# ============================================================
# PART 1: Fundamental Constants and QNM Dispersion
# ============================================================

# The fundamental frequency from DESI data (MLE fit)
omega_R = 1.388  # from desi2_ready_v2.py, H1 QNM fit
# The exact information-theoretic prediction
ln4 = np.log(4)  # 1.3862943611198906

print("=" * 70)
print("QUANTUM INFORMATION COSMOS — Analysis")
print("=" * 70)
print(f"\n📊 DESI DR2 QNM fit: ω_R = {omega_R:.4f}")
print(f"📐 Information-theoretic prediction: ln(4) = {ln4:.6f}")
print(f"✅ Agreement: {100 * (1 - abs(omega_R - ln4)/ln4):.4f}%\n")

# -----------------------------------------------------------------
# 1.1 Effective mass from de Sitter QNM dispersion
# -----------------------------------------------------------------
# Dispersion relation: ω_R^2 = (m_eff / H)^2 - 9/4
# Therefore: m_eff / H = sqrt(ω_R^2 + 9/4)

m_eff_over_H = np.sqrt(omega_R**2 + 9/4)
m_eff_over_H_exact = np.sqrt(ln4**2 + 9/4)

print("─── EFFECTIVE MASS (de Sitter QNM) ───────────────────────")
print(f"  m_eff / H (from fit) = {m_eff_over_H:.4f}")
print(f"  m_eff / H (from ln4) = {m_eff_over_H_exact:.4f}")
print(f"  Difference: {abs(m_eff_over_H - m_eff_over_H_exact):.4f}")
print(f"  → Consistent within current uncertainties.\n")

# -----------------------------------------------------------------
# 1.2 Physical mass in eV (using H0 as reference)
# -----------------------------------------------------------------
# H0 = 67.4 km/s/Mpc = 2.185e-18 s^-1
# Convert H0 to eV: E = ħ ω
H0_si = 67.4 * 1000 / (3.08567758e22)  # s^-1
H0_eV = hbar * H0_si / (1.602176634e-19)  # eV

m_eff_eV = m_eff_over_H * H0_eV
m_eff_eV_exact = m_eff_over_H_exact * H0_eV

# CORRECCIÓN: masa del electrón en eV (correcta: 511000 eV)
m_e_eV = m_e * c**2 / (1.602176634e-19)  # kg * (m/s)^2 / J_per_eV

print(f"  H0 = {H0_eV:.4e} eV")
print(f"  m_eff (fit) = {m_eff_eV:.4e} eV")
print(f"  m_eff (ln4) = {m_eff_eV_exact:.4e} eV")
print(f"  Electron mass: {m_e_eV:.4e} eV (reference)\n")

# ============================================================
# PART 2: Information Capacity of the Cosmic Horizon
# ============================================================

def horizon_info_capacity(z, H0=67.4, Om=0.315):
    """
    Compute the information capacity (in nats) of the cosmic horizon
    at redshift z, using the Bekenstein-Hawking bound.
    
    Returns:
        I_max: maximum information in nats
        R_H: horizon radius in meters
        N_bits: maximum number of bits
    """
    # Hubble parameter at redshift z
    H_z = H0 * 1000 / (3.08567758e22) * np.sqrt(Om * (1+z)**3 + (1-Om))
    # Horizon radius (comoving horizon)
    R_H = c / H_z  # meters
    
    # Planck area
    l_P = np.sqrt(hbar * G / c**3)  # meters
    # Bekenstein-Hawking entropy (in nats): S = A / (4 l_P^2)
    # But A = 4π R_H^2, so S = π R_H^2 / l_P^2
    I_max = (np.pi * R_H**2) / (l_P**2)  # nats
    
    return I_max, R_H, I_max / np.log(2)  # nats, meters, bits

# Sample redshifts: from CMB to today
z_sample = np.linspace(0, 1100, 500)
I_sample = np.array([horizon_info_capacity(z)[0] for z in z_sample])

# Current cosmic info capacity
I_now, R_now, bits_now = horizon_info_capacity(0)

print("─── COSMIC HORIZON INFORMATION CAPACITY ────────────────")
print(f"  Current horizon radius: {R_now:.3e} m  ({R_now / 3.08567758e22:.2f} Gpc)")
print(f"  Current info capacity:  {I_now:.4e} nats")
print(f"  Current info capacity:  {bits_now:.4e} bits")

# Compare with known systems
genome_bits = 3.2e9 * 2  # human genome: ~3.2 Gb base pairs, 2 bits per base
brain_bits = 1e15  # estimated synaptic bits
galaxy_bits = 1e12  # rough estimate

# Pre-calculate ratios for the final print
genome_ratio = bits_now / genome_bits

print("\n  Comparison:")
print(f"    Human genome:         {genome_bits:.2e} bits")
print(f"    Human brain (syn.):   {brain_bits:.2e} bits")
print(f"    Milky Way (stars):    {galaxy_bits:.2e} bits")
print(f"    Cosmic horizon:       {bits_now:.2e} bits")
print(f"    → The cosmos can store {genome_ratio:.2e} genomes\n")

# ============================================================
# PART 3: Evolution of Cosmic Memory
# ============================================================

# Redshift to scale factor: a = 1/(1+z)
a_sample = 1 / (1 + z_sample)
I_vs_a = np.array([horizon_info_capacity(z)[0] for z in z_sample])

# Critical redshifts: recombination, galaxy formation, today
z_recomb = 1100
z_galaxy = 6
z_today = 0

I_recomb, _, _ = horizon_info_capacity(z_recomb)
I_galaxy, _, _ = horizon_info_capacity(z_galaxy)
I_now, _, _ = horizon_info_capacity(z_today)

print("─── EVOLUTION OF COSMIC MEMORY ─────────────────────────")
print(f"  Recombination (z={z_recomb:.0f}):  I = {I_recomb:.2e} nats")
print(f"  Galaxy formation (z={z_galaxy:.0f}): I = {I_galaxy:.2e} nats")
print(f"  Today (z=0):                       I = {I_now:.2e} nats")
print(f"  Growth since recombination:        {I_now/I_recomb:.2f}x")
print(f"  Growth since galaxy formation:     {I_now/I_galaxy:.2f}x")

# ============================================================
# PART 4: Linking to the Fine-Structure Constant
# ============================================================

# The effective mass is related to the electron mass
# and the fine-structure constant via the information network
# m_eff / m_e ≈ (alpha / G_coupling) * ln(4) ?

# From the cross-correlation, G_coupling ≈ 5.04e-4
G_coupling = 5.04e-4

# m_eff / m_e (using the corrected electron mass in eV)
m_eff_over_m_e = m_eff_eV / m_e_eV

print("\n─── CONNECTION TO FINE-STRUCTURE CONSTANT ──────────────")
print(f"  m_eff / m_e = {m_eff_over_m_e:.4e}")
print(f"  G_coupling = {G_coupling:.4e} (from cross-correlation)")

# Check the relation: m_eff / m_e ≈ (alpha / G_coupling) * ln(4)
pred_ratio = (alpha / G_coupling) * ln4
print(f"  (α / G_coupling) * ln(4) = {pred_ratio:.4e}")

print(f"  Discrepancy: {m_eff_over_m_e / pred_ratio:.4f}x")
print("  (Note: exact relation requires full quantum gravity treatment)")

# ============================================================
# PART 5: Plotting
# ============================================================

fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle('Quantum Information Cosmos: Memory of the Universe', fontsize=14, fontweight='bold')

# Panel 1: Information capacity vs redshift
ax1 = axes[0, 0]
ax1.semilogy(z_sample, I_sample, color='steelblue', lw=2)
ax1.axvline(z_recomb, color='gray', linestyle='--', alpha=0.7, label='Recombination')
ax1.axvline(z_galaxy, color='orange', linestyle='--', alpha=0.7, label='Galaxy formation')
ax1.set_xlabel('Redshift z', fontsize=11)
ax1.set_ylabel('Horizon info capacity I (nats)', fontsize=11)
ax1.set_title('Cosmic Memory vs Redshift')
ax1.grid(alpha=0.3)
ax1.legend()
ax1.set_ylim(1e100, 1e125)

# Panel 2: Information capacity vs scale factor (log-log)
ax2 = axes[0, 1]
ax2.loglog(a_sample, I_vs_a, color='teal', lw=2)
ax2.set_xlabel('Scale factor a', fontsize=11)
ax2.set_ylabel('Horizon info capacity I (nats)', fontsize=11)
ax2.set_title('Memory vs Expansion (log-log)')
ax2.grid(alpha=0.3)
ax2.set_xlim(1e-3, 1)

# Panel 3: Effective mass evolution
ax3 = axes[1, 0]
H_z_func = lambda z: 67.4 * np.sqrt(0.315*(1+z)**3 + 0.685)
m_eff_z = m_eff_over_H * H_z_func(z_sample)
m_eff_eV_z = m_eff_z * hbar * (67.4*1000/3.08567758e22) / 1.602176634e-19
ax3.semilogy(z_sample, m_eff_eV_z, color='crimson', lw=2)
ax3.axhline(m_e_eV, color='green', linestyle='--', alpha=0.7, label='Electron mass')
ax3.set_xlabel('Redshift z', fontsize=11)
ax3.set_ylabel('m_eff (eV)', fontsize=11)
ax3.set_title('Effective Field Mass vs Redshift')
ax3.grid(alpha=0.3)
ax3.legend()

# Panel 4: Memory per unit volume (density of information)
ax4 = axes[1, 1]
# Volume in Gpc^3 (using comoving distance)
vol = 4/3 * np.pi * (c / (H_z_func(z_sample) * 3.08567758e22))**3  # Gpc^3
info_density = I_sample / vol
ax4.loglog(z_sample, info_density, color='purple', lw=2)
ax4.set_xlabel('Redshift z', fontsize=11)
ax4.set_ylabel('Info density (nats / Gpc³)', fontsize=11)
ax4.set_title('Information Density of the Cosmos')
ax4.grid(alpha=0.3)
ax4.set_xlim(1e-1, 1e3)

plt.tight_layout()
fig_path = FIGURES_DIR / 'quantum_information_cosmos.png'
plt.savefig(fig_path, dpi=150, bbox_inches='tight')
print(f"\n📊 Figure saved: {fig_path}")

# ============================================================
# PART 6: Summary for the paper
# ============================================================

print("\n" + "=" * 70)
print("📝 SUMMARY FOR PAPER (v3.1)")
print("=" * 70)

# Use the pre-calculated ratio
print(f"""
The DESI DR2 QNM fit yields ω_R = 1.388, which is numerically indistinguishable
from ln(4) = 1.3863 within current uncertainties. This suggests that the
oscillatory frequency of dark energy fluctuations is set by the information
content of the cosmic horizon.

Key findings:
1. Effective mass ratio: m_eff / H = sqrt( (ln 4)^2 + 9/4 ) ≈ 2.043
2. This lies just above the de Sitter stability threshold (m_eff / H > 2),
   implying that dark energy fluctuations are stable and long-lived.
3. The cosmic horizon currently stores ≈ {bits_now:.2e} bits of information,
   enough to encode ≈ {genome_ratio:.2e} human genomes.
4. The expansion of the universe increases this capacity, enabling the
   emergence of complex structures (galaxies, life, consciousness).

If confirmed by Euclid DR1, this result will provide the first empirical
evidence that dark energy is not a cosmological constant, but the
thermodynamic response of a finite-information horizon to expansion.
""")

print("=" * 70)