#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ou_bao_stochastic_test.py
=========================
Stochastic Dark Energy Network Test: OU and QNM Kernels vs DESI BAO

Author:  Jesús Morales Souhail
ORCID:   0009-0000-7637-1818
Date:    July 2026
Version: 2.1 (corrected, reviewer-ready)

DATA SOURCE (current):
  DESI DR1 BAO, arXiv:2404.03000, Table 3
  Isotropic alpha measurements, 7 redshift bins.

  To update to DESI DR2 (arXiv:2503.14738):
    Replace z_eff, alpha, sigma arrays in the DATA SECTION below.
    Verify that S_z values match the DR2 fiducial cosmology YAML.

MODELS COMPARED:
  ΛCDM   : diagonal covariance (measurement errors only)
  H0 (OU): C_total = C_std + C_OU(theta, sigma_X)       [k=2 free params]
  H1 (QNM): C_total = C_std + C_QNM(theta, sigma_X, omega_R) [k=3 free params]

  OU kernel:  C(Δx) = σ_X² · exp(-θ · Δx)
  QNM kernel: C(Δx) = σ_X² · exp(-θ · Δx) · cos(ω_R · Δx)

  Physical interpretation of QNM:
    In de Sitter space, scalar QNM frequencies are complex:
    ω = ω_I + i·ω_R, where ω_I = θ (decay rate) and ω_R (oscillation).
    The dispersion relation: ω_R² = (m_eff/H)² - 9/4
    If m_eff/H < 3/2: ω_R = 0 → OU limit (H0).
    If m_eff/H > 3/2: ω_R > 0 → oscillatory kernel (H1).
    If θ ≈ 0 in fit: undamped oscillation (unphysical for de Sitter QNM).

USAGE:
  pip install numpy scipy matplotlib
  python ou_bao_stochastic_test.py

OUTPUTS:
  - Console: model comparison table, fit parameters, lag correlations
  - figures/test_desi_QNM.png: 4-panel diagnostic figure

REFERENCES:
  [1] DESI DR1 BAO: arXiv:2404.03000 (2024)
  [2] DESI DR2 BAO: arXiv:2503.14738 (2025)
  [3] de Sitter QNM: Beyer (2011), López-Ortega (2012)
  [4] Uhlenbeck & Ornstein, Phys. Rev. 36, 823 (1930)
"""

import numpy as np
from scipy.linalg import cholesky, solve_triangular
from scipy.optimize import minimize
from scipy.integrate import quad
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import warnings
import os

warnings.filterwarnings('ignore')

# ============================================================
# SECTION 1: DATA
# Source: DESI DR1, arXiv:2404.03000, Table 3
# Replace these arrays with DR2 values (arXiv:2503.14738) when available.
# ============================================================
DATA_SOURCE = "DESI DR1 (arXiv:2404.03000)"

z_eff  = np.array([0.295, 0.510, 0.706, 0.934, 1.321, 1.484, 2.330])
alpha  = np.array([0.9857, 0.9911, 0.9749, 0.9886, 0.9911, 1.0032, 0.9971])
sigma  = np.array([0.0093, 0.0077, 0.0067, 0.0046, 0.0071, 0.0153, 0.0082])

# Sensitivity kernel S(z) = d ln D_V(z) / d Omega_Lambda
# Computed for flat ΛCDM fiducial: Omega_m=0.315, H0=67.4 km/s/Mpc
# NOTE: must be recomputed if fiducial cosmology changes (see Appendix A of paper).
S_z = np.array([-0.284, -0.462, -0.595, -0.719, -0.870, -0.917, -1.070])

# ============================================================
# SECTION 2: BASE OU PARAMETERS (calibrated, Section 3.2 of paper)
# ============================================================
THETA_BASE = 1.2        # mean-reversion rate (calibrated)
SIGMA_OU   = 2.31e-2    # sigma_X (calibrated from DESI DR1 anchor bin)
SIGMA_X2   = SIGMA_OU**2

# Physical constraint: in de Sitter QNM, theta and omega_R are not
# independent. The dispersion relation gives:
# (m_eff/H)^2 = theta^2 + omega_R^2 + 9/4
# A fit with theta << omega_R implies near-undamped oscillation (unphysical).
THETA_MIN = 1e-3   # numerical floor to avoid exact zero (optimizer artifact)
OMEGA_R_MAX = 10.0  # upper bound for omega_R scan

# ============================================================
# SECTION 3: KERNEL FUNCTIONS
# ============================================================

def kernel_OU(delta_x, theta, sigma_X2):
    """
    Ornstein-Uhlenbeck kernel (H0).
    C(Δx) = σ_X² · exp(-θ · Δx)
    Limit: omega_R → 0 (field mass m_eff >= 3H/2 in de Sitter).
    """
    return sigma_X2 * np.exp(-theta * delta_x)


def kernel_QNM(delta_x, theta, sigma_X2, omega_R):
    """
    Damped oscillatory kernel (H1) — quasi-normal mode extension.
    C(Δx) = σ_X² · exp(-θ · Δx) · cos(ω_R · Δx)

    Physical note: theta = ω_I (imaginary part), omega_R = real part.
    For a scalar field in de Sitter: ω_R² = (m_eff/H)² - 9/4.
    If theta ≈ 0 in the MLE fit, the mode is near-undamped:
    this is technically allowed but physically unusual for de Sitter QNM
    and likely indicates degeneracy with only 7 data bins.
    """
    return sigma_X2 * np.exp(-theta * delta_x) * np.cos(omega_R * delta_x)


def build_cov_total(z_arr, sigma_arr, S_arr, theta, sigma_X2, omega_R=0.0):
    """
    Total covariance: C_total = C_inst + C_signal

    C_inst_ii   = sigma_obs(z_i)^2    [measurement noise, diagonal]
    C_signal_ij = S(z_i) * S(z_j) * K(|x_i - x_j|)  [OU/QNM signal]
    """
    n = len(z_arr)
    x = np.log(1.0 + z_arr)

    C_inst   = np.diag(sigma_arr**2)
    C_signal = np.zeros((n, n))

    for i in range(n):
        for j in range(n):
            dx = abs(x[i] - x[j])
            if omega_R == 0.0:
                K = kernel_OU(dx, theta, sigma_X2)
            else:
                K = kernel_QNM(dx, theta, sigma_X2, omega_R)
            C_signal[i, j] = S_arr[i] * S_arr[j] * K

    return C_inst + C_signal

# ============================================================
# SECTION 4: LIKELIHOOD
# ============================================================

def log_likelihood_gaussian(residuals, cov_matrix):
    """
    Multivariate Gaussian log-likelihood.
    logL = -0.5 * [r^T C^{-1} r + ln|C| + n*ln(2π)]
    Returns -inf if C is not positive definite.
    """
    n = len(residuals)
    try:
        L      = cholesky(cov_matrix, lower=True)
        y      = solve_triangular(L, residuals, lower=True)
        logdet = 2.0 * np.sum(np.log(np.diag(L)))
        return -0.5 * (np.dot(y, y) + logdet + n * np.log(2.0 * np.pi))
    except np.linalg.LinAlgError:
        return -np.inf

# ============================================================
# SECTION 5: MODEL FITTING (MLE)
# ============================================================

def fit_OU(residuals, z_arr, sigma_arr, S_arr):
    """
    H0: OU pure. Free parameters: theta, sigma_X.
    Returns dict with best-fit params, logL, AIC, BIC.
    """
    def neg_logL(params):
        theta_, log_s2 = params
        if theta_ <= THETA_MIN or log_s2 < -20.0 or log_s2 > -1.0:
            return 1e10
        C = build_cov_total(z_arr, sigma_arr, S_arr,
                            theta_, np.exp(log_s2), omega_R=0.0)
        return -log_likelihood_gaussian(residuals, C)

    x0  = [1.2, np.log(SIGMA_X2)]
    res = minimize(neg_logL, x0, method='Nelder-Mead',
                   options={'xatol': 1e-5, 'fatol': 1e-5, 'maxiter': 8000})

    theta_opt = res.x[0]
    s2_opt    = np.exp(res.x[1])
    logL      = -res.fun
    k, n      = 2, len(residuals)

    return {
        'theta': theta_opt, 'sigma_X': np.sqrt(s2_opt),
        'sigma_X2': s2_opt, 'omega_R': 0.0,
        'logL': logL, 'k': k,
        'AIC': 2*k - 2*logL,
        'BIC': k*np.log(n) - 2*logL,
    }


def fit_QNM(residuals, z_arr, sigma_arr, S_arr):
    """
    H1: QNM damped oscillatory. Free parameters: theta, sigma_X, omega_R.
    Grid search over omega_R to avoid local minima.
    Physical constraint: theta >= THETA_MIN (no exact undamped modes).
    """
    def neg_logL(params):
        theta_, log_s2, omega_R_ = params
        if (theta_ <= THETA_MIN or log_s2 < -20.0 or
                log_s2 > -1.0 or omega_R_ < 0.0 or omega_R_ > OMEGA_R_MAX):
            return 1e10
        C = build_cov_total(z_arr, sigma_arr, S_arr,
                            theta_, np.exp(log_s2), omega_R=omega_R_)
        return -log_likelihood_gaussian(residuals, C)

    best_res, best_val = None, np.inf
    for om0 in [0.5, 1.0, 1.4, 2.0, 3.0, 5.0]:
        x0  = [1.2, np.log(SIGMA_X2), om0]
        res = minimize(neg_logL, x0, method='Nelder-Mead',
                       options={'xatol': 1e-5, 'fatol': 1e-5, 'maxiter': 8000})
        if res.fun < best_val:
            best_val, best_res = res.fun, res

    theta_opt   = best_res.x[0]
    s2_opt      = np.exp(best_res.x[1])
    omega_R_opt = best_res.x[2]
    logL        = -best_res.fun
    k, n        = 3, len(residuals)

    result = {
        'theta': theta_opt, 'sigma_X': np.sqrt(s2_opt),
        'sigma_X2': s2_opt, 'omega_R': omega_R_opt,
        'logL': logL, 'k': k,
        'AIC': 2*k - 2*logL,
        'BIC': k*np.log(n) - 2*logL,
    }

    # Physical sanity check
    result['theta_near_zero'] = (theta_opt < 0.05)
    if result['theta_near_zero']:
        result['warning'] = (
            f"θ = {theta_opt:.4f} ≈ 0: near-undamped oscillation. "
            "Likely a numerical artifact with N=7 bins, not a physical QNM. "
            "The dispersion relation gives m_eff/H = "
            f"{np.sqrt(omega_R_opt**2 + 9/4):.3f}, "
            "but the near-zero decay rate is unphysical for de Sitter QNM."
        )

    # de Sitter dispersion: m_eff / H from best-fit omega_R
    result['m_eff_over_H'] = np.sqrt(omega_R_opt**2 + theta_opt**2 + 9.0/4.0)

    return result

# ============================================================
# SECTION 6: VALIDATION — Linearity of Δη vs Δx
# ============================================================

def test_linearity_conformal_time(z_arr, H0=67.4, Om=0.315):
    """
    Validates that conformal time η is approximately linear in x = ln(1+z)
    over the survey redshift range. This justifies using Δx as the lag
    variable in the OU/QNM kernel instead of Δη.

    r > 0.999: projection C(Δη) ≈ C(Δx) is valid.
    r ∈ [0.99, 0.999]: acceptable, include 2nd-order correction in v3.1.
    r < 0.99: projection invalid, requires explicit Jacobian.

    Note: η is computed from z=0 (today), not from the CMB.
    This is appropriate for the DESI observational range (z < 2.5).
    """
    c_kms = 299792.458

    def H_z(z_):
        return H0 * np.sqrt(Om * (1.0 + z_)**3 + (1.0 - Om))

    eta_arr, x_arr = [], []
    for z in z_arr:
        eta_val, _ = quad(lambda zp: c_kms / H_z(zp), 0.0, z)
        eta_arr.append(eta_val)
        x_arr.append(np.log(1.0 + z))

    eta_arr = np.array(eta_arr)
    x_arr   = np.array(x_arr)

    # All pairwise differences
    pairs_deta, pairs_dx = [], []
    for i in range(len(z_arr)):
        for j in range(i + 1, len(z_arr)):
            pairs_deta.append(abs(eta_arr[i] - eta_arr[j]))
            pairs_dx.append(abs(x_arr[i]   - x_arr[j]))

    pairs_deta = np.array(pairs_deta)
    pairs_dx   = np.array(pairs_dx)
    r, _       = pearsonr(pairs_dx, pairs_deta)

    return r, pairs_deta, pairs_dx, x_arr, eta_arr

# ============================================================
# SECTION 7: LAG CORRELATIONS (whitened residuals)
# ============================================================

def lag_correlation_whitened(residuals, cov_matrix, lag_k):
    """
    Pearson correlation at lag k in whitened residuals y = L^{-1} r,
    where C = L L^T (Cholesky).

    Under H_null (ΛCDM): all lags should be ~0.
    Under H_OU:  lags decay as exp(-theta * mean_Δx_k).
    Under H_QNM: lags follow exp(-theta * mean_Δx_k) * cos(omega_R * mean_Δx_k).

    With N=7 bins, σ(ρ) ≈ 1/sqrt(N-3) ≈ 0.5, so 95% CI ≈ ±1.0.
    No individual lag is significant at N=7.
    """
    if lag_k >= len(residuals):
        return np.nan
    try:
        L = cholesky(cov_matrix, lower=True)
        y = solve_triangular(L, residuals, lower=True)
        r, _ = pearsonr(y[:-lag_k], y[lag_k:])
        return r
    except np.linalg.LinAlgError:
        return np.nan

# ============================================================
# SECTION 8: MAIN COMPUTATION
# ============================================================

residuals = alpha - 1.0

# ΛCDM baseline (diagonal covariance)
C_LCDM    = np.diag(sigma**2)
logL_LCDM = log_likelihood_gaussian(residuals, C_LCDM)

# OU calibrated (fixed parameters from paper Section 3.2)
C_OU_cal    = build_cov_total(z_eff, sigma, S_z, THETA_BASE, SIGMA_X2, 0.0)
logL_OU_cal = log_likelihood_gaussian(residuals, C_OU_cal)

# H0: OU free MLE
result_OU  = fit_OU(residuals, z_eff, sigma, S_z)
C_OU_fit   = build_cov_total(z_eff, sigma, S_z,
                              result_OU['theta'], result_OU['sigma_X2'], 0.0)

# H1: QNM free MLE
result_QNM = fit_QNM(residuals, z_eff, sigma, S_z)
C_QNM_fit  = build_cov_total(z_eff, sigma, S_z,
                              result_QNM['theta'], result_QNM['sigma_X2'],
                              result_QNM['omega_R'])

# Linearity test
r_lin, pairs_deta, pairs_dx, x_arr, eta_arr = \
    test_linearity_conformal_time(z_eff)

# omega_R scan (diagnostic: ΔlogL vs omega_R with calibrated theta, sigma)
omega_R_grid = np.linspace(0.0, OMEGA_R_MAX, 100)
logL_scan    = np.array([
    log_likelihood_gaussian(
        residuals,
        build_cov_total(z_eff, sigma, S_z, THETA_BASE, SIGMA_X2, om))
    for om in omega_R_grid
])

# ============================================================
# SECTION 9: PRINT RESULTS
# ============================================================

SEP = "=" * 70

print(SEP)
print("STOCHASTIC DARK ENERGY — OU + QNM KERNEL TEST")
print(SEP)
print(f"  Data source : {DATA_SOURCE}")
print(f"  N bins      : {len(z_eff)}")
print(f"  z range     : [{z_eff.min():.3f}, {z_eff.max():.3f}]")
print(f"  WARNING: Results are PRELIMINARY with N=7 bins.")
print(f"  95% CI on lag correlations: ≈ ±1.0 (non-significant).")
print()

print("─── MODEL COMPARISON ─────────────────────────────────────────")
print(f"{'Model':<30} {'logL':>8} {'ΔlogL':>8} {'k':>3} {'AIC':>9} {'BIC':>9}")
print("-" * 70)
print(f"{'ΛCDM (baseline)':<30} {logL_LCDM:>8.3f} {'0.000':>8} {'0':>3} {'ref':>9} {'ref':>9}")
print(f"{'OU calibrated (θ=1.2)':<30} {logL_OU_cal:>8.3f} "
      f"{logL_OU_cal-logL_LCDM:>+8.3f} {'—':>3} {'—':>9} {'—':>9}")
print(f"{'H0: OU free MLE':<30} {result_OU['logL']:>8.3f} "
      f"{result_OU['logL']-logL_LCDM:>+8.3f} {result_OU['k']:>3} "
      f"{result_OU['AIC']:>9.3f} {result_OU['BIC']:>9.3f}")
print(f"{'H1: QNM free MLE':<30} {result_QNM['logL']:>8.3f} "
      f"{result_QNM['logL']-logL_LCDM:>+8.3f} {result_QNM['k']:>3} "
      f"{result_QNM['AIC']:>9.3f} {result_QNM['BIC']:>9.3f}")
print()

print("─── BEST-FIT PARAMETERS ──────────────────────────────────────")
print(f"  H0 (OU):  θ = {result_OU['theta']:.4f},  "
      f"σ_X = {result_OU['sigma_X']:.5f},  ω_R = 0 (fixed)")
print(f"  H1 (QNM): θ = {result_QNM['theta']:.4f},  "
      f"σ_X = {result_QNM['sigma_X']:.5f},  "
      f"ω_R = {result_QNM['omega_R']:.4f}")
print(f"  H1 m_eff/H (de Sitter dispersion) = {result_QNM['m_eff_over_H']:.4f}")
print()

# Physical warning for near-zero theta
if result_QNM.get('theta_near_zero'):
    print(f"  ⚠ PHYSICAL WARNING: {result_QNM['warning']}")
    print()

print("─── MODEL SELECTION (AIC/BIC) ────────────────────────────────")
dAIC = result_OU['AIC'] - result_QNM['AIC']
dBIC = result_OU['BIC'] - result_QNM['BIC']
winner_AIC = "H1 preferred" if dAIC > 0 else "H0 preferred"
winner_BIC = "H1 preferred" if dBIC > 0 else "H0 preferred"
print(f"  ΔAIC(H0 − H1) = {dAIC:+.3f}  → {winner_AIC}  (|Δ| > 2 notable, > 6 strong)")
print(f"  ΔBIC(H0 − H1) = {dBIC:+.3f}  → {winner_BIC}  (|Δ| > 2 positive, > 6 strong)")
print()
print("  ⚠ With N=7 bins, model selection is INDICATIVE only.")
print("  ⚠ Decisive test: >20 bins (Euclid DR1, expected H2 2026).")
print()

print("─── QNM CONSISTENCY: θ/ω_R RATIO ────────────────────────────")
if result_QNM['omega_R'] > 0.05:
    ratio = result_QNM['theta'] / result_QNM['omega_R']
    print(f"  θ/ω_R = {result_QNM['theta']:.4f} / "
          f"{result_QNM['omega_R']:.4f} = {ratio:.4f}")
    print(f"  (This ratio = ω_I/ω_R constrains m_eff/H in de Sitter QNM.)")
    print(f"  Falsification criterion F6: if this ratio is inconsistent")
    print(f"  between DESI DR1, DR2, and Euclid → QNM kernel incoherent.")
else:
    print(f"  ω_R ≈ {result_QNM['omega_R']:.4f}: fit converged near OU limit.")
    print(f"  H1 is numerically degenerate with H0 at this precision.")
print()

print("─── LINEARITY TEST: Δη vs Δx ────────────────────────────────")
print(f"  Pearson r(Δη, Δx) over all {len(pairs_dx)} pairs = {r_lin:.6f}")
if r_lin > 0.999:
    print(f"  ✓ Excellent linearity. Projection C(Δη) ≈ C(Δx) is VALID "
          f"for z ∈ [{z_eff.min():.2f}, {z_eff.max():.2f}].")
elif r_lin > 0.99:
    print(f"  ≈ Good linearity (r < 0.999). Consider 2nd-order correction "
          f"in future versions.")
else:
    print(f"  ✗ Non-linear. Projection C(Δη) → C(Δx) requires "
          f"explicit Jacobian. Results may be biased.")
print()

print("─── LAG CORRELATIONS (whitened residuals) ────────────────────")
print(f"  Note: 95% CI ≈ ±{1.96/np.sqrt(len(z_eff)-3):.2f} with N={len(z_eff)} bins. "
      f"No lag is individually significant.")
print()
x_bins = np.log(1.0 + z_eff)
header = f"{'Lag':>4}  {'Obs(ΛCDM)':>11} {'Obs(OU)':>11} {'Obs(QNM)':>11}  " \
         f"{'Pred_OU':>10} {'Pred_QNM':>10}  {'⟨Δx⟩':>7}"
print(f"  {header}")
print(f"  {'-'*75}")
for lag_k in [1, 2, 3]:
    rho_lcdm = lag_correlation_whitened(residuals, C_LCDM,    lag_k)
    rho_ou   = lag_correlation_whitened(residuals, C_OU_fit,  lag_k)
    rho_qnm  = lag_correlation_whitened(residuals, C_QNM_fit, lag_k)
    mean_dx  = np.mean([abs(x_bins[i+lag_k] - x_bins[i])
                        for i in range(len(z_eff) - lag_k)])
    pred_ou  = np.exp(-result_OU['theta']  * mean_dx)
    pred_qnm = (np.exp(-result_QNM['theta'] * mean_dx)
                * np.cos(result_QNM['omega_R'] * mean_dx))
    print(f"  {lag_k:>4}  {rho_lcdm:>+11.4f} {rho_ou:>+11.4f} {rho_qnm:>+11.4f}  "
          f"{pred_ou:>+10.4f} {pred_qnm:>+10.4f}  {mean_dx:>7.4f}")
print()

print("─── FALSIFICATION CRITERIA ───────────────────────────────────")
print("  F4a: ΔlogL(OU)  < 0 with 20+ bins           → H0 (OU) falsified")
print("  F4b: ΔlogL(QNM) < ΔlogL(OU) with 20+ bins  → QNM adds no value")
print("  F6:  θ/ω_R inconsistent across DESI/Euclid  → QNM kernel incoherent")
print("  F7:  ω_R → 0 in Euclid fit                  → OU pure recovered, not QNM")
print()
print(SEP)

# ============================================================
# SECTION 10: FIGURES (4 panels)
# ============================================================

fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle(
    f'Stochastic Dark Energy: OU + QNM Kernel Test — {DATA_SOURCE}',
    fontsize=13, fontweight='bold'
)

# --- Panel 1: BAO residuals with model uncertainty bands ---
ax1 = axes[0, 0]
ax1.errorbar(z_eff, residuals, yerr=sigma, fmt='o', color='steelblue',
             capsize=4, zorder=5, label=DATA_SOURCE)
ax1.axhline(0, color='gray', linestyle='--', alpha=0.5, label='ΛCDM')
sig_ou  = np.array([abs(S_z[i]) * result_OU['sigma_X']  for i in range(len(z_eff))])
sig_qnm = np.array([abs(S_z[i]) * result_QNM['sigma_X'] for i in range(len(z_eff))])
ax1.fill_between(z_eff, -sig_ou, sig_ou, alpha=0.20, color='orange',
                 label=f"OU ±1σ (θ={result_OU['theta']:.2f})")
ax1.fill_between(z_eff, -sig_qnm, sig_qnm, alpha=0.15, color='crimson',
                 label=f"QNM ±1σ (θ={result_QNM['theta']:.2f}, "
                       f"ω_R={result_QNM['omega_R']:.2f})")
ax1.set_xlabel('z_eff'); ax1.set_ylabel('α − 1')
ax1.set_title('BAO Residuals + Model Uncertainty Bands')
ax1.legend(fontsize=8); ax1.grid(alpha=0.3)

# --- Panel 2: ω_R scan (ΔlogL vs omega_R) ---
ax2 = axes[0, 1]
ax2.plot(omega_R_grid, logL_scan - logL_LCDM, color='crimson', lw=2,
         label='QNM scan (θ=1.2, σ=calibrated)')
ax2.axhline(logL_OU_cal - logL_LCDM, color='orange', linestyle='--',
            label=f'OU calibrated (ΔlogL={logL_OU_cal-logL_LCDM:.2f})')
ax2.axvline(result_QNM['omega_R'], color='crimson', linestyle=':',
            alpha=0.8, label=f"MLE ω_R = {result_QNM['omega_R']:.3f}")
ax2.axhline(0, color='gray', alpha=0.4, label='ΛCDM')
ax2.set_xlabel('ω_R  (QNM oscillation frequency in e-folds)')
ax2.set_ylabel('ΔlogL vs ΛCDM')
ax2.set_title('ω_R Scan (calibrated θ, σ — diagnostic only)')
ax2.legend(fontsize=8); ax2.grid(alpha=0.3)
# Rayleigh limit annotation
from matplotlib.patches import FancyArrowPatch
omega_R_min_DESI = 2 * np.pi / 0.944
ax2.axvline(omega_R_min_DESI, color='navy', linestyle='-.', alpha=0.6,
            label=f'Rayleigh limit DESI\n(ω_R,min={omega_R_min_DESI:.1f})')
ax2.legend(fontsize=7)

# --- Panel 3: Linearity test Δη vs Δx ---
ax3 = axes[1, 0]
ax3.scatter(pairs_dx, pairs_deta / pairs_deta.max(), color='teal',
            alpha=0.7, s=60, zorder=5,
            label=f'DESI pairs (r={r_lin:.5f})')
p_fit = np.polyfit(pairs_dx, pairs_deta / pairs_deta.max(), 1)
x_line = np.linspace(0, pairs_dx.max() * 1.05, 100)
ax3.plot(x_line, np.polyval(p_fit, x_line), 'k--', alpha=0.5,
         label='Linear fit')
ax3.set_xlabel('Δx = |ln(1+z_i) − ln(1+z_j)|')
ax3.set_ylabel('Δη / max(Δη)  [normalized conformal time]')
ax3.set_title(f'Projection Validity: Δη vs Δx\n'
              f'(r={r_lin:.5f}; valid if r > 0.999)')
ax3.legend(fontsize=8); ax3.grid(alpha=0.3)

# --- Panel 4: Kernel shapes OU vs QNM ---
ax4 = axes[1, 1]
dx_plot   = np.linspace(0.0, 2.5, 400)
K_ou_plot = np.exp(-result_OU['theta']  * dx_plot)
K_qnm_plot = (np.exp(-result_QNM['theta'] * dx_plot)
              * np.cos(result_QNM['omega_R'] * dx_plot))
K_env     = np.exp(-result_QNM['theta'] * dx_plot)
ax4.plot(dx_plot, K_ou_plot,  color='orange', lw=2,
         label=f"OU: exp(−{result_OU['theta']:.2f}·Δx)")
ax4.plot(dx_plot, K_qnm_plot, color='crimson', lw=2,
         label=f"QNM: exp(−{result_QNM['theta']:.2f}·Δx)·cos({result_QNM['omega_R']:.2f}·Δx)")
ax4.plot(dx_plot,  K_env, 'k:', alpha=0.35, lw=1.5, label='QNM envelope')
ax4.plot(dx_plot, -K_env, 'k:', alpha=0.35, lw=1.5)
ax4.axhline(0, color='gray', alpha=0.4)
for x_b in x_bins:
    ax4.axvline(x_b - x_bins[0], color='steelblue', alpha=0.12, lw=0.8)
ax4.set_xlabel('Δx = |ln(1+z_i) − ln(1+z_j)|')
ax4.set_ylabel('C(Δx) / σ_X²')
ax4.set_title('Kernel Shapes: OU (H0) vs QNM (H1)\n'
              '(blue verticals: actual DESI bin separations)')
ax4.legend(fontsize=8); ax4.grid(alpha=0.3); ax4.set_ylim(-1.25, 1.25)

plt.tight_layout()
os.makedirs('plots', exist_ok=True)
plt.savefig('figures/test_desi_QNM.png', dpi=150, bbox_inches='tight')
print("Figure saved: figures/test_desi_QNM.png")