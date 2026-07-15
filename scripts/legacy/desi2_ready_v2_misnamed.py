"""
Test de Red Estocástica para DESI DR2
Fecha: Febrero 2026

INSTRUCCIONES PARA DESI DR2:
1. Descarga los alphas BAO de la Tabla principal de arXiv:2503.14738
2. Sustituye los arrays z_eff, alpha, sigma con los valores de DR2
3. Ejecuta: python desi2_ready_v2.py

NOVEDADES (extensión QNM):
- Kernel oscilatorio amortiguado: C(Δx) = σ² e^{-θΔx} cos(ω_R Δx)
- Motivación: modos cuasi-normales del horizonte de de Sitter
- OU puro es el límite omega_R → 0 (masa efectiva > 3H/2)
- H0 (OU puro) vs H1 (oscilatorio) comparados con AIC/BIC
- Test de linealidad Δη vs Δx para validar proyección temporal
"""

import numpy as np
from scipy.linalg import cholesky, solve_triangular
from scipy.optimize import minimize_scalar, minimize
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# DATOS - SUSTITUIR CON DESI DR2 CUANDO ESTÉN DISPONIBLES
# ============================================================
# Actualmente: DESI DR1 (arXiv:2404.03000, Tabla 3)
# Para DR2: reemplaza estos 3 arrays con valores de arXiv:2503.14738

z_eff  = np.array([0.295, 0.510, 0.706, 0.934, 1.321, 1.484, 2.330])
alpha  = np.array([0.9857, 0.9911, 0.9749, 0.9886, 0.9911, 1.0032, 0.9971])
sigma  = np.array([0.0093, 0.0077, 0.0067, 0.0046, 0.0071, 0.0153, 0.0082])

# ============================================================
# KERNEL DE SENSIBILIDAD BAO - S(z) = ∂ ln D_V / ∂Ω_Λ
# Calculado con cosmología base plana (Ω_m=0.315, H0=67.4)
# ⚠ Debe actualizarse con fiducial exacta de DESI (archivo YAML)
# ============================================================
S_z = np.array([-0.284, -0.462, -0.595, -0.719, -0.870, -0.917, -1.070])

# ============================================================
# PARÁMETROS MODELO OU BASE
# ============================================================
theta_base  = 1.2
sigma_OU    = 2.31e-2   # σ_X = desviación estándar de δΩ_Λ
sigma_X2    = sigma_OU**2

# ============================================================
# KERNELS DE CORRELACIÓN
# ============================================================

def kernel_OU(delta_x, theta, sigma_X2):
    """
    Kernel Ornstein-Uhlenbeck puro (H0).
    C(Δx) = σ_X² · exp(-θ · Δx)
    Límite: omega_R = 0 (masa efectiva > 3H/2 en de Sitter)
    """
    return sigma_X2 * np.exp(-theta * delta_x)

def kernel_QNM(delta_x, theta, sigma_X2, omega_R):
    """
    Kernel oscilatorio amortiguado (H1) - extensión QNM.
    C(Δx) = σ_X² · exp(-θ · Δx) · cos(omega_R · Δx)
    
    Motivación física:
    - Los modos cuasi-normales del horizonte de de Sitter tienen frecuencia
      compleja ω = ω_I + i·ω_R (Beyer 2011, López-Ortega 2012)
    - Para campo escalar con masa m_eff < 3H/2: ω_R real y != 0
      → correlaciones oscilatorias en el tiempo
    - Para m_eff >= 3H/2: ω_R → 0, recupera OU puro (H0)
    - theta ≡ ω_I (tasa de decaimiento), omega_R (frecuencia de oscilación)
    - La relación theta/omega_R = ω_I/ω_R depende de m_eff/H
      → restricción entre parámetros, no son independientes
    
    Proyección a BAO:
    C_α(zᵢ, zⱼ) = S(zᵢ) · S(zⱼ) · C(|ln(1+zᵢ) - ln(1+zⱼ)|)
    
    Validez: verificar linealidad Δη vs Δx en rango z < 3 (ver test abajo)
    """
    return sigma_X2 * np.exp(-theta * delta_x) * np.cos(omega_R * delta_x)

def build_cov_total(z_arr, sigma_arr, S_arr, theta, sigma_X2, omega_R=0.0):
    """
    Construye covarianza total: C_total = C_inst + C_kernel
    
    C_inst = diag(sigma_obs²)   [incertidumbre instrumental]
    C_kernel_ij = S(zᵢ)·S(zⱼ)·K(|Δx_ij|)   [señal estocástica]
    
    Si omega_R=0 → OU puro (H0)
    Si omega_R>0 → oscilatorio QNM (H1)
    """
    n = len(z_arr)
    x = np.log(1 + z_arr)  # variable logarítmica
    
    C_inst = np.diag(sigma_arr**2)
    C_signal = np.zeros((n, n))
    
    for i in range(n):
        for j in range(n):
            delta_x = abs(x[i] - x[j])
            if omega_R == 0.0:
                K = kernel_OU(delta_x, theta, sigma_X2)
            else:
                K = kernel_QNM(delta_x, theta, sigma_X2, omega_R)
            C_signal[i, j] = S_arr[i] * S_arr[j] * K
    
    return C_inst + C_signal

def log_likelihood(residuals, cov_matrix):
    """
    log-likelihood Gaussiano multivariado.
    logL = -0.5 * [r^T C^{-1} r + ln|C| + n·ln(2π)]
    """
    n = len(residuals)
    try:
        L = cholesky(cov_matrix, lower=True)
        y = solve_triangular(L, residuals, lower=True)
        log_det = 2 * np.sum(np.log(np.diag(L)))
        logL = -0.5 * (np.dot(y, y) + log_det + n * np.log(2 * np.pi))
        return logL
    except Exception:
        return -np.inf

# ============================================================
# TEST DE LINEALIDAD Δη vs Δx
# Valida que la proyección kernel QNM es válida en z < 3
# ============================================================

def test_linearity_deta_dx(z_arr, H0=67.4, Om=0.315):
    """
    Verifica la linealidad entre tiempo conforme η y x = ln(1+z).
    
    En de Sitter exacto: Δη ∝ e^{-x}  (no lineal)
    En el rango observacional z < 3 con Λ dominante tardío:
    la relación se aproxima a lineal → el kernel en Δη ≈ kernel en Δx
    
    Este test numérico justifica (o invalida) la proyección C(Δη) → C(Δx).
    """
    from scipy.integrate import quad
    
    c = 2.998e5  # km/s
    
    def H(z_):
        return H0 * np.sqrt(Om * (1+z_)**3 + (1 - Om))
    
    # tiempo conforme η(z) = ∫ dz / H(z) · c  [en unidades Mpc/km·s]
    eta = []
    x_arr = []
    for z in z_arr:
        eta_val, _ = quad(lambda zp: c / H(zp), 0, z)
        eta.append(eta_val)
        x_arr.append(np.log(1 + z))
    
    eta = np.array(eta)
    x_arr = np.array(x_arr)
    
    # Normalizar para comparar forma
    eta_norm = (eta - eta.min()) / (eta.max() - eta.min())
    x_norm   = (x_arr - x_arr.min()) / (x_arr.max() - x_arr.min())
    
    # Pearson r entre Δη y Δx para todos los pares
    pairs_deta = []
    pairs_dx   = []
    for i in range(len(z_arr)):
        for j in range(i+1, len(z_arr)):
            pairs_deta.append(abs(eta[i] - eta[j]))
            pairs_dx.append(abs(x_arr[i] - x_arr[j]))
    
    pairs_deta = np.array(pairs_deta)
    pairs_dx   = np.array(pairs_dx)
    r = np.corrcoef(pairs_deta, pairs_dx)[0, 1]
    
    return r, pairs_deta, pairs_dx, x_arr, eta

# ============================================================
# AJUSTE LIBRE DE PARÁMETROS
# ============================================================

def fit_model(residuals, z_arr, sigma_arr, S_arr, model='OU'):
    """
    Ajuste por MLE de los parámetros del modelo.
    H0 (OU):  parámetros libres = (theta, sigma_X2)     → k=2
    H1 (QNM): parámetros libres = (theta, sigma_X2, omega_R) → k=3
    """
    def neg_logL_OU(params):
        theta_, log_s2 = params
        if theta_ <= 0 or log_s2 < -20 or log_s2 > 0:
            return 1e10
        s2 = np.exp(log_s2)
        C = build_cov_total(z_arr, sigma_arr, S_arr, theta_, s2, omega_R=0.0)
        return -log_likelihood(residuals, C)
    
    def neg_logL_QNM(params):
        theta_, log_s2, omega_R_ = params
        if theta_ <= 0 or log_s2 < -20 or log_s2 > 0 or omega_R_ < 0:
            return 1e10
        s2 = np.exp(log_s2)
        C = build_cov_total(z_arr, sigma_arr, S_arr, theta_, s2, omega_R=omega_R_)
        return -log_likelihood(residuals, C)
    
    if model == 'OU':
        x0 = [1.2, np.log(sigma_X2)]
        res = minimize(neg_logL_OU, x0, method='Nelder-Mead',
                       options={'xatol': 1e-4, 'fatol': 1e-4, 'maxiter': 5000})
        theta_opt = res.x[0]
        s2_opt    = np.exp(res.x[1])
        omega_R_opt = 0.0
        logL_opt  = -res.fun
        k = 2
    else:  # QNM
        best_res = None
        best_val = np.inf
        # Grid search en omega_R para evitar mínimos locales
        for om0 in [0.5, 1.0, 2.0, 3.0, 5.0]:
            x0 = [1.2, np.log(sigma_X2), om0]
            res = minimize(neg_logL_QNM, x0, method='Nelder-Mead',
                           options={'xatol': 1e-4, 'fatol': 1e-4, 'maxiter': 5000})
            if res.fun < best_val:
                best_val = res.fun
                best_res = res
        theta_opt   = best_res.x[0]
        s2_opt      = np.exp(best_res.x[1])
        omega_R_opt = best_res.x[2]
        logL_opt    = -best_res.fun
        k = 3
    
    n    = len(residuals)
    AIC  = 2*k - 2*logL_opt
    BIC  = k*np.log(n) - 2*logL_opt
    
    return {
        'theta': theta_opt,
        'sigma_X': np.sqrt(s2_opt),
        'sigma_X2': s2_opt,
        'omega_R': omega_R_opt,
        'logL': logL_opt,
        'k': k,
        'AIC': AIC,
        'BIC': BIC,
    }

# ============================================================
# CÁLCULO PRINCIPAL
# ============================================================

residuals = alpha - 1.0

# --- ΛCDM baseline ---
C_LCDM = np.diag(sigma**2)
logL_LCDM = log_likelihood(residuals, C_LCDM)

# --- H0: OU puro con parámetros calibrados ---
C_OU_calib = build_cov_total(z_eff, sigma, S_z, theta_base, sigma_X2, omega_R=0.0)
logL_OU_calib = log_likelihood(residuals, C_OU_calib)

# --- H0: OU puro ajustado libre ---
fit_OU = fit_model(residuals, z_eff, sigma, S_z, model='OU')
C_OU_fit = build_cov_total(z_eff, sigma, S_z, fit_OU['theta'], fit_OU['sigma_X2'], omega_R=0.0)
logL_OU_fit = log_likelihood(residuals, C_OU_fit)

# --- H1: QNM oscilatorio ajustado libre ---
fit_QNM = fit_model(residuals, z_eff, sigma, S_z, model='QNM')
C_QNM_fit = build_cov_total(z_eff, sigma, S_z, fit_QNM['theta'], fit_QNM['sigma_X2'], fit_QNM['omega_R'])
logL_QNM_fit = log_likelihood(residuals, C_QNM_fit)

# --- Test de linealidad Δη vs Δx ---
r_linear, pairs_deta, pairs_dx, x_arr, eta_arr = test_linearity_deta_dx(z_eff)

# ============================================================
# BARRIDO DE omega_R (diagnóstico)
# ============================================================
omega_R_grid = np.linspace(0, 8, 80)
logL_scan = []
for om in omega_R_grid:
    C_tmp = build_cov_total(z_eff, sigma, S_z, theta_base, sigma_X2, omega_R=om)
    logL_scan.append(log_likelihood(residuals, C_tmp))
logL_scan = np.array(logL_scan)

# ============================================================
# CORRELACIONES DE LAG (covarianza completa)
# ============================================================
def lag_corr_full(residuals, cov_matrix, k):
    """Correlación de lag usando whitening con covarianza completa."""
    if k >= len(residuals):
        return np.nan
    try:
        L = cholesky(cov_matrix, lower=True)
        y = solve_triangular(L, residuals, lower=True)
        return np.corrcoef(y[:-k], y[k:])[0, 1]
    except Exception:
        return np.nan

# ============================================================
# PRINT RESULTADOS
# ============================================================
print("=" * 65)
print("TEST RED ESTOCÁSTICA - EXTENSIÓN QNM (v2)")
print("=" * 65)
print(f"N bins: {len(z_eff)}")
print()

print("─── COMPARACIÓN DE MODELOS ───────────────────────────────────")
print(f"{'Modelo':<30} {'logL':>8} {'ΔlogL':>8} {'k':>3} {'AIC':>8} {'BIC':>8}")
print("-" * 65)
print(f"{'ΛCDM (baseline)':<30} {logL_LCDM:>8.3f} {'—':>8} {'0':>3} {'—':>8} {'—':>8}")

dlogL_OU_c = logL_OU_calib - logL_LCDM
print(f"{'OU calibrado (θ=1.2)':<30} {logL_OU_calib:>8.3f} {dlogL_OU_c:>8.3f} {'2':>3} {'—':>8} {'—':>8}")

dlogL_OU_f = fit_OU['logL'] - logL_LCDM
print(f"{'OU ajustado libre (H0)':<30} {fit_OU['logL']:>8.3f} {dlogL_OU_f:>8.3f} {fit_OU['k']:>3} {fit_OU['AIC']:>8.3f} {fit_OU['BIC']:>8.3f}")

dlogL_QNM = fit_QNM['logL'] - logL_LCDM
print(f"{'QNM oscilatorio (H1)':<30} {fit_QNM['logL']:>8.3f} {dlogL_QNM:>8.3f} {fit_QNM['k']:>3} {fit_QNM['AIC']:>8.3f} {fit_QNM['BIC']:>8.3f}")

print()
print("─── PARÁMETROS AJUSTADOS ─────────────────────────────────────")
print(f"  H0 (OU puro):  θ = {fit_OU['theta']:.3f},  σ_X = {fit_OU['sigma_X']:.4f},  ω_R = 0 (fijo)")
print(f"  H1 (QNM):      θ = {fit_QNM['theta']:.3f},  σ_X = {fit_QNM['sigma_X']:.4f},  ω_R = {fit_QNM['omega_R']:.3f}")
print()

print("─── SELECCIÓN DE MODELO (AIC/BIC) ────────────────────────────")
dAIC = fit_OU['AIC'] - fit_QNM['AIC']
dBIC = fit_OU['BIC'] - fit_QNM['BIC']
print(f"  ΔAIC (H0 - H1) = {dAIC:+.3f}  {'→ H1 preferido' if dAIC > 0 else '→ H0 preferido'} (>2 = notable, >6 = fuerte)")
print(f"  ΔBIC (H0 - H1) = {dBIC:+.3f}  {'→ H1 preferido' if dBIC > 0 else '→ H0 preferido'} (>2 = positivo, >6 = fuerte)")
print()
print("  ⚠ Con 7 bins, selección de modelo es indicativa.")
print("  ⚠ Test decisivo requiere 20+ bins (Euclid DR1, oct 2026).")
print()

print("─── RELACIÓN θ/ω_R (test de consistencia QNM) ───────────────")
if fit_QNM['omega_R'] > 0.01:
    ratio = fit_QNM['theta'] / fit_QNM['omega_R']
    print(f"  θ/ω_R = {fit_QNM['theta']:.3f} / {fit_QNM['omega_R']:.3f} = {ratio:.3f}")
    print(f"  Este ratio = ω_I/ω_R fija m_eff/H via espectro QNM de de Sitter.")
    print(f"  Si QNM es el mecanismo real, este ratio debe ser consistente")
    print(f"  entre datasets independientes (DESI DR2, Euclid) → falsable.")
else:
    print(f"  ω_R ≈ 0: ajuste prefiere OU puro (masa efectiva >= 3H/2 en de Sitter)")
print()

print("─── TEST LINEALIDAD Δη vs Δx ─────────────────────────────────")
print(f"  Correlación Pearson(Δη, Δx) para todos los pares = {r_linear:.6f}")
if r_linear > 0.999:
    print(f"  ✓ Linealidad excelente en z ∈ [{z_eff.min():.3f}, {z_eff.max():.3f}]")
    print(f"    → Proyección kernel QNM: C(Δη) ≈ C(Δx) VÁLIDA en este rango")
elif r_linear > 0.99:
    print(f"  ~ Linealidad buena pero no perfecta. Considerar corrección de orden 2.")
else:
    print(f"  ✗ No lineal. La proyección C(Δη) → C(Δx) requiere jacobiano explícito.")
print()

print("─── CORRELACIONES DE LAG ─────────────────────────────────────")
print(f"{'Lag':<6} {'Obs(ΛCDM)':>10} {'Obs(OU)':>10} {'Obs(QNM)':>10} {'Pred OU':>10} {'Pred QNM':>10}")
print("-" * 60)
x_bins = np.log(1 + z_eff)
lags_lnz = [abs(x_bins[i+1] - x_bins[i]) for i in range(len(z_eff)-1)]
for k in [1, 2, 3]:
    rho_LCDM = lag_corr_full(residuals, C_LCDM, k)
    rho_OU_c = lag_corr_full(residuals, C_OU_calib, k)
    rho_QNM  = lag_corr_full(residuals, C_QNM_fit, k)
    # predicciones teóricas promediando sobre pares
    mean_dx = np.mean([abs(x_bins[i+k] - x_bins[i]) for i in range(len(z_eff)-k)])
    pred_OU  = np.exp(-fit_OU['theta'] * mean_dx)
    pred_QNM = np.exp(-fit_QNM['theta'] * mean_dx) * np.cos(fit_QNM['omega_R'] * mean_dx)
    print(f"  {k:<4} {rho_LCDM:>10.4f} {rho_OU_c:>10.4f} {rho_QNM:>10.4f} {pred_OU:>10.4f} {pred_QNM:>10.4f}")
print("  IC95% ≈ ±0.8 con 7 bins. Resultados indicativos, no significativos.")
print()

print("─── CRITERIOS DE FALSACIÓN (actualizados) ────────────────────")
print("  F4a: ΔlogL(OU) < 0 con 20+ bins → OU muere")
print("  F4b: ΔlogL(QNM) < ΔlogL(OU) con 20+ bins → extensión QNM no aporta")
print("  F6 (nuevo): θ/ω_R inconsistente entre DESI y Euclid → QNM incoherente")
print("  F7 (nuevo): ω_R → 0 en ajuste Euclid → OU puro recuperado (no QNM)")
print()
print("=" * 65)

# ============================================================
# FIGURA (4 paneles)
# ============================================================
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle('Red Estocástica + Extensión QNM — DESI BAO', fontsize=13, fontweight='bold')

# Panel 1: Residuos con bandas de ambos modelos
ax1 = axes[0, 0]
z_fine = np.linspace(z_eff.min(), z_eff.max(), 200)
ax1.errorbar(z_eff, residuals, yerr=sigma, fmt='o', color='steelblue',
             capsize=4, label='DESI DR1', zorder=5)
ax1.axhline(0, color='gray', linestyle='--', alpha=0.5, label='ΛCDM')
# Banda OU
sigma_floor_OU = np.sqrt([S_z[i]**2 * fit_OU['sigma_X2'] for i in range(len(z_eff))])
ax1.fill_between(z_eff, -sigma_floor_OU, sigma_floor_OU,
                 alpha=0.2, color='orange', label=f'OU ±1σ (θ={fit_OU["theta"]:.2f})')
# Banda QNM (σ efectiva = σ_X · |S(z)| · |cos oscila|^{1/2} aproximado)
sigma_floor_QNM = np.sqrt([S_z[i]**2 * fit_QNM['sigma_X2'] for i in range(len(z_eff))])
ax1.fill_between(z_eff, -sigma_floor_QNM, sigma_floor_QNM,
                 alpha=0.15, color='crimson', label=f'QNM ±1σ (θ={fit_QNM["theta"]:.2f}, ω_R={fit_QNM["omega_R"]:.2f})')
ax1.set_xlabel('z_eff')
ax1.set_ylabel('α - 1')
ax1.set_title('Residuos BAO')
ax1.legend(fontsize=8)
ax1.grid(alpha=0.3)

# Panel 2: Barrido de omega_R
ax2 = axes[0, 1]
ax2.plot(omega_R_grid, logL_scan - logL_LCDM, color='crimson', lw=2)
ax2.axhline(logL_OU_calib - logL_LCDM, color='orange', linestyle='--', label=f'OU calibrado (ΔlogL={logL_OU_calib-logL_LCDM:.2f})')
ax2.axvline(fit_QNM['omega_R'], color='crimson', linestyle=':', alpha=0.7,
            label=f'ω_R óptimo = {fit_QNM["omega_R"]:.2f}')
ax2.axhline(0, color='gray', linestyle='-', alpha=0.3, label='ΛCDM')
ax2.set_xlabel('ω_R (frecuencia oscilación QNM)')
ax2.set_ylabel('ΔlogL vs ΛCDM')
ax2.set_title('Barrido de ω_R (θ=1.2, σ calibrado)')
ax2.legend(fontsize=8)
ax2.grid(alpha=0.3)

# Panel 3: Test linealidad Δη vs Δx
ax3 = axes[1, 0]
ax3.scatter(pairs_dx, pairs_deta / pairs_deta.max(), color='teal', alpha=0.7,
            s=60, label=f'Pares DESI (r={r_linear:.5f})')
x_line = np.linspace(0, max(pairs_dx)*1.05, 100)
# Fit lineal
p = np.polyfit(pairs_dx, pairs_deta / pairs_deta.max(), 1)
ax3.plot(x_line, np.polyval(p, x_line), 'k--', alpha=0.5, label='Ajuste lineal')
ax3.set_xlabel('Δx = |ln(1+zᵢ) - ln(1+zⱼ)| (variable BAO)')
ax3.set_ylabel('Δη / max(Δη)  (tiempo conforme, normalizado)')
ax3.set_title(f'Validación proyección: Δη vs Δx\n(r={r_linear:.5f} → proyección válida si r>0.999)')
ax3.legend(fontsize=8)
ax3.grid(alpha=0.3)

# Panel 4: Kernels OU vs QNM (visualización teórica)
ax4 = axes[1, 1]
dx_plot = np.linspace(0, 2.5, 300)
K_OU  = np.exp(-fit_OU['theta']  * dx_plot)
K_QNM = np.exp(-fit_QNM['theta'] * dx_plot) * np.cos(fit_QNM['omega_R'] * dx_plot)
K_envelope = np.exp(-fit_QNM['theta'] * dx_plot)
ax4.plot(dx_plot, K_OU,  color='orange', lw=2, label=f'OU puro: e^(-{fit_OU["theta"]:.2f}Δx)')
ax4.plot(dx_plot, K_QNM, color='crimson', lw=2,
         label=f'QNM: e^(-{fit_QNM["theta"]:.2f}Δx)·cos({fit_QNM["omega_R"]:.2f}Δx)')
ax4.plot(dx_plot,  K_envelope, 'k:', alpha=0.4, lw=1.5, label='Envolvente QNM')
ax4.plot(dx_plot, -K_envelope, 'k:', alpha=0.4, lw=1.5)
ax4.axhline(0, color='gray', alpha=0.4)
# Marcar posiciones de los bins
for x_b in np.log(1 + z_eff):
    ax4.axvline(x_b - np.log(1 + z_eff[0]), color='steelblue', alpha=0.15, lw=0.8)
ax4.set_xlabel('Δx = |ln(1+zᵢ) - ln(1+zⱼ)|')
ax4.set_ylabel('Kernel de correlación C(Δx) / σ_X²')
ax4.set_title('Kernels: OU puro (H0) vs QNM oscilatorio (H1)\n(líneas verticales: separaciones reales entre bins)')
ax4.legend(fontsize=8)
ax4.grid(alpha=0.3)
ax4.set_ylim(-1.2, 1.2)

plt.tight_layout()
import os
os.makedirs('plots', exist_ok=True)
plt.savefig('plots/test_desi_QNM.png', dpi=150, bbox_inches='tight')
print("\nFigura guardada: plots/test_desi_QNM.png")