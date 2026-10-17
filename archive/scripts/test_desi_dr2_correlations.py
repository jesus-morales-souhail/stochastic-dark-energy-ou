#!/usr/bin/env python3
# Lag-correlation test on DESI DR2 BAO residuals (OU kernel)

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import cholesky, solve_triangular
from scipy import stats
import sys
import os

class bcolors:
    HEADER = '\033[95m'; OKBLUE = '\033[94m'; OKGREEN = '\033[92m'
    WARNING = '\033[93m'; FAIL = '\033[91m'; ENDC = '\033[0m'; BOLD = '\033[1m'

print(f"{bcolors.HEADER}{'='*70}{bcolors.ENDC}")
print(f"{bcolors.BOLD}BAO residual lag correlations — DESI DR2{bcolors.ENDC}")
print(f"{bcolors.BOLD}   Ornstein–Uhlenbeck stochastic model{bcolors.ENDC}")
print(f"{bcolors.HEADER}{'='*70}{bcolors.ENDC}")

# ------------------------------------------------------------
# 1. CARGAR DATOS
# ------------------------------------------------------------
# Intentamos diferentes nombres de archivo posibles
archivos_data = [
    "desi_dr2_bao_alphas.txt",  # Formato: z, alpha, error
    "DESI_DR2_BAO.txt",         # Formato alternativo [citation:2]
    "../cobaya_packages_dr2/data/bao_data/desi_dr2_bao.npz"  # Si se instaló con Cobaya
]

data_cargados = False
for archivo in archivos_data:
    if os.path.exists(archivo):
        print(f"{bcolors.OKGREEN}Loading: {archivo}{bcolors.ENDC}")
        if archivo.endswith('.npz'):
            datos = np.load(archivo)
            z = datos['z_eff']
            alpha = datos['alpha']
            error = datos['alpha_err']
            # La covarianza completa puede estar en datos['cov']
            if 'cov' in datos:
                cov = datos['cov']
            else:
                cov = np.diag(error**2)
        else:
            # Asumimos formato texto: columnas z, alpha, error
            data = np.loadtxt(archivo)
            z = data[:, 0]
            alpha = data[:, 1]
            error = data[:, 2] if data.shape[1] > 2 else None
            cov = np.diag(error**2) if error is not None else None
        data_cargados = True
        break

if not data_cargados:
    print(f"{bcolors.FAIL}No DESI DR2 data files found.{bcolors.ENDC}")
    print(f"   Ejecuta primero: python obtener_datos_desi_dr2.py")
    sys.exit(1)

# ------------------------------------------------------------
# 2. CARGAR COVARIANZA (si no se cargó antes)
# ------------------------------------------------------------
if cov is None:
    archivos_cov = [
        "desi_dr2_bao_cov.txt",
        "DESI_DR2_BAO_Cov.txt",  # [citation:2]
        "desi_dr2_bao_covariance.txt"
    ]
    for archivo in archivos_cov:
        if os.path.exists(archivo):
            print(f"Loading covariance: {archivo}")
            cov = np.loadtxt(archivo)
            break
    
    if cov is None:
        print(f"{bcolors.WARNING}No covariance file; using diagonal errors.{bcolors.ENDC}")
        if error is not None:
            cov = np.diag(error**2)
        else:
            print(f"{bcolors.FAIL}No error column; aborting.{bcolors.ENDC}")
            sys.exit(1)

# ------------------------------------------------------------
# 3. RESIDUOS Y WHITENING
# ------------------------------------------------------------
# Para DESI, alpha ya está normalizado: alpha=1 para ΛCDM fiducial [citation:7]
residuals = alpha - 1.0

print(f"\n{bcolors.BOLD} Estadísticas básicas:{bcolors.ENDC}")
print(f"   Número de bins: {len(z)}")
print(f"   Redshifts: {z}")
print(f"   Alphas: {alpha}")
print(f"   Residuos (alpha-1): {np.round(residuals, 4)}")

# Whitening con la matriz de covarianza completa
try:
    L = cholesky(cov, lower=True)
    y = solve_triangular(L, residuals, lower=True)
    print(f"{bcolors.OKGREEN}Whitening with full covariance complete.{bcolors.ENDC}")
except np.linalg.LinAlgError:
    print(f"{bcolors.WARNING}Cholesky failed; adding regularisation...{bcolors.ENDC}")
    cov_reg = cov + np.eye(len(cov)) * 1e-6
    L = cholesky(cov_reg, lower=True)
    y = solve_triangular(L, residuals, lower=True)

# ------------------------------------------------------------
# 4. CORRELACIONES LAG
# ------------------------------------------------------------
def lag_correlation(y, k):
    """Calcula correlación de Pearson para lag k"""
    if k >= len(y):
        return np.nan
    a = y[:-k]
    b = y[k:]
    return np.corrcoef(a, b)[0, 1]

def lag_correlation_with_ci(y, k, alpha=0.05):
    """Calcula correlación e intervalo de confianza aproximado"""
    n = len(y) - k
    if n < 4:
        return np.nan, (np.nan, np.nan)
    
    rho = lag_correlation(y, k)
    # Transformada de Fisher para IC
    z_r = np.arctanh(np.clip(rho, -0.9999, 0.9999))
    se = 1.0 / np.sqrt(n - 3)
    z_crit = stats.norm.ppf(1 - alpha/2)
    
    ci_low = np.tanh(z_r - z_crit * se)
    ci_high = np.tanh(z_r + z_crit * se)
    
    return rho, (ci_low, ci_high)

print(f"\n{bcolors.BOLD}Results:{bcolors.ENDC}")
print("-" * 50)
print(f"{'Lag':<6} {'Correlación':<12} {'IC 95%':<30} {'n_efectivo':<10}")
print("-" * 50)

resultados = []
for k in [1, 2, 3, 4, 5]:  # DR2 tiene más bins, podemos mirar más lags
    rho, (ci_low, ci_high) = lag_correlation_with_ci(y, k)
    n_eff = len(y) - k
    if not np.isnan(rho):
        print(f"  {k:<5} {rho:<12.4f} [{ci_low:<6.4f}, {ci_high:<6.4f}]  {n_eff:<10}")
        resultados.append((k, rho, ci_low, ci_high))

print("-" * 50)

# ------------------------------------------------------------
# 5. COMPARACIÓN CON MODELO OU
# ------------------------------------------------------------
print(f"\n{bcolors.BOLD} PREDICCIÓN DEL MODELO OU (Red Estocástica):{bcolors.ENDC}")
print(f"   Lag 1: 0.78")
print(f"   Lag 2: 0.62")
print(f"   Lag 3: 0.49")

print(f"\n{bcolors.BOLD} INTERPRETACIÓN:{bcolors.ENDC}")
print(f"   • Si las correlaciones observadas están DENTRO del IC 95% de la predicción OU,")
print(f"     el modelo es consistente con los datos.")
print(f"   • Si están FUERA del IC 95% y son consistentemente positivas, el modelo gana puntos.")
print(f"   • Si son cercanas a 0 (ruido blanco), has puesto un límite superior a las fluctuaciones.")
print(f"   • Con DR2, los intervalos serán mucho más estrechos que con DR1.")

# ------------------------------------------------------------
# 6. VISUALIZACIÓN
# ------------------------------------------------------------
fig, axes = plt.subplots(2, 1, figsize=(10, 8))

# Gráfico 1: Residuos whitened
axes[0].plot(z, y, 'o-', color='blue', markersize=8, linewidth=2)
axes[0].axhline(0, color='gray', linestyle='--', alpha=0.7)
axes[0].set_xlabel('Redshift z', fontsize=12)
axes[0].set_ylabel('Residuos whitened', fontsize=12)
axes[0].set_title('DESI DR2: Residuos de BAO (whitened)', fontsize=14)
axes[0].grid(True, alpha=0.3)
axes[0].fill_between(z, -2, 2, color='gray', alpha=0.1, label='Zona de ruido blanco (2σ)')
axes[0].legend()

# Gráfico 2: Correlaciones
lags_plot = [r[0] for r in resultados]
corrs_plot = [r[1] for r in resultados]
err_low = [r[1] - r[2] for r in resultados]
err_high = [r[3] - r[1] for r in resultados]
err = [err_low, err_high]

axes[1].errorbar(lags_plot, corrs_plot, yerr=err, fmt='o', color='red', 
                 markersize=8, capsize=5, capthick=2, label='Observado')
axes[1].axhline(0, color='gray', linestyle='-', alpha=0.5)
axes[1].axhline(0.78, color='green', linestyle='--', alpha=0.7, label='OU Lag 1')
axes[1].axhline(0.62, color='orange', linestyle='--', alpha=0.7, label='OU Lag 2')
axes[1].axhline(0.49, color='purple', linestyle='--', alpha=0.7, label='OU Lag 3')
axes[1].set_xlabel('Lag (k)', fontsize=12)
axes[1].set_ylabel('Correlación', fontsize=12)
axes[1].set_title('Correlaciones lag en residuos BAO', fontsize=14)
axes[1].grid(True, alpha=0.3)
axes[1].legend()
axes[1].set_ylim(-1.1, 1.1)

plt.tight_layout()
import os
os.makedirs('plots', exist_ok=True)
plt.savefig('plots/test_desi_dr2_resultado.png', dpi=150)
print(f"\n{bcolors.OKGREEN}Figure saved: plots/test_desi_dr2_resultado.png{bcolors.ENDC}")

# ------------------------------------------------------------
# 7. Summary
# ------------------------------------------------------------
print(f"\n{bcolors.HEADER}{'='*70}{bcolors.ENDC}")
print(f"{bcolors.BOLD}Summary{bcolors.ENDC}")
print(f"{bcolors.HEADER}{'='*70}{bcolors.ENDC}")

if len(z) >= 12:
    print(f"With {len(z)} bins the lag test has limited statistical power.")
    
    # Evaluar si las correlaciones son consistentes con OU
    if resultados and resultados[0][1] > 0.5:
        print(f"   {bcolors.OKGREEN} Las correlaciones son positivas y altas (>0.5).{bcolors.ENDC}")
        print(f"   Esto APOYA la hipótesis de la Red Estocástica.")
    elif resultados and resultados[0][1] < 0.2 and abs(resultados[0][1]) < 0.2:
        print(f"   {bcolors.WARNING}Lag correlations are consistent with zero (<0.2).{bcolors.ENDC}")
        print(f"   Esto PONE UN LÍMITE SUPERIOR a las fluctuaciones estocásticas.")
        print(f"   Es un resultado publicable como 'constraint'.")
    else:
        print(f"   {bcolors.WARNING}Intermediate result; more data needed.{bcolors.ENDC}")
else:
    print(f"With {len(z)} bins, confidence intervals remain wide.")
    print(f"   Los resultados son preliminares hasta tener la muestra completa de DR2.")

print(f"\n{bcolors.BOLD} Para publicar:{bcolors.ENDC}")
print(f"   • Si el resultado apoya el modelo: 'Evidence for stochastic dark energy from DESI DR2'")
print(f"   • Si el resultado pone un límite: 'Constraints on stochastic dark energy from DESI DR2'")
print(f"\n{bcolors.HEADER}{'='*70}{bcolors.ENDC}")