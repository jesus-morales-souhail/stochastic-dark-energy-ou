import numpy as np
from scipy.optimize import minimize
from scipy.integrate import quad
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# DATOS REALES DE DESI DR2
# ============================================================
z_eff = np.array([0.295, 0.510, 0.706, 0.934, 1.321, 1.484, 2.330])
alpha_obs = np.array([1.0030, 0.9947, 1.0016, 0.9960, 1.0020, 0.9963, 1.0008])
sigma_obs = np.array([0.0097, 0.0072, 0.0057, 0.0049, 0.0063, 0.0088, 0.0120])

# Parámetros fiduciales
H0 = 67.4
Om = 0.315
Ol = 1.0 - Om
c_kms = 299792.458

# ============================================================
# FUNCIONES COSMOLÓGICAS ESTÁNDAR
# ============================================================
def E2_z(z, Om=0.315, Ol=0.685):
    """E(z)^2 = Ωm(1+z)^3 + ΩΛ(z)"""
    return Om * (1+z)**3 + Ol

def H_z(z, H0=67.4, Om=0.315, Ol=0.685):
    return H0 * np.sqrt(E2_z(z, Om, Ol))

def DM_standard(z, H0=67.4, Om=0.315, Ol=0.685):
    """Distancia comóvil transversal estándar (Mpc)."""
    integrand = lambda zp: 1.0 / H_z(zp, H0, Om, Ol)
    return c_kms * quad(integrand, 0, z)[0]

def DH_standard(z, H0=67.4, Om=0.315, Ol=0.685):
    return c_kms / H_z(z, H0, Om, Ol)

def DV_standard(z, H0=67.4, Om=0.315, Ol=0.685):
    dm = DM_standard(z, H0, Om, Ol)
    dh = DH_standard(z, H0, Om, Ol)
    return (dm**2 * dh * z) ** (1.0/3.0)

# ============================================================
# EoS EFECTIVA Y DENSIDAD DE ENERGÍA OSCURA MODIFICADA
# ============================================================
def t_z(z, H0=67.4, Om=0.315):
    """Tiempo conforme desde z=0 hasta z."""
    integrand = lambda zp: 1.0 / ((1+zp) * H_z(zp, H0, Om))
    return quad(integrand, 0, z)[0]

def Omega_DE_z(z, w0, wa, sigma, theta, Om=0.315):
    """
    Densidad de energía oscura relativa a la crítica.
    Ω_DE(z) = (1 - Ωm) * exp(3 * ∫_0^z [1 + w_eff(z')]/(1+z') dz')
    """
    def integrand(zp):
        w_CPL = w0 + wa * zp / (1 + zp)
        tz = t_z(zp)
        Omega_background = (1 - Om) * (1+zp)**3 / E2_z(zp, Om)
        if theta < 1e-6:
            correction = sigma**2 * tz / Omega_background
        else:
            correction = (sigma**2 / (2 * theta * Omega_background)) * (1 - np.exp(-2 * theta * tz))
        w_eff_val = w_CPL + correction
        # Restricción física: w_eff >= -1
        if w_eff_val < -1:
            w_eff_val = -1
        return (1 + w_eff_val) / (1 + zp)
    
    exponent = 3 * quad(integrand, 0, z, limit=100)[0]
    return (1 - Om) * np.exp(exponent)

def E2_eff(z, w0, wa, sigma, theta, Om=0.315):
    """E(z)^2 con Ω_DE(z) variable."""
    return Om * (1+z)**3 + Omega_DE_z(z, w0, wa, sigma, theta, Om)

def H_eff(z, w0, wa, sigma, theta, H0=67.4, Om=0.315):
    return H0 * np.sqrt(E2_eff(z, w0, wa, sigma, theta, Om))

def DM_eff(z, w0, wa, sigma, theta, H0=67.4, Om=0.315):
    integrand = lambda zp: 1.0 / H_eff(zp, w0, wa, sigma, theta, H0, Om)
    return c_kms * quad(integrand, 0, z, limit=100)[0]

def DH_eff(z, w0, wa, sigma, theta, H0=67.4, Om=0.315):
    return c_kms / H_eff(z, w0, wa, sigma, theta, H0, Om)

def DV_eff(z, w0, wa, sigma, theta, H0=67.4, Om=0.315):
    dm = DM_eff(z, w0, wa, sigma, theta, H0, Om)
    dh = DH_eff(z, w0, wa, sigma, theta, H0, Om)
    return (dm**2 * dh * z) ** (1.0/3.0)

# ============================================================
# VALORES FIDUCIALES PARA DESI DR2 (estándar)
# ============================================================
# Usamos las distancias estándar con Ωm=0.315, Ol=0.685
DV_fid = DV_standard(0.295)
DM_fid = np.array([DM_standard(z) for z in [0.510, 0.706, 0.934, 1.321, 1.484, 2.330]])
DH_fid = np.array([DH_standard(z) for z in [0.510, 0.706, 0.934, 1.321, 1.484, 2.330]])
z_fid_DM = np.array([0.510, 0.706, 0.934, 1.321, 1.484, 2.330])

def alpha_iso_pred(z_idx, w0, wa, sigma, theta):
    z = z_eff[z_idx]
    if z_idx == 0:
        DV_pred = DV_eff(z, w0, wa, sigma, theta)
        return DV_pred / DV_fid
    else:
        i = z_idx - 1
        DM_pred = DM_eff(z, w0, wa, sigma, theta)
        DH_pred = DH_eff(z, w0, wa, sigma, theta)
        return (DM_pred / DM_fid[i])**(2/3) * (DH_pred / DH_fid[i])**(1/3)

# ============================================================
# VEROSIMILITUD
# ============================================================
def log_likelihood(params):
    w0, wa, sigma, theta = params
    if sigma < 0 or theta < 0:
        return -1e10
    try:
        alpha_pred = np.array([alpha_iso_pred(i, w0, wa, sigma, theta) for i in range(7)])
        residuals = alpha_obs - alpha_pred
        chi2 = np.sum(residuals**2 / sigma_obs**2)
        return -0.5 * chi2
    except:
        return -1e10

def neg_log_likelihood(params):
    return -log_likelihood(params)

# ============================================================
# AJUSTE MLE
# ============================================================
print("="*60)
print("AJUSTE MLE CON EoS EFECTIVA (VERSIÓN CORREGIDA)")
print("="*60)

# Valores iniciales: CPL best-fit de DESI DR2
x0 = [-0.87, -0.41, 1e-4, 0.001]
bounds = [(-2, 0), (-2, 2), (0, 0.01), (1e-4, 10)]

result = minimize(neg_log_likelihood, x0, method='L-BFGS-B', bounds=bounds, options={'maxiter': 1000})
w0_opt, wa_opt, sigma_opt, theta_opt = result.x

print(f"w0   = {w0_opt:.4f}")
print(f"wa   = {wa_opt:.4f}")
print(f"σ    = {sigma_opt:.2e}")
print(f"θ    = {theta_opt:.4f}")
print(f"logL = {-result.fun:.3f}")

# ============================================================
# AJUSTE CPL PURO
# ============================================================
def log_likelihood_CPL(params):
    w0, wa = params
    try:
        alpha_pred = np.array([alpha_iso_pred(i, w0, wa, 0.0, 0.001) for i in range(7)])
        residuals = alpha_obs - alpha_pred
        chi2 = np.sum(residuals**2 / sigma_obs**2)
        return -0.5 * chi2
    except:
        return -1e10

result_CPL = minimize(lambda p: -log_likelihood_CPL(p), [-0.87, -0.41],
                      method='L-BFGS-B', bounds=[(-2, 0), (-2, 2)])
w0_CPL, wa_CPL = result_CPL.x
logL_CPL = -result_CPL.fun

print("\n" + "="*60)
print("AJUSTE CPL PURO")
print("="*60)
print(f"w0_CPL   = {w0_CPL:.4f}")
print(f"wa_CPL   = {wa_CPL:.4f}")
print(f"logL_CPL = {logL_CPL:.3f}")

# ============================================================
# COMPARACIÓN
# ============================================================
print("\n" + "="*60)
print("COMPARACIÓN DE MODELOS")
print("="*60)
delta_logL = logL_CPL - (-result.fun)
print(f"ΔlogL (CPL - EoS_eff) = {delta_logL:.3f}")
if delta_logL > 2:
    print("CPL es preferido sobre EoS_eff (evidencia fuerte).")
elif delta_logL > 0.5:
    print("CPL es ligeramente preferido sobre EoS_eff.")
else:
    print("No hay evidencia a favor de ninguno de los dos modelos.")

# ============================================================
# PREDICCIÓN PARA EUCLID
# ============================================================
def alpha_iso_pred_z(z, w0, wa, sigma, theta):
    # Interpolamos los valores fiduciales de DM/DH a z
    DM_fid_interp = np.interp(z, z_fid_DM, DM_fid)
    DH_fid_interp = np.interp(z, z_fid_DM, DH_fid)
    DM_pred = DM_eff(z, w0, wa, sigma, theta)
    DH_pred = DH_eff(z, w0, wa, sigma, theta)
    return (DM_pred / DM_fid_interp)**(2/3) * (DH_pred / DH_fid_interp)**(1/3)

z_euclid = np.linspace(0.9, 1.8, 20)
try:
    alpha_euclid = np.array([alpha_iso_pred_z(z, w0_opt, wa_opt, sigma_opt, theta_opt) for z in z_euclid])
    alpha_euclid_CPL = np.array([alpha_iso_pred_z(z, w0_CPL, wa_CPL, 0.0, 0.001) for z in z_euclid])
    
    print("\n" + "="*60)
    print("PREDICCIÓN PARA EUCLID DR1")
    print("="*60)
    print(f"z_min = {z_euclid[0]:.2f}, z_max = {z_euclid[-1]:.2f}, N_bins = {len(z_euclid)}")
    print(f"alpha_iso medio (EoS_eff) = {np.mean(alpha_euclid):.4f}")
    print(f"alpha_iso medio (CPL)     = {np.mean(alpha_euclid_CPL):.4f}")
    print(f"Diferencia media = {np.mean(alpha_euclid - alpha_euclid_CPL):.6f}")
except Exception as e:
    print(f"\n⚠️ Error en la predicción: {e}")

print("\n" + "="*60)
print("FIN DEL ANÁLISIS")
print("="*60)