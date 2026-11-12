"""Paso 2: matrices B y A en p0, comprobación de A por diferencias finitas, ln Z(sigma_X) y cota 1.92."""
import json, time, math
import numpy as np
from scipy.linalg import cholesky, solve_triangular, cho_solve
from modelo import Modelo

t0 = time.time()
m = Modelo()
p0 = np.load("salidas/p0.npy")
out = {"p0": p0.tolist()}
N = m.N
y0, partes = m.vector_modelo(p0, partes=True)
r = m.d - y0
E, DM, DH = partes["E"], partes["DM"], partes["DH"]

# ------------------------------------------------------------------ B (1595 x 3)
pasos = m.cte["finite_difference_steps_B"]
B = np.zeros((m.n_dat, 3))
for j, (nombre, eps) in enumerate([("h", pasos["h"]), ("omega_m", pasos["omega_m"]), ("M", pasos["M"])]):
    pp, pm = p0.copy(), p0.copy()
    pp[j] += eps; pm[j] -= eps
    B[:, j] = (m.vector_modelo(pp) - m.vector_modelo(pm)) / (2 * eps)
print("B columna M: SN todos 1 ->", np.allclose(B[15:, 2], 1.0, atol=1e-9), " resto 0 ->", np.all(B[:15, 2] == 0))
print("B fila omega_m:", B[14].tolist())
out["B_col_M_sn_max_dev_de_1"] = float(np.max(np.abs(B[15:, 2] - 1.0)))
out["B_col_M_no_sn_max_abs"] = float(np.max(np.abs(B[:15, 2])))
out["B_fila_omega_m"] = B[14].tolist()

# ------------------------------------------------------------------ A analítica (1595 x N), convenio B
du = m.du
g = 0.5 * (1.0 + m.z) * E ** (-3.0)            # 1/2 (1+z_m) E_m^-3
# pesos trapecio acumulado: W_j = 1/2 (du_{j-1} + du_j) para 0<j<k ; w^{(k)}_k = 1/2 du_{k-1} ; w^{(k)}_0 = 1/2 du_0
W_int = np.zeros(N)
W_int[1:-1] = 0.5 * (du[:-1] + du[1:])
w_last = np.zeros(N); w_last[1:] = 0.5 * du          # peso del nodo k en D_M(u_k)
# dDM[k, j] = -g_j * w^{(k)}_j  (j>=1, j<=k) ; dDM[k, 0] = sum_{1<=m<=k} g_m w^{(k)}_m
dDM = np.zeros((N, N))
kk, jj = np.tril_indices(N, k=-1)                    # j < k
dDM[kk, jj] = -g[jj] * W_int[jj]
idx = np.arange(N)
dDM[idx, idx] = -g * w_last                          # j == k
# columna 0: w^{(k)}_0 g_0 - sum_{m<=k} w_m g_m  -> con signo: dDM[k,0] = -(w0 g0 - sum) = sum_{1<=m<=k} ...
dDM[:, 0] = -dDM[:, 1:].sum(axis=1)
dDH = np.zeros((N, N))
dDH[idx, idx] = -0.5 * E ** (-3.0)
dDH[:, 0] -= dDH[idx, idx]                          # - delta_0j
# dlnD
dlnDM = dDM / np.where(DM[:, None] > 0, DM[:, None], 1.0)
dlnDM[0, :] = 0.0
dlnDH = dDH / DH[:, None]
A = np.zeros((m.n_dat, N))
for i, (k, tipo) in enumerate(zip(m.k_bao, m.tipo_bao)):
    if tipo == "DM":
        dl = dlnDM[k]
    elif tipo == "DH":
        dl = dlnDH[k]
    else:
        dl = (2.0 / 3.0) * dlnDM[k] + (1.0 / 3.0) * dlnDH[k]
    A[i] = y0[i] * dl
kf = m.k_final
A[13] = -y0[13] * dlnDM[kf]                 # lectura literal: ∂ ln D~_M(u_final)
A[14] = 0.0
k, t = m.k_sn, m.t_sn
dDM_sn = (1.0 - t)[:, None] * dDM[k] + t[:, None] * dDM[k + 1]
A[15:] = (5.0 / math.log(10.0)) * dDM_sn / partes["DM_sn"][:, None]
sumas = A.sum(axis=1)
print("máx |suma de fila de A| =", np.max(np.abs(sumas)))
out["A_max_abs_suma_filas"] = float(np.max(np.abs(sumas)))
out["A_suma_filas_ok_1e-12"] = bool(np.max(np.abs(sumas)) < 1e-12)
del dDM, dDH, dlnDM, dlnDH

# ------------------------------------------------------------------ comprobación de A por diferencias finitas
eps = 1e-4
rng = np.random.default_rng(20261112)
nodos_todos = np.sort(rng.choice(N, size=20, replace=False))
nodos_bajo = np.sort(rng.choice(np.where(m.z <= 2.33)[0], size=20, replace=False))
# variante theta* "completa" (incluyendo el término de corrección en z*) para diagnóstico
def columna_fd(j):
    X = np.zeros(N)
    if j == 0:
        X[1:] = -eps
    else:
        X[j] = eps
    yp = m.vector_modelo(p0, X)
    ym = m.vector_modelo(p0, -X)
    return (yp - ym) / (2 * eps)

def comparar(nodos, umbral_rel):
    res = []
    peor = 0.0
    for j in nodos:
        col_fd = columna_fd(j)
        col_an = A[:, j]
        esc = np.max(np.abs(col_an))
        mask = np.abs(col_an) > umbral_rel * esc if esc > 0 else np.zeros(m.n_dat, bool)
        if mask.sum() == 0:
            res.append(dict(nodo=int(j), z=float(m.z[j]), n_filas=0)); continue
        rel = np.abs(col_fd[mask] - col_an[mask]) / np.abs(col_an[mask])
        abs_err = np.max(np.abs(col_fd - col_an))
        i_peor = np.where(mask)[0][np.argmax(rel)]
        res.append(dict(nodo=int(j), z=float(m.z[j]), n_filas=int(mask.sum()), max_rel=float(rel.max()),
                        fila_peor=int(i_peor), A_an_peor=float(col_an[i_peor]), A_fd_peor=float(col_fd[i_peor]),
                        max_abs=float(abs_err), max_abs_col=float(esc), ok_1pct=bool(rel.max() < 0.01)))
        peor = max(peor, rel.max())
    return res, peor

res_a, peor_a = comparar(nodos_todos, 1e-6)
res_b, peor_b = comparar(nodos_bajo, 1e-6)
print("FD (20 nodos de toda la malla): peor rel =", peor_a)
for d in res_a: print("  ", d)
print("FD (20 nodos con z<=2.33): peor rel =", peor_b)
for d in res_b: print("  ", d)
out["fd_nodos_toda_malla"] = res_a; out["fd_peor_toda_malla"] = float(peor_a)
out["fd_nodos_z_le_2.33"] = res_b; out["fd_peor_z_le_2.33"] = float(peor_b)
out["fd_umbral_relativo_fila_no_nula"] = 1e-6
out["fd_eps"] = eps
# diagnóstico extra: fila theta* con fórmula completa vs literal en nodos altos
j_alto = [j for j in nodos_todos if m.z[j] > 2.33]
diag = []
for j in j_alto[:5]:
    col_fd = columna_fd(j)
    diag.append(dict(nodo=int(j), z=float(m.z[j]), A_theta_an=float(A[13, j]), A_theta_fd=float(col_fd[13])))
out["fd_diag_theta_nodos_altos"] = diag
print("diag theta* nodos altos:", diag)

# ------------------------------------------------------------------ prior OU y lnZ
K = np.exp(-np.abs(m.u[:, None] - m.u[None, :]))   # theta = 1
AKAt = A @ K @ A.T
AKAt = 0.5 * (AKAt + AKAt.T)
np.save("salidas/A.npy", A); np.save("salidas/B.npy", B); np.save("salidas/AKAt.npy", AKAt)

def lnZ(sig):
    S = m.C + sig ** 2 * AKAt
    L = cholesky(S, lower=True)
    Sr = cho_solve((L, True), r)
    SB = cho_solve((L, True), B)
    F = B.T @ SB
    b = B.T @ Sr
    LF = cholesky(F, lower=True)
    bFb = float(b @ cho_solve((LF, True), b))
    logdetS = 2.0 * np.sum(np.log(np.diag(L)))
    logdetF = 2.0 * np.sum(np.log(np.diag(LF)))
    rSr = float(r @ Sr)
    val = -0.5 * (rSr - bFb + logdetS + logdetF)
    return val, dict(rSr=rSr, bFb=bFb, logdetS=float(logdetS), logdetF=float(logdetF),
                     chi2_min_lineal=rSr - bFb, dp=cho_solve((LF, True), b).tolist())

malla = list(m.cte["sigma_X_grid"])
ext = list(m.cte["sigma_X_grid_extension_if_max_at_last_node"])
tabla = {}
def evaluar(malla):
    for s in malla:
        if s not in tabla:
            v, det = lnZ(s); tabla[s] = (v, det)
            print(f"sigma_X={s:<5}  lnZ={v:.6f}  lnZ-lnZ0={v - tabla[0.0][0]:.6f}  chi2min_lin={det['chi2_min_lineal']:.4f} logdetS={det['logdetS']:.4f} logdetF={det['logdetF']:.6f}")
evaluar(malla)
extendida = False
sigmas = sorted(tabla)
vals = [tabla[s][0] - tabla[0.0][0] for s in sigmas]
if int(np.argmax(vals)) == len(sigmas) - 1:
    extendida = True
    print("máximo en el último nodo: extendiendo la malla")
    evaluar(ext)
    sigmas = sorted(tabla); vals = [tabla[s][0] - tabla[0.0][0] for s in sigmas]
out["malla_extendida"] = extendida
out["lnZ_rel"] = [dict(sigma_X=s, lnZ_menos_lnZ0=float(v), lnZ_abs=float(tabla[s][0]), **tabla[s][1]) for s, v in zip(sigmas, vals)]
out["lnZ_todos_finitos"] = bool(np.all(np.isfinite(vals)))
out["lnZ0_menos_lnZ0"] = float(tabla[0.0][0] - tabla[0.0][0])

# ------------------------------------------------------------------ regla 1.92
k0 = int(np.argmax(vals))
t = vals[k0] - 1.92
cota = None; k_cruce = None
for kk in range(k0, len(sigmas) - 1):
    if vals[kk] >= t > vals[kk + 1]:
        k_cruce = kk; break
if k_cruce is None:
    cota_txt = "sin cota en la malla"
elif k_cruce == 0:
    cota_txt = "cota < 0.01"
else:
    sk, sk1 = sigmas[k_cruce], sigmas[k_cruce + 1]
    cota = math.exp(math.log(sk) + (vals[k_cruce] - t) / (vals[k_cruce] - vals[k_cruce + 1]) * (math.log(sk1) - math.log(sk)))
    cota_txt = f"{cota:.4g}"
print(f"máximo global: k0={k0} sigma_X={sigmas[k0]} v={vals[k0]:.6f}; umbral t={t:.6f}; k cruce={k_cruce}; cota = {cota_txt}")
out["maximo"] = dict(k0=k0, sigma_X=sigmas[k0], v=float(vals[k0]))
out["umbral_t"] = float(t); out["k_cruce"] = k_cruce; out["cota"] = cota; out["cota_txt"] = cota_txt
out["tiempo_s"] = time.time() - t0
with open("salidas/paso2.json", "w") as f:
    json.dump(out, f, indent=1)
print("tiempo", out["tiempo_s"])
