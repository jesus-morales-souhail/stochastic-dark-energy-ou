"""Paso 1: comprobaciones de datos/malla y ajuste del punto fiducial p0."""
import json, sys, time, platform
import numpy as np, scipy
from scipy.optimize import least_squares
from modelo import Modelo

t0 = time.time()
m = Modelo()
out = {}
out["versiones"] = dict(python=sys.version, numpy=np.__version__, scipy=scipy.__version__,
                        plataforma=platform.platform())
print("python", sys.version.split()[0], "numpy", np.__version__, "scipy", scipy.__version__)

# --- comprobaciones de datos y malla
print("SN totales en .dat:", m.sn["n_total"], " seleccionadas:", m.n_sn)
print("u_max =", repr(m.u_max), " (constantes grid_u_max =", repr(m.cte["grid_u_max"]), ")")
print("n lineal =", m.n_lin, " N nodos tras añadir BAO =", m.N, " (esperado", m.cte["grid_n_nodes"], ")")
print("asimetría máx C_SN =", m.sn["cov_asim"])
print("pasos du: min", m.du.min(), "max", m.du.max())
print("índices nodos BAO:", m.k_bao.tolist())
print("longitud vector datos:", m.n_dat)
out["n_sn_total"] = int(m.sn["n_total"]); out["n_sn_sel"] = int(m.n_sn)
out["u_max"] = m.u_max; out["n_lineal"] = int(m.n_lin); out["N_nodos"] = int(m.N)
out["N_esperado"] = int(m.cte["grid_n_nodes"]); out["n_datos"] = int(m.n_dat)
out["asimetria_C_SN"] = m.sn["cov_asim"]
out["du_min"] = float(m.du.min()); out["du_max"] = float(m.du.max())
out["k_bao"] = m.k_bao.tolist()
out["zHD_sn_min"] = float(m.sn["zHD"].min()); out["zHD_sn_max"] = float(m.sn["zHD"].max())

# --- punto de arranque
ff = m.cte["fiducial_fit"]
h0, om0 = ff["start_h"], ff["start_omega_m"]
M0 = m.M_analitico(h0, om0)
chi2_ini = m.chi2((h0, om0, M0))
print(f"arranque: h={h0} om={om0} M0={M0:.9f}  chi2={chi2_ini:.6f}")
out["arranque"] = dict(h=h0, omega_m=om0, M=M0, chi2=chi2_ini)

# --- ajuste
lb = [ff["bounds_h"][0], ff["bounds_omega_m"][0], ff["bounds_M"][0]]
ub = [ff["bounds_h"][1], ff["bounds_omega_m"][1], ff["bounds_M"][1]]

def ajustar(x0, gtol=None, tag=""):
    kw = dict(method="trf", bounds=(lb, ub), x_scale=ff["x_scale"], ftol=1e-10, xtol=1e-10, max_nfev=200)
    if gtol is not None:
        kw["gtol"] = gtol
    res = least_squares(lambda p: m.residuo_blanqueado(p), x0, **kw)
    print(f"[{tag}] status={res.status} nfev={res.nfev} msg='{res.message}'")
    print(f"[{tag}] p = {res.x.tolist()}  chi2 = {2*res.cost:.6f}")
    return res

historial = []
res = ajustar([h0, om0, M0], tag="ajuste 1")
historial.append(dict(status=int(res.status), nfev=int(res.nfev), message=res.message,
                      p=res.x.tolist(), chi2=float(2 * res.cost)))
n_rep = 0
while res.status == 0 and n_rep < 10:   # parada por max_nfev: repetir desde el punto alcanzado
    n_rep += 1
    res = ajustar(res.x, tag=f"repetición {n_rep}")
    historial.append(dict(status=int(res.status), nfev=int(res.nfev), message=res.message,
                          p=res.x.tolist(), chi2=float(2 * res.cost)))
p0 = res.x.copy()
chi2_0 = m.chi2(p0)
print(f"p0 = h={p0[0]:.6f}  omega_m={p0[1]:.6f}  M={p0[2]:.6f}   chi2(p0)={chi2_0:.4f}")
out["ajuste_historial"] = historial
out["ajuste_repeticiones_por_max_nfev"] = n_rep
out["p0"] = dict(h=float(p0[0]), omega_m=float(p0[1]), M=float(p0[2]))
out["chi2_p0"] = chi2_0

# --- sensibilidad a gtol (no fijado en la especificación)
res_g = ajustar(p0, gtol=1e-14, tag="sens. gtol=1e-14 desde p0")
out["sensibilidad_gtol"] = dict(p=res_g.x.tolist(), chi2=float(2 * res_g.cost), status=int(res_g.status),
                                dif_p=(res_g.x - p0).tolist())
# --- sensibilidad al punto de arranque
res_b = ajustar([0.70, 0.15, M0], tag="sens. arranque (0.70,0.15)")
out["sensibilidad_arranque"] = dict(p=res_b.x.tolist(), chi2=float(2 * res_b.cost), dif_p=(res_b.x - p0).tolist())

# --- descomposición de chi2 en p0
y, partes = m.vector_modelo(p0, partes=True)
r = m.d - y
from scipy.linalg import solve_triangular
w13 = solve_triangular(np.linalg.cholesky(m.C_bao), r[:13], lower=True)
out["chi2_bloques"] = dict(bao=float(w13 @ w13), theta=float(r[13] ** 2 / m.C[13, 13]),
                           omega_m=float(r[14] ** 2 / m.C[14, 14]),
                           sn=float(chi2_0 - w13 @ w13 - r[13] ** 2 / m.C[13, 13] - r[14] ** 2 / m.C[14, 14]))
out["modelo_p0"] = dict(y_bao=partes["y_bao"].tolist(), y_theta=partes["y_th"], rd=partes["rd"], rs=partes["rs"],
                        z_star=partes["zs"], DM_final=float(partes["DM"][-1]), DM_star=float(partes["DM_star"]))
print("chi2 por bloques:", out["chi2_bloques"])
print("modelo en p0:", out["modelo_p0"])
out["tiempo_s"] = time.time() - t0
with open("salidas/paso1.json", "w") as f:
    json.dump(out, f, indent=1)
np.save("salidas/p0.npy", p0)
print("tiempo", out["tiempo_s"])
