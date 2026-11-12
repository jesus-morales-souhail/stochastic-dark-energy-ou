"""
Reproducción a ciegas: modelo exacto (convenio B), datos, malla y utilidades.
Escrito desde cero siguiendo únicamente ESPECIFICACION.md.
"""
import json
import math
import os
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.linalg import cholesky, solve_triangular

AQUI = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(AQUI)


# ---------------------------------------------------------------- constantes
def cargar_constantes():
    with open(os.path.join(RAIZ, "tablas", "constantes.json")) as f:
        return json.load(f)


# ---------------------------------------------------------------- tablas r_d, r_s, z*
class Tablas:
    """Splines cúbicos (not-a-knot, por defecto de scipy) en omega_m."""

    def __init__(self):
        t = np.loadtxt(os.path.join(RAIZ, "tablas", "rd_rs_zstar_vs_omega_m.csv"),
                       delimiter=",", comments="#")
        self.om = t[:, 0]
        self.rd = CubicSpline(t[:, 0], t[:, 1])      # Mpc
        self.rs = CubicSpline(t[:, 0], t[:, 2])      # Mpc
        self.zs = CubicSpline(t[:, 0], t[:, 3])
        self.n_nodos = t.shape[0]


# ---------------------------------------------------------------- malla en u
def construir_malla(cte):
    z_ref = cte["z_star_reference"]
    u_max = math.log(1.0 + z_ref)
    du = cte["grid_step_du"]
    n = math.ceil(u_max / du) + 1
    u = np.linspace(0.0, u_max, n)
    z_bao = sorted({row["z"] for row in cte["bao_rows"]})
    u = np.concatenate([u, np.log1p(np.array(z_bao))])
    u = np.unique(u)          # ordena y elimina duplicados exactos
    return u, u_max, n, z_bao


# ---------------------------------------------------------------- datos
def cargar_bao(cte):
    d = os.path.join(RAIZ, "datos", "desi_dr2")
    filas = []
    with open(os.path.join(d, "desi_gaussian_bao_ALL_GCcomb_mean.txt")) as f:
        for linea in f:
            if linea.startswith("#") or not linea.strip():
                continue
            z, val, nombre = linea.split()
            filas.append((float(z), float(val), nombre))
    cov = np.loadtxt(os.path.join(d, "desi_gaussian_bao_ALL_GCcomb_cov.txt"))
    cov = 0.5 * (cov + cov.T)
    z = np.array([r[0] for r in filas])
    val = np.array([r[1] for r in filas])
    tipo = [r[2].split("_")[0]  for r in filas]   # DV / DM / DH
    # comprobación de coherencia con constantes.json
    for i, row in enumerate(cte["bao_rows"]):
        assert abs(row["z"] - z[i]) < 1e-9, (i, row, z[i])
        assert row["kind"] == tipo[i], (i, row, tipo[i])
    return z, val, tipo, cov


def cargar_pantheon():
    d = os.path.join(RAIZ, "datos", "pantheon_plus")
    with open(os.path.join(d, "Pantheon+SH0ES.dat")) as f:
        cab = f.readline().split()
        filas = [l.split() for l in f if l.strip()]
    ic = {c: i for i, c in enumerate(cab)}
    n_tot = len(filas)
    zHD = np.array([float(r[ic["zHD"]]) for r in filas])
    zHEL = np.array([float(r[ic["zHEL"]]) for r in filas])
    mb = np.array([float(r[ic["m_b_corr"]]) for r in filas])
    iscal = np.array([int(float(r[ic["IS_CALIBRATOR"]])) for r in filas])
    sel = np.where((zHD > 0.01) & (iscal == 0))[0]
    with open(os.path.join(d, "Pantheon+SH0ES_STAT+SYS.cov")) as f:
        ncov = int(f.readline())
        cov = np.fromfile(f, sep="\n", dtype=float)
    assert ncov == n_tot, (ncov, n_tot)
    assert cov.size == ncov * ncov, cov.size
    cov = cov.reshape(ncov, ncov)
    cov = cov[np.ix_(sel, sel)]
    return dict(n_total=n_tot, sel=sel, zHD=zHD[sel], zHEL=zHEL[sel],
                mb=mb[sel], cov=cov, cov_asim=float(np.max(np.abs(cov - cov.T))))


# ---------------------------------------------------------------- modelo exacto
class Modelo:
    def __init__(self):
        self.cte = cargar_constantes()
        c = self.cte
        self.c_km_s = c["c_km_s"]
        self.Or_h2 = c["Omega_r_h2"]
        self.z_ref = c["z_star_reference"]
        self.CAL = c["theta_calibration_CAL"]
        self.tab = Tablas()
        self.u, self.u_max, self.n_lin, self.z_bao = construir_malla(c)
        self.N = self.u.size
        self.z = np.expm1(self.u)
        self.du = np.diff(self.u)
        # BAO
        self.z_bao_f, self.d_bao, self.tipo_bao, self.C_bao = cargar_bao(c)
        self.k_bao = np.array([self._indice_nodo(math.log1p(z)) for z in self.z_bao_f])
        # SN
        self.sn = cargar_pantheon()
        self.n_sn = self.sn["sel"].size
        u_hd = np.log1p(self.sn["zHD"])
        k = np.searchsorted(self.u, u_hd, side="right") - 1
        k = np.clip(k, 0, self.N - 2)
        self.k_sn = k
        self.t_sn = (u_hd - self.u[k]) / (self.u[k + 1] - self.u[k])
        self.u_sn = u_hd
        self.log_1pzhel = np.log10(1.0 + self.sn["zHEL"])
        # vector de datos y covarianza total
        self.d = np.concatenate([self.d_bao, [c["theta_star_100_planck"], c["omega_m_planck"]], self.sn["mb"]])
        self.n_dat = self.d.size
        C = np.zeros((self.n_dat, self.n_dat))
        C[:13, :13] = self.C_bao
        C[13, 13] = c["sigma_theta_star_100"] ** 2
        C[14, 14] = c["sigma_omega_m_planck"] ** 2
        C[15:, 15:] = self.sn["cov"]
        self.C = C
        self.L = cholesky(C, lower=True)          # C = L L^T
        self.L_sn = cholesky(self.sn["cov"], lower=True)
        self.k_final = self.N - 1

    def _indice_nodo(self, u):
        k = int(np.argmin(np.abs(self.u - u)))
        assert abs(self.u[k] - u) < 1e-12, (u, self.u[k])
        return k

    # ---- fondo
    def E2(self, h, om, X=None):
        Om = om / h ** 2
        Or = self.Or_h2 / h ** 2
        z1 = 1.0 + self.z
        E2 = Om * z1 ** 3 + Or * z1 ** 4 + (1.0 - Om - Or)
        if X is not None:
            E2 = E2 + (X - X[0])
        return E2

    def distancias(self, h, om, X=None):
        """Devuelve E, D~_M, D~_H en todos los nodos (adimensionales)."""
        E2 = self.E2(h, om, X)
        if np.any(E2 <= 0):
            raise ValueError("E^2 <= 0 en algún nodo: historia no admisible")
        E = np.sqrt(E2)
        f = (1.0 + self.z) / E
        DM = np.zeros(self.N)
        DM[1:] = np.cumsum(0.5 * (f[:-1] + f[1:]) * self.du)
        DH = 1.0 / E
        return E, DM, DH

    def vector_modelo(self, p, X=None, partes=False):
        h, om, M = p
        E, DM, DH = self.distancias(h, om, X)
        cH = self.c_km_s / (100.0 * h)
        # BAO
        rd = float(self.tab.rd(om))
        y_bao = np.empty(13)
        for i, (k, tipo) in enumerate(zip(self.k_bao, self.tipo_bao)):
            if tipo == "DM":
                Dt = DM[k]
            elif tipo == "DH":
                Dt = DH[k]
            elif tipo == "DV":
                Dt = (self.z[k] * DM[k] ** 2 * DH[k]) ** (1.0 / 3.0)
            else:
                raise ValueError(tipo)
            y_bao[i] = cH / rd * Dt
        # theta*
        kf = self.k_final
        zs = float(self.tab.zs(om))
        rs = float(self.tab.rs(om))
        DM_star = DM[kf] + DH[kf] * (zs - self.z_ref)
        y_th = 100.0 * rs / (cH * DM_star) / self.CAL
        # omega_m
        y_om = om
        # SN
        k, t = self.k_sn, self.t_sn
        DM_sn = (1.0 - t) * DM[k] + t * DM[k + 1]
        y_sn = 5.0 * (self.log_1pzhel + np.log10(DM_sn)) + M
        y = np.concatenate([y_bao, [y_th, y_om], y_sn])
        if partes:
            return y, dict(E=E, DM=DM, DH=DH, DM_star=DM_star, DM_sn=DM_sn, rd=rd, rs=rs, zs=zs, y_bao=y_bao, y_th=y_th)
        return y

    # ---- chi2
    def residuo_blanqueado(self, p, X=None):
        r = self.d - self.vector_modelo(p, X)
        return solve_triangular(self.L, r, lower=True)

    def chi2(self, p, X=None):
        w = self.residuo_blanqueado(p, X)
        return float(w @ w)

    def M_analitico(self, h, om):
        """Mínimo analítico de chi2 en M (solo filas SN) para (h, om), X = 0."""
        y0 = self.vector_modelo((h, om, 0.0))
        r_sn = self.sn["mb"] - y0[15:]
        uno = np.ones(self.n_sn)
        Li_r = solve_triangular(self.L_sn, r_sn, lower=True)
        Li_1 = solve_triangular(self.L_sn, uno, lower=True)
        return float((Li_1 @ Li_r) / (Li_1 @ Li_1))
