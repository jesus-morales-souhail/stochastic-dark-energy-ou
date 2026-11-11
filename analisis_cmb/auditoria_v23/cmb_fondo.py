#!/usr/bin/env python3
"""
Exploration only (CMB_FONDO.md): does the compressed CMB bound what DESI alone cannot? Full linear-response OU kernel with
the background calibrated by Planck. Nothing in results/ or papers is changed.

Data vector (15 rows): the 13 DESI DR2 numbers (official covariance) + 100 theta_* + omega_m of Planck, with the values,
errors and the theta_* calibration of analisis_cmb/escalon/escalon.py, treatment (B) (imported from there).
Background: flat LCDM (h, omega_m) with radiation, r_d(omega_m) and r_s(z_*; omega_m) as in escalon.py.
Residual: stationary OU X(z) in x = ln(1+z), full linear response of verificacion.py (delta E^2 = X(z) - X(0)(1+z)^3,
flatness with Omega_L fixed, D_M integrated), extended to the theta_* row: delta ln theta_* = -delta ln D_M(z_*), X
integrated from z = 0 to z_*. Grid: uniform in u = ln(1+z), step 2e-3, from 0 to z_* (3500 nodes; trapezoid; data z on the
nearest node, |du| <= 1e-3); the z > 5 part of the theta_* response is measured (variance fraction). The omega_m row has no
response to X (omega_m is a parameter of the early universe; the residual lives in the dark-energy sector). Residual
covariance sigma_X^2 W K W^T on the 15 rows, K_ab = exp(-theta |u_a - u_b|), computed in O(N) per theta with the recursion of
the exponential kernel on the uniform grid.
Variant "Lambda absorbs X(0)": delta E^2 = X(z) - X(0) (the constant part of a DE residual is Lambda; the matter density is
omega_m). Then a constant X has no effect at all.
Background treatments: profiled (least squares over (h, omega_m) at every (sigma_X, theta)) and marginalised (linearised
around the LCDM best fit of each data set, flat prior), as verificacion.py / perfil_13.py.
Bounds on sigma_X with Delta lnL = -1.92 and -1.35, theta in {0.001, 0.1, 0.3, 1, 3, 10} and theta free; data sets: BAO
only, BAO + theta_*, BAO + omega_m, BAO + CMB (both). Check: BAO only, marginalised, must reproduce the full-kernel bounds of
VERIFICACION.md (0.475 at theta = 1).

    ../venv/bin/python cmb_fondo.py     # writes cmb_fondo.{json,txt}
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, least_squares, minimize_scalar
from scipy.signal import lfilter

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "escalon"))
import escalon as E  # noqa: E402  (data, Planck numbers, r_d / r_s splines, theta_* calibration, distances)
import verificacion as V  # noqa: E402  (bound(), SIG, P13.max_theta)

C_KMS = E.C_KMS
OBS13, COV13 = E.OBS, E.COV
TH_P, SIG_TH, WM_P, SIG_WM, OMR, ZSTAR = E.TH_P, E.SIG_TH, E.WM_P, E.SIG_WM, E.OMR, E.ZSTAR
Z13, KINDS13 = E.Z, E.CU.KINDS
THETAS = (0.001, 0.1, 0.3, 1.0, 3.0, 10.0)
THRS = {"192": -1.92, "135": -1.35}
OM_FID, H_FID = 1.0 - V.P13.OL_FID, E.H0_P / 100  # kernel fiducial: the S-table convention (Omega_L = 0.6847), Planck h
DU = 2e-3
Y15 = np.concatenate([OBS13, [TH_P, WM_P]])
CMEAS15 = np.zeros((15, 15))
CMEAS15[:13, :13] = COV13
CMEAS15[13, 13], CMEAS15[14, 14] = SIG_TH**2, SIG_WM**2
SETS = {"BAO": list(range(13)), "BAO+theta": list(range(14)), "BAO+omega_m": list(range(13)) + [14], "BAO+CMB": list(range(15))}


# ----------------------------------------------------------------------------------------------------------------------
# background model
# ----------------------------------------------------------------------------------------------------------------------
def model15(h: float, wm: float) -> np.ndarray:
    om, orad = wm / h**2, OMR / h**2
    d, _, _ = E.dist_lo(om, orad, 1.0, 1.0)
    bao = C_KMS / (100.0 * h * float(E.RD_SPL(wm))) * d
    return np.concatenate([bao, [E.theta_of_B(h, wm, 0.0, None), wm]])


# ----------------------------------------------------------------------------------------------------------------------
# full linear response on a uniform grid in u = ln(1+z), up to z_*
# ----------------------------------------------------------------------------------------------------------------------
class Response:
    def __init__(self, lambda_absorbs_x0: bool = False):
        u_star = np.log1p(ZSTAR)
        N = int(np.ceil(u_star / DU)) + 1
        self.u = np.linspace(0.0, u_star, N)
        self.z = np.expm1(self.u)
        self.du = self.u[1] - self.u[0]
        om, orad = OM_FID, OMR / H_FID**2
        E2 = om * (1 + self.z) ** 3 + orad * (1 + self.z) ** 4 + (1 - om - orad)
        Es = np.sqrt(E2)
        # trapezoid weights for int_0^{z_k} f dz = int_0^{u_k} f e^u du, cumulative
        base = self.du * (1 + self.z)  # dz = e^u du
        self.idx = [int(np.argmin(np.abs(self.u - np.log1p(z)))) for z in Z13]
        self.idx_star = N - 1
        self.max_du_data = float(max(abs(self.u[i] - np.log1p(z)) for i, z in zip(self.idx, Z13)))
        rows = []

        def cumw(k):
            w = np.zeros(N)
            w[: k + 1] = base[: k + 1]
            w[0] *= 0.5
            w[k] *= 0.5
            return w

        def wM(k):
            cw = cumw(k)
            dm = float(np.sum(cw / Es))
            w = -cw / (2 * Es**3) / dm
            x0 = float(np.sum(cw * (1 + self.z) ** 3 / (2 * Es**3))) / dm
            return w, x0, dm

        def wH(k):
            w = np.zeros(N)
            w[k] = -1.0 / (2 * E2[k])
            return w, (1 + self.z[k]) ** 3 / (2 * E2[k])

        for k, kind in zip(self.idx, KINDS13):
            m, m0, _ = wM(k)
            hh, h0 = wH(k)
            if kind == "DM":
                w, x0 = m, m0
            elif kind == "DH":
                w, x0 = hh, h0
            else:
                w, x0 = (2 * m + hh) / 3, (2 * m0 + h0) / 3
            rows.append((w, x0))
        m, m0, self.dm_star = wM(self.idx_star)
        rows.append((-m, -m0))  # delta ln theta_* = -delta ln D_M(z_*)
        rows.append((np.zeros(N), 0.0))  # omega_m row: no response
        self.Wx = np.array([r[0] for r in rows])  # response to X(z'), z' > 0 (local + integrated)
        x0 = np.array([r[1] for r in rows])
        self.W = self.Wx.copy()
        if lambda_absorbs_x0:
            self.W[:, 0] += -self.Wx.sum(axis=1)  # delta E^2 = X(z) - X(0)
        else:
            self.W[:, 0] += x0  # delta E^2 = X(z) - X(0)(1+z)^3
        self.coef_x0 = self.W[:, 0].copy()
        self.sumW = self.W.sum(axis=1)
        self.hi = self.z > 5.0

    def Kw(self, theta: float, W: np.ndarray) -> np.ndarray:
        """(K W^T)^T for the exponential kernel on the uniform grid: forward + backward recursions, O(N)."""
        rho = np.exp(-theta * self.du)
        F = lfilter([1.0], [1.0, -rho], W, axis=1)
        G = lfilter([1.0], [1.0, -rho], W[:, ::-1], axis=1)[:, ::-1] - W
        return F + G

    def M(self, theta: float, scale: np.ndarray) -> np.ndarray:
        Ws = self.W * scale[:, None]
        return Ws @ self.Kw(theta, Ws).T

    def theta_row_hi_fraction(self, theta: float) -> float:
        w = self.W[13]
        tot = float(w @ self.Kw(theta, w[None, :])[0])
        wh = np.where(self.hi, w, 0.0)
        hi = float(wh @ self.Kw(theta, wh[None, :])[0])
        return hi / tot


# ----------------------------------------------------------------------------------------------------------------------
# likelihoods per data set
# ----------------------------------------------------------------------------------------------------------------------
class Fit:
    def __init__(self, rows: list[int], R: Response):
        self.rows = np.array(rows)
        self.y = Y15[self.rows]
        self.C = CMEAS15[np.ix_(self.rows, self.rows)]
        self.R = R
        self.n = len(rows)
        # LCDM best fit
        r = least_squares(lambda p: self.white(self.C, self.y - model15(*p)[self.rows]), [0.68, 0.143], bounds=([0.5, 0.08], [0.9, 0.25]),
                          x_scale=[0.01, 0.001], ftol=1e-12, xtol=1e-12)
        self.p0 = r.x.copy()
        self.chi2_0 = float(2 * r.cost)
        m0 = model15(*self.p0)
        self.scale = m0  # raw-unit scaling of the response rows (delta y = y delta ln y); omega_m row: 0 response anyway
        # linear templates for (c)
        T = np.zeros((15, 2))
        for j, eps in ((0, 1e-4), (1, 1e-5)):
            p = self.p0.copy()
            p[j] += eps
            mp = model15(*p)
            p[j] -= 2 * eps
            mm = model15(*p)
            T[:, j] = (mp - mm) / (2 * eps)
        self.T = T[self.rows]
        self.r0 = self.y - m0[self.rows]
        self._cache = {}

    @staticmethod
    def white(C, r):
        L = np.linalg.cholesky(C)
        return np.linalg.solve(L, r)

    def cov(self, sigma: float, theta: float) -> np.ndarray:
        key = round(float(theta), 12)
        if key not in self._cache:
            if len(self._cache) > 5000:
                self._cache.clear()
            self._cache[key] = self.R.M(theta, self.scale)[np.ix_(self.rows, self.rows)]
        return self.C + sigma**2 * self._cache[key]

    @staticmethod
    def _ll(L, r):
        return float(-0.5 * (r @ np.linalg.solve(L.T, np.linalg.solve(L, r)) + 2 * np.sum(np.log(np.diag(L))) + len(r) * np.log(2 * np.pi)))

    def ll_b(self, sigma: float, theta: float) -> float:
        """profiled over (h, omega_m)"""
        L = np.linalg.cholesky(self.cov(sigma, theta))
        res = least_squares(lambda p: np.linalg.solve(L, self.y - model15(*p)[self.rows]), self.p0, bounds=([0.5, 0.08], [0.9, 0.25]),
                            x_scale=[0.01, 0.001], ftol=1e-10, xtol=1e-10)
        return self._ll(L, self.y - model15(*res.x)[self.rows])

    def ll_c(self, sigma: float, theta: float) -> float:
        """marginalised (linearised, flat prior)"""
        L = np.linalg.cholesky(self.cov(sigma, theta))
        solve = lambda b: np.linalg.solve(L.T, np.linalg.solve(L, b))
        cr, cT = solve(self.r0), solve(self.T)
        F = self.T.T @ cT
        b = self.T.T @ cr
        return float(-0.5 * (self.r0 @ cr - b @ np.linalg.solve(F, b) + 2 * np.sum(np.log(np.diag(L))) + np.linalg.slogdet(F)[1]
                             + self.n * np.log(2 * np.pi)))


def bounds_for(fit: Fit, treat: str) -> dict:
    fun = fit.ll_b if treat == "b" else fit.ll_c
    out = {}
    for th in THETAS + ("free",):
        g = (lambda s, f=fun: V.P13.max_theta(f, s)[0]) if th == "free" else (lambda s, f=fun, t=th: f(s, t))
        ll = np.array([g(float(s)) for s in V.SIG])
        rec = {"sigma_at_max": float(V.SIG[int(np.argmax(ll))]), "dlnL_max_vs_sigma0": float(ll.max() - ll[0])}
        for name, thr in THRS.items():
            rec[name] = V.bound(g, thr)["sigma_95"]
        out[f"theta_{th}" if th != "free" else "theta_free"] = rec
    return out


def main() -> None:
    t0 = time.time()
    R = Response()
    Rl = Response(lambda_absorbs_x0=True)
    out = {"grid": {"nodes": int(len(R.u)), "du": R.du, "max_du_data": R.max_du_data, "dm_star_dimless": R.dm_star},
           "theta_row_z_gt_5_variance_fraction": {f"theta_{t:g}": R.theta_row_hi_fraction(t) for t in (0.001, 0.1, 1.0, 10.0)},
           "response": {"rows": E.CU.G["labels"] + ["100theta_star", "omega_m"], "sum_W": R.sumW.tolist(), "coef_X0": R.coef_x0.tolist(),
                        "sum_W_lambda_variant": Rl.sumW.tolist(), "coef_X0_lambda_variant": Rl.coef_x0.tolist()},
           "planck": {"theta100": TH_P, "sigma": SIG_TH, "omega_m": WM_P, "sigma_omega_m": SIG_WM, "CAL": E.CAL}}
    # check: the S-table / s_type convention
    S13 = np.array([V.P13.s_type(k, z) for k, z in zip(KINDS13, Z13)])
    out["checks"] = {"max_abs_sumW_plus_s_type_13": float(np.max(np.abs(R.sumW[:13] + S13))),
                     "theta_row_sumW": float(R.sumW[13]), "theta_row_coef_X0": float(R.coef_x0[13])}
    print(json.dumps({k: out[k] for k in ("grid", "theta_row_z_gt_5_variance_fraction", "checks")}, indent=1), flush=True)
    out["fits"] = {}
    for name, rows in SETS.items():
        fit = Fit(rows, R)
        rec = {"lcdm": {"H0": 100 * fit.p0[0], "omega_m": float(fit.p0[1]), "chi2": fit.chi2_0}}
        for treat in ("b", "c"):
            t1 = time.time()
            rec["profiled" if treat == "b" else "marginalised"] = bounds_for(fit, treat)
            print(f"{name} {treat}: {json.dumps(rec['profiled' if treat == 'b' else 'marginalised'])}  [{time.time() - t1:.0f} s]", flush=True)
        out["fits"][name] = rec
    fitl = Fit(SETS["BAO+CMB"], Rl)
    out["fits"]["BAO+CMB|Lambda_absorbs_X0"] = {"lcdm": {"H0": 100 * fitl.p0[0], "omega_m": float(fitl.p0[1]), "chi2": fitl.chi2_0},
                                                "profiled": bounds_for(fitl, "b"), "marginalised": bounds_for(fitl, "c")}
    fitl2 = Fit(SETS["BAO"], Rl)
    out["fits"]["BAO|Lambda_absorbs_X0"] = {"marginalised": bounds_for(fitl2, "c")}
    out["runtime_s"] = time.time() - t0
    (HERE / "cmb_fondo.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    fmt = lambda v: f"{v:.4g}" if v is not None else "sin cota"
    keys = [f"theta_{t}" for t in THETAS] + ["theta_free"]
    L = [f"malla: {out['grid']['nodes']} nodos, du = {R.du:.4g}, |du| datos <= {R.max_du_data:.2e}; fraccion z > 5 de Var(dln theta*): "
         + ", ".join(f"theta={k[6:]}: {v:.2e}" for k, v in out["theta_row_z_gt_5_variance_fraction"].items()),
         f"comprobaciones: {json.dumps(out['checks'])}", ""]
    for name, rec in out["fits"].items():
        if "lcdm" in rec:
            L.append(f"{name}: LCDM H0 = {rec['lcdm']['H0']:.2f}, omega_m = {rec['lcdm']['omega_m']:.5f}, chi2 = {rec['lcdm']['chi2']:.2f}")
        for treat in ("profiled", "marginalised"):
            if treat not in rec:
                continue
            for thr in ("192", "135"):
                L.append(f"  {treat:13s} {thr}: " + "  ".join(f"{k[6:]}={fmt(rec[treat][k][thr])}" for k in keys))
            L.append("  " + " " * 13 + "sigma_at_max: " + "  ".join(f"{k[6:]}={rec[treat][k]['sigma_at_max']:.2g}({rec[treat][k]['dlnL_max_vs_sigma0']:+.2f})" for k in keys))
    (HERE / "cmb_fondo.txt").write_text("\n".join(L) + "\n", encoding="utf-8")
    print("\n".join(L))


if __name__ == "__main__":
    main()
