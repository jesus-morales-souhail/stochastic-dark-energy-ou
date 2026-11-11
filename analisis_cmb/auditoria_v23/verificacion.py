#!/usr/bin/env python3
"""
Exploration only (VERIFICACION.md): three checks prompted by an external audit of PREPRINT v2.2.1. Nothing in results/ or
papers is changed.

1. Tachyonic rank-1 mode (C4). scripts/tachyonic_rank1_mle.py builds u_i = sigma_0 S_i exp(+t_i/t_c) with t = lookback time
   (larger in the past), sigma_0 = 1e-61 fixed, and reports as "representative" the scan point closest to Delta lnL = -11.35
   (line 97). Here: (i) the scan is reproduced and the neighbours of that point listed; (ii) the amplitude ratio between the
   last and the first bin is shown; then the mode that grows in cosmic time, u_i = A_0 S_i exp(-t_i/t_c) (A_0 = amplitude today,
   in alpha units), is profiled over t_c and the 95 % limit on A_0 is computed (Delta lnL = -1.92 and -1.35), per t_c and with
   t_c profiled; compared with the OU bound at theta -> 0 (constant mode, t_c -> infinity).
2. Full linear response of the OU kernel. The pipeline uses delta alpha_i = S_i X(z_i) (local). Full response, flat, H0 fixed,
   fiducial Omega_L = 0.6847 (the S_z convention): delta E^2(z) = X(z) - X(0)(1+z)^3, delta ln D_H = -delta E^2/(2E^2),
   delta ln D_M(z) = -(1/Dhat_M) int_0^z delta E^2/(2E^3) dz', delta ln D_V = (2/3) delta ln D_M + (1/3) delta ln D_H.
   Written as delta ln D_i = sum_a W_ia X(z_a) on a Simpson grid (nodes at z = 0 and the data z), the OU covariance propagates
   as C_ij = sigma_X^2 (W K W^T)_ij, K_ab = exp(-theta |x_a - x_b|), x = ln(1+z). Check: sum_a W_ia = d ln D_i / d Omega_L
   (= -S_i) and C -> S S^T sigma_X^2 at theta -> 0. Bounds (fixed / profiled / marginalised background, perfil_fondo.py and
   perfil_13.py machinery) at theta = 0.001, 0.1, 1, 3, 10 and theta free, with the local and the full kernel, 7 alpha and 13 numbers.
3. Threshold. Main bounds (7 alpha, local kernel) with Delta lnL = -1.92 and -1.35; coverage by injection: 2000 universes
   alpha = 1 + S X + n, n ~ N(0, C_7x7) (seed 20261105), X an OU draw at (sigma_t, theta_t) with sigma_t = each bound and
   theta_t in {0.001, 1}; coverage = fraction of universes whose limit (same treatment, theta fixed at theta_t or free) is
   >= sigma_t, i.e. lnL_prof(sigma_t) - max lnL_prof >= threshold. Also the limits of noise-only universes (median, 16-84 %).

    ../venv/bin/python verificacion.py     # writes verificacion.{json,txt}
"""

from __future__ import annotations

import json
import sys
import time
from functools import lru_cache
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq, minimize_scalar

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE.parent / "fondo"))
sys.path.insert(0, str(ROOT / "scripts"))
import perfil_13 as P13  # noqa: E402
import perfil_fondo as PF  # noqa: E402
from desi_dr2_data import load_alpha_dv, load_gaussian_bao_vector, logL_gaussian, measurement_cov  # noqa: E402

Z7, ALPHA7, S7, CM7, X7 = PF.Z, PF.ALPHA, PF.S, PF.CM, PF.X
SIG = PF.SIG
THR_192, THR_135 = -1.92, -1.35
OM_FID = 1.0 - P13.OL_FID  # 0.3153, the fiducial of the S_z table
THETAS = (0.001, 0.1, 1.0, 3.0, 10.0)
N_UNI, SEED_NOISE = 2000, 20261105
C_KMS = 299792.458


# ----------------------------------------------------------------------------------------------------------------------
# generic bound with a threshold argument (perfil_fondo.bound with THR as parameter)
# ----------------------------------------------------------------------------------------------------------------------
def bound(g, thr: float) -> dict:
    ll = np.array([g(float(s)) for s in SIG])
    llmax = float(ll.max())
    rel = ll - llmax
    out = {"sigma_at_max": float(SIG[int(np.argmax(ll))])}
    j = int(np.nonzero(rel >= thr)[0][-1])
    if j + 1 < len(SIG):
        out.update(sigma_95=float(brentq(lambda s: g(float(s)) - llmax - thr, SIG[j], SIG[j + 1], xtol=1e-10 * SIG[j])), crossing="inside grid")
        return out
    ext = np.geomspace(SIG[-1], 1.0, 25)[1:]
    lle = np.array([g(float(s)) for s in ext])
    k = np.nonzero(lle - llmax < thr)[0]
    if len(k):
        k = int(k[0])
        lo = SIG[-1] if k == 0 else ext[k - 1]
        out.update(sigma_95=float(brentq(lambda s: g(float(s)) - llmax - thr, lo, ext[k], xtol=1e-10 * lo)), crossing="extended to 1")
    else:
        out.update(sigma_95=None, crossing="no crossing up to sigma_X = 1")
    return out


# ----------------------------------------------------------------------------------------------------------------------
# 1. tachyonic
# ----------------------------------------------------------------------------------------------------------------------
H0_GYR = 67.4 * 1.022712165e-3


def lookback_gyr(z: float, om: float = 0.315) -> float:
    return quad(lambda zp: 1.0 / ((1.0 + zp) * np.sqrt(om * (1 + zp) ** 3 + 1 - om)), 0.0, z, epsabs=1e-10)[0] / H0_GYR


def tachyonic() -> dict:
    d = load_alpha_dv(prefer_file=True, use_full_cov=True)
    z, alpha, S = d["z"], d["alpha"], d["S_z"]
    C = measurement_cov(d)
    t = np.array([lookback_gyr(float(zi)) for zi in z])
    r = alpha - 1.0
    ll0 = logL_gaussian(r, C)
    # (i) reproduce the scan of the script (sigma_0 = 1e-61, exp(+t/t_c))
    tcs = np.logspace(-4, 4, 500)
    rows = []
    for tc in tcs:
        u = 1e-61 * S * np.exp(np.clip(t / tc, -700, 700))
        ll = logL_gaussian(r, C + np.outer(u, u)) if np.max(np.abs(u)) <= 1e12 else -np.inf
        rows.append((tc, ll - ll0 if np.isfinite(ll) else np.nan, float(np.max(np.abs(u)))))
    rows = np.array(rows)
    fin = rows[np.isfinite(rows[:, 1])]
    act = fin[fin[:, 2] > 1e-6]
    i_near = int(np.argmin(np.abs(act[:, 1] + 11.35)))
    i_best = int(np.argmax(act[:, 1]))
    near = act[i_near]
    j = int(np.nonzero(fin[:, 0] == near[0])[0][0])
    neigh = [{"tc_Gyr": float(fin[k, 0]), "dlnL": float(fin[k, 1]), "max_amp_alpha": float(fin[k, 2])} for k in range(max(j - 3, 0), min(j + 4, len(fin)))]
    # (ii) direction in time
    ratio_past_today = float(np.exp((t[-1] - t[0]) / near[0]))
    # correct sign: growing in cosmic time, amplitude today A_0 (alpha units)
    def ll_grow(A0: float, tc: float) -> float:
        u = A0 * S * np.exp(-t / tc)
        return logL_gaussian(r, C + np.outer(u, u))

    def ll_prof_tc(A0: float, tc_min: float = 0.03) -> tuple[float, float]:
        best, best_tc = -1e30, np.nan
        for tc in np.geomspace(tc_min, 1e3, 120):
            v = ll_grow(A0, tc)
            if v > best:
                best, best_tc = v, tc
        return best, float(best_tc)

    out = {"ll_LCDM": ll0, "scan_near_minus_11.35": {"tc_Gyr": float(near[0]), "dlnL": float(near[1]), "max_amp_alpha": float(near[2])},
           "scan_best_active": {"tc_Gyr": float(act[i_best, 0]), "dlnL": float(act[i_best, 1]), "max_amp_alpha": float(act[i_best, 2])},
           "scan_neighbours": neigh, "lookback_Gyr": t.tolist(),
           "amp_ratio_z2.33_over_z0.295_at_that_tc": ratio_past_today,
           "tc_where_sigma0_exp_tmax_tc_equals_0.03": float(t[-1] / np.log(0.03 / 1e-61)),
           "growing_mode": {"per_tc": {}, "profiled_over_tc": {}}}
    sv = float(1 / np.sqrt(S @ np.linalg.solve(C, S)))
    out["ou_theta0_reference"] = {"sigma_v": sv, "occam_bound_192": float(np.sqrt(np.exp(3.84) - 1) * sv),
                                  "occam_bound_135": float(np.sqrt(np.exp(2.70) - 1) * sv), "C2_bound": 0.02633}
    for tc in (0.1, 0.3, 1.0, 3.0, 10.0, 100.0, np.inf):
        g = lambda a, tc=tc: ll_grow(a, tc) if np.isfinite(tc) else logL_gaussian(r, C + a**2 * np.outer(S, S))
        rec = {}
        for name, thr in (("192", THR_192), ("135", THR_135)):
            b = bound(g, thr)
            rec[name] = b["sigma_95"]
        rec["dlnL_max"] = float(max(g(float(s)) for s in SIG) - ll0)
        rec["amp_at_z2.33_over_A0"] = float(np.exp(-t[-1] / tc)) if np.isfinite(tc) else 1.0
        out["growing_mode"]["per_tc"][f"tc_{tc:g}"] = rec
    for tc_min, lab in ((0.03, "tc_free"), (1.0, "tc_ge_1Gyr")):
        gp = lambda a, m=tc_min: ll_prof_tc(a, m)[0]
        rec = {"dlnL_max": float(max(gp(float(s)) for s in SIG) - ll0), "dlnL_at_A0_1": float(gp(1.0) - ll0)}
        for name, thr in (("192", THR_192), ("135", THR_135)):
            b = bound(gp, thr)
            rec[name] = {"A0_95": b["sigma_95"], "crossing": b["crossing"],
                         "tc_at_bound": ll_prof_tc(b["sigma_95"], tc_min)[1] if b["sigma_95"] is not None else None}
        out["growing_mode"]["profiled_over_tc"][lab] = rec
    return out


# ----------------------------------------------------------------------------------------------------------------------
# 2. full linear response
# ----------------------------------------------------------------------------------------------------------------------
NSUB = 128


class Response:
    """W (n_obs x N) with delta ln D_i = sum_a W_ia X(z_a); nodes: z = 0 and Simpson segments between the data z."""

    def __init__(self, z_rows, kinds, om: float = OM_FID):
        zu = np.unique(np.asarray(z_rows, float))
        edges = np.concatenate([[0.0], zu])
        nodes, wts, seg_end = [0.0], [], {}
        for a, b in zip(edges[:-1], edges[1:]):
            zz = np.linspace(a, b, NSUB + 1)
            w = np.ones(NSUB + 1)
            w[1:-1:2], w[2:-1:2] = 4.0, 2.0
            w *= (b - a) / NSUB / 3.0
            nodes.extend(zz[1:].tolist())
            wts.append(w)
            seg_end[float(b)] = len(nodes) - 1
        self.z = np.array(nodes)
        self.x = np.log1p(self.z)
        E2 = om * (1 + self.z) ** 3 + (1 - om)
        E = np.sqrt(E2)
        N = len(self.z)
        # cumulative Simpson weights for int_0^{z_k} ... at each unique z (segment end)
        cumw = {}
        acc = np.zeros(N)
        start = 0
        for k, (b, w) in enumerate(zip(zu, wts)):
            idx = np.arange(start, start + NSUB + 1)
            acc = acc.copy()
            acc[idx] += w
            cumw[float(b)] = acc.copy()
            start += NSUB
        self.W = np.zeros((len(z_rows), N))
        self.sum_check = []
        for i, (zi, k) in enumerate(zip(z_rows, kinds)):
            zi = float(zi)
            j = seg_end[zi]
            dm = float(np.sum(cumw[zi] / E))
            wH = np.zeros(N)
            wH[j] += -1.0 / (2 * E2[j])
            wH[0] += (1 + zi) ** 3 / (2 * E2[j])
            wM = -cumw[zi] / (2 * E**3) / dm
            wM[0] += float(np.sum(cumw[zi] * (1 + self.z) ** 3 / (2 * E**3))) / dm
            self.W[i] = {"DH": wH, "DM": wM, "DV": (2 * wM + wH) / 3}[k]
        self.dx = np.abs(self.x[:, None] - self.x[None, :])
        self._cache = {}

    def M(self, theta: float) -> np.ndarray:
        key = round(float(theta), 12)
        if key not in self._cache:
            if len(self._cache) > 4000:
                self._cache.clear()
            K = np.exp(-theta * self.dx)
            self._cache[key] = self.W @ K @ self.W.T
        return self._cache[key]


def treatments(D: P13.Data, kern, thetas=THETAS, thr: float = THR_192, theta_free: bool = True) -> dict:
    """Bounds (a) fixed, (b) profiled, (c) marginalised with a given kernel function kern(sigma, theta) -> covariance."""
    D.kern = kern
    out = {}
    for name, fun in (("a_fixed", D.ll_a), ("b_profiled", D.ll_b), ("c_marginalised", D.ll_c)):
        d = {}
        if theta_free:
            d["theta_free"] = bound(lambda s, f=fun: P13.max_theta(f, s)[0], thr)["sigma_95"]
        for th in thetas:
            d[f"theta_{th:g}"] = bound(lambda s, f=fun, t=th: f(s, t), thr)["sigma_95"]
        out[name] = d
    return out


def data13() -> tuple[P13.Data, np.ndarray, list, np.ndarray]:
    g = load_gaussian_bao_vector()
    kinds = [lab[:2] for lab in g["labels"]]
    z = np.asarray(g["z_rows"], float)
    obs, cov = np.asarray(g["mean"], float), np.asarray(g["cov"], float)
    pred = np.array([C_KMS / (100 * PF.HRD_F) * P13.dhat(k, PF.OM_F, zz) for k, zz in zip(kinds, z)])
    S13 = np.array([P13.s_type(k, zz) for k, zz in zip(kinds, z)])
    return P13.Data(z, kinds, obs / pred, cov / np.outer(pred, pred), S13), z, kinds, S13


def full_kernel() -> dict:
    out = {}
    D7 = P13.Data(Z7, ["DV"] * 7, ALPHA7, CM7, S7)
    D13, z13, kinds13, S13 = data13()
    R7, R13 = Response(Z7, ["DV"] * 7), Response(z13, kinds13)
    out["checks"] = {
        "7alpha": {"sum_W": R7.W.sum(axis=1).tolist(), "minus_S_table": (-S7).tolist(),
                   "max_abs_diff": float(np.max(np.abs(R7.W.sum(axis=1) + S7))),
                   "max_rel_diff_cov_theta0.001": float(np.max(np.abs(R7.M(0.001) / np.outer(S7, S7) - 1)))},
        "13": {"max_abs_diff_sum_W_vs_minus_s_type": float(np.max(np.abs(R13.W.sum(axis=1) + S13))),
               "max_rel_diff_cov_theta0.001": float(np.max(np.abs(R13.M(0.001) / np.outer(S13, S13) - 1)))},
        "nodes": int(len(R7.z))}
    # variance ratio of the full response at each bin against the local one, theta = 1 and 10
    for th in (1.0, 3.0, 10.0):
        out["checks"][f"diag_ratio_full_over_local_7alpha_theta_{th:g}"] = (np.diag(R7.M(th)) / S7**2).tolist()
    for lab, D, R in (("7alpha", D7, R7), ("13", D13, R13)):
        loc = lambda s, t, D=D: D.cov + np.outer(D.S, D.S) * s**2 * np.exp(-t * D.dx)
        ful = lambda s, t, R=R, D=D: D.cov + s**2 * R.M(t)
        out[lab] = {"local": treatments(D, loc), "full": treatments(D, ful)}
        out[lab]["ratio_full_over_local"] = {tr: {k: (None if out[lab]["full"][tr][k] is None or out[lab]["local"][tr][k] is None
                                                      else out[lab]["full"][tr][k] / out[lab]["local"][tr][k])
                                                  for k in out[lab]["local"][tr]} for tr in out[lab]["local"]}
        print(lab, json.dumps(out[lab]), flush=True)
    return out


# ----------------------------------------------------------------------------------------------------------------------
# 3. threshold and coverage (7 alpha, local kernel)
# ----------------------------------------------------------------------------------------------------------------------
def kern7(s, t):
    return CM7 + np.outer(S7, S7) * s**2 * np.exp(-t * np.abs(X7[:, None] - X7[None, :]))


def ll_fixed(y, s, t):
    return logL_gaussian(y - 1.0, kern7(s, t))


def prof_free(y, s):
    return P13.max_theta(lambda ss, tt: ll_fixed(y, ss, tt), s)[0]


_U = {}


def _worker(i: int):
    n = _U["noise"][i]
    res = {}
    # noise-only limits
    y = 1.0 + n
    for lab, g in (("theta_free", lambda s: prof_free(y, s)), ("theta_0.001", lambda s: ll_fixed(y, s, 0.001)), ("theta_1", lambda s: ll_fixed(y, s, 1.0))):
        ll = np.array([g(float(s)) for s in SIG])
        llmax = ll.max()
        for name, thr in (("192", THR_192), ("135", THR_135)):
            res[f"noise|{lab}|{name}"] = bound(g, thr)["sigma_95"] or 1.0
    # injected: coverage indicator = lnL_prof(sigma_t) - max_sigma lnL_prof >= thr
    for case, (th_t, s_t, treat) in _U["cases"].items():
        X = _U["L_inj"][case] @ _U["gauss"][i]
        y = 1.0 + S7 * X + n
        g = (lambda s: prof_free(y, s)) if treat == "free" else (lambda s, t=th_t: ll_fixed(y, s, t))
        ll = np.array([g(float(s)) for s in SIG] + [g(float(s_t))])
        res[f"cov|{case}"] = float(ll[-1] - ll[:-1].max())
    return res


def coverage(bounds_main: dict) -> dict:
    rng = np.random.default_rng(SEED_NOISE)
    noise = rng.standard_normal((N_UNI, 7)) @ np.linalg.cholesky(CM7).T
    gauss = np.random.default_rng(SEED_NOISE + 1).standard_normal((N_UNI, 7))
    cases, L_inj = {}, {}
    for th_t, key in ((0.001, "theta_0.001"), (1.0, "theta_1")):
        for name in ("192", "135"):
            for treat in ("fixed", "free"):
                s_t = bounds_main[name]["a_fixed"][key if treat == "fixed" else "theta_free"]
                case = f"{key}|{name}|{treat}"
                cases[case] = (th_t, s_t, treat)
                K = s_t**2 * np.exp(-th_t * np.abs(X7[:, None] - X7[None, :]))
                L_inj[case] = np.linalg.cholesky(K + 1e-12 * s_t**2 * np.eye(7))
    _U.update(noise=noise, gauss=gauss, cases=cases, L_inj=L_inj)
    with Pool() as pool:
        R = pool.map(_worker, range(N_UNI), chunksize=20)
    out = {"n": N_UNI, "cases": {}, "noise_only_limits": {}}
    for case, (th_t, s_t, treat) in cases.items():
        d = np.array([r[f"cov|{case}"] for r in R])
        thr = THR_192 if "|192|" in case else THR_135
        c = float(np.mean(d >= thr))
        out["cases"][case] = {"theta_true": th_t, "sigma_true": s_t, "treatment": treat, "threshold": thr,
                              "coverage": c, "coverage_err": float(np.sqrt(c * (1 - c) / N_UNI))}
    for key in R[0]:
        if key.startswith("noise|"):
            v = np.array([r[key] for r in R])
            out["noise_only_limits"][key[6:]] = {"median": float(np.median(v)), "p16": float(np.percentile(v, 16)), "p84": float(np.percentile(v, 84)),
                                                 "frac_at_extension_edge": float(np.mean(v >= 1.0))}
    return out


def main() -> None:
    t0 = time.time()
    out = {"tachyonic": tachyonic()}
    print("tachyonic:", json.dumps(out["tachyonic"], indent=1), flush=True)
    out["full_kernel"] = full_kernel()
    print(f"[{time.time() - t0:.0f} s]", flush=True)
    D7 = P13.Data(Z7, ["DV"] * 7, ALPHA7, CM7, S7)
    loc = lambda s, t: D7.cov + np.outer(D7.S, D7.S) * s**2 * np.exp(-t * D7.dx)
    out["threshold"] = {"192": treatments(D7, loc, thetas=(0.001, 1.0, 10.0), thr=THR_192),
                        "135": treatments(D7, loc, thetas=(0.001, 1.0, 10.0), thr=THR_135)}
    print("threshold:", json.dumps(out["threshold"]), flush=True)
    out["coverage"] = coverage(out["threshold"])
    print("coverage:", json.dumps(out["coverage"], indent=1), flush=True)
    out["runtime_s"] = time.time() - t0
    (HERE / "verificacion.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    summary(out)


def summary(out: dict) -> None:
    L = []
    T = out["tachyonic"]
    L.append(f"1. taquionico: punto 'representativo' del barrido: tc = {T['scan_near_minus_11.35']['tc_Gyr']:.4e} Gyr, dlnL = {T['scan_near_minus_11.35']['dlnL']:.3f}, "
             f"amplitud max = {T['scan_near_minus_11.35']['max_amp_alpha']:.3g} (en alpha); mejor activo dlnL = {T['scan_best_active']['dlnL']:.2e}")
    L.append("   vecinos: " + "; ".join(f"tc={n['tc_Gyr']:.4f}: {n['dlnL']:+.3f} (amp {n['max_amp_alpha']:.2g})" for n in T["scan_neighbours"]))
    L.append(f"   cociente de amplitud z=2.33 / z=0.295 en ese tc: {T['amp_ratio_z2.33_over_z0.295_at_that_tc']:.3g} (crece hacia el pasado)")
    fz = lambda v: f"{v:.4g}" if v is not None else "sin cota"
    L.append("   modo que crece en tiempo cosmico, cota sobre A_0 (alpha) al 95 % (-1.92 / -1.35): " + "; ".join(
        f"tc={k[3:]}: {fz(v['192'])} / {fz(v['135'])}" for k, v in T["growing_mode"]["per_tc"].items()))
    for lab, P in T["growing_mode"]["profiled_over_tc"].items():
        f = lambda v: f"{v:.4g}" if v is not None else "sin cota"
        L.append(f"   perfilando tc ({lab}): {f(P['192']['A0_95'])} (-1.92, tc en la cota = {P['192']['tc_at_bound']}) / {f(P['135']['A0_95'])} (-1.35); "
                 f"max dlnL = {P['dlnL_max']:.2e}; dlnL en A0 = 1: {P['dlnL_at_A0_1']:.3f}")
    L.append(f"   OU theta->0: {T['ou_theta0_reference']['C2_bound']} (Occam {T['ou_theta0_reference']['occam_bound_192']:.4g} / {T['ou_theta0_reference']['occam_bound_135']:.4g})")
    F = out["full_kernel"]
    L.append(f"2. nucleo completo: nodos {F['checks']['nodes']}; 7alpha max|sum W + S| = {F['checks']['7alpha']['max_abs_diff']:.2e}, "
             f"cov rel diff theta=0.001: {F['checks']['7alpha']['max_rel_diff_cov_theta0.001']:.2e}; 13: {F['checks']['13']['max_abs_diff_sum_W_vs_minus_s_type']:.2e}, {F['checks']['13']['max_rel_diff_cov_theta0.001']:.2e}")
    keys = ["theta_free"] + [f"theta_{t:g}" for t in THETAS]
    for lab in ("7alpha", "13"):
        L.append(f"   {lab}: {'tratamiento':16s} " + " ".join(f"{k:>22s}" for k in keys))
        for tr in ("a_fixed", "b_profiled", "c_marginalised"):
            fmt = lambda v: f"{v:.4g}" if v is not None else "sin cota"
            L.append(f"   {'':8s} {tr:16s} " + " ".join(f"{fmt(F[lab]['local'][tr][k]):>7s}/{fmt(F[lab]['full'][tr][k]):>7s}"
                                                       f"({(F[lab]['ratio_full_over_local'][tr][k] or float('nan')):.2f})" for k in keys))
    L.append("3. umbral (7 alpha, local): " + json.dumps(out["threshold"]))
    for case, c in out["coverage"]["cases"].items():
        L.append(f"   cobertura {case}: sigma_t = {c['sigma_true']:.4g}, cobertura = {c['coverage']:.4f} +- {c['coverage_err']:.4f}")
    for k, v in out["coverage"]["noise_only_limits"].items():
        L.append(f"   limite con solo ruido {k}: mediana {v['median']:.4g} (16-84 %: {v['p16']:.4g}-{v['p84']:.4g})")
    (HERE / "verificacion.txt").write_text("\n".join(L) + "\n", encoding="utf-8")
    print("\n".join(L))


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--txt-only":
        summary(json.loads((HERE / "verificacion.json").read_text()))
    else:
        main()
