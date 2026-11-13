#!/usr/bin/env python3
"""
JACKKNIFE_BINS.md: does the preference for sigma_X > 0 (maximum of ln Z(sigma_X) - ln Z(0), linear analytic marginal,
convention B) live in a single DESI DR2 bin? Context: Kim, Mota & Tamosiunas, arXiv:2607.28918 v3 (22-10-2026), "An E-value
Reanalysis of DESI's Dynamical Dark Energy Signal" (w0wa, fixed Planck background: their e-value drops from 33.97 to 0.49
without LRG2). Our model and analysis are different (OU residual, refitted fiducial).

Set-up: that of marginal_guiada / respuesta_v3A (exacto grid G0, 3506 nodes with the DESI z exact; data = 13 DESI DR2 numbers
with the official 13x13 covariance + 100 theta_* + omega_m of Planck [+ Pantheon+ 1580 SN, STAT+SYS]; flat priors on
(h, omega_m[, M]); prior x ~ N(0, sigma^2 K_theta) with K the OU kernel on the grid). The response A is built with the
DEFINITION of respuesta_v3A.build_A_v3 (same code, exact-Jacobian SN rows, rows sum to zero) at the fiducial p0 refitted to
the kept rows with the ORIGINAL optimiser of each combination:
  - without SN: exacto.fit (least_squares trf, bounds [0.4, 1.0] x [0.05, 0.30], x_scale [0.01, 0.001], ftol = xtol = 1e-10,
    max_nfev 300, start (0.682, 0.1407)), on (h, omega_m);
  - with Pantheon+: exacto_sn.fit_combo (least_squares trf, bounds [0.4, 1.0] x [0.05, 0.30] x [-50, 50], x_scale
    [0.01, 0.001, 0.01], ftol = xtol = 1e-10, max_nfev 200, start (0.682, 0.1407, M analytic)), on (h, omega_m, M).
Removing a bin = removing its rows from the data vector and taking the principal submatrix of the covariance (of C and of
S = C + sigma^2 A K A^T): this is the marginal distribution of the kept data, NOT the conditional (Schur complement).
ln Z_lin(sigma) on the grid of maximos_globales (perfil_fondo.SIG up to 0.05, then geomspace(0.05, 1, 40)), ln Z(0) at
sigma = 0 exactly; maximum (value, node) over the grid; bound by a drop of 1.92 (and 1.35) from the global maximum, crossing
found with brentq between grid nodes (maximos_globales.bound_from). Checks at import: the all-bins A with Pantheon+ equals
respuesta_v3A.build_A_v3() bit for bit; the all-bins p0 equals marginal_guiada.P0 (with SN) and exacto_sn.P0 (without).

    ../venv/bin/python -I jackknife_bins.py paso1      # all bins only -> jackknife_bins_paso1.json
    ../venv/bin/python -I jackknife_bins.py todo       # all bins + 7 jackknife cases -> jackknife_bins.json
Nothing existing is modified.
"""
from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.linalg import cho_factor, solve_triangular
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import exacto_sn as S  # noqa: E402
import marginal_guiada as MG  # noqa: E402
import respuesta_v3A as RV  # noqa: E402
import maximos_globales as MGL  # noqa: E402  (GRID, SIG, bound_from; its main() is not run)

MODE = sys.argv[1] if len(sys.argv) > 1 else "paso1"
G, X, E = S.G0, S.X, S.E
SN = S.SN["all"]
NSN, NROWS = SN["n"], MG.NROWS
BINS = {"BGS": [0], "LRG1": [1, 2], "LRG2": [3, 4], "LRG3+ELG1": [5, 6], "ELG2": [7, 8], "QSO": [9, 10], "Lya": [11, 12]}
THETAS = (0.3, 1.0, 3.0, 10.0)
GRID = MGL.GRID
THRS = {"192": -1.92, "135": -1.35}
D = MG.D
CFULL = MG.C
LABELS = list(E.CU.G["labels"])
COV13 = X.COV13


# ----------------------------------------------------------------------------------------------------------------------
# fiducial refits (original optimisers, rows restricted)
# ----------------------------------------------------------------------------------------------------------------------
def fit_no_sn(keep_bao: list[int]) -> tuple[float, np.ndarray]:
    """exacto.fit on the 15-row set with only the BAO rows `keep_bao` (plus theta_*, omega_m)."""
    L = np.linalg.cholesky(COV13[np.ix_(keep_bao, keep_bao)])
    kb = np.array(keep_bao)

    def resid(p):  # exacto.resid, rows restricted (same model function exacto.model, same solves)
        v, m = X.model(p[0], p[1], np.zeros(X.N))
        if v is None:
            return np.full(len(kb) + 2, 1e3 * (1.0 - m))
        r = X.Y15 - v
        return np.concatenate([np.linalg.solve(L, r[kb]), [r[13] / X.SIG_TH, r[14] / X.SIG_WM]])

    r = least_squares(resid, np.asarray((0.682, 0.1407)), bounds=X.BOUNDS_P, method="trf", x_scale=[0.01, 0.001], ftol=1e-10, xtol=1e-10, max_nfev=300)
    return float(2.0 * r.cost), np.array([r.x[0], r.x[1], 0.0])


def fit_with_sn(keep_bao: list[int], x0=(0.682, 0.1407)) -> tuple[float, np.ndarray]:
    """exacto_sn.fit_combo(with_cmb=True) with only the BAO rows `keep_bao`."""
    L = np.linalg.cholesky(COV13[np.ix_(keep_bao, keep_bao)])
    kb = np.array(keep_bao)
    hist, _ = S.history(G, x0[1] / x0[0] ** 2, S.OMR / x0[0] ** 2, np.zeros(G.N))
    M0 = 0.0
    if hist is not None:
        rr = SN["mb"] - S.mu_sn(G, SN, hist[0])
        M0 = float((SN["Ci1"] @ rr) / SN["e"])

    def resid(p):
        h, wm, M = p
        hist, m = S.history(G, wm / h**2, S.OMR / h**2, np.zeros(G.N))
        nr = len(kb) + 2 + NSN
        if hist is None:
            return np.full(nr, 1e3 * (1.0 - m))
        dm_u, inv = hist
        r = X.Y15 - S.vec15(G, h, wm, dm_u, inv)
        rs = SN["mb"] - S.mu_sn(G, SN, dm_u) - M
        return np.concatenate([solve_triangular(L, r[kb], lower=True, check_finite=False), [r[13] / X.SIG_TH, r[14] / X.SIG_WM],
                               solve_triangular(SN["L"], rs, lower=True, check_finite=False)])  # as exacto_sn.resid_combo

    r = least_squares(resid, np.array([x0[0], x0[1], M0]), bounds=([0.4, 0.05, -50.0], [1.0, 0.30, 50.0]), method="trf",
                      x_scale=[0.01, 0.001, 0.01], ftol=1e-10, xtol=1e-10, max_nfev=200)
    return float(2.0 * r.cost), np.array(r.x)


# ----------------------------------------------------------------------------------------------------------------------
# linear pieces at a fiducial (build_A_v3 definition with p0, B0 injected)
# ----------------------------------------------------------------------------------------------------------------------
def model_at(p):
    return MG.exact_model(p[0], p[1], p[2], np.zeros((1, G.N)))[0][0]


def build_A_v3_at(p, B0full):
    """Same code as respuesta_v3A.build_A_v3, with (h, omega_m) and B0 given instead of marginal_guiada's."""
    E2, dm, dDM, Wdm, WdH = RV._ingredients(G, p[0], p[1], MG.OMR)
    rows = []
    for i, (k, kind) in enumerate(zip(G.IDX, X.KINDS13)):
        rows.append(Wdm[k] if kind == "DM" else (WdH[i] if kind == "DH" else (2 * Wdm[k] + WdH[i]) / 3))
    A = np.zeros((NROWS, G.N))
    A[:13] = np.array(rows) * B0full[:13, None]
    A[13] = -Wdm[-1] * B0full[13]
    A[15:] = RV._sn_rows_v3(dDM, dm, MG._k, MG._t, MG.MAG)
    return A


def build_B_at(p, npar):
    B = np.zeros((NROWS, 3))
    for j, eps in enumerate((1e-4, 1e-5, 1e-3)):
        pp, pm = np.array(p, float), np.array(p, float)
        pp[j] += eps
        pm[j] -= eps
        B[:, j] = (model_at(pp) - model_at(pm)) / (2 * eps)
    return B[:, :npar]


def profile_case(rows, p, Bm, A, K_by_theta):
    """ln Z_lin on GRID and at sigma = 0 for the kept rows; maximum and bounds per theta."""
    rows = np.array(rows)
    B0 = model_at(p)
    r = (D - B0)[rows]
    C = CFULL[np.ix_(rows, rows)]
    Ar, Br = A[rows], Bm[rows]
    out = {}
    for theta in THETAS:
        AK = Ar @ K_by_theta[theta]
        Msig = AK @ Ar.T
        Msig = 0.5 * (Msig + Msig.T)

        def lz(s):
            return MG.lnZ_lin(cho_factor(C + s * s * Msig, lower=True), Br, r)

        lz0 = lz(0.0)
        ll = np.array([lz(float(s)) for s in GRID])
        k = int(np.argmax(ll))
        rec = {"lnZ0": lz0, "lnZ_rel": (ll - lz0).tolist(), "max_lnZ_minus_lnZ0": float(ll[k] - lz0), "sigma_at_max": float(GRID[k]),
               "max_class": "sigma->0" if k == 0 else ("interior" if k < len(GRID) - 1 else "upper edge (1.0)"),
               "max_lnZ_minus_lnZ_grid0": float(ll[k] - ll[0])}
        for name, thr in THRS.items():
            rec[f"bound_{name}"] = MGL.bound_from(ll, lz, thr, False)
        out[f"theta_{theta:g}"] = rec
    return out


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    t0 = time.time()
    # ---- covariance structure of the 13 DESI numbers
    s = np.sqrt(np.diag(COV13))
    corr = COV13 / np.outer(s, s)
    offblock = np.ones((13, 13), bool)
    for b in BINS.values():
        for i in b:
            for j in b:
                offblock[i, j] = False
    cov_info = {"labels": LABELS, "bins": BINS, "z": X.Z13.tolist(), "within_bin_corr_DM_DH": {n: float(corr[b[0], b[1]]) for n, b in BINS.items() if len(b) == 2},
                "max_abs_cov_off_block": float(np.max(np.abs(COV13[offblock]))), "max_abs_corr_off_block": float(np.max(np.abs(corr[offblock]))),
                "n_off_block_nonzero": int(np.sum(COV13[offblock] != 0.0)), "sigma_over_obs_pct": (100 * s / X.OBS13).tolist()}
    # ---- checks of the set-up against the existing modules
    checks = {}
    c0_sn, p_sn = fit_with_sn(list(range(13)))
    checks["p0_with_sn_vs_marginal_guiada_P0_max_abs"] = float(np.max(np.abs(p_sn - MG.P0)))
    checks["chi2_with_sn_vs_marginal_guiada"] = [c0_sn, MG._c0]
    B0_sn = model_at(p_sn)
    A_sn = build_A_v3_at(p_sn, B0_sn)
    checks["A_all_bins_with_sn_vs_build_A_v3_max_abs"] = float(np.max(np.abs(A_sn - RV.build_A_v3())))
    checks["A_code_identity_at_marginal_guiada_P0_max_abs"] = float(np.max(np.abs(build_A_v3_at(MG.P0, MG.B0) - RV.build_A_v3())))
    checks["A_rows_sum_max_abs"] = float(np.max(np.abs(A_sn.sum(axis=1))))
    c0_no, p_no = fit_no_sn(list(range(13)))
    checks["p0_no_sn_vs_exacto_sn_P0_max_abs"] = float(np.max(np.abs(p_no[:2] - np.array(S.P0))))
    checks["chi2_no_sn"] = c0_no
    print(json.dumps(checks, indent=1), flush=True)
    assert checks["A_code_identity_at_marginal_guiada_P0_max_abs"] == 0.0, "build_A_v3_at must be build_A_v3 bit for bit at the same p0"
    assert checks["A_all_bins_with_sn_vs_build_A_v3_max_abs"] == 0.0, "A with all bins must equal build_A_v3 bit for bit"
    assert checks["p0_with_sn_vs_marginal_guiada_P0_max_abs"] == 0.0 and checks["p0_no_sn_vs_exacto_sn_P0_max_abs"] < 1e-12, "refits must reproduce the originals"
    # ---- kernels
    K = {th: MG.kernel(th) for th in THETAS}
    print(f"kernels: {time.time() - t0:.0f} s", flush=True)
    cases = {"todos": list(range(13))}
    if MODE == "todo":
        for n, b in BINS.items():
            cases[f"sin_{n}"] = [i for i in range(13) if i not in b]
    out = {"mode": MODE, "grid": {"n": len(GRID), "sigma_min": float(GRID[0]), "sigma_max": float(GRID[-1]), "n_below_0.05": int(len(MGL.SIG))},
           "thetas": THETAS, "covariance_13": cov_info, "checks": checks, "combos": {}}
    for combo in ("sinSN", "Pantheon+"):
        out["combos"][combo] = {}
        for case, keep in cases.items():
            t1 = time.time()
            if combo == "sinSN":
                c0, p = fit_no_sn(keep)
                rows = keep + [13, 14]
                npar = 2
            else:
                c0, p = fit_with_sn(keep)
                rows = keep + [13, 14] + list(range(15, NROWS))
                npar = 3
            B0 = model_at(p)
            A = build_A_v3_at(p, B0)
            Bm = build_B_at(p, npar)
            rec = {"rows_kept_bao": keep, "n_rows": len(rows), "chi2_lcdm": c0, "p0": {"h": float(p[0]), "omega_m": float(p[1]), "M": float(p[2])},
                   "thetas": profile_case(rows, p, Bm, A, K)}
            out["combos"][combo][case] = rec
            print(f"{combo:9s} {case:14s} chi2={c0:8.3f} h={p[0]:.5f} wm={p[1]:.5f}  " + "  ".join(
                f"th{th:g}: max {rec['thetas'][f'theta_{th:g}']['max_lnZ_minus_lnZ0']:+.3f}@{rec['thetas'][f'theta_{th:g}']['sigma_at_max']:.3g} "
                f"b192={rec['thetas'][f'theta_{th:g}']['bound_192']}" for th in THETAS) + f"  [{time.time() - t1:.0f} s]", flush=True)
    # ---- references for step 1
    mgj = json.load(open(HERE / "maximos_globales.json"))
    psj = json.load(open(HERE / "precision_sn.json"))["table"]
    v3a = json.load(open(HERE / "exacto_sn_lineal_v3A_maximos_v3A.json"))["combos"]["BAO+CMB+SN"]
    pl = json.load(open(HERE / "positividad_lin_v3A.json"))
    refs = {"sinSN_linear_maximos_globales (cmb_fondo grid, fiducial OM_FID)": {f"theta_{t:g}": {k: mgj["cmb_fondo"]["BAO+CMB|Lambda_absorbs_X0"]["marginalised"][f"theta_{float(t)}"][k]
                                                                               for k in ("sigma_at_max_grid", "lnL_max_minus_sigma0", "192_global", "135_global")} for t in THETAS},
            "sinSN_exact_MC_precision_sn (2000 paths)": {f"theta_{t:g}": {k: psj[f"theta_{t:g}|BAO+CMB"][k] for k in ("lnL_max_at_sigma", "lnL_max", "marg_192", "marg_135")} for t in THETAS},
            "Pantheon+_linear_exacto_sn_lineal_v3A_maximos (cmb_fondo grid, fiducial OM_FID)": {f"theta_{t:g}": {k: v3a[f"theta_{t:g}"][k] for k in ("sigma_at_max_grid", "lnL_max_minus_sigma0", "192_global", "135_global")} for t in THETAS},
            "Pantheon+_linear_positividad_lin_v3A (build_A_v3, grid G, its own sigma cells)": {f"theta_{t:g}": {"bound_lin_192": pl["thetas"][f"theta_{t:g}"]["bound_lin_192"], "bound_lin_135": pl["thetas"][f"theta_{t:g}"]["bound_lin_135"],
                                                                                                  "max_lnZ_lin_rel": float(max(pl["thetas"][f"theta_{t:g}"]["lnZ_lin_rel"])), "max_lin_at_sigma": pl["thetas"][f"theta_{t:g}"]["max_lin_at_sigma"]} for t in THETAS}}
    out["references"] = refs
    # ---- jackknife summary
    if MODE == "todo":
        summ = {}
        for combo in out["combos"]:
            summ[combo] = {}
            for th in THETAS:
                key = f"theta_{th:g}"
                base = out["combos"][combo]["todos"]["thetas"][key]
                d = {}
                for n in BINS:
                    r = out["combos"][combo][f"sin_{n}"]["thetas"][key]
                    d[n] = {"max": r["max_lnZ_minus_lnZ0"], "d_max": r["max_lnZ_minus_lnZ0"] - base["max_lnZ_minus_lnZ0"], "sigma_at_max": r["sigma_at_max"],
                            "bound_192": r["bound_192"], "d_bound_192_rel": (r["bound_192"] / base["bound_192"] - 1) if (r["bound_192"] and base["bound_192"]) else None}
                summ[combo][key] = {"todos": {"max": base["max_lnZ_minus_lnZ0"], "sigma_at_max": base["sigma_at_max"], "bound_192": base["bound_192"]}, "sin_bin": d,
                                    "bin_moving_max_most": max(d, key=lambda n: abs(d[n]["d_max"])),
                                    "bin_moving_bound_most": max((n for n in d if d[n]["d_bound_192_rel"] is not None), key=lambda n: abs(d[n]["d_bound_192_rel"]), default=None)}
        out["summary"] = summ
        s1 = summ["sinSN"]["theta_1"]
        others = {n: v["d_max"] for n, v in s1["sin_bin"].items() if n != "LRG2"}
        p14 = {n: v["d_max"] for n, v in summ["Pantheon+"]["theta_1"]["sin_bin"].items()}
        out["predictions"] = {"P12": {"statement": "sin SN, theta=1: sin LRG2 el maximo queda por debajo de la mitad del de todos los bins",
                                      "max_all": s1["todos"]["max"], "max_sin_LRG2": s1["sin_bin"]["LRG2"]["max"], "holds": bool(s1["sin_bin"]["LRG2"]["max"] < 0.5 * s1["todos"]["max"])},
                              "P13": {"statement": "sin SN, theta=1: quitar cualquier otro bin cambia el maximo menos de +-0.5", "d_max_others": others,
                                      "holds": bool(all(abs(v) < 0.5 for v in others.values()))},
                              "P14": {"statement": "con Pantheon+, theta=1: quitar cualquier bin cambia el maximo menos de +-0.5", "d_max": p14,
                                      "holds": bool(all(abs(v) < 0.5 for v in p14.values()))}}
        print(json.dumps(out["predictions"], indent=1), flush=True)
    out["runtime_s"] = time.time() - t0
    name = "jackknife_bins.json" if MODE == "todo" else "jackknife_bins_paso1.json"
    (HERE / name).write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")
    print("escrito", name, f"{out['runtime_s']:.0f} s")


if __name__ == "__main__":
    main()
