#!/usr/bin/env python3
"""
_v3A versions of exacto_sn.linear_bounds (EXACTO_SN linear bounds, three combinations) and of the exacto_sn part of
maximos_globales.py (global-maximum bounds on the extended grid), with the SN rows of the response built with the definition of
respuesta_v3A.build_A_v3 applied on their own grid (cmb_fondo uniform grid, fiducial CF.OM_FID):
    row_SN = (5/ln 10) [(1 - t) dD~_M(u_k) + t dD~_M(u_{k+1})] / D~_M(u_HD),  D~_M(u_HD) = (1 - t) D~_M(u_k) + t D~_M(u_{k+1}),
then column 0 <- -sum (convention B). The old rows interpolated the per-node d ln D~_M (errata 129). Everything else is copied
unchanged from exacto_sn.linear_bounds (l. 272-363) and maximos_globales.main (l. 90-178). The old functions and json files are not
modified. Mode 'old' reproduces the originals (check); mode 'v3A' uses the corrected rows.

    OMP_NUM_THREADS=1  ../venv/bin/python -I exacto_sn_lineal_v3A.py exacto_sn old|v3A     (exacto_sn ran with 1 thread)
    OMP_NUM_THREADS=16 ../venv/bin/python -I exacto_sn_lineal_v3A.py maximos  old|v3A     (maximos_globales ran with the default)
Writes exacto_sn_lineal_v3A_<which>_<mode>.json.
"""
import json, sys, time
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cmb_fondo as CF  # noqa: E402
import exacto_sn as S  # noqa: E402
import maximos_globales as MGL  # noqa: E402  (GRID, SIG, THRS, bound_from, entry; its main() is not run)

which, mode = sys.argv[1], sys.argv[2]
tag = ("_" + sys.argv[3]) if len(sys.argv) > 3 else ""   # optional output tag (e.g. the BLAS thread count)
t0 = time.time()
R = CF.Response(lambda_absorbs_x0=True)
u, z = R.u, R.z
om, orad = CF.OM_FID, S.OMR / CF.H_FID**2
E2 = om * (1 + z) ** 3 + orad * (1 + z) ** 4 + (1 - om - orad)
Es = np.sqrt(E2)
base = R.du * (1 + z)
N = len(u)
sn = S.SN["all"]

def cw_at(k):
    cw = np.zeros(N); cw[: k + 1] = base[: k + 1]; cw[0] *= 0.5; cw[k] *= 0.5
    return cw

def wM_at(k):                      # old: per-node d ln D~_M
    cw = cw_at(k); dm = float(np.sum(cw / Es))
    return -cw / (2 * Es**3) / dm

def dDM_at(k):                     # v3A ingredients: d D~_M(u_k) (un-normalised) and D~_M(u_k)
    cw = cw_at(k)
    return -cw / (2 * Es**3), float(np.sum(cw / Es))

Wsn = np.zeros((sn["n"], N))
for i, uu in enumerate(sn["uhd"]):
    k = int(np.searchsorted(u, uu)) - 1
    t = (uu - u[k]) / (u[k + 1] - u[k])
    if mode == "old":
        row = (1 - t) * wM_at(k) + t * wM_at(k + 1)
    else:
        d0, m0_ = dDM_at(k); d1, m1_ = dDM_at(k + 1)
        row = ((1 - t) * d0 + t * d1) / ((1 - t) * m0_ + t * m1_)
    row[0] -= row.sum()
    Wsn[i] = S.MAG * row
W = np.vstack([R.W[:15], Wsn])
fit15 = CF.Fit(CF.SETS["BAO+CMB"], R)
scale = np.concatenate([fit15.scale[:15], np.ones(sn["n"])])
Ws = W * scale[:, None]
y = np.concatenate([S.X.Y15, sn["mb"]])
Cm = np.zeros((15 + sn["n"], 15 + sn["n"]))
Cm[:15, :15] = S.X.CMEAS15 if hasattr(S.X, "CMEAS15") else CF.CMEAS15
Cm[15:, 15:] = np.linalg.inv(sn["Ci"])
rows_of = {"SN": list(range(15, 15 + sn["n"])), "BAO+SN": list(range(13)) + list(range(15, 15 + sn["n"])), "BAO+CMB+SN": list(range(15 + sn["n"]))}
G0 = S.G0
dX0 = np.zeros(G0.N)

def model_full(h, wm, M):
    hist, _ = S.history(G0, wm / h**2, S.OMR / h**2, dX0)
    return np.concatenate([S.vec15(G0, h, wm, *hist), S.mu_sn(G0, sn, hist[0]) + M])

out = {"which": which, "mode": mode, "sn_row_definition": "old: interpolated per-node d ln D~_M" if mode == "old" else "v3A: d ln of the interpolated D~_M (build_A_v3 definition on the cmb_fondo grid)", "combos": {}}
for combo, rows in rows_of.items():
    rows = np.array(rows)
    if combo == "SN":
        c0, om_b = S.fit_sn_only(G0, dX0, sn)
        h_b, wm_b = S.H_SN, om_b * S.H_SN**2
        hist, _ = S.history(G0, om_b, S.OMR / S.H_SN**2, dX0)
        rr = sn["mb"] - S.mu_sn(G0, sn, hist[0])
        M_b = float((sn["Ci1"] @ rr) / sn["e"])
        pars = ["Omega_m", "M"]
    else:
        c0, h_b, wm_b, M_b = S.fit_combo(G0, dX0, sn, combo == "BAO+CMB+SN")
        pars = ["h", "omega_m", "M"]
    m0 = model_full(h_b, wm_b, M_b)
    T = []
    for pname in pars:
        eps = {"h": 1e-4, "omega_m": 1e-5, "Omega_m": 1e-4, "M": 1e-3}[pname]
        def mp(sign):
            if pname == "h": return model_full(h_b + sign * eps, wm_b, M_b)
            if pname == "omega_m": return model_full(h_b, wm_b + sign * eps, M_b)
            if pname == "Omega_m": return model_full(h_b, (om_b + sign * eps) * S.H_SN**2, M_b)
            return model_full(h_b, wm_b, M_b + sign * eps)
        T.append((mp(+1) - mp(-1)) / (2 * eps))
    T = np.array(T).T[rows]
    r0 = (y - m0)[rows]
    C = Cm[np.ix_(rows, rows)]
    Wr = Ws[rows]
    rec = {"lcdm_chi2": c0, "params_best": {"h": h_b, "omega_m": wm_b, "M": M_b}}
    for theta in S.THETAS:
        Kw = R.Kw(theta, Wr)
        Mth = Wr @ Kw.T
        def ll_c(s, Mth=Mth):
            L = np.linalg.cholesky(C + s * s * Mth)
            solve = lambda b: np.linalg.solve(L.T, np.linalg.solve(L, b))
            cr, cT = solve(r0), solve(T)
            F = T.T @ cT
            b = T.T @ cr
            return float(-0.5 * (r0 @ cr - b @ np.linalg.solve(F, b) + 2 * np.sum(np.log(np.diag(L))) + np.linalg.slogdet(F)[1]))
        if which == "exacto_sn":
            rec[f"theta_{theta:g}"] = {name: CF.V.bound(ll_c, thr)["sigma_95"] for name, thr in CF.THRS.items()}
        else:
            old = json.load(open(HERE / "exacto_sn.json"))["linear"][combo][f"theta_{theta:g}"]
            rec[f"theta_{theta:g}"] = MGL.entry(ll_c, old["192"], old["135"])
        print(which, mode, combo, theta, json.dumps(rec[f"theta_{theta:g}"], default=float)[:160], flush=True)
    out["combos"][combo] = rec
out["seconds"] = time.time() - t0
(HERE / f"exacto_sn_lineal_v3A_{which}_{mode}{tag}.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")
print("done", which, mode, f"{out['seconds']:.0f} s")
