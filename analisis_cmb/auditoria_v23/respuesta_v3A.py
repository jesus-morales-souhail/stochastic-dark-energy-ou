#!/usr/bin/env python3
"""
Corrected linear response (v3A): the SN rows of A are the exact Jacobian of the forward model used by the exact likelihood
(exacto_sn.mu_sn: D~_M interpolated linearly in u between the two bracketing nodes). Definition fixed in REPRODUCCION.md,
'Puntos a fijar en la especificacion v3':
    dy_SN/dX_j = (5/ln 10) * [(1 - t) dD~_M(u_k)/dX_j + t dD~_M(u_{k+1})/dX_j] / D~_M(u_HD),  D~_M(u_HD) = (1 - t) D~_M(u_k) + t D~_M(u_{k+1}).
The old build_A (marginal_guiada.py:147, errata 129) interpolated the per-node d ln D~_M instead. Everything else (BAO rows, theta_*
row, omega_m row, convention B: rows sum to zero) is unchanged and is computed with the same code paths. The old functions are not
modified: this module provides build_A_v3() for the marginal_guiada set-up (Pantheon+, grid G) and make_linear_v3(sn) for the
per-sample set-up of muestras_sn.make_linear.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
from scipy.linalg import cho_factor
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import marginal_guiada as MG  # noqa: E402
import muestras_sn as MS  # noqa: E402

_MAKE_LINEAR_ORIG = MS.make_linear  # captured at import, so that monkeypatching MS.make_linear with make_linear_v3 does not recurse


def _ingredients(G, h, wm, OMR):
    om, orad = wm / h**2, OMR / h**2
    E2 = om * G.p3 + orad * G.p4 + (1 - om - orad)
    Es = np.sqrt(E2)
    N = G.N
    wint = np.zeros(N)
    wint[1:] += 0.5 * G.DUK
    wint[:-1] += 0.5 * G.DUK
    g = G.p1 / (2 * Es**3)
    f = G.p1 / Es
    dm = np.concatenate([[0.0], np.cumsum(0.5 * (f[:-1] + f[1:]) * G.DUK)])
    T = np.tril(np.tile(wint, (N, 1)), -1)
    T[np.arange(1, N), np.arange(1, N)] = 0.5 * G.DUK
    T[0, 0] = 0.0
    dDM = -(T * g[None, :])            # d D~_M(u_k) / d X_j before convention B
    dDM[:, 0] -= dDM.sum(axis=1)       # convention B: Lambda absorbs X(0) -> rows sum to zero
    Wdm = dDM / np.where(dm > 0, dm, 1.0)[:, None]   # d ln D~_M(u_k) / d X_j (per node), as build_A
    WdH = np.zeros((len(G.IDX), N))
    for i, k in enumerate(G.IDX):
        WdH[i, k] = -1.0 / (2 * E2[k])
        WdH[i, 0] += 1.0 / (2 * E2[k])
    return E2, dm, dDM, Wdm, WdH


def _sn_rows_v3(dDM, dm, k, t, MAG):
    DM_sn = (1 - t) * dm[k] + t * dm[k + 1]
    return MAG * ((1 - t)[:, None] * dDM[k] + t[:, None] * dDM[k + 1]) / DM_sn[:, None]


def build_A_v3():
    """As marginal_guiada.build_A but with the exact-Jacobian SN rows."""
    G, X, P0, B0 = MG.G, MG.X, MG.P0, MG.B0
    E2, dm, dDM, Wdm, WdH = _ingredients(G, P0[0], P0[1], MG.OMR)
    rows = []
    for i, (k, kind) in enumerate(zip(G.IDX, X.KINDS13)):
        rows.append(Wdm[k] if kind == "DM" else (WdH[i] if kind == "DH" else (2 * Wdm[k] + WdH[i]) / 3))
    A = np.zeros((MG.NROWS, G.N))
    A[:13] = np.array(rows) * B0[:13, None]
    A[13] = -Wdm[-1] * B0[13]
    A[15:] = _sn_rows_v3(dDM, dm, MG._k, MG._t, MG.MAG)
    return A


def make_linear_v3(sn: dict) -> dict:
    """As muestras_sn.make_linear (same P0, B0, Bm, C) but with the exact-Jacobian SN rows in A."""
    L = _MAKE_LINEAR_ORIG(sn)   # unchanged original function; its A is replaced below
    G0, X, MAG = MS.G0, MS.X, MS.MAG
    h, wm = L["P0"][0], L["P0"][1]
    E2, dm, dDM, Wdm, WdH = _ingredients(G0, h, wm, MS.OMR)
    rows = []
    for i, (k, kind) in enumerate(zip(G0.IDX, X.KINDS13)):
        rows.append(Wdm[k] if kind == "DM" else (WdH[i] if kind == "DH" else (2 * Wdm[k] + WdH[i]) / 3))
    A = np.zeros_like(L["A"])
    A[:13] = np.array(rows) * L["B0"][:13, None]
    A[13] = -Wdm[-1] * L["B0"][13]
    A[15:] = _sn_rows_v3(dDM, dm, sn["_k"], sn["_t"], MAG)
    L["A"] = A
    L["A_rows_sum_max_abs"] = float(np.max(np.abs(A.sum(axis=1))))
    L["A_version"] = "v3A (exact-Jacobian SN rows)"
    return L


if __name__ == "__main__":
    import json, os
    A3 = build_A_v3(); A_old = MG.build_A()
    Arep = np.load(os.path.expanduser("~/Proyectos/01_Fisica_y_Cosmologia/reproduccion_ciega/mi_codigo/salidas/A.npy"))
    jc = json.load(open(HERE / "jacobiano_sn_check.json"))
    conv = [r for r in jc["per_node"] if r["converged"]]
    # jacobian check: recompute the numerical columns (h/10) for the converged nodes and compare with A3
    import jacobiano_sn_check as J  # noqa: E402  (recomputes; read-only)
    errs = []
    for r in conv:
        j = r["node"]; c2 = J.col(j, 1e-4); scale = np.max(np.abs(c2)); mask = np.abs(c2) > 1e-6 * scale
        errs.append(float(np.max(np.abs(A3[15:, j][mask] - c2[mask]) / np.abs(c2[mask]))))
    chk = {"max_abs_diff_vs_reproducer_Anpy_all_rows": float(np.max(np.abs(A3 - Arep))), "max_abs_diff_vs_reproducer_Anpy_SN_rows": float(np.max(np.abs(A3[15:] - Arep[15:]))),
           "max_rel_err_vs_numerical_jacobian_converged_nodes": float(max(errs)), "n_converged_nodes": len(conv), "rows_sum_max_abs": float(np.max(np.abs(A3.sum(axis=1)))),
           "max_abs_diff_rows_0_15_vs_old_build_A": float(np.max(np.abs(A3[:15] - A_old[:15]))), "max_abs_diff_SN_rows_vs_old_build_A": float(np.max(np.abs(A3[15:] - A_old[15:]))),
           "passes": bool(float(np.max(np.abs(A3 - Arep))) < 1e-8 and max(errs) < 2e-6)}
    (HERE / "respuesta_v3A_check.json").write_text(json.dumps(chk, indent=1), encoding="utf-8")
    print(json.dumps(chk, indent=1))
