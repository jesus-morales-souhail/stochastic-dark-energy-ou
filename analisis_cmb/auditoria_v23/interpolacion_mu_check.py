#!/usr/bin/env python3
"""
P10: in the exact forward model at p0 and X = 0, mu of the 1580 Pantheon+ SN by linear interpolation of D~_M in u on the grid
(exacto_sn.py:147, np.interp) versus direct integration of (1+z)/E(u) du from 0 to u_HD of each SN (scipy quad, LCDM E at p0).
Writes interpolacion_mu_check.json.
"""
import json, sys
from pathlib import Path
import numpy as np
from scipy.integrate import quad
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import marginal_guiada as MG  # noqa: E402
import exacto_sn as S  # noqa: E402
G, SN, P0 = MG.G, MG.SN, MG.P0
h, wm = P0[0], P0[1]; om, orad = wm / h**2, MG.OMR / h**2
E2 = om * G.p3 + orad * G.p4 + (1 - om - orad)
f = G.p1 / np.sqrt(E2)
dm_u = np.concatenate([[0.0], np.cumsum(0.5 * (f[:-1] + f[1:]) * G.DUK)])
mu_interp = S.mu_sn(G, SN, dm_u)  # + M, dropped (same in both)
def integrand(u):
    z = np.expm1(u)
    return (1 + z) / np.sqrt(om * (1 + z) ** 3 + orad * (1 + z) ** 4 + (1 - om - orad))
dm_direct = np.array([quad(integrand, 0.0, u, epsabs=0, epsrel=1e-12, limit=200)[0] for u in SN["uhd"]])
mu_direct = 5.0 * np.log10((1.0 + SN["zhel"]) * dm_direct)
d = mu_interp - mu_direct
i = int(np.argmax(np.abs(d)))
# also: trapezoid on the grid evaluated exactly at the node (no interpolation error): the node values themselves vs quad
dm_nodes_q = np.array([quad(integrand, 0.0, u, epsabs=0, epsrel=1e-12, limit=200)[0] for u in G.U[[G.IDX[0], G.IDX[3], G.IDX[6]]]])
out = {"n_sn": int(SN["n"]), "p0": {"h": float(h), "omega_m": float(wm)}, "max_abs_dmu_mag": float(np.abs(d).max()), "median_abs_dmu_mag": float(np.median(np.abs(d))),
       "p95_abs_dmu_mag": float(np.percentile(np.abs(d), 95)), "worst_case": {"zHD": float(SN["zhd"][i]), "dmu_mag": float(d[i]), "u_offset_from_node_in_du": float((SN["uhd"][i] - G.U[MG._k[i]]) / (G.U[MG._k[i] + 1] - G.U[MG._k[i]]))},
       "trapezoid_grid_vs_quad_at_three_BAO_nodes_rel": [float(dm_u[k] / q - 1) for k, q in zip([G.IDX[0], G.IDX[3], G.IDX[6]], dm_nodes_q)],
       "P10_max_lt_1e-3_mag": bool(np.abs(d).max() < 1e-3), "definition": "mu_interp - mu_direct, M dropped; direct = quad of (1+z)/E(u) in u, rel tol 1e-12"}
(HERE / "interpolacion_mu_check.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
print(json.dumps(out, indent=1))
