#!/usr/bin/env python3
"""
P1 of the blind reproduction: which SN row of A is the exact Jacobian of the forward model used by the exact likelihood?
Forward model: marginal_guiada.exact_model (identical SN treatment to exacto_sn.mu_sn: D~_M interpolated LINEARLY in u = ln(1+z) between
the two bracketing grid nodes, then mu = 5 log10((1+z_hel) D~_M) + M). Numerical Jacobian of the SN rows (mu, magnitudes) with respect
to X_j at X = 0 and p0, convention B (X enters as X - X(0): node j >= 1 perturbed by +-h; node 0 = all other nodes perturbed by -+h),
central differences with h = 1e-3 and h/10 = 1e-4; convergence criterion (fixed before running): on non-negligible elements
(|A_ij| > 1e-6 max_i |A_ij| of the column) relative difference between h and h/10 < 1e-6, or absolute < 1e-12 if the element is
below 1e-6; the h/10 result is used. Compared with (i) our A (marginal_guiada.build_A, SN rows = MAG * linear interpolation in u of the
per-node d ln D~_M) and (ii) the reproducer's A (reproduccion_ciega/mi_codigo/salidas/A.npy, read only; its SN rows = MAG * d D~_M(u_HD)
/ D~_M(u_HD) with d D~_M interpolated, paso2_lineal_lnZ.py line 65), and (iii) that formula rebuilt from our own arrays as a check.
All three in the same units (d mu / d X_j, magnitudes per unit X), same 1580 rows in file order, same p0 (ours; the reproducer's p0
differs by ~1e-9). Writes jacobiano_sn_check.json. Nothing in the repository is modified.
"""
import json, os, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import marginal_guiada as MG  # noqa: E402

G, SN, P0 = MG.G, MG.SN, MG.P0
N = G.N; NSN = SN["n"]
A_ours = MG.build_A()
A_rep_path = Path(os.path.expanduser("~/Proyectos/01_Fisica_y_Cosmologia/reproduccion_ciega/mi_codigo/salidas/A.npy"))
A_rep = np.load(A_rep_path)  # read only
# (iii) reproducer's SN-row formula rebuilt from our arrays: d D~_M per node from build_A's ingredients
h, wm = P0[0], P0[1]
om, orad = wm / h**2, MG.OMR / h**2
E2 = om * G.p3 + orad * G.p4 + (1 - om - orad); Es = np.sqrt(E2)
wint = np.zeros(N); wint[1:] += 0.5 * G.DUK; wint[:-1] += 0.5 * G.DUK
g = G.p1 / (2 * Es**3); f = G.p1 / Es
dm = np.concatenate([[0.0], np.cumsum(0.5 * (f[:-1] + f[1:]) * G.DUK)])
T = np.tril(np.tile(wint, (N, 1)), -1); T[np.arange(1, N), np.arange(1, N)] = 0.5 * G.DUK; T[0, 0] = 0.0
dDM = -(T * g[None, :]); dDM[:, 0] -= dDM.sum(axis=1)  # d D~_M(u_k)/d X_j, convention B
k_, t_ = MG._k, MG._t
DM_sn = (1 - t_) * dm[k_] + t_ * dm[k_ + 1]
A_rep_rebuilt_sn = MG.MAG * ((1 - t_)[:, None] * dDM[k_] + t_[:, None] * dDM[k_ + 1]) / DM_sn[:, None]
# ---- node sets
rep = json.load(open(HERE / "reproduccion_ciega/mi_codigo/salidas/paso2.json"))
nodes_prev = [r["nodo"] for r in rep["fd_nodos_z_le_2.33"]]
order = np.argsort(SN["zhd"]); low, mid, high = order[:10], order[NSN // 2 - 5: NSN // 2 + 5], order[-10:]
brack = sorted(set(int(v) for s in np.concatenate([low, mid, high]) for v in (k_[s], k_[s] + 1)))
nodes = sorted(set(nodes_prev) | set(brack) | {0})
# ---- numerical Jacobian of the SN rows
def model_sn(dX):
    out, bad = MG.exact_model(P0[0], P0[1], P0[2], dX[None, :])
    assert not bad[0]
    return out[0, 15:]
def col(j, hstep):
    dX = np.zeros(N)
    if j == 0:
        dX[1:] = -hstep; yp = model_sn(dX); dX[1:] = +hstep; ym = model_sn(dX)
    else:
        dX[j] = +hstep; yp = model_sn(dX); dX[j] = -hstep; ym = model_sn(dX)
    return (yp - ym) / (2 * hstep)
res = []; conv_all = True; worst = 0.0
for j in nodes:
    c1, c2 = col(j, 1e-3), col(j, 1e-4)
    scale = np.max(np.abs(c2)); mask = np.abs(c2) > 1e-6 * scale if scale > 0 else np.zeros(NSN, bool)
    big = mask & (np.abs(c2) > 1e-6); small = mask & ~big
    rel = np.abs(c1[big] - c2[big]) / np.abs(c2[big]) if big.any() else np.zeros(0)
    ab = np.abs(c1[small] - c2[small]) if small.any() else np.zeros(0)
    conv = bool((rel < 1e-6).all() and (ab < 1e-12).all())
    conv_all &= conv; worst = max(worst, float(rel.max()) if rel.size else 0.0)
    def err(Acol):
        d = np.abs(Acol[mask] - c2[mask]) / np.abs(c2[mask])
        return {"max": float(d.max()), "median": float(np.median(d))} if mask.any() else {"max": None, "median": None}
    res.append({"node": int(j), "z": float(G.Z[j]), "n_nonneg": int(mask.sum()), "converged": conv, "max_rel_diff_h_vs_h10": float(rel.max()) if rel.size else 0.0,
                "max_abs_diff_small": float(ab.max()) if ab.size else 0.0, "max_abs_jac": float(scale),
                "ours": err(A_ours[15:, j]), "reproducer_Anpy": err(A_rep[15:, j]), "reproducer_formula_rebuilt": err(A_rep_rebuilt_sn[:, j])})
agg = {}
for key in ("ours", "reproducer_Anpy", "reproducer_formula_rebuilt"):
    mx = [r[key]["max"] for r in res if r[key]["max"] is not None]; md = [r[key]["median"] for r in res if r[key]["median"] is not None]
    agg[key] = {"max_over_nodes": float(max(mx)), "median_of_node_medians": float(np.median(md)), "median_over_all_elements_proxy": float(np.median(md))}
out = {"forward_model": "marginal_guiada.exact_model (= exacto_sn.mu_sn): D~_M linear in u between bracketing nodes; mu = 5 log10((1+z_hel) D~_M) + M",
       "units": "all Jacobians in d mu / d X_j (magnitudes per unit X); convention B; 1580 SN in file order; p0 = marginal_guiada.P0",
       "h_steps": [1e-3, 1e-4], "criterion": "non-negligible: |A_ij| > 1e-6 * column max; converged if rel diff h vs h/10 < 1e-6 (elements > 1e-6) or abs < 1e-12 (smaller); h/10 used",
       "n_nodes": len(nodes), "nodes": nodes, "node_sets": {"previous_check_z_le_2.33": nodes_prev, "bracketing_10_lowest_10_middle_10_highest_z_SN": brack, "node_0": [0]},
       "converged_all": conv_all, "worst_rel_diff_h_vs_h10": worst, "per_node": res, "aggregate": agg,
       "P5": {"ours_max_rel_err": agg["ours"]["max_over_nodes"], "ours_lt_1e-5": agg["ours"]["max_over_nodes"] < 1e-5, "ours_better_than_reproducer": agg["ours"]["max_over_nodes"] < agg["reproducer_Anpy"]["max_over_nodes"]},
       "check_reproducer_Anpy_vs_rebuilt_formula_max_abs": float(np.max(np.abs(A_rep[15:, nodes] - A_rep_rebuilt_sn[:, nodes]))),
       "check_reproducer_Anpy_vs_ours_rows_0_15_max_abs": float(np.max(np.abs(A_rep[:15] - A_ours[:15])))}
(HERE / "jacobiano_sn_check.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
print("nodes", len(nodes), "converged", conv_all, "worst h vs h/10", worst)
print(json.dumps(agg, indent=1)); print("P5", out["P5"]); print("rep A.npy vs rebuilt:", out["check_reproducer_Anpy_vs_rebuilt_formula_max_abs"], "| rep rows 0-15 vs ours:", out["check_reproducer_Anpy_vs_ours_rows_0_15_max_abs"])
for r in res[:6] + res[-3:]: print(r["node"], f"z={r['z']:.4f}", "conv", r["converged"], f"ours {r['ours']['max']:.2e}/{r['ours']['median']:.2e}", f"rep {r['reproducer_Anpy']['max']:.2e}/{r['reproducer_Anpy']['median']:.2e}")
