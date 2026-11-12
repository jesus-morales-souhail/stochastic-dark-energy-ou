#!/usr/bin/env python3
"""
Side-by-side comparison of the blind reproduction (reproduccion_ciega/mi_codigo/salidas/paso2.json, paso1.json; copies under
auditoria_v23/reproduccion_ciega/) with our linear marginal for OU theta = 1 (positividad_lin.json and marginal_guiada.json:
thetas.theta_1.lnZ_lin_rel on the same 10 sigma_X nodes, same normalisation ln Z(sigma) - ln Z(0), flat priors on (h, omega_m, M),
same constants). Plus one diagnostic, not a new result: our ln Z recomputed with the Pantheon+ covariance taken directly from the
file (exacto_sn.load_sn stores Ci = inv(C) and marginal_guiada rebuilds C = inv(Ci); the reproducer never inverts twice).
Writes comparacion_reproduccion.json. Nothing in the repository results is modified.
"""
import json, sys
from pathlib import Path
import numpy as np
from scipy.linalg import cho_factor
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import marginal_guiada as MG  # noqa: E402
import precision_sn as PS  # noqa: E402
import exacto_sn as S  # noqa: E402

rep1 = json.load(open(HERE / "reproduccion_ciega/mi_codigo/salidas/paso1.json")); rep2 = json.load(open(HERE / "reproduccion_ciega/mi_codigo/salidas/paso2.json"))
pl = json.load(open(HERE / "positividad_lin.json"))["thetas"]["theta_1"]; mg = json.load(open(HERE / "marginal_guiada.json"))
sig_r = [c["sigma_X"] for c in rep2["lnZ_rel"]]; v_r = [c["lnZ_menos_lnZ0"] for c in rep2["lnZ_rel"]]
assert sig_r == pl["sigmas"] == mg["thetas"]["theta_1"]["sigmas"]
v_pl = pl["lnZ_lin_rel"]; v_mg = mg["thetas"]["theta_1"]["lnZ_lin_rel"]
nodes = [{"sigma_X": s, "ours_positividad_lin": a, "ours_marginal_guiada": m, "reproducer": r, "diff_reproducer_minus_ours": r - a, "abs_diff_lt_1e-4": abs(r - a) < 1e-4} for s, a, m, r in zip(sig_r, v_pl, v_mg, v_r)]
b_ours = pl["bound_lin_192"]; b_rep = rep2["cota"]
out = {"what_is_compared": "ln Z(sigma_X) - ln Z(0), analytic linear marginal, OU theta = 1, convention B, BAO + theta_* + omega_m + Pantheon+, flat priors on (h, omega_m, M); same 10 nodes in the same order; same constants (constantes.json exported from our code)",
       "reproducer_json_precision": {"lnZ_rel": "full double precision (paso2.json lnZ_rel[*].lnZ_menos_lnZ0)", "cota": "full double precision (paso2.json cota)", "p0": "full double precision (paso1.json p0)"},
       "nodes": nodes, "P1_all_nodes_within_1e-4": all(n["abs_diff_lt_1e-4"] for n in nodes), "P1_n_nodes_within_1e-4": int(sum(n["abs_diff_lt_1e-4"] for n in nodes)),
       "P1_max_abs_diff": float(max(abs(n["diff_reproducer_minus_ours"]) for n in nodes)),
       "bound": {"ours_positividad_lin": b_ours, "ours_marginal_guiada": mg["thetas"]["theta_1"]["lin_192"], "reproducer": b_rep, "rel_diff": (b_rep - b_ours) / b_ours, "P2_rel_diff_lt_1e-4": abs(b_rep - b_ours) / b_ours < 1e-4},
       "fiducial": {"ours": mg["checks"]["p0"], "reproducer": {**rep1["p0"], "chi2": rep1["chi2_p0"]}},
       "max_lnZ": {"ours_sigma": pl["max_lin_at_sigma"], "reproducer_sigma": rep2["maximo"]["sigma_X"]}}
# ---- diagnostic: covariance taken directly from the file (single read, no double inversion)
dat, cov = S.SN_DIR / "Pantheon+SH0ES.dat", S.SN_DIR / "Pantheon+SH0ES_STAT+SYS.cov"
with open(dat) as f: head = f.readline().split()
raw = np.genfromtxt(dat, skip_header=1, dtype=None, encoding="ascii", names=head)
n = int(open(cov).readline()); Cfull = np.loadtxt(cov, skiprows=1).reshape(n, n)
idx = np.nonzero((raw["zHD"] > 0.01) & (raw["IS_CALIBRATOR"] == 0))[0]
Cs = Cfull[np.ix_(idx, idx)]
C2 = MG.C.copy(); C2[15:, 15:] = Cs
Bm, A = MG.build_B(), MG.build_A()
K = MG.kernel(1.0); Msig = (A @ K) @ A.T; Msig = 0.5 * (Msig + Msig.T)
r = MG.D - MG.B0
def prof(C):
    return [MG.lnZ_lin(cho_factor(C + s * s * Msig, lower=True), Bm, r) for s in sig_r]
z_std = np.array(prof(MG.C)); z_dir = np.array(prof(C2))
z_std -= z_std[0]; z_dir -= z_dir[0]
sig = np.array(sig_r); m = sig > 0
out["diagnostic_direct_covariance"] = {"note": "our ln Z recomputed in this script with C_SN read once from the file (Cs) instead of inv(Ci); the standard pipeline value recomputed here is also given to show the script reproduces positividad_lin.json",
    "max_abs_diff_C_direct_vs_inv_inv": float(np.max(np.abs(C2[15:, 15:] - MG.C[15:, 15:]))),
    "lnZ_rel_standard_recomputed": z_std.tolist(), "lnZ_rel_direct_covariance": z_dir.tolist(),
    "max_abs_diff_recomputed_vs_positividad_lin": float(np.max(np.abs(z_std - np.array(v_pl)))),
    "max_abs_diff_direct_vs_reproducer": float(np.max(np.abs(z_dir - np.array(v_r)))), "max_abs_diff_standard_vs_reproducer": float(np.max(np.abs(z_std - np.array(v_r)))),
    "bound_standard_recomputed": PS.lnl_bound(sig[m], z_std[m], 1.92), "bound_direct_covariance": PS.lnl_bound(sig[m], z_dir[m], 1.92)}
(HERE / "comparacion_reproduccion.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
for nd in nodes: print(f"{nd['sigma_X']:<5g} ours {nd['ours_positividad_lin']:+.9f} rep {nd['reproducer']:+.9f} diff {nd['diff_reproducer_minus_ours']:+.2e} {'ok' if nd['abs_diff_lt_1e-4'] else 'NO'}")
print("P1", out["P1_all_nodes_within_1e-4"], out["P1_n_nodes_within_1e-4"], "max", out["P1_max_abs_diff"]); print("bound", out["bound"])
print("diag", json.dumps({k: v for k, v in out["diagnostic_direct_covariance"].items() if "lnZ_rel" not in k}, indent=1))
