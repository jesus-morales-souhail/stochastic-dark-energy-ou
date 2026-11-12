#!/usr/bin/env python3
"""
Sensitivity of ln Z(sigma_X) - ln Z(0) and of the -1.92 bound (OU theta = 1, convention B, BAO + theta_* + omega_m + Pantheon+) to the
SN rows of A, everything else identical (same data, p0, B, prior, analytic marginal, grid, rule; Msig = A K A' without symmetrisation,
as positividad.py; run with OMP_NUM_THREADS=16, the BLAS thread count with which positividad_lin.json was produced, so that variant (a) reproduces the reference bit for bit): (a) our A (build_A); (b) our A with the SN rows replaced by the reproducer's formula rebuilt from our arrays
(d ln of the interpolated D~_M); (c) the reproducer's full A.npy (read only). Writes sensibilidad_fila_sn.json.
"""
import json, os, sys
from pathlib import Path
import numpy as np
from scipy.linalg import cho_factor
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import marginal_guiada as MG  # noqa: E402
import precision_sn as PS  # noqa: E402
import jacobiano_sn_check as J  # noqa: E402  (runs the check again; provides A_ours, A_rep, A_rep_rebuilt_sn)

SIG = [0.0, 0.01, 0.03, 0.05, 0.07, 0.1, 0.15, 0.2, 0.3, 0.4]
REF = 0.27194769579562567
Bm = MG.build_B(); K = MG.kernel(1.0); r = MG.D - MG.B0
A_b = J.A_ours.copy(); A_b[15:] = J.A_rep_rebuilt_sn
def profile(A):
    AK = A @ K; Msig = AK @ A.T
    z = np.array([MG.lnZ_lin(cho_factor(MG.C + s * s * Msig, lower=True), Bm, r) for s in SIG]); z -= z[0]
    sig = np.array(SIG); m = sig > 0
    return z, PS.lnl_bound(sig[m], z[m], 1.92)
out = {"sigmas": SIG, "reference_bound": REF, "variants": {}}
for name, A in (("a_ours_build_A", J.A_ours), ("b_ours_with_reproducer_SN_rows", A_b), ("c_reproducer_A_npy", J.A_rep)):
    z, b = profile(A)
    out["variants"][name] = {"lnZ_rel": z.tolist(), "bound_192": b, "rel_diff_vs_reference": (b - REF) / REF, "max_lnZ_at_sigma": float(SIG[int(np.argmax(z))])}
    print(name, "bound", repr(b), "rel", (b - REF) / REF)
out["a_equals_reference_exactly"] = bool(out["variants"]["a_ours_build_A"]["bound_192"] == REF)
rep = json.load(open(HERE / "reproduccion_ciega/mi_codigo/salidas/paso2.json"))
vr = np.array([c["lnZ_menos_lnZ0"] for c in rep["lnZ_rel"]])
for name in out["variants"]:
    out["variants"][name]["max_abs_diff_lnZ_vs_reproducer"] = float(np.max(np.abs(np.array(out["variants"][name]["lnZ_rel"]) - vr)))
out["P6"] = {"winner_by_criterion": "the A closest to the numerical Jacobian (jacobiano_sn_check.json)", "rel_change_b": out["variants"]["b_ours_with_reproducer_SN_rows"]["rel_diff_vs_reference"], "rel_change_c": out["variants"]["c_reproducer_A_npy"]["rel_diff_vs_reference"]}
(HERE / "sensibilidad_fila_sn.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
print("a == REF exactly:", out["a_equals_reference_exactly"]); print({k: v["max_abs_diff_lnZ_vs_reproducer"] for k, v in out["variants"].items()})
