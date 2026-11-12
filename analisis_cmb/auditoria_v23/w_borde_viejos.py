#!/usr/bin/env python3
"""
w at the table edge a = 1e-3 in the old CMB runs (realizaciones.py / malla, picos.py, rapido.py), regenerated with their seeds
and their grids. CAMB extrapolates w = w(1e-3) below the table; with w < 1/3 the dark-energy density grows backwards slower
than radiation (rho_DE ~ a^-3(1+w) vs a^-4) and cannot dominate before recombination: no edge artefact.
Writes auditoria_v23/w_borde_viejos.json.
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import ou_cmb as M  # noqa: E402

lna_f, tab_idx, main = M.production_grids()
h = lna_f[1] - lna_f[0]
lna_tab = lna_f[tab_idx]
om_l = M.omega_l0()
Z200 = np.random.default_rng(M.SEED).standard_normal((200, len(lna_f)))
RES = HERE.parent / "resultados" / "takahashi"
out = {"a_edge": float(np.exp(lna_tab[main][0])), "cases": []}
unit = {}


def case(label, theta, sigma, delta, n):
    if theta not in unit:
        unit[theta] = M.ou_uniform(Z200, h, theta, 1.0)
    xs = M.smooth_uniform(sigma * unit[theta][:n], h, delta)
    rho, w, bad = M.rho_w_table(lna_tab, main, xs[:, tab_idx], om_l)
    w_edge = w[:, 0]
    rec = {"label": label, "theta": theta, "sigma_X": sigma, "Delta": delta, "n": n, "n_rho_le_0": int(bad.sum()),
           "w_edge_max": float(np.nanmax(w_edge)), "w_edge_min": float(np.nanmin(w_edge)),
           "abs_wp1_max_table": float(np.nanmax(np.abs(w + 1))), "w_edge_lt_one_third_all": bool(np.nanmax(w_edge) < 1 / 3)}
    out["cases"].append(rec)
    print(f"{label:48s} theta={theta:<6g} sigma={sigma:.3e} Delta={delta:<6g} n={n:3d}  w_edge max={rec['w_edge_max']:+.4f} min={rec['w_edge_min']:+.4f}  |w+1|max tabla={rec['abs_wp1_max_table']:.3f}", flush=True)


# realizaciones.py section 3.5 (and the 2.506e-2 node) and the malla: cases from the result folders
for folder in sorted(list(RES.glob("sigma_*")) + list((RES / "malla").glob("sigma_*")), key=lambda p: (p.parent.name, float(p.name[6:]))):
    meta = json.loads((folder / "meta.json").read_text())
    for f in sorted(folder.glob("observables_theta*_Delta*.npz")):
        tag = f.stem[len("observables_"):]
        theta = float(tag.split("_")[0][5:])
        delta = float(tag.split("_")[1][5:])
        case(f"{folder.parent.name}/{folder.name}/{tag}", theta, meta["sigma_X"], delta, meta["n_real"])
# picos.py: first 60, Delta 0.05
pf = json.loads((HERE.parent / "fondo" / "perfil_fondo.json").read_text())["treatments"]["c_marginalised"]["theta_1"]["sigma_95"]
for lab, th, sg in (("picos/theta1_C2", 1.0, 0.0263282756994459), ("picos/theta0.001_C2", 0.001, 0.0263282756994459), ("picos/theta1_marg", 1.0, pf)):
    case(lab, th, sg, 0.05, 60)
# rapido.py: CASES x sigmas, first 60
pf10 = json.loads((HERE.parent / "fondo" / "perfil_fondo.json").read_text())["treatments"]
for th, dl in ((3.0, 0.05), (10.0, 0.0125), (10.0, 0.05), (10.0, 0.025), (3.0, 0.025)):
    for sg in (pf10["a_fixed"]["theta_10"]["sigma_95"], pf10["c_marginalised"]["theta_10"]["sigma_95"]):
        case(f"rapido/theta{th:g}_Delta{dl:g}", th, sg, dl, 60)
out["w_edge_max_over_all_cases"] = max(c["w_edge_max"] for c in out["cases"])
out["all_below_one_third"] = all(c["w_edge_lt_one_third_all"] for c in out["cases"])
(HERE / "w_borde_viejos.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
print("max w_edge over all cases:", out["w_edge_max_over_all_cases"], "all < 1/3:", out["all_below_one_third"])
