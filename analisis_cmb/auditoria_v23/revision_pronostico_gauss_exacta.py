#!/usr/bin/env python3
"""
REVISION, diagnostico: los fallos de Rayleigh/Kuiper en revision_pronostico_gaussianidad.py ¿son no gaussianidad o el desplazamiento
del centro (media de A estimada con error Monte Carlo ~0.0045 sigma, y la propia media de B)? Se repiten los tests de angulo en B
blanqueando con la media y covarianza EXACTAS del modelo linealizado (affine_map) y, ademas, en 20 mitades B nuevas con otras
semillas (blanqueadas con la A original y con la exacta). Escribe revision_pronostico_gauss_exacta.json. Nada existente se modifica.
"""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import revision_pronostico_reproduccion as R  # noqa: E402
import revision_pronostico_gaussianidad as G  # noqa: E402

out = {}
for ic, combo in enumerate(("sinSN", "Pantheon+")):
    st = R.setup(combo); op = R.case_ops(combo, st, 0.6)
    for ih, hyp in enumerate(("H0", "H1")):
        a, Cx, _ = R.affine_map(hyp, st, op, 0.6)
        Lx = np.linalg.inv(np.linalg.cholesky(Cx))
        rA = R.simulate(hyp, st, op, 0.6, R.N_UNI // 2, np.random.default_rng([R.SEED_BASE, ic, 0, ih, 0]))[0]
        rB = R.simulate(hyp, st, op, 0.6, R.N_UNI // 2, np.random.default_rng([R.SEED_BASE, ic, 0, ih, 1]))[0]
        LA = np.linalg.inv(np.linalg.cholesky(np.cov(rA.T)))
        tx = G.tests_2d((rB - a) @ Lx.T)
        shift = (rA.mean(0) - a) @ Lx.T
        rec = {"B_whitened_with_exact": {k: tx[k] for k in ("angle_rayleigh", "angle_kuiper", "radius_ks_vs_chi2dof", "fails")},
               "A_mean_minus_exact_whitened": shift.tolist(), "A_mean_offset_norm_over_sigma": float(np.hypot(*shift)),
               "rayleigh_p_new_B_with_A_center": [], "rayleigh_p_new_B_with_exact_center": []}
        for k in range(20):
            rBk = R.simulate(hyp, st, op, 0.6, R.N_UNI // 2, np.random.default_rng([R.SEED_BASE, ic, 0, ih, 50 + k, 4242]))[0]
            rec["rayleigh_p_new_B_with_A_center"].append(G.rayleigh(np.mod(np.arctan2(*((rBk - rA.mean(0)) @ LA.T)[:, ::-1].T), 2 * np.pi))[2])
            rec["rayleigh_p_new_B_with_exact_center"].append(G.rayleigh(np.mod(np.arctan2(*((rBk - a) @ Lx.T)[:, ::-1].T), 2 * np.pi))[2])
        out[f"{combo}/{hyp}"] = rec
        print(f"{combo} {hyp}: B blanqueada con la exacta: Rayleigh p={tx['angle_rayleigh']['p']:.3f} Kuiper p={tx['angle_kuiper']['p']:.3f} KS radio p={tx['radius_ks_vs_chi2dof']['p']:.3f} fails={tx['fails']}; "
              f"|media A - exacta|/sigma={rec['A_mean_offset_norm_over_sigma']:.4f}; 20 B nuevas: Rayleigh p (centro A) min={min(rec['rayleigh_p_new_B_with_A_center']):.3f} mediana={np.median(rec['rayleigh_p_new_B_with_A_center']):.3f} | (centro exacto) min={min(rec['rayleigh_p_new_B_with_exact_center']):.3f} mediana={np.median(rec['rayleigh_p_new_B_with_exact_center']):.3f}", flush=True)
(HERE / "revision_pronostico_gauss_exacta.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")
print("escrito revision_pronostico_gauss_exacta.json")
