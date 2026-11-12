#!/usr/bin/env python3
"""FAMILIAS.md, comprobación del suelo de P en los escalones: masa del prior sin cambio sobre los datos frente a P(chi2_min <= chi2_obs) con sigma_X grande. Sin ajustes: los chi2 vienen de familias_chi2.npz y los desfases se regeneran de las semillas. Escribe familias_suelo_escalones.json."""
import sys, json, numpy as np
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import familias as F
raw = np.load(HERE / "familias_chi2.npz"); d = json.load(open(HERE / "familias.json")); obs = d["chi2_obs"]
u_bao = np.log1p(2.33); u_sn = np.log1p(F.SN["zhd"].max()); u_data = max(u_bao, u_sn); umax = F.X._u[-1]
print(f"u_max BAO = {u_bao:.3f}, u_max SN = {u_sn:.3f} (z = {F.SN['zhd'].max():.3f}), u(z*) = {umax:.3f}")
out = {}
for l in F.ELLS:
    rec = {}
    for s in (0.3, 0.5, 0.7, 1.0):
        rng = np.random.default_rng(F.cell_seed("Steps", s, l))
        phi, V = F.steps_draw(l, F.N_REAL, rng, umax)
        u1 = l - phi                                  # first step
        arr = raw[f"Steps__l_{l:g}__sigma_{s:g}"]; ok = arr[:, 0] > 0; chi = arr[:, 2]
        # consistency: drho at z=2.33 == 0 iff u1 > u_bao
        dz233 = arr[:, 3 + len(F.Z_EVAL) - 1]
        assert np.array_equal(dz233 == 0, u1 > u_bao - 1e-12), "inconsistency with stored drho"
        surv = ok & (chi <= obs)
        f_data = float(np.mean(u1 > u_data)); f_all = float(np.mean(u1 > umax)); f_z6 = float(np.mean(u1 > 2.0))
        P = float(np.mean(chi[ok] <= obs))
        cond = lambda m: float(np.mean(chi[ok & m] <= obs)) if (ok & m).sum() else None
        rec[f"sigma_{s:g}"] = {"P": P, "f_first_segment_covers_BAO_SN": f_data, "f_no_step_in_range": f_all, "f_first_step_beyond_z6.4": f_z6,
                               "P_given_covers_BAO_SN": cond(u1 > u_data), "P_given_no_step": cond(u1 > umax), "P_given_step_in_(1.2,2)": cond((u1 > u_data) & (u1 <= 2.0)),
                               "P_given_step_beyond_2": cond(u1 > 2.0), "P_given_step_below_1.2": cond(u1 <= u_data),
                               "frac_survivors_with_first_segment_covering_BAO_SN": float(np.mean(u1[surv] > u_data)) if surv.any() else None}
        print(f"l={l:<4g} s={s:<4g} P={P:.3f} | cubre BAO+SN {f_data:.3f} (P|cubre {rec[f'sigma_{s:g}']['P_given_covers_BAO_SN']}) | sin escalon {f_all:.3f} (P| {rec[f'sigma_{s:g}']['P_given_no_step']}) | salto z>6.4 {f_z6:.3f} (P| {rec[f'sigma_{s:g}']['P_given_step_beyond_2']}) | P|salto en (1.2,2) {rec[f'sigma_{s:g}']['P_given_step_in_(1.2,2)']} | P|salto z<2.3 {rec[f'sigma_{s:g}']['P_given_step_below_1.2']} | surv que cubren {rec[f'sigma_{s:g}']['frac_survivors_with_first_segment_covering_BAO_SN']}")
    out[f"l_{l:g}"] = rec
json.dump({"u_max_BAO": u_bao, "u_max_SN": u_sn, "u_zstar": umax, "by_l": out}, open(HERE / "familias_suelo_escalones.json", "w"), indent=1)

# --- desglose por z del primer escalón (l = 1 y 3)
J = json.load(open(HERE / "familias_suelo_escalones.json"))
bins = [(0.0, 0.5), (0.5, 1.0), (1.0, 1.5), (1.5, 2.33)]
for l in (1.0, 3.0):
    for s in (0.5, 1.0):
        rng = np.random.default_rng(F.cell_seed("Steps", s, l)); phi, V = F.steps_draw(l, F.N_REAL, rng, umax); z1 = np.expm1(l - phi)
        arr = raw[f"Steps__l_{l:g}__sigma_{s:g}"]; ok = arr[:, 0] > 0; chi = arr[:, 2]; surv = ok & (chi <= obs)
        rec = {}
        for a, b in bins:
            m = (z1 > a) & (z1 <= b)
            rec[f"z1_in_({a},{b}]"] = {"frac": float(m.mean()), "P_given": float(np.mean(chi[ok & m] <= obs)) if (ok & m).any() else None, "frac_of_survivors": float(np.mean(m[surv])) if surv.any() else None}
        J["by_l"][f"l_{l:g}"][f"sigma_{s:g}"]["by_first_step_z"] = rec
        print(l, s, {k: (round(v["frac"], 3), None if v["P_given"] is None else round(v["P_given"], 3), None if v["frac_of_survivors"] is None else round(v["frac_of_survivors"], 3)) for k, v in rec.items()})
json.dump(J, open(HERE / "familias_suelo_escalones.json", "w"), indent=1)
