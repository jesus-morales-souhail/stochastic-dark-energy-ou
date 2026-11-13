#!/usr/bin/env python3
"""P11 table: exacto_sn.linear_bounds and maximos_globales (exacto_sn part) with the old SN rows (reproduction) and the corrected ones
(_v3A), the A35 ratios (SN alone / BAO alone, BAO alone from maximos_globales cmb_fondo, unaffected) and the A26 figure (inyeccion
real_check.marg_linear old / v3A). Writes CORRECCION_FILA_SN_cierre_v3A.json."""
import json
from pathlib import Path
HERE = Path(__file__).resolve().parent
J = lambda p: json.load(open(HERE / p))
e_old, e_new = J("exacto_sn.json")["linear"], J("exacto_sn_lineal_v3A_exacto_sn_v3A.json")["combos"]
e_rep = J("exacto_sn_lineal_v3A_exacto_sn_old.json")["combos"]
m_old, m_new, m_rep = J("maximos_globales.json")["exacto_sn"], J("exacto_sn_lineal_v3A_maximos_v3A.json")["combos"], J("exacto_sn_lineal_v3A_maximos_old_omp16.json")["combos"]
m_rep1 = J("exacto_sn_lineal_v3A_maximos_old_omp1.json")["combos"]
bao = J("maximos_globales.json")["cmb_fondo"]["BAO|Lambda_absorbs_X0"]["marginalised"]
rel = lambda a, b: None if (a is None or b is None or a in (-1.0,) or a == 0) else (b - a) / a
rows = []; worst = 0.0; lnz_worst = 0.0; node_changes = []
for c in e_old:
    for th in ("theta_0.1", "theta_0.3", "theta_1", "theta_3", "theta_10"):
        for n in ("192", "135"):
            a, b = e_old[c][th][n], e_new[c][th][n]; r = rel(a, b)
            rows.append({"source": "exacto_sn.linear_bounds", "cell": f"{c} {th} {n}", "old": a, "new": b, "rel": r, "repro_old_abs": abs(e_rep[c][th][n] - a) if a is not None and e_rep[c][th][n] is not None else 0.0})
            if r is not None: worst = max(worst, abs(r))
for c in m_old:
    for th in m_old[c]:
        a, b = m_old[c][th], m_new[c][th]
        for n in ("192_global", "135_global", "192_clipped"):
            r = rel(a[n], b[n]); rows.append({"source": "maximos_globales exacto_sn", "cell": f"{c} {th} {n}", "old": a[n], "new": b[n], "rel": r,
                                               "repro_old_abs_16h": abs(m_rep[c][th][n] - a[n]) if a[n] and m_rep[c][th][n] else 0.0, "repro_old_abs_1h": abs(m_rep1[c][th][n] - a[n]) if a[n] and m_rep1[c][th][n] else 0.0})
            if r is not None: worst = max(worst, abs(r))
        rows.append({"source": "maximos_globales exacto_sn", "cell": f"{c} {th} máximo global: nodo σ / ln Z − ln Z(0)", "old": [a["sigma_at_max_grid"], a["lnL_max_minus_sigma0"]], "new": [b["sigma_at_max_grid"], b["lnL_max_minus_sigma0"]],
                     "same_node": a["sigma_at_max_grid"] == b["sigma_at_max_grid"], "dlnZ": b["lnL_max_minus_sigma0"] - a["lnL_max_minus_sigma0"]})
        node_changes.append(a["sigma_at_max_grid"] == b["sigma_at_max_grid"]); lnz_worst = max(lnz_worst, abs(b["lnL_max_minus_sigma0"] - a["lnL_max_minus_sigma0"]))
# A35 ratios
a35 = {}
for th, thk in (("theta_1", "theta_1.0"), ("theta_3", "theta_3.0"), ("theta_10", "theta_10.0")):
    sn_o, sn_n = m_old["SN"][th]["192_global"], m_new["SN"][th]["192_global"]; b_ = bao[thk]["192_global"]
    a35[th] = {"SN_old": sn_o, "SN_new": sn_n, "BAO_only": b_, "ratio_old": b_ / sn_o, "ratio_new": b_ / sn_n, "SN_quoted_old": round(sn_o, 3), "SN_quoted_new": round(sn_n, 3), "ratio_quoted_old": round(b_ / sn_o, 2), "ratio_quoted_new": round(b_ / sn_n, 2)}
# A26: inyeccion real_check.marg_linear
i_old, i_new = J("inyeccion.json")["real_check"], J("inyeccion_v3A.json")["real_check"]
a26 = {"marg_linear_old": i_old["marg_linear"], "marg_linear_new": i_new["marg_linear"], "quoted_old": round(i_old["marg_linear"], 4), "quoted_new": round(i_new["marg_linear"], 4), "ref_exact": i_old["ref"]["marg_exact"], "pct_old": round(100 * (i_old["marg_linear"] / 0.287 - 1)), "pct_new": round(100 * (i_new["marg_linear"] / 0.287 - 1))}
quoted_changes = [r["cell"] for r in rows if r.get("rel") is not None and round(r["old"], 3) != round(r["new"], 3)]
P11 = {"SN_alone_quoted_unchanged": all(a35[t]["SN_quoted_old"] == a35[t]["SN_quoted_new"] for t in a35), "ratios_quoted_unchanged": all(a35[t]["ratio_quoted_old"] == a35[t]["ratio_quoted_new"] for t in a35),
       "A26_quoted_unchanged": a26["quoted_old"] == a26["quoted_new"] and a26["pct_old"] == a26["pct_new"], "all_maxima_same_node": all(node_changes), "max_abs_dlnZ_at_max": lnz_worst, "max_abs_rel_change_any_bound": worst,
       "quoted_3dec_changes_anywhere": quoted_changes}
P11["holds"] = bool(P11["SN_alone_quoted_unchanged"] and P11["ratios_quoted_unchanged"] and P11["A26_quoted_unchanged"] and P11["all_maxima_same_node"] and lnz_worst < 1e-3)
out = {"rows": rows, "A35": a35, "A26": a26, "P11": P11, "reproduction_old_rows": {"exacto_sn_max_abs": max(r.get("repro_old_abs", 0) for r in rows), "maximos_16h_max_abs": max(r.get("repro_old_abs_16h", 0) for r in rows), "maximos_1h_max_abs": max(r.get("repro_old_abs_1h", 0) for r in rows)}}
(HERE / "CORRECCION_FILA_SN_cierre_v3A.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
print(json.dumps(P11, indent=1)); print(json.dumps(a35, indent=1)); print(json.dumps(a26, indent=1)); print("repro:", out["reproduction_old_rows"])
