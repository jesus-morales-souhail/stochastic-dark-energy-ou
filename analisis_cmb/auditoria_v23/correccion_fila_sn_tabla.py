#!/usr/bin/env python3
"""Builds CORRECCION_FILA_SN_v3A.json: (1) reproduction check of the old-A runs against the original JSONs; (2) old/new/relative-change
table for every A-dependent figure of positividad_lin, familias (linear), muestras_sn (linear) and inyeccion; (3) P7-P9."""
import json
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; R = HERE / "v3A_runs"
def J(p): return json.load(open(p))
out = {"reproduction_old_A": {}, "tables": {}, "predictions": {}}
def rel(a, b):
    if a is None or b is None: return None
    return (b - a) / a if a != 0 else (b - a)
def maxdiff(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float); return float(np.max(np.abs(x - y)))
# ---- positividad_lin
o, oldr, new = J(HERE / "positividad_lin.json"), J(R / "positividad_lin_old/positividad_lin.json"), J(R / "positividad_lin_new/positividad_lin.json")
rep = {}; rows = []
for k in o["thetas"]:
    rep[k] = {"bound_lin_192": abs(oldr["thetas"][k]["bound_lin_192"] - o["thetas"][k]["bound_lin_192"]), "bound_pos_lin_192": abs(oldr["thetas"][k]["bound_pos_lin_192"] - o["thetas"][k]["bound_pos_lin_192"]), "lnZ_lin_rel": maxdiff(oldr["thetas"][k]["lnZ_lin_rel"], o["thetas"][k]["lnZ_lin_rel"])}
    for q, lab in (("bound_lin_192", "lineal −1.92"), ("bound_lin_135", "lineal −1.35"), ("bound_pos_lin_192", "con ρ_DE > 0, −1.92")):
        a, b = o["thetas"][k][q], new["thetas"][k][q]; rows.append({"cell": f"{k} {lab}", "old": a, "new": b, "rel": rel(a, b)})
    a, b = o["thetas"][k]["max_lin_at_sigma"], new["thetas"][k]["max_lin_at_sigma"]; rows.append({"cell": f"{k} σ del máximo", "old": a, "new": b, "rel": rel(a, b)})
    a, b = max(o["thetas"][k]["lnZ_lin_rel"]), max(new["thetas"][k]["lnZ_lin_rel"]); rows.append({"cell": f"{k} máx ln Z − ln Z(0)", "old": a, "new": b, "abs": b - a})
    a, b = o["thetas"][k]["frac_neg_at_lin_ceiling"], new["thetas"][k]["frac_neg_at_lin_ceiling"]; rows.append({"cell": f"{k} ρ_DE ≤ 0 en el techo lineal", "old": a, "new": b, "abs": b - a})
if "prediction_3" in o and "prediction_3" in new:
    rows.append({"cell": "P_prior(ρ_DE>0 | σ = cota pos θ = 1)", "old": o["prediction_3"]["P_prior_pos"], "new": new["prediction_3"]["P_prior_pos"], "abs": new["prediction_3"]["P_prior_pos"] - o["prediction_3"]["P_prior_pos"]})
out["reproduction_old_A"]["positividad_lin"] = {"max_abs_diff": max(max(v.values()) for v in rep.values()), "per_theta": rep}
out["tables"]["positividad_lin"] = rows
b8 = [r for r in rows if "−1.92" in r["cell"] and "−1.35" not in r["cell"]]
out["predictions"]["P7"] = {"n_bounds": len(b8), "max_abs_rel_change": max(abs(r["rel"]) for r in b8), "holds_lt_1e-3": max(abs(r["rel"]) for r in b8) < 1e-3}
# ---- familias
o, oldr, new = J(HERE / "familias.json")["linear"], J(R / "familias_old/familias_linear.json"), J(R / "familias_new/familias_linear.json")
rep = {}; rows = []; mx = 0.0; mxl = 0.0
for k in o:
    rep[k] = {"192": (abs(oldr[k]["192"] - o[k]["192"]) if o[k]["192"] is not None and oldr[k]["192"] is not None else (0.0 if o[k]["192"] == oldr[k]["192"] else float("nan"))), "lnZ_rel": maxdiff(oldr[k]["lnZ_rel"], o[k]["lnZ_rel"])}
    for q in ("192", "135"):
        a, b = o[k][q], new[k][q]; r_ = rel(a, b) if (a not in (None, -1.0) and b not in (None, -1.0)) else None
        rows.append({"cell": f"{k} {q}", "old": a, "new": b, "rel": r_}); mx = max(mx, abs(r_)) if r_ is not None else mx
    a, b = o[k]["lnZ_max_rel"], new[k]["lnZ_max_rel"]; rows.append({"cell": f"{k} máx ln Z", "old": a, "new": b, "abs": b - a}); mxl = max(mxl, abs(b - a))
out["reproduction_old_A"]["familias_linear"] = {"max_abs_diff": max(max(v.values()) for v in rep.values()), "per_cell": rep}
out["tables"]["familias_linear"] = rows; fam_mx, fam_mxl = mx, mxl
# ---- muestras_sn
o, oldr, new = J(HERE / "muestras_sn.json")["linear"], J(R / "muestras_sn_old/muestras_sn_linear.json"), J(R / "muestras_sn_new/muestras_sn_linear.json")
rep = {}; rows = []; mx = 0.0; mxl = 0.0
for k in o:
    rep[k] = {"192": (abs(oldr[k]["192"] - o[k]["192"]) if o[k]["192"] is not None and oldr[k]["192"] is not None else 0.0), "lnZ_rel": maxdiff(oldr[k]["lnZ_rel"], o[k]["lnZ_rel"])}
    for q in ("192", "135"):
        a, b = o[k][q], new[k][q]; r_ = rel(a, b) if (a not in (None, -1.0) and b not in (None, -1.0)) else None
        rows.append({"cell": f"{k} {q}", "old": a, "new": b, "rel": r_}); mx = max(mx, abs(r_)) if r_ is not None else mx
    a, b = o[k]["lnZ_max_rel"], new[k]["lnZ_max_rel"]; rows.append({"cell": f"{k} máx ln Z", "old": a, "new": b, "abs": b - a}); mxl = max(mxl, abs(b - a))
out["reproduction_old_A"]["muestras_sn_linear"] = {"max_abs_diff": max(max(v.values()) for v in rep.values()), "per_cell": rep}
out["tables"]["muestras_sn_linear"] = rows
out["predictions"]["P8"] = {"familias_max_abs_rel_bound_change": fam_mx, "familias_max_abs_lnZmax_change": fam_mxl, "muestras_max_abs_rel_bound_change": mx, "muestras_max_abs_lnZmax_change": mxl,
                            "holds": bool(max(fam_mx, mx) < 1e-2 and max(fam_mxl, mxl) < 0.05)}
# ---- inyeccion
o, oldr, new = J(HERE / "inyeccion.json"), J(R / "inyeccion_old/inyeccion.json"), J(R / "inyeccion_new/inyeccion.json")
rep = {"real_check_marg": abs(oldr["real_check"]["marg_linear"] - o["real_check"]["marg_linear"]), "real_check_5pct": abs(oldr["real_check"]["5pct_linear"] - o["real_check"]["5pct_linear"]), "delta": abs(oldr["delta"] - o["delta"])}
for s in o["scenarios"]:
    rep[f"{s}_marg_median"] = abs(oldr["scenarios"][s]["marg"]["median"] - o["scenarios"][s]["marg"]["median"]); rep[f"{s}_b5_median"] = abs(oldr["scenarios"][s]["5pct"]["median"] - o["scenarios"][s]["5pct"]["median"])
    rep[f"{s}_frac_excl0"] = abs(oldr["scenarios"][s]["frac_excl0"] - o["scenarios"][s]["frac_excl0"])
    if o["scenarios"][s].get("coverage_marg") is not None: rep[f"{s}_coverage_marg"] = abs(oldr["scenarios"][s]["coverage_marg"] - o["scenarios"][s]["coverage_marg"])
rep["exact_check_ratio_min"] = abs(oldr["exact_check"]["ratio_min"] - o["exact_check"]["ratio_min"])
out["reproduction_old_A"]["inyeccion"] = {"max_abs_diff": max(rep.values()), "items": rep}
rows = []; cov_changes = []; excl = None
for s in o["scenarios"]:
    so, sn_ = o["scenarios"][s], new["scenarios"][s]
    for q, lab in (("marg", "marginal −1.92 mediana"), ("5pct", "5 % mediana")):
        a, b = so[q]["median"], sn_[q]["median"]; rows.append({"cell": f"{s} {lab}", "old": a, "new": b, "rel": rel(a, b)})
    if so.get("coverage_marg") is not None:
        a, b = so["coverage_marg"], sn_["coverage_marg"]; rows.append({"cell": f"{s} cobertura marginal", "old": a, "new": b, "abs": b - a, "universes": round((b - a) * 200)}); cov_changes.append(abs(round((b - a) * 200)))
        a, b = so["coverage_5pct"], sn_["coverage_5pct"]; rows.append({"cell": f"{s} cobertura 5 %", "old": a, "new": b, "abs": b - a})
    a, b = so["frac_excl0"], sn_["frac_excl0"]; rows.append({"cell": f"{s} σ_X = 0 excluido", "old": a, "new": b, "abs": b - a})
    if s == "C": excl = (a, b)
    a, b = so["ratio_vs_A"]["marg"], sn_["ratio_vs_A"]["marg"]
    if a is not None and b is not None: rows.append({"cell": f"{s} marginal /A", "old": a, "new": b, "rel": rel(a, b)})
    a, b = so["dlnL_max"]["median"], sn_["dlnL_max"]["median"]; rows.append({"cell": f"{s} Δln Z máx mediana", "old": a, "new": b, "abs": b - a})
for q in ("marg_linear", "5pct_linear", "dlnL_max_linear"):
    a, b = o["real_check"][q], new["real_check"][q]; rows.append({"cell": f"datos reales {q}", "old": a, "new": b, "rel": rel(a, b)})
rows.append({"cell": "comprobación exacta: cociente mín/máx", "old": [o["exact_check"]["ratio_min"], o["exact_check"]["ratio_max"]], "new": [new["exact_check"]["ratio_min"], new["exact_check"]["ratio_max"]]})
rows.append({"cell": "δ del escenario C", "old": o["delta"], "new": new["delta"], "rel": rel(o["delta"], new["delta"])})
out["tables"]["inyeccion"] = rows
out["predictions"]["P9"] = {"max_coverage_change_universes": max(cov_changes), "C_frac_excl0_old_new": excl, "C_frac_excl0_change_points": 100 * (excl[1] - excl[0]), "holds": bool(max(cov_changes) <= 1 and abs(100 * (excl[1] - excl[0])) <= 2)}
out["predictions"]["P10"] = J(HERE / "interpolacion_mu_check.json")["P10_max_lt_1e-3_mag"]
(HERE / "CORRECCION_FILA_SN_v3A.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
print("reproduction old A, max abs diff per target:", {k: v["max_abs_diff"] for k, v in out["reproduction_old_A"].items()})
print(json.dumps(out["predictions"], indent=1))
for t in ("positividad_lin", "inyeccion"):
    for r in out["tables"][t]: print(t, r["cell"], r["old"], r["new"], r.get("rel", r.get("abs")))
