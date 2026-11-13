#!/usr/bin/env python3
"""
REVISION 2 (comprobación de pronostico_dr3_lrg2_v2 / PRONOSTICO_DR3_REGISTRO.md). Nada existente se modifica.
1. Densidades exactas de v2 (a, Sigma) frente a mi affine_map de la primera revisión, en los 12 casos y en el barrido (42 puntos).
2. Reproducción de la decisión registrada (densidades exactas + mitad B con las semillas del borrador) con mi simulador.
3. c0, c1, potencia de alta precisión muestreando directamente de N(a, Sigma) (4e6 universos por hipótesis): la distribución de Λ
   está determinada por (a0, Σ0, a1, Σ1), así que esto da el valor "sin error MC" con el que comparar los registrados ± MC.
4. Cotejo numérico de las tablas de PRONOSTICO_DR3_REGISTRO.md con pronostico_dr3_lrg2_v2.json (tolerancia: medio dígito del md).
5. Efecto de la tolerancia de la regla c: C3 de LRG2 dentro de lo permitido (ρ² 0.54 / 0.66 por componente, correlación −0.395 → −0.49)
   pero calibración con ρ²_k·C_kk (media 0.6), como hace el script; tasas reales de error bajo el anidamiento matricial exacto
   n3 = C3 C⁻¹ n2 + e, Cov(e) = C3 − C3 C⁻¹ C3, y lo que daría calibrar con ese anidamiento.
6. Camino ρ² por bin de v2 frente a una generación propia con R = diag(ρ²_fila).
Escribe revision2_pronostico_comprobacion.json.
"""
from __future__ import annotations
import json, re, sys
from pathlib import Path
import numpy as np
from scipy.stats import multivariate_normal
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import revision_pronostico_reproduccion as R   # noqa: E402  (mi primera revisión: setup, case_ops, affine_map, simulate, decision)
import pronostico_dr3_lrg2_v2 as V2            # noqa: E402  (no ejecuta nada al importar)
P = R.P
D, COV13, LRG2, KEEP, C22 = R.D, R.COV13, R.LRG2, R.KEEP, R.C22
J2JSON = json.load(open(HERE / "pronostico_dr3_lrg2_v2.json"))
MD = (HERE / "PRONOSTICO_DR3_REGISTRO.md").read_text(encoding="utf-8")
out = {}


def lam_of(dens):
    g0 = multivariate_normal(mean=dens["H0"][0], cov=dens["H0"][1]); g1 = multivariate_normal(mean=dens["H1"][0], cov=dens["H1"][1])
    return lambda r: g1.logpdf(r) - g0.logpdf(r)


def hp_decision(dens, n=4_000_000, seed=777):
    """c0, c1, power sampling r3 directly from the exact Gaussians (the affine model makes this exact)."""
    rng = np.random.default_rng(seed)
    lam = lam_of(dens)
    L0 = lam(rng.multivariate_normal(dens["H0"][0], dens["H0"][1], n))
    L1 = lam(rng.multivariate_normal(dens["H1"][0], dens["H1"][1], n))
    c0, c1 = float(np.percentile(L0, 95)), float(np.percentile(L1, 5))
    return {"c0": c0, "c1": c1, "power": float(np.mean(L1 > c0)), "P_lt_c1_H0": float(np.mean(L0 < c1)), "P_gt3_H1": float(np.mean(L1 > 3)), "P_ltm3_H0": float(np.mean(L0 < -3)), "n": n}


# ---------------------------------------------------------------------------------------------- 1–3: main cases
sec = {}
for ic, combo in enumerate(("sinSN", "Pantheon+")):
    st = R.setup(combo); st2 = V2.setup(combo)
    for ir, rho2 in enumerate(P.RHO2S):
        key = f"rho2_{rho2:g}"
        jc = J2JSON["combos"][combo]["cases"][key]
        op = R.case_ops(combo, st, rho2); op2 = V2.case_ops(st2, rho2)
        rec = {"V3k_v2_vs_mine": float(np.max(np.abs(op2["V3k"] - op["V3k"])))}
        dens_mine, dens_v2 = {}, {}
        for hyp in ("H0", "H1"):
            a, Sig, lin = R.affine_map(hyp, st, op, rho2)
            a2, Sig2 = np.array(jc["exact"][hyp]["a"]), np.array(jc["exact"][hyp]["Sigma"])
            dens_mine[hyp], dens_v2[hyp] = (a, Sig), (a2, Sig2)
            rec[hyp] = {"a_max_abs_diff": float(np.max(np.abs(a - a2))), "Sigma_max_abs_diff": float(np.max(np.abs(Sig - Sig2))),
                        "a_diff_over_sigma": float(np.max(np.abs(a - a2) / np.sqrt(np.diag(Sig)))), "E_d2_v2": jc["exact"][hyp]["E_d2"]}
        # registered decision reproduced with my simulator and the draft's B seeds
        lam = lam_of(dens_v2)
        B = {hyp: R.simulate(hyp, st, op, rho2, R.N_UNI // 2, np.random.default_rng([R.SEED_BASE, ic, ir, ih, 1])) for ih, hyp in enumerate(("H0", "H1"))}
        dec = R.decision(lam(B["H0"][0]), lam(B["H1"][0]))
        jd = jc["decision"]
        rec["registered_decision_repro_diff"] = {"c0": dec["c0"] - jd["c0"], "c1": dec["c1"] - jd["c1"], "power": dec["power"] - jd["power"],
                                                 "med_d2_H0": float(np.median(B["H0"][1])) - jd["med_d2_H0"], "med_d2_H1": float(np.median(B["H1"][1])) - jd["med_d2_H1"]}
        hp = hp_decision(dens_v2)
        e = jc["mc_error"]["stats"]
        rec["high_precision"] = hp
        rec["registered_minus_hp_over_mc_std"] = {k: (jd[k] - hp[k]) / e[k]["std"] for k in ("c0", "c1", "power", "P_lt_c1_H0")}
        sec[f"{combo}/{key}"] = rec
        print(f"{combo} {key}: a/Σ diff máx {max(rec[h]['a_max_abs_diff'] for h in ('H0','H1')):.1e}/{max(rec[h]['Sigma_max_abs_diff'] for h in ('H0','H1')):.1e}; "
              f"repro decisión Δ c0={rec['registered_decision_repro_diff']['c0']:.1e} power={rec['registered_decision_repro_diff']['power']:.1e}; "
              f"alta precisión c0={hp['c0']:.4f} c1={hp['c1']:.4f} power={hp['power']:.4f} | registrado {jd['c0']:.3f}±{e['c0']['std']:.3f} {jd['c1']:.3f}±{e['c1']['std']:.3f} {jd['power']:.4f}±{e['power']['std']:.4f} "
              f"(desv. en σ_MC: {np.round(list(rec['registered_minus_hp_over_mc_std'].values()), 2).tolist()})", flush=True)
    # sweep
    sw = {}
    for k, rho2 in enumerate(V2.RHO2_GRID):
        js = J2JSON["combos"][combo]["rho2_sweep"][f"{rho2:.2f}"]
        op = R.case_ops(combo, st, rho2)
        dens = {hyp: R.affine_map(hyp, st, op, rho2)[:2] for hyp in ("H0", "H1")}
        case, _, dens_v2s, _, _ = V2.calibrate_case(st2, rho2, {hyp: [R.SEED_BASE, ic, V2.SEED_GRID + k, ih, 1] for ih, hyp in enumerate(("H0", "H1"))})
        hp = hp_decision(dens, n=1_000_000, seed=1000 + k)
        sw[f"{rho2:.2f}"] = {"a_Sigma_max_abs_diff_vs_v2": float(max(np.max(np.abs(dens[h][0] - dens_v2s[h][0])) for h in dens) + max(np.max(np.abs(dens[h][1] - dens_v2s[h][1])) for h in dens)),
                            "v2_calibrate_case_repro_diff": {q: case["decision"][q] - js[q] for q in ("c0", "c1", "power")},
                            "hp": hp, "registered": {q: js[q] for q in ("c0", "c1", "power")}, "reg_minus_hp": {q: js[q] - hp[q] for q in ("c0", "c1", "power")}}
    sec[f"{combo}/sweep"] = sw
    print(f"{combo} barrido: máx |a,Σ v2 − mías| {max(v['a_Sigma_max_abs_diff_vs_v2'] for v in sw.values()):.1e}; máx |repro| {max(abs(x) for v in sw.values() for x in v['v2_calibrate_case_repro_diff'].values()):.1e}; "
          f"registrado − alta precisión: c0 máx {max(abs(v['reg_minus_hp']['c0']) for v in sw.values()):.3f}, potencia máx {max(abs(v['reg_minus_hp']['power']) for v in sw.values()):.4f}", flush=True)
out["1_3_densidades_y_decisiones"] = sec

# ---------------------------------------------------------------------------------------------- 4: md tables vs json
NUM = re.compile(r"[-+]?\d+\.\d+(?:[eE][-+]?\d+)?")
# table-aware parser: tables identified by a substring of their header line
TABLES = {}
cur = None
for line in MD.splitlines():
    if line.startswith("|---"):
        continue
    if line.startswith("|"):
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if cur is None or cells[0] in ("combinación", "ρ² común"):
            cur = line; TABLES[cur] = {}
            continue
        TABLES[cur][tuple(cells[:3])] = line
    else:
        cur = None


def table(header_sub):
    for h, rows in TABLES.items():
        if header_sub in h:
            return rows
    raise KeyError(header_sub)


def cmp(header_sub, key, expected):
    rows = table(header_sub)
    line = None
    for k, l in rows.items():
        if k[:len(key)] == tuple(key):
            line = l; break
    rec = {"table": header_sub, "row": key, "found": line is not None}
    if line is None:
        return rec
    s = re.sub(r"\(Δ[^)]*\)", "", line)
    toks = NUM.findall(s)
    nums = [float(t) for t in toks]; decs = [len(t.split(".")[1].split("e")[0].split("E")[0]) for t in toks]
    rec["n_md"], rec["n_json"] = len(nums), len(expected)
    rec["mismatches"] = [{"i": i, "md": m, "json": x} for i, (m, dcs, x) in enumerate(zip(nums, decs, expected)) if abs(m - x) > 0.5 * 10.0 ** (-dcs) + 1e-9]
    rec["count_mismatch"] = len(nums) != len(expected)
    return rec


tab = []
ORIGJ = json.load(open(HERE / "pronostico_dr3_lrg2.json"))
for combo in ("sinSN", "Pantheon+"):
    jc = J2JSON["combos"][combo]
    for rho2 in P.RHO2S:
        c = jc["cases"][f"rho2_{rho2:g}"]; d = c["decision"]; e = c["mc_error"]["stats"]; o = ORIGJ["combos"][combo]["cases"][f"rho2_{rho2:g}"]["decision"]
        tab.append(cmp("| c0 | c1 |", [combo, f"{rho2:g}"], [rho2, d["c0"], e["c0"]["std"], d["c1"], e["c1"]["std"], d["power"], e["power"]["std"], d["P_lt_c1_H0"], e["P_lt_c1_H0"]["std"],
                                                           d["inconclusive_H0"], e["inconclusive_H0"]["std"], d["inconclusive_H1"], e["inconclusive_H1"]["std"]]))
        r = c["draft_reproduction"]
        tab.append(cmp("| borrador", [combo, f"{rho2:g}"], [rho2, o["c0_pct95_H0"], o["c1_pct5_H1"], o["power_P_Lambda_gt_c0_given_H1"], r["c0"], r["c1"], r["power"], d["c0"], d["c1"], d["power"], e["c0"]["std"], e["c1"]["std"], e["power"]["std"]]))
        for hyp in ("H0", "H1"):
            sv = c["sim_vs_exact"][hyp]; ex = c["exact"][hyp]
            tab.append(cmp("a (media exacta)", [combo, f"{rho2:g}", hyp], [rho2] + ex["a"] + list(np.array(ex["Sigma"]).ravel()) + sv["AB_mean_minus_exact_over_sigma"] + [100 * v for v in sv["A_var_over_exact_minus_1"]] + sv["draft_A_mean_minus_exact_over_sigma"]))
            tab.append(cmp("d²_DR3 pct", [combo, f"{rho2:g}", hyp], [rho2] + d["d2_pct_5_50_95"][hyp] + [d[f"med_d2_{hyp}"], e[f"med_d2_{hyp}"]["std"], ex["E_d2"]]))
            tab.append(cmp("Λ pct", [combo, f"{rho2:g}", hyp], [rho2] + d["Lambda_pct_5_50_95"][hyp] + [d[f"P_gt3_{hyp}"], e[f"P_gt3_{hyp}"]["std"], d[f"P_ltm3_{hyp}"], e[f"P_ltm3_{hyp}"]["std"]]))
for rho2 in V2.RHO2_GRID:
    a, b = J2JSON["combos"]["sinSN"]["rho2_sweep"][f"{rho2:.2f}"], J2JSON["combos"]["Pantheon+"]["rho2_sweep"][f"{rho2:.2f}"]
    tab.append(cmp("ρ² común", [f"{rho2:.2f}"], [rho2, a["c0"], a["c1"], a["power"], a["P_lt_c1_H0"], a["inconclusive_H0"], a["inconclusive_H1"], a["med_d2_H0"], a["med_d2_H1"],
                                              b["c0"], b["c1"], b["power"], b["med_d2_H0"], b["med_d2_H1"]]))
pr = J2JSON["predictions"]
cal = json.load(open(HERE / "pronostico_dr3_lrg2_v2_calibracion_falso.json")); res = json.load(open(HERE / "pronostico_dr3_lrg2_v2_resultado_falso.json"))
text_checks = {"P18_power_0.656±0.003": [pr["P18"]["power"], pr["P18"]["mc_std"]], "P18_Pantheon_0.685": pr["P18"]["power_Pantheon+"], "P19": [pr["P19"]["median_H0"], pr["P19"]["median_H1"], pr["P19"]["mc_std"]],
               "falso_c0_c1_sinSN_(+1.013/-1.599)": [cal["combos"]["sinSN"]["c0"], cal["combos"]["sinSN"]["c1"]], "falso_c0_c1_P+_(+0.978/-1.573)": [cal["combos"]["Pantheon+"]["c0"], cal["combos"]["Pantheon+"]["c1"]],
               "falso_r3_d2_Lambda_sinSN_(-0.182,-0.620;10.38;+1.127)": [res["combos"]["sinSN"]["r3"], res["combos"]["sinSN"]["d2_DR3"], res["combos"]["sinSN"]["Lambda"]], "falso_Lambda_P+_(1.169)": res["combos"]["Pantheon+"]["Lambda"]}
bad_rows = [t for t in tab if (not t.get("found")) or t.get("mismatches") or t.get("count_mismatch")]
out["4_tablas_md_vs_json"] = {"n_rows_checked": len(tab), "n_rows_with_problems": len(bad_rows), "problems": bad_rows, "text_checks": text_checks,
                             "max_AB_mean_over_sigma_all_cases": float(max(abs(x) for cb in J2JSON["combos"].values() for c in cb["cases"].values() for h in c["sim_vs_exact"].values() for x in h["AB_mean_minus_exact_over_sigma"])),
                             "max_A_var_rel_all_cases": float(max(abs(x) for cb in J2JSON["combos"].values() for c in cb["cases"].values() for h in c["sim_vs_exact"].values() for x in h["A_var_over_exact_minus_1"])),
                             "max_Ed2_sim_minus_exact": float(max(abs(h["d2_mean_B_vs_exact_E_d2"][0] - h["d2_mean_B_vs_exact_E_d2"][1]) for cb in J2JSON["combos"].values() for c in cb["cases"].values() for h in c["sim_vs_exact"].values()))}
print(f"tablas md vs json: {len(tab)} filas comprobadas, {len(bad_rows)} con problemas: {json.dumps(bad_rows)[:2000]}")
print("text_checks:", json.dumps(text_checks))
print("máx |media AB − exacta|/σ:", out["4_tablas_md_vs_json"]["max_AB_mean_over_sigma_all_cases"], "máx var rel:", out["4_tablas_md_vs_json"]["max_A_var_rel_all_cases"], "máx |E d2 sim − exacta|:", out["4_tablas_md_vs_json"]["max_Ed2_sim_minus_exact"])


# ---------------------------------------------------------------------------------------------- 5–6: general (matrix) nesting
def gen_affine_general(hyp, st, Rm, Le, cov3):
    """r3 = a + M z under the nesting n3 = Rm n2 + e (Cov(e) = Le Le^T), refit with cov3 (linearised). Returns (a, Sigma, V3k)."""
    npar, J, J2 = st["npar"], st["J"], st["J2"]
    _, _, S3, CJ3, _ = P.fisher(st.get("combo", "sinSN"), KEEP, st["p2"], cov3)
    G3 = S3 @ CJ3.T
    G3_bao, const_rest = G3[:, :len(KEEP)], G3[:, len(KEEP):] @ (st["d_rest"] - st["m2_rest"])
    V3k = cov3[np.ix_(LRG2, LRG2)] + J2 @ S3 @ J2.T
    Lp = st["L_full"] if hyp == "H0" else st["L2"]
    p_mean = st["p_full"][:npar] if hyp == "H0" else st["p2"][:npar]
    nz = npar + 13 + (2 if hyp == "H1" else 0)

    def r3_of(z):
        p_true = p_mean + Lp @ z[:npar]
        m_true = st["m2_bao"] + J[:13] @ (p_true - st["p2"][:npar])
        x = m_true.copy()
        if hyp == "H1":
            x[LRG2] += D[LRG2] - m_true[LRG2] + st["L22"] @ z[npar + 13:]
        d3 = x + Rm @ (D[:13] - x) + Le @ z[npar:npar + 13]
        dp3 = G3_bao @ (d3[KEEP] - st["m2_bao"][KEEP]) + const_rest
        return d3[LRG2] - st["m2_bao"][LRG2] - J2 @ dp3
    a = r3_of(np.zeros(nz)); M = np.column_stack([r3_of(np.eye(nz)[k]) - a for k in range(nz)])
    return a, M @ M.T, V3k


st = R.setup("sinSN"); st["combo"] = "sinSN"; st2 = V2.setup("sinSN")
# 6: per-bin rho2 path of v2 vs my general generator with R = diag
rb = [0.5, 0.6, 0.6, 0.7, 0.6, 0.6, 0.6]
rr = V2.rho2_rows(rb)
cov3 = COV13 * np.outer(np.sqrt(rr), np.sqrt(rr))
Le = np.linalg.cholesky(cov3 - cov3 @ np.linalg.solve(COV13, cov3))     # = diag(sqrt(rr(1-rr))) C diag(...) for proportional blocks
op2 = V2.case_ops(st2, rb)
sec6 = {}
for hyp in ("H0", "H1"):
    a2, S2v, _, _ = V2.exact_density(hyp, st2, op2)
    a, Sg, _ = gen_affine_general(hyp, st, np.diag(rr), Le, cov3)
    sec6[hyp] = {"a_max_abs_diff": float(np.max(np.abs(a - a2))), "Sigma_max_abs_diff": float(np.max(np.abs(Sg - S2v)))}
out["6_rho2_por_bin_v2_vs_general"] = {"rho2_bins": rb, **sec6}
print("ρ² por bin, v2 frente a general:", sec6)

# 5: rule-c tolerance effect
c11, c22 = COV13[3, 3], COV13[4, 4]
cases5 = {}
for label, (rDM, rDH, corr3) in {"dentro_tolerancia_0.54_0.66_corr-0.49": (0.54, 0.66, -0.49), "dentro_tolerancia_0.66_0.54_corr-0.30": (0.66, 0.54, -0.30), "solo_corr-0.49": (0.6, 0.6, -0.49), "solo_spread_0.54_0.66": (0.54, 0.66, -0.3953), "tol_estricta_0.585_0.615_corr-0.425": (0.585, 0.615, -0.425), "tol_estricta_0.615_0.585_corr-0.365": (0.615, 0.585, -0.365)}.items():
    C3 = 0.6 * COV13
    C3[3, 3], C3[4, 4] = rDM * c11, rDH * c22
    C3[3, 4] = C3[4, 3] = corr3 * np.sqrt(C3[3, 3] * C3[4, 4])
    rc = V2.rupture_checks(C3)
    Rm = C3 @ np.linalg.inv(COV13)
    Ce = C3 - C3 @ np.linalg.solve(COV13, C3)
    ev = np.linalg.eigvalsh(0.5 * (Ce + Ce.T))
    Le = np.linalg.cholesky(Ce + 1e-15 * np.eye(13)) if ev.min() > -1e-12 else None
    rec = {"rupture_status": rc["status"], "broken": rc["broken"], "rho2_per_bin_LRG2": rc["rho2_per_bin"]["LRG2"], "spread_LRG2": rc["component_rel_spread"]["LRG2"], "corr_change": rc["corr_DM_DH_LRG2"]["change"], "Cov_e_min_eig": float(ev.min())}
    if Le is not None:
        rho2_bins = [rc["rho2_per_bin"][n] for n in V2.BIN_NAMES]
        cal, opc, dens_cal, lam_cal, _ = V2.calibrate_case(st2, rho2_bins, {hyp: [R.SEED_BASE, 0, V2.SEED_EVAL, ih, 1] for ih, hyp in enumerate(("H0", "H1"))})
        c0, c1 = cal["decision"]["c0"], cal["decision"]["c1"]
        true = {hyp: gen_affine_general(hyp, st, Rm, Le, C3)[:2] for hyp in ("H0", "H1")}
        rng = np.random.default_rng(55)
        n = 2_000_000
        L0 = lam_cal(rng.multivariate_normal(true["H0"][0], true["H0"][1], n)); L1 = lam_cal(rng.multivariate_normal(true["H1"][0], true["H1"][1], n))
        rec["script_calibration"] = {"c0": c0, "c1": c1, "nominal_power": cal["decision"]["power"]}
        rec["true_rates_with_script_calibration"] = {"false_alarm_P_gt_c0_H0": float(np.mean(L0 > c0)), "power_P_gt_c0_H1": float(np.mean(L1 > c0)), "P_lt_c1_H1": float(np.mean(L1 < c1)), "P_lt_c1_H0": float(np.mean(L0 < c1))}
        # what calibrating with the matrix nesting would give (same B size as the script, but from exact Gaussians)
        hp = hp_decision(true, n=n, seed=56)
        rec["matrix_nesting_calibration"] = {k: hp[k] for k in ("c0", "c1", "power", "P_lt_c1_H0")}
        rec["density_shift"] = {hyp: {"mean_shift_over_sigma": ((true[hyp][0] - dens_cal[hyp][0]) / np.sqrt(np.diag(dens_cal[hyp][1]))).tolist(), "var_ratio": (np.diag(true[hyp][1]) / np.diag(dens_cal[hyp][1])).tolist()} for hyp in true}
    cases5[label] = rec
    print(f"regla c, {label}: estado {rc['status']} {rc['broken']}; " + (f"calibración script c0={rec['script_calibration']['c0']:.3f} potencia nominal {rec['script_calibration']['nominal_power']:.3f}; reales: falsa alarma {rec['true_rates_with_script_calibration']['false_alarm_P_gt_c0_H0']:.4f}, potencia {rec['true_rates_with_script_calibration']['power_P_gt_c0_H1']:.4f}, P(Λ<c1|H1) {rec['true_rates_with_script_calibration']['P_lt_c1_H1']:.4f}; anidamiento matricial: c0={rec['matrix_nesting_calibration']['c0']:.3f} potencia {rec['matrix_nesting_calibration']['power']:.3f}" if Le is not None else "Cov(e) no PSD"))
out["5_regla_c_tolerancia"] = cases5
(HERE / "revision2_pronostico_comprobacion.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")
print("escrito revision2_pronostico_comprobacion.json")
