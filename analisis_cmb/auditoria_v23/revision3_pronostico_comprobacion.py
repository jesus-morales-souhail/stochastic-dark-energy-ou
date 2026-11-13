#!/usr/bin/env python3
"""
REVISION 3: comprobación de pronostico_dr3_lrg2_v3 / v4 y sus registros. Nada existente se modifica; escribe revision3_pronostico_comprobacion.json.
1. Anidamiento matricial: identidades Cov(n3, n2) = C3, Cov(n3) = C3; equivalencia Cov(e) ⪰ 0 ⟺ C2 − C3 ⪰ 0 (aleatorio); R y Fe de v3.nesting.
2. v3: densidades frente a mi affine_map (revisión 1) y umbrales de alta precisión frente a mi muestreo independiente.
3. v4: modelo en z arbitrario (nodos frente a la malla; interpolación de D_M frente a integral directa; D_H cerrado frente a la malla); ajustador en z2
   frente al de v3; densidades a dz = 0 frente a v3 (en unidades absolutas); unidades fraccionales (a_frac · m = a_abs).
4. Conversión de μ (foco): relación exacta μ3 = μ2 ⊙ m(p_true, z3)/m(p_true, z2); término omitido μ2 ⊙ [ρ(p_true) − ρ(p2)] sobre 10⁵ muestras del
   posterior de (p, μ) con el modelo exacto; efecto en c1 y potencia con la expresión completa (ρ exacto por universo) frente a la linealizada, muestras
   pareadas, dz = ±0.01, en unidades del error MC.
5. Denominador m(p3, z3) (evaluación) frente a m(p2, z3) (densidades): efecto en r3, Λ en la prueba falsa z0712 y en universos simulados.
6. Tablas de los dos registros frente a sus json.
"""
from __future__ import annotations
import json, re, sys, time
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.linalg import cho_solve, cho_factor
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import revision_pronostico_reproduccion as R    # noqa: E402
import pronostico_dr3_lrg2_v3 as V3             # noqa: E402
import pronostico_dr3_lrg2_v4 as V4             # noqa: E402
P, JB = V3.P, V3.JB
D, COV13, LRG2, KEEP, C22 = V3.D, V3.COV13, V3.LRG2, V3.KEEP, V3.C22
J3S = json.load(open(HERE / "pronostico_dr3_lrg2_v3.json")); J4S = json.load(open(HERE / "pronostico_dr3_lrg2_v4.json"))
MD3 = (HERE / "PRONOSTICO_DR3_REGISTRO_v3.md").read_text(encoding="utf-8"); MD4 = (HERE / "PRONOSTICO_DR3_REGISTRO_v4.md").read_text(encoding="utf-8")
out = {}
t0 = time.time()

# ------------------------------------------------------------------ 1. nesting
rng = np.random.default_rng(3)
sec = {"random_checks": []}
for k in range(200):
    n = rng.integers(1, 4)
    A = rng.standard_normal((n, n)); C2 = A @ A.T + 0.1 * np.eye(n)
    B = rng.standard_normal((n, n)); C3 = C2 - (B @ B.T) * rng.uniform(0, 0.3) if k % 2 == 0 else rng.uniform(0.2, 1.5) * C2 + 0.05 * (B @ B.T)
    C3 = 0.5 * (C3 + C3.T)
    if np.linalg.eigvalsh(C3).min() <= 0:
        continue
    Rm = C3 @ np.linalg.inv(C2); Ce = C3 - C3 @ np.linalg.solve(C2, C3)
    e1 = np.linalg.eigvalsh(0.5 * (Ce + Ce.T)).min(); e2 = np.linalg.eigvalsh(C2 - C3).min()
    sec["random_checks"].append({"Cov_n3n2_eq_C3": float(np.max(np.abs(Rm @ C2 - C3))), "Cov_n3_eq_C3": float(np.max(np.abs(Rm @ C2 @ Rm.T + Ce - C3))),
                                 "min_eig_Cov_e": float(e1), "min_eig_C2_minus_C3": float(e2), "same_sign": bool((e1 >= -1e-12) == (e2 >= -1e-12))})
sec["identities_max_abs"] = float(max(max(c["Cov_n3n2_eq_C3"], c["Cov_n3_eq_C3"]) for c in sec["random_checks"]))
sec["psd_equivalence_holds_in_all"] = bool(all(c["same_sign"] for c in sec["random_checks"]))
sec["n_random"] = len(sec["random_checks"]); sec.pop("random_checks")
rep, Rm, Fe = V3.nesting(COV13, 0.6 * COV13)
sec["v3_nesting_rho2_0.6"] = {"R_minus_0.6I_max": float(np.max(np.abs(Rm - 0.6 * np.eye(13)))), "FeFeT_minus_0.24C2_max": float(np.max(np.abs(Fe @ Fe.T - 0.24 * COV13))), "status": rep["status"]}
# a non-proportional C3 that passes: R C2 = C3 and Fe Fe^T = C3 - C3 C2^-1 C3
C3n = 0.6 * COV13.copy(); C3n[3, 3] *= 0.9; C3n[4, 4] *= 1.1; C3n[3, 4] = C3n[4, 3] = -0.45 * np.sqrt(C3n[3, 3] * C3n[4, 4])
rep, Rm, Fe = V3.nesting(COV13, C3n)
sec["v3_nesting_nonproportional"] = {"status": rep["status"], "RC2_minus_C3_max": float(np.max(np.abs(Rm @ COV13 - C3n))), "FeFeT_minus_Cove_max": float(np.max(np.abs(Fe @ Fe.T - (C3n - C3n @ np.linalg.solve(COV13, C3n)))))}
# MC check of the nesting with that C3: Cov(n3), Cov(n3, n2)
z = rng.standard_normal((400000, 13)); n2 = z @ np.linalg.cholesky(COV13).T; n3 = n2 @ Rm.T + rng.standard_normal((400000, 13)) @ Fe.T
sec["v3_nesting_nonproportional"]["MC_Cov_n3_rel_err_max"] = float(np.max(np.abs(np.cov(n3.T) - C3n) / np.sqrt(np.outer(np.diag(C3n), np.diag(C3n)))))
sec["v3_nesting_nonproportional"]["MC_Cov_n3n2_rel_err_max"] = float(np.max(np.abs((n3 - n3.mean(0)).T @ (n2 - n2.mean(0)) / (len(n2) - 1) - C3n) / np.sqrt(np.outer(np.diag(C3n), np.diag(COV13)))))
# tolerance behaviour: rounding negative vs real negative
Cr = COV13.copy(); Cr3 = 0.6 * COV13; Cr3[4, 4] = COV13[4, 4] * (1 + 1e-12)   # C2 - C3 with a -1e-12-level eigenvalue in LRG2? (D_H var slightly > DR2's)
rep_r, _, _ = V3.nesting(COV13, Cr3)
sec["tolerance_test_DH_var_1+1e-12_of_DR2"] = {"status": rep_r["status"], "LRG2_C2_minus_C3_real_negative": rep_r["bins"]["LRG2"]["C2_minus_C3"]["real_negative"], "rounding_negative": rep_r["bins"]["LRG2"]["C2_minus_C3"]["rounding_negative"]}
out["1_anidamiento"] = sec
print("1 anidamiento:", json.dumps(sec)[:900], flush=True)

# ------------------------------------------------------------------ 2. v3 densities and hp vs mine
sec = {}
for ic, combo in enumerate(("sinSN", "Pantheon+")):
    st = R.setup(combo); st3 = V3.setup(combo)
    for ir, rho2 in enumerate(P.RHO2S):
        jc = J3S["combos"][combo]["cases"][f"rho2_{rho2:g}"]
        op = R.case_ops(combo, st, rho2)
        rec = {}
        dens = {}
        for hyp in ("H0", "H1"):
            a, Sg, _ = R.affine_map(hyp, st, op, rho2)
            dens[hyp] = (a, Sg)
            rec[hyp] = {"a_diff": float(np.max(np.abs(a - np.array(jc["exact"][hyp]["a"])))), "Sigma_diff": float(np.max(np.abs(Sg - np.array(jc["exact"][hyp]["Sigma"]))))}
        # my own independent high-precision sample (different seed): 8e6 (2e7 in the principal case), block errors, combined-sigma deviation
        g = {h: V3.Gauss(*dens[h]) for h in dens}
        rg = np.random.default_rng([99, ic, ir])
        nblk = 20 if (combo == "sinSN" and rho2 == 0.6) else 8
        L0 = np.concatenate([g["H1"].logpdf(x) - g["H0"].logpdf(x) for x in (g["H0"].sample(1_000_000, rg) for _ in range(nblk))])
        L1 = np.concatenate([g["H1"].logpdf(x) - g["H0"].logpdf(x) for x in (g["H1"].sample(1_000_000, rg) for _ in range(nblk))])
        c0, c1 = float(np.percentile(L0, 95)), float(np.percentile(L1, 5)); pw = float(np.mean(L1 > c0))
        bl = [(float(np.percentile(L0[k*1_000_000:(k+1)*1_000_000], 95)), float(np.percentile(L1[k*1_000_000:(k+1)*1_000_000], 5))) for k in range(nblk)]
        pwb = [float(np.mean(L1[k*1_000_000:(k+1)*1_000_000] > bl[k][0])) for k in range(nblk)]
        my_err = {"c0": float(np.std([b[0] for b in bl], ddof=1) / np.sqrt(nblk)), "c1": float(np.std([b[1] for b in bl], ddof=1) / np.sqrt(nblk)), "power": float(np.std(pwb, ddof=1) / np.sqrt(nblk))}
        h = jc["hp"]
        comb = {k: float(np.hypot(my_err[k], h["mc_err"][k])) for k in my_err}
        rec["hp_mine_vs_v3"] = {"N_mine": nblk * 1_000_000, "c0": [c0, h["c0"], my_err["c0"], h["mc_err"]["c0"]], "c1": [c1, h["c1"], my_err["c1"], h["mc_err"]["c1"]], "power": [pw, h["power"], my_err["power"], h["mc_err"]["power"]],
                               "dev_in_combined_sigma": {"c0": abs(c0 - h["c0"]) / comb["c0"], "c1": abs(c1 - h["c1"]) / comb["c1"], "power": abs(pw - h["power"]) / comb["power"]}}
        sec[f"{combo}/rho2_{rho2:g}"] = rec
        print(f"2 {combo} {rho2}: a/Σ diff {max(rec[h]['a_diff'] for h in ('H0','H1')):.1e}/{max(rec[h]['Sigma_diff'] for h in ('H0','H1')):.1e}; hp mine {c0:.4f}±{my_err['c0']:.4f}/{c1:.4f}±{my_err['c1']:.4f}/{pw:.4f}±{my_err['power']:.4f} vs v3 {h['c0']:.4f}±{h['mc_err']['c0']:.4f}/{h['c1']:.4f}±{h['mc_err']['c1']:.4f}/{h['power']:.4f}±{h['mc_err']['power']:.4f} (dev σ comb.: {np.round(list(rec['hp_mine_vs_v3']['dev_in_combined_sigma'].values()), 1).tolist()})", flush=True)
out["2_v3_densidades_umbrales"] = sec

# ------------------------------------------------------------------ 3. v4 model at arbitrary z, fitter, densities at dz = 0, units
sec = {}
st3 = V3.setup("sinSN"); p2 = st3["p2"]; Z2 = V4.Z2
m_nodes = V4.model_bao(p2, Z2); m_grid = JB.model_at(p2)[:13]
sec["model_bao_at_DR2_nodes_vs_grid_rel_max"] = float(np.max(np.abs(m_nodes / m_grid - 1)))
h, wm = p2[0], p2[1]; om, orad = wm / h**2, V3.P.S.OMR / h**2
E = lambda zz: np.sqrt(om * (1 + zz)**3 + orad * (1 + zz)**4 + (1 - om - orad))
pref = V3.P.S.C_KMS / (100.0 * h * float(V3.P.S.E.RD_SPL(wm)))
rel = []
for zt in (0.696, 0.701, 0.711, 0.716, 0.52, 0.944, 2.34):
    zz = Z2.copy(); zz[3] = zt; zz[4] = zt
    mb = V4.model_bao(p2, zz)
    dm_q = quad(lambda x: 1 / E(x), 0, zt, epsabs=0, epsrel=1e-12)[0]
    rel.append({"z": zt, "DM_interp_vs_quad_rel": float(mb[3] / (pref * dm_q) - 1), "DH_closed_vs_formula_rel": float(mb[4] / (pref / E(zt)) - 1)})
sec["model_bao_off_node"] = rel
sec["max_abs_rel_DM_interp_vs_quad"] = float(max(abs(r["DM_interp_vs_quad_rel"]) for r in rel))
# fitter at z2 vs v3/original fitter
for combo in ("sinSN", "Pantheon+"):
    chi_a, pa = P.fit_c(combo, KEEP, D, 0.6 * COV13); chi_b, pb = V4.fit_cz(combo, KEEP, D, 0.6 * COV13, Z2)
    sec[f"fit_cz_vs_fit_c_{combo}"] = {"p_max_abs_diff": float(np.max(np.abs(pa - pb))), "chi2_diff": chi_b - chi_a}
# densities dz = 0 vs v3 (absolute units) and fractional units check
for combo in ("sinSN", "Pantheon+"):
    st = V3.setup(combo)
    _, op4 = V4.case_ops(st, 0.6 * COV13, Z2)
    _, op3 = V3.case_ops(st, 0.6 * COV13)
    m2L = op4["m2_bao3"][LRG2]
    r = {"m2_LRG2_at_p2": m2L.tolist()}
    for hyp in ("H0", "H1"):
        a4, S4, _ = V4.exact_density(hyp, st, op4); a3, S3, _ = V3.exact_density(hyp, st, op3)
        r[hyp] = {"a_frac_times_m_minus_a3_abs": float(np.max(np.abs(a4 * m2L - a3))), "Sigma_frac_scaled_minus_Sigma3": float(np.max(np.abs(S4 * np.outer(m2L, m2L) - S3)))}
    r["V3k_frac_scaled_minus_V3k_abs"] = float(np.max(np.abs(op4["V3k"] * np.outer(m2L, m2L) - op3["V3k"])))
    sec[f"dz0_vs_v3_{combo}"] = r
out["3_v4_modelo_unidades"] = sec
print("3 v4:", json.dumps({k: v for k, v in sec.items() if k != "model_bao_off_node"}), flush=True)

# ------------------------------------------------------------------ 4. mu conversion (focus)
st = V3.setup("sinSN"); npar = 2
L2T = st["LF2"]
rng = np.random.default_rng(2024)
N_POST = 100_000
z_p = rng.standard_normal((N_POST, npar)); z_mu = rng.standard_normal((N_POST, 2))
p_true = V3.p_draw(st, "H1", z_p)                                  # N(p2, Sigma2)
sec = {"formula": "H1 fraccional: dato_L = m_L(p, z)(1 + mu_f) + ruido en los dos releases => mu2 = m_L(p_true, z2) mu_f, mu3 = m_L(p_true, z3) mu_f, "
                  "luego mu3_i = mu2_i * m_i(p_true, z3) / m_i(p_true, z2), componente a componente (D_M, D_H). v4 usa el cociente en p2: "
                  "termino omitido Delta_i = mu2_i [rho_i(p_true) - rho_i(p2)], rho_i(p) = m_i(p, z3)/m_i(p, z2)."}
for dz in (-0.01, 0.01):
    z3 = Z2.copy(); z3[LRG2] += dz
    _, op4 = V4.case_ops(st, 0.6 * COV13, z3)
    ratio_p2 = op4["ratio"]
    # exact model at p_true (non-linear) for the ratio and for mu2's mean; chunked
    t1 = time.time()
    ratio_exact = np.empty((N_POST, 2)); m2L_exact = np.empty((N_POST, 2))
    for i in range(N_POST):
        pt = np.array([p_true[i, 0], p_true[i, 1], 0.0])
        mz2 = V4.model_bao(pt, Z2)[LRG2]; mz3 = V4.model_bao(pt, z3)[LRG2]
        ratio_exact[i] = mz3 / mz2; m2L_exact[i] = mz2
    mu2 = D[LRG2][None, :] - m2L_exact + z_mu @ st["L22"].T           # posterior of mu2 given p_true (exact model)
    delta = mu2 * (ratio_exact - ratio_p2[None, :])                   # omitted term, absolute units at z3
    # sigma of r3 (absolute units at z3) from the registered fractional densities at this dz
    key = f"dz_{dz:+g}"
    a_frac = np.array(J4S["P24"][key]["a_fractional"]["H1"]); m3 = np.array(J4S["P24"][key]["m2_bao3_LRG2"])
    # Sigma not stored for P24 in the json: compute with v4
    a4, S4, _ = V4.exact_density("H1", st, op4)
    sig_abs = np.sqrt(np.diag(S4)) * m3
    dfrac = delta / m3[None, :]
    rec = {"z3_LRG2": float(Z2[3] + dz), "ratio_at_p2": ratio_p2.tolist(), "ratio_exact_pct_5_50_95": np.percentile(ratio_exact, [5, 50, 95], axis=0).T.tolist(),
           "ratio_minus_ratio_p2_abs_max": float(np.max(np.abs(ratio_exact - ratio_p2))), "sigma_r3_abs_units": sig_abs.tolist(), "sigma_r3_fractional": (np.sqrt(np.diag(S4))).tolist(),
           "omitted_term_abs": {"median_abs": np.median(np.abs(delta), axis=0).tolist(), "pct99_abs": np.percentile(np.abs(delta), 99, axis=0).tolist(), "max_abs": np.max(np.abs(delta), axis=0).tolist()},
           "omitted_term_over_sigma": {"median": np.median(np.abs(delta) / sig_abs, axis=0).tolist(), "pct99": np.percentile(np.abs(delta) / sig_abs, 99, axis=0).tolist(), "max": np.max(np.abs(delta) / sig_abs, axis=0).tolist()},
           "a_fractional_json_vs_recomputed_max_abs": float(np.max(np.abs(a4 - a_frac))), "seconds_exact_model": time.time() - t1}
    # Effect on c1 and power: paired simulation, linearised ratio vs exact ratio in the mu3 term. r3 affine part from gen_r3; the extra term
    # enters d3 at the LRG2 rows and propagates through the (linear) refit: r3_extra = s3 * (E2 delta)_L - s3 * J3_L G3 (E2 delta)_keep = s3 * delta (LRG2 rows are not in KEEP).
    # ratio_exact(p) is needed per universe: use a first-order expansion in p (finite differences at p2) checked against the exact values above.
    eps = np.array([1e-4, 1e-5])
    dr = np.zeros((2, npar))
    for jj in range(npar):
        pp, pm = np.array(p2, float), np.array(p2, float); pp[jj] += eps[jj]; pm[jj] -= eps[jj]
        rp = V4.model_bao(pp, z3)[LRG2] / V4.model_bao(pp, Z2)[LRG2]; rm = V4.model_bao(pm, z3)[LRG2] / V4.model_bao(pm, Z2)[LRG2]
        dr[:, jj] = (rp - rm) / (2 * eps[jj])
    ratio_lin1 = ratio_p2[None, :] + (p_true - p2[:npar]) @ dr.T
    rec["ratio_first_order_vs_exact_max_abs"] = float(np.max(np.abs(ratio_lin1 - ratio_exact)))
    dens = {hyp: V3.Gauss(*V4.exact_density(hyp, st, op4)[:2]) for hyp in ("H0", "H1")}
    lam = lambda r: dens["H1"].logpdf(r) - dens["H0"].logpdf(r)
    NB = 4_000_000; CH = 500_000
    rs = np.random.default_rng([20261114, 0, 4, 1, 7, 77])
    L_lin, L_ext = [], []
    for _ in range(NB // CH):
        zp = rs.standard_normal((CH, npar)); ze = rs.standard_normal((CH, 13)); zm = rs.standard_normal((CH, 2))
        r_lin = V4.gen_r3("H1", st, op4, zp, ze, zm)
        pt = V3.p_draw(st, "H1", zp); m_true2 = st["m2_bao"][None, :] + (pt - p2[:npar]) @ st["J"][:13].T
        mu2b = D[LRG2][None, :] - m_true2[:, LRG2] + zm @ st["L22"].T
        extra = mu2b * ((pt - p2[:npar]) @ dr.T) * op4["s3"][None, :]       # first-order omitted term, fractional units
        L_lin.append(lam(r_lin)); L_ext.append(lam(r_lin + extra))
    L_lin, L_ext = np.concatenate(L_lin), np.concatenate(L_ext)
    r0 = np.random.default_rng([20261114, 0, 4, 0, 7, 77])
    L0 = np.concatenate([lam(dens["H0"].sample(CH, r0)) for _ in range(NB // CH)])
    c0 = float(np.percentile(L0, 95))
    c1_lin, c1_ext = float(np.percentile(L_lin, 5)), float(np.percentile(L_ext, 5))
    pw_lin, pw_ext = float(np.mean(L_lin > c0)), float(np.mean(L_ext > c0))
    nb = 8; blk = NB // nb
    c1_err = float(np.std([np.percentile(L_lin[k * blk:(k + 1) * blk], 5) for k in range(nb)], ddof=1) / np.sqrt(nb))
    pw_err = float(np.std([np.mean(L_lin[k * blk:(k + 1) * blk] > c0) for k in range(nb)], ddof=1) / np.sqrt(nb))
    rec["effect_on_thresholds"] = {"c0_H0": c0, "c1_linearised": c1_lin, "c1_exact_ratio": c1_ext, "d_c1": c1_ext - c1_lin, "mc_err_c1": c1_err, "d_c1_over_mc": (c1_ext - c1_lin) / c1_err,
                                   "power_linearised": pw_lin, "power_exact_ratio": pw_ext, "d_power": pw_ext - pw_lin, "mc_err_power": pw_err, "d_power_over_mc": (pw_ext - pw_lin) / pw_err,
                                   "Lambda_max_abs_change_per_universe": float(np.max(np.abs(L_ext - L_lin))), "Lambda_rms_change": float(np.sqrt(np.mean((L_ext - L_lin)**2))),
                                   "json_P24": {k: J4S["P24"][key][k] for k in ("c0", "c1", "power")}}
    sec[key] = rec
    print(f"4 dz={dz:+}: ratio p2 {np.round(ratio_p2, 5).tolist()}, |ratio−ratio_p2| máx {rec['ratio_minus_ratio_p2_abs_max']:.2e}; término omitido mediana {np.round(rec['omitted_term_abs']['median_abs'], 7).tolist()} máx {np.round(rec['omitted_term_abs']['max_abs'], 6).tolist()} abs; "
          f"en σ: mediana {np.round(rec['omitted_term_over_sigma']['median'], 6).tolist()} máx {np.round(rec['omitted_term_over_sigma']['max'], 5).tolist()}; "
          f"Δc1 {rec['effect_on_thresholds']['d_c1']:+.5f} ({rec['effect_on_thresholds']['d_c1_over_mc']:+.2f} σ_MC), Δpotencia {rec['effect_on_thresholds']['d_power']:+.6f} ({rec['effect_on_thresholds']['d_power_over_mc']:+.2f} σ_MC), |ΔΛ| máx {rec['effect_on_thresholds']['Lambda_max_abs_change_per_universe']:.2e}", flush=True)
out["4_conversion_mu"] = sec
print("4 extra:", {k: (v["ratio_first_order_vs_exact_max_abs"], v["a_fractional_json_vs_recomputed_max_abs"], v["seconds_exact_model"]) for k, v in sec.items() if k.startswith("dz")})

# ------------------------------------------------------------------ 5. denominator m(p3, z3) vs m(p2, z3)
sec = {}
resf = json.load(open(HERE / "pronostico_dr3_lrg2_v4_resultado_falso_z0712.json")); calf = json.load(open(HERE / "pronostico_dr3_lrg2_v4_calibracion_falso_z0712.json"))
for combo in ("sinSN",):
    c = resf["combos"][combo]; cal = calf["combos"][combo]
    r3_abs = np.array(c["r3_absolute"]); m3 = np.array(c["m3_LRG2"]); m2 = np.array(cal["m2_bao3_LRG2"])
    dens = {h: V3.Gauss(np.array(cal["exact_densities_fractional"][h]["a"]), np.array(cal["exact_densities_fractional"][h]["Sigma"])) for h in ("H0", "H1")}
    lam = lambda r: float(dens["H1"].logpdf(r[None, :])[0] - dens["H0"].logpdf(r[None, :])[0])
    L_p3, L_p2 = lam(r3_abs / m3), lam(r3_abs / m2)
    sig = np.sqrt(np.diag(dens["H0"].Sig))
    sec[f"falso_z0712_{combo}"] = {"m3_over_m2_minus_1": (m3 / m2 - 1).tolist(), "r3_frac_p3_denominator": (r3_abs / m3).tolist(), "r3_frac_p2_denominator": (r3_abs / m2).tolist(),
                                  "diff_r3_over_sigma": ((r3_abs / m3 - r3_abs / m2) / sig).tolist(), "Lambda_p3_denominator": L_p3, "Lambda_p2_denominator": L_p2, "d_Lambda": L_p3 - L_p2, "Lambda_json": c["Lambda"]}
# over simulated universes (H0 and H1, dz = 0.006): linearised dp3 -> m(p3) = m(p2) + J3 dp3
st = V3.setup("sinSN"); z3 = Z2.copy(); z3[LRG2] += 0.006
_, op4 = V4.case_ops(st, 0.6 * COV13, z3)
dens = {hyp: V3.Gauss(*V4.exact_density(hyp, st, op4)[:2]) for hyp in ("H0", "H1")}
lam = lambda r: dens["H1"].logpdf(r) - dens["H0"].logpdf(r)
rs = np.random.default_rng(5)
for hyp in ("H0", "H1"):
    n = 400_000
    zp = rs.standard_normal((n, 2)); ze = rs.standard_normal((n, 13)); zm = rs.standard_normal((n, 2)) if hyp == "H1" else None
    pt = V3.p_draw(st, hyp, zp); dp = pt - st["p2"][:2]
    m_true2 = st["m2_bao"][None, :] + dp @ st["J"][:13].T; m_true3 = op4["m2_bao3"][None, :] + dp @ op4["J3"][:13].T
    x2, x3 = m_true2.copy(), m_true3.copy()
    if hyp == "H1":
        mu2 = D[LRG2][None, :] - m_true2[:, LRG2] + zm @ st["L22"].T; x2[:, LRG2] += mu2; x3[:, LRG2] += mu2 * op4["ratio"][None, :]
    n2 = D[:13][None, :] - x2; d3 = x3 + n2 @ op4["R"].T + ze @ op4["Fe"].T
    dp3 = (d3[:, KEEP] - op4["m2_bao3"][KEEP][None, :]) @ op4["G3_bao"].T + op4["const_rest"][None, :]
    m3L = op4["m2_bao3"][LRG2][None, :] + dp3 @ op4["J3"][LRG2].T
    r_abs = d3[:, LRG2] - m3L
    r_p2, r_p3 = r_abs * op4["s3"][None, :], r_abs / m3L
    Lp2, Lp3 = lam(r_p2), lam(r_p3)
    sig = np.sqrt(np.diag(dens[hyp].Sig))
    sec[f"simulated_{hyp}_dz0.006"] = {"m3_over_m2_minus_1_pct_5_50_95": np.percentile(m3L / op4["m2_bao3"][LRG2][None, :] - 1, [5, 50, 95], axis=0).T.tolist(),
                                      "dr3_over_sigma_pct_50_99_max": [np.percentile(np.abs(r_p3 - r_p2) / sig, [50, 99], axis=0).T.tolist(), (np.max(np.abs(r_p3 - r_p2), axis=0) / sig).tolist()],
                                      "dLambda_pct_50_99_max": [float(np.percentile(np.abs(Lp3 - Lp2), 50)), float(np.percentile(np.abs(Lp3 - Lp2), 99)), float(np.max(np.abs(Lp3 - Lp2)))]}
out["5_denominador"] = sec
print("5 denominador:", json.dumps(sec)[:1500], flush=True)

# ------------------------------------------------------------------ 6. tables vs json
NUM = re.compile(r"[-+]?\d+\.\d+(?:[eE][-+]?\d+)?")


def parse_tables(MD):
    T, cur = {}, None
    for line in MD.splitlines():
        if line.startswith("|---"):
            continue
        if line.startswith("|"):
            cells = [c.strip() for c in line.strip().strip("|").split("|")]
            if cur is None:
                cur = line; T[cur] = {}; continue
            T[cur][tuple(cells[:3])] = line
        else:
            cur = None
    return T


def cmp(T, header_sub, key, expected, label=""):
    rows = next((r for h, r in T.items() if header_sub in h), None)
    rec = {"table": header_sub, "row": key, "found": rows is not None}
    if rows is None:
        return rec
    line = next((l for k, l in rows.items() if k[:len(key)] == tuple(key)), None)
    rec["found"] = line is not None
    if line is None:
        return rec
    toks = NUM.findall(re.sub(r"\(Δ[^)]*\)", "", line))
    nums = [float(t) for t in toks]; decs = [len(t.split(".")[1].split("e")[0].split("E")[0]) for t in toks]
    rec["n_md"], rec["n_json"] = len(nums), len(expected)
    rec["mismatches"] = [{"i": i, "md": m, "json": x} for i, (m, dcs, x) in enumerate(zip(nums, decs, expected)) if abs(m - x) > 0.5 * 10.0 ** (-dcs) * (1 if abs(m) < 1e-6 or dcs < 10 else 1) + 1e-9 and not (abs(m) < 1e-10 and abs(x) < 1e-10)]
    rec["count_mismatch"] = len(nums) != len(expected)
    return rec


T3 = parse_tables(MD3); T4 = parse_tables(MD4)
tab = []
for combo in ("sinSN", "Pantheon+"):
    for rho2 in P.RHO2S:
        c = J3S["combos"][combo]["cases"][f"rho2_{rho2:g}"]; h = c["hp"]; e = h["mc_err"]; v2 = c["vs_v2"]
        tab.append(cmp(T3, "| c0 | c1 |", [combo, f"{rho2:g}"], [rho2, h["c0"], e["c0"], h["c1"], e["c1"], h["power"], e["power"], h["P_lt_c1_H0"], e["P_lt_c1_H0"], h["inconclusive_H0"], e["inconclusive_H0"], h["inconclusive_H1"], e["inconclusive_H1"]]))
        tab.append(cmp(T3, "frente a v2", [combo, f"{rho2:g}"], [rho2, v2["v2_B_half"]["c0"], v2["v2_mc_err"]["c0"], v2["v2_B_half"]["c1"], v2["v2_mc_err"]["c1"], v2["v2_B_half"]["power"], v2["v2_mc_err"]["power"],
                                                              h["c0"], e["c0"], h["c1"], e["c1"], h["power"], e["power"], max(abs(x) for x in v2["deviation_in_v2_sigma"].values())]))
        for hyp in ("H0", "H1"):
            tab.append(cmp(T3, "d²_DR3 pct", [combo, f"{rho2:g}", hyp], [rho2] + h["d2_pct_5_50_95"][hyp] + [h[f"med_d2_{hyp}"], e[f"med_d2_{hyp}"], c["exact"][hyp]["E_d2"]]))
            tab.append(cmp(T3, "Λ pct", [combo, f"{rho2:g}", hyp], [rho2] + h["Lambda_pct_5_50_95"][hyp] + [h[f"P_gt3_{hyp}"], e[f"P_gt3_{hyp}"], h[f"P_ltm3_{hyp}"], e[f"P_ltm3_{hyp}"]]))
for name in V3.BIN_NAMES:
    b = J3S["combos"]["sinSN"]["cases"]["rho2_0.6"]["nesting"]["bins"][name]
    tab.append(cmp(T3, "| bin | autovalores C2", [name], b["C2"]["eigenvalues"] + b["C3"]["eigenvalues"] + b["C2_minus_C3"]["eigenvalues"] + b["Cov_e"]["eigenvalues"] + [b["rho2_eff_trace"]]))
for rho2 in V3.RHO2_GRID:
    a, b = J3S["combos"]["sinSN"]["rho2_sweep"][f"{rho2:.2f}"], J3S["combos"]["Pantheon+"]["rho2_sweep"][f"{rho2:.2f}"]
    tab.append(cmp(T3, "ρ² común", [f"{rho2:.2f}"], [rho2, a["c0"], a["mc_err"]["c0"], a["c1"], a["mc_err"]["c1"], a["power"], a["mc_err"]["power"], a["P_lt_c1_H0"], a["inconclusive_H0"], a["inconclusive_H1"], a["med_d2_H0"], a["med_d2_H1"],
                                                 b["c0"], b["c1"], b["power"], b["mc_err"]["power"], b["med_d2_H0"], b["med_d2_H1"]]))
for combo in ("sinSN", "Pantheon+"):
    for rho2 in P.RHO2S:
        v = J4S["combos"][combo]["cases"][f"rho2_{rho2:g}"]["vs_v3"]
        tab.append(cmp(T4, "máx \\|Δa\\| (unidades absolutas)", [combo, f"{rho2:g}"], [rho2, v["max_abs_diff_a_absolute_units"], v["max_abs_diff_Sigma_absolute_units"], v["diff_hp"]["c0"], v["diff_hp"]["c1"], v["diff_hp"]["power"], v["diff_hp"]["med_d2_H0"], v["diff_hp"]["med_d2_H1"]]))
for dz in V4.SHIFTS:
    q = J4S["P24"][f"dz_{dz:+g}"]
    tab.append(cmp(T4, "Δz_eff de LRG2", [f"{dz:+g} ({0.706 + dz:.3f})"], [dz, q["z_LRG2"], q["c0"], q["mc_err"]["c0"], q["c1"], q["mc_err"]["c1"], q["power"], q["mc_err"]["power"], q["d_c0"], q["d_c1"], q["d_power"]]))
q0 = J4S["P24"]["dz_0"]
tab.append(cmp(T4, "Δz_eff de LRG2", ["0 (0.706)"], [0.706, q0["c0"], q0["mc_err"]["c0"], q0["c1"], q0["mc_err"]["c1"], q0["power"], q0["mc_err"]["power"]]))
bad = [t for t in tab if (not t.get("found")) or t.get("mismatches") or t.get("count_mismatch")]
text4 = {"m2_LRG2_p2_(17.5603,20.0907)": J4S["combos"]["sinSN"]["cases"]["rho2_0.6"]["m2_bao3_LRG2"], "falso_z0712_c0_c1_power_(1.0260,-1.5918,0.6582)": [calf["combos"]["sinSN"][k] for k in ("c0", "c1", "power")],
         "falso_z0712_r3_d2_Lambda": [resf["combos"]["sinSN"]["r3_fractional"], resf["combos"]["sinSN"]["d2_DR3"], resf["combos"]["sinSN"]["Lambda"]], "v3_falso_c0_c1_(1.0260,-1.5921)": [json.load(open(HERE / "pronostico_dr3_lrg2_v3_calibracion_falso.json"))["combos"]["sinSN"][k] for k in ("c0", "c1")],
         "P24_json": J4S["predictions"]["P24"], "v3_predictions": J3S["predictions"]}
out["6_tablas"] = {"n_rows": len(tab), "n_bad": len(bad), "bad": bad, "text_checks": text4}
print(f"6 tablas: {len(tab)} filas, {len(bad)} con problemas: {json.dumps(bad)[:1500]}")
print("text checks:", json.dumps(text4)[:800])
out["runtime_s"] = time.time() - t0
(HERE / "revision3_pronostico_comprobacion.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")
print("escrito revision3_pronostico_comprobacion.json", f"{out['runtime_s']:.0f} s")
