#!/usr/bin/env python3
"""
REVISION independiente de pronostico_dr3_lrg2.py (no modifica nada existente; escribe revision_pronostico_reproduccion.json).

1. Reproduce el caso rho2 = 0.6, sinSN (y Pantheon+ como extra) con las semillas del json, reimplementando la generacion
   linealizada a partir de las ecuaciones del borrador (NO llamando a la funcion `simulate` del original, que esta dentro de main),
   y compara r3_mean_A, r3_cov_A, percentiles de d2 (B), c0, c1, potencia con los del json.
2. Calcula la media y la covarianza EXACTAS de r3 en el modelo linealizado (r3 es afin en (z_p, z_e, z_mu)) y las compara con
   las simuladas (A) y con las "formulas de referencia" del borrador. Esto da una referencia analitica que si deberia coincidir
   con las simulaciones (al error Monte Carlo), cosa que la referencia del borrador no hace.
3. Con las densidades exactas (no estimadas) calcula c0, c1 y potencia en B y, con 20 reejecuciones con otras semillas, el error
   Monte Carlo de c0, c1, potencia y medianas.

    ../venv/bin/python -I revision_pronostico_reproduccion.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.stats import chi2 as chi2_dist, multivariate_normal

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import pronostico_dr3_lrg2 as P  # noqa: E402  (fit_c, fisher, constants; its main() is not run)

PB, JB = P.PB, P.JB
D, COV13, LRG2, KEEP, C22 = P.D, P.COV13, P.LRG2, P.KEEP, P.C22
SEED_BASE, N_UNI, N_VAL = P.SEED_BASE, P.N_UNI, P.N_VAL
ORIG = json.load(open(HERE / "pronostico_dr3_lrg2.json"))


def setup(combo):
    """DR2 constants exactly as main() of the original computes them (same calls, same order)."""
    npar = P.npar_of(combo)
    chi2_full, p_full = PB.fit_d(combo, list(range(13)), D)
    Jf, Ff, Sf, _, _ = P.fisher(combo, list(range(13)), p_full, COV13)
    chi2_2, p2 = PB.fit_d(combo, KEEP, D)
    J, F2, S2, CJ2, rows2 = P.fisher(combo, KEEP, p2, COV13)
    m2 = JB.model_at(p2)
    r2 = (D - m2)[LRG2]
    J2 = J[LRG2]
    Q = J2 @ S2 @ J2.T
    V2 = C22 + Q
    return dict(npar=npar, p_full=p_full, Sf=Sf, p2=p2, S2=S2, J=J, J2=J2, m2=m2, r2=r2, Q=Q, V2=V2, rows2=rows2,
                m2_bao=m2[:13], m2_rest=m2[rows2[len(KEEP):]], d_rest=D[rows2[len(KEEP):]],
                L_full=np.linalg.cholesky(Sf), L2=np.linalg.cholesky(S2), L22=np.linalg.cholesky(C22), L13=np.linalg.cholesky(COV13))


def case_ops(combo, st, rho2):
    cov3 = rho2 * COV13
    J3f, F3, S3, CJ3, _ = P.fisher(combo, KEEP, st["p2"], cov3)
    G3 = S3 @ CJ3.T
    G3_bao, G3_rest = G3[:, :len(KEEP)], G3[:, len(KEEP):]
    const_rest = G3_rest @ (st["d_rest"] - st["m2_rest"])
    V3k = rho2 * C22 + st["J2"] @ S3 @ st["J2"].T
    return dict(cov3=cov3, S3=S3, G3_bao=G3_bao, const_rest=const_rest, V3k=V3k, V3ki=np.linalg.inv(V3k))


def affine_map(hyp, st, op, rho2):
    """r3 = a + M z with z = (z_p, z_e[, z_mu]) standard normal: exact mean a and covariance M M^T of the linearised r3.
    Written from the equations of the draft (not from the original code)."""
    npar, J, J2 = st["npar"], st["J"], st["J2"]
    Jb = J[:13]
    Lp = st["L_full"] if hyp == "H0" else st["L2"]
    p_mean = st["p_full"][:npar] if hyp == "H0" else st["p2"][:npar]
    nz = npar + 13 + (2 if hyp == "H1" else 0)
    # p_true = p_mean + Lp z_p ; dp = p_true - p2 ; m_true = m2 + Jb dp
    # H0: d3 = m_true + rho2 (D - m_true) + e = (1 - rho2) m_true + rho2 D + e
    # H1: mu = D_L - m_true_L + L22 z_mu ; n2 = D - m_true - E2 mu ; d3 = m_true + E2 mu + rho2 n2 + e
    #     -> d3 = (1 - rho2)(m_true + E2 mu) + rho2 D + e
    # refit: dp3 = G3_bao (d3_keep - m2_keep) + const_rest ; r3 = d3_L - m2_L - J2 dp3
    def d3_of(z):
        z_p = z[:npar]
        z_e = z[npar:npar + 13]
        p_true = p_mean + Lp @ z_p
        m_true = st["m2_bao"] + Jb @ (p_true - st["p2"][:npar])
        e = np.sqrt(rho2 * (1 - rho2)) * (st["L13"] @ z_e)
        if hyp == "H1":
            mu = D[LRG2] - m_true[LRG2] + st["L22"] @ z[npar + 13:]
            x = m_true.copy()
            x[LRG2] += mu
            return (1 - rho2) * x + rho2 * D[:13] + e
        return (1 - rho2) * m_true + rho2 * D[:13] + e

    def r3_of(z):
        d3 = d3_of(z)
        dp3 = op["G3_bao"] @ (d3[KEEP] - st["m2_bao"][KEEP]) + op["const_rest"]
        return d3[LRG2] - st["m2_bao"][LRG2] - J2 @ dp3

    a = r3_of(np.zeros(nz))
    M = np.column_stack([r3_of(np.eye(nz)[k]) - a for k in range(nz)])
    # check linearity: r3(2 e_k) - a = 2 (r3(e_k) - a)
    lin_err = max(float(np.max(np.abs(r3_of(2 * np.eye(nz)[k]) - a - 2 * M[:, k]))) for k in range(nz))
    return a, M @ M.T, lin_err


def simulate(hyp, st, op, rho2, n, rng):
    """Re-implementation of the linearised generation + refit from the draft's equations (same draw order as the original:
    z_p (n x npar), z_e (n x 13), z_mu (n x 2) if H1), so the same seed must give the same universes."""
    npar = st["npar"]
    z_p = rng.standard_normal((n, npar))
    z_e = rng.standard_normal((n, 13))
    z_mu = rng.standard_normal((n, 2)) if hyp == "H1" else None
    p_true = (st["p_full"][:npar] + z_p @ st["L_full"].T) if hyp == "H0" else (st["p2"][:npar] + z_p @ st["L2"].T)
    m_true = st["m2_bao"][None, :] + (p_true - st["p2"][:npar]) @ st["J"][:13].T
    e = np.sqrt(rho2 * (1 - rho2)) * (z_e @ st["L13"].T)
    x = m_true.copy()
    if hyp == "H1":
        mu = D[LRG2][None, :] - m_true[:, LRG2] + z_mu @ st["L22"].T
        x[:, LRG2] += mu
    d3 = (1 - rho2) * x + rho2 * D[:13][None, :] + e          # algebraically identical to m_true + E2 mu + rho2 n2 + e
    dp3 = (d3[:, KEEP] - st["m2_bao"][KEEP][None, :]) @ op["G3_bao"].T + op["const_rest"][None, :]
    r3 = d3[:, LRG2] - st["m2_bao"][LRG2][None, :] - dp3 @ st["J2"].T
    d2 = np.einsum("ni,ij,nj->n", r3, op["V3ki"], r3)
    return r3, d2


def decision(lam0, lam1):
    c0, c1 = float(np.percentile(lam0, 95)), float(np.percentile(lam1, 5))
    return dict(c0=c0, c1=c1, c1_lt_c0=bool(c1 < c0), power=float(np.mean(lam1 > c0)), P_lt_c1_H0=float(np.mean(lam0 < c1)))


def main():
    t0 = time.time()
    out = {"note": "revision independiente; nada existente modificado", "cases": {}}
    for ic, combo in enumerate(("sinSN", "Pantheon+")):
        st = setup(combo)
        o2 = ORIG["combos"][combo]["DR2"]
        out["cases"][combo] = {"DR2_check": {"r2_max_abs_diff": float(np.max(np.abs(st["r2"] - np.array(o2["r2"])))),
                                             "p2_max_abs_diff": float(np.max(np.abs(st["p2"] - np.array(o2["p2_without_LRG2"])))),
                                             "V2_max_abs_diff": float(np.max(np.abs(st["V2"] - np.array(o2["V2"]))))}, "rho2": {}}
        for ir, rho2 in enumerate(P.RHO2S):
            if combo == "Pantheon+" and rho2 != 0.6:
                continue
            oc = ORIG["combos"][combo]["cases"][f"rho2_{rho2:g}"]
            op = case_ops(combo, st, rho2)
            rec = {"V3k_max_abs_diff_vs_json": float(np.max(np.abs(op["V3k"] - np.array(oc["V3k_linearised"])))), "hyp": {}}
            halves = {}
            exact = {}
            for ih, hyp in enumerate(("H0", "H1")):
                oh = oc["hyp"][hyp]
                H = {}
                for half, code in (("A", 0), ("B", 1)):
                    seed = [SEED_BASE, ic, ir, ih, code]
                    assert seed == oh["seeds"][half], (seed, oh["seeds"][half])
                    halves[(hyp, half)] = simulate(hyp, st, op, rho2, N_UNI // 2, np.random.default_rng(seed))
                rA, dA = halves[(hyp, "A")]
                rB, dB = halves[(hyp, "B")]
                H["r3_mean_A_max_abs_diff_vs_json"] = float(np.max(np.abs(rA.mean(axis=0) - np.array(oh["r3_mean_A"]))))
                H["r3_cov_A_max_abs_diff_vs_json"] = float(np.max(np.abs(np.cov(rA.T) - np.array(oh["r3_cov_A"]))))
                pB = np.percentile(dB, [5, 50, 95])
                H["d2_pct_B"] = pB.tolist()
                H["d2_pct_B_max_abs_diff_vs_json"] = float(np.max(np.abs(pB - np.array(oh["d2_percentiles_5_50_95_B"]))))
                # exact linear-Gaussian mean/cov
                a, Cx, lin_err = affine_map(hyp, st, op, rho2)
                exact[hyp] = (a, Cx)
                sd = np.sqrt(np.diag(Cx))
                rf = oc["reference"]["H0"] if hyp == "H0" else oc["reference"]["H1_marginalised_flat_mu"]
                H["exact_linear"] = {"mean": a.tolist(), "cov": Cx.tolist(), "linearity_check_max_abs": lin_err,
                                     "sim_A_mean_minus_exact_over_sigma": ((rA.mean(axis=0) - a) / sd).tolist(),
                                     "mc_sigma_of_mean_over_sigma": float(1 / np.sqrt(N_UNI // 2)),
                                     "sim_A_cov_diag_over_exact_minus_1": (np.diag(np.cov(rA.T)) / np.diag(Cx) - 1).tolist(),
                                     "mc_rel_sigma_of_var": float(np.sqrt(2 / (N_UNI // 2))),
                                     "draft_reference_mean_minus_exact_over_sigma": ((np.array(rf["mean"]) - a) / sd).tolist(),
                                     "draft_reference_cov_diag_over_exact_minus_1": (np.diag(np.array(rf["cov"])) / np.diag(Cx) - 1).tolist(),
                                     "exact_d2_of_mean": float(a @ op["V3ki"] @ a), "exact_E_d2": float(a @ op["V3ki"] @ a + np.trace(op["V3ki"] @ Cx))}
                H["d2_mean_B"] = float(dB.mean())
                rec["hyp"][hyp] = H
            # Lambda as in the original: Gaussian densities from A, evaluated in B
            dens = {h: multivariate_normal(mean=halves[(h, "A")][0].mean(axis=0), cov=np.cov(halves[(h, "A")][0].T)) for h in ("H0", "H1")}
            lam = {h: dens["H1"].logpdf(halves[(h, "B")][0]) - dens["H0"].logpdf(halves[(h, "B")][0]) for h in ("H0", "H1")}
            dec = decision(lam["H0"], lam["H1"])
            od = oc["decision"]
            dec["diff_vs_json"] = {"c0": dec["c0"] - od["c0_pct95_H0"], "c1": dec["c1"] - od["c1_pct5_H1"], "power": dec["power"] - od["power_P_Lambda_gt_c0_given_H1"]}
            dec["Lambda_pct_5_50_95"] = {h: np.percentile(lam[h], [5, 50, 95]).tolist() for h in lam}
            rec["decision_repro"] = dec
            # Lambda with the EXACT densities (no estimation), evaluated in B
            dens_x = {h: multivariate_normal(mean=exact[h][0], cov=exact[h][1]) for h in exact}
            lam_x = {h: dens_x["H1"].logpdf(halves[(h, "B")][0]) - dens_x["H0"].logpdf(halves[(h, "B")][0]) for h in ("H0", "H1")}
            rec["decision_exact_densities_B"] = decision(lam_x["H0"], lam_x["H1"])
            # Monte Carlo error: 20 re-runs with other seeds, full pipeline (A estimates, B evaluates), N_UNI/2 each
            reps = []
            for k in range(20):
                hv = {}
                for ih, hyp in enumerate(("H0", "H1")):
                    for half, code in (("A", 10 + 2 * k), ("B", 11 + 2 * k)):
                        hv[(hyp, half)] = simulate(hyp, st, op, rho2, N_UNI // 2, np.random.default_rng([SEED_BASE, ic, ir, ih, code, 999]))
                dn = {h: multivariate_normal(mean=hv[(h, "A")][0].mean(axis=0), cov=np.cov(hv[(h, "A")][0].T)) for h in ("H0", "H1")}
                lm = {h: dn["H1"].logpdf(hv[(h, "B")][0]) - dn["H0"].logpdf(hv[(h, "B")][0]) for h in ("H0", "H1")}
                d = decision(lm["H0"], lm["H1"])
                d["med_d2_H0"], d["med_d2_H1"] = float(np.median(hv[("H0", "B")][1])), float(np.median(hv[("H1", "B")][1]))
                reps.append(d)
            rec["mc_error_20_reruns"] = {k: {"mean": float(np.mean([r[k] for r in reps])), "std": float(np.std([r[k] for r in reps], ddof=1))}
                                         for k in ("c0", "c1", "power", "med_d2_H0", "med_d2_H1")}
            out["cases"][combo]["rho2"][f"rho2_{rho2:g}"] = rec
            print(f"{combo} rho2={rho2}: repro diffs c0={dec['diff_vs_json']['c0']:.2e} c1={dec['diff_vs_json']['c1']:.2e} power={dec['diff_vs_json']['power']:.2e}; "
                  f"meanA diff {rec['hyp']['H0']['r3_mean_A_max_abs_diff_vs_json']:.1e}/{rec['hyp']['H1']['r3_mean_A_max_abs_diff_vs_json']:.1e}; "
                  f"d2pct diff {rec['hyp']['H0']['d2_pct_B_max_abs_diff_vs_json']:.1e}/{rec['hyp']['H1']['d2_pct_B_max_abs_diff_vs_json']:.1e}", flush=True)
            for hyp in ("H0", "H1"):
                x = rec["hyp"][hyp]["exact_linear"]
                print(f"   {hyp}: sim-exact mean/sigma {np.round(x['sim_A_mean_minus_exact_over_sigma'], 4).tolist()} (MC 1 sigma {x['mc_sigma_of_mean_over_sigma']:.4f}); "
                      f"var ratio-1 {np.round(x['sim_A_cov_diag_over_exact_minus_1'], 4).tolist()} (MC {x['mc_rel_sigma_of_var']:.4f}); "
                      f"draft ref - exact mean/sigma {np.round(x['draft_reference_mean_minus_exact_over_sigma'], 3).tolist()}, var ratio-1 {np.round(x['draft_reference_cov_diag_over_exact_minus_1'], 3).tolist()}; linearity {x['linearity_check_max_abs']:.1e}")
            print(f"   exact-density decision: {rec['decision_exact_densities_B']}")
            print(f"   MC error (20 reruns): { {k: (round(v['mean'], 4), round(v['std'], 4)) for k, v in rec['mc_error_20_reruns'].items()} }")
    out["runtime_s"] = time.time() - t0
    (HERE / "revision_pronostico_reproduccion.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")
    print("escrito revision_pronostico_reproduccion.json", f"{out['runtime_s']:.0f} s")


if __name__ == "__main__":
    main()
