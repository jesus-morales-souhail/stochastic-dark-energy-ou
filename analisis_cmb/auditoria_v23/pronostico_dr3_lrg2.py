#!/usr/bin/env python3
"""
PRONOSTICO_DR3_BORRADOR.md (draft, not registered): forecast of the LRG2 residual in DESI DR3 under two hypotheses, using that DR3
contains DR2 (nested data). Set-up and fitters of predictivo_bins.py / jackknife_bins.py (imported, not modified).

ASSUMPTIONS of the forecast (not guaranteed properties of the future analysis):
- Nesting: for the 7 BAO bins, n3 = rho2 n2 + e, e independent of n2, Cov(e) = rho2 (1 - rho2) C, so Cov(n3) = rho2 C and
  Cov(n3, n2) = rho2 C; C = DR2 block covariance. External data (theta_*, omega_m, Pantheon+) unchanged. rho2 = 0.6 (main),
  0.5 and 0.72 (sensitivity), the same in every bin.
- H0: LCDM, no own deviation in LRG2. H1: LCDM + a fixed deviation mu (2 components) only in LRG2, flat unbounded prior; mu is
  physical and does NOT shrink between DR2 and DR3; only the noise is transmitted with the factor rho2.
- Residuals r = data_LRG2 - LCDM prediction fitted without LRG2, same procedure as PREDICTIVO_BINS (optimiser, priors, Fisher
  F = J^T C^-1 J); in DR3 the refit uses the other six DR3 bins (which also change) with covariance rho2 C for BAO.
Generation (linearised around the DR2 fit without LRG2, p2; validated case by case against the full computation):
  H0: p_true ~ N(p_full, Sigma_full) (Gaussian approximation of the LCDM posterior with all DR2 data); n2 = d2 - m(p_true) on the
      13 BAO rows; d3 = m(p_true) + rho2 n2 + e.
  H1: (p_true, mu) from the joint posterior with flat unbounded prior on mu: p_true ~ N(p2, Sigma2), mu | p_true ~ N(d2_LRG2 -
      m_LRG2(p_true), C_LRG2); n2 = d2 - m(p_true) - E2 mu; d3 = m(p_true) + E2 mu + rho2 n2 + e (mu fixed, noise transmitted).
  Refit without LRG2 on d3 with C3 (rho2 C on BAO, original on externals): p3; r3 = d3_LRG2 - m_LRG2(p3); Sigma3 = Fisher^-1 of
  that refit; V3 = rho2 C_LRG2 + J3 Sigma3 J3^T; d2_DR3 = r3^T V3^-1 r3.
100000 universes per (hypothesis, combination, rho2), halves A (densities) and B (Gaussianity, thresholds, power), seeds recorded.
Lambda(r3) = ln p(r3 | DR2, H1) - ln p(r3 | DR2, H0) with the A densities; c0 = 95th percentile of Lambda under H0 (B), c1 = 5th
under H1 (B); rule only if c1 < c0.

    ../venv/bin/python -I pronostico_dr3_lrg2.py        # writes pronostico_dr3_lrg2.json
Nothing existing is modified.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.linalg import solve_triangular
from scipy.optimize import least_squares
from scipy.stats import chi2 as chi2_dist, gaussian_kde, multivariate_normal

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import predictivo_bins as PB  # noqa: E402  (fitters with data argument, model, set-up; its main() is not run)

JB = PB.JB
S, X, G = PB.S, PB.X, PB.G
SN, NSN, NROWS = PB.SN, PB.NSN, PB.NROWS
D, CFULL, COV13 = PB.D, PB.CFULL, PB.COV13
LRG2 = [3, 4]
KEEP = [i for i in range(13) if i not in LRG2]
RHO2S = (0.6, 0.5, 0.72)
N_UNI, N_VAL = 100000, 20
SEED_BASE = 20261114
COMBOS = ("sinSN", "Pantheon+")
C22 = COV13[np.ix_(LRG2, LRG2)]


# ----------------------------------------------------------------------------------------------------------------------
# fitters and Fisher with the BAO covariance as an argument (same optimiser settings as predictivo_bins / jackknife_bins)
# ----------------------------------------------------------------------------------------------------------------------
def fit_c(combo, keep, d, cov13):
    L = np.linalg.cholesky(cov13[np.ix_(keep, keep)])
    kb = np.array(keep)
    y15 = d[:15]
    if combo == "sinSN":
        def resid(p):
            v, m = X.model(p[0], p[1], np.zeros(X.N))
            if v is None:
                return np.full(len(kb) + 2, 1e3 * (1.0 - m))
            r = y15 - v
            return np.concatenate([np.linalg.solve(L, r[kb]), [r[13] / X.SIG_TH, r[14] / X.SIG_WM]])
        r = least_squares(resid, np.asarray((0.682, 0.1407)), bounds=X.BOUNDS_P, method="trf", x_scale=[0.01, 0.001], ftol=1e-10, xtol=1e-10, max_nfev=300)
        return float(2.0 * r.cost), np.array([r.x[0], r.x[1], 0.0])
    mb = d[15:]
    x0 = (0.682, 0.1407)
    hist, _ = S.history(G, x0[1] / x0[0] ** 2, S.OMR / x0[0] ** 2, np.zeros(G.N))
    M0 = float((SN["Ci1"] @ (mb - S.mu_sn(G, SN, hist[0]))) / SN["e"]) if hist is not None else 0.0

    def resid(p):
        h, wm, M = p
        hist, m = S.history(G, wm / h**2, S.OMR / h**2, np.zeros(G.N))
        if hist is None:
            return np.full(len(kb) + 2 + NSN, 1e3 * (1.0 - m))
        dm_u, inv = hist
        r = y15 - S.vec15(G, h, wm, dm_u, inv)
        rs = mb - S.mu_sn(G, SN, dm_u) - M
        return np.concatenate([solve_triangular(L, r[kb], lower=True, check_finite=False), [r[13] / X.SIG_TH, r[14] / X.SIG_WM],
                               solve_triangular(SN["L"], rs, lower=True, check_finite=False)])
    r = least_squares(resid, np.array([x0[0], x0[1], M0]), bounds=([0.4, 0.05, -50.0], [1.0, 0.30, 50.0]), method="trf",
                      x_scale=[0.01, 0.001, 0.01], ftol=1e-10, xtol=1e-10, max_nfev=200)
    return float(2.0 * r.cost), np.array(r.x)


def npar_of(combo):
    return 2 if combo == "sinSN" else 3


def rows_of(combo, keep):
    return list(keep) + [13, 14] + (list(range(15, NROWS)) if combo == "Pantheon+" else [])


def cinv_apply(combo, keep, cov13, V):
    """C^-1 V for rows ordered [keep, 13, 14, SN...] with BAO covariance cov13 (block-diagonal C)."""
    nb = len(keep)
    out = [np.linalg.solve(cov13[np.ix_(keep, keep)], V[:nb]), np.array([1 / X.SIG_TH**2, 1 / X.SIG_WM**2])[:, None] * V[nb:nb + 2]]
    if combo == "Pantheon+":
        out.append(SN["Ci"] @ V[nb + 2:])
    return np.concatenate(out)


def fisher(combo, keep, p, cov13):
    """J (NROWS x npar, finite differences of the exact model at p), F = J^T C^-1 J on rows [keep, 13, 14, SN], Sigma = F^-1."""
    J = JB.build_B_at(p, npar_of(combo))
    rows = np.array(rows_of(combo, keep))
    CJ = cinv_apply(combo, keep, cov13, J[rows])
    F = J[rows].T @ CJ
    return J, F, np.linalg.inv(F), CJ, rows


def eig(M):
    return np.linalg.eigvalsh(0.5 * (M + M.T)).tolist()


# ----------------------------------------------------------------------------------------------------------------------
def main():
    t0 = time.time()
    out = {"assumptions": "see module docstring", "rho2": RHO2S, "N_universes_per_case": N_UNI, "halves": "A = first N/2 (densities), B = second N/2 (Gaussianity, thresholds, power), different seeds",
           "seed_base": SEED_BASE, "combos": {}}
    for ic, combo in enumerate(COMBOS):
        t1 = time.time()
        npar = npar_of(combo)
        rec = {}
        # ---- DR2: full fit, fit without LRG2, Fisher, residual
        chi2_full, p_full = PB.fit_d(combo, list(range(13)), D)
        chk_fit = fit_c(combo, list(range(13)), D, COV13)
        rec["check_fit_c_vs_predictivo_fit_d"] = {"p0_max_abs_diff": float(np.max(np.abs(chk_fit[1] - p_full))), "chi2_diff": chk_fit[0] - chi2_full}
        Jf, Ff, Sf, _, _ = fisher(combo, list(range(13)), p_full, COV13)
        chi2_2, p2 = PB.fit_d(combo, KEEP, D)
        J, F2, S2, CJ2, rows2 = fisher(combo, KEEP, p2, COV13)
        m2 = JB.model_at(p2)
        r2 = (D - m2)[LRG2]
        J2 = J[LRG2]
        Q = J2 @ S2 @ J2.T
        V2 = C22 + Q
        d2_dr2 = float(r2 @ np.linalg.solve(V2, r2))
        pbj = json.load(open(HERE / "predictivo_bins.json"))["combos"][combo]["bins"]["LRG2"]
        rec["DR2"] = {"p_full": p_full.tolist(), "chi2_full": chi2_full, "Sigma_full": Sf.tolist(), "p2_without_LRG2": p2.tolist(), "chi2_without_LRG2": chi2_2,
                      "Sigma2": S2.tolist(), "r2": r2.tolist(), "Q": Q.tolist(), "V2": V2.tolist(), "d2_DR2": d2_dr2, "p_local_DR2": float(chi2_dist.sf(d2_dr2, 2)),
                      "check_vs_predictivo_bins": {"d2_diff": d2_dr2 - pbj["d2"], "r_max_abs_diff": float(np.max(np.abs(r2 - np.array(pbj["r"])))), "p0_max_abs_diff": float(np.max(np.abs(p2 - np.array(pbj["p0_minus_k"]))))}}
        # ---- H1 posterior checks (joint Fisher of (p, mu) at the joint MLE (p2, mu_hat = r2), all rows of the combination)
        rows_all = np.array(rows_of(combo, list(range(13))))
        J_all = J[rows_all]
        CJ_all = cinv_apply(combo, list(range(13)), COV13, J_all)
        Fpp = J_all.T @ CJ_all
        C22i = np.linalg.inv(C22)
        Fpm = J2.T @ C22i
        Fj = np.block([[Fpp, Fpm], [Fpm.T, C22i]])
        Sj = np.linalg.inv(Fj)
        chk = {"1_C_LRG2_invertible": {"eigenvalues": eig(C22), "det": float(np.linalg.det(C22)), "ok": bool(np.all(np.linalg.eigvalsh(C22) > 0))},
               "2_mu_block_of_Fisher_invertible": {"block": "C_LRG2^-1", "eigenvalues": eig(C22i), "ok": bool(np.all(np.linalg.eigvalsh(C22i) > 0))},
               "3_joint_Fisher_pd_inverse_finite": {"eigenvalues_F_joint": eig(Fj), "ok": bool(np.all(np.linalg.eigvalsh(0.5 * (Fj + Fj.T)) > 0) and np.all(np.isfinite(Sj))),
                                                    "max_abs_F_Sinv_minus_I": float(np.max(np.abs(Fj @ Sj - np.eye(len(Fj)))))},
               "4_eigenvalues": {"F_joint": eig(Fj), "Sigma_joint_posterior": eig(Sj), "Sigma_joint_p_block_vs_Sigma2_max_abs_diff": float(np.max(np.abs(Sj[:npar, :npar] - S2))),
                                 "Sigma_joint_mu_block_vs_V2_max_abs_diff": float(np.max(np.abs(Sj[npar:, npar:] - V2)))},
               "5_posterior_of_mu_proper": {"requires": "posterior of p with the other six bins + externals proper: F2 = J^T C^-1 J on those rows positive definite",
                                            "eigenvalues_F2": eig(F2), "ok": bool(np.all(np.linalg.eigvalsh(F2) > 0)), "note": "marginalising mu with a flat prior, the LRG2 rows contribute a constant independent of p"}}
        chk["all_ok"] = bool(all(chk[k]["ok"] for k in ("1_C_LRG2_invertible", "2_mu_block_of_Fisher_invertible", "3_joint_Fisher_pd_inverse_finite", "5_posterior_of_mu_proper")))
        rec["H1_posterior_checks"] = chk
        print(f"{combo}: DR2 d2={d2_dr2:.3f} (predictivo {pbj['d2']:.3f}); H1 checks ok={chk['all_ok']}; eig F_joint={np.array(eig(Fj))}", flush=True)
        assert chk["all_ok"], "PARA: comprobación del posterior de H1 fallida"
        # constants of the linearised refit
        m2_bao, m2_rest = m2[:13], m2[rows2[len(KEEP):]]
        d_rest = D[rows2[len(KEEP):]]
        L_full = np.linalg.cholesky(Sf)
        L2 = np.linalg.cholesky(S2)
        L22 = np.linalg.cholesky(C22)
        L13 = np.linalg.cholesky(COV13)
        rec["cases"] = {}
        for ir, rho2 in enumerate(RHO2S):
            cov3 = rho2 * COV13
            J3f, F3, S3, CJ3, _ = fisher(combo, KEEP, p2, cov3)            # linearised: J at p2, C3 -> Sigma3 constant
            G3 = S3 @ CJ3.T                                                 # npar x n_rows(-2):  p3 - p2 = G3 (d3_(-2) - m_(-2)(p2))
            G3_bao, G3_rest = G3[:, :len(KEEP)], G3[:, len(KEEP):]
            const_rest = G3_rest @ (d_rest - m2_rest)                       # externals unchanged
            V3k = rho2 * C22 + J2 @ S3 @ J2.T
            V3ki = np.linalg.inv(V3k)
            # reference formulas (step 1), not used for Lambda
            V3 = rho2 * C22 + Q
            C32 = rho2 * C22 + Q
            K = C32 @ np.linalg.inv(V2)
            ref = {"H0": {"mean": (K @ r2).tolist(), "cov": (V3 - K @ V2 @ K.T).tolist(), "K": K.tolist(), "V3": V3.tolist()},
                   "H1_conditional_mu_known": {"mean_formula": "rho2 r2 + (1 - rho2) mu", "cov": (rho2 * (1 - rho2) * C22).tolist()},
                   "H1_marginalised_flat_mu": {"mean": r2.tolist(), "cov": ((1 - rho2) * C22).tolist(),
                                               "integration_check": {"mean_diff_max_abs": float(np.max(np.abs(rho2 * r2 + (1 - rho2) * r2 - r2))),
                                                                     "cov_diff_max_abs": float(np.max(np.abs(rho2 * (1 - rho2) * C22 + (1 - rho2) ** 2 * C22 - (1 - rho2) * C22)))}},
                   "H1_with_refit": "r3 - r2 ~ (rho2 - 1) n2 + e - J2 (p3 - p2)", "compared_with_simulations": "H1_marginalised_flat_mu"}
            case = {"rho2": rho2, "Sigma3_linearised": S3.tolist(), "V3k_linearised": V3k.tolist(), "reference": ref, "hyp": {}}

            def simulate(hyp, n, rng):
                """Linearised generation + refit; returns r3 (n x 2), d2 (n), and the draws (for the validation)."""
                z_p = rng.standard_normal((n, npar))
                z_e = rng.standard_normal((n, 13))
                z_mu = rng.standard_normal((n, 2)) if hyp == "H1" else None
                p_true = (p_full[:npar] + z_p @ L_full.T) if hyp == "H0" else (p2[:npar] + z_p @ L2.T)
                dp = p_true - p2[:npar]
                m_true_bao = m2_bao[None, :] + dp @ J[:13].T
                e = np.sqrt(rho2 * (1 - rho2)) * (z_e @ L13.T)
                n2 = D[:13][None, :] - m_true_bao
                mu = None
                if hyp == "H1":
                    mu = (D[LRG2][None, :] - m_true_bao[:, LRG2]) + z_mu @ L22.T
                    n2 = n2.copy()
                    n2[:, LRG2] -= mu
                d3 = m_true_bao + rho2 * n2 + e
                if hyp == "H1":
                    d3[:, LRG2] += mu
                dp3 = (d3[:, KEEP] - m2_bao[KEEP][None, :]) @ G3_bao.T + const_rest[None, :]
                r3 = d3[:, LRG2] - m2_bao[LRG2][None, :] - dp3 @ J2.T
                d2 = np.einsum("ni,ij,nj->n", r3, V3ki, r3)
                return r3, d2, {"z_p": z_p, "z_e": z_e, "z_mu": z_mu, "p_true": p_true, "mu": mu, "e": e}

            for ih, hyp in enumerate(("H0", "H1")):
                H = {"seeds": {}}
                halves = {}
                for half in ("A", "B"):
                    seed = [SEED_BASE, ic, ir, ih, 0 if half == "A" else 1]
                    H["seeds"][half] = seed
                    r3, d2, _ = simulate(hyp, N_UNI // 2, np.random.default_rng(seed))
                    halves[half] = (r3, d2)
                rA, dA = halves["A"]
                rB, dB = halves["B"]
                H["r3_mean_A"], H["r3_cov_A"] = rA.mean(axis=0).tolist(), np.cov(rA.T).tolist()
                H["d2_percentiles_5_50_95_B"] = np.percentile(dB, [5, 50, 95]).tolist()
                H["d2_mean_B"] = float(dB.mean())
                # Gaussianity in B with A's mean/cov
                mu_A, cov_A = rA.mean(axis=0), np.cov(rA.T)
                dm = np.einsum("ni,ij,nj->n", rB - mu_A, np.linalg.inv(cov_A), rB - mu_A)
                pcts = np.percentile(dm, [5, 50, 95])
                ref_q = chi2_dist.ppf([0.05, 0.5, 0.95], 2)
                frac95 = float(np.mean(dm > ref_q[2]))
                gauss_ok = bool(np.all(np.abs(pcts / ref_q - 1) < 0.10) and abs(frac95 - 0.05) < 0.005)
                H["gaussianity_B"] = {"mahalanobis_pct_5_50_95": pcts.tolist(), "chi2_2_quantiles": ref_q.tolist(), "rel_diff": (pcts / ref_q - 1).tolist(), "frac_above_q95": frac95, "gaussian_ok": gauss_ok}
                # comparison with the reference formula
                rf = ref["H0"] if hyp == "H0" else ref["H1_marginalised_flat_mu"]
                H["vs_reference"] = {"mean_diff": (mu_A - np.array(rf["mean"])).tolist(), "mean_diff_over_sigma": ((mu_A - np.array(rf["mean"])) / np.sqrt(np.diag(cov_A))).tolist(),
                                     "cov_diff": (cov_A - np.array(rf["cov"])).tolist(), "cov_rel_diff_diag": (np.diag(cov_A) / np.diag(np.array(rf["cov"])) - 1).tolist()}
                H["_A"], H["_B"] = (rA, dA), (rB, dB)
                case["hyp"][hyp] = H
            # densities (A) and Lambda (B)
            dens = {}
            for hyp in ("H0", "H1"):
                rA = case["hyp"][hyp]["_A"][0]
                if case["hyp"][hyp]["gaussianity_B"]["gaussian_ok"]:
                    dens[hyp] = ("gaussian", multivariate_normal(mean=rA.mean(axis=0), cov=np.cov(rA.T)))
                else:
                    dens[hyp] = ("kde", gaussian_kde(rA.T))
            def lam(r):
                return (dens["H1"][1].logpdf(r) if dens["H1"][0] == "gaussian" else dens["H1"][1].logpdf(r.T)) - (dens["H0"][1].logpdf(r) if dens["H0"][0] == "gaussian" else dens["H0"][1].logpdf(r.T))
            lamB = {hyp: lam(case["hyp"][hyp]["_B"][0]) for hyp in ("H0", "H1")}
            c0, c1 = float(np.percentile(lamB["H0"], 95)), float(np.percentile(lamB["H1"], 5))
            dec = {"density_type": {h: dens[h][0] for h in dens}, "c0_pct95_H0": c0, "c1_pct5_H1": c1, "c1_lt_c0": bool(c1 < c0),
                   "Lambda_percentiles_5_50_95": {h: np.percentile(lamB[h], [5, 50, 95]).tolist() for h in lamB},
                   "frac_abs_Lambda_gt_3": {h: float(np.mean(np.abs(lamB[h]) > 3)) for h in lamB}, "frac_Lambda_gt_3": {h: float(np.mean(lamB[h] > 3)) for h in lamB},
                   "frac_Lambda_lt_minus3": {h: float(np.mean(lamB[h] < -3)) for h in lamB}}
            if c1 < c0:
                dec.update({"power_P_Lambda_gt_c0_given_H1": float(np.mean(lamB["H1"] > c0)), "P_Lambda_lt_c1_given_H0": float(np.mean(lamB["H0"] < c1)),
                            "inconclusive_given_H0": float(np.mean((lamB["H0"] >= c1) & (lamB["H0"] <= c0))), "inconclusive_given_H1": float(np.mean((lamB["H1"] >= c1) & (lamB["H1"] <= c0)))})
            else:
                dec["note"] = "las regiones se solapan: el procedimiento no permite la regla de tres resultados con esos errores"
            case["decision"] = dec
            # validation of the linearised computation: N_VAL universes per hypothesis, full generation (exact model) + full refit
            case["validation"] = {}
            for ih, hyp in enumerate(("H0", "H1")):
                seed = [SEED_BASE, ic, ir, ih, 2]
                r3l, d2l, dr = simulate(hyp, N_VAL, np.random.default_rng(seed))
                diffs = []
                for u in range(N_VAL):
                    p_t = np.array([dr["p_true"][u][0], dr["p_true"][u][1], dr["p_true"][u][2] if npar == 3 else 0.0])
                    m_true = JB.model_at(p_t)
                    n2 = D[:13] - m_true[:13]
                    d3 = D.copy()
                    if hyp == "H1":
                        n2 = n2.copy()
                        n2[LRG2] -= dr["mu"][u]
                    d3[:13] = m_true[:13] + rho2 * n2 + dr["e"][u]
                    if hyp == "H1":
                        d3[LRG2] += dr["mu"][u]
                    _, p3 = fit_c(combo, KEEP, d3, cov3)
                    J3, _, S3f, _, _ = fisher(combo, KEEP, p3, cov3)
                    r3 = (d3 - JB.model_at(p3))[LRG2]
                    V = rho2 * C22 + J3[LRG2] @ S3f @ J3[LRG2].T
                    d2f = float(r3 @ np.linalg.solve(V, r3))
                    diffs.append({"universe": u, "d2_full": d2f, "d2_lin": float(d2l[u]), "abs_diff": abs(d2f - d2l[u]), "r3_max_abs_diff": float(np.max(np.abs(r3 - r3l[u])))})
                worst = max(diffs, key=lambda x: x["abs_diff"])
                case["validation"][hyp] = {"seed": seed, "n": N_VAL, "max_abs_diff_d2": worst["abs_diff"], "at_universe": worst["universe"], "max_r3_abs_diff": max(x["r3_max_abs_diff"] for x in diffs),
                                           "passes_0.1": bool(worst["abs_diff"] <= 0.1), "all": diffs}
                print(f"  {combo} rho2={rho2} {hyp}: d2 pct(5/50/95)={np.round(case['hyp'][hyp]['d2_percentiles_5_50_95_B'], 2).tolist()} gauss_ok={case['hyp'][hyp]['gaussianity_B']['gaussian_ok']} "
                      f"validation max|dd2|={worst['abs_diff']:.3e} (u {worst['universe']})", flush=True)
                assert worst["abs_diff"] <= 0.1, "PARA: linealizado frente a completo difiere en más de 0.1"
            print(f"  {combo} rho2={rho2}: c0={c0:.3f} c1={c1:.3f} c1<c0={c1 < c0} power={dec.get('power_P_Lambda_gt_c0_given_H1')} dens={dec['density_type']}", flush=True)
            for hyp in ("H0", "H1"):
                case["hyp"][hyp].pop("_A"), case["hyp"][hyp].pop("_B")
            rec["cases"][f"rho2_{rho2:g}"] = case
        rec["seconds"] = time.time() - t1
        out["combos"][combo] = rec
    main_case = out["combos"]["sinSN"]["cases"]["rho2_0.6"]
    pw = main_case["decision"].get("power_P_Lambda_gt_c0_given_H1")
    med0, med1 = main_case["hyp"]["H0"]["d2_percentiles_5_50_95_B"][1], main_case["hyp"]["H1"]["d2_percentiles_5_50_95_B"][1]
    out["predictions"] = {"P18": {"statement": "rho2 = 0.6, falsa alarma 5 %: potencia para distinguir H0 de H1 en LRG2 entre 0.5 y 0.75 (evaluada en sin SN)", "power_sinSN": pw,
                                  "power_Pantheon+": out["combos"]["Pantheon+"]["cases"]["rho2_0.6"]["decision"].get("power_P_Lambda_gt_c0_given_H1"), "holds": bool(pw is not None and 0.5 <= pw <= 0.75)},
                          "P19": {"statement": "rho2 = 0.6: mediana de d2_DR3 de LRG2 entre 3 y 6 bajo H0 y entre 9 y 15 bajo H1 (evaluada en sin SN)", "median_H0_sinSN": med0, "median_H1_sinSN": med1,
                                  "holds": bool(3 <= med0 <= 6 and 9 <= med1 <= 15)}}
    out["runtime_s"] = time.time() - t0
    (HERE / "pronostico_dr3_lrg2.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")
    print(json.dumps(out["predictions"], indent=1))
    print("escrito pronostico_dr3_lrg2.json", f"{out['runtime_s']:.0f} s")


if __name__ == "__main__":
    main()
