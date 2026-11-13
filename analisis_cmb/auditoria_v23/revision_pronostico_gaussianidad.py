#!/usr/bin/env python3
"""
REVISION: gaussianidad bidimensional de r3 (mitad B, blanqueada con media y covarianza de la mitad A), rho2 = 0.6, H0 y H1,
sinSN y Pantheon+. Tests: Rayleigh y Kuiper (angulo uniforme), KS del radio frente a chi_2 (Rayleigh(1)), radio por sectores de
angulo (Kruskal-Wallis + KS por sector + chi2 de independencia en malla 8 x 5), asimetria y curtosis por componente, correlacion.
Si algun test relevante falla (p < 0.01 en Rayleigh/Kuiper/KS radio/Kruskal/independencia, o |asimetria| > 0.05 o |exceso curtosis| > 0.1),
recalcula c0, c1 y potencia con KDE gaussiano cuyo ancho de banda se elige SOLO en A (validacion cruzada por verosimilitud dentro de A);
B solo calibra umbrales y evalua potencia. Escribe revision_pronostico_gaussianidad.json. Nada existente se modifica.

    ../venv/bin/python -I revision_pronostico_gaussianidad.py [--kde]   (--kde fuerza el calculo KDE aunque no falle nada)
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.stats import chi2 as chi2_dist, kstest, kruskal, skew, kurtosis, chi2_contingency, multivariate_normal

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import pronostico_dr3_lrg2 as P  # noqa: E402
import revision_pronostico_reproduccion as R  # noqa: E402  (setup, case_ops, simulate, decision; its main() is not run)

FORCE_KDE = "--kde" in sys.argv
ORIG = json.load(open(HERE / "pronostico_dr3_lrg2.json"))


def kuiper(theta):
    """Kuiper V for angles in [0, 2pi) against uniform; asymptotic p-value (Stephens 1970)."""
    u = np.sort(theta / (2 * np.pi))
    n = len(u)
    i = np.arange(1, n + 1)
    dplus = np.max(i / n - u)
    dminus = np.max(u - (i - 1) / n)
    V = dplus + dminus
    lam = (np.sqrt(n) + 0.155 + 0.24 / np.sqrt(n)) * V
    j = np.arange(1, 101)
    p = float(np.sum(2 * (4 * j**2 * lam**2 - 1) * np.exp(-2 * j**2 * lam**2)))
    return float(V), min(1.0, max(0.0, p))


def rayleigh(theta):
    n = len(theta)
    Rbar = np.hypot(np.mean(np.cos(theta)), np.mean(np.sin(theta)))
    Z = n * Rbar**2
    p = np.exp(-Z) * (1 + (2 * Z - Z**2) / (4 * n) - (24 * Z - 132 * Z**2 + 76 * Z**3 - 9 * Z**4) / (288 * n**2))
    return float(Rbar), float(Z), float(min(1.0, max(0.0, p)))


def tests_2d(wB):
    n = len(wB)
    r = np.hypot(wB[:, 0], wB[:, 1])
    th = np.mod(np.arctan2(wB[:, 1], wB[:, 0]), 2 * np.pi)
    Rbar, Z, p_ray = rayleigh(th)
    V, p_kui = kuiper(th)
    rayleigh_cdf = lambda x: 1 - np.exp(-x**2 / 2)
    ks = kstest(r, rayleigh_cdf)
    sect = np.floor(th / (2 * np.pi) * 8).astype(int)
    groups = [r[sect == k] for k in range(8)]
    kw = kruskal(*groups)
    ks_sect = [kstest(g, rayleigh_cdf) for g in groups]
    med_sect = [float(np.median(g)) for g in groups]
    q = np.quantile(r, [0.2, 0.4, 0.6, 0.8])
    rbin = np.searchsorted(q, r)
    tab = np.zeros((8, 5), int)
    np.add.at(tab, (sect, rbin), 1)
    chi2c = chi2_contingency(tab)
    sk = [float(skew(wB[:, k])) for k in range(2)]
    ku = [float(kurtosis(wB[:, k])) for k in range(2)]
    corr = float(np.corrcoef(wB.T)[0, 1])
    var = np.var(wB, axis=0, ddof=1).tolist()
    mean = wB.mean(axis=0).tolist()
    se_sk, se_ku, se_corr = np.sqrt(6 / n), np.sqrt(24 / n), 1 / np.sqrt(n)
    out = {"n": int(n),
           "angle_rayleigh": {"Rbar": Rbar, "Z": Z, "p": p_ray}, "angle_kuiper": {"V": V, "p": p_kui},
           "radius_ks_vs_chi2dof": {"D": float(ks.statistic), "p": float(ks.pvalue)},
           "radius_by_8_sectors": {"kruskal_H": float(kw.statistic), "kruskal_p": float(kw.pvalue), "median_per_sector": med_sect, "median_chi2dof_ref": float(np.sqrt(2 * np.log(2))),
                                   "ks_p_per_sector": [float(t.pvalue) for t in ks_sect], "min_ks_p_sector": float(min(t.pvalue for t in ks_sect))},
           "independence_chi2_8x5": {"chi2": float(chi2c[0]), "dof": int(chi2c[2]), "p": float(chi2c[1])},
           "components": {"mean": mean, "var": var, "skew": sk, "skew_se": float(se_sk), "excess_kurtosis": ku, "kurtosis_se": float(se_ku), "corr": corr, "corr_se": float(se_corr),
                          "skew_over_se": (np.array(sk) / se_sk).tolist(), "kurt_over_se": (np.array(ku) / se_ku).tolist(), "corr_over_se": corr / se_corr}}
    fails = []
    for k, p in (("rayleigh", p_ray), ("kuiper", p_kui), ("ks_radius", ks.pvalue), ("kruskal_sectors", kw.pvalue), ("independence", chi2c[1])):
        if p < 0.01:
            fails.append(k)
    if max(abs(s) for s in sk) > 0.05:
        fails.append("skew")
    if max(abs(k) for k in ku) > 0.1:
        fails.append("kurtosis")
    out["fails"] = fails
    return out


def kde_logpdf(train, bw, x, chunk=2000):
    """Gaussian KDE with bandwidth matrix bw^2 * Cov(train) (Scott-type scaling), log density at x."""
    cov = np.cov(train.T)
    Li = np.linalg.inv(np.linalg.cholesky(cov))
    tw = train @ Li.T / bw
    xw = x @ Li.T / bw
    lognorm = -np.log(len(train)) - np.log(2 * np.pi) - 2 * np.log(bw) - 0.5 * np.log(np.linalg.det(cov))
    out = np.empty(len(x))
    for i in range(0, len(x), chunk):
        xx = xw[i:i + chunk]
        d2 = (xx**2).sum(1)[:, None] + (tw**2).sum(1)[None, :] - 2 * xx @ tw.T
        m = (-0.5 * d2).max(axis=1, keepdims=True)
        out[i:i + chunk] = (m[:, 0] + np.log(np.exp(-0.5 * d2 - m).sum(axis=1))) + lognorm
    return out


def cv_bandwidth(train, rng, n_tr=8000, n_te=4000):
    """Likelihood cross-validation INSIDE A: fit on one random subset of A, score on another disjoint subset of A."""
    idx = rng.permutation(len(train))
    tr, te = train[idx[:n_tr]], train[idx[n_tr:n_tr + n_te]]
    grid = np.geomspace(0.05, 1.0, 25)
    scores = [float(np.mean(kde_logpdf(tr, b, te))) for b in grid]
    best = grid[int(np.argmax(scores))]
    # the chosen bandwidth is for n_tr points; rescale to the full A by Scott's n^(-1/6) law
    return float(best * (n_tr / len(train)) ** (1 / 6)), {"grid": grid.tolist(), "scores": scores, "best_on_subset": float(best)}


def main():
    t0 = time.time()
    out = {"rho2": 0.6, "cases": {}}
    for ic, combo in enumerate(("sinSN", "Pantheon+")):
        st = R.setup(combo)
        rho2, ir = 0.6, 0
        op = R.case_ops(combo, st, rho2)
        halves = {}
        rec = {}
        for ih, hyp in enumerate(("H0", "H1")):
            for half, code in (("A", 0), ("B", 1)):
                halves[(hyp, half)] = R.simulate(hyp, st, op, rho2, R.N_UNI // 2, np.random.default_rng([R.SEED_BASE, ic, ir, ih, code]))[0]
            rA, rB = halves[(hyp, "A")], halves[(hyp, "B")]
            muA, covA = rA.mean(axis=0), np.cov(rA.T)
            Li = np.linalg.inv(np.linalg.cholesky(covA))
            wB = (rB - muA) @ Li.T
            rec[hyp] = tests_2d(wB)
            # the original's Mahalanobis check, recomputed
            dm = (wB**2).sum(1)
            rec[hyp]["mahalanobis_pct_5_50_95"] = np.percentile(dm, [5, 50, 95]).tolist()
            rec[hyp]["mahalanobis_pct_json"] = ORIG["combos"][combo]["cases"]["rho2_0.6"]["hyp"][hyp]["gaussianity_B"]["mahalanobis_pct_5_50_95"]
            t = rec[hyp]
            print(f"{combo} {hyp}: Rayleigh p={t['angle_rayleigh']['p']:.3f} Kuiper p={t['angle_kuiper']['p']:.3f} KS radio p={t['radius_ks_vs_chi2dof']['p']:.3f} "
                  f"Kruskal p={t['radius_by_8_sectors']['kruskal_p']:.3f} (min KS sector p={t['radius_by_8_sectors']['min_ks_p_sector']:.3f}) indep p={t['independence_chi2_8x5']['p']:.3f} "
                  f"skew={np.round(t['components']['skew'], 4).tolist()} (/se {np.round(t['components']['skew_over_se'], 2).tolist()}) "
                  f"kurt={np.round(t['components']['excess_kurtosis'], 4).tolist()} (/se {np.round(t['components']['kurt_over_se'], 2).tolist()}) corr={t['components']['corr']:.4f} (/se {t['components']['corr_over_se']:.2f}) fails={t['fails']}", flush=True)
        any_fail = bool(rec["H0"]["fails"] or rec["H1"]["fails"])
        rec["any_fail"] = any_fail
        if any_fail or FORCE_KDE:
            rng = np.random.default_rng([R.SEED_BASE, ic, 7777])
            bws, dens = {}, {}
            for hyp in ("H0", "H1"):
                bw, info = cv_bandwidth(halves[(hyp, "A")], rng)
                bws[hyp] = {"bandwidth_rel_to_cov": bw, "cv": info}
            lam = {}
            for hyp in ("H0", "H1"):
                xB = halves[(hyp, "B")]
                lam[hyp] = kde_logpdf(halves[("H1", "A")], bws["H1"]["bandwidth_rel_to_cov"], xB) - kde_logpdf(halves[("H0", "A")], bws["H0"]["bandwidth_rel_to_cov"], xB)
            dec = R.decision(lam["H0"], lam["H1"])
            # gaussian decision for the same B, for the comparison
            dg = {h: multivariate_normal(mean=halves[(h, "A")].mean(axis=0), cov=np.cov(halves[(h, "A")].T)) for h in ("H0", "H1")}
            lg = {h: dg["H1"].logpdf(halves[(h, "B")]) - dg["H0"].logpdf(halves[(h, "B")]) for h in ("H0", "H1")}
            decg = R.decision(lg["H0"], lg["H1"])
            rec["kde"] = {"bandwidths": bws, "decision_kde": dec, "decision_gaussian_same_B": decg,
                          "diff_kde_minus_gaussian": {k: dec[k] - decg[k] for k in ("c0", "c1", "power")}}
            print(f"   KDE (bw {bws['H0']['bandwidth_rel_to_cov']:.3f}/{bws['H1']['bandwidth_rel_to_cov']:.3f}): c0={dec['c0']:.3f} c1={dec['c1']:.3f} power={dec['power']:.4f} | gaussian c0={decg['c0']:.3f} c1={decg['c1']:.3f} power={decg['power']:.4f}", flush=True)
        out["cases"][combo] = rec
    out["runtime_s"] = time.time() - t0
    (HERE / "revision_pronostico_gaussianidad.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")
    print("escrito revision_pronostico_gaussianidad.json", f"{out['runtime_s']:.0f} s")


if __name__ == "__main__":
    main()
