#!/usr/bin/env python3
"""
EBOSS_LRG.md: is the DESI DR2 LRG2 deviation (z_eff = 0.706; p local 0.028 under LCDM, PREDICTIVO_BINS.md) compatible with, or contradicted by,
the SDSS/eBOSS DR16 LRG BAO-only measurement at z_eff = 0.698 (Bautista et al. 2021, arXiv:2007.08993; Gil-Marin et al. 2020)? eBOSS cannot
confirm LRG2 (larger errors, shared volume). The expectation under each hypothesis is computed and committed BEFORE the eBOSS values are read.

Model (Gaussian approximation; fractional deviations from the LCDM prediction with p_ref = fit without DESI LRG2 of PREDICTIVO_BINS, "sinSN"
and "Pantheon+" separately; model evaluated at the z_eff of each measurement, 0.706 for DESI and 0.698 for eBOSS, with
pronostico_dr3_lrg2_v4.model_bao):
    delta_D = mu + J_D q + n_D,   delta_E = mu + J_E q + n_E,
mu = physical deviation (0 under H0; under H1 the same fractional deviation at z = 0.706 and 0.698: a physical assumption of the test, not a
consequence of the redshifts being close), q = p_true - p_ref ~ N(0, Sigma) with Sigma of the fit without DESI LRG2, n_D, n_E measurement errors
(fractional covariances C_D, C_E), cross covariance C_x = rho_x C_D^1/2 C_E^1/2 (symmetric square roots), rho_x = 0, 0.3, 0.5 (a hypothetical
family for sensitivity, not a measure of the shared volume; no published cross covariance exists).
J is the Jacobian of the fractional RESIDUAL: J = -d ln model / d p at p_ref, same definition for DESI and eBOSS (central differences).
Predictive conditionals:
    H0: V_D = C_D + J_D Sigma J_D^T, V_E = C_E + J_E Sigma J_E^T, V_DE = C_x + J_D Sigma J_E^T;
        delta_E | delta_D, H0 ~ N(V_DE^T V_D^-1 delta_D, V_E - V_DE^T V_D^-1 V_DE).
    H1 (mu free, flat prior): delta_E - delta_D = (J_E - J_D) q + (n_E - n_D); marginalising the common shift mu with a flat prior, observing
        delta_D carries no information on q (mu absorbs delta_D), so delta_E | delta_D, H1 ~ N(delta_D, C_Delta),
        C_Delta = C_D + C_E - C_x - C_x^T + (J_E - J_D) Sigma (J_E - J_D)^T.
Modes (order is mandatory):
    potencia  -> reads ONLY the metadata of the eBOSS data file (z_eff and row kinds: columns 0 and 2, never column 1) and the covariance file;
                 checks; 100000 simulated delta_E per conditional; power; eboss_lrg_potencia.json
    resultado -> (only after the commit of potencia) reads the eBOSS values; observed p local, pulls, signs, percentiles, compatibility,
                 Lambda_eBOSS; eboss_lrg_resultado.json
    ../venv/bin/python -I eboss_lrg_analisis.py potencia | resultado
Nothing existing is modified.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
from scipy.linalg import cho_factor, cho_solve
from scipy.stats import chi2 as chi2_dist

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import pronostico_dr3_lrg2_v4 as V4  # noqa: E402  (model_bao at arbitrary z; set-up; not modified)

D, COV13, LRG2, Z2 = V4.D, V4.COV13, V4.LRG2, V4.Z2
ROOT = HERE.parent.parent
BAO_DIR = ROOT / "data/desi_dr2_local/dr2_data/dr2-bao-zenodo/cosmology_chains/bao_data"
EBOSS_DAT, EBOSS_COV = BAO_DIR / "sdss_DR16_LRG_BAO_DMDH.dat", BAO_DIR / "sdss_DR16_LRG_BAO_DMDH_covtot.txt"
RHOS = (0.0, 0.3, 0.5)
COMBOS = ("sinSN", "Pantheon+")
N_SIM, N_CHECK, SEED = 100000, 2000, 20261116
EPS = (1e-4, 1e-5, 1e-3)


def sha256_file(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def read_metadata():
    """z_eff and row kinds of the eBOSS data file WITHOUT reading the values (column 1 is discarded unparsed)."""
    z, kinds = [], []
    with open(EBOSS_DAT) as f:
        for line in f:
            parts = line.split()
            if len(parts) == 3:
                z.append(float(parts[0]))
                kinds.append(parts[2])
    assert kinds == ["DM_over_rs", "DH_over_rs"] and len(set(z)) == 1, "formato inesperado del fichero de eBOSS"
    return z[0], kinds


def sqrtm_sym(C):
    w, U = np.linalg.eigh(0.5 * (C + C.T))
    assert np.all(w > 0)
    return (U * np.sqrt(w)) @ U.T


def pd(M):
    Ms = 0.5 * (M + M.T)
    w = np.linalg.eigvalsh(Ms)
    try:
        cho_factor(Ms, lower=True)
        ok = True
    except np.linalg.LinAlgError:
        ok = False
    return {"eigenvalues": w.tolist(), "pd": bool(ok and np.all(w > 0))}


def model_pair(p, z):
    """(D_M/r_d, D_H/r_d) of the LRG2 rows evaluated at redshift z (same code path as pronostico_dr3_lrg2_v4.model_bao)."""
    zz = Z2.copy()
    zz[LRG2] = z
    return V4.model_bao(p, zz)[LRG2]


def J_resid(p, z, npar):
    """J = -d ln model / d p (2 x npar), central differences."""
    J = np.zeros((2, npar))
    for jj in range(npar):
        pp, pm = np.array(p, float), np.array(p, float)
        pp[jj] += EPS[jj]
        pm[jj] -= EPS[jj]
        J[:, jj] = -(np.log(model_pair(pp, z)) - np.log(model_pair(pm, z))) / (2 * EPS[jj])
    return J


def setup(combo, C_eboss, z_eboss):
    pb = json.load(open(HERE / "predictivo_bins.json"))["combos"][combo]["bins"]["LRG2"]
    p_ref = np.array(pb["p0_minus_k"])
    Sigma = np.array(pb["Sigma"])
    npar = Sigma.shape[0]
    m_D, m_E = model_pair(p_ref, Z2[3]), model_pair(p_ref, z_eboss)
    delta_D = (D[LRG2] - m_D) / m_D
    assert float(np.max(np.abs(D[LRG2] - m_D - np.array(pb["r"])))) < 1e-9, "r de DESI LRG2 no reproduce PREDICTIVO_BINS"
    C_D = COV13[np.ix_(LRG2, LRG2)] / np.outer(m_D, m_D)
    C_E = C_eboss / np.outer(m_E, m_E)
    J_D, J_E = J_resid(p_ref, Z2[3], npar), J_resid(p_ref, z_eboss, npar)
    return dict(combo=combo, p_ref=p_ref, Sigma=Sigma, npar=npar, m_D=m_D, m_E=m_E, delta_D=delta_D, C_D=C_D, C_E=C_E, J_D=J_D, J_E=J_E, sD=sqrtm_sym(C_D), sE=sqrtm_sym(C_E))


def conditionals(st, rho):
    C_x = rho * st["sD"] @ st["sE"]
    S = st["Sigma"]
    V_D = st["C_D"] + st["J_D"] @ S @ st["J_D"].T
    V_E = st["C_E"] + st["J_E"] @ S @ st["J_E"].T
    V_DE = C_x + st["J_D"] @ S @ st["J_E"].T
    K = cho_solve(cho_factor(V_D, lower=True), V_DE).T                   # V_DE^T V_D^-1
    S0 = V_E - K @ V_DE
    S0 = 0.5 * (S0 + S0.T)
    dJ = st["J_E"] - st["J_D"]
    C_Delta = st["C_D"] + st["C_E"] - C_x - C_x.T + dJ @ S @ dJ.T
    C_Delta = 0.5 * (C_Delta + C_Delta.T)
    joint = np.block([[st["C_D"], C_x], [C_x.T, st["C_E"]]])
    checks = {"joint_meas_cov_pd": pd(joint), "C_Delta_pd": pd(C_Delta), "S_H0_pd": pd(S0), "V_E_pd": pd(V_E)}
    checks["all_ok"] = bool(all(checks[k]["pd"] for k in checks))
    return dict(rho=rho, C_x=C_x, V_D=V_D, V_E=V_E, V_DE=V_DE, m_H0=K @ st["delta_D"], S_H0=S0, m_H1=st["delta_D"].copy(), S_H1=C_Delta, checks=checks)


def sample(mean, cov, n, rng):
    return mean + rng.standard_normal((n, 2)) @ np.linalg.cholesky(cov).T


def d2_of(x, V):
    f = cho_factor(V, lower=True)
    x = np.atleast_2d(x)
    return np.einsum("ni,ni->n", x, cho_solve(f, x.T).T)


def logpdf(x, mean, cov):
    f = cho_factor(cov, lower=True)
    r = x - mean
    return float(-0.5 * (r @ cho_solve(f, r) + 2 * np.sum(np.log(np.diag(f[0]))) + 2 * np.log(2 * np.pi)))


def git_commit_of(path):
    try:
        return subprocess.run(["git", "log", "-1", "--format=%H", "--", str(path)], capture_output=True, text=True, cwd=HERE, check=True).stdout.strip() or None
    except Exception:
        return None


def potencia():
    t0 = time.time()
    z_eboss, kinds = read_metadata()
    C_eboss = np.loadtxt(EBOSS_COV)
    assert C_eboss.shape == (2, 2)
    out = {"mode": "potencia", "note": "se leen solo los metadatos (z_eff, tipo de fila) y la covarianza de eBOSS; los valores observados no se leen en este paso",
           "eboss_dat_file": str(EBOSS_DAT.relative_to(ROOT)), "eboss_dat_sha256": sha256_file(EBOSS_DAT), "eboss_cov_file": str(EBOSS_COV.relative_to(ROOT)), "eboss_cov_sha256": sha256_file(EBOSS_COV),
           "eboss_cov": C_eboss.tolist(), "eboss_sigma": np.sqrt(np.diag(C_eboss)).tolist(), "eboss_corr_DM_DH": float(C_eboss[0, 1] / np.sqrt(C_eboss[0, 0] * C_eboss[1, 1])), "z_eboss": z_eboss, "row_kinds": kinds,
           "rhos": RHOS, "N_sim": N_SIM, "N_check": N_CHECK, "seed": SEED, "combos": {}}
    rng = np.random.default_rng(SEED)
    for combo in COMBOS:
        st = setup(combo, C_eboss, z_eboss)
        rec = {"p_ref": st["p_ref"].tolist(), "Sigma": st["Sigma"].tolist(), "m_D": st["m_D"].tolist(), "m_E": st["m_E"].tolist(), "delta_D": st["delta_D"].tolist(),
               "C_D": st["C_D"].tolist(), "C_E": st["C_E"].tolist(), "J_D": st["J_D"].tolist(), "J_E": st["J_E"].tolist(), "rho": {}}
        for rho in RHOS:
            c = conditionals(st, rho)
            r = {k: (c[k].tolist() if isinstance(c[k], np.ndarray) else c[k]) for k in ("C_x", "V_D", "V_E", "V_DE", "m_H0", "S_H0", "m_H1", "S_H1", "checks")}
            assert c["checks"]["all_ok"], f"PARA: comprobación de matrices fallida ({combo}, rho_x = {rho})"
            chk = {}
            for hyp in ("H0", "H1"):
                x = sample(c[f"m_{hyp}"], c[f"S_{hyp}"], N_CHECK, rng)
                sd = np.sqrt(np.diag(c[f"S_{hyp}"]))
                dm = (x.mean(axis=0) - c[f"m_{hyp}"]) / sd
                dv = np.diag(np.cov(x.T)) / np.diag(c[f"S_{hyp}"]) - 1
                chk[hyp] = {"mean_diff_over_sigma": dm.tolist(), "mc_sigma_mean": float(1 / np.sqrt(N_CHECK)), "var_rel_diff": dv.tolist(), "mc_sigma_var": float(np.sqrt(2 / N_CHECK)),
                            "ok": bool(np.all(np.abs(dm) < 4 / np.sqrt(N_CHECK)) and np.all(np.abs(dv) < 4 * np.sqrt(2 / N_CHECK)))}
            assert chk["H0"]["ok"] and chk["H1"]["ok"], "PARA: la muestra pequeña no reproduce la media/covarianza analítica"
            r["check_small_sample"] = chk
            sims = {}
            for hyp in ("H0", "H1"):
                x = sample(c[f"m_{hyp}"], c[f"S_{hyp}"], N_SIM, rng)
                d2 = d2_of(x, c["V_E"])
                p = chi2_dist.sf(d2, 2)
                sims[hyp] = {"P_p_lt_0.05": float(np.mean(p < 0.05)), "P_p_lt_0.01": float(np.mean(p < 0.01)), "d2_pct_5_50_95": np.percentile(d2, [5, 50, 95]).tolist(), "d2_mean": float(d2.mean()),
                             "d2_sorted_every_100": np.sort(d2)[::100].tolist()}
            r["sims"] = sims
            r["power_P_p_lt_0.05_given_H1"], r["P_p_lt_0.05_given_H0"] = sims["H1"]["P_p_lt_0.05"], sims["H0"]["P_p_lt_0.05"]
            rec["rho"][f"rho_{rho:g}"] = r
            print(f"{combo} rho_x={rho}: checks ok; P(p<0.05|H1)={sims['H1']['P_p_lt_0.05']:.4f} P(p<0.05|H0)={sims['H0']['P_p_lt_0.05']:.4f}; d2 median H0 {sims['H0']['d2_pct_5_50_95'][1]:.2f} H1 {sims['H1']['d2_pct_5_50_95'][1]:.2f}", flush=True)
        out["combos"][combo] = rec
    p20 = out["combos"]["sinSN"]["rho"]["rho_0"]["power_P_p_lt_0.05_given_H1"]
    out["P20"] = {"statement": "Bajo H1 (desvío de LRG2 real, con el tamaño y la incertidumbre estimados en DESI), la probabilidad de que eBOSS dé p local < 0.05 frente a ΛCDM es < 0.4. (Se evalúa con ρ_x = 0; los demás ρ_x son sensibilidad, no nuevas predicciones.)",
                  "power_sinSN_rho0": p20, "power_Pantheon+_rho0": out["combos"]["Pantheon+"]["rho"]["rho_0"]["power_P_p_lt_0.05_given_H1"], "holds": bool(p20 < 0.4)}
    out["runtime_s"] = time.time() - t0
    (HERE / "eboss_lrg_potencia.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")
    print(json.dumps(out["P20"], indent=1))
    print("escrito eboss_lrg_potencia.json")


def resultado():
    pot = json.load(open(HERE / "eboss_lrg_potencia.json"))
    commit = git_commit_of(HERE / "eboss_lrg_potencia.json")
    assert commit, "PARA: eboss_lrg_potencia.json no está commiteado todavía"
    assert sha256_file(EBOSS_COV) == pot["eboss_cov_sha256"] and sha256_file(EBOSS_DAT) == pot["eboss_dat_sha256"]
    C_eboss = np.array(pot["eboss_cov"])
    z_eboss = pot["z_eboss"]
    raw = np.genfromtxt(EBOSS_DAT, dtype=None, encoding="ascii")        # now the values are read
    D_E = np.array([float(r[1]) for r in raw])
    out = {"mode": "resultado", "potencia_commit": commit, "eboss_dat_sha256": sha256_file(EBOSS_DAT), "eboss_values_DM_DH_over_rd": D_E.tolist(), "z_eboss": z_eboss, "combos": {}}
    for combo in COMBOS:
        st = setup(combo, C_eboss, z_eboss)
        delta_E = (D_E - st["m_E"]) / st["m_E"]
        rec = {"m_E": st["m_E"].tolist(), "delta_E": delta_E.tolist(), "delta_D": st["delta_D"].tolist(), "rho": {}}
        for rho in RHOS:
            c = conditionals(st, rho)
            d2E = float(d2_of(delta_E, c["V_E"])[0])
            pulls = (delta_E / np.sqrt(np.diag(c["V_E"]))).tolist()          # sigmas of V_E (eBOSS vs LCDM on its own), not of the conditional
            Delta = st["delta_D"] - delta_E
            d2c = float(Delta @ cho_solve(cho_factor(c["S_H1"], lower=True), Delta))
            lam = logpdf(delta_E, c["m_H1"], c["S_H1"]) - logpdf(delta_E, c["m_H0"], c["S_H0"])
            sims = pot["combos"][combo]["rho"][f"rho_{rho:g}"]["sims"]
            pct = {h: float(np.mean(np.array(sims[h]["d2_sorted_every_100"]) <= d2E)) for h in ("H0", "H1")}
            r = {"d2_E": d2E, "p_local": float(chi2_dist.sf(d2E, 2)), "pull_marginal_V_E": pulls, "signs": [int(np.sign(v)) for v in delta_E], "same_signs_as_DESI_LRG2": bool(np.all(np.sign(delta_E) == np.sign(st["delta_D"]))),
                 "percentile_d2E_in_sims": pct, "compat": {"Delta": Delta.tolist(), "d2": d2c, "p": float(chi2_dist.sf(d2c, 2))}, "Lambda_eBOSS": lam,
                 "logp_H1": logpdf(delta_E, c["m_H1"], c["S_H1"]), "logp_H0": logpdf(delta_E, c["m_H0"], c["S_H0"])}
            rec["rho"][f"rho_{rho:g}"] = r
            print(f"{combo} rho_x={rho}: p_local={r['p_local']:.4f} pulls={np.round(pulls, 2).tolist()} pct H0/H1={pct['H0']:.3f}/{pct['H1']:.3f} compat p={r['compat']['p']:.4f} Lambda={lam:+.3f}", flush=True)
        out["combos"][combo] = rec
    s0 = out["combos"]["sinSN"]["rho"]["rho_0"]
    out["P21"] = {"statement": "El p local observado de eBOSS frente a ΛCDM es > 0.05.", "p_local_sinSN": s0["p_local"], "p_local_Pantheon+": out["combos"]["Pantheon+"]["rho"]["rho_0"]["p_local"], "holds": bool(s0["p_local"] > 0.05)}
    out["P22"] = {"statement": "Con ρ_x = 0, DESI LRG2 y eBOSS son compatibles entre sí (p > 0.05 para su diferencia).", "p_compat_sinSN": s0["compat"]["p"], "p_compat_Pantheon+": out["combos"]["Pantheon+"]["rho"]["rho_0"]["compat"]["p"], "holds": bool(s0["compat"]["p"] > 0.05)}
    (HERE / "eboss_lrg_resultado.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")
    print(json.dumps({"P21": out["P21"], "P22": out["P22"]}, indent=1))
    print("escrito eboss_lrg_resultado.json")


if __name__ == "__main__":
    {"potencia": potencia, "resultado": resultado}[sys.argv[1]]()
