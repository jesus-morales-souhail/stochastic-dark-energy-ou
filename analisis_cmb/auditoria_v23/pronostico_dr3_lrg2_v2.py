#!/usr/bin/env python3
"""
PRONOSTICO_DR3_REGISTRO.md: version 2 of the DR3 forecast of the LRG2 residual, with the changes asked by the independent review
(REVISION_PRONOSTICO_DR3.md, commit 99f879d). Imports the set-up and fitters of pronostico_dr3_lrg2.py (not modified; its main()
is not run). Same assumptions, hypotheses, thresholds (5 % / 5 %, |Lambda| > 3), P18 and P19 as the draft. Differences:
  - EXACT densities of the linearised model: r3 = a + M z with z standard normal (z_p, z_e[, z_mu]) -> N(a, M M^T) under H0 and H1;
    Lambda uses these densities. The simulations are a check (compared with the exact moments and with the draft's numbers).
  - c0 and c1 calibrated with the SAME half-B simulations as the draft (same seeds). Monte Carlo errors of c0, c1, power, P(Lambda < c1 | H0),
    inconclusive fractions, medians of d2 and P(|Lambda| > 3) from 20 re-runs of B with new seeds (recorded).
  - rho2 per bin (vector of 7) besides the common value; table of power against a common rho2 from 0.40 to 0.80 (step 0.02).
  - "DR3 evaluation" mode (not used yet): `calibrar` reads ONLY the DR3 covariance, computes rho2_k, applies the rupture rules, recalibrates
    c0 and c1 with the frozen seeds and writes the calibration with its sha256; `evaluar` then reads the DR3 values and computes Lambda
    once (refuses to overwrite a result). Tested with fake data (`falso`: DR2 rescaled), never with real DR3 data.

    ../venv/bin/python -I pronostico_dr3_lrg2_v2.py                       # registration run -> pronostico_dr3_lrg2_v2.json
    ../venv/bin/python -I pronostico_dr3_lrg2_v2.py falso                 # writes fake DR3 files (DR2 rescaled) for the test
    ../venv/bin/python -I pronostico_dr3_lrg2_v2.py calibrar COV.json TAG # -> pronostico_dr3_lrg2_v2_calibracion_TAG.json (+ .sha256)
    ../venv/bin/python -I pronostico_dr3_lrg2_v2.py evaluar COV.json VALORES.json TAG   # -> pronostico_dr3_lrg2_v2_resultado_TAG.json
DR3 input files (JSON): {"labels": [13 labels in the DR2 order], "cov": 13x13} and {"labels": [...], "mean": [13]}.
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
from scipy.stats import chi2 as chi2_dist, multivariate_normal

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import pronostico_dr3_lrg2 as P  # noqa: E402

PB, JB = P.PB, P.JB
D, COV13, LRG2, KEEP, C22 = P.D, P.COV13, P.LRG2, P.KEEP, P.C22
SEED_BASE, N_UNI = P.SEED_BASE, P.N_UNI
NB = N_UNI // 2
RHO2S = P.RHO2S
COMBOS = ("sinSN", "Pantheon+")
BINS = JB.BINS
BIN_NAMES = list(BINS)
LABELS = JB.LABELS
RHO2_GRID = [round(0.40 + 0.02 * k, 2) for k in range(21)]
RHO2_LRG2_RANGE = (0.40, 0.80)
N_REP = 20
LAMBDA_STRONG = 3.0
SEED_EVAL = 500          # seeds of the DR3 evaluation mode: [SEED_BASE, ic, SEED_EVAL, ih, 1] (frozen here, before DR3)
SEED_GRID = 100          # seeds of the rho2 sweep: [SEED_BASE, ic, SEED_GRID + k, ih, 1]
SEED_REP = 2000          # seeds of the MC-error re-runs: [SEED_BASE, ic, ir, ih, 1, SEED_REP + rep]


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def git_head():
    try:
        return subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True, cwd=HERE, check=True).stdout.strip()
    except Exception:
        return None


def rho2_rows(rho2):
    """Common value or per-bin vector (7, in BIN_NAMES order) -> per-row vector (13)."""
    if np.isscalar(rho2):
        return np.full(13, float(rho2))
    r = np.zeros(13)
    for v, (name, rows) in zip(rho2, BINS.items()):
        r[rows] = v
    return r


# ----------------------------------------------------------------------------------------------------------------------
# DR2 set-up (same calls as the draft) and per-case operators
# ----------------------------------------------------------------------------------------------------------------------
def setup(combo):
    npar = P.npar_of(combo)
    chi2_full, p_full = PB.fit_d(combo, list(range(13)), D)
    _, _, Sf, _, _ = P.fisher(combo, list(range(13)), p_full, COV13)
    chi2_2, p2 = PB.fit_d(combo, KEEP, D)
    J, F2, S2, _, rows2 = P.fisher(combo, KEEP, p2, COV13)
    m2 = JB.model_at(p2)
    r2 = (D - m2)[LRG2]
    J2 = J[LRG2]
    Q = J2 @ S2 @ J2.T
    return dict(combo=combo, npar=npar, p_full=p_full, chi2_full=chi2_full, Sf=Sf, p2=p2, chi2_2=chi2_2, S2=S2, J=J, J2=J2, m2=m2, r2=r2, Q=Q, V2=C22 + Q,
                m2_bao=m2[:13], m2_rest=m2[rows2[len(KEEP):]], d_rest=D[rows2[len(KEEP):]],
                L_full=np.linalg.cholesky(Sf), L2=np.linalg.cholesky(S2), L22=np.linalg.cholesky(C22), L13=np.linalg.cholesky(COV13))


def case_ops(st, rho2):
    rr = rho2_rows(rho2)
    s = np.sqrt(rr)
    cov3 = COV13 * np.outer(s, s)                       # block k: rho2_k C_kk
    _, _, S3, CJ3, _ = P.fisher(st["combo"], KEEP, st["p2"], cov3)
    G3 = S3 @ CJ3.T
    const_rest = G3[:, len(KEEP):] @ (st["d_rest"] - st["m2_rest"])
    V3k = cov3[np.ix_(LRG2, LRG2)] + st["J2"] @ S3 @ st["J2"].T
    Le = np.diag(np.sqrt(rr * (1 - rr))) @ st["L13"]     # Cholesky factor of Cov(e) = diag(sqrt(rho2(1-rho2))) C diag(...)
    return dict(rho2_rows=rr, cov3=cov3, S3=S3, G3_bao=G3[:, :len(KEEP)], const_rest=const_rest, V3k=V3k, V3ki=np.linalg.inv(V3k), Le=Le)


def gen_r3(hyp, st, op, z_p, z_e, z_mu):
    """Linearised generation + refit, vectorised over universes (rows). Same equations and draw order as the draft."""
    npar, rr = st["npar"], op["rho2_rows"]
    p_true = (st["p_full"][:npar] + z_p @ st["L_full"].T) if hyp == "H0" else (st["p2"][:npar] + z_p @ st["L2"].T)
    m_true = st["m2_bao"][None, :] + (p_true - st["p2"][:npar]) @ st["J"][:13].T
    e = z_e @ op["Le"].T
    x = m_true.copy()
    if hyp == "H1":
        mu = D[LRG2][None, :] - m_true[:, LRG2] + z_mu @ st["L22"].T
        x[:, LRG2] += mu                                   # x = m_true + E2 mu
    n2 = D[:13][None, :] - x                               # n2 = D - m_true (- E2 mu)
    d3 = x + rr[None, :] * n2 + e                          # d3 = m_true (+ E2 mu) + rho2 n2 + e
    dp3 = (d3[:, KEEP] - st["m2_bao"][KEEP][None, :]) @ op["G3_bao"].T + op["const_rest"][None, :]
    return d3[:, LRG2] - st["m2_bao"][LRG2][None, :] - dp3 @ st["J2"].T


def simulate(hyp, st, op, n, rng):
    npar = st["npar"]
    z_p = rng.standard_normal((n, npar))
    z_e = rng.standard_normal((n, 13))
    z_mu = rng.standard_normal((n, 2)) if hyp == "H1" else None
    r3 = gen_r3(hyp, st, op, z_p, z_e, z_mu)
    return r3, np.einsum("ni,ij,nj->n", r3, op["V3ki"], r3)


def exact_density(hyp, st, op):
    """r3 is affine in z = (z_p, z_e[, z_mu]): a = r3(0), M_k = r3(e_k) - a. Returns a, Sigma = M M^T, linearity check."""
    npar = st["npar"]
    nz = npar + 13 + (2 if hyp == "H1" else 0)

    def r3_of(z):
        z = np.atleast_2d(z)
        return gen_r3(hyp, st, op, z[:, :npar], z[:, npar:npar + 13], z[:, npar + 13:] if hyp == "H1" else None)[0]

    a = r3_of(np.zeros(nz))
    M = np.column_stack([r3_of(np.eye(nz)[k]) - a for k in range(nz)])
    lin = max(float(np.max(np.abs(r3_of(2 * np.eye(nz)[k]) - a - 2 * M[:, k]))) for k in range(nz))
    Sig = M @ M.T
    return a, Sig, lin, {"a": a.tolist(), "Sigma": Sig.tolist(), "E_d2": float(a @ op["V3ki"] @ a + np.trace(op["V3ki"] @ Sig)), "linearity_check_max_abs": lin}


def lam_fun(dens):
    g0, g1 = multivariate_normal(mean=dens["H0"][0], cov=dens["H0"][1]), multivariate_normal(mean=dens["H1"][0], cov=dens["H1"][1])
    return lambda r: g1.logpdf(r) - g0.logpdf(r)


def decide(lam0, lam1, d20=None, d21=None):
    c0, c1 = float(np.percentile(lam0, 95)), float(np.percentile(lam1, 5))
    d = {"c0": c0, "c1": c1, "c1_lt_c0": bool(c1 < c0), "power": float(np.mean(lam1 > c0)), "P_lt_c1_H0": float(np.mean(lam0 < c1)),
         "inconclusive_H0": float(np.mean((lam0 >= c1) & (lam0 <= c0))), "inconclusive_H1": float(np.mean((lam1 >= c1) & (lam1 <= c0))),
         "P_gt3_H0": float(np.mean(lam0 > LAMBDA_STRONG)), "P_gt3_H1": float(np.mean(lam1 > LAMBDA_STRONG)), "P_ltm3_H0": float(np.mean(lam0 < -LAMBDA_STRONG)),
         "P_ltm3_H1": float(np.mean(lam1 < -LAMBDA_STRONG))}
    if d20 is not None:
        d["med_d2_H0"], d["med_d2_H1"] = float(np.median(d20)), float(np.median(d21))
    return d


def calibrate_case(st, rho2, seeds):
    """Exact densities + c0/c1 from half-B simulations with the given seeds {H0: seed, H1: seed}; returns record and the B samples."""
    op = case_ops(st, rho2)
    dens, exact = {}, {}
    for hyp in ("H0", "H1"):
        a, Sig, lin, rec = exact_density(hyp, st, op)
        dens[hyp], exact[hyp] = (a, Sig), rec
    lam = lam_fun(dens)
    B = {hyp: simulate(hyp, st, op, NB, np.random.default_rng(seeds[hyp])) for hyp in ("H0", "H1")}
    lamB = {hyp: lam(B[hyp][0]) for hyp in B}
    dec = decide(lamB["H0"], lamB["H1"], B["H0"][1], B["H1"][1])
    dec["Lambda_pct_5_50_95"] = {h: np.percentile(lamB[h], [5, 50, 95]).tolist() for h in lamB}
    dec["d2_pct_5_50_95"] = {h: np.percentile(B[h][1], [5, 50, 95]).tolist() for h in B}
    return {"rho2_rows": op["rho2_rows"].tolist(), "V3k": op["V3k"].tolist(), "Sigma3": op["S3"].tolist(), "exact": exact, "seeds_B": seeds, "decision": dec}, op, dens, lam, B


# ----------------------------------------------------------------------------------------------------------------------
# registration run
# ----------------------------------------------------------------------------------------------------------------------
def registration():
    t0 = time.time()
    orig = json.load(open(HERE / "pronostico_dr3_lrg2.json"))
    out = {"version": "v2 (cambios de REVISION_PRONOSTICO_DR3.md)", "script_sha256": sha256_file(__file__), "git_head_at_run": git_head(), "N_universes_per_case": N_UNI, "N_B": NB,
           "seeds": {"B_same_as_draft": "[SEED_BASE, i_comb, i_rho2, i_hyp, 1]", "mc_reruns": f"[SEED_BASE, i_comb, i_rho2, i_hyp, 1, {SEED_REP} + rep], rep = 0..{N_REP - 1}",
                     "rho2_sweep": f"[SEED_BASE, i_comb, {SEED_GRID} + k, i_hyp, 1], k = index in rho2_grid", "dr3_evaluation": f"[SEED_BASE, i_comb, {SEED_EVAL}, i_hyp, 1]"},
           "rho2_grid": RHO2_GRID, "rho2_LRG2_range": RHO2_LRG2_RANGE, "combos": {}}
    for ic, combo in enumerate(COMBOS):
        st = setup(combo)
        o2 = orig["combos"][combo]["DR2"]
        rec = {"DR2": {"p_full": st["p_full"].tolist(), "p2": st["p2"].tolist(), "r2": st["r2"].tolist(), "Q": st["Q"].tolist(), "V2": st["V2"].tolist(),
                       "check_vs_draft_max_abs": float(max(np.max(np.abs(st["r2"] - np.array(o2["r2"]))), np.max(np.abs(st["p2"] - np.array(o2["p2_without_LRG2"]))), np.max(np.abs(st["V2"] - np.array(o2["V2"])))))},
               "cases": {}, "rho2_sweep": {}}
        for ir, rho2 in enumerate(RHO2S):
            oc = orig["combos"][combo]["cases"][f"rho2_{rho2:g}"]
            seeds = {hyp: [SEED_BASE, ic, ir, ih, 1] for ih, hyp in enumerate(("H0", "H1"))}
            for hyp in seeds:
                assert seeds[hyp] == oc["hyp"][hyp]["seeds"]["B"], "semilla B distinta de la del borrador"
            case, op, dens, lam, B = calibrate_case(st, rho2, seeds)
            case["V3k_vs_draft_max_abs"] = float(np.max(np.abs(op["V3k"] - np.array(oc["V3k_linearised"]))))
            # reproduction of the draft's estimator (A-estimated densities) with the draft's A seeds
            A = {hyp: simulate(hyp, st, op, NB, np.random.default_rng([SEED_BASE, ic, ir, ih, 0])) for ih, hyp in enumerate(("H0", "H1"))}
            densA = {hyp: (A[hyp][0].mean(axis=0), np.cov(A[hyp][0].T)) for hyp in A}
            lamA = lam_fun(densA)
            decA = decide(lamA(B["H0"][0]), lamA(B["H1"][0]), B["H0"][1], B["H1"][1])
            od = oc["decision"]
            case["draft_reproduction"] = {"c0": decA["c0"], "c1": decA["c1"], "power": decA["power"], "diff_vs_draft_json": {"c0": decA["c0"] - od["c0_pct95_H0"], "c1": decA["c1"] - od["c1_pct5_H1"], "power": decA["power"] - od["power_P_Lambda_gt_c0_given_H1"]},
                                          "d2_pct_diff_vs_draft": {h: float(np.max(np.abs(np.array(case["decision"]["d2_pct_5_50_95"][h]) - np.array(oc["hyp"][h]["d2_percentiles_5_50_95_B"])))) for h in ("H0", "H1")}}
            # simulations (A and B) against the exact moments
            case["sim_vs_exact"] = {}
            for hyp in ("H0", "H1"):
                a, Sig = dens[hyp]
                sd = np.sqrt(np.diag(Sig))
                rAB = np.vstack([A[hyp][0], B[hyp][0]])
                case["sim_vs_exact"][hyp] = {"A_mean_minus_exact_over_sigma": ((A[hyp][0].mean(axis=0) - a) / sd).tolist(), "B_mean_minus_exact_over_sigma": ((B[hyp][0].mean(axis=0) - a) / sd).tolist(),
                                             "AB_mean_minus_exact_over_sigma": ((rAB.mean(axis=0) - a) / sd).tolist(), "mc_sigma_mean_over_sigma_A_or_B": float(1 / np.sqrt(NB)),
                                             "A_var_over_exact_minus_1": (np.diag(np.cov(A[hyp][0].T)) / np.diag(Sig) - 1).tolist(), "B_var_over_exact_minus_1": (np.diag(np.cov(B[hyp][0].T)) / np.diag(Sig) - 1).tolist(),
                                             "mc_rel_sigma_var": float(np.sqrt(2 / NB)), "draft_A_mean_minus_exact_over_sigma": ((np.array(oc["hyp"][hyp]["r3_mean_A"]) - a) / sd).tolist(),
                                             "d2_mean_B_vs_exact_E_d2": [float(B[hyp][1].mean()), case["exact"][hyp]["E_d2"]]}
            # Monte Carlo errors: N_REP re-runs of B with new seeds
            reps = []
            seeds_rep = []
            for rep in range(N_REP):
                sr = {hyp: [SEED_BASE, ic, ir, ih, 1, SEED_REP + rep] for ih, hyp in enumerate(("H0", "H1"))}
                seeds_rep.append(sr)
                Br = {hyp: simulate(hyp, st, op, NB, np.random.default_rng(sr[hyp])) for hyp in sr}
                reps.append(decide(lam(Br["H0"][0]), lam(Br["H1"][0]), Br["H0"][1], Br["H1"][1]))
            keys = [k for k in reps[0] if k != "c1_lt_c0"]
            case["mc_error"] = {"n_reruns": N_REP, "seeds": seeds_rep, "stats": {k: {"mean": float(np.mean([r[k] for r in reps])), "std": float(np.std([r[k] for r in reps], ddof=1))} for k in keys},
                                "c1_lt_c0_in_all_reruns": bool(all(r["c1_lt_c0"] for r in reps))}
            rec["cases"][f"rho2_{rho2:g}"] = case
            d = case["decision"]
            e = case["mc_error"]["stats"]
            print(f"{combo} rho2={rho2}: exact c0={d['c0']:.3f}±{e['c0']['std']:.3f} c1={d['c1']:.3f}±{e['c1']['std']:.3f} power={d['power']:.4f}±{e['power']['std']:.4f} "
                  f"(draft {od['c0_pct95_H0']:.3f} / {od['c1_pct5_H1']:.3f} / {od['power_P_Lambda_gt_c0_given_H1']:.4f}; repro diffs {case['draft_reproduction']['diff_vs_draft_json']}) "
                  f"med d2 {d['med_d2_H0']:.2f}±{e['med_d2_H0']['std']:.2f} / {d['med_d2_H1']:.2f}±{e['med_d2_H1']['std']:.2f}; E_d2 exact {case['exact']['H0']['E_d2']:.3f} / {case['exact']['H1']['E_d2']:.3f}", flush=True)
        # power against a common rho2, continuous
        for k, rho2 in enumerate(RHO2_GRID):
            seeds = {hyp: [SEED_BASE, ic, SEED_GRID + k, ih, 1] for ih, hyp in enumerate(("H0", "H1"))}
            case, op, dens, lam, B = calibrate_case(st, rho2, seeds)
            d = case["decision"]
            rec["rho2_sweep"][f"{rho2:.2f}"] = {"seeds_B": seeds, "c0": d["c0"], "c1": d["c1"], "c1_lt_c0": d["c1_lt_c0"], "power": d["power"], "P_lt_c1_H0": d["P_lt_c1_H0"],
                                                "inconclusive_H0": d["inconclusive_H0"], "inconclusive_H1": d["inconclusive_H1"], "med_d2_H0": d["med_d2_H0"], "med_d2_H1": d["med_d2_H1"],
                                                "E_d2_exact": {h: case["exact"][h]["E_d2"] for h in ("H0", "H1")}}
        print(f"{combo} sweep: " + " ".join(f"{r}:{v['power']:.3f}" for r, v in rec["rho2_sweep"].items()), flush=True)
        out["combos"][combo] = rec
    m = out["combos"]["sinSN"]["cases"]["rho2_0.6"]
    pw, e = m["decision"]["power"], m["mc_error"]["stats"]
    out["predictions"] = {"P18": {"statement": "rho2 = 0.6, falsa alarma 5 %: potencia entre 0.5 y 0.75 (sin SN)", "power": pw, "mc_std": e["power"]["std"], "holds": bool(0.5 <= pw <= 0.75),
                                  "power_Pantheon+": out["combos"]["Pantheon+"]["cases"]["rho2_0.6"]["decision"]["power"]},
                          "P19": {"statement": "rho2 = 0.6: mediana de d2_DR3 entre 3 y 6 bajo H0 y entre 9 y 15 bajo H1 (sin SN)", "median_H0": m["decision"]["med_d2_H0"], "median_H1": m["decision"]["med_d2_H1"],
                                  "mc_std": [e["med_d2_H0"]["std"], e["med_d2_H1"]["std"]], "holds": bool(3 <= m["decision"]["med_d2_H0"] <= 6 and 9 <= m["decision"]["med_d2_H1"] <= 15)}}
    out["runtime_s"] = time.time() - t0
    (HERE / "pronostico_dr3_lrg2_v2.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")
    print(json.dumps(out["predictions"], indent=1))
    print("escrito pronostico_dr3_lrg2_v2.json", f"{out['runtime_s']:.0f} s")


# ----------------------------------------------------------------------------------------------------------------------
# DR3 evaluation mode: calibrar (covariance only) -> evaluar (values, once)
# ----------------------------------------------------------------------------------------------------------------------
def read_dr3_cov(path):
    j = json.load(open(path))
    assert list(j["labels"]) == LABELS, "etiquetas o orden distintos de DR2: no evaluable sin conversión declarada"
    C3 = np.array(j["cov"], float)
    assert C3.shape == (13, 13)
    return C3


def rupture_checks(C3):
    """Operational rho2_k and the rupture rules of PRONOSTICO_DR3_REGISTRO.md (b, c)."""
    rc = np.diag(C3) / np.diag(COV13)
    offblock = np.ones((13, 13), bool)
    for rows in BINS.values():
        for i in rows:
            for jj in rows:
                offblock[i, jj] = False
    rep = {"rho2_per_component": rc.tolist(), "rho2_per_bin": {}, "component_rel_spread": {}, "broken": [], "max_abs_offblock_cov3": float(np.max(np.abs(C3[offblock])))}
    if rep["max_abs_offblock_cov3"] != 0.0:
        rep["broken"].append("correlaciones entre bins distintas de cero en la covarianza de DR3")
    for name, rows in BINS.items():
        v = rc[rows]
        rep["rho2_per_bin"][name] = float(v.mean())
        spread = float((v.max() - v.min()) / v.mean()) if len(rows) > 1 else 0.0
        rep["component_rel_spread"][name] = spread
        if spread > 0.20:
            rep["broken"].append(f"componentes de {name} con rho2 que difieren más del 20 % ({spread:.2f})")
    corr2 = COV13[3, 4] / np.sqrt(COV13[3, 3] * COV13[4, 4])
    corr3 = C3[3, 4] / np.sqrt(C3[3, 3] * C3[4, 4])
    rep["corr_DM_DH_LRG2"] = {"DR2": float(corr2), "DR3": float(corr3), "change": float(corr3 - corr2)}
    if abs(corr3 - corr2) > 0.1:
        rep["broken"].append(f"correlación D_M–D_H de LRG2 cambia {corr3 - corr2:+.3f} (> 0.1)")
    r2l = rep["rho2_per_bin"]["LRG2"]
    if not (RHO2_LRG2_RANGE[0] <= r2l <= RHO2_LRG2_RANGE[1]):
        rep["broken"].append(f"rho2_LRG2 = {r2l:.3f} fuera de [{RHO2_LRG2_RANGE[0]}, {RHO2_LRG2_RANGE[1]}]")
    rep["status"] = "modelo roto" if rep["broken"] else "ok"
    return rep


def calibrar(cov_path, tag):
    cov_sha = sha256_file(cov_path)
    C3 = read_dr3_cov(cov_path)
    rep = rupture_checks(C3)
    cal = {"mode": "calibrar", "tag": tag, "cov_file": str(cov_path), "cov_sha256": cov_sha, "script_sha256": sha256_file(__file__), "git_head_at_run": git_head(), "rupture_checks": rep, "combos": {}}
    if rep["status"] == "ok":
        rho2_bins = [rep["rho2_per_bin"][n] for n in BIN_NAMES]
        for ic, combo in enumerate(COMBOS):
            st = setup(combo)
            seeds = {hyp: [SEED_BASE, ic, SEED_EVAL, ih, 1] for ih, hyp in enumerate(("H0", "H1"))}
            case, op, dens, lam, B = calibrate_case(st, rho2_bins, seeds)
            cal["combos"][combo] = {"rho2_bins": rho2_bins, "seeds_B": seeds, "exact_densities": case["exact"], "V3k": case["V3k"], "c0": case["decision"]["c0"], "c1": case["decision"]["c1"],
                                    "c1_lt_c0": case["decision"]["c1_lt_c0"], "power": case["decision"]["power"], "P_lt_c1_H0": case["decision"]["P_lt_c1_H0"], "p2": st["p2"].tolist(), "r2": st["r2"].tolist()}
    outp = HERE / f"pronostico_dr3_lrg2_v2_calibracion_{tag}.json"
    assert not outp.exists(), f"ya existe {outp.name}: la calibración no se sobrescribe"
    outp.write_text(json.dumps(cal, indent=1, default=float), encoding="utf-8")
    sha = sha256_file(outp)
    (HERE / f"pronostico_dr3_lrg2_v2_calibracion_{tag}.sha256").write_text(f"{sha}  {outp.name}\n", encoding="utf-8")
    print(json.dumps({"status": rep["status"], "broken": rep["broken"], "rho2_per_bin": rep["rho2_per_bin"], "calibracion": {c: {k: cal["combos"][c][k] for k in ("c0", "c1", "c1_lt_c0", "power")} for c in cal["combos"]}, "sha256": sha}, indent=1))
    return outp


def evaluar(cov_path, val_path, tag):
    calp = HERE / f"pronostico_dr3_lrg2_v2_calibracion_{tag}.json"
    resp = HERE / f"pronostico_dr3_lrg2_v2_resultado_{tag}.json"
    assert calp.exists(), "primero calibrar"
    assert not resp.exists(), f"ya existe {resp.name}: Λ se calcula una sola vez"
    sha_rec = (HERE / f"pronostico_dr3_lrg2_v2_calibracion_{tag}.sha256").read_text().split()[0]
    assert sha256_file(calp) == sha_rec, "la calibración no coincide con su sha256 registrado"
    cal = json.load(open(calp))
    assert cal["cov_sha256"] == sha256_file(cov_path), "la covarianza no es la usada en la calibración"
    assert cal["script_sha256"] == sha256_file(__file__), "el script ha cambiado desde la calibración"
    res = {"mode": "evaluar", "tag": tag, "calibration_sha256": sha_rec, "values_file": str(val_path), "values_sha256": sha256_file(val_path), "status": cal["rupture_checks"]["status"], "combos": {}}
    if cal["rupture_checks"]["status"] != "ok":
        res["result"] = "modelo roto"
    else:
        C3 = read_dr3_cov(cov_path)
        j = json.load(open(val_path))
        assert list(j["labels"]) == LABELS
        d3 = D.copy()
        d3[:13] = np.array(j["mean"], float)
        for combo in COMBOS:
            c = cal["combos"][combo]
            _, p3 = P.fit_c(combo, KEEP, d3, C3)                      # full refit without LRG2 with the DR3 covariance
            J3, _, S3, _, _ = P.fisher(combo, KEEP, p3, C3)
            r3 = (d3 - JB.model_at(p3))[LRG2]
            V = C3[np.ix_(LRG2, LRG2)] + J3[LRG2] @ S3 @ J3[LRG2].T
            dens = {h: (np.array(c["exact_densities"][h]["a"]), np.array(c["exact_densities"][h]["Sigma"])) for h in ("H0", "H1")}
            L = float(lam_fun(dens)(r3))
            verdict = "favorece H1" if L > c["c0"] else ("favorece H0" if L < c["c1"] else "no concluyente")
            res["combos"][combo] = {"p3": p3.tolist(), "r3": r3.tolist(), "d2_DR3": float(r3 @ np.linalg.solve(V, r3)), "Lambda": L, "c0": c["c0"], "c1": c["c1"], "verdict": verdict,
                                    "strong": bool(abs(L) > LAMBDA_STRONG)}
        res["result"] = {c: res["combos"][c]["verdict"] for c in COMBOS}
    resp.write_text(json.dumps(res, indent=1, default=float), encoding="utf-8")
    print(json.dumps({k: res[k] for k in ("status", "result")} | {"combos": {c: {k: v[k] for k in ("r3", "d2_DR3", "Lambda", "verdict")} for c, v in res["combos"].items()}}, indent=1))
    return resp


def falso():
    """Fake DR3 files from DR2 (test of the evaluation mode; never real data): cov = 0.6 C (ok), cov = 0.3 C (rho2 out of range), values = DR2."""
    for tag, f in (("falso", 0.6), ("falso_roto", 0.3)):
        p = HERE / f"pronostico_dr3_lrg2_v2_{tag}_dr3_cov.json"
        assert not p.exists()
        p.write_text(json.dumps({"note": f"DATOS FALSOS: covarianza de DR2 × {f}; solo para probar el modo evaluación", "labels": LABELS, "cov": (f * COV13).tolist()}, indent=1), encoding="utf-8")
    p = HERE / "pronostico_dr3_lrg2_v2_falso_dr3_valores.json"
    assert not p.exists()
    p.write_text(json.dumps({"note": "DATOS FALSOS: valores de DR2 sin cambio; solo para probar el modo evaluación", "labels": LABELS, "mean": D[:13].tolist()}, indent=1), encoding="utf-8")
    print("escritos ficheros falsos")


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "registro"
    if mode == "registro":
        registration()
    elif mode == "falso":
        falso()
    elif mode == "calibrar":
        calibrar(Path(sys.argv[2]), sys.argv[3])
    elif mode == "evaluar":
        evaluar(Path(sys.argv[2]), Path(sys.argv[3]), sys.argv[4])
    else:
        raise SystemExit("modo desconocido")
