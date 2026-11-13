#!/usr/bin/env python3
"""
PRONOSTICO_DR3_REGISTRO_v3.md: version 3 of the DR3 forecast of the LRG2 residual, with the corrections of the second review
(REVISION2_PRONOSTICO_DR3.md, commit c7df2d8). Imports the set-up of pronostico_dr3_lrg2_v2.py / pronostico_dr3_lrg2.py (not modified).
Unchanged: P18, P19, H0, H1, the nesting idea, the decision rule (c0 = pct 95 of Lambda | H0, c1 = pct 5 of Lambda | H1, three outcomes).
New:
  a. Matrix nesting for the evaluation: n3 = R n2 + e, R = C3 C2^-1 (block per bin), e independent of n2, Cov(e) = C3 - C3 C2^-1 C3, with the
     DR3 covariance C3 as published. Then Cov(n3, n2) = R C2 = C3 and Cov(n3) = R C2 R^T + Cov(e) = C3 C2^-1 C3 + C3 - C3 C2^-1 C3 = C3.
     With C3 = rho2 C2: R = rho2 I, Cov(e) = rho2 (1 - rho2) C2 (the previous case). Checks per bin (eigenvalues recorded; tolerance: in a
     matrix that must be PSD, lambda < -1e-10 lambda_max is a real failure, -1e-10 lambda_max <= lambda < 0 is rounding): C2, C3 symmetric and
     positive definite (Cholesky), C2 - C3 PSD, Cov(e) PSD. Any real failure -> "modelo roto".
  b. Exact densities of r3 under H0 / H1 with that C3 (affine map r3 = a + M z). For C3 = rho2 C2 they must equal those of v2 (checked).
  c. High-precision thresholds: c0, c1 from 1e7 samples of N(a, Sigma) under H0 and H1 (seeds recorded), MC error from 10 blocks; same method
     in the main cases and in the rho2 sweep.
  d. Range: rho2_LRG2,eff = tr(C3_LRG2) / tr(C2_LRG2) must be in [0.40, 0.80]. The 20 % / 0.1 tolerances of v2 are removed.
  e. Rule b with tolerance: inter-bin correlation coefficients of DR3 with |r| >= 1e-6 -> "modelo roto".
  f. Identification of bins in DR3 by label: same kind and tracer, z_eff within +-0.01 of DR2's (LRG2: LRG with z in 0.706 +- 0.01).
  g. Single principal output: combination "sinSN"; "Pantheon+" is computed and reported as secondary.
  h. Evaluation mode: `calibrar COV.json TAG` reads only the covariance; `evaluar COV.json VALORES.json TAG` reads the values once. The sha256
     checks of the script and of the calibration are internal consistency controls, not the commitment (the commitment is the public commit order).
Linear algebra: Cholesky factorisations and linear solves; no explicit inverses (the only inverse-like operation is a Cholesky solve).

    ../venv/bin/python -I pronostico_dr3_lrg2_v3.py                        # registration run -> pronostico_dr3_lrg2_v3.json (+ .sha256 of the script)
    ../venv/bin/python -I pronostico_dr3_lrg2_v3.py falso                  # fake DR3 files (DATOS FALSOS)
    ../venv/bin/python -I pronostico_dr3_lrg2_v3.py calibrar COV.json TAG
    ../venv/bin/python -I pronostico_dr3_lrg2_v3.py evaluar COV.json VALORES.json TAG
DR3 input files (JSON): {"labels": [...], "cov": NxN} and {"labels": [...], "mean": [N]}, labels "KIND_TRACER_zZ" as in DR2.
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
from scipy.linalg import cho_factor, cho_solve, solve_triangular

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import pronostico_dr3_lrg2_v2 as V2  # noqa: E402

P, JB = V2.P, V2.JB
D, COV13, LRG2, KEEP, C22 = V2.D, V2.COV13, V2.LRG2, V2.KEEP, V2.C22
BINS, BIN_NAMES, LABELS = V2.BINS, V2.BIN_NAMES, V2.LABELS
SEED_BASE = V2.SEED_BASE
RHO2S = V2.RHO2S
RHO2_GRID = V2.RHO2_GRID
COMBOS = ("sinSN", "Pantheon+")
PRINCIPAL = "sinSN"
N_HP = 10_000_000
N_BLOCKS = 10
CHUNK = 1_000_000
EIG_TOL = 1e-10
CORR_TOL = 1e-6
Z_TOL = 0.01
RHO2_RANGE = (0.40, 0.80)
LAMBDA_STRONG = 3.0
SEED_HP = 7          # [SEED_BASE, i_comb, case_index, i_hyp, SEED_HP]; case_index = i_rho2 (main), 100 + k (sweep), 500 (evaluation)


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


def parse_label(lab):
    kind, rest = lab.split("_", 1)
    tracer, z = rest.rsplit("_z", 1)
    return kind, tracer, float(z)


# ----------------------------------------------------------------------------------------------------------------------
# matrix nesting: checks, R = C3 C2^-1, factor of Cov(e)
# ----------------------------------------------------------------------------------------------------------------------
def psd_check(M, name):
    Ms = 0.5 * (M + M.T)
    ev = np.linalg.eigvalsh(Ms)
    lmax = float(max(ev.max(), 0.0))
    thr = -EIG_TOL * lmax
    real_neg = [float(v) for v in ev if v < thr]
    rounding = [float(v) for v in ev if thr <= v < 0]
    return {"matrix": name, "eigenvalues": ev.tolist(), "max_abs_asym": float(np.max(np.abs(M - M.T))), "tolerance_neg": thr, "real_negative": real_neg, "rounding_negative": rounding,
            "ok": len(real_neg) == 0}


def pd_check(M, name):
    rep = psd_check(M, name)
    try:
        cho_factor(0.5 * (M + M.T), lower=True)
        rep["cholesky_ok"] = True
    except np.linalg.LinAlgError:
        rep["cholesky_ok"] = False
    rep["ok"] = bool(rep["ok"] and rep["cholesky_ok"] and len(rep["rounding_negative"]) == 0)
    return rep


def nesting(C2, C3):
    """Per-bin checks and operators: R (13x13 block), Fe (13x13 factor of Cov(e), Cov(e) = Fe Fe^T). Returns report, R, Fe."""
    rep = {"bins": {}, "status": "ok", "broken": []}
    R = np.zeros((13, 13))
    Fe = np.zeros((13, 13))
    for name, rows in BINS.items():
        ix = np.ix_(rows, rows)
        A2, A3 = C2[ix], C3[ix]
        b = {"C2": pd_check(A2, "C2"), "C3": pd_check(A3, "C3"), "C2_minus_C3": psd_check(A2 - A3, "C2 - C3")}
        if b["C2"]["ok"] and b["C3"]["ok"]:
            f2 = cho_factor(0.5 * (A2 + A2.T), lower=True)
            Rb = cho_solve(f2, A3).T                       # (C2^-1 C3)^T = C3 C2^-1 (C3 symmetric)
            Ce = A3 - Rb @ A3                              # C3 - C3 C2^-1 C3
            Ce = 0.5 * (Ce + Ce.T)
            b["Cov_e"] = psd_check(Ce, "Cov(e) = C3 - C3 C2^-1 C3")
            w, U = np.linalg.eigh(Ce)
            w = np.where(w < 0, 0.0, w)                    # rounding negatives -> 0 (real negatives already flagged)
            R[ix] = Rb
            Fe[ix] = U * np.sqrt(w)[None, :]
        else:
            b["Cov_e"] = {"ok": False, "note": "C2 or C3 not positive definite"}
        b["rho2_eff_trace"] = float(np.trace(A3) / np.trace(A2))
        bad = [k for k in ("C2", "C3", "C2_minus_C3", "Cov_e") if not b[k]["ok"]]
        if bad:
            rep["broken"].append(f"{name}: " + ", ".join(bad))
        rep["bins"][name] = b
    if rep["broken"]:
        rep["status"] = "modelo roto"
    return rep, R, Fe


def interbin_check(C3):
    s = np.sqrt(np.diag(C3))
    corr = C3 / np.outer(s, s)
    off = np.ones((13, 13), bool)
    for rows in BINS.values():
        for i in rows:
            for j in rows:
                off[i, j] = False
    m = float(np.max(np.abs(corr[off])))
    return {"max_abs_interbin_corr": m, "tolerance": CORR_TOL, "ok": m < CORR_TOL}


# ----------------------------------------------------------------------------------------------------------------------
# DR2 set-up and Fisher with Cholesky
# ----------------------------------------------------------------------------------------------------------------------
def fisher_c(combo, keep, p, cov13):
    """J (NROWS x npar), Cholesky factor of F = J^T C^-1 J on rows [keep, 13, 14, SN], C^-1 J (via solves)."""
    J = JB.build_B_at(p, P.npar_of(combo))
    rows = np.array(P.rows_of(combo, keep))
    CJ = P.cinv_apply(combo, keep, cov13, J[rows])
    F = J[rows].T @ CJ
    return J, cho_factor(0.5 * (F + F.T), lower=True), CJ, rows


def setup(combo):
    npar = P.npar_of(combo)
    chi2_full, p_full = P.PB.fit_d(combo, list(range(13)), D)
    _, Ff, _, _ = fisher_c(combo, list(range(13)), p_full, COV13)
    chi2_2, p2 = P.PB.fit_d(combo, KEEP, D)
    J, F2f, _, rows2 = fisher_c(combo, KEEP, p2, COV13)
    m2 = JB.model_at(p2)
    r2 = (D - m2)[LRG2]
    J2 = J[LRG2]
    Q = J2 @ cho_solve(F2f, J2.T)
    # factors of the posterior covariances without explicit inverses: Sigma = F^-1 = L^-T L^-1 -> p = p_hat + L^-T z
    return dict(combo=combo, npar=npar, p_full=p_full, chi2_full=chi2_full, LFfull=Ff[0], p2=p2, chi2_2=chi2_2, LF2=F2f[0], J=J, J2=J2, m2=m2, r2=r2, Q=Q,
                V2=C22 + Q, m2_bao=m2[:13], m2_rest=m2[rows2[len(KEEP):]], d_rest=D[rows2[len(KEEP):]], L22=np.linalg.cholesky(C22))


def case_ops(st, C3):
    rep, R, Fe = nesting(COV13, C3)
    if rep["status"] != "ok":
        return rep, None
    J, F3f, CJ3, _ = fisher_c(st["combo"], KEEP, st["p2"], C3)
    G3 = cho_solve(F3f, CJ3.T)                             # F3^-1 J_(-2)^T C3^-1
    const_rest = G3[:, len(KEEP):] @ (st["d_rest"] - st["m2_rest"])
    V3k = C3[np.ix_(LRG2, LRG2)] + st["J2"] @ cho_solve(F3f, st["J2"].T)
    V3k = 0.5 * (V3k + V3k.T)
    return rep, dict(C3=C3, R=R, Fe=Fe, G3_bao=G3[:, :len(KEEP)], const_rest=const_rest, V3k=V3k, V3kf=cho_factor(V3k, lower=True), LF3=F3f[0])


def p_draw(st, hyp, z_p):
    """p_true = p_hat + L^-T z (Sigma = F^-1 = L^-T L^-1)."""
    L = st["LFfull"] if hyp == "H0" else st["LF2"]
    base = st["p_full"] if hyp == "H0" else st["p2"]
    return base[: st["npar"]][None, :] + solve_triangular(L.T, z_p.T, lower=False).T


def gen_r3(hyp, st, op, z_p, z_e, z_mu):
    """Linearised generation (matrix nesting) + linearised refit, vectorised over universes."""
    p_true = p_draw(st, hyp, z_p)
    m_true = st["m2_bao"][None, :] + (p_true - st["p2"][: st["npar"]]) @ st["J"][:13].T
    x = m_true.copy()
    if hyp == "H1":
        mu = D[LRG2][None, :] - m_true[:, LRG2] + z_mu @ st["L22"].T
        x[:, LRG2] += mu
    n2 = D[:13][None, :] - x
    d3 = x + n2 @ op["R"].T + z_e @ op["Fe"].T               # d3 = m_true (+ E2 mu) + R n2 + e
    dp3 = (d3[:, KEEP] - st["m2_bao"][KEEP][None, :]) @ op["G3_bao"].T + op["const_rest"][None, :]
    return d3[:, LRG2] - st["m2_bao"][LRG2][None, :] - dp3 @ st["J2"].T


def exact_density(hyp, st, op):
    npar = st["npar"]
    nz = npar + 13 + (2 if hyp == "H1" else 0)

    def r3_of(z):
        z = np.atleast_2d(z)
        return gen_r3(hyp, st, op, z[:, :npar], z[:, npar:npar + 13], z[:, npar + 13:] if hyp == "H1" else None)[0]

    a = r3_of(np.zeros(nz))
    M = np.column_stack([r3_of(np.eye(nz)[k]) - a for k in range(nz)])
    lin = max(float(np.max(np.abs(r3_of(2 * np.eye(nz)[k]) - a - 2 * M[:, k]))) for k in range(nz))
    Sig = M @ M.T
    Sig = 0.5 * (Sig + Sig.T)
    return a, Sig, lin


class Gauss:
    def __init__(self, a, Sig):
        self.a, self.Sig = a, Sig
        self.L = np.linalg.cholesky(Sig)
        self.logdet = 2 * np.sum(np.log(np.diag(self.L)))

    def logpdf(self, r):
        y = solve_triangular(self.L, (r - self.a).T, lower=True)
        return -0.5 * (np.sum(y * y, axis=0) + self.logdet + 2 * np.log(2 * np.pi))

    def sample(self, n, rng):
        return self.a + rng.standard_normal((n, 2)) @ self.L.T


def hp_thresholds(dens, V3kf, seeds):
    """c0, c1, power... from N_HP exact samples per hypothesis, errors from N_BLOCKS blocks; d2 quantiles too."""
    lam, d2 = {}, {}
    for hyp in ("H0", "H1"):
        rng = np.random.default_rng(seeds[hyp])
        L, Q = [], []
        for _ in range(N_HP // CHUNK):
            r = dens[hyp].sample(CHUNK, rng)
            L.append(dens["H1"].logpdf(r) - dens["H0"].logpdf(r))
            y = cho_solve(V3kf, r.T)
            Q.append(np.sum(r.T * y, axis=0))
        lam[hyp], d2[hyp] = np.concatenate(L), np.concatenate(Q)

    def stats(l0, l1, q0, q1):
        c0, c1 = float(np.percentile(l0, 95)), float(np.percentile(l1, 5))
        return {"c0": c0, "c1": c1, "c1_lt_c0": bool(c1 < c0), "power": float(np.mean(l1 > c0)), "P_lt_c1_H0": float(np.mean(l0 < c1)),
                "inconclusive_H0": float(np.mean((l0 >= c1) & (l0 <= c0))), "inconclusive_H1": float(np.mean((l1 >= c1) & (l1 <= c0))),
                "P_gt3_H0": float(np.mean(l0 > LAMBDA_STRONG)), "P_gt3_H1": float(np.mean(l1 > LAMBDA_STRONG)), "P_ltm3_H0": float(np.mean(l0 < -LAMBDA_STRONG)), "P_ltm3_H1": float(np.mean(l1 < -LAMBDA_STRONG)),
                "med_d2_H0": float(np.median(q0)), "med_d2_H1": float(np.median(q1))}

    full = stats(lam["H0"], lam["H1"], d2["H0"], d2["H1"])
    nb = N_HP // N_BLOCKS
    blocks = [stats(lam["H0"][k * nb:(k + 1) * nb], lam["H1"][k * nb:(k + 1) * nb], d2["H0"][k * nb:(k + 1) * nb], d2["H1"][k * nb:(k + 1) * nb]) for k in range(N_BLOCKS)]
    err = {k: float(np.std([b[k] for b in blocks], ddof=1) / np.sqrt(N_BLOCKS)) for k in full if k != "c1_lt_c0"}
    full["mc_err"] = err
    full["Lambda_pct_5_50_95"] = {h: np.percentile(lam[h], [5, 50, 95]).tolist() for h in lam}
    full["d2_pct_5_50_95"] = {h: np.percentile(d2[h], [5, 50, 95]).tolist() for h in d2}
    full["N"], full["n_blocks"], full["seeds"] = N_HP, N_BLOCKS, seeds
    return full


def calibrate(st, C3, seeds):
    rep, op = case_ops(st, C3)
    if op is None:
        return rep, None, None
    dens, exact = {}, {}
    for hyp in ("H0", "H1"):
        a, Sig, lin = exact_density(hyp, st, op)
        dens[hyp] = Gauss(a, Sig)
        exact[hyp] = {"a": a.tolist(), "Sigma": Sig.tolist(), "linearity_check_max_abs": lin, "E_d2": float(a @ cho_solve(op["V3kf"], a) + np.trace(cho_solve(op["V3kf"], Sig)))}
    hp = hp_thresholds(dens, op["V3kf"], seeds)
    return rep, {"nesting": rep, "V3k": op["V3k"].tolist(), "exact": exact, "hp": hp}, (op, dens)


# ----------------------------------------------------------------------------------------------------------------------
# registration run
# ----------------------------------------------------------------------------------------------------------------------
def registration():
    t0 = time.time()
    v2 = json.load(open(HERE / "pronostico_dr3_lrg2_v2.json"))
    out = {"version": "v3 (cambios de REVISION2_PRONOSTICO_DR3.md)", "script_sha256": sha256_file(__file__), "git_head_at_run": git_head(), "principal_combo": PRINCIPAL,
           "N_hp": N_HP, "n_blocks": N_BLOCKS, "eig_tol": EIG_TOL, "corr_tol": CORR_TOL, "z_tol": Z_TOL, "rho2_LRG2_range": RHO2_RANGE,
           "seeds": {"hp": f"[SEED_BASE, i_comb, case_index, i_hyp, {SEED_HP}]: case_index = i_rho2 (0.6, 0.5, 0.72), 100 + k (sweep), 500 (DR3 evaluation)"},
           "nesting_derivation": "n3 = R n2 + e, R = C3 C2^-1 (block per bin), e ⊥ n2, Cov(e) = C3 - C3 C2^-1 C3. Cov(n3, n2) = R C2 = C3; "
                                 "Cov(n3) = R C2 R^T + Cov(e) = C3 C2^-1 C3 + C3 - C3 C2^-1 C3 = C3. With C3 = rho2 C2: R = rho2 I, Cov(e) = rho2 (1 - rho2) C2.",
           "rho2_grid": RHO2_GRID, "combos": {}}
    for ic, combo in enumerate(COMBOS):
        st = setup(combo)
        o2 = v2["combos"][combo]["DR2"]
        rec = {"DR2": {"p_full": st["p_full"].tolist(), "p2": st["p2"].tolist(), "r2": st["r2"].tolist(), "Q": st["Q"].tolist(), "V2": st["V2"].tolist(),
                       "check_vs_v2_max_abs": float(max(np.max(np.abs(st["r2"] - np.array(o2["r2"]))), np.max(np.abs(st["p2"] - np.array(o2["p2"]))), np.max(np.abs(st["Q"] - np.array(o2["Q"])))))},
               "cases": {}, "rho2_sweep": {}}
        for ir, rho2 in enumerate(RHO2S):
            seeds = {hyp: [SEED_BASE, ic, ir, ih, SEED_HP] for ih, hyp in enumerate(("H0", "H1"))}
            rep, case, _ = calibrate(st, rho2 * COV13, seeds)
            assert case is not None, "PARA: anidamiento roto con C3 = rho2 C2"
            vc = v2["combos"][combo]["cases"][f"rho2_{rho2:g}"]
            case["vs_v2"] = {"max_abs_diff_a": max(float(np.max(np.abs(np.array(case["exact"][h]["a"]) - np.array(vc["exact"][h]["a"])))) for h in ("H0", "H1")),
                             "max_abs_diff_Sigma": max(float(np.max(np.abs(np.array(case["exact"][h]["Sigma"]) - np.array(vc["exact"][h]["Sigma"])))) for h in ("H0", "H1")),
                             "max_abs_diff_V3k": float(np.max(np.abs(np.array(case["V3k"]) - np.array(vc["V3k"])))),
                             "v2_B_half": {k: vc["decision"][k] for k in ("c0", "c1", "power", "med_d2_H0", "med_d2_H1")}, "v2_mc_err": {k: vc["mc_error"]["stats"][k]["std"] for k in ("c0", "c1", "power", "med_d2_H0", "med_d2_H1")}}
            case["vs_v2"]["deviation_in_v2_sigma"] = {k: float((case["hp"][k] - vc["decision"][k]) / vc["mc_error"]["stats"][k]["std"]) for k in ("c0", "c1", "power", "med_d2_H0", "med_d2_H1")}
            rec["cases"][f"rho2_{rho2:g}"] = case
            h = case["hp"]
            print(f"{combo} rho2={rho2}: a/Sigma vs v2 {case['vs_v2']['max_abs_diff_a']:.1e}/{case['vs_v2']['max_abs_diff_Sigma']:.1e}; hp c0={h['c0']:.4f}±{h['mc_err']['c0']:.4f} c1={h['c1']:.4f}±{h['mc_err']['c1']:.4f} "
                  f"power={h['power']:.4f}±{h['mc_err']['power']:.4f} med d2 {h['med_d2_H0']:.3f}/{h['med_d2_H1']:.3f} (v2 B: {vc['decision']['c0']:.3f}/{vc['decision']['c1']:.3f}/{vc['decision']['power']:.4f})  [{time.time() - t0:.0f} s]", flush=True)
        for k, rho2 in enumerate(RHO2_GRID):
            seeds = {hyp: [SEED_BASE, ic, 100 + k, ih, SEED_HP] for ih, hyp in enumerate(("H0", "H1"))}
            _, case, _ = calibrate(st, rho2 * COV13, seeds)
            h = case["hp"]
            rec["rho2_sweep"][f"{rho2:.2f}"] = {k2: h[k2] for k2 in ("c0", "c1", "c1_lt_c0", "power", "P_lt_c1_H0", "inconclusive_H0", "inconclusive_H1", "med_d2_H0", "med_d2_H1")} | {"mc_err": h["mc_err"], "seeds": seeds}
        print(f"{combo} sweep: " + " ".join(f"{r}:{v['power']:.4f}" for r, v in rec["rho2_sweep"].items()) + f"  [{time.time() - t0:.0f} s]", flush=True)
        out["combos"][combo] = rec
    m = out["combos"][PRINCIPAL]["cases"]["rho2_0.6"]["hp"]
    out["predictions"] = {"P18": {"statement": "rho2 = 0.6, falsa alarma 5 %: potencia entre 0.5 y 0.75 (sin SN)", "power": m["power"], "mc_err": m["mc_err"]["power"], "holds": bool(0.5 <= m["power"] <= 0.75)},
                          "P19": {"statement": "rho2 = 0.6: mediana de d2_DR3 entre 3 y 6 bajo H0 y entre 9 y 15 bajo H1 (sin SN)", "median_H0": m["med_d2_H0"], "median_H1": m["med_d2_H1"],
                                  "holds": bool(3 <= m["med_d2_H0"] <= 6 and 9 <= m["med_d2_H1"] <= 15)},
                          "P23": {"statement": "Con los umbrales de alta precisión, la potencia del caso principal (rho2 = 0.6, sin SN) queda en 0.656 ± 0.010.", "power": m["power"], "mc_err": m["mc_err"]["power"],
                                  "c0": m["c0"], "c1": m["c1"], "holds": bool(abs(m["power"] - 0.656) <= 0.010)}}
    out["runtime_s"] = time.time() - t0
    (HERE / "pronostico_dr3_lrg2_v3.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")
    (HERE / "pronostico_dr3_lrg2_v3.sha256").write_text(f"{out['script_sha256']}  pronostico_dr3_lrg2_v3.py\n", encoding="utf-8")
    print(json.dumps(out["predictions"], indent=1))
    print("escrito pronostico_dr3_lrg2_v3.json", f"{out['runtime_s']:.0f} s")


# ----------------------------------------------------------------------------------------------------------------------
# DR3 evaluation mode
# ----------------------------------------------------------------------------------------------------------------------
def map_labels(labels3):
    """DR2 row -> DR3 row: same kind and tracer, |z3 - z2| <= Z_TOL, unique. Returns index array (13) and report, or None."""
    parsed3 = [parse_label(l) for l in labels3]
    idx, rep = [], {"matches": {}, "unmatched": [], "ambiguous": [], "z_shifts": {}, "ignored_dr3_rows": []}
    for lab in LABELS:
        k2, t2, z2 = parse_label(lab)
        cand = [i for i, (k3, t3, z3) in enumerate(parsed3) if k3 == k2 and t3 == t2 and abs(z3 - z2) <= Z_TOL]
        if len(cand) == 1:
            idx.append(cand[0])
            rep["matches"][lab] = labels3[cand[0]]
            if parsed3[cand[0]][2] != z2:
                rep["z_shifts"][lab] = parsed3[cand[0]][2] - z2
        elif not cand:
            rep["unmatched"].append(lab)
        else:
            rep["ambiguous"].append(lab)
    rep["ignored_dr3_rows"] = [labels3[i] for i in range(len(labels3)) if i not in idx]
    rep["ok"] = not rep["unmatched"] and not rep["ambiguous"]
    return (np.array(idx) if rep["ok"] else None), rep


def self_check():
    reg = HERE / "pronostico_dr3_lrg2_v3.sha256"
    cur = sha256_file(__file__)
    if reg.exists():
        return {"registered": reg.read_text().split()[0], "current": cur, "match": reg.read_text().split()[0] == cur}
    return {"registered": None, "current": cur, "match": None}


def calibrar(cov_path, tag):
    j = json.load(open(cov_path))
    idx, lrep = map_labels(list(j["labels"]))
    cal = {"mode": "calibrar", "tag": tag, "cov_file": str(cov_path), "cov_sha256": sha256_file(cov_path), "script_check": self_check(), "git_head_at_run": git_head(), "labels": lrep,
           "principal_combo": PRINCIPAL, "combos": {}}
    if idx is None:
        cal["status"] = "no evaluable"
        cal["reason"] = "bins de DR2 sin homólogo único en DR3 (mismo tipo y trazador, z_eff ± 0.01)"
    else:
        Cfull = np.array(j["cov"], float)
        C3 = Cfull[np.ix_(idx, idx)]
        cal["interbin"] = interbin_check(C3)
        rep, R, Fe = nesting(COV13, C3)
        cal["nesting"] = rep
        rho2_eff = {n: rep["bins"][n]["rho2_eff_trace"] for n in BIN_NAMES}
        cal["rho2_eff_trace"] = rho2_eff
        broken = list(rep["broken"])
        if not cal["interbin"]["ok"]:
            broken.append(f"correlación entre bins |r| máx {cal['interbin']['max_abs_interbin_corr']:.2e} ≥ {CORR_TOL}")
        if not (RHO2_RANGE[0] <= rho2_eff["LRG2"] <= RHO2_RANGE[1]):
            broken.append(f"rho2_LRG2,eff = {rho2_eff['LRG2']:.3f} fuera de [{RHO2_RANGE[0]}, {RHO2_RANGE[1]}]")
        cal["broken"] = broken
        cal["status"] = "modelo roto" if broken else "ok"
        if cal["status"] == "ok":
            for ic, combo in enumerate(COMBOS):
                st = setup(combo)
                seeds = {hyp: [SEED_BASE, ic, 500, ih, SEED_HP] for ih, hyp in enumerate(("H0", "H1"))}
                _, case, _ = calibrate(st, C3, seeds)
                cal["combos"][combo] = {"seeds": seeds, "exact_densities": case["exact"], "V3k": case["V3k"], "c0": case["hp"]["c0"], "c1": case["hp"]["c1"], "c1_lt_c0": case["hp"]["c1_lt_c0"],
                                        "power": case["hp"]["power"], "P_lt_c1_H0": case["hp"]["P_lt_c1_H0"], "mc_err": case["hp"]["mc_err"], "p2": st["p2"].tolist(), "r2": st["r2"].tolist()}
    outp = HERE / f"pronostico_dr3_lrg2_v3_calibracion_{tag}.json"
    assert not outp.exists(), f"ya existe {outp.name}"
    outp.write_text(json.dumps(cal, indent=1, default=float), encoding="utf-8")
    sha = sha256_file(outp)
    (HERE / f"pronostico_dr3_lrg2_v3_calibracion_{tag}.sha256").write_text(f"{sha}  {outp.name}\n", encoding="utf-8")
    summary = {"status": cal["status"], "broken": cal.get("broken", cal.get("reason")), "rho2_eff_trace": cal.get("rho2_eff_trace"), "script_check": cal["script_check"],
               "calibracion": {c: {k: cal["combos"][c][k] for k in ("c0", "c1", "c1_lt_c0", "power")} for c in cal["combos"]}, "sha256": sha}
    print(json.dumps(summary, indent=1))
    return outp


def evaluar(cov_path, val_path, tag):
    calp = HERE / f"pronostico_dr3_lrg2_v3_calibracion_{tag}.json"
    resp = HERE / f"pronostico_dr3_lrg2_v3_resultado_{tag}.json"
    assert calp.exists(), "primero calibrar"
    assert not resp.exists(), f"ya existe {resp.name}: Λ se calcula una sola vez (control de coherencia; el compromiso es el orden de los commits públicos)"
    sha_rec = (HERE / f"pronostico_dr3_lrg2_v3_calibracion_{tag}.sha256").read_text().split()[0]
    assert sha256_file(calp) == sha_rec, "la calibración no coincide con su sha256 registrado"
    cal = json.load(open(calp))
    assert cal["cov_sha256"] == sha256_file(cov_path), "la covarianza no es la usada en la calibración"
    sc = self_check()
    assert sc["match"] is not False, "el script no coincide con el sha256 registrado"
    res = {"mode": "evaluar", "tag": tag, "calibration_sha256": sha_rec, "values_file": str(val_path), "values_sha256": sha256_file(val_path), "script_check": sc, "status": cal["status"],
           "principal_combo": PRINCIPAL, "combos": {}}
    if cal["status"] != "ok":
        res["result_principal"] = cal["status"]
    else:
        j = json.load(open(cov_path))
        idx, _ = map_labels(list(j["labels"]))
        C3 = np.array(j["cov"], float)[np.ix_(idx, idx)]
        v = json.load(open(val_path))
        idxv, lrep = map_labels(list(v["labels"]))
        assert idxv is not None and list(idxv) == list(idx), "las etiquetas del fichero de valores no coinciden con las de la covarianza"
        d3 = D.copy()
        d3[:13] = np.array(v["mean"], float)[idxv]
        for combo in COMBOS:
            c = cal["combos"][combo]
            _, p3 = P.fit_c(combo, KEEP, d3, C3)                      # full refit without LRG2 with the DR3 covariance
            J3, F3f, _, _ = fisher_c(combo, KEEP, p3, C3)
            r3 = (d3 - JB.model_at(p3))[LRG2]
            V = C3[np.ix_(LRG2, LRG2)] + J3[LRG2] @ cho_solve(F3f, J3[LRG2].T)
            dens = {h: Gauss(np.array(c["exact_densities"][h]["a"]), np.array(c["exact_densities"][h]["Sigma"])) for h in ("H0", "H1")}
            L = float(dens["H1"].logpdf(r3[None, :])[0] - dens["H0"].logpdf(r3[None, :])[0])
            verdict = "favorece H1" if L > c["c0"] else ("favorece H0" if L < c["c1"] else "no concluyente")
            res["combos"][combo] = {"p3": p3.tolist(), "r3": r3.tolist(), "d2_DR3": float(r3 @ cho_solve(cho_factor(V, lower=True), r3)), "Lambda": L, "c0": c["c0"], "c1": c["c1"], "verdict": verdict,
                                    "strong_abs_Lambda_gt_3": bool(abs(L) > LAMBDA_STRONG)}
        res["result_principal"] = res["combos"][PRINCIPAL]["verdict"]
        res["result_secondary_Pantheon+"] = res["combos"]["Pantheon+"]["verdict"]
    resp.write_text(json.dumps(res, indent=1, default=float), encoding="utf-8")
    print(json.dumps({k: res[k] for k in res if k in ("status", "result_principal", "result_secondary_Pantheon+")} | {"combos": {c: {k: v[k] for k in ("r3", "d2_DR3", "Lambda", "verdict")} for c, v in res["combos"].items()}}, indent=1))
    return resp


def falso():
    """Fake DR3 files (DATOS FALSOS) for the test of the evaluation mode: (1) C3 = 0.6 C2; (2) one direction less precise than DR2 (LRG2 D_H
    variance x 2); (3) inter-bin correlation 1e-3 between LRG1 and LRG2; (4) inter-bin correlation 1e-9 (below tolerance); values = DR2."""
    C = 0.6 * COV13
    C2 = C.copy()
    C2[4, 4] = 2.0 * COV13[4, 4]
    C3 = C.copy()
    s = np.sqrt(np.diag(C3))
    for i, jj in ((1, 3), (2, 4)):
        C3[i, jj] = C3[jj, i] = 1e-3 * s[i] * s[jj]
    C4 = C.copy()
    C4[1, 3] = C4[3, 1] = 1e-9 * s[1] * s[3]
    for tag, M, note in (("falso", C, "0.6 × C_DR2"), ("falso_direccion", C2, "0.6 × C_DR2 con Var(D_H de LRG2) = 2 × la de DR2 (una dirección menos precisa)"),
                         ("falso_corr", C3, "0.6 × C_DR2 con correlación 1e-3 entre LRG1 y LRG2"), ("falso_corr_tol", C4, "0.6 × C_DR2 con correlación 1e-9 entre LRG1 y LRG2 (bajo la tolerancia)")):
        p = HERE / f"pronostico_dr3_lrg2_v3_{tag}_dr3_cov.json"
        assert not p.exists()
        p.write_text(json.dumps({"note": f"DATOS FALSOS: {note}; solo para probar el modo evaluación", "labels": LABELS, "cov": M.tolist()}, indent=1), encoding="utf-8")
    p = HERE / "pronostico_dr3_lrg2_v3_falso_dr3_valores.json"
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
