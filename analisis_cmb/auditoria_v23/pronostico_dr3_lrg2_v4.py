#!/usr/bin/env python3
"""
PRONOSTICO_DR3_REGISTRO_v4.md: version 4 of the DR3 forecast of the LRG2 residual. Correction limited to the treatment of the redshift
(imports pronostico_dr3_lrg2_v3.py, not modified). Unchanged: P18, P19, P23, H0, H1, the matrix nesting, the decision rule, thresholds.

Redshifts: the DR2 conditioning data keep their z_eff (the DR2 fits, the DR2 noise n2 = d2 - m(p, z_DR2) and the DR2 residual of LRG2 at
z = 0.706 are untouched). The DR3 predictions (model and Jacobian) are evaluated at the z_eff published by DR3 for each bin, so that a shift of
z_eff is absorbed by the prediction, not by the residual. Model at an arbitrary z: D_M by linear interpolation in u = ln(1+z) of the same grid
history (error < 1e-6 relative, as for the SN), D_H = 1/E(z) in closed form, D_V = (z D_M^2 D_H)^(1/3); at the DR2 nodes it reproduces the grid
model to 2e-16.
Generation (linearised around p2): n2 = d2 - m(p_true, z_DR2) [- E2 mu2 under H1]; d3 = m(p_true, z_DR3) [+ E2 mu3] + R n2 + e, with
mu3 = mu2 * m(p2, z_DR3)_LRG2 / m(p2, z_DR2)_LRG2 (H1, new assumption: the persistent deviation mu is FRACTIONAL and the same at z_eff DR2 and
z_eff DR3 within +-0.01; the scaling is linearised at p2). The DR3 refit without LRG2, V3,k and the exact densities use the model at z_DR3.
Residuals are FRACTIONAL: r3 = (d3 - m(p3, z_DR3)) / m(z_DR3) (denominator m(p2, z_DR3) in the linearised densities, m(p3, z_DR3) in the
evaluation; the difference is second order). Lambda, c0, c1, power and d2 are invariant under this rescaling, so with z_DR3 = z_DR2 v4 must
reproduce v3 to rounding (checked). If |z_eff DR3 - z_eff DR2| > 0.01 for LRG2: "no evaluable"; for the other six bins the shift only moves the
prediction and is recorded.
Evaluation mode, reading order: (1) z_eff (from the labels) and covariance of DR3 -> calibration with sha256; (2) only then the measured values.

    ../venv/bin/python -I pronostico_dr3_lrg2_v4.py                        # registration: v4 = v3 check (6 cases) + P24 table -> pronostico_dr3_lrg2_v4.json
    ../venv/bin/python -I pronostico_dr3_lrg2_v4.py falso                  # fake DR3 files (DATOS FALSOS): z_LRG2 = 0.712 and 0.720
    ../venv/bin/python -I pronostico_dr3_lrg2_v4.py calibrar COV.json TAG
    ../venv/bin/python -I pronostico_dr3_lrg2_v4.py evaluar COV.json VALORES.json TAG
Nothing existing is modified.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.linalg import cho_factor, cho_solve, solve_triangular
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import pronostico_dr3_lrg2_v3 as V3  # noqa: E402

P, JB = V3.P, V3.JB
S, X, G = P.S, P.X, P.G
SN, NSN, NROWS = P.SN, P.NSN, P.NROWS
D, COV13, LRG2, KEEP, C22 = V3.D, V3.COV13, V3.LRG2, V3.KEEP, V3.C22
BINS, BIN_NAMES, LABELS = V3.BINS, V3.BIN_NAMES, V3.LABELS
SEED_BASE, SEED_HP, RHO2S, COMBOS, PRINCIPAL = V3.SEED_BASE, V3.SEED_HP, V3.RHO2S, V3.COMBOS, V3.PRINCIPAL
Z_TOL, RHO2_RANGE, LAMBDA_STRONG = V3.Z_TOL, V3.RHO2_RANGE, V3.LAMBDA_STRONG
Z2 = np.array(X.Z13, float)
SHIFTS = (-0.01, -0.005, 0.005, 0.01)
sha256_file, git_head, parse_label, map_labels, Gauss, hp_thresholds, nesting, interbin_check = V3.sha256_file, V3.git_head, V3.parse_label, V3.map_labels, V3.Gauss, V3.hp_thresholds, V3.nesting, V3.interbin_check


# ----------------------------------------------------------------------------------------------------------------------
# BAO model at arbitrary redshifts (same grid history; D_M interpolated in u; D_H closed form)
# ----------------------------------------------------------------------------------------------------------------------
def model_bao(p, z):
    h, wm = p[0], p[1]
    om, orad = wm / h**2, S.OMR / h**2
    hist, m = S.history(G, om, orad, np.zeros(G.N))
    if hist is None:
        return None
    dm_u, _ = hist
    u = np.log1p(z)
    dm = np.interp(u, G.U, dm_u)
    inv = 1.0 / np.sqrt(om * (1 + z) ** 3 + orad * (1 + z) ** 4 + (1 - om - orad))
    d = np.where(X.IS_DM, dm, np.where(X.IS_DH, inv, np.cbrt(z * dm * dm * inv)))
    return S.C_KMS / (100.0 * h * float(S.E.RD_SPL(wm))) * d


def model_full_z(p, z):
    m = JB.model_at(p)
    m[:13] = model_bao(p, z)
    return m


def build_J_z(p, z, npar):
    """J (NROWS x npar): externals and SN rows from JB.build_B_at (independent of the BAO z); BAO rows by central differences of model_bao at z."""
    J = JB.build_B_at(p, npar)
    for jj, eps in enumerate((1e-4, 1e-5, 1e-3)[:npar]):
        pp, pm = np.array(p, float), np.array(p, float)
        pp[jj] += eps
        pm[jj] -= eps
        J[:13, jj] = (model_bao(pp, z) - model_bao(pm, z)) / (2 * eps)
    return J


def fisher_cz(combo, keep, p, cov13, z):
    J = build_J_z(p, z, P.npar_of(combo))
    rows = np.array(P.rows_of(combo, keep))
    CJ = P.cinv_apply(combo, keep, cov13, J[rows])
    F = J[rows].T @ CJ
    return J, cho_factor(0.5 * (F + F.T), lower=True), CJ, rows


def fit_cz(combo, keep, d, cov13, z):
    """As pronostico_dr3_lrg2.fit_c (same optimiser settings) with the BAO model evaluated at z (DR3 redshifts)."""
    L = np.linalg.cholesky(cov13[np.ix_(keep, keep)])
    kb = np.array(keep)
    y15 = d[:15]
    if combo == "sinSN":
        def resid(p):
            v = model_bao(p, z)
            if v is None:
                return np.full(len(kb) + 2, 1e3)
            r = y15[:13] - v
            return np.concatenate([np.linalg.solve(L, r[kb]), [(y15[13] - JB.model_at(np.array([p[0], p[1], 0.0]))[13]) / X.SIG_TH, (y15[14] - p[1]) / X.SIG_WM]])
        r = least_squares(resid, np.asarray((0.682, 0.1407)), bounds=X.BOUNDS_P, method="trf", x_scale=[0.01, 0.001], ftol=1e-10, xtol=1e-10, max_nfev=300)
        return float(2.0 * r.cost), np.array([r.x[0], r.x[1], 0.0])
    mb = d[15:]
    x0 = (0.682, 0.1407)
    hist, _ = S.history(G, x0[1] / x0[0] ** 2, S.OMR / x0[0] ** 2, np.zeros(G.N))
    M0 = float((SN["Ci1"] @ (mb - S.mu_sn(G, SN, hist[0]))) / SN["e"]) if hist is not None else 0.0

    def resid(p):
        m = model_full_z(p, z)
        if not np.all(np.isfinite(m)):
            return np.full(len(kb) + 2 + NSN, 1e3)
        r = np.concatenate([y15, mb]) - m
        return np.concatenate([solve_triangular(L, r[kb], lower=True, check_finite=False), [r[13] / X.SIG_TH, r[14] / X.SIG_WM],
                               solve_triangular(SN["L"], r[15:], lower=True, check_finite=False)])
    r = least_squares(resid, np.array([x0[0], x0[1], M0]), bounds=([0.4, 0.05, -50.0], [1.0, 0.30, 50.0]), method="trf",
                      x_scale=[0.01, 0.001, 0.01], ftol=1e-10, xtol=1e-10, max_nfev=200)
    return float(2.0 * r.cost), np.array(r.x)


# ----------------------------------------------------------------------------------------------------------------------
# linearised generation with DR2 conditioning redshifts and DR3 prediction redshifts
# ----------------------------------------------------------------------------------------------------------------------
def case_ops(st, C3, z3):
    rep, R, Fe = nesting(COV13, C3)
    if rep["status"] != "ok":
        return rep, None
    m2_bao3 = model_bao(st["p2"], z3)
    J3, F3f, CJ3, _ = fisher_cz(st["combo"], KEEP, st["p2"], C3, z3)
    G3 = cho_solve(F3f, CJ3.T)
    const_rest = G3[:, len(KEEP):] @ (st["d_rest"] - st["m2_rest"])
    s3 = 1.0 / m2_bao3[LRG2]                                          # fractional scaling of the LRG2 residual
    V3k_abs = C3[np.ix_(LRG2, LRG2)] + J3[LRG2] @ cho_solve(F3f, J3[LRG2].T)
    V3k = (s3[:, None] * V3k_abs * s3[None, :])
    V3k = 0.5 * (V3k + V3k.T)
    ratio = m2_bao3[LRG2] / st["m2_bao"][LRG2]                        # fractional mu: mu3 = ratio * mu2 (linearised at p2)
    return rep, dict(C3=C3, z3=z3, R=R, Fe=Fe, m2_bao3=m2_bao3, J3=J3, G3_bao=G3[:, :len(KEEP)], const_rest=const_rest, s3=s3, ratio=ratio, V3k=V3k, V3kf=cho_factor(V3k, lower=True))


def gen_r3(hyp, st, op, z_p, z_e, z_mu):
    p_true = V3.p_draw(st, hyp, z_p)
    dp = p_true - st["p2"][: st["npar"]]
    m_true2 = st["m2_bao"][None, :] + dp @ st["J"][:13].T              # DR2 conditioning redshifts
    m_true3 = op["m2_bao3"][None, :] + dp @ op["J3"][:13].T            # DR3 prediction redshifts
    x2, x3 = m_true2.copy(), m_true3.copy()
    if hyp == "H1":
        mu2 = D[LRG2][None, :] - m_true2[:, LRG2] + z_mu @ st["L22"].T
        x2[:, LRG2] += mu2
        x3[:, LRG2] += mu2 * op["ratio"][None, :]
    n2 = D[:13][None, :] - x2
    d3 = x3 + n2 @ op["R"].T + z_e @ op["Fe"].T
    dp3 = (d3[:, KEEP] - op["m2_bao3"][KEEP][None, :]) @ op["G3_bao"].T + op["const_rest"][None, :]
    r3_abs = d3[:, LRG2] - op["m2_bao3"][LRG2][None, :] - dp3 @ op["J3"][LRG2].T
    return r3_abs * op["s3"][None, :]


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
    return a, 0.5 * (Sig + Sig.T), lin


def calibrate(st, C3, z3, seeds):
    rep, op = case_ops(st, C3, z3)
    if op is None:
        return rep, None
    dens, exact = {}, {}
    for hyp in ("H0", "H1"):
        a, Sig, lin = exact_density(hyp, st, op)
        dens[hyp] = Gauss(a, Sig)
        exact[hyp] = {"a": a.tolist(), "Sigma": Sig.tolist(), "linearity_check_max_abs": lin, "E_d2": float(a @ cho_solve(op["V3kf"], a) + np.trace(cho_solve(op["V3kf"], Sig)))}
    hp = hp_thresholds(dens, op["V3kf"], seeds)
    return rep, {"nesting": rep, "z3": np.asarray(z3).tolist(), "z_shifts": (np.asarray(z3) - Z2).tolist(), "m2_bao3_LRG2": op["m2_bao3"][LRG2].tolist(), "V3k_fractional": op["V3k"].tolist(),
                 "exact_fractional": exact, "hp": hp}


# ----------------------------------------------------------------------------------------------------------------------
def registration():
    t0 = time.time()
    v3 = json.load(open(HERE / "pronostico_dr3_lrg2_v3.json"))
    out = {"version": "v4 (corrección del redshift: predicción de DR3 en z_eff de DR3, residuos fraccionales)", "script_sha256": sha256_file(__file__), "git_head_at_run": git_head(),
           "assumption_H1_fractional": "mu persistente fraccional, el mismo en z_eff DR2 y z_eff DR3 dentro de ±0.01; mu3 = mu2 · m(p2, z3)/m(p2, z2) (linealizado en p2)",
           "shifts_P24": SHIFTS, "seeds": {"reproduction": "[SEED_BASE, i_comb, i_rho2, i_hyp, 7] (as v3)", "P24": "[SEED_BASE, 0, 600 + i_shift, i_hyp, 7]"}, "combos": {}, "P24": {}}
    for ic, combo in enumerate(COMBOS):
        st = V3.setup(combo)
        rec = {"cases": {}}
        for ir, rho2 in enumerate(RHO2S):
            seeds = {hyp: [SEED_BASE, ic, ir, ih, SEED_HP] for ih, hyp in enumerate(("H0", "H1"))}
            _, case = calibrate(st, rho2 * COV13, Z2, seeds)
            assert case is not None
            vc = v3["combos"][combo]["cases"][f"rho2_{rho2:g}"]
            s3 = 1.0 / np.array(case["m2_bao3_LRG2"])
            da = max(float(np.max(np.abs(np.array(case["exact_fractional"][h]["a"]) / s3 - np.array(vc["exact"][h]["a"])))) for h in ("H0", "H1"))
            dS = max(float(np.max(np.abs(np.array(case["exact_fractional"][h]["Sigma"]) / np.outer(s3, s3) - np.array(vc["exact"][h]["Sigma"])))) for h in ("H0", "H1"))
            dh = {k: case["hp"][k] - vc["hp"][k] for k in ("c0", "c1", "power", "P_lt_c1_H0", "med_d2_H0", "med_d2_H1")}
            case["vs_v3"] = {"max_abs_diff_a_absolute_units": da, "max_abs_diff_Sigma_absolute_units": dS, "diff_hp": dh, "max_abs_diff_hp": float(max(abs(v) for v in dh.values()))}
            rec["cases"][f"rho2_{rho2:g}"] = case
            print(f"{combo} rho2={rho2} dz=0: a/Sigma vs v3 {da:.1e}/{dS:.1e}; c0/c1/power diff vs v3 {dh['c0']:.1e}/{dh['c1']:.1e}/{dh['power']:.1e}  [{time.time() - t0:.0f} s]", flush=True)
        out["combos"][combo] = rec
    # P24: LRG2 z_eff shifted, principal case (sinSN, rho2 = 0.6)
    st = V3.setup(PRINCIPAL)
    base = out["combos"][PRINCIPAL]["cases"]["rho2_0.6"]["hp"]
    for i, dz in enumerate(SHIFTS):
        z3 = Z2.copy()
        z3[LRG2] += dz
        seeds = {hyp: [SEED_BASE, 0, 600 + i, ih, SEED_HP] for ih, hyp in enumerate(("H0", "H1"))}
        _, case = calibrate(st, 0.6 * COV13, z3, seeds)
        h = case["hp"]
        out["P24"][f"dz_{dz:+g}"] = {"z_LRG2": float(Z2[3] + dz), "seeds": seeds, "c0": h["c0"], "c1": h["c1"], "power": h["power"], "mc_err": {k: h["mc_err"][k] for k in ("c0", "c1", "power")},
                                     "d_c0": h["c0"] - base["c0"], "d_c1": h["c1"] - base["c1"], "d_power": h["power"] - base["power"], "med_d2_H0": h["med_d2_H0"], "med_d2_H1": h["med_d2_H1"],
                                     "m2_bao3_LRG2": case["m2_bao3_LRG2"], "a_fractional": {hh: case["exact_fractional"][hh]["a"] for hh in ("H0", "H1")}}
        print(f"P24 dz={dz:+g}: c0 {h['c0']:.4f} ({h['c0'] - base['c0']:+.4f}) c1 {h['c1']:.4f} ({h['c1'] - base['c1']:+.4f}) power {h['power']:.4f} ({h['power'] - base['power']:+.4f})  [{time.time() - t0:.0f} s]", flush=True)
    out["P24"]["dz_0"] = {"c0": base["c0"], "c1": base["c1"], "power": base["power"], "mc_err": {k: base["mc_err"][k] for k in ("c0", "c1", "power")}}
    out["predictions"] = {"P24": {"statement": "Desplazar el z_eff de LRG2 en ±0.01 (recalculando el modelo y el jacobiano en el nuevo z) cambia c0 y c1 menos de 0.05 y la potencia del caso principal menos de 0.02.",
                                  "max_abs_d_c0": max(abs(out["P24"][f"dz_{dz:+g}"]["d_c0"]) for dz in SHIFTS), "max_abs_d_c1": max(abs(out["P24"][f"dz_{dz:+g}"]["d_c1"]) for dz in SHIFTS),
                                  "max_abs_d_power": max(abs(out["P24"][f"dz_{dz:+g}"]["d_power"]) for dz in SHIFTS)}}
    out["predictions"]["P24"]["holds"] = bool(out["predictions"]["P24"]["max_abs_d_c0"] < 0.05 and out["predictions"]["P24"]["max_abs_d_c1"] < 0.05 and out["predictions"]["P24"]["max_abs_d_power"] < 0.02)
    out["runtime_s"] = time.time() - t0
    (HERE / "pronostico_dr3_lrg2_v4.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")
    (HERE / "pronostico_dr3_lrg2_v4.sha256").write_text(f"{out['script_sha256']}  pronostico_dr3_lrg2_v4.py\n", encoding="utf-8")
    print(json.dumps(out["predictions"], indent=1))
    print("escrito pronostico_dr3_lrg2_v4.json", f"{out['runtime_s']:.0f} s")


# ----------------------------------------------------------------------------------------------------------------------
# DR3 evaluation mode: (1) z_eff + covariance -> calibration; (2) values
# ----------------------------------------------------------------------------------------------------------------------
def self_check():
    reg = HERE / "pronostico_dr3_lrg2_v4.sha256"
    cur = sha256_file(__file__)
    return {"registered": reg.read_text().split()[0] if reg.exists() else None, "current": cur, "match": (reg.read_text().split()[0] == cur) if reg.exists() else None}


def calibrar(cov_path, tag):
    j = json.load(open(cov_path))
    labels3 = list(j["labels"])                       # only the labels (z_eff) and the covariance are read here, never "mean"
    idx, lrep = map_labels(labels3)
    cal = {"mode": "calibrar", "tag": tag, "cov_file": str(cov_path), "cov_sha256": sha256_file(cov_path), "script_check": self_check(), "git_head_at_run": git_head(), "labels": lrep,
           "keys_read": ["labels", "cov"], "principal_combo": PRINCIPAL, "combos": {}}
    if idx is None:
        cal["status"] = "no evaluable"
        cal["reason"] = "bins de DR2 sin homólogo único en DR3 (mismo tipo y trazador, z_eff ± 0.01); LRG2 fuera de 0.706 ± 0.01 cae aquí"
    else:
        z3 = np.array([parse_label(labels3[i])[2] for i in idx])
        C3 = np.array(j["cov"], float)[np.ix_(idx, idx)]
        cal["z3"], cal["z_shifts_per_row"] = z3.tolist(), (z3 - Z2).tolist()
        cal["interbin"] = interbin_check(C3)
        rep, _, _ = nesting(COV13, C3)
        cal["nesting"] = rep
        rho2_eff = {n: rep["bins"][n]["rho2_eff_trace"] for n in BIN_NAMES}
        cal["rho2_eff_trace"] = rho2_eff
        broken = list(rep["broken"])
        if not cal["interbin"]["ok"]:
            broken.append(f"correlación entre bins |r| máx {cal['interbin']['max_abs_interbin_corr']:.2e} ≥ {V3.CORR_TOL}")
        if not (RHO2_RANGE[0] <= rho2_eff["LRG2"] <= RHO2_RANGE[1]):
            broken.append(f"rho2_LRG2,eff = {rho2_eff['LRG2']:.3f} fuera de [{RHO2_RANGE[0]}, {RHO2_RANGE[1]}]")
        cal["broken"] = broken
        cal["status"] = "modelo roto" if broken else "ok"
        if cal["status"] == "ok":
            for ic, combo in enumerate(COMBOS):
                st = V3.setup(combo)
                seeds = {hyp: [SEED_BASE, ic, 500, ih, SEED_HP] for ih, hyp in enumerate(("H0", "H1"))}
                _, case = calibrate(st, C3, z3, seeds)
                cal["combos"][combo] = {"seeds": seeds, "exact_densities_fractional": case["exact_fractional"], "V3k_fractional": case["V3k_fractional"], "m2_bao3_LRG2": case["m2_bao3_LRG2"],
                                        "c0": case["hp"]["c0"], "c1": case["hp"]["c1"], "c1_lt_c0": case["hp"]["c1_lt_c0"], "power": case["hp"]["power"], "P_lt_c1_H0": case["hp"]["P_lt_c1_H0"],
                                        "mc_err": case["hp"]["mc_err"], "p2": st["p2"].tolist(), "r2_absolute_z_DR2": st["r2"].tolist()}
    outp = HERE / f"pronostico_dr3_lrg2_v4_calibracion_{tag}.json"
    assert not outp.exists(), f"ya existe {outp.name}"
    outp.write_text(json.dumps(cal, indent=1, default=float), encoding="utf-8")
    sha = sha256_file(outp)
    (HERE / f"pronostico_dr3_lrg2_v4_calibracion_{tag}.sha256").write_text(f"{sha}  {outp.name}\n", encoding="utf-8")
    print(json.dumps({"status": cal["status"], "reason_or_broken": cal.get("broken", cal.get("reason")), "z_shifts": cal.get("z_shifts_per_row"), "labels": {k: lrep[k] for k in ("unmatched", "z_shifts", "ignored_dr3_rows")},
                      "calibracion": {c: {k: cal["combos"][c][k] for k in ("c0", "c1", "c1_lt_c0", "power")} for c in cal["combos"]}, "sha256": sha}, indent=1))
    return outp


def evaluar(cov_path, val_path, tag):
    calp = HERE / f"pronostico_dr3_lrg2_v4_calibracion_{tag}.json"
    resp = HERE / f"pronostico_dr3_lrg2_v4_resultado_{tag}.json"
    assert calp.exists(), "primero calibrar"
    assert not resp.exists(), f"ya existe {resp.name}: Λ se calcula una sola vez (control de coherencia; el compromiso es el orden de los commits públicos)"
    sha_rec = (HERE / f"pronostico_dr3_lrg2_v4_calibracion_{tag}.sha256").read_text().split()[0]
    assert sha256_file(calp) == sha_rec, "la calibración no coincide con su sha256 registrado"
    cal = json.load(open(calp))
    assert cal["cov_sha256"] == sha256_file(cov_path), "la covarianza no es la usada en la calibración"
    sc = self_check()
    assert sc["match"] is not False, "el script no coincide con el sha256 registrado"
    res = {"mode": "evaluar", "tag": tag, "calibration_sha256": sha_rec, "values_file": str(val_path), "values_sha256": sha256_file(val_path), "script_check": sc, "status": cal["status"], "principal_combo": PRINCIPAL, "combos": {}}
    if cal["status"] != "ok":
        res["result_principal"] = cal["status"]
    else:
        j = json.load(open(cov_path))
        idx, _ = map_labels(list(j["labels"]))
        C3 = np.array(j["cov"], float)[np.ix_(idx, idx)]
        z3 = np.array(cal["z3"])
        v = json.load(open(val_path))
        idxv, _ = map_labels(list(v["labels"]))
        assert idxv is not None and list(idxv) == list(idx), "las etiquetas del fichero de valores no coinciden con las de la covarianza"
        d3 = D.copy()
        d3[:13] = np.array(v["mean"], float)[idxv]
        for combo in COMBOS:
            c = cal["combos"][combo]
            _, p3 = fit_cz(combo, KEEP, d3, C3, z3)                     # full refit without LRG2, DR3 covariance, model at z_DR3
            J3, F3f, _, _ = fisher_cz(combo, KEEP, p3, C3, z3)
            m3 = model_bao(p3, z3)
            r3_abs = d3[LRG2] - m3[LRG2]
            r3 = r3_abs / m3[LRG2]                                       # fractional, denominator m(p3, z_DR3)
            V_abs = C3[np.ix_(LRG2, LRG2)] + J3[LRG2] @ cho_solve(F3f, J3[LRG2].T)
            V = V_abs / np.outer(m3[LRG2], m3[LRG2])
            dens = {h: Gauss(np.array(c["exact_densities_fractional"][h]["a"]), np.array(c["exact_densities_fractional"][h]["Sigma"])) for h in ("H0", "H1")}
            L = float(dens["H1"].logpdf(r3[None, :])[0] - dens["H0"].logpdf(r3[None, :])[0])
            verdict = "favorece H1" if L > c["c0"] else ("favorece H0" if L < c["c1"] else "no concluyente")
            res["combos"][combo] = {"p3": p3.tolist(), "z3_LRG2": z3[LRG2].tolist(), "m3_LRG2": m3[LRG2].tolist(), "r3_absolute": r3_abs.tolist(), "r3_fractional": r3.tolist(),
                                    "d2_DR3": float(r3 @ cho_solve(cho_factor(0.5 * (V + V.T), lower=True), r3)), "Lambda": L, "c0": c["c0"], "c1": c["c1"], "verdict": verdict, "strong_abs_Lambda_gt_3": bool(abs(L) > LAMBDA_STRONG)}
        res["result_principal"] = res["combos"][PRINCIPAL]["verdict"]
        res["result_secondary_Pantheon+"] = res["combos"]["Pantheon+"]["verdict"]
    resp.write_text(json.dumps(res, indent=1, default=float), encoding="utf-8")
    print(json.dumps({k: res[k] for k in res if k in ("status", "result_principal", "result_secondary_Pantheon+")} | {"combos": {c: {k: v[k] for k in ("z3_LRG2", "r3_fractional", "d2_DR3", "Lambda", "verdict")} for c, v in res["combos"].items()}}, indent=1))
    return resp


def falso():
    """Fake DR3 files (DATOS FALSOS): z_eff of LRG2 = 0.712 (evaluable) and 0.720 (no evaluable), C3 = 0.6 C2, values = DR2."""
    for tag, z in (("falso_z0712", 0.712), ("falso_z0720", 0.720)):
        labels = [l.replace("_z0.706", f"_z{z:.3f}") for l in LABELS]
        p = HERE / f"pronostico_dr3_lrg2_v4_{tag}_dr3_cov.json"
        assert not p.exists()
        p.write_text(json.dumps({"note": f"DATOS FALSOS: covarianza 0.6 × C_DR2 con z_eff de LRG2 = {z}; solo para probar el modo evaluación", "labels": labels, "cov": (0.6 * COV13).tolist()}, indent=1), encoding="utf-8")
        p = HERE / f"pronostico_dr3_lrg2_v4_{tag}_dr3_valores.json"
        assert not p.exists()
        p.write_text(json.dumps({"note": f"DATOS FALSOS: valores de DR2 sin cambio con z_eff de LRG2 = {z}; solo para probar el modo evaluación", "labels": labels, "mean": D[:13].tolist()}, indent=1), encoding="utf-8")
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
