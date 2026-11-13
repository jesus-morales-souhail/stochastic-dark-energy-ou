#!/usr/bin/env python3
"""
PRONOSTICO_DR3_REGISTRO_v5.md: version 5 of the DR3 forecast of the LRG2 residual. Formally replaces v3 and v4 for the registration (they stay as
history). Imports pronostico_dr3_lrg2_v4.py (not modified) and applies the corrections of the third review (REVISION3_PRONOSTICO_DR3.md, commit
53bcb70). Unchanged: P18, P19, P23, P24, H0, H1, the matrix nesting, the thresholds, the decision rule, the simulations (none is redone).
  a. Inclusive boundary of rule f with an explicit tolerance: a DR3 bin is the homologue of a DR2 bin if it has the same kind and tracer and
     |z_DR3 - z_DR2| <= 0.01 + 1e-9 (floating point: 0.716 - 0.706 = 0.010000000000000009).
  b. Self-checks that fail closed: if the .sha256 file of the script or of the calibration is missing, or the hash does not match, `calibrar`
     and `evaluar` refuse with a declared output ("control de coherencia fallido"), never with an exception. They remain consistency controls,
     not the commitment (the commitment is the public commit order).
  c. Pairing by identity in `evaluar`: each row of the values file is paired with the calibration by tracer, quantity (D_M/r_d or D_H/r_d) and
     z_eff (|dz| <= 1e-9: the values must be at the z_eff that was calibrated, from the same release), not by index. Values at another z,
     missing rows or ambiguous rows -> "no evaluable"; reordered rows -> evaluable; extra rows (new bins) -> ignored (rule e).
  d. Frozen denominator: the fractional residual and its covariance in `evaluar` use m(p2, z_DR3) stored in the calibration (the same as in
     the densities), not m(p3, z_DR3); the Lambda with the v4 denominator is also reported for comparison.
  e. With the DR2 z_eff, the v5 calibration reproduces the v4 one to rounding (checked in mode registro).
Modes:
    ../venv/bin/python -I pronostico_dr3_lrg2_v5.py registro                  # boundary tests, v5 = v4 check -> pronostico_dr3_lrg2_v5.json (+ .sha256 of the script)
    ../venv/bin/python -I pronostico_dr3_lrg2_v5.py calibrar COV.json TAG     # reads only labels (z_eff) and cov
    ../venv/bin/python -I pronostico_dr3_lrg2_v5.py evaluar COV.json VALORES.json TAG
    ../venv/bin/python -I pronostico_dr3_lrg2_v5.py falso                     # extra fake files (DATOS FALSOS) for the identity-pairing tests
Environment variable PRONOSTICO_V5_SHA256_FILE (tests only) points the script self-check to another registry file.
Nothing existing is modified.
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

import numpy as np
from scipy.linalg import cho_factor, cho_solve

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import pronostico_dr3_lrg2_v4 as V4  # noqa: E402

V3, P, JB = V4.V3, V4.P, V4.JB
D, COV13, LRG2, KEEP, Z2 = V4.D, V4.COV13, V4.LRG2, V4.KEEP, V4.Z2
BINS, BIN_NAMES, LABELS = V4.BINS, V4.BIN_NAMES, V4.LABELS
SEED_BASE, SEED_HP, COMBOS, PRINCIPAL = V4.SEED_BASE, V4.SEED_HP, V4.COMBOS, V4.PRINCIPAL
RHO2_RANGE, LAMBDA_STRONG = V4.RHO2_RANGE, V4.LAMBDA_STRONG
Z_TOL, Z_EPS, Z_PAIR_EPS = 0.01, 1e-9, 1e-9
sha256_file, git_head, parse_label, nesting, interbin_check = V4.sha256_file, V4.git_head, V4.parse_label, V4.nesting, V4.interbin_check
FAIL = "control de coherencia fallido"


# ----------------------------------------------------------------------------------------------------------------------
# rule f with inclusive boundary; identity pairing
# ----------------------------------------------------------------------------------------------------------------------
def map_labels(labels3):
    """DR2 row -> DR3 row: same kind and tracer, |z3 - z2| <= Z_TOL + Z_EPS, unique."""
    parsed3 = [parse_label(l) for l in labels3]
    idx, rep = [], {"tolerance": Z_TOL + Z_EPS, "matches": {}, "unmatched": [], "ambiguous": [], "z_shifts": {}, "ignored_dr3_rows": []}
    for lab in LABELS:
        k2, t2, z2 = parse_label(lab)
        cand = [i for i, (k3, t3, z3) in enumerate(parsed3) if k3 == k2 and t3 == t2 and abs(z3 - z2) <= Z_TOL + Z_EPS]
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


def pair_values(labels_cal, labels_val):
    """Calibrated row (kind, tracer, z_cal) -> values row with the same kind and tracer and |z - z_cal| <= Z_PAIR_EPS."""
    pv = [parse_label(l) for l in labels_val]
    idx, rep = [], {"pair_tolerance_z": Z_PAIR_EPS, "pairs": {}, "missing_or_other_z": [], "ambiguous": [], "ignored_value_rows": []}
    for lab in labels_cal:
        k, t, z = parse_label(lab)
        cand = [i for i, (kv, tv, zv) in enumerate(pv) if kv == k and tv == t and abs(zv - z) <= Z_PAIR_EPS]
        if len(cand) == 1:
            idx.append(cand[0])
            rep["pairs"][lab] = labels_val[cand[0]]
        elif not cand:
            rep["missing_or_other_z"].append(lab)
        else:
            rep["ambiguous"].append(lab)
    rep["ignored_value_rows"] = [labels_val[i] for i in range(len(labels_val)) if i not in idx]
    rep["ok"] = not rep["missing_or_other_z"] and not rep["ambiguous"]
    return (np.array(idx) if rep["ok"] else None), rep


def self_check():
    reg = Path(os.environ.get("PRONOSTICO_V5_SHA256_FILE", HERE / "pronostico_dr3_lrg2_v5.sha256"))
    cur = sha256_file(__file__)
    if not reg.exists():
        return {"registered": None, "current": cur, "match": None, "status": FAIL, "reason": f"falta el fichero {reg.name}"}
    registered = reg.read_text().split()[0]
    ok = registered == cur
    return {"registered": registered, "current": cur, "match": ok, "status": "ok" if ok else FAIL, "reason": None if ok else "el sha256 del script no coincide con el registrado"}


def write_json(path, obj):
    path.write_text(json.dumps(obj, indent=1, default=float), encoding="utf-8")


# ----------------------------------------------------------------------------------------------------------------------
# DR3 evaluation mode
# ----------------------------------------------------------------------------------------------------------------------
def calibrar(cov_path, tag):
    outp = HERE / f"pronostico_dr3_lrg2_v5_calibracion_{tag}.json"
    assert not outp.exists(), f"ya existe {outp.name}"
    cal = {"mode": "calibrar", "version": "v5", "tag": tag, "cov_file": str(cov_path), "cov_sha256": sha256_file(cov_path), "script_check": self_check(), "git_head_at_run": git_head(),
           "keys_read": ["labels", "cov"], "principal_combo": PRINCIPAL, "combos": {}}
    if cal["script_check"]["status"] != "ok":
        cal["status"] = FAIL
        cal["reason"] = cal["script_check"]["reason"]
    else:
        j = json.load(open(cov_path))
        labels3 = list(j["labels"])
        idx, lrep = map_labels(labels3)
        cal["labels"] = lrep
        if idx is None:
            cal["status"] = "no evaluable"
            cal["reason"] = "bins de DR2 sin homólogo único en DR3 (mismo tipo y trazador, |Δz| ≤ 0.01 + 1e-9)"
        else:
            z3 = np.array([parse_label(labels3[i])[2] for i in idx])
            C3 = np.array(j["cov"], float)[np.ix_(idx, idx)]
            cal["labels_calibrated"] = [labels3[i] for i in idx]
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
                    _, case = V4.calibrate(st, C3, z3, seeds)
                    cal["combos"][combo] = {"seeds": seeds, "exact_densities_fractional": case["exact_fractional"], "V3k_fractional": case["V3k_fractional"], "m2_bao3_LRG2": case["m2_bao3_LRG2"],
                                            "c0": case["hp"]["c0"], "c1": case["hp"]["c1"], "c1_lt_c0": case["hp"]["c1_lt_c0"], "power": case["hp"]["power"], "P_lt_c1_H0": case["hp"]["P_lt_c1_H0"],
                                            "mc_err": case["hp"]["mc_err"], "p2": st["p2"].tolist(), "r2_absolute_z_DR2": st["r2"].tolist()}
    write_json(outp, cal)
    sha = sha256_file(outp)
    (HERE / f"pronostico_dr3_lrg2_v5_calibracion_{tag}.sha256").write_text(f"{sha}  {outp.name}\n", encoding="utf-8")
    print(json.dumps({"status": cal["status"], "reason_or_broken": cal.get("broken", cal.get("reason")), "script_check": cal["script_check"]["status"], "z_shifts": cal.get("z_shifts_per_row"),
                      "labels": ({k: cal["labels"][k] for k in ("unmatched", "ambiguous", "z_shifts", "ignored_dr3_rows")} if "labels" in cal else None),
                      "calibracion": {c: {k: cal["combos"][c][k] for k in ("c0", "c1", "c1_lt_c0", "power")} for c in cal["combos"]}, "sha256": sha}, indent=1))
    return outp


def evaluar(cov_path, val_path, tag):
    resp = HERE / f"pronostico_dr3_lrg2_v5_resultado_{tag}.json"
    assert not resp.exists(), f"ya existe {resp.name}: Λ se calcula una sola vez (control de coherencia; el compromiso es el orden de los commits públicos)"
    calp = HERE / f"pronostico_dr3_lrg2_v5_calibracion_{tag}.json"
    shap = HERE / f"pronostico_dr3_lrg2_v5_calibracion_{tag}.sha256"
    res = {"mode": "evaluar", "version": "v5", "tag": tag, "values_file": str(val_path), "values_sha256": sha256_file(val_path), "script_check": self_check(), "principal_combo": PRINCIPAL, "combos": {}}

    def declare(status, reason):
        res["status"], res["result_principal"], res["reason"] = status, status, reason
        write_json(resp, res)
        print(json.dumps({"status": status, "result_principal": status, "reason": reason}, indent=1))
        return resp

    if res["script_check"]["status"] != "ok":
        return declare(FAIL, res["script_check"]["reason"])
    if not calp.exists():
        return declare(FAIL, "falta la calibración")
    if not shap.exists():
        return declare(FAIL, "falta el .sha256 de la calibración")
    sha_rec = shap.read_text().split()[0]
    res["calibration_sha256"] = sha_rec
    if sha256_file(calp) != sha_rec:
        return declare(FAIL, "la calibración no coincide con su sha256 registrado")
    cal = json.load(open(calp))
    if cal.get("script_check", {}).get("status") != "ok" or cal["script_check"]["current"] != res["script_check"]["current"]:
        return declare(FAIL, "la calibración se hizo con un script distinto o con el control fallido")
    if cal["cov_sha256"] != sha256_file(cov_path):
        return declare(FAIL, "la covarianza no es la usada en la calibración")
    if cal["status"] != "ok":
        return declare(cal["status"], cal.get("reason", cal.get("broken")))
    # identity pairing of the values rows with the calibrated rows
    v = json.load(open(val_path))
    idxv, prep = pair_values(cal["labels_calibrated"], list(v["labels"]))
    res["pairing"] = prep
    if idxv is None:
        return declare("no evaluable", "filas de valores sin pareja en la calibración (trazador, cantidad y z_eff idénticos): " + ", ".join(prep["missing_or_other_z"] + prep["ambiguous"]))
    j = json.load(open(cov_path))
    idxc, _ = map_labels(list(j["labels"]))
    C3 = np.array(j["cov"], float)[np.ix_(idxc, idxc)]
    z3 = np.array(cal["z3"])
    d3 = D.copy()
    d3[:13] = np.array(v["mean"], float)[idxv]
    for combo in COMBOS:
        c = cal["combos"][combo]
        _, p3 = V4.fit_cz(combo, KEEP, d3, C3, z3)
        J3, F3f, _, _ = V4.fisher_cz(combo, KEEP, p3, C3, z3)
        m3 = V4.model_bao(p3, z3)
        m2 = np.array(c["m2_bao3_LRG2"])                                 # frozen denominator (same as in the densities)
        r3_abs = d3[LRG2] - m3[LRG2]
        r3 = r3_abs / m2
        V_abs = C3[np.ix_(LRG2, LRG2)] + J3[LRG2] @ cho_solve(F3f, J3[LRG2].T)
        V = V_abs / np.outer(m2, m2)
        dens = {h: V3.Gauss(np.array(c["exact_densities_fractional"][h]["a"]), np.array(c["exact_densities_fractional"][h]["Sigma"])) for h in ("H0", "H1")}
        lam = float(dens["H1"].logpdf(r3[None, :])[0] - dens["H0"].logpdf(r3[None, :])[0])
        r3_v4 = r3_abs / m3[LRG2]
        lam_v4 = float(dens["H1"].logpdf(r3_v4[None, :])[0] - dens["H0"].logpdf(r3_v4[None, :])[0])
        verdict = "favorece H1" if lam > c["c0"] else ("favorece H0" if lam < c["c1"] else "no concluyente")
        res["combos"][combo] = {"p3": p3.tolist(), "z3_LRG2": z3[LRG2].tolist(), "m3_LRG2_p3": m3[LRG2].tolist(), "m2_LRG2_frozen": m2.tolist(), "r3_absolute": r3_abs.tolist(), "r3_fractional": r3.tolist(),
                                "d2_DR3": float(r3 @ cho_solve(cho_factor(0.5 * (V + V.T), lower=True), r3)), "Lambda": lam, "Lambda_v4_denominator_m_p3": lam_v4, "dLambda_vs_v4": lam - lam_v4,
                                "c0": c["c0"], "c1": c["c1"], "verdict": verdict, "strong_abs_Lambda_gt_3": bool(abs(lam) > LAMBDA_STRONG)}
    res["status"] = "ok"
    res["result_principal"] = res["combos"][PRINCIPAL]["verdict"]
    res["result_secondary_Pantheon+"] = res["combos"]["Pantheon+"]["verdict"]
    write_json(resp, res)
    print(json.dumps({"status": "ok", "result_principal": res["result_principal"], "result_secondary_Pantheon+": res["result_secondary_Pantheon+"], "ignored_value_rows": prep["ignored_value_rows"],
                      "combos": {c: {k: vv[k] for k in ("z3_LRG2", "r3_fractional", "d2_DR3", "Lambda", "Lambda_v4_denominator_m_p3", "verdict")} for c, vv in res["combos"].items()}}, indent=1))
    return resp


# ----------------------------------------------------------------------------------------------------------------------
def registration():
    t0 = time.time()
    v4 = json.load(open(HERE / "pronostico_dr3_lrg2_v4.json"))
    out = {"version": "v5 (correcciones de REVISION3_PRONOSTICO_DR3.md; sustituye a v3 y v4 para el registro)", "script_sha256": sha256_file(__file__), "git_head_at_run": git_head(),
           "z_tolerance_rule_f": Z_TOL + Z_EPS, "z_pair_tolerance": Z_PAIR_EPS, "boundary_tests": [], "v5_vs_v4": {}}
    # a. boundary tests of rule f (labels of DR2 with one z replaced)
    for name, base_lab, zs in (("LRG2", "DM_LRG_z0.706", (0.696, 0.6965, 0.706, 0.7155, 0.716, 0.6959, 0.7161, 0.720)), ("LRG1", "DM_LRG_z0.510", (0.500, 0.520, 0.4999, 0.5201))):
        z0 = parse_label(base_lab)[2]
        for z in zs:
            labels = [l.replace(f"_z{z0:.3f}", f"_z{z:g}") if base_lab.split("_z")[0].split("_", 1)[1] in l and f"_z{z0:.3f}" in l else l for l in LABELS]
            idx, rep = map_labels(labels)
            out["boundary_tests"].append({"bin": name, "z": z, "dz": z - z0, "evaluable": bool(rep["ok"]), "unmatched": rep["unmatched"], "ambiguous": rep["ambiguous"]})
            print(f"frontera {name} z={z:g} dz={z - z0:+.4f}: {'evaluable' if rep['ok'] else 'no evaluable'}", flush=True)
    # e. v5 = v4 (rho2 = 0.6, DR2 redshifts, same seeds as v4 registration)
    for ic, combo in enumerate(COMBOS):
        st = V3.setup(combo)
        seeds = {hyp: [SEED_BASE, ic, 0, ih, SEED_HP] for ih, hyp in enumerate(("H0", "H1"))}
        _, case = V4.calibrate(st, 0.6 * COV13, Z2, seeds)
        vc = v4["combos"][combo]["cases"]["rho2_0.6"]
        out["v5_vs_v4"][combo] = {"d_c0": case["hp"]["c0"] - vc["hp"]["c0"], "d_c1": case["hp"]["c1"] - vc["hp"]["c1"], "d_power": case["hp"]["power"] - vc["hp"]["power"],
                                  "max_abs_d_a": max(float(np.max(np.abs(np.array(case["exact_fractional"][h]["a"]) - np.array(vc["exact_fractional"][h]["a"])))) for h in ("H0", "H1")),
                                  "max_abs_d_Sigma": max(float(np.max(np.abs(np.array(case["exact_fractional"][h]["Sigma"]) - np.array(vc["exact_fractional"][h]["Sigma"])))) for h in ("H0", "H1"))}
        print(f"v5 vs v4 {combo}: {out['v5_vs_v4'][combo]}", flush=True)
    out["runtime_s"] = time.time() - t0
    write_json(HERE / "pronostico_dr3_lrg2_v5.json", out)
    (HERE / "pronostico_dr3_lrg2_v5.sha256").write_text(f"{out['script_sha256']}  pronostico_dr3_lrg2_v5.py\n", encoding="utf-8")
    print("escrito pronostico_dr3_lrg2_v5.json", f"{out['runtime_s']:.0f} s")


def falso():
    """Extra fake files (DATOS FALSOS) for the identity-pairing tests: reordered rows; a new bin in front; values labelled at z0.706 exist already
    (pronostico_dr3_lrg2_v3_falso_dr3_valores.json) for the 'other z' test; a wrong .sha256 registry for the altered-hash test."""
    lab712 = [l.replace("_z0.706", "_z0.712") for l in LABELS]
    vals = D[:13].tolist()
    perm = list(range(13))[::-1]
    p = HERE / "pronostico_dr3_lrg2_v5_falso_reordenado_dr3_valores.json"
    assert not p.exists()
    p.write_text(json.dumps({"note": "DATOS FALSOS: valores de DR2 (z_eff de LRG2 = 0.712) con las filas en orden inverso; solo para probar el emparejamiento", "labels": [lab712[i] for i in perm], "mean": [vals[i] for i in perm]}, indent=1), encoding="utf-8")
    p = HERE / "pronostico_dr3_lrg2_v5_falso_binnuevo_dr3_valores.json"
    assert not p.exists()
    p.write_text(json.dumps({"note": "DATOS FALSOS: valores de DR2 (z_eff de LRG2 = 0.712) con un bin nuevo delante (DV_ELG2b_z1.100); solo para probar el emparejamiento", "labels": ["DV_ELG2b_z1.100"] + lab712, "mean": [30.0] + vals}, indent=1), encoding="utf-8")
    p = HERE / "pronostico_dr3_lrg2_v5_falso_sha256_alterado.sha256"
    assert not p.exists()
    p.write_text("0" * 64 + "  pronostico_dr3_lrg2_v5.py\n", encoding="utf-8")
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
