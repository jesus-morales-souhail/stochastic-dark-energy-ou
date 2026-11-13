#!/usr/bin/env python3
"""
REVISION 3: intentos de romper el modo evaluación de pronostico_dr3_lrg2_v4.py (y una prueba de v3). Módulos cargados desde su ruta original con
HERE redirigido a un directorio de trabajo fuera del repositorio (argumento 1). Nada se escribe ni se borra en auditoria_v23 salvo este json.
    ../venv/bin/python -I revision3_pronostico_romper.py /ruta/de/trabajo  -> revision3_pronostico_romper.json
"""
from __future__ import annotations
import hashlib, importlib.util, json, sys, shutil
from pathlib import Path
import numpy as np
AUD = Path(__file__).resolve().parent
WORK = Path(sys.argv[1]).resolve(); WORK.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(AUD))
import pronostico_dr3_lrg2_v3 as V3   # noqa: E402  (HERE of v3 stays in auditoria; only used for its functions)
spec = importlib.util.spec_from_file_location("v4mod", AUD / "pronostico_dr3_lrg2_v4.py")
V4 = importlib.util.module_from_spec(spec); spec.loader.exec_module(V4)
V4.HERE = WORK
LABELS, COV13, D = V4.LABELS, V4.COV13, V4.D
res = {"work_dir": str(WORK), "attempts": {}}


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def attempt(name, fn, expect):
    try:
        r = fn()
        res["attempts"][name] = {"expected": expect, "outcome": "EJECUTÓ sin negarse", "detail": r}
    except AssertionError as e:
        res["attempts"][name] = {"expected": expect, "outcome": "se negó (AssertionError)", "detail": str(e)}
    except Exception as e:
        res["attempts"][name] = {"expected": expect, "outcome": "error distinto (no es una salida declarada)", "detail": repr(e)[:300]}
    print(f"{name}: {res['attempts'][name]['outcome']} — {str(res['attempts'][name]['detail'])[:200]}", flush=True)


def wjson(name, obj):
    p = WORK / name; p.write_text(json.dumps(obj, indent=1), encoding="utf-8"); return p


def lab(z_lrg2=0.706, extra=None, rename=None):
    L = [l.replace("_z0.706", f"_z{z_lrg2:.3f}") for l in LABELS]
    if rename:
        L = [rename.get(l, l) for l in L]
    return L + (extra or [])


def calsum(tag):
    c = json.loads((WORK / f"pronostico_dr3_lrg2_v4_calibracion_{tag}.json").read_text())
    return {"status": c["status"], "z_shifts": c.get("z_shifts_per_row"), "reason": c.get("reason"), "broken": c.get("broken"), "script_check": c["script_check"], "c0": {k: v["c0"] for k, v in c["combos"].items()}}


def ressum(p):
    r = json.load(open(p)); return {k: r[k] for k in r if k in ("status", "result_principal", "result_secondary_Pantheon+", "script_check")} | {"L": {c: v["Lambda"] for c, v in r["combos"].items()}}


cov06 = wjson("cov_z0706.json", {"labels": lab(), "cov": (0.6 * COV13).tolist()})
cov12 = wjson("cov_z0712.json", {"labels": lab(0.712), "cov": (0.6 * COV13).tolist()})
val06 = wjson("val_z0706.json", {"labels": lab(), "mean": D[:13].tolist()})
val12 = wjson("val_z0712.json", {"labels": lab(0.712), "mean": D[:13].tolist()})
# 0. the script's own sha256 file is NOT in WORK: self_check -> registered None
res["attempts"]["0_self_check_sin_fichero_sha256"] = {"outcome": "informativo", "detail": V4.self_check()}
print("0:", V4.self_check())
# 1. evaluar sin calibrar
attempt("1_evaluar_sin_calibrar", lambda: V4.evaluar(cov12, val12, "a1"), "negarse")
# 2. calibrar con z 0.712 y evaluar con valores etiquetados z 0.706 (z distinto entre ficheros) — sin fichero .sha256 del script en WORK
attempt("2a_calibrar_z0712", lambda: calsum("a2") if V4.calibrar(cov12, "a2") else None, "ok, shift +0.006")
attempt("2b_evaluar_valores_z0706_con_cov_z0712_y_sin_sha256_del_script", lambda: ressum(V4.evaluar(cov12, val06, "a2")), "¿se niega? (etiquetas con z distinto; self_check None)")
# 3. now with the registered .sha256 copied into WORK: self_check True; then a wrong one
shutil.copy(AUD / "pronostico_dr3_lrg2_v4.sha256", WORK / "pronostico_dr3_lrg2_v4.sha256")
res["attempts"]["3_self_check_con_sha256_registrado"] = {"outcome": "informativo", "detail": V4.self_check()}
attempt("3a_calibrar_z0712_bis", lambda: calsum("a3") if V4.calibrar(cov12, "a3") else None, "ok")
attempt("3b_evaluar_normal", lambda: ressum(V4.evaluar(cov12, val12, "a3")), "ok una vez")
attempt("3c_evaluar_otra_vez", lambda: V4.evaluar(cov12, val12, "a3"), "negarse")
(WORK / "pronostico_dr3_lrg2_v4.sha256").write_text("0" * 64 + "  pronostico_dr3_lrg2_v4.py\n")
attempt("3d_calibrar_con_sha256_registrado_incorrecto", lambda: calsum("a3d") if V4.calibrar(cov12, "a3d") else None, "calibrar no comprueba (solo registra match False)")
attempt("3e_evaluar_con_sha256_registrado_incorrecto", lambda: V4.evaluar(cov12, val12, "a3d"), "negarse")
shutil.copy(AUD / "pronostico_dr3_lrg2_v4.sha256", WORK / "pronostico_dr3_lrg2_v4.sha256")
# 4. values file with rows reordered (same labels) -> declared outcome?
perm = list(range(13))[::-1]
valperm = wjson("val_perm.json", {"labels": [lab(0.712)[i] for i in perm], "mean": [D[:13][i] for i in perm]})
attempt("4a_calibrar", lambda: calsum("a4") if V4.calibrar(cov12, "a4") else None, "ok")
attempt("4b_evaluar_valores_reordenados", lambda: ressum(V4.evaluar(cov12, valperm, "a4")), "¿se niega con AssertionError (no es 'no evaluable')?")
# 5. duplicate LRG rows within tolerance (0.706 and 0.712) -> ambiguous -> no evaluable
Cd = np.zeros((15, 15)); Cd[:13, :13] = 0.6 * COV13; Cd[13:, 13:] = 0.6 * COV13[np.ix_([3, 4], [3, 4])]
covdup = wjson("cov_dup.json", {"labels": lab(0.706, extra=["DM_LRG_z0.712", "DH_LRG_z0.712"]), "cov": Cd.tolist()})
attempt("5_calibrar_LRG_duplicado_0.706_y_0.712", lambda: calsum("a5") if V4.calibrar(covdup, "a5") else None, "no evaluable (ambiguo)")
# 6. boundary of the tolerance: LRG2 at 0.716 (|dz| = 0.01, inclusive) and 0.7161
for z, tag in ((0.716, "a6a"), (0.7161, "a6b")):
    c = wjson(f"cov_z{z}.json", {"labels": lab(z), "cov": (0.6 * COV13).tolist()})
    attempt(f"6_calibrar_LRG2_z{z}", lambda c=c, tag=tag: calsum(tag) if V4.calibrar(c, tag) else None, "0.716 ok (inclusive); 0.7161 no evaluable")
# 7. shift in another bin only (LRG1 0.510 -> 0.520) -> ok with the shift recorded
covl1 = wjson("cov_lrg1.json", {"labels": [l.replace("_z0.510", "_z0.520") for l in lab()], "cov": (0.6 * COV13).tolist()})
attempt("7_calibrar_LRG1_desplazado_+0.01", lambda: calsum("a7") if V4.calibrar(covl1, "a7") else None, "ok, shift registrado en LRG1")
# 8. tracer renamed (LRG -> LRG2) -> unmatched -> no evaluable
covren = wjson("cov_ren.json", {"labels": lab(rename={"DM_LRG_z0.706": "DM_LRG2_z0.706", "DH_LRG_z0.706": "DH_LRG2_z0.706"}), "cov": (0.6 * COV13).tolist()})
attempt("8_calibrar_trazador_renombrado", lambda: calsum("a8") if V4.calibrar(covren, "a8") else None, "no evaluable")
# 9. covariance file that also contains 'mean': calibrar must not read it (keys_read) and the sha covers the whole file
covmean = wjson("cov_with_mean.json", {"labels": lab(0.712), "cov": (0.6 * COV13).tolist(), "mean": (D[:13] * 1.01).tolist()})
attempt("9_calibrar_fichero_con_mean", lambda: json.loads((WORK / "pronostico_dr3_lrg2_v4_calibracion_a9.json").read_text())["keys_read"] if V4.calibrar(covmean, "a9") else None, "ok; keys_read = labels, cov")
# 10. boundary: LRG2 at 0.715 (ok), 0.716 (|dz| = 0.01: text says inclusive), 0.696, and LRG1 shifted +0.009 (shift recorded)
for z, tag in ((0.715, "b10a"), (0.696, "b10b")):
    c = wjson(f"cov_z{z}.json", {"labels": lab(z), "cov": (0.6 * COV13).tolist()})
    attempt(f"10_calibrar_LRG2_z{z}", lambda c=c, tag=tag: calsum(tag) if V4.calibrar(c, tag) else None, "0.715 ok; 0.696 debería ser ok (±0.01 inclusive)")
covl1b = wjson("cov_lrg1b.json", {"labels": [l.replace("_z0.510", "_z0.519") for l in lab()], "cov": (0.6 * COV13).tolist()})
attempt("10c_calibrar_LRG1_desplazado_+0.009", lambda: calsum("b10c") if V4.calibrar(covl1b, "b10c") else None, "ok, shift +0.009 registrado en LRG1")
res["attempts"]["10d_aritmetica"] = {"outcome": "informativo", "detail": {"0.716-0.706": 0.716 - 0.706, "0.706-0.696": 0.706 - 0.696, "0.52-0.51": 0.52 - 0.51, "<=0.01": [0.716 - 0.706 <= 0.01, 0.706 - 0.696 <= 0.01]}}
# 11. rewrite a valid calibration (a4) + its .sha256 (known hole) -> runs with the altered c0
def tamper():
    p = WORK / "pronostico_dr3_lrg2_v4_calibracion_a4.json"; c = json.loads(p.read_text()); c["combos"]["sinSN"]["c0"] = 9.0; p.write_text(json.dumps(c, indent=1))
    (WORK / "pronostico_dr3_lrg2_v4_calibracion_a4.sha256").write_text(f"{sha(p)}  x\n"); return "c0 = 9 y sha reescrito"
attempt("11a_manipular_calibracion_y_sha", tamper, "—")
attempt("11b_evaluar_tras_manipular", lambda: ressum(V4.evaluar(cov12, val12, "a4")), "se ejecuta con c0 = 9 (agujero conocido; compromiso = commits públicos)")
# 12. values with extra new bins -> positions differ -> ?
valx = wjson("val_extra.json", {"labels": ["DM_LRG_z0.300"] + lab(0.712), "mean": [10.0] + D[:13].tolist()})
attempt("12a_calibrar", lambda: calsum("b12") if V4.calibrar(cov12, "b12") else None, "ok")
attempt("12_evaluar_valores_con_bin_nuevo_delante", lambda: ressum(V4.evaluar(cov12, valx, "b12")), "¿AssertionError por índices distintos aunque las etiquetas homólogas existan?")
res["work_files"] = sorted(p.name for p in WORK.iterdir())
(AUD / "revision3_pronostico_romper.json").write_text(json.dumps(res, indent=1, default=str), encoding="utf-8")
print("escrito revision3_pronostico_romper.json")
