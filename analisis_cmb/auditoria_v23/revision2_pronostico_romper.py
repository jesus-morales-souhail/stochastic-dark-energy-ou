#!/usr/bin/env python3
"""
REVISION 2: intentos de romper el modo evaluación DR3 de pronostico_dr3_lrg2_v2.py. El módulo se carga desde su ruta original (así
script_sha256 es el registrado) pero con HERE redirigido a un directorio de trabajo aparte (argumento 1), de modo que calibrar/evaluar
escriben ALLÍ y no en auditoria_v23. Nada existente se modifica; no se borra nada (los intentos que exigirían borrar se describen, no se hacen).
    ../venv/bin/python -I revision2_pronostico_romper.py /ruta/de/trabajo   -> revision2_pronostico_romper.json (en auditoria_v23)
"""
from __future__ import annotations
import hashlib, importlib.util, json, sys, traceback
from pathlib import Path
import numpy as np
AUD = Path(__file__).resolve().parent
WORK = Path(sys.argv[1]).resolve()
WORK.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(AUD))
spec = importlib.util.spec_from_file_location("v2mod", AUD / "pronostico_dr3_lrg2_v2.py")
V2 = importlib.util.module_from_spec(spec); spec.loader.exec_module(V2)
V2.HERE = WORK                      # all outputs of calibrar/evaluar/falso go to WORK
LABELS, COV13, D = V2.LABELS, V2.COV13, V2.D
res = {"work_dir": str(WORK), "labels_DR2": LABELS, "attempts": {}}


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def attempt(name, fn, expect):
    try:
        r = fn()
        res["attempts"][name] = {"expected": expect, "outcome": "EJECUTÓ sin negarse", "detail": r}
    except AssertionError as e:
        res["attempts"][name] = {"expected": expect, "outcome": "se negó", "detail": str(e)}
    except Exception as e:
        res["attempts"][name] = {"expected": expect, "outcome": "error distinto", "detail": repr(e)}
    print(f"{name}: {res['attempts'][name]['outcome']} — {str(res['attempts'][name]['detail'])[:160]}")


def wjson(name, obj):
    p = WORK / name; p.write_text(json.dumps(obj, indent=1), encoding="utf-8"); return p


cov_ok = wjson("cov_ok.json", {"labels": LABELS, "cov": (0.6 * COV13).tolist()})
cov_other = wjson("cov_other.json", {"labels": LABELS, "cov": (0.62 * COV13).tolist()})
vals = wjson("vals.json", {"labels": LABELS, "mean": D[:13].tolist()})
vals2 = wjson("vals2.json", {"labels": LABELS, "mean": (D[:13] * 1.001).tolist()})

# 1. evaluate without calibrating
attempt("1_evaluar_sin_calibrar", lambda: V2.evaluar(cov_ok, vals, "t1"), "negarse: primero calibrar")
# 2. calibrate, then evaluate normally, then evaluate again (same tag)
attempt("2a_calibrar", lambda: str(V2.calibrar(cov_ok, "t2")), "ok")
attempt("2b_evaluar", lambda: json.load(open(V2.evaluar(cov_ok, vals, "t2")))["result"], "ok, una vez")
attempt("2c_evaluar_otra_vez_mismo_tag", lambda: V2.evaluar(cov_ok, vals, "t2"), "negarse")
attempt("2d_evaluar_otra_vez_otros_valores_mismo_tag", lambda: V2.evaluar(cov_ok, vals2, "t2"), "negarse")
# 3. evaluate with a covariance different from the calibrated one
attempt("3a_calibrar", lambda: str(V2.calibrar(cov_ok, "t3")), "ok")
attempt("3b_evaluar_con_otra_cov", lambda: V2.evaluar(cov_other, vals, "t3"), "negarse: covarianza distinta")
# 4. tamper with the calibration after calibrating (change c0), without touching the .sha256
attempt("4a_calibrar", lambda: str(V2.calibrar(cov_ok, "t4")), "ok")
def tamper4():
    p = WORK / "pronostico_dr3_lrg2_v2_calibracion_t4.json"; c = json.loads(p.read_text()); c["combos"]["sinSN"]["c0"] = 5.0; p.write_text(json.dumps(c, indent=1)); return "c0 cambiado a 5.0"
attempt("4b_manipular_calibracion", tamper4, "—")
attempt("4c_evaluar_tras_manipular", lambda: V2.evaluar(cov_ok, vals, "t4"), "negarse: sha256 de la calibración")
# 5. same, but also rewrite the .sha256 file (the attacker controls the whole directory)
def tamper5():
    p = WORK / "pronostico_dr3_lrg2_v2_calibracion_t4.json"; (WORK / "pronostico_dr3_lrg2_v2_calibracion_t4.sha256").write_text(f"{sha(p)}  {p.name}\n"); return "sha256 reescrito"
attempt("5a_reescribir_sha256", tamper5, "—")
attempt("5b_evaluar_con_sha_reescrito", lambda: json.load(open(V2.evaluar(cov_ok, vals, "t4")))["combos"]["sinSN"]["verdict"], "se ejecuta (el .sha256 no es un compromiso externo)")
# 6. a second evaluation of the same data under another tag (the 'once' is per tag)
attempt("6a_calibrar_tag_nuevo", lambda: str(V2.calibrar(cov_ok, "t6")), "ok")
attempt("6b_evaluar_mismos_datos_otro_tag", lambda: json.load(open(V2.evaluar(cov_ok, vals2, "t6")))["result"], "se ejecuta (otra etiqueta)")
# 7. labels permuted / off-block tiny value / rho2 out of range / corr change / component spread
perm = LABELS[1:] + LABELS[:1]
cov_perm = wjson("cov_perm.json", {"labels": perm, "cov": (0.6 * COV13).tolist()})
attempt("7a_calibrar_etiquetas_permutadas", lambda: str(V2.calibrar(cov_perm, "t7a")), "negarse: etiquetas")
Ct = 0.6 * COV13; Ct[0, 1] = Ct[1, 0] = 1e-15
attempt("7b_calibrar_fuera_de_bloque_1e-15", lambda: json.loads((WORK / "pronostico_dr3_lrg2_v2_calibracion_t7b.json").read_text())["rupture_checks"]["status"] if V2.calibrar(wjson("cov_off.json", {"labels": LABELS, "cov": Ct.tolist()}), "t7b") else None, "modelo roto por un 1e-15 (umbral != 0)")
Cs = 0.6 * COV13; Cs[3, 3] *= 0.9; Cs[4, 4] *= 1.1
attempt("7c_spread_LRG2_20pct_exacto", lambda: json.loads((WORK / "pronostico_dr3_lrg2_v2_calibracion_t7c.json").read_text())["rupture_checks"] if V2.calibrar(wjson("cov_spread.json", {"labels": LABELS, "cov": Cs.tolist()}), "t7c") else None, "límite de la regla (0.9/1.1 -> spread 0.2)")
# 8. modify a copy of the script (SEED_EVAL) and run calibrar+evaluar with it: do the checks notice?
src = (AUD / "pronostico_dr3_lrg2_v2.py").read_text(encoding="utf-8").replace("SEED_EVAL = 500", "SEED_EVAL = 501")
mod_path = WORK / "pronostico_dr3_lrg2_v2_MODIFICADO.py"; mod_path.write_text(src, encoding="utf-8")
spec2 = importlib.util.spec_from_file_location("v2mod2", mod_path)
try:
    V2b = importlib.util.module_from_spec(spec2); spec2.loader.exec_module(V2b); V2b.HERE = WORK
    attempt("8a_calibrar_con_script_modificado", lambda: json.loads((WORK / "pronostico_dr3_lrg2_v2_calibracion_t8.json").read_text())["script_sha256"] if V2b.calibrar(cov_ok, "t8") else None, "se ejecuta; registra OTRO sha256 (no el registrado add948…)")
    attempt("8b_evaluar_con_script_modificado", lambda: json.load(open(V2b.evaluar(cov_ok, vals, "t8")))["result"], "se ejecuta: el script no se compara con el sha registrado")
    attempt("8c_evaluar_t8_con_script_original", lambda: V2.evaluar(cov_ok, vals, "t8b") if (WORK / "pronostico_dr3_lrg2_v2_calibracion_t8b.json").write_text((WORK / "pronostico_dr3_lrg2_v2_calibracion_t8.json").read_text()) and (WORK / "pronostico_dr3_lrg2_v2_calibracion_t8b.sha256").write_text(f"{sha(WORK / 'pronostico_dr3_lrg2_v2_calibracion_t8b.json')}  x\n") else None, "negarse: el script ha cambiado")
except Exception as e:
    res["attempts"]["8_script_modificado"] = {"outcome": "error", "detail": repr(e), "tb": traceback.format_exc()[-800:]}
# 9. does calibrar read anything but the covariance? (static: it opens only cov_path) — recorded by inspection
res["attempts"]["9_calibrar_solo_covarianza"] = {"outcome": "por inspección del código: calibrar abre solo cov_path (read_dr3_cov); evaluar abre cov, valores y calibración", "detail": ""}
res["attempts"]["10_borrar_resultado"] = {"outcome": "no probado (prohibido borrar): borrar pronostico_dr3_lrg2_v2_resultado_TAG.json permitiría reevaluar; el script no puede impedirlo", "detail": ""}
res["work_files"] = sorted(p.name for p in WORK.iterdir())
(AUD / "revision2_pronostico_romper.json").write_text(json.dumps(res, indent=1, default=str), encoding="utf-8")
print("escrito revision2_pronostico_romper.json; ficheros de trabajo en", WORK)
