#!/usr/bin/env python3
"""
Our tabulated r_d(omega_m), r_s(z_*; omega_m) and z_*(omega_m) (escalon.RD_SPL, RS_SPL, zstar_of, evaluated exactly as the code
does) against CAMB, with the fixed parameters of everpresent/rd_planck.py (ou_cmb.PLANCK18: H0 = 67.36, ombh2 = 0.02237,
mnu = 0.06 eV, tau = 0.0544, ns = 0.9649, As = 2.1e-9, flat, DarkEnergyPPF with w = -1; CAMB default accuracy, T_CMB and Y_He
defaults, N_eff default) varying only omch2 = omega_m - ombh2 - omnuh2. Writes tablas_camb_check.json. No table is modified.
"""
import json, sys, platform
from pathlib import Path
import numpy as np
import camb
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent))
import exacto_sn as S  # noqa: E402
import ou_cmb as M  # noqa: E402
X, E = S.X, S.E
p = M.PLANCK18
base = M.make_params(halofit="lineal")
omnuh2 = float(base.omnuh2)
def run(wm):
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=p["H0"], ombh2=p["ombh2"], omch2=wm - p["ombh2"] - omnuh2, tau=p["tau"], mnu=p["mnu"], omk=0.0)
    pars.InitPower.set_params(As=p["As"], ns=p["ns"])
    from camb.dark_energy import DarkEnergyPPF
    de = DarkEnergyPPF(); de.set_params(w=-1.0, wa=0.0); pars.DarkEnergy = de
    d = camb.get_background(pars).get_derived_params()
    return pars, {"rdrag": float(d["rdrag"]), "zdrag": float(d["zdrag"]), "zstar": float(d["zstar"]), "rstar": float(d["rstar"]), "thetastar100": float(d["thetastar"])}
wm_list = [0.12, 0.13, 0.14, 0.1409233, 0.15, 0.16, float(X.WM_P)] + [float(w) for w in np.arange(0.12, 0.1601, 0.005)]
wm_list = sorted(set(round(w, 7) for w in wm_list))
rows = []
for wm in wm_list:
    pars, d = run(wm)
    t = {"rd": float(E.RD_SPL(wm)), "rs": float(E.RS_SPL(wm)), "zstar": float(E.zstar_of(wm))}
    rows.append({"omega_m": wm, "omch2": float(pars.omch2), "camb": d, "table": t,
                 "rel_diff_table_minus_camb": {"rd": t["rd"] / d["rdrag"] - 1, "rs_zstar": t["rs"] / d["rstar"] - 1, "zstar": t["zstar"] / d["zstar"] - 1}})
pars0, d0 = run(float(X.WM_P))
fixed = {"H0": float(pars0.H0), "ombh2": float(pars0.ombh2), "omnuh2": omnuh2, "mnu_eV": p["mnu"], "tau": p["tau"], "ns": p["ns"], "As": p["As"], "omk": 0.0,
         "TCMB_K": float(pars0.TCMB), "YHe": float(pars0.YHe), "N_eff_massless": float(pars0.num_nu_massless), "num_nu_massive": int(pars0.num_nu_massive), "nu_mass_eigenstates": int(pars0.nu_mass_eigenstates),
         "dark_energy": "DarkEnergyPPF w = -1, wa = 0", "accuracy": {"AccuracyBoost": float(pars0.Accuracy.AccuracyBoost), "lAccuracyBoost": float(pars0.Accuracy.lAccuracyBoost), "lSampleBoost": float(pars0.Accuracy.lSampleBoost), "note": "CAMB defaults (get_background); same call as everpresent/rd_planck.py"},
         "recomb": str(type(pars0.Recomb).__name__), "nonlinear": "none"}
rng = [r for r in rows if 0.12 - 1e-9 <= r["omega_m"] <= 0.16 + 1e-9]
summary = {"planck_point": next(r for r in rows if abs(r["omega_m"] - float(X.WM_P)) < 1e-9)["rel_diff_table_minus_camb"],
           "max_abs_rel_diff_in_0.12_0.16": {k: float(max(abs(r["rel_diff_table_minus_camb"][k]) for r in rng)) for k in ("rd", "rs_zstar", "zstar")},
           "max_abs_rel_diff_in_0.13_0.15": {k: float(max(abs(r["rel_diff_table_minus_camb"][k]) for r in rng if 0.13 - 1e-9 <= r["omega_m"] <= 0.15 + 1e-9)) for k in ("rd", "rs_zstar", "zstar")},
           "reference_values_rd_planck_json": {"rdrag_Mpc": S.E.RDJ["rdrag_Mpc"], "zstar": S.E.RDJ["zstar"], "rs_zstar_Mpc": S.E.RDJ["rs_zstar_Mpc"], "omega_m_incl_nu": S.E.RDJ["omega_m_incl_nu"]},
           "camb_at_planck_here": d0}
summary["P3"] = bool(abs(summary["planck_point"]["rd"]) < 0.003 and summary["max_abs_rel_diff_in_0.12_0.16"]["rd"] < 0.01)
summary["P4"] = bool(summary["max_abs_rel_diff_in_0.12_0.16"]["zstar"] < 0.005)
out = {"versions": {"camb": camb.__version__, "python": platform.python_version(), "numpy": np.__version__}, "fixed_parameters": fixed, "rows": rows, "summary": summary,
       "table_origin": "escalon.py: RD_SPL = spline of NORM * int_{z_d(wm)}^{inf} c_s/H dz (no dark energy), z_d(wm) = z_d^CAMB * EH98(wm)/EH98(wm_Planck), NORM = r_d^CAMB / integral at wm_Planck; RS_SPL same integral from z_*(wm) with the same NORM; z_*(wm) = z_*^CAMB * HS96(wm)/HS96(wm_Planck); CAMB values from everpresent/resultados/rd_planck.json"}
(HERE / "tablas_camb_check.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
print("camb", camb.__version__, "omnuh2", omnuh2); print(json.dumps(fixed, indent=1))
for r in rows: print(f"wm {r['omega_m']:<9g} rd tab {r['table']['rd']:.4f} camb {r['camb']['rdrag']:.4f} ({100*r['rel_diff_table_minus_camb']['rd']:+.3f} %) | rs* tab {r['table']['rs']:.4f} camb {r['camb']['rstar']:.4f} ({100*r['rel_diff_table_minus_camb']['rs_zstar']:+.3f} %) | z* tab {r['table']['zstar']:.3f} camb {r['camb']['zstar']:.3f} ({100*r['rel_diff_table_minus_camb']['zstar']:+.4f} %)")
print(json.dumps(summary, indent=1))
