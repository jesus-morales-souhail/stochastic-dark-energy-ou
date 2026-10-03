#!/usr/bin/env python3
"""
Statistics for one resultados/<halofit>/sigma_<s>/ folder -> resumen.json, dchi2_*.txt, tabla.md.

Delta chi^2_r = chi^2(OU_r) - chi^2(LCDM) on Planck bins (TT: upper error if model > data,
lower otherwise; phiphi: bandpower Error column). Every curve has H0 re-solved so that
100 theta_* equals the LCDM reference. Realisations that failed in an observable are excluded
from that observable's percentiles and counted in the table.

Usage: ./venv/bin/python estadisticos.py resultados/takahashi/sigma_2.506e-02
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

import planck_datos as P

THETAS = ("0.001", "1", "20")
DELTAS = ("0.05", "0.10")

# Planck 2018 VI (arXiv:1807.06209) Table 2, TT,TE,EE+lowE+lensing
H0_PLANCK, H0_PLANCK_ERR = 67.36, 0.54
# DESI DR2 (arXiv:2503.14738) Sec. "Cosmological constraints in the LCDM model", DESI+BBN (LCDM).
# Comparison only: the BAO already entered the sigma_X fit.
H0_DESI_BBN, H0_DESI_BBN_ERR = 68.51, 0.58


def targets():
    tt, lens = P.load_tt(), P.load_lensing()
    return {"TT": (tt["Dl"], tt["err_lo"], tt["err_hi"]), "pp": (lens["pp"], lens["err"], lens["err"])}


def chi2(model, data, elo, ehi):
    d = np.asarray(model) - data
    return np.sum((d / np.where(d > 0, ehi, elo)) ** 2, axis=-1)


def band_gap_edge(vals, ref, elo, ehi):
    lo, hi = np.percentile(vals, [16, 84], axis=0)
    gap = np.where(lo > ref, (lo - ref) / ehi, np.where(hi < ref, (ref - hi) / elo, 0.0))
    edge = np.maximum((hi - ref) / ehi, (ref - lo) / elo)
    return float(np.max(gap)), float(np.max(edge)), lo, hi


def resumen(folder: Path) -> dict:
    folder = Path(folder)
    meta = json.loads((folder / "meta.json").read_text())
    refs = meta["references"]
    tg = targets()
    lcdm = refs["LCDM"]
    out = {"meta": meta, "H0_planck": [H0_PLANCK, H0_PLANCK_ERR], "H0_desi_bbn": [H0_DESI_BBN, H0_DESI_BBN_ERR],
           "references": {}, "cases": {}, "delta_sensitivity": {}}

    chi_ref = {o: {k: float(chi2(v[o], *tg[o])) for k, v in refs.items()} for o in tg}
    for k, v in refs.items():
        out["references"][k] = {"H0": v["H0"], "thetastar100_at_H0planck": v["thetastar_H0planck"],
                                "fail": v["fail"],
                                **{f"dchi2_{o}_vs_LCDM": chi_ref[o][k] - chi_ref[o]["LCDM"] for o in tg}}

    bands = {}
    for th in THETAS:
        for de in DELTAS:
            tag = f"theta{th}_Delta{de}"
            f = folder / f"observables_{tag}.npz"
            if not f.exists():
                continue
            d = np.load(f)
            fl = json.loads((folder / f"fallos_{tag}.json").read_text())["per_realisation"]
            case = {"n_rho_le_0": int(d["rho_bad"].sum()), "obs": {}}
            H0 = d["H0"]
            okH = np.isfinite(H0)
            pc = np.percentile(H0[okH], [2.5, 16, 50, 84, 97.5]).tolist()
            case["H0"] = {
                "n_fail": int((~okH).sum()) - case["n_rho_le_0"],
                "p2.5_16_50_84_97.5": pc,
                "frac_outside_planck_2sigma": float(np.mean(np.abs(H0[okH] - H0_PLANCK) > 2 * H0_PLANCK_ERR)),
                "frac_outside_desi_bbn_2sigma": float(np.mean(np.abs(H0[okH] - H0_DESI_BBN) > 2 * H0_DESI_BBN_ERR)),
                "thetastar100_at_H0planck_p16_50_84":
                    np.percentile(d["thetastar_H0planck"][np.isfinite(d["thetastar_H0planck"])], [16, 50, 84]).tolist(),
            }
            bands[(th, de, "H0")] = (pc[1], pc[3])
            for o in tg:
                v = d[o]
                ok = np.all(np.isfinite(v), axis=1)
                n_camb = sum(1 for m in fl.values() if o in m and not m[o].startswith("silent"))
                n_nan = sum(1 for m in fl.values() if o in m and m[o].startswith("silent"))
                dchi = np.full(len(v), np.nan)
                dchi[ok] = chi2(v[ok], *tg[o]) - chi_ref[o]["LCDM"]
                np.savetxt(folder / f"dchi2_{o}_{tag}.txt", dchi,
                           header=f"Delta chi2_r = chi2(OU_r) - chi2(LCDM), {o}, {tag}, sigma_X={meta['sigma_X']:.6e}, "
                                  f"H0 re-solved for theta_*; NaN = failure (see fallos_{tag}.json)")
                x = dchi[ok]
                gap, edge, lo, hi = band_gap_edge(v[ok], np.asarray(lcdm[o]), tg[o][1], tg[o][2])
                bands[(th, de, o)] = (lo, hi)
                chi_r = chi2(v[ok], *tg[o])
                case["obs"][o] = {
                    "n_valid": int(ok.sum()), "n_fail_camb": n_camb, "n_fail_silent_nan": n_nan,
                    "dchi2_p50": float(np.percentile(x, 50)), "dchi2_p95": float(np.percentile(x, 95)),
                    "frac_dchi2_gt_4": float(np.mean(x > 4)),
                    "band68_gap_vs_LCDM_sigma": gap, "band68_edge_vs_LCDM_sigma": edge,
                    "dchi2_vs_CPL_repo_p95": float(np.percentile(chi_r - chi_ref[o]["CPL_repo"], 95)),
                    "dchi2_vs_CPL_joint_p95": float(np.percentile(chi_r - chi_ref[o]["CPL_joint"], 95)),
                }
            out["cases"][tag] = case

    sig = {"TT": 0.5 * (tg["TT"][1] + tg["TT"][2]), "pp": tg["pp"][1], "H0": np.array([H0_PLANCK_ERR])}
    for th in THETAS:
        a, b = f"theta{th}_Delta0.05", f"theta{th}_Delta0.10"
        if a not in out["cases"] or b not in out["cases"]:
            continue
        e = {}
        for o in ("TT", "pp", "H0"):
            la, ha = (np.atleast_1d(x) for x in bands[(th, "0.05", o)])
            lb, hb = (np.atleast_1d(x) for x in bands[(th, "0.10", o)])
            e[o] = {"band68_diff_sigma": float(max(np.max(np.abs(la - lb) / sig[o]), np.max(np.abs(ha - hb) / sig[o])))}
            if o != "H0":
                e[o]["dchi2_p95_diff"] = out["cases"][a]["obs"][o]["dchi2_p95"] - out["cases"][b]["obs"][o]["dchi2_p95"]
        out["delta_sensitivity"][f"theta{th}"] = e

    (folder / "resumen.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    (folder / "tabla.md").write_text(tabla(out), encoding="utf-8")
    return out


def g(x, nd=3):
    return f"{x:.{nd}g}"


def tabla(s: dict) -> str:
    m = s["meta"]
    L = [f"# σ_X = {m['sigma_X']:.4e} — {m['n_real']} realizaciones, semilla {m['seed']}, lensing no lineal: {m['halofit']}",
         "",
         f"H0 de cada curva resuelto para 100θ_* = {m['thetastar100_target_LCDM']:.6f} (ΛCDM de referencia), ω_b, ω_c, mnu fijos, plano.",
         "Δχ² = χ²(curva) − χ²(ΛCDM) frente a datos Planck. Fuente: `resumen.json`, `dchi2_*.txt`, `fallos_*.json` de esta carpeta.",
         "",
         "## Curvas de referencia", "",
         "| curva | H0 resuelto | 100θ_* con H0 = 67.36 | Δχ² TT 2≤l<30 | Δχ² φφ L 8–400 | fallos |",
         "|---|---|---|---|---|---|"]
    for k, v in s["references"].items():
        L.append(f"| {k} | {v['H0']:.3f} | {v['thetastar100_at_H0planck']:.5f} | {g(v['dchi2_TT_vs_LCDM'])} | "
                 f"{g(v['dchi2_pp_vs_LCDM'])} | {len(v['fail']) or 0} |")
    L += ["", "Ficheros Planck: TT `COM_PowerSpect_CMB-TT-full_R3.01.txt`; φφ `smicadx12_Dec5_ftl_mv2_ndclpp_p_teb_consext8_bandpowers.dat`.",
          "", "## OU: TT y φφ (Δχ² frente a ΛCDM)", "",
          "| obs. | θ | Δ | ρ≤0 | fallos (CAMB+NaN) | válidas | Δχ² p50 | Δχ² p95 | frac >4 | hueco banda68 [σ] | borde banda68 [σ] | p95 vs CPL repo | p95 vs CPL joint | nota |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for o, name in (("TT", "TT l<30"), ("pp", "φφ")):
        for tag, c in s["cases"].items():
            th, de = tag.split("_")[0][5:], tag.split("Delta")[1]
            e = c["obs"][o]
            note = {"0.001": "θ de DESI a esta σ_X (borde de malla)", "20": "dominado por el suavizado (1/θ = Δ)"}.get(th, "")
            L.append(f"| {name} | {th} | {de} | {c['n_rho_le_0']} | {e['n_fail_camb']}+{e['n_fail_silent_nan']} | "
                     f"{e['n_valid']} | {g(e['dchi2_p50'])} | {g(e['dchi2_p95'])} | {e['frac_dchi2_gt_4']:.3f} | "
                     f"{g(e['band68_gap_vs_LCDM_sigma'])} | {g(e['band68_edge_vs_LCDM_sigma'])} | "
                     f"{g(e['dchi2_vs_CPL_repo_p95'])} | {g(e['dchi2_vs_CPL_joint_p95'])} | {note} |")
    L += ["", f"## OU: H0 resuelto (Planck {H0_PLANCK} ± {H0_PLANCK_ERR}; DESI+BBN {H0_DESI_BBN} ± {H0_DESI_BBN_ERR}, sólo comparación)", "",
          "| θ | Δ | fallos H0 | p2.5 | p16 | p50 | p84 | p97.5 | frac fuera Planck ±2σ | frac fuera DESI+BBN ±2σ | 100θ_* con H0=67.36 [p16, p50, p84] |",
          "|---|---|---|---|---|---|---|---|---|---|---|"]
    for tag, c in s["cases"].items():
        th, de = tag.split("_")[0][5:], tag.split("Delta")[1]
        h = c["H0"]
        p = h["p2.5_16_50_84_97.5"]
        t = h["thetastar100_at_H0planck_p16_50_84"]
        L.append(f"| {th} | {de} | {h['n_fail']} | " + " | ".join(f"{x:.2f}" for x in p)
                 + f" | {h['frac_outside_planck_2sigma']:.3f} | {h['frac_outside_desi_bbn_2sigma']:.3f} | "
                   f"[{t[0]:.5f}, {t[1]:.5f}, {t[2]:.5f}] |")
    L += ["", "## Sensibilidad a Δ (0.05 frente a 0.10)", "",
          "| θ | obs. | Δ banda68 [σ_Planck] | Δ p95 de Δχ² |", "|---|---|---|---|"]
    for th, e in s["delta_sensitivity"].items():
        for o in ("TT", "pp", "H0"):
            L.append(f"| {th[5:]} | {o} | {g(e[o]['band68_diff_sigma'])} | {g(e[o]['dchi2_p95_diff']) if 'dchi2_p95_diff' in e[o] else '—'} |")
    L.append("")
    L.append("σ de H0 en la sensibilidad a Δ: error de Planck (0.54).")
    return "\n".join(L) + "\n"


if __name__ == "__main__":
    print(tabla(resumen(Path(sys.argv[1]))))
