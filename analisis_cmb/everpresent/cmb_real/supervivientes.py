#!/usr/bin/env python3
"""
Exploration only (EXPLORACION.md, point 1): every Everpresent realisation that passes A or B under any
criterion of EVERPRESENT_vs_DESI.md 3a, 3d or 3e, with its index/seed and what of its history is stored.

Reads (does not write) ../resultados/: bao_m2.json, robustez.json, rd_real_m2.npz, rd_real_m1_alpha*.npz,
bao_m1_alpha*.npz. Recomputes the 3d Model 1 chi^2 (omega_m free, r_d by formula) with robustez.py's
functions, and regenerates the Model 2 Omega_hat on its 401 + 57 sampled points with the same seeds as
bao_sim.py / rd_real.py, to report the sign of Omega_de there.

Writes cmb_real/supervivientes.json.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
EV = HERE.parent
sys.path.insert(0, str(EV))

import bao_sim as B  # noqa: E402
import robustez as RB  # noqa: E402

RES = EV / "resultados"
A_THR, B_THR = B.CHI2_LCDM, B.CHI2_LCDM + 4.0


def model2() -> list[dict]:
    rob = json.loads((RES / "robustez.json").read_text())
    passB_3d = {d["index"]: d for d in rob["m2_passing_B"]}  # exact re-minimisation, 3d
    rr = np.load(RES / "rd_real_m2.npz")
    chi3a = np.load(RES / "bao_m2_bestfit.npz")["chi2"]
    crit = {
        "3a_fixed": chi3a,
        "3e_a_rdreal_wm_fixed": rr["chi_a"],
        "3e_b_rdreal_wm_free": rr["chi_b"],
        "3e_c_rdreal_wm_2sigma_nowindow": rr["chi_c"],
    }
    idx = set(np.nonzero(chi3a < B_THR)[0]) | set(passB_3d)
    for k in ("chi_a", "chi_b", "chi_c"):
        idx |= set(np.nonzero(rr[k] < B_THR)[0])
    idx = sorted(int(i) for i in idx)

    # the same Omega_hat as bao_sim.py (401 points) and rd_real.py (57 points z in [1000, 1e5])
    mu, alpha = 0.9804, 0.8824
    L = B.chol_m2(mu)
    Z = np.random.default_rng(B.SEED).standard_normal((401, 20000))
    Ohat = alpha * L @ Z
    lam_rec = B.LAM_ALL[400]
    zw = np.geomspace(1000.0, 5000.0, 10)
    lam_new = np.unique(np.concatenate([np.arange(-np.log(1e5), -np.log(1000.0) + 1e-9, 0.1), -np.log1p(zw)]))
    k = np.exp(-mu * (lam_new - lam_rec) ** 2 / 2)
    Knn = np.exp(-mu * (lam_new[:, None] - lam_new[None, :]) ** 2 / 2) - np.outer(k, k)
    w, V = np.linalg.eigh(Knn)
    Ls = V * np.sqrt(np.clip(w, 0, None))[None, :]
    Zn = np.random.default_rng(B.SEED + 1).standard_normal((len(lam_new), 20000))
    Onew = alpha * (k[:, None] * (Ohat[400][None, :] / alpha) + Ls @ Zn)

    out = []
    for j in idx:
        low = np.tanh(Ohat[:400, j])          # z in [0, 3]
        high = np.tanh(Onew[:, j])            # z in [1000, 1e5]
        rec = float(np.tanh(Ohat[400, j]))
        row = {"index": j, "Ode0": float(low[-1]), "Ode_rec": rec,
               "Ode_min_z0_3": float(low.min()), "Ode_max_z0_3": float(low.max()),
               "Ode_min_z1e3_1e5": float(high.min()), "Ode_max_z1e3_1e5": float(high.max()),
               "positive_on_sampled_points": bool(low.min() > 0 and high.min() > 0 and rec > 0),
               "window_proxy": bool(rr["proxy_win"][j]), "passes": {}}
        for name, chi in crit.items():
            c = float(chi[j])
            if c < B_THR:
                row["passes"][name] = {"chi2": c, "A": c <= A_THR, "B": True}
        if j in passB_3d:
            d = passB_3d[j]
            row["passes"]["3d_wm_free"] = {"chi2": d["chi2"], "A": d["chi2"] <= A_THR, "B": True,
                                           "omega_m": d["omega_m"], "proxy_rec": d["cmb_proxy"]}
        if "3e_b_rdreal_wm_free" in row["passes"]:
            row["passes"]["3e_b_rdreal_wm_free"]["omega_m"] = float(rr["wm_b"][j])
        if "3e_c_rdreal_wm_2sigma_nowindow" in row["passes"]:
            row["passes"]["3e_c_rdreal_wm_2sigma_nowindow"]["omega_m"] = float(rr["wm_c"][j])
        out.append(row)
    return out


def model1(alpha: float) -> list[dict]:
    R = np.load(RES / f"bao_m1_alpha{alpha}.npz")
    rr = np.load(RES / f"rd_real_m1_alpha{alpha}.npz")
    t = RB.theory13(R["dmh"], R["eh"], R["H0"], R["z_unique"])
    _, chi_3d, wm_3d = RB.chi2_free_wm(t)
    crit = {"3c_fixed": R["chi2"], "3d_wm_free": chi_3d, "3e_a_rdreal_wm_fixed": rr["chi_a"],
            "3e_b_rdreal_wm_free": rr["chi_b"], "3e_c_rdreal_wm_2sigma_nowindow": rr["chi_c"]}
    # 3d is scaled-only for Model 1 (no exact re-minimisation): include runs within 1 of B as borderline
    pos = set()
    for name, chi in crit.items():
        thr = B_THR + (1.0 if name == "3d_wm_free" else 0.0)
        pos |= set(np.nonzero(chi < thr)[0])
    zw = R["z_win"]
    inwin = (zw >= 1000.0) & (zw <= 5000.0)
    out = []
    for p in sorted(int(q) for q in pos):
        row = {"npz_row": p, "seed_chunk": int(R["seed_chunk"][p]), "idx_in_chunk": int(R["idx"][p]),
               "chunk": int(R["seed_chunk"][p] - B.SEED), "H0_wm0": float(R["H0"][p]), "OmL0": float(R["OmL0"][p]),
               "Ode_rec": float(R["Ode_rec"][p]),
               "Ode_min_window_1e3_1e4": float(np.nanmin(R["win_ode"][p])),
               "Ode_max_window_1e3_1e4": float(np.nanmax(R["win_ode"][p])),
               "window_proxy": bool(rr["proxy_win"][p]), "passes": {}}
        for name, chi in crit.items():
            c = float(chi[p])
            if c < B_THR:
                row["passes"][name] = {"chi2": c, "A": c <= A_THR, "B": True}
            elif name == "3d_wm_free" and c < B_THR + 1.0:
                row["passes"][name] = {"chi2": c, "A": False, "B": False, "borderline_pm1": True}
        if "3d_wm_free" in row["passes"]:
            row["passes"]["3d_wm_free"]["omega_m"] = float(wm_3d[p])
        for name, key in (("3e_b_rdreal_wm_free", "wm_b"), ("3e_c_rdreal_wm_2sigma_nowindow", "wm_c")):
            if name in row["passes"]:
                row["passes"][name]["omega_m"] = float(rr[key][p])
        out.append(row)
    return out


def main() -> None:
    m2 = model2()
    res = {"criteria": {"A": A_THR, "B_lt": B_THR},
           "model2": m2,
           "model2_counts": {
               "n_union": len(m2),
               "n_positive_on_sampled_points": sum(r["positive_on_sampled_points"] for r in m2),
               "n_Ode0_negative": sum(r["Ode0"] < 0 for r in m2)},
           "model1_alpha0.01": model1(0.01), "model1_alpha0.005": model1(0.005)}
    (HERE / "supervivientes.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    print(json.dumps(res["model2_counts"]))
    for r in m2:
        print(r["index"], f"Ode0 {r['Ode0']:+.2f} rec {r['Ode_rec']:+.2f} min03 {r['Ode_min_z0_3']:+.2f} "
              f"minhi {r['Ode_min_z1e3_1e5']:+.2f} pos {r['positive_on_sampled_points']}",
              {k: round(v["chi2"], 2) for k, v in r["passes"].items()})
    for a in ("0.01", "0.005"):
        for r in res[f"model1_alpha{a}"]:
            print(a, r["npz_row"], r["chunk"], r["idx_in_chunk"], f"H0 {r['H0_wm0']:.1f} OmL0 {r['OmL0']:.2f}",
                  {k: round(v["chi2"], 2) for k, v in r["passes"].items()})


if __name__ == "__main__":
    main()
