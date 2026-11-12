#!/usr/bin/env python3
"""
Sweep (BARRIDO_MAXIMOS.md): every json of results/ and analisis_cmb/ that stores a profile or a bound measured as a drop from a
maximum. For each entry: grid (min, max), where the maximum is, class (interior / physical lower edge sigma_X = 0 / artificial
upper edge = failure), and, for failures or edge-of-theta cases, the bound recomputed with the grid extended twice (the maximum
must not move in the second extension). Nothing historical is modified. Writes barrido_maximos.json.
Pieces:
  A. OU profile of the 7 alpha (results/profile_sigma_x, perfil_fondo (a), reloj clocks, desi_dr2_real_bao, everpresent T0):
     sigma grid 1e-6..5e-2, theta (1e-3, 20); theta at the bound = 1e-3 (edge) -> theta range extended to (1e-5, 200) and (1e-7, 2000).
  B. perfil_fondo (b), (c) and perfil_13 (a, b, c) (+ cuesta's theta = 0.3, 3): sigma grid extended to 1 for the maximum; theta
     range extended for theta free.
  C. verificacion.py full-response kernel (7 alpha and 13, local and full, a/b/c, theta fixed and free): sigma extended to 1.
  D. perfil_sn: sigma extension already stored in the json (monotone decrease beyond 0.05); theta free at the bound = 1e-3 -> extended.
  E. everpresent Z2 (Zwane, mu free): mu at the bound = 0.1 (edge of [0.1, 10]) -> mu range extended to [0.01, 10] and [0.001, 10].
  F. escalon cells with the best fit or the 95 % upper edge at f = 3: f extended to 6 and 12.
  G. read-only classification of the rest (fisher equaliser, precision/convergencia/marginal_guiada, exacto P-crossings, tachyonic,
     ou_qnm, joint_w0wa, malla, gof, cmb_fondo and exacto_sn from maximos_globales.json).
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, minimize_scalar

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "fondo"))
sys.path.insert(0, str(HERE.parent / "reloj"))
sys.path.insert(0, str(HERE.parent / "escalon"))
sys.path.insert(0, str(HERE.parent / "everpresent"))
import perfil_13 as P13  # noqa: E402
import perfil_fondo as PF  # noqa: E402
import profile_sigma_x_desi as PS  # noqa: E402
import verificacion as V  # noqa: E402

SIG = PF.SIG
EXT = np.geomspace(SIG[-1], 1.0, 40)[1:]
GRID = np.concatenate([SIG, EXT])
THR = -1.92
TH_RANGES = {"base (1e-3, 20)": (1e-3, 20.0), "ext1 (1e-5, 200)": (1e-5, 200.0), "ext2 (1e-7, 2000)": (1e-7, 2000.0)}
OUT = {"rows": [], "details": {}}


def max_theta_ext(fun, sigma, th_min, th_max):
    best_ll, best_th = -1e30, 1.0
    brs = [(th_min, 1e-3), (1e-3, 0.5), (0.1, 3.0), (1.0, 15.0), (15.0, th_max)] if th_min < 1e-3 else [(1e-3, 0.5), (0.1, 3.0), (1.0, 15.0)]
    if th_max > 20 and th_min >= 1e-3:
        brs.append((15.0, th_max))
    for br in brs:
        if br[0] >= br[1]:
            continue
        r = minimize_scalar(lambda th: -fun(sigma, th), bounds=br, method="bounded", options={"xatol": 1e-5})
        if -r.fun > best_ll:
            best_ll, best_th = -float(r.fun), float(r.x)
    for th in np.geomspace(th_min, th_max, 60):
        v = fun(sigma, th)
        if v > best_ll:
            best_ll, best_th = v, float(th)
    return best_ll, best_th


def bound_global(g, thr=THR, grid=GRID):
    """Profile on grid; max over the whole grid; bound = crossing of thr below the max, outward from the max."""
    ll = np.array([g(float(s)) for s in grid])
    k0 = int(np.argmax(ll))
    llmax = float(ll[k0])
    below = np.nonzero((ll - llmax < thr) & (np.arange(len(ll)) > k0))[0]
    b = float(brentq(lambda s: g(float(s)) - llmax - thr, grid[below[0] - 1], grid[below[0]], xtol=1e-10 * grid[below[0] - 1])) if len(below) else None
    cls = "sigma->0" if k0 == 0 else ("interior" if k0 < len(grid) - 1 else "upper edge (artificial)")
    return {"sigma_at_max": float(grid[k0]), "lnL_max_rel_sigma0": float(llmax - ll[0]), "max_class": cls,
            "clipped_max_at_0.05": bool(int(np.argmax(ll[: len(SIG)])) == len(SIG) - 1 and k0 >= len(SIG) - 1), "bound": b,
            "bound_class": ("none up to 1" if b is None else ("beyond 0.05 (extension)" if b > SIG[-1] else "inside the historical grid"))}


def row(file, quantity, grid, max_at, cls, hist, recomputed=None, note=""):
    OUT["rows"].append({"file": file, "quantity": quantity, "grid": grid, "max_at": max_at, "class": cls, "historical": hist, "recomputed": recomputed, "note": note})


def main():
    t0 = time.time()
    # ------------------------------------------------------------------ A. OU 7-alpha profile, theta range
    ref = json.loads((ROOT / "results/profile_sigma_x/profile_sigma_x.json").read_text())
    detA = {}
    for lab, (tmin, tmax) in TH_RANGES.items():
        g = lambda s, tmin=tmin, tmax=tmax: max_theta_ext(PF.ll_a, s, tmin, tmax)[0]
        ll = np.array([g(float(s)) for s in SIG])
        llmax = ll.max()
        j = int(np.nonzero(ll - llmax >= THR)[0][-1])
        b = float(brentq(lambda s: g(float(s)) - llmax - THR, SIG[j], SIG[j + 1], xtol=1e-10 * SIG[j]))
        th_b = max_theta_ext(PF.ll_a, b, tmin, tmax)[1]
        detA[lab] = {"bound": b, "theta_at_bound": th_b, "sigma_at_max": float(SIG[int(np.argmax(ll))]), "theta_at_max": max_theta_ext(PF.ll_a, float(SIG[int(np.argmax(ll))]), tmin, tmax)[1]}
        print("A", lab, detA[lab], f"[{time.time() - t0:.0f} s]", flush=True)
    # fully correlated limit (rank-1): the theta -> 0 supremum
    g0 = lambda s: PF.ll_a(s, 1e-12)
    ll0 = np.array([g0(float(s)) for s in SIG]); j = int(np.nonzero(ll0 - ll0.max() >= THR)[0][-1])
    detA["theta->0 (1e-12)"] = {"bound": float(brentq(lambda s: g0(float(s)) - ll0.max() - THR, SIG[j], SIG[j + 1], xtol=1e-10 * SIG[j]))}
    OUT["details"]["A_profile_sigma_x_theta_range"] = detA
    b0, b1, b2 = detA["base (1e-3, 20)"]["bound"], detA["ext1 (1e-5, 200)"]["bound"], detA["ext2 (1e-7, 2000)"]["bound"]
    note = f"theta at the bound = 1e-3 (lower edge); extending theta: {b0:.5f} -> {b1:.5f} ({100*(b1/b0-1):+.2f} %) -> {b2:.5f} ({100*(b2/b0-1):+.2f} %); theta->0 supremum {detA['theta->0 (1e-12)']['bound']:.5f}; max at sigma = 1e-6 where theta is flat (theta_best = 20 there is irrelevant)"
    row("results/profile_sigma_x/profile_sigma_x.json", "sigma_X 95 % (theta free)", "sigma 1e-6..5e-2 (48), theta (1e-3, 20)", "sigma = 1e-6", "sigma->0; theta at lower edge of its range", ref["sigma_X_95CL_upper_profile"], {"ext1": b1, "ext2": b2}, note)
    for f, q, h in (("results/desi_dr2_real_bao/desi_dr2_real_bao.json", "sigma_X_95_profile_upper_approx (old grid node)", 0.02175257993422772),
                    ("analisis_cmb/everpresent/resultados/profile_desi_everpresent.json", "T0_ou_theta_profiled upper95_refined", 0.0263282756994459)):
        row(f, q, "sigma 1e-6..5e-2, theta (1e-3, 20)", "sigma = 1e-6", "sigma->0; theta at lower edge", h, {"ext1": b1, "ext2": b2}, "same likelihood as profile_sigma_x; same theta extension")
    # clocks (reloj)
    import perfil_reloj as PR
    rel = json.loads((HERE.parent / "reloj/perfil_reloj.json").read_text())
    x_save = PS.x.copy()
    detR = {}
    for name, x in PR.clocks(PS.z_eff).items():
        PS.x = x
        res = {}
        for lab, (tmin, tmax) in TH_RANGES.items():
            g = lambda s, tmin=tmin, tmax=tmax: max_theta_ext(PF.ll_a, s, tmin, tmax)[0]
            ll = np.array([g(float(s)) for s in SIG]); llmax = ll.max(); j = int(np.nonzero(ll - llmax >= THR)[0][-1])
            res[lab] = float(brentq(lambda s: g(float(s)) - llmax - THR, SIG[j], SIG[j + 1], xtol=1e-10 * SIG[j]))
        detR[name] = res
        h = rel["clocks"][name]["sigma_X_95"]
        row("analisis_cmb/reloj/perfil_reloj.json", f"clock {name}: sigma_X 95 % (theta free)", "sigma 1e-6..5e-2, theta (1e-3, 20)", "sigma = 1e-6", "sigma->0; theta at lower edge", h,
            {"ext1": res["ext1 (1e-5, 200)"], "ext2": res["ext2 (1e-7, 2000)"]}, f"base reproduces {res['base (1e-3, 20)']:.5f}; ext changes {100*(res['ext2 (1e-7, 2000)']/h-1):+.2f} %")
        print("R", name, res, flush=True)
    PS.x = x_save
    OUT["details"]["A_reloj"] = detR
    # ------------------------------------------------------------------ B. perfil_fondo / perfil_13 / cuesta: sigma extension + theta extension
    pf = json.loads((HERE.parent / "fondo/perfil_fondo.json").read_text())["treatments"]
    p13 = json.loads((HERE.parent / "fondo/perfil_13.json").read_text())
    D7 = P13.Data(PF.Z, ["DV"] * 7, PF.ALPHA, PF.CM, PF.S)
    D13, _, _, _ = V.data13()
    detB = {}
    for lab, D, hist in (("perfil_fondo (7 alpha)", D7, pf), ("perfil_13 (13 numbers)", D13, p13["thirteen"])):
        detB[lab] = {}
        for tr, fun in (("a_fixed", D.ll_a), ("b_profiled", D.ll_b), ("c_marginalised", D.ll_c)):
            for th in ("free", 0.001, 0.1, 0.3, 1.0, 3.0, 10.0):
                key = "theta_free" if th == "free" else f"theta_{th:g}"
                h = hist[tr].get(key, {}).get("sigma_95", "n/a") if isinstance(hist[tr].get(key, {}), dict) else "n/a"
                if th == "free":
                    rec = {}
                    for rl, (tmin, tmax) in TH_RANGES.items():
                        g = lambda s, f=fun, tmin=tmin, tmax=tmax: max_theta_ext(f, s, tmin, tmax)[0]
                        rec[rl] = bound_global(g)
                    base = rec["base (1e-3, 20)"]
                    rc = {k: v["bound"] for k, v in rec.items()}
                    note = f"theta free: base {rc['base (1e-3, 20)']} -> ext1 {rc['ext1 (1e-5, 200)']} -> ext2 {rc['ext2 (1e-7, 2000)']}"
                else:
                    base = bound_global(lambda s, f=fun, t=th: f(s, t))
                    rc, note = {"global_max": base["bound"]}, ""
                detB[lab][f"{tr}|{key}"] = rec if th == "free" else base
                cls = base["max_class"] + ("; theta at lower edge" if th == "free" else "")
                row(f"analisis_cmb/fondo/{'perfil_fondo' if '7' in lab else 'perfil_13'}.json" + (" (+ cuesta bounds_13)" if (th in (0.3, 3.0) and '13' in lab) else ""),
                    f"{tr} {key}", "sigma 1e-6..5e-2 (ext. to 1 for the crossing), theta (1e-3, 20)", f"sigma = {base['sigma_at_max']:.3g}", cls, h, rc, note)
            print("B", lab, tr, f"[{time.time() - t0:.0f} s]", flush=True)
    OUT["details"]["B_fondo_13"] = detB
    # ------------------------------------------------------------------ C. verificacion full kernel
    ver = json.loads((HERE / "verificacion.json").read_text())["full_kernel"]
    R7, R13 = V.Response(V.Z7, ["DV"] * 7), V.Response(D13.z, D13.kinds)
    detC = {}
    for lab, D, R in (("7alpha", D7, R7), ("13", D13, R13)):
        for kern_lab, kern in (("local", lambda s, t, D=D: D.cov + np.outer(D.S, D.S) * s**2 * np.exp(-t * D.dx)), ("full", lambda s, t, R=R, D=D: D.cov + s**2 * R.M(t))):
            D.kern = kern
            for tr, fun in (("a_fixed", D.ll_a), ("b_profiled", D.ll_b), ("c_marginalised", D.ll_c)):
                for th in (0.001, 0.1, 1.0, 3.0, 10.0, "free"):
                    key = "theta_free" if th == "free" else f"theta_{th:g}"
                    g = (lambda s, f=fun: V.P13.max_theta(f, s)[0]) if th == "free" else (lambda s, f=fun, t=th: f(s, t))
                    b = bound_global(g)
                    h = ver[lab][kern_lab][tr][key]
                    detC[f"{lab}|{kern_lab}|{tr}|{key}"] = b
                    row("analisis_cmb/auditoria_v23/verificacion.json", f"{lab} {kern_lab} {tr} {key}", "sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20)", f"sigma = {b['sigma_at_max']:.3g}",
                        b["max_class"] + ("; theta at lower edge" if th == "free" else ""), h, {"global_max": b["bound"]}, "")
            print("C", lab, kern_lab, f"[{time.time() - t0:.0f} s]", flush=True)
    OUT["details"]["C_verificacion"] = detC
    # ------------------------------------------------------------------ D. SN: theta range at the free-theta bound
    sys.path.insert(0, str(HERE.parent / "sn"))
    import perfil_sn as PSN
    snj = json.loads((HERE.parent / "sn/perfil_sn.json").read_text())
    detD = {}
    for lab, (tmin, tmax) in TH_RANGES.items():
        g = lambda s, tmin=tmin, tmax=tmax: max_theta_ext(lambda sg, th: PSN.max_over_om(sg, th)[0], s, tmin, tmax)[0]
        # profile only near the historical bound 0.1166: coarse grid 0.05..0.3 then brentq; max from the json (sigma_at_max 0.00998, ll_max)
        llmax = g(0.009979836078583298)
        grid = np.geomspace(0.05, 0.3, 12)
        ll = np.array([g(float(s)) for s in grid])
        k = int(np.nonzero(ll - llmax < THR)[0][0])
        b = float(brentq(lambda s: g(float(s)) - llmax - THR, grid[k - 1] if k > 0 else 0.03, grid[k], xtol=1e-6))
        detD[lab] = {"bound": b, "lnL_max": llmax}
        print("D", lab, detD[lab], f"[{time.time() - t0:.0f} s]", flush=True)
    OUT["details"]["D_sn"] = detD
    h = snj["free_theta"]["sigma_95"]
    row("analisis_cmb/sn/perfil_sn.json", "free_theta sigma_95 (Omega_m profiled)", "sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20)", f"sigma = {snj['free_theta']['sigma_at_max']:.4f}", "interior (lnL +0.026); theta at the bound = 1e-3 (lower edge)", h,
        {"ext1": detD["ext1 (1e-5, 200)"]["bound"], "ext2": detD["ext2 (1e-7, 2000)"]["bound"]}, f"base reproduces {detD['base (1e-3, 20)']['bound']:.4f}; the stored extension profile decreases monotonically beyond 0.05 (max interior)")
    for th, v in snj["fixed_theta"].items():
        row("analisis_cmb/sn/perfil_sn.json", f"fixed theta = {th}: sigma_95", "sigma 1e-6..5e-2 (ext. to 1)", f"sigma = {v['sigma_at_max']:.3g}", "sigma->0" if v["sigma_at_max"] < 1e-5 else "interior", v["sigma_95"], None, "stored extension profile decreases monotonically beyond 0.05")
    # ------------------------------------------------------------------ E. everpresent Z2: mu range
    # Zwane kernel of everpresent/desi_dr2_data.add_zwane_kernel (its module clashes with scripts/desi_dr2_data): inline
    _lam = -np.log1p(PS.z_eff)
    _d2 = (_lam[:, None] - _lam[None, :]) ** 2
    _SS = np.outer(PS.S_z, PS.S_z)
    ll_zw = lambda mu, a: PS.logL_gaussian(PS.residuals, PS.C_MEAS + _SS * a**2 * np.exp(-mu * _d2 / 2.0))
    GRID_ZW = np.geomspace(1e-6, 2.0, 80)
    ev = json.loads((HERE.parent / "everpresent/resultados/profile_desi_everpresent.json").read_text())
    detE = {}
    for lab, mu_lo in (("base [0.1, 10]", 0.1), ("ext1 [0.01, 10]", 0.01), ("ext2 [0.001, 10]", 0.001)):
        def prof_mu(a, mu_lo=mu_lo):
            lg = np.linspace(np.log(mu_lo), np.log(10.0), 61)
            lls = [ll_zw(np.exp(v), a) for v in lg]
            k = int(np.argmax(lls))
            lo, hi = lg[max(k - 1, 0)], lg[min(k + 1, len(lg) - 1)]
            r = minimize_scalar(lambda v: -ll_zw(np.exp(v), a), bounds=(lo, hi), method="bounded", options={"xatol": 1e-6})
            return (-r.fun, float(np.exp(r.x))) if -r.fun >= lls[k] else (lls[k], float(np.exp(lg[k])))
        grid = GRID_ZW
        ll = np.array([prof_mu(float(a))[0] for a in grid]); llmax = ll.max(); j = int(np.nonzero(ll - llmax >= THR)[0][-1])
        b = float(brentq(lambda a: prof_mu(float(a))[0] - llmax - THR, grid[j], grid[j + 1], xtol=1e-8))
        detE[lab] = {"bound": b, "mu_at_bound": prof_mu(b)[1], "amp_at_max": float(grid[int(np.argmax(ll))])}
        print("E", lab, detE[lab], f"[{time.time() - t0:.0f} s]", flush=True)
    OUT["details"]["E_zwane_mu"] = detE
    row("analisis_cmb/everpresent/resultados/profile_desi_everpresent.json", "Z2_zwane_mu_free upper95_refined (alpha)", "alpha 1e-6..2 (80), mu [0.1, 10]", "alpha = 1e-6", "sigma->0 (alpha); mu at the bound = 0.1 (lower edge of its range)",
        ev["Z2_zwane_mu_free"]["upper95_refined"], {"ext1": detE["ext1 [0.01, 10]"]["bound"], "ext2": detE["ext2 [0.001, 10]"]["bound"]}, f"base reproduces {detE['base [0.1, 10]']['bound']:.5f}")
    for k in ("T1_ou_theta1", "T2_ou_theta_fit", "Z1_zwane_mu_fixed"):
        v = ev[k]
        row("analisis_cmb/everpresent/resultados/profile_desi_everpresent.json", f"{k} upper95_refined", "amp 1e-6..0.05 or 2 (Z1)", f"amp = {v['best_amp']:g}", "sigma->0; nuisance interior or fixed", v["upper95_refined"], None, f"nuisance at bound {v['nuisance_at_upper95']}")
    # ------------------------------------------------------------------ F. escalon: f edges
    import escalon as E
    esc = json.loads((HERE.parent / "escalon/escalon.json").read_text())
    detF = {}
    for treat, cells in (("A", ["zt_1.6|w_0.3", "zt_2|w_0.3"]), ("B", ["zt_3|w_0.05", "zt_3|w_0.3"])):
        for cname in cells:
            cell = [c for c in E.CELLS if c.name == cname][0]
            fun = E.profile_fun(treat, cell, E.OBS, E.TH_P if treat == "B" else None, E.WM_P if treat == "B" else None)
            c0 = esc["real"][treat]["chi2_0"]
            rec = {}
            for lab, fmax in (("base [-2, 3]", 3.0), ("ext1 [-2, 6]", 6.0), ("ext2 [-2, 12]", 12.0)):
                fg = np.round(np.arange(-2.0, fmax + 1e-9, 0.1), 3)
                vals = np.array([fun(f) for f in fg])
                fb, cb = E.refine(fun, fg, vals)
                iv = E.interval95(fun, fg, vals, fb, cb)
                rec[lab] = {"f_best": fb, "dchi2": cb - c0, "hi95": iv["hi"]["f"], "hi_flag": iv["hi"]["flag"], "n_H2_le_0": int(np.sum(~np.isfinite(vals)))}
            detF[f"{treat}|{cname}"] = rec
            h = esc["real"][treat]["cells"][cname]
            row("analisis_cmb/escalon/escalon.json", f"({treat}) {cname}: f_best, dchi2, 95 % upper", "f in [-2, 3] step 0.05", f"f = {h['f_best']:.3f}", "upper edge of the f grid (artificial)" if h["f_best"] >= 2.99 else "interior; 95 % upper interval at the edge",
                {"f_best": h["f_best"], "dchi2": h["dchi2"], "hi95": h["interval95"]["hi"]["f"]}, {k: {"f_best": v["f_best"], "dchi2": v["dchi2"], "hi95": v["hi95"]} for k, v in rec.items()},
                "global p-value of ESCALON.md not recomputed with the extended f range")
            print("F", treat, cname, rec, f"[{time.time() - t0:.0f} s]", flush=True)
    OUT["details"]["F_escalon"] = detF
    # ------------------------------------------------------------------ G. read-only
    mg = json.loads((HERE / "maximos_globales.json").read_text())
    for name, trs in mg["cmb_fondo"].items():
        for tr, ths in trs.items():
            for key, r in ths.items():
                fail = r["max_class"] != "sigma->0" and r["clipped_max_at_grid_edge_0.05"]
                row("analisis_cmb/auditoria_v23/cmb_fondo.json", f"{name} {tr} {key}", "sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20)", f"sigma = {r['sigma_at_max_grid']:.3g}",
                    ("upper edge of the historical max grid 0.05 (failure)" if fail else r["max_class"]) + ("; theta at lower edge" if key == "theta_free" and not fail else ""),
                    r["192_historical"], {"global_max_192": r["192_global"], "global_max_135": r["135_global"]}, f"change {r['192_change_pct']:+.1f} %" if r["192_change_pct"] is not None else "")
    for combo, ths in mg["exacto_sn"].items():
        for key, r in ths.items():
            fail = r["max_class"] != "sigma->0" and r["clipped_max_at_grid_edge_0.05"]
            row("analisis_cmb/auditoria_v23/exacto_sn.json", f"linear {combo} {key}", "sigma 1e-6..5e-2 (ext. to 1)", f"sigma = {r['sigma_at_max_grid']:.3g}", "upper edge of the historical max grid 0.05 (failure)" if fail else r["max_class"],
                r["192_historical"], {"global_max_192": r["192_global"], "global_max_135": r["135_global"]}, f"change {r['192_change_pct']:+.1f} %" if r["192_change_pct"] is not None else "")
    fi = json.loads((HERE.parent / "fondo/fisher_diseno.json").read_text())["c_equaliser"]
    for k, v in fi.items():
        row("analisis_cmb/fondo/fisher_diseno.json", f"c_equaliser {k}: A_95 (40 periods)", "A in [0, 1] (4001)", f"A_at_max <= {max(v['A_at_max']):.4f}", "sigma->0 / interior near 0", f"A_95 max {max(v['A_95']):.3f}", None, "no A_95 at the edge 1")
    for f, q, g, m, c, h, n in (
        ("analisis_cmb/auditoria_v23/precision_sn.json", "marginal MC bounds (4 combos x 5 theta)", "sigma 0.01..0.7 (12)", "lnL max at sigma = 0 (SN) or 0.03-0.2 (others)", "interior / sigma->0; 'sin cota (<= 0.7)' entries are explicit (upper edge = no bound, not a max)", "see table", "none"),
        ("analisis_cmb/auditoria_v23/convergencia_sn.json", "marginal MC bounds vs N", "sigma 0.03..0.4 (0.7 for theta <= 0.3)", "lnL max at sigma = 0.03-0.4 (interior)", "interior", "see table", "none"),
        ("analisis_cmb/auditoria_v23/marginal_guiada.json", "exact guided marginal bounds", "sigma 0..0.4 (0.7 for theta = 0.3)", "lnZ max at sigma = 0.03-0.2 (interior)", "interior", "see table", "none"),
        ("analisis_cmb/auditoria_v23/exacto.json / exacto_sn.json", "5 % bounds (P crossings)", "sigma 0.05..0.7 / 0.01..0.7", "not a drop from a maximum", "n/a (P crossing; 'sin cota (<= 0.7)' explicit)", "see tables", "none"),
        ("results/tachyonic_rank1/SUMMARY.json", "scan dlnL(t_c)", "t_c 1e-4..1e4 (500)", "plateau dlnL = 0 for t_c >= 0.09 up to the upper edge", "physical edge (no growth = LCDM; amplitude -> 0); redone in VERIFICACION.md", "-12.7 'representative' (not a bound)", "none"),
        ("results/ou_qnm_mle/ou_qnm_mle.json", "MLE (sigma_X, theta, omega_R)", "ln sigma_X^2 in [-20, -1], theta >= 1e-3, omega_R [0, 10]", "sigma_X at the lower bound 4.54e-5", "sigma->0 (physical lower edge); theta, omega_R interior", "no bound (MLE)", "none"),
        ("results/joint_w0wa_sigma/joint_w0wa_sigma.json", "joint fits", "sigma_X free", "sigma_X = 0 at the best fit", "sigma->0", "no bound", "none"),
        ("analisis_cmb/resultados/takahashi/malla/malla.json", "bound = largest sigma with p95(dchi2) < 4", "sigma 1e-4..2.506e-2 (10)", "p95 < 4 on the whole grid", "no bound on the grid (explicit); not a drop from a maximum", "none", "none"),
        ("analisis_cmb/cuesta/gof.json", "GOF bound P(chi2 <= 10.27) = 5 %", "sigma 0.003..0.15 (16)", "not a drop from a maximum", "n/a (P crossing)", "see GOF.md", "none")):
        row(f, q, g, m, c, h, None, n)
    OUT["runtime_s"] = time.time() - t0
    (HERE / "barrido_maximos.json").write_text(json.dumps(OUT, indent=1, default=float), encoding="utf-8")
    print("rows:", len(OUT["rows"]), f"[{time.time() - t0:.0f} s]")


if __name__ == "__main__":
    main()
