#!/usr/bin/env python3
"""
Exploration only (MODOS.md): do the few modes of X that the data measure well explain how much each family is bounded in
FAMILIAS.md? Nothing in results/ or papers changes.

Phase 1 (`--fase 1`; does not open familias.json nor FAMILIAS.md):
a. Fisher matrix on X with the canonical linear response A of marginal_guiada.py (BAO + theta_* + omega_m + SN, convention B,
   exacto grid, 1595 rows): F = A_b' P A_b with P = C^-1 - C^-1 B (B' C^-1 B)^-1 B' C^-1 (h, omega_m, M marginalised), and
   X piecewise constant on nb uniform bins in u = ln(1+z) up to z = 3 plus one bin up to z_* (A_b = A M, M the node->bin
   indicator). nb = 20 and 40.
b. Normalisation with D = diag(bin widths in u): F_norm = D^-1/2 F D^-1/2 and, per family, K_norm = D^1/2 K_b D^1/2 with K_b the
   covariance of the bin averages of X (node means of K on the grid). Always the pair F_norm / K_norm. Check tr(F K_b) = tr(F_norm K_norm).
c. Eigenvalues lambda_i and modes e_i of F_norm; sigma(a_i) = lambda_i^-1/2. The constant direction (D^1/2 1, exactly null in
   convention B) is removed; any other null direction (lambda < 1e-10 lambda_max) is reported.
d. The same with P = C^-1 (h, omega_m, M fixed), for prediction 8.
e. Per family (OU, Matern 3/2, Gaussian, triangular kernel with Gaussian prior) and l in FAMILIAS, sigma_X = 1:
   v_i = e_i' K_norm e_i, S_k = sum_{i<=k} lambda_i v_i, k = 1, 2, 3, all.
f. Fixed formula sigma_pred = c / sqrt(S_k). Two constants calibrated on (OU, l = 1) with k = 2 and frozen: c_5 against 0.169
   (5 % bound) and c_lin against 0.272 (linear -1.92 bound with the canonical A, MARGINAL_GUIADA.md). Same c for every k, cell
   and family (the per-k recalibrated constants are reported separately, as information only).
g. 20 vs 40 bins: ratios of the first lambda_i; if two lambda differ by less than 20 % the subspace they span is compared through
   principal angles (modes mapped to the grid, L2(du) inner product).
h. modos_pred.json is written; its sha256 goes to MODOS.md.
Phase 2 (`--fase 2`, after the phase-1 commit): opens familias.json, compares sigma_pred with the 5 % and linear -1.92 bounds
cell by cell (ratio, median |ln ratio| per k), including the real Steps, draws the first three modes against z
(modos_modos.png; the sign of each mode is arbitrary) and writes modos.json.

    ../venv/bin/python modos.py --fase 1
    ../venv/bin/python modos.py --fase 2
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import cho_solve, subspace_angles

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import marginal_guiada as MG  # noqa: E402

U = MG.G.U
N = len(U)
USTAR = float(U[-1])
U3 = float(np.log(4.0))  # z = 3
ELLS = (0.1, 0.3, 1.0, 3.0, 10.0)
FAMS = ("OU", "Matern32", "Gauss", "Triangular")
NBS = (20, 40)
KS = (1, 2, 3)
CAL = {"5pct": 0.169, "lin192": 0.272}
K_CAL = 2
NB_CAL = 20
NULL_TOL = 1e-10
# trapezoid weights on the grid (L2(du) inner product for grid functions)
WGRID = np.zeros(N)
WGRID[1:] += 0.5 * np.diff(U)
WGRID[:-1] += 0.5 * np.diff(U)


def kfun(fam: str, r, l: float):
    r = np.abs(np.asarray(r, dtype=float))
    if fam == "OU":
        return np.exp(-r / l)
    if fam == "Matern32":
        a = np.sqrt(3.0) * r / l
        return (1.0 + a) * np.exp(-a)
    if fam == "Gauss":
        return np.exp(-(r * r) / (2.0 * l * l))
    if fam == "Triangular":
        return np.maximum(1.0 - r / l, 0.0)
    raise ValueError(fam)


class Bins:
    def __init__(self, nb: int):
        self.nb = nb
        self.edges = np.concatenate([np.linspace(0.0, U3, nb + 1), [USTAR]])
        self.n = nb + 1
        idx = np.searchsorted(self.edges, U, side="right") - 1
        self.idx = np.clip(idx, 0, self.n - 1)
        self.widths = np.diff(self.edges)
        self.D = np.diag(self.widths)
        self.Dm12 = np.diag(self.widths**-0.5)
        self.Dp12 = np.diag(self.widths**0.5)
        self.M = np.zeros((N, self.n))
        self.M[np.arange(N), self.idx] = 1.0
        self.counts = self.M.sum(axis=0)
        assert self.counts.min() > 0
        self.Mavg = self.M / self.counts[None, :]
        self.u_c = 0.5 * (self.edges[:-1] + self.edges[1:])
        self.z_lo, self.z_hi, self.z_c = np.expm1(self.edges[:-1]), np.expm1(self.edges[1:]), np.expm1(self.u_c)


def fisher(A: np.ndarray, Bm: np.ndarray, b: Bins, marginalise: bool) -> np.ndarray:
    Ab = A @ b.M  # (NROWS, n)
    CiA = cho_solve(MG.CF, Ab)
    if marginalise:
        CiB = cho_solve(MG.CF, Bm)
        G = Bm.T @ CiB
        CiA = CiA - CiB @ np.linalg.solve(G, Bm.T @ CiA)
    F = Ab.T @ CiA
    return 0.5 * (F + F.T)


def modes(F: np.ndarray, b: Bins) -> dict:
    Fn = b.Dm12 @ F @ b.Dm12
    w, V = np.linalg.eigh(Fn)
    order = np.argsort(w)[::-1]
    w, V = w[order], V[:, order]
    lam_max = float(w[0])
    null = np.abs(w) < NULL_TOL * lam_max
    n_const = np.sqrt(b.widths)
    n_const /= np.linalg.norm(n_const)
    Vn = V[:, null]
    const_in_null = float(np.linalg.norm(Vn.T @ n_const)) if null.any() else 0.0
    keep = ~null
    return {"Fn": Fn, "lam": w[keep], "E": V[:, keep], "n_null": int(null.sum()), "null_lams": w[null].tolist(), "const_overlap_with_null_span": const_in_null,
            "Fn_const_residual": float(np.linalg.norm(Fn @ n_const) / lam_max), "lam_max": lam_max}


def kernel_bins(fam: str, l: float, b: Bins) -> np.ndarray:
    K = kfun(fam, U[:, None] - U[None, :], l)
    return b.Mavg.T @ K @ b.Mavg


def mode_on_grid(e: np.ndarray, b: Bins) -> np.ndarray:
    """Normalised-coordinate mode e -> amplitude a(u) = D^-1/2 e on the grid (piecewise constant)."""
    a = e / np.sqrt(b.widths)
    return a[b.idx]


def describe_mode(e: np.ndarray, b: Bins) -> dict:
    a = e / np.sqrt(b.widths)  # X-space amplitude per bin
    k = int(np.argmax(np.abs(e)))
    w2 = e**2 / np.sum(e**2)
    half = 0.5 * np.abs(a).max()
    above = np.abs(a) >= half
    # width: z range of the bins at >= half the maximum amplitude (connected run containing the maximum)
    ka = int(np.argmax(np.abs(a)))
    lo = ka
    while lo - 1 >= 0 and above[lo - 1]:
        lo -= 1
    hi = ka
    while hi + 1 < b.n and above[hi + 1]:
        hi += 1
    sign_changes = int(np.sum(np.diff(np.sign(a[a != 0])) != 0))
    return {"bin_max_abs_e": k, "z_range_bin_max_abs_e": [float(b.z_lo[k]), float(b.z_hi[k])], "bin_max_abs_a": ka,
            "z_range_bin_max_abs_a": [float(b.z_lo[ka]), float(b.z_hi[ka])], "z_half_max_range_a": [float(b.z_lo[lo]), float(b.z_hi[hi])],
            "sign_at_max": int(np.sign(a[ka])), "sign_changes": sign_changes, "frac_weight_e2_above_z0.5": float(np.sum(w2[b.z_c > 0.5])),
            "frac_weight_e2_last_bin": float(w2[-1]), "e": e.tolist(), "a_per_bin": a.tolist()}


def phase1() -> dict:
    A, Bm = MG.build_A(), MG.build_B()
    out = {"ells": ELLS, "families": FAMS, "nbs": NBS, "u_z3": U3, "u_zstar": USTAR, "null_tol": NULL_TOL, "calibration": {"cell": "OU|l_1", "k": K_CAL, "nb": NB_CAL, "targets": CAL},
           "A_rows_sum_max_abs": float(np.max(np.abs(A.sum(axis=1)))), "bins": {}}
    R = {}
    for nb in NBS:
        b = Bins(nb)
        F_m = fisher(A, Bm, b, True)
        F_f = fisher(A, Bm, b, False)
        mm, mf = modes(F_m, b), modes(F_f, b)
        rec = {"n_bins": b.n, "edges_u": b.edges.tolist(), "z_edges": np.expm1(b.edges).tolist(), "widths_u": b.widths.tolist(),
               "marginalised": {"lambda": mm["lam"].tolist(), "sigma_a": (mm["lam"] ** -0.5).tolist(), "n_null": mm["n_null"], "null_lambdas": mm["null_lams"],
                                "const_overlap_with_null_span": mm["const_overlap_with_null_span"], "Fn_const_residual_rel": mm["Fn_const_residual"],
                                "n_modes_sigma_lt_0.1": int(np.sum(mm["lam"] ** -0.5 < 0.1)), "modes_first3": [describe_mode(mm["E"][:, i], b) for i in range(3)]},
               "fixed": {"lambda": mf["lam"].tolist(), "sigma_a": (mf["lam"] ** -0.5).tolist(), "n_null": mf["n_null"], "null_lambdas": mf["null_lams"],
                         "const_overlap_with_null_span": mf["const_overlap_with_null_span"], "n_modes_sigma_lt_0.1": int(np.sum(mf["lam"] ** -0.5 < 0.1)),
                         "modes_first3": [describe_mode(mf["E"][:, i], b) for i in range(3)]}}
        # prediction 8: information of the first mode with and without marginalisation
        e1f = mf["E"][:, 0]
        e1m = mm["E"][:, 0]
        rec["marginalisation_effect"] = {"lambda1_marg": float(mm["lam"][0]), "lambda1_fixed": float(mf["lam"][0]), "reduction_lambda1": float(1 - mm["lam"][0] / mf["lam"][0]),
                                         "info_along_e1_fixed_under_marg": float(e1f @ mm["Fn"] @ e1f), "reduction_along_e1_fixed": float(1 - (e1f @ mm["Fn"] @ e1f) / mf["lam"][0]),
                                         "overlap_e1_marg_e1_fixed": float(abs(e1m @ e1f)), "trace_marg_over_fixed": float(np.trace(mm["Fn"]) / np.trace(mf["Fn"]))}
        # families: S_k
        fam_rec = {}
        for fam in FAMS:
            for l in ELLS:
                Kb = kernel_bins(fam, l, b)
                Kn = b.Dp12 @ Kb @ b.Dp12
                tr1, tr2 = float(np.trace(F_m @ Kb)), float(np.trace(mm["Fn"] @ Kn))
                v = np.einsum("ij,jk,ki->i", mm["E"].T, Kn, mm["E"])
                lv = mm["lam"] * v
                S = {f"k_{k}": float(np.sum(lv[:k])) for k in KS}
                S["k_all"] = float(np.sum(lv))
                fam_rec[f"{fam}|l_{l:g}"] = {"v": v.tolist(), "lambda_v": lv.tolist(), "S": S, "trace_F_Kb": tr1, "trace_Fn_Kn": tr2, "trace_rel_diff": abs(tr1 - tr2) / abs(tr1),
                                             "S_all_over_trace": S["k_all"] / tr1, "frac_S_first3": float(np.sum(lv[:3]) / np.sum(lv))}
        rec["families"] = fam_rec
        R[nb] = (b, mm, mf, rec)
        out["bins"][f"nb_{nb}"] = rec
    # constants, frozen (OU, l = 1, k = 2, nb = 20)
    S_cal = out["bins"][f"nb_{NB_CAL}"]["families"]["OU|l_1"]["S"][f"k_{K_CAL}"]
    c = {name: tgt * np.sqrt(S_cal) for name, tgt in CAL.items()}
    out["calibration"]["S_k2_OU_l1"] = S_cal
    out["calibration"]["c"] = c
    out["calibration"]["c_per_k_info_only"] = {name: {kk: tgt * np.sqrt(out["bins"][f"nb_{NB_CAL}"]["families"]["OU|l_1"]["S"][kk]) for kk in list(f"k_{k}" for k in KS) + ["k_all"]} for name, tgt in CAL.items()}
    for nb in NBS:
        for key, fr in out["bins"][f"nb_{nb}"]["families"].items():
            fr["sigma_pred"] = {name: {kk: float(c[name] / np.sqrt(S)) for kk, S in fr["S"].items()} for name in CAL}
    # 20 vs 40 bins
    b20, m20, _, _ = R[20]
    b40, m40, _, _ = R[40]
    lam20, lam40 = m20["lam"], m40["lam"]
    cmp = {"lambda_ratio_40_over_20_first6": (lam40[:6] / lam20[:6]).tolist(), "max_abs_rel_diff_first3": float(np.max(np.abs(lam40[:3] / lam20[:3] - 1)))}
    # mode overlaps in L2(du): map to the grid
    G20 = np.array([mode_on_grid(m20["E"][:, i], b20) for i in range(6)]).T
    G40 = np.array([mode_on_grid(m40["E"][:, i], b40) for i in range(6)]).T
    Wh = np.sqrt(WGRID)[:, None]
    Q20, Q40 = G20 * Wh, G40 * Wh
    Q20 /= np.linalg.norm(Q20, axis=0)
    Q40 /= np.linalg.norm(Q40, axis=0)
    cmp["mode_overlap_one_to_one_first6"] = [float(abs(Q20[:, i] @ Q40[:, i])) for i in range(6)]
    # near-degenerate groups (consecutive lambdas within 20 %) -> principal angles of the spanned subspaces
    groups, g = [], [0]
    for i in range(1, 6):
        if abs(lam20[i] / lam20[i - 1] - 1) < 0.2 or abs(lam40[i] / lam40[i - 1] - 1) < 0.2:
            g.append(i)
        else:
            groups.append(g)
            g = [i]
    groups.append(g)
    cmp["groups_first6"] = groups
    cmp["principal_angles_deg_by_group"] = [np.degrees(subspace_angles(Q20[:, gg], Q40[:, gg])).tolist() for gg in groups]
    out["compare_20_40"] = cmp
    out["stop_conditions"] = {"first3_lambda_differ_gt_20pct": bool(cmp["max_abs_rel_diff_first3"] > 0.2),
                              "more_than_one_null_direction": {f"nb_{nb}": bool(out["bins"][f"nb_{nb}"]["marginalised"]["n_null"] > 1) for nb in NBS}}
    return out


def _default(o):
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(str(type(o)))


def summary1(out: dict) -> None:
    for nb in NBS:
        r = out["bins"][f"nb_{nb}"]
        m = r["marginalised"]
        print(f"nb={nb}: n_null(marg)={m['n_null']} null lams={m['null_lambdas']} const overlap={m['const_overlap_with_null_span']:.6f}; "
              f"lambda first6 = {[round(x, 3) for x in m['lambda'][:6]]}; sigma_a first6 = {[round(x, 4) for x in m['sigma_a'][:6]]}; n(sigma<0.1)={m['n_modes_sigma_lt_0.1']}")
        f = r["fixed"]
        print(f"        fixed: n_null={f['n_null']} lambda first6 = {[round(x, 3) for x in f['lambda'][:6]]}; n(sigma<0.1)={f['n_modes_sigma_lt_0.1']}; marg effect {json.dumps({k: round(v, 3) for k, v in r['marginalisation_effect'].items()})}")
        for i, md in enumerate(m["modes_first3"]):
            print(f"        mode {i + 1}: max|e| in z {md['z_range_bin_max_abs_e']}, max|a| in z {md['z_range_bin_max_abs_a']}, half-max z {md['z_half_max_range_a']}, sign changes {md['sign_changes']}, weight above z=0.5: {md['frac_weight_e2_above_z0.5']:.3f}, last bin {md['frac_weight_e2_last_bin']:.3f}")
        print("        tr check max rel diff:", max(v["trace_rel_diff"] for v in r["families"].values()))
    print("calibration:", json.dumps(out["calibration"], default=_default))
    print("20 vs 40:", json.dumps(out["compare_20_40"], default=_default))
    print("stop:", out["stop_conditions"])
    r = out["bins"]["nb_20"]["families"]
    print("sigma_pred (nb=20, c frozen): family|l : 5pct k1/k2/k3/all ; lin k1/k2/k3/all ; frac S first3")
    for key, fr in r.items():
        sp = fr["sigma_pred"]
        print(f"   {key:18s} {[round(sp['5pct'][k], 3) for k in ('k_1', 'k_2', 'k_3', 'k_all')]}  {[round(sp['lin192'][k], 3) for k in ('k_1', 'k_2', 'k_3', 'k_all')]}  {fr['frac_S_first3']:.3f}")


# ----------------------------------------------------------------------------------------------------------------------
def phase2() -> dict:
    pred = json.loads((HERE / "modos_pred.json").read_text())
    fam = json.loads((HERE / "familias.json").read_text())
    out = {"modos_pred_sha256": hashlib.sha256((HERE / "modos_pred.json").read_bytes()).hexdigest(), "cells": {}, "summary": {}}
    lin_ou1 = fam["linear"]["OU|l_1"]["192"]
    out["OU_l1_linear_in_familias"] = {"value": lin_ou1, "calibration_target": CAL["lin192"], "equal": abs(lin_ou1 - CAL["lin192"]) < 5e-4}
    out["OU_l1_5pct_in_familias"] = {"value": fam["bounds"]["OU|l_1"]["exact_5pct"], "calibration_target": CAL["5pct"]}
    kmap = {"OU": "OU", "Matern32": "Matern32", "Gauss": "Gauss", "Steps": "Triangular", "StepsGauss": "Triangular"}
    for nb in NBS:
        fr = pred["bins"][f"nb_{nb}"]["families"]
        cells = {}
        for ffam, kf in kmap.items():
            for l in ELLS:
                key = f"{ffam}|l_{l:g}"
                if key not in fam["bounds"]:
                    continue
                b5 = fam["bounds"][key]["exact_5pct"]
                lin = fam["linear"].get(key if ffam != "StepsGauss" else f"Steps|l_{l:g}", {}).get("192")
                sp = fr[f"{kf}|l_{l:g}"]["sigma_pred"]
                rec = {"kernel": kf, "bound_5pct": b5, "bound_lin192": lin, "sigma_pred_5pct": sp["5pct"], "sigma_pred_lin": sp["lin192"], "is_calibration": key == "OU|l_1"}
                for kk in ("k_1", "k_2", "k_3", "k_all"):
                    rec[f"ratio_5pct_{kk}"] = (sp["5pct"][kk] / b5) if isinstance(b5, float) and b5 > 0 else None
                    rec[f"ratio_lin_{kk}"] = (sp["lin192"][kk] / lin) if isinstance(lin, float) and lin > 0 else None
                cells[key] = rec
        out["cells"][f"nb_{nb}"] = cells
        S = {}
        for rule in ("5pct", "lin"):
            for kk in ("k_1", "k_2", "k_3", "k_all"):
                for group, fams in (("OU_Matern_Gauss", ("OU", "Matern32", "Gauss")), ("StepsGauss", ("StepsGauss",)), ("Steps", ("Steps",))):
                    rat = [c[f"ratio_{rule}_{kk}"] for k2, c in cells.items() if k2.split("|")[0] in fams and not c["is_calibration"] and c[f"ratio_{rule}_{kk}"] is not None]
                    if not rat:
                        continue
                    la = np.abs(np.log(rat))
                    S[f"{rule}|{kk}|{group}"] = {"n": len(rat), "ratios": rat, "median_abs_ln_ratio": float(np.median(la)), "max_abs_ln_ratio": float(la.max()),
                                                 "frac_within_1.3": float(np.mean(la <= np.log(1.3))), "frac_within_1.2": float(np.mean(la <= np.log(1.2))), "frac_within_1.5": float(np.mean(la <= np.log(1.5)))}
        # unbounded cells: sigma_pred > 1?
        unb = {}
        for key, c in cells.items():
            if c["bound_5pct"] is None:
                unb[key + "|5pct"] = {kk: c["sigma_pred_5pct"][kk] for kk in ("k_1", "k_2", "k_3", "k_all")}
            if c["bound_lin192"] is None:
                unb[key + "|lin"] = {kk: c["sigma_pred_lin"][kk] for kk in ("k_1", "k_2", "k_3", "k_all")}
        out["summary"][f"nb_{nb}"] = {"stats": S, "unbounded_cells_sigma_pred": unb}
    # figure: first three modes (marginalised, nb = 20 solid, nb = 40 dashed), X-space amplitude a(u) = D^-1/2 e against z
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(3, 1, figsize=(7.5, 8.5), sharex=True)
    for i, ax in enumerate(axes):
        for nb, ls in ((20, "-"), (40, "--")):
            r = pred["bins"][f"nb_{nb}"]
            md = r["marginalised"]["modes_first3"][i]
            z_edges = np.array(r["z_edges"])
            a = np.array(md["a_per_bin"])
            # sign convention for the plot only: positive at the bin of max |a| of the 20-bin mode
            if nb == 20:
                s = np.sign(a[int(np.argmax(np.abs(a)))]) or 1.0
                a20_grid = a
            else:
                # align the sign with the 20-bin mode through the L2 overlap
                b20, b40 = Bins(20), Bins(40)
                s = np.sign(np.sum(mode_on_grid(np.array(a20_grid) * np.sqrt(b20.widths), b20) * mode_on_grid(a * np.sqrt(b40.widths), b40) * WGRID)) or 1.0
            zz = np.repeat(z_edges, 2)[1:-1]
            aa = np.repeat(s * a, 2)
            ax.plot(zz, aa, ls, lw=1.6 if nb == 20 else 1.0, label=f"{nb} + 1 tramos" + ("" if nb == 20 else " (último tramo hasta z*)"))
        lam = pred["bins"]["nb_20"]["marginalised"]["lambda"][i]
        ax.axhline(0, color="0.6", lw=0.6)
        ax.set_xscale("log")
        ax.set_xlim(0.03, 1200)
        ax.set_ylabel(f"modo {i + 1}: a(z)")
        ax.text(0.02, 0.92, f"λ_{i + 1} = {lam:.3g}, σ(a_{i + 1}) = {lam ** -0.5:.3g}", transform=ax.transAxes, fontsize=9, va="top")
        ax.axvline(0.5, color="0.75", lw=0.6, ls=":")
    axes[0].legend(fontsize=8, loc="upper right")
    axes[-1].set_xlabel("z (escala logarítmica; el último tramo va de z = 3 a z* = 1090)")
    fig.suptitle("Modos de Fisher de X(u) (BAO + θ* + ω_m + SN, convenio B, h, ω_m, M marginalizados)\nEl signo de cada modo es arbitrario; a(u) = D^{-1/2} e por tramo", fontsize=10)
    fig.tight_layout()
    fig.savefig(HERE / "modos_modos.png", dpi=130)
    out["figure"] = "modos_modos.png"
    return out


def summary2(out: dict) -> None:
    print("OU l=1 in familias.json:", out["OU_l1_linear_in_familias"], out["OU_l1_5pct_in_familias"])
    for nb in NBS:
        print(f"nb={nb}")
        for key, c in out["cells"][f"nb_{nb}"].items():
            f = lambda v: "   -  " if v is None else f"{v:6.3f}"
            print(f"   {key:18s} 5%: {f(c['bound_5pct'])} pred k1/k2/k3/all {[round(c['sigma_pred_5pct'][k], 3) for k in ('k_1', 'k_2', 'k_3', 'k_all')]} ratio k2 {f(c['ratio_5pct_k_2'])} | lin: {f(c['bound_lin192'])} pred {[round(c['sigma_pred_lin'][k], 3) for k in ('k_1', 'k_2', 'k_3', 'k_all')]} ratio k2 {f(c['ratio_lin_k_2'])}{'  [cal]' if c['is_calibration'] else ''}")
        for key, s in out["summary"][f"nb_{nb}"]["stats"].items():
            print(f"   {key:28s} n={s['n']:2d} median|ln r|={s['median_abs_ln_ratio']:.3f} max={s['max_abs_ln_ratio']:.3f} within 1.2/1.3/1.5: {s['frac_within_1.2']:.2f}/{s['frac_within_1.3']:.2f}/{s['frac_within_1.5']:.2f}")
        print("   unbounded:", json.dumps({k: {kk: round(v, 3) for kk, v in d.items()} for k, d in out["summary"][f"nb_{nb}"]["unbounded_cells_sigma_pred"].items()}))


if __name__ == "__main__":
    fase = sys.argv[sys.argv.index("--fase") + 1] if "--fase" in sys.argv else "1"
    if fase == "1":
        out = phase1()
        (HERE / "modos_pred.json").write_text(json.dumps(out, indent=1, default=_default), encoding="utf-8")
        summary1(out)
        print("sha256 modos_pred.json:", hashlib.sha256((HERE / "modos_pred.json").read_bytes()).hexdigest())
    else:
        out = phase2()
        (HERE / "modos.json").write_text(json.dumps(out, indent=1, default=_default), encoding="utf-8")
        summary2(out)
