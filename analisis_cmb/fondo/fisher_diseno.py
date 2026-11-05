#!/usr/bin/env python3
"""
Exploration only (FISHER_DISENO.md): where to put the next snapshot, by Fisher information on (sigma_X, theta).

Data model, as in perfil_13.py: rows (z, type) in alpha units, alpha = obs / pred(DESI DR2 best-fit flat LCDM: Omega_m,f =
0.29740, h r_d,f = 101.538 Mpc), covariance C (official 13x13 rescaled), residual delta alpha_i = S_i X(z_i) with S by type,
OU kernel K_ij = S_i S_j sigma_X^2 exp(-theta |x_i - x_j|), x = ln(1 + z). Background (Omega_m, h r_d) marginalised with a flat
prior through the linear templates T = [d alpha/d Omega_m, alpha] at the reference: the marginal likelihood has the projected
precision W = S^-1 - S^-1 T (T' S^-1 T)^-1 T' S^-1 (Sigma = C + K), and the Fisher matrix of the kernel parameters is
F_ab = 1/2 tr(W dSigma/da W dSigma/db), a, b in (sigma_X, theta), evaluated at sigma_X = 1.5e-2 and theta in {1, 0.3}.
"F(sigma_X)" and "F(theta)" below are the marginal informations 1/(F^-1)_aa (the other parameter marginalised); the diagonal
elements are also saved.
(a) New snapshot: one extra D_V/r_d with 1% error, uncorrelated, at z in 40 points of [0.05, 4]; relative gain F_new/F_13 - 1.
(b) Survey: DESI's total volume in N equal bins in ln(1 + z) between z = 0.1 and 2.3; each bin gives D_M/r_d and D_H/r_d with
    fractional errors eps_M sqrt(N/7), eps_H sqrt(N/7) and the mean DESI D_M-D_H correlation, eps_M, eps_H = mean fractional
    errors of DESI's six anisotropic bins (so N = 7 reproduces DESI's mean errors); N in {3, 5, 7, 10, 15, 20, 30, 50, 100}.
(c) Equaliser: harmonic residual X = A sin(2 pi ln a / P + phi), ln a = -x; delta alpha_i = S_i X_i; 95% bound on A versus P in
    40 log points of [0.05, 5], phase phi marginalised numerically (uniform, 180 points) and background marginalised (W with
    Sigma = C, no OU kernel); bound = upper crossing of Delta lnL = -1.92 from the maximum over A in [0, 1] (the convention of
    FONDO.md); "sin cota" if no crossing up to A = 1. With the real 13 DESI numbers, with the expected DESI data (alpha = 1),
    and with the expected data of the N = 20 survey of (b).

Writes analisis_cmb/fondo/fisher_diseno.{json,txt,png}; nothing in results/, figures/ or papers.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import perfil_13 as P13  # noqa: E402
import perfil_fondo as PF  # noqa: E402
from desi_dr2_data import load_gaussian_bao_vector  # noqa: E402

OM_F, HRD_F, C_KMS = PF.OM_F, PF.HRD_F, 299792.458
SIGMA_FID, THETAS = 1.5e-2, (1.0, 0.3)
Z_NEW = np.linspace(0.05, 4.0, 40)
NS = (3, 5, 7, 10, 15, 20, 30, 50, 100)
PERIODS = np.geomspace(0.05, 5.0, 40)
PHI = np.linspace(0, 2 * np.pi, 180, endpoint=False)
A_GRID = np.linspace(0, 1, 4001)
THR = -1.92


class Config:
    def __init__(self, z, kinds, cov, alpha=None):
        self.z, self.kinds, self.cov = np.asarray(z, float), list(kinds), np.asarray(cov, float)
        self.x = np.log1p(self.z)
        self.dx = np.abs(self.x[:, None] - self.x[None, :])
        self.S = np.array([P13.s_type(k, zz) for k, zz in zip(self.kinds, self.z)])
        e = 1e-4
        f0 = np.array([P13.dhat(k, OM_F, zz) for k, zz in zip(self.kinds, self.z)])
        df = np.array([(P13.dhat(k, OM_F + e, zz) - P13.dhat(k, OM_F - e, zz)) / (2 * e) for k, zz in zip(self.kinds, self.z)]) / f0
        self.T = np.column_stack([df, np.ones(len(self.z))])
        self.r0 = np.zeros(len(self.z)) if alpha is None else np.asarray(alpha, float) - 1.0

    @staticmethod
    def W(sigma_mat, T):
        Si = np.linalg.inv(sigma_mat)
        A = Si @ T
        return Si - A @ np.linalg.solve(T.T @ A, A.T)

    def fisher(self, sigma, theta):
        E = np.exp(-theta * self.dx)
        SS = np.outer(self.S, self.S)
        Sig = self.cov + SS * sigma**2 * E
        dS, dT = 2 * sigma * SS * E, -(sigma**2) * SS * E * self.dx
        W = self.W(Sig, self.T)
        a, b = W @ dS, W @ dT
        F = 0.5 * np.array([[np.trace(a @ a), np.trace(a @ b)], [np.trace(b @ a), np.trace(b @ b)]])
        Fi = np.linalg.inv(F)
        return {"F_ss": F[0, 0], "F_tt": F[1, 1], "F_st": F[0, 1], "Fm_s": 1 / Fi[0, 0], "Fm_t": 1 / Fi[1, 1],
                "corr": F[0, 1] / np.sqrt(F[0, 0] * F[1, 1])}

    def harmonic_bound(self, P):
        """95% bound on A for X = A sin(-2 pi x / P + phi), phi and background marginalised."""
        W0 = self.W(self.cov, self.T)
        U = self.S[None, :] * np.sin(-2 * np.pi * self.x[None, :] / P + PHI[:, None])  # n_phi x n
        # chi2(A, phi) = (r0 - A u)' W0 (r0 - A u) = c0 - 2 A c1(phi) + A^2 c2(phi)
        c0 = self.r0 @ W0 @ self.r0
        c1 = U @ (W0 @ self.r0)
        c2 = np.einsum("pi,ij,pj->p", U, W0, U)
        lnL = lambda A: -0.5 * c0 + np.log(np.maximum(np.mean(np.exp(A * c1 - 0.5 * A * A * c2)), 1e-300))
        ll = np.array([lnL(A) for A in A_GRID])
        j = int(np.argmax(ll))
        rel = ll - ll[j]
        k = np.nonzero((rel < THR) & (np.arange(len(A_GRID)) > j))[0]
        if len(k) == 0:
            return None, float(A_GRID[j])
        k = int(k[0])
        return float(brentq(lambda A: lnL(A) - ll[j] - THR, A_GRID[k - 1], A_GRID[k], xtol=1e-9)), float(A_GRID[j])


def desi_config() -> tuple[Config, dict]:
    g = load_gaussian_bao_vector()
    kinds = [lab[:2] for lab in g["labels"]]
    z = np.asarray(g["z_rows"], float)
    pred = np.array([C_KMS / (100 * HRD_F) * P13.dhat(k, OM_F, zz) for k, zz in zip(kinds, z)])
    cov = g["cov"] / np.outer(pred, pred)
    alpha = g["mean"] / pred
    err = np.sqrt(np.diag(cov))
    iM = [i for i, k in enumerate(kinds) if k == "DM"]
    iH = [i for i, k in enumerate(kinds) if k == "DH"]
    rho = []
    for i in iM:
        j = [jj for jj in iH if z[jj] == z[i]][0]
        rho.append(cov[i, j] / (err[i] * err[j]))
    info = {"eps_M": float(err[iM].mean()), "eps_H": float(err[iH].mean()), "rho_MH": float(np.mean(rho))}
    return Config(z, kinds, cov, alpha), info


def survey_config(N: int, info: dict) -> Config:
    edges = np.linspace(np.log1p(0.1), np.log1p(2.3), N + 1)
    zc = np.expm1(0.5 * (edges[:-1] + edges[1:]))
    z, kinds = np.repeat(zc, 2), ["DM", "DH"] * N
    eM, eH = info["eps_M"] * np.sqrt(N / 7), info["eps_H"] * np.sqrt(N / 7)
    cov = np.zeros((2 * N, 2 * N))
    for i in range(N):
        cov[2 * i, 2 * i], cov[2 * i + 1, 2 * i + 1] = eM**2, eH**2
        cov[2 * i, 2 * i + 1] = cov[2 * i + 1, 2 * i] = info["rho_MH"] * eM * eH
    return Config(z, kinds, cov)


def main() -> None:
    D, info = desi_config()
    base = {f"theta_{th:g}": D.fisher(SIGMA_FID, th) for th in THETAS}
    # (a) new snapshot
    a = {}
    for th in THETAS:
        gs, gt, cs, ct = [], [], [], []
        for zn in Z_NEW:
            Cn = np.zeros((14, 14))
            Cn[:13, :13], Cn[13, 13] = D.cov, 1e-4
            cfg = Config(np.append(D.z, zn), D.kinds + ["DV"], Cn)
            f = cfg.fisher(SIGMA_FID, th)
            gs.append(f["Fm_s"] / base[f"theta_{th:g}"]["Fm_s"] - 1)
            gt.append(f["Fm_t"] / base[f"theta_{th:g}"]["Fm_t"] - 1)
            cs.append(f["F_ss"] / base[f"theta_{th:g}"]["F_ss"] - 1)
            ct.append(f["F_tt"] / base[f"theta_{th:g}"]["F_tt"] - 1)
        gs, gt = np.array(gs), np.array(gt)
        a[f"theta_{th:g}"] = {"gain_sigma": gs.tolist(), "gain_theta": gt.tolist(),
                              "gain_sigma_conditional": cs, "gain_theta_conditional": ct,
                              "gain_at_z2_sigma": float(np.interp(2.0, Z_NEW, gs)), "gain_at_z2_theta": float(np.interp(2.0, Z_NEW, gt)),
                              "z_opt_sigma": float(Z_NEW[int(np.argmax(gs))]), "z_opt_theta": float(Z_NEW[int(np.argmax(gt))]),
                              "gain_max_sigma": float(gs.max()), "gain_max_theta": float(gt.max())}
    # (b) survey
    b = {}
    for th in THETAS + (3.0, 10.0):
        fs, ft, cs, ct = [], [], [], []
        for N in NS:
            f = survey_config(N, info).fisher(SIGMA_FID, th)
            fs.append(f["Fm_s"]), ft.append(f["Fm_t"]), cs.append(f["F_ss"]), ct.append(f["F_tt"])
        fs, ft, cs, ct = np.array(fs), np.array(ft), np.array(cs), np.array(ct)
        n90 = int(NS[int(np.nonzero(ft >= 0.9 * ft.max())[0][0])])
        if th not in THETAS:
            base[f"theta_{th:g}"] = D.fisher(SIGMA_FID, th)
        b[f"theta_{th:g}"] = {"F_sigma": fs.tolist(), "F_theta": ft.tolist(), "N_90_theta": n90, "extra": th not in THETAS,
                              "F_sigma_cond_rel_N7": (cs / cs[NS.index(7)]).tolist(), "F_theta_cond_rel_N7": (ct / ct[NS.index(7)]).tolist(),
                              "F_sigma_rel_N7": (fs / fs[NS.index(7)]).tolist(), "F_theta_rel_N7": (ft / ft[NS.index(7)]).tolist(),
                              "desi13_F_sigma": base[f"theta_{th:g}"]["Fm_s"], "desi13_F_theta": base[f"theta_{th:g}"]["Fm_t"]}
    # (c) equaliser
    D_exp = Config(D.z, D.kinds, D.cov)
    S20 = survey_config(20, info)
    c = {}
    for name, cfg in (("desi_real", D), ("desi_expected", D_exp), ("survey_N20", S20)):
        bounds, amax = zip(*(cfg.harmonic_bound(P) for P in PERIODS))
        bl = [bb if bb is not None else None for bb in bounds]
        vals = np.array([bb if bb is not None else np.nan for bb in bounds])
        c[name] = {"A_95": bl, "A_at_max": list(amax), "P_best": float(PERIODS[int(np.nanargmin(vals))]),
                   "A_95_best": float(np.nanmin(vals)), "n_sin_cota": int(np.sum(np.isnan(vals)))}
    out = {"reference": {"Omega_m": OM_F, "h_rd_Mpc": HRD_F}, "sigma_fid": SIGMA_FID, "thetas": THETAS,
           "desi_errors": info, "fisher_desi13": base, "z_new": Z_NEW.tolist(), "a_new_snapshot": a, "N": NS, "b_survey": b,
           "periods": PERIODS.tolist(), "c_equaliser": c}
    (HERE / "fisher_diseno.json").write_text(json.dumps(out, indent=1), encoding="utf-8")

    # figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
    for th, ls in zip(THETAS, ("-", "--")):
        k = f"theta_{th:g}"
        ax[0].plot(Z_NEW, 100 * np.array(a[k]["gain_sigma"]), ls, color="C0", label=f"F(σ_X), θ = {th:g}")
        ax[0].plot(Z_NEW, 100 * np.array(a[k]["gain_theta"]), ls, color="C3", label=f"F(θ), θ = {th:g}")
        ax[1].plot(NS, b[k]["F_sigma_rel_N7"], ls, color="C0", marker="o", ms=3, label=f"F(σ_X), θ = {th:g}")
        ax[1].plot(NS, b[k]["F_theta_rel_N7"], ls, color="C3", marker="o", ms=3, label=f"F(θ), θ = {th:g}")
    ax[0].set(xlabel="z de la foto nueva (D_V/r_d, error 1 %)", ylabel="ganancia de información [%]", title="(a) foto nueva")
    ax[1].set(xscale="log", xlabel="N tramos (volumen de DESI)", ylabel="F / F(N = 7)", title="(b) cartografiado")
    ax[1].axhline(1, color="k", lw=0.5)
    for name, lab, st in (("desi_real", "DESI, 13 números reales", "-"), ("desi_expected", "DESI, datos esperados", ":"),
                          ("survey_N20", "cartografiado N = 20, esperado", "--")):
        v = np.array([x if x is not None else np.nan for x in c[name]["A_95"]])
        ax[2].plot(PERIODS, v, st, label=lab)
    ax[2].set(xscale="log", yscale="log", xlabel="P (periodo en ln a)", ylabel="A al 95 %", title="(c) ecualizador")
    for aa in ax:
        aa.legend(fontsize=8), aa.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(HERE / "fisher_diseno.png", dpi=130)

    lines = [f"Fisher en sigma_X = {SIGMA_FID}, fondo marginalizado; F = informacion marginal 1/(F^-1)_aa",
             f"errores medios DESI: eps_M = {info['eps_M']:.4f}, eps_H = {info['eps_H']:.4f}, rho_MH = {info['rho_MH']:.3f}"]
    for th in THETAS:
        k = f"theta_{th:g}"
        f = base[k]
        lines += [f"theta = {th:g}: DESI 13: F(sigma) = {f['Fm_s']:.4g} (1/sqrt = {1/np.sqrt(f['Fm_s']):.4f}), F(theta) = {f['Fm_t']:.4g} "
                  f"(1/sqrt = {1/np.sqrt(f['Fm_t']):.3f}), corr = {f['corr']:+.3f}",
                  f"  (a) z optimo: sigma {a[k]['z_opt_sigma']:.2f} (+{100*a[k]['gain_max_sigma']:.1f} %), theta {a[k]['z_opt_theta']:.2f} "
                  f"(+{100*a[k]['gain_max_theta']:.1f} %)",
                  f"      ganancia en z = 2: sigma +{100*a[k]['gain_at_z2_sigma']:.0f} %, theta +{100*a[k]['gain_at_z2_theta']:.0f} %; "
                  f"condicional (diagonal) max: sigma +{100*max(a[k]['gain_sigma_conditional']):.0f} %, theta +{100*max(a[k]['gain_theta_conditional']):.0f} %",
                  "  (b) N: " + " ".join(f"{N:>6d}" for N in NS),
                  "      F(sigma)/F7: " + " ".join(f"{v:6.3f}" for v in b[k]["F_sigma_rel_N7"]),
                  "      F(theta)/F7: " + " ".join(f"{v:6.3f}" for v in b[k]["F_theta_rel_N7"]) + f"   N_90 = {b[k]['N_90_theta']}",
                  "      cond. F_tt/F7: " + " ".join(f"{v:6.3f}" for v in b[k]["F_theta_cond_rel_N7"])]
    for th in (3.0, 10.0):
        k = f"theta_{th:g}"
        lines += [f"  (b) extra theta = {th:g}: F(theta)/F7: " + " ".join(f"{v:6.3f}" for v in b[k]["F_theta_rel_N7"]) + f"   N_90 = {b[k]['N_90_theta']}",
                  f"      F(sigma)/F7: " + " ".join(f"{v:6.3f}" for v in b[k]["F_sigma_rel_N7"])]
    for name in ("desi_real", "desi_expected", "survey_N20"):
        lines.append(f"  (c) {name}: mejor P = {c[name]['P_best']:.3f}, A_95 = {c[name]['A_95_best']:.4f}, sin cota en {c[name]['n_sin_cota']} P")
        lines.append("      P:    " + " ".join(f"{P:7.3f}" for P in PERIODS))
        lines.append("      A_95: " + " ".join(f"{v:7.4f}" if v is not None else "   ----" for v in c[name]["A_95"]))
    (HERE / "fisher_diseno.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
