#!/usr/bin/env python3
"""
Robustness of the simulation test (bao_sim.py) against "fixed parameters": omega_m free per realisation.

For a fixed realisation the dark-energy history Omega_L(z) is invariant under omega_m -> c omega_m
(Das I Sec. 4.4, Model 1 with S_0 = 0; by construction in Model 2), H0 scales as sqrt(omega_m), so the
dimensionless distances D_M H0/c and H0/H(z) are invariant. With r_d(omega_m) from Aubourg et al. 2015
(omega_b = 0.02237, omega_nu = 0.000645 fixed; normalised to CAMB at Planck, see rd_planck.json) every
D/r_d scales by the same factor
        s(omega_m) = sqrt(omega_m0 / omega_m) * r_d(omega_m0) / r_d(omega_m),
so chi^2(s) = s^2 a - 2 s b + c is minimised exactly over omega_m in [0.10, 0.20] by clipping
s* = b/a to [s(0.20), s(0.10)]. Radiation breaks the invariance at the 1e-3 level in H^2 at z < 2.4;
the size of that error is checked for Model 2 by recomputing chi^2 exactly (check_exact_m2).

References:
  A  : chi^2 <= 29.84 (LCDM-Planck, oscilador.json)               -- the generous criterion
  B  : chi^2 - 29.84 < 4
  A' : chi^2 <= chi^2_min of LCDM fitted to DESI DR2 alone        -- the demanding criterion
       (Omega_m and H0 r_d free; this script)

Writes resultados/robustez.json.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize_scalar

import bao_sim as B

HERE = Path(__file__).resolve().parent
OUT = HERE / "resultados"
C = B.C
OMB = 0.02237
RDJ = json.loads((OUT / "rd_planck.json").read_text())
OMNU = RDJ["rd_aubourg"]["omega_nu"]
OM0 = B.OM
RD0 = B.R_D
WM_RANGE = (0.10, 0.20)


def rd_aubourg(omega_b, omega_cb, omega_nu):
    """Aubourg et al. 2015 (arXiv:1411.1074), Sec. II, eq. labelled eqn:rd (same as rd_planck.py)."""
    return 55.154 * np.exp(-72.3 * (omega_nu + 0.0006) ** 2) / (np.asarray(omega_cb) ** 0.25351 * omega_b**0.12807)


def scale(wm):
    wm = np.asarray(wm, dtype=float)
    return np.sqrt(OM0 / wm) * rd_aubourg(OMB, OM0 - OMNU, OMNU) / rd_aubourg(OMB, wm - OMNU, OMNU)


WM_GRID = np.linspace(*WM_RANGE, 20001)
S_GRID = scale(WM_GRID)  # decreasing in omega_m
S_LO, S_HI = S_GRID[-1], S_GRID[0]


def theory13(dmh7: np.ndarray, eh7: np.ndarray, H0: np.ndarray, zu: np.ndarray) -> np.ndarray:
    """(n, 7) dimensionless distances -> (13, n) D/r_d at omega_m0 and r_d = RD0."""
    iz = np.searchsorted(zu, B.ZR)
    DM = dmh7[:, iz].T * C / H0[None, :]
    DH = eh7[:, iz].T * C / H0[None, :]
    DV = (B.ZR[:, None] * DM**2 * DH) ** (1 / 3)
    K = np.array(B.KIND)[:, None]
    return np.where(K == "DM", DM, np.where(K == "DH", DH, DV)) / RD0


def chi2_free_wm(t: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """t: (13, n). Returns chi2 at omega_m0, chi2 minimised over omega_m in range, omega_m*."""
    m = B.BAO["mean"][:, None]
    Ci = B.CINV
    a = np.einsum("in,ij,jn->n", t, Ci, t)
    b = np.einsum("in,ij,jn->n", t, Ci, np.broadcast_to(m, t.shape))
    c = float(B.BAO["mean"] @ Ci @ B.BAO["mean"])
    chi_fix = a - 2 * b + c
    s = np.clip(b / a, S_LO, S_HI)
    chi_opt = s * s * a - 2 * s * b + c
    wm = np.interp(s, S_GRID[::-1], WM_GRID[::-1])
    return chi_fix, chi_opt, wm


def fr(chi: np.ndarray, ref_a: float, mask=None) -> dict:
    chi = chi if mask is None else chi[mask]
    A = chi <= B.CHI2_LCDM
    Bc = chi - B.CHI2_LCDM < 4
    Ap = chi <= ref_a
    return {"n": int(len(chi)), "nA": int(A.sum()), "fA": float(A.mean()), "nB": int(Bc.sum()), "fB": float(Bc.mean()),
            "nA_prime": int(Ap.sum()), "fA_prime": float(Ap.mean()),
            "chi2_min": float(chi.min()), "chi2_p5": float(np.percentile(chi, 5)), "chi2_p50": float(np.percentile(chi, 50))}


def lcdm_desi() -> dict:
    """Flat LCDM fitted to DESI DR2 alone: Omega_m and H0 r_d free (radiation fixed at h = 0.6736)."""
    orad = B.ORAD / 0.6736**2

    def dimless(om):
        E = lambda z: np.sqrt(om * (1 + z) ** 3 + orad * (1 + z) ** 4 + 1 - om - orad)
        out = []
        for z, k in zip(B.ZR, B.KIND):
            dm = quad(lambda x: 1 / E(x), 0, z, epsrel=1e-10)[0]
            dh = 1 / E(z)
            out.append({"DM": dm, "DH": dh, "DV": (z * dm * dm * dh) ** (1 / 3)}[k])
        return np.array(out)  # times c/(H0 r_d)

    def chi_om(om):
        t = dimless(om)
        a = t @ B.CINV @ t
        b = t @ B.CINV @ B.BAO["mean"]
        amp = b / a
        return float(B.BAO["mean"] @ B.CINV @ B.BAO["mean"] - b * b / a), amp

    r = minimize_scalar(lambda om: chi_om(om)[0], bounds=(0.2, 0.4), method="bounded", options={"xatol": 1e-7})
    chi, amp = chi_om(r.x)
    return {"Omega_m": float(r.x), "H0_rd_kms": float(C / amp), "chi2_min": chi,
            "note": "amp = c/(H0 r_d) fitted analytically; 11 dof (13 points, 2 parameters)"}


def main() -> None:
    out = {"omega_m_range": WM_RANGE, "omega_b": OMB, "omega_nu": OMNU, "omega_m0": OM0, "rd0_camb": RD0,
           "rd_formula_rel_diff_vs_camb": RDJ["rd_aubourg"]["rel_diff_vs_camb"],
           "s_range": [float(S_LO), float(S_HI)]}
    L = lcdm_desi()
    out["lcdm_desi"] = L
    ref = L["chi2_min"]
    print(f"LCDM solo DESI: chi2_min {ref:.3f}  Omega_m {L['Omega_m']:.4f}  H0 r_d {L['H0_rd_kms']:.1f} km/s;"
          f"  s en [{S_LO:.4f}, {S_HI:.4f}] para omega_m en {WM_RANGE}")

    # ---------------- Model 2: regenerate the same 20000 realisations (same seed) and get distances
    mu, alpha = 0.9804, 0.8824
    L2 = B.chol_m2(mu)
    rng = np.random.default_rng(B.SEED)
    Z = rng.standard_normal((401, 20000))
    Ohat = alpha * L2 @ Z
    chi_saved = np.load(OUT / "bao_m2_bestfit.npz")["chi2"]
    Ode = np.tanh(Ohat)
    zg = np.exp(-B.LAM) - 1.0
    H = 100.0 * np.sqrt((OM0 * (1 + zg[:, None]) ** 3 + B.ORAD * (1 + zg[:, None]) ** 4) / (1.0 - Ode[:400]))
    f = C * np.exp(-B.LAM)[:, None] / H
    cum = np.concatenate([np.zeros((1, f.shape[1])), np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(B.LAM)[:, None], axis=0)])
    DMg = cum[-1][None, :] - cum
    zu = np.unique(B.ZR)
    lam_u = -np.log1p(zu)
    k = np.clip(np.searchsorted(B.LAM, lam_u) - 1, 0, len(B.LAM) - 2)
    w = ((lam_u - B.LAM[k]) / (B.LAM[k + 1] - B.LAM[k]))[:, None]
    DMu = (1 - w) * DMg[k] + w * DMg[k + 1]
    Hu = (1 - w) * H[k] + w * H[k + 1]
    H0 = H[-1]
    dmh = (DMu * H0[None, :] / C).T
    eh = (H0[None, :] / Hu).T
    t = theory13(dmh, eh, H0, zu)
    chi_fix, chi_opt, wm = chi2_free_wm(t)
    proxy = np.abs(Ode[400]) < 0.25
    out["m2_check_regenerated_vs_saved_max_abs"] = float(np.max(np.abs(chi_fix - chi_saved)))
    out["m2_fixed"] = fr(chi_fix, ref)
    # Exact treatment (radiation not rescaled, r_d by formula): the scaled chi^2 has errors up to ~1 for
    # the best realisations (radiation/matter ~1e-3 at z ~ 2.3, DESI distances at ~0.5 %), so every
    # realisation with scaled chi2_opt < 29.84 + 4 + 2 is re-minimised exactly over omega_m.
    om_keep, rd_keep = B.OM, B.R_D

    def chi_exact(jj, w_):
        B.OM = w_
        B.R_D = RD0 * rd_aubourg(OMB, w_ - OMNU, OMNU) / rd_aubourg(OMB, OM0 - OMNU, OMNU)
        ce, _, _ = B.m2_chi2_fast(Ohat[:, [jj]])
        return float(ce[0])

    cand = np.nonzero(chi_opt < B.CHI2_LCDM + 6.0)[0]
    errs = []
    chi_ex, wm_ex = chi_opt.copy(), wm.copy()
    for jj in cand:
        r_ = minimize_scalar(lambda w_: chi_exact(jj, w_), bounds=WM_RANGE, method="bounded", options={"xatol": 1e-6})
        errs.append(float(r_.fun - chi_opt[jj]))
        chi_ex[jj], wm_ex[jj] = r_.fun, r_.x
    B.OM, B.R_D = om_keep, rd_keep
    out["check_exact_m2"] = {"n_candidates_reminimised": int(len(cand)), "max_abs_dchi2_exact_minus_scaled": float(np.max(np.abs(errs))),
                             "median_dchi2": float(np.median(errs))}
    out["m2_free_wm_scaled_only"] = fr(chi_opt, ref)
    chi_opt, wm = chi_ex, wm_ex
    out["m2_free_wm"] = fr(chi_opt, ref)
    out["m2_free_wm_cmb"] = fr(chi_opt, ref, proxy)
    out["m2_free_wm_cmb"]["n_proxy"] = int(proxy.sum())
    jb = int(np.argmin(chi_opt))
    out["m2_best"] = {"index": jb, "chi2": float(chi_opt[jb]), "omega_m": float(wm[jb]),
                      "H0": float(H0[jb] * np.sqrt(wm[jb] / OM0)), "Ode0": float(Ode[399, jb]), "Ode_rec": float(Ode[400, jb])}
    pj = np.nonzero(proxy)[0][int(np.argmin(chi_opt[proxy]))]
    out["m2_best_with_cmb"] = {"index": int(pj), "chi2": float(chi_opt[pj]), "omega_m": float(wm[pj]),
                               "H0": float(H0[pj] * np.sqrt(wm[pj] / OM0)), "Ode0": float(Ode[399, pj]), "Ode_rec": float(Ode[400, pj])}
    out["m2_joint_over_all"] = {k: int(((chi_opt <= th) & proxy).sum()) for k, th in
                                (("nA_and_cmb", B.CHI2_LCDM), ("nB_and_cmb", B.CHI2_LCDM + 4), ("nAprime_and_cmb", ref))}
    pas = np.nonzero(chi_opt - B.CHI2_LCDM < 4)[0]
    out["m2_passing_B"] = [{"index": int(q), "chi2": float(chi_opt[q]), "omega_m": float(wm[q]), "cmb_proxy": bool(proxy[q]),
                            "H0": float(H0[q] * np.sqrt(wm[q] / OM0))} for q in pas[np.argsort(chi_opt[pas])]]
    within1 = np.abs(wm / OM0 - 1) <= 0.01
    out["m2_within_1pct_planck_wm"] = {"nA": int(((chi_opt <= B.CHI2_LCDM) & within1).sum()),
                                      "nB": int(((chi_opt - B.CHI2_LCDM < 4) & within1).sum()),
                                      "nA_and_cmb": int(((chi_opt <= B.CHI2_LCDM) & within1 & proxy).sum()),
                                      "nB_and_cmb": int(((chi_opt - B.CHI2_LCDM < 4) & within1 & proxy).sum())}
    out["m2_frac_at_wm_bound"] = {"low_0.10": float(np.mean(wm <= WM_RANGE[0] + 1e-5)),
                                  "high_0.20": float(np.mean(wm >= WM_RANGE[1] - 1e-5))}
    print("M2 fijo   :", out["m2_fixed"])
    print("M2 wm libre:", out["m2_free_wm"])
    print("M2 wm libre + CMB:", out["m2_free_wm_cmb"])
    print("M2 mejor con CMB:", out["m2_best_with_cmb"], " conjunta sobre 20000:", out["m2_joint_over_all"])
    print("M2 mejor:", out["m2_best"], " wm en borde:", out["m2_frac_at_wm_bound"],
          " check exacto:", out["check_exact_m2"], " regenerado vs guardado:", out["m2_check_regenerated_vs_saved_max_abs"])

    # ---------------- Model 1: needs the dimensionless distances saved by bao_sim.py m1
    for alpha in (0.01, 0.005):
        p = OUT / f"bao_m1_alpha{alpha}.npz"
        R = np.load(p)
        if "dmh" not in R:
            out[f"m1_alpha{alpha}"] = "npz without distances"
            print(f"M1 alpha {alpha}: el npz no tiene distancias")
            continue
        t = theory13(R["dmh"], R["eh"], R["H0"], R["z_unique"])
        chi_fix, chi_opt, wm = chi2_free_wm(t)
        proxy = np.abs(R["Ode_rec"]) < 0.25
        n_seeds = json.loads((OUT / f"bao_m1_alpha{alpha}.json").read_text())["n_seeds"]
        d = {"n_seeds": n_seeds, "check_fixed_vs_saved_max_abs": float(np.max(np.abs(chi_fix - R["chi2"]))),
             "fixed": fr(chi_fix, ref), "free_wm": fr(chi_opt, ref), "free_wm_cmb": fr(chi_opt, ref, proxy)}
        j = int(np.argmin(chi_opt))
        d["best"] = {"chi2": float(chi_opt[j]), "omega_m": float(wm[j]), "H0": float(R["H0"][j] * np.sqrt(wm[j] / OM0)),
                     "OmL0": float(R["OmL0"][j]), "Ode_rec": float(R["Ode_rec"][j]),
                     "seed": int(R["seed_chunk"][j]), "index": int(R["idx"][j])}
        # scaled chi^2 only (no exact re-minimisation possible without re-running): count runs within
        # +-1 of each threshold (the scaled-vs-exact error in Model 2 is <~1, check_exact_m2)
        d["borderline_pm1"] = {k: int((np.abs(chi_opt - th) < 1.0).sum()) for k, th in
                               (("A", B.CHI2_LCDM), ("B", B.CHI2_LCDM + 4), ("A_prime", ref))}
        d["joint_over_all"] = {k: int(((chi_opt <= th) & proxy).sum()) for k, th in
                               (("nA_and_cmb", B.CHI2_LCDM), ("nB_and_cmb", B.CHI2_LCDM + 4), ("nAprime_and_cmb", ref))}
        pas = np.nonzero(chi_opt - B.CHI2_LCDM < 5)[0]
        d["lowest_chi2_runs"] = [{"chi2": float(chi_opt[q]), "omega_m": float(wm[q]), "cmb_proxy": bool(proxy[q]),
                                  "H0": float(R["H0"][q] * np.sqrt(wm[q] / OM0))} for q in pas[np.argsort(chi_opt[pas])]]
        d["frac_at_wm_bound"] = {"low_0.10": float(np.mean(wm <= WM_RANGE[0] + 1e-6)),
                                 "high_0.20": float(np.mean(wm >= WM_RANGE[1] - 1e-6))}
        H0s = R["H0"] * np.sqrt(wm / OM0)
        d["H0_after_p16_p50_p84"] = [float(v) for v in np.percentile(H0s, [16, 50, 84])]
        out[f"m1_alpha{alpha}"] = d
        print(f"M1 alpha {alpha}:", json.dumps(d))
    (OUT / "robustez.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    print("Escrito resultados/robustez.json")


if __name__ == "__main__":
    main()
