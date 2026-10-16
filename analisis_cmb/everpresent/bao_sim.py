#!/usr/bin/env python3
"""
Everpresent Lambda against DESI DR2 BAO by simulation (statistic of Das II: fraction of realisations
compatible with the data).

Data: DESI DR2 Gaussian BAO, 13 distances D/r_d, full 13x13 covariance, as in oscilador/
(loader: local copy desi_dr2_data.load_gaussian_bao_vector). r_d fixed to Planck 2018 LCDM
(147.103 Mpc, resultados/rd_planck.json) for every realisation: an assumption, equivalent to
conditioning on realisations that fit the CMB acoustic scale. Reference chi^2(LCDM-Planck) = 29.84
(oscilador/resultados/oscilador.json, CAMB background); the analytic distances used here give 29.80
for the same LCDM (rd_planck.json).

Criteria per realisation: (A) chi^2 <= 29.84;  (B) chi^2 - 29.84 < 4.
CMB proxy (Zwane's criterion): |Omega_de(z = 1100)| < 0.25.

  m2        Zwane Model 2 at (alpha, mu) = (0.8824, 0.9804): 20000 realisations, exact Cholesky on 400
            points in lambda = ln a in [-ln 4, 0] plus lambda_rec = -ln 1101; Omega_de = tanh(Omega_hat).
            H^2(z) = 100^2 [om (1+z)^3 + orad (1+z)^4] / (1 - Omega_de(z)).
            Then the alpha grid {0.02, 0.05, 0.1, 0.2, 0.4, 0.8824}, 2000 realisations each, same normals.
  m1 ALPHA NSEEDS
            Das Model 1 (ahmed_sim.run_das, n = 1e5, start a = 1e-4, S_0 = 0, K = 8 pi / 3) with
            omega_m = 0.14301 (Planck incl. nu, consistent with r_d); chunks of 20000 seeds with
            seed = 20261016 + chunk index. H0 from the simulation.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

from desi_dr2_data import load_gaussian_bao_vector

HERE = Path(__file__).resolve().parent
OUT = HERE / "resultados"
C = 299792.458
RD = json.loads((OUT / "rd_planck.json").read_text())
R_D = RD["rdrag_Mpc"]
OM = RD["omega_m_incl_nu"]
ORAD = RD["omega_r"]
CHI2_LCDM = json.loads((HERE.parent / "oscilador" / "resultados" / "oscilador.json").read_text())["lcdm"]["chi2_bao"]
BAO = load_gaussian_bao_vector()
CINV = np.linalg.inv(BAO["cov"])
ZR = BAO["z_rows"]
KIND = [lab[:2] for lab in BAO["labels"]]
SEED = 20261016
Z_REC = 1100.0


def chi2_from(DM: np.ndarray, DH: np.ndarray) -> np.ndarray:
    """DM, DH: (13, n) in Mpc at the 13 rows. Returns chi^2 per realisation."""
    DV = (ZR[:, None] * DM**2 * DH) ** (1.0 / 3.0)
    th = np.where(np.array(KIND)[:, None] == "DM", DM, np.where(np.array(KIND)[:, None] == "DH", DH, DV)) / R_D
    d = th - BAO["mean"][:, None]
    return np.einsum("in,ij,jn->n", d, CINV, d)


def fractions(chi2: np.ndarray, proxy: np.ndarray | None = None) -> dict:
    A = chi2 <= CHI2_LCDM
    B = chi2 - CHI2_LCDM < 4.0
    out = {"n": int(len(chi2)), "chi2_p5": float(np.percentile(chi2, 5)), "chi2_p50": float(np.percentile(chi2, 50)),
           "chi2_p95": float(np.percentile(chi2, 95)), "chi2_min": float(np.min(chi2)),
           "nA": int(A.sum()), "fA": float(A.mean()), "nB": int(B.sum()), "fB": float(B.mean())}
    if proxy is not None:
        out.update(n_cmb=int(proxy.sum()), f_cmb=float(proxy.mean()), nA_cmb=int((A & proxy).sum()),
                   fA_cmb=float((A & proxy).mean()), nB_cmb=int((B & proxy).sum()), fB_cmb=float((B & proxy).mean()))
    return out


# ------------------------------------------------------------------ Model 2 (Zwane)
LAM = np.linspace(-np.log(4.0), 0.0, 400)
LAM_ALL = np.append(LAM, -np.log(1.0 + Z_REC))


def chol_m2(mu: float, jitter: float = 1e-10) -> np.ndarray:
    d = LAM_ALL[:, None] - LAM_ALL[None, :]
    K = np.exp(-mu * d**2 / 2.0)
    return np.linalg.cholesky(K + jitter * np.eye(len(K)))


def m2_chi2(Ohat: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Ohat: (401, n). Returns chi2, H0, Omega_de(z_rec)."""
    Ode = np.tanh(Ohat)
    zg = np.exp(-LAM) - 1.0
    H = 100.0 * np.sqrt((OM * (1 + zg[:, None]) ** 3 + ORAD * (1 + zg[:, None]) ** 4) / (1.0 - Ode[:400]))
    f = C * np.exp(-LAM)[:, None] / H  # dD_M/dlambda (with sign: D_M = int_lambda^0 f)
    cum = np.concatenate([np.zeros((1, f.shape[1])), np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(LAM)[:, None], axis=0)])
    DM_grid = cum[-1][None, :] - cum  # int from lambda to 0
    lam_r = -np.log1p(ZR)
    DM = np.array([[np.interp(l, LAM, DM_grid[:, j]) for j in range(f.shape[1])] for l in lam_r])
    Hz = np.array([[np.interp(l, LAM, H[:, j]) for j in range(f.shape[1])] for l in lam_r])
    return chi2_from(DM, C / Hz), H[-1], Ode[400]


def m2_chi2_fast(Ohat: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Vectorised interpolation (linear in lambda) of m2_chi2."""
    Ode = np.tanh(Ohat)
    zg = np.exp(-LAM) - 1.0
    H = 100.0 * np.sqrt((OM * (1 + zg[:, None]) ** 3 + ORAD * (1 + zg[:, None]) ** 4) / (1.0 - Ode[:400]))
    f = C * np.exp(-LAM)[:, None] / H
    cum = np.concatenate([np.zeros((1, f.shape[1])), np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(LAM)[:, None], axis=0)])
    DM_grid = cum[-1][None, :] - cum
    lam_r = -np.log1p(ZR)
    k = np.clip(np.searchsorted(LAM, lam_r) - 1, 0, len(LAM) - 2)
    w = ((lam_r - LAM[k]) / (LAM[k + 1] - LAM[k]))[:, None]
    DM = (1 - w) * DM_grid[k] + w * DM_grid[k + 1]
    Hz = (1 - w) * H[k] + w * H[k + 1]
    return chi2_from(DM, C / Hz), H[-1], Ode[400]


def main_m2() -> None:
    out = {"assumptions": {"r_d_Mpc": R_D, "omega_m": OM, "omega_r": ORAD, "chi2_lcdm_ref": CHI2_LCDM,
                           "cholesky_jitter": 1e-10, "lambda_grid": "400 points in [-ln 4, 0] + -ln 1101",
                           "H0": "from H^2(0) = 100^2 (om + orad)/(1 - Omega_de(0))"}}
    mu, alpha = 0.9804, 0.8824
    L = chol_m2(mu)
    rng = np.random.default_rng(SEED)
    # check of the cross-correlation rec <-> lambda = -ln 4 (independence)
    out["corr_rec_vs_z3"] = float(np.exp(-mu * (LAM_ALL[400] - LAM[0]) ** 2 / 2))
    # direct check of fast vs explicit interpolation on 50 realisations
    zt = rng.standard_normal((401, 50))
    c1, _, _ = m2_chi2(alpha * L @ zt)
    c2, _, _ = m2_chi2_fast(alpha * L @ zt)
    out["check_interp_max_abs_diff"] = float(np.max(np.abs(c1 - c2)))

    # pipeline checks: LCDM-Planck history through the same code, and Omega_de = 0
    h2 = 0.6736**2
    zall = np.exp(-LAM_ALL) - 1.0
    ode_l = (h2 - OM - ORAD) / (OM * (1 + zall) ** 3 + ORAD * (1 + zall) ** 4 + h2 - OM - ORAD)
    c_l, h_l, _ = m2_chi2_fast(np.arctanh(ode_l)[:, None])
    c_0, h_0, _ = m2_chi2_fast(np.zeros((401, 1)))
    out["check_lcdm_history"] = {"chi2": float(c_l[0]), "H0": float(h_l[0])}
    out["check_no_dark_energy"] = {"chi2": float(c_0[0]), "H0": float(h_0[0])}
    print(f"checks: LCDM por este pipeline chi2 {c_l[0]:.3f} (H0 {h_l[0]:.2f}); Omega_de = 0: chi2 {c_0[0]:.0f} (H0 {h_0[0]:.2f})")

    rng = np.random.default_rng(SEED)
    Z = rng.standard_normal((401, 20000))
    Ohat = alpha * L @ Z
    chi2, H0, Orec = m2_chi2_fast(Ohat)
    for j in np.nonzero(chi2 - CHI2_LCDM < 4.0)[0]:
        out.setdefault("passing_realisations", []).append(
            {"index": int(j), "chi2": float(chi2[j]), "H0": float(H0[j]), "Ode0": float(np.tanh(Ohat[399, j])),
             "Ode_z1": float(np.tanh(np.interp(-np.log(2.0), LAM, Ohat[:400, j]))), "Ode_rec": float(Orec[j])})
    proxy = np.abs(Orec) < 0.25
    out["best_fit_zwane"] = {"alpha": alpha, "mu": mu, **fractions(chi2, proxy),
                             "H0_p16_p50_p84": [float(v) for v in np.percentile(H0, [16, 50, 84])]}
    np.savez_compressed(OUT / "bao_m2_bestfit.npz", chi2=chi2, H0=H0, Ode_rec=Orec)
    r = out["best_fit_zwane"]
    print(f"M2 (0.8824, 0.9804), 20000: chi2 p5/p50/p95 {r['chi2_p5']:.1f}/{r['chi2_p50']:.1f}/{r['chi2_p95']:.1f}"
          f" min {r['chi2_min']:.2f}; A {r['nA']} ({r['fA']:.4%}); B {r['nB']} ({r['fB']:.4%}); CMB-proxy {r['n_cmb']}"
          f" ({r['f_cmb']:.3%}); A&CMB {r['nA_cmb']}; B&CMB {r['nB_cmb']} ({r['fB_cmb']:.4%})")

    grid = [0.02, 0.05, 0.1, 0.2, 0.4, 0.8824]
    Z2 = Z[:, :2000]
    LZ = L @ Z2
    out["alpha_grid"] = {}
    chis = {}
    for a in grid:
        c, h0, orc = m2_chi2_fast(a * LZ)
        chis[f"alpha_{a}"] = c
        out["alpha_grid"][str(a)] = fractions(c, np.abs(orc) < 0.25)
        g = out["alpha_grid"][str(a)]
        print(f"  alpha {a:<6}: chi2 p50 {g['chi2_p50']:8.1f}  fA {g['fA']:.4f}  fB {g['fB']:.4f}  f_cmb {g['f_cmb']:.3f}")
    np.savez_compressed(OUT / "bao_m2_alphagrid.npz", **chis)
    for key in ("fA", "fB"):
        f = np.array([out["alpha_grid"][str(a)][key] for a in grid])
        for thr in (0.05, 0.001):
            j = np.nonzero(f < thr)[0]
            if not len(j):
                out[f"alpha_{key}_below_{thr}"] = None
                continue
            j = j[0]
            if j == 0:
                out[f"alpha_{key}_below_{thr}"] = {"below_at_first_grid_point": grid[0]}
                continue
            # log-log interpolation between the bracketing grid points (a zero count is floored at 0.5/2000)
            la = np.interp(np.log(thr), [np.log(max(f[j], 0.5 / 2000)), np.log(f[j - 1])],
                           [np.log(grid[j]), np.log(grid[j - 1])])
            out[f"alpha_{key}_below_{thr}"] = {"bracket": [grid[j - 1], grid[j]], "loglog_interp": float(np.exp(la))}
            print(f"  {key} < {thr}: entre alpha {grid[j - 1]} y {grid[j]} (interp log-log {np.exp(la):.3g})")
    (OUT / "bao_m2.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    print("Escrito resultados/bao_m2.json, bao_m2_bestfit.npz, bao_m2_alphagrid.npz")


# ------------------------------------------------------------------ Model 1 (Das)
def main_m1(alpha: float, nseeds: int, chunk: int = 20000) -> None:
    from ahmed_sim import run_das

    zu = np.unique(ZR)
    z_obs = np.append(zu, Z_REC)
    rows = {k: [] for k in ("chi2", "H0", "OmL0", "Ode_rec", "seed_chunk", "idx", "dmh", "eh")}
    n_done = n_alive = 0
    t0 = time.time()
    for c in range(int(np.ceil(nseeds / chunk))):
        m = min(chunk, nseeds - n_done)
        r = run_das(alpha, m, seed=SEED + c, om=OM, orad=ORAD, n=100000, i_ini=10, z_obs=z_obs)
        ok = r["alive"]
        n_done += m
        n_alive += int(ok.sum())
        a = r["obs_a"][:, None]
        tau, u = r["obs_tau"][:, ok], r["obs_u"][:, ok]
        H = 100.0 * np.sqrt(OM / a**3 + ORAD / a**4 + u)
        DM_u = (C / 100.0) * (r["tau0"][ok][None, :] - tau[:-1])
        DH_u = C / H[:-1]
        iz = np.searchsorted(zu, ZR)
        chi2 = chi2_from(DM_u[iz], DH_u[iz])
        rows["chi2"].append(chi2)
        rows["H0"].append(r["H0"][ok])
        rows["OmL0"].append(r["OmL0"][ok])
        rows["Ode_rec"].append(u[-1] / (H[-1] / 100.0) ** 2)
        rows["seed_chunk"].append(np.full(ok.sum(), SEED + c))
        rows["idx"].append(np.nonzero(ok)[0])
        # dimensionless distances at the unique DESI z (invariant under omega_m rescaling, Das I Sec. 4.4):
        # D_M H0 / c and H0 / H(z); shape (n_complete, 7)
        H0c = r["H0"][ok][None, :]
        rows["dmh"].append((DM_u * H0c / C).T)
        rows["eh"].append((H0c / H[:-1]).T)
        print(f"  chunk {c}: {m} semillas, completan {ok.sum()}, acumulado {n_done} en {time.time() - t0:.0f} s", flush=True)
    R = {k: np.concatenate(v) for k, v in rows.items()}
    R["z_unique"] = zu
    np.savez_compressed(OUT / f"bao_m1_alpha{alpha}.npz", **R)
    proxy = np.abs(R["Ode_rec"]) < 0.25
    out = {"alpha": alpha, "n_seeds": n_done, "n_complete": n_alive, "frac_complete": n_alive / n_done,
           "omega_m": OM, "n_steps": 100000, "a_start": 1e-4, "seconds": time.time() - t0,
           "unconditional_over_seeds": {"nA": int((R["chi2"] <= CHI2_LCDM).sum()),
                                        "nB": int((R["chi2"] - CHI2_LCDM < 4).sum()),
                                        "fA": float((R["chi2"] <= CHI2_LCDM).sum() / n_done),
                                        "fB": float((R["chi2"] - CHI2_LCDM < 4).sum() / n_done)},
           "over_completed": fractions(R["chi2"], proxy),
           "H0_complete_p16_p50_p84": [float(v) for v in np.percentile(R["H0"], [16, 50, 84])]}
    win = (R["H0"] >= 63) & (R["H0"] <= 77)
    out["H0_window_63_77"] = {"n": int(win.sum())}
    if win.sum() >= 20:
        out["H0_window_63_77"].update(fractions(R["chi2"][win], proxy[win]))
    else:
        out["H0_window_63_77"]["note"] = "fewer than 20 runs: count only"
        out["H0_window_63_77"]["chi2_values"] = [float(v) for v in np.sort(R["chi2"][win])]
    j = int(np.argmin(R["chi2"]))
    out["best"] = {"chi2": float(R["chi2"][j]), "H0": float(R["H0"][j]), "OmL0": float(R["OmL0"][j]),
                   "Ode_rec": float(R["Ode_rec"][j]), "seed": int(R["seed_chunk"][j]), "index": int(R["idx"][j])}
    (OUT / f"bao_m1_alpha{alpha}.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps({k: v for k, v in out.items() if k != "H0_window_63_77" or True}, indent=1)[:3000])


if __name__ == "__main__":
    if sys.argv[1] == "m2":
        main_m2()
    else:
        main_m1(float(sys.argv[2]), int(sys.argv[3]))
