#!/usr/bin/env python3
"""
PREDICTIVO_BINS.md: how anomalous is each DESI DR2 bin under LCDM (sigma_X = 0)? Leave-one-bin-out predictive check with a
look-elsewhere correction over the 7 bins by LCDM simulations. Set-up and fitters of jackknife_bins.py (imported, not modified):
exacto grid G0, data = 13 DESI DR2 numbers (official 13x13 covariance) + 100 theta_* + omega_m of Planck [+ Pantheon+, 1580 SN,
STAT+SYS]; combinations "sinSN" (parameters h, omega_m) and "Pantheon+" (h, omega_m, M).

Step 0: the noise of the excluded bin must be uncorrelated with every row used in the fit (other BAO bins, theta_*, omega_m, SN):
        max |C[bin rows, other rows]| is recorded and must be exactly 0 (otherwise V_k would need an extra term).
Step 1 (observed, full non-linear refit): for each bin k, LCDM fit to all rows but bin k with the ORIGINAL optimiser (exacto.fit /
        exacto_sn.fit_combo settings, see jackknife_bins.py) -> p0_(-k); Fisher F = J^T C^-1 J with J the finite-difference Jacobian
        of the exact model (eps 1e-4, 1e-5, 1e-3 in h, omega_m, M) on the kept rows and C their (block-diagonal) covariance;
        Sigma_(-k) = F^-1; chi^2 = r^T C^-1 r with M an explicit parameter (SN combination). Prediction of bin k: r_k = d_k -
        m_k(p0_(-k)); V_k = C_kk + J_k Sigma_(-k) J_k^T; d^2_k = r_k^T V_k^-1 r_k; p_local = chi2.sf(d^2, n_k); marginal pulls
        r_i / sqrt(V_ii) (sign kept; they ignore the correlation between the two components of a bin).
Step 2 (look-elsewhere): N = 20000 LCDM data sets = model at the full-data p0 + Gaussian noise with the full covariance of the
        combination; in each, the 7 leave-one-bin-out refits are LINEARISED (Gauss-Newton around the full p0: r_k = n_k -
        J_k Sigma_(-k) J_(-k)^T C_(-k)^-1 n_(-k), everything linear in the noise); p_min per simulation;
        p_global = (N_ext + 1) / (N + 1) with N_ext = #{p_min_sim <= p_min_obs}; MC error sqrt(p (1 - p) / N).
        Validation: in the first 20 simulations the 7 full non-linear refits are repeated and |d^2_lin - d^2_full| is recorded
        (max, simulation, bin); the run stops if any exceeds 0.1.
Step 3 (descriptive): signs of the residuals of LRG1 and LRG3+ELG1 against LRG2 in D_M/r_d and D_H/r_d.

    ../venv/bin/python -I predictivo_bins.py        # writes predictivo_bins.json
Nothing existing is modified.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.linalg import solve_triangular
from scipy.optimize import least_squares
from scipy.stats import chi2 as chi2_dist

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import jackknife_bins as JB  # noqa: E402  (set-up, fitters, model_at, build_B_at; its main() is not run)

S, X, G, MG = JB.S, JB.X, JB.G, JB.MG
SN, NSN, NROWS = JB.SN, JB.NSN, JB.NROWS
BINS, D, CFULL, LABELS = JB.BINS, JB.D, JB.CFULL, JB.LABELS
COV13 = JB.COV13
ZEFF = {"BGS": 0.295, "LRG1": 0.510, "LRG2": 0.706, "LRG3+ELG1": 0.934, "ELG2": 1.321, "QSO": 1.484, "Lya": 2.33}
N_SIM, N_VAL, SEED = 20000, 20, 20261113
EPS = (1e-4, 1e-5, 1e-3)


# ----------------------------------------------------------------------------------------------------------------------
# fitters with the data vector as an argument (same code paths and optimiser settings as jackknife_bins.fit_*)
# ----------------------------------------------------------------------------------------------------------------------
def fit_no_sn_d(keep_bao, y15):
    L = np.linalg.cholesky(COV13[np.ix_(keep_bao, keep_bao)])
    kb = np.array(keep_bao)

    def resid(p):
        v, m = X.model(p[0], p[1], np.zeros(X.N))
        if v is None:
            return np.full(len(kb) + 2, 1e3 * (1.0 - m))
        r = y15 - v
        return np.concatenate([np.linalg.solve(L, r[kb]), [r[13] / X.SIG_TH, r[14] / X.SIG_WM]])

    r = least_squares(resid, np.asarray((0.682, 0.1407)), bounds=X.BOUNDS_P, method="trf", x_scale=[0.01, 0.001], ftol=1e-10, xtol=1e-10, max_nfev=300)
    return float(2.0 * r.cost), np.array([r.x[0], r.x[1], 0.0])


def fit_with_sn_d(keep_bao, y15, mb, x0=(0.682, 0.1407)):
    L = np.linalg.cholesky(COV13[np.ix_(keep_bao, keep_bao)])
    kb = np.array(keep_bao)
    hist, _ = S.history(G, x0[1] / x0[0] ** 2, S.OMR / x0[0] ** 2, np.zeros(G.N))
    M0 = 0.0
    if hist is not None:
        rr = mb - S.mu_sn(G, SN, hist[0])
        M0 = float((SN["Ci1"] @ rr) / SN["e"])

    def resid(p):
        h, wm, M = p
        hist, m = S.history(G, wm / h**2, S.OMR / h**2, np.zeros(G.N))
        nr = len(kb) + 2 + NSN
        if hist is None:
            return np.full(nr, 1e3 * (1.0 - m))
        dm_u, inv = hist
        r = y15 - S.vec15(G, h, wm, dm_u, inv)
        rs = mb - S.mu_sn(G, SN, dm_u) - M
        return np.concatenate([solve_triangular(L, r[kb], lower=True, check_finite=False), [r[13] / X.SIG_TH, r[14] / X.SIG_WM],
                               solve_triangular(SN["L"], rs, lower=True, check_finite=False)])

    r = least_squares(resid, np.array([x0[0], x0[1], M0]), bounds=([0.4, 0.05, -50.0], [1.0, 0.30, 50.0]), method="trf",
                      x_scale=[0.01, 0.001, 0.01], ftol=1e-10, xtol=1e-10, max_nfev=200)
    return float(2.0 * r.cost), np.array(r.x)


def fit_d(combo, keep_bao, d):
    if combo == "sinSN":
        return fit_no_sn_d(keep_bao, d[:15])
    return fit_with_sn_d(keep_bao, d[:15], d[15:])


# ----------------------------------------------------------------------------------------------------------------------
# block-diagonal inverse covariance of a row set, Fisher, leave-one-out prediction
# ----------------------------------------------------------------------------------------------------------------------
def cinv_blocks(keep_bao, with_sn):
    """C^-1 of rows keep_bao + [13, 14] (+ SN) as blocks (bao inverse, cmb diag, SN inverse)."""
    return {"bao": np.linalg.inv(COV13[np.ix_(keep_bao, keep_bao)]), "cmb": np.array([1 / X.SIG_TH**2, 1 / X.SIG_WM**2]), "sn": SN["Ci"] if with_sn else None}


def apply_cinv(blocks, keep_bao, V):
    """C^-1 V for V with rows ordered as [keep_bao, 13, 14, SN...]."""
    nb = len(keep_bao)
    out = [blocks["bao"] @ V[:nb], blocks["cmb"][:, None] * V[nb:nb + 2]]
    if blocks["sn"] is not None:
        out.append(blocks["sn"] @ V[nb + 2:])
    return np.concatenate(out)


def rows_of(combo, keep_bao):
    return list(keep_bao) + [13, 14] + (list(range(15, NROWS)) if combo == "Pantheon+" else [])


def leave_one_out_obs(combo, name, d):
    """Full non-linear refit without bin `name` on data vector d; prediction of the bin; returns record (and p0)."""
    b = BINS[name]
    keep = [i for i in range(13) if i not in b]
    npar = 2 if combo == "sinSN" else 3
    chi2m, p = fit_d(combo, keep, d)
    rows = np.array(rows_of(combo, keep))
    J = JB.build_B_at(p, npar)              # NROWS x npar, finite differences of the exact model at p0_(-k)
    blocks = cinv_blocks(keep, combo == "Pantheon+")
    CJ = apply_cinv(blocks, keep, J[rows])
    F = J[rows].T @ CJ
    Sigma = np.linalg.inv(F)
    m = JB.model_at(p)
    r = (d - m)[b]
    Jk = J[b]
    V = COV13[np.ix_(b, b)] + Jk @ Sigma @ Jk.T
    d2 = float(r @ np.linalg.solve(V, r))
    return {"rows": b, "n": len(b), "chi2_minus_k": chi2m, "p0_minus_k": p.tolist(), "Sigma": Sigma.tolist(), "r": r.tolist(), "V": V.tolist(),
            "pull_marginal": (r / np.sqrt(np.diag(V))).tolist(), "d2": d2, "p_local": float(chi2_dist.sf(d2, len(b)))}


def linear_operators(combo, p0):
    """For each bin: L_k (n_k x NROWS) with r_k = L_k (d - m0) under the linearised leave-one-out refit, and V_k."""
    npar = 2 if combo == "sinSN" else 3
    J = JB.build_B_at(p0, npar)
    ops = {}
    for name, b in BINS.items():
        keep = [i for i in range(13) if i not in b]
        rows = np.array(rows_of(combo, keep))
        blocks = cinv_blocks(keep, combo == "Pantheon+")
        CJ = apply_cinv(blocks, keep, J[rows])          # C_(-k)^-1 J_(-k)
        Sigma = np.linalg.inv(J[rows].T @ CJ)
        Lk = np.zeros((len(b), NROWS))
        Lk[np.arange(len(b)), b] = 1.0
        Lk[:, rows] -= J[b] @ Sigma @ CJ.T                # r_k = n_k - J_k Sigma J_(-k)^T C^-1 n_(-k)
        V = COV13[np.ix_(b, b)] + J[b] @ Sigma @ J[b].T
        ops[name] = {"L": Lk, "V": V, "Vinv": np.linalg.inv(V), "n": len(b)}
    return ops


def noise(combo, n, rng):
    """Gaussian noise with the full covariance of the combination (block diagonal: 13x13, theta_*, omega_m, SN)."""
    L13 = np.linalg.cholesky(COV13)
    cols = NROWS if combo == "Pantheon+" else 15
    out = np.zeros((n, cols))
    out[:, :13] = rng.standard_normal((n, 13)) @ L13.T
    out[:, 13] = X.SIG_TH * rng.standard_normal(n)
    out[:, 14] = X.SIG_WM * rng.standard_normal(n)
    if combo == "Pantheon+":
        out[:, 15:] = rng.standard_normal((n, NSN)) @ SN["L"].T
    return out


def main():
    t0 = time.time()
    out = {"N_sim": N_SIM, "N_validation": N_VAL, "seed": SEED, "fisher": "F = J^T C^-1 J, J finite-difference Jacobian of the exact model (eps 1e-4, 1e-5, 1e-3) on the kept rows at the minimum, "
            "C block-diagonal covariance of the kept rows; chi^2 = r^T C^-1 r, M explicit parameter in the SN combination; same in the 7 cases and the full data",
           "step0": {}, "combos": {}}
    # ---- step 0: noise of each bin uncorrelated with everything else used
    for combo in ("sinSN", "Pantheon+"):
        rec = {}
        for name, b in BINS.items():
            keep = [i for i in range(13) if i not in b]
            others = rows_of(combo, keep)
            rec[name] = float(np.max(np.abs(CFULL[np.ix_(b, others)])))
        out["step0"][combo] = {"max_abs_cov_bin_vs_fit_rows": rec, "all_zero": bool(all(v == 0.0 for v in rec.values()))}
    print("paso 0:", json.dumps(out["step0"]), flush=True)
    assert all(out["step0"][c]["all_zero"] for c in out["step0"]), "PARA: correlación entre el bin excluido y datos del ajuste"
    # ---- check: the data-argument fitters reproduce jackknife_bins' fitters on the real data
    chk = {}
    for combo in ("sinSN", "Pantheon+"):
        a = JB.fit_no_sn(list(range(13))) if combo == "sinSN" else JB.fit_with_sn(list(range(13)))
        bb = fit_d(combo, list(range(13)), D)
        chk[combo] = {"chi2_diff": a[0] - bb[0], "p0_max_abs_diff": float(np.max(np.abs(a[1] - bb[1])))}
    out["check_fitters_vs_jackknife_bins"] = chk
    print("fitters:", chk, flush=True)
    assert all(v["p0_max_abs_diff"] == 0.0 for v in chk.values())
    rng = np.random.default_rng(SEED)
    for combo in ("sinSN", "Pantheon+"):
        t1 = time.time()
        rec = {}
        chi2f, p0 = fit_d(combo, list(range(13)), D)
        rec["p0_full"] = p0.tolist()
        rec["chi2_full"] = chi2f
        # step 1: observed
        rec["bins"] = {}
        for name in BINS:
            rec["bins"][name] = leave_one_out_obs(combo, name, D)
            rec["bins"][name]["z_eff"] = ZEFF[name]
            r = rec["bins"][name]
            print(f"{combo:9s} sin {name:10s} chi2={r['chi2_minus_k']:9.3f} d2={r['d2']:.3f} p={r['p_local']:.4f} pulls={np.round(r['pull_marginal'], 2).tolist()}", flush=True)
        p_loc_obs = {n: rec["bins"][n]["p_local"] for n in BINS}
        rec["p_min_obs"] = min(p_loc_obs.values())
        rec["bin_p_min_obs"] = min(p_loc_obs, key=p_loc_obs.get)
        # linearised operators around the full p0 (also applied to the observed data as a by-product)
        ops = linear_operators(combo, p0)
        m0 = JB.model_at(p0)
        dobs = D - m0
        cols = NROWS if combo == "Pantheon+" else 15
        rec["linearised_obs"] = {}
        for name, o in ops.items():
            rk = o["L"][:, :cols] @ dobs[:cols]
            d2 = float(rk @ o["Vinv"] @ rk)
            # internal check: Cov(r_k) = L C L^T must equal V_k
            LCL = o["L"][:, :cols] @ CFULL[:cols, :cols] @ o["L"][:, :cols].T
            rec["linearised_obs"][name] = {"d2": d2, "p_local": float(chi2_dist.sf(d2, o["n"])), "max_abs_LCLt_minus_V": float(np.max(np.abs(LCL - o["V"])))}
        # step 2: simulations
        nz = noise(combo, N_SIM, rng)
        d2_sim = np.zeros((N_SIM, len(BINS)))
        for j, (name, o) in enumerate(ops.items()):
            rk = nz @ o["L"][:, :cols].T                      # N x n_k
            d2_sim[:, j] = np.einsum("ni,ij,nj->n", rk, o["Vinv"], rk)
        p_sim = np.column_stack([chi2_dist.sf(d2_sim[:, j], o["n"]) for j, o in enumerate(ops.values())])
        p_min_sim = p_sim.min(axis=1)
        n_ext = int(np.sum(p_min_sim <= rec["p_min_obs"]))
        pg = (n_ext + 1) / (N_SIM + 1)
        rec["sims"] = {"N": N_SIM, "N_ext": n_ext, "p_global": pg, "mc_err": float(np.sqrt(pg * (1 - pg) / N_SIM)),
                       "p_local_mc_check": {n: float(np.mean(p_sim[:, j] <= p_loc_obs[n])) for j, n in enumerate(BINS)},
                       "d2_sim_mean_per_bin": d2_sim.mean(axis=0).tolist(), "p_min_sim_quantiles_1_5_50": np.percentile(p_min_sim, [1, 5, 50]).tolist(),
                       "which_bin_min_counts": {n: int(np.sum(p_sim.argmin(axis=1) == j)) for j, n in enumerate(BINS)}}
        # validation: full refits in the first N_VAL simulations
        diffs = []
        for s in range(N_VAL):
            d_s = np.concatenate([m0[:cols] + nz[s], m0[cols:]]) if combo == "sinSN" else m0 + nz[s]
            for j, name in enumerate(BINS):
                full = leave_one_out_obs(combo, name, d_s)
                diffs.append({"sim": s, "bin": name, "d2_full": full["d2"], "d2_lin": float(d2_sim[s, j]), "abs_diff": abs(full["d2"] - d2_sim[s, j])})
        worst = max(diffs, key=lambda x: x["abs_diff"])
        rec["validation"] = {"n_sims": N_VAL, "max_abs_diff_d2": worst["abs_diff"], "at_sim": worst["sim"], "at_bin": worst["bin"],
                             "max_abs_diff_d2_per_bin": {n: max(x["abs_diff"] for x in diffs if x["bin"] == n) for n in BINS}, "all": diffs, "passes_0.1": bool(worst["abs_diff"] <= 0.1)}
        print(f"{combo}: p_min_obs={rec['p_min_obs']:.4f} ({rec['bin_p_min_obs']}) N_ext={n_ext} p_global={pg:.4f}±{rec['sims']['mc_err']:.4f}; validation max |d2_lin-d2_full|={worst['abs_diff']:.2e} (sim {worst['sim']}, {worst['bin']})  [{time.time() - t1:.0f} s]", flush=True)
        out["combos"][combo] = rec
        assert worst["abs_diff"] <= 0.1, "PARA: el reajuste linealizado difiere del completo en más de 0.1 en d2"
    # step 3: signs of neighbours vs LRG2 (observed, full refit)
    st3 = {}
    for combo in ("sinSN", "Pantheon+"):
        bins = out["combos"][combo]["bins"]
        sg = {n: {"DM": float(np.sign(bins[n]["r"][0])), "DH": float(np.sign(bins[n]["r"][1]))} for n in ("LRG1", "LRG2", "LRG3+ELG1")}
        st3[combo] = {"signs": sg, "residuals": {n: bins[n]["r"] for n in ("LRG1", "LRG2", "LRG3+ELG1")},
                      "LRG1_same_sign_as_LRG2": {k: sg["LRG1"][k] == sg["LRG2"][k] for k in ("DM", "DH")},
                      "LRG3+ELG1_same_sign_as_LRG2": {k: sg["LRG3+ELG1"][k] == sg["LRG2"][k] for k in ("DM", "DH")}}
    out["step3"] = st3
    b = out["combos"]["sinSN"]["bins"]
    out["predictions"] = {"P15": {"statement": "sin SN, p local de LRG2 entre 0.003 y 0.05", "p_local_LRG2": b["LRG2"]["p_local"], "holds": bool(0.003 <= b["LRG2"]["p_local"] <= 0.05)},
                          "P16": {"statement": "sin SN, p global (peor bin de 7, calibrado con simulaciones LCDM) > 0.05", "p_global": out["combos"]["sinSN"]["sims"]["p_global"], "holds": bool(out["combos"]["sinSN"]["sims"]["p_global"] > 0.05)},
                          "P17": {"statement": "sin SN, p local de BGS > 0.05", "p_local_BGS": b["BGS"]["p_local"], "holds": bool(b["BGS"]["p_local"] > 0.05)}}
    out["runtime_s"] = time.time() - t0
    (HERE / "predictivo_bins.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")
    print(json.dumps(out["predictions"], indent=1))
    print("escrito predictivo_bins.json", f"{out['runtime_s']:.0f} s")


if __name__ == "__main__":
    main()
