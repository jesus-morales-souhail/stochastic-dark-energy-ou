#!/usr/bin/env python3
"""
Checks added to MARGINAL_GUIADA.md after the main run (no figure of marginal_guiada.json is changed):
(a) Matheron sampling in the cell (theta = 1, sigma_X = 0.3, the cell of the exact bound): mean and covariance of the x samples
    against the analytic linear posterior (dp marginalised), on a subset of nodes; maximum relative errors.
(b) Linear bound at theta = 1 with the response A rebuilt at the S-table fiducial (Omega_m = 0.3153) instead of p0, to see how much
    of the 9-13 % difference with cmb_fondo.py / exacto_sn.py at theta <= 1 is the fiducial.
Writes marginal_guiada_check.json.
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import cho_factor, cho_solve

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import marginal_guiada as M  # noqa: E402

theta, sigma, N = 1.0, 0.3, M.N_SAMP
Bm, A = M.build_B(), M.build_A()
K = M.kernel(theta)
AK = A @ K
Msig = AK @ A.T
Msig = 0.5 * (Msig + Msig.T)
r = M.D - M.B0
S = M.C + sigma**2 * Msig
Sf = cho_factor(S, lower=True)
SB = cho_solve(Sf, Bm)
Sp = np.linalg.inv(Bm.T @ SB)
dp_hat = Sp @ (Bm.T @ cho_solve(Sf, r))
# analytic marginal posterior of x (dp marginalised): mean mu = sigma^2 K A' S^-1 (r - B dp_hat);
# cov = sigma^2 K - sigma^4 K A' S^-1 A K + sigma^4 K A' S^-1 B Sp B' S^-1 A K
G = sigma**2 * cho_solve(Sf, AK).T  # (N_grid x 1595) = sigma^2 K A' S^-1
mu = G @ (r - Bm @ dp_hat)
sub = np.unique(np.concatenate([np.arange(0, M.G.N, 50), M.G.IDX, [M.G.N - 1]]))
Gs = G[sub]
cov = sigma**2 * K[np.ix_(sub, sub)] - Gs @ (sigma**2 * AK[:, sub]) + Gs @ Bm @ Sp @ Bm.T @ Gs.T
# samples with Matheron (fresh seed)
rng = np.random.default_rng(M.SEED + 1)
dp = dp_hat[None, :] + rng.standard_normal((N, 3)) @ np.linalg.cholesky(Sp).T
x0 = M.S.ou_paths(M.G, sigma, theta, N, rng)
n0 = rng.standard_normal((N, M.NROWS)) @ M.LC.T
R = r[None, :] - dp @ Bm.T - x0 @ A.T - n0
x = x0 + sigma**2 * cho_solve(Sf, R.T).T @ AK
xs = x[:, sub]
m_emp = xs.mean(axis=0)
c_emp = np.cov(xs.T)
sd = np.sqrt(np.diag(cov))
err_mean = np.abs(m_emp - mu[sub]) / sd  # in units of the posterior std
err_mean_rel = np.abs(m_emp - mu[sub]) / np.maximum(np.abs(mu[sub]), 1e-12)
err_cov_diag = np.abs(np.diag(c_emp) / np.diag(cov) - 1)
err_cov_frob = np.linalg.norm(c_emp - cov) / np.linalg.norm(cov)
corr_an = cov / np.outer(sd, sd)
corr_emp = c_emp / np.outer(np.sqrt(np.diag(c_emp)), np.sqrt(np.diag(c_emp)))
out = {"cell": {"theta": theta, "sigma_X": sigma, "n_samples": N, "n_nodes_checked": int(len(sub))},
       "mean": {"max_abs_err_in_posterior_std": float(err_mean.max()), "rms_err_in_posterior_std": float(np.sqrt(np.mean(err_mean**2))),
                "expected_rms_mc": float(1 / np.sqrt(N)), "max_rel_err_vs_mean_where_mean_gt_0.01": float(err_mean_rel[np.abs(mu[sub]) > 0.01].max()),
                "posterior_mean_range": [float(mu.min()), float(mu.max())], "posterior_std_range": [float(sd.min()), float(sd.max())]},
       "cov": {"max_rel_err_diag": float(err_cov_diag.max()), "rms_rel_err_diag": float(np.sqrt(np.mean(err_cov_diag**2))), "expected_rms_mc_diag": float(np.sqrt(2 / N)),
               "frobenius_rel_err": float(err_cov_frob), "max_abs_err_corr": float(np.max(np.abs(corr_emp - corr_an)))}}
print(json.dumps(out, indent=1), flush=True)
# (b) fiducial variant of A at Omega_m = 0.3153 (S-table), theta = 1 and 0.3
res_fid = {}
P0_saved = M.P0.copy()
h_fid = M.P0[0]
wm_fid = (1.0 - M.V.P13.OL_FID if hasattr(M, "V") else 0.3153) * h_fid**2
M.P0 = np.array([h_fid, wm_fid, P0_saved[2]])
A_fid = M.build_A()
M.P0 = P0_saved
for th in (0.3, 1.0, 3.0, 10.0):
    Kt = M.kernel(th)
    out_t = {}
    for lab, AA in (("A_p0", A), ("A_fid_0.3153", A_fid)):
        Mt = AA @ Kt @ AA.T
        Mt = 0.5 * (Mt + Mt.T)
        sig = [0.01, 0.03, 0.05, 0.07, 0.1, 0.15, 0.2, 0.3, 0.4] + ([0.5, 0.7] if th == 0.3 else [])
        lz = [M.lnZ_lin(cho_factor(M.C + s * s * Mt, lower=True), Bm, r) for s in sig]
        lz0 = M.lnZ_lin(cho_factor(M.C, lower=True), Bm, r)
        out_t[lab] = M.bounds_from(np.array(sig), np.array(lz) - lz0, 1.92)
    res_fid[f"theta_{th:g}"] = out_t
    print(th, out_t, flush=True)
out["fiducial_variant"] = {"note": "linear bound -1.92 with A built at p0 (canonical) and with A built at the S-table fiducial Omega_m = 0.3153 (same grid)",
                           "omega_m_fid_used": float(wm_fid), "Omega_m_fid": float(wm_fid / h_fid**2), "Omega_m_p0": float(P0_saved[1] / P0_saved[0] ** 2), "bounds": res_fid}
(HERE / "marginal_guiada_check.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
