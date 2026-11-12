#!/usr/bin/env python3
"""
Check of MARGINAL_GUIADA.md: the guided importance sampling repeated at theta = 1 and 0.3 with the A of cmb_fondo.py /
exacto_sn.py (uniform grid, nearest-node BAO rows, S-table fiducial, scaling at the BAO + CMB fit) as the guide instead of the
canonical A. The identity Z_exact = Z_lin E_q[w] holds for any linear proposal, so the exact bound must come out the same
(0.287 and 0.500) within the bootstrap error; only ln Z_lin, the correction ln <w> and the ESS may differ.
The path lives on the uniform grid of cmb_fondo.py (3499 nodes, the same uniform grid as exacto.py's before the 7 DESI nodes
are inserted); for the exact chi^2 it is interpolated linearly onto the exacto grid (only the 7 inserted nodes are new).
Writes marginal_guiada_cf.json.
"""
import json
import sys
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.linalg import cho_factor, cho_solve

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import marginal_guiada as M  # noqa: E402
import cmb_fondo as CF  # noqa: E402
import exacto_sn as S  # noqa: E402
import precision_sn as PS  # noqa: E402

t0 = time.time()
R = CF.Response(lambda_absorbs_x0=True)
u, z = R.u, R.z
assert np.allclose(u, S.X._u), "uniform grids differ"
om, orad = CF.OM_FID, M.OMR / CF.H_FID**2
E2 = om * (1 + z) ** 3 + orad * (1 + z) ** 4 + (1 - om - orad)
Es = np.sqrt(E2)
base = R.du * (1 + z)
N = len(u)
sn = M.SN


def wM_at(k):
    cw = np.zeros(N)
    cw[: k + 1] = base[: k + 1]
    cw[0] *= 0.5
    cw[k] *= 0.5
    dm = float(np.sum(cw / Es))
    return -cw / (2 * Es**3) / dm


Wsn = np.zeros((sn["n"], N))
for i, uu in enumerate(sn["uhd"]):
    k = int(np.searchsorted(u, uu)) - 1
    t = (uu - u[k]) / (u[k + 1] - u[k])
    row = (1 - t) * wM_at(k) + t * wM_at(k + 1)
    row[0] -= row.sum()
    Wsn[i] = M.MAG * row
fit15 = CF.Fit(CF.SETS["BAO+CMB"], R)
scale = np.concatenate([fit15.scale[:15], np.ones(sn["n"])])
A_cf = np.vstack([R.W[:15], Wsn]) * scale[:, None]
print("A_cf built", A_cf.shape, "rows sum max", float(np.max(np.abs(A_cf.sum(axis=1)))), f"[{time.time() - t0:.0f} s]", flush=True)
Bm = M.build_B()
Gu = SimpleNamespace(U=u, DUK=np.diff(u), N=N)
# interpolation of the uniform-grid path onto the exacto grid: indices and weights
ku = np.clip(np.searchsorted(u, M.G.U) - 1, 0, N - 2)
tu = np.clip((M.G.U - u[ku]) / (u[ku + 1] - u[ku]), 0.0, 1.0)


def to_g0(x):
    return (1 - tu)[None, :] * x[:, ku] + tu[None, :] * x[:, ku + 1]


def cell(theta, sigma, K, AK, Msig, rng):
    r = M.D - M.B0
    Sf = cho_factor(M.C + sigma**2 * Msig, lower=True)
    lz = M.lnZ_lin(Sf, Bm, r)
    SB = cho_solve(Sf, Bm)
    Sp = np.linalg.inv(Bm.T @ SB)
    dp_hat = Sp @ (Bm.T @ cho_solve(Sf, r))
    dp = dp_hat[None, :] + rng.standard_normal((M.N_SAMP, 3)) @ np.linalg.cholesky(Sp).T
    x0 = S.ou_paths(Gu, sigma, theta, M.N_SAMP, rng) if sigma > 0 else np.zeros((M.N_SAMP, N))
    n0 = rng.standard_normal((M.N_SAMP, M.NROWS)) @ M.LC.T
    Rr = r[None, :] - dp @ Bm.T - x0 @ A_cf.T - n0
    x = x0 + sigma**2 * cho_solve(Sf, Rr.T).T @ AK if sigma > 0 else x0
    chi_lin = M.chi2_of(M.B0[None, :] + dp @ Bm.T + x @ A_cf.T)
    p = M.P0[None, :] + dp
    xg = to_g0(x)
    m_ex, bad = M.exact_model(p[:, 0], p[:, 1], p[:, 2], xg - xg[:, [0]])
    chi_ex = np.where(bad, np.inf, M.chi2_of(np.where(bad[:, None], 0.0, m_ex)))
    w = np.where(bad, 0.0, np.exp(np.minimum(-0.5 * (chi_ex - chi_lin), 700)))
    idx = rng.integers(0, M.N_SAMP, size=(M.N_BOOT, M.N_SAMP))
    return {"sigma": sigma, "lnZ_lin": lz, "ln_mean_w": float(np.log(max(w.mean(), 1e-300))), "ln_mean_w_boot": np.log(np.maximum(w[idx].mean(axis=1), 1e-300)).tolist(),
            "ess": float(w.sum() ** 2 / np.sum(w**2)), "w_max_frac": float(w.max() / w.sum()), "frac_H2_le_0": float(bad.mean())}


ref = json.load(open(HERE / "marginal_guiada.json"))["thetas"]
out = {"guide": "A of cmb_fondo.py / exacto_sn.py (uniform grid, S-table fiducial, nearest-node BAO rows, scaling at the BAO+CMB fit)", "thetas": {}}
rng = np.random.default_rng(M.SEED + 7)
for theta in (1.0, 0.3):
    t1 = time.time()
    K = np.exp(-theta * np.abs(u[:, None] - u[None, :]))
    AK = A_cf @ K
    Msig = AK @ A_cf.T
    Msig = 0.5 * (Msig + Msig.T)
    sigs = (0.0, 0.01) + M.SIG_MAIN + ((0.5, 0.7) if theta == 0.3 else ())
    cells = [cell(theta, s, K, AK, Msig, rng) for s in sigs]
    sig = np.array([c["sigma"] for c in cells])
    lz = np.array([c["lnZ_lin"] for c in cells])
    lnw = np.array([c["ln_mean_w"] for c in cells])
    lz_lin_rel, lz_ex_rel = lz - lz[0], (lz + lnw) - (lz[0] + lnw[0])
    rec = {"sigmas": sig.tolist(), "lnZ_lin_rel": lz_lin_rel.tolist(), "ln_mean_w": lnw.tolist(), "lnZ_exact_rel": lz_ex_rel.tolist(),
           "ess": [c["ess"] for c in cells], "w_max_frac": [c["w_max_frac"] for c in cells], "frac_H2_le_0": [c["frac_H2_le_0"] for c in cells],
           "lin_192": M.bounds_from(sig[1:], lz_lin_rel[1:], 1.92), "exact_192": M.bounds_from(sig[1:], lz_ex_rel[1:], 1.92), "exact_135": M.bounds_from(sig[1:], lz_ex_rel[1:], 1.35)}
    bb = []
    for b in range(M.N_BOOT):
        lnw_b = np.array([c["ln_mean_w_boot"][b] for c in cells])
        bb.append(M.bounds_from(sig[1:], ((lz + lnw_b) - (lz[0] + lnw_b[0]))[1:], 1.92))
    bb = np.array(bb)
    fin = np.isfinite(bb) & (bb > 0)
    rec["exact_192_boot_p16_p84"] = [float(np.percentile(bb[fin], 16)), float(np.percentile(bb[fin], 84))]
    k = int(np.argmin(np.abs(np.log(sig[1:] / rec["exact_192"])))) + 1
    rec["at_exact_192"] = {"sigma_cell": float(sig[k]), "ess": cells[k]["ess"], "ess_frac_N": cells[k]["ess"] / M.N_SAMP, "correction": float(lnw[k] - lnw[0]), "w_max_frac": cells[k]["w_max_frac"]}
    rr = ref[f"theta_{theta:g}"]
    rec["canonical"] = {"exact_192": rr["exact_192"], "boot_p16_p84": rr["exact_192_boot_p16_p84"], "lin_192": rr["lin_192"], "ess_at_bound": rr["at_exact_192"]["ess"], "correction": rr["at_exact_192"]["correction"]}
    lo = min(rec["exact_192_boot_p16_p84"][0], rr["exact_192_boot_p16_p84"][0]); hi = max(rec["exact_192_boot_p16_p84"][1], rr["exact_192_boot_p16_p84"][1])
    rec["agree_within_bootstrap"] = bool(abs(rec["exact_192"] - rr["exact_192"]) <= (rec["exact_192_boot_p16_p84"][1] - rec["exact_192_boot_p16_p84"][0]) / 2 + (rr["exact_192_boot_p16_p84"][1] - rr["exact_192_boot_p16_p84"][0]) / 2)
    rec["seconds"] = time.time() - t1
    out["thetas"][f"theta_{theta:g}"] = rec
    print(f"theta {theta:g} [{rec['seconds']:.0f} s]: lin192(A_cf)={rec['lin_192']:.3f} exact192={rec['exact_192']:.3f} {rec['exact_192_boot_p16_p84']} canonical {rr['exact_192']:.3f} {rr['exact_192_boot_p16_p84']} "
          f"| ESS at bound {rec['at_exact_192']['ess']:.0f} ({100 * rec['at_exact_192']['ess_frac_N']:.1f} %), correction {rec['at_exact_192']['correction']:+.3f} (canonical {rr['at_exact_192']['correction']:+.3f}) | agree: {rec['agree_within_bootstrap']}", flush=True)
    for c, a, b2 in zip(cells, lz_lin_rel, lz_ex_rel):
        print(f"   sigma={c['sigma']:<5g} lnZlin={a:+8.3f} ln<w>={c['ln_mean_w']:+7.3f} lnZex={b2:+8.3f} ess={c['ess']:7.0f} wmax={c['w_max_frac']:.3f} H2<=0={c['frac_H2_le_0']:.4f}", flush=True)
(HERE / "marginal_guiada_cf.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
