#!/usr/bin/env python3
"""
Exploration only (VALVULA.md): how much the tanh "valve" of Model 2 does. bao_sim.py and rd_real.py use
Omega_de = tanh(Omega_hat), so 1 - Omega_de is in (0, 2) and H^2 >= (1/2) 100^2 [w_m (1+z)^3 + w_r (1+z)^4]: H^2 <= 0
cannot happen (Model 1 can, and "crashes"). Without the valve, Omega_de = Omega_hat; a realisation is discarded
("crash") if Omega_hat >= 1 at any stored point.

  (a) fraction of the 20000 realisations with |Omega_hat| > 1 at lambda_rec, and anywhere in lambda in [-ln 4, 0]
      (also anywhere on the 458 stored points, and the split Omega_hat >= 1 / Omega_hat <= -1);
  (b) for the 65 Model 2 survivors of supervivientes.json (the 6 Model 1 survivors have no valve) and the 5 of
      criterion (c): Omega_hat(z_rec), Omega_de(z_rec) = tanh, and the factor in H at z_rec relative to no dark energy,
      1/sqrt(1 - Omega_de), with and without the valve;
  (c) criterion (c) of rd_real.py without the valve: rd_real.model2() run with tanh replaced by the identity (its own
      code, its npz redirected to a scratch directory), after discarding the crashed realisations.
Checks: the regenerated Omega_hat reproduces bao_m2.json (chi^2 through bao_sim.m2_chi2_fast, counts A, B, proxy),
and rd_real.model2() run through the same wrapper with the real tanh reproduces resultados/rd_real_m2.npz.

Writes cmb_real/valvula.json; nothing in resultados/.
"""

from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
EV = HERE.parent
sys.path.insert(0, str(EV))

import bao_sim as B  # noqa: E402
import rd_real as RR  # noqa: E402

RES = EV / "resultados"
C5 = (1540, 8163, 11072, 12909, 16616)


class _Numpy:
    """numpy with tanh replaced (identity for 'no valve'); everything else delegated."""

    def __init__(self, valve: bool):
        self.valve = valve

    def __getattr__(self, k):
        return getattr(np, k)

    def tanh(self, x):
        return np.tanh(x) if self.valve else np.asarray(x)


def rd_real_model2(valve: bool, tmp: Path) -> tuple[dict, dict]:
    """rd_real.model2() with or without the valve; returns (its json dict, its npz arrays)."""
    np_keep, out_keep = RR.np, RR.OUT
    RR.np, RR.OUT = _Numpy(valve), tmp
    try:
        with np.errstate(invalid="ignore", divide="ignore"):
            d = RR.model2()
    finally:
        RR.np, RR.OUT = np_keep, out_keep
    z = np.load(tmp / "rd_real_m2.npz")
    return d, {k: z[k] for k in z.files}


def factor(o):
    o = np.asarray(o, dtype=float)
    return np.where(o < 1.0, 1.0 / np.sqrt(np.clip(1.0 - o, 1e-300, None)), np.nan)


def main() -> None:
    out = {}
    # ---------------- regenerate Omega_hat and check against bao_m2.json
    L = B.chol_m2(0.9804)
    Ohat = 0.8824 * L @ np.random.default_rng(B.SEED).standard_normal((401, 20000))
    chi2, H0, Orec = B.m2_chi2_fast(Ohat)
    ref = json.loads((RES / "bao_m2.json").read_text())["best_fit_zwane"]
    saved = np.load(RES / "bao_m2_bestfit.npz")["chi2"]
    proxy = np.abs(Orec) < 0.25
    out["check_bao_m2"] = {"max_abs_dchi2_vs_npz": float(np.max(np.abs(chi2 - saved))),
                           "nA": int((chi2 <= B.CHI2_LCDM).sum()), "nA_ref": ref["nA"],
                           "nB": int((chi2 - B.CHI2_LCDM < 4).sum()), "nB_ref": ref["nB"],
                           "n_cmb": int(proxy.sum()), "n_cmb_ref": ref["n_cmb"]}
    # the 57 points of rd_real.py (z in [1e3, 1e5]), same seed and construction as there
    mu, alpha = 0.9804, 0.8824
    lam_rec = B.LAM_ALL[400]
    zw = np.geomspace(1000.0, 5000.0, 10)
    lam_new = np.unique(np.concatenate([np.arange(-np.log(1e5), -np.log(1000.0) + 1e-9, 0.1), -np.log1p(zw)]))
    k = np.exp(-mu * (lam_new - lam_rec) ** 2 / 2)
    Knn = np.exp(-mu * (lam_new[:, None] - lam_new[None, :]) ** 2 / 2) - np.outer(k, k)
    w, V = np.linalg.eigh(Knn)
    Ls = V * np.sqrt(np.clip(w, 0, None))[None, :]
    Zn = np.random.default_rng(B.SEED + 1).standard_normal((len(lam_new), 20000))
    Onew = alpha * (k[:, None] * (Ohat[400][None, :] / alpha) + Ls @ Zn)
    allp = np.vstack([Ohat, Onew])

    # ---------------- (a)
    def frac(m):
        return {"n": int(m.sum()), "frac": float(m.mean())}

    out["a"] = {
        "abs_gt_1_at_rec": frac(np.abs(Ohat[400]) > 1), "ge_1_at_rec": frac(Ohat[400] >= 1),
        "le_-1_at_rec": frac(Ohat[400] <= -1),
        "abs_gt_1_somewhere_z0_3": frac(np.any(np.abs(Ohat[:400]) > 1, axis=0)),
        "ge_1_somewhere_z0_3": frac(np.any(Ohat[:400] >= 1, axis=0)),
        "abs_gt_1_somewhere_458": frac(np.any(np.abs(allp) > 1, axis=0)),
        "ge_1_somewhere_458_crash": frac(np.any(allp >= 1, axis=0)),
        "ge_1_somewhere_401_z0_3_and_rec": frac(np.any(Ohat >= 1, axis=0)),
        "sigma_Omega_hat": 0.8824, "P_abs_gt_1_gaussian_single_point": float(__import__("math").erfc(1 / 0.8824 / 2**0.5))}
    crash = np.any(allp >= 1, axis=0)

    # ---------------- (b)
    S = json.loads((HERE / "supervivientes.json").read_text())
    rows = []
    for sv in S["model2"]:
        j = sv["index"]
        oh = float(Ohat[400, j])
        rows.append({"index": j, "Omega_hat_rec": oh, "Omega_de_rec_valve": float(np.tanh(oh)),
                     "factor_H_rec_valve": float(factor(np.tanh(oh))), "factor_H_rec_no_valve": float(factor(oh)),
                     "crash_at_rec": bool(oh >= 1), "crash_anywhere": bool(crash[j]), "criterion_c": j in C5,
                     "max_Omega_hat": float(allp[:, j].max())})
    out["b_survivors_model2"] = rows
    out["b_criterion_c"] = [r for r in rows if r["criterion_c"]]

    # ---------------- (c)
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        d_v, z_v = rd_real_model2(True, tmp)
        ref_npz = np.load(RES / "rd_real_m2.npz")
        out["check_rd_real_with_valve"] = {k: bool(np.array_equal(z_v[k], ref_npz[k], equal_nan=True)) for k in ref_npz.files}
        d_n, z_n = rd_real_model2(False, tmp)
    keep = ~crash
    A, Bt = B.CHI2_LCDM, B.CHI2_LCDM + 4

    def count(z, mask_extra):
        cc = z["chi_c"]
        ok = keep & np.isfinite(cc) & mask_extra
        return {"nA": int((ok & (cc <= A)).sum()), "nB": int((ok & (cc < Bt)).sum()),
                "chi2_min": float(np.min(np.where(ok, cc, np.inf))) if ok.any() else None,
                "passing": [{"index": int(i), "chi_c_no_valve": float(cc[i]), "chi_c_valve": float(z_v["chi_c"][i]),
                             "Omega_hat_rec": float(Ohat[400, i]), "window_no_valve": bool(z["proxy_win"][i]),
                             "in_65_survivors": int(i) in surv65} for i in np.nonzero(ok & (cc < Bt))[0]]}

    ones = np.ones(20000, bool)
    surv65 = {r_["index"] for r_ in rows}
    out["c"] = {"discarded_crash": int(crash.sum()), "discarded_frac": float(crash.mean()),
                "remaining": int(keep.sum()),
                "no_valve": {"with_window": count(z_n, z_n["proxy_win"].astype(bool)), "without_window": count(z_n, ones),
                             "n_window": int((keep & z_n["proxy_win"].astype(bool)).sum())},
                "valve_reference": {"with_window": {"nA": d_v["c_rd_real_window_wm_2sigma"]["nA"], "nB": d_v["c_rd_real_window_wm_2sigma"]["nB"]},
                                    "without_window": {"nA": d_v["c_without_window"]["nA"], "nB": d_v["c_without_window"]["nB"]},
                                    "n_window": d_v["proxy"]["n_window"]},
                "c5_no_valve": {str(j): {"crash": bool(crash[j]), "chi_c": float(z_n["chi_c"][j]),
                                         "chi_c_valve": float(z_v["chi_c"][j]), "window": bool(z_n["proxy_win"][j])} for j in C5}}
    (HERE / "valvula.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    print(json.dumps({k: v for k, v in out.items() if k != "b_survivors_model2"}, indent=1))


if __name__ == "__main__":
    main()
