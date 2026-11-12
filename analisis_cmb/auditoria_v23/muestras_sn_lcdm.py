#!/usr/bin/env python3
"""
Closure of MUESTRAS_SN.md, point 2b: chi^2_min of flat LCDM (X = 0) for each SN sample combined with BAO + theta_* + omega_m,
its degrees of freedom and p-value, the split of chi^2_min between the 15 rows and the SN block, and the tension statistic
Delta chi^2 = chi^2_min(joint) - chi^2_min(15 rows alone) - chi^2_min(SN alone), compared with the loosening of the 5 % rule and
of the linear rule read from muestras_sn.json. No bound is recomputed: the joint best fits are the ones stored in muestras_sn.json
(lcdm) and the 15-row best fit is exacto_sn.P0; everything here is an evaluation of chi^2 at those points plus the SN-only
one-parameter profile already used for the Omega_m check.

    ../venv/bin/python -I muestras_sn_lcdm.py      # writes muestras_sn_lcdm.json
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.stats import chi2 as chi2_dist

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import exacto_sn as S  # noqa: E402
import muestras_sn as MS  # noqa: E402

G0 = S.G0


def main() -> None:
    res = json.loads((HERE / "muestras_sn.json").read_text(encoding="utf-8"))
    samples = MS.load_samples()
    zero = np.zeros(G0.N)
    out = {}
    # 15 rows alone (BAO + theta_* + omega_m), 2 parameters (h, omega_m): exacto_sn.P0, re-minimised here as a check
    sn_dummy = samples["union3"]

    def chi2_15(p):
        r = S.resid_combo([p[0], p[1], 0.0], G0, zero, sn_dummy, True)
        return float(np.sum(r[:15] ** 2))

    m15 = minimize(chi2_15, np.array(S.P0), method="Nelder-Mead", options={"xatol": 1e-8, "fatol": 1e-10})
    c15_min = float(m15.fun)
    out["rows15_alone"] = {"chi2_min": c15_min, "h": float(m15.x[0]), "omega_m": float(m15.x[1]), "chi2_at_exacto_P0": chi2_15(S.P0),
                           "n": 15, "n_par": 2, "dof": 13, "p": float(chi2_dist.sf(c15_min, 13))}
    for samp, sn in samples.items():
        lc = res["lcdm"][samp]
        p = [lc["H0"] / 100.0, lc["omega_m"], lc["M"]]
        r = S.resid_combo(p, G0, zero, sn, True)
        c_tot, c_15, c_sn = float(np.sum(r**2)), float(np.sum(r[:15] ** 2)), float(np.sum(r[15:] ** 2))
        n = 15 + sn["n"]
        dof = n - 3
        c_sn_min, om_sn = S.fit_sn_only(G0, zero, sn)
        d_ten = c_tot - c15_min - c_sn_min
        rat = {k: v for k, v in res["ratios_vs_pantheon"].items() if k.startswith(samp + "|")}
        r5 = [v["5pct"] for v in rat.values() if v.get("5pct") is not None]
        r192 = [v["192"] for v in rat.values() if v.get("192") is not None]
        lnLmax = [v["lnZ_max_rel"] for k, v in res["linear"].items() if k.startswith(samp + "|")]
        out[samp] = {
            "n_sn": sn["n"], "n_rows": n, "n_par": 3, "dof": dof,
            "chi2_min_joint": c_tot, "chi2_obs_json": res["chi2_obs"][samp], "chi2_joint_minus_dof": c_tot - dof,
            "z_score": (c_tot - dof) / np.sqrt(2.0 * dof), "p_value": float(chi2_dist.sf(c_tot, dof)),
            "split_at_joint_fit": {"rows15": c_15, "sn": c_sn},
            "sn_alone": {"chi2_min": c_sn_min, "Omega_m": om_sn, "dof": sn["n"] - 2, "p_value": float(chi2_dist.sf(c_sn_min, sn["n"] - 2))},
            "tension": {"delta_chi2": d_ten, "dof": 1, "p_value": float(chi2_dist.sf(d_ten, 1)), "sigma_equiv": float(np.sqrt(max(d_ten, 0.0)))},
            "loosening_vs_pantheon": {"5pct": {"n": len(r5), "min": min(r5) if r5 else None, "median": float(np.median(r5)) if r5 else None, "max": max(r5) if r5 else None},
                                      "192": {"n": len(r192), "min": min(r192) if r192 else None, "median": float(np.median(r192)) if r192 else None, "max": max(r192) if r192 else None}},
            "lnL_max_rel": {"min": min(lnLmax), "max": max(lnLmax), "n_cells": len(lnLmax), "n_above_1.92": int(sum(v > 1.92 for v in lnLmax))},
        }
    # orderings
    order_ten = sorted(MS.SAMPLES, key=lambda s: out[s]["tension"]["delta_chi2"])
    order_p = sorted(MS.SAMPLES, key=lambda s: -out[s]["p_value"])
    order_5 = sorted(MS.SAMPLES, key=lambda s: out[s]["loosening_vs_pantheon"]["5pct"]["median"] if s != "pantheon" else 1.0)
    order_192 = sorted(MS.SAMPLES, key=lambda s: out[s]["loosening_vs_pantheon"]["192"]["median"] if s != "pantheon" else 1.0)
    out["orderings"] = {"by_tension_delta_chi2": order_ten, "by_joint_p_value_desc": order_p, "by_5pct_loosening_median": order_5,
                        "by_192_loosening_median": order_192,
                        "prediction_5pct_follows_tension": order_5 == order_ten, "prediction_192_follows_tension": order_192 == order_ten,
                        "prediction_5pct_follows_joint_p": order_5 == order_p, "prediction_192_follows_joint_p": order_192 == order_p}
    (HERE / "muestras_sn_lcdm.json").write_text(json.dumps(out, indent=1, default=float), encoding="utf-8")
    print(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
