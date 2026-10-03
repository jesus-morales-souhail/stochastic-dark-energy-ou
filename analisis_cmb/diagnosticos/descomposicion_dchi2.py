#!/usr/bin/env python3
"""
Split Delta chi^2_r = chi^2(OU_r) - chi^2(LCDM) into
    linear    = 2 sum_b r_b delta_b / sigma_b^2
    quadratic = sum_b delta_b^2 / sigma_b^2
with r = LCDM - data, delta = OU_r - LCDM. TT: sigma side taken from the sign of the OU residual
(as in estadisticos.py); the split is exact only where that side equals LCDM's side, so the
remainder (Delta chi^2 - linear - quadratic) is reported.
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
import estadisticos as E

D = HERE / "resultados" / "takahashi" / "sigma_2.506e-02"
TAG = "theta1_Delta0.05"
d = np.load(D / f"observables_{TAG}.npz")
lcdm = json.loads((D / "meta.json").read_text())["references"]["LCDM"]
tg = E.targets()
out = {"source": f"{D.relative_to(HERE)}/observables_{TAG}.npz", "obs": {}}
lines = [f"Delta chi2 split, {TAG}, sigma_X=2.506e-2 (source {out['source']})"]
for o in ("TT", "pp"):
    data, elo, ehi = tg[o]
    L = np.asarray(lcdm[o])
    m = d[o][np.all(np.isfinite(d[o]), axis=1)]
    r = L - data
    delta = m - L
    s = np.where((m - data) > 0, ehi, elo)          # side used for the OU chi2
    lin = np.sum(2 * r * delta / s**2, axis=1)
    quad = np.sum(delta**2 / s**2, axis=1)
    dchi = E.chi2(m, data, elo, ehi) - E.chi2(L, data, elo, ehi)
    rem = dchi - lin - quad
    flip = np.mean(np.sign(m - data) != np.sign(L - data))
    e = {"p95_linear": float(np.percentile(lin, 95)), "p95_quadratic": float(np.percentile(quad, 95)),
         "p50_linear": float(np.percentile(lin, 50)), "p50_quadratic": float(np.percentile(quad, 50)),
         "p95_dchi2": float(np.percentile(dchi, 95)), "max_abs_remainder": float(np.max(np.abs(rem))),
         "fraction_bins_side_flip": float(flip)}
    out["obs"][o] = e
    lines.append(f"{o}: p95 linear={e['p95_linear']:.4g} p95 quadratic={e['p95_quadratic']:.4g} "
                 f"p95 dchi2={e['p95_dchi2']:.4g} | p50 linear={e['p50_linear']:.3g} p50 quad={e['p50_quadratic']:.3g} "
                 f"| max|remainder|={e['max_abs_remainder']:.2g} side flips={flip:.3f}")
(D / "descomposicion_dchi2.json").write_text(json.dumps(out, indent=2))
(D / "descomposicion_dchi2.txt").write_text("\n".join(lines) + "\n")
print("\n".join(lines))
