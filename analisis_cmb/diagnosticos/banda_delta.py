#!/usr/bin/env python3
"""68 % band difference between Delta=0.05 and Delta=0.10 per observable and theta, in sigma_Planck.
Rows with NaN in an observable (CAMB failure or silent NaN) are excluded and counted."""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import planck_datos as P
D = Path(__file__).resolve().parents[1] / "resultados" / (sys.argv[1] if len(sys.argv) > 1 else "sigma_2.506e-02")
tt, lens = P.load_tt(), P.load_lensing()
sig = {"TT": 0.5 * (tt["err_lo"] + tt["err_hi"]), "pp": lens["err"], "thetastar": np.array([P.THETASTAR100_ERR])}
lines = ["obs theta n_valid(D0.05) n_valid(D0.10) max|dp16|/sigma max|dp84|/sigma"]
for th in ("0.001", "1", "20"):
    a = np.load(D / f"observables_theta{th}_Delta0.05.npz")
    b = np.load(D / f"observables_theta{th}_Delta0.10.npz")
    for o in ("TT", "pp", "thetastar"):
        va = a[o][~np.isnan(a[o]).any(1)]
        vb = b[o][~np.isnan(b[o]).any(1)]
        pa, pb = np.percentile(va, [16, 84], axis=0), np.percentile(vb, [16, 84], axis=0)
        d16 = np.max(np.abs(pa[0] - pb[0]) / sig[o]); d84 = np.max(np.abs(pa[1] - pb[1]) / sig[o])
        lines.append(f"{o} {th} {len(va)} {len(vb)} {d16:.3f} {d84:.3f}")
(D / "banda_delta.txt").write_text("\n".join(lines) + "\nTT sigma = mean of asymmetric errors\n")
print("\n".join(lines))
