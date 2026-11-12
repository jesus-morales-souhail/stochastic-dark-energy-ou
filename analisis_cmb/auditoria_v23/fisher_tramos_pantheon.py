#!/usr/bin/env python3
"""
Fisher error of Omega_m with supernovae alone (Pantheon+, the 1580 SN of exacto_sn.load_sn: zHD > 0.01, no calibrators), h fixed,
M marginalised analytically: (i) the 1580 individual SN with the full STAT+SYS covariance; (ii) the same SN compressed by GLS into
the 22 redshift bins of inyeccion.py (scenario D: edges at the geometric midpoints of the Union3 nodes, first edge 0.01, last the
largest zHD; T = (E' W E)^-1 E' W, C_b = (E' W E)^-1). Fiducial: the LCDM best fit p0 = (h, omega_m, M) of BAO + theta_* + omega_m +
Pantheon+ used by inyeccion.py (muestras_sn.make_linear), convention B, X = 0. Derivatives: central finite differences of the exact
distance modulus with respect to omega_m at fixed h (step 1e-5, as make_linear), and exactly 1 with respect to M.
Fisher F = B' W B on (omega_m, M); marginalising M analytically with a flat prior gives sigma^2(omega_m) = 1 / (F_ww - F_wM^2 / F_MM),
which is the (omega_m, omega_m) element of F^-1. sigma(Omega_m) = sigma(omega_m) / h^2 (h fixed). Writes fisher_tramos_pantheon.json.

    ../venv/bin/python -I fisher_tramos_pantheon.py
"""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import inyeccion as I  # noqa: E402  (make_full, make_binned; its main() is not run)
import exacto_sn as S  # noqa: E402

full = I.make_full(); binned = I.make_binned(full)
sn = full["sn"]; h0, wm0, M0 = [float(v) for v in full["P0"]]
assert sn is S.SN["all"] and sn["n"] == 1580
def fisher(B, C):
    W = np.linalg.inv(C); F = B.T @ W @ B
    s2_marg = 1.0 / (F[0, 0] - F[0, 1] ** 2 / F[1, 1])
    s2_inv = float(np.linalg.inv(F)[0, 0])
    return F, s2_marg, s2_inv
B_full = full["Bm"][15:, 1:]  # columns: omega_m, M (h fixed)
F_full, s2m_full, s2i_full = fisher(B_full, full["C"][15:, 15:])
B_bin = binned["Bm"][15:, 1:]
F_bin, s2m_bin, s2i_bin = fisher(B_bin, binned["C"][15:, 15:])
out = {
 "sample": "Pantheon+ (exacto_sn.load_sn 'all': zHD > 0.01, IS_CALIBRATOR == 0; 1580 SN; data files analisis_cmb/sn/datos/Pantheon+SH0ES.dat and Pantheon+SH0ES_STAT+SYS.cov, sha256 in auditoria_v23/datos/FUENTES.json)",
 "n_sn": int(sn["n"]), "covariance": "STAT+SYS 1701x1701 restricted to the 1580 rows (C_SN = inv(sn['Ci']), as muestras_sn.make_linear)",
 "fiducial": {"h": h0, "omega_m": wm0, "M": M0, "Omega_m": wm0 / h0**2, "source": "LCDM best fit of BAO + theta_* + omega_m + Pantheon+ (exacto_sn.fit_combo, M explicit), convention B, X = 0; the point used by inyeccion.py"},
 "fixed": "h fixed at the fiducial; Omega_r h^2 = OMR fixed; X = 0",
 "derivatives": "central finite differences of the exact mu(z) with respect to omega_m at fixed h (step 1e-5; Bm[:, 1] of make_linear); d mu / dM = 1 exactly",
 "fisher_definition": "F = B' W B with B = [d mu/d omega_m, d mu/dM] (n x 2), W = C_SN^-1; M marginalised analytically with a flat prior: sigma^2(omega_m) = 1/(F_ww - F_wM^2/F_MM) = (F^-1)_ww; sigma(Omega_m) = sigma(omega_m)/h^2",
 "bins": {"definition": "inyeccion.make_binned: edges = [0.01, geometric midpoints of the 22 Union3 nodes (datos/union3/mu_mat_union3_cosmo=2_mu.fits, first row), max(zHD, last node) + 1e-6]; bins with no SN dropped; GLS compression T = (E'WE)^-1 E'W with E the membership matrix and W = C_SN^-1; C_b = (E'WE)^-1; the model and derivatives are compressed with the same T",
          "n_bins": int(binned["nb"]), "edges": binned["edges"], "counts": binned["counts"]},
 "individual_1580": {"F": F_full.tolist(), "sigma2_omega_m_marg": s2m_full, "sigma2_omega_m_inv": s2i_full, "sigma_omega_m": float(np.sqrt(s2m_full)), "sigma_Omega_m": float(np.sqrt(s2m_full)) / h0**2},
 "binned_22": {"F": F_bin.tolist(), "sigma2_omega_m_marg": s2m_bin, "sigma2_omega_m_inv": s2i_bin, "sigma_omega_m": float(np.sqrt(s2m_bin)), "sigma_Omega_m": float(np.sqrt(s2m_bin)) / h0**2},
}
out["ratio_binned_over_individual"] = out["binned_22"]["sigma_Omega_m"] / out["individual_1580"]["sigma_Omega_m"]
out["values_quoted_in_INYECCION_md"] = {"individual": 0.01745, "binned": 0.01754, "note": "computed in-session on 2026-11-12 and only written in the note; this json is the first stored calculation"}
(HERE / "fisher_tramos_pantheon.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
print(json.dumps({k: out[k] for k in ("fiducial", "individual_1580", "binned_22", "ratio_binned_over_individual")}, indent=1))
print("bins", out["bins"]["n_bins"], out["bins"]["counts"])
