#!/usr/bin/env python3
"""
Deterministic oscillator: X(ln a) = A sin(2 pi ln a / P + phi), same pipeline as the OU analysis
(reading (a), zero fixed at a = 1, H0 re-solved for 100 theta_* = LCDM, halofit takahashi,
same Planck files), 60 curves: A in {0.005, 0.01, 0.025}, P in {0.1, 0.2, 0.5, 1, 2}, phi in {0, pi/2, pi, 3pi/2}.

    rho_DE(a)/rho_L0 = 1 + (X(a) - X(1)) / Omega_L0
    w(a) = -1 - (1/3) d ln rho_DE / d ln a   (analytic derivative; X is smooth, no Delta smoothing)

Per curve: Delta chi^2 (TT 2<=l<30, phiphi L 8-400) vs LCDM, H0, and Delta chi^2 vs LCDM against
DESI DR2 Gaussian BAO (13 distances D/r_d, full 13x13 covariance) loaded with the repo's
scripts/desi_dr2_data.load_gaussian_bao_vector (unmodified). The repo's distance functions
(joint_w0wa_sigma_desi.DM/DH/DV, eos_efectiva.DM_cpl/...) only accept CPL (w0, wa) with
H0 = 67.4, Om = 0.315 and no r_d, so distances here come from CAMB (D_M, c/H, D_V, rdrag)
for the same cosmology as the CMB part (H0 re-solved).

Outputs: resultados/{oscilador.json, oscilador.csv, tabla.md}
"""

from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "2")

import json
import multiprocessing as mp
import sys
from itertools import product
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
CMB = HERE.parent
REPO = CMB.parent
sys.path.insert(0, str(CMB))
OUT = HERE / "resultados"

AS = (0.005, 0.01, 0.025)
PS = (0.1, 0.2, 0.5, 1.0, 2.0)
PHIS = (0.0, 0.5 * np.pi, np.pi, 1.5 * np.pi)
PHI_LABEL = {0.0: "0", 0.5 * np.pi: "π/2", np.pi: "π", 1.5 * np.pi: "3π/2"}
N_TAB = 2000
C_KMS = 299792.458
N_WORKERS = 8

_CTX = {}


def _init(target, om_l, bao):
    import planck_datos as P

    _CTX.update(target=target, om_l=om_l, bao=bao, lens=P.load_lensing())


def w_table(A, P, phi, om_l):
    lna = np.linspace(np.log(1e-3), 0.0, N_TAB)
    x = A * np.sin(2 * np.pi * lna / P + phi)
    dx = A * (2 * np.pi / P) * np.cos(2 * np.pi * lna / P + phi)
    rho = 1.0 + (x - A * np.sin(phi)) / om_l
    if np.any(rho <= 0):
        return None, None
    return np.exp(lna), -1.0 - (dx / om_l) / rho / 3.0


def bao_vector(res, bao):
    """Theory for the 13 DESI DR2 Gaussian BAO quantities (D/r_d), from CAMB results."""
    rd = res.get_derived_params()["rdrag"]
    out = []
    for z, lab in zip(bao["z_rows"], bao["labels"]):
        dm = res.comoving_radial_distance(z)
        dh = C_KMS / res.hubble_parameter(z)
        q = lab[:2]
        out.append({"DM": dm, "DH": dh, "DV": (z * dm * dm * dh) ** (1 / 3)}[q] / rd)
    return np.array(out)


def _curve(task):
    import ou_cmb as M
    import planck_datos as P

    kind = task[0]
    if kind == "cpl":
        def make(h):
            return M.make_params(w0=task[1], wa=task[2], H0=h, halofit="takahashi")
    else:
        a, w = w_table(*task[1:4], _CTX["om_l"])
        if a is None:
            return {"fail": "rho_DE <= 0"}

        def make(h):
            return M.make_params(a=a, w=w, H0=h, halofit="takahashi")
    out = {"fail": None}
    try:
        out["H0"], out["thetastar_solved"], _ = M.solve_H0_for_thetastar(_CTX["target"], make)
        pars = make(out["H0"])
        r = M.run(pars, keep_results=True)
        out["TT"] = r["DlTT"][2:30]
        out["pp"] = P.bin_pp(r["Clpp"], _CTX["lens"])
        out["bao"] = bao_vector(r["results"], _CTX["bao"])
        del r["results"]
        for k in ("TT", "pp", "bao"):
            if not np.all(np.isfinite(out[k])):
                out["fail"] = f"silent NaN in {k}"
    except Exception as exc:  # noqa: BLE001
        out["fail"] = f"{type(exc).__name__}: {exc}".replace("\n", " ")
    return out


def main():
    import ou_cmb as M
    import estadisticos as E

    sys.path.insert(0, str(REPO / "scripts"))
    from desi_dr2_data import load_gaussian_bao_vector  # repo loader, unmodified

    bao = load_gaussian_bao_vector()
    bao_small = {"z_rows": bao["z_rows"], "labels": bao["labels"]}
    cinv = np.linalg.inv(bao["cov"])
    om_l = M.omega_l0()
    target = M.thetastar_background(M.make_params(halofit="takahashi"))
    tg = E.targets()

    grid = list(product(AS, PS, PHIS))
    tasks = [("cpl", -1.0, 0.0), ("cpl", -0.42, -1.75)] + [("osc", A, P, phi) for A, P, phi in grid]
    with mp.get_context("spawn").Pool(N_WORKERS, initializer=_init, initargs=(target, om_l, bao_small)) as pool:
        res = pool.map(_curve, tasks)
    lcdm, desi = res[0], res[1]

    def chi2_bao(v):
        d = v - bao["mean"]
        return float(d @ cinv @ d)

    def dchis(r):
        return {"TT": float(E.chi2(r["TT"], *tg["TT"]) - E.chi2(lcdm["TT"], *tg["TT"])),
                "pp": float(E.chi2(r["pp"], *tg["pp"]) - E.chi2(lcdm["pp"], *tg["pp"])),
                "bao": chi2_bao(r["bao"]) - chi2_bao(lcdm["bao"])}

    rows = []
    for (A, P, phi), r in zip(grid, res[2:]):
        row = {"A": A, "P": P, "phi": phi, "phi_label": PHI_LABEL[phi], "fail": r["fail"]}
        if r["fail"] is None:
            d = dchis(r)
            row.update(H0=r["H0"], dchi2_bao=d["bao"], dchi2_TT=d["TT"], dchi2_pp=d["pp"],
                       permitida=bool(d["bao"] < 4), distinguible_cmb=bool(max(d["TT"], d["pp"]) > 4))
        rows.append(row)

    summary = {
        "bao_source": bao["source"], "bao_cov_source": bao["cov_source"],
        "bao_loader": "scripts/desi_dr2_data.py::load_gaussian_bao_vector (unmodified)",
        "bao_distances": "CAMB comoving_radial_distance, c/hubble_parameter, D_V, rdrag (repo distance functions are CPL-only)",
        "thetastar100_target": target, "Omega_L0": om_l,
        "lcdm": {"H0": lcdm["H0"], "chi2_bao": chi2_bao(lcdm["bao"])},
        "desi_dr2_cmb_check": {"H0": desi["H0"], **dchis(desi)},
        "rows": rows,
        "n_fail": sum(r["fail"] is not None for r in rows),
        "n_permitida": sum(r.get("permitida", False) for r in rows),
        "n_distinguible": sum(r.get("distinguible_cmb", False) for r in rows),
        "n_permitida_y_distinguible": sum(r.get("permitida", False) and r.get("distinguible_cmb", False) for r in rows),
    }
    OUT.mkdir(exist_ok=True)
    (OUT / "oscilador.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    hdr = "A,P,phi,dchi2_bao,permitida,dchi2_TT,dchi2_pp,H0,distinguible_cmb,fail"
    csv = [hdr] + [f"{r['A']},{r['P']},{r['phi_label']},{r.get('dchi2_bao', '')},{r.get('permitida', '')},"
                   f"{r.get('dchi2_TT', '')},{r.get('dchi2_pp', '')},{r.get('H0', '')},{r.get('distinguible_cmb', '')},"
                   f"{r['fail'] or ''}" for r in rows]
    (OUT / "oscilador.csv").write_text("\n".join(csv) + "\n", encoding="utf-8")

    L = ["# Oscilador X = A sin(2π ln a / P + φ) — Δχ² frente a ΛCDM", "",
         f"Fuente: `resultados/oscilador.json`, `resultados/oscilador.csv`. ΛCDM: H0 = {lcdm['H0']:.3f}, "
         f"χ²_BAO = {summary['lcdm']['chi2_bao']:.2f} (13 puntos). Comprobación DESI DR2+CMB (w0=−0.42, wa=−1.75): "
         f"H0 = {desi['H0']:.2f}, Δχ²_BAO = {summary['desi_dr2_cmb_check']['bao']:.2f}.", "",
         "| A | P | φ | Δχ²_BAO | permitida | Δχ²_TT | Δχ²_φφ | H0 | distinguible en CMB |",
         "|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        if r["fail"]:
            L.append(f"| {r['A']} | {r['P']} | {r['phi_label']} | fallo: {r['fail']} | | | | | |")
            continue
        L.append(f"| {r['A']} | {r['P']} | {r['phi_label']} | {r['dchi2_bao']:.3g} | {'sí' if r['permitida'] else 'no'} | "
                 f"{r['dchi2_TT']:.3g} | {r['dchi2_pp']:.3g} | {r['H0']:.2f} | {'sí' if r['distinguible_cmb'] else 'no'} |")
    L += ["", f"Fallos: {summary['n_fail']}. Permitidas: {summary['n_permitida']}/60. Distinguibles en CMB: "
              f"{summary['n_distinguible']}/60. Permitidas y distinguibles: {summary['n_permitida_y_distinguible']}."]
    (OUT / "tabla.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print("\n".join(L))


if __name__ == "__main__":
    main()
