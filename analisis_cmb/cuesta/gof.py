#!/usr/bin/env python3
"""
Exploration only (GOF.md): a bound on sigma_X from the goodness of fit of LCDM, with the set-up of cuesta.py.

For theta in {0.1, 0.3, 1, 3, 10} and sigma_X in geomspace(3e-3, 0.15, 16), N_UNI universes per cell exactly as in cuesta.py
(13 DESI DR2 numbers, DESI best-fit LCDM background, OU residual S_i X(z_i) with S by type, official 13x13 noise with the same
N_UNI draws in every cell). Each universe is fitted with flat LCDM only (Omega_m, h r_d; cuesta.fit_lcdm) and chi^2_min is kept.
Goodness-of-fit bound: the sigma_X at which P(chi^2_min <= CHI2_REAL | sigma_X, theta) = 0.05, interpolated linearly in
ln sigma_X; CHI2_REAL = 10.27 is the chi^2 of LCDM on the real 13 numbers (cuesta.json). The mean chi^2_min per cell is kept too.
Alongside: the 95% bounds with the background fixed (a), profiled (b) and marginalised (c) at the same five theta, with the
perfil_13.py machinery on the same 13 numbers (TRECE.md gives them at theta = 0.1, 1, 10).

Writes analisis_cmb/cuesta/gof.{json,txt} and gof_chi2.npz; nothing in results/, figures/ or papers.
"""

from __future__ import annotations

import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cuesta as C  # noqa: E402  (data, reference background, universes(), fit_lcdm(), bounds machinery)

THETAS = (0.1, 0.3, 1.0, 3.0, 10.0)
SIGMAS = np.geomspace(3e-3, 0.15, 16)
N_UNI = C.N_UNI
SEED_OU = 9_000_000
P_LEVEL = 0.05


def chi2_lcdm(y: np.ndarray) -> float:
    return C.fit_lcdm(y)[0]


def bounds_three() -> dict:
    pred = C.PRED_REF
    D = C.P13.Data(C.Z, C.KINDS, C.OBS / pred, C.COV / np.outer(pred, pred), C.S13)
    out = {}
    for th in THETAS:
        out[f"theta_{th:g}"] = {"a_fixed": C.PF.bound(lambda s, t=th: D.ll_a(s, t))["sigma_95"],
                                "b_profiled": C.PF.bound(lambda s, t=th: D.ll_b(s, t))["sigma_95"],
                                "c_marginalised": C.PF.bound(lambda s, t=th: D.ll_c(s, t))["sigma_95"]}
    return out


def gof_bound(p: np.ndarray) -> float | None:
    """sigma_X where P(chi2 <= CHI2_REAL) crosses P_LEVEL, linear in ln sigma (P is decreasing in sigma)."""
    ls = np.log(SIGMAS)
    below = np.nonzero(p <= P_LEVEL)[0]
    if len(below) == 0:
        return None
    k = int(below[0])
    if k == 0:
        return float(SIGMAS[0])
    f = (p[k - 1] - P_LEVEL) / (p[k - 1] - p[k])
    return float(np.exp(ls[k - 1] + f * (ls[k] - ls[k - 1])))


def main() -> None:
    t0 = time.time()
    chi2_real = C.fit_both(C.OBS)[1]
    print(f"chi2 real de LCDM: {chi2_real:.3f}", flush=True)
    bnd = bounds_three()
    print("cotas (13 numeros):", json.dumps(bnd), flush=True)
    noise = np.random.default_rng(C.SEED_NOISE).standard_normal((N_UNI, 13)) @ C.L_COV.T
    chi2 = np.zeros((len(THETAS), len(SIGMAS), N_UNI))
    with Pool() as pool:
        ctrl = np.array(pool.map(chi2_lcdm, list(C.universes(0.0, 0.0, noise, 0)), chunksize=50))
        for i, th in enumerate(THETAS):
            for j, s in enumerate(SIGMAS):
                ys = C.universes(float(s), th, noise, SEED_OU + 100 * i + j)
                chi2[i, j] = pool.map(chi2_lcdm, list(ys), chunksize=50)
            print(f"theta = {th:g} hecho [{time.time() - t0:.0f} s]", flush=True)
    p = np.mean(chi2 <= chi2_real, axis=2)
    mean = chi2.mean(axis=2)
    p_ctrl, mean_ctrl = float(np.mean(ctrl <= chi2_real)), float(ctrl.mean())
    cells = {}
    for i, th in enumerate(THETAS):
        g = gof_bound(p[i])
        cells[f"theta_{th:g}"] = {"gof_bound": g, "P_chi2_le_real": p[i].tolist(), "mean_chi2": mean[i].tolist(),
                                  "mean_chi2_at_gof_bound": float(np.interp(np.log(g), np.log(SIGMAS), mean[i])) if g else None,
                                  **bnd[f"theta_{th:g}"]}
    np.savez_compressed(HERE / "gof_chi2.npz", thetas=np.array(THETAS), sigmas=SIGMAS, chi2=chi2.astype(np.float32), control=ctrl)
    out = {"chi2_real_lcdm": chi2_real, "p_level": P_LEVEL, "n_universes": N_UNI, "thetas": THETAS, "sigmas": SIGMAS.tolist(),
           "control": {"P_chi2_le_real": p_ctrl, "mean_chi2": mean_ctrl}, "seeds": {"noise": C.SEED_NOISE, "ou_base": SEED_OU},
           "cells": cells, "runtime_s": time.time() - t0}
    (HERE / "gof.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    fmt = lambda v: f"{v:.4g}" if v else "sin cota"
    lines = [f"chi2 real de LCDM con los 13 numeros: {chi2_real:.2f}; control sigma_X = 0: P(chi2 <= real) = {p_ctrl:.3f}, "
             f"<chi2> = {mean_ctrl:.2f}; cota: P(chi2_min <= {chi2_real:.2f}) = {P_LEVEL}", "",
             f"{'theta':>6s} {'gof':>9s} {'fijo':>9s} {'perfilado':>9s} {'marginal':>9s} {'<chi2> en gof':>13s}"]
    for th in THETAS:
        c = cells[f"theta_{th:g}"]
        lines.append(f"{th:6g} {fmt(c['gof_bound']):>9s} {fmt(c['a_fixed']):>9s} {fmt(c['b_profiled']):>9s} "
                     f"{fmt(c['c_marginalised']):>9s} {c['mean_chi2_at_gof_bound']:13.2f}")
    lines += ["", "P(chi2_min <= real) y <chi2> por celda (filas theta, columnas sigma_X):",
              "sigma_X: " + " ".join(f"{s:7.4f}" for s in SIGMAS)]
    for i, th in enumerate(THETAS):
        lines.append(f"P  {th:4g}: " + " ".join(f"{v:7.3f}" for v in p[i]))
        lines.append(f"<chi2>{th:4g}: " + " ".join(f"{v:7.2f}" for v in mean[i]))
    (HERE / "gof.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
