#!/usr/bin/env python3
"""
Exploration only (MUESTRAS_SN.md): do the per-family bounds of FAMILIAS.md survive a change of supernova sample? Nothing in
results/ or papers changes.

Set-up of familias.py unchanged (convention B, exact history on the exacto grid, BAO + theta_* + omega_m + SN, discard if H^2 <= 0
at the LCDM best fit of the 15 rows, P(chi^2_min <= chi^2_obs) = 5 % interpolated in ln sigma_X, same seeds per cell), changing
only the SN sample. Samples, each read as numbers only from its own folder datos/<sample>/ (provenance in FUENTE.json):
- pantheon: Pantheon+ as analisis_cmb/sn (exacto_sn.load_sn; reference).
- des5yr:   DES-SN5YR Hubble diagram (tag 1.3): MU, MUERR_FINAL, zHD, zHEL; C = diag(MUERR_FINAL^2) + STAT+SYS.txt.gz (systematics).
- dovekie:  DES-Dovekie (main): MU, zHD, zHEL; C = inverse of the released inverse covariance (STAT+SYS.npz, upper triangle).
- union3:   Union3/UNITY1.5 compressed distances: 22 nodes (z, mu) and their inverse covariance; zHEL = z.
The absolute offset (M) is a fit parameter in fit_combo (equivalent to the analytic marginalisation) and marginalised analytically
in the linear marginal (template) and in the SN-only LCDM fit (chi2_sn).
Cells: OU theta in {0.3, 1, 3, 10} (l = 1/theta), Matern 3/2 and Gaussian l in {0.3, 1}; sigma_X in SIGMAS; 2000 paths per cell,
the SAME paths for the four samples (seed per cell as in familias.py / exacto_sn.py). Per sample and cell: 5 % bound with bootstrap,
linear marginalised bounds (-1.92, -1.35) with the response A rebuilt for that sample (marginal_guiada.build_A logic), global
maximum of ln L and its sigma_X, Delta rho(0.1) at the bounds (sigma_b sqrt(2 (1 - K(ln 1.1)))).
Checks: (a) pantheon reproduces the FAMILIAS / EXACTO_SN cells; (b) flat LCDM with each sample alone: Omega_m against the
published value (FUENTE.json), stop if it differs by more than 0.01 (`--check` runs only this).

    ../venv/bin/python -I muestras_sn.py --check
    ../venv/bin/python -I muestras_sn.py            # writes muestras_sn.{json,txt}
"""

from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import gzip
import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import exacto as X  # noqa: E402
import exacto_sn as S  # noqa: E402
import familias as F  # noqa: E402
import marginal_guiada as MG  # noqa: E402
import precision_sn as PS  # noqa: E402

G0 = S.G0
OMR, ZSTAR, C_KMS, MAG = S.OMR, S.ZSTAR, S.C_KMS, S.MAG
SAMPLES = ("pantheon", "des5yr", "dovekie", "union3")
CELLS = [("OU", 1.0 / 0.3), ("OU", 1.0), ("OU", 1.0 / 3.0), ("OU", 0.1), ("Matern32", 0.3), ("Matern32", 1.0), ("Gauss", 0.3), ("Gauss", 1.0)]
THETA_OF = {("OU", 1.0 / 0.3): 0.3, ("OU", 1.0): 1.0, ("OU", 1.0 / 3.0): 3.0, ("OU", 0.1): 10.0}
SIGMAS = (0.03, 0.05, 0.07, 0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 1.0)
N_REAL, N_BOOT, SEED_M = S.N_REAL, 1000, 20261117
U01 = float(np.log1p(0.1))
GRID_LIN = F.GRID_LIN
THRS = F.THRS
DATA = HERE / "datos"
PUB = json.loads((DATA / "FUENTES.json").read_text(encoding="utf-8"))


def label(fam: str, l: float) -> str:
    return f"OU|theta_{THETA_OF[(fam, l)]:g}" if fam == "OU" else f"{fam}|l_{l:g}"


# ----------------------------------------------------------------------------------------------------------------------
# samples
# ----------------------------------------------------------------------------------------------------------------------
def make_sn(zhd, zhel, mu, cov) -> dict:
    zhd, zhel, mu = (np.asarray(v, dtype=float) for v in (zhd, zhel, mu))
    cov = 0.5 * (cov + cov.T)
    L = np.linalg.cholesky(cov)
    Ci = np.linalg.inv(cov)
    one = np.ones(len(mu))
    return {"zhd": zhd, "zhel": zhel, "mb": mu, "L": L, "Ci": Ci, "Ci1": Ci @ one, "e": float(one @ Ci @ one), "n": len(mu), "uhd": np.log1p(zhd)}


def load_des5yr() -> dict:
    d = DATA / "des_sn5yr"
    rows = [l.strip().split(",") for l in open(d / "DES-SN5YR_HD.csv", encoding="ascii") if l.strip()]
    head, body = rows[0], rows[1:]
    col = {k: i for i, k in enumerate(head)}
    get = lambda k: np.array([float(r[col[k]]) for r in body])
    with gzip.open(d / "STAT+SYS.txt.gz", "rt") as f:
        n = int(f.readline())
        sys_ = np.loadtxt(f).reshape(n, n)
    assert n == len(body) == 1829
    cov = sys_ + np.diag(get("MUERR_FINAL") ** 2)
    return make_sn(get("zHD"), get("zHEL"), get("MU"), cov)


def load_dovekie() -> dict:
    d = DATA / "dovekie"
    names, body = None, []
    for l in open(d / "DES-Dovekie_HD.csv", encoding="ascii"):
        p = l.split()
        if not p:
            continue
        if p[0] == "VARNAMES:":
            names = p[1:]
        elif p[0] == "SN:":
            body.append(p[1:])
    col = {k: i for i, k in enumerate(names)}
    get = lambda k: np.array([float(r[col[k]]) for r in body])
    z = np.load(d / "STAT+SYS.npz", allow_pickle=False)
    n = int(np.asarray(z["nsn"]).ravel()[0])
    assert n == len(body) == 1820
    inv = np.zeros((n, n))
    inv[np.triu_indices(n)] = z["cov"].astype(float)
    il = np.tril_indices(n, -1)
    inv[il] = inv.T[il]
    cov = np.linalg.inv(inv)
    return make_sn(get("zHD"), get("zHEL"), get("MU"), cov)


def load_union3() -> dict:
    b = (DATA / "union3" / "mu_mat_union3_cosmo=2_mu.fits").read_bytes()
    hdr = b[:2880].decode("ascii")
    assert "BITPIX  =                  -64" in hdr and "NAXIS1  =                   23" in hdr and "NAXIS2  =                   23" in hdr
    end = b.find(b"END     ")
    start = (end // 2880 + 1) * 2880
    a = np.frombuffer(b[start:start + 23 * 23 * 8], dtype=">f8").reshape(23, 23).astype(float)
    z, mu, icov = a[0, 1:], a[1:, 0], a[1:, 1:]
    return make_sn(z, z, mu, np.linalg.inv(0.5 * (icov + icov.T)))


def load_samples() -> dict:
    return {"pantheon": S.SN["all"], "des5yr": load_des5yr(), "dovekie": load_dovekie(), "union3": load_union3()}


# ----------------------------------------------------------------------------------------------------------------------
# linear machinery per sample (marginal_guiada logic with the sample as argument)
# ----------------------------------------------------------------------------------------------------------------------
def exact_model(sn, h, wm, M, dX):
    h, wm, M = np.atleast_1d(h), np.atleast_1d(wm), np.atleast_1d(M)
    om, orad = (wm / h**2)[:, None], (OMR / h**2)[:, None]
    e2 = om * G0.p3 + orad * G0.p4 + (1.0 - om - orad) + dX
    bad = e2.min(axis=1) <= 0.0
    e2 = np.where(e2 > 0, e2, 1.0)
    inv = 1.0 / np.sqrt(e2)
    f = G0.p1 * inv
    dm_u = np.concatenate([np.zeros((len(h), 1)), np.cumsum(0.5 * (f[:, :-1] + f[:, 1:]) * G0.DUK, axis=1)], axis=1)
    dm, dh = dm_u[:, G0.IDX], inv[:, G0.IDX]
    d = np.where(X.IS_DM, dm, np.where(X.IS_DH, dh, np.cbrt(X.Z13 * dm * dm * dh)))
    bao = (C_KMS / (100.0 * h * S.E.RD_SPL(wm)))[:, None] * d
    dm_star = C_KMS / (100.0 * h) * (dm_u[:, -1] + inv[:, -1] * (np.array([S.E.zstar_of(w) for w in wm]) - ZSTAR))
    th = 100.0 * S.E.RS_SPL(wm) / dm_star / S.E.CAL
    k = sn["_k"]
    t = sn["_t"]
    mu = 5.0 * np.log10((1.0 + sn["zhel"])[None, :] * ((1 - t)[None, :] * dm_u[:, k] + t[None, :] * dm_u[:, k + 1])) + M[:, None]
    out = np.concatenate([bao, th[:, None], wm[:, None], mu], axis=1)
    out[bad] = np.nan
    return out, bad


def make_linear(sn: dict) -> dict:
    sn["_k"] = np.searchsorted(G0.U, sn["uhd"]) - 1
    sn["_t"] = (sn["uhd"] - G0.U[sn["_k"]]) / (G0.U[sn["_k"] + 1] - G0.U[sn["_k"]])
    nrows = 15 + sn["n"]
    D = np.concatenate([X.Y15, sn["mb"]])
    C = np.zeros((nrows, nrows))
    C[:13, :13] = X.COV13
    C[13, 13], C[14, 14] = X.SIG_TH**2, X.SIG_WM**2
    C[15:, 15:] = np.linalg.inv(sn["Ci"])
    c0, h0, wm0, M0 = S.fit_combo(G0, np.zeros(G0.N), sn, True)
    P0 = np.array([h0, wm0, M0])
    B0 = exact_model(sn, P0[0], P0[1], P0[2], np.zeros((1, G0.N)))[0][0]
    Bm = np.zeros((nrows, 3))
    for j, eps in enumerate((1e-4, 1e-5, 1e-3)):
        pp, pm = P0.copy(), P0.copy()
        pp[j] += eps
        pm[j] -= eps
        Bm[:, j] = (exact_model(sn, *pp, np.zeros((1, G0.N)))[0][0] - exact_model(sn, *pm, np.zeros((1, G0.N)))[0][0]) / (2 * eps)
    # full linear response (convention B, Lambda absorbs X(0)), raw units
    h, wm = P0[0], P0[1]
    om, orad = wm / h**2, OMR / h**2
    E2 = om * G0.p3 + orad * G0.p4 + (1 - om - orad)
    Es = np.sqrt(E2)
    N = G0.N
    wint = np.zeros(N)
    wint[1:] += 0.5 * G0.DUK
    wint[:-1] += 0.5 * G0.DUK
    g = G0.p1 / (2 * Es**3)
    f = G0.p1 / Es
    dm = np.concatenate([[0.0], np.cumsum(0.5 * (f[:-1] + f[1:]) * G0.DUK)])
    T = np.tril(np.tile(wint, (N, 1)), -1)
    T[np.arange(1, N), np.arange(1, N)] = 0.5 * G0.DUK
    T[0, 0] = 0.0
    Wdm = -(T * g[None, :]) / np.where(dm > 0, dm, 1.0)[:, None]
    Wdm[:, 0] -= Wdm.sum(axis=1)
    WdH = np.zeros((len(G0.IDX), N))
    for i, k in enumerate(G0.IDX):
        WdH[i, k] = -1.0 / (2 * E2[k])
        WdH[i, 0] += 1.0 / (2 * E2[k])
    rows = []
    for i, (k, kind) in enumerate(zip(G0.IDX, X.KINDS13)):
        rows.append(Wdm[k] if kind == "DM" else (WdH[i] if kind == "DH" else (2 * Wdm[k] + WdH[i]) / 3))
    A = np.zeros((nrows, N))
    A[:13] = np.array(rows) * B0[:13, None]
    A[13] = -Wdm[-1] * B0[13]
    A[15:] = MAG * ((1 - sn["_t"])[:, None] * Wdm[sn["_k"]] + sn["_t"][:, None] * Wdm[sn["_k"] + 1])
    return {"D": D, "C": C, "CF": cho_factor(C, lower=True), "P0": P0, "B0": B0, "Bm": Bm, "A": A, "chi2_lcdm": c0, "A_rows_sum_max_abs": float(np.max(np.abs(A.sum(axis=1))))}


LIN = {}


def _linear(args):
    samp, fam, l = args
    t0 = time.time()
    L = LIN[samp]
    A, Bm = L["A"], L["Bm"]
    K = F.kernel_matrix(fam, G0, l)
    Msig = (A @ K) @ A.T
    Msig = 0.5 * (Msig + Msig.T)
    r = L["D"] - L["B0"]

    def lnz(s):
        return MG.lnZ_lin(cho_factor(L["C"] + s * s * Msig, lower=True), Bm, r)

    z0 = lnz(0.0)
    prof = np.array([lnz(float(s)) - z0 for s in GRID_LIN])
    vals = np.concatenate([[0.0], prof])
    nodes = np.concatenate([[0.0], GRID_LIN])
    k0 = int(np.argmax(vals))
    rec = {"sample": samp, "cell": label(fam, l), "lnZ_rel": prof.tolist(), "sigma_at_max": float(nodes[k0]), "lnZ_max_rel": float(vals[k0]),
           "max_class": "sigma->0" if k0 == 0 else ("upper edge (1.0)" if k0 == len(nodes) - 1 else "interior")}
    for name, drop in THRS.items():
        thr = vals[k0] - drop
        below = np.nonzero((vals < thr) & (np.arange(len(vals)) > k0))[0]
        if len(below) == 0:
            rec[name] = None
        else:
            k = int(below[0])
            lo = nodes[k - 1]
            rec[name] = -1.0 if lo == 0.0 else float(brentq(lambda s: lnz(s) - z0 - thr, lo, nodes[k], xtol=1e-10 * lo))
    rec["seconds"] = time.time() - t0
    return rec


# ----------------------------------------------------------------------------------------------------------------------
# cells: same paths for all samples
# ----------------------------------------------------------------------------------------------------------------------
SN_ALL = {}


def cell_seed(fam: str, sigma: float, l: float):
    if l in F.ELLS and fam in F.FAMILIES:
        return F.cell_seed(fam, sigma, l)
    theta = THETA_OF[(fam, l)]
    if sigma in X.SIGMAS or sigma in (0.01, 0.02, 0.03):
        return S.cell_seed(sigma, theta)
    return [SEED_M, 0, CELLS.index((fam, l)), SIGMAS.index(sigma)]


def _cell(args):
    fam, l, sigma = args
    t0 = time.time()
    rng = np.random.default_rng(cell_seed(fam, sigma, l))
    Xp, _ = F.sample(fam, G0, sigma, l, N_REAL, rng)
    dX = Xp - Xp[:, [0]]
    h0, wm0 = S.P0
    om0, or0 = wm0 / h0**2, OMR / h0**2
    ok = (om0 * G0.p3 + or0 * G0.p4 + (1 - om0 - or0) + dX).min(axis=1) > 0.0
    chi = np.full((N_REAL, len(SAMPLES)), np.nan)
    for i in np.nonzero(ok)[0]:
        for j, samp in enumerate(SAMPLES):
            chi[i, j] = S.fit_combo(G0, dX[i], SN_ALL[samp], True)[0]
    return {"fam": fam, "l": l, "sigma": sigma, "ok": ok, "chi": chi, "seconds": time.time() - t0}


def crossing_boot(P_boot):
    out = np.full(P_boot.shape[0], np.nan)
    for b in range(P_boot.shape[0]):
        c = S.crossing(SIGMAS, P_boot[b])
        out[b] = -1.0 if isinstance(c, str) else (np.nan if c is None else c)
    return out


# ----------------------------------------------------------------------------------------------------------------------
def check_omega(samples: dict) -> dict:
    out = {}
    zero = np.zeros(G0.N)
    for samp, sn in samples.items():
        c, om = S.fit_sn_only(G0, zero, sn)
        pub = PUB[{"pantheon": "pantheon_plus", "des5yr": "des_sn5yr"}.get(samp, samp)]["Omega_m_publicado"]
        out[samp] = {"n_sn": sn["n"], "Omega_m_fit": om, "chi2_rel": c, "Omega_m_published": pub["valor"], "published_error": pub["error"], "source": pub["fuente"],
                     "diff": om - pub["valor"], "within_0.01": abs(om - pub["valor"]) <= 0.01}
    return out


def main(check_only: bool) -> None:
    t0 = time.time()
    samples = load_samples()
    SN_ALL.update(samples)
    chk = check_omega(samples)
    print(json.dumps(chk, indent=1), flush=True)
    if check_only:
        (HERE / "muestras_sn_check.json").write_text(json.dumps(chk, indent=1), encoding="utf-8")
        return
    if not all(v["within_0.01"] for v in chk.values()):
        print("Omega_m check failed; stopping before any bound.", flush=True)
        return
    out = {"samples": SAMPLES, "cells": [label(f, l) for f, l in CELLS], "sigmas": SIGMAS, "n_real": N_REAL, "omega_check": chk, "chi2_obs": {}, "lcdm": {}}
    zero = np.zeros(G0.N)
    for samp, sn in samples.items():
        c, h, wm, M = S.fit_combo(G0, zero, sn, True)
        out["chi2_obs"][samp] = c
        out["lcdm"][samp] = {"H0": 100 * h, "omega_m": wm, "M": M, "n_sn": sn["n"]}
    print("chi2_obs:", out["chi2_obs"], flush=True)
    # cells
    jobs = [(fam, l, s) for fam, l in CELLS for s in SIGMAS]
    jobs.sort(key=lambda j: (j[0] != "Gauss", j[1]))
    res = []
    with Pool() as pool:
        for i, r in enumerate(pool.imap_unordered(_cell, jobs, chunksize=1)):
            res.append(r)
            if i % 16 == 0:
                print(f"  celda {i + 1}/{len(jobs)} ({time.time() - t0:.0f} s) {r['fam']} l={r['l']:.3g} s={r['sigma']} {r['seconds']:.0f} s", flush=True)
    print(f"celdas: {time.time() - t0:.0f} s", flush=True)
    raw = {}
    out["P"] = {}
    for r in res:
        key = f"{label(r['fam'], r['l'])}|sigma_{r['sigma']:g}"
        raw[key.replace("|", "__")] = np.column_stack([r["ok"], r["chi"]])
        ok, chi = r["ok"], r["chi"]
        rec = {"n_ok": int(ok.sum()), "frac_H2_le_0": float(1 - ok.mean())}
        for j, samp in enumerate(SAMPLES):
            p = float(np.mean(chi[ok, j] <= out["chi2_obs"][samp]))
            rec[f"P_{samp}"] = p
        out["P"][key] = rec
    np.savez_compressed(HERE / "muestras_sn_chi2.npz", **raw, columns=np.array(["ok"] + [f"chi2_{s}" for s in SAMPLES]))
    # 5 % bounds with bootstrap
    rng_b = np.random.default_rng([SEED_M, 1])
    out["bounds"] = {}
    for fam, l in CELLS:
        lab = label(fam, l)
        for j, samp in enumerate(SAMPLES):
            P = np.array([out["P"][f"{lab}|sigma_{s:g}"][f"P_{samp}"] for s in SIGMAS])
            b5 = S.crossing(SIGMAS, P)
            Pb = np.empty((N_BOOT, len(SIGMAS)))
            for js, s in enumerate(SIGMAS):
                arr = raw[f"{lab}|sigma_{s:g}".replace("|", "__")]
                okm = arr[:, 0] > 0
                ind = (arr[okm, 1 + j] <= out["chi2_obs"][samp]).astype(float)
                idx = rng_b.integers(0, len(ind), size=(N_BOOT, len(ind)))
                Pb[:, js] = ind[idx].mean(axis=1)
            bb = crossing_boot(Pb)
            fin = np.isfinite(bb) & (bb > 0)
            out["bounds"][f"{samp}|{lab}"] = {"P_vs_sigma": P.tolist(), "exact_5pct": b5, "boot_p16_p84": [float(np.percentile(bb[fin], 16)), float(np.percentile(bb[fin], 84))] if fin.sum() > 10 else None,
                                             "boot_frac_no_bound": float(np.mean(np.isnan(bb)))}
    # linear bounds per sample (A rebuilt)
    for samp, sn in samples.items():
        LIN[samp] = make_linear(sn)
        print(f"  lineal {samp}: A filas suman {LIN[samp]['A_rows_sum_max_abs']:.1e}; chi2 LCDM {LIN[samp]['chi2_lcdm']:.3f}", flush=True)
    ljobs = [(samp, fam, l) for samp in SAMPLES for fam, l in CELLS]
    with Pool(16) as pool:
        lin = pool.map(_linear, ljobs, chunksize=1)
    out["linear"] = {f"{r['sample']}|{r['cell']}": r for r in lin}
    print(f"lineal: {time.time() - t0:.0f} s", flush=True)
    # Delta rho(0.1) at the bounds, ratios, checks
    out["drho01"] = {}
    for samp in SAMPLES:
        for fam, l in CELLS:
            lab = label(fam, l)
            kk = float(F.kfun(fam, U01, l))
            b5 = out["bounds"][f"{samp}|{lab}"]["exact_5pct"]
            l192 = out["linear"][f"{samp}|{lab}"]["192"]
            out["drho01"][f"{samp}|{lab}"] = {"K_u01": kk, "5pct": (b5 * np.sqrt(2 * (1 - kk))) if isinstance(b5, float) and b5 > 0 else None,
                                              "192": (l192 * np.sqrt(2 * (1 - kk))) if isinstance(l192, float) and l192 > 0 else None}
    out["ratios_vs_pantheon"] = {}
    for samp in SAMPLES[1:]:
        for fam, l in CELLS:
            lab = label(fam, l)
            rec = {}
            for rule, get in (("5pct", lambda s: out["bounds"][f"{s}|{lab}"]["exact_5pct"]), ("192", lambda s: out["linear"][f"{s}|{lab}"]["192"]), ("135", lambda s: out["linear"][f"{s}|{lab}"]["135"])):
                a, b = get(samp), get("pantheon")
                rec[rule] = (a / b) if isinstance(a, float) and isinstance(b, float) and a > 0 and b > 0 else None
            out["ratios_vs_pantheon"][f"{samp}|{lab}"] = rec
    out["family_ratios_l1"] = {}
    for samp in SAMPLES:
        rec = {}
        for fam in ("Matern32", "Gauss"):
            for rule, get in (("5pct", lambda s, lb: out["bounds"][f"{s}|{lb}"]["exact_5pct"]), ("192", lambda s, lb: out["linear"][f"{s}|{lb}"]["192"])):
                a, b = get(samp, f"{fam}|l_1"), get(samp, "OU|theta_1")
                rec[f"{fam}/OU|{rule}"] = (a / b) if isinstance(a, float) and isinstance(b, float) and a > 0 and b > 0 else None
        out["family_ratios_l1"][samp] = rec
    out["drho01_spread"] = {}
    for samp in SAMPLES:
        for rule in ("5pct", "192"):
            vals = [v[rule] for k, v in out["drho01"].items() if k.startswith(samp + "|") and v[rule] is not None]
            out["drho01_spread"][f"{samp}|{rule}"] = {"n": len(vals), "min": min(vals) if vals else None, "max": max(vals) if vals else None, "ratio": (max(vals) / min(vals)) if vals else None}
    # (a) reproduction of FAMILIAS / EXACTO_SN with pantheon
    fam_json = json.loads((HERE / "familias.json").read_text())
    rep = {}
    for fam, l in CELLS:
        lab = label(fam, l)
        if l in F.ELLS and fam in F.FAMILIES:
            k = f"{fam}|l_{l:g}"
            rep[lab] = {"familias_5pct": fam_json["bounds"][k]["exact_5pct"], "here_5pct": out["bounds"][f"pantheon|{lab}"]["exact_5pct"],
                        "familias_192": fam_json["linear"][k]["192"], "here_192": out["linear"][f"pantheon|{lab}"]["192"]}
            for r_ in ("5pct", "192"):
                a, b = rep[lab][f"here_{r_}"], rep[lab][f"familias_{r_}"]
                rep[lab][f"rel_diff_{r_}"] = (a / b - 1) if isinstance(a, float) and isinstance(b, float) else ("both none" if a is None and b is None else "mismatch")
    exsn = json.loads((HERE / "exacto_sn.json").read_text())
    for fam, l in CELLS:
        if fam == "OU" and l not in F.ELLS:
            th = THETA_OF[(fam, l)]
            ref = exsn["bounds"][f"theta_{th:g}"]["BAO+CMB+SN"]["P_vs_sigma"]
            ref_s = exsn["sigmas"]
            d = {}
            for s, p in zip(ref_s, ref):
                if s in SIGMAS:
                    d[f"sigma_{s:g}"] = {"exacto_sn": p, "here": out["P"][f"{label(fam, l)}|sigma_{s:g}"]["P_pantheon"]}
            rep[label(fam, l)] = {"P_shared_cells_vs_exacto_sn": d, "max_abs_dP": max(abs(v["exacto_sn"] - v["here"]) for v in d.values())}
    out["reproduction_pantheon"] = rep
    out["shared_sne"] = {"pantheon_used_vs_des5yr_by_CID": 335, "of_which_low_z_by_name": 190, "pantheon_used_vs_dovekie_by_CID": 339, "des5yr_vs_dovekie": 1718,
                         "note": "contados en la preparacion (CID de los ficheros oficiales); Union3 no publica identificadores"}
    out["runtime_s"] = time.time() - t0
    (HERE / "muestras_sn.json").write_text(json.dumps(out, indent=1, default=F._json_default), encoding="utf-8")
    summary(out)


def summary(out: dict) -> None:
    f5 = lambda v: "sin cota" if v is None else (f"< {v[1:]}" if isinstance(v, str) else f"{v:.3f}")
    L = [f"chi2_obs por muestra: {json.dumps({k: round(v, 2) for k, v in out['chi2_obs'].items()})}; n_SN: {json.dumps({k: v['n_sn'] for k, v in out['lcdm'].items()})}",
         "Omega_m SN solas frente al publicado: " + json.dumps({k: (round(v['Omega_m_fit'], 4), v['Omega_m_published'], round(v['diff'], 4)) for k, v in out["omega_check"].items()}), ""]
    L.append("Cotas: celda | muestra | 5 % (boot 16-84) | lineal -1.92 | -1.35 | max ln L (sigma) | cociente/Pantheon+ 5 % | /Pantheon+ -1.92 | Delta rho(0.1) 5 % | -1.92")
    for lab in out["cells"]:
        for samp in out["samples"]:
            b = out["bounds"][f"{samp}|{lab}"]
            li = out["linear"][f"{samp}|{lab}"]
            rt = out["ratios_vs_pantheon"].get(f"{samp}|{lab}", {})
            dr = out["drho01"][f"{samp}|{lab}"]
            bp = b["boot_p16_p84"]
            L.append(f"   {lab:16s} {samp:9s} {f5(b['exact_5pct']):>8s} {('(' + f'{bp[0]:.3f}-{bp[1]:.3f}' + ')') if bp else '':>15s} | {PS.fmt(li['192']):>8s} | {PS.fmt(li['135']):>8s} | {li['lnZ_max_rel']:+.2f} ({li['sigma_at_max']:.3g}, {li['max_class']}) | "
                     f"{(f'{rt.get('5pct'):.2f}' if rt.get('5pct') else '  -  '):>5s} | {(f'{rt.get('192'):.2f}' if rt.get('192') else '  -  '):>5s} | {(f'{dr['5pct']:.3f}' if dr['5pct'] else '  -  '):>6s} | {(f'{dr['192']:.3f}' if dr['192'] else '  -  '):>6s}")
        L.append("")
    L.append("Cocientes entre familias en l = 1: " + json.dumps({s: {k: (round(v, 3) if v else None) for k, v in r.items()} for s, r in out["family_ratios_l1"].items()}))
    L.append("Dispersion de Delta rho(0.1) entre celdas con cota: " + json.dumps({k: {kk: (round(vv, 3) if isinstance(vv, float) else vv) for kk, vv in v.items()} for k, v in out["drho01_spread"].items()}))
    L.append("Reproduccion con Pantheon+: " + json.dumps(out["reproduction_pantheon"], default=F._json_default)[:1500])
    L.append(f"tiempo: {out['runtime_s']:.0f} s")
    (HERE / "muestras_sn.txt").write_text("\n".join(L) + "\n", encoding="utf-8")
    print("\n".join(L))


if __name__ == "__main__":
    main(check_only="--check" in sys.argv)
