#!/usr/bin/env python3
"""
OU residual (same kernel as the DESI BAO fit) -> rho_DE(a) -> w(a) -> CAMB (PPF).

Definitions fixed in NOTAS_ANALISIS_CMB.md:
  - X(x), x = ln(1+z) = -ln a: stationary Gaussian process with
    Cov = sigma_X^2 exp(-theta |dx|)  (scripts/desi_dr2_data.py::add_ou_kernel).
  - Reading (a): delta ln rho_DE = X / Omega_L0, normalised at a = 1:
        rho_DE(a) / rho_L0 = 1 + (X_s(a) - X_s(1)) / Omega_L0
    with X_s the Gaussian-smoothed X (std Delta in ln a).
  - The stationary process is padded beyond a = 1 (ln a > 0, purely numerical) so the
    smoothing kernel is symmetric at a = 1. The far edge a = 1e-3 stays truncated.
  - w(a) = -1 - (1/3) d ln rho_DE / d ln a (finite differences on the smoothed curve).
  - Cosmology fixed to Planck 2018 base LCDM for every curve.

Does not import or modify the BAO fitting scripts.
"""

from __future__ import annotations

import numpy as np
import camb
from camb.dark_energy import DarkEnergyPPF

SEED = 20261003
SIGMA_X = 0.025063155506929828  # results/profile_sigma_x/profile_sigma_x.json
A_MIN = 1e-3
PAD = 0.4  # ln a padding beyond a = 1: 4 x Delta_max (Delta_max = 0.10)

PLANCK18 = dict(H0=67.36, ombh2=0.02237, omch2=0.1200, tau=0.0544, ns=0.9649, As=2.1e-9, mnu=0.06)
LMAX = 2500


def lna_grid(n: int) -> np.ndarray:
    """Uniform grid in ln a from ln(1e-3) to 0 (increasing a)."""
    return np.linspace(np.log(A_MIN), 0.0, n)


def padded_grid(n: int, pad: float = PAD) -> tuple[np.ndarray, np.ndarray]:
    """Main n-point grid plus points with the same spacing up to ln a = pad. Returns (lna, main_mask)."""
    main = lna_grid(n)
    h = main[1] - main[0]
    extra = np.arange(1, int(np.ceil(pad / h)) + 1) * h
    lna = np.concatenate([main, extra])
    mask = np.zeros(len(lna), bool)
    mask[:n] = True
    return lna, mask


def ou_paths(x: np.ndarray, theta: float, sigma_x: float, z: np.ndarray) -> np.ndarray:
    """
    Exact draws of the stationary GP with Cov = sigma^2 exp(-theta|dx|) at points x
    (any order), from standard normals z of shape (n_real, len(x)).

    The process is Markov, so the Cholesky factor of this covariance (points sorted) is
    exactly this AR(1) recursion; works on irregular grids and avoids the ill-conditioned
    Cholesky at theta -> 0. X is linear in sigma_x for fixed z (common random numbers).
    """
    z = np.atleast_2d(z)
    order = np.argsort(x)
    xs = x[order]
    zs = z[:, order]
    out = np.empty_like(zs)
    out[:, 0] = sigma_x * zs[:, 0]
    rho = np.exp(-theta * np.diff(xs))
    amp = sigma_x * np.sqrt(1.0 - rho * rho)
    for k in range(1, len(xs)):
        out[:, k] = rho[k - 1] * out[:, k - 1] + amp[k - 1] * zs[:, k]
    res = np.empty_like(out)
    res[:, order] = out
    return res


def smoothing_matrix(lna: np.ndarray, delta: float) -> np.ndarray:
    """Gaussian kernel (std = delta in ln a), renormalised over the grid at the edges."""
    d = lna[:, None] - lna[None, :]
    w = np.exp(-0.5 * (d / delta) ** 2)
    return w / w.sum(axis=1, keepdims=True)


def rho_and_w(lna: np.ndarray, mask: np.ndarray, x_s: np.ndarray, omega_l0: float):
    """
    rho_DE/rho_L0 and w on the main grid (a <= 1) from smoothed X on the padded grid.
    x_s: (n_real, len(lna)). The derivative is taken on the padded grid, so a = 1 is interior.
    """
    x_s = np.atleast_2d(x_s)
    i1 = np.flatnonzero(mask)[-1]  # a = 1
    rho = 1.0 + (x_s - x_s[:, [i1]]) / omega_l0
    if np.any(rho <= 0):
        raise ValueError("rho_DE <= 0 somewhere: reading (a) breaks for this draw")
    w = -1.0 - np.gradient(np.log(rho), lna, axis=1) / 3.0
    return rho[:, mask], w[:, mask]


def make_params(a: np.ndarray | None = None, w: np.ndarray | None = None,
                w0: float | None = None, wa: float | None = None,
                H0: float | None = None, halofit: str | None = None) -> camb.CAMBparams:
    """Planck 2018 base LCDM parameters (H0 overridable), flat, DarkEnergyPPF.
    halofit: None = CAMB default (mead2020); 'takahashi'; 'lineal' = no non-linear correction."""
    p = PLANCK18
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=p["H0"] if H0 is None else H0, ombh2=p["ombh2"], omch2=p["omch2"],
                       tau=p["tau"], mnu=p["mnu"], omk=0.0)
    if halofit is not None and halofit != "lineal":
        pars.NonLinearModel.set_params(halofit_version=halofit)
    pars.InitPower.set_params(As=p["As"], ns=p["ns"])
    pars.set_for_lmax(LMAX, lens_potential_accuracy=1)
    de = DarkEnergyPPF()
    if a is not None:
        de.set_w_a_table(a, w)
    elif w0 is not None:
        de.set_params(w=w0, wa=wa)
    pars.DarkEnergy = de
    if halofit == "lineal":
        pars.NonLinear = camb.model.NonLinear_none
    return pars


def run(pars: camb.CAMBparams, keep_results: bool = True) -> dict:
    res = camb.get_results(pars)
    tt = res.get_cmb_power_spectra(pars, CMB_unit="muK", spectra=["total"])["total"][:, 0]
    pp = res.get_lens_potential_cls(lmax=LMAX)[:, 0]  # [L(L+1)]^2 C_L^phiphi / 2pi
    der = res.get_derived_params()
    out = {
        "DlTT": tt,
        "Clpp": pp,
        "thetastar100": der["thetastar"],
        "cosmomc_theta100": 100.0 * res.cosmomc_theta(),
    }
    if keep_results:
        out["results"] = res
    return out


def omega_l0() -> float:
    """Omega_DE today in the Planck 2018 LCDM baseline (mnu = 0.06 included)."""
    return camb.get_background(make_params()).omega_de


# ---------------------------------------------------------------------------
# Production scheme (NOTAS: "Convergencia de malla"): the path lives on a fine uniform
# grid with step h_tab / K_PATH; CAMB gets a 2000-point table subsampled from it.
# ---------------------------------------------------------------------------

N_TAB = 2000
K_PATH = 16


def production_grids(n_tab: int = N_TAB, k: int = K_PATH, pad: float = PAD):
    """Fine path grid (ln a, step h_tab/k) incl. padding, table indices into it, table main mask."""
    l0 = np.log(A_MIN)
    h_tab = -l0 / (n_tab - 1)
    n_pad_tab = int(np.ceil(pad / h_tab))
    n_fine = (n_tab - 1 + n_pad_tab) * k + 1
    lna_fine = l0 + (h_tab / k) * np.arange(n_fine)
    tab_idx = np.arange(0, n_fine, k)
    main = np.arange(len(tab_idx)) < n_tab
    return lna_fine, tab_idx, main


def ou_uniform(z: np.ndarray, h: float, theta: float, sigma_x: float) -> np.ndarray:
    """Exact stationary OU on a uniform grid (rows = realisations); linear in sigma_x for fixed z."""
    from scipy.signal import lfilter

    z = np.atleast_2d(z)
    r = np.exp(-theta * h)
    x0 = sigma_x * z[:, 0]
    y, _ = lfilter([sigma_x * np.sqrt(1.0 - r * r)], [1.0, -r], z[:, 1:], axis=1, zi=(r * x0)[:, None])
    return np.concatenate([x0[:, None], y], axis=1)


def smooth_uniform(x: np.ndarray, h: float, delta: float) -> np.ndarray:
    """Gaussian smoothing (std delta in ln a) along axis 1, renormalised at the edges."""
    from scipy.ndimage import gaussian_filter1d

    s = delta / h
    num = gaussian_filter1d(x, s, axis=-1, mode="constant", truncate=6)
    den = gaussian_filter1d(np.ones(x.shape[-1]), s, mode="constant", truncate=6)
    return num / den


def rho_w_table(lna_tab: np.ndarray, main: np.ndarray, xs_tab: np.ndarray, omega_l0: float):
    """
    Like rho_and_w but does not raise: returns (rho, w, bad) with bad[r] = True when
    rho_DE <= 0 anywhere for realisation r (reading (a) breaks); w is NaN for those rows.
    """
    xs_tab = np.atleast_2d(xs_tab)
    i1 = np.flatnonzero(main)[-1]
    rho = 1.0 + (xs_tab - xs_tab[:, [i1]]) / omega_l0
    bad = np.any(rho <= 0, axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        w = -1.0 - np.gradient(np.log(np.where(rho > 0, rho, np.nan)), lna_tab, axis=1) / 3.0
    w[bad] = np.nan
    return rho[:, main], w[:, main], bad


def thetastar_background(pars: camb.CAMBparams) -> float:
    """100 theta_* from the background + thermal history only (no perturbations)."""
    return camb.get_background(pars).get_derived_params()["thetastar"]


def solve_H0_for_thetastar(target100: float, make, lo: float = 40.0, hi: float = 100.0,
                           tol: float = 1e-5, max_iter: int = 60) -> tuple[float, float, int]:
    """
    Bisection on H0 (omega_b, omega_c, mnu fixed; flat, Omega_DE = 1 - Omega_m - Omega_r)
    so that 100 theta_*(H0) = target100. make(H0) -> CAMBparams. 100 theta_* increases with H0.
    Returns (H0, achieved 100 theta_*, iterations); raises ValueError if not bracketed.
    """
    f_lo = thetastar_background(make(lo)) - target100
    f_hi = thetastar_background(make(hi)) - target100
    if f_lo * f_hi > 0:
        raise ValueError(f"theta_* not bracketed in H0=[{lo},{hi}]: f={f_lo:.3g},{f_hi:.3g}")
    for it in range(1, max_iter + 1):
        mid = 0.5 * (lo + hi)
        f_mid = thetastar_background(make(mid)) - target100
        if f_lo * f_mid <= 0:
            hi = mid
        else:
            lo, f_lo = mid, f_mid
        if hi - lo < tol:
            break
    h0 = 0.5 * (lo + hi)
    return h0, thetastar_background(make(h0)), it
