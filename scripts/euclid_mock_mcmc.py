#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Euclid-style mock BAO MCMC for vacuum-relaxation protocol.

Samples {w0, wa, theta, sigma_X} under C = C_std + C_OU on a synthetic
Euclid-like BAO vector. Runs three pre-registered scenarios:
  null   (F cup E0): pure measurement noise
  E1     : detectable amplitude, weak mean-reversion
  E2     : detectable amplitude + O(1) mean-reversion

Author: Jesús Morales Souhail
ORCID:  0009-0000-7637-1818
"""

from __future__ import annotations

import os
import json
import numpy as np
from scipy.integrate import cumulative_trapezoid, trapezoid
from scipy.linalg import cholesky, solve_triangular
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    import emcee
except ImportError as e:
    raise SystemExit("emcee is required: pip install emcee") from e

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "euclid_mcmc")
FIG = os.path.join(ROOT, "figures")
os.makedirs(OUT, exist_ok=True)
os.makedirs(FIG, exist_ok=True)

# ---------------------------------------------------------------------------
# Fiducial cosmology (repository convention)
# ---------------------------------------------------------------------------
H0 = 67.4
OM_FID = 0.315
C_KMS = 299792.458
THETA_MIN = 1e-3
THETA_MAX = 10.0
LOG_SX_MIN = np.log(1e-6)
LOG_SX_MAX = np.log(1e-2)

RNG = np.random.default_rng(42)


# ---------------------------------------------------------------------------
# Background: flat CPL
# ---------------------------------------------------------------------------
def w_cpl(z, w0, wa):
    return w0 + wa * z / (1.0 + z)


def E_z(z, Om, w0, wa):
    z_arr = np.atleast_1d(np.asarray(z, dtype=float))
    # rho_DE / rho_DE0 = (1+z)^{3(1+w0+wa)} exp(-3 wa z/(1+z))
    f = (1.0 + z_arr) ** (3.0 * (1.0 + w0 + wa)) * np.exp(
        -3.0 * wa * z_arr / (1.0 + z_arr)
    )
    out = np.sqrt(Om * (1.0 + z_arr) ** 3 + (1.0 - Om) * f)
    if np.ndim(z) == 0:
        return float(out[0])
    return out


def chi_of_z(z, Om, w0, wa, ngrid=800):
    """Comoving distance chi = int_0^z c dz'/H(z') in Mpc (c/H0 units * c/H0)."""
    z = float(z)
    if z <= 0:
        return 0.0
    zs = np.linspace(0.0, z, ngrid)
    invH = 1.0 / (H0 * E_z(zs, Om, w0, wa))
    return float(trapezoid(C_KMS * invH, zs))


def D_V(z, Om, w0, wa):
    """Isotropic BAO distance proxy D_V(z)."""
    z = float(z)
    chi = chi_of_z(z, Om, w0, wa)
    DM = chi  # flat
    DH = C_KMS / (H0 * float(E_z(z, Om, w0, wa)))
    return (z * DH * DM**2) ** (1.0 / 3.0)


def alpha_model(z_arr, w0, wa, Om=OM_FID, w0_fid=-1.0, wa_fid=0.0):
    """alpha = D_V(model) / D_V(fiducial LCDM)."""
    out = np.empty(len(z_arr))
    for i, z in enumerate(z_arr):
        out[i] = D_V(z, Om, w0, wa) / D_V(z, Om, w0_fid, wa_fid)
    return out


def S_of_z(z_arr, Om=OM_FID, eps=1e-3):
    """Numerical S(z) = d ln D_V / d Omega_L along flat direction Om + Ol = 1."""
    S = np.empty(len(z_arr))
    for i, z in enumerate(z_arr):
        Ol0 = 1.0 - Om
        # vary Ol keeping flat: Om' = 1 - (Ol0+eps)
        dv_p = D_V(z, 1.0 - (Ol0 + eps), -1.0, 0.0)
        dv_m = D_V(z, 1.0 - (Ol0 - eps), -1.0, 0.0)
        # d ln DV / d Ol ≈ (DV+ - DV-) / (2 eps DV)
        dv0 = D_V(z, Om, -1.0, 0.0)
        S[i] = (dv_p - dv_m) / (2.0 * eps * dv0)
    return S


# ---------------------------------------------------------------------------
# OU covariance (same structure as repository BAO pipeline)
# ---------------------------------------------------------------------------
def build_cov_OU(z_arr, sigma_arr, S_arr, theta, sigma_X):
    n = len(z_arr)
    x = np.log(1.0 + z_arr)
    C = np.diag(sigma_arr**2)
    s2 = sigma_X**2
    for i in range(n):
        for j in range(n):
            dx = abs(x[i] - x[j])
            C[i, j] += S_arr[i] * S_arr[j] * s2 * np.exp(-theta * dx)
    # numerical jitter
    C += np.eye(n) * 1e-12
    return C


def log_like_gaussian(resid, C):
    try:
        L = cholesky(C, lower=True)
        y = solve_triangular(L, resid, lower=True)
        logdet = 2.0 * np.sum(np.log(np.diag(L)))
        n = len(resid)
        return -0.5 * (np.dot(y, y) + logdet + n * np.log(2.0 * np.pi))
    except np.linalg.LinAlgError:
        return -np.inf


# ---------------------------------------------------------------------------
# Mock Euclid-like data vector
# ---------------------------------------------------------------------------
def make_euclid_grid(n_bins=24, zmin=0.9, zmax=1.8):
    z = np.linspace(zmin, zmax, n_bins)
    # Percent-level BAO errors, slightly better at mid-z (schematic, not official forecast)
    sig = 0.008 + 0.004 * ((z - 1.3) / 0.5) ** 2
    sig = np.clip(sig, 0.005, 0.015)
    S = S_of_z(z)
    return z, sig, S


def draw_mock(z, sig, S, true):
    """Draw alpha from N(alpha_model, C_std + C_OU)."""
    mean = alpha_model(z, true["w0"], true["wa"])
    C = build_cov_OU(z, sig, S, true["theta"], true["sigma_X"])
    return RNG.multivariate_normal(mean, C)


# ---------------------------------------------------------------------------
# Priors and posterior (protocol §3)
# ---------------------------------------------------------------------------
# theta in log space for sampling stability: log10(theta)
# sigma_X in log space: log(sigma_X)

def log_prior(p):
    w0, wa, log10_theta, log_sx = p
    if not (-1.5 <= w0 <= -0.5):
        return -np.inf
    if not (-2.0 <= wa <= 1.0):
        return -np.inf
    theta = 10.0 ** log10_theta
    if not (THETA_MIN <= theta <= THETA_MAX):
        return -np.inf
    if not (LOG_SX_MIN <= log_sx <= LOG_SX_MAX):
        return -np.inf
    # flat in log10(theta) and log(sigma_X) inside bounds → constant
    return 0.0


def log_prob(p, z, sig, S, data):
    lp = log_prior(p)
    if not np.isfinite(lp):
        return -np.inf
    w0, wa, log10_theta, log_sx = p
    theta = 10.0 ** log10_theta
    sx = np.exp(log_sx)
    model = alpha_model(z, w0, wa)
    resid = data - model
    C = build_cov_OU(z, sig, S, theta, sx)
    return lp + log_like_gaussian(resid, C)



def run_mcmc(z, sig, S, data, nwalkers=32, nsteps=2000, discard=700, thin=4,
             fixed_background=None):
    """If fixed_background=(w0,wa), sample only (log10_theta, log_sx)."""
    if fixed_background is None:
        ndim = 4

        def log_prob_local(p):
            return log_prob(p, z, sig, S, data)

        p0 = []
        for _ in range(nwalkers):
            p0.append([
                -1.0 + 0.03 * RNG.normal(),
                0.0 + 0.05 * RNG.normal(),
                np.log10(0.3) + 0.25 * RNG.normal(),
                np.log(0.008) + 0.4 * RNG.normal(),
            ])
        p0 = np.array(p0)
        p0[:, 0] = np.clip(p0[:, 0], -1.49, -0.51)
        p0[:, 1] = np.clip(p0[:, 1], -1.9, 0.9)
        p0[:, 2] = np.clip(p0[:, 2], np.log10(THETA_MIN) + 1e-3, np.log10(THETA_MAX) - 1e-3)
        p0[:, 3] = np.clip(p0[:, 3], LOG_SX_MIN + 1e-3, LOG_SX_MAX - 1e-3)
    else:
        ndim = 2
        w0f, waf = fixed_background

        def log_prob_local(p):
            log10_theta, log_sx = p
            theta = 10.0 ** log10_theta
            sx = np.exp(log_sx)
            if not (THETA_MIN <= theta <= THETA_MAX):
                return -np.inf
            if not (LOG_SX_MIN <= log_sx <= LOG_SX_MAX):
                return -np.inf
            model = alpha_model(z, w0f, waf)
            resid = data - model
            C = build_cov_OU(z, sig, S, theta, sx)
            return log_like_gaussian(resid, C)

        p0 = []
        for _ in range(nwalkers):
            p0.append([
                np.log10(0.5) + 0.3 * RNG.normal(),
                np.log(0.008) + 0.4 * RNG.normal(),
            ])
        p0 = np.array(p0)
        p0[:, 0] = np.clip(p0[:, 0], np.log10(THETA_MIN) + 1e-3, np.log10(THETA_MAX) - 1e-3)
        p0[:, 1] = np.clip(p0[:, 1], LOG_SX_MIN + 1e-3, LOG_SX_MAX - 1e-3)

    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob_local)
    sampler.run_mcmc(p0, nsteps, progress=False)
    chain = sampler.get_chain(discard=discard, thin=thin, flat=True)
    logp = sampler.get_log_prob(discard=discard, thin=thin, flat=True)
    return chain, logp, sampler


def summarize_chain(chain, logp, true, label, fixed_background=False):
    if fixed_background:
        theta = 10.0 ** chain[:, 0]
        sx = np.exp(chain[:, 1])
        w0 = np.full_like(theta, true["w0"])
        wa = np.full_like(theta, true["wa"])
    else:
        w0, wa = chain[:, 0], chain[:, 1]
        theta = 10.0 ** chain[:, 2]
        sx = np.exp(chain[:, 3])

    def q(a):
        return {
            "median": float(np.median(a)),
            "p16": float(np.percentile(a, 16)),
            "p84": float(np.percentile(a, 84)),
            "p02": float(np.percentile(a, 2.5)),
            "p97": float(np.percentile(a, 97.5)),
        }

    sx_med = float(np.median(sx))
    sx_p02 = float(np.percentile(sx, 2.5))
    sx_p97 = float(np.percentile(sx, 97.5))
    sx_95u = float(np.percentile(sx, 95))
    th_p02 = float(np.percentile(theta, 2.5))
    th_med = float(np.median(theta))

    if sx_95u < 3e-3 and sx_med < 2e-3:
        region = "F_union_E0"
    elif sx_p02 > 3e-3 and th_p02 > 0.4:
        region = "E2_mean_reversion"
    elif sx_p02 > 3e-3:
        region = "E1_amplitude"
    else:
        region = "borderline_or_weak"

    return {
        "scenario": label,
        "fixed_background": fixed_background,
        "true": true,
        "region_decision": region,
        "w0": q(w0),
        "wa": q(wa),
        "theta": q(theta),
        "sigma_X": q(sx),
        "sigma_X_95_upper": sx_95u,
        "theta_95_lower": float(np.percentile(theta, 5)),
        "max_log_prob": float(np.max(logp)),
        "n_samples": int(chain.shape[0]),
    }, {"w0": w0, "wa": wa, "theta": theta, "sigma_X": sx, "logp": logp}


def corner_like_plot(samples, true, label, path):
    names = ["w0", "wa", "theta", "sigma_X"]
    truths = [true["w0"], true["wa"], true["theta"], true["sigma_X"]]
    fig, axes = plt.subplots(2, 2, figsize=(8, 7))
    axes = axes.ravel()
    for ax, name, truth in zip(axes, names, truths):
        ax.hist(samples[name], bins=40, color="C0", alpha=0.8, density=True)
        ax.axvline(truth, color="k", ls="--", lw=1.2, label="truth")
        ax.axvline(np.median(samples[name]), color="C3", ls="-", lw=1, label="median")
        ax.set_xlabel(name)
        ax.legend(fontsize=7)
    fig.suptitle(f"Euclid mock MCMC — {label}")
    fig.tight_layout()
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def main():
    print("Building Euclid-like grid...")
    z, sig, S = make_euclid_grid(n_bins=24)
    print(f"  N={len(z)}, z=[{z[0]:.2f},{z[-1]:.2f}], |S|~[{np.abs(S).min():.3f},{np.abs(S).max():.3f}]")

    # SNR-honest truths: C_OU comparable to C_std when sigma_X ~ 0.01
    scenarios = {
        "null_F_E0": {"w0": -1.0, "wa": 0.0, "theta": 0.5, "sigma_X": 1e-6},
        "E1_amplitude": {"w0": -1.0, "wa": 0.0, "theta": 0.1, "sigma_X": 0.012},
        "E2_mean_reversion": {"w0": -1.0, "wa": 0.0, "theta": 1.5, "sigma_X": 0.015},
    }

    all_summaries = {}
    nwalkers, nsteps, discard, thin = 28, 1600, 500, 3

    for name, true in scenarios.items():
        print(f"\n=== Scenario: {name} (fixed background first) ===")
        data = draw_mock(z, sig, S, true)
        np.save(os.path.join(OUT, f"mock_alpha_{name}.npy"), data)

        # Tier A: fixed LCDM background — validates OU kernel recovery
        chain_f, logp_f, _ = run_mcmc(
            z, sig, S, data,
            nwalkers=nwalkers, nsteps=nsteps, discard=discard, thin=thin,
            fixed_background=(true["w0"], true["wa"]),
        )
        # expand chain to 4-col form for storage: pad w0,wa
        chain_pad = np.column_stack([
            np.full(len(chain_f), true["w0"]),
            np.full(len(chain_f), true["wa"]),
            chain_f[:, 0],
            chain_f[:, 1],
        ])
        summary_f, samples_f = summarize_chain(
            chain_f, logp_f, true, name + "_fixedBG", fixed_background=True
        )
        all_summaries[name + "_fixedBG"] = summary_f
        np.savez(os.path.join(OUT, f"chain_{name}_fixedBG.npz"), chain=chain_pad, logp=logp_f)
        corner_like_plot(samples_f, true, name + " (fixed BG)", os.path.join(FIG, f"euclid_mcmc_{name}_fixedBG.png"))
        corner_like_plot(samples_f, true, name + " (fixed BG)", os.path.join(OUT, f"euclid_mcmc_{name}_fixedBG.png"))
        print(f"  fixedBG decision: {summary_f['region_decision']}")
        print(f"  sigma_X = {summary_f['sigma_X']['median']:.4f} [{summary_f['sigma_X']['p02']:.4f},{summary_f['sigma_X']['p97']:.4f}]")
        print(f"  theta   = {summary_f['theta']['median']:.3f} [{summary_f['theta']['p02']:.3f},{summary_f['theta']['p97']:.3f}]")

        # Tier B: free {w0,wa,theta,sigma_X}
        print(f"=== Scenario: {name} (free background) ===")
        chain, logp, _ = run_mcmc(
            z, sig, S, data,
            nwalkers=nwalkers, nsteps=nsteps, discard=discard, thin=thin,
            fixed_background=None,
        )
        summary, samples = summarize_chain(chain, logp, true, name + "_freeBG", fixed_background=False)
        all_summaries[name + "_freeBG"] = summary
        np.savez(os.path.join(OUT, f"chain_{name}_freeBG.npz"), chain=chain, logp=logp)
        corner_like_plot(samples, true, name + " (free BG)", os.path.join(FIG, f"euclid_mcmc_{name}_freeBG.png"))
        corner_like_plot(samples, true, name + " (free BG)", os.path.join(OUT, f"euclid_mcmc_{name}_freeBG.png"))
        print(f"  freeBG decision: {summary['region_decision']}")
        print(f"  sigma_X = {summary['sigma_X']['median']:.4f} [{summary['sigma_X']['p02']:.4f},{summary['sigma_X']['p97']:.4f}]")
        print(f"  theta   = {summary['theta']['median']:.3f} [{summary['theta']['p02']:.3f},{summary['theta']['p97']:.3f}]")
        print(f"  w0      = {summary['w0']['median']:.3f} [{summary['w0']['p02']:.3f},{summary['w0']['p97']:.3f}]")

    with open(os.path.join(OUT, "mcmc_summaries.json"), "w", encoding="utf-8") as f:
        json.dump(all_summaries, f, indent=2)

    lines = [
        "Euclid mock MCMC — vacuum-relaxation protocol",
        "=" * 50,
        f"N_bins={len(z)}, z in [{z[0]:.2f},{z[-1]:.2f}]",
        f"walkers={nwalkers}, steps={nsteps}, discard={discard}, thin={thin}",
        "Tier A: fixed (w0,wa)=truth. Tier B: free {w0,wa,theta,sigma_X}.",
        "Note: sigma_X ~ 0.01 is used for SNR-honest recovery under ~1% BAO errors;",
        "literature targets 1e-5..1.5e-4 require multi-survey power beyond this mock.",
        "",
    ]
    for name, s in all_summaries.items():
        lines.append(f"[{name}] decision={s['region_decision']}")
        lines.append(f"  truth={s['true']}")
        lines.append(
            f"  sigma_X={s['sigma_X']['median']:.4g} "
            f"[{s['sigma_X']['p02']:.4g},{s['sigma_X']['p97']:.4g}]"
        )
        lines.append(
            f"  theta={s['theta']['median']:.3g} "
            f"[{s['theta']['p02']:.3g},{s['theta']['p97']:.3g}]"
        )
        lines.append("")
    text = "\n".join(lines)
    print("\n" + text)
    with open(os.path.join(OUT, "mcmc_summary.txt"), "w", encoding="utf-8") as f:
        f.write(text + "\n")
    print("Wrote", OUT)


if __name__ == "__main__":
    main()
