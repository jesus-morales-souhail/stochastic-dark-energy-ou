#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Joint Euclid-like BAO + Type Ia SN mock MCMC for the vacuum-relaxation protocol.

Samples {w0, wa, theta, sigma_X} with:
  - BAO: C = C_std + C_OU(theta, sigma_X), Euclid-forecast-style errors
  - SNe: diagonal distance-modulus likelihood (Pantheon+-like binned), shared {w0, wa}

Scenarios: null (F∪E0), E1 (amplitude), E2 (mean-reversion).
Reports fixed-background and free-background tiers for BAO-only vs joint.

Author: Jesús Morales Souhail
ORCID:  0009-0000-7637-1818
"""

from __future__ import annotations

import os
import json
import numpy as np
from scipy.integrate import trapezoid
from scipy.linalg import cholesky, solve_triangular
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    import emcee
except ImportError as e:
    raise SystemExit("emcee is required: pip install emcee") from e

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "euclid_joint_mcmc")
FIG = os.path.join(ROOT, "figures")
os.makedirs(OUT, exist_ok=True)
os.makedirs(FIG, exist_ok=True)

H0 = 67.4
OM_FID = 0.315
C_KMS = 299792.458
THETA_MIN, THETA_MAX = 1e-3, 10.0
LOG_SX_MIN, LOG_SX_MAX = np.log(1e-6), np.log(3e-2)
RNG = np.random.default_rng(7)


# ---------------------------------------------------------------------------
# Cosmology
# ---------------------------------------------------------------------------
def E_z(z, Om, w0, wa):
    z_arr = np.atleast_1d(np.asarray(z, dtype=float))
    f = (1.0 + z_arr) ** (3.0 * (1.0 + w0 + wa)) * np.exp(
        -3.0 * wa * z_arr / (1.0 + z_arr)
    )
    out = np.sqrt(Om * (1.0 + z_arr) ** 3 + (1.0 - Om) * f)
    if np.ndim(z) == 0:
        return float(out[0])
    return out


def chi_of_z(z, Om, w0, wa, ngrid=600):
    z = float(z)
    if z <= 0:
        return 0.0
    zs = np.linspace(0.0, z, ngrid)
    invH = 1.0 / (H0 * E_z(zs, Om, w0, wa))
    return float(trapezoid(C_KMS * invH, zs))


def D_V(z, Om, w0, wa):
    chi = chi_of_z(z, Om, w0, wa)
    DH = C_KMS / (H0 * E_z(z, Om, w0, wa))
    return (float(z) * DH * chi**2) ** (1.0 / 3.0)


def d_L(z, Om, w0, wa):
    return (1.0 + float(z)) * chi_of_z(z, Om, w0, wa)


def mu_model(z, Om, w0, wa):
    # distance modulus; absolute calibration absorbed into relative fit via M free? 
    # Here we fix absolute scale to fiducial H0 (relative cosmology only).
    return 5.0 * np.log10(max(d_L(z, Om, w0, wa), 1e-6)) + 25.0


def alpha_model(z_arr, w0, wa, Om=OM_FID):
    out = np.empty(len(z_arr))
    for i, z in enumerate(z_arr):
        out[i] = D_V(z, Om, w0, wa) / D_V(z, Om, -1.0, 0.0)
    return out


def S_of_z(z_arr, Om=OM_FID, eps=1e-3):
    S = np.empty(len(z_arr))
    for i, z in enumerate(z_arr):
        Ol0 = 1.0 - Om
        dv_p = D_V(z, 1.0 - (Ol0 + eps), -1.0, 0.0)
        dv_m = D_V(z, 1.0 - (Ol0 - eps), -1.0, 0.0)
        dv0 = D_V(z, Om, -1.0, 0.0)
        S[i] = (dv_p - dv_m) / (2.0 * eps * dv0)
    return S


# ---------------------------------------------------------------------------
# Data vectors
# ---------------------------------------------------------------------------
def make_euclid_bao(n_bins=20, mode="optimistic"):
    """
    Euclid spectroscopic BAO-like vector.

    mode='optimistic': ~0.3–0.6% on alpha (schematic of high-SNR spectroscopic forecasts)
    mode='conservative': ~0.6–1.2% (closer to earlier internal mocks)
    """
    # Slightly wider than DR1 photo window: spectroscopic clustering 0.9–1.8
    z = np.linspace(0.9, 1.8, n_bins)
    if mode == "optimistic":
        # Better mid-range, slightly worse at edges
        sig = 0.0035 + 0.0025 * ((z - 1.35) / 0.45) ** 2
        sig = np.clip(sig, 0.003, 0.007)
    else:
        sig = 0.007 + 0.004 * ((z - 1.3) / 0.5) ** 2
        sig = np.clip(sig, 0.005, 0.014)
    S = S_of_z(z)
    return z, sig, S


def make_sn_bins(n_bins=25, zmin=0.1, zmax=1.5):
    """
    Compressed SN Ia distance-modulus bins (Pantheon+-like schematic).
    Effective per-bin errors ~ 0.02–0.04 mag after compression (not raw single-SN 0.15).
    """
    z = np.linspace(zmin, zmax, n_bins)
    # more SNe at low z -> slightly smaller error
    sig_mu = 0.025 + 0.02 * (z / zmax)
    return z, sig_mu


def build_cov_OU(z_arr, sigma_arr, S_arr, theta, sigma_X):
    n = len(z_arr)
    x = np.log(1.0 + z_arr)
    C = np.diag(sigma_arr.astype(float) ** 2)
    s2 = float(sigma_X) ** 2
    for i in range(n):
        for j in range(n):
            C[i, j] += S_arr[i] * S_arr[j] * s2 * np.exp(-theta * abs(x[i] - x[j]))
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


def draw_bao(z, sig, S, true):
    mean = alpha_model(z, true["w0"], true["wa"])
    C = build_cov_OU(z, sig, S, true["theta"], true["sigma_X"])
    return RNG.multivariate_normal(mean, C)


def draw_sn(z_sn, sig_mu, true):
    mu = np.array([mu_model(z, OM_FID, true["w0"], true["wa"]) for z in z_sn])
    return mu + RNG.normal(0.0, sig_mu)


# ---------------------------------------------------------------------------
# Posterior
# ---------------------------------------------------------------------------
def log_prior(p):
    w0, wa, log10_th, log_sx = p
    if not (-1.5 <= w0 <= -0.5):
        return -np.inf
    if not (-2.0 <= wa <= 1.0):
        return -np.inf
    th = 10.0 ** log10_th
    if not (THETA_MIN <= th <= THETA_MAX):
        return -np.inf
    if not (LOG_SX_MIN <= log_sx <= LOG_SX_MAX):
        return -np.inf
    return 0.0


def log_prob_joint(p, bao, sn, use_sn=True):
    lp = log_prior(p)
    if not np.isfinite(lp):
        return -np.inf
    w0, wa, log10_th, log_sx = p
    th = 10.0 ** log10_th
    sx = np.exp(log_sx)

    z, sig, S, data_a = bao
    model_a = alpha_model(z, w0, wa)
    resid_a = data_a - model_a
    C = build_cov_OU(z, sig, S, th, sx)
    ll = log_like_gaussian(resid_a, C)
    if not np.isfinite(ll):
        return -np.inf

    if use_sn:
        z_sn, sig_mu, data_mu = sn
        model_mu = np.array([mu_model(zz, OM_FID, w0, wa) for zz in z_sn])
        # absolute offset degeneracy with H0: subtract mean residual (relative SN)
        r = data_mu - model_mu
        r = r - np.mean(r)
        ll_sn = -0.5 * np.sum((r / sig_mu) ** 2 + np.log(2.0 * np.pi * sig_mu**2))
        ll += ll_sn
    return lp + ll


def run_mcmc(bao, sn, use_sn, nwalkers=28, nsteps=1400, discard=450, thin=3):
    ndim = 4
    p0 = np.array(
        [
            [
                -1.0 + 0.04 * RNG.normal(),
                0.0 + 0.08 * RNG.normal(),
                np.log10(0.4) + 0.25 * RNG.normal(),
                np.log(0.008) + 0.35 * RNG.normal(),
            ]
            for _ in range(nwalkers)
        ]
    )
    p0[:, 0] = np.clip(p0[:, 0], -1.49, -0.51)
    p0[:, 1] = np.clip(p0[:, 1], -1.9, 0.9)
    p0[:, 2] = np.clip(p0[:, 2], np.log10(THETA_MIN) + 1e-3, np.log10(THETA_MAX) - 1e-3)
    p0[:, 3] = np.clip(p0[:, 3], LOG_SX_MIN + 1e-3, LOG_SX_MAX - 1e-3)

    sampler = emcee.EnsembleSampler(
        nwalkers, ndim, log_prob_joint, args=(bao, sn, use_sn)
    )
    sampler.run_mcmc(p0, nsteps, progress=False)
    chain = sampler.get_chain(discard=discard, thin=thin, flat=True)
    logp = sampler.get_log_prob(discard=discard, thin=thin, flat=True)
    return chain, logp


def summarize(chain, true, label):
    w0, wa = chain[:, 0], chain[:, 1]
    th = 10.0 ** chain[:, 2]
    sx = np.exp(chain[:, 3])

    def q(a):
        return {
            "median": float(np.median(a)),
            "p02": float(np.percentile(a, 2.5)),
            "p16": float(np.percentile(a, 16)),
            "p84": float(np.percentile(a, 84)),
            "p97": float(np.percentile(a, 97.5)),
        }

    sx_med, sx_p02, sx_95u = np.median(sx), np.percentile(sx, 2.5), np.percentile(sx, 95)
    th_med, th_p02 = np.median(th), np.percentile(th, 2.5)

    if sx_95u < 3e-3 and sx_med < 2e-3:
        region = "F_union_E0"
    elif sx_p02 > 3e-3 and th_p02 > 0.4:
        region = "E2_mean_reversion"
    elif sx_p02 > 3e-3:
        region = "E1_amplitude"
    else:
        region = "borderline_or_weak"

    return {
        "label": label,
        "true": true,
        "region": region,
        "w0": q(w0),
        "wa": q(wa),
        "theta": q(th),
        "sigma_X": q(sx),
        "sigma_X_95_upper": float(sx_95u),
    }


def panel_plot(chain, true, title, path):
    th = 10.0 ** chain[:, 2]
    sx = np.exp(chain[:, 3])
    fig, axes = plt.subplots(1, 4, figsize=(12, 3.2))
    for ax, arr, name, tv in zip(
        axes,
        [chain[:, 0], chain[:, 1], th, sx],
        ["w0", "wa", "theta", "sigma_X"],
        [true["w0"], true["wa"], true["theta"], true["sigma_X"]],
    ):
        ax.hist(arr, bins=35, density=True, color="C0", alpha=0.8)
        ax.axvline(tv, color="k", ls="--", lw=1, label="truth")
        ax.axvline(np.median(arr), color="C3", lw=1, label="median")
        ax.set_xlabel(name)
        ax.legend(fontsize=6)
    axes[0].set_ylabel("density")
    fig.suptitle(title, fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def main():
    print("=== Joint Euclid BAO + SN mock MCMC ===")
    z_b, sig_b, S = make_euclid_bao(n_bins=20, mode="optimistic")
    z_sn, sig_mu = make_sn_bins(n_bins=25)
    print(
        f"BAO: N={len(z_b)}, sig~[{sig_b.min():.4f},{sig_b.max():.4f}] (optimistic forecast-style)"
    )
    print(f"SN:  N={len(z_sn)}, sig_mu~[{sig_mu.min():.3f},{sig_mu.max():.3f}] mag")

    scenarios = {
        "null_F_E0": {"w0": -1.0, "wa": 0.0, "theta": 0.5, "sigma_X": 1e-6},
        "E1_amplitude": {"w0": -1.0, "wa": 0.0, "theta": 0.1, "sigma_X": 0.008},
        "E2_mean_reversion": {"w0": -1.0, "wa": 0.0, "theta": 1.5, "sigma_X": 0.012},
    }

    summaries = {}
    nwalkers, nsteps = 28, 1500

    for name, true in scenarios.items():
        print(f"\n--- {name} ---")
        data_a = draw_bao(z_b, sig_b, S, true)
        data_mu = draw_sn(z_sn, sig_mu, true)
        bao = (z_b, sig_b, S, data_a)
        sn = (z_sn, sig_mu, data_mu)
        np.savez(
            os.path.join(OUT, f"mock_{name}.npz"),
            alpha=data_a,
            mu=data_mu,
            z_bao=z_b,
            z_sn=z_sn,
            true=np.array([true["w0"], true["wa"], true["theta"], true["sigma_X"]]),
        )

        # BAO-only
        print("  BAO-only MCMC...")
        chain_b, _ = run_mcmc(bao, sn, use_sn=False, nwalkers=nwalkers, nsteps=nsteps)
        s_b = summarize(chain_b, true, name + "_BAO_only")
        summaries[s_b["label"]] = s_b
        panel_plot(
            chain_b,
            true,
            f"{name} — BAO only (optimistic Euclid)",
            os.path.join(FIG, f"joint_mcmc_{name}_BAO.png"),
        )
        panel_plot(
            chain_b,
            true,
            f"{name} — BAO only",
            os.path.join(OUT, f"joint_mcmc_{name}_BAO.png"),
        )
        print(
            f"  BAO-only: region={s_b['region']} "
            f"sx={s_b['sigma_X']['median']:.4f} "
            f"th={s_b['theta']['median']:.3f} "
            f"w0={s_b['w0']['median']:.3f}"
        )

        # Joint BAO+SN
        print("  BAO+SN MCMC...")
        chain_j, _ = run_mcmc(bao, sn, use_sn=True, nwalkers=nwalkers, nsteps=nsteps)
        s_j = summarize(chain_j, true, name + "_BAO_SN")
        summaries[s_j["label"]] = s_j
        panel_plot(
            chain_j,
            true,
            f"{name} — BAO + SN joint",
            os.path.join(FIG, f"joint_mcmc_{name}_BAO_SN.png"),
        )
        panel_plot(
            chain_j,
            true,
            f"{name} — BAO + SN joint",
            os.path.join(OUT, f"joint_mcmc_{name}_BAO_SN.png"),
        )
        print(
            f"  BAO+SN:   region={s_j['region']} "
            f"sx={s_j['sigma_X']['median']:.4f} "
            f"th={s_j['theta']['median']:.3f} "
            f"w0={s_j['w0']['median']:.3f}"
        )

    # Comparison figure: w0 width BAO vs joint
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.4))
    for ax, name in zip(axes, scenarios.keys()):
        # reload from last chains via re-summary files only — re-plot from stored if needed
        # use summaries for text; for hist need chains — re-load mock and skip, show text bars
        sb = summaries[name + "_BAO_only"]
        sj = summaries[name + "_BAO_SN"]
        ax.bar(
            [0, 1],
            [sb["w0"]["p84"] - sb["w0"]["p16"], sj["w0"]["p84"] - sj["w0"]["p16"]],
            color=["C0", "C1"],
            tick_label=["BAO", "BAO+SN"],
        )
        ax.set_title(name.replace("_", " "))
        ax.set_ylabel(r"68% width on $w_0$")
    fig.suptitle("Background tightening: BAO only vs BAO+SN", y=1.03)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "joint_mcmc_w0_width_compare.png"), dpi=140, bbox_inches="tight")
    fig.savefig(os.path.join(OUT, "joint_mcmc_w0_width_compare.png"), dpi=140, bbox_inches="tight")
    plt.close()

    with open(os.path.join(OUT, "joint_summaries.json"), "w", encoding="utf-8") as f:
        json.dump(summaries, f, indent=2)

    lines = [
        "Joint Euclid-like BAO + SN mock MCMC",
        "=" * 50,
        "BAO: optimistic forecast-style errors (~0.3–0.7% on alpha).",
        "SN:  25 compressed mu bins, relative fit (mean residual removed).",
        "Parameters: {w0, wa, theta, sigma_X} with OU kernel on BAO only.",
        "",
    ]
    for k, s in summaries.items():
        lines.append(f"[{k}] region={s['region']}")
        lines.append(f"  true={s['true']}")
        lines.append(
            f"  sigma_X={s['sigma_X']['median']:.4g} "
            f"[{s['sigma_X']['p02']:.4g},{s['sigma_X']['p97']:.4g}]"
        )
        lines.append(
            f"  theta={s['theta']['median']:.3g} "
            f"[{s['theta']['p02']:.3g},{s['theta']['p97']:.3g}]"
        )
        lines.append(
            f"  w0={s['w0']['median']:.3f} "
            f"[{s['w0']['p02']:.3f},{s['w0']['p97']:.3f}]  "
            f"wa={s['wa']['median']:.3f} "
            f"[{s['wa']['p02']:.3f},{s['wa']['p97']:.3f}]"
        )
        lines.append("")
    lines.append(
        "Interpretation: SN primarily tightens {w0,wa}; OU parameters still "
        "require BAO covariance SNR. Joint analysis reduces false absorption "
        "of OU signal into smooth DE when SNR is sufficient."
    )
    text = "\n".join(lines)
    print("\n" + text)
    with open(os.path.join(OUT, "joint_summary.txt"), "w", encoding="utf-8") as f:
        f.write(text + "\n")
    print("Wrote", OUT)


if __name__ == "__main__":
    main()
