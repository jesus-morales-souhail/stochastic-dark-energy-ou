#!/usr/bin/env python3
"""
Reproduction of Ahmed, Dodelson, Greene, Sorkin (astro-ph/0209274), Sec. III.

Their eq. (5):  S_{i+1} = S_i + alpha xi_{i+1} sqrt(dN_i),  rho_L = S/N,
                dN_i = V(t_{i+1}) - V(t_i)               (eq. 4),
                V(t) = 4pi/3 int_0^t dt' a'^3 [int_t'^t dt''/a'']^3   (eq. 3, past light cone),
Friedmann only (eq. 6), radiation + matter + rho_L, S_0 = 0 at a = 1e-32. A run dies when
rho_tot <= 0.

Two unit conventions for "N = V" and rho_L:
  (R)  8 pi G = hbar = c = 1:  H^2 = rho / 3,      N = V / l_red^4
  (NR) G = hbar = c = 1:       H^2 = 8 pi rho / 3, N = V / l_P^4
Writing time in 1/H0, lengths in c/H0, Vh = V H0^4/c^4 and Sh = sum xi sqrt(dVh), every
power of (H0 t_unit) cancels and the dark-energy term in E^2 = H^2/H0^2 is
        u = K alpha Sh / Vh,     K = 1/3 (R)  or  8 pi / 3 (NR).
So the two conventions differ only by K (a factor 8 pi), and no 10^240 numbers appear.

V is integrated without cancellation through J_k(t) = int_0^t dt' a'^3 (eta - eta')^k:
        dJ_0/dt = a^3,  dJ_k/dt = k J_{k-1} / a,  V = 4pi/3 J_3,
with J_k(a_i) from pure radiation (a^2 proportional to t):
        J_k = a^{5+k} / C^{k+1} * 4! k! / (5+k)!,  C = sqrt(Omega_r).
Steps are uniform in ln a (h ~ 0.01); rho_L is held fixed within a step (RK4 for the J's)
and updated at the end of it, as in eq. (5). Background: Planck 2018 (n_count.py).
"""

from __future__ import annotations

import json
from math import factorial
from pathlib import Path

import numpy as np

from n_count import OMM, OMR

HERE = Path(__file__).resolve().parent
OUT = HERE / "resultados"
SEED = 20261016
A_INI = 1e-32
N_REAL = 200
H_STEP = 0.01
K = {"R": 1.0 / 3.0, "NR": 8.0 * np.pi / 3.0}
SIGMA_X = 0.025063155506929828  # results/profile_sigma_x/profile_sigma_x.json
Z_DESI = (0.3, 2.3)
Z_REC_MAX = 3.0  # record Omega_L(z) for z <= Z_REC_MAX


def run(alpha: float, conv: str, n_real: int = N_REAL, h: float = H_STEP, seed: int = SEED, a_ini: float = A_INI) -> dict:
    lna0 = np.log(a_ini)
    nsteps = int(round(-lna0 / h))
    h = -lna0 / nsteps
    k = K[conv]
    rng = np.random.default_rng(seed)
    xi = rng.standard_normal((nsteps, n_real))

    C = np.sqrt(OMR)
    a = a_ini
    J = np.array([a ** (5 + j) / C ** (j + 1) * factorial(4) * factorial(j) / factorial(5 + j) for j in range(4)])
    J = np.repeat(J[:, None], n_real, axis=1)
    V = 4 * np.pi / 3 * J[3]
    S = np.zeros(n_real)
    u = np.zeros(n_real)
    alive = np.ones(n_real, bool)
    a_death = np.full(n_real, np.nan)

    n_rec = int(np.ceil(np.log(1 + Z_REC_MAX) / h)) + 1
    rec_lna = lna0 + h * np.arange(nsteps + 1)[-n_rec:]
    rec = np.full((n_rec, n_real), np.nan)
    i_rec0 = nsteps + 1 - n_rec

    def rhs(lna, J, u):
        aa = np.exp(lna)
        E2 = OMR / aa**4 + OMM / aa**3 + u
        bad = E2 <= 0
        H = np.sqrt(np.where(bad, np.nan, E2))
        d = np.empty_like(J)
        d[0] = aa**3 / H
        for j in (1, 2, 3):
            d[j] = j * J[j - 1] / (aa * H)
        return d, bad

    for i in range(nsteps):
        lna = lna0 + i * h
        k1, b1 = rhs(lna, J, u)
        k2, b2 = rhs(lna + h / 2, J + h / 2 * k1, u)
        k3, b3 = rhs(lna + h / 2, J + h / 2 * k2, u)
        k4, b4 = rhs(lna + h, J + h * k3, u)
        Jn = J + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        died = alive & (b1 | b2 | b3 | b4 | ~np.all(np.isfinite(Jn), axis=0))
        Vn = 4 * np.pi / 3 * Jn[3]
        dV = np.where(died | ~alive, 0.0, Vn - V)
        Sn = S + xi[i] * np.sqrt(np.maximum(dV, 0.0))
        un = k * alpha * Sn / Vn
        an = np.exp(lna + h)
        E2n = OMR / an**4 + OMM / an**3 + un
        died |= alive & (E2n <= 0)
        a_death[died] = an
        alive &= ~died
        J = np.where(alive, Jn, J)
        V = np.where(alive, Vn, V)
        S = np.where(alive, Sn, S)
        u = np.where(alive, un, u)
        if i + 1 >= i_rec0:
            rec[i + 1 - i_rec0] = np.where(alive, un / E2n, np.nan)

    return {"alive": alive, "a_death": a_death, "rec_lna": rec_lna, "rec": rec, "OmL0": rec[-1], "nsteps": nsteps}


# ---------------------------------------------------------------------------------------------
# Das, Nasiri, Yazdi (arXiv:2304.03819, "Das I"), Model 1, followed literally. Kept in parallel
# to run() above (Ahmed-style, steps in ln a from a = 1e-32), which is not modified.
#   Lambda_i = 8 pi G S_i / V_i,  S_i = S_{i-1} + alpha sqrt(V_i - V_{i-1}) xi_i / G   (4.11)-(4.12)
#     => Lambda_i = 8 pi alpha sum_j sqrt(dV_j) xi_j / V_i                              (3.16), c = 1
#   H^2 = 8 pi G (rho_m + rho_r)/3 + Lambda/3                                           (4.2)
#   crash when H^2 < 0 (Sec. 4.3); steps uniform in a, a_i = i/n (Sec. 4.2);
#   start z ~ 1e4 with the volume up to then from matter + radiation only (Sec. 5.2),
#   S_0 = 0; n = 10000, i_initial = 1 (the values quoted in Sec. 5.3 and in Das II Sec. 4.1);
#   Omega_m h^2 = 0.14 (Sec. 5.1).
# In units H = 100 km/s/Mpc, lengths c/(100 km/s/Mpc):  H^2 = om a^-3 + orad a^-4 + u,
#   u = Lambda/3 = (8 pi alpha / 3) Sh / Vh  -> K_Das = 8 pi / 3 (= K["NR"] above).
# ---------------------------------------------------------------------------------------------
OMR_H2 = OMR * (0.6736**2)  # omega_r of n_count (photons + 3.046 massless nu)
K_DAS = 8.0 * np.pi / 3.0


def _j_init_mr(a0: float, om: float, orad: float) -> np.ndarray:
    """J_k(a0) = int_0^a0 a^3 (eta(a0) - eta(a))^k dt for a matter + radiation universe, no Lambda."""
    from scipy.integrate import quad

    H = lambda a: np.sqrt(om / a**3 + orad / a**4)
    eta = lambda a: quad(lambda x: 1.0 / (x * x * H(x)), 0.0, a, epsabs=0, epsrel=1e-12, limit=200)[0]
    e0 = eta(a0)
    return np.array([quad(lambda a: a**2 / H(a) * (e0 - eta(a)) ** k, 0.0, a0, epsabs=0, epsrel=1e-10, limit=200)[0]
                     for k in range(4)])


def _tau_init_mr(a0: float, om: float, orad: float) -> float:
    from scipy.integrate import quad

    return quad(lambda x: 1.0 / (x * x * np.sqrt(om / x**3 + orad / x**4)), 0.0, a0, epsabs=0, epsrel=1e-12)[0]


def run_das(alpha: float, n_real: int, seed: int = SEED, om: float = 0.14, orad: float = OMR_H2, n: int = 10000,
            i_ini: int = 1, z_rec: float | None = None, sub_rel: float = 0.02, z_obs=None,
            rd_wm=None, z_win=None) -> dict:
    """Model 1 as in Das I Sec. 4.2. Within each a-step Lambda is fixed and the J_k are integrated
    with RK4 substeps of relative size <= sub_rel (stands in for their Romberg integration).
    Optional: records u(a) for a >= 1/(1+z_rec); records conformal time tau and u at the grid
    point nearest to each z in z_obs (tau in units c/(100 km/s/Mpc)).
    rd_wm: omega_m values for which r_d is accumulated in the run from a_start to a_d(omega_m), with
    H_k^2 = w_k a^-3 + orad a^-4 + u w_k / om (Das I Sec. 4.4 scaling, radiation not rescaled);
    Simpson per step with the step's Lambda. Returns the part z < 1/a_start - 1 (no normalisation).
    z_win: z values where Omega_de = u / H^2 (at om) is stored."""
    rng = np.random.default_rng(seed)
    a_grid = np.arange(i_ini, n + 1) / n
    J0 = np.append(_j_init_mr(a_grid[0], om, orad), _tau_init_mr(a_grid[0], om, orad))  # row 4 = tau
    J = np.repeat(J0[:, None], n_real, axis=1)
    V = 4 * np.pi / 3 * J[3]
    S = np.zeros(n_real)
    u = np.zeros(n_real)
    alive = np.ones(n_real, bool)
    a_death = np.full(n_real, np.nan)
    i_rec = int(np.searchsorted(a_grid, 1.0 / (1.0 + z_rec))) if z_rec is not None else len(a_grid)
    rec_a = a_grid[i_rec:]
    rec_u = np.full((len(rec_a), n_real), np.nan)
    z_obs = np.atleast_1d(z_obs) if z_obs is not None else np.zeros(0)
    i_obs = np.array([int(np.argmin(np.abs(a_grid - 1.0 / (1.0 + z)))) for z in z_obs], dtype=int)
    obs_a = a_grid[i_obs] if len(i_obs) else np.zeros(0)
    obs_tau = np.full((len(i_obs), n_real), np.nan)
    obs_u = np.full((len(i_obs), n_real), np.nan)
    if rd_wm is not None:
        import rdcalc

        rd_wm = np.asarray(rd_wm, dtype=float)
        a_d = 1.0 / (1.0 + rdcalc.zd(rd_wm))
        rd_acc = np.zeros((len(rd_wm), n_real))
        a_max_rd = float(a_d.max())
    z_win = np.atleast_1d(z_win) if z_win is not None else np.zeros(0)
    i_win = np.array([int(np.argmin(np.abs(a_grid - 1.0 / (1.0 + z)))) for z in z_win], dtype=int)
    win_ode = np.full((len(i_win), n_real), np.nan)

    def rhs(a, J, u):
        H2 = om / a**3 + orad / a**4 + u
        bad = H2 <= 0
        H = np.sqrt(np.where(bad, np.nan, H2))
        d = np.empty_like(J)
        d[0] = a**2 / H
        for j in (1, 2, 3):
            d[j] = j * J[j - 1] / (a * a * H)
        d[4] = 1.0 / (a * a * H)
        return d, bad

    for i in range(1, len(a_grid)):
        a0, a1 = a_grid[i - 1], a_grid[i]
        if rd_wm is not None and a0 < a_max_rd:
            # r_d integrand c_s / (a^2 H_k) over [a0, min(a1, a_d_k)], Lambda of this step (u)
            b = np.minimum(a1, a_d)[:, None]
            ok_k = (a0 < a_d)[:, None]
            acc = np.zeros_like(rd_acc)
            am = 0.5 * (a0 + b)
            for aa, wgt in ((np.full_like(b, a0), 1.0), (am, 4.0), (b, 1.0)):
                H2 = rd_wm[:, None] / aa**3 + orad / aa**4 + u[None, :] * rd_wm[:, None] / om
                H = 100.0 * np.sqrt(np.where(H2 > 0, H2, np.nan))
                acc += wgt * rdcalc.cs(1.0 / aa - 1.0) / (aa * aa * H)
            rd_acc += np.where(ok_k & alive[None, :], acc * (b - a0) / 6.0, 0.0)
        m = max(1, int(np.ceil((a1 - a0) / (sub_rel * a0))))
        h = (a1 - a0) / m
        Jn = J.copy()
        bad = np.zeros(n_real, bool)
        for k in range(m):
            a = a0 + k * h
            k1, b1 = rhs(a, Jn, u)
            k2, b2 = rhs(a + h / 2, Jn + h / 2 * k1, u)
            k3, b3 = rhs(a + h / 2, Jn + h / 2 * k2, u)
            k4, b4 = rhs(a + h, Jn + h * k3, u)
            Jn = Jn + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
            bad |= b1 | b2 | b3 | b4
        xi = rng.standard_normal(n_real)
        died = alive & (bad | ~np.all(np.isfinite(Jn), axis=0))
        Vn = 4 * np.pi / 3 * Jn[3]
        dV = np.where(alive & ~died, Vn - V, 0.0)
        Sn = S + xi * np.sqrt(np.maximum(dV, 0.0))
        un = K_DAS * alpha * Sn / Vn
        died |= alive & (om / a1**3 + orad / a1**4 + un <= 0)
        a_death[died] = a1
        alive &= ~died
        J = np.where(alive, Jn, J)
        V = np.where(alive, Vn, V)
        S = np.where(alive, Sn, S)
        u = np.where(alive, un, u)
        if i >= i_rec:
            rec_u[i - i_rec] = np.where(alive, un, np.nan)
        for k in np.nonzero(i_obs == i)[0]:
            obs_tau[k] = J[4]
            obs_u[k] = u
        for k in np.nonzero(i_win == i)[0]:
            win_ode[k] = np.where(alive, u / (om / a1**3 + orad / a1**4 + u), np.nan)
    H0 = 100.0 * np.sqrt(om + orad + u)
    return {"alive": alive, "a_death": a_death, "rec_a": rec_a, "rec_u": rec_u, "om": om, "orad": orad,
            "H0": np.where(alive, H0, np.nan), "OmL0": np.where(alive, u / (om + orad + u), np.nan), "V0hat": V,
            "tau0": J[4], "obs_a": obs_a, "obs_tau": obs_tau, "obs_u": obs_u,
            "rd_part": rd_acc if rd_wm is not None else None, "win_ode": win_ode}


def main_das() -> None:
    from scipy.integrate import quad

    OUT.mkdir(exist_ok=True)
    out = {"source": "Das, Nasiri, Yazdi arXiv:2304.03819 Model 1", "seed": SEED, "omega_m": 0.14, "omega_r": OMR_H2,
           "n": 10000, "i_initial": 1, "a_start": 1e-4, "S0": 0.0}

    # Check H0^4 V for a universe without dark energy (Das I Sec. 5.1: 0.033)
    om = 0.14
    H = lambda a: np.sqrt(om / a**3 + OMR_H2 / a**4)
    J = _j_init_mr(1.0, om, OMR_H2)
    h_own = np.sqrt(om + OMR_H2)
    v_own = 4 * np.pi / 3 * J[3] * h_own**4
    eds = 36 * np.pi * 3 * factorial(8) * factorial(3) / factorial(12) * (2 / 3) ** 4
    # our earlier number: same m + r universe, V in units of the fiducial H0 = 67.36 (n_count.py)
    v_fid = 4 * np.pi / 3 * _j_init_mr(1.0, OMM, OMR)[3]
    out["check_V"] = {"H0own4_V_matter_radiation_om0.14": float(v_own), "H0_own_kms": float(100 * h_own),
                      "EdS_analytic": float(eds), "earlier_0.3435_in_fiducial_H0_units": float(v_fid),
                      "earlier_times_(H0own/H0fid)^4": float(v_fid * (OMM + OMR) ** 2)}
    print(f"H0^4 V sin energia oscura (om = 0.14, m + r): {v_own:.4f}  (H0 propio {100 * h_own:.2f} km/s/Mpc);"
          f" EdS analitico {eds:.4f}")
    print(f"  el 'V sin Lambda' anterior: {v_fid:.4f} en unidades de H0 fiducial 67.36; x (H0 propio/H0 fid)^4 ="
          f" {v_fid * (OMM + OMR) ** 2:.4f}")

    # Typical |Lambda|/(3H^2) = 8 pi alpha / (3 H^2 sqrt(V)) in radiation and matter domination (no Lambda)
    rd = np.sqrt(4 * np.pi / 3 * factorial(4) * factorial(3) / factorial(8))  # H^2 sqrt(V), pure radiation
    md = np.sqrt(eds)  # H0^2 sqrt(V), Einstein-de Sitter
    out["typical_amplitude"] = {"H2sqrtV_RD": float(rd), "H2sqrtV_MD": float(md),
                                "OmL_per_alpha_RD": float(K_DAS / rd), "OmL_per_alpha_MD": float(K_DAS / md)}
    print(f"H^2 sqrt(V): RD {rd:.4f}, MD {md:.4f};  |Lambda|/(3H^2) = {K_DAS / rd:.1f} alpha (RD), {K_DAS / md:.1f} alpha (MD)")

    das_ref = {0.005: "9990/10000", 0.01: "3228/10000", 0.012: "535/10000", 0.02: "<10/10000"}
    out["fractions"] = {}
    for al in (0.005, 0.01, 0.012, 0.02):
        r = run_das(al, 2000)
        k = int(r["alive"].sum())
        p = k / 2000
        err = np.sqrt(p * (1 - p) / 2000)
        h0 = r["H0"][r["alive"]]
        o0 = r["OmL0"][r["alive"]]
        out["fractions"][str(al)] = {"complete": k, "n": 2000, "frac": p, "binom_err": err, "das": das_ref[al],
                                     "H0_mean": float(np.mean(h0)) if k else None, "H0_std": float(np.std(h0)) if k else None,
                                     "OmL0_std": float(np.std(o0)) if k else None}
        print(f"  alpha {al}: completan {k}/2000 = {p:.4f} +- {err:.4f}   (Das I: {das_ref[al]})"
              + (f"   H0 {np.mean(h0):.1f} +- {np.std(h0):.1f}  sigma(OmL0) {np.std(o0):.3f}" if k else ""))
    print("Dependencia con la discretizacion (mismo arranque a = 1e-4):")
    out["discretization"] = {}
    for n, i_ini, nr in ((10000, 1, 2000), (100000, 10, 1000), (1000000, 100, 300)):
        for al in (0.005, 0.01, 0.012):
            r = run_das(al, nr, n=n, i_ini=i_ini)
            k = int(r["alive"].sum())
            h0 = r["H0"][r["alive"]]
            out["discretization"][f"n{n}_alpha{al}"] = {"complete": k, "n_seeds": nr, "frac": k / nr,
                                                       "H0_mean": float(np.mean(h0)), "H0_std": float(np.std(h0))}
            print(f"  n={n:>7} alpha={al}: {k}/{nr} = {k / nr:.4f}  H0 {np.mean(h0):.1f} +- {np.std(h0):.1f}")
    (OUT / "das_model1.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    print("Escrito resultados/das_model1.json")


def stats(res: dict) -> dict:
    ok = res["alive"]
    om0 = res["OmL0"][ok]
    out = {"n": int(len(ok)), "n_complete": int(ok.sum()), "frac_complete": float(ok.mean())}
    if ok.sum() >= 2:
        p16, p50, p84 = np.percentile(om0, [16, 50, 84])
        out.update(OmL0_p16=float(p16), OmL0_p50=float(p50), OmL0_p84=float(p84), OmL0_std=float(np.std(om0, ddof=1)),
                   OmL0_mean_abs=float(np.mean(np.abs(om0))))
        # Omega_L(z) on a grid uniform in ln(1+z), z in DESI range
        x_rec = -res["rec_lna"]  # ln(1+z) at recorded steps (ascending as index decreases)
        order = np.argsort(x_rec)
        xs, R = x_rec[order], res["rec"][order][:, ok]
        xg = np.linspace(np.log(1 + Z_DESI[0]), np.log(1 + Z_DESI[1]), 50)
        Og = np.array([np.interp(xg, xs, R[:, j]) for j in range(R.shape[1])]).T  # (50, n_ok)
        dO = Og - om0[None, :]
        per_z = np.std(dO, axis=1, ddof=1)
        o05 = np.array([np.interp(np.log(1.5), xs, R[:, j]) for j in range(R.shape[1])])
        o2 = np.array([np.interp(np.log(3.0), xs, R[:, j]) for j in range(R.shape[1])])
        corr = float(np.corrcoef(o05, o2)[0, 1])
        out.update(
            dOmL_DESI_std_pooled=float(np.std(dO, ddof=1)),
            dOmL_DESI_std_perz_min=float(per_z.min()),
            dOmL_DESI_std_perz_max=float(per_z.max()),
            corr_z05_z2=corr,
            theta_equiv=float(-np.log(corr) / np.log(2.0)) if corr > 0 else None,  # dx = ln(3/1.5)
            ratio_to_sigmaX=float(np.std(dO, ddof=1) / SIGMA_X),
        )
    return out


def process(res: dict, z_eff: np.ndarray) -> dict:
    """Normalisation-independent shape of the process, over runs that reach a = 1.

    (i)  r = sigma(Omega_L(z) - Omega_L(0), pooled over 0.3 < z < 2.3) / sigma(Omega_L(0))
    (ii) 7x7 correlation of Omega_L at the DESI z_eff, and the OU theta that best fits its
         off-diagonal entries with exp(-theta |d ln(1+z)|) (least squares).
    """
    from scipy.optimize import minimize_scalar

    s = stats(res)
    ok = res["alive"]
    x_rec = -res["rec_lna"]
    order = np.argsort(x_rec)
    xs, R = x_rec[order], res["rec"][order][:, ok]
    xz = np.log1p(z_eff)
    Oz = np.array([np.interp(xz, xs, R[:, j]) for j in range(R.shape[1])])  # (n_ok, 7)
    xg = np.linspace(np.log(1 + Z_DESI[0]), np.log(1 + Z_DESI[1]), 50)
    Og = np.array([np.interp(xg, xs, R[:, j]) for j in range(R.shape[1])])
    Cz = np.corrcoef(Oz.T)
    dx = np.abs(xz[:, None] - xz[None, :])
    iu = np.triu_indices(len(xz), 1)
    fit = minimize_scalar(lambda th: np.sum((Cz[iu] - np.exp(-th * dx[iu])) ** 2), bounds=(1e-3, 20), method="bounded")
    return {
        "n_complete": s["n_complete"],
        "sigma_OmL0": s["OmL0_std"],
        "sigma_dOmL_DESI": s["dOmL_DESI_std_pooled"],
        "r": s["dOmL_DESI_std_pooled"] / s["OmL0_std"],
        "r_level": float(np.std(Og, ddof=1) / s["OmL0_std"]),  # sigma(Omega_L(z)) itself, not minus today
        "z_eff": [float(z) for z in z_eff],
        "corr_7x7": Cz.tolist(),
        "theta_fit": float(fit.x),
        "theta_fit_rms_resid": float(np.sqrt(fit.fun / len(iu[0]))),
    }


def main() -> None:
    from desi_dr2_data import load_alpha_dv

    OUT.mkdir(exist_ok=True)
    summary = {"seed": SEED, "a_ini": A_INI, "n_real": N_REAL, "h_lna": H_STEP, "conventions": {}, "h_check": {},
               "f_scan": {}, "process": {}}

    print(f"Etapa 3a: alpha = 0.02, {N_REAL} realizaciones, seed {SEED}")
    for h in (H_STEP, 0.003):
        for conv in ("R", "NR"):
            r = run(0.02, conv, h=h)
            s = stats(r)
            ad = r["a_death"][~r["alive"]]
            if len(ad):
                s["a_death_p16_p50_p84"] = [float(v) for v in np.percentile(ad, [16, 50, 84])]
            (summary["conventions"] if h == H_STEP else summary["h_check"])[conv] = s
            line = f"  h={h:<5} {conv:>2}: completan {s['n_complete']}/{s['n']}"
            if "OmL0_std" in s:
                line += (f"  Omega_L(a=1) p16 {s['OmL0_p16']:+.4f} p50 {s['OmL0_p50']:+.4f} "
                         f"p84 {s['OmL0_p84']:+.4f} std {s['OmL0_std']:.4f}")
            if len(ad):
                line += f"  a_muerte p16/p50/p84 {np.percentile(ad, [16, 50, 84])}"
            print(line)
    for a_ini in (1.0 / 1100,):
        for conv in ("R", "NR"):
            s = stats(run(0.02, conv, a_ini=a_ini))
            summary["conventions"][f"{conv}_from_a{a_ini:.3g}"] = s
            print(f"  desde a={a_ini:.3g} {conv:>2}: completan {s['n_complete']}/{s['n']}"
                  + (f"  std {s['OmL0_std']:.4f}" if "OmL0_std" in s else ""))

    print("\nBarrido diagnostico: K = f * (1/3), alpha = 0.02")
    for f in (1, 2, 3, 4, 5, 6, 6.5, 7, 7.5, 8, 12, 8 * np.pi):
        K["f"] = K["R"] * f
        s = stats(run(0.02, "f"))
        summary["f_scan"][f"{f:.4g}"] = s
        print(f"  f={f:6.3g}: completan {s['n_complete']:3d}/200"
              + (f"  std {s['OmL0_std']:.3f}  <|OmL0|> {s['OmL0_mean_abs']:.3f}" if "OmL0_std" in s else ""))

    z_eff = load_alpha_dv(prefer_file=True, use_full_cov=True)["z"]
    print("\nForma del proceso (normalizacion-independiente)")
    for lab, conv, al, f in (("R_alpha0.02", "R", 0.02, None), ("R_alpha0.01", "R", 0.01, None), ("f7_alpha0.02", "f", 0.02, 7.0)):
        if f is not None:
            K["f"] = K["R"] * f
        p = process(run(al, conv), z_eff)
        summary["process"][lab] = p
        print(f"  {lab}: n={p['n_complete']}  sigma(OmL0)={p['sigma_OmL0']:.4f}  sigma(dOmL DESI)={p['sigma_dOmL_DESI']:.4f}"
              f"  r={p['r']:.4f}  r_level={p['r_level']:.4f}  theta_fit={p['theta_fit']:.4f} (rms resid {p['theta_fit_rms_resid']:.3f})")
    print("  corr 7x7 (R_alpha0.02), z_eff =", np.round(z_eff, 3))
    print(np.array2string(np.array(summary["process"]["R_alpha0.02"]["corr_7x7"]), precision=3, suppress_small=True))
    (OUT / "ahmed_sim.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print("\nEscrito resultados/ahmed_sim.json")


if __name__ == "__main__":
    import sys

    main_das() if len(sys.argv) > 1 and sys.argv[1] == "das" else main()
