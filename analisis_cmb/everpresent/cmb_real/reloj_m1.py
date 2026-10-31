#!/usr/bin/env python3
"""
Exploration only (RELOJ_M1.md): does the survival of Model 1 (Das I Sec. 4.2, ahmed_sim.run_das) depend on how the
steps are distributed? The walk S_i = S_{i-1} + alpha sqrt(dV_i) xi_i is a Wiener process in V; Lambda is held fixed
within a step and H^2 > 0 is checked at every RK4 substep, as in run_das.

Clocks (same start a = 1e-4, S_0 = 0, K = 8 pi / 3, alpha = 0.01, omega_m = 0.14301, omega_r of ahmed_sim):
  a   steps uniform in a: the grid of run_das, a = arange(i_ini, n + 1) / n
  lna steps uniform in ln a, same number of steps
  V   steps uniform in the simulation's own 4-volume: dV = (V_ref(a=1) - V(a_start)) / n, with V_ref the 4-volume of
      the same universe without dark energy; each realisation advances from its own state with the local power law
      V ~ a^p, a1 = a (1 + dV/V)^(1/p) (approximately uniform in V; the actual dV per step is recorded), and steps
      are clipped to land on the DESI redshifts and on a = 1 (not on z_rec, which would add an early step).
Random numbers: realisation i uses np.random.default_rng(20261016 + i); its k-th step uses its k-th normal, in every
clock. (Section 1 of EVERPRESENT_vs_DESI.md used one generator for all columns, default_rng(20261016), and
omega_m = 0.14; that table is reproduced separately with run_das itself.)
BAO: chi^2 with r_d from each run (accumulated in the run for the 25 omega_m of rd_real.py: Simpson per step in the
grid clocks, as run_das; per RK4 substep in the V clock, whose early steps are long) and omega_m within +-1.6 % of
Planck: rd_real.minimise, as criterion (c) without window of EVERPRESENT_vs_DESI.md 3e.

Checks: the generalised loop with the run_das grid and run_das's generator reproduces run_das (alive, H0, conformal
time, r_d); run_das reproduces the 13.2 % (n = 1e4) and 27.5 % (n = 1e5) of section 1 with that section's seeds.

    python3 reloj_m1.py            # writes cmb_real/reloj_m1.json
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
EV = HERE.parent
sys.path.insert(0, str(EV))

import ahmed_sim as AS  # noqa: E402
import bao_sim as B  # noqa: E402
import rd_real as RR  # noqa: E402
import rdcalc as RC  # noqa: E402
import robustez as RB  # noqa: E402

ALPHA, OM, ORAD = 0.01, B.OM, B.ORAD
A_START = 1e-4
N_REAL = 2000
SEED0 = 20261016
WM_RD = RR.WM_RD
ZU = np.unique(B.ZR)
A_OBS = 1.0 / (1.0 + np.append(ZU, B.Z_REC))  # z order, as run_das (z_rec last; it is not used in the chi^2)
A_EQ = ORAD / OM


class Streams:
    """One generator per realisation (seed SEED0 + i); next() returns the k-th normal of every stream."""

    def __init__(self, n: int, seed0: int = SEED0, block: int = 2000):
        self.g = [np.random.default_rng(seed0 + i) for i in range(n)]
        self.block, self.k, self.buf = block, 0, None

    def next(self) -> np.ndarray:
        if self.buf is None or self.k == self.block:
            self.buf = np.stack([g.standard_normal(self.block) for g in self.g], axis=1)
            self.k = 0
        self.k += 1
        return self.buf[self.k - 1]


class Single:
    """run_das's generator: one default_rng(seed), n normals per step."""

    def __init__(self, n: int, seed: int):
        self.r, self.n = np.random.default_rng(seed), n

    def next(self) -> np.ndarray:
        return self.r.standard_normal(self.n)


def rhs(a, J, u, om, orad):
    H2 = om / a**3 + orad / a**4 + u
    bad = H2 <= 0
    H = np.sqrt(np.where(bad, np.nan, H2))
    d = np.empty_like(J)
    d[0] = a**2 / H
    for j in (1, 2, 3):
        d[j] = j * J[j - 1] / (a * a * H)
    d[4] = 1.0 / (a * a * H)
    return d, bad


def simpson_rd(a0, b, u, om, orad, a_d):
    """r_d integrand c_s/(a^2 H_k) on [a0, min(b, a_d_k)] for the 25 omega_m (Das I Sec. 4.4 scaling), as run_das."""
    bb = np.minimum(b[None, :], a_d[:, None])
    ok = a0[None, :] < a_d[:, None]
    acc = np.zeros((len(a_d), len(u)))
    am = 0.5 * (a0[None, :] + bb)
    for aa, wgt in ((np.broadcast_to(a0[None, :], bb.shape), 1.0), (am, 4.0), (bb, 1.0)):
        H2 = WM_RD[:, None] / aa**3 + orad / aa**4 + u[None, :] * WM_RD[:, None] / om
        H = 100.0 * np.sqrt(np.where(H2 > 0, H2, np.nan))
        acc += wgt * RC.cs(1.0 / aa - 1.0) / (aa * aa * H)
    return np.where(ok, acc * (bb - a0[None, :]) / 6.0, 0.0)


def run_grid(a_grid: np.ndarray, gen, n_real: int, om: float = OM, orad: float = ORAD, sub_rel: float = 0.02) -> dict:
    """run_das's loop on an arbitrary common grid (observations at the nearest grid point, as run_das)."""
    J0 = np.append(AS._j_init_mr(a_grid[0], om, orad), AS._tau_init_mr(a_grid[0], om, orad))
    J = np.repeat(J0[:, None], n_real, axis=1)
    V = 4 * np.pi / 3 * J[3]
    S, u = np.zeros(n_real), np.zeros(n_real)
    alive = np.ones(n_real, bool)
    i_obs = np.array([int(np.argmin(np.abs(a_grid - x))) for x in A_OBS])
    obs_tau = np.full((len(A_OBS), n_real), np.nan)
    obs_u = np.full((len(A_OBS), n_real), np.nan)
    a_d = 1.0 / (1.0 + RC.zd(WM_RD))
    rd_acc = np.zeros((len(WM_RD), n_real))
    n_before_eq = int(np.sum(a_grid[1:] <= A_EQ))
    for i in range(1, len(a_grid)):
        a0, a1 = a_grid[i - 1], a_grid[i]
        if a0 < a_d.max():
            rd_acc += np.where(alive[None, :], simpson_rd(np.full(n_real, a0), np.full(n_real, a1), u, om, orad, a_d), 0.0)
        m = max(1, int(np.ceil((a1 - a0) / (sub_rel * a0))))
        h = (a1 - a0) / m
        Jn = J.copy()
        bad = np.zeros(n_real, bool)
        for k in range(m):
            a = a0 + k * h
            k1, b1 = rhs(a, Jn, u, om, orad)
            k2, b2 = rhs(a + h / 2, Jn + h / 2 * k1, u, om, orad)
            k3, b3 = rhs(a + h / 2, Jn + h / 2 * k2, u, om, orad)
            k4, b4 = rhs(a + h, Jn + h * k3, u, om, orad)
            Jn = Jn + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
            bad |= b1 | b2 | b3 | b4
        xi = gen.next()
        died = alive & (bad | ~np.all(np.isfinite(Jn), axis=0))
        Vn = 4 * np.pi / 3 * Jn[3]
        dV = np.where(alive & ~died, Vn - V, 0.0)
        Sn = S + xi * np.sqrt(np.maximum(dV, 0.0))
        un = AS.K_DAS * ALPHA * Sn / Vn
        died |= alive & (om / a1**3 + orad / a1**4 + un <= 0)
        alive &= ~died
        J = np.where(alive, Jn, J)
        V = np.where(alive, Vn, V)
        S = np.where(alive, Sn, S)
        u = np.where(alive, un, u)
        for kk in np.nonzero(i_obs == i)[0]:
            obs_tau[kk] = J[4]
            obs_u[kk] = u
    return {"alive": alive, "u": u, "tau0": J[4], "obs_a": a_grid[i_obs], "obs_tau": obs_tau, "obs_u": obs_u,
            "rd_part": rd_acc, "steps": len(a_grid) - 1, "steps_before_eq": np.full(n_real, n_before_eq)}


def run_v(n: int, gen, n_real: int, om: float = OM, orad: float = ORAD, sub_rel: float = 0.02) -> dict:
    """Steps of (approximately) constant dV in each realisation's own 4-volume."""
    J0 = np.append(AS._j_init_mr(A_START, om, orad), AS._tau_init_mr(A_START, om, orad))
    v_start = 4 * np.pi / 3 * J0[3]
    v_ref1 = 4 * np.pi / 3 * AS._j_init_mr(1.0, om, orad)[3]
    dV_target = (v_ref1 - v_start) / n
    J = np.repeat(J0[:, None], n_real, axis=1)
    a = np.full(n_real, A_START)
    V = np.full(n_real, v_start)
    S, u = np.zeros(n_real), np.zeros(n_real)
    alive = np.ones(n_real, bool)
    done = np.zeros(n_real, bool)
    obs_tau = np.full((len(A_OBS), n_real), np.nan)
    obs_u = np.full((len(A_OBS), n_real), np.nan)
    stops_a = np.append(np.sort(1.0 / (1.0 + ZU)), 1.0)  # DESI redshifts (ascending a) and a = 1; not z_rec
    stop_obs = [int(np.argmin(np.abs(A_OBS - x))) for x in stops_a[:-1]]  # index of each stop in A_OBS (z order)
    next_obs = np.zeros(n_real, int)
    a_d = 1.0 / (1.0 + RC.zd(WM_RD))
    rd_acc = np.zeros((len(WM_RD), n_real))
    steps = np.zeros(n_real, int)
    steps_eq = np.zeros(n_real, int)
    dv_rel = []
    while True:
        act = alive & ~done
        if not act.any():
            break
        H = np.sqrt(np.maximum(om / a**3 + orad / a**4 + u, 1e-300))
        dVda = 4 * np.pi / 3 * 3 * J[2] / (a * a * H)
        p = a * dVda / V  # local power law V ~ a^p
        a1 = a * (1.0 + dV_target / V) ** (1.0 / p)
        stop = stops_a[np.minimum(next_obs, len(stops_a) - 1)]
        a1 = np.where(act, np.minimum(a1, stop), a)
        da = a1 - a
        m = max(1, int(np.ceil(np.max(np.where(act, da / (sub_rel * a), 0.0)))))
        h = da / m
        Jn = J.copy()
        bad = np.zeros(n_real, bool)
        for k in range(m):
            aa = a + k * h
            if np.any(act & (aa < a_d.max())):
                rd_acc += np.where(act[None, :], simpson_rd(aa, aa + h, u, om, orad, a_d), 0.0)
            k1, b1 = rhs(aa, Jn, u, om, orad)
            k2, b2 = rhs(aa + h / 2, Jn + h / 2 * k1, u, om, orad)
            k3, b3 = rhs(aa + h / 2, Jn + h / 2 * k2, u, om, orad)
            k4, b4 = rhs(aa + h, Jn + h * k3, u, om, orad)
            Jn = np.where(act, Jn + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4), Jn)
            bad |= act & (b1 | b2 | b3 | b4)
        xi = gen.next()
        died = act & (bad | ~np.all(np.isfinite(Jn), axis=0))
        Vn = 4 * np.pi / 3 * Jn[3]
        dV = np.where(act & ~died, Vn - V, 0.0)
        if act.sum():
            dv_rel.append(float(np.median(dV[act & ~died] / dV_target)) if (act & ~died).any() else np.nan)
        Sn = S + xi * np.sqrt(np.maximum(dV, 0.0))
        un = AS.K_DAS * ALPHA * Sn / Vn
        died |= act & (om / a1**3 + orad / a1**4 + un <= 0)
        alive &= ~died
        ok = act & alive
        J = np.where(ok, Jn, J)
        V = np.where(ok, Vn, V)
        S = np.where(ok, Sn, S)
        u = np.where(ok, un, u)
        a = np.where(ok, a1, a)
        steps += ok
        steps_eq += ok & (a1 <= A_EQ)
        hit = ok & (next_obs < len(stop_obs)) & (a >= stops_a[np.minimum(next_obs, len(stops_a) - 1)] - 1e-15)
        for c in np.nonzero(hit)[0]:
            obs_tau[stop_obs[next_obs[c]], c] = J[4, c]
            obs_u[stop_obs[next_obs[c]], c] = u[c]
        next_obs = next_obs + hit
        done |= ok & (a >= 1.0 - 1e-15)
    return {"alive": alive, "u": u, "tau0": J[4], "obs_a": A_OBS, "obs_tau": obs_tau, "obs_u": obs_u, "rd_part": rd_acc,
            "steps": steps, "steps_before_eq": steps_eq, "dV_target": dV_target,
            "dV_over_target_median_p5_p95": [float(np.nanmedian(dv_rel)), float(np.nanpercentile(dv_rel, 5)),
                                             float(np.nanpercentile(dv_rel, 95))]}


def summary(r: dict, om: float = OM, orad: float = ORAD) -> dict:
    ok = r["alive"]
    H0 = 100.0 * np.sqrt(om + orad + r["u"])
    OmL = r["u"] / (om + orad + r["u"])
    out = {"n_real": int(len(ok)), "n_alive": int(ok.sum()), "frac_alive": float(ok.mean()),
           "steps_median_survivors": float(np.median(np.broadcast_to(r["steps"], ok.shape)[ok])),
           "steps_before_eq_median_survivors": float(np.median(np.broadcast_to(r["steps_before_eq"], ok.shape)[ok])),
           "binomial_sd": float(np.sqrt(ok.mean() * (1 - ok.mean()) / len(ok)))}
    if ok.sum() >= 2:
        out.update(H0_p16_p50_p84=[float(v) for v in np.percentile(H0[ok], [16, 50, 84])], H0_std=float(np.std(H0[ok])),
                   OmL0_p16_p50_p84=[float(v) for v in np.percentile(OmL[ok], [16, 50, 84])], OmL0_std=float(np.std(OmL[ok])))
        # BAO chi^2 with own r_d, omega_m within +-1.6 % (rd_real.model1 recipe)
        a = r["obs_a"][:, None]
        tau, u = r["obs_tau"][:, ok], r["obs_u"][:, ok]
        H = 100.0 * np.sqrt(om / a**3 + orad / a**4 + u)
        DM = (B.C / 100.0) * (r["tau0"][ok][None, :] - tau[:-1])
        H0o = H0[ok][None, :]
        dmh, eh = (DM * H0o / B.C).T, (H0o / H[:-1]).T
        tail = np.array([RC.rd_standard(w, 1.0 / A_START - 1.0) for w in WM_RD])
        rd = RC.NORM * (r["rd_part"][:, ok].T + tail[None, :])
        t_rd0 = RB.theory13(dmh, eh, H0[ok], ZU)
        i0 = int(np.argmin(np.abs(WM_RD - OM)))
        t0 = t_rd0 * (RC.RD_CAMB / rd[:, i0])[None, :]
        chi_c, wm_c = RR.minimise(t0, rd, RR.WM_C)
        j = int(np.argmin(chi_c))
        out["bao_chi2_c_min"] = float(chi_c[j])
        out["bao_chi2_c_min_omega_m"] = float(wm_c[j])
        out["bao_n_below_B"] = int(np.sum(chi_c < B.CHI2_LCDM + 4))
        out["bao_chi2_c_p50"] = float(np.median(chi_c[np.isfinite(chi_c)]))
    return out


def grid(clock: str, n: int) -> np.ndarray:
    i_ini = int(round(n * A_START))
    if clock == "a":
        return np.arange(i_ini, n + 1) / n
    return np.exp(np.linspace(np.log(A_START), 0.0, n - i_ini + 1))


def main() -> None:
    out = {"alpha": ALPHA, "omega_m": OM, "omega_r": ORAD, "a_start": A_START, "n_real": N_REAL,
           "seeds": f"default_rng({SEED0} + i), i = 0..{N_REAL - 1}", "a_eq": A_EQ}
    # check 1: section 1 of EVERPRESENT_vs_DESI.md with its own seeds and omega_m = 0.14 (run_das itself)
    t0 = time.time()
    s1 = {}
    for n, i_ini, nr in ((10000, 1, 2000), (100000, 10, 1000)):
        r = AS.run_das(ALPHA, nr, n=n, i_ini=i_ini)
        s1[f"n{n}"] = {"alive": int(r["alive"].sum()), "n": nr, "frac": float(r["alive"].mean())}
    out["check_section1_run_das"] = s1
    print("seccion 1:", s1, f"{time.time() - t0:.0f} s", flush=True)
    # check 2: the generalised loop with run_das's grid and generator reproduces run_das (n = 1e4, 2000 columns)
    ref = AS.run_das(ALPHA, N_REAL, seed=SEED0, om=OM, orad=ORAD, n=10000, i_ini=1, z_obs=np.append(ZU, B.Z_REC),
                     rd_wm=WM_RD)
    mine = run_grid(grid("a", 10000), Single(N_REAL, SEED0), N_REAL)
    ok = ref["alive"]
    out["check_run_grid_vs_run_das"] = {
        "alive_equal": bool(np.array_equal(ok, mine["alive"])),
        "u_max_abs_diff_alive": float(np.max(np.abs(ref["H0"][ok] - 100 * np.sqrt(OM + ORAD + mine["u"][ok])))),
        "tau_obs_max_rel_diff": float(np.nanmax(np.abs(mine["obs_tau"][:, ok] / ref["obs_tau"][:, ok] - 1))),
        "rd_part_max_rel_diff": float(np.nanmax(np.abs(mine["rd_part"][:, ok] / ref["rd_part"][:, ok] - 1)))}
    print("run_grid vs run_das:", out["check_run_grid_vs_run_das"], flush=True)
    # the three clocks
    out["clocks"] = {}
    alive = {}
    for n in (10000, 100000):
        for clock in ("a", "lna", "V"):
            t0 = time.time()
            gen = Streams(N_REAL)
            r = run_v(n, gen, N_REAL) if clock == "V" else run_grid(grid(clock, n), gen, N_REAL)
            s = summary(r)
            if clock == "V":
                s["dV_over_target_median_p5_p95"] = r["dV_over_target_median_p5_p95"]
            s["seconds"] = time.time() - t0
            alive[f"{clock}_n{n}"] = r["alive"]
            out["clocks"][f"{clock}_n{n}"] = s
            print(clock, n, json.dumps(s), flush=True)
            (HERE / "reloj_m1.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    # paired comparison (same seeds): survivors in common
    keys = list(alive)
    out["overlap_survivors"] = {f"{x}|{y}": int(np.sum(alive[x] & alive[y])) for i, x in enumerate(keys) for y in keys[i + 1:]}
    (HERE / "reloj_m1.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    print(json.dumps(out["overlap_survivors"]))


if __name__ == "__main__":
    main()
