#!/usr/bin/env python3
"""
Exploration only (EXPLORACION.md, point 1): can a single Model 1 run of bao_sim.py m1 be re-simulated
from (alpha, chunk, idx) alone, and does it reproduce what the npz stored?

bao_sim.main_m1 calls ahmed_sim.run_das(alpha, 20000, seed=SEED + chunk, ...). Inside run_das each step
draws rng.standard_normal(n_real) for the whole chunk, and the columns evolve independently (the RK4
substep count depends only on the a-grid). So run idx is reproduced by running run_das with n_real = 1
and an rng that, at each call, draws the full 20000 normals of the chunk and returns column idx.
z_rec = 1/a_start - 1 makes run_das keep the full u(a) history (rec_a, rec_u), which bao_sim did not save.

Usage: python3 resim_m1.py ALPHA NPZ_ROW [NPZ_ROW ...]
Writes cmb_real/historias_m1/alpha<ALPHA>_row<NPZ_ROW>.npz (a, u, H0, ...) and prints the comparison.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
EV = HERE.parent
sys.path.insert(0, str(EV))

import ahmed_sim as AS  # noqa: E402
import bao_sim as B  # noqa: E402

CHUNK = 20000


class _Column:
    """rng stand-in: each standard_normal(n) call draws the whole chunk and returns column idx."""

    def __init__(self, seed: int, idx: int, width: int = CHUNK, factory=np.random.default_rng):
        self.g = factory(seed)
        self.idx, self.width = idx, width

    def standard_normal(self, n):
        return self.g.standard_normal(self.width)[[self.idx]]


def resim(alpha: float, chunk: int, idx: int) -> dict:
    real = np.random.default_rng
    AS.np.random.default_rng = lambda seed: _Column(seed, idx, factory=real)  # only run_das calls it here
    try:
        zu = np.unique(B.ZR)
        r = AS.run_das(alpha, 1, seed=B.SEED + chunk, om=B.OM, orad=B.ORAD, n=100000, i_ini=10,
                       z_obs=np.append(zu, B.Z_REC), z_rec=1.0 / 1e-4 - 1.0,
                       rd_wm=np.unique(np.append(np.geomspace(0.10, 0.20, 24), B.OM)),
                       z_win=np.geomspace(1000.0, 1.0e4, 25))
    finally:
        AS.np.random.default_rng = real
    a = r["obs_a"][:, None]
    tau, u = r["obs_tau"], r["obs_u"]
    H = 100.0 * np.sqrt(B.OM / a**3 + B.ORAD / a**4 + u)
    C = B.C
    DM = C / 100.0 * (r["tau0"][None, :] - tau[:-1])
    iz = np.searchsorted(zu, B.ZR)
    chi2 = B.chi2_from(DM[iz], (C / H[:-1])[iz])
    return {"r": r, "chi2": float(chi2[0]), "H0": float(r["H0"][0]), "OmL0": float(r["OmL0"][0]),
            "Ode_rec": float(u[-1, 0] / (H[-1, 0] / 100.0) ** 2), "alive": bool(r["alive"][0])}


def main() -> None:
    alpha = float(sys.argv[1])
    R = np.load(EV / "resultados" / f"bao_m1_alpha{alpha}.npz")
    (HERE / "historias_m1").mkdir(exist_ok=True)
    for row in map(int, sys.argv[2:]):
        chunk, idx = int(R["seed_chunk"][row] - B.SEED), int(R["idx"][row])
        t0 = time.time()
        s = resim(alpha, chunk, idx)
        dt = time.time() - t0
        r = s["r"]
        cmp = {k: (s[k], float(R[k][row])) for k in ("chi2", "H0", "OmL0", "Ode_rec")}
        dwin = float(np.nanmax(np.abs(r["win_ode"][:, 0] - R["win_ode"][row])))
        drd = float(np.max(np.abs(r["rd_part"][:, 0] - R["rd_part"][row])))
        print(f"alpha {alpha} row {row} (chunk {chunk}, idx {idx}): alive {s['alive']}, {dt:.0f} s")
        for k, (x, y) in cmp.items():
            print(f"   {k:8s} resim {x:.10g}  npz {y:.10g}  diff {x - y:.2e}")
        print(f"   win_ode max|diff| {dwin:.2e}   rd_part max|diff| {drd:.2e}")
        np.savez_compressed(HERE / "historias_m1" / f"alpha{alpha}_row{row}.npz", a=r["rec_a"], u=r["rec_u"][:, 0],
                            om=B.OM, orad=B.ORAD, H0=s["H0"], chi2=s["chi2"], chunk=chunk, idx=idx, seconds=dt)


if __name__ == "__main__":
    main()
