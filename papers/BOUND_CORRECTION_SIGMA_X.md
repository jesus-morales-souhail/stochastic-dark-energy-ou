# Bound correction: DESI residual ceiling

**Date:** August 2026  
**Source of truth:** `scripts/profile_sigma_x_desi.py` → `results/profile_sigma_x/`

## What changed

| | Former (working target) | Corrected (profile 95% CL) |
|--|-------------------------|----------------------------|
| Number | $\sigma_X < 1.5\times 10^{-4}$ | $\sigma_X < 2.5\times 10^{-2}$ |
| Exact grid | — | $2.506\times 10^{-2}$ |
| Criterion | not equal to free-$\theta$ profile | $\Delta\ln\mathcal{L}\ge -1.92$ (1 dof) from profile max |
| Data | 7 isotropic $\alpha$ | same |
| Cov | earlier diagonal / aspirational | full Gaussian BAO $13\times13$ projected to $\alpha$ |

## What did **not** change

- Best-fit residual: $\sigma_X\to 0$ (null still preferred).
- Rank-1 tachyonic exclusion when active.
- Euclid target $10^{-5}$ for amplification gap $\sim 10^{56}$.

## Why the number moved

With only 7 compressed BAO $\alpha$ bins and free mean-reversion $\theta$, the free-$\theta$ profile likelihood is weakly constraining on amplitude. The former $1.5\times 10^{-4}$ was carried as a **working target** / tighter programme number; it was **not** the profile 95% CL of this script. Aligning claim C2 to the script avoids overstating precision.

## Downstream

- Amplification gap to DESI ceiling: $\sim 10^{59}$ (was $\sim 10^{57}$).
- Three-gate G1 headroom: $\mathrm{ceiling}/\sigma_{\mathrm{free}}\approx 2.5\times 10^{-2}/8.5\times 10^{-5}\approx 294$ (was $\sim 1.76$).
- Slip upper bound from residual ceiling is weaker (larger ceiling) — still far below Maus/Sakr if one uses the *measured* null preference rather than the upper edge alone.

