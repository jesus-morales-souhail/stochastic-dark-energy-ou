# amplification routes 1–3 — numerical verdict

| | |
|:--|:--|
| **Author** | Jesús Morales Souhail · [github.com/jesus-morales-souhail](https://github.com/jesus-morales-souhail) |
| **Analysis** | Corpus July 2026 |
| **Parent note** | [`papers/amplification-gap.md`](../../papers/amplification-gap.md) §6.3 |
| **Scripts** | `scripts/amplification/route{1,2,3}_*.py`, `run_all_routes.py` |
| **Status** | Machine-checked hygiene freeze of route scans. **Not a DESI detection claim.** |

I freeze here the **machine-checked** conclusions of the amplification-route scans.
I do **not** invent a successful microphysical amplifier. This file is a numerical lock-in for theory bottlenecks — **not** a DESI claim that vacuum noise was seen.

---

## Route 1 — local / meso causal-set seed

**Files:** `route1_N_eff_required.csv`, `route1_N_eff_grid.csv` (1800 rows)

| Target | $\sigma_{0,\mathrm{eff}}$ | $N_{\mathrm{eff}}$ |
|:-------|:--------------------------|:-------------------|
| DESI 95% CL ceiling | $2.633\times 10^{-2}$ | $1.44\times 10^{3}$ (`route1_N_eff_required.csv` line 3) |
| Window | $1.0\times 10^{-4}$ | $1.00\times 10^{8}$ |
| Euclid residual scale | $1.0\times 10^{-5}$ | $1.00\times 10^{10}$ |

Global horizon count $N_{\mathrm{BH}}\sim 10^{122}$ is larger by $\sim 10^{112}$–$10^{119}$ for these $N_{\mathrm{eff}}$ ($10^{112.0}$ Euclid, line 2; $10^{114.0}$ window, line 4; $10^{118.8}$ DESI, line 3).

**Verdict:** Only route that reaches target **amplitudes** — by **redefining** $\sigma_0$, not by multiplying $10^{-61}$.

---

## Route 2 — late horizon-exit OU toy

**File:** `route2_horizon_exit_scan.csv` (24 configs × 4000 paths)

| Check | Result |
|:------|:-------|
| $\mathrm{rms}/\sigma$ | $0.36$–$0.97$($\mathcal{O}(1)$) |
| Freeze vs restore ratio | $1.000$ (no freeze bonus at $\Delta x\sim\mathcal{O}(1)$) |
| Residual scaling | $\propto\sigma$ |

**Verdict:** Freeze-out **preserves** seed; it does **not** amplify. Cannot close $10^{56}$ alone.

---

## Route 3 — double-well Langevin

**File:** `route3_avalanche_scan.csv` (288 jobs × 20 000 paths)
Grid: $a\in\{1,10\}$, $b\in[10^{-6},0.1]$, $\sigma\in[10^{-8},10^{-3}]$, $\Delta x=1$.

| Check | Result |
|:------|:-------|
| Gain $\mathrm{p95}/\sigma$ | $1.93$–$2.08$ (median $\approx 1.97$) |
| $\mathrm{p95}\ge 10^{-5}$ with $\sigma\le 10^{-8}$ | **0 / 288** |
| BAO-safe $\mathrm{p95}<2.633\times 10^{-2}$ | **288 / 288** (max p95 = $2.08\times 10^{-3}$, `route3_avalanche_scan.csv` line 289: $a=10$, $b=0.1$, $\sigma=10^{-3}$) |
| BAO-unsafe | **0 / 288** |

**BAO-safe (this scan):** the whole grid ($\sigma\le 10^{-3}$) is BAO-safe. *(Updated 2026-10-17 with the DESI bound $2.633\times10^{-2}$ from `results/profile_sigma_x/profile_sigma_x.json`; the former 216/288, 72/288 and $\sigma\lesssim 4\times10^{-5}$ were computed with the old $1.5\times10^{-4}$ target.)*

**Verdict:** Soft-potential avalanche is **not** a seed→Euclid free amplification. Large residuals track large $\sigma$.

---

## One-line programme status

**No free amplification confirmed numerically on Routes 2–3; Route 1 is a seed-redefinition card with a causal-set price.**
Default: publish **limits + exclusions + hypothesis map**.
I do not claim detection of Planck-scale vacuum noise.

*End of numerical freeze.*

