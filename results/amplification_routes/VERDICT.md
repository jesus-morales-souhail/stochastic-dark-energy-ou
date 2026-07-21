# Act III routes 1–3 — numerical verdict

| | |
|:--|:--|
| **Analysis** | Corpus July 2026 |
| **Parent note** | [`papers/amplification-gap.md`](../../papers/amplification-gap.md) §6.3 |
| **Scripts** | `scripts/amplification/route{1,2,3}_*.py`, `run_all_routes.py` |

This file freezes the **machine-checked** conclusions of the amplification-route scans. 
It does **not** invent a successful microphysical amplifier.

---

## Route 1 — local / meso causal-set seed

**Files:** `route1_N_eff_required.csv`, `route1_N_eff_grid.csv` (1800 rows)

| Target | $\sigma_{0,\mathrm{eff}}$ | $N_{\mathrm{eff}}$ |
|:-------|:--------------------------|:-------------------|
| DESI 95% CL ceiling | $1.5\times 10^{-4}$ | $4.44\times 10^{7}$ |
| Window | $1.0\times 10^{-4}$ | $1.00\times 10^{8}$ |
| Euclid residual scale | $1.0\times 10^{-5}$ | $1.00\times 10^{10}$ |

Global horizon count $N_{\mathrm{BH}}\sim 10^{122}$ is larger by $\sim 10^{112}$–$10^{114}$ for these $N_{\mathrm{eff}}$.

**Verdict:** Only route that reaches target **amplitudes** — by **redefining** $\sigma_0$, not by multiplying $10^{-61}$.

---

## Route 2 — late horizon-exit OU toy

**File:** `route2_horizon_exit_scan.csv` (24 configs × 4000 paths)

| Check | Result |
|:------|:-------|
| $\mathrm{rms}/\sigma$ | $0.36$–$0.97$ ($\mathcal{O}(1)$) |
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
| BAO-safe $\mathrm{p95}<1.5\times 10^{-4}$ | **216 / 288** |
| BAO-unsafe | **72 / 288**, all with $\sigma\gtrsim 1.2\times 10^{-4}$ |

**BAO-safe rule of thumb (this scan):** $\sigma\lesssim 4\times 10^{-5}$.

**Verdict:** Soft-potential avalanche is **not** a seed→Euclid free amplification. Large residuals track large $\sigma$.

---

## One-line programme status

**No free amplification confirmed numerically on Routes 2–3; Route 1 is a seed-redefinition card with a causal-set price.** 
Default: publish **limits + exclusions + hypothesis map**.
