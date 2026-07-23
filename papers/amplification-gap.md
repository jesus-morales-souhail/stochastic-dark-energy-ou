# Amplification gap: Sorkin seed to BAO residual

**Author:** Jesús Morales Souhail  
**Date:** July 2026  
**Status:** Act III closure note — theory bottleneck, not a detection claim  
**Filename:** `papers/amplification-gap.md`  
**Former name:** `papers/amplification-no-free-lunch.md` ([redirect stub](amplification-no-free-lunch.md))  

**Related:**  
`papers/fundamental-vs-emergent-vacuum-relaxation.md`,  
`notes/desqueezing-relaxation-vacuum-fluctuations-note.md`,  
`papers/anisotropic-slip-option0.md`,  
`papers/scale-operator-experiment-map.md`,  
`papers/resume.txt`,  
`scripts/amplifier_audit.py`,  
`results/amplification_routes/VERDICT.md`

---

## 1. The only scientifically open front (after Acts I–II, IV–V)

Acts I–II deliver a **null residual** under the OU/QNM kernel and a **geometric candidate** (SDiff) for isotropic smoothness. 
Acts IV–V close **wrong operators** and **wrong scales** (slip amplitude starvation, Born/band-limit, tesseract \(B_4\) seal, wavefront \(T\neq\) OU).

The corpus already states the remaining hole:

> *This note does not … identify a microphysical amplifier from \(10^{-61}\) to \(10^{-5}\).* 
> — `fundamental-vs-emergent-vacuum-relaxation.md`

This note **fixes that gap as a quantitative problem** and audits candidate amplifiers. It does **not** invent a successful mechanism.

---

## 2. Quantified bottleneck

| Quantity | Value | Role |
|----------|-------|------|
| Sorkin / Bekenstein seed | \(\sigma_0 \sim 1.18\times 10^{-61}\) (\(H_0=67.4\)) | Motivational UV amplitude (holographic \(d=2\)) |
| DESI DR2 BAO residual bound | \(\sigma_X < 1.5\times 10^{-4}\) (95% CL) | **Measured** effective late-time amplitude (OU kernel) |
| Euclid-scale residual target | \(\sim 10^{-5}\) | **Aspirational** path-integrated residual of interest |
| Minimum kick for detection (small \(\theta\)) | \(A_0^{\min}\sim 10^{-5}\) | From \(\sigma_{\rm res}=A_0 e^{-\theta\Delta x}\) |

### 2.1 Two gaps (do not mix labels)

These are **different questions** with answers that differ by about **one decade**. Never paste “\(\sim 10^{56}\)” onto a calculation that used the DESI ceiling.

| Gap | Definition | Exact OOM (\(\sigma_0=1.18\times 10^{-61}\)) | Soft \(r\) needed (\(e^{2r}\sigma_0=\sigma_{\mathrm{tgt}}\)) |
|:----|:-----------|:---------------------------------------------|:-----------------------------------------------------------|
| **\(G_{\mathrm{Euclid}}\)** | \(10^{-5}/\sigma_0\) | \(\mathbf{8.47\times 10^{55}\approx 10^{56}}\) | \(r\approx 64.4\) |
| **\(G_{\mathrm{DESI}}\)** | \(1.5\times 10^{-4}/\sigma_0\) | \(\mathbf{1.27\times 10^{57}\approx 10^{57}}\) | \(r\approx 65.7\) |

**Label rule (mandatory):**

- “\(\sim 10^{56}\)” **only** for the **Euclid-scale target** \(\sigma\sim 10^{-5}\) (or \(A_0^{\min}\sim 10^{-5}\)).  
- “\(\sim 10^{57}\)” for the **measured DESI ceiling** \(\sigma_X=1.5\times 10^{-4}\).  
- Density form of the same ratio: \(\lvert\delta\rho\rvert_{\max}/(\sigma_0\rho_\Lambda)=\sigma_{\mathrm{tgt}}/\sigma_0\) — so DESI density gap is \(\sim 10^{57}\), not \(10^{56}\).

**What DESI actually measured** is a **null residual bound**, not a new \(\rho_\Lambda\). Mean vacuum density \(\rho_\Lambda\sim\Omega_{\Lambda 0}\rho_{\mathrm{crit}}\) is imported from standard cosmology (e.g. Planck-class \(H_0,\Omega_\Lambda\)); this corpus constrains **fractional residual amplitude** \(\sigma_X\) on top of that background.

DESI constrains **effective** amplitude after any amplification — **not** \(\sigma_0\) itself.  
Without a physical map \(\sigma_0\to A_0\), “Euclid will see vacuum noise” is **not** a theorem of this repository.

Runnable audit: `python scripts/amplifier_audit.py` · exact ratios: `python scripts/gap_two_targets.py`.  
Headlines discipline: [`HONEST_HEADLINES.md`](HONEST_HEADLINES.md).

---

## 3. Amplifier audit (classes vs closed bounds)

| Amplifier | Schematic gain | Verdict in this corpus |
|-----------|----------------|------------------------|
| **(a)** Coherent / tachyonic growth | Exponential in \(t/t_c\) | **Excluded** in-repo with correct rank-1 covariance (\(\Delta\chi^2\sim +23\)) — `quantum-fluid-instabilities-desi-dr2.md` |
| **(b)** Desqueezing, \(r=1.5\) | \(e^{2r}\approx 20\sim 10^{1}\) | **Orders short** (~55 decades missing to \(10^{56}\)) |
| **(c)** RMS accumulate \(\sqrt{N}\) | \(\sqrt{N}\sim 10^{61}\) | Seed is **already** \(1/\sqrt{N}\); free \(\sqrt{N}\) boost **double-counts** unless growth is **coherent** → reduces to (a) |
| **(d)** Inflation-style freeze-out (~60 e-folds) | \(e^{60}\sim 10^{26}\) | Only \(\mathcal{O}(10^{26})\); and the natural seed is \(H/M_{\mathrm{Pl}}\sim 10^{-5}\), **not** \(\sigma_0\sim 10^{-61}\). Mapping inflationary noise onto Sorkin is a **new ansatz** |
| **(e)** Phase-transition jump | Discrete feature | Tension with **smooth** BAO residual structure if large |

**No audited channel** maps bare Sorkin \(\sigma_0\) to Euclid/DESI residual scales without:

1. already being **excluded** (a), 
2. delivering only \(\mathcal{O}(1)\)–\(\mathcal{O}(10^{26})\) (b,d), 
3. reintroducing a **free large \(A_0\)** (renaming the problem), or 
4. **abandoning** the Sorkin seed definition (d).

Every closed channel that looked like a “shortcut” — anisotropic slip, casual GRB phase stories, pupil/tesseract optics — **inherits this same gap**. That is why they are amplitude-starved even when the **operator** is right (Option 0).

---

## 4. The only candidate that deserves a serious *new* theory exam: (d)

### 4.0 Concrete literature (revalidated): Gordon & Wands (2005)

**Reference:** C. Gordon & D. Wands, Phys. Rev. D **71**, 123505 (2005); arXiv:[astro-ph/0504132](https://arxiv.org/abs/astro-ph/0504132).  
**Theory-repo expansion:** [measurable-stochastic-vacuum `inflation-spectator-seed-gordon-wands.md`](https://github.com/jesus-morales-souhail/measurable-stochastic-vacuum/blob/main/papers/inflation-spectator-seed-gordon-wands.md).

They **do not** start from Sorkin \(\sigma_0\sim 10^{-61}\). Light spectator: \(\mathcal{P}_Q^{1/2}=H_{\mathrm{inf}}/(2\pi)\). For *their* large isocurvature / quadrupole target, frozen evolution conflicts with tensor bounds; they require \(\delta Q_f/\delta Q_i>45\) (Mexican-hat radial roll) — **factor \(\sim 45\)**, not \(10^{56}\).

**Do not** equate that “tachyonic” roll with the **excluded** coherent GPE/tachyonic residual on BAO in this corpus (different seed, different observable, different math).  
**Do not** claim it closes DESI residual detectability without a new seed→BAO map.  
**Do** treat it as a respected, structurally distinct open door: missing factor \(\mathcal{O}(10^{1}\)–\(10^{2})\) + unmeasured \(H_{\mathrm{inf}}\), not free \(G\sim 10^{56}\) on holographic counting.

### 4.1 Why (d) is special

Inflation is the **only** standard paradigm known to convert sub-horizon quantum fluctuations into large-scale classical curvature perturbations with a huge effective boost over many e-folds. Schematic amplitude after freeze-out is set by \(H/M_{\mathrm{Pl}}\) (and slow-roll factors), not by \(1/\sqrt{N_{\mathrm{BH,today}}}\).

### 4.2 Why (d) is **not** a free amplification for this repo

| Requirement | Status |
|-------------|--------|
| Gain \(\sim 10^{56}\) from \(\sigma_0\) alone | **Fails** — \(e^{60}\sim 10^{26}\) is ~30 decades short of \(10^{56}\) if one insists on starting from \(\sigma_0\) |
| Use \(H/M_{\mathrm{Pl}}\sim 10^{-5}\) as the seed | **Changes the microphysical story** — no longer pure late-time Sorkin Poisson counting |
| Late-time “horizon exit” analogue for a **homogeneous DE residual** \(X(x)\) | **Not derived** — inflation applies to **perturbation modes** freezing outside the Hubble radius, not automatically to a background OU residual on \(x=\ln a\) |
| Compatibility with BAO smoothness | Any large residual after freeze-out must still satisfy \(\sigma_X < 1.5\times 10^{-4}\) |

**Honest statement:** 
(d) is worth a **separate theory paper** only if one (i) redefines the seed, (ii) derives a late-time freeze-out map for the DE residual, and (iii) re-runs BAO constraints on that map. It is **not** a patch that rescues \(\sigma_0\sim 10^{-61}\) inside the current OU likelihood.

### 4.3 What not to do

- Do not claim desqueezing \(r\sim\mathcal{O}(1)\) closes the gap. 
- Do not claim slip or GRBs will “see” Sorkin seeds without \(A_0\). 
- Do not smuggle \(\sqrt{N}\) twice. 
- Do not treat free \(A_0\sim 10^{-5}\) as microphysics.

---

## 5. How Act III sits in the five-act narrative

| Act | Content | Amplification |
|-----|---------|----------------|
| I | BAO null + model kills | Bounds **effective** amplitude |
| II | SDiff geometry | Explains **isotropic** smallness without a particle |
| **III** | **This note** | **No free amplification** from \(\sigma_0\) to \(10^{-5}\) |
| IV | Slip Option 0 | Right operator; still needs \(A_0\) |
| V | Method no-gos | Wrong scale/operator closed |

**Programme status if (d) fails as a late-time DE residual theory:**

- The repository remains what it already is: **honest limits, exclusions, and geometric interpretation**. 
- Amplification stays **formally open** as a theory problem. 
- That is a successful scientific posture, not a failed discovery.

---

## 6. Three *allowed* working hypotheses (each with a price)

Linear gain \(\times 10^{56}\) on white noise is dead. The only coherent ways to **change the problem** (not multiply noise by magic) are:

| Route | Mechanism | How it “reaches” \(10^{-5}\) | Price / open problem | Script |
|-------|-----------|------------------------------|----------------------|--------|
| **1. Scale-invariant / local seed** | Sorkin count uses \(N_{\mathrm{eff}}\ll N_{\mathrm{BH}}\) (local / meso correlation) | Raise **starting** \(\sigma_{0,\mathrm{eff}}=1/\sqrt{N_{\mathrm{eff}}}\) | Redefine causal-set boundary counting | `scripts/amplification/route1_local_causal_set_seed.py` |
| **2. Late horizon exit** | \(\theta(x)\to 0\) after mode exits horizon during late acceleration; freeze residual | Stretch / freeze **before** restoration kills the kick | Must derive DE-mode freeze-out (≠ 60 inflation e-folds; late \(\Delta x=\mathcal{O}(1)\)) | `scripts/amplification/route2_late_horizon_exit.py` |
| **3. Nonlinear avalanche** | Tiny seed is only a **trigger**; large jump from multi-well / threshold dynamics paid by background potential | Burst amplitude set by potential, not \(e^{2r}\) | New \(V(X)\); must stay BAO-smooth (not always-large) | `scripts/amplification/route3_nonlinear_avalanche.py` |

**Batch runner:** `python scripts/amplification/run_all_routes.py` (add `--heavy` for denser grids). 
**Outputs:** `results/amplification_routes/*.csv`.

### 6.1 Reading of Route 2 vs classical inflation

Late-time acceleration supplies \(\Delta x=\ln(a_0/a_{\mathrm{exit}})\sim\mathcal{O}(1)\), **not** \(\sim 60\). 
Even with perfect freeze, a bare \(\sigma_0\sim 10^{-61}\) does **not** become \(10^{-5}\) from e-fold stretch alone. Route 2 only becomes powerful if combined with Route 1 (larger seed) or a **non-Sorkin** mode amplitude \(\sim H/M_{\mathrm{Pl}}\).

### 6.2 Decision tree

```
Need detection-scale residual without free A0?
 ├─ Redefine seed (Route 1)? → justify N_eff; re-run BAO with new prior
 ├─ Freeze-out (Route 2)? → derive θ(k,a); note Δx_late ≪ 60
 ├─ Avalanche (Route 3)? → derive V(X); scan BAO-safe rare-burst region
 └─ None → programme = honest limits + exclusions (success)
```

**Recommended default for publication of the present corpus:**  
Act III closed as an **amplification gap on linear amplifiers** (no free gain factor across \(\sim 10^{56}\)); Routes 1–3 listed as **the only open theory cards**, each with a price.

### 6.3 Numerical scan results (HEAD — machine tables)

**Data:** `results/amplification_routes/` 
`route1_N_eff_required.csv`, `route1_N_eff_grid.csv` (1800 rows), 
`route2_horizon_exit_scan.csv` (24 configs × 4000 paths), 
`route3_avalanche_scan.csv` (288 jobs × 20 000 paths). 
**Compact verdict file:** `results/amplification_routes/VERDICT.md`.

#### Route 1 — redefining \(N_{\mathrm{eff}}\) (not a gain factor)

Pure Poisson counting \(\sigma_{0,\mathrm{eff}}=1/\sqrt{N_{\mathrm{eff}}}\):

| Target residual | \(\sigma_{0,\mathrm{eff}}\) | \(N_{\mathrm{eff}}\) required | vs \(N_{\mathrm{BH}}\sim 10^{122}\) |
|-----------------|-----------------------------|-------------------------------|-------------------------------------|
| DESI 95% ceiling | \(1.5\times 10^{-4}\) | \(\mathbf{4.44\times 10^{7}}\) | smaller by \(\sim 10^{114}\) |
| Intermediate window | \(10^{-4}\) | \(1.00\times 10^{8}\) | \(\sim 10^{114}\) |
| Euclid-scale residual | \(10^{-5}\) | \(\mathbf{1.00\times 10^{10}}\) | smaller by \(\sim 10^{112}\) |

Implied local scale if \(N_{\mathrm{eff}}=N_{\mathrm{BH}}(L/R_H)^p\) (order-of-magnitude only):

| Target | \(p=2\) (area) \(L\) | \(p=3\) (volume) \(L\) | \(p=4\) (4-vol) \(L\) |
|--------|----------------------|------------------------|------------------------|
| DESI \(N_{\mathrm{eff}}\sim 4.4\times 10^{7}\) | \(\sim 10^{-31} \mathrm{m}\) | \(\sim 10^{-12} \mathrm{m}\) | \(\sim\mathrm{mm}\) |
| Euclid \(N_{\mathrm{eff}}\sim 10^{10}\) | \(\sim 10^{-30} \mathrm{m}\) | \(\sim 10^{-11} \mathrm{m}\) | \(\sim\mathrm{cm}\) |

**Verdict R1:** The **only** card that reaches DESI/Euclid residual *amplitudes* without linear gain. **Price:** abandon global horizon counting; justify a meso-scale causal-set / correlation volume. That is a **new theory claim**, not a free amplification on \(\sigma_0\sim 10^{-61}\).

#### Route 2 — late freeze-out toy (96k paths total)

SDE \(dX=-\theta(x)X dx+\sigma dW\) with \(\theta\) dropping after a late \(x_{\mathrm{exit}}\) (\(\Delta x\sim\mathcal{O}(1)\)).

| Metric | Result |
|--------|--------|
| \(\mathrm{rms}(X_f)/\sigma\) | \(\sim 0.36\)–\(0.97\) (always \(\mathcal{O}(1)\)) |
| Freeze vs restore (\(\theta_{\mathrm{super}}=0\) vs \(10^{-4}\)) | ratio \(\mathbf{1.000}\) at fixed \(\theta_{\mathrm{sub}}\) |
| Scaling | residual \(\propto\sigma\) (no decade gain) |
| Smallest \(\sigma\) in scan | \(10^{-6}\) (float-safe); bare \(10^{-61}\) not representable |

**Verdict R2:** Late-time freeze with \(\Delta x\sim\mathcal{O}(1)\) **does not amplify**. It at best *preserves* an already large seed. Without Route 1 or a non-Sorkin \(H/M_{\mathrm{Pl}}\) seed, Route 2 **cannot** close the \(10^{56}\) gap. Matches §6.1.

#### Route 3 — double-well Langevin avalanche (288 × 20k)

Potential \(V=\tfrac14 a X^4-\tfrac12 b X^2\), \(a\in\{1,10\}\), \(b\in[10^{-6},0.1]\), \(\sigma\in[10^{-8},10^{-3}]\), \(\Delta x=1\).

| Metric | Result |
|--------|--------|
| Gain \(\mathrm{p95}(|X|)/\sigma\) | **min 1.93 · median 1.97 · max 2.08** |
| Linearized growth ceiling \(e^{b\Delta x}\) at \(b_{\max}=0.1\) | \(\approx 1.11\) (scan not in deep-avalanche regime) |
| Jobs with \(\mathrm{p95}\ge 10^{-5}\) and \(\sigma\le 10^{-8}\) | **0 / 288** |
| Jobs with \(\mathrm{p95}\ge 10^{-5}\) and \(\sigma < 10^{-6}\) | **0 / 288** |
| BAO-safe (\(\mathrm{p95}<1.5\times 10^{-4}\)) | **216 / 288** |
| BAO-unsafe (\(\mathrm{p95}\ge 1.5\times 10^{-4}\)) | **72 / 288** — **all** have \(\sigma\gtrsim 1.2\times 10^{-4}\) |
| Euclid window \(10^{-5}\le\mathrm{p95}<1.5\times 10^{-4}\) | 72 jobs — residual still \(\sim\sigma\), not seed-triggered |

**BAO-safe rule of thumb in this scan:** keep drive noise \(\sigma\lesssim 4\times 10^{-5}\) (p95 stays below DESI ceiling for all \((a,b)\) tested). 
**BAO-unsafe:** \(\sigma\gtrsim 1.2\times 10^{-4}\) → p95 exceeds DESI for **every** \((a,b)\).

**Verdict R3:** In the soft-potential, late-\(\Delta x=\mathcal{O}(1)\) window scanned here, the “avalanche” **does not** convert a tiny seed into a DESI/Euclid residual. Gain remains \(\mathcal{O}(1)\). Large residuals require **large \(\sigma\)** (renaming free \(A_0\)) or a **new** regime \(b \Delta x\gg 1\) / multi-e-fold threshold dynamics **not** present in this grid — that would be a separate theory paper, not a free amplification inside the present corpus.

#### Act III decision after the scan

```
Linear amplifiers (a–c,e) → closed (audit + exclusions)
Route 2 freeze alone → numerically dead as amplifier (gain ~1)
Route 3 soft double-well → no seed→Euclid trigger (gain ~2; BAO set by σ)
Route 1 local N_eff → only card that hits target amplitude
 (by redefining the seed; price: causal-set theory)
Default publication posture → limits + exclusions + hypothesis map
 (success without detection claim)
```

---

## 7. Key equations (corpus-fixed)

Path residual:

\[
\sigma_{\rm res}(x)=A_0 e^{-\theta\Delta x}.
\]

Stationary OU variance (when \(\theta>0\)):

\[
\mathrm{Var}(X)=\frac{\sigma^2}{2\theta}.
\]

Desqueezing half-life (open systems):

\[
t_{1/2}=\frac{\ln 2}{\gamma},\qquad \gamma\leftrightarrow\theta H_0\ \text{(map, not identity of units)}.
\]

Bare seed motivation:

\[
\sigma_0\sim\frac{1}{\sqrt{N}},\quad N\sim 10^{122}.
\]

---

## 8. Bottom line

**Amplification gap (bottom line):** every channel that could make vacuum noise *visible* either fails gain by many decades, is already excluded, or replaces \(\sigma_0\) with a different microphysics. There is no free gain factor that closes \(\sim 10^{56}\) on linear amplifiers alone.
The only amplifier class worth a *serious new theory exam* is an **inflation-style freeze-out (d)** — and only after abandoning a pure Sorkin late-time seed and deriving a DE-residual map. Until then, the project’s strength is **what it already closed**, not a promised detection.

**After Routes 1–3 machine scans (§6.3):** Route 2 is numerically dead as an amplifier; Route 3 soft double-well yields gain \(\sim 2\) and BAO bounds set by \(\sigma\); Route 1 alone hits target amplitudes by redefining \(N_{\mathrm{eff}}\). Publication default remains **honest limits + exclusions + open hypothesis cards**.

---

*End of Act III amplification-gap note.*
