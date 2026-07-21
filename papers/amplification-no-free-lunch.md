# Amplification: no free lunch

**Author:** Jesús Morales Souhail  
**Date:** July 2026  
**Status:** Act III closure note — theory bottleneck, not a detection claim  
**Related:**  
`papers/fundamental-vs-emergent-vacuum-relaxation.md`,  
`notes/desqueezing-relaxation-vacuum-fluctuations-note.md`,  
`papers/anisotropic-slip-option0.md`,  
`papers/scale-operator-experiment-map.md`,  
`papers/resume.txt`,  
`scripts/amplifier_audit.py`

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
| Sorkin / Bekenstein seed | \(\sigma_0 \sim 10^{-61} \sim 1/\sqrt{N}\), \(N\sim 10^{122}\) | Motivational UV amplitude |
| DESI DR2 BAO residual bound | \(\sigma_X < 1.5\times 10^{-4}\) (95% CL) | Effective late-time amplitude (OU kernel) |
| Euclid-scale residual target | \(\sim 10^{-5}\) | Path-integrated residual of interest |
| Minimum kick for detection (small \(\theta\)) | \(A_0^{\min}\sim 10^{-5}\) | From \(\sigma_{\rm res}=A_0 e^{-\theta\Delta x}\) |

**Gap (order of magnitude):**

\[
\frac{A_0^{\min}}{\sigma_0}\sim \frac{10^{-5}}{10^{-61}} = 10^{56}
\qquad\text{(to}\ \sim 10^{57}\ \text{depending on \(\theta\Delta x\))}.
\]

DESI constrains **effective** amplitude after any amplification — **not** \(\sigma_0\) itself.  
Without a physical map \(\sigma_0\to A_0\), “Euclid will see vacuum noise” is **not** a theorem of this repository.

Runnable audit: `python scripts/amplifier_audit.py`.

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

### 4.1 Why (d) is special

Inflation is the **only** standard paradigm known to convert sub-horizon quantum fluctuations into large-scale classical curvature perturbations with a huge effective boost over many e-folds. Schematic amplitude after freeze-out is set by \(H/M_{\mathrm{Pl}}\) (and slow-roll factors), not by \(1/\sqrt{N_{\mathrm{BH,today}}}\).

### 4.2 Why (d) is **not** a free lunch for this repo

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
| **III** | **This note** | **No free lunch** from \(\sigma_0\) to \(10^{-5}\) |
| IV | Slip Option 0 | Right operator; still needs \(A_0\) |
| V | Method no-gos | Wrong scale/operator closed |

**Programme status if (d) fails as a late-time DE residual theory:**

- The repository remains what it already is: **honest limits, exclusions, and geometric interpretation**.  
- Amplification stays **formally open** as a theory problem.  
- That is a successful scientific posture, not a failed discovery.

---

## 6. Three *allowed* working hypotheses (not free lunches)

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
   ├─ Redefine seed (Route 1)?  → justify N_eff; re-run BAO with new prior
   ├─ Freeze-out (Route 2)?     → derive θ(k,a); note Δx_late ≪ 60
   ├─ Avalanche (Route 3)?      → derive V(X); scan BAO-safe rare-burst region
   └─ None                      → programme = honest limits + exclusions (success)
```

**Recommended default for publication of the present corpus:**  
Act III closed as **no free lunch on linear amplifiers**; Routes 1–3 listed as **the only open theory cards**, each with a price.

---

## 7. Key equations (corpus-fixed)

Path residual:

\[
\sigma_{\rm res}(x)=A_0\,e^{-\theta\Delta x}.
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

**There is no free lunch:** every channel that could make vacuum noise *visible* either fails gain by many decades, is already excluded, or replaces \(\sigma_0\) with a different microphysics.  
The only amplifier class worth a *serious new theory exam* is an **inflation-style freeze-out (d)** — and only after abandoning a pure Sorkin late-time seed and deriving a DE-residual map. Until then, the project’s strength is **what it already closed**, not a promised detection.

---

*End of Act III amplification note.*
