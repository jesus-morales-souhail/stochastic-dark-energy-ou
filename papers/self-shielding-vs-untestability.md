# Self-shielding vs untestability: how the vacuum-noise mechanism hides

**Author:** Jesús Morales Souhail 
**Date:** July 2026 
**Status:** Conceptual note (honest status of the programme) 
**Repository:** https://github.com/jesus-morales-souhail/stochastic-dark-energy-ou 
**Related:** `papers/anisotropic-slip-option0.md`, `papers/sdiff-fundamental-vs-emergent.md`, `papers/fundamental-vs-emergent-vacuum-relaxation.md`, `papers/resume.txt`

---

## 1. The question

Does the SDiff / unimodular + stochastic residual mechanism **shield itself** so that it cannot be tested?

**Short answer:** It does not by itself rewrite the definition of a test. It **does** place its cleanest signatures at amplitudes and sectors where **current** probes are blind **unless** there is large amplification or a different (anisotropic) operator. That is **structural smallness**, not a free pass forever — and not a license for hype.

---

## 2. Three layers of “shielding” (precise, not mystical)

### 2.1 Conceptual layer — what DE *is* in this frame

In unimodular / SDiff-oriented narratives, vacuum pieces \(\propto g_{\mu\nu}\) are projected out of the local sourcing of curvature. Dark energy can be read as a **global / integration-constant** sector rather than a new local particle species.

**What this shields against:** naive “look for a new DE particle in the lab” searches. 
**What it does not shield against:** background \(w(z)\) tests, residual BAO smoothness, growth/slip if a non-isotropic operator exists, or external premises (e.g. spacetime discreteness via non-gravitational channels).

### 2.2 Amplitude layer — numbers from this repo

Public DESI DR2 BAO residual analysis in this repository:

\[
\sigma_X < 1.5\times 10^{-4}\quad(95\%~\mathrm{CL}).
\]

MLE drives \(\sigma_X\to 0\); stochastic extensions are not preferred (AIC penalty). 
See `papers/resume.txt`, `papers/stochastic-dark-energy-desi-dr2.md`.

**What this means:** any **isotropic residual** of the size of the bare Sorkin-scale seed (\(\sim 10^{-61}\)) is **hopelessly** below BAO. Even the **effective** bound \(10^{-4}\) is already a null for OU/QNM residuals on the kernel used.

**Shielding?** Only in the sense that **small predicted effects are hard to see**. That is ordinary experimental physics, not a conspiracy.

### 2.3 Operator layer — SDiff kills only isotropic vacuum stress

SDiff-type protection targets \(T_{\mu\nu}\propto g_{\mu\nu}\). Anisotropic stress is the structural **gap** (`papers/anisotropic-slip-option0.md`).

Gravitational slip (Maus et al. arXiv:2505.20656):

\[
\gamma=\frac{\Phi}{\Psi},\qquad \gamma=1.17\pm 0.11\ \text{(GR: }1\text{)},
\]

sensitivity \(\mathcal{O}(0.1)\) on \(|\gamma-1|\). 
Order-of-magnitude from \(\sigma_X\sim 10^{-4}\), \(f=1\):

\[
|\eta-1|\sim 2\times 10^{-4},
\]

still **orders of magnitude** below current DESI-era slip sensitivity and below naive Euclid few-percent forecasts **without amplification**.

**Shielding?** The isotropic channel is protected; the anisotropic crack is real **but amplitude-starved** — same amplification problem as desqueezing.

---

## 3. “The question the mechanism blinds”

The bad question:

> “Does the vacuum tremble?” (generic, isotropic, microscopic)

is **almost designed** to be invisible to BAO residual kernels once SDiff (or extreme freezing) is in place and amplitudes are Planck-small.

The good questions (already in the programme’s DNA):

| Question | Status |
|----------|--------|
| Is **this** OU residual preferred by DESI BAO? | **No** (null + upper limit) |
| Is **this** tachyonic growing mode preferred? | **No** (strongly disfavoured) |
| Can **this** \(\sigma_X\) produce observable slip without \(A_0\)? | **No** (Option 0 scaling) |
| Is there a **physical** amplifier from \(\sim 10^{-61}\) to \(\sim 10^{-4}\)? | **Open** — central theory problem |
| Does the framework **predict** \(H_0\) or \(S_8\) shifts? | **Not derived** in the BAO pipeline; do not claim it |

---

## 4. Paradox without cynicism

| Cynical framing | Honest framing |
|-----------------|----------------|
| “Compatible because untestable forever” | “Compatible **today** because amplitudes/operators sit below present kernels” |
| “Shield means we can never be wrong” | “Shield means **generic** isotropic noise is the wrong target; kill **specific** models” |
| “Euclid will see Sorkin seeds” | “Euclid can tighten \(\sigma_X\) and test **amplified** or **non-isotropic** extensions if defined” |

**Publishable without a detection:** upper limits, model exclusion, Option 0 honesty, improved covariance, Euclid readiness.

---

## 5. What not to do

- Do not invent tesseract / superoscillation shortcuts for cosmology. 
- Do not equate \(\ln 4 \approx \omega_R\) (archived script) with a DESI discovery — the production QNM fit drives \(\omega_R\to 0\). 
- Do not claim slip will “save” visibility without an amplification model. 
- Do not build home-made Boltzmann codes as the first move (Option 1).

---

## 6. Bottom line

The mechanism **blinds the wrong question** (generic isotropic vacuum tremor at Planck scale) and **leaves open the right ones** (specific models, tighter limits, physical amplification, carefully defined anisotropic operators). 

That is not a castle of cards immune to data: **DESI already killed the unamplified OU detection claim and the coherent tachyonic growth model.** The shield is not perfect; it is **selective**.

---

*End of note.*
