# Option 0: Anisotropic stress, gravitational slip, and the SDiff “gap”

**Author:** Jesús Morales Souhail 
**Date:** July 2026 
**Status:** One-page research note (literature + order-of-magnitude) — **not** a Boltzmann implementation 
**Repository:** https://github.com/jesus-morales-souhail/stochastic-dark-energy-ou 
**Related notes:** `sdiff-fundamental-vs-emergent.md`, `notes/desqueezing-relaxation-vacuum-fluctuations-note.md`, `papers/resume.txt`

---

## 0. Purpose

The BAO / OU pipeline in this repository constrains **background residual amplitude** only:

\[
\sigma_X < 1.5\times 10^{-4}\quad(95\%~\mathrm{CL}),
\]

from public DESI DR2 BAO summary statistics (`papers/resume.txt`). It does **not** compute \(\Phi\), \(\Psi\), or the gravitational slip \(\eta=\Phi/\Psi\).

**Option 0** asks, *before* any CLASS/CAMB work:

1. What published DESI-era numbers exist for slip / growth / MG?
2. Can they be compared, even roughly, to \(\sigma_X\) if a fraction \(f\) of the noise is anisotropic?
3. Is the anisotropic channel a shortcut around the **amplification** problem already discussed in the desqueezing note?

**Conclusion (preview):** the anisotropic “gap” is a *real structural crack* in SDiff protection, and \(\eta\) is the right *kind* of observable — but with \(\sigma_X\sim 10^{-4}\) it inherits the same amplitude problem as the isotropic residual. It is an **additional channel that turns on with amplification**, not a free amplification.

---

## 1. What this repository actually does (and does not)

| Present | Absent |
|--------|--------|
| CPL background + OU/QNM residual kernel on BAO distances | Boltzmann hierarchy (CLASS/CAMB) |
| Joint / profile MLE for \(\{w_0,w_a,\sigma_X\}\) | Poisson / slip equations for \(\Phi,\Psi\) |
| Euclid **BAO** mock forecasts | \(\eta(a,k)\) likelihoods |
| Unimodular / SDiff narrative notes | Full \(\eta(a,k)\) MCMC |
| `scripts/slip_bridge.py` (scaling map \(\sigma_X\to\|\gamma-1\|\)) | hi_class / MGCAMB Boltzmann |

Scripts such as `ou_bao_likelihood.py`, `joint_w0wa_sigma_desi.py`, `profile_sigma_x_desi.py`, `eos_efectiva.py`, and `euclid_*` operate on **distances** (\(D_V\), \(D_M\), \(D_H\)) obtained from \(\int dz/H(z)\). That is by design: the primary DESI DR2 product used here is **BAO-only residual smoothness**.

**Implication:** the “grieta” (anisotropic leakage past SDiff) is **quantified as amplitude-starved** via `slip_bridge.py`, not tested as a new MCMC. Option 0 remains literature + scaling.

---

## 2. Structural “gap” in SDiff protection (why \(\eta\) is the right idea)

Volume-preserving / unimodular structure projects out stress-energy contributions of the form

\[
T_{\mu\nu} = V(x) g_{\mu\nu}
\]

(the **isotropic** vacuum-like piece). That is a strong algebraic filter, not a complete eraser of all vacuum-sector noise.

If the underlying source is **granular** (Planck-scale discreteness, Sorkin–Poisson-type fluctuations), there is no *a priori* reason for the stress to be purely isotropic at every point. A shear (anisotropic) piece is generically expected. That piece is **not** of the form \(V g_{\mu\nu}\), so the trace-free / SDiff projection **need not kill it**.

**Observational handle:** gravitational slip

\[
\eta \equiv \frac{\Phi}{\Psi}.
\]

- In GR without anisotropic stress: \(\Phi=\Psi\) \(\Rightarrow\) \(\eta=1\).
- Anisotropic stress in a dark sector can yield \(\eta\neq 1\).
- Measurement: weak lensing (sensitive to \(\Phi+\Psi\)) combined with galaxy dynamics / RSD (sensitive to \(\Psi\)).

This is a **different channel** from BAO residual kernels and from tensor modes (GW): imperfect fluids with anisotropic stress can affect the two potentials without the same signature in tensors (standard lore; see e.g. imperfect-fluid / DE anisotropic-stress literature).

**Honesty (degeneracy):** \(\eta\neq 1\) does **not** uniquely prove “emergent SDiff leakage.” There is a well-known degeneracy with other modified-gravity sectors (classic review: Clifton, Ferreira, Padilla & Skordis). A detection would be a **real clue**, not an automatic identification of the mechanism.

---

## 3. Published / forecast numbers (use carefully)

Numbers below are **indicative** from public DESI-era MG / growth literature. Before any citation in a paper, open the original work and copy **definition + dataset + fiducial**.

### 3.1 Full-shape / MG-style parameters (DESI clustering era)

Analyses of DESI clustering / full-shape modelling often report deviations of the form \(\mu_0\), \(\Sigma_0\) (and related \(\gamma_0\) in some MG expansions), with GR-consistent central values near zero. Example scale of reported errors (order of magnitude from DESI 2024 MG full-shape style results in the public literature):

| Parameter (schematic) | Typical scale of constraint | Comment |
|----------------------|----------------------------|---------|
| \(\mu_0\) | \(\sim -0.2 \pm 0.1\) | ~1–2σ of zero in some combinations — check paper |
| \(\Sigma_0\) | \(\sim 0 \pm 0.06\) | largely GR-like |
| \(\gamma_0\) (MG expansion, **not** the growth index \(f=\Omega_m^\gamma\)) | \(\sim -0.16 \pm 0.07\) | ~2σ of zero in some fits — **do not** call this “no deviation” without quoting σ |

**Do not** lump these into a single slogan “GR preferred, no story.” State each with its own significance.

### 3.2 The number \(\gamma = 1.17 \pm 0.11\) — **now verified (wrong name in casual tables)**

**Paper (opened, not just the abstract blurb):** 
Maus et al., *A joint analysis of 3D clustering and galaxy × CMB-lensing cross-correlations with DESI DR1 galaxies*, 
[arXiv:2505.20656](https://arxiv.org/abs/2505.20656).

**What \(\gamma\) is in *that* paper (Appendix E):**

\[
\gamma \equiv \frac{\Phi}{\Psi} ,
\qquad
\Phi_\gamma \equiv \frac{\Phi+\Psi}{2} \equiv \Psi \frac{1+\gamma}{2} .
\]

- This is the **gravitational slip** (Newtonian potential \(\Psi\) vs lensing/Weyl-related \(\Phi\)). 
- **GR limit: \(\gamma = 1\)** (not \(0.55\)). 
- **Result:** \(\gamma = 1.17 \pm 0.11\), mild \(\sim 1.5\sigma\) from GR; authors note it may be a **projection effect** (Appendix E). 
- Probe: RSD (3D clustering) + galaxy×CMB-lensing (Planck PR4 / ACT DR6) with DESI DR1 BGS+LRG.

**What it is *not*:**

- **Not** the growth-rate index \(f(z)=\Omega_m(z)^\gamma\) of structure growth (that book-value is \(\gamma\simeq 0.55\) in GR/ΛCDM). 
- Confusing the two was the earlier table error: a slip \(\gamma\sim 1.17\) is *near* GR’s \(1.0\); a growth index \(\gamma\sim 1.17\) would be a catastrophe. Once the definition is read, the number is **internally consistent**.

**For our Option 0 note:** this is actually *better* than a random growth-index: it is a **published DESI-era slip constraint of the type we care about** (RSD + lensing). 
With \(\sigma(\gamma)\sim 0.11\), current sensitivity to \(|\eta-1|\) (or \(|\gamma-1|\)) is still \(\mathcal{O}(0.1)\), far above the \(\sim 10^{-4}\) scale of the \(f=1\), \(\sigma_X\sim 10^{-4}\) estimate below.

### 3.3 Forecasts on slip \(\eta\) (or equivalent \(\gamma=\Phi/\Psi\))

Model-independent forecasts combining Euclid-like and DESI-like surveys for anisotropic stress / slip typically quote **percent-level** reach on \(\eta\) (few percent, not \(10^{-4}\)). Exact \(\sigma(\eta)\) depends on survey assumptions; use as \(\mathcal{O}(0.03\text{–}0.05)\) order-of-magnitude, not a hard law.

### 3.4 Review to read cover-to-cover (Option 0 reading list)

1. **Dark Energy After DESI DR2: Observational Status, Reconstructions, and Physical Models** (2026 review literature) — section on \(\eta(a,k)=\Phi/\Psi\), growth index, RSD + weak lensing. 
2. Clifton, Ferreira, Padilla & Skordis — modified gravity review (degeneracies: anisotropic stress vs other MG). 
3. DESI collaboration papers on full-shape / MG parameters (\(\mu,\Sigma\)) with DR1/DR2-era clustering. 
4. Euclid + DESI-like forecasts for anisotropic stress / slip. 
5. Plaza, León & Kraiselburd (2025) — unimodular + DESI DR2 (different science goal: holographic DE / \(H_0\)); useful for **pipeline craft**, not for our \(\sigma_X\) map.

---

## 4. Order-of-magnitude map \(\sigma_X \rightarrow |\eta-1|\)

**Runnable bridge:** `python scripts/slip_bridge.py` 
(sub-horizon anisotropy equation; not a Boltzmann code).

Let \(\varepsilon\in[0,1]\) (also written \(f\) below) be the fraction of the DE-sector residual stress that is **anisotropic** (shear), and let \(\sigma_X\) be the fractional residual amplitude constrained by BAO. With \(\pi_T=\varepsilon \sigma_X \rho_X\) and \(\delta_m\sim 1\):

\[
|\gamma-1|
 = 
2 \varepsilon \sigma_X \frac{\rho_X}{\rho_m |\delta_m|}.
\]

**Machine check (best case \(\varepsilon=1\), \(\sigma_X=1.5\times 10^{-4}\), \(\delta_m=1\)):**

| \(z\) | \(\rho_X/\rho_m\) | max \(\|\gamma-1\|\) | gap vs Maus \(\sim 0.17\) | gap vs Sakr floor \(\sim 0.05\) |
|-------|-------------------|----------------------|---------------------------|--------------------------------|
| 0.5 | 0.644 | \(1.93\times 10^{-4}\) | \(\sim 880\times\) | \(\sim 260\times\) |
| 1.0 | 0.272 | \(8.16\times 10^{-5}\) | \(\sim 2\times 10^{3}\) | \(\sim 610\times\) |
| 1.5 | 0.139 | \(4.18\times 10^{-5}\) | \(\sim 4\times 10^{3}\) | \(\sim 1.2\times 10^{3}\) |

Compare:

| Instrument / status | Rough sensitivity to \(\eta-1\) | vs \(\sim 10^{-4}\) (best case) |
|---------------------|----------------------------------|----------------------------------|
| Current DESI-era slip (Maus \(\gamma=1.17\pm 0.11\)) | \(\mathcal{O}(0.1)\) | **far above** prediction |
| Euclid-like forecasts (Sakr et al.) | \(\mathcal{O}(0.03\text{–}0.05)\) | **still above** by \(\sim 10^{2}\)–\(10^{3}\) |

**Interpretation:** even the most optimistic anisotropic fraction does **not** make the gap visible at current DESI precision, nor automatically at Euclid, **unless** there is **amplification** of the bare residual amplitude (or a shear sector **decoupled** from the BAO-bounded isotropic \(\sigma_X\)).

That word — **amplification** — is the same one that appears in Act III (`amplification-gap.md`) when mapping Sorkin-scale seeds (\(\sim 10^{-61}\)) up to BAO-relevant \(\sigma_X\) (factors of order \(A_0/\sigma_0\sim 10^{56}\)). The anisotropic channel **does not evade** that problem; it **inherits** it.

| \(\varepsilon\) | \(\|\gamma-1\|\) at \(z=0.5\), \(\sigma_X=1.5\times 10^{-4}\) |
|-----------------|------------------------------------------------------------------|
| 1 | \(1.93\times 10^{-4}\) |
| 0.1 | \(1.93\times 10^{-5}\) |
| 0.01 | \(1.93\times 10^{-6}\) |

---

## 5. What Option 0 does *and* does not conclude

**Does conclude:**

1. DESI **BAO residual** bounds (\(\sigma_X\)) do **not** by themselves kill an anisotropic gap — they bound a different operator (isotropic residual covariance on distances). 
2. Published growth/MG analyses with DESI-like data are **GR-compatible at the \(\mathcal{O}(0.1)\) level**, which is still **too coarse** to see \(|\eta-1|\sim 10^{-4}\). 
3. \(\eta\) remains the **conceptually correct** crack to watch **if** a large amplification of vacuum residual ever becomes physical. 
4. **Do not** build homemade Boltzmann codes (Option 1) as the next step. If code is needed later: **hi_class / MGCAMB / EFTCAMB** (Option 2), with a **pre-specified** target amplitude from the table above.

**Does not conclude:**

- That Euclid will “easily” see the gap without amplification. 
- That \(\eta\neq 1\) would uniquely select emergent SDiff. 
- That any \(\gamma\sim 1.17\) number can be cited without checking the paper’s definition.

---

## 6. Checklist (one afternoon → one week)

- [ ] Open the DESI DR2 dark-energy review section on \(\eta\) / \(\gamma\) / RSD+lensing; copy **exact** definitions and tables. 
- [ ] For each MG/growth paper: write one line “parameter = …, dataset = …, GR limit = …, \(n\sigma\) = …”. 
- [ ] Fill the \(f\) table with your own \(\sigma_X^{\mathrm{eff}}\) (use \(1.5\times 10^{-4}\) and \(5\times 10^{-5}\) as bookends from `resume.txt`). 
- [ ] Decision gate: if no path to \(|\eta-1|\gtrsim 10^{-2}\) without huge \(A_0\), **pause** Boltzmann work; focus on amplification physics (desqueezing / open systems / non-Sorkin seeds). 
- [ ] If a path exists (e.g. large \(f\) + amplification model), only then open **hi_class** or **MGCAMB** with a single anisotropic-stress parameter and DESI lensing+clustering likelihoods already in the community.

---

## 7. Closing sentence (for the project narrative)

**This anisotropic / gravitational-slip channel is subject to the same amplification problem as the rest of the project; it is not a route that avoids amplification, but an additional route that would activate together with it.**

That statement links this Option 0 note to the desqueezing / Sorkin-seed discussion and prevents “easy slip with Euclid” from becoming a new form of hype.

---

## 8. References (minimal, to open first)

- Morales Souhail, J. — this repository: `papers/resume.txt`, `papers/stochastic-dark-energy-desi-dr2.md`, `papers/sdiff-fundamental-vs-emergent.md`. 
- DESI Collaboration BAO DR2 public products (arXiv:2503.14738 and companion papers). 
- Clifton, T., Ferreira, P. G., Padilla, A., & Skordis, C. — modified gravity review (degeneracies). 
- DESI full-shape / MG analyses reporting \(\mu,\Sigma\) (and any \(\gamma_0\) in their MG basis). 
- Euclid + DESI-like forecasts for anisotropic stress / \(\eta\). 
- Plaza, León & Kraiselburd (2025) — unimodular + DESI DR2 (different science goal; pipeline craft only).

---

**Runnable check:** `python scripts/slip_bridge.py`

*End of Option 0 note. No Boltzmann code was modified or added.*
