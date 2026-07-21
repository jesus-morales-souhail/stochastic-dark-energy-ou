# Wavefront shaping vs the OU residual: what transfers and what does not

**Author:** Jesús Morales Souhail  
**Date:** July 2026  
**Status:** Boundary note (optics of disordered media ≠ BAO OU model)  
**Related:**  
`papers/self-shielding-triple-barrier.md`,  
`papers/optics-ou-analogies-and-limits.md`,  
`papers/anisotropic-slip-option0.md`,  
`papers/scale-operator-experiment-map.md`,  
`papers/stochastic-dark-energy-desi-dr2.md`

---

## 1. Verdict first

| Statement | Verdict |
|-----------|---------|
| Wavefront shaping, optical time reversal, and Wigner–Smith operators are real lab physics | **True** |
| The OU process in **this repo** is the same object as a lab transmission matrix \(T\) of a disordered slab | **False** |
| Characterizing a **fixed** \(T\) lets disorder help focus light | **True** (lab, quasi-static medium) |
| Therefore we can “use the vacuum’s disorder as a tool” by conjugating phase through the OU model | **Not supported** |
| Expanding Maxwell on a stochastic metric is a possible *formal* calculation | **Possible in principle** |
| That calculation is already “well-defined and closed” for DESI \(\sigma_X\) and GRB data | **False** — wrong identification of fields and scales |
| “Principal modes immune to DE fluctuations” as an observational channel today | **Speculation**, not a derived pipeline |

**Do not put wavefront-shaping-of-the-vacuum as a core result of this repository.**  
**Do put:** honest scale/operator separation and, if ever, a carefully scoped GR+Maxwell calculation with its own amplitude assumptions.

---

## 2. What the repo’s OU process actually is

From the production analysis (`ou_bao_*`, `papers/stochastic-dark-energy-desi-dr2.md`):

\[
\mathrm{d}X = -\theta X\,\mathrm{d}x + \sigma\,\mathrm{d}W_x,
\qquad x=\ln a,
\qquad
X \sim \frac{\delta\rho_\Lambda}{\rho_{\Lambda,0}},
\]

with residual BAO covariance

\[
(C_{\mathrm{OU}})_{ij}
= S(z_i)S(z_j)\,\frac{\sigma_X^2}{2\theta}\,e^{-\theta|x_i-x_j|}
\quad\text{(or equivalent normalization conventions).}
\]

| Property | Lab disordered medium \(T\) | Repo OU |
|----------|----------------------------|---------|
| State variable | EM field on input/output channels | Scalar residual \(X(x)\) on **cosmic expansion history** |
| Correlation | Spatial (and maybe temporal) structure of scatterer | In **e-fold time** \(x=\ln a\) between BAO redshifts |
| Measurement | Interferometry / SLM / cameras | Distance ratios \(D_V,D_M,D_H\) |
| Control | Shape input wavefront | **No** control of “input modes” of the cosmos |
| Time reversal | Send conjugated field **back through same realization** of \(T\) | No second pass of the same metric realization available for GRBs |

Identifying

\[
\langle\delta g_{\mu\nu}(x)\,\delta g_{\alpha\beta}(x')\rangle
\stackrel{?}{=}
\sigma_X^2\,e^{-|x-x'|/\xi}
\times\text{(polarization structure)}
\]

with the BAO residual \(\sigma_X\) is an **extra theoretical step** that the BAO pipeline **does not** perform. The published bound \(\sigma_X < 1.5\times 10^{-4}\) is **not** a free amplitude for arbitrary metric shear felt by photons.

---

## 3. Three optical ideas — transfer filter

### 3.1 Wavefront shaping (complex media)

**Lab truth:** For a **quasi-static** linear medium, the input–output map is a matrix \(T\). Measuring \(T\) (or optimizing intensity without full \(T\)) allows focusing through opacity.

**Cosmology transfer:**  
- There is no operational “measure \(T\) of the vacuum” with an SLM.  
- GRB photons make **one way** trips through an expanding universe; you do not reconfigure the “input” at the source with a measured \(T\).  
- Fluctuations of DE residual on BAO are constrained **after** integration against \(S(z)\), not as a scattering matrix on Fourier modes of \(E\).

**Useful analogy only:** “If the disorder is unknown and uncontrolled, you cannot shape through it.” That supports **characterization**, not a free focusing tool for DE.

### 3.2 Optical time reversal

**Lab truth:** For a **time-independent** passive linear medium, phase conjugation undoes multipath.

**Cosmology transfer:**  
- A photon from a GRB does not reverse through the **same** stochastic realization of \(\delta g_{\mu\nu}(x)\).  
- The background expands; the “medium” is not a fixed lab slab.  
- Stationarity of an OU process in \(x=\ln a\) is **not** the same as time-reversal symmetry of a fixed scatterer for EM waves.

**Verdict:** Time reversal is **not** “valid for the repo” in the lab sense. Gaussianity ≠ optical TR.

### 3.3 Wigner–Smith / principal modes

**Lab truth:** \(Q=-i T^{-1}\partial_\omega T\) defines delay operators; eigenmodes can minimize delay spread in multi-path systems.

**Cosmology transfer:**  
- Requires a well-defined frequency-dependent EM transmission operator through the medium.  
- That operator is **not** the BAO OU kernel.  
- Claiming “modes immune to \(\sigma_X\)” without deriving \(T(\omega)\) from Maxwell + metric is **undeclared power**.

---

## 4. On the proposed perturbative \(T\)

A legitimate research program (outside the current BAO code) would look like:

1. Specify the **metric** model: which components of \(h_{\mu\nu}\) fluctuate, gauge, correlation length in **physical** coordinates (not just \(x=\ln a\)).  
2. Propagate Maxwell (or geometric optics / WKB) on FLRW + \(h_{\mu\nu}\).  
3. Expand the scattering / phase operator in the amplitude of \(h\).  
4. **Map** that amplitude to observables (arrival-time variance, polarization rotation, etc.).  
5. **Separately** argue how (if at all) that amplitude relates to the BAO residual \(\sigma_X\).

Until (1)–(5) are done, the expansion

\[
T(k,k')=\delta(k-k')+\sigma_X\,\tilde T_1+\sigma_X^2\,\tilde T_2+\cdots
\]

with \(\sigma_X\) from DESI BAO is **not** a theorem of this repository.

**SDiff remark (partially right):** Maxwell in metric form involves \(\sqrt{-g}\) and \(g^{\mu\nu}\). Trace / conformal pieces and shear enter **differently**; isotropic vacuum-like stress is special in gravity, but **photon propagation** still needs an explicit calculation — one cannot simply quote “SDiff ⇒ photons see nothing” or “shear ⇒ GRB phase delay” without the expansion.

**Amplitude remark:** Even if a linear map \(\delta\phi\sim \mathcal{O}(\sigma_h)\) exists, using \(\sigma_h\sim 10^{-4}\) is already optimistic relative to BAO residual meaning; Sorkin-scale seeds are far smaller. Order-of-magnitude GRB constraints typically kill large Lorentz-violating or large stochastic refractive models; they do **not** automatically test the BAO OU residual without a bridge.

---

## 5. Reframing that *is* allowed

| Bad reframe | Better reframe |
|-------------|----------------|
| “Chaos focuses DE for us like an SLM” | “Uncontrolled stochasticity without a measured \(T\) cannot be inverted” |
| “OU is reversible ⇒ time-reverse the cosmos” | “OU is a Gaussian model for residual correlations in \(x\); TR of EM is a different structure” |
| “Principal modes are immune channels to observe DE” | “If we ever derive a delay operator from \(h_{\mu\nu}\), its eigenmodes are a calculational tool — not an immunity badge” |
| “Grieta = anisotropic \(T\) measured by GRB” | “Grieta = anisotropic stress may source slip \(\Phi/\Psi\) (gravity sector); photon channel needs its own map and amplitude” |

The **good** methodological move from disordered-media optics is:

> **Characterise the actual operator that couples to your probe.**  
> Do not fight a cartoon enemy; do not claim control you do not have.

That supports Option 0 (slip for anisotropic stress) and BAO residual kernels — **not** vacuum wavefront shaping.

---

## 6. Answers to the four open questions (as posed)

| Question | Answer |
|----------|--------|
| Analogue of \(T\) for the quantum vacuum? | Only after defining a **metric** noise model + Maxwell map. Not automatic from BAO OU. |
| Does repo OU give analytic \(T\)? | **No.** It gives \(C_{ij}\) for distance residuals. |
| Is optical TR valid for the repo? | **Not in the lab sense.** No second pass through the same realization. |
| Where is the grieta? | In **gravity**: isotropic \(V g_{\mu\nu}\) vs shear → slip (Option 0). Photon \(T\) is a **different** project with its own honesty requirements. |

---

## 7. What goes where

| Item | Git `stochastic-dark-energy-ou` | Separate optics sandbox / PC |
|------|--------------------------------|------------------------------|
| BAO OU / \(\sigma_X\) / model kills | **Yes** | — |
| Option 0 slip + data pack | **Yes** | — |
| No-go tesseract / triple barrier | **Yes** | — |
| This boundary note | **Yes** | — |
| Full Maxwell + stochastic metric code | Only if scoped as **new theory module** with explicit assumptions | Prefer separate folder first |
| Lab wavefront shaping / SLM experiments | **No** | Optics lab / other project |
| “Focus vacuum with conjugated \(E_{\mathrm{in}}\)” as a result | **No** | Fantasy until (1)–(5) exist |

---

## 8. Recommended next step (if any calculation)

**Not:** “derive \(T\) for the repo” as if \(T\) already exists.  

**Instead, if pursued at all:**

1. Write a one-page **specification**: which \(h_{\mu\nu}\), what correlator in **physical** \((t,\mathbf{x})\), what observable (e.g. \(\langle(\Delta t)^2\rangle\) for two energies).  
2. Compute geometric-optics delay to linear order in \(h\).  
3. Insert an **independent** amplitude prior (not silently equal to BAO \(\sigma_X\)).  
4. Compare order-of-magnitude to published GRB / LIV bounds.  
5. Decide if the bridge to BAO \(\sigma_X\) is even possible.

Until then: **do not** open a “principal transmission modes of dark energy” work package.

---

## 9. One-sentence closure

**Disordered-media optics teaches characterization and control of a fixed linear map \(T\); the OU model constrains a scalar residual on the expansion history — different object, different operator, different experiment. Using chaos as a tool requires measuring that map; we have not measured a vacuum \(T\), and the repo does not define one.**

---

*End of boundary note.*
