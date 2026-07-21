# No-go note: Tesseract phase cavities, superoscillations, and the Born rule

**Author:** Jesús Morales Souhail 
**Date:** July 2026 
**Status:** Negative result / conceptual no-go (optics + probability) 
**Repository:** https://github.com/jesus-morales-souhail/stochastic-dark-energy-ou 
**Related:** `papers/anisotropic-slip-option0.md`, archived `quantum_information_cosmos.py` (local_archive only)

---

## 1. Claim under examination

A proposed scheme asserts that:

1. Optical wave equations are time-reversal symmetric, so one can inject a **phase-conjugated / “negative diffraction”** wavefront.
2. The **8 cubic cells of a 4D hypercube (tesseract)** supply a phase cavity / Coxeter \(B_4\) symmetry that controls \((k_x,k_y,k_z,k_t)\) better than ordinary optics.
3. The result is a **superoscillatory focus** smaller than the Airy disk, with the physical medium “undoing” pre-distortion.
4. With a “sufficiently perfect phase,” a **single photon** will land in the sub-diffraction core rather than the side lobes.

Items (1) and (3) partially touch real physics. Items (2) and (4) do not survive as stated.

---

## 2. What is real

### 2.1 Time-reversal / phase conjugation

Maxwell’s equations in a passive, linear, lossless medium are time-reversal invariant in the standard sense. **Optical phase conjugation** and adaptive optics are established techniques: pre-distort a wavefront so free-space propagation partially undoes aberration. This does **not** require 4D geometry.

### 2.2 Superoscillations (Berry & Popescu and related work)

A band-limited function can, on a finite interval, oscillate faster than its highest Fourier component. This is **real** and experimentally demonstrated. It does **not** violate Maxwell: the spectrum remains band-limited.

### 2.3 The energy / probability tax

A standard structural fact of superoscillations: the more extreme the local superoscillation (tighter sub-wavelength focus), the **smaller** the fraction of \(L^2\)-norm in the core relative to the side lobes. Schematically:

\[
\frac{E_{\mathrm{core}}}{E_{\mathrm{total}}} \ll 1
\qquad\text{(often exponentially small in the “how super” parameter).}
\]

For a single photon, the Born rule identifies \(|\psi|^2\) with detection probability. If the classical energy fraction in the core is \(\varepsilon\), the single-shot probability of detection in the core is of order \(\varepsilon\), not \(\sim 1\).

A toy numerical demonstration lives in:

`scripts/superoscillation_energy_cost_demo.py`

---

## 3. Where the argument breaks

### 3.1 Born is not a technical bug

**False conclusion:** “0.01% in the core means we almost have it; perfect the phase further and the photon will go to the centre.”

**Correct conclusion:** A **perfect** phase for a given superoscillatory design **already produces** that \(\varepsilon : (1-\varepsilon)\) split. Further “perfection” of the same design does not convert \(\varepsilon\to 1\). Asking the photon to “always hit the centre” is asking the Born measure to stop being a probability measure.

Numbers like “0.01% / 99.99%” in informal texts are order-of-magnitude illustrations, not a free parameter one dials with geometry. The qualitative statement stands: **side lobes carry almost all the norm**.

### 3.2 Tesseract / 8 cubes / \(B_4\): geometry ≠ undeclared control power

**True:** A tesseract has **8 cubic cells** as 3D boundary facets. The hyperoctahedral / Coxeter group \(B_4\) is a real symmetry group of the 4-cube (\(|B_4|=4!\cdot 2^4=384\) signed permutations of \(\mathbb{R}^4\)).

**Sealed count (optical phase space, first order):** 
On the paraxial phase space \((x,p_x,y,p_y)\) with symplectic form
\[
J=\mathrm{blockdiag}(J_2,J_2),\qquad
J_2=\begin{pmatrix}0&1\\-1&0\end{pmatrix},
\]
exactly **32** of the 384 signed permutations satisfy \(MJM^{\mathsf T}=J\) 
(check: \(4\times 4\times 2=32\) — symplectic rotations in each plane plus exchange of the two planes with compatible signs). 

Those 32 generate a joint **commutant of dimension 1** inside \(\mathfrak{sp}(4,\mathbb{R})\) (\(\dim=10\)). 
Runnable verification: `scripts/b4_symplectic_count.py`.

**Reading of the seal:** as a *symmetry constraint*, \(B_4\) **removes** 9 of 10 first-order optical freedoms; it is **never** an amplifier of channels. The “8” of the 4-cube is **polytope boundary combinatorics**, not the rank of a Maxwell-derived operator.

**What the device actually is (Maxwell → ABCD):** 
Paraxial Helmholtz \(\to\) Fresnel propagators \(P_d\) (metaplectic) + thin phase screens \(L_\Phi\) (e.g. SLM). 
Composition \(M=P_{d_n}L_{\Phi_{n-1}}\cdots L_{\Phi_1}P_{d_1}\). 
If every \(\Phi\) is quadratic, \(M\in\mathrm{Sp}(4,\mathbb{R})\) (Collins kernel); image condition \(B_{\mathrm{tot}}=0\) recovers Descartes. 
Controllable DOF of a phase-only SLM: \(N_{\mathrm{pix}}\) (e.g. \(1920\times 1080\)) times spectral bins if used — **hardware**, not 8/16/24/32 from the 4-cube.

**Not shown in tesseract proposals:**

- A map from “8 projected cubes” \(\to\) a physical field \(E(\mathbf{r},t)\) obeying Maxwell that *extends* \(\mathrm{Sp}(4,\mathbb{R})\) or \(N_{\mathrm{pix}}\). 
- Any advantage over standard phase conjugation / holography / Fourier optics.

This is the same **pattern of undeclared power** as in the archived script `quantum_information_cosmos.py`: take a clean mathematical object (\(\ln 4\), \(B_4\)) and **declare** it sets a physical frequency or optical control law without derivation.

### 3.3 Decoherence is secondary to Born

Thermal index fluctuations and inelastic scattering do destroy phase conjugation in real media. That is true — and **orthogonal** to the main no-go: even in a perfect deterministic medium with perfect phase, superoscillation still parks almost all probability in the side lobes.

---

## 4. Formal no-gos

### 4.1 Optics + quantum detection (Born)

**Assumptions**

1. Propagation is linear, band-limited (effective \(|k|\le k_{\max}\)).
2. Detection of one photon is described by a POVM / Born rule from a normalized field mode \(\psi\).
3. The “sub-Airy core” is a region \(C\) with \(\int_C |\psi|^2 = \varepsilon \ll 1\) for any superoscillatory design that meaningfully beats the diffraction scale on \(C\).

**Conclusion**

No choice of deterministic pre-phase (tesseract-inspired or otherwise) that yields such a \(\psi\) can make the single-shot detection probability in \(C\) approach 1 without either:

- abandoning band-limitation / Maxwell in free space, or 
- post-selecting / discarding trials (which reintroduces the energy cost as failed shots), or 
- changing the definition of “success” (e.g. classical multi-photon intensity peaks without single-photon certainty).

### 4.2 First-order device symmetry (\(B_4\) vs \(\mathrm{Sp}(4,\mathbb{R})\))

**Assumptions:** paraxial phase-space maps in \(\mathrm{Sp}(4,\mathbb{R})\); candidate “tesseract symmetry” = signed-permutation action of \(B_4\) on \(\mathbb{R}^4\).

**Facts (machine-checked):** 32 of 384 elements of \(B_4\) are optical-symplectic; their commutant in \(\mathfrak{sp}(4,\mathbb{R})\) has dimension **1**.

**Conclusion:** Imposing \(B_4\) as a symmetry **restricts** first-order devices (\(10\to 1\) free directions in the Lie algebra, up to that commutant). It cannot be a derivation of **extra** independent optical channels equal to the number of cells/vertices of the 4-cube. Any claim “8 cells → 8 channels” confuses **polytope combinatorics** with **operator rank**.

---

## 5. What remains scientifically useful

| Idea | Status |
|------|--------|
| Phase conjugation / adaptive optics | Real engineering |
| Superoscillations | Real math/physics; energy tax real |
| Sub-diffraction classical intensity features | Possible with high side-lobe cost |
| Deterministic single-photon sub-Airy “hit every time” | **No-go** under Born + band-limit |
| Tesseract as superior phase engine | **Unsupported claim** |

A clean PC-scale project (no ad hoc free-parameter tuning): run `scripts/superoscillation_energy_cost_demo.py` and report \(E_{\mathrm{core}}/E_{\mathrm{total}}\) vs bandwidth — **confirming** the tax, not bypassing it.

---

## 6. Relation to the cosmology programme

The cosmology repo’s **publishable** content is BAO residual limits, model exclusion (e.g. tachyonic fluid), and Option 0 on anisotropic slip (which **inherits** the amplification problem). Optical tesseract schemes are **out of scope** and should not be mixed into DESI/SDiff papers.

The same **methodological** standard applies: no elegant object (\(\ln 4\), \(B_4\), SDiff) gets observational power without a derived map to a measurable and an honest amplitude estimate.

---

*End of no-go note.*
