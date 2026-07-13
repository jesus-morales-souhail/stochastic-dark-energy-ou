# SDiff as Fundamental Symmetry vs Emergent Feature: Implications for the Cosmological Vacuum

**Author:** Jesús Morales Souhail  
**Date:** July 2026  
**Repository:** https://github.com/jesus-morales-souhail/stochastic-dark-energy-ou

---

## 1. Introduction

The reduction of the spacetime gauge group from the full diffeomorphism group \(\mathrm{Diff}(M)\) to the volume-preserving subgroup \(\mathrm{SDiff}(M)\) is the defining feature of unimodular gravity. This reduction has direct consequences for the coupling of vacuum energy to spacetime curvature.

Two logically distinct interpretations exist:

- **Interpretation A**: \(\mathrm{SDiff}(M)\) is a *fundamental* symmetry imposed at the level of the action.
- **Interpretation B**: \(\mathrm{SDiff}(M)\) is an *emergent* low-energy symmetry arising from the quantum entanglement structure of spacetime.

Both interpretations are compatible with current data (\(\sigma_X < 1.5 \times 10^{-4}\) from DESI DR2), but they make different predictions for future observations.

---

## 2. Mathematical Distinction: Diff(M) vs SDiff(M)

| Property                      | \(\mathrm{Diff}(M)\)                          | \(\mathrm{SDiff}(M)\)                              |
|-------------------------------|-----------------------------------------------|----------------------------------------------------|
| Definition                    | All smooth reparametrizations                 | Volume-preserving diffeomorphisms only             |
| Free functions                | 4 arbitrary functions \(\xi^\mu(x)\)          | 3 arbitrary functions with \(\partial_\mu \xi^\mu = 0\) |
| Determinant of the metric     | Can vary locally (\(\sqrt{-g}\) dynamical)    | Fixed: \(\sqrt{-g} =\) constant                    |
| Consequence for vacuum energy | Local vacuum fluctuations can source curvature| Vacuum contributions of the form \(\rho g_{\mu\nu}\) are algebraically projected out |

In unimodular gravity the restriction \(\sqrt{-g} =\) constant is imposed (via Lagrange multiplier or directly in the action), reducing the gauge group to \(\mathrm{SDiff}(M)\).

---

## 3. Two Interpretations

### 3.1 Interpretation A: SDiff as a Fundamental Symmetry

\(\mathrm{SDiff}(M)\) is imposed from the beginning as a fundamental principle of the theory.

**Advantages:**
- Minimal modification of General Relativity.
- The suppression of vacuum fluctuations (\(\sigma_X \to 0\)) follows directly from the algebraic structure (trace-free projection).
- No additional dynamical mechanism is required.

**Limitations:**
- The choice of \(\mathrm{SDiff}(M)\) over \(\mathrm{Diff}(M)\) remains a postulate without deeper justification from first principles.
- Does not predict the specific value of \(\Lambda\); it only explains why it does not fluctuate.

### 3.2 Interpretation B: SDiff as an Emergent Symmetry

\(\mathrm{SDiff}(M)\) emerges at low energies as a consequence of the quantum entanglement structure of spacetime (holographic principle, Ryu-Takayanagi, etc.).

**Advantages:**
- Opens the possibility that the value of \(\Lambda\) could be determined by the entanglement structure.
- Allows for small residual fluctuations at scales below current sensitivity.
- Connects unimodular gravity with approaches to quantum gravity based on information and entanglement.

**Limitations:**
- Highly speculative. No complete quantum gravity framework currently demonstrates that \(\mathrm{SDiff}(M)\) emerges from entanglement.
- Predicts possible small stochastic signals that have not yet been observed.

---

## 4. Observational Discrimination

Current and near-future data can distinguish between the two interpretations:

| Observable                    | Fundamental SDiff                          | Emergent SDiff                              | Discriminating Power |
|-------------------------------|--------------------------------------------|---------------------------------------------|----------------------|
| Stochastic amplitude \(\sigma_X\) | Exactly zero (or below numerical floor)   | Small but non-zero residual possible        | High (Euclid DR1)    |
| Euclid DR1 (2026)             | No detection of stochastic signal          | Possible detection of small \(\sigma_X\)    | High                 |
| Primordial gravitational waves| No scalar mode                             | Possible scalar mode at high energies       | Medium               |
| DESI DR2 (current)            | Compatible                                 | Compatible                                  | None (already passed)|

**Key test:** If Euclid DR1 detects a stochastic signal with \(0 < \sigma_X \lesssim 10^{-4}\), the fundamental interpretation would be strongly disfavored. If no signal is found, both remain viable, with the fundamental option being more minimal.

---

## 5. Current Status and Outlook

From the perspective of existing data and theoretical economy:

- **Fundamental SDiff** is currently the more economical and consistent framework. It explains the observed smoothness of the vacuum (\(\sigma_X < 1.5 \times 10^{-4}\)) without additional assumptions.
- **Emergent SDiff** remains an interesting possibility but requires further development within a quantum gravity framework.

The reduction of the gauge group to volume-preserving diffeomorphisms is the only known structural mechanism that naturally protects the cosmological vacuum from local quantum fluctuations while remaining compatible with all current observations.

---

## 6. Open Questions

Three major questions emerge from this framework:

1. **Topology**: Is the total four-volume \(\mathcal{V}_4\) a fixed global observable determined from the Big Bang, or can it be influenced by quantum dynamics?
2. **Thermodynamics of the gauge group**: Is \(\mathrm{SDiff}(M)\) a fundamental symmetry or does it emerge from maximal entanglement at low energies?
3. **Structure of spacetime**: Do the DESI limits on \(\sigma_X\) rule out models with significant local metric fluctuations (e.g., spacetime foam or discrete causal sets with large variance)?

Future data from Euclid DR1 will provide the first quantitative test capable of distinguishing these possibilities.

---

*This note complements the observational analysis presented in `stochastic-dark-energy-desi-dr2.md` and `quantum-fluid-instabilities-desi-dr2.md`.*
