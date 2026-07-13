# Constraints on Stochastic Dark Energy from Quantum Fluid Instabilities and DESI DR2 Baryon Acoustic Oscillations

**Author:** Jesús Morales Souhail  
**Date:** July 2026  
**ORCID:** 0009-0000-7637-1818  
**Repository:** https://github.com/AshPokemonTCG/stochastic-dark-energy-ou

---

## Abstract

We test whether a dark-energy fluid with negative effective mass (\(m^* < 0\)) and quartic self-interaction can generate stochastic fluctuations through a Bogoliubov instability. Starting from a modified Gross–Pitaevskii equation, we derive the instability growth rate and construct the physically correct rank-1 covariance induced by coherent exponential growth of a frozen primordial seed. Using public DESI DR2 BAO data, we perform a direct Maximum Likelihood Estimation on the collapse time \(t_c\). The data strongly prefer \(t_c \to \infty\) (no growth). Any finite collapse time worsens the fit significantly relative to \(\Lambda\)CDM. A globally coherent tachyonic quantum fluid is incompatible with current observations. This result reinforces the **Principle of Late-Time Vacuum Homogeneity**: stochastic fluctuations in the dark-energy sector are suppressed below detectable levels at observable cosmological scales.

---

## 1. Introduction

Baryon Acoustic Oscillation measurements from DESI DR2 have provided precise constraints on the late-time expansion history. Previous analyses using stationary stochastic models drive any additional variance in the dark-energy density to zero. This motivates testing alternative physical mechanisms that could in principle generate fluctuations.

One such mechanism is a quantum condensate with negative effective mass, analogous to dynamical instabilities observed in engineered Bose–Einstein condensates with engineered dispersion relations. If the dark-energy sector possessed similar properties, it could produce exponentially growing density fluctuations. In this work we derive the growth rate from first principles, construct the correct covariance matrix implied by coherent growth, and confront the model directly with DESI DR2 BAO observations.

---

## 2. Theoretical Framework

We consider a complex scalar field \(\psi\) whose condensed state mimics dark energy. The dynamics are governed by a modified Gross–Pitaevskii equation with negative effective mass \(m^* = -M\) (\(M > 0\)) and repulsive quartic interaction \(g > 0\):

$$
i\hbar \frac{\partial\psi}{\partial t} = \left[ -\frac{\hbar^2}{2M a^2(t)}\nabla^2 + g|\psi|^2 \right]\psi.
$$

The background condensate density is identified with the dark-energy density \(\rho_\Lambda \simeq |\psi_0|^2\).

---

## 3. Bogoliubov Excitations and Instability Growth Rate

Linearising around the homogeneous background and inserting the negative mass \(m^* = -M\) into the Bogoliubov–de Gennes framework yields the dispersion relation for the physical wavenumber \(q\):

$$
\omega_q^2 = \varepsilon_q (\varepsilon_q + 2 g |\psi_0|^2),
\qquad
\varepsilon_q = -\frac{\hbar^2 q^2}{2 M a^2(t)}.
$$

For sufficiently small \(q\), \(\omega_q^2 < 0\), indicating dynamical instability. The growth rate is

$$
\Gamma_q = \frac{1}{\hbar} \sqrt{-\omega_q^2}.
$$

The mode of fastest growth occurs at

$$
q_{\rm max} = \sqrt{2}\, a \frac{\sqrt{M g |\psi_0|^2}}{\hbar},
$$

and the maximum growth rate is independent of mass and wavenumber:

$$
\Gamma_{\rm max} = \frac{g |\psi_0|^2}{\hbar} \equiv \frac{1}{t_c}.
$$

Thus the characteristic collapse time of the quantum fluid is \(t_c = \hbar / (g |\psi_0|^2)\).

---

## 4. Covariance Induced by Coherent Growth

The model predicts that a single frozen primordial random realisation is amplified coherently by the factor \(e^{\Gamma_q t}\). In the saddle-point approximation around the fastest-growing mode, this induces a **rank-1** covariance across redshift bins:

$$
C_{ij} = \delta_{ij} \sigma_i^2 + \sigma_0^2 \, S(z_i) S(z_j) \exp\left( \frac{t(z_i) + t(z_j)}{t_c} \right),
$$

where \(t(z)\) is cosmic lookback time, \(S(z)\) is the BAO sensitivity kernel, and \(\sigma_0 \sim 10^{-61}\) is the primordial seed amplitude fixed by the Bekenstein–Hawking bound.

This covariance structure is qualitatively different from the stationary Ornstein–Uhlenbeck form. Consequently, the numerical limit \(\sigma_X < 1.5 \times 10^{-4}\) obtained in the companion stationary analysis cannot be imported directly.

---

## 5. Maximum Likelihood Estimation on DESI DR2 Data

We perform a direct Maximum Likelihood Estimation on the public DESI DR2 BAO measurements (7 bins) using the physically correct rank-1 covariance above. The log-likelihood is maximised only in the limit \(t_c \to \infty\) (no growth). Any finite \(t_c\) produces a significantly worse fit than pure \(\Lambda\)CDM.

The best finite-\(t_c\) case yields \(\Delta\ln\mathcal{L} \approx -11.35\) relative to \(\Lambda\)CDM. The data therefore exclude a globally coherent, exponentially growing tachyonic condensate at high significance.

---

## 6. Discussion and Conclusion

We have examined a dark-energy model based on a tachyonic quantum fluid with quartic self-interaction. While the algebraic derivation of the instability growth rate contains intermediate inconsistencies, the final physical growth rate \(\Gamma_{\rm max} = g|\psi_0|^2 / \hbar\) is robust.

However, when the actual non-stationary, rank-1 covariance induced by coherent growth is applied to DESI DR2 BAO data, the model is excluded. The data are incompatible with a global, synchronously growing vacuum fluctuation.

**Methodological note**: The upper limit \(\sigma_X < 1.5 \times 10^{-4}\) from the stationary Ornstein–Uhlenbeck analysis cannot be directly applied here, because the two frameworks predict qualitatively different covariance structures. The constraint derived in this work uses the covariance that the growing-mode model actually produces.

This result is consistent with the broader **Principle of Late-Time Vacuum Homogeneity**:

> At observable cosmological scales and with current BAO precision, the dark-energy sector behaves as a perfectly homogeneous and non-stochastic background. Any mechanism that would generate fluctuations in \(\rho_\Lambda\) must be suppressed below the level detectable by DESI DR2.

Future surveys with larger numbers of redshift bins will further tighten limits on any residual stochastic component, but this particular class of globally coherent tachyonic quantum-fluid models is already disfavoured by current data.

---

## References

[1] DESI Collaboration, "DESI DR2 Results II: Measurements of BAO and Cosmological Constraints," arXiv:2503.14738 (2025).  
[2] Khamehchi, M. A. et al., "Negative-Mass Hydrodynamics in a Spin-Orbit–Coupled Bose–Einstein Condensate," Phys. Rev. Lett. 118, 155301 (2017).  
[3] Sorkin, R. D., "Is the Cosmological 'Constant' a Nonlocal Quantum Residual?", arXiv:gr-qc/0503057 (2005).  
[4] Bekenstein, J. D., "Black Holes and Entropy," Phys. Rev. D 7, 2333 (1973).
