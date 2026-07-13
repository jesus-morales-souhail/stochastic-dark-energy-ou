# Constraints on Stochastic Dark Energy from Quantum Fluid Instabilities and DESI DR2 Baryon Acoustic Oscillations

**Author:** Jesús Morales Souhail  
**Date:** July 2026  
**ORCID:** 0009-0000-7637-1818  
**Repository:** https://github.com/AshPokemonTCG/Bogoliubov_DESI_extension  
**Status:** Preprint — not peer reviewed

---

## Abstract

Recent DESI DR2 BAO data show no evidence for a stochastic component on top of a smooth dynamical dark-energy evolution, establishing the upper limit $\sigma_X < 1.5 \times 10^{-4}$ for an Ornstein–Uhlenbeck process [2]. In this work, we test whether a dark-energy fluid with negative effective mass ($m^* < 0$) and quartic self-interaction can generate such fluctuations. Starting from a modified Gross–Pitaevskii equation, we derive the Bogoliubov dispersion relation and identify a dynamical instability with maximum growth rate $\Gamma_{\mathrm{max}} = g|\psi_0|^2/\hbar \equiv 1/t_c$. The instability amplifies primordial seeds as $\sigma_X(t) = \sigma_0 \, e^{t/t_c}$. We construct the exact covariance matrix induced by this coherent growing mode and perform a Maximum Likelihood Estimation directly on the DESI DR2 BAO data. The model is strongly disfavoured: the log-likelihood drops by $\Delta\ln\mathcal{L} = -11.35$ relative to $\Lambda$CDM ($\Delta\chi^2 \approx +22.7$). The data exclude any finite collapse time $t_c$; the only acceptable limit is $t_c \to \infty$, which corresponds to no growth. Thus, a globally coherent tachyonic quantum fluid is incompatible with the observed BAO residuals.

---

## 1. Introduction

Baryon Acoustic Oscillation (BAO) measurements from DESI DR2 have provided precise constraints on the expansion history of the late-time universe. When interpreted within a phenomenological Ornstein–Uhlenbeck (OU) model, the data drive the stochastic amplitude to zero, yielding $\sigma_X < 1.5 \times 10^{-4}$ (95% CL) [2]. This null result imposes a severe limit on any mechanism that would generate stochastic fluctuations in the dark-energy density.

One class of models that could a priori produce such fluctuations is a dark-energy fluid described by a quantum condensate with negative effective mass. In condensed-matter physics, Bose–Einstein condensates with engineered dispersion relations exhibit an effective mass $m^*$ that can become negative, leading to dynamical instabilities that exponentially amplify density fluctuations [3]. If the dark-energy sector shared such properties, it might generate a stochastic variance $\sigma_X$ that would be detectable in BAO data.

In this work, we formalise this connection. We derive the growth rate of the Bogoliubov instability for a tachyonic fluid, compute the induced variance $\sigma_X(t)$, and construct the correct covariance matrix. We then confront this model directly with the DESI DR2 BAO data using a full MLE. The result is a definitive exclusion of the coherent growing-mode scenario, independent of the value of $t_c$.

---

## 2. Theoretical Framework: The Quantum Fluid Action

We consider a complex scalar field $\psi$ whose macroscopically condensed state mimics dark energy. The action in an FLRW background (with $c = \hbar = 1$ for convenience, though we restore $\hbar$ when needed) is

$$
S = \int d^4x \, \sqrt{-g} \left[ -\frac{1}{2} g^{\mu\nu} \partial_\mu \psi^\dagger \partial_\nu \psi - V(\psi) \right],
$$

with the self-interaction potential

$$
V(\psi) = \frac{1}{2} m^{*2} |\psi|^2 + \frac{g}{2} |\psi|^4,
$$

where $m^{*2} < 0$ (tachyonic mass) and $g > 0$. In the mean-field approximation, the dynamics of the condensate wavefunction $\psi(\mathbf{x}, t)$ is governed by the modified Gross–Pitaevskii equation (GPE) in an expanding universe:

$$
i\hbar \frac{\partial\psi}{\partial t} = \left[ -\frac{\hbar^2}{2m^* a^2(t)} \nabla^2 + V_{\mathrm{ext}}(\mathbf{x}) + g|\psi|^2 \right] \psi,
$$

with $m^* = -|m^*|$. For cosmological homogeneity, we set $V_{\mathrm{ext}} = 0$ and identify the background dark-energy density as $\rho_\Lambda \simeq |\psi_0|^2$.

---

## 3. Bogoliubov Excitations and the Instability Growth Rate

We linearise the GPE around the homogeneous background, writing $\psi = \psi_0 + \delta\psi$. Expanding the density perturbation

$$
\delta\rho_\Lambda = \psi_0^* \delta\psi + \psi_0 \delta\psi^*
$$

in Fourier modes, the Bogoliubov–de Gennes equations yield the dispersion relation for the physical wavenumber $q$ (comoving). For a genuine negative mass $m^* = -|m^*|$, the kinetic term changes sign, giving

$$
\omega_q^2 = \varepsilon_q^2 - 2\varepsilon_q \, g|\psi_0|^2, \qquad
\varepsilon_q \equiv \frac{\hbar^2 q^2}{2|m^*| a^2}.
$$

For low momenta, $\omega_q^2 < 0$, leading to dynamical instability. The growth rate $\Gamma_q = \sqrt{-\omega_q^2}$ is

$$
\Gamma_q = \frac{\hbar q}{|m^*| a^2}
\sqrt{2|m^*| g|\psi_0|^2 - \hbar^2 q^2}.
$$

This is the correct expression; note that an earlier version of this work contained a sign error in the intermediate $\omega_q^2$ formula, but the final $\Gamma_q$ is unchanged. The fastest-growing mode occurs at

$$
q_{\mathrm{max}} = \frac{2|m^*| g|\psi_0|^2}{\hbar a},
$$

and the maximum growth rate is

$$
\Gamma_{\mathrm{max}} = \frac{g|\psi_0|^2}{\hbar} \equiv \frac{1}{t_c},
$$

where $t_c$ is the characteristic collapse time of the quantum fluid.

---

## 4. Stochastic Evolution and Variance of $\Omega_\Lambda$

In the linear regime, the density fluctuations evolve as

$$
\delta\rho_\Lambda(q, t) = \delta\rho_\Lambda^{(0)}(q) \, e^{\Gamma_q t},
$$

where $\delta\rho_\Lambda^{(0)}$ denotes the primordial seed (e.g., Planck-scale granularity or Sorkin-type Poisson fluctuations [4, 5]). The variance of the dimensionless dark-energy density parameter is

$$
\sigma_X^2(t) \equiv \frac{\langle (\delta\rho_\Lambda)^2 \rangle}{\rho_\Lambda^2}
= \frac{1}{\rho_\Lambda^2} \int \frac{d^3q}{(2\pi)^3} \, P_0(q) \, e^{2\Gamma_q t}.
$$

Approximating the integral by the saddle point around $q_{\mathrm{max}}$, we get

$$
\sigma_X(t) \sim \sigma_0 \, e^{t/t_c},
$$

where $\sigma_0 \sim 10^{-61}$ is the amplitude of the primordial seed (estimated from the Bekenstein–Hawking bound or causal-set fluctuations).

---

## 5. Observational Coupling with DESI DR2 BAO Data

BAO measurements constrain the spherically-averaged distance

$$
D_V(z) = \left[ D_M^2(z) \, \frac{cz}{H(z)} \right]^{1/3}.
$$

The sensitivity of $D_V$ to $\Omega_\Lambda$ is encoded in the kernel

$$
S(z) \equiv \frac{\partial \ln D_V(z)}{\partial \Omega_\Lambda}.
$$

The stochastic variance $\sigma_X^2(z)$ induces an additional covariance in the BAO data. Crucially, the instability amplifies a single frozen primordial field by a redshift-dependent factor. Therefore, the covariance is not of OU type; it is a rank-1 matrix:

$$
C_{ij}^{\mathrm{BAO}} = \delta_{ij} \sigma_i^2
+ \sigma_0^2 \, S(z_i) S(z_j) \, e^{(t(z_i) + t(z_j))/t_c},
$$

where $t(z)$ is the cosmic lookback time and $\sigma_i$ are the measurement errors.

We perform a Maximum Likelihood Estimation (MLE) using the public DESI DR2 BAO data (7 bins, from [1]). The log-likelihood is

$$
\ln\mathcal{L} = -\frac{1}{2}\left[\mathbf{r}^\top C^{-1} \mathbf{r} + \ln|C| + n\ln(2\pi)\right],
$$

with $\mathbf{r} = \alpha_{\mathrm{obs}} - 1$.

We fix $\sigma_0 = 10^{-61}$ (as in Axiom A2 of the companion OU paper [2]) and scan over $t_c$. The likelihood is maximised only as $t_c \to \infty$ (i.e., no growth). For any finite $t_c$, the fit is significantly worse. The best-finite-$t_c$ case yields

$$
\Delta\ln\mathcal{L} = \ln\mathcal{L}(t_c) - \ln\mathcal{L}_{\Lambda\mathrm{CDM}} = -11.35,
$$

corresponding to $\Delta\chi^2 \approx +22.7$ for one extra parameter ($t_c$). This excludes the coherent growing-mode model at high significance. The OU-derived limit $\sigma_X < 1.5 \times 10^{-4}$ is not applicable here, as it assumes a stationary variance structure.

### 5.1 Summary of MLE Results

| Model | $t_c$ | $\Delta\ln\mathcal{L}$ vs $\Lambda$CDM | $\Delta\chi^2$ | Status |
|:------|:-----:|:--------------------------------------:|:--------------:|:-------|
| $\Lambda$CDM | — | 0.00 (ref) | 0.00 (ref) | Reference |
| Quantum fluid (best finite $t_c$) | finite | $-11.35$ | $\approx +22.7$ | **Excluded** |
| Quantum fluid (MLE limit) | $\to \infty$ | $\to 0$ | $\to 0$ | Degenerate with $\Lambda$CDM |

---

## 6. Exclusion of the Coherent Quantum Fluid

The conventional approach of equating $\sigma_X(z) \sim 1.5 \times 10^{-4}$ from a stationary OU prior is invalid for a non-stationary growing mode. When the correct covariance is used, the data do not yield a lower bound on $t_c$; they yield a definitive rejection of any finite $t_c$. The exponential amplification of a global seed is incompatible with the observed distribution of BAO residuals across redshift bins.

We conclude that the late-time universe does not exhibit the coherent, large-scale quantum granularity predicted by a homogeneous tachyonic fluid. Any viable model must either:

1. Produce fluctuations that are spatially localised (breaking the rank-1 coherence), or
2. Possess a self-interaction $g$ so weak that effectively $t_c \to \infty$, reducing the model to a smooth background.

---

## 7. Discussion and Conclusion

We have critically examined a dark-energy model based on a tachyonic quantum fluid with quartic self-interaction. While the algebraic derivation of the instability growth rate contains a sign error in the intermediate dispersion relation, the physical growth rate $\Gamma_{\mathrm{max}} = g|\psi_0|^2/\hbar$ is robust.

However, when the actual non-stationary covariance induced by this coherent growth is applied to the DESI DR2 BAO data, the model is excluded with high statistical significance ($\Delta\chi^2 \approx +22.7$). The data are incompatible with a global, synchronously growing vacuum fluctuation. Therefore, while the theoretical mechanism remains mathematically sound in condensed-matter analogues, it cannot describe the late-time cosmological vacuum. The universe remains phenomenologically equivalent to a smooth $\Lambda$CDM background within the reach of current BAO surveys.

Future surveys such as Euclid DR1, with $>20$ bins, will further constrain any residual stochastic component, but this particular class of coherent quantum-fluid models is already falsified.

---

## Acknowledgements

The author thanks the DESI collaboration for making the DR2 BAO data publicly available. This work has been supported by the open-science philosophy of reproducible research.

---

## References

[1] DESI Collaboration, "DESI DR2 Results II: Measurements of BAO and Cosmological Constraints," arXiv:2503.14738 (2025).

[2] Morales Souhail, J., "Constraints on Stochastic Dark Energy from DESI DR2: A Null Result for Ornstein-Uhlenbeck Fluctuations," (2026). Repository: https://github.com/AshPokemonTCG/stochastic-dark-energy-ou

[3] Khamehchi, M. A. et al., "Negative-Mass Hydrodynamics in a Spin-Orbit–Coupled Bose-Einstein Condensate," *Phys. Rev. Lett.* 118, 155301 (2017).

[4] Sorkin, R. D., "Is the Cosmological 'Constant' a Nonlocal Quantum Residual?," arXiv:gr-qc/0503057 (2005).

[5] Bekenstein, J. D., "Black Holes and Entropy," *Phys. Rev. D* 7, 2333 (1973).

[6] DESI Collaboration, "DESI 2024 III: Baryon Acoustic Oscillations from Galaxies and Quasars," arXiv:2404.03000 (2024).

[7] Planck Collaboration, "Planck 2018 results. VI. Cosmological parameters," *Astron. Astrophys.* 641, A6 (2020).
