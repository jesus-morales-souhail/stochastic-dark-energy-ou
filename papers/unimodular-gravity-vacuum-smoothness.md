# The Smoothness of the Vacuum as a Geometric Necessity: Unimodular Gravity and the Null Result of DESI DR2

**Author:** Jesús Morales Souhail  
**Date:** July 2026  
**ORCID:** 0009-0000-7637-1818  
**Repository:** https://github.com/jesus-morales-souhail/stochastic-dark-energy-ou/  
**Status:** Preprint — not peer reviewed  
**Type:** Theoretical synthesis

---

## Abstract

Recent DESI DR2 BAO data establish a strict upper limit on stochastic fluctuations of the dark-energy density: $\sigma_X < 1.5 \times 10^{-4}$ (95% CL) [1]. This empirical fact forces a logical dichotomy: the dark energy is either (i) a pure geometric constant with no local degrees of freedom, or (ii) a dynamical scalar field so ultralight and Hubble-frozen that it is observationally indistinguishable from a constant. We then ask: what symmetry of nature can enforce this rigidity in the presence of local quantum fluctuations of matter? We examine the two candidate symmetries — conformal (Weyl) invariance and volume-preserving diffeomorphisms (SDiff). Classical Weyl invariance is broken by the conformal anomaly, which reintroduces local vacuum fluctuations. In contrast, unimodular gravity — whose gauge group is SDiff — cancels any local vacuum fluctuation $V(x) g_{\mu\nu}$ identically through the trace-free Einstein equations. We show that the cosmological constant emerges as a global integration constant,

$$
\Lambda = \frac{1}{4V_4} \int_M d^4x \, \sqrt{-g} \left(R + 8\pi G\, T^{\lambda}{}_{\lambda}\right),
$$

which is a consistency condition, not a predictive equation. The value of $\Lambda$ remains a free global parameter, but its fluctuations are rigorously zero. We conclude that the observed smoothness of the vacuum is not an accident; it is a geometric necessity that singles out unimodular gravity as the only known structure compatible with the data. Future surveys, such as Euclid DR1, will provide the decisive test.

---

## 1. Introduction

The $\Lambda$CDM model has been the standard framework of cosmology for decades. Its central feature — a cosmological constant with no dynamics and no fluctuations — has always been treated as a phenomenological placeholder, waiting to be explained by a more fundamental theory. The expectation has been that quantum gravity would eventually provide a dynamical origin for $\Lambda$ and perhaps even predict small, observable fluctuations.

That expectation has not been fulfilled.

The DESI DR2 BAO data, when analysed within a phenomenological Ornstein-Uhlenbeck (OU) framework, show no evidence for stochastic fluctuations [1]. The Maximum Likelihood Estimation drives the amplitude to zero:

$$
\sigma_X \to 0, \qquad \omega_R \to 0,
$$

establishing a conservative upper limit:

$$
\sigma_X < 1.5 \times 10^{-4}
$$

**(95% CL).**

Moreover, a specific class of models — those based on a tachyonic quantum fluid with negative effective mass — is excluded with high significance [2]:

$$
\Delta\chi^2 \approx +22.7,
$$

meaning that the data reject a coherent growing-mode scenario entirely.

These results are not a "non-discovery." They are a structural constraint on the space of possible theories. In this paper, we take this empirical limit and trace its logical consequences to their ultimate conclusion: the choice of the fundamental symmetry that protects the vacuum from local quantum fluctuations. We show that unimodular gravity — the theory whose gauge group is the group of volume-preserving diffeomorphisms (SDiff) — is the only known structure that can guarantee $\sigma_X = 0$ in a mathematically consistent way.

---

## 2. The Observational Constraint

We define the fractional fluctuation of the dark-energy density as:

$$
X(x) \equiv \frac{\delta\rho_\Lambda(x)}{\bar{\rho}_\Lambda},
$$

with variance $\sigma_X^2 \equiv \langle X^2 \rangle$. The DESI DR2 BAO data, combined with the sensitivity kernel $S(z) = \partial \ln D_V / \partial \Omega_\Lambda$, yield [1]:

$$
\sigma_X < 1.5 \times 10^{-4}
$$

**(95% CL).**

This is a direct, empirical statement: the dark-energy density is homogeneous to better than one part in ten thousand across the observable universe.

---

## 3. The Logical Dichotomy

Under the classical theory of fields and Einstein's equations, only two structures can account for this level of homogeneity:

### 3.1 Static geometry: $\Lambda$ as an integration constant

If dark energy is a pure cosmological constant, its energy-momentum tensor is

$$
T_{\mu\nu} = -\Lambda \, g_{\mu\nu}.
$$

By the Bianchi identity, $\nabla_\mu T^{\mu\nu} = 0$ implies $\partial_\mu \Lambda = 0$. There are no local degrees of freedom. Consequently,

$$
\sigma_X = 0 \quad \text{exactly}.
$$

### 3.2 A dynamical field frozen by Hubble friction

If dark energy is a scalar field $\phi$ with energy density

$$
\rho_\phi = \frac{1}{2}\dot{\phi}^2 + V(\phi),
$$

then small fluctuations $\delta\phi$ obey the wave equation in an expanding universe:

$$
\delta\ddot{\phi} + 3H\,\delta\dot{\phi} + \left(\frac{k^2}{a^2} + m_\phi^2\right)\delta\phi = 0.
$$

If $m_\phi \ll H_0$ (where $H_0 \sim 10^{-33}\ \text{eV}$), the friction term $3H\,\delta\dot{\phi}$ forces the mode into the "frozen" regime, effectively erasing inhomogeneities. The limits derived in [1] give:

$$
m_\phi \lesssim 10^{-5}\ \text{eV}, \qquad \left|\frac{V'}{V}\right| M_{\mathrm{Pl}} \lesssim 10^{-2}.
$$

In this regime, the field is observationally indistinguishable from a cosmological constant. Both alternatives lead to the same phenomenological outcome: $\sigma_X$ is zero or below detection.

---

## 4. The Two Symmetry Candidates

To distinguish which of the two alternatives is more fundamental, we must ask: what symmetry of nature prohibits local quantum fluctuations of matter from coupling to the cosmic volume element? Two candidates have been proposed.

### 4.1 Candidate A: Conformal (Weyl) invariance

Weyl invariance requires the action to be invariant under local rescalings of the metric: $g_{\mu\nu}(x) \to \Omega^2(x)\, g_{\mu\nu}(x)$. Classically, this forces the trace of the energy-momentum tensor to vanish: $T^{\mu}{}_{\mu} = 0$.

**Problem:** At the quantum level, Weyl invariance is broken by the conformal anomaly:

$$
\langle T^{\mu}{}_{\mu} \rangle_{\mathrm{anomaly}} \propto R^2 + C_{\mu\nu\rho\sigma} C^{\mu\nu\rho\sigma}.
$$

Since the local curvature $R(x)$ varies in space and time due to structure formation, the anomaly reintroduces local vacuum fluctuations. Thus, Weyl invariance cannot guarantee $\sigma_X = 0$ in a universe with matter.

### 4.2 Candidate B: Volume-preserving diffeomorphisms (SDiff)

In unimodular gravity, the gauge group is reduced from the full diffeomorphism group $\mathrm{Diff}(M)$ to the subgroup that preserves the four-volume element:

$$
\sqrt{-g} = \text{constant}.
$$

The resulting field equation is the trace-free part of Einstein's equation:

$$
R_{\mu\nu} - \frac{1}{4} R\, g_{\mu\nu} = 8\pi G \left(T_{\mu\nu} - \frac{1}{4} T^{\lambda}{}_{\lambda}\, g_{\mu\nu}\right).
$$

Now, any local vacuum fluctuation of matter has the form

$$
T_{\mu\nu}^{\mathrm{vac}}(x) = V(x)\, g_{\mu\nu}.
$$

Substituting this into the trace-free equation:

$$
V(x)\, g_{\mu\nu} - \frac{1}{4}(4V(x))\, g_{\mu\nu} = 0 \quad \text{identically}.
$$

**Conclusion:** In unimodular gravity, the trace-free projection cancels any local vacuum fluctuation, regardless of its magnitude or inhomogeneity. This is a purely algebraic cancellation, independent of the details of the quantum fields.

---

## 5. The Cosmological Constant as a Global Integration Constant

To understand how $\Lambda$ emerges in unimodular gravity, we take the divergence of the trace-free field equation. Using the Bianchi identity and the local conservation of energy-momentum, one obtains:

$$
\nabla_\nu \left(R + 8\pi G\, T^{\lambda}{}_{\lambda}\right) = 0.
$$

This implies that $R + 8\pi G\, T^{\lambda}{}_{\lambda}$ is a constant throughout spacetime:

$$
R + 8\pi G\, T^{\lambda}{}_{\lambda} = 4\Lambda,
$$

where $\Lambda$ is a global integration constant. Substituting this back into the trace-free equation recovers the standard Einstein equations:

$$
G_{\mu\nu} + \Lambda g_{\mu\nu} = 8\pi G\, T_{\mu\nu}.
$$

The value of $\Lambda$ is given by the global average:

$$
\Lambda = \frac{1}{4V_4} \int_M d^4x \, \sqrt{-g} \left(R + 8\pi G\, T^{\lambda}{}_{\lambda}\right),
$$

where $V_4 = \int d^4x \sqrt{-g}$ is the total four-volume of spacetime.

It is worth emphasizing the precise scope of this result. The unimodular structure ensures that local vacuum energy contributions of the form \(T_{\mu\nu} \propto g_{\mu\nu}\) are projected out by the trace-free Einstein equations. This provides a natural mechanism for suppressing inhomogeneous fluctuations in the vacuum energy density. 

However, this protection does not extend to the homogeneous background evolution. A time-dependent dark energy component with equation of state \(w(z) \neq -1\) remains compatible with unimodular gravity, as long as fluctuations around this background are sufficiently suppressed. In this sense, unimodular gravity constrains the *fluctuations* of the vacuum but does not, by itself, require the background energy density to be strictly constant. For a more detailed discussion on this distinction, see `sdiff-fundamental-vs-emergent.md`

---

## 6. The Meaning of the Equation: A Constraint, Not a Prediction

The equation above is not a predictive equation. It is a consistency condition — a global constraint that the universe must satisfy. If we attempted to use it to calculate $\Lambda$ deductively, we would fall into a circular argument:

1. To compute $V_4$ and the integral of $R$, we need to know the global metric $g_{\mu\nu}$.
2. To know $g_{\mu\nu}$, we must solve the Einstein equations.
3. To solve the Einstein equations, we need to know $\Lambda$ in advance.

Thus, the equation does not determine the numerical value of $\Lambda$. It only tells us that, whatever its value, it is a global constant — fixed by the initial or boundary conditions of the universe — and that local fluctuations of matter cannot affect it.

The smallness of $\Lambda$ ($\sim 10^{-122}$ in Planck units) is not a fine-tuning problem in this framework. It is a consequence of the vastness of the universe: the denominator $V_4$ is enormous, so the average of any finite quantity over it tends to zero. The discrepancy between the quantum prediction and the observed value is not an inconsistency; it is the signature that the quantum fluctuations are local, while $\Lambda$ is global.

---

## 7. The Cancellation of Local Vacuum Fluctuations ($\sigma_X = 0$)

The reason unimodular gravity is the natural candidate to explain the extreme homogeneity of the universe lies in the algebraic structure of its matter source. If we model local quantum vacuum fluctuations as a stochastic contribution of the form

$$
T_{\mu\nu}^{\mathrm{vac}}(x) = V(x)\, g_{\mu\nu},
$$

then its coupling to the metric is processed through the trace-free tensor:

$$
T_{\mu\nu}^{\mathrm{vac}} - \frac{1}{4} T^{\lambda}{}_{\lambda}{}^{\mathrm{vac}}\, g_{\mu\nu}
= V(x)\, g_{\mu\nu} - \frac{1}{4}(4V(x))\, g_{\mu\nu} \equiv 0.
$$

**Consequence:** Any local fluctuation of the quantum vacuum — regardless of its origin, energy scale, or stochastic behaviour — is identically cancelled at the algebraic level before it can curve spacetime. The metric is completely immune to this local quantum noise.

This decouples the macroscopic value of the observed cosmological constant from the catastrophic 120-order-of-magnitude problem of the quantum zero-point energy, fixing the variance of the dark-energy density exactly to:

$$
\sigma_X = 0.
$$

---

## 8. Observational Prospects: The Test of Euclid DR1

The hypothesis that unimodular gravity is the true regulatory framework behind cosmic homogeneity is highly testable in the medium term through empirical exclusion:

| Test | Prediction of Unimodular Gravity | What would falsify it |
|:-----|:--------------------------------|:---------------------|
| Euclid DR1 BAO (2026) | $\sigma_X = 0$ exactly (within experimental uncertainty) | Any statistically significant deviation ($\sigma_X > 0$) |
| Gravitational waves | Only two tensor polarisations (same as GR) | Detection of a scalar polarisation mode |
| ISW effect | No late-time ISW anomaly from dark-energy clustering | Detection of a residual ISW signal from dark-energy perturbations |

If Euclid DR1 detects any stochastic fluctuation in the dark-energy density above the noise, unimodular gravity would be ruled out. If it confirms the null result (and tightens the limit to $\sigma_X \lesssim 10^{-5}$), the case for unimodular gravity would become compelling.

---

## 9. Conclusion

We have traced the empirical limit $\sigma_X < 1.5 \times 10^{-4}$ to its ultimate logical conclusion. The only known structure that can guarantee the absolute smoothness of the vacuum in the presence of local quantum fluctuations is unimodular gravity — a theory whose gauge group is the group of volume-preserving diffeomorphisms $\mathrm{SDiff}(M)$.

This conclusion does not rest on speculative models. It follows from:

- The observational constraint that $\sigma_X$ is below $10^{-4}$.
- The logical dichotomy that only a pure constant or a frozen field can satisfy this bound.
- The mathematical fact that Weyl invariance is broken by the conformal anomaly, while SDiff cancels local vacuum fluctuations identically.

The cosmological constant $\Lambda$ emerges as a global integration constant, fixed by the initial or boundary conditions of the universe. Its numerical value remains a mystery, but its homogeneity is rigorously guaranteed: $\sigma_X = 0$. The equation that defines $\Lambda$ is a consistency condition, not a predictive equation — a constraint that the universe must satisfy, not a formula from which we can derive its value.

We therefore propose that the **Principle of Vacuum Smoothness** — the requirement that $\sigma_X = 0$ exactly — is not an ad-hoc assumption but a necessary condition for any consistent theory of gravity coupled to quantum matter. Unimodular gravity satisfies this principle naturally.

Future surveys, such as Euclid DR1 (expected 2026), may further tighten the limit on $\sigma_X$ or even detect a residual signal. But if the smoothness persists, the theoretical community must take seriously the possibility that the universe is geometrically unimodular at its foundation.

---

## Acknowledgements

The author thanks the DESI collaboration for making the DR2 BAO data publicly available. This work has been supported by the open-science philosophy of reproducible research.

---

## References

[1] Morales Souhail, J., "Constraints on Stochastic Dark Energy from DESI DR2: A Null Result for Ornstein-Uhlenbeck Fluctuations," (2026).

[2] Morales Souhail, J., "Constraints on Stochastic Dark Energy from Quantum Fluid Instabilities and DESI DR2 Baryon Acoustic Oscillations," (2026).

[3] DESI Collaboration, "DESI DR2 Results II: Measurements of BAO and Cosmological Constraints," arXiv:2503.14738 (2025).

[4] Einstein, A., "The Principle of General Relativity," 1915.

[5] Wald, R. M., *General Relativity*, University of Chicago Press (1984).

[6] van der Bij, J. J., van Dam, H., & Ng, Y. J., "The Exchange of Massless Spin-2 Particles," *Physica A* 116, 307 (1982).

[7] Ellis, J., "The Trace Anomaly and Its Implications," 1977.

[8] Ashtekar, A. et al., "BMS supertranslations and Weinberg's soft graviton theorem," *JHEP* 2015, 152 (2015).

[9] Weinberg, S., "The Cosmological Constant Problem," *Rev. Mod. Phys.* 61, 1 (1989).

[10] Carroll, S. M., "Quintessence and the Rest of the World," *Phys. Rev. Lett.* 81, 3067 (1998).
