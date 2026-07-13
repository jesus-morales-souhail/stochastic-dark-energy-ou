# The Smoothness of the Vacuum and the Unimodular Structure of Spacetime

**Author:** Jesús Morales Souhail  
**Date:** July 2026  
**ORCID:** 0009-0000-7637-1818  
**Repository:** https://github.com/jesus-morales-souhail/stochastic-dark-energy-ou
**Status:** Preprint — not peer reviewed  
**Type:** Theoretical synthesis

---

## Abstract

Recent DESI DR2 BAO data establish a strict upper limit on stochastic fluctuations of the dark-energy density: $\sigma_X < 1.5 \times 10^{-4}$ (95% CL). This empirical fact forces a logical dichotomy: the dark energy is either (i) a pure geometric constant ($\Lambda$) with no local degrees of freedom, or (ii) a dynamical scalar field so ultralight and Hubble-frozen that it is observationally indistinguishable from a constant. We then ask: which fundamental symmetry of nature can enforce this rigidity in the presence of local quantum fluctuations of matter? We examine the two candidate symmetries: conformal (Weyl) invariance and volume-preserving diffeomorphisms (SDiff). Classical Weyl invariance is broken by the conformal anomaly, which reintroduces local vacuum fluctuations. In contrast, unimodular gravity — whose gauge group is SDiff — cancels any local vacuum fluctuation $V(x) g_{\mu\nu}$ identically through the trace-free Einstein equations. We conclude that the observed smoothness of the vacuum is not an accident; it is a geometric necessity that singles out unimodular gravity (or any theory with SDiff as a fundamental gauge symmetry) as the only known structure compatible with the data. The value of $\Lambda$ remains a free integration constant, but its fluctuations are rigorously zero.

---

## 1. Introduction

The $\Lambda$CDM model has been the standard framework of cosmology for decades. Its central feature — a cosmological constant with no dynamics and no fluctuations — has always been treated as a phenomenological placeholder, waiting to be explained by a more fundamental theory. The expectation has been that quantum gravity would eventually provide a dynamical origin for $\Lambda$ and perhaps even predict small, observable fluctuations.

That expectation has not been fulfilled.

The DESI DR2 BAO data, when analysed within a phenomenological Ornstein-Uhlenbeck (OU) framework, show no evidence for stochastic fluctuations [1]. The Maximum Likelihood Estimation drives the amplitude to zero:

$$
\sigma_X \to 0, \qquad \omega_R \to 0,
$$

establishing a conservative upper limit [1]:

$$
\sigma_X < 1.5 \times 10^{-4} \quad (95\%\ \text{CL}).
$$

Moreover, a specific class of models — those based on a tachyonic quantum fluid with negative effective mass — is excluded with high significance [2]:

$$
\Delta\chi^2 \approx +22.7,
$$

meaning that the data reject a coherent growing-mode scenario entirely.

These results are not a "non-discovery." They are a structural constraint on the space of possible theories. In this paper, we take this empirical limit and trace its logical consequences to their ultimate conclusion: the choice of the fundamental symmetry that protects the vacuum from local quantum fluctuations.

---

## 2. The Observational Constraint

We define the fractional fluctuation of the dark-energy density as:

$$
X(x) \equiv \frac{\delta\rho_\Lambda(x)}{\bar{\rho}_\Lambda},
$$

with variance $\sigma_X^2 \equiv \langle X^2 \rangle$. The DESI DR2 BAO data, combined with the sensitivity kernel $S(z) = \partial \ln D_V / \partial \Omega_\Lambda$, yield [1]:

$$
\sigma_X < 1.5 \times 10^{-4} \quad (95\%\ \text{CL}).
$$

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

## 5. Unimodular Gravity as the Logical Consequence

| Symmetry | Classical cancellation | Survives quantum anomaly? | Verdict |
|:---------|:----------------------|:--------------------------|:--------|
| Weyl invariance | Yes, trace-free | No, anomaly reintroduces fluctuations | **Rejected** |
| SDiff (unimodular) | Yes, trace-free projection | Yes, algebraic cancellation independent of matter | **Accepted** |

The mathematical structure of unimodular gravity provides the only known symmetry that can absolutely shield the cosmic volume element from local quantum noise. The cosmological constant $\Lambda$ emerges as an integration constant from the conservation law:

$$
\nabla_\mu \left(R + 8\pi G\, T^{\lambda}{}_{\lambda}\right) = 0
\quad \Longrightarrow \quad
R + 8\pi G\, T^{\lambda}{}_{\lambda} = 4\Lambda = \text{global constant}.
$$

Thus, $\Lambda$ is a single number for the entire universe. It has no local dynamics, no propagating degrees of freedom, and no fluctuations: $\sigma_X = 0$ exactly.

---

## 6. Discussion: What Remains Open

Unimodular gravity explains why $\Lambda$ does not fluctuate, but it does not predict its numerical value. The observed value $\Lambda \sim 10^{-122} M_{\mathrm{Pl}}^4$ remains a free integration constant, fixed by the initial conditions or the global topology of the universe.

This is not a flaw; it is a feature. The problem of the cosmological constant is now split into two separate questions:

1. **Why is $\Lambda$ so small?** → This is a question about the global state of the universe (initial conditions, topology, or a yet-unknown selection principle).

2. **Why does $\Lambda$ not fluctuate?** → This is answered by unimodular gravity (SDiff symmetry).

The first question is open. The second is closed by the data and by the mathematical consistency of unimodular gravity.

---

## 7. Conclusions

We have traced the empirical limit $\sigma_X < 1.5 \times 10^{-4}$ to its ultimate logical conclusion. The only known structure that can guarantee the absolute smoothness of the vacuum in the presence of local quantum fluctuations is unimodular gravity — i.e., a theory whose gauge group is the group of volume-preserving diffeomorphisms $\mathrm{SDiff}(M)$.

This conclusion does not rest on speculative models. It follows from:

- The observational constraint that $\sigma_X$ is below $10^{-4}$.
- The logical dichotomy that only a pure constant or a frozen field can satisfy this bound.
- The mathematical fact that Weyl invariance is broken by the conformal anomaly, while SDiff cancels local vacuum fluctuations identically.

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
