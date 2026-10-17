# The Smoothness of the Vacuum as a Geometric Necessity: Unimodular Gravity and the Null Result of DESI DR2

**Author:** Jesús Morales Souhail
**Date:** July 2026
**ORCID:** 0009-0000-7637-1818
**Repository:** https://github.com/jesus-morales-souhail/stochastic-dark-energy-ou/
**Status:** Preprint
**Type:** Theoretical synthesis

---

## Abstract

Recent DESI DR2 BAO data give a 95% CL profile upper limit (over $\theta$) on stochastic fluctuations of the dark-energy density: $\sigma_X < 2.6\times 10^{-2}$ [1]. I take this empirical fact as forcing a logical dichotomy: dark energy is either (i) a pure geometric constant with no local degrees of freedom, or (ii) a dynamical scalar field so ultralight and Hubble-frozen that it is observationally indistinguishable from a constant. I then ask: what symmetry of nature can enforce this rigidity in the presence of local quantum fluctuations of matter? I examine the two candidate symmetries — conformal (Weyl) invariance and volume-preserving diffeomorphisms (SDiff). Classical Weyl invariance is broken by the conformal anomaly, which reintroduces local vacuum fluctuations. In contrast, unimodular gravity — whose gauge group is SDiff — cancels any local vacuum fluctuation $V(x) g_{\mu\nu}$ identically through the trace-free Einstein equations. I show that the cosmological constant emerges as a global integration constant,


$$
\Lambda = \frac{1}{4V_4} \int_M d^4x \sqrt{-g} \left(R + 8\pi G T^{\lambda}{}_{\lambda}\right),
$$


which is a consistency condition, not a predictive equation. The value of $\Lambda$ remains a free global parameter, but its local fluctuations of the form $V(x)g_{\mu\nu}$ are cancelled identically. The same framework with causal-set discreteness predicts order-one fluctuations of the global $\Lambda$ (Everpresent Λ, Ahmed et al. 2004); which case holds is decided by data (Section 7). Future surveys, such as Euclid DR1, will provide the decisive test.

---

## 1. Introduction

The $\Lambda$CDM model has been the standard framework of cosmology for decades. Its central feature — a cosmological constant with no dynamics and no fluctuations — has always been treated as a phenomenological placeholder, waiting for a more fundamental theory. The usual expectation has been that quantum gravity would eventually give $\Lambda$ a dynamical origin and perhaps even predict small, observable fluctuations.

That expectation has not been fulfilled.

When I analyse the DESI DR2 BAO data within a phenomenological Ornstein–Uhlenbeck (OU) framework, I find no evidence for stochastic fluctuations [1]. The best fit drives the amplitude to zero:


$$
\sigma_X \to 0, \qquad \omega_R \to 0,
$$


and the profile likelihood over $\theta$ gives the 95% CL upper limit:


$$
\sigma_X < 2.6\times 10^{-2}
$$


**(95% CL).**

I also exclude a specific class of models — those based on a tachyonic quantum fluid with negative effective mass — with high significance [2]:


$$
\Delta\chi^2 \approx +25
$$
so the data reject a coherent growing-mode scenario entirely.

I do not read these results as a mere “non-discovery.” They are a structural constraint on the space of possible theories. In this paper I take the empirical limit and follow its logical consequences to the end: the choice of fundamental symmetry that protects the vacuum from local quantum fluctuations. I show that unimodular gravity — the theory whose gauge group is the group of volume-preserving diffeomorphisms (SDiff) — is the only known structure that can guarantee $\sigma_X = 0$ in a mathematically consistent way.

---

## 2. The Observational Constraint

I define the dark-energy density fluctuation in units of today's critical density as:


$$
X(x) \equiv \frac{\delta\rho_{\rm DE}(x)}{\rho_{{\rm crit},0}} = \Omega_{\Lambda 0}\,\frac{\delta\rho_\Lambda(x)}{\bar{\rho}_\Lambda},
$$


with variance $\sigma_X^2 \equiv \langle X^2 \rangle$ ($= \sigma^2/2\theta$ for the OU process; equivalently, rms $\delta\rho_\Lambda/\bar{\rho}_\Lambda = \sigma_X/\Omega_{\Lambda 0}$, below $3.84\times 10^{-2}$ for the bound that follows). The DESI DR2 BAO data, combined with the sensitivity kernel $S(z) = \partial \ln D_V / \partial \Omega_\Lambda$, yield [1]:


$$
\sigma_X < 2.6\times 10^{-2}
$$


**(95% CL).**

This is a direct empirical statement: the dark-energy density is homogeneous at the level of the residual bound (profile 95% CL $\sigma_X<2.6\times 10^{-2}$; null preferred) across the observable universe.

The Planck 2018 CMB does not add a bound at this level: at σ_X = 2.633×10⁻² the residual is indistinguishable from ΛCDM in TT (ℓ < 30) and lensing (`analisis_cmb/RESULTADO.md`).

---

## 3. The Logical Dichotomy

Under the classical theory of fields and Einstein’s equations, only two structures can account for this level of homogeneity.

### 3.1 Static geometry: $\Lambda$ as an integration constant

If dark energy is a pure cosmological constant, its energy-momentum tensor is


$$
T_{\mu\nu} = -\Lambda g_{\mu\nu}.
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
\delta\ddot{\phi} + 3H \delta\dot{\phi} + \left(\frac{k^2}{a^2} + m_\phi^2\right)\delta\phi = 0.
$$


If $m_\phi \ll H_0$ (where $H_0 \sim 10^{-33}\ \text{eV}$), the friction term $3H \delta\dot{\phi}$ forces the mode into the “frozen” regime and effectively erases inhomogeneities. A scalar-field reading would translate $\sigma_X$ into a bound on $m_\phi$; it is not derived here.


In this regime the field is observationally indistinguishable from a cosmological constant. Both alternatives lead to the same phenomenological outcome: $\sigma_X$ is zero or below detection.

---

## 4. The Two Symmetry Candidates

To decide which of the two alternatives is more fundamental, I ask: what symmetry of nature prohibits local quantum fluctuations of matter from coupling to the cosmic volume element? Two candidates have been proposed.

### 4.1 Candidate A: Conformal (Weyl) invariance

Weyl invariance requires the action to be invariant under local rescalings of the metric: $g_{\mu\nu}(x) \to \Omega^2(x) g_{\mu\nu}(x)$. Classically this forces the trace of the energy-momentum tensor to vanish: $T^{\mu}{}_{\mu} = 0$.

**Problem:** At the quantum level, Weyl invariance is broken by the conformal anomaly:


$$
\langle T^{\mu}{}_{\mu} \rangle_{\mathrm{anomaly}} \propto R^2 + C_{\mu\nu\rho\sigma} C^{\mu\nu\rho\sigma}.
$$


Since the local curvature $R(x)$ varies in space and time because of structure formation, the anomaly reintroduces local vacuum fluctuations. Weyl invariance therefore cannot guarantee $\sigma_X = 0$ in a universe with matter.

### 4.2 Candidate B: Volume-preserving diffeomorphisms (SDiff)

In unimodular gravity the gauge group is reduced from the full diffeomorphism group $\mathrm{Diff}(M)$ to the subgroup that preserves the four-volume element:


$$
\sqrt{-g} = \text{constant}.
$$


The resulting field equation is the trace-free part of Einstein’s equation:


$$
R_{\mu\nu} - \frac{1}{4} R g_{\mu\nu} = 8\pi G \left(T_{\mu\nu} - \frac{1}{4} T^{\lambda}{}_{\lambda} g_{\mu\nu}\right).
$$


Any local vacuum fluctuation of matter has the form


$$
T_{\mu\nu}^{\mathrm{vac}}(x) = V(x) g_{\mu\nu}.
$$


Substituting into the trace-free equation:


$$
V(x) g_{\mu\nu} - \frac{1}{4}(4V(x)) g_{\mu\nu} = 0 \quad \text{identically}.
$$


**Conclusion:** In unimodular gravity the trace-free projection cancels any local vacuum fluctuation, regardless of its magnitude or inhomogeneity. This is a purely algebraic cancellation, independent of the details of the quantum fields.

---

## 5. The Cosmological Constant as a Global Integration Constant

To see how $\Lambda$ emerges in unimodular gravity, I take the divergence of the trace-free field equation. Using the Bianchi identity and local conservation of energy-momentum, one obtains:


$$
\nabla_\nu \left(R + 8\pi G T^{\lambda}{}_{\lambda}\right) = 0.
$$


This implies that $R + 8\pi G T^{\lambda}{}_{\lambda}$ is a constant throughout spacetime:


$$
R + 8\pi G T^{\lambda}{}_{\lambda} = 4\Lambda,
$$


where $\Lambda$ is a global integration constant. Substituting this back into the trace-free equation recovers the standard Einstein equations:


$$
G_{\mu\nu} + \Lambda g_{\mu\nu} = 8\pi G T_{\mu\nu}.
$$


The value of $\Lambda$ is given by the global average:


$$
\Lambda = \frac{1}{4V_4} \int_M d^4x \sqrt{-g} \left(R + 8\pi G T^{\lambda}{}_{\lambda}\right),
$$


where $V_4 = \int d^4x \sqrt{-g}$ is the total four-volume of spacetime.

I want to be careful about the precise scope of this result. The unimodular structure ensures that local vacuum energy contributions of the form $T_{\mu\nu} \propto g_{\mu\nu}$ are projected out by the trace-free Einstein equations. That gives a natural mechanism for suppressing inhomogeneous fluctuations in the vacuum energy density.

This protection does not, however, extend to the homogeneous background evolution. A time-dependent dark energy component with equation of state $w(z) \neq -1$ remains compatible with unimodular gravity, as long as fluctuations around that background are sufficiently suppressed. In this sense unimodular gravity constrains the *fluctuations* of the vacuum but does not, by itself, require the background energy density to be strictly constant. I discuss this distinction in more detail in `sdiff-fundamental-vs-emergent.md`.

---

## 6. The Meaning of the Equation: A Constraint, Not a Prediction

The equation above is not a predictive equation. It is a consistency condition — a global constraint that the universe must satisfy. If one tried to use it to calculate $\Lambda$ deductively, one would fall into a circular argument:

1. To compute $V_4$ and the integral of $R$, one needs the global metric $g_{\mu\nu}$.
2. To know $g_{\mu\nu}$, one must solve the Einstein equations.
3. To solve the Einstein equations, one needs $\Lambda$ in advance.

So the equation does not determine the numerical value of $\Lambda$. It only tells us that, whatever its value, it is a global constant — fixed by the initial or boundary conditions of the universe — and that local fluctuations of matter cannot affect it.

The smallness of $\Lambda$($\sim 10^{-122}$ in Planck units) is not a fine-tuning problem in this framework. It is a consequence of the vastness of the universe: the denominator $V_4$ is enormous, so the average of any finite quantity over it tends to zero. The discrepancy between the quantum prediction and the observed value is not an inconsistency; it is the signature that the quantum fluctuations are local, while $\Lambda$ is global.

The same four-volume, built from discrete elements, is the origin of the Everpresent Λ fluctuation δΛ ∼ V_4^{-1/2} (Ahmed et al. 2004; Sorkin 2007); see the end of Section 7.

---

## 7. The Cancellation of Local Vacuum Fluctuations ($\sigma_X = 0$)

The reason I take unimodular gravity as the natural candidate for the extreme homogeneity of the universe is the algebraic structure of its matter source. If I model local quantum vacuum fluctuations as a stochastic contribution of the form


$$
T_{\mu\nu}^{\mathrm{vac}}(x) = V(x) g_{\mu\nu},
$$


then its coupling to the metric is processed through the trace-free tensor:


$$
T_{\mu\nu}^{\mathrm{vac}} - \frac{1}{4} T^{\lambda}{}_{\lambda}{}^{\mathrm{vac}} g_{\mu\nu} = V(x) g_{\mu\nu} - \frac{1}{4}(4V(x)) g_{\mu\nu} \equiv 0.
$$


**Consequence:** Any local fluctuation of the quantum vacuum — regardless of its origin, energy scale, or stochastic behaviour — is identically cancelled at the algebraic level before it can curve spacetime. The metric is completely immune to this local quantum noise.

This decouples the macroscopic value of the observed cosmological constant from the catastrophic 120-order-of-magnitude problem of the quantum zero-point energy, and fixes the variance of the dark-energy density exactly to:


$$
\sigma_X = 0.
$$

**Relation to Everpresent Λ.** Ahmed, Dodelson, Greene and Sorkin (2004) start from the same unimodular framework — Λ conjugate to the four-volume V_4 — and, adding causal-set discreteness, reach the opposite conclusion: Λ fluctuates, δΛ ∼ V_4^{-1/2}, at the level of the critical density in every epoch (Everpresent Λ; Zwane et al. 2018). The two statements concern different objects. The cancellation in this section acts on the local term V(x) g_μν of the trace-free equations; it says nothing about the global Λ conjugate to V_4, which is the object whose fluctuations Everpresent Λ predicts. Which picture describes our universe is decided empirically, not by this algebra. The order-one version, in its two published implementations, does not fit DESI DR2 BAO together with the CMB: once r_d is computed for each realization, dark energy is small around recombination and ω_m is within ±2σ of Planck, no realization fits as well as ΛCDM (0/20000 and 0/500000; `analisis_cmb/everpresent/EVERPRESENT_vs_DESI.md`). A small stationary version is bounded by σ_X < 2.6×10⁻². In this sense σ_X = 0 is the prediction of the local cancellation for local fluctuations, not a theorem about Λ.


---

## 8. Observational Prospects: The Test of the Full Euclid DR1 Release (expected mid-2027)

The hypothesis that unimodular gravity is the true regulatory framework behind cosmic homogeneity is highly testable in the medium term through empirical exclusion:

| Test | Prediction of Unimodular Gravity | What would falsify it |
|:-----|:--------------------------------|:---------------------|
| Euclid DR1 BAO (full release, expected mid-2027) | $\sigma_X = 0$ exactly (within experimental uncertainty) | Any statistically significant deviation ($\sigma_X > 0$) |
| Gravitational waves | Only two tensor polarisations (same as GR) | Detection of a scalar polarisation mode |
| ISW effect | No late-time ISW anomaly from dark-energy clustering | Detection of a residual ISW signal from dark-energy perturbations |

*Note on the ISW row: at the DESI ceiling, Planck 2018 TT at ℓ < 30 and lensing do not distinguish the residual from ΛCDM (p95 Δχ² ≤ 0.067), so with current precision this row cannot falsify anything.*

If Euclid DR1 detects any stochastic fluctuation in the dark-energy density above the noise, unimodular gravity would be ruled out. If it confirms the null result (and tightens the limit to $\sigma_X \lesssim 10^{-5}$), the case for unimodular gravity would become compelling.

---

## 9. Conclusion

I have traced the empirical limit $\sigma_X < 2.6\times 10^{-2}$ to its logical conclusion. The only known structure that can guarantee the absolute smoothness of the vacuum in the presence of local quantum fluctuations is unimodular gravity — a theory whose gauge group is the group of volume-preserving diffeomorphisms $\mathrm{SDiff}(M)$.

This conclusion does not rest on speculative models. It follows from:

- the observational constraint that $\sigma_X$ is below $2.6\times 10^{-2}$ (95% CL);
- the logical dichotomy that only a pure constant or a frozen field can satisfy this bound;
- the mathematical fact that Weyl invariance is broken by the conformal anomaly, while SDiff cancels local vacuum fluctuations identically.

The cosmological constant $\Lambda$ emerges as a global integration constant, fixed by the initial or boundary conditions of the universe. Its numerical value remains a mystery, but its homogeneity is rigorously guaranteed: $\sigma_X = 0$. The equation that defines $\Lambda$ is a consistency condition, not a predictive equation — a constraint that the universe must satisfy, not a formula from which we can derive its value.

I therefore propose that the **Principle of Vacuum Smoothness** — the requirement that $\sigma_X = 0$ exactly — is not an ad-hoc assumption but a necessary condition for any consistent theory of gravity coupled to quantum matter. Unimodular gravity satisfies this principle naturally.

Future surveys, such as the full Euclid DR1 release with cosmology products (expected mid-2027), may further tighten the limit on $\sigma_X$ or even detect a residual signal. If the smoothness persists, I think the theoretical community must take seriously the possibility that the universe is geometrically unimodular at its foundation.

---

## Acknowledgements

I thank the DESI collaboration for making the DR2 BAO data publicly available.

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

