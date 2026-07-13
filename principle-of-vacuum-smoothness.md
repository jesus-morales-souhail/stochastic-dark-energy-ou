# The Principle of Vacuum Smoothness: What Quantum Gravity Must Explain

**Author:** Jesús Morales Souhail  
**Date:** July 2026  
**ORCID:** 0009-0000-7637-1818  
**Repository:** https://github.com/AshPokemonTCG/stochastic-dark-energy-ou  
**Status:** Preprint — not peer reviewed  
**Type:** Theoretical / Methodological Essay

---

## Abstract

Recent BAO measurements from DESI DR2 establish a strict upper limit on stochastic fluctuations of the dark-energy density: $\sigma_X < 1.5 \times 10^{-4}$ (95% CL). When combined with the exclusion of coherent quantum-fluid models ($\Delta\chi^2 \approx +22.7$), these results suggest that the late-time vacuum is smooth to an extraordinary degree. This paper argues that this observational fact should be elevated from a mere empirical limit to a methodological principle: any viable theory of quantum gravity must predict $\sigma_X = 0$ (or a suppression indistinguishable from zero) as a consequence of its fundamental symmetries, not as a fine-tuned adjustment. We explore the implications of this principle for three classes of quantum-gravity models — unimodular gravity, asymptotic symmetries (BMS/holography), and self-tuning mechanisms — and identify the key theoretical question that remains open: what symmetry of nature forbids the vacuum from fluctuating? We conclude that the path toward a consistent quantum gravity must begin by accepting the rigidity of the cosmological vacuum as a foundational constraint.

---

## 1. Introduction: The Observational Fact

The $\Lambda$CDM model has been the standard framework for cosmology for decades. Its central feature — a cosmological constant with no dynamics and no fluctuations — has always been treated as a phenomenological placeholder, waiting to be explained by a more fundamental theory. The expectation has been that quantum gravity would eventually provide a dynamical origin for $\Lambda$ and perhaps even predict small, observable fluctuations.

That expectation has not been fulfilled.

The DESI DR2 BAO data, analysed within a phenomenological Ornstein-Uhlenbeck (OU) framework, show no evidence for stochastic fluctuations [1]. The Maximum Likelihood Estimation drives the amplitude to zero:

$$
\sigma_X \to 0, \qquad \omega_R \to 0,
$$

establishing a conservative upper limit:

$$
\sigma_X < 1.5 \times 10^{-4} \quad (95\%\ \text{CL}).
$$

Moreover, a specific class of models — those based on a tachyonic quantum fluid with negative effective mass — is excluded with high significance:

$$
\Delta\chi^2 \approx +22.7,
$$

meaning that the data reject a coherent growing-mode scenario entirely. The only viable limit for such models is $t_c \to \infty$, i.e., no growth, which reduces them to a smooth background indistinguishable from $\Lambda$CDM [2].

These results are not a "non-discovery." They are a structural constraint on the space of possible theories.

---

## 2. The Principle: From Limit to Axiom

In standard scientific practice, an empirical upper limit is simply that: a limit. It does not prove that the quantity is exactly zero; it only bounds it. However, in the case of $\sigma_X$, the limit is so small ($10^{-4}$) and the models that predict larger values are so severely penalised, that a methodological shift is warranted.

We propose the following principle:

> **Principle of Vacuum Smoothness (PVS):**  
> Any viable theory of quantum gravity must predict $\sigma_X = 0$ (or a suppression so strong that it is observationally indistinguishable from zero) as a consequence of its fundamental structure, not as a result of parameter tuning.

This principle is not a claim that $\sigma_X$ is exactly zero. It is a criterion for theory selection: theories that cannot guarantee this level of smoothness are, for all practical purposes, incompatible with the observed universe.

The PVS inverts the standard "top-down" methodology of quantum gravity. Instead of constructing a theory at the Planck scale and then checking whether its infrared limit produces fluctuations, we start from the observed infrared smoothness and demand that any fundamental theory reproduce it automatically.

---

## 3. The Filter: Excluded Classes of Models

The PVS acts as a filter on the space of quantum-gravity models. The following classes are excluded:

### 3.1. Models with intrinsic local fluctuations

Any theory in which the vacuum energy density is a local quantum field (e.g., quintessence, $k$-essence, or any scalar field with a non-trivial potential) will generate quantum fluctuations. Even if the background is smooth, the field's quantum nature implies $\langle \delta\phi^2 \rangle \neq 0$, which propagates to $\sigma_X > 0$. Unless the field is extremely massive (which would make it irrelevant for dark energy) or has a strong shift symmetry that suppresses fluctuations, such models are incompatible with the PVS.

### 3.2. Models with coherent growing modes

The tachyonic quantum fluid analysed in [2] is a prime example. Its instability amplifies primordial seeds as

$$
\sigma_X(t) = \sigma_0 \, e^{t/t_c}.
$$

Even with the smallest possible seed ($\sigma_0 \sim 10^{-61}$), the growth would exceed the observed limit unless $t_c > 10^8$ years — which, as shown in [2], is not sufficient to save the model. The model is excluded outright.

### 3.3. Models with stochastic memory (OU-type processes)

The Ornstein-Uhlenbeck process, while not producing coherent growth, generates redshift correlations that would be detectable in BAO residuals. The fact that the MLE drives $\sigma_X \to 0$ and that the predicted lag correlations are absent (the observed pattern being an algebraic artefact of whitening with $N = 7$ bins) means that any model with a similar memory kernel is strongly disfavoured [1].

---

## 4. Surviving Candidates: What Must a Theory of QG Look Like?

If the PVS is adopted, the number of viable quantum-gravity candidates shrinks dramatically. The surviving theories must, by construction, have no fluctuating vacuum energy. We consider three possibilities.

### 4.1. Unimodular Gravity

In unimodular gravity, the metric determinant is constrained: $\sqrt{-g} = \text{constant}$. This reduces the gauge group from $\mathrm{Diff}(M)$ to $\mathrm{SDiff}(M)$ (volume-preserving diffeomorphisms). As a consequence, the cosmological constant $\Lambda$ appears not as a source term but as an integration constant of the field equations.

**Advantage:** By construction, $\Lambda$ is a global constant. It has no local dynamics, no propagating degrees of freedom, and no fluctuations. The PVS is satisfied exactly: $\sigma_X = 0$.

**Challenge:** The theory does not predict the value of $\Lambda$. It only says that it does not fluctuate. The numerical value $10^{-122}$ remains unexplained. Moreover, the quantum formulation of unimodular gravity is not fully understood; the path integral over $\mathrm{SDiff}(M)$ is technically more complex than standard GR.

**Verdict:** Unimodular gravity is the most natural fit for the PVS, but it is a phenomenological refuge, not a complete theory. It solves the fluctuation problem but leaves the cosmological constant problem unsolved.

### 4.2. Asymptotic Symmetries (BMS) and Holography

The Bondi-Metzner-Sachs (BMS) group describes the asymptotic symmetries of asymptotically flat spacetimes. In recent years, it has been suggested that BMS supertranslations may play a role in the infrared structure of gravity and in the holographic encoding of bulk degrees of freedom.

**Advantage:** If the vacuum energy is encoded in the boundary data (e.g., as a charge of the asymptotic symmetry group), then fluctuations of the bulk vacuum would correspond to violations of the boundary conditions, which are prohibited. The PVS could be a consequence of holographic consistency.

**Challenge:** This is speculative. There is no rigorous derivation that connects BMS symmetries to a fluctuating bulk vacuum. The connection to de Sitter space (where the universe actually lives) is even more tenuous, as BMS symmetries are primarily defined for asymptotically flat spacetimes.

**Verdict:** Promising but not yet a theory. The PVS provides a motivation to pursue this line of research, but it is not a solution.

### 4.3. Self-Tuning Mechanisms with Trans-Planckian Suppression

In modified gravity models (e.g., DGP, galileons, or chameleon-like screening), the vacuum energy can be dynamically cancelled by a non-linear coupling between the scalar field and the curvature.

**Advantage:** These models can, in principle, cancel large vacuum energies and leave a small, smooth remnant.

**Challenge:** The limit derived in Appendix H of [1] imposes a severe constraint. If the dark-energy field is coupled to the electromagnetic sector, the coupling constant must satisfy:

$$
|\zeta| \lesssim 1.8 \times 10^{-2},
$$

which implies a cutoff scale:

$$
M_{\mathrm{cut}} \equiv \frac{M_{\mathrm{Pl}}}{|\zeta|} \gtrsim 55 \, M_{\mathrm{Pl}}.
$$

This places the new physics at trans-Planckian scales, where the effective field theory breaks down. A UV-complete theory would be required, and none exists without introducing ghosts or instabilities.

**Verdict:** Self-tuning models are effectively excluded by the PVS unless a consistent UV-completion is found. The cutoff scale being above the Planck scale suggests that the mechanism is either unnatural or requires a symmetry we do not yet understand.

---

## 5. The Open Question: What Symmetry Forbids Fluctuations?

The PVS shifts the central question of quantum gravity. It is no longer:

> "Why is the cosmological constant so small?"

but rather:

> "What symmetry of nature forbids the vacuum from fluctuating?"

Three candidate symmetries are worth exploring:

### A. Scale invariance (conformal symmetry)

If the vacuum state is exactly scale-invariant, then $\langle T_{\mu\nu} \rangle = 0$ identically, and there are no fluctuations. This would require that the cosmological constant is generated by a spontaneous breaking of scale invariance, but once broken, the fluctuations would be suppressed by the Goldstone mechanism. This is the idea behind some models of "conformal cosmology."

### B. Shift symmetry (axionic protection)

If the dark-energy field is an axion-like particle with a shift symmetry $\phi \to \phi + c$ (which is broken only by a periodic potential), then quantum corrections to the potential are suppressed. The mass of the field is protected, and fluctuations are negligible. This is the standard argument for quintessence with $m_\phi \lesssim 10^{-5}\ \text{eV}$.

**However:** The PVS requires more than a small mass. It requires that $\sigma_X$ be suppressed to better than $10^{-4}$. Even for an axion, quantum fluctuations of the field itself would generically produce $\langle \delta\phi^2 \rangle > 0$, which would propagate to $\sigma_X > 0$. To satisfy the PVS, the field would need to be frozen to its vacuum expectation value with no local excitations — which is tantamount to saying that it is not a field at all.

### C. Asymptotic safety and the fixed-point structure

In asymptotically safe gravity, the cosmological constant is a relevant coupling that flows to a fixed point in the ultraviolet. If the fixed point is such that $\Lambda = 0$ exactly (or very small) in the infrared, then the vacuum energy could be naturally small. But would it be smooth? The fluctuations of the metric itself would still be present in the quantum theory. The PVS would require that these metric fluctuations do not couple to the vacuum energy in a way that produces observable $\sigma_X$. This is a non-trivial requirement that has not been demonstrated.

---

## 6. Conclusion: The Path Forward

The PVS is not a theory. It is a constraint on theories. It tells us what any successful theory of quantum gravity must not do: it must not produce measurable fluctuations of the cosmological vacuum.

This constraint is remarkably strong. It excludes:

- All models with local quantum fields as the source of dark energy.
- All models with coherent instabilities that amplify fluctuations.
- All stochastic models with memory kernels that produce redshift correlations.
- Most self-tuning mechanisms (due to the trans-Planckian cutoff).

The surviving candidates are few:

- Unimodular gravity (which provides a phenomenological explanation but no prediction of $\Lambda$).
- Holographic or asymptotic-symmetry models (which are promising but not yet developed).
- UV-complete theories with symmetries that protect the vacuum (which we have not yet identified).

The PVS reframes the problem of quantum gravity. It tells us that the correct theory must have a symmetry that makes the vacuum rigid — unresponsive to local perturbations, incapable of sustaining excitations, and devoid of memory. The search for that symmetry is now the central task of theoretical cosmology.

---

## Acknowledgements

The author thanks the DESI collaboration for making the DR2 BAO data publicly available. This work has been supported by the open-science philosophy of reproducible research. The ideas presented here are an interpretation of the empirical results obtained in [1, 2] and do not constitute a new observational claim.

---

## References

[1] Morales Souhail, J., "Constraints on Stochastic Dark Energy from DESI DR2: A Null Result for Ornstein-Uhlenbeck Fluctuations," arXiv:xxxx.xxxxx (2026).

[2] Morales Souhail, J., "Constraints on Stochastic Dark Energy from Quantum Fluid Instabilities and DESI DR2 Baryon Acoustic Oscillations," arXiv:yyyy.yyyyy (2026).

[3] DESI Collaboration, "DESI DR2 Results II: Measurements of BAO and Cosmological Constraints," arXiv:2503.14738 (2025).

[4] Bekenstein, J. D., "Black Holes and Entropy," *Phys. Rev. D* 7, 2333 (1973).

[5] Sorkin, R. D., "Is the Cosmological 'Constant' a Nonlocal Quantum Residual?," arXiv:gr-qc/0503057 (2005).

[6] Weinberg, S., "The Cosmological Constant Problem," *Rev. Mod. Phys.* 61, 1 (1989).

[7] Carroll, S. M., "Quintessence and the Rest of the World," *Phys. Rev. Lett.* 81, 3067 (1998).

[8] Ashtekar, A. et al., "BMS supertranslations and Weinberg's soft graviton theorem," *JHEP* 2015, 152 (2015).

[9] Reuter, M., "Nonperturbative Evolution Equation for Quantum Gravity," *Phys. Rev. D* 57, 971 (1998).

[10] Dvali, G., Gabadadze, G., & Porrati, M., "4D Gravity on a Brane in 5D Minkowski Space," *Phys. Lett. B* 485, 208 (2000).

---

> **Note on the title:** The term "smoothness" is used here in the sense of absence of stochastic fluctuations, not in the sense of differentiability of the metric. The vacuum is smooth because it does not fluctuate, not because it is mathematically differentiable.