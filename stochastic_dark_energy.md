# Stochastic Dark Energy from a Finite Information Network: BAO Covariance, Cross-Correlation Evidence, and the Case for a Refined Kernel

**Author:** Jesús Morales Souhail  
**Date:** July 2026  
**Version:** v3.0 (Revised)

---

## Abstract

We present a phenomenological framework in which the cosmological constant is not strictly constant but includes a small stochastic component motivated by the finite information content of the observable universe (Bekenstein–Hawking bound) and the root‑√N fluctuation scaling suggested by unimodular gravity and causal‑set arguments. We model late‑time fluctuations in Ω_Λ as an Ornstein–Uhlenbeck (OU) process in the logarithmic scale factor x = ln a, with mean‑reversion rate θ and diffusion amplitude σ.

We derive three observational tests using public DESI and Pantheon+ data:

1. **BAO likelihood with OU covariance.** The OU model improves the fit over flat ΛCDM with Δlog𝓛 = +6.75 for DESI DR2 BAO data (7 bins), yielding θ = 0.765 and σ_X = 0.018. AIC and BIC favour the OU extension. **With 7 bins, this result is preliminary.**

2. **Angular cross‑correlation (δ_g × δΩ_Λ).** Using DESI DR1 LRG galaxy density maps and Pantheon+ Hubble residuals, we find r_cross = 0.1673 ± 0.0613 (2.73σ excess). **This is preliminary: only 67 overlapping pixels and systematics uncontrolled.**

3. **Redshift lag correlations in BAO residuals.** The OU kernel predicts ρ_1 ≈ 0.78, ρ_2 ≈ 0.62, ρ_3 ≈ 0.49. DESI DR2 shows mixed results: lag-1 positive (consistent), lags 2–3 negative (inconsistent). **With N=7 bins and 95% CI ≈ ±1.0, no lag is individually significant.**

**Conclusion:** The finite‑information stochastic framework is **not falsified** by current data. The OU likelihood improves the fit, but the model remains degenerate with only 7 bins. The decisive test requires >20 bins from Euclid DR1 (expected H2 2026).

---

## 1. Introduction

The ΛCDM model has been remarkably successful in describing a wide range of cosmological observations. However, recent BAO measurements from DESI, when combined with CMB and supernova data, have shown a statistical preference for dynamical dark energy in the w_0 > -1, w_a < 0 quadrant, with tensions ranging from 2.6σ to 4.2σ depending on dataset combinations [1, 2]. While not yet a discovery, these tensions motivate exploring minimal, falsifiable extensions that do not require new classical fields.

This work proposes one such extension. The starting point is the Bekenstein–Hawking entropy bound, which implies that the observable universe has a finite number of effective degrees of freedom (N ~ 10^122 in Planck units). In unimodular gravity, the cosmological constant Λ appears as an integration constant conjugate to spacetime four‑volume. If spacetime is fundamentally discrete (as in causal‑set theory), Poisson fluctuations in the number of elements N induce residual fluctuations δΛ ~ 1/√N. This is not a derivation of the observed value of Λ, but it provides a plausible order‑of‑magnitude motivation for a small stochastic component.

Rather than attempting a full quantum‑gravity derivation, we work at the phenomenological level and model the late‑time fluctuations of Ω_Λ as an Ornstein–Uhlenbeck process in ln a. The key observable consequence is an additive covariance term in BAO distance measurements: a redshift‑dependent "precision floor" that cannot be beaten by increasing sample size alone. We test this prediction against public DESI DR1/DR2 BAO products and Pantheon+ supernova data, and we outline the analysis protocol for the upcoming Euclid DR1.

---

## 2. Axiomatic Foundation

### Axiom A1: Finite Information Bound
The maximum entropy S contained in a region with horizon area A satisfies the Bekenstein–Hawking bound:

$$S \leq \frac{A}{4 G \hbar}.$$

For the observable universe, A ~ 10^122 (in Planck units), which implies a finite effective Hilbert‑space dimension. Continuum field descriptions are therefore effective coarse‑grainings, not fundamental.

### Axiom A2: Stochastic Λ from Discreteness (Sorkin Mechanism)
In unimodular gravity, Λ can be interpreted as a constant of integration conjugate to the spacetime four‑volume. If spacetime consists of N discrete elements with N = V / L_P^4, Poisson fluctuations yield

$$\delta \Lambda \sim \frac{1}{\sqrt{N}}.$$

With N ~ 10^122, this gives δΛ ~ 10^-61 in Planck units. This motivates the existence of a small, non‑zero stochastic component, though we do not claim a full derivation of ρ_Λ from this argument.

### Axiom A3: Effective Stochastic Dynamics (OU Closure)
Define X(x) ≡ δΩ_Λ(x), where x = ln a is the logarithmic scale factor. We model X as an Ornstein–Uhlenbeck process:

$$dX = -\theta \, X \, dx + \sigma \, dW_x,$$

with stationary variance

$$\operatorname{Var}(X) = \frac{\sigma^2}{2\theta}.$$

The OU process captures finite memory and yields analytic redshift correlations. It is adopted as a minimal phenomenological closure.

### Axiom A4: Late‑Time Activation (Degenerate Parameter)

To preserve early‑universe constraints (CMB, BBN), a smooth activation factor g(z) is conceptually introduced:

$$\sigma_{\rm eff}(z) = \sigma \, g(z), \qquad g(z) = \frac{1}{1 + \exp[-(x - x_*) / \Delta]},$$

where x = ln(1/(1+z)). However, for the redshift range probed by DESI (z ≲ 2.3) and Euclid (z ≲ 2.0), g(z) is essentially unconstrained by data. We set z_* = 1.5 (illustrative), with the understanding that this parameter is degenerate with σ and cannot be independently calibrated. For all numerical results in Sections 4–6, g(z) ≈ 1 in the observed range, meaning the effective amplitude σ_X absorbs any early‑time suppression.

---

## 3. BAO Sensitivity Kernel and Precision Floor

### 3.1 Definition of S(z)
We use the standard isotropic BAO distance proxy:

$$D_V(z) = \left[ D_M(z)^2 \, \frac{c z}{H(z)} \right]^{1/3},$$

Define the sensitivity kernel to Ω_Λ:

$$S(z) \equiv \frac{\partial \ln D_V(z)}{\partial \Omega_\Lambda}.$$

Linear propagation yields an induced BAO scatter (the "precision floor"):

$$\sigma_{\alpha,\rm floor}(z) \approx |S(z)| \, \sigma_{\Omega_\Lambda}(z).$$

### 3.2 Calibration with DESI DR1 (Illustrative Upper Limit)

To illustrate the magnitude of the precision floor, we allocate a conservative fraction f_net = 0.15 of the observed variance at the anchor bin z_eff = 0.71 to an irreducible network component. **This is a plausible upper limit for current DESI precision; the actual fraction is unconstrained and is precisely what the data will test.** Using σ_α,obs(0.71) ≈ 0.012:

$$\sigma_{\alpha,\rm floor}(0.71) = \sqrt{0.15} \times 0.012 = 4.65 \times 10^{-3}.$$

Using the kernel value S(0.71) = -0.312, we infer σ_Ω_Λ(0.71) ≈ 1.49 × 10^-2. For θ = 1.2, this gives σ ≈ 2.31 × 10^-2. **These parameters are illustrative only; the actual fit (Section 4.2) treats θ and σ_X as free parameters.**

### 3.3 Precision Floor Implementation (Discrete Kernel)

For the numerical results in Sections 4–6, we use the discrete sensitivity kernel S(z). This yields:

| z_eff | Tracer       | \|S(z)\|/\|S(0.706)\| | σ_α,floor |
|:---:|:---|:---:|:---:|
| 0.295 | BGS          | 0.477                 | 2.22 × 10^-3  |
| 0.510 | LRG1         | 0.777                 | 3.61 × 10^-3  |
| **0.706** | **LRG2** | **1.000**             | **4.65 × 10^-3** |
| **0.934** | **LRG3+ELG1** | **1.208**       | **5.62 × 10^-3** |
| 1.321 | ELG2         | 1.462                 | 6.80 × 10^-3  |
| 1.484 | QSO          | 1.541                 | 7.17 × 10^-3  |
| 2.330 | Lyα          | 1.798                 | 8.36 × 10^-3  |

> **Note:** An exact integral susceptibility kernel χ(z, z') is presented in Appendix C for theoretical completeness. It is **not** used in this version but provides a roadmap for v3.1.

---

## 4. Test 1: BAO Likelihood with OU Covariance

### 4.1 Methodology
The standard BAO analysis assumes Gaussian residuals with diagonal covariance. Our model adds an OU‑induced component to the total covariance:

$$C_{\rm total} = C_{\rm std} + C_{\rm OU},$$

where $$(C_{\rm OU})_{ij} = S(z_i) S(z_j) \operatorname{Cov}[X(x_i), X(x_j)]$$ and $$\operatorname{Cov}[X(x_i), X(x_j)] = \frac{\sigma^2}{2\theta} \exp[-\theta |x_i - x_j|].$$

### 4.2 Results

| Model               | θ | σ_X | Δlog𝓛 (vs ΛCDM) | AIC      | BIC      |
|:--------------------|:---:|:---:|:---:|:---:|:---:|
| ΛCDM                | —  | —   | 0.00 | ref  | ref  |
| OU pure (H0)        | 0.765 | 0.018 | **+6.75** | −38.73 | −38.84 |
| Oscillatory (QNM)   | 0.000 | 0.019 | +7.00 | −37.24 | −37.40 |

**Interpretation:** With only 7 bins, this result is **preliminary**. The OU covariance structure provides a better fit than white noise, but this does not confirm the physical reality of stochastic dark energy; systematics or different deterministic models could mimic this pattern. The lag‑correlation test (Section 6) is the critical discriminator.

---

## 5. Test 2: Angular Cross‑Correlation δ_g × δΩ_Λ

### 5.1 Data and Method
We construct a HEALPix map (N_side = 32) of galaxy overdensities from DESI DR1 LRG catalogs (NGC: 1.476M, SGC: 662k objects, 0.4 < z < 0.8). As a proxy for δΩ_Λ, we use Pantheon+ Hubble residuals Δμ = μ_obs - μ_ΛCDM.

### 5.2 Results (Preliminary)

Combining NGC and SGC:

$$r_{\rm cross} = 0.1673 \pm 0.0613 \quad (Z \approx 2.73\sigma \text{ excess}).$$

> **Important caveat:** With only **67 overlapping pixels**, the statistical power is limited. **Systematic effects (Galactic dust, imaging systematics) have not been controlled.** This is a preliminary motivation for DESI DR2 analysis with imaging weights, **not a confirmed detection.**

---

## 6. Test 3: Redshift Lag Correlations in BAO Residuals (The Critical Test)

### 6.1 The Smoking‑Gun Prediction

The OU process predicts positive lag correlations. For a DESI‑like 9-bin grid, the predictions are[^1]:

$$\rho_1 \approx 0.78, \quad \rho_2 \approx 0.62, \quad \rho_3 \approx 0.49.$$

[^1]: Generated by `predicted_lag_correlations()` in `ou_bao_likelihood.py` with θ=0.765, σ=0.018, using z = {0.20, 0.35, 0.51, 0.65, 0.80, 0.95, 1.10, 1.25, 1.40}.

### 6.2 Results from DESI DR2

| Lag | DR2 (7 bins) | OU Prediction (θ=1.2) | 95% CI |
|:---:|:---:|:---:|:---:|
| 1   | **+0.37** | +0.83 | ±1.0 |
| 2   | -0.32 | +0.85 | ±1.0 |
| 3   | -0.81 | +0.85 | ±1.0 |

### 6.3 Interpretation

For N = 7 bins, the standard error is σ_ρ ≈ 1/√(N-3) ≈ 0.5, giving 95% CI ≈ ±1.0. **None of the measured lags are individually significant at 95% confidence.** The data are consistent with white noise at current resolution.

This apparent inconsistency with the positive BAO likelihood (Δlog𝓛 = +6.75) is resolved by noting that the global covariance structure can improve fit even if individual lags are insignificant. **This ambiguity requires >20 bins to resolve.**

---

## 7. Unified Discussion

| Test | Result | Status |
|:---|:---|:---|
| **BAO likelihood** | OU improves fit (Δlog𝓛 = +6.75) | Preliminary, 7 bins |
| **Cross-correlation** | r = 0.167 (2.73σ), 67 pixels | Preliminary, systematics uncontrolled |
| **Lag correlations** | Mixed; none individually significant | Not conclusive with 7 bins |

**Current status:** The stochastic framework is **not falsified** but **not confirmed**. It remains viable and falsifiable. Higher-z data will decide.

---

## 8. Falsification Criteria (Data‑Kills‑Model)

**Current status of falsification.** With only 7 BAO bins, the model space is highly degenerate. DESI DR2 alone cannot definitively falsify or confirm the stochastic framework. The criteria below apply to future data (Euclid DR1, DESI DR3). **Presently, none are met—meaning the model is not ruled out, but it is not confirmed.**

| Criterion | Condition for Exclusion | Current Status |
|:---|:---|:---|
| **F1** (Variance floor absent) | σ²_α,obs < σ²_α,floor in multiple bins | **Alive** |
| **F2** (w_a → 0 at >5σ) | No room for stochastic component | **Alive** |
| **F3** (ISW exclusion) | σ_Ω_Λ incompatible with CMB-LSS | **Pending** |
| **F4** (Bayesian evidence) | ln K > 5 favouring ΛCDM | **Not computed** |
| **F5** (Lag correlations) | All lags ≤ 0 with 20+ bins | **Partial (7 bins)** |

---

## 9. Near‑Term Observational Program: Euclid DR1

Euclid Data Release 1 (expected H2 2026) will provide >20 redshift bins, resolving the current degeneracies with unprecedented statistical power. The analysis pipeline is fully ready.

**Notably, Euclid DR1's narrower redshift baseline (z ∈ [0.9, 1.8]) yields a higher Rayleigh frequency limit (ω_R,min ≈ 16.2) than DESI (6.66), meaning Euclid cannot geometrically resolve intermediate-frequency quasi-normal mode (QNM) oscillations that DESI could potentially detect (Appendix D). The decisive contribution of Euclid will be statistical power through >20 bins, not oscillatory frequency resolution.** Upon release, we will re-run the lag-correlation test and precision-floor measurement with this enhanced sample.

---

## 10. Conclusion

The finite‑information stochastic framework is **not falsified** by DESI DR2 data. The OU likelihood improves the fit and no criterion decisively excludes the model. However, **with only 7 BAO bins, the framework is not yet confirmed**. The pure OU kernel may require refinement (redshift‑dependent θ, oscillatory components). The decisive test requires Euclid DR1 (>20 bins, expected H2 2026).

---

## Appendix A — BAO Sensitivity Kernel S(z): Numerical Implementation

```python
import numpy as np

C_KMS = 299792.458  # km/s

def E_z(z, Om, Ol):
    return np.sqrt(Om * (1 + z)**3 + Ol)

def H_z(z, H0, Om, Ol):
    return H0 * E_z(z, Om, Ol)

def D_M(z, H0, Om, Ol, n=4096):
    zs = np.linspace(0.0, z, n)
    integrand = 1.0 / H_z(zs, H0, Om, Ol)
    chi = np.trapz(integrand, zs)
    return C_KMS * chi

def D_V(z, H0, Om, Ol):
    DM = D_M(z, H0, Om, Ol)
    Hz = H_z(z, H0, Om, Ol)
    return (DM * DM * (C_KMS * z) / Hz) ** (1.0 / 3.0)

def kernel_S_iso(z, H0=67.4, Ol0=0.685, delta=1e-4):
    Ol_p = Ol0 + delta
    Ol_m = Ol0 - delta
    Om_p = 1.0 - Ol_p
    Om_m = 1.0 - Ol_m
    DV_p = D_V(z, H0, Om_p, Ol_p)
    DV_m = D_V(z, H0, Om_m, Ol_m)
    return (np.log(DV_p) - np.log(DV_m)) / (2.0 * delta)

if __name__ == "__main__":
    for z in [0.51, 0.71, 0.93, 1.50]:
        print(f"z={z:.2f}, S(z)={kernel_S_iso(z):+.4f}")
```

---

## Appendix B — Reproducibility and Data Access

All scripts are publicly available:
- `ou_bao_likelihood.py`: BAO likelihood and lag correlations
- `cross_correlation_DESI.py`: Cross-correlation analysis

Data: DESI DR1/DR2 BAO at data.desi.lbl.gov, Pantheon+ at github.com/PantheonPlusSH0ES

---

## Appendix C — Theoretical Outlook: Integral Susceptibility χ(z, z') (for v3.1)

The following derivation is **not** used in v3.0 but provided for future reference.

Starting from the perturbed Hubble parameter:

$$H^2(z) = H_0^2 \left[ \Omega_m (1+z)^3 + \Omega_\Lambda + X(z) \right],$$

one obtains (to first order in X):

$$\delta \ln D_V(z) = -\frac{c}{3 H_0 D_{M,0}(z)} \int_0^z \frac{X(z')}{E_0(z')^3} dz' - \frac{1}{6 E_0(z)^2} X(z).$$

Defining χ(z, z') via δ ln D_V(z) = ∫_0^z χ(z, z') X(z') dz':

$$\chi(z, z') = -\frac{c}{3 H_0 D_{M,0}(z) E_0(z')^3} \Theta(z - z') - \frac{1}{6 E_0(z)^2} \delta_D(z - z').$$

Full numerical implementation deferred to v3.1. Current version uses discrete S(z).

---

## Appendix D — Geometric Rayleigh Resolution Limit (Survey-Independent)

This is a **purely geometric consequence** of survey redshift coverage, independent of any dark‑energy model.

For redshift range Δx = ln(1+z_max) - ln(1+z_min), the minimum resolvable frequency is:

$$\omega_{R,\min} = \frac{2\pi}{\Delta x}.$$

Any oscillation with ω_R < ω_R,min produces <1 visible cycle and is indistinguishable from a monotonic (pure OU) trend.

**DESI DR2:** z ∈ [0.295, 2.330] → Δx ≈ 0.944 → ω_R,min ≈ 6.66

**Euclid DR1:** z ∈ [0.9, 1.8] → Δx ≈ 0.388 → ω_R,min ≈ 16.2

**Implication:** If the true kernel is QNM with ω_R ~ 8, DESI can detect it but Euclid cannot (geometric limit, not noise).

---

## References

[1] DESI Collaboration, arXiv:2404.03000 (2024).

[2] DESI Collaboration, arXiv:2503.14738 (2025).

[3] Bekenstein, J. D., Phys. Rev. D 7, 2333 (1973).

[4] Sorkin, R. D., arXiv:gr-qc/0503057 (2005).

[5] Scolnic, D. et al., ApJ 938, 113 (2022).

[6] Uhlenbeck, G. E. & Ornstein, L. S., Phys. Rev. 36, 823 (1930).
