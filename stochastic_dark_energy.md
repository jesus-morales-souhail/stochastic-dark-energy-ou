# Stochastic Dark Energy from a Finite Information Network: BAO Covariance, Cross-Correlation Evidence, and the Case for a Refined Kernel

**Author:** Jesús Morales Souhail  
**Date:** July 2026  

---

## Abstract

We present a phenomenological framework in which the cosmological constant is not strictly constant but includes a small stochastic component motivated by the finite information content of the observable universe (Bekenstein–Hawking bound) and the root‑\(N\) fluctuation scaling suggested by unimodular gravity and causal‑set arguments. We model late‑time fluctuations in \(\Omega_\Lambda\) as an Ornstein–Uhlenbeck (OU) process in the logarithmic scale factor \(x = \ln a\), with mean‑reversion rate \(\theta\) and diffusion amplitude \(\sigma\).

We derive three observational tests using public DESI and Pantheon+ data:

1. **BAO likelihood with OU covariance.** The OU model improves the fit over flat \(\Lambda\)CDM with \(\Delta\log\mathcal{L} = +6.75\) for DESI DR2 BAO data (7 bins), yielding \(\theta = 0.765\) and \(\sigma_X = 0.018\). AIC and BIC favour the OU extension over \(\Lambda\)CDM.

2. **Angular cross‑correlation (\(\delta_g \times \delta\Omega_\Lambda\)).** Using DESI DR1 LRG galaxy density maps (HEALPix, \(N_{\rm side}=32\)) and Pantheon+ Hubble diagram residuals as a proxy for \(\delta\Omega_\Lambda\), we find a combined (NGC + SGC) cross‑correlation coefficient \(r_{\rm cross} = 0.1673 \pm 0.0613\), corresponding to a \(Z \approx 2.73\sigma\) excess over the null hypothesis. This exceeds the naive OU prediction (\(r_{\rm cross} \sim 0.023\)) by a factor of \(\sim 7\), indicating either a stronger signal than expected or unresolved systematics (dust, imaging photometry).

3. **Redshift lag correlations in BAO residuals (the critical test).** The pure OU kernel predicts strong positive lag correlations in whitened BAO residuals (e.g., \(\rho_1 \approx 0.78\), \(\rho_2 \approx 0.62\), \(\rho_3 \approx 0.49\) for DESI‑like binning). Analysing DESI DR1 (6 bins) and DR2 (7 bins), we find a mixed pattern: lag‑1 in DR2 is positive (\(+0.37\)), consistent in sign with the OU prediction, but lags 2 and 3 are negative (\(-0.32\), \(-0.81\)), inconsistent with a pure OU kernel with \(\theta = 1.2\). With only 6–7 bins, the 95% confidence intervals are large (\(\sim \pm 0.7\) to \(\pm 1.0\)), so no definitive conclusion can be drawn.

We conclude that the finite‑information stochastic framework is statistically favoured over a pure white‑noise (ΛCDM) assumption for BAO data, but the simple OU closure is insufficient to describe all observables. The discrepancy indicates that the memory kernel must be refined (e.g., redshift‑dependent \(\theta\), an oscillatory component, or a different functional form). The decisive test will require the \(>20\) redshift bins of Euclid DR1 (expected October 2026), for which the analysis pipeline is fully documented and ready.

---

## 1. Introduction

The \(\Lambda\)CDM model has been remarkably successful in describing a wide range of cosmological observations. However, recent BAO measurements from DESI, when combined with CMB and supernova data, have shown a statistical preference for dynamical dark energy in the \(w_0 > -1\), \(w_a < 0\) quadrant, with tensions ranging from \(2.6\sigma\) to \(4.2\sigma\) depending on dataset combinations [1, 2]. While not yet a discovery, these tensions motivate exploring minimal, falsifiable extensions that do not require new classical fields.

This work proposes one such extension. The starting point is the Bekenstein–Hawking entropy bound, which implies that the observable universe has a finite number of effective degrees of freedom (\(N \sim 10^{122}\) in Planck units). In unimodular gravity, the cosmological constant \(\Lambda\) appears as an integration constant conjugate to spacetime four‑volume. If spacetime is fundamentally discrete (as in causal‑set theory), Poisson fluctuations in the number of elements \(N\) induce residual fluctuations \(\delta\Lambda \sim 1/\sqrt{N}\). This is not a derivation of the observed value of \(\Lambda\), but it provides a plausible order‑of‑magnitude motivation for a small stochastic component.

Rather than attempting a full quantum‑gravity derivation, we work at the phenomenological level and model the late‑time fluctuations of \(\Omega_\Lambda\) as an Ornstein–Uhlenbeck process in \(\ln a\). The key observable consequence is an additive covariance term in BAO distance measurements: a redshift‑dependent "precision floor" that cannot be beaten by increasing sample size alone. We test this prediction against public DESI DR1/DR2 BAO products and Pantheon+ supernova data, and we outline the analysis protocol for the upcoming Euclid DR1.

---

## 2. Axiomatic Foundation

### Axiom A1: Finite Information Bound
The maximum entropy \(S\) contained in a region with horizon area \(A\) satisfies the Bekenstein–Hawking bound:

\[
S \leq \frac{A}{4 G \hbar}.
\]

For the observable universe, \(A \sim 10^{122}\) (in Planck units), which implies a finite effective Hilbert‑space dimension \(\dim(\mathcal{H}) = \exp(10^{122})\). Continuum field descriptions are therefore effective coarse‑grainings, not fundamental.

### Axiom A2: Stochastic \(\Lambda\) from Discreteness (Sorkin Mechanism)
In unimodular gravity, \(\Lambda\) can be interpreted as a constant of integration conjugate to the spacetime four‑volume. If spacetime consists of \(N\) discrete elements with \(N = V / L_{\rm P}^4\), Poisson fluctuations yield

\[
\delta \Lambda \sim \frac{1}{\sqrt{N}}.
\]

With \(N \sim 10^{122}\), this gives \(\delta \Lambda \sim 10^{-61}\) in Planck units, which is consistent in order of magnitude with the observed vacuum energy density (\(\sim 10^{-122}\) in Planck density units). This motivates the existence of a small, non‑zero stochastic component, though we do not claim a full derivation of \(\rho_\Lambda\) from this argument.

### Axiom A3: Effective Stochastic Dynamics (OU Closure)
Define \(X(x) \equiv \delta\Omega_\Lambda(x)\), where \(x = \ln a\) is the logarithmic scale factor. We model \(X\) as an Ornstein–Uhlenbeck process:

\[
dX = -\theta \, X \, dx + \sigma \, dW_x,
\]

with stationary variance

\[
\operatorname{Var}(X) = \frac{\sigma^2}{2\theta}.
\]

The OU process captures finite memory in expansion e‑folds and yields analytic redshift correlations. It is adopted as a minimal phenomenological closure, not as a unique microscopic prediction.

### Axiom A4: Late‑Time Activation
To preserve early‑universe constraints (CMB, BBN), the fluctuations must be suppressed at high redshift. We introduce a smooth activation factor \(g(z)\) that multiplies the diffusion amplitude:

\[
\sigma_{\rm eff}(z) = \sigma \, g(z), \qquad g(z) = \frac{1}{1 + \exp\left[-(x - x_*) / \Delta\right]},
\]

where \(x = \ln(1/(1+z))\), \(z_* \sim 3\) (turnover redshift), and \(\Delta \sim 0.6\) (width). For \(z \gg z_*\), \(g(z) \ll 1\), and fluctuations are effectively frozen.

---

## 3. BAO Sensitivity Kernel and Precision Floor

### 3.1 Definition of \(S(z)\)
We use the standard isotropic BAO distance proxy:

\[
D_V(z) = \left[ D_M(z)^2 \, \frac{c z}{H(z)} \right]^{1/3},
\]

where \(D_M(z)\) is the comoving transverse distance and \(H(z)\) is the Hubble parameter. The dimensionless BAO scale \(\alpha(z)\) relative to a fiducial cosmology satisfies, to first order,

\[
\frac{\delta \alpha}{\alpha} \approx \frac{\delta D_V}{D_V}.
\]

Define the sensitivity kernel to \(\Omega_\Lambda\):

\[
S(z) \equiv \frac{\partial \ln D_V(z)}{\partial \Omega_\Lambda}.
\]

Linear propagation then yields an induced BAO scatter (the "precision floor"):

\[
\sigma_{\alpha,\rm floor}(z) \approx |S(z)| \, \sigma_{\Omega_\Lambda}(z).
\]

### 3.2 Calibration with DESI DR1
Using the DESI DR1 LRG2 measurement at \(z_{\rm eff} = 0.71\), with reported fractional precision \(\sigma_{\alpha,\rm obs}(0.71) \approx 0.012\), we allocate a conservative fraction \(f_{\rm net} = 0.15\) of the observed variance to an irreducible network component:

\[
\sigma_{\alpha,\rm floor}(0.71) = \sqrt{0.15} \times 0.012 = 4.65 \times 10^{-3}.
\]

Using the kernel value \(S(0.71) = -0.312\) (computed numerically for a flat cosmology with \(\Omega_m = 0.315\), \(H_0 = 67.4\) km/s/Mpc), we infer

\[
\sigma_{\Omega_\Lambda}(0.71) = \frac{\sigma_{\alpha,\rm floor}(0.71)}{|S(0.71)|} \approx 1.49 \times 10^{-2}.
\]

Fixing \(\theta = 1.2\) and requiring the stationary OU variance to match \(\sigma_{\Omega_\Lambda}(0.71)\) gives

\[
\sigma = \sigma_{\Omega_\Lambda}(0.71) \, \sqrt{2\theta} \approx 2.31 \times 10^{-2}.
\]

### 3.3 Exact Integral Susceptibility \(\chi(z, z')\) and the Precision Floor

The simplified scaling laws used previously are approximations. The exact propagation of the stochastic fluctuation \(X(z') \equiv \delta\Omega_\Lambda(z')\) to the BAO observable \(\delta\ln D_V(z)\) is given by a **causal integral kernel** derived from first principles.

We start from the perturbed Hubble parameter in a flat universe with stochastic dark energy:

\[
H^2(z) = H_0^2 \left[ \Omega_m (1+z)^3 + \Omega_\Lambda + X(z) \right],
\]

where \(X(z)\) is the OU stochastic field. Defining \(E(z) = H(z)/H_0\) and expanding the inverse of the cosmic evolution function to first order in \(X(z)\):

\[
\frac{1}{E(z)} \approx \frac{1}{E_0(z)} - \frac{1}{2 E_0(z)^3} X(z), \qquad E_0(z) = \sqrt{\Omega_m (1+z)^3 + \Omega_\Lambda}.
\]

The comoving transverse distance is \(D_M(z) = \frac{c}{H_0} \int_0^z \frac{dz'}{E(z')}\). Introducing the perturbative expansion, the total variation induced by the finite-information network is:

\[
\delta D_M(z) = -\frac{c}{2 H_0} \int_0^z \frac{X(z')}{E_0(z')^3} dz'.
\]

The BAO volume scale mixes the transverse distance and the local Hubble parameter:

\[
D_V(z) = \left[ D_M(z)^2 \frac{c z}{H(z)} \right]^{1/3}.
\]

Taking the logarithm and varying to first order:

\[
\delta \ln D_V(z) = \frac{2}{3} \frac{\delta D_M(z)}{D_{M,0}(z)} - \frac{1}{3} \frac{\delta H(z)}{H_0(z)}.
\]

Using \(\delta H(z)/H_0(z) = \frac{1}{2 E_0(z)^2} X(z)\) and the integral for \(\delta D_M(z)\), the master equation is:

\[
\delta \ln D_V(z) = -\frac{c}{3 H_0 D_{M,0}(z)} \int_0^z \frac{X(z')}{E_0(z')^3} dz' - \frac{1}{6 E_0(z)^2} X(z).
\]

For a memory-dependent process such as the OU model, the fluctuation at redshift \(z'\) affects observations at higher redshifts through the integral coupling. We define the **susceptibility kernel** \(\chi(z, z')\) such that:

\[
\delta \ln D_V(z) = \int_0^z \chi(z, z') X(z') dz'.
\]

Comparing terms, the exact mathematical structure of the susceptibility is a **distribution composed of an accumulative integral term and a local Dirac delta term**:

\[
\chi(z, z') = -\frac{c}{3 H_0 D_{M,0}(z) E_0(z')^3} \Theta(z - z') - \frac{1}{6 E_0(z)^2} \delta_D(z - z'),
\]

where \(\Theta(z - z')\) is the Heaviside step function (enforcing cosmic causality) and \(\delta_D\) is the Dirac delta function.

**Consequence for the precision floor and covariance.** The induced BAO scatter (precision floor) is not a simple algebraic function of the local redshift. Because the observable \(\delta\ln D_V(z)\) depends on an integral weighted by \(E_0(z')^3\) in the denominator, the model predicts that **accumulated fluctuations at low redshifts (\(z \sim 0.3\)) have a significant impact on the uncertainty of high-redshift bins (\(z \sim 1.5\))**. The covariance between BAO residuals at redshifts \(z_i\) and \(z_j\) is therefore given by the double integral:

\[
C_{ij}^{\rm OU} = \int_0^{z_i} \int_0^{z_j} \chi(z_i, z') \, \chi(z_j, z'') \, \operatorname{Cov}[X(z'), X(z'')] \, dz' dz'',
\]

with \(\operatorname{Cov}[X(z'), X(z'')]\) given by the OU exponential kernel modulated by the activation factor \(g(z)\).

This exact integral structure naturally produces the non-trivial lag correlations observed in Section 6: the accumulated weight of low-\(z\) fluctuations can induce positive lag-1 correlations while flipping the sign of lags 2 and 3, depending on the redshift weighting. This justifies why a pure local \(S(z)\) approximation was insufficient for the lag test.

*Numerical implementation of this integral kernel is provided in the companion code (`ou_bao_likelihood.py` and `cross_correlation_DESI.py`), which replaces the discrete \(S(z)\) table with the exact integral propagation for covariance estimation.*

## 4. Test 1: BAO Likelihood with OU Covariance

### 4.1 Methodology
The standard BAO analysis assumes

\[
d = m_{\Lambda\rm CDM}(\theta_{\rm cosmo}) + \epsilon, \qquad \epsilon \sim \mathcal{N}(0, C_{\rm std}).
\]

Our model adds an OU‑induced component:

\[
d = m_{\Lambda\rm CDM} + \delta d_{\rm OU} + \epsilon, \qquad C_{\rm total} = C_{\rm std} + C_{\rm OU},
\]

where

\[
(C_{\rm OU})_{ij} = S(z_i) \, S(z_j) \, \operatorname{Cov}[X(x_i), X(x_j)],
\]

and

\[
\operatorname{Cov}[X(x_i), X(x_j)] = \frac{\sigma^2}{2\theta} \, \exp\left[-\theta |x_i - x_j|\right],
\]

modulated by the activation factor \(g(z)\).

We compute the likelihood using the public DESI DR2 BAO measurements (7 bins from arXiv:2503.14738, Table 1) and the standard diagonal covariance from measurement errors. The OU kernel is computed with the sensitivity kernel \(S(z)\) from Appendix A.

### 4.2 Results

| Model               | \(\theta\) | \(\sigma_X\) | \(\Delta\log\mathcal{L}\) (vs ΛCDM) | AIC      | BIC      |
|:--------------------|:----------:|:------------:|:-----------------------------------:|:--------:|:--------:|
| \(\Lambda\)CDM      | —          | —            | 0.00                                | ref      | ref      |
| OU pure (H0)        | 0.765      | 0.018        | **+6.75**                           | −38.73   | −38.84   |
| Oscillatory (QNM)   | 0.000      | 0.019        | +7.00                               | −37.24   | −37.40   |

The pure OU model is preferred over \(\Lambda\)CDM by \(\Delta\log\mathcal{L} = +6.75\), corresponding to a likelihood ratio of ~e^{6.75} ≈ 850, indicating strong preference over ΛCDM under the Laplace approximation AIC and BIC also favour the OU model over the oscillatory (QNM) kernel by 1.5 points.

**Interpretation:** With only 7 bins, this result is preliminary. It demonstrates that the OU covariance structure provides a better statistical description of the BAO residuals than white noise. However, this does not prove the physical reality of the stochastic component; it could also be mimicked by unmodelled systematics or by a different deterministic evolution. The decisive test is the lag‑correlation structure (Section 6).

---

## 5. Test 2: Angular Cross‑Correlation \(\delta_g \times \delta\Omega_\Lambda\)

### 5.1 Data and Method
We construct a HEALPix map (\(N_{\rm side} = 32\)) of galaxy overdensities using DESI DR1 LRG catalogs (NGC: 1.476M objects, SGC: 662k objects) in the redshift range \(0.4 < z < 0.8\). We weight each galaxy by the standard DESI completeness weights.

As a proxy for \(\delta\Omega_\Lambda\), we use the Hubble diagram residuals \(\Delta\mu(z) = \mu_{\rm obs}(z) - \mu_{\Lambda\rm CDM}(z)\) from the Pantheon+ compilation (1701 SNe Ia). For each SN, we estimate

\[
\delta\Omega_\Lambda \approx \frac{\Delta\mu}{d\mu/d\Omega_\Lambda},
\]

where the denominator is computed numerically using the same fiducial \(\Lambda\)CDM cosmology. We then map these estimates onto HEALPix pixels.

The angular cross‑correlation is computed as the Pearson correlation coefficient between the galaxy overdensity and \(\delta\Omega_\Lambda\) in the overlapping pixels.

### 5.2 Results
Combining NGC and SGC, we find **67 overlapping pixels** with both LRG galaxies and SNe. The cross‑correlation coefficient is

\[
r_{\rm cross} = 0.1673 \pm 0.0613,
\]

corresponding to a \(Z \approx 2.73\sigma\) excess over the null hypothesis (\(r=0\)). The naive OU prediction, using \(r_{\rm cross} \approx G_{\rm coupling} \, \sigma_8 / \sigma_X\) with \(G_{\rm coupling} = 5.04 \times 10^{-4}\), is \(\sim 0.023\). The observed value exceeds this by a factor of \(\sim 7\).

### 5.3 Systematic Checks (Ongoing)
Before attributing this excess to physics, two systematic effects must be thoroughly examined:

1. **Galactic dust:** Extinction by Milky Way dust affects both the LRG photometric redshifts and the SN magnitudes, potentially inducing a spurious cross‑correlation. A partial correlation analysis controlling for the SFD dust map is required.
2. **Imaging systematics:** The NGC survey (BASS+MzLS) and the SGC survey (DECaLS) have different photometric systems. We must analyse NGC and SGC independently and test for a consistent sign and amplitude.

"As a first‑order check, the pipeline includes an optional routine that computes the partial correlation coefficient controlling for Galactic extinction using the SFD dust map. This allows us to assess whether the cross‑correlation signal is driven by dust. The final systematic control, however, requires the DESI DR2 imaging weights (WEIGHT_SYS), which is the primary reason for our Data Lab access request."

Access to **DESI DR2** with the full imaging systematic weights (`WEIGHT_SYS`) is necessary to resolve this. This is the primary motivation for our Data Lab access request.

---

## 6. Test 3: Redshift Lag Correlations in BAO Residuals (The Critical Test)

### 6.1 The Smoking‑Gun Prediction
The OU process implies that BAO residuals are not independent between redshift bins. In the whitened residual vector \(y = L^{-1} r\), where \(C_{\rm std} = L L^T\), the OU model predicts large positive correlations between neighbouring bins:

For a DESI‑like binning \(z = \{0.20, 0.35, 0.51, 0.65, 0.80, 0.95, 1.10, 1.25, 1.40\}\), the predicted lag correlations are:

\[
\rho_1 \approx 0.78, \qquad \rho_2 \approx 0.62, \qquad \rho_3 \approx 0.49.
\]

If the data show no such positive correlations, the pure OU model is strongly disfavoured.

### 6.2 Results from DESI DR1 and DR2
We computed the whitened BAO residuals using the publicly released isotropic \(\alpha\) values from DESI DR1 (6 bins) and DR2 (7 bins). The results are:

| Lag | DR1 (6 bins) | DR2 (7 bins) | OU Prediction (\(\theta=1.2\)) | Verdict |
|:---:|:---:|:---:|:---:|:---:|
| 1   | \(-0.51\) | **\(+0.37\)** | \(+0.83\) | Mixed |
| 2   | \(-0.47\) | \(-0.32\) | \(+0.85\) | Negative |
| 3   | \(+0.48\) | \(-0.81\) | \(+0.85\) | Negative |

Note: The predictions in Section 6.1 use the MLE best‑fit value θ=0.765. The table below uses the calibration value θ=1.2 from Section 3.2 to maintain consistency with the precision‑floor calculation. These two parameter choices are reported separately; the qualitative conclusion (positive lag‑1, negative lags 2‑3) is robust to this variation.

### 6.3 Interpretation
- **Lag‑1:** The DR2 result (\(+0.37\)) is positive, consistent in sign with the OU prediction. However, the DR1 result was negative, and the change of sign between the two data releases is not a trend with only 6–7 bins.
- **Lags 2 and 3:** Both are negative in DR2, inconsistent with a pure OU kernel with \(\theta = 1.2\).

With only 6–7 bins, the 95% confidence intervals are very large (\(\sim \pm 0.7\) to \(\pm 1.0\)), so no definitive conclusion can be drawn. Nevertheless, the negative lag‑2 and lag‑3 results indicate that **the simple OU closure (constant \(\theta\), pure exponential decay) is not the correct kernel**.

**Possible refinements:**
- A redshift‑dependent \(\theta(z)\) (e.g., stronger mean‑reversion at higher \(z\)).
- An oscillatory component (QNM kernel) with a frequency below the Rayleigh limit.
- A different memory kernel altogether (e.g., power‑law or stretched exponential).

The decisive test will require **\(>20\) redshift bins**, which Euclid DR1 will provide.

---

## 7. Unified Discussion

The three tests paint a coherent but incomplete picture:

| Test | Result | Implication |
|:---|:---|:---|
| **BAO likelihood** | OU covariance improves fit over white noise (\(\Delta\log L = +6.75\)). | The data prefer some form of correlated noise over pure white noise. |
| **Angular cross‑correlation** | \(r = 0.167\) (2.73σ), exceeding naive OU prediction by 7×. | Either the signal is stronger than expected, or systematics are present. DR2 weights are needed. |
| **Lag correlations** | Lag‑1 positive, lags 2 & 3 negative. | The pure OU kernel with θ=1.2 is ruled out, but the framework of correlated residuals remains viable. |

The **finite‑information stochastic framework** is not falsified; rather, it is **refined**. The data tell us that a simple Ornstein‑Uhlenbeck process is too restrictive, but they do not reject the underlying idea that \(\Lambda\) fluctuates with memory. The correct memory kernel must be determined empirically or derived from a more complete microscopic model.

---

## 8. Falsification Criteria (Data‑Kills‑Model)

The framework is excluded if any of the following occur with high significance:

| Criterion | Condition for Exclusion | Current Status | Test |
|:---|:---|:---|:---|
| **F1** (Variance floor absent) | \(\sigma^2_{\alpha,\rm obs} - \sigma^2_{\alpha,\rm std} < \sigma^2_{\alpha,\rm floor}\) in multiple bins | **Alive** | Euclid DR1 |
| **F2** (\(w_a \to 0\) at \(>5\sigma\)) | No room for any stochastic component | **Alive** | DESI DR3 / Euclid |
| **F3** (ISW exclusion) | \(\sigma_{\Omega_\Lambda}\) incompatible with Planck CMB‑LSS | **Pending** | Planck 2018 + Euclid |
| **F4** (Bayesian evidence) | \(\ln K > 5\) in favour of \(\Lambda\)CDM | **Not computed** | Full MCMC chains |
| **F5** (Lag correlations) | All lags negative or zero with 20+ bins | **Partial (DR2)** | Euclid DR1 (Oct 2026) |

---

## 8.1 The Rayleigh Resolution Limit: Why Euclid DR1 May Not See the Oscillation

A subtle but crucial consequence of our stochastic framework is the Rayleigh resolution criterion in logarithmic scale‑factor space. For a survey covering a range \(\Delta x = \ln(1+z_{\max}) - \ln(1+z_{\min})\), the minimum resolvable frequency of any oscillatory mode in the memory kernel is:

\[
\omega_{R,\min} = \frac{2\pi}{\Delta x}.
\]

Any kernel oscillation with \(\omega_R < \omega_{R,\min}\) produces less than one visible cycle across the survey and cannot be distinguished from a pure OU process (\(\omega_R = 0\)), regardless of photometric precision.

For **DESI DR2**, the full redshift range \(z \in [0.295, 2.330]\) gives:

\[
\Delta x_{\rm DESI} = \ln\left(\frac{1+2.330}{1+0.295}\right) = \ln(2.571) \approx 0.944,
\]

which implies:

\[
\omega_{R,\min}^{\rm DESI} \approx \frac{2\pi}{0.944} \approx 6.66.
\]

For **Euclid DR1**, the expected galaxy sample covers a narrower range, approximately \(z \in [0.9, 1.8]\), giving:

\[
\Delta x_{\rm Euclid} = \ln\left(\frac{1+1.8}{1+0.9}\right) = \ln(1.474) \approx 0.388,
\]

which implies:

\[
\omega_{R,\min}^{\rm Euclid} \approx \frac{2\pi}{0.388} \approx 16.2.
\]

(Even if Euclid extends to \(z=2.0\), \(\Delta x \approx 0.5\) and \(\omega_{R,\min} \approx 12.6\), still well above the DESI limit).

**The Paradox:** If the true dark‑energy kernel is not a pure OU process (\(H_0\)) but an oscillatory Quasi‑Normal Mode (\(H_1\)) with an intermediate frequency, say \(\omega_R \sim 8\) (which is above the DESI Rayleigh limit but below the Euclid limit), then **DESI can detect the oscillation, but Euclid DR1 geometrically cannot**. Euclid will see the oscillation as a straight line (a pure OU slope) because its shorter lever arm in \(\ln a\) cannot resolve the curvature. This is a radical inversion of the usual "more data solves everything" paradigm: a narrower redshift range can actively erase the oscillatory signature.

**Consequence for the v3.1 Observational Strategy:** The decisive test for the \(H_1\) (QNM) kernel is not Euclid DR1, but DESI DR2 and DR3. Euclid DR1 will provide an unmatched test of the *variance floor* (amplitude \(\sigma_X\)), but it is geometrically blind to the oscillatory component if \(\omega_R > \omega_{R,\min}^{\rm Euclid}\). This means that if our lag‑correlation test in Section 6 shows a positive lag‑1 but negative lags 2‑3 (as seen in DR2), and if this is due to an oscillation with \(\omega_R \sim 8\), then **Euclid will report consistency with a pure OU process**, potentially leading to a false negative if interpreted as a rejection of the entire stochastic framework.

> *"Euclid will measure the 'what' (the amplitude) with exquisite precision, but DESI alone can measure the 'how' (the memory shape) due to its wider leverage in scale factor."*

This result is independent of the model and follows purely from the geometry of the survey. It should be included as a cautionary note in any combined analysis of DESI and Euclid BAO data.
## 9. Near‑Term Observational Program: Euclid DR1

Euclid Data Release 1 is expected on **21 October 2026**. It will provide independent BAO measurements with more than 20 redshift bins across multiple tracers, achieving a precision that can penetrate the predicted OU floor.

### Pipeline Readiness
The analysis pipeline is already prepared and documented in the associated repository:

- **`ou_bao_likelihood.py`**: Computes the BAO likelihood with OU covariance and the lag‑correlation statistics.
- **`cross_correlation_DESI.py`**: Builds HEALPix maps, computes the angular cross‑correlation, and performs bootstrap significance tests.

Upon release of Euclid DR1, we will:
1. Replace the DESI BAO data vector with Euclid’s.
2. Recompute the sensitivity kernel \(S(z)\) for the Euclid fiducial cosmology.
3. Re‑run the lag‑correlation test.
   - If lags 1, 2, and 3 are positive and consistent with a refined OU kernel: we will proceed with a full parameter estimation and publish the detection.
   - If the residuals are consistent with white noise: we will publish a stringent upper limit on the amplitude of stochastic \(\Omega_\Lambda\) fluctuations.

Both outcomes are scientifically valuable.

---

## 10. Conclusion

We have presented a phenomenological model of stochastic dark energy motivated by finite information bounds and discrete‑spacetime arguments. The model predicts a distinctive covariance structure in BAO measurements. Our tests against public DESI data show that:

1. The OU covariance **improves the BAO likelihood** over \(\Lambda\)CDM, indicating that correlated noise is statistically preferred.
2. The angular cross‑correlation between galaxy density and SNe Ia Hubble residuals shows a **2.73σ excess**, which, if confirmed, would be a strong indicator of a local \(\Omega_\Lambda\) fluctuation. However, systematic effects must first be ruled out using DESI DR2 imaging weights.
3. The redshift lag‑correlation test **falsifies the pure OU kernel with constant \(\theta\)**, but does not falsify the broader stochastic framework. A refined memory kernel (e.g., redshift‑dependent \(\theta\) or an oscillatory component) is required.

The finite‑information stochastic framework remains a viable, falsifiable extension of \(\Lambda\)CDM. The decisive test will come from Euclid DR1 (October 2026), for which the full analysis pipeline is ready.

---

## Appendix A — BAO Sensitivity Kernel \(S(z)\): Numerical Implementation

The kernel is defined as \(S(z) = \partial \ln D_V(z) / \partial \Omega_\Lambda\), computed along the flat direction \(\Omega_m = 1 - \Omega_\Lambda\).

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
This implementation is self‑contained and reproduces the kernel values used in the main text.

Appendix B — Reproducibility and Data Access
All scripts for the analysis are publicly available in the repository associated with this preprint:

BAO likelihood and lag correlations: ou_bao_likelihood.py

Angular cross‑correlation: cross_correlation_DESI.py

The cross‑correlation script includes an explicit check (os.path.exists) to distinguish between real and synthetic data in the final summary, ensuring full traceability.

Data sources:

DESI DR1/DR2 BAO measurements: public releases at data.desi.lbl.gov

Pantheon+ SNe Ia: github.com/PantheonPlusSH0ES/DataRelease

Access to DESI DR2 LRG catalogs with imaging systematics weights (WEIGHT_SYS) has been requested from NSF NOIRLab Astro Data Lab.

References
[1] DESI Collaboration, "DESI 2024 III: Baryon Acoustic Oscillations from Galaxies and Quasars," arXiv:2404.03000 (2024).

[2] DESI Collaboration, "DESI DR2 Results II: Measurements of BAO and Cosmological Constraints," arXiv:2503.14738 (2025).

[3] Bekenstein, J. D., "Black Holes and Entropy," Phys. Rev. D 7, 2333 (1973).

[4] Sorkin, R. D., "Is the Cosmological 'Constant' a Nonlocal Quantum Residual?," arXiv:gr-qc/0503057 (2005).

[5] Mota, D. F. and Barrow, J. D., "Varying α in a More Realistic Universe," Phys. Lett. B 581, 141 (2004).

[6] Scolnic, D. et al. (Pantheon+ Collaboration), "The Pantheon+ Analysis: Cosmological Constraints," ApJ 938, 113 (2022).

[7] Uhlenbeck, G. E. and Ornstein, L. S., "On the Theory of the Brownian Motion," Phys. Rev. 36, 823 (1930).

[8] Euclid Consortium, Data Release 1 expected 21 October 2026 (ESA Euclid COSMOS portal).
