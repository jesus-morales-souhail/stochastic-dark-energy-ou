# Constraints on Stochastic Dark Energy from DESI DR2 BAO: A Null Result and the $10^{56}$ Amplification Bottleneck

| | |
|:--|:--|
| Author | Jesús Morales Souhail |
| ORCID | [0009-0000-7637-1818](https://orcid.org/0009-0000-7637-1818) |
| Email | jmskjym@gmail.com |
| Affiliation | Independent researcher |
| Date | v2, October 2026 (v1: July 2026) |
| Status | Preprint |
| Code | [stochastic-dark-energy-ou](https://github.com/jesus-morales-souhail/stochastic-dark-energy-ou) |
| License | MIT (code); CC BY 4.0 (text) |

---

## Abstract

I test whether public DESI DR2 baryon acoustic oscillation (BAO) summary statistics require a stationary stochastic component in the late-time dark-energy sector, on top of a smooth background. Fluctuations are modelled as an Ornstein–Uhlenbeck (OU) process in logarithmic scale factor $x=\ln a$, with a quasi-normal-mode (QNM) oscillatory extension. Using a Gaussian likelihood with an additive residual covariance kernel, the best fit drives the stochastic amplitude to zero, and the profile likelihood over $\theta$ ($\Delta\ln\mathcal{L}=-1.92$, full measurement covariance) gives the 95% confidence upper limit

$$\sigma_X < 2.6\times 10^{-2}$$

under the stated phenomenological kernel and public BAO summary statistics (official Gaussian BAO covariance projected to isotropic α in the baseline pipeline). Nested extensions that free $\{w_0,w_a,\sigma_X,\theta\}$ are not preferred over a smooth background ($\Delta\mathrm{AIC}\approx +4$).

Separately, a coherent tachyonic growth model with the correct rank-1 covariance induced by a single growing mode is strongly disfavoured: $\Delta\ln\mathcal{L}\approx -12.7$ relative to $\Lambda$CDM ($\Delta\chi^2\approx +25$). The two published implementations of Everpresent Λ, a zero-mean Λ with order-one fluctuations, do not fit the same BAO data together with the CMB: with r_d computed for each realization, small dark energy around recombination and ω_m within ±2σ of Planck, none fits as well as ΛCDM (0/20000 and 0/500000 realizations).

At the DESI ceiling the OU residual is invisible to Planck 2018 (TT at ℓ < 30 and lensing; 95th percentile of Δχ² ≤ 0.067), so the CMB gives no leverage on σ_X, and its only imprint, a scatter σ(H0) ≈ 0.48 km/s/Mpc, is 12 times smaller than the H0 tension. A percent-level harmonic residual is likewise allowed by DESI and invisible in the CMB (§3.5).

I then quantify an amplification bottleneck: a Sorkin–Bekenstein Poisson seed with horizon-area counting, $\sigma_0\sim 10^{-61}$ (four-volume counting gives order-one fluctuations instead — Everpresent Λ, a different model) lies $\sim 10^{56}$ below Euclid-scale residual amplitudes ($\sim 10^{-5}$) and $\sim 10^{59}$ below the DESI residual ceiling measured here ($2.6\times 10^{-2}$). Audited linear amplifiers are excluded, short by many decades, or redefine the seed. Machine scans of three allowed non-linear cards show that late-time freeze-out alone yields unit gain, a soft double-well avalanche yields gain $\sim 2$, and only a local effective count $N_{\mathrm{eff}}\sim 1.4\times 10^{3}$ (DESI) / $10^{10}$ (Euclid) reaches target amplitudes—by redefining the microphysical seed, not by free amplification.

Gravitational slip $\gamma=\Phi/\Psi$ is the correct operator for anisotropic leakage past volume-preserving (SDiff) structure, and with the BAO-bounded residual it is marginal, not starved: full anisotropy gives at most $\lvert\gamma-1\rvert\approx 5\times 10^{-2}$ at $z=0.5$, about $3.4\times$ below the current measurement (Maus et al.) and at the level of the Sakr et al. forecast floor ($\approx 1.0\times$), falling to $2.4\times$ and $4.7\times$ below that floor at $z=1$ and $1.5$. The σ_X → slip map is an order-of-magnitude estimate: it holds for δ_m ≈ 1 (quasi-linear scales) and does not extrapolate to δ_m ≪ 1; ε < 1 or z ≳ 1 push it below reach.

I claim limits, exclusions, and a map of open hypotheses—not a detection of vacuum discreteness or a new dark-energy fluid.


---

## 1. Introduction

### 1.1 Motivation

DESI BAO measurements have sharpened late-time distance constraints and, with CMB and supernovae, have renewed interest in dynamical dark energy (DESI Collaboration, arXiv:2503.14738 and companions). Independently of whether the *background* equation of state departs from $w=-1$, a separate question is whether *residuals* about a smooth expansion history require a stochastic dark-energy component.

A finite information bound on the observable universe ($N\sim 10^{122}$ Planck units) has long motivated order-of-magnitude vacuum fluctuations $\delta\Lambda\sim 1/\sqrt{N}$ (Sorkin-type arguments in unimodular / causal-set settings). That seed is many orders of magnitude below BAO precision. The programme I follow is therefore:

1. **Data:** constrain the effective residual amplitude under a minimal, falsifiable kernel.
2. **Models:** exclude coherent amplification mechanisms with the correct covariance.
3. **Theory hygiene:** state which amplifiers are closed and which open cards remain—without smuggling free parameters as “microphysics.”

Everpresent Λ (Ahmed et al. 2004; Zwane et al. 2018), in which the same counting applied to the four-volume gives order-one fluctuations, is a different model, tested directly in §4.2.

### 1.2 Scope and non-claims

This paper and the associated public repository address (1)–(3) using public BAO summary statistics and open analysis code. I do **not** claim:

- detection of Planck-scale vacuum noise;
- that Euclid will see unamplified Sorkin seeds;
- a unique identification of unimodular / SDiff gravity;
- peer-reviewed status of this preprint;
- results from optical laboratory metaphors (those live in a separate exploratory repository; see Section 8).

### 1.3 Reading guide

| Section | Content |
|:--------|:--------|
| 2 | Data and likelihood |
| 3 | Stationary OU / QNM null, $\sigma_X$ limit and CMB consistency (§3.5) |
| 4 | Exclusions: coherent tachyonic growth (§4.1) and Everpresent Λ (§4.2) |
| 5 | Amplification bottleneck and route scans |
| 6 | Slip / SDiff structural remark (amplitude honesty) |
| 7 | Limitations |
| 8 | Software, repositories, and claim boundaries |
| 9 | Conclusions |

---

## 2. Data and method

### 2.1 Dataset

I use **public DESI DR2 BAO summary statistics** (arXiv:2503.14738): isotropic BAO $\alpha$ measurements in seven effective redshift bins, with published measurement uncertainties. The baseline pipeline adopts the official DESI Gaussian BAO 13×13 covariance projected to isotropic α (block-diagonal across bins) in the residual likelihood. Off-diagonal covariances between bins are absent in the public product; they are not required to establish the *null preference* of the stochastic extension under the kernels tested.

**Important separation:** DESI analyses that combine BAO with CMB+SN often prefer dynamical $w_0$–$w_a$. The reference for $\alpha=1$ in my BAO-only residual analysis is the fiducial $\Lambda$CDM cosmology of the published $\alpha$ products ($w_0=-1$, $w_a=0$). As a check, a CPL fit to the same data lands at $w_0\approx -1.001$, $w_a\approx +0.014$ (`results/eos_cpl_desi_dr2/`, projected covariance, figure6 α), indistinguishable from that reference. The stochastic question is about *residuals*, not a re-fit of the full multi-probe DE equation of state.

**Secondary dataset (Section 3.5 only).** Planck 2018 TT at 2 ≤ ℓ < 30 and the Planck lensing bandpowers (L = 8–400), with cosmological parameters fixed to Planck 2018 except H0, which is solved for each curve to keep θ_* at its ΛCDM value. It is used as a consistency check, not in the σ_X likelihood.

### 2.2 Residual process

Let $X(x)\equiv\delta\rho_{\mathrm{DE}}(x)/\rho_{\mathrm{crit},0}$ ($=\delta\Omega_\Lambda$ in today's units) with $x=\ln a$, so that $\delta\ln\rho_{\mathrm{DE}}=X/\Omega_{\Lambda0}$. The stationary OU process is


$$
\mathrm{d}X=-\theta X\,\mathrm{d}x+\sigma\,\mathrm{d}W_x,\quad \sigma_X^2\equiv\mathrm{Var}(X)=\frac{\sigma^{2}}{2\theta}\quad(\theta>0).
$$


I work with the effective residual amplitude $\sigma_X$ (the stationary standard deviation of $X$, not the diffusion $\sigma$), which enters the BAO residual covariance through the sensitivity kernel $S_i\equiv S(z_i)=\partial\ln D_V(z_i)/\partial\Omega_\Lambda$ (flat), $\delta\alpha_i=S_iX(z_i)$. For bins $i,j$ with logarithmic separation $\Delta x_{ij}$,


$$
C^{\mathrm{OU}}_{ij}=S_iS_j\,\sigma_X^2 e^{-\theta\Delta x_{ij}},\quad C^{\mathrm{QNM}}_{ij}=S_iS_j\,\sigma_X^2 e^{-\theta\Delta x_{ij}}\cos(\omega_R\Delta x_{ij}).
$$


The total covariance is $C_{\mathrm{tot}}=C_{\mathrm{meas}}+C_{\mathrm{stoch}}$. The Gaussian log-likelihood is standard:


$$
-2\ln\mathcal{L}=\ln\det C_{\mathrm{tot}}+\mathbf{r}^{\mathrm{T}}C_{\mathrm{tot}}^{-1}\mathbf{r}+\mathrm{const},
$$


with residual vector $\mathbf{r}$ relative to the smooth background prediction.

### 2.3 Nested models and information criteria

| Model | Stochastic content | Free stochastic params |
|:------|:-------------------|:-----------------------|
| $\Lambda$CDM baseline | $ C_{\mathrm{meas}}$ only | — |
| H0 (OU) | $C_{\mathrm{OU}}(\theta,\sigma_X)$ | 2 |
| H1 (QNM) | $C_{\mathrm{QNM}}(\theta,\sigma_X,\omega_R)$ | 3 |
| Nested CPL+$(\sigma,\theta)$ | background free + residual | see Section 3 |

AIC/BIC differences use the usual penalties on the number of free parameters.

### 2.4 Software

Primary scripts (repository root):

| Script | Role |
|:-------|:-----|
| `scripts/ou_bao_likelihood.py` | OU / QNM residual MLE |
| `scripts/eos_efectiva.py` | CPL background + nested $(\sigma,\theta)$ |
| `scripts/profile_sigma_x_desi.py` | Profile likelihood for $\sigma_X$ |
| `scripts/joint_w0wa_sigma_desi.py` | Joint $\{w_0,w_a,\sigma_X\}$ |
| `scripts/gpe/` | Tachyonic / GPE exclusion pipeline |
| `scripts/amplifier_audit.py` | Amplification audit |
| `scripts/amplification/route{1,2,3}_*.py` | Non-linear cards |
| `scripts/slip_bridge.py` | map $\sigma_X\to\lvert\gamma-1\rvert$ (order of magnitude) |
| `analisis_cmb/realizaciones.py`, `analisis_cmb/estadisticos.py` | CMB consistency: OU realizations through CAMB, $\Delta\chi^2$ against Planck 2018 (§3.5) |
| `analisis_cmb/oscilador/oscilador.py` | Harmonic residual against DESI BAO and Planck 2018 (§3.5) |
| `analisis_cmb/everpresent/bao_sim.py`, `rd_real.py`, `robustez.py` | Everpresent Λ realizations against DESI DR2 BAO (§4.2) |

---

## 3. Results I — Stationary stochastic null

### 3.1 Maximum-likelihood point

Under a fixed smooth background and public BAO summary statistics, with the official DESI measurement covariance in all three likelihoods, the MLE drives the stochastic amplitude to the lower bound of the search range, $\sigma_X = 4.54\times 10^{-5}$ ($\ln\sigma_X^2 = -20$):

| Model | $\theta$ | $\sigma_X$ | $\omega_R$ | $\ln\mathcal{L}$ | $\Delta\ln\mathcal{L}$ vs $\Lambda$CDM | AIC | BIC |
|:------|:-----------|:-------------|:-------------|:------------------|:-------------------------------------------|:----|:----|
| $\Lambda$CDM | — | — | — | $24.8996$ | $0$ (ref) | $-49.80$ | $-49.80$ |
| OU (H0) | $266$ (unconstrained) | $4.54\times 10^{-5}$ (lower bound) | — | $24.8995$ | $-4.0\times 10^{-5}$ | $-45.80$ | $-45.91$ |
| QNM (H1) | $0.043$ (unconstrained) | $4.54\times 10^{-5}$ (lower bound) | $6.93$ (unconstrained) | $24.8995$ | $-3.8\times 10^{-5}$ | $-43.80$ | $-43.96$ |

Relative to $\Lambda$CDM, $\Delta$AIC $= +4.00$ (OU) and $+6.00$ (QNM), $\Delta$BIC $= +3.89$ and $+5.84$: the extra parameters bring no gain in likelihood. With $\sigma_X$ at the bound, $\theta$ and $\omega_R$ are not constrained. Source: [`results/ou_qnm_mle/ou_qnm_mle.json`](../results/ou_qnm_mle/ou_qnm_mle.json) (`scripts/ou_bao_likelihood.py`).

An earlier diagonal-covariance run (old run, [`archive/results/ou_bao_desi_dr2_run.txt`](../archive/results/ou_bao_desi_dr2_run.txt)) gave $\theta = 0.001$, $\sigma_X = 5\times 10^{-5}$ and AIC $-50.03$ / $-48.03$ (OU / QNM).

The data do not require a stationary stochastic residual. The preferred stochastic amplitude is consistent with zero within the resolution of the likelihood surface; values at the bound reflect the search range, not a detection.

### 3.2 Working upper limit

From the residual analysis (profile / likelihood ratio under the OU kernel as implemented in the public pipeline): $\sigma_X < 2.6\times 10^{-2}$ (95% CL). The profile is shown in Fig. 1 ([`figures/profile_sigma_x.png`](../figures/profile_sigma_x.png)).


This is an upper limit on an effective residual amplitude, not a measurement of a Planck-scale seed.

### 3.3 Nested background freedom

When $\{w_0,w_a\}$ are freed together with residual parameters, the stochastic extension is not preferred ($\Delta\mathrm{AIC}\approx +4$). Marginalising over a flexible smooth background drives $\sigma_X$ toward zero: extra residual variance is absorbed by background freedom when the number of BAO bins is small ($N_{\mathrm{bin}}=7$). Fig. 2 ([`figures/joint_w0wa_sigma.png`](../figures/joint_w0wa_sigma.png)) shows the joint fit of $\{w_0, w_a, \sigma_X\}$.

### 3.4 Claim statement (BAO residuals)

> Under the phenomenological OU/QNM residual kernels and public DESI DR2 BAO summary statistics, there is no statistical evidence for stationary stochastic dark-energy fluctuations. The data prefer smooth evolution; $\sigma_X$ is bounded.

### 3.5 CMB consistency (Planck 2018)

I map the OU residual onto a single expansion history, δ ln ρ_DE = X/Ω_Λ0 with X smoothed on Δ = 0.05 in ln a, convert it to w(a) and run CAMB (PPF) with ω_b, ω_c and Σm_ν fixed to Planck 2018 and H0 solved so that 100θ_* = 1.041198. For 200 realizations at σ_X = 2.633×10⁻²:

| θ | Δχ²_TT (ℓ<30) p50 / p95 | Δχ²_φφ p50 / p95 | fraction Δχ² > 4 | H0 p16 / p50 / p84 [km/s/Mpc] |
|:--|:--|:--|:--|:--|
| 0.001 | −1.5×10⁻⁴ / 2.7×10⁻³ | −1.2×10⁻⁴ / 2.3×10⁻³ | 0 | 67.34 / 67.36 / 67.38 |
| 1 | −2.5×10⁻³ / 0.067 | −2.8×10⁻³ / 0.054 | 0 | 66.92 / 67.36 / 67.88 |

On a grid of ten values from 2.506×10⁻² down to 10⁻⁴ the 95th percentile stays below 4 everywhere and scales linearly with σ_X, so Planck does not bound σ_X below the DESI ceiling. The only CMB imprint is a scatter in H0, σ(H0) = 0.48 km/s/Mpc at θ = 1 and 0.020 at θ = 0.001, i.e. 12 and 280 times smaller than the Planck–SH0ES difference (5.68 km/s/Mpc): the residual does not explain the H0 tension.

A deterministic harmonic residual, X = A sin(2π ln a/P + φ) with A ∈ {0.005, 0.01, 0.025}, P ∈ {0.1, …, 2} and four phases, is allowed by DESI DR2 BAO (13 distances, full covariance, Δχ²_BAO < 4) in 58 of 60 cases, and none of the 60 is distinguishable in the CMB (largest Δχ² among the allowed: 0.070 in TT, 0.63 in lensing).

Figure and data: `analisis_cmb/RESULTADO.md`, `analisis_cmb/oscilador/`.

---

## 4. Results II — Model exclusions

### 4.1 Coherent tachyonic growth

A distinct hypothesis is that a dark-energy quantum fluid with negative effective mass undergoes Bogoliubov instability, yielding a coherent growing mode $\sigma_X(t)=\sigma_0 e^{t/t_c}$. The induced BAO residual covariance is **rank-1**, not OU-stationary. Using that covariance in an MLE on the same BAO summary statistics yields

$$
\Delta\ln\mathcal{L}=\ln\mathcal{L}(t_c)-\ln\mathcal{L}_{\Lambda\mathrm{CDM}}\approx -12.7\quad(\Delta\chi^2\approx +25)\ \text{at representative active }t_c;\ \text{preferred }t_c\to\infty.
$$

Any finite collapse / growth time $t_c$ is disfavoured; the acceptable limit is $t_c\to\infty$ (no growth).

**Claim (model kill):**

> A globally coherent tachyonic growing mode with rank-1 residual covariance is incompatible with DESI DR2 BAO residuals under the pipeline stated in the repository.

This exclusion is logically independent of the stationary OU bound: the covariance structure is different, and the correct structure must be used.

### 4.2 Everpresent Λ tested directly

The seed σ_0 ∼ 10⁻⁶¹ uses N ∼ 10¹²² from horizon-area (Bekenstein–Hawking) counting. With the four-volume counting of Everpresent Λ, δΛ tracks H² and the relative amplitude is of order one with no amplifier; the ∼10⁵⁶ gap is specific to the area-counting choice. Zwane et al. (2018) find that Everpresent Λ fits the CMB, BAO and SN data of that time about as well as ΛCDM and eases the low-redshift H0 and Lyα BAO tensions. The order-one alternative to the area-counting seed has been tested directly. In its two published implementations (Zwane et al. 2018, Model 2; Das et al. 2023, Model 1), computing r_d with each realization's own history, requiring small dark energy over z ∈ [1000, 5000] and ω_m within ±2σ of Planck, none of the realizations fits DESI DR2 BAO as well as ΛCDM, nor within Δχ² < 4: 0/20000 for Model 2 (< 1.5×10⁻⁴ at 95 %) and 0/300000 and 0/200000 for Model 1 with α = 0.01 and 0.005 (< 1.0×10⁻⁵, < 1.5×10⁻⁵). Without the CMB conditions, rare realizations pass, all with ω_m far from Planck or order-one dark energy at recombination. As small fluctuations around ΛCDM, DESI DR2 bounds the Zwane et al. amplitude to α < 1.85×10⁻² (μ fixed). These results hold for those phenomenological implementations, with the CMB only through a proxy (|Ω_de| < 0.25), not for the idea of Everpresent Λ in general.

| Implementation | Realizations | Parameters fixed to Planck (r_d, ω_m): pass A / B | Final criterion (c): pass A / B | 95 % upper limit, (c) |
|:--|--:|:--|:--|:--|
| Model 2 (Zwane et al. 2018; α = 0.8824, μ = 0.9804) | 20000 | 1 / 1 (fails the CMB proxy) | 0 / 0 | < 1.5×10⁻⁴ |
| Model 1 (Das et al. 2023), α = 0.01 | 300000 | 0 / 0 (best Δχ² = +17.3) | 0 / 0 | < 1.0×10⁻⁵ |
| Model 1, α = 0.005 | 200000 | 0 / 0 (best Δχ² = +30.5) | 0 / 0 | < 1.5×10⁻⁵ |
| Final criterion (c), all realizations | 520000 | — | **0 / 0** | — |

Conditions of criterion (c): r_d computed from each realization's own expansion history; small dark energy around recombination, $\lvert\Omega_{\rm de}\rvert < 0.25$ over $z\in[1000, 5000]$ (CMB proxy); $\omega_m$ within ±1.6 % (±2σ) of Planck 2018. Data: DESI DR2 BAO, 13 distances with full covariance. Pass A: χ² ≤ 29.84 (ΛCDM with Planck parameters); pass B: Δχ² < 4. Source: [`analisis_cmb/everpresent/EVERPRESENT_vs_DESI.md`](../analisis_cmb/everpresent/EVERPRESENT_vs_DESI.md) §3a, §3c and §3e.

---

## 5. Results III — Amplification bottleneck

### 5.1 Quantified gap (two targets — do not mix)

| Quantity | Value |
|:---------|:------|
| Motivational UV seed (Sorkin / Bekenstein) | $\sigma_0\approx 1.18\times 10^{-61}$ ($H_0=67.4$) |
| DESI residual bound (this work, measured) | $\sigma_X<2.6\times 10^{-2}$ (95% CL)|
| Euclid-scale residual of interest (aspirational) | $\sim 10^{-5}$ |
| $G_{\mathrm{Euclid}}=10^{-5}/\sigma_0$ | $8.5\times 10^{55}\approx 10^{56}$ |
| $G_{\mathrm{DESI}}=2.6\times 10^{-2}/\sigma_0$ | $2.24\times 10^{59}\approx 10^{59}$ |

The headline “$\sim 10^{56}$” always means the **Euclid-scale** target. Using the DESI ceiling in a density ratio and labelling it $10^{56}$ is incorrect by a factor $\sim 2\times 10^{3}$. Both gaps remain unbridgeable by soft linear amplifiers ($r\sim 64$–$68$).

DESI constrains the **effective** residual amplitude after any amplification—not $\sigma_0$, and not a new value of $\rho_\Lambda$. Without a physical map $\sigma_0\to A_0$, “Euclid will see vacuum noise” is not a theorem of this analysis.

### 5.2 Linear amplifier audit

| Class | Schematic gain | Status |
|:------|:---------------|:-------|
| Coherent / tachyonic growth | exponential | Excluded (Section 4) |
| Open-system desqueezing ($r\sim 1.5$) | $e^{2r}\sim 20$ | $\sim 55$ decades short |
| Naive $\sqrt{N}$ boost | $\sim 10^{61}$ | Double-counts unless coherent (→ excluded) |
| Inflation-style $e^{60}$ | $\sim 10^{26}$ | Still short from $\sigma_0$; natural seed $\sim H/M_{\mathrm{Pl}}$ |
| Large phase jump | discrete feature | Tensions with smooth BAO residuals if large |

### 5.3 Three non-linear cards (machine scans)

| Route | Idea | Numerical verdict |
|:------|:-----|:------------------|
| **1. Local $N_{\mathrm{eff}}$** | $\sigma_{0,\mathrm{eff}}=1/\sqrt{N_{\mathrm{eff}}}$, $N_{\mathrm{eff}}\ll N_{\mathrm{BH}}$ | $N_{\mathrm{eff}}\sim 1.44\times 10^{3}$ (DESI), $1.00\times 10^{10}$ (Euclid). Only amplitude card—redefines the seed. |
| **2. Late freeze-out** | $\theta\to 0$ after late exit, $\Delta x\sim\mathcal{O}(1)$ | Freeze/restore gain $=1.000$; $\mathrm{rms}/\sigma=\mathcal{O}(1)$. Dead as amplifier alone. |
| **3. Soft avalanche** | double-well Langevin trigger | Gain $\mathrm{p95}/\sigma\sim 2$; 0/288 jobs with $\mathrm{p95}\ge 10^{-5}$ and $\sigma\le 10^{-8}$. All 288 jobs BAO-safe (max p95 $=2.1\times 10^{-3}$). No free gain. |

Full tables: [`results/amplification_routes/VERDICT.md`](../results/amplification_routes/VERDICT.md) and [`papers/amplification-gap.md`](../papers/amplification-gap.md) §6.3.

### 5.4 Claim statement (amplification)

> **Amplification gap:** $\sigma_0\to 10^{-5}$ is $\sim 10^{56}$ (Euclid target); $\sigma_0\to 2.6\times 10^{-2}$ is $\sim 10^{59}$ (DESI ceiling); linear free amplifiers fail both.

---

## 6. Structural remarks — SDiff and gravitational slip

### 6.1 Geometry as interpretation, not detection

Volume-preserving diffeomorphisms / unimodular structure project out isotropic vacuum stress of the form $T_{\mu\nu}=V(x)g_{\mu\nu}$. That is a candidate explanation of late-time residual smoothness, not a particle detection and not a substitute for the BAO likelihood.

The same unimodular starting point, combined with causal-set discreteness, predicts order-one fluctuations of the global Λ (Everpresent Λ; Ahmed et al. 2004). The projection above removes local terms of the form $V(x)g_{\mu\nu}$; it does not act on the global Λ conjugate to the four-volume.

### 6.2 Option 0: anisotropic gap is real in principle, marginal in amplitude

Anisotropic stress is not of the form $Vg_{\mu\nu}$. Sub-horizon slip satisfies (Newtonian gauge)


$$
\lvert\gamma-1\rvert=\left| \frac{\Phi-\Psi}{\Psi}\right| =2\varepsilon\sigma_X\frac{\rho_{\mathrm{crit},0}}{\rho_m\lvert\delta_m\rvert},
$$


with $\varepsilon\in[0,1]$ the anisotropic fraction of the residual, $\pi_T=\varepsilon\sigma_X\rho_{\mathrm{crit},0}$ (since $X\equiv\delta\rho_{\mathrm{DE}}/\rho_{\mathrm{crit},0}$). With $\sigma_X=2.633\times 10^{-2}$, $\varepsilon=1$, $\delta_m=1$:

| $z$ | max $\lvert\gamma-1\rvert$ | vs Maus $\lvert\gamma-1\rvert\sim 0.17$ | vs forecast floor $\sim 0.05$ |
|:------|:-----------------------------|:------------------------------------------|:--------------------------------|
| 0.5 | $4.95\times 10^{-2}$ | $\sim 3.4\times$ short | $\sim 1.0\times$ (at the floor) |
| 1.0 | $2.09\times 10^{-2}$ | $\sim 8.1\times$ short | $\sim 2.4\times$ short |
| 1.5 | $1.07\times 10^{-2}$ | $\sim 16\times$ short | $\sim 4.7\times$ short |

References: Maus et al., arXiv:2505.20656; Sakr et al., arXiv:2501.07477; runnable map: [`scripts/slip_bridge.py`](../scripts/slip_bridge.py).

**Claim (Slip / Option 0):**

> Slip is the correct *kind* of operator for anisotropic leakage past SDiff. With the BAO-bounded residual amplitude it is marginal, not starved: at ε = 1 it reaches $\lvert\gamma-1\rvert\approx 5\times 10^{-2}$ at $z=0.5$, level with the forecast floor ($\approx 1.0\times$) and $\approx 3.4\times$ below the current measurement. It is not a demonstrated shortcut. The σ_X → slip map is an order-of-magnitude estimate: it holds for δ_m ≈ 1 (quasi-linear scales) and does not extrapolate to δ_m ≪ 1; ε < 1 or z ≳ 1 push it below reach.

I do not implement a homemade Boltzmann hierarchy. Future work, if amplitudes become reachable, should use community codes (e.g. hi_class / MGCAMB) with a pre-specified target. (The CMB consistency check of §3.5 uses CAMB with the residual mapped to a tabulated w(a) in the PPF module; it is not a Boltzmann treatment of anisotropic stress.)

---

## 7. Limitations

1. **Summary statistics:** baseline residual analysis uses published BAO $\alpha$ with the official Gaussian BAO covariance projected to $\alpha$ (block-diagonal across bins). Cross-bin covariances, if released, may change the numerical factor in $\sigma_X$ while preserving the null preference if residuals remain consistent with noise.
2. **Kernel dependence:** bounds are for the OU/QNM residual kernels stated; other residual models need their own covariance and MLE.
3. **Background degeneracy:** with only seven BAO bins, flexible $w(z)$ can absorb residual variance (Section 3.3). Multi-probe fits are a separate analysis. The CMB check of §3.5 is a consistency test with parameters fixed except H0, not a multi-probe fit.
4. **Cross-correlation (DESI galaxies × SN residuals):** an exploratory angular correlation test exists (`archive/scripts/cross_correlation_DESI.py`) but is preliminary (systematics not fully controlled; DR1 weights). It is not part of the primary claim set of this preprint.
5. **Sorkin seed:** used as motivational UV scale (horizon-area counting) for the amplification discussion, not as a fitted parameter. See §4.2 for Everpresent Λ.
6. **Independent research:** this work has not undergone peer review at the time of writing.
7. **CMB check:** parameters other than H0 fixed to Planck 2018 (no MCMC); smoothing scale Δ introduced for the check; reading (a), δ ln ρ_DE = X/Ω_Λ0 (reading (b) would scale the amplitude by 1.46); halofit with tabulated w(a) in PPF not verified.
8. **Everpresent Λ test:** the results of §4.2 hold for the two published phenomenological implementations, with the CMB only through a proxy ($\lvert\Omega_{\rm de}\rvert < 0.25$ over $z\in[1000,5000]$), not for the idea of Everpresent Λ in general.

---

## 8. Software, repositories, and claim boundaries

### 8.1 Primary scientific repository

**https://github.com/jesus-morales-souhail/stochastic-dark-energy-ou**

Contains the BAO pipelines, amplification audits, slip bridge, figures, and technical notes supporting Sections 2–6.

### 8.2 Exploratory repository (separate)

Method-hygiene digressions, optical analogies, and other non-cosmological scale–operator checks are not part of the peer-review claim set:

**https://github.com/jesus-morales-souhail/stochastic-de-exploratory-notes**

Theory (if the null is not the end of the story)—minimal model as simple as $\Lambda$, residual scale uniqueness—lives here:

**https://github.com/jesus-morales-souhail/measurable-stochastic-vacuum**  
Start with the README there, then `papers/r1_kernel/NOTE_uniqueness_residual_grain.md` and `papers/core/VERIFIED_RESULTS.md` (`pytest -q`).

That theory repo does not alter the empirical claims of this preprint; it uses them only as a posteriori tests.

| In primary repo (claims) | In exploratory repo (not claims) |
|:-------------------------|:---------------------------------|
| BAO OU/QNM limits | Car–drone pupil diffraction demo |
| Tachyonic exclusion | Tesseract /$B_{4}$ optical combinatorics |
| Amplification audit + routes 1–3 | Wavefront-shaping vs OU vacuum essays |
| Slip amplitude map | Superoscillation energy-tax demos |
| SDiff as interpretation | Maxwell→device→$M$ lab notes |

**Rule:** if a result does not constrain cosmology at the correct scale and operator, it does not enter the abstract or Sections 3–6.

### 8.3 Reproducibility

```bash
git clone https://github.com/jesus-morales-souhail/stochastic-dark-energy-ou.git
cd stochastic-dark-energy-ou
pip install -r requirements.txt
python scripts/ou_bao_likelihood.py
python scripts/profile_sigma_x_desi.py
python scripts/amplifier_audit.py
python scripts/amplification/run_all_routes.py
python scripts/slip_bridge.py
# CMB consistency: requires CAMB (pip install camb) and the Planck files listed in analisis_cmb/README.md
cd analisis_cmb && python realizaciones.py --halofit takahashi && python estadisticos.py resultados/takahashi/sigma_2.633e-02
```

---

## 9. Conclusions

1. **Null residual:** DESI DR2 BAO summary statistics do not favour stationary OU/QNM stochastic dark energy; $\sigma_X<2.6\times 10^{-2}$ (95% CL) under the stated kernel.
2. **Model kill:** coherent tachyonic growth with rank-1 covariance is excluded ($\Delta\chi^2\approx +25$).
3. **Everpresent Λ:** its two published implementations do not fit DESI DR2 BAO together with the CMB conditions (0/20000 and 0/500000 realizations; §4.2).
4. **CMB consistency:** at the DESI ceiling the OU residual is invisible to Planck 2018 (p95 Δχ² ≤ 0.067) and does not explain the H0 tension (σ(H0) ≈ 0.48 km/s/Mpc); a percent-level harmonic residual is equally allowed by DESI (58/60) and invisible in the CMB.
5. **Amplification gap:** $\sigma_0\to 10^{-5}$ is $\sim 10^{56}$ (Euclid target); $\sigma_0\to 2.6\times 10^{-2}$ is $\sim 10^{59}$ (DESI ceiling); linear free amplifiers fail both.
6. **Slip honesty:** anisotropic leakage past SDiff is a real structural crack in principle and a marginal channel in practice — at best level with forecast floors ($\approx 1.0\times$ the Sakr et al. floor, $\approx 3.4\times$ below Maus et al., at $z=0.5$). The σ_X → slip map is an order-of-magnitude estimate: it holds for δ_m ≈ 1 (quasi-linear scales) and does not extrapolate to δ_m ≪ 1; ε < 1 or z ≳ 1 push it below reach.
7. **Scientific posture:** the strength of this programme is what it closes—limits, exclusions, and an honest map of open theory prices—not a promised detection.

Future data (Euclid BAO multi-bin products; full DESI covariances) can tighten $\sigma_X$ with the same residual likelihood. The CMB gives no additional leverage at this amplitude (Section 3.5); BAO remains the route. Future theory work must pay the price of Routes 1 or 3 (or an inflation-style seed redefinition)—not invent linear gain on white noise.

---

## Version history

- **v1 (July 2026):** σ_X bound quoted as the grid node 2.5×10⁻²; derived numbers (amplification gaps, N_eff, slip, appendix bounds) computed with an earlier working target of 1.5×10⁻⁴; CPL fits with diagonal covariance.
- **v2 (October 2026):** bound = root of the profile, 2.6×10⁻² (2.633×10⁻²); all derived numbers regenerated; CPL fits with the official DESI covariance; CMB consistency and harmonic residual (§3.5) and Everpresent Λ (§4.2) added; slip map corrected (π_T = εσ_X ρ_crit,0).
- Audit and before/after records: [`analisis_cmb/AUDITORIA.md`](../analisis_cmb/AUDITORIA.md); blocks [A](../analisis_cmb/bloque_A/COMPARACION.md), [B1](../analisis_cmb/bloque_B/COMPARACION_B1.md), [B2](../analisis_cmb/bloque_B/COMPARACION_B2.md), [C](../analisis_cmb/bloque_C/COMPARACION_C.md), [D](../analisis_cmb/bloque_D/COMPARACION_D.md), [E](../analisis_cmb/bloque_E/COMPARACION_E.md).

---

## Acknowledgments

This work uses publicly released DESI DR2 BAO products and community literature on gravitational slip and forecasts (Maus et al.; Sakr et al.). Any errors are mine alone.

---

## Data availability

- DESI DR2 BAO: [data.desi.lbl.gov](https://data.desi.lbl.gov/public/) · [arXiv:2503.14738](https://arxiv.org/abs/2503.14738)
- Analysis code: [stochastic-dark-energy-ou](https://github.com/jesus-morales-souhail/stochastic-dark-energy-ou)
- Numerical amplification tables: [`results/amplification_routes/`](../results/amplification_routes/)

---

## References

1. DESI Collaboration, *DESI DR2 Results II: Measurements of Baryon Acoustic Oscillations*, [arXiv:2503.14738](https://arxiv.org/abs/2503.14738).
2. Maus et al., *Joint 3D clustering and galaxy × CMB-lensing with DESI DR1*, [arXiv:2505.20656](https://arxiv.org/abs/2505.20656).
3. Sakr, Zheng & Casas, forecasts on anisotropic stress /$\eta$, [arXiv:2501.07477](https://arxiv.org/abs/2501.07477).
4. Ma & Bertschinger, *Cosmological Perturbation Theory in the Synchronous and Conformal Newtonian Gauges*, Astrophys. J. Suppl. (1995), [arXiv:astro-ph/9506072](https://arxiv.org/abs/astro-ph/9506072).
5. Clifton, Ferreira, Padilla & Skordis, *Modified Gravity and Cosmology*, Phys. Rep. (2012).
6. Ahmed, M., Dodelson, S., Greene, P. B., & Sorkin, R. D., "Everpresent Λ," Phys. Rev. D 69, 103523 (2004), arXiv:astro-ph/0209274.
7. Zwane, N., Afshordi, N., & Sorkin, R. D., "Cosmological tests of Everpresent Λ," Class. Quantum Grav. (2018), doi:10.1088/1361-6382/aadc36, arXiv:1703.06265.
8. Sorkin, R. D., "Is the cosmological 'constant' a nonlocal quantum residue of discreteness of the causal set type?," AIP Conf. Proc. 957, 142 (2007), arXiv:0710.1675.
9. Morales Souhail, J., this repository and companion notes (October 2026).

