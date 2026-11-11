# Constraints on Stochastic Dark Energy from DESI DR2 BAO: A Null Result and the $10^{56}$ Amplification Bottleneck

| | |
|:--|:--|
| Author | Jesús Morales Souhail |
| ORCID | [0009-0000-7637-1818](https://orcid.org/0009-0000-7637-1818) |
| Email | jmskjym@gmail.com |
| Affiliation | Independent researcher |
| Date | v2.2.1, November 2026 (v2.2: November 2026; v2: October 2026; v1: July 2026) |
| Status | Preprint |
| Code | [stochastic-dark-energy-ou](https://github.com/jesus-morales-souhail/stochastic-dark-energy-ou) |
| License | MIT (code); CC BY 4.0 (text) |

---

## Abstract

I test whether public DESI DR2 baryon acoustic oscillation (BAO) summary statistics require a stationary stochastic component in the late-time dark-energy sector, on top of a smooth background. Fluctuations are modelled as an Ornstein–Uhlenbeck (OU) process in logarithmic scale factor $x=\ln a$, with a quasi-normal-mode (QNM) oscillatory extension. Using a Gaussian likelihood with an additive residual covariance kernel, the best fit drives the stochastic amplitude to zero, and the profile likelihood over $\theta$ ($\Delta\ln\mathcal{L}=-1.92$, full measurement covariance) gives the 95% confidence upper limit

$$\sigma_X < 2.6\times 10^{-2}$$

under the stated phenomenological kernel and public BAO summary statistics (official Gaussian BAO covariance projected to isotropic α in the baseline pipeline). This number is reached in the fully correlated limit of the kernel, which is a shift of $\Omega_m$: it measures consistency with the DESI best-fit background and disappears if $\Omega_m$ and $h r_\mathrm{d}$ are marginalised. A residual that varies between redshifts (θ ≥ 1) is bounded at $\sigma_X < 1.5\times 10^{-2}$ with the background fixed or profiled, and $< 4.6\times 10^{-2}$ with it marginalised. Nested extensions that free $\{w_0,w_a,\sigma_X,\theta\}$ are not preferred over a smooth background ($\Delta\mathrm{AIC}\approx +4$).

Separately, a coherent tachyonic growth model with the correct rank-1 covariance induced by a single growing mode is strongly disfavoured: $\Delta\ln\mathcal{L}\approx -12.7$ relative to $\Lambda$CDM ($\Delta\chi^2\approx +25$). The two published implementations of Everpresent Λ, a zero-mean Λ with order-one fluctuations, do not fit the same BAO data together with the CMB: with r_d computed for each realization, small dark energy around recombination and ω_m within ±2σ of Planck, none fits as well as ΛCDM (0/20000 and 0/500000 realizations). With the acoustic angle θ* computed for each realization in place of the recombination condition, Model 1 stays at 0 (≥ 21σ), while Model 2 admits an expected fraction ≤ 2×10⁻⁵, all with the expansion rate at recombination off by ≥ 11 % against a ≈ 0.6 % Planck sensitivity (N_eff analogy); the full CMB likelihood with a sign-changing dark-energy density was not run.

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

**Important separation:** DESI analyses that combine BAO with CMB+SN often prefer dynamical $w_0$–$w_a$. The reference for $\alpha=1$ in my BAO-only residual analysis is the one of the published $\alpha$ products: the DESI DR2 best-fit flat $\Lambda$CDM ($\Omega_m = 0.2975$, $h r_\mathrm{d} = 101.54$ Mpc; arXiv:2503.14738, Sec. VI, eq. `constraints_lcdm_desi`), which I recover from the data themselves to $3\times 10^{-4}$ in α ([`analisis_cmb/fondo/FONDO.md`](../analisis_cmb/fondo/FONDO.md)); it is not the Planck cosmology. The CPL check of `results/eos_cpl_desi_dr2/` ($w_0\approx -1.001$, $w_a\approx +0.014$) computed its α with a Planck background ($\Omega_m = 0.315$, $H_0 = 67.4$) and with $\Omega_m$ and $h r_\mathrm{d}$ fixed, against data that are relative to the DESI best fit: at $w_0=-1$, $w_a=0$ it reproduces the reference by construction, so only a departure of shape could show up, and its closeness to $\Lambda$CDM is not an independent result. For comparison, DESI DR2 BAO alone with all parameters free gives $w_0=-0.48^{+0.35}_{-0.17}$, $w_a<-1.34$ (68 %), prior-dependent and only 1.7σ from $\Lambda$CDM (arXiv:2503.14738, Sec. VII A, eq. `w0wa_DESI`). The stochastic question is about *residuals*, not a re-fit of the full multi-probe DE equation of state.

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

From the residual analysis (profile / likelihood ratio under the OU kernel as implemented in the public pipeline): $\sigma_X < 2.6\times 10^{-2}$ (95% CL). The profile is shown in Fig. 1 ([`figures/profile_sigma_x.png`](../figures/profile_sigma_x.png)). At the limit the profile sits at the lower edge of the θ search (θ = 10⁻³, fully correlated bins); the θ → 0 supremum is 2.644 × 10⁻² (+0.4 %). Redefining the stationarity clock as cosmic time, conformal time or z instead of ln a moves the limit by less than 0.7 % ([`analisis_cmb/reloj/RELOJ.md`](../analisis_cmb/reloj/RELOJ.md)): it is set by the rank-1, fully correlated limit of the kernel and does not depend on the clock.

**What the limit measures** ([`analisis_cmb/fondo/FONDO.md`](../analisis_cmb/fondo/FONDO.md)). The fully correlated mode (θ → 0) is, to first order, a change of $\Omega_m$: $S_i$ and $\partial\ln\alpha_i/\partial\Omega_m$ have cosine 0.99996. Because the α are referred to the DESI best fit, the data place this mode at $0.02\sigma$, and the limit is the $\ln\det$ (Occam) term of the profile, $\sqrt{e^{3.84}-1}\,\sigma_v = 6.75\,\sigma_v$ with $\sigma_v = (S^{\mathrm T}C^{-1}S)^{-1/2} = 0.0039$, i.e. $2.644\times 10^{-2}$ at θ → 0: consistency with the reference background, not a measurement of a residual. 95 % limits on $\sigma_X$ by treatment of the background ($\Omega_m$, $h r_\mathrm{d}$):

| background | θ free | θ = 0.001 | θ = 0.1 | θ = 1 | θ = 10 |
|:--|:--|:--|:--|:--|:--|
| fixed (DESI best fit) | $2.63\times 10^{-2}$ | $2.63\times 10^{-2}$ | $2.10\times 10^{-2}$ | $1.44\times 10^{-2}$ | $1.29\times 10^{-2}$ |
| profiled | $2.63\times 10^{-2}$ | $2.63\times 10^{-2}$ | $2.13\times 10^{-2}$ | $1.50\times 10^{-2}$ | $1.33\times 10^{-2}$ |
| marginalised (flat prior) | none | none | 0.140 | $4.6\times 10^{-2}$ | $2.2\times 10^{-2}$ |

The θ → 0 part of the limit disappears when the background is marginalised. The bound on a stochastic residual that varies between redshifts (θ ≳ 1) is $1.3$–$1.5\times 10^{-2}$ with the background fixed or profiled, and $2.2$–$4.6\times 10^{-2}$ marginalised. Pantheon+ supernovae with the same kernel bound $\sigma_X$ only at 0.06–0.12, 2.4–4.5 times weaker than DESI, with the same θ → 0 mode degenerate with $\Omega_m$ ([`analisis_cmb/sn/SN.md`](../analisis_cmb/sn/SN.md)). Using the 13 published numbers ($D_M$ and $D_H$ separately, $D_V$ for BGS) instead of the 7 isotropic α tightens every limit by a factor 1.11–1.15 with the background fixed or profiled ($2.3\times 10^{-2}$ at θ → 0, $1.3\times 10^{-2}$ at θ = 1) and by a factor 1.14–1.25 marginalised ($3.8\times 10^{-2}$ at θ = 1), without changing its meaning; the baseline pipeline keeps the 7 α ([`analisis_cmb/fondo/TRECE.md`](../analisis_cmb/fondo/TRECE.md)). Two independent checks follow. (i) The goodness of fit of $\Lambda$CDM to the 13 numbers ($\chi^2 = 10.27$ for 11 d.o.f.) limits the extra scatter a residual may add: requiring $P(\chi^2 \le 10.27) \ge 5\,\%$ gives $\sigma_X < 4.4\times 10^{-2}$ (θ = 1), $2.8\times 10^{-2}$ (θ = 3) and $2.2\times 10^{-2}$ (θ = 10), 11–25 % above the marginalised limits, which it reproduces without using the kernel; at the fixed-background limits the $\chi^2$ does not move (P from 0.49 to 0.45–0.19), confirming that those limits measure consistency with the reference, not scatter ([`analisis_cmb/cuesta/GOF.md`](../analisis_cmb/cuesta/GOF.md)). (ii) A harmonic residual $A\sin(2\pi\ln a/P + \phi)$ with phase and background marginalised is bounded at $A < 0.013$–$0.028$ for $P \le 1$ (median 0.017), the bound growing ∝ P beyond, where only the slope survives the background ([`analisis_cmb/fondo/FISHER_DISENO.md`](../analisis_cmb/fondo/FISHER_DISENO.md)).


This is an upper limit on an effective residual amplitude, not a measurement of a Planck-scale seed.

### 3.3 Nested background freedom

When $\{w_0,w_a\}$ are freed together with residual parameters, the stochastic extension is not preferred ($\Delta\mathrm{AIC}\approx +4$). Marginalising over a flexible smooth background drives $\sigma_X$ toward zero: extra residual variance is absorbed by background freedom when the number of BAO bins is small ($N_{\mathrm{bin}}=7$). Fig. 2 ([`figures/joint_w0wa_sigma.png`](../figures/joint_w0wa_sigma.png)) shows the joint fit of $\{w_0, w_a, \sigma_X\}$. A stochastic residual does not mimic the $w_0$–$w_a$ preference of DESI DR2 ($\Delta\chi^2 = -4.7$, 1.7σ with BAO alone): in 2000 simulated $\Lambda$CDM universes per case with an OU residual at the limits, the best-fit $(w_0, w_a)$ falls in the DESI quadrant ($w_0 > -1$, $w_a < 0$) half of the time, as without residual; the residual only inflates $\Delta\chi^2$ because the fit covariance ignores it, raising the fraction of universes with $\Delta\chi^2 \ge 4.7$ in that quadrant from 4.75 % (noise alone) to 6–18 % at the fixed-background limits (θ = 0.1–3), and to 5–9 % (against 3.75 %) once the $\Lambda$CDM fit is also required to remain acceptable ($\chi^2 \le 19.7$ for 11 d.o.f.) ([`analisis_cmb/cuesta/CUESTA.md`](../analisis_cmb/cuesta/CUESTA.md)).

### 3.4 Claim statement (BAO residuals)

> Under the phenomenological OU/QNM residual kernels and public DESI DR2 BAO summary statistics, there is no statistical evidence for stationary stochastic dark-energy fluctuations. The data prefer smooth evolution; $\sigma_X$ is bounded.

### 3.5 CMB consistency (Planck 2018)

I map the OU residual onto a single expansion history, δ ln ρ_DE = X/Ω_Λ0 with X smoothed on Δ = 0.05 in ln a, convert it to w(a) and run CAMB (PPF) with ω_b, ω_c and Σm_ν fixed to Planck 2018 and H0 solved so that 100θ_* = 1.041198. For 200 realizations at σ_X = 2.633×10⁻²:

| θ | Δχ²_TT (ℓ<30) p50 / p95 | Δχ²_φφ p50 / p95 | fraction Δχ² > 4 | H0 p16 / p50 / p84 [km/s/Mpc] |
|:--|:--|:--|:--|:--|
| 0.001 | −1.5×10⁻⁴ / 2.7×10⁻³ | −1.2×10⁻⁴ / 2.3×10⁻³ | 0 | 67.34 / 67.36 / 67.38 |
| 1 | −2.5×10⁻³ / 0.067 | −2.8×10⁻³ / 0.054 | 0 | 66.92 / 67.36 / 67.88 |

On a grid of ten values from 2.506×10⁻² down to 10⁻⁴ the 95th percentile stays below 4 everywhere and scales linearly with σ_X, so Planck does not bound σ_X below the DESI ceiling. The only CMB imprint is a scatter in H0, σ(H0) = 0.48 km/s/Mpc at θ = 1 and 0.020 at θ = 0.001, i.e. 12 and 280 times smaller than the Planck–SH0ES difference (5.68 km/s/Mpc): the residual does not explain the H0 tension. At the marginalised θ = 1 ceiling of §3.2 ($4.611\times 10^{-2}$, ×1.751), this linear scaling gives a p95 of $0.067\to 0.117$ (TT, θ = 1): still invisible.

A deterministic harmonic residual, X = A sin(2π ln a/P + φ) with A ∈ {0.005, 0.01, 0.025}, P ∈ {0.1, …, 2} and four phases, is allowed by DESI DR2 BAO (13 distances, full covariance, Δχ²_BAO < 4) in 58 of 60 cases, and none of the 60 is distinguishable in the CMB (largest Δχ² among the allowed: 0.070 in TT, 0.63 in lensing).

At high ℓ the residual is equally invisible: for the first 60 of these realizations, the Gaussian distance to $\Lambda$CDM in TT, TE and EE over ℓ ∈ [30, 2500] (cosmic variance plus Planck-like noise, $f_\mathrm{sky}=0.6$) has p95 $\Delta\chi^2 = 3.5\times 10^{-3}$ (θ = 1) and $5.3\times 10^{-6}$ (θ = 0.001); it comes from lensing in the acoustic troughs, not from peak positions (which θ_* fixes), and grows as $\sigma_X^2$: $1.1\times 10^{-2}$ at the marginalised θ = 1 ceiling ([`analisis_cmb/cmb/PICOS.md`](../analisis_cmb/cmb/PICOS.md)).

Fast residuals (θ = 3 and 10 at their limits, $1.3$ and $2.2\times 10^{-2}$) are equally invisible: p95 $\Delta\chi^2 = 0.066$–$0.12$ in TT ℓ < 30 + lensing and $\le 2.5\times 10^{-3}$ in the peaks with linear lensing, less than at θ = 1, because the CMB integrates over the memory 1/θ. The implied equation of state oscillates around −1 with an amplitude that depends on the smoothing Δ (p95 |w + 1| = 0.07 at θ = 3, Δ = 0.05; 0.43 at θ = 10, Δ = 0.0125), because an OU path is not differentiable, while no CMB observable depends on Δ (≤ 4 % between Δ = 0.05 and 0.0125). The Takahashi halofit, which uses w(a) explicitly, produces a spurious growth of $C_L^{\phi\phi}$ for such w (p95 $\Delta\chi^2$ up to 0.65 at θ = 10, Δ = 0.0125, with P(k, z) unchanged); the physical reading is linear lensing, and at θ = 1 the effect is negligible ([`analisis_cmb/cmb/RAPIDO.md`](../analisis_cmb/cmb/RAPIDO.md)).

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

The seed σ_0 ∼ 10⁻⁶¹ uses N ∼ 10¹²² from horizon-area (Bekenstein–Hawking) counting. With the four-volume counting of Everpresent Λ, δΛ tracks H² and the relative amplitude is of order one with no amplifier; the ∼10⁵⁶ gap is specific to the area-counting choice. Zwane et al. (2018) find that Everpresent Λ fits the CMB, BAO and SN data of that time about as well as ΛCDM and eases the low-redshift H0 and Lyα BAO tensions. The order-one alternative to the area-counting seed has been tested directly. In its two published implementations (Zwane et al. 2018, Model 2; Das et al. 2023, Model 1), computing r_d with each realization's own history, requiring small dark energy over z ∈ [1000, 5000] and ω_m within ±2σ of Planck, none of the realizations fits DESI DR2 BAO as well as ΛCDM, nor within Δχ² < 4: 0/20000 for Model 2 (< 1.5×10⁻⁴ at 95 %) and 0/300000 and 0/200000 for Model 1 with α = 0.01 and 0.005 (< 1.0×10⁻⁵, < 1.5×10⁻⁵). Counts are over seeds, with realizations that reach H² ≤ 0 counted as not reproducing the universe; the completed fraction itself is a discretization artefact (for α = 0.01, 27.6 % with Das et al.'s uniform-in-a steps in the production run, 32.3 % in their paper, ≈ 16 % in the continuum limit with steps in ln a, where a √Δ extrapolation and the Broadie–Glasserman–Kou continuity correction agree to 0.5 points), while the survivors' H0 and Ω_Λ(0) and the BAO counts do not depend on the stepping ([`analisis_cmb/everpresent/cmb_real/RELOJ_M1.md`](../analisis_cmb/everpresent/cmb_real/RELOJ_M1.md)). Without the CMB conditions, rare realizations pass, all with ω_m far from Planck or order-one dark energy at recombination. As small fluctuations around ΛCDM, DESI DR2 bounds the Zwane et al. amplitude to α < 1.85×10⁻² (μ fixed). These results hold for those phenomenological implementations, not for the idea of Everpresent Λ in general. Replacing the proxy by the acoustic angle θ* computed from each realization's own history (the full Planck likelihood is not available because CAMB and CLASS cannot follow a dark-energy density that changes sign, as 70 of the 71 BAO survivors do), Model 1 still gives 0 (closest at 21σ). For Model 2 the stretch z ∈ (3, 1000), unconstrained by BAO and sampled here (100 draws per realization), can tune θ*: 5/20000 realizations pass criterion (c) with θ* within 5σ of Planck for some draws (expected fraction 8.5×10⁻⁶; ≤ 2×10⁻⁵ within 10σ), all with |Ω_de(z ≈ 1100)| ≥ 0.26, i.e. H at recombination off from ΛCDM by ≥ 11 % (≥ 6 % in every draw that brings any survivor within 10σ of θ*), whereas Planck's sensitivity to the expansion rate there is ≈ 0.6 % (ΔN_eff = 0.34 at 95 %; a generous background analogy, not a formal bound). The Model 2 exclusion therefore holds with the recombination condition; the full CMB likelihood with a sign-changing dark-energy density was not run. Model 2 maps Ω_de = tanh Ω̂ (Zwane et al. 2018, eq. 2.10), which keeps H² > 0 by construction; |Ω̂| exceeds 1 somewhere in 87.5 % of the 20000 realizations. Removing the map, which is our modification and not a published model, 58.7 % of the realizations reach H² ≤ 0, two of the five θ*-survivors are discarded and the other three fit BAO worse (χ² from 9–28 to 63–1209), while the count with the recombination condition stays at 0: the exclusion does not rest on the map ([`analisis_cmb/everpresent/cmb_real/VALVULA.md`](../analisis_cmb/everpresent/cmb_real/VALVULA.md)).

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

The headline “$\sim 10^{56}$” always means the **Euclid-scale** target. Using the DESI ceiling in a density ratio and labelling it $10^{56}$ is incorrect by a factor $\sim 2\times 10^{3}$. Both gaps remain unbridgeable by soft linear amplifiers ($r\sim 64$–$68$). Against the marginalised θ = 1 ceiling ($4.611\times 10^{-2}$) the DESI gap becomes $G_{\mathrm{DESI}}=3.92\times 10^{59}$ (×1.751); the conclusion is unchanged.

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
| **1. Local $N_{\mathrm{eff}}$** | $\sigma_{0,\mathrm{eff}}=1/\sqrt{N_{\mathrm{eff}}}$, $N_{\mathrm{eff}}\ll N_{\mathrm{BH}}$ | $N_{\mathrm{eff}}\sim 1.44\times 10^{3}$ (DESI), $1.00\times 10^{10}$ (Euclid). Only amplitude card—redefines the seed. $N_{\mathrm{eff}}\propto\sigma_X^{-2}$: at the marginalised θ = 1 ceiling ($4.611\times 10^{-2}$), $\approx 470$ (DESI). |
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

> Slip is the correct *kind* of operator for anisotropic leakage past SDiff. With the BAO-bounded residual amplitude it is marginal, not starved: at ε = 1 it reaches $\lvert\gamma-1\rvert\approx 5\times 10^{-2}$ at $z=0.5$, level with the forecast floor ($\approx 1.0\times$) and $\approx 3.4\times$ below the current measurement. It is not a demonstrated shortcut. The σ_X → slip map is an order-of-magnitude estimate: it holds for δ_m ≈ 1 (quasi-linear scales) and does not extrapolate to δ_m ≪ 1; ε < 1 or z ≳ 1 push it below reach. At the marginalised θ = 1 ceiling ($4.611\times 10^{-2}$, ×1.751) the slip at $z=0.5$ becomes $8.67\times 10^{-2}$, $\approx 1.7\times$ the Sakr et al. (Euclid-like) forecast floor: within reach.

I do not implement a homemade Boltzmann hierarchy. Future work, if amplitudes become reachable, should use community codes (e.g. hi_class / MGCAMB) with a pre-specified target. (The CMB consistency check of §3.5 uses CAMB with the residual mapped to a tabulated w(a) in the PPF module; it is not a Boltzmann treatment of anisotropic stress.)

---

## 7. Limitations

1. **Summary statistics:** baseline residual analysis uses published BAO $\alpha$ with the official Gaussian BAO covariance projected to $\alpha$ (block-diagonal across bins). Cross-bin covariances, if released, may change the numerical factor in $\sigma_X$ while preserving the null preference if residuals remain consistent with noise. Using the full 13-number vector tightens the limits by a factor 1.11–1.25 (§3.2).
2. **Kernel dependence and the fully correlated limit:** bounds are for the OU/QNM residual kernels stated; other residual models need their own covariance and MLE. The headline limit is reached in the fully correlated limit (θ → 0), which is a background parameter, not a residual: a shift of $\Omega_m$ relative to the DESI best fit. That is why it does not depend on the stationarity clock (ln a, t, η, z within ±0.7 %), why DESI DR2 does not resolve θ, and why it vanishes when $\Omega_m$ and $h r_\mathrm{d}$ are marginalised (§3.2); with only seven BAO bins, flexible $w(z)$ can likewise absorb residual variance (Section 3.3). Multi-probe fits are a separate analysis. The CMB check of §3.5 is a consistency test with parameters fixed except H0, not a multi-probe fit.
3. **Cross-correlation (DESI galaxies × SN residuals):** an exploratory angular correlation test exists (`archive/scripts/cross_correlation_DESI.py`) but is preliminary (systematics not fully controlled; DR1 weights). It is not part of the primary claim set of this preprint.
4. **Sorkin seed:** used as motivational UV scale (horizon-area counting) for the amplification discussion, not as a fitted parameter. See §4.2 for Everpresent Λ.
5. **Independent research:** this work has not undergone peer review at the time of writing.
6. **CMB check:** parameters other than H0 fixed to Planck 2018 (no MCMC); smoothing scale Δ introduced for the check; reading (a), δ ln ρ_DE = X/Ω_Λ0 (reading (b) would scale the amplitude by 1.46); the Takahashi halofit misbehaves for a rapidly oscillating w(a) ([`analisis_cmb/cmb/RAPIDO.md`](../analisis_cmb/cmb/RAPIDO.md)); results at θ ≥ 3 use linear lensing.
7. **Everpresent Λ test:** the results of §4.2 hold for the two published phenomenological implementations, not for the idea of Everpresent Λ in general. The CMB enters as θ* computed from each realization's history, together with the proxy $\lvert\Omega_{\rm de}\rvert < 0.25$ over $z\in[1000,5000]$, not as the full likelihood, because CAMB and CLASS do not admit a dark-energy density that changes sign, which 70 of the 71 BAO survivors and the Zwane et al. best fit itself have; the stretch $z\in(3,1000)$ of Model 2 is sampled (N = 100 draws per realization). The comparison of H at recombination with Planck's N_eff sensitivity is a background analogy, not a formal bound; the full likelihood with signed ρ_DE (a modified CAMB) was not run.

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

1. **Null residual:** DESI DR2 BAO summary statistics do not favour stationary OU/QNM stochastic dark energy; $\sigma_X<2.6\times 10^{-2}$ (95% CL) under the stated kernel. That limit is set by the fully correlated mode, i.e. consistency with the DESI best-fit background; for a residual varying between redshifts (θ ≥ 1), $\sigma_X<1.5\times 10^{-2}$ with the background fixed or profiled and $<4.6\times 10^{-2}$ marginalised (§3.2).
2. **Model kill:** coherent tachyonic growth with rank-1 covariance is excluded ($\Delta\chi^2\approx +25$).
3. **Everpresent Λ:** its two published implementations do not fit DESI DR2 BAO together with the CMB conditions (0/20000 and 0/500000 realizations; §4.2). Model 1: 0 also with θ* computed for each realization (≥ 21σ). Model 2: 0 with the recombination condition; θ* alone admits an expected fraction ≤ 2×10⁻⁵, all with H at recombination off by ≥ 11 % against a ≈ 0.6 % Planck sensitivity (N_eff analogy); the full likelihood with signed ρ_DE was not run.
4. **CMB consistency:** at the DESI ceiling the OU residual is invisible to Planck 2018 (p95 Δχ² ≤ 0.067) and does not explain the H0 tension (σ(H0) ≈ 0.48 km/s/Mpc); a percent-level harmonic residual is equally allowed by DESI (58/60) and invisible in the CMB.
5. **Amplification gap:** $\sigma_0\to 10^{-5}$ is $\sim 10^{56}$ (Euclid target); $\sigma_0\to 2.6\times 10^{-2}$ is $\sim 10^{59}$ (DESI ceiling); linear free amplifiers fail both.
6. **Slip honesty:** anisotropic leakage past SDiff is a real structural crack in principle and a marginal channel in practice — at best level with forecast floors ($\approx 1.0\times$ the Sakr et al. floor, $\approx 3.4\times$ below Maus et al., at $z=0.5$). The σ_X → slip map is an order-of-magnitude estimate: it holds for δ_m ≈ 1 (quasi-linear scales) and does not extrapolate to δ_m ≪ 1; ε < 1 or z ≳ 1 push it below reach.
7. **Scientific posture:** the strength of this programme is what it closes—limits, exclusions, and an honest map of open theory prices—not a promised detection.

Future data (Euclid BAO multi-bin products; full DESI covariances) can tighten $\sigma_X$ with the same residual likelihood. The CMB gives no additional leverage at this amplitude (Section 3.5); BAO remains the route. Future theory work must pay the price of Routes 1 or 3 (or an inflation-style seed redefinition)—not invent linear gain on white noise. Two questions follow from §3.2. Standard sirens measure $D_L$ without the $M_B$ calibration of supernovae, but their present precision is far from the effect at the limit: GW170817 gave $D_L = 40^{+8}_{-14}$ Mpc (Abbott et al. 2017, PRL 119, 161101, arXiv:1710.05832), and the GWTC-3 dark sirens $H_0 = 68^{+8}_{-6}$ km/s/Mpc with the galaxy-catalogue method and a fixed BBH population (Abbott et al. 2023, ApJ 949, 76, arXiv:2111.03604), i.e. ≈ 10–35 % per measurement, whereas the fully correlated residual at the limit changes $D_L$ by 0.2 % at $z = 0.1$ to 2.2 % at $z = 2$ ($S(z)\times 2.633\times 10^{-2}$; [`analisis_cmb/sn/SN.md`](../analisis_cmb/sn/SN.md)); like supernovae, they would not by themselves break the degeneracy of that mode with $\Omega_m$. The growth rate is the velocity-side signature: at the DESI limit the residual shifts $f\sigma_8$ by ≈ 3.8 % at $z = 0.5$, 3.5 times below the DESI DR1 full-shape errors (arXiv:2411.12021, Table 9; cosmology in arXiv:2411.12022), with $\Omega_m$ fixed by the CMB (the constant mode is otherwise degenerate with $\Omega_m$, as for the distances); it becomes marginal at the per-bin precision forecast for DESI Y5, 1.0–2.3 % per $\Delta z = 0.1$ bin at $0.75 \le z \le 1.25$ ($k_{\max} = 0.2$–$0.1\,h$ Mpc$^{-1}$; arXiv:1611.00036, Table 2.3), where the shift is 0.7–2.8 σ per bin, while Euclid preparation VII (arXiv:1910.09273) projects the per-bin $f\sigma_8$ onto cosmological parameters without quoting them ([`analisis_cmb/fondo/CRECIMIENTO.md`](../analisis_cmb/fondo/CRECIMIENTO.md)). Survey design: with the background marginalised, the Fisher information on $(\sigma_X, θ)$ is set by the survey volume, not by the number of bins (N = 7 already holds 90 % of the information on θ for θ ≤ 1; finer binning pays only for θ ≳ 3), and a new BAO point at z > 2.5 adds more than any within the DESI range (+114 % at z = 3 and +273 % at z = 4 for θ = 1, against ≤ 15 % anywhere in 0.3 < z < 2.3): the next measurement belongs to z = 2.5–5 surveys (DESI-II, Spec-S5) rather than to the Euclid range ([`analisis_cmb/fondo/FISHER_DISENO.md`](../analisis_cmb/fondo/FISHER_DISENO.md)). Strong-lens time delays (TDCOSMO: H0 to 2.4 % from six lenses with parametric mass profiles, Wong et al. 2020, arXiv:1907.04869; 8 % from seven lenses, or 5 % adding 33 SLACS lenses, once the mass-sheet degeneracy is constrained by kinematics alone, Birrer et al. 2020, arXiv:2007.02941) and redshift drift (Loeb 1998, arXiv:astro-ph/9802122; ELT-ANDES) would measure the expansion with neither ruler nor candle; the effect of the residual at the limit is the same 0.2 % (z = 0.1) to 2.2 % (z = 2) on distances as for the sirens. Turning the headline limit into a bound on stochasticity requires resolving θ, and that needs a residual near the present limits: with DESI DR2 and the background fixed, separating θ = 1 from θ = 10⁻³ needs $\sigma_X \geq 2.2\times 10^{-2}$, above the θ = 1 limit of $1.44\times 10^{-2}$; with the background marginalised a narrow window opens ($2.8$–$4.6\times 10^{-2}$ if the truth is θ = 1); a Euclid-like survey (13 bins at 0.5 %) needs $\sigma_X \geq 1.5\times 10^{-2}$, i.e. it resolves θ only if the residual sits at the present ceiling ([`analisis_cmb/fondo/RESOLVER_THETA.md`](../analisis_cmb/fondo/RESOLVER_THETA.md)).

---

## Version history

- **v1 (July 2026):** σ_X bound quoted as the grid node 2.5×10⁻²; derived numbers (amplification gaps, N_eff, slip, appendix bounds) computed with an earlier working target of 1.5×10⁻⁴; CPL fits with diagonal covariance.
- **v2 (October 2026):** bound = root of the profile, 2.6×10⁻² (2.633×10⁻²); all derived numbers regenerated; CPL fits with the official DESI covariance; CMB consistency and harmonic residual (§3.5) and Everpresent Λ (§4.2) added; slip map corrected (π_T = εσ_X ρ_crit,0).
- **v2.1 (November 2026):** meaning of the σ_X limit stated: the θ → 0 part is the background mode relative to the DESI best fit (Occam term, no bound with the background marginalised), θ ≥ 1 limits by background treatment (§3.2, Pantheon+ cross-check); reference of the α named and the CPL check qualified (§2.1); high-ℓ CMB check (§3.5); derived numbers at the marginalised θ = 1 ceiling (§3.5, §5, §6.2); stationarity clock (§3.2), tanh map of Model 2 and convergence of Model 1 survival (§4.2), θ* per realization (§4.2); outlook on standard sirens and on resolving θ (§9). Before/after: block [F](../analisis_cmb/bloque_F/COMPARACION_F.md).
- **v2.2 (November 2026):** goodness-of-fit and harmonic-residual checks of the limits (§3.2); a stochastic residual does not mimic the $w_0$–$w_a$ preference of DESI DR2 (§3.3); fast residuals θ = 3, 10 in the CMB and the Takahashi halofit caveat (§3.5, limitation 6); survey design by Fisher, strong-lens time delays and redshift drift (§9); evidence of C2, C8 and C9 extended. Before/after: block [G](../analisis_cmb/bloque_G/COMPARACION_G.md).
- **v2.2.1 (November 2026):** Zenodo deposit (`.zenodo.json`, `CITATION.cff`); final proofreading of the manuscript files; no change to the analysis or the numbers.
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
9. Morales Souhail, J., this repository and companion notes (2026).
10. Abbott, B. P., et al. (LIGO Scientific and Virgo Collaborations), "GW170817: Observation of Gravitational Waves from a Binary Neutron Star Inspiral," Phys. Rev. Lett. 119, 161101 (2017), [arXiv:1710.05832](https://arxiv.org/abs/1710.05832).
11. Abbott, R., et al. (LIGO, Virgo and KAGRA Collaborations), "Constraints on the Cosmic Expansion History from GWTC-3," Astrophys. J. 949, 76 (2023), [arXiv:2111.03604](https://arxiv.org/abs/2111.03604).
12. Brout, D., et al., "The Pantheon+ Analysis: Cosmological Constraints," Astrophys. J. 938, 110 (2022), [arXiv:2202.04077](https://arxiv.org/abs/2202.04077).
13. Wong, K. C., et al. (H0LiCOW), "H0LiCOW XIII. A 2.4% measurement of $H_0$ from lensed quasars: 5.3σ tension between early- and late-Universe probes," Mon. Not. R. Astron. Soc. 498, 1420 (2020), [arXiv:1907.04869](https://arxiv.org/abs/1907.04869).
14. Birrer, S., et al. (TDCOSMO), "TDCOSMO IV: Hierarchical time-delay cosmography — joint inference of the Hubble constant and galaxy density profiles," Astron. Astrophys. 643, A165 (2020), [arXiv:2007.02941](https://arxiv.org/abs/2007.02941).
15. Loeb, A., "Direct Measurement of Cosmological Parameters from the Cosmic Deceleration of Extragalactic Objects," Astrophys. J. 499, L111 (1998), [arXiv:astro-ph/9802122](https://arxiv.org/abs/astro-ph/9802122).

