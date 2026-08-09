# Euclid BAO analysis for vacuum smoothness and mean-reversion

Jesús Morales Souhail · github.com/jesus-morales-souhail · July 2026

Related notes: `fundamental-vs-emergent-vacuum-relaxation.md`, `sdiff-fundamental-vs-emergent.md`, `stochastic-dark-energy-desi-dr2.md`.

---

## Abstract

I set out a Euclid-scale BAO analysis for residual stochastic fluctuations of dark energy, extending the DESI DR2 working limit $\sigma_X < 1.5\times 10^{-4}$ (95% CL). The fit uses the parameter vector $\{w_0,w_a,\theta,\sigma_X\}$, with optional $\omega_R$ for a nested QNM kernel. Path residuals of the form $\sigma_{\rm res}=A_0 e^{-\theta\Delta x}$ define operational regions (F, E0–E3). For $\theta\sim\mathcal{O}(1)$, amplitude-based damping can be visible when $A_0$ exceeds the measurement noise, while lag–shape tests with a Pearson estimator alone are weak unless $N$ is large or the full OU kernel is used in the likelihood. A BAO-only null does not separate fundamental geometric silence (F) from an emergent but unobservable seed (E0); a residual in $10^{-5}\lesssim\sigma_X\lesssim 10^{-4}$ would favour E1-type effective physics.

---

## 1. What is constrained

Three layers enter the interpretation:

| Layer | Quantity | Role |
|-------|----------|------|
| Geometry | Local $T_{\mu\nu}=V(x)g_{\mu\nu}$ projected by SDiff | Exact silence if fundamental (region F) |
| Microscopic seed | Sorkin / Bekenstein–Hawking $\sigma_0\sim 10^{-61}$ | Always $\ll$ BAO noise |
| Effective BAO residual | $\sigma_{\rm res}=A_0 e^{-\theta\Delta x}$ after the expansion path | What DESI and Euclid constrain |

Euclid therefore constrains the effective residual amplitude and, when the data allow, the shape of the OU kernel. It does not measure $\sigma_0$ directly. A deep null remains compatible with both F and E0.

---

## 2. Parameter vector

Baseline Euclid BAO-focused model:


$$
\Theta = \{w_0, w_a, \theta, \sigma_X\}.
$$


| Parameter | Role |
|-----------|------|
| $w_0,w_a$ | Smooth CPL background (breaks DE–noise degeneracy) |
| $\theta$ | OU mean-reversion in $x=\ln a$; $\Gamma_{\rm phys}=\theta H$ |
| $\sigma_X$ | Amplitude scale of the additive OU covariance |

Equivalent reparameterizations (Jacobian handled as usual):

- $\{\theta, A_0\}$ with $\sigma_{\rm res}(z)\sim A_0 e^{-\theta\Delta x(z)}$ for a single frozen kick;
- $\{\Gamma_0, A_0\}$ with $\Gamma_0=\theta H_0$ at $z=0$.

Additive kernel as in this repository:


$$
(C_{\rm OU})_{ij}=S(z_i)S(z_j) \sigma_X^{2} e^{-\theta\lvert x_i-x_j\rvert},
$$


plus the survey covariance $C_{\rm std}$.

Optional nested extension:


$$
\Theta_+ = \{w_0, w_a, \theta, \sigma_X, \omega_R\}
$$


for a damped oscillatory (QNM) kernel. If $\omega_R\to 0$, the model reduces to pure OU.

---

## 3. Priors

Default priors for a reproducible pipeline:

| Parameter | Prior | Rationale |
|-----------|--------|-----------|
| $w_0$ | Uniform $[-1.5, -0.5]$, or Gaussian about an external CPL with width $\sim 0.1$ | Does not force $w=-1$ |
| $w_a$ | Uniform $[-2, 1]$, or Gaussian of width $\sim 0.3$–$0.5$ | Covers the DESI-preferred quadrant |
| $\theta$ | Log-uniform on $[10^{-3}, 10]$ | From numerical floor to strong damping |
| $\sigma_X$ | Log-uniform on $[10^{-6}, 10^{-2}]$, or half-Gaussian at 0 with scale $10^{-4}$ | Spans Euclid targets and the DESI limit |
| $\omega_R$ (optional) | Uniform $[0, 20]$ with $\theta\ge 10^{-3}$ | Avoids undamped modes |

Physical constraints: $\theta > 0$, $\sigma_X \ge 0$; for QNM, $\theta \ge \theta_{\min}\approx 10^{-3}$. Flatness $\Omega_m+\Omega_{\rm DE}=1$ is assumed in the BAO kernel unless $\Omega_m$ is sampled jointly.

A prior peaked at $\sigma_X\sim 10^{-61}$ is inappropriate: that scale is a UV seed, not the effective BAO parameter.

---

## 4. Likelihood and model comparison

Gaussian residual likelihood (as in `scripts/ou_bao_stochastic_test.py`):


$$
-2\ln\mathcal{L} = \mathbf{r}^{\top} C^{-1}\mathbf{r}+\ln\det C+\mathrm{const}, \qquad C=C_{\rm std}+C_{\rm OU}(\theta,\sigma_X).
$$


For Euclid, I replace the DESI 7-bin vector by the Euclid BAO vector ($\gtrsim 20$ bins when available) and recompute $S(z)$.

| Question | Statistic |
|----------|-----------|
| Need for a noise component | Nested CPL vs CPL+$\{\theta,\sigma_X\}$: $\Delta\mathrm{AIC}$, $\Delta\mathrm{BIC}$, or Bayes factor |
| Amplitude | 95% upper limit or interval on $\sigma_X$ with free $\{w_0,w_a,\theta\}$ |
| Mean-reversion | Marginal $p(\theta\mid\mathrm{data})$ from the full kernel |
| QNM | Nested OU vs QNM on $\omega_R$; require $\omega_R$ away from 0 and $\theta$ not at the floor |

I reserve detection language for cases where the Bayes factor (or $\Delta\mathrm{AIC}$) and the $\sigma_X$ posterior both support a non-zero amplitude.

---

## 5. Operational regions

With residual scale of order $\sigma_X$:

| Region | Criterion | Interpretation |
|--------|-----------|----------------|
| **F** | $\sigma_X$ consistent with 0; no shape signal | Compatible with geometric silence |
| **E0** | Observationally as F; theory assumes a Sorkin-only seed | Emergent but invisible |
| **E1** | $10^{-5}\lesssim \sigma_X \lesssim 1.5\times 10^{-4}$, $\theta$ small | Emergent amplitude above the Poisson seed |
| **E2** | $\theta$ constrained away from 0 with reduced amplitude relative to no-damping | Detectable mean-reversion |
| **E3** | $\sigma_X \gtrsim 1.5\times 10^{-4}$ with free background | Tension with the DESI null |

Euclid alone cannot split F from E0. A null should be reported as **F $\cup$ E0**.

---

## 6. Detectability of $\theta\sim\mathcal{O}(1)$

### 6.1 Amplitude-based damping

If a no-damping model would over-predict the residual and the absolute uncertainty on $\sigma_X$ is $s$,


$$
A_0\bigl(1-e^{-\theta\Delta x}\bigr)\gtrsim s \qquad\Rightarrow\qquad \theta \gtrsim -\frac{1}{\Delta x}\ln\Bigl(1-\frac{s}{A_0}\Bigr) \quad(s<A_0).
$$


Examples for a DESI-wide path $\Delta x=0.94$:

| $s$ | $A_0$ | Min $\theta$ (order of magnitude) |
|------:|--------:|------------------------------------:|
| $10^{-5}$ | $3\times 10^{-5}$ | $\sim 0.4$ |
| $10^{-5}$ | $10^{-4}$ | $\sim 0.1$ |
| $3\times 10^{-5}$ | $10^{-4}$ | $\sim 0.4$ |
| $5\times 10^{-5}$ | $1.5\times 10^{-4}$ | $\sim 0.4$ |

On a narrower Euclid DR1 path ($\Delta x\sim 0.39$), the same amplitude criterion requires larger $\theta$ by a factor $\sim 0.94/0.39\sim 2.4$. Thus $\theta\sim\mathcal{O}(1)$ is relevant when $A_0$ is a few times the noise $s$.

### 6.2 Lag correlations

For lag-1 Pearson correlations with error $\sim 1/\sqrt{N-3}$:

| Survey (schematic) | $N$ | $\sigma_\rho$ | Min $\theta$ (order) |
|--------------------|------:|----------------:|-----------------------:|
| DESI-like | 7 | $\sim 0.5$ | $\sim 4$ (unusable) |
| Euclid-like | 20 | $\sim 0.24$ | $\gtrsim 5$ |
| Euclid-like | 40 | $\sim 0.16$ | $\gtrsim 7$ |

With $N=7$, no lag is individually significant in the present DESI analysis. Lag tables remain useful diagnostics; E2 claims should rest on the joint posterior of $\theta$ from $C_{\rm OU}$.

---

## 7. Outcome scenarios

**Marginal $\sigma_X$.** Posterior near $10^{-5}$–few $\times 10^{-5}$, $\theta$ unconstrained. Consistent with weak E1; systematics and free $\{w_0,w_a\}$ required before a firm claim.

**Strong null.** Upper limit $\sigma_X\ll 10^{-5}$, $\theta$ prior-like. Compatible with F $\cup$ E0.

**Detection with shape.** $\sigma_X$ in the E1 window and $p(\theta\mid\mathrm{data})$ excludes $\theta\to 0$. Strongest BAO-only support for dynamical relaxation (E1/E2).

**Apparent detection that vanishes when $\{w_0,w_a\}$ are free.** Background degeneracy; not an E1 detection. This mode already appears in the DESI BAO-only analysis of this repository.

---

## 8. Multi-probe extension

If Euclid BAO leaves F/E0 degenerate, a natural joint analysis includes:

| Probe | Contribution |
|-------|----------------|
| Pantheon+ / future SN | Background $\{w_0,w_a\}$ |
| Growth $f\sigma_8$ | Smooth DE versus clustering |
| ISW / CMB lensing | DE perturbations if present |
| Imaging systematics weights | Required for angular cross-correlations |

I report $\sigma_X$ and $\theta$ conditional on free smooth DE.

---

## 9. Analysis sequence

1. Construct the Euclid BAO data vector, $C_{\rm std}$, and $S(z)$.
2. Sample $\{w_0,w_a,\theta,\sigma_X\}$ with the priors of §3.
3. Report the marginal interval or upper limit on $\sigma_X$, the posterior on $\theta$, and nested model comparison versus CPL-only.
4. Use lag correlations only as diagnostics.
5. Assign F $\cup$ E0 / E1 / E2 / E3 using §§5–7.
6. Compare fixed versus free $\{w_0,w_a\}$; optionally nest QNM.

---

## 10. Numerical products

| File | Content |
|------|---------|
| `results/euclid_protocol/euclid_forecast_grid.csv` | Residual forecasts over $(\theta,A_0)$ |
| `figures/euclid_A0_theta_plane.png` | Contours in $(A_0,\theta)$ |
| `figures/euclid_theta_shape_detectability.png` | Lag-shape $\theta$ thresholds vs $N$ |
| `scripts/desqueezing/euclid_protocol_forecasts.py` | Regenerates grids and figures |

```bash
python scripts/desqueezing/euclid_protocol_forecasts.py
```

---

## 11. Summary

- Fit $\{w_0,w_a,\theta,\sigma_X\}$ with weakly informative log priors on $\theta$ and $\sigma_X$.
- E1 is the BAO window in which Euclid can favour residual amplitude above a pure Poisson seed.
- E2 requires damping large enough to exceed amplitude noise when $A_0>s$; Pearson lags alone are usually insufficient.
- F versus E0 remains a theoretical distinction under a deep null.

---

## 12. DESI real BAO residual test (repository)

Reference implementation: `scripts/desi_dr2_real_bao_test.py` (emcee). Euclid-like vector of 24 bins ($z\in[0.9,1.8]$), OU kernel as above; samples $\{\theta,\sigma_X\}$ at fixed background or the full $\{w_0,w_a,\theta,\sigma_X\}$.

With percent-level BAO uncertainties ($\sigma_\alpha\sim 0.5\%$–$1.2\%$) and $S(z)\sim\mathcal{O}(1)$, the OU term competes with $C_{\rm std}$ only for


$$
\sigma_X \sim \mathrm{few}\times 10^{-3}\ \text{to}\ 10^{-2}.
$$


The scientific window $10^{-5}$–$1.5\times 10^{-4}$ therefore sits below the single-survey SNR of this simplified forecast. The DESI working limit remains a phenomenological upper bound from multi-bin MLE behaviour, not a high-SNR detection scale in a simplified forecast vector with $\sim 1\%$ errors.

| Scenario (truth) | Fixed $ (w_0,w_a)$ | Free $ (w_0,w_a)$ |
|------------------|---------------------|---------------------|
| Null $\sigma_X=10^{-6}$ | Prior-like posterior | Same; background widens further |
| E1 $\sigma_X=0.012,\ \theta=0.1$ | Partial amplitude recovery; $\theta$ weak | Signal largely absorbed by background freedom |
| E2 $\sigma_X=0.015,\ \theta=1.5$ | Amplitude recovered in order of magnitude; $\theta$ poorly pinned | Recovery degrades with free $w_0,w_a$ |

I free $\{w_0,w_a\}$ before interpreting E1/E2. Products: `results/euclid_mcmc/`, `figures/euclid_mcmc_*_{fixed,free}BG.png`.

```bash
pip install emcee
python scripts/desi_dr2_real_bao_test.py
```

---

## 13. Joint BAO on DESI (real)

Implementation: `scripts/euclid_joint_bao_sne_mcmc.py`.

- **BAO:** 20 bins, $z\in[0.9,1.8]$, forecast-style $\sigma_\alpha\sim 0.3\%$–$0.7\%$.
- **SN:** 25 compressed distance-modulus bins (relative fit), shared $\{w_0,w_a\}$.
- **Noise sector:** OU kernel on BAO only.

| Scenario | Probe | $\sigma_X$ (illustrative) | $\theta$ | Comment |
|----------|-------|----------------------------:|-----------:|---------|
| Null | BAO | $\sim 2\times 10^{-4}$ (upper, prior-like) | unconstrained | No false high-$\sigma_X$ peak |
| Null | BAO+SN | similar $\sigma_X$; tighter $w_0$ | unconstrained | SN helps the background |
| E1 ($\sigma_X=0.008,\theta=0.1$) | BAO / joint | weak recovery | weak | Edge of optimistic SNR |
| E2 ($\sigma_X=0.012,\theta=1.5$) | BAO | $\sim 0.016\ [0.000,0.029]$ | poorly pinned | Amplitude recovered at right order |
| E2 | BAO+SN | $\sim 0.007\ [0.000,0.027]$ | median higher | SN reduces DE–noise trade-off |

Optimistic Euclid errors move the practical recovery threshold toward $\sigma_X\sim 1\%$, but $10^{-5}$–$10^{-4}$ remains demanding. SN data primarily tighten $\{w_0,w_a\}$. Products: `results/euclid_joint_mcmc/`, `figures/joint_mcmc_*.png`.

```bash
python scripts/joint_w0wa_sigma_desi.py
```

---

## References

[1] Morales Souhail, J., DESI OU/QNM analysis, this repository.
[2] Morales Souhail, J., Fundamental vs emergent vacuum relaxation, this repository.
[3] DESI Collaboration, arXiv:2503.14738 (2025).
[4] Euclid Collaboration, DR1 BAO forecasts (as available at analysis time).

