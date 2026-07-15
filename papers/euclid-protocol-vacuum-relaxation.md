# A Minimal Euclid Protocol for Fundamental versus Emergent Vacuum Relaxation

**Author:** Jesús Morales Souhail  
**Date:** July 2026  
**ORCID:** [0009-0000-7637-1818](https://orcid.org/0009-0000-7637-1818)  
**Repository:** https://github.com/jesus-morales-souhail/stochastic-dark-energy-ou  
**Status:** Analysis protocol note — not peer reviewed  
**Companions:** `fundamental-vs-emergent-vacuum-relaxation.md`, `sdiff-fundamental-vs-emergent.md`, `stochastic-dark-energy-desi-dr2.md`

---

## Abstract

We specify a minimal but sharp analysis protocol for Euclid-scale BAO (and near-term joint) tests of vacuum smoothness. Building on the DESI null \(\sigma_X < 1.5\times 10^{-4}\) (95% CL) and the path-integrated residual \(\sigma_{\rm res}=A_0 e^{-\theta\Delta x}\), we define the parameter vector, priors, decision regions (F / E0 / E1 / E2 / E3), and outcome scenarios. We refine region **E2** (\(\theta\sim\mathcal{O}(1)\)): amplitude-based damping is potentially visible when \(A_0\) exceeds the measurement noise, whereas lag-shape detection with a simple Pearson estimator is weak unless \(N\) is large or a full OU-kernel likelihood is used. Euclid alone does not separate fundamental SDiff (**F**) from emergent Sorkin-only silence (**E0**); a detection in \(10^{-5}\lesssim\sigma_X\lesssim 10^{-4}\) selects **E1**-type physics.

---

## 1. Logic of the test (what is actually being constrained)

Three layers must not be confused:

| Layer | Object | Status |
|-------|--------|--------|
| Geometric | Local \(T_{\mu\nu}=V(x)g_{\mu\nu}\) projected out by SDiff | **F** if exact |
| Micro seed | Sorkin / Bekenstein–Hawking \(\sigma_0\sim 10^{-61}\) | Always \(\ll\) BAO noise |
| Effective BAO | Residual amplitude after expansion path | What DESI/Euclid measure |

**Honest target:** Euclid tests the **effective** residual \(\sigma_{\rm res}=A_0 e^{-\theta\Delta x}\) (and, if the data allow, the shape of the OU kernel). It does **not** directly measure \(\sigma_0\), and a null cannot distinguish **F** from **E0**.

---

## 2. Parameter vector

### 2.1 Recommended baseline (Euclid BAO-focused)

$$
\Theta = \{w_0,\, w_a,\, \theta,\, \sigma_X\}.
$$

| Parameter | Role |
|-----------|------|
| \(w_0,w_a\) | Smooth CPL background (breaks DE / noise degeneracy) |
| \(\theta\) | OU mean-reversion in \(x=\ln a\); \(\Gamma_{\rm phys}=\theta H\) |
| \(\sigma_X\) | Effective stationary / residual amplitude scale of the OU kernel (same role as \(\sqrt{\mathrm{Var}(X)}\) or the \(A_0\)-like overall scale of \(C_{\rm OU}\)) |

Equivalent reparameterizations (do not change physics if Jacobian is handled):

- \(\{\theta,\, A_0\}\) with \(\sigma_{\rm res}(z)\sim A_0 e^{-\theta\Delta x(z)}\) for a single frozen kick;  
- \(\{\Gamma_0,\, A_0\}\) with \(\Gamma_0=\theta H_0\) at \(z=0\).

**Baseline recommendation:** fit \(\{w_0,w_a,\theta,\sigma_X\}\) with the additive kernel already used in this repository,

$$
(C_{\rm OU})_{ij}=S(z_i)S(z_j)\,\sigma_X^2\,e^{-\theta\lvert x_i-x_j\rvert},
$$

plus the survey covariance \(C_{\rm std}\).

Optional extension (only if data support it):

$$
\Theta_+ = \{w_0,\, w_a,\, \theta,\, \sigma_X,\, \omega_R\}
$$

for the damped oscillatory (QNM) kernel. Prefer nested comparison: if \(\omega_R\to 0\), recover OU.

---

## 3. Priors (deliberately conservative)

Priors must not smuggle a detection. Suggested defaults for a public, reproducible pipeline:

| Parameter | Prior | Rationale |
|-----------|--------|-----------|
| \(w_0\) | Uniform \([-1.5,\,-0.5]\) or Gaussian centered on external CPL with width \(\sim 0.1\) | Wide enough not to force \(w=-1\) |
| \(w_a\) | Uniform \([-2,\,1]\) or Gaussian width \(\sim 0.3\)–\(0.5\) | Covers DESI-preferred quadrant without hard walls at 0 |
| \(\theta\) | Log-uniform on \([10^{-3},\,10]\) | From MLE numerical floor to strongly damped e-fold scales |
| \(\sigma_X\) | Log-uniform on \([10^{-6},\,10^{-2}]\) or half-Gaussian at 0 with scale \(10^{-4}\) | Straddles Euclid target and DESI limit; does not force a detection |
| \(\omega_R\) (optional) | Uniform \([0,\,20]\) with \(\theta\ge 10^{-3}\) | Avoids undamped oscillations (unphysical QNM) |

**Hard constraints (physics):**

- \(\theta > 0\), \(\sigma_X \ge 0\);  
- if QNM: \(\theta \ge \theta_{\min}\approx 10^{-3}\) (same floor as this repository’s optimizer);  
- background flatness \(\Omega_m+\Omega_{\rm DE}=1\) as in the BAO kernel construction unless jointly sampling \(\Omega_m\).

**Do not** place a prior on \(\sigma_X\) peaked at \(10^{-61}\): that is a UV seed, not the effective BAO parameter.

---

## 4. Likelihood and statistics

### 4.1 Likelihood

Gaussian BAO residual model (as implemented in `scripts/ou_bao_stochastic_test.py`):

$$
-2\ln\mathcal{L}
=
\mathbf{r}^\top C^{-1}\mathbf{r}+\ln\det C+\mathrm{const},
\qquad
C=C_{\rm std}+C_{\rm OU}(\theta,\sigma_X).
$$

For Euclid, replace the 7-bin DESI vector by the Euclid BAO vector (\(\gtrsim 20\) bins when available) and recompute \(S(z)\).

### 4.2 Decision statistics (pre-registered)

| Question | Statistic | Prefer |
|----------|-----------|--------|
| Is a noise component needed? | Nested \(\Lambda\)CDM/CPL vs CPL+\(\{\theta,\sigma_X\}\): \(\Delta\mathrm{AIC}\), \(\Delta\mathrm{BIC}\), or Bayes factor \(B_{10}\) | Detection language only if \(B_{10}\gtrsim 3\)–\(5\) (or \(\Delta\mathrm{AIC}\lesssim -6\)) **and** posterior on \(\sigma_X\) away from 0 |
| Amplitude only | 95% upper limit on \(\sigma_X\) with free \(\{w_0,w_a,\theta\}\) | Profile or marginal posterior |
| Shape / mean-reversion | Marginal posterior \(p(\theta\mid\mathrm{data})\); compare to prior | Prefer full-kernel likelihood over single-lag Pearson |
| QNM | Nested OU vs QNM on \(\omega_R\) | Require \(\omega_R\) posterior away from 0 **and** \(\theta\) not at floor |

---

## 5. Region map (operational)

Using the DESI-path residual scale \(\sigma_{\rm res}\sim\sigma_X\) (same order as the OU amplitude entering \(C_{\rm OU}\)):

| Region | Operational criterion (posterior) | Interpretation |
|--------|-----------------------------------|----------------|
| **F** | \(\sigma_X\) consistent with 0 at Euclid precision; no shape signal | Compatible with geometric silence |
| **E0** | Same observationally as F; theory assumes Sorkin-only seed | Emergent but invisible |
| **E1** | \(10^{-5}\lesssim \sigma_X \lesssim 1.5\times 10^{-4}\) with \(\theta\) small (posterior mass at \(\theta\lesssim 0.3\)) | Emergent + effective amplification |
| **E2** | \(\theta\) constrained away from 0 **and** amplitude reduced relative to a no-damping model | Detectable mean-reversion / path damping |
| **E3** | \(\sigma_X \gtrsim 1.5\times 10^{-4}\) at high posterior probability with free background | Tension with DESI null |

**Euclid alone cannot split F from E0.** Label a null as “F \(\cup\) E0”.

---

## 6. Region E2 refined: when is \(\theta\sim\mathcal{O}(1)\) detectable?

Two different notions of “detectable damping” must be separated.

### 6.1 Amplitude-based damping (primary for E2)

Suppose the effective kick \(A_0\) is large enough that a no-damping model would over-predict the residual. With absolute uncertainty \(s\) on \(\sigma_X\),

$$
A_0\bigl(1-e^{-\theta\Delta x}\bigr)\gtrsim s
\qquad\Rightarrow\qquad
\theta \gtrsim -\frac{1}{\Delta x}\ln\Bigl(1-\frac{s}{A_0}\Bigr)
\quad(s<A_0).
$$

**Numerical examples** (\(\Delta x=0.94\) DESI-wide path):

| \(s\) (abs. unc.) | \(A_0\) | Min \(\theta\) (order of magnitude) |
|------------------:|--------:|------------------------------------:|
| \(10^{-5}\) | \(3\times 10^{-5}\) | \(\sim 0.4\) |
| \(10^{-5}\) | \(10^{-4}\) | \(\sim 0.1\) |
| \(3\times 10^{-5}\) | \(10^{-4}\) | \(\sim 0.4\) |
| \(5\times 10^{-5}\) | \(1.5\times 10^{-4}\) | \(\sim 0.4\) |

On a **narrower** Euclid DR1 path (\(\Delta x\sim 0.39\), \(z\sim 0.9\)–\(1.8\)), the same amplitude criterion requires **larger** \(\theta\) by a factor \(\sim 0.94/0.39\sim 2.4\).

**Reading for E2:** \(\theta\sim\mathcal{O}(1)\) is in the right ballpark for amplitude-based damping **if** \(A_0\) is a few times the noise \(s\). If \(A_0\sim s\), damping cannot be seen from amplitude alone.

### 6.2 Lag-shape detection (secondary; often weak)

For lag-1 Pearson correlations with error \(\sim 1/\sqrt{N-3}\), requiring \(1-e^{-\theta\Delta x_{\rm lag}}\gtrsim k\,\sigma_\rho\) gives:

| Survey (schematic) | \(N\) | \(\sigma_\rho\) | Min \(\theta\) (order, \(k=1\)) |
|--------------------|------:|----------------:|--------------------------------:|
| DESI-like | 7 | \(\sim 0.5\) | \(\sim 4\) (effectively **unusable**) |
| Euclid-like | 20 | \(\sim 0.24\) | \(\gtrsim 5\) (wide path) — still large |
| Euclid-like | 40 | \(\sim 0.16\) | \(\gtrsim 7\) |

This matches the existing DESI lag analysis: with \(N=7\), no lag is individually significant. **Falsifier #3 (positive OU-like lags) is therefore a high bar** unless:

- \(N\) is large **and** the full OU kernel is fit in the likelihood (not single-lag Pearson alone), or  
- external data (SNe residuals, multi-tracer BAO) tighten the covariance shape.

**Protocol rule:** report lag diagnostics, but base E2 claims on the **joint posterior of \(\theta\) from \(C_{\rm OU}\)**, not on a single lag \(p\)-value.

---

## 7. Pre-registered outcome scenarios

### 7.1 Marginal detection in \(\sigma_X\)

- Posterior: \(\sigma_X\) peaks near \(10^{-5}\)–few\(\times 10^{-5}\), 95% interval excludes 0 weakly.  
- \(\theta\) unconstrained or piled at the lower edge.  
- **Decision:** favor **E1** over F/E0 at low-to-moderate significance; **do not** claim fundamental SDiff is ruled out at high confidence. Demand systematics control and joint \(\{w_0,w_a\}\).

### 7.2 Strong null

- 95% (or 99%) upper limit \(\sigma_X\ll 10^{-5}\).  
- \(\theta\) returns to prior.  
- **Decision:** data compatible with **F \(\cup\) E0**. Emergent scenarios without amplification become less attractive but not eliminated if one accepts \(A_0\sim\sigma_0\).

### 7.3 Detection + shape (\(\theta\) away from 0)

- \(\sigma_X\) in the E1 window **and** \(p(\theta\mid\mathrm{data})\) excludes \(\theta\to 0\) (e.g. 95% lower bound \(\theta\gtrsim 0.1\)–\(0.3\) depending on \(A_0\) and \(s\)).  
- Lag diagnostics consistent with \(e^{-\theta\Delta x}\) (supporting, not decisive alone).  
- **Decision:** strongest empirical support for **dynamical relaxation (E1/E2)**. This is the most damaging outcome for pure algebraic silence (**F**) among BAO-only tests.

### 7.4 Apparent detection that vanishes when \(\{w_0,w_a\}\) are free

- **Decision:** background degeneracy; do **not** claim E1. This failure mode is already emphasized in the DESI OU analysis of this repository.

---

## 8. Minimal joint upgrade (when BAO alone is not enough)

If Euclid BAO still leaves F/E0 degenerate, the next logical layer is a **pre-declared joint**:

| Probe | What it adds |
|-------|----------------|
| Pantheon+ / future SN | Background \(\{w_0,w_a\}\) and residual cross-checks |
| Redshift-space growth \(f\sigma_8\) | Smooth DE vs clustering; consistency of background |
| ISW / CMB lensing cross | Optional; sensitive to DE perturbations if present |
| Imaging systematics weights | Mandatory before claiming angular cross-correlations |

The vacuum-relaxation protocol itself stays the same: always quote \(\sigma_X\) and \(\theta\) **conditional on free smooth DE**.

---

## 9. One-page checklist (copy into an analysis notebook)

1. Build Euclid BAO data vector + \(C_{\rm std}\) + \(S(z)\).  
2. Sample \(\{w_0,w_a,\theta,\sigma_X\}\) with priors in §3.  
3. Report:  
   - marginal 95% upper limit / interval on \(\sigma_X\);  
   - marginal posterior on \(\theta\);  
   - \(\Delta\mathrm{AIC}/\Delta\mathrm{BIC}\) or Bayes factor vs CPL-only;  
   - lag table as diagnostic only.  
4. Assign region F∪E0 / E1 / E2 / E3 using §5–§7.  
5. Stress test: fix vs free \(\{w_0,w_a\}\); optional QNM nest.  
6. No Sorkin prior on \(\sigma_X\); discuss amplification only in interpretation.

---

## 10. Numerical products

| File | Content |
|------|---------|
| `results/euclid_protocol/euclid_forecast_grid.csv` | Residual forecasts over \((\theta,A_0)\) |
| `figures/euclid_A0_theta_plane.png` | Detection / DESI-tension contours in \((A_0,\theta)\) |
| `figures/euclid_theta_shape_detectability.png` | Rough lag-shape \(\theta\) thresholds vs \(N\) |
| `scripts/desqueezing/euclid_protocol_forecasts.py` | Regenerates grids and figures |

```bash
python scripts/desqueezing/euclid_protocol_forecasts.py
```

---

## 11. Bottom line

- Fit \(\{w_0,w_a,\theta,\sigma_X\}\) with log-ish priors on \(\theta,\sigma_X\).  
- **E1** is the only BAO-scale region where Euclid can positively favor emergence-with-amplification.  
- **E2** needs \(\theta\) large enough that damping exceeds amplitude noise when \(A_0>s\); lag Pearson alone is usually insufficient.  
- **F vs E0** remains a theoretical split under a deep null; only new physics (amplification) or new observables (shape + multi-probe) break it.

That is the sharpest Euclid-facing protocol the present framework supports without overselling.

---



---

## 12. Mock MCMC validation (this repository)

A reference implementation lives in `scripts/euclid_mock_mcmc.py` (emcee). It builds a 24-bin Euclid-like BAO vector (\(z\in[0.9,1.8]\)), uses the repository OU kernel, and samples either \(\{\theta,\sigma_X\}\) at fixed background or the full \(\{w_0,w_a,\theta,\sigma_X\}\).

### 12.1 Sensitivity reality check

With percent-level BAO uncertainties (\(\sigma_\alpha\sim 0.5\%\text{--}1.2\%\)) and \(S(z)\sim\mathcal{O}(1)\), the OU contribution \(S_i S_j\sigma_X^2\) only competes with \(C_{\rm std}\) when

$$
\sigma_X \sim \mathrm{few}\times 10^{-3}\ \text{to}\ 10^{-2}.
$$

The scientific window \(10^{-5}\text{--}1.5\times 10^{-4}\) therefore sits **below the single-mock SNR** of this simplified forecast. The published DESI working limit remains a phenomenological upper bound from multi-bin MLE behaviour, not a claim that \(\sigma_X\sim 10^{-4}\) is a high-SNR detection scale in a 24-bin mock with \(\sim 1\%\) errors.

### 12.2 What the mocks show (illustrative run)

| Scenario (truth) | Fixed \((w_0,w_a)\) | Free \((w_0,w_a)\) |
|------------------|---------------------|---------------------|
| Null \(\sigma_X=10^{-6}\) | Posterior prior-like; no false high-\(\sigma_X\) claim | Same; background widens further |
| E1 \(\sigma_X=0.012,\ \theta=0.1\) | Amplitude only partially recovered; \(\theta\) weakly constrained | Signal largely absorbed by background freedom |
| E2 \(\sigma_X=0.015,\ \theta=1.5\) | Amplitude recovered at the right order; \(\theta\) still poorly pinned (kernel degeneracy) | Recovery degrades once \(w_0,w_a\) are free |

**Lesson for the protocol:** always free \(\{w_0,w_a\}\) before claiming E1/E2; fixed-background fits can look artificially optimistic. Lag-shape constraints remain secondary to the full-kernel posterior on \(\theta\).

### 12.3 Products

| Path | Content |
|------|---------|
| `scripts/euclid_mock_mcmc.py` | Mock generation + emcee sampler |
| `results/euclid_mcmc/mcmc_summary.txt` | Latest run summary |
| `results/euclid_mcmc/mcmc_summaries.json` | Machine-readable posteriors |
| `figures/euclid_mcmc_*_fixedBG.png` | 1D posterior panels (fixed background) |
| `figures/euclid_mcmc_*_freeBG.png` | 1D posterior panels (free background) |

```bash
pip install emcee
python scripts/euclid_mock_mcmc.py
```


## References

[1] Morales Souhail, J., DESI OU/QNM analysis, this repository.  
[2] Morales Souhail, J., Fundamental vs emergent vacuum relaxation, this repository.  
[3] DESI Collaboration, arXiv:2503.14738 (2025).  
[4] Euclid Collaboration, DR1 BAO forecasts (as available at analysis time).
