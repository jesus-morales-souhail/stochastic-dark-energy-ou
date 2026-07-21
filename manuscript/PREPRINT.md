# Constraints on Stochastic Dark Energy from DESI DR2 BAO: A Null Result and the \(10^{56}\) Amplification Bottleneck

**Author:** Jesús Morales Souhail  
**ORCID:** [0009-0000-7637-1818](https://orcid.org/0009-0000-7637-1818)  
**Email:** jmskjym@gmail.com  
**Affiliation:** Independent researcher  
**Date:** July 2026  
**Status:** Preprint — **not peer reviewed**  
**Code & data pipeline:** https://github.com/jesus-morales-souhail/stochastic-dark-energy-ou  
**License:** Code MIT · Text CC BY 4.0  

---

## Abstract

We test whether public DESI DR2 baryon acoustic oscillation (BAO) summary statistics require a *stationary stochastic* component in the late-time dark-energy sector, on top of a smooth background. Fluctuations are modelled as an Ornstein–Uhlenbeck (OU) process in logarithmic scale factor \(x=\ln a\), with a quasi-normal-mode (QNM) oscillatory extension. Using a Gaussian likelihood with an additive residual covariance kernel, maximum-likelihood estimation drives the stochastic amplitude to the numerical floor. We report the working 95% confidence upper limit

\[
\sigma_X < 1.5\times 10^{-4}
\]

under the stated phenomenological kernel and public BAO summary statistics (diagonal measurement errors in the baseline pipeline). Nested extensions that free \(\{w_0,w_a,\sigma_X,\theta\}\) are not preferred over a smooth background (\(\Delta\mathrm{AIC}\approx +4\)).

Separately, a *coherent* tachyonic growth model with the correct rank-1 covariance induced by a single growing mode is strongly disfavoured: \(\Delta\ln\mathcal{L}\approx -11.35\) relative to \(\Lambda\)CDM (\(\Delta\chi^2\approx +22.7\)).

We then quantify an *amplification bottleneck*: a Sorkin–Bekenstein Poisson seed \(\sigma_0\sim 10^{-61}\) lies \(\sim 10^{56}\) below residual amplitudes of interest for Euclid-scale BAO (\(\sim 10^{-5}\)). Audited linear amplifiers are excluded, short by many decades, or redefine the seed. Machine scans of three allowed non-linear cards show that late-time freeze-out alone yields unit gain, a soft double-well avalanche yields gain \(\sim 2\), and only a *local* effective count \(N_{\mathrm{eff}}\sim 4\times 10^{7}\) (DESI) / \(10^{10}\) (Euclid) reaches target amplitudes—by redefining the microphysical seed, not by free lunch.

Gravitational slip \(\gamma=\Phi/\Psi\) is the correct *operator* for anisotropic leakage past volume-preserving (SDiff) structure, but inherits the same amplitude starvation: even full anisotropy of the BAO-bounded residual gives \(|\gamma-1|\sim 10^{-4}\), \(\mathcal{O}(10^{2}\)–\(10^{3})\) below current and forecast floors.

**We claim limits, exclusions, and a map of open hypotheses—not a detection of vacuum discreteness or a new dark-energy fluid.**

---

## 1. Introduction

### 1.1 Motivation

DESI BAO measurements have sharpened late-time distance constraints and, in combination with CMB and supernovae, have stimulated renewed interest in dynamical dark energy (DESI Collaboration, arXiv:2503.14738 and companions). Independently of whether the *background* equation of state departs from \(w=-1\), a logically distinct question is whether *residuals* about a smooth expansion history require a stochastic dark-energy component.

A finite information bound on the observable universe (\(N\sim 10^{122}\) Planck units) has long motivated order-of-magnitude vacuum fluctuations \(\delta\Lambda\sim 1/\sqrt{N}\) (Sorkin-type arguments in unimodular / causal-set settings). That seed is many orders of magnitude below BAO precision. The scientific programme is therefore:

1. **Data:** constrain the *effective* residual amplitude under a minimal, falsifiable kernel.  
2. **Models:** exclude coherent amplification mechanisms with the *correct* covariance.  
3. **Theory hygiene:** state which amplifiers are closed and which open cards remain—without smuggling free parameters as “microphysics.”

### 1.2 Scope and non-claims

This paper and the associated public repository address (1)–(3) using public BAO summary statistics and open analysis code. We do **not** claim:

- detection of Planck-scale vacuum noise;  
- that Euclid will see unamplified Sorkin seeds;  
- a unique identification of unimodular / SDiff gravity;  
- peer-reviewed status of this preprint;  
- results from optical laboratory metaphors (those live in a separate exploratory repository; see §8).

### 1.3 Reading guide

| Section | Content |
|---------|---------|
| §2 | Data and likelihood |
| §3 | Stationary OU / QNM null and \(\sigma_X\) limit |
| §4 | Exclusion of coherent tachyonic growth |
| §5 | Amplification bottleneck and route scans |
| §6 | Slip / SDiff structural remark (amplitude honesty) |
| §7 | Limitations |
| §8 | Software, repositories, and claim boundaries |
| §9 | Conclusions |

---

## 2. Data and method

### 2.1 Dataset

We use **public DESI DR2 BAO summary statistics** (arXiv:2503.14738): isotropic BAO \(\alpha\) measurements in seven effective redshift bins, with published measurement uncertainties. The baseline pipeline analysed here adopts **diagonal** measurement errors in the residual likelihood. Full off-diagonal BAO covariances are a natural extension for tightening the limit (§7); they are not required to establish the *null preference* of the stochastic extension under the kernels tested.

**Important separation:** DESI analyses that combine BAO with CMB+SN often prefer dynamical \(w_0\)–\(w_a\). Our BAO-*only* residual analysis uses a CPL background near \(\Lambda\)CDM (\(w_0\approx -0.99\), \(w_a\approx -0.02\)) as the smooth reference for \(\alpha=1\) under the fiducial cosmology of the published \(\alpha\) products. The stochastic question is about *residuals*, not a re-fit of the full multi-probe DE equation of state.

### 2.2 Residual process

Let \(X(x)\equiv \delta\Omega_X(x)\) with \(x=\ln a\). The stationary OU process is

\[
\mathrm{d}X = -\theta X\,\mathrm{d}x + \sigma\,\mathrm{d}W_x,
\qquad
\mathrm{Var}(X)=\frac{\sigma^2}{2\theta}\quad (\theta>0).
\]

We work with an effective residual amplitude \(\sigma_X\) entering the BAO residual covariance. For bins \(i,j\) with logarithmic separation \(\Delta x_{ij}\),

\[
C^{\mathrm{OU}}_{ij}
=
\sigma_X^2\,\mathrm{e}^{-\theta\,\Delta x_{ij}},
\qquad
C^{\mathrm{QNM}}_{ij}
=
\sigma_X^2\,\mathrm{e}^{-\theta\,\Delta x_{ij}}\cos(\omega_R\,\Delta x_{ij}).
\]

The total covariance is \(C_{\mathrm{tot}}=C_{\mathrm{meas}}+C_{\mathrm{stoch}}\). The Gaussian log-likelihood is standard:

\[
-2\ln\mathcal{L}
=
\ln\det C_{\mathrm{tot}}
+
\mathbf{r}^{\mathsf T}C_{\mathrm{tot}}^{-1}\mathbf{r}
+
\mathrm{const},
\]

with residual vector \(\mathbf{r}\) relative to the smooth background prediction.

### 2.3 Nested models and information criteria

| Model | Stochastic content | Free stochastic params |
|-------|--------------------|------------------------|
| \(\Lambda\)CDM baseline | \(C_{\mathrm{meas}}\) only | — |
| H0 (OU) | \(C_{\mathrm{OU}}(\theta,\sigma_X)\) | 2 |
| H1 (QNM) | \(C_{\mathrm{QNM}}(\theta,\sigma_X,\omega_R)\) | 3 |
| Nested CPL+\((\sigma,\theta)\) | background free + residual | see §3 |

AIC/BIC differences use the usual penalties on the number of free parameters.

### 2.4 Software

Primary scripts (repository root):

| Script | Role |
|--------|------|
| `scripts/ou_bao_likelihood.py` | OU / QNM residual MLE |
| `scripts/eos_efectiva.py` | CPL background + nested \((\sigma,\theta)\) |
| `scripts/profile_sigma_x_desi.py` | Profile likelihood for \(\sigma_X\) |
| `scripts/joint_w0wa_sigma_desi.py` | Joint \(\{w_0,w_a,\sigma_X\}\) |
| `scripts/gpe/` | Tachyonic / GPE exclusion pipeline |
| `scripts/amplifier_audit.py` | Amplification audit |
| `scripts/amplification/route{1,2,3}_*.py` | Non-linear cards |
| `scripts/slip_bridge.py` | \(\sigma_X\to|\gamma-1|\) OOM map |

---

## 3. Results I — Stationary stochastic null

### 3.1 Maximum-likelihood point

Under fixed smooth background and public BAO summary statistics, the MLE drives the stochastic amplitude to the numerical floor:

| Model | \(\theta\) | \(\sigma_X\) | \(\omega_R\) | \(\Delta\ln\mathcal{L}\) vs \(\Lambda\)CDM | AIC | BIC |
|-------|------------|--------------|--------------|---------------------------------------------|-----|-----|
| OU (H0) | \(0.001\) (boundary) | \(5\times 10^{-5}\) | — | \(0.00\) | \(-50.03\) | \(-50.14\) |
| QNM (H1) | \(0.001\) | \(5\times 10^{-5}\) | \(0.00\) | \(0.00\) | \(-48.03\) | \(-48.19\) |

**Interpretation:** the data do not require a stationary stochastic residual. The preferred stochastic amplitude is consistent with zero within the resolution of the likelihood surface; reported floor values reflect optimizer bounds, not a detection.

### 3.2 Working upper limit

From the residual analysis (profile / likelihood ratio under the OU kernel as implemented in the public pipeline):

\[
\boxed{\sigma_X < 1.5\times 10^{-4}\quad (95\%~\mathrm{CL})}
\]

This is an **upper limit on an effective residual amplitude**, not a measurement of a Planck-scale seed.

### 3.3 Nested background freedom

When \(\{w_0,w_a\}\) are freed together with residual parameters, the stochastic extension is **not preferred** (\(\Delta\mathrm{AIC}\approx +4\)). Marginalising over a flexible smooth background drives \(\sigma_X\) toward zero: extra residual variance is absorbed by background freedom when the number of BAO bins is small (\(N_{\mathrm{bin}}=7\)).

### 3.4 Claim statement (Act I)

> Under the phenomenological OU/QNM residual kernels and public DESI DR2 BAO summary statistics, there is **no statistical evidence** for stationary stochastic dark-energy fluctuations. The data prefer smooth evolution; \(\sigma_X\) is bounded.

---

## 4. Results II — Coherent tachyonic growth excluded

A distinct hypothesis is that a dark-energy quantum fluid with negative effective mass undergoes Bogoliubov instability, yielding a *coherent* growing mode \(\sigma_X(t)=\sigma_0\,e^{t/t_c}\). The induced BAO residual covariance is **rank-1**, not OU-stationary. Using that covariance in an MLE on the same BAO summary statistics yields

\[
\Delta\ln\mathcal{L}
=
\ln\mathcal{L}(t_c)-\ln\mathcal{L}_{\Lambda\mathrm{CDM}}
=
-11.35
\qquad
\bigl(\Delta\chi^2\approx +22.7\bigr).
\]

Any finite collapse / growth time \(t_c\) is disfavoured; the acceptable limit is \(t_c\to\infty\) (no growth).

**Claim statement (model kill):**

> A globally coherent tachyonic growing mode with rank-1 residual covariance is **incompatible** with DESI DR2 BAO residuals under the pipeline stated in the repository.

This exclusion is logically independent of the stationary OU bound: the covariance structure is different, and the correct structure must be used.

---

## 5. Results III — Amplification bottleneck (no free lunch)

### 5.1 Quantified gap

| Quantity | Value |
|----------|-------|
| Motivational UV seed (Sorkin / Bekenstein counting) | \(\sigma_0\sim 10^{-61}\) |
| DESI residual bound (this work) | \(\sigma_X<1.5\times 10^{-4}\) (95% CL) |
| Euclid-scale residual of interest | \(\sim 10^{-5}\) |
| Required gain \(\sigma_{\mathrm{target}}/\sigma_0\) | \(\sim 10^{56}\)–\(10^{57}\) |

DESI constrains the **effective** amplitude after any amplification. Without a physical map \(\sigma_0\to A_0\), “Euclid will see vacuum noise” is not a theorem of this analysis.

### 5.2 Linear amplifier audit

| Class | Schematic gain | Status |
|-------|----------------|--------|
| Coherent / tachyonic growth | exponential | **Excluded** (§4) |
| Open-system desqueezing (\(r\sim 1.5\)) | \(e^{2r}\sim 20\) | \(\sim 55\) decades short |
| Naive \(\sqrt{N}\) boost | \(\sim 10^{61}\) | Double-counts: seed already \(\propto 1/\sqrt{N}\) unless coherent → reduces to excluded class |
| Inflation-style \(e^{60}\) | \(\sim 10^{26}\) | Still short if forced from \(\sigma_0\); natural inflationary seed is \(\sim H/M_{\mathrm{Pl}}\), **not** \(10^{-61}\) |
| Large phase jump | discrete feature | Tensions with smooth BAO residuals if large |

### 5.3 Three non-linear cards (machine scans)

| Route | Idea | Numerical verdict |
|-------|------|-------------------|
| **1. Local \(N_{\mathrm{eff}}\)** | \(\sigma_{0,\mathrm{eff}}=1/\sqrt{N_{\mathrm{eff}}}\) with \(N_{\mathrm{eff}}\ll N_{\mathrm{BH}}\) | \(N_{\mathrm{eff}}\sim 4.44\times 10^{7}\) (DESI), \(1.00\times 10^{10}\) (Euclid). **Only card that reaches target amplitude**—by **redefining the seed**. |
| **2. Late freeze-out** | \(\theta\to 0\) after late horizon exit, \(\Delta x\sim\mathcal{O}(1)\) | Freeze/restore gain \(=1.000\); \(\mathrm{rms}/\sigma=\mathcal{O}(1)\). **Dead as amplifier alone.** |
| **3. Soft avalanche** | double-well Langevin trigger | Gain \(\mathrm{p95}/\sigma\sim 2\); **0/288** jobs with \(\mathrm{p95}\ge 10^{-5}\) and \(\sigma\le 10^{-8}\). BAO-safe if \(\sigma\lesssim 4\times 10^{-5}\) in the scanned grid. **No free lunch.** |

Full tables: `results/amplification_routes/VERDICT.md` and `papers/amplification-no-free-lunch.md` §6.3 in the code repository.

### 5.4 Claim statement (Act III)

> There is **no free lunch** for linear amplification of a pure Sorkin seed to BAO/Euclid residual scales. Publication of the present corpus is consistent with **honest limits + exclusions + a short map of open hypothesis cards**, each with an explicit theoretical price.

---

## 6. Structural remarks — SDiff and gravitational slip

### 6.1 Geometry as interpretation, not detection

Volume-preserving diffeomorphisms / unimodular structure project out isotropic vacuum stress of the form \(T_{\mu\nu}=V(x)\,g_{\mu\nu}\). That is a *candidate explanation* of late-time residual smoothness, not a particle detection and not a substitute for the BAO likelihood.

### 6.2 Option 0: anisotropic gap is real in principle, starved in amplitude

Anisotropic stress is **not** of the form \(V g_{\mu\nu}\). Sub-horizon slip satisfies (Newtonian gauge)

\[
|\gamma-1|
=
\left|\frac{\Phi-\Psi}{\Psi}\right|
=
2\,\varepsilon\,\sigma_X\,\frac{\rho_X}{\rho_m\,|\delta_m|},
\]

with \(\varepsilon\in[0,1]\) the anisotropic fraction of the residual. With \(\sigma_X=1.5\times 10^{-4}\), \(\varepsilon=1\), \(\delta_m=1\):

| \(z\) | max \(|\gamma-1|\) | vs Maus \(\gamma=1.17\pm 0.11\) (\(|\gamma-1|\sim 0.17\)) | vs Sakr-like floor \(\sim 0.05\) |
|-------|---------------------|----------------------------------------------------------|----------------------------------|
| 0.5 | \(1.93\times 10^{-4}\) | \(\sim 880\times\) short | \(\sim 260\times\) short |
| 1.0 | \(8.16\times 10^{-5}\) | \(\sim 2\times 10^{3}\) | \(\sim 610\times\) |
| 1.5 | \(4.18\times 10^{-5}\) | \(\sim 4\times 10^{3}\) | \(\sim 1.2\times 10^{3}\) |

(Maus et al., arXiv:2505.20656; forecasts Sakr et al., arXiv:2501.07477; runnable map: `scripts/slip_bridge.py`.)

**Claim statement (Act IV):**

> Slip is the correct *kind* of operator for anisotropic leakage past SDiff, but with the BAO-bounded residual amplitude it is **not** an observational shortcut. It inherits the amplification problem.

We do **not** implement a homemade Boltzmann hierarchy. Future work, if amplitudes become reachable, should use community codes (e.g. hi_class / MGCAMB) with a pre-specified target.

---

## 7. Limitations (required for peer review)

1. **Summary statistics:** baseline residual analysis uses published BAO \(\alpha\) and diagonal errors. Full DESI covariance matrices may change the numerical factor in \(\sigma_X\) while preserving the null preference if residuals remain consistent with noise.  
2. **Kernel dependence:** bounds are for the OU/QNM residual kernels stated; other residual models need their own covariance and MLE.  
3. **Background degeneracy:** with only seven BAO bins, flexible \(w(z)\) can absorb residual variance (§3.3). Multi-probe fits are a separate analysis.  
4. **Cross-correlation (DESI galaxies × SN residuals):** an exploratory angular correlation test exists in the repository (`scripts/cross_correlation_DESI.py`) but is **preliminary** (systematics not fully controlled; DR1 weights). It is **not** part of the primary claim set of this preprint.  
5. **Sorkin seed:** used as *motivational UV scale* for the amplification discussion, not as a fitted parameter.  
6. **Independent research:** this work has not undergone peer review at the time of writing.

---

## 8. Software, repositories, and claim boundaries

### 8.1 Primary scientific repository

**https://github.com/jesus-morales-souhail/stochastic-dark-energy-ou**

Contains the BAO pipelines, amplification audits, slip bridge, figures, and technical notes supporting §§2–6.

### 8.2 Exploratory / speculative repository (separate)

Method-hygiene digressions, optical analogies, band-limit / superoscillation demos, and other **non-cosmological** scale–operator checks are **not** part of the peer-review claim set. They are maintained separately:

**https://github.com/jesus-morales-souhail/stochastic-de-exploratory-notes**  
*(local path if not yet pushed: `~/Proyectos/01_Fisica_y_Cosmologia/stochastic-de-exploratory-notes/`)*

| In primary repo (claims) | In exploratory repo (not claims) |
|--------------------------|----------------------------------|
| BAO OU/QNM limits | Car–drone pupil diffraction pedagogy |
| Tachyonic exclusion | Tesseract / \(B_4\) optical combinatorics |
| Amplification audit + routes 1–3 | Wavefront-shaping vs OU vacuum essays |
| Slip amplitude map | Superoscillation energy-tax demos |
| SDiff as *interpretation* | Maxwell→device→\(M\) lab notes |

**Rule:** if a result does not constrain cosmology at the correct scale and operator, it does not enter the abstract or §3–6.

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
```

---

## 9. Conclusions

1. **Null residual:** DESI DR2 BAO summary statistics do not favour stationary OU/QNM stochastic dark energy; \(\sigma_X < 1.5\times 10^{-4}\) (95% CL) under the stated kernel.  
2. **Model kill:** coherent tachyonic growth with rank-1 covariance is excluded (\(\Delta\chi^2\approx +22.7\)).  
3. **No free lunch:** the gap from \(\sigma_0\sim 10^{-61}\) to residual targets \(\sim 10^{-5}\) is \(\sim 10^{56}\); audited linear amplifiers fail; freeze-out and soft avalanche scans do not provide free gain; only redefining \(N_{\mathrm{eff}}\) hits target amplitudes.  
4. **Slip honesty:** anisotropic leakage past SDiff is a real structural crack in principle and an amplitude-starved channel in practice.  
5. **Scientific posture:** the strength of this programme is **what it closes**—limits, exclusions, and an honest map of open theory prices—not a promised detection.

Future data (Euclid BAO multi-bin products; full DESI covariances) can tighten \(\sigma_X\) with the same residual likelihood. Future *theory* work must pay the price of Routes 1 or 3 (or an inflation-style seed redefinition)—not invent linear gain on white noise.

---

## Acknowledgments

This work uses publicly released DESI DR2 BAO products and community literature on gravitational slip and forecasts (Maus et al.; Sakr et al.). Any errors are the author’s alone.

---

## Data availability

- DESI DR2 BAO: https://data.desi.lbl.gov/public/ · arXiv:2503.14738  
- Analysis code: https://github.com/jesus-morales-souhail/stochastic-dark-energy-ou  
- Numerical amplification tables: `results/amplification_routes/`  

---

## References (minimal set for the claim spine)

1. DESI Collaboration, *DESI DR2 Results II: Measurements of Baryon Acoustic Oscillations*, arXiv:2503.14738.  
2. Maus et al., *Joint 3D clustering and galaxy × CMB-lensing with DESI DR1*, arXiv:2505.20656.  
3. Sakr, Zheng & Casas, forecasts on anisotropic stress / \(\eta\), arXiv:2501.07477.  
4. Ma & Bertschinger, *Cosmological Perturbation Theory in the Synchronous and Conformal Newtonian Gauges*, Astrophys. J. Suppl. (1995), arXiv:astro-ph/9506072.  
5. Clifton, Ferreira, Padilla & Skordis, *Modified Gravity and Cosmology*, Phys. Rep. (2012).  
6. Sorkin, causal-set / unimodular vacuum fluctuation arguments (see repository notes for curated pointers).  
7. Morales Souhail, J., this repository and companion notes (July 2026).  

*End of preprint manuscript.*
