# Constraints on Stochastic Dark Energy from Quantum Fluid Instabilities and DESI DR2 Baryon Acoustic Oscillations

Jesús Morales Souhail · github.com/jesus-morales-souhail · July 2026

ORCID: [0009-0000-7637-1818](https://orcid.org/0009-0000-7637-1818)  
Repository: https://github.com/jesus-morales-souhail/stochastic-dark-energy-ou

---

## Abstract

Recent DESI DR2 BAO data show no evidence for a stationary stochastic residual on top of a smooth background. Under the free-$\theta$ OU profile on seven isotropic $\alpha$ bins with full measurement covariance, the 95% CL upper limit is $\sigma_X < 2.5\times 10^{-2}$ (exact grid $2.506\times 10^{-2}$); the best fit is $\sigma_X\to 0$. Here I test a different hypothesis: a dark-energy fluid with negative effective mass ($m^{\ast}<0$) and quartic self-interaction. From a modified Gross–Pitaevskii equation I derive the Bogoliubov growth rate $\Gamma_{\max}=g|\psi_0|^2/\hbar\equiv 1/t_c$ and the induced rank-1 residual covariance. An MLE on the real DESI DR2 isotropic $\alpha$ vector shows that the likelihood peaks at no growth ($t_c\to\infty$, $\Delta\ln\mathcal{L}\to 0$). When growth is large enough that the rank-1 term matters, the model is excluded at $\Delta\ln\mathcal{L}\approx -12.7$ ($\Delta\chi^{2}\approx +25$) for a representative active $t_c$. A globally coherent tachyonic quantum fluid is incompatible with the observed BAO residuals under this pipeline.

---

## 1. Introduction

Baryon Acoustic Oscillation (BAO) measurements from DESI DR2 constrain the late-time expansion history. In the companion OU analysis in this repository, a stationary residual amplitude is driven toward zero. The free-$\theta$ profile 95% CL ceiling is $\sigma_X < 2.5\times 10^{-2}$ (7 bins, full meas. cov.) [2]; the null ($\sigma_X\to 0$) is preferred. That limits mechanisms that would source stochastic fluctuations in the dark-energy density.

One class of models that could, *a priori*, produce such fluctuations is a dark-energy fluid described as a quantum condensate with negative effective mass. In condensed-matter physics, Bose–Einstein condensates with engineered dispersion can have $m^{\ast}<0$ and undergo dynamical instabilities that amplify density fluctuations [3]. If the dark-energy sector shared that structure, it might generate a residual $\sigma_X$ visible in BAO data.

In this note I formalise that connection: Bogoliubov growth rate, induced variance, **rank-1** residual covariance, and an MLE on the real DESI DR2 isotropic $\alpha$ vector. The coherent growing-mode scenario is excluded whenever growth is active at the seed amplitude used in this programme; the preferred limit is $t_c\to\infty$ (no growth).

---

## 2. Theoretical Framework: The Quantum Fluid Action

I consider a complex scalar field $\psi$ whose macroscopically condensed state mimics dark energy. The action in an FLRW background is

$$
\mathcal{S}=\int\mathrm{d}^{4}x\,\sqrt{-g}\left[-\frac{1}{2}g^{\mu\nu}\partial_{\mu}\psi^{\dagger}\partial_{\nu}\psi-V(\psi)\right],
$$

with the self-interaction potential

$$
V(\psi)=\frac{1}{2}(m^{\ast})^{2}|\psi|^{2}+\frac{g}{2}|\psi|^{4},
$$

where $(m^{\ast})^{2}<0$ (tachyonic mass term) and $g>0$. In the mean-field approximation, the condensate wavefunction $\psi(\mathbf{x},t)$ obeys a modified Gross–Pitaevskii equation (GPE) in an expanding universe. Restoring $\hbar$ and writing the inertial mass as $m^{\ast}=-|m^{\ast}|$,

$$
\mathrm{i}\hbar\frac{\partial\psi}{\partial t}=\left[-\frac{\hbar^{2}}{2m^{\ast}a^{2}(t)}\nabla^{2}+V_{\mathrm{ext}}(\mathbf{x})+g|\psi|^{2}\right]\psi.
$$

For cosmological homogeneity I set $V_{\mathrm{ext}}=0$ and identify the background dark-energy density scale with the condensate density, $\rho_{\Lambda}\simeq|\psi_{0}|^{2}$ (mean-field bookkeeping; not a full stress-tensor derivation).

---

## 3. Bogoliubov Excitations and the Instability Growth Rate

I linearise around a homogeneous background, $\psi=\psi_{0}+\delta\psi$, and expand the density perturbation $\delta\rho_{\Lambda}=\psi_{0}^{\ast}\delta\psi+\psi_{0}\delta\psi^{\ast}$ in comoving Fourier modes $q$. With $m^{\ast}=-|m^{\ast}|$ the kinetic term flips sign relative to the usual positive-mass GPE. Defining the kinetic energy scale

$$
\varepsilon_{q}\equiv\frac{\hbar^{2}q^{2}}{2|m^{\ast}|a^{2}(t)},
$$

the Bogoliubov frequency (as an energy) satisfies

$$
\omega_{q}^{2}=\varepsilon_{q}^{2}-2\varepsilon_{q}\,g|\psi_{0}|^{2}.
$$

For $\varepsilon_{q}<2g|\psi_{0}|^{2}$ one has $\omega_{q}^{2}<0$ and a dynamical instability. The growth rate is

$$
\Gamma_{q}=\frac{1}{\hbar}\sqrt{\varepsilon_{q}\bigl(2g|\psi_{0}|^{2}-\varepsilon_{q}\bigr)}.
$$

Equivalently, for $\varepsilon_{q}<2g|\psi_{0}|^{2}$,

$$
\Gamma_{q}=\frac{q}{a^{2}\sqrt{2|m^{\ast}|}}\sqrt{2g|\psi_{0}|^{2}-\frac{\hbar^{2}q^{2}}{2|m^{\ast}|a^{2}}}.
$$

The fastest-growing mode sits at $\varepsilon_{q}=g|\psi_{0}|^{2}$, i.e.

$$
q_{\max}=\frac{a}{\hbar}\sqrt{2|m^{\ast}|\,g|\psi_{0}|^{2}},
$$

with maximum growth rate

$$
\Gamma_{\max}=\frac{g|\psi_{0}|^{2}}{\hbar}\equiv\frac{1}{t_{c}}.
$$

Here $t_{c}$ is the characteristic collapse / growth time of the quantum fluid. An older draft wrote $q_{\max}\propto 1/a$ instead of $\propto a$; that was bookkeeping. $\Gamma_{\max}$ is unchanged.

---

## 4. Stochastic Evolution and Variance of $\Omega_{\Lambda}$

In the linear regime,

$$
\delta\rho_{\Lambda}(q,t)=\delta\rho_{\Lambda}^{(0)}(q)\,e^{\Gamma_{q}t},
$$

where $\delta\rho_{\Lambda}^{(0)}$ is a primordial seed (Planck-scale granularity or Sorkin-type Poisson fluctuations [4,5]). The dimensionless residual variance is

$$
\sigma_{X}^{2}(t)\equiv\frac{\langle(\delta\rho_{\Lambda})^{2}\rangle}{\rho_{\Lambda}^{2}}=\frac{1}{\rho_{\Lambda}^{2}}\int\frac{\mathrm{d}^{3}q}{(2\pi)^{3}}P_{0}(q)\,e^{2\Gamma_{q}t}.
$$

A saddle-point estimate around $q_{\max}$ gives the working form used below,

$$
\sigma_{X}(t)\sim\sigma_{0}\,e^{t/t_{c}},
$$

with seed amplitude $\sigma_{0}\sim 10^{-61}$ (Bekenstein–Hawking / causal-set counting, same order used in the amplification-gap notes).

---

## 5. Observational Coupling with DESI DR2 BAO Data

BAO summary statistics constrain the isotropic dilation $\alpha(z)\equiv D_{V}(z)/r_{s}$ relative to a fiducial cosmology. The sensitivity of $D_{V}$ to $\Omega_{\Lambda}$ is

$$
S(z)\equiv\frac{\partial\ln D_{V}(z)}{\partial\Omega_{\Lambda}}.
$$

A single coherent growing mode does **not** produce an OU-stationary residual covariance. It produces a **rank-1** update. With residual $\mathbf{r}=\boldsymbol{\alpha}_{\mathrm{obs}}-\mathbf{1}$ and lookback time $t(z)$,

$$
C_{ij}=C^{\mathrm{meas}}_{ij}+\sigma_{0}^{2}S(z_{i})S(z_{j})\exp\bigl[(t(z_{i})+t(z_{j}))/t_{c}\bigr],
$$

where $C^{\mathrm{meas}}$ is the official Gaussian BAO measurement covariance projected to isotropic $\alpha$ (7×7; block-diagonal across bins in the public product). Equivalently, $C=C^{\mathrm{meas}}+\mathbf{u}\mathbf{u}^{\mathsf{T}}$ with $u_{i}=\sigma_{0}S_{i}e^{t_{i}/t_{c}}$.

The Gaussian log-likelihood is

$$
\ln\mathcal{L}=-\frac{1}{2}\left[\mathbf{r}^{\mathsf{T}}C^{-1}\mathbf{r}+\ln|C|+n\ln(2\pi)\right].
$$

### 5.1 Real DESI DR2 numbers (this repository)

I run the MLE on the real isotropic $\alpha$ vector from the public DESI DR2 Zenodo pack (`figure6/DESI_DR2_alpha_DV_over_rs.txt`; loader `scripts/desi_dr2_data.py`). Pipeline: `scripts/tachyonic_rank1_mle.py`. Output: `results/tachyonic_rank1/SUMMARY.txt`.

| Quantity | Value |
|:---------|:------|
| Bins $n$ | 7 |
| Data file | DESI DR2 $\alpha$ (real, not mock) |
| Measurement cov | full Gaussian BAO $13\times 13$ projected to $\alpha$ |
| $\ln\mathcal{L}_{\Lambda\mathrm{CDM}}$ | $24.8996$ |
| $\chi^{2}_{\Lambda\mathrm{CDM}}$ | $4.49$ |
| Seed $\sigma_{0}$ | $10^{-61}$ |
| No-growth limit $t_{c}\to\infty$ | $\Delta\ln\mathcal{L}\to 0$ (preferred) |
| Representative *active* $t_{c}$ | $\approx 0.074\,\mathrm{Gyr}$ |
| $\Delta\ln\mathcal{L}$ at that $t_{c}$ | $\approx -12.68$ |
| $\Delta\chi^{2}\approx -2\Delta\ln\mathcal{L}$ | $\approx +25.4$ |
| Deeper in the active window | $\Delta\ln\mathcal{L}$ can reach $\lesssim -30$ |

**Reading of the scan.** With $\sigma_{0}=10^{-61}$, most large $t_{c}$ leave the rank-1 term invisible (recovering $\Lambda$CDM). When $t_{c}$ is short enough that $\sigma_{0}e^{t/t_{c}}$ is no longer negligible compared with the BAO errors, the log-likelihood drops by $\mathcal{O}(10)$ or more. The data therefore do **not** prefer a finite growth time: the acceptable limit is $t_{c}\to\infty$ (no coherent growth).

The stationary OU profile ceiling $\sigma_{X}<2.5\times 10^{-2}$ is a **different** operator and must not be mixed into this rank-1 test.

---

## 6. Exclusion of the Coherent Quantum Fluid

Equating a stationary OU prior to this non-stationary growing mode is invalid. With the rank-1 covariance above, the data do not return a useful lower bound on $t_{c}$; they reject the region of parameter space where coherent growth is large enough to matter, and they are happy with $t_{c}\to\infty$.

I conclude that the late-time BAO residuals do not show the coherent, large-scale quantum granularity of a homogeneous tachyonic fluid. Any viable model must either (i) break the global rank-1 coherence (localised fluctuations), or (ii) keep $g$ so small that $t_{c}\to\infty$, i.e. reduce to a smooth background under this probe.

---

## 7. Discussion and Conclusion

I examined a dark-energy model based on a tachyonic quantum fluid with quartic self-interaction. The Bogoliubov analysis gives a clean maximum growth rate $\Gamma_{\max}=g|\psi_{0}|^{2}/\hbar$. The observationally relevant object is not that formula alone but the **rank-1 residual covariance** it induces on BAO $\alpha$ measurements.

On public DESI DR2 data (real $\alpha$ vector, 7 bins), that model is excluded whenever growth is active at the Sorkin seed scale used in this programme. The smooth measurement-cov $\Lambda$CDM residual description remains preferred. Condensed-matter analogues of negative-mass hydrodynamics remain interesting laboratory physics; they are not a description of late-time cosmological vacuum residuals at DESI precision under this pipeline.

Future surveys (e.g. Euclid) will tighten residual tests further, but this coherent rank-1 class is already ruled out as a DESI-scale explanation of BAO residuals.

---

## Acknowledgements

I thank the DESI collaboration for making the DR2 BAO data publicly available.

---

## References

[1] DESI Collaboration, "DESI DR2 Results II: Measurements of BAO and Cosmological Constraints," arXiv:2503.14738 (2025).

[2] Morales Souhail, J., this repository: `papers/stochastic-dark-energy-desi-dr2.md`, `scripts/profile_sigma_x_desi.py` (2026). Profile 95% CL $\sigma_X<2.5\times 10^{-2}$; null preferred.

[3] Khamehchi, M. A. et al., "Negative-Mass Hydrodynamics in a Spin-Orbit–Coupled Bose–Einstein Condensate," Phys. Rev. Lett. **118**, 155301 (2017).

[4] Sorkin, R. D., "Is the Cosmological 'Constant' a Nonlocal Quantum Residual?," arXiv:gr-qc/0503057 (2005).

[5] Bekenstein, J. D., "Black Holes and Entropy," Phys. Rev. D **7**, 2333 (1973).

[6] DESI Collaboration, "DESI 2024 III: Baryon Acoustic Oscillations from Galaxies and Quasars," arXiv:2404.03000 (2024).

[7] Planck Collaboration, "Planck 2018 results. VI. Cosmological parameters," Astron. Astrophys. **641**, A6 (2020).

---

## Addendum: Principle of Late-Time Vacuum Homogeneity

This exclusion is consistent with the broader Principle of Late-Time Vacuum Homogeneity: at observable cosmological scales and with current BAO precision, the dark-energy sector behaves as a homogeneous, non-stochastic background under the tests in this repository. See `principle-of-vacuum-smoothness.md`.

## Reproducibility

```bash
python scripts/tachyonic_rank1_mle.py
# writes results/tachyonic_rank1/SUMMARY.txt and SUMMARY.json
```
