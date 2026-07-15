# Desqueezing Relaxation Time from Open Quantum Systems and its Mapping to Late-Time Vacuum Fluctuations in Stochastic Dark Energy

**Author:** Jesús Morales Souhail  
**Date:** July 2026  
**ORCID:** 0009-0000-7637-1818  
**Repository:** https://github.com/jesus-morales-souhail/stochastic-dark-energy-ou  
**Companion numerics:** open quantum simulations (`QuTiP`), project `simulaciones_cuanticas`  
**Status:** Research note — not peer reviewed  
**Type:** Self-contained synthesis note (~paper-length)

---

## Abstract

We present a self-contained research note that closes the loop between (i) a controlled open-quantum-system simulation of **squeezed-state desqueezing** under a thermal Lindblad bath and (ii) the phenomenological **Ornstein–Uhlenbeck (OU)** description of late-time vacuum fluctuations used in our DESI DR2 stochastic dark-energy analysis. The laboratory sector yields a clean half-life for the anomalous correlator,

$$
t_{1/2}\!\left(\lvert\langle a^2\rangle\rvert\right)=\frac{\ln 2}{\gamma},
$$

independent of the bath occupation \(n_{\rm th}\) over the scanned range. The cosmological sector identifies the physical damping rate with the repository continuity structure,

$$
\Gamma_{\rm phys}(z)=\theta\,H(z),
$$

so that at the present epoch \(\gamma\leftrightarrow\theta H_0\). Using **only** amplitudes and limits already fixed in the stochastic-dark-energy-ou corpus — in particular the Sorkin–Bekenstein seed \(\sigma_0\sim 10^{-61}\) and the DESI bound \(\sigma_X<1.5\times 10^{-4}\) (95% CL) — we show that pure Poisson vacuum noise remains unobservable at BAO precision. Euclid DR1 can probe **emergent** vacuum smoothness only if an effective kick amplitude \(A_0\gtrsim 10^{-5}\) arises after amplification or from non-Sorkin sources. We discuss implications for volume-preserving diffeomorphisms (SDiff) as a fundamental versus emergent symmetry, and we state explicitly what the null result does and does not constrain.

---

## 1. Introduction

### 1.1 The observational fact

Public DESI DR2 baryon acoustic oscillation (BAO) measurements, analysed within an OU / quasi-normal-mode (QNM) covariance framework, drive the amplitude of stationary stochastic fluctuations of the dark-energy density to the numerical floor. The working 95% confidence upper limit is

$$
\sigma_X < 1.5\times 10^{-4},
$$

where \(X\equiv\delta\Omega_\Lambda/\bar{\Omega}_\Lambda\) (or an equivalent fractional density contrast in the notation of the companion papers). This is not merely a non-detection: it is a structural constraint on how “noisy” the late-time vacuum is allowed to be at BAO scales.

### 1.2 The theoretical fork

Two interpretations of vacuum protection compete inside the same dataset:

1. **Fundamental SDiff.** The gauge group of gravity is reduced from \(\mathrm{Diff}(M)\) to volume-preserving diffeomorphisms \(\mathrm{SDiff}(M)\) (unimodular gravity). Local vacuum stresses of the form \(T_{\mu\nu}=V(x)g_{\mu\nu}\) are algebraically projected out of the trace-free field equations. Then \(\sigma_X=0\) exactly.
2. **Emergent SDiff.** The same protection arises only as a low-energy, information-theoretic effective symmetry. Residual fluctuations and a **finite relaxation time** after a vacuum perturbation are then logically expected.

### 1.3 Why open quantum systems

A squeezed vacuum is a controlled, pure-state perturbation away from the ground state. A thermal Lindblad bath restores a thermal fixed point. The half-life of the anomalous correlator \(\lvert\langle a^2\rangle\rvert\) is a precise operational definition of “how long non-smoothness lasts.” Mapping that half-life onto the OU mean-reversion parameter \(\theta\) converts a laboratory law into a cosmological correlation time.

### 1.4 Scope and honesty about amplitudes

This note uses:

- the QuTiP desqueezing law measured in this project;
- the OU continuity structure and numerical anchors of the stochastic-dark-energy-ou repository;
- **no** free “illustrative” kick amplitudes of order \(10^{-2}\).

The bare UV seed is the repository value \(\sigma_0\sim 10^{-61}\).

---

## 2. Open-system setup and the desqueezing law

### 2.1 Model

We consider a single damped harmonic oscillator in the Markovian thermal Lindblad equation

$$
\dot\rho=-i[H,\rho]
+\gamma(n_{\rm th}+1)\,\mathcal{D}[a]\rho
+\gamma n_{\rm th}\,\mathcal{D}[a^\dagger]\rho,
$$

with

$$
H=\omega a^\dagger a,
\qquad
\mathcal{D}[L]\rho=L\rho L^\dagger-\tfrac12\{L^\dagger L,\rho\}.
$$

The initial state is a squeezed vacuum,

$$
\lvert\psi_0\rangle=S(r)\lvert 0\rangle,
\qquad
S(r)=\exp\!\left[\tfrac{r}{2}\big(a^2-(a^\dagger)^2\big)\right],
$$

with fiducial parameters validated against Hilbert-space truncation:

$$
N=150,\qquad r=1.5,\qquad \omega=2.4
$$

(in simulation units). The thermal occupation \(n_{\rm th}\) is scanned independently of \(\gamma\).

### 2.2 Observables

We track:

- mean occupation \(\langle n\rangle=\langle a^\dagger a\rangle\);
- position variance \(\mathrm{Var}(q)\) with \(q=(a+a^\dagger)/\sqrt{2}\);
- the anomalous correlator \(\lvert\langle a^2\rangle\rvert\).

The first two quantities oscillate because free evolution rotates the squeeze ellipse in the \((q,p)\) plane. The modulus \(\lvert\langle a^2\rangle\rvert\) isolates the **loss of squeezing** without that kinematic rotation.

### 2.3 Result

Across \(\gamma\in\{0.6,1.2,2.4,3.6,4.8\}\) and \(n_{\rm th}\in[0,1]\),

$$
\boxed{
t_{1/2}\!\left(\lvert\langle a^2\rangle\rvert\right)
=\frac{\ln 2}{\gamma}
\qquad\text{with}\qquad
\gamma\,t_{1/2}=\ln 2
}
$$

to better than \(0.05\%\) relative accuracy. **Independence of \(n_{\rm th}\)** is a central result: the bath temperature sets the late-time floor

$$
\mathrm{Var}(q)_{\rm ss}=n_{\rm th}+\tfrac12,
$$

but not the half-life of the anomalous correlator.

Thus the laboratory export is simply

$$
t_{\rm relax}\propto\frac{1}{\gamma}.
$$

---

## 3. Stochastic dark energy as an OU process (repository structure)

### 3.1 Axioms (compressed)

From the main OU analysis:

- **A1 (finite information):** the observable universe has effective Hilbert-space dimension set by the Bekenstein–Hawking scale, \(N\sim 10^{122}\).
- **A2 (Sorkin seed):** discreteness induces Poisson fluctuations \(\delta\Lambda\sim 1/\sqrt{N}\), i.e.

$$
\sigma_0\sim 10^{-61}.
$$

- **A3 (OU closure):** fractional fluctuations \(X(x)\) in logarithmic scale factor \(x=\ln a\) obey

$$
dX=-\theta X\,dx+\sigma\,dW_x,
$$

with stationary variance \(\mathrm{Var}(X)=\sigma^2/(2\theta)\) when \(\theta>0\).

### 3.2 Continuity form and physical rate

Writing the OU process in cosmic time via \(dx=H\,dt\),

$$
\frac{dX}{dt}=-\theta H(t)\,X+\xi(t).
$$

Therefore

$$
\boxed{\Gamma_{\rm phys}(t)=\theta\,H(t)=\theta\,H(z).}
$$

This is not an extra postulate: it is the chain rule applied to the repository definition of \(\theta\).

### 3.3 BAO observable consequence

Induced BAO covariance between redshift bins takes the schematic form

$$
(C_{\rm OU})_{ij}=S(z_i)S(z_j)\,\mathrm{Cov}[X(x_i),X(x_j)],
$$

with

$$
\mathrm{Cov}[X(x_i),X(x_j)]=\frac{\sigma^2}{2\theta}\,e^{-\theta\lvert x_i-x_j\rvert},
$$

and sensitivity kernel \(S(z)=\partial\ln D_V/\partial\Omega_\Lambda\) computed on the flat \(\Lambda\)CDM fiducial

$$
H_0=67.4\,\mathrm{km\,s^{-1}\,Mpc^{-1}},\qquad\Omega_m=0.315.
$$

The DESI path width in the repository is

$$
\Delta x_{\rm DESI}=0.94
\quad(z\in[0.295,2.33]).
$$

### 3.4 Empirical anchors from the repository

| Quantity | Value | Status |
|----------|------:|--------|
| \(\sigma_X\) (95% CL) | \(<1.5\times 10^{-4}\) | DESI DR2 working limit |
| MLE floor | \(\theta\sim 10^{-3}\), \(\sigma_X\sim 5\times 10^{-5}\) | numerical boundary |
| Superseded script calibration | \(\theta=1.2\), \(\sigma_X\simeq 0.023\) | **not** preferred after free MLE |
| Euclid target / window | \(\sim 10^{-5}\); detection \(0<\sigma_X\lesssim 10^{-4}\) | forecast / discrimination |

---

## 4. The \(\gamma\leftrightarrow\theta H(z)\) dictionary

### 4.1 Identification

The Lindblad rate \(\gamma\) has dimension of inverse time. The natural identification with the OU sector is

$$
\gamma\;\longleftrightarrow\;\Gamma_{\rm phys}(z)=\theta\,H(z).
$$

At the present epoch,

$$
\gamma=\theta H_0
\qquad\Rightarrow\qquad
t_{1/2}=\frac{\ln 2}{\theta H_0}.
$$

Along a finite expansion history, the survival of a frozen kick of amplitude \(A_0\) is

$$
\sigma_{\rm res}(x)=A_0\,e^{-\theta\Delta x},
\qquad
\Delta x=\ln\!\frac{1+z_{\rm past}}{1+z_{\rm now}}.
$$

### 4.2 What “\(\gamma\sim H_0\)” really means

One should **not** impose \(\gamma=H_0\) by hand. That equality holds if and only if \(\theta\sim 1\). In particular, for a half-life equal to the Hubble time \(t_{H_0}=H_0^{-1}\simeq 14.51\,\mathrm{Gyr}\),

$$
\theta=\ln 2\simeq 0.693,
\qquad
\gamma=(\ln 2)\,H_0.
$$

For a half-life equal to the \(\Lambda\)CDM age \(t_0\simeq 13.8\,\mathrm{Gyr}\) at the paper fiducial, \(\theta\simeq 0.73\).

### 4.3 Repository \(\theta\) values translated to \(t_{1/2}\)

| Case | \(\theta\) | \(\gamma=\theta H_0\) [Gyr\(^{-1}\)] | \(t_{1/2}\) [Gyr] | \(\Delta x_{1/2}=\ln 2/\theta\) |
|------|----------:|-------------------------------------:|------------------:|-------------------------------:|
| MLE floor | \(10^{-3}\) | \(6.9\times 10^{-5}\) | \(\sim 1.0\times 10^{4}\) | \(693\) |
| Superseded \(\theta=1.2\) | \(1.2\) | \(8.3\times 10^{-2}\) | \(8.4\) | \(0.58\) |
| Half-decay on DESI path | \(0.737\) | \(5.1\times 10^{-2}\) | \(13.6\) | \(0.94\) |
| Hubble half-life | \(0.693\) | \(4.8\times 10^{-2}\) | \(14.5\) | \(1\) |

**Key reading:** the MLE floor implies **essentially no damping** across the DESI path,

$$
e^{-10^{-3}\times 0.94}\simeq 0.999.
$$

Hubble-scale restoration of smoothness requires \(\theta\sim\mathcal{O}(1)\).

### 4.4 Page and scrambling scales (\(S\sim 10^{122}\))

| Scale | \(T\) [Gyr] | \(\theta=\ln 2/(T H_0)\) |
|-------|------------:|-------------------------:|
| \(t_H\ln S\) (scrambling) | \(\sim 4\times 10^{3}\) | \(\sim 2.5\times 10^{-3}\) |
| \(t_H\sqrt{S}\) | \(\sim 10^{62}\) | \(\sim 10^{-62}\) |
| \(t_H S\) | \(\sim 10^{123}\) | \(\sim 10^{-122}\) |

Page-time half-lives correspond to mean-reversion rates **invisible** on BAO e-fold paths. Scrambling \(\theta\) sits near the MLE numerical floor, not near \(\theta\sim 1\).

---

## 5. Residuals under the repository seed only

### 5.1 The calculation that must not be cheated

With the Sorkin seed fixed by the repository,

$$
\sigma_{\rm res}=\sigma_0\,e^{-\theta\Delta x},
\qquad
\sigma_0\sim 10^{-61}.
$$

For every \(\theta\) in Section 4, and for both \(\Delta x=0.94\) (DESI) and \(\Delta x\simeq 7\) (since recombination),

$$
\sigma_{\rm res}\sim 10^{-61}\ll 10^{-5}\ll 1.5\times 10^{-4}.
$$

**Pure Poisson discreteness is unobservable at DESI and Euclid precision.** This is the first-principles conclusion.

### 5.2 What DESI actually bounds

The limit \(\sigma_X<1.5\times 10^{-4}\) constrains the **effective late-time amplitude** — after any amplification, non-linear growth, or non-Sorkin noise — not the bare \(\sigma_0\).

### 5.3 Amplification required for detectability

Inverting \(\sigma_{\rm res}=A_0 e^{-\theta\Delta x}\),

$$
A_0=\mathrm{target}\times e^{+\theta\Delta x}.
$$

To reach Euclid’s \(\sim 10^{-5}\) on the DESI path with \(\theta\sim 10^{-3}\) already needs \(A_0\sim 10^{-5}\), i.e.

$$
\frac{A_0}{\sigma_0}\sim 10^{56}.
$$

Hubble-scale damping on a recombination-length path still requires \(A_0/\sigma_0\sim 10^{56}\)–\(10^{58}\) depending on \(\theta\). Any **observable** residual under emergent SDiff therefore demands physics beyond bare Sorkin Poisson noise.

---

## 6. Implications for fundamental vs emergent SDiff

| Question | Fundamental SDiff | Emergent SDiff |
|----------|-------------------|----------------|
| Local \(V(x)g_{\mu\nu}\) stresses | Cancelled algebraically | Cancelled only effectively |
| Bare \(\sigma_0\sim 10^{-61}\) | Irrelevant to geometry | Present but unobservable |
| Finite desqueezing / relaxation time | No local bath coupling | Yes: \(t_{1/2}=\ln 2/(\theta H)\) |
| DESI null | Fully expected | Fully expected unless amplified noise exists |
| Euclid detection \(0<\sigma_X\lesssim 10^{-4}\) | Strongly disfavours pure fundamental protection **if** the signal is real and vacuum-sector | Favours effective residual noise with \(A_0\gtrsim 10^{-5}\) |
| Euclid deeper null | Strengthens fundamental (or emergent with tiny \(A_0\)) | Still allowed |

The open-system analogy **demonstrates** that finite relaxation is dynamically natural once a bath couples to a perturbed vacuum. It does **not** prove that cosmology realises such a bath. It supplies the scaling law any emergent Markovian mechanism must obey.

---

## 7. What Euclid DR1 can and cannot see (emergent scenario)

### 7.1 Cannot see (with repo seed only)

- Bare Sorkin fluctuations \(\sigma_0\sim 10^{-61}\).
- Page-time restoration dynamics (\(\theta\sim 10^{-62}\) and smaller).
- Distinguishing fundamental SDiff from “emergent but Sorkin-only” smoothness: both predict a null at \(\sigma_X\sim 10^{-5}\).

### 7.2 Can see (if nature supplies larger effective \(A_0\))

- A residual in the window \(0<\sigma_X\lesssim 10^{-4}\), as anticipated in the SDiff fundamental-vs-emergent note.
- Path-dependent damping if \(\theta\sim\mathcal{O}(0.1\)–\(1)\) and the effective kick is not negligible: then the BAO covariance kernel shape (not only the amplitude) becomes informative once \(N_{\rm bins}\gtrsim 20\).

### 7.3 Operational forecast

| Scenario | Euclid expectation |
|----------|--------------------|
| Fundamental SDiff | Null at \(\sigma_X\sim 10^{-5}\) and below |
| Emergent + Sorkin only | Null (same as fundamental for practical purposes) |
| Emergent + amplified / non-Sorkin \(A_0\gtrsim 10^{-5}\) | Possible detection; compare to DESI limit \(1.5\times 10^{-4}\) |
| Emergent + \(\theta\sim\mathcal{O}(1)\) + moderate \(A_0\) | Amplitude reduced along \(\Delta x\); joint \(\{\theta,\sigma_X\}\) fit becomes meaningful |

---

## 8. Connection to the broader research programme

This note sits downstream of:

1. the DESI DR2 OU/QNM null result and falsification table (F1–F5);
2. the Principle of Late-Time Vacuum Homogeneity;
3. the unimodular / SDiff geometric argument against local vacuum noise;
4. the exclusion of a globally coherent tachyonic quantum-fluid instability;
5. the open-system desqueezing calibration reported here.

Upstream, it supplies a dynamical language — relaxation time, bath rate, residual amplitude — for discussing whether vacuum smoothness is a hard symmetry or an attractor.

---

## 9. Limitations

1. **Single-mode Markovian bath.** Cosmological “baths” need not be Lindblad; non-Markovian or expansion-driven noise may modify the \(\ln 2/\gamma\) law.
2. **Dictionary, not derivation.** \(\gamma\leftrightarrow\theta H\) is a dimensional and structural identification, not a microscopic quantum-gravity derivation of \(\theta\).
3. **Background degeneracy.** As emphasised in the OU papers, \(\sigma_X\) and smooth \(\{w_0,w_a\}\) remain partially degenerate with few BAO bins.
4. **No claim of detection.** All cosmological statements are upper limits or forecasts.

---

## 10. Conclusions

1. Open-system desqueezing yields a robust law \(t_{1/2}=\ln 2/\gamma\), independent of \(n_{\rm th}\) for \(\lvert\langle a^2\rangle\rvert\).
2. Repository OU dynamics imply \(\Gamma_{\rm phys}(z)=\theta H(z)\); the lab–cosmos map is \(\gamma\leftrightarrow\theta H_0\) at \(z=0\).
3. Cosmic-age relaxation requires \(\theta\sim\mathcal{O}(1)\); the DESI MLE floor \(\theta\sim 10^{-3}\) implies almost no path damping.
4. With the repository seed \(\sigma_0\sim 10^{-61}\), residuals stay unobservable: **Euclid cannot see pure Sorkin noise**.
5. A Euclid detection in \(0<\sigma_X\lesssim 10^{-4}\) would require effective amplitudes far above \(\sigma_0\), and would favour **emergent** residual vacuum noise over pure fundamental silence — but only in that amplified sector.
6. A persistent null is compatible with fundamental SDiff and with emergent SDiff that never amplifies the seed.

The synthesis is therefore both constructive and restrictive: finite desqueezing is natural in open quantum systems; cosmological visibility of that logic is not automatic once the repository amplitudes are respected.

---

## Acknowledgements

The DESI collaboration is thanked for public BAO data products. Numerical simulations use QuTiP. This note consolidates results from the stochastic dark energy and vacuum-smoothness research programme.

---

## References

[1] Morales Souhail, J., “Constraints on Stochastic Dark Energy from DESI DR2: A Null Result for Ornstein-Uhlenbeck Fluctuations,” (2026).  
[2] Morales Souhail, J., “Principle of Late-Time Vacuum Homogeneity,” (2026).  
[3] Morales Souhail, J., “SDiff as Fundamental Symmetry vs Emergent Feature,” (2026).  
[4] Morales Souhail, J., “The Smoothness of the Vacuum as a Geometric Necessity,” (2026).  
[5] Morales Souhail, J., “Quantum Fluid Instabilities and DESI DR2,” (2026).  
[6] DESI Collaboration, “DESI DR2 Results II,” arXiv:2503.14738 (2025).  
[7] Bekenstein, J. D., Phys. Rev. D 7, 2333 (1973).  
[8] Sorkin, R. D., arXiv:gr-qc/0503057 (2005).  
[9] Uhlenbeck, G. E. & Ornstein, L. S., Phys. Rev. 36, 823 (1930).  
[10] Johansson, J. R., Nation, P. D. & Nori, F., Comp. Phys. Comm. 183, 1760 (2012); 184, 1234 (2013).  
[11] Walls, D. F. & Milburn, G. J., *Quantum Optics*, Springer.  
[12] Breuer, H.-P. & Petruccione, F., *The Theory of Open Quantum Systems*, Oxford.

---

## Appendix A — Numerical file map

| Path | Content |
|------|---------|
| `code/desqueezing_relax_time.py` | QuTiP scan \(t_{1/2}(\gamma,n_{\rm th})\) |
| `code/cosmological_mapping_from_repo.py` | First-principles Gyr / Sorkin / Euclid tables |
| `code/termico_vs_squeezed_recomendado.py` | Thermal vs squeezed baselines |
| `results/desqueezing_relax_time/` | Figures and CSV of the half-life law |
| `results/cosmological_mapping_from_repo/` | Repo-anchored mapping tables |
| `docs/MEMORIA_REPO_stochastic-dark-energy-ou.md` | Index of the GitHub science corpus |
| `docs/mapping_tables_from_repo.md` | Compact tables |
| `docs/desqueezing-relaxation-vacuum-fluctuations-note.md` | This note |

## Appendix B — Constants used

$$
\begin{aligned}
H_0&=67.4\ \mathrm{km\,s^{-1}\,Mpc^{-1}},\\
t_{H_0}&\simeq 14.51\ \mathrm{Gyr},\\
\Omega_m&=0.315,\\
\Delta x_{\rm DESI}&=0.94,\\
N&\sim 10^{122},\\
\sigma_0&\sim 10^{-61},\\
\sigma_X^{\rm DESI}&<1.5\times 10^{-4},\\
\sigma_X^{\rm Euclid\ target}&\sim 10^{-5}.
\end{aligned}
$$

## Appendix C — One-line summary for readers in a hurry

**Desqueezing dies as \(e^{-\gamma t}\); cosmology writes \(\gamma=\theta H\); Sorkin’s \(10^{-61}\) never reaches Euclid; only an amplified effective \(A_0\) could.**
