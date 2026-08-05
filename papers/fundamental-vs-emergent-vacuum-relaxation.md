# Fundamental versus Emergent Vacuum Smoothness: Path-Integrated Relaxation, Amplification Bounds, and Euclid Discriminators

**Author:** Jesús Morales Souhail
**Date:** July 2026
**ORCID:** [0009-0000-7637-1818](https://orcid.org/0009-0000-7637-1818)
**Repository:** https://github.com/jesus-morales-souhail/stochastic-dark-energy-ou
**Status:** Technical note
**Related:** `sdiff-fundamental-vs-emergent.md`, `desqueezing-relaxation-vacuum-fluctuations-note.md`, `stochastic-dark-energy-desi-dr2.md`

---

## Abstract

My DESI DR2 BAO analyses in this repository drive stationary stochastic fluctuations of dark energy to a null working limit $\sigma_X < 1.5\times 10^{-4}$ (95% CL). The open question is whether late-time vacuum smoothness is an exact geometric law (fundamental volume-preserving diffeomorphisms, SDiff) or an attractor with finite relaxation time (emergent SDiff). Combining the Ornstein–Uhlenbeck continuity structure $\Gamma_{\rm phys}(z)=\theta H(z)$ with the open-system desqueezing law $t_{1/2}=\ln 2/\gamma$, I derive path-integrated residuals $\sigma_{\rm res}=A_0 e^{-\theta\Delta x}$, estimate the kick amplitude needed for Euclid-scale detectability, and show that a pure Sorkin–Bekenstein seed $\sigma_0\sim 10^{-61}$ remains unobservable for cosmologically relevant $\theta$. I partition theory space into regions F and E0–E3 with corresponding Euclid and lag-correlation tests.

---

## 1. The question

The empirical situation is not “how large is $\sigma_X$?” — the data push it down — but:

> Is the observed late-time vacuum homogeneity a **fundamental** property of the gauge structure of spacetime, or the **end state of a relaxation process** with a finite characteristic time?

If fundamental, local vacuum stresses of the form $T_{\mu\nu}=V(x)g_{\mu\nu}$ are algebraically projected out (unimodular / SDiff). If emergent, a perturbation of the vacuum should decay with a finite $\tau$, and residual fluctuations may survive depending on the kick amplitude and the mean-reversion rate.

---

## 2. Minimal dynamical structure

### 2.1 OU sector (this repository)

I take fractional fluctuations $X\equiv\delta\Omega_\Lambda$ in $x=\ln a$:


dX=-\theta X dx+\sigma dW_x.


With $dx=H dt$,


\Gamma_{\rm phys}(z)=\theta H(z).


For the homogeneous mean (noise off), a frozen kick of amplitude $A_0$ damps as


\boxed{\sigma_{\rm res}(x)=A_0 e^{-\theta\Delta x}, \qquad \Delta x=\ln\frac{1+z_{\rm past}}{1+z_{\rm now}}.}


### 2.2 Laboratory desqueezing

Open-system simulations (QuTiP Lindblad, $N=150$, $r=1.5$) yield


\boxed{t_{1/2}\!\left(\lvert\langle a^2\rangle\rvert\right)=\frac{\ln 2}{\gamma},}


independent of thermal occupation $n_{\rm th}$ for this metric. The present-epoch dictionary is


\gamma \longleftrightarrow \theta H_0 \qquad\Rightarrow\qquad t_{1/2}=\frac{\ln 2}{\theta H_0}.


Cosmic-age half-life requires $\theta\sim\mathcal{O}(1)$. The DESI MLE numerical floor $\theta\sim 10^{-3}$ implies essentially no damping on the DESI path $\Delta x\simeq 0.94$.

---

## 3. Quantitative results (fiducial $H_0=67.4$, $\Omega_m=0.315$)

| Quantity | Value |
|----------|------:|
| Age of Universe $t_0$ | $13.80 \mathrm{Gyr}$ |
| Hubble time $t_{H_0}$ | $14.51 \mathrm{Gyr}$ |
| DESI path $\Delta x$ | $0.94$ |
| Post-recombination $\Delta x=\ln(1+z_{\rm rec})$ | $\simeq 6.99$ |
| DESI working limit | $\sigma_X<1.5\times 10^{-4}$ (95% CL) |
| Euclid target scale | $\sim 10^{-5}$ |
| Sorkin seed | $\sigma_0\sim 10^{-61}$ |

### 3.1 Sorkin-only corollary

For any $\theta\ge 0$ and any path with $\Delta x=\mathcal{O}(1\text{–}10)$,


\sigma_{\rm res}\le\sigma_0\sim 10^{-61}\ll 10^{-5}.


**A pure Poisson discreteness seed is invisible to DESI and to Euclid.**
DESI therefore constrains the **effective** late-time amplitude after amplification or non-Sorkin coupling — not $\sigma_0$ itself.

### 3.2 Amplification to reach Euclid on the DESI path

Minimum kick $A_0^{\min}=\sigma_{\rm target} e^{+\theta\Delta x}$:

| $\theta$ | $A_0^{\min}$ for $\sigma=10^{-5}$ | $A_0^{\min}/\sigma_0$ |
|-----------:|--------------------------------------:|------------------------:|
| $10^{-3}$ | $\sim 1.0\times 10^{-5}$ | $\sim 10^{56}$ |
| $0.1$ | $\sim 1.1\times 10^{-5}$ | $\sim 10^{56}$ |
| $\ln 2\simeq 0.69$ | $\sim 1.9\times 10^{-5}$ | $\sim 10^{56}$ |
| $1.2$ | $\sim 3.1\times 10^{-5}$ | $\sim 10^{56}$ |
| $5$ | $\sim 1.1\times 10^{-3}$ | $\sim 10^{58}$ |

Any **observable** residual under an emergent picture requires physics beyond bare Sorkin noise (or an equivalent effective $A_0$).

### 3.3 Physical half-life (Mapping A)


t_{1/2}=\frac{\ln 2}{\theta H_0}.


| $\theta$ | $t_{1/2}$ |
|-----------:|------------:|
| $10^{-3}$ | $\sim 10^{4} \mathrm{Gyr}$ |
| $\ln 2$ | $\sim t_{H_0}$ |
| $\sim 0.73$ | $\sim t_0$ |
| $1.2$ | $\sim 8.4 \mathrm{Gyr}$ |

---

## 4. Theory-space partition

| Region | Definition | DESI status | Euclid expectation |
|--------|------------|-------------|--------------------|
| **F** | Fundamental SDiff: geometric $A_0^{\rm geo}=0$ for local $Vg_{\mu\nu}$ | Compatible | Deeper null |
| **E0** | Emergent + Sorkin-only seed | Compatible | Null (practically identical to F) |
| **E1** | Emergent + $A_0\in[10^{-5},1.5\times 10^{-4}]$, small $\theta$ | Compatible if residual $<1.5\times 10^{-4}$ | **Possible detection** |
| **E2** | Emergent + $\theta\sim\mathcal{O}(1)$ + moderate $A_0$ | Compatible if damped below DESI | Residual near/below $10^{-5}$ after long paths |
| **E3** | Small $\theta$ + $A_0>1.5\times 10^{-4}$ | **Tension with DESI null** | Already disfavored |

Euclid alone cannot separate **F** from **E0**. Separation requires either a residual detection (favouring E1-type physics) or an independent prior on $A_0$.

---

## 5. Observational tests

1. **Euclid residual** $0<\sigma_X\lesssim 10^{-4}$ in a joint fit with free $\{w_0,w_a\}$: disfavours pure **F** (if systematics are controlled).
2. **Euclid null** $\sigma_X\ll 10^{-5}$ with $\gtrsim 20$ bins: forces $A_0 e^{-\theta\Delta x}$ small; compatible with **F** and **E0**.
3. **Lag structure** in whitened BAO residuals matching an OU kernel: supports dynamical mean-reversion rather than pure algebraic projection.
4. **Open-system map:** a Markovian bath mapped to cosmology should obey $t_{1/2}=\ln 2/\gamma$ with $\gamma\leftrightarrow\theta H$.

I do not derive SDiff from quantum gravity in this note, nor do I identify a microphysical amplifier from $10^{-61}$ to $10^{-5}$. I organise the fundamental/emergent alternatives into quantitative regions and tests, subject to a full joint analysis with free background parameters and the survey covariance.

---

## 6. Numerical products

| File | Content |
|------|---------|
| `results/sdiff_discrimination/sdiff_discrimination_grid.csv` | Grid over $\theta$ and $A_0$ scenarios |
| `results/sdiff_discrimination/summary.txt` | Numerical summary |
| `figures/sdiff_residual_vs_theta.png` | Residual vs $\theta$ |
| `figures/sdiff_A0_min_for_detection.png` | Minimum $A_0$ for detection |
| `figures/sdiff_t_half_vs_theta.png` | Half-life vs $\theta$ |

```bash
python scripts/desqueezing/sdiff_discrimination.py
```

---

## 7. Conclusion

The DESI null raises the status of vacuum smoothness as a quantitative question:

- **Fundamental SDiff** predicts geometric silence at BAO scales.
- **Emergent SDiff** allows finite relaxation, with open-system scaling $t_{1/2}=\ln 2/\gamma$ and residual $A_0 e^{-\theta\Delta x}$.
- A pure Sorkin seed remains unobservable at DESI/Euclid precision.
- A residual in $10^{-5}$–$10^{-4}$ would indicate effective physics well above that seed.

I specify a Euclid-facing fit with free $\{w_0,w_a,\theta,\sigma_X\}$ in `euclid-protocol-vacuum-relaxation.md`.

---

## References

[1] Morales Souhail, J., “Constraints on Stochastic Dark Energy from DESI DR2,” this repository (`papers/stochastic-dark-energy-desi-dr2.md`).
[2] Morales Souhail, J., “SDiff as Fundamental Symmetry vs Emergent Feature,” this repository.
[3] Morales Souhail, J., “Desqueezing Relaxation Time…,” this repository (`notes/`).
[4] DESI Collaboration, arXiv:2503.14738 (2025).
[5] Sorkin, R. D., arXiv:gr-qc/0503057 (2005).
[6] Bekenstein, J. D., Phys. Rev. D 7, 2333 (1973).
[7] Uhlenbeck, G. E. & Ornstein, L. S., Phys. Rev. 36, 823 (1930).
