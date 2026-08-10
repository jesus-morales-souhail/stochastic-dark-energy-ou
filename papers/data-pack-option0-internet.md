# Data pack (internet): scale/operator probes for Option 0 and closed tests

Jesús Morales Souhail · github.com/jesus-morales-souhail · July 2026

I compiled this pack so that `scale-operator-experiment-map.md` and `anisotropic-slip-option0.md` have citable external numbers next to the repository numbers. Every row has an arXiv or DOI where I could find one, and I state the definition I am using.

---

## A) This repository (internal)

| Quantity | Value | Source in repo |
|----------|-------|----------------|
| $\sigma_X$ (95% CL upper limit) | $< 2.5\times 10^{-2}$ | `papers/resume.txt`, DESI DR2 BAO residual OU |
| Best-fit behaviour | $\sigma_X \to 0$, $\omega_R \to 0$ | MLE; $\Delta\mathrm{AIC}=+4$ vs ΛCDM for free OU |
| Tachyonic coherent growth | Excluded ($\Delta\chi^2\sim +25$) | `papers/quantum-fluid-instabilities-desi-dr2.md` |
| Motivational seed | $\sigma_0\sim 10^{-61}\sim 1/\sqrt{N}$, $N\sim 10^{122}$ | Bekenstein–Hawking / Sorkin (motivational, not a detection) |
| BAO public data used | DESI DR2 summary stats | arXiv:[2503.14738](https://arxiv.org/abs/2503.14738) |

**Operator in repo:** residual kernel on BAO distances $D_V,D_M,D_H$.
**Not in repo:** $\Phi,\Psi,\eta$, CLASS/CAMB.

---

## B) Gravitational slip (correct operator for the “gap”)

### B1. Maus et al. (DESI DR1 × CMB lensing) — verified primary

| Field | Value |
|-------|--------|
| Paper | Maus, White, Sailer et al., *A joint analysis of 3D clustering and galaxy × CMB-lensing cross-correlations with DESI DR1 galaxies* |
| arXiv | [2505.20656](https://arxiv.org/abs/2505.20656) (v3 Oct 2025); JCAP 11 (2025) 077 |
| Definition | $\gamma = \Phi/\Psi$, $\Phi_\gamma=(\Phi+\Psi)/2=\Psi(1+\gamma)/2$(Appendix E) |
| GR value | $\gamma = 1.0$ |
| Result | $\gamma = 1.17 \pm 0.11$(~1.5σ from GR; authors note possible projection effects) |
| Data | DESI DR1 BGS+LRG full-shape + recon; Planck PR4 + ACT DR6 κ; photometric Legacy LRGs |
| Headline structure | $\sigma_8=0.803\pm 0.017$, $\Omega_m=0.3037\pm 0.0069$, $S_8=0.808\pm 0.017$ |

I do not confuse this with the growth index $f=\Omega_m^\gamma$(GR $\gamma\simeq 0.55$).

### B2. Forecast anisotropic stress $\eta$ (Euclid-like + DESI-like)

| Field | Value |
|-------|--------|
| Paper | Sakr, Zheng, Casas, *Model-independent forecasts for the cosmological anisotropic stress* |
| arXiv | [2501.07477](https://arxiv.org/abs/2501.07477) |
| Definition | Effective anisotropic stress / slip $\eta$ (model-independent observables from clustering + lensing) |
| Forecast (Euclid-like photo ± DESI-like spectro) | $\eta$ constant: **~5%**; $z$-only: **&lt;10%** average; free $z,k$: **at least ~30%** |

Earlier forecast I keep on the list: Amendola et al., arXiv:[1311.4765](https://arxiv.org/abs/1311.4765).

---

## C) Modified gravity full-shape (DESI 2024 / DR1 era)

| Field | Value |
|-------|--------|
| Paper | DESI Collaboration / Ishak et al., *Modified Gravity Constraints from the Full Shape Modeling of Clustering Measurements from DESI 2024* |
| arXiv | [2411.12026](https://arxiv.org/abs/2411.12026) (companion to full-shape cosmology [2411.12022](https://arxiv.org/abs/2411.12022)) |
| Parameters | $\mu_0$ (clustering of massive particles), $\Sigma_0$ (lensing-related) |
| DESI alone (example from 2411.12022 abstract) | $\mu_0 = 0.11^{+0.45}_{-0.54}$ (consistent with GR zero) |
| DESI + CMB + DESY3 | $\mu_0 = 0.04 \pm 0.22$, $\Sigma_0 = 0.044 \pm 0.047$(GR-compatible) |

These are not the same symbol as Maus slip $\gamma=\Phi/\Psi$, but they are the right scale of structure-growth / MG operators (RSD + lensing), unlike pupil diffraction.

---

## D) DESI DR2 BAO (background; this repo’s data class)

| Field | Value |
|-------|--------|
| Paper | DESI Collaboration, *DESI DR2 Results II: BAO and Cosmological Constraints* |
| arXiv | [2503.14738](https://arxiv.org/abs/2503.14738) |
| Sample | &gt;14M galaxies and quasars, 3 years |
| ΛCDM | BAO consistent with smooth distance–redshift; mild tension with CMB parameters (~2.3σ class) |
| $w_0 w_a$ | Preference for evolving DE in combinations with CMB ± SN (significance sample-dependent, often ~3σ class) |

This is the correct scale+operator family for $\sigma_X$ residual tests on distances — already used in this repository.

---

## E) Order-of-magnitude comparison (repo + internet)

Using $\sigma_X \sim 10^{-4}$, $\Omega_{\mathrm{DE}}/\Omega_m\sim 2$, $f=1$:


$$
\lvert \eta-1 \rvert\sim \mathcal{O}\!\left(f\cdot\frac{\Omega_{\mathrm{DE}}}{\Omega_m}\cdot\sigma_X\right)\sim 2\times 10^{-4}.
$$


| Probe | Rough reach on $\lvert \eta-1 \rvert$ or $\lvert \gamma-1 \rvert$ | vs $2\times 10^{-4}$ |
|-------|-----------------------------------------------|------------------------|
| Maus et al. slip $\sigma(\gamma)\sim 0.11$ | $\mathcal{O}(0.1)$ | ~500× larger than prediction |
| Sakr et al. forecast (const $\eta$) | $\sim 5\% = 0.05$ | ~250× larger |
| Sakr free $z,k$ | $\gtrsim 30\%$ | still ≫ |

**Conclusion (unchanged):** the external literature confirms Option 0. The correct operator exists and is measured or forecast; the amplitude still requires amplification before a detection story makes sense.

---

## F) Closed wrong scale experiment (local)

| Item | Value |
|------|--------|
| Script | `scripts/car_drone_pupil_newton_einstein.py` |
| $v$ | 120 km/h $\Rightarrow \beta\sim 1.1\times 10^{-7}$, $\gamma-1\sim 6\times 10^{-15}$ |
| Pupil | $D=1 \mathrm{mm}$, $\theta_{\mathrm{Airy}}\sim 6.7\times 10^{-4}$ rad |
| Verdict | **Wrong scale + wrong operator** for cosmic expansion |

---

## G) Reading / download checklist

- [x] Maus et al. 2505.20656 — slip definition + number
- [x] Sakr et al. 2501.07477 — $\eta$ forecasts
- [x] DESI full-shape MG 2411.12026 / 2411.12022 — $\mu_0,\Sigma_0$
- [x] DESI DR2 BAO 2503.14738 — background distances
- [ ] Full PDF of any “Dark Energy After DESI DR2” review if I want a single secondary survey (optional)
- [ ] Plaza et al. unimodular + DESI (pipeline craft only)

---

## H) Links (quick)

| Topic | URL |
|-------|-----|
| Slip DESI DR1 | https://arxiv.org/abs/2505.20656 |
| $\eta$ forecasts | https://arxiv.org/abs/2501.07477 |
| DESI MG full-shape | https://arxiv.org/abs/2411.12026 |
| DESI FS cosmology | https://arxiv.org/abs/2411.12022 |
| DESI DR2 BAO | https://arxiv.org/abs/2503.14738 |
| This repo | https://github.com/jesus-morales-souhail/stochastic-dark-energy-ou |

---

Data pack only. No Boltzmann run. Definitions are mandatory when I cite a number.

