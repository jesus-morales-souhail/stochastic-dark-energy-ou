# Scale / operator map: closed results and correct next experiments

**Author:** Jesús Morales Souhail 
**Date:** July 2026 
**Status:** Experimental logic map (closes the car–drone–pupil result; points to correct probes) 
**Related:** 
`papers/car-drone-pupil-newton-einstein.md`, 
`scripts/car_drone_pupil_newton_einstein.py`, 
`papers/anisotropic-slip-option0.md`, 
`papers/self-shielding-vs-untestability.md`, 
`papers/no-go-superoscillation-tesseract.md`, 
`papers/resume.txt`

---

## 0. Rule

A measurement is only informative if **both** are right:

1. **Scale** — length, time, energy, redshift range where the effect is \(\mathcal{O}(1)\) or at least above noise. 
2. **Operator** — which field / combination of observables couples to the mechanism (e.g. isotropic residual on \(D_V\) vs slip \(\Phi/\Psi\) vs free-photon travel time).

Wrong scale **or** wrong operator ⇒ null result that **does not** falsify the deep premise; it only falsifies “this effect is large here.”

---

## 1. Closed result: car + drone @ 120 km/h + pupil diffraction

| Item | Value |
|------|--------|
| Setup | Car and drone co-moving at \(120 \mathrm{km/h}\); \(D=1 \mathrm{mm}\); \(\lambda=550 \mathrm{nm}\) |
| Newton / Galileo | Relative \(v=0\) ⇒ **identical** Airy pattern inside and for the drone |
| Einstein (SR) | Relative \(\beta=0\) ⇒ same; even vs road, \(\beta/\theta_{\mathrm{Airy}}\sim 10^{-4}\) |
| Hubble stretch | \(\Delta L\sim 10^{-16} \mathrm{m}\) per metre per minute |
| Tesseract / 4D optical claim | **Not present** in any consistent calculation |

### Verdict (closed)

| Question asked | Verdict |
|----------------|---------|
| Does expansion deform the Airy pattern in a pupil at highway speed? | **No** (wrong **scale**) |
| Does co-moving motion create a “hypercube projection” of the diffraction pattern? | **No** (wrong **operator** / non-derived geometry) |
| What was actually tested? | Local EM diffraction + inertial SR at \(\beta\sim 10^{-7}\) |

**Status: CLOSED.** Documented in `car-drone-pupil-newton-einstein.md` and the runnable script.

**One-line moral:** local EM + inertial physics **shield** the experiment from cosmology — same *kind* of “wrong question” as asking BAO for unamplified Planck-scale isotropic noise.

---

## 2. Experiment matrix: other scales, correct operators

### 2.1 Cosmology / vacuum (this repository)

| Scientific question | Wrong setup | Correct scale | Correct operator / data | Status in this project |
|---------------------|-------------|---------------|-------------------------|------------------------|
| Unamplified isotropic vacuum noise (Sorkin \(\sigma_0\sim 10^{-61}\)) | Pupil, lab interferometer, single BAO bin “wiggle hunt” | Cosmological path length + enormous amplification if any | Residual kernel on distances; still needs \(A_0\) | **Not detectable**; null on effective \(\sigma_X\) |
| Effective late-time residual amplitude | Claiming detection without Occam | DESI BAO multi-bin, public covariances | \(\sigma_X\) on BAO residual covariance (OU/QNM) | **Done:** \(\sigma_X < 1.5\times 10^{-4}\) (95% CL) |
| Coherent tachyonic growth of DE fluid | Stationary OU proxy for growing mode | Same BAO, correct rank-1 covariance | Growing-mode \(\sigma_X(t)\) | **Done:** model excluded |
| SDiff “gap” (anisotropic stress) | Isotropic BAO residual only | RSD + weak lensing / CMB lensing | Slip \(\gamma=\Phi/\Psi\) or \(\mu,\Sigma\) | **Option 0:** Maus+ \(\gamma=1.17\pm 0.11\) (GR\(=1\)); \(|\eta-1|\sim 10^{-4}\) still below reach without amplification |
| Background \(w(z)\) | Stochastic kernel alone | DESI+CMB+SN | CPL \(w_0,w_a\) | Independent of SDiff smoothness |
| Physical amplification \(10^{-61}\to 10^{-4}\) | Hand-waving \(10^{56}\) factor | Lab open systems + derived map | Desqueezing \(t_{1/2}=\ln 2/\gamma\), open-system models | **Open theory problem** |
| Spacetime discreteness (premise, not shield) | Gravity-only channels if shield applies | Astrophysical baselines, high-\(E\) photons | Energy-dependent time of flight (GRB / γ-ray), etc. | Literature-driven; not BAO residual |

### 2.2 Local / lab (PC-scale, no ad hoc free-parameter tuning)

| Scientific question | Wrong setup | Correct scale | Correct operator | Tool on this PC |
|---------------------|-------------|---------------|------------------|-----------------|
| “Free work” from ideal current / levitation | Ignore Lenz / 3rd law | Force, power, temperature of SC | Mechanical reaction or magnetic drag \(P_{\mathrm{load}}=-P_{\mathrm{drag}}\) | Analytical + optional FEM later |
| Energy balance under forced work | “Efficiency 100% useful” claim | Molecular / nm–μm MD box | Work vs heat to thermostat (1st law check) | **GROMACS 2026.2 + CUDA** (already installed) |
| Superoscillation “beats Airy for free” | “Perfect phase ⇒ photon always in core” | Band-limited 1D/2D field | \(E_{\mathrm{core}}/E_{\mathrm{total}}\) under Born | `scripts/superoscillation_energy_cost_demo.py` (**done**) |
| Tesseract optical engine | Assign \(B_4\) power without Maxwell map | — | — | **No-go** (`no-go-superoscillation-tesseract.md`) |

### 2.3 Optics / relativity (extensions of the car experiment)

| Scientific question | Wrong setup | Correct scale | Correct operator | Experiment |
|---------------------|-------------|---------------|------------------|------------|
| Inertial diffraction | Expect Hubble in the eye | \(\beta\sim 10^{-7}\) highway | Airy vs aberration | **Closed** (car–drone) |
| Relativistic aberration of diffraction | 120 km/h | \(v\sim 10^2 \mathrm{km/s}\) or lab \(\beta\) with precision metrology | Angular pattern + Doppler | Not required for cosmology programme |
| Equivalence principle (accelerating car) | Constant \(v\) | Acceleration \(a\), light-crossing time of apparatus | Effective \(g\) tilt of rays / clocks | Optional add-on (not expansion) |
| Cosmic expansion with light | Pupil in car | Gpc path, \(z\sim 0.1\text{–}2\) | Redshift, BAO \(D_M,D_H\), SN \(\mu(z)\) | DESI / Euclid (this repo) |

---

## 3. Ordered programme (what to do next)

### Already closed
1. Car–drone–pupil @ 120 km/h → **wrong scale / wrong operator** for expansion and for “4D diffraction.” 
2. Superoscillation energy tax demo → Born peaje real. 
3. Tesseract optical claim → no-go (undeclared power). 
4. DESI BAO OU residual → \(\sigma_X\) upper limit; no detection. 
5. Tachyonic coherent growth → excluded with correct covariance. 
6. Option 0 slip → Maus et al. definition verified; amplitude still kills easy detection.

### Next experiments (correct scale + operator)

| Priority | Experiment | Scale | Operator | Why |
|----------|------------|-------|----------|-----|
| **P0** | Tighten \(\sigma_X\) with full DESI cov / future DR | BAO multi-bin | Residual kernel | Publishable null science |
| **P0** | More **concrete models** with sharp BAO predictions | Same | Model-specific cov | Already works (GPE case) |
| **P1** | Amplification physics (open systems) | Lab + theory map | \(\gamma\leftrightarrow\theta H_0\), \(A_0\) | Only way unamplified seeds matter |
| **P1** | Option 0 reading complete (slip papers) | DESI RSD+lensing | \(\gamma=\Phi/\Psi\), \(\mu,\Sigma\) | Cite, don’t invent |
| **P2** | hi_class / MGCAMB **only if** \(A_0\) + \(f\) make \(|\eta-1|\) reachable | Perturbations | Anisotropic stress | Not before amplitude honesty |
| **P2** | GROMACS “work vs heat” micro-demo | nm–μm, ns–μs | 1st law balance | Local “shield + extraction” pedagogy |
| **Avoid** | Pupil / car / tesseract as cosmology tests | Human scale | Diffraction | Closed wrong |

---

## 4. Decision rule (one sentence)

> If the predicted signal at the **correct operator** is still \(\ll\) instrument noise after all honest amplitudes, the next step is **theory of amplification or a new operator** — not a more clever geometry at the wrong scale.

---

## 5. Commands

```bash
# Closed local experiment
python scripts/car_drone_pupil_newton_einstein.py

# Born / energy tax (optics)
python scripts/superoscillation_energy_cost_demo.py
```

---

*End of map. Car–drone–pupil: CLOSED as wrong scale/operator. Next work: correct scales and operators only.*
