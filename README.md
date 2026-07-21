# stochastic-dark-energy-ou

**Author:** Jesús Morales Souhail  
**ORCID:** [0009-0000-7637-1818](https://orcid.org/0009-0000-7637-1818)  
**Contact:** jmskjym@gmail.com  
**Status:** Independent research (preprint drafts — **not** peer reviewed) · July 2026  

Independent analysis of **stochastic fluctuations in the dark-energy sector** with **public DESI DR2 BAO** products: code, numerical results, and a **coherent set of notes** on what the data constrain, what geometry may protect, and which speculative paths are closed.

---

## Narrative arc (how the papers fit together)

The project asks one empirical question and two structural ones:

1. **Empirical:** Do DESI DR2 BAO residual statistics require a stationary stochastic (OU/QNM) component on top of a smooth background?  
2. **Structural:** Can volume-preserving diffeomorphisms (SDiff / unimodular gravity) explain **late-time vacuum smoothness**?  
3. **Method:** What is the **correct scale and operator** for each claim — and what is numerology?

### Act I — Data (BAO residuals)

| Result | Detail |
|--------|--------|
| CPL background (BAO-only) | \(w_0 \approx -0.99\), \(w_a \approx -0.02\) (near ΛCDM) |
| OU / QNM residuals | Amplitude driven to the numerical floor when the background is free |
| Working upper limit | **\(\sigma_X < 1.5\times 10^{-4}\) (95% CL)** under the phenomenological OU kernel |
| Tachyonic quantum fluid | Coherent growing mode **excluded** (\(\Delta\ln\mathcal{L}\approx -11.35\) vs ΛCDM) |

**Primary papers:** `stochastic-dark-energy-desi-dr2.md`, `quantum-fluid-instabilities-desi-dr2.md`, `resume.txt`, `sensitivity_kernel_table.md`.

### Act II — Geometry of smoothness (SDiff)

Unimodular / SDiff structure projects out local vacuum stress of the form \(T_{\mu\nu}\propto g_{\mu\nu}\). That is a candidate explanation of **why residuals are small**, not a detection of a new particle.

**Primary papers:** `principle-of-vacuum-smoothness.md`, `smoothness-of-the-vacuum-unimodular.md`, `unimodular-gravity-vacuum-smoothness.md`, `sdiff-fundamental-vs-emergent.md`.

### Act III — Open systems and amplification

A pure Sorkin–Bekenstein seed \(\sigma_0\sim 10^{-61}\) sits far below BAO reach. Open-system **desqueezing** supplies a dynamical language for finite relaxation (\(t_{1/2}=\ln 2/\gamma\)) and makes the **amplification** problem explicit: without a physical map from seed to \(A_0\gtrsim 10^{-5}\), Euclid-scale residual detections are not automatic.

**Primary notes:** `notes/desqueezing-relaxation-vacuum-fluctuations-note.md`, `fundamental-vs-emergent-vacuum-relaxation.md`, `euclid-protocol-vacuum-relaxation.md`.

### Act IV — Option 0: the anisotropic “gap” (no homemade Boltzmann)

SDiff does **not** cancel shear. Gravitational slip \(\gamma=\Phi/\Psi\) (Maus et al., arXiv:2505.20656: \(\gamma=1.17\pm 0.11\), GR\(=1\)) is the right **kind** of operator. Order-of-magnitude: even \(f=1\) and \(\sigma_X\sim 10^{-4}\) give \(|\eta-1|\sim 10^{-4}\), far below current \(\mathcal{O}(0.1)\) sensitivity. **The gap inherits the amplification problem; it is not a shortcut.**

**Primary papers:** `anisotropic-slip-option0.md`, `data-pack-option0-internet.md`.

### Act V — Closed wrong paths (scale / operator hygiene)

These notes keep the corpus honest and prevent mixing lab metaphors into DESI claims:

| Closed claim | Why closed | Document |
|--------------|------------|----------|
| Highway pupil / “tesseract diffraction” tests expansion | Wrong **scale** and **operator** | `car-drone-pupil-newton-einstein.md` |
| Phase tricks force single-photon sub-Airy certainty | Born + band-limit | `no-go-superoscillation-tesseract.md`, `self-shielding-triple-barrier.md` |
| Tesseract/\(B_4\) as channel amplifier | 32/384 symplectic; commutant dim 1 in \(\mathfrak{sp}(4,\mathbb{R})\) — symmetry *kills* DOF | `no-go-superoscillation-tesseract.md` §3.2–4.2, `scripts/b4_symplectic_count.py` |
| Lab wavefront \(T\) = BAO OU vacuum | Different objects | `wavefront-shaping-vs-ou-vacuum.md` |
| \(\ln 4\equiv\omega_R\) or \(B_4\) as undeclared physics | Pattern hygiene | `pattern-undeclared-physical-power.md` |
| \(f_{\mathrm{cosmo}}\equiv 1/\theta\) as optical identity | Analogy only | `optics-ou-analogies-and-limits.md` |

**Master map:** `scale-operator-experiment-map.md`  
**Shielding without cynicism:** `self-shielding-vs-untestability.md`

### What this programme **does not** claim

- Detection of unamplified Planck-scale vacuum noise.  
- That Euclid will “easily” see Sorkin seeds or slip at \(10^{-4}\) without amplification.  
- That optical tesseracts, SLMs, or wavefront shaping of the cosmic vacuum are part of the DESI result.  
- Peer-reviewed status.

### What it **does** claim

- Under the stated OU/QNM residual kernel and public BAO summary statistics, the data **prefer smooth** evolution; \(\sigma_X\) is bounded.  
- A coherent tachyonic growth model with rank-1 covariance is **disfavoured**.  
- SDiff is a **structural candidate** for isotropic smoothness, with a **shear gap** that still needs amplitude honesty.  
- Method: **correct scale + correct operator**, or no claim.

---

## Reading order (cohesive path)

1. `papers/resume.txt` — numbers only  
2. `papers/stochastic-dark-energy-desi-dr2.md` — main BAO analysis  
3. `papers/quantum-fluid-instabilities-desi-dr2.md` — model kill  
4. `papers/principle-of-vacuum-smoothness.md` + `sdiff-fundamental-vs-emergent.md` — geometry  
5. `notes/desqueezing-relaxation-vacuum-fluctuations-note.md` — amplification  
6. `papers/anisotropic-slip-option0.md` + `data-pack-option0-internet.md` — slip / Option 0  
7. `papers/scale-operator-experiment-map.md` — what is closed vs next  
8. Boundary notes (V) as needed when analogies appear  

---

## Repository layout

```
papers/       Scientific notes (English) — full narrative above
scripts/      Analysis code
  desqueezing/   QuTiP Lindblad desqueezing + cosmological mapping
  gpe/           Gross–Pitaevskii / Bogoliubov numerics
figures/      Figures
results/      Numerical outputs (BAO fits, forecasts, desqueezing)
notes/        Desqueezing synthesis and mapping tables
```

Detailed file index: **`STRUCTURE.md`**.

### Separated project (not this repo)

Lab **transmission-matrix** wavefront shaping lives here (do not mix into DESI papers):

`~/Proyectos/04_Optica_medios_complejos/`

---

## Quick start

```bash
git clone https://github.com/jesus-morales-souhail/stochastic-dark-energy-ou.git
cd stochastic-dark-energy-ou
python -m venv venv && source venv/bin/activate   # optional
pip install -r requirements.txt

# BAO OU / QNM residual kernel (DESI DR2 summary statistics)
python scripts/ou_bao_likelihood.py

# CPL background (+ optional nested σ, θ extension)
python scripts/eos_efectiva.py

# Desqueezing relaxation scan (requires qutip)
python scripts/desqueezing/desqueezing_relax_time.py

# Cosmological mapping tables
python scripts/desqueezing/cosmological_mapping.py

# Euclid forecast MCMC (mock likelihood; requires emcee)
python scripts/euclid_mock_mcmc.py
python scripts/euclid_joint_bao_sne_mcmc.py

# Pedagogy / hygiene demos (not cosmology likelihoods)
python scripts/superoscillation_energy_cost_demo.py
python scripts/car_drone_pupil_newton_einstein.py
```

**Dependencies:** `numpy`, `scipy`, `matplotlib`, `astropy`, `healpy`, `qutip` (see `requirements.txt`).

---

## Data sources (public)

- DESI DR2 BAO: [data.desi.lbl.gov](https://data.desi.lbl.gov/public/) · arXiv:[2503.14738](https://arxiv.org/abs/2503.14738)  
- Slip / structure: Maus et al. arXiv:[2505.20656](https://arxiv.org/abs/2505.20656)  
- \(\eta\) forecasts: Sakr et al. arXiv:[2501.07477](https://arxiv.org/abs/2501.07477)  
- Pantheon+ SN Ia (cross-correlation pipeline only)  

Compiled links: `papers/data-pack-option0-internet.md`.

---

## Reproducibility notes

- Core BAO OU/QNM scripts use **published BAO summary statistics** (redshift bins, \(\alpha\), \(\sigma\)) and do not require bulk catalog access.  
- Cross-correlation scripts may request large public FITS products; see script headers.  
- Large local archives and obsolete drafts (e.g. \(\ln 4\) narrative scripts) stay **off** the public tree (`.gitignore` / `local_archive/`).

---

## Citation

```text
Morales Souhail, J. (2026). stochastic-dark-energy-ou (GitHub).
https://github.com/jesus-morales-souhail/stochastic-dark-energy-ou
ORCID: 0009-0000-7637-1818
```

---

## License

Code: MIT (unless otherwise noted in individual files).  
Text: CC BY 4.0 for author-written notes.

---

## Contact

**Jesús Morales Souhail** · jmskjym@gmail.com · [ORCID](https://orcid.org/0009-0000-7637-1818)
