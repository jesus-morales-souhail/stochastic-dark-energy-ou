# stochastic-dark-energy-ou

**Author:** Jesús Morales Souhail  
**ORCID:** [0009-0000-7637-1818](https://orcid.org/0009-0000-7637-1818)  
**Contact:** jmskjym@gmail.com  
**Status:** Independent research (preprint drafts — not peer reviewed) · July 2026  

Independent analysis of **stochastic fluctuations in the dark-energy sector** using **public DESI DR2 BAO** products. The repository contains analysis code, numerical results, and draft scientific notes on:

1. Stationary Ornstein–Uhlenbeck (OU) and quasi-normal-mode (QNM) covariance models  
2. Non-stationary rank-1 covariance from a tachyonic quantum-fluid (Bogoliubov) instability  
3. Geometric protection of vacuum smoothness (volume-preserving diffeomorphisms / unimodular gravity)  
4. Open-system **desqueezing** as a dynamical analogue of finite vacuum relaxation  

Cosmological constraints use **public DESI BAO products** (and Pantheon+ where noted). Core BAO likelihood scripts run from published summary statistics; no proprietary catalogs are required.

---

## Main findings (summary)

| Result | Detail |
|--------|--------|
| BAO-only CPL background | **w₀ ≈ −0.99, w_a ≈ −0.02** (near flat ΛCDM; `results/eos_cpl_desi_dr2/`) |
| OU / QNM residuals | Stochastic amplitude at the numerical floor when the background is free |
| Working upper limit | **σ_X < 1.5 × 10⁻⁴ (95% CL)** under the phenomenological OU kernel |
| Tachyonic quantum fluid | Globally coherent growing mode **excluded** (\(\Delta\ln\mathcal{L}\approx -11.35\) vs \(\Lambda\)CDM) |
| Vacuum smoothness | Interpreted via SDiff / unimodular structure; Euclid DR1 as a future discriminator |
| Desqueezing (open systems) | \(t_{1/2}(\lvert\langle a^2\rangle\rvert)=\ln 2/\gamma\); Sorkin-scale seeds \(\sim 10^{-61}\) lie far below BAO sensitivity |

Full narrative and caveats: see `papers/` and `papers/resume.txt`.

---

## Repository layout

```
papers/       English scientific notes (primary reading)
scripts/      Runnable analysis code
  desqueezing/   QuTiP Lindblad desqueezing + cosmological mapping
  gpe/           Gross–Pitaevskii / Bogoliubov numerics
figures/      PNG figures referenced in the notes
results/      Numerical tables and figure data for desqueezing
notes/        Desqueezing synthesis note and mapping tables
```

---

## Papers (`papers/`)

| File | Description |
|------|-------------|
| [`stochastic-dark-energy-desi-dr2.md`](papers/stochastic-dark-energy-desi-dr2.md) | Main OU/QNM analysis vs DESI DR2 |
| [`quantum-fluid-instabilities-desi-dr2.md`](papers/quantum-fluid-instabilities-desi-dr2.md) | Tachyonic quantum fluid + rank-1 covariance (full version) |
| [`principle-of-vacuum-smoothness.md`](papers/principle-of-vacuum-smoothness.md) | Principle of Late-Time Vacuum Homogeneity |
| [`smoothness-of-the-vacuum-unimodular.md`](papers/smoothness-of-the-vacuum-unimodular.md) | Unimodular gravity and vacuum smoothness |
| [`unimodular-gravity-vacuum-smoothness.md`](papers/unimodular-gravity-vacuum-smoothness.md) | Extended geometric synthesis |
| [`sdiff-fundamental-vs-emergent.md`](papers/sdiff-fundamental-vs-emergent.md) | Fundamental vs emergent SDiff; Euclid tests |
| [`sensitivity_kernel_table.md`](papers/sensitivity_kernel_table.md) | BAO sensitivity kernel \(S(z)\) |
| [`resume.txt`](papers/resume.txt) | Compact numerical summary |

Supporting note:

| File | Description |
|------|-------------|
| [`notes/desqueezing-relaxation-vacuum-fluctuations-note.md`](notes/desqueezing-relaxation-vacuum-fluctuations-note.md) | Open-system half-life and map to \(\theta H(z)\) |
| [`papers/fundamental-vs-emergent-vacuum-relaxation.md`](papers/fundamental-vs-emergent-vacuum-relaxation.md) | Path-integrated residuals; SDiff regions F/E0–E3 |
| [`papers/euclid-protocol-vacuum-relaxation.md`](papers/euclid-protocol-vacuum-relaxation.md) | Euclid BAO fit: parameters, priors, and region assignment |

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
```

**Dependencies:** `numpy`, `scipy`, `matplotlib`, `astropy`, `healpy`, `qutip` (see `requirements.txt`).

---

## Data sources (public)

- DESI DR2 BAO: [data.desi.lbl.gov](https://data.desi.lbl.gov/public/) · arXiv:2503.14738  
- Pantheon+ SN Ia (cross-correlation pipeline only)

---

## Reproducibility notes

- Core BAO OU/QNM scripts use **published BAO summary statistics** (redshift bins, \(\alpha\), \(\sigma\)) and do not require bulk catalog access.
- Cross-correlation scripts may request large public FITS products; see script headers.
- Large local archives and draft Word histories are kept off-Git (see `.gitignore`).

---

## Citation

If you use this code or notes, please cite the repository and ORCID:

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

For questions about the analysis or access to public DESI products used here:  
**Jesús Morales Souhail** · jmskjym@gmail.com · [ORCID](https://orcid.org/0009-0000-7637-1818)
