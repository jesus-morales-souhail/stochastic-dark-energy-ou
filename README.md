# Stochastic Dark Energy: Constraints from DESI DR2

**Author:** Jesús Morales Souhail  
**Contact:** jmskjym@gmail.com  
**ORCID:** 0009-0000-7637-1818
**Repository:** https://github.com/AshPokemonTCG/stochastic-dark-energy-ou

This repository contains the numerical tools, statistical analysis scripts, and
technical documentation to test the **Stochastic Dark Energy (SDE)** hypothesis
against public DESI DR2 BAO data. The SDE framework postulates that the
cosmological constant is not strictly static, but experiences small stochastic
fluctuations motivated by the finite information content of the observable
universe (Bekenstein–Hawking bound) and the ~1/√N fluctuation scaling suggested
by unimodular gravity and causal-set arguments.

Fluctuations in Ω_Λ are modeled as an **Ornstein–Uhlenbeck (OU) process**
in x = ln a, with mean-reversion rate θ and diffusion amplitude σ. A
**quasi-normal mode (QNM) extension** with an oscillatory kernel is also
tested.

---

## Key Results (DESI DR2, July 2026)

### Test 1 — BAO Likelihood with OU Covariance

Using DESI DR2 BAO data (7 bins, arXiv:2503.14738), the maximum likelihood
estimation (MLE) forces the stochastic amplitude to zero:

| Model | θ | σ_X | ω_R | ΔlogL vs ΛCDM | AIC | BIC |
|---|---|---|---|---|---|---|
| **ΛCDM** | — | — | — | 0.00 (ref) | ref | ref |
| **H0: OU free MLE** | 0.001 | 5 × 10⁻⁵ | 0 (fixed) | 0.00 | -50.03 | -50.14 |
| **H1: QNM free MLE** | 0.001 | 5 × 10⁻⁵ | 0.00 | 0.00 | -48.03 | -48.19 |

**Conclusion:** There is **no evidence** for stochastic fluctuations in the
DESI DR2 BAO data. The data are fully consistent with smooth CPL evolution
(\(w_0 \approx -0.87, w_a \approx -0.41\)) plus instrumental noise.

**Strict upper limit (95% CL):** \(\sigma_X < 1.5 \times 10^{-4}\).

### Test 2 — Angular Cross-Correlation (δ_g × δΩ_Λ)

DESI DR1 LRG (1,476,135 + 662,000 objects) × Pantheon+ residuals (1,701 SNe):
r_cross = 0.1673 ± 0.0613 (Z ≈ 2.73σ, 67 overlapping pixels)
OU prediction: r_pred ~ 0.023

text

**Preliminary.** Imaging systematics not yet controlled. DESI DR2
with WEIGHT_SYS column required for systematic-controlled reanalysis.

### Test 3 — Lag Correlations in BAO Residuals

| Lag | DR2 (7 bins) | OU Prediction | 95% CI | Significant? |
|---|---|---|---|---|
| 1 | +0.37 | +0.83 | ±1.0 | No |
| 2 | −0.32 | +0.85 | ±1.0 | No |
| 3 | −0.81 | +0.85 | ±1.0 | No |

None individually significant at 7 bins. Decisive test: >20 bins (Euclid DR1,
H2 2026).

---

## Requirements and Installation

Python 3.8+

```bash
pip install -r requirements.txt
requirements.txt:

text
numpy>=1.24
scipy>=1.10
astropy>=5.3
healpy>=1.16
matplotlib>=3.7
Verify:

bash
python -c "import numpy, scipy, astropy, healpy, matplotlib; print('All OK')"
Quickstart
Test 1 — BAO likelihood (no catalog needed)
Uses hardcoded public DESI DR2 BAO numbers from arXiv:2503.14738:

bash
python code/ou_bao_likelihood.py
Test 2 — Cross-correlation (requires catalogs)
bash
# Download DESI DR1 LRG (~200 MB total)
wget -O data/LRG_clustering_N.fits \
  https://data.desi.lbl.gov/public/edr/vac/edr/lss/v2.0/LSScats/LRG_clustering_N.fits
wget -O data/LRG_clustering_S.fits \
  https://data.desi.lbl.gov/public/edr/vac/edr/lss/v2.0/LSScats/LRG_clustering_S.fits

# Download Pantheon+ (~1 MB)
wget -O data/PantheonPlus.dat \
  "https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/Pantheon+SH0ES.dat"

# Run pipeline
python code/cross_correlation_DESI.py \
    --lrg_n data/LRG_clustering_N.fits \
    --lrg_s data/LRG_clustering_S.fits \
    --sne   data/PantheonPlus.dat
If catalogs are absent, the script runs in synthetic validation mode.

Data Sources
Dataset	Objects	Source
DESI DR1 LRG NGC	1,476,135	data.desi.lbl.gov/public/edr/
DESI DR1 LRG SGC	662,000	data.desi.lbl.gov/public/edr/
DESI DR1/DR2 BAO	7 bins	arXiv:2404.03000, arXiv:2503.14738
Pantheon+ SNe Ia	1,701	github.com/PantheonPlusSH0ES
DESI DR2 LRG catalogs with imaging systematic weights (WEIGHT_SYS) have been
requested from NSF NOIRLab Astro Data Lab for the cross-correlation
systematic analysis.

Falsification Criteria
Criterion	Condition	Status
F1 — Variance floor	σ²_obs < σ²_floor in multiple bins	Alive
F2 — w_a frozen	w_a → 0 at >5σ	Alive
F3 — ISW	σ_Ω_Λ incompatible with CMB-LSS	Pending
F4 — Bayes	ln K > 5 favouring ΛCDM	Not computed
F5 — Lag	All lags ≤ 0 with 20+ bins	Partial (7 bins)
Citation
bibtex
@misc{moralesssouhail2026constraints,
  author = {Morales Souhail, Jesús},
  title  = {Constraints on Stochastic Dark Energy from DESI DR2},
  year   = {2026},
  url    = {https://github.com/AshPokemonTCG/stochastic-dark-energy-ou},
}
License
Copyright (c) 2026 Jesús Morales Souhail

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

References
DESI Collaboration, arXiv:2404.03000 (2024)

DESI Collaboration, arXiv:2503.14738 (2025)

Scolnic et al. (Pantheon+), ApJ 938, 113 (2022)

Bekenstein, Phys. Rev. D 7, 2333 (1973)

Sorkin, arXiv:gr-qc/0503057 (2005)

Uhlenbeck & Ornstein, Phys. Rev. 36, 823 (1930)
