# Stochastic Dark Energy: Finite Information Network Analysis

**Author:** Jesús Morales Souhail  
**Contact:** jmskjym@gmail.com  
**Version:** v3.0 (July 2026) | **Status:** Preprint — not peer reviewed

This repository contains the numerical tools, statistical analysis scripts, and
technical documentation to evaluate the **Stochastic Dark Energy (SDE)**
framework. This phenomenological model postulates that the cosmological constant
is not strictly static, but experiences stochastic fluctuations motivated by the
finite information content of the observable universe (Bekenstein–Hawking bound)
and the ~1/√N fluctuation scaling suggested by unimodular gravity and causal-set
arguments.

Fluctuations in Ω_Λ are modeled as an **Ornstein–Uhlenbeck (OU) process**
in x = ln a, with mean-reversion rate θ and diffusion amplitude σ.

---

## Repository Structure

```
stochastic-dark-energy-ou/
├── paper/
│   └── stochastic_dark_energy_v3.md   ← full preprint (v3.0)
├── code/
│   ├── ou_bao_likelihood.py           ← BAO likelihood + lag-correlation test
│   └── cross_correlation_DESI.py      ← HEALPix cross-correlation pipeline
├── docs/
│   ├── sensitivity_kernel_table.md    ← S(z) for all DESI tracers
│   └── INSTRUCCIONES_EXPERIMENTO.md   ← full protocol
├── results/
│   └── RESUMEN_COMPLETO_FASES.txt     ← verified results (DR1)
├── requirements.txt
└── README.md
```

---

## Key Results (July 2026)

### Test 1 — BAO Likelihood with OU Covariance

| Model | θ | σ_X | ΔlogL vs ΛCDM | AIC | BIC |
|---|---|---|---|---|---|
| **OU pure** | 0.765 | 0.018 | **+6.75** | −38.73 | −38.84 |
| QNM (oscillatory) | 0.000 | 0.019 | +7.00 | −37.24 | −37.40 |
| ΛCDM | — | — | 0 | ref | ref |

AIC and BIC favour OU over ΛCDM. **Preliminary with 7 bins.**

### Test 2 — Angular Cross-Correlation (δ_g × δΩ_Λ)

DESI DR1 LRG (1,476,135 + 662,000 objects) × Pantheon+ residuals (1,701 SNe):

```
r_cross = 0.1673 ± 0.0613    (Z ≈ 2.73σ,  67 overlapping pixels)
OU prediction:  r_pred ~ 0.023
```

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
```

**`requirements.txt`:**
```
numpy>=1.24
scipy>=1.10
astropy>=5.3
healpy>=1.16
matplotlib>=3.7
```

Verify:
```bash
python -c "import numpy, scipy, astropy, healpy, matplotlib; print('All OK')"
```

---

## Quickstart

### Test 1 — BAO likelihood (no catalog needed)

Uses hardcoded public DESI DR2 BAO numbers from arXiv:2503.14738:

```bash
python code/ou_bao_likelihood.py
```

### Test 2 — Cross-correlation (requires catalogs)

```bash
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
```

If catalogs are absent, the script runs in **synthetic validation mode**.

---

## Data Sources

| Dataset | Objects | Source |
|---|---|---|
| DESI DR1 LRG NGC | 1,476,135 | data.desi.lbl.gov/public/edr/ |
| DESI DR1 LRG SGC | 662,000 | data.desi.lbl.gov/public/edr/ |
| DESI DR1/DR2 BAO | 7 bins | arXiv:2404.03000, arXiv:2503.14738 |
| Pantheon+ SNe Ia | 1,701 | github.com/PantheonPlusSH0ES |

DESI DR2 LRG catalogs with imaging systematic weights (WEIGHT_SYS) have been
requested from **NSF NOIRLab Astro Data Lab** for the cross-correlation
systematic analysis.

---

## Falsification Criteria

| Criterion | Condition | Status |
|---|---|---|
| F1 — Variance floor | σ²_obs < σ²_floor in multiple bins | Alive |
| F2 — w_a frozen | w_a → 0 at >5σ | Alive |
| F3 — ISW | σ_Ω_Λ incompatible with CMB-LSS | Pending |
| F4 — Bayes | ln K > 5 favouring ΛCDM | Not computed |
| F5 — Lag | All lags ≤ 0 with 20+ bins | Partial (7 bins) |

---

## Citation

```bibtex
@misc{moralesssouhail2026,
  author = {Morales Souhail, Jesús},
  title  = {Stochastic Dark Energy from a Finite Information Network},
  year   = {2026},
  url    = {https://github.com/[username]/stochastic-dark-energy-ou},
  note   = {Preprint v3.0}
}
```

---

## License

Code: MIT | Text: CC BY 4.0

---

## References

1. DESI Collaboration, arXiv:2404.03000 (2024)
2. DESI Collaboration, arXiv:2503.14738 (2025)
3. Scolnic et al. (Pantheon+), ApJ 938, 113 (2022)
4. Bekenstein, Phys. Rev. D 7, 2333 (1973)
5. Sorkin, arXiv:gr-qc/0503057 (2005)
6. Uhlenbeck & Ornstein, Phys. Rev. 36, 823 (1930)
