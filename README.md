Stochastic Dark Energy from a Finite Information Network: BAO Covariance Signatures, Cross-Correlation Evidence, and Falsifiable Predictions
# Stochastic Dark Energy: Ornstein-Uhlenbeck Model vs. DESI BAO

**Author:** Jesús Morales Souhail  
**Contact:** jmskjym@gmail.com  
**ORCID:** [tu ORCID]  
**Repository:** https://github.com/[tu-usuario]/stochastic-dark-energy-ou

---

## What this is

A phenomenological test of whether dark energy fluctuations follow an
Ornstein–Uhlenbeck (OU) stochastic process, calibrated against public DESI DR1
BAO measurements. The model produces a specific, falsifiable signature:
**exponentially-decaying redshift correlations in BAO residuals**, expressed as
an additive covariance component \(C_{\rm OU}\) that can be tested immediately
against public DESI data products.

This is not a claim that dark energy *is* stochastic. It is a model with a
concrete prediction that can be confirmed or ruled out with existing data.

---

## Key results (reproducible from public data)

### BAO likelihood test

Using DESI DR2 BAO measurements (7 redshift bins, \(z = 0.295\)–\(2.33\)):

| Model               | \(\theta\) | \(\sigma_X\) | \(\Delta\log\mathcal{L}\) vs ΛCDM | AIC      | BIC      |
|:--------------------|:----------:|:------------:|:-----------------------------------:|:--------:|:--------:|
| **H0: OU pure**     | 0.765      | 0.018        | **+6.75**                           | −38.73   | −38.84   |
| H1: QNM (oscillatory)| 0.000     | 0.019        | +7.00                               | −37.24   | −37.40   |
| ΛCDM baseline       | —          | —            | 0.00                                | ref      | ref      |

→ **AIC/BIC prefer H0 (pure OU) over H1 (QNM) by 1.5 points.**  
→ With only 7 bins, all results are *preliminary*. The decisive test requires 20+ bins (Euclid DR1, October 2026, or full DESI DR2 with systematic weights).

### Cross‑correlation test (galaxy density × Hubble residuals)

Preliminary angular cross‑correlation between DESI DR1 LRG maps and
Pantheon+ SNe Ia residuals as a proxy for \(\delta\Omega_\Lambda\):

| Dataset             | Objects       | Overlap pixels | \(r_{\rm cross}\) | Z‑score |
|:--------------------|:--------------|:--------------:|:-----------------:|:-------:|
| NGC + SGC combined  | 1,278,805 LRG | 67             | \(0.1673 \pm 0.0613\) | **2.73σ** |

Naive OU prediction: \(r_{\rm cross} \sim 0.023\). The observed value exceeds this by a factor of \(\sim 7\).

**Honest assessment:** Two systematic explanations are not yet ruled out:

1. **Galactic dust:** extinction affects both the LRG photometric redshifts and the SN magnitudes.
2. **NGC/SGC imaging systematics:** BASS+MzLS vs DECaLS photometric systems.

**This is why we are requesting DESI DR2 access:** larger overlap, proper imaging systematic weights (`WEIGHT_SYS`), and the lag‑correlation test on whitened BAO residuals described in the paper.

### Precision floor and the diagnostic bin \(z = 0.934\)

Using the exact sensitivity kernel \(S(z) \equiv \partial \ln D_V / \partial \Omega_\Lambda\), the predicted irreducible BAO scatter ("precision floor") is:

\[
\sigma_{\alpha,\rm floor}(z) = \sigma_{\alpha,\rm floor}(0.706) \times \frac{|S(z)|}{|S(0.706)|}, \qquad \sigma_{\alpha,\rm floor}(0.706) = 4.65\times 10^{-3}.
\]

| \(z_{\rm eff}\) | Tracer       | \(|S(z)|/|S(0.706)|\) | \(\sigma_{\alpha,\rm floor}\) |
|:---------------:|:-------------|:---------------------:|:-----------------------------:|
| 0.295           | BGS          | 0.477                 | \(2.22 \times 10^{-3}\)       |
| 0.510           | LRG1         | 0.777                 | \(3.61 \times 10^{-3}\)       |
| **0.706**       | **LRG2**     | **1.000**             | **\(4.65 \times 10^{-3}\)**   |
| **0.934**       | **LRG3+ELG1**| **1.208**             | **\(5.62 \times 10^{-3}\)**   |
| 1.321           | ELG2         | 1.462                 | \(6.80 \times 10^{-3}\)       |
| 1.484           | QSO          | 1.541                 | \(7.17 \times 10^{-3}\)       |
| 2.330           | Ly\(\alpha\) | 1.798                 | \(8.36 \times 10^{-3}\)       |

**Diagnostic bin \(z = 0.934\) (LRG3+ELG1):** the predicted OU floor exceeds the DESI DR2 observational error (\(\sigma_{\alpha,\rm obs}=0.0049\)) by a margin that translates into a \(\sim 2.8\sigma\) sensitivity excess. This is the cleanest test available with current DESI DR2 data.

---

## The falsifiable predictions

### 1. Lag correlations in BAO residuals (smoking gun)

For DESI‑like redshift binning, the OU model predicts strong positive lag correlations in whitened BAO residuals:
Lag-1: ρ₁ ≈ 0.78
Lag-2: ρ₂ ≈ 0.62
Lag-3: ρ₃ ≈ 0.49

text

Under ΛCDM (no OU), these correlations should be consistent with zero.  
The test requires only public BAO data vectors and covariance matrices already released by DESI — **no proprietary data needed for this test**.

### 2. Rayleigh resolution limit for Euclid DR1

A geometric result independent of the model: the minimum resolvable frequency in logarithmic scale‑factor space is

\[
\omega_{R,\min} = \frac{2\pi}{\Delta x}, \qquad \Delta x = \ln\frac{1+z_{\max}}{1+z_{\min}}.
\]

| Survey       | \(\Delta x\) | \(\omega_{R,\min}\) | What it can test |
|:-------------|:------------:|:-------------------:|:-----------------|
| DESI DR2     | 0.944        | 6.66                | Oscillatory modes with \(\omega_R > 6.66\) |
| Euclid DR1   | ~0.39–0.50   | ~12.6–16.2          | Only modes with \(\omega_R > 12.6\) |

**If the true kernel is an oscillatory QNM with \(\omega_R \sim 8\), DESI can detect it, but Euclid DR1 geometrically cannot.** Euclid will see the oscillation as a straight line (pure OU slope) because its shorter lever arm in \(\ln a\) cannot resolve the curvature. This is a radical inversion of the usual "more data solves everything" paradigm.

---

## Repository structure
README.md ← this file
paper/
stochastic_dark_energy_v3.0.md ← full preprint (English)
code/
cross_correlation_DESI.py ← HEALPix pipeline (galaxy × SNe Ia)
ou_bao_likelihood.py ← BAO likelihood with C_OU covariance
docs/
sensitivity_kernel_table.md ← S(z) values for all DESI tracers
INSTRUCCIONES_EXPERIMENTO.md ← full experimental protocol (optional)
results/
RESUMEN_COMPLETO_FASES.txt ← verified numerical results (DR1)

text

---

## Data sources (all public)

| Dataset                     | Size        | Source |
|:----------------------------|:------------|:-------|
| DESI DR1 LRG NGC            | 136.56 MB   | data.desi.lbl.gov/public/edr/ |
| DESI DR1 LRG SGC            | 61.3 MB     | data.desi.lbl.gov/public/edr/ |
| DESI DR1 BAO measurements   | —           | arXiv:2404.03000, Table 1 |
| DESI DR2 BAO measurements   | —           | arXiv:2503.14738 |
| Pantheon+ SNe Ia            | 1,701 obj.  | github.com/PantheonPlusSH0ES |

**Requested from Astro Data Lab (DESI DR2):**

```sql
SELECT RA, DEC, Z, WEIGHT, WEIGHT_SYS, PHOTSYS
FROM desi_dr2.lss.lrg_clustering
WHERE Z BETWEEN 0.4 AND 0.8 AND ZWARN = 0
Installation and reproduction
bash
git clone https://github.com/[tu-usuario]/stochastic-dark-energy-ou
cd stochastic-dark-energy-ou
pip install numpy scipy astropy healpy matplotlib

# Run BAO likelihood test (no catalog needed — uses public DR2 numbers)
python code/ou_bao_likelihood.py

# Run cross-correlation pipeline (requires LRG + Pantheon+ data)
python code/cross_correlation_DESI.py
All numerical results in results/ are fully reproducible from public data.

References
DESI Collaboration, "DESI 2024 III: Baryon Acoustic Oscillations from Galaxies and Quasars," arXiv:2404.03000 (2024).

DESI Collaboration, "DESI DR2 Results II: Measurements of BAO and Cosmological Constraints," arXiv:2503.14738 (2025).

Scolnic, D. et al. (Pantheon+ Collaboration), "The Pantheon+ Analysis: Cosmological Constraints," ApJ 938, 113 (2022).

Uhlenbeck, G. E. and Ornstein, L. S., "On the Theory of the Brownian Motion," Phys. Rev. 36, 823 (1930).

Bekenstein, J. D., "Black Holes and Entropy," Phys. Rev. D 7, 2333 (1973).

Sorkin, R. D., "Is the Cosmological 'Constant' a Nonlocal Quantum Residual?," arXiv:gr-qc/0503057 (2005).

Acknowledgments
This work used public data products from the Dark Energy Spectroscopic Instrument (DESI) collaboration and the Pantheon+ supernova compilation.
NSF NOIRLab Astro Data Lab access is requested for the DESI DR2 systematic analysis described in Section 6 of the paper.
