# stochastic-dark-energy-ou

Ornstein-Uhlenbeck model for dark energy fluctuations tested against DESI DR2 BAO data, together with a quantum fluid instability analysis.

## Overview

This repository contains research on stochastic components in dark energy using public DESI DR2 Baryon Acoustic Oscillation data.

**Main scientific results:**
- Strong null result for stationary stochastic fluctuations modeled as an Ornstein–Uhlenbeck process.
- Exclusion of globally coherent tachyonic quantum fluid instabilities.
- Proposal of the **Principle of Late-Time Vacuum Homogeneity**: at observable cosmological scales, the dark energy sector behaves as a perfectly smooth, non-stochastic background.

## Main Documents

| File                                | Description |
|-------------------------------------|-----------|
| `stochastic_dark_energy.md`         | Main paper. Ornstein–Uhlenbeck + QNM analysis with DESI DR2. Includes corrected statistical tests and the Principle of Vacuum Homogeneity. |
| `quantum-fluid-instabilities-desi-dr2.md`                          | Companion paper testing a tachyonic quantum fluid model with Bogoliubov instability using the correct rank-1 growing covariance. |
| `principle-of-vacuum-smoothness.md` | Short document stating the proposed principle that future theoretical models should start from σ_X = 0. |

## Analysis Code

- `ou_bao_likelihood.py` — Main analysis pipeline (BAO likelihood, MLE fitting, lag correlations, and sensitivity kernel).

## Key Findings

- When the cosmological background is allowed to vary, stationary stochastic models are not favored by DESI DR2 data.
- A globally coherent growing tachyonic mode is strongly disfavored.
- The data support treating **σ_X = 0** as a working principle for model building in dark energy and quantum gravity.

## Data Sources

- DESI DR2 BAO: https://data.desi.lbl.gov/public/
- Pantheon+ Supernovae: https://github.com/PantheonPlusSH0ES/DataRelease

## Status

Preprint (not peer-reviewed). Actively developed.

## Author

Jesús Morales Souhail  
ORCID: [0009-0000-7637-1818](https://orcid.org/0009-0000-7637-1818)
