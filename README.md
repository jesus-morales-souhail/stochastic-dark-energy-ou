# stochastic-dark-energy-ou

Testing Ornstein-Uhlenbeck fluctuations in dark energy against DESI DR2 BAO data, plus a check for tachyonic instabilities in a quantum fluid model.

## What this repo contains

This is an independent analysis using public DESI DR2 Baryon Acoustic Oscillation measurements. The goal was to see whether the data prefer any kind of stochastic component on top of dynamical dark energy.

The main findings so far are:

- No sign of stationary stochastic fluctuations (Ornstein–Uhlenbeck type). When the background cosmology is allowed to vary, the amplitude of fluctuations is driven to zero.
- A globally coherent growing tachyonic mode (from a simple quantum fluid model) is strongly disfavored by the data.
- This leads to a working hypothesis: at the redshifts and scales currently probed by BAO, the dark energy sector appears remarkably smooth. I call this the **Principle of Late-Time Vacuum Homogeneity**.

## Main documents

| File | Content |
|------|---------|
| `stochastic-dark-energy-desi-dr2.md` | Main analysis: OU + QNM models fitted to DESI DR2 |
| `quantum-fluid-instabilities-desi-dr2.md` | Test of tachyonic instability with rank-1 growing covariance |
| `principle-of-vacuum-smoothness.md` | Short note on the proposed principle |

## Code

- `ou_bao_stochastic_test.py` — Main likelihood pipeline (BAO covariance with OU/QNM kernels, MLE fitting)

## Data

- DESI DR2 BAO: https://data.desi.lbl.gov/public/
- Pantheon+ supernovae (for cross-checks)

## Status

This is work in progress (July 2026). Not peer-reviewed. The analysis is being extended for future datasets like Euclid DR1.

## Author

Jesús Morales Souhail  
ORCID: [0009-0000-7637-1818](https://orcid.org/0009-0000-7637-1818)
