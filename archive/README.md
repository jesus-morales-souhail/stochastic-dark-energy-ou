# archive/ — retired material

Everything here was moved on 2026-10-17 (git history is preserved; use `git log --follow <file>`). It is kept for the record. **It is not maintained, it may not run, and none of it supports the preprint's claims.** Each `.py` and `.md` file carries a two-line header with its previous path and the reason; for the other files the reason is listed below.

Paths below are relative to `archive/`; the previous path is the same without the `archive/` prefix unless stated.

## Scripts

| File | Reason |
|:--|:--|
| `scripts/ou_bao_stochastic_test.py` | DR1 variant of the OU/QNM likelihood with hard-coded data and diagonal covariance. The preprint uses `scripts/ou_bao_likelihood.py` (DR2, full covariance). |
| `scripts/cross_correlation_DESI.py` | Preliminary DESI galaxies × SN residual cross-correlation (non-claim N5; PREPRINT §6 mentions it as preliminary). Paths depend on the working directory. |
| `scripts/euclid_mock_mcmc.py` | Deprecated stub of the synthetic Euclid MCMC (N9); it only forwarded to `scripts/desi_dr2_real_bao_test.py`, which no longer sits next to it. |
| `scripts/euclid_joint_bao_sne_mcmc.py` | Same, for the synthetic Euclid BAO+SN joint MCMC. |
| `scripts/test_desi_dr2_correlations.py` | Stand-alone DR2 correlation test, not used by the preprint; writes to `plots/`. |
| `analisis_cmb/diagnosticos/obsoleto/realizaciones_H0fijo_mead2020.py` | Obsolete CMB diagnostic (fixed \(H_0\), HMCode mead2020); superseded by `analisis_cmb/realizaciones.py`. |
| `analisis_cmb/diagnosticos/obsoleto/tabla_H0fijo.py` | Table for that obsolete run. |

## Results

| File | Reason |
|:--|:--|
| `results/euclid_mcmc/` (13 files: `RETIRED.md`, 10 PNG, `mcmc_summaries.json`, `mcmc_summary.txt`) | Synthetic Euclid-like BAO vector, not data (non-claim N9). |
| `results/euclid_joint_mcmc/` (10 files: `RETIRED.md`, 7 PNG, `joint_summaries.json`, `joint_summary.txt`) | Synthetic Euclid-like BAO+SN vector, not data (N9). |
| `results/ou_bao_desi_dr2_run.txt` | Console log of an old OU/QNM run with diagonal covariance (θ = 0.001, σ_X = 5e-5); its log-likelihood does not match the full-covariance pipeline. |

## Figures (previously in `figures/`)

None of these is cited by the preprint, and no current script produces them.

| File | Origin |
|:--|:--|
| `figures/Figure_1.png`, `figures/Figure_2.png` | Uploaded by hand (commit `500f05b`, 2026-07-13). |
| `figures/falsos_positivos_detectores.png` | Uploaded by hand (`826dd3e`, 2026-07-06). |
| `figures/exclusion_plot.png` | Uploaded by hand (`f247fcb`, 2026-07-05). The function that would draw it in `scripts/ou_bao_likelihood.py` is never called. |
| `figures/test_desi_QNM.png` | Written only by the archived DR1 `ou_bao_stochastic_test.py`. |
| `figures/joint_mcmc_E1_amplitude_BAO.png`, `figures/joint_mcmc_E2_mean_reversion_BAO.png`, `figures/joint_mcmc_null_F_E0_BAO.png` | Synthetic Euclid joint MCMC (`470619a`, 2026-07-15); a different version from the files of the same name in `results/euclid_joint_mcmc/`. |

The 14 figures in `figures/` that were byte-identical copies of `results/euclid_mcmc/` and `results/euclid_joint_mcmc/` were removed, not archived; the originals are in `results/` above.
