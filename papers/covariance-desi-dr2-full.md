# DESI DR2 measurement covariance in this repository

Jesús Morales Souhail · [github.com/jesus-morales-souhail](https://github.com/jesus-morales-souhail)

## What “full covariance” means here

DESI DR2 Results II (arXiv:[2503.14738](https://arxiv.org/abs/2503.14738)) and the public Cobaya-style pack
(`CobayaSampler/bao_data`, also in the Zenodo DR2 supplementary) ship a **Gaussian BAO** mean vector and
a **13×13 covariance** for compressed distances:

| Index | Observable |
|------:|------------|
| 0 | $D_V/r_s$ (BGS, $z\sim 0.295$) |
| 1–2 | $D_M/r_s$, $D_H/r_s$ (LRG, $z\sim 0.51$) |
| 3–4 | LRG $z\sim 0.71$ |
| 5–6 | LRG+ELG $z\sim 0.93$ |
| 7–8 | ELG $z\sim 1.32$ |
| 9–10 | QSO $z\sim 1.48$ |
| 11–12 | Lyα $D_H$, $D_M$ ($z\sim 2.33$; order DH then DM in the mean file) |

**Structure of $C_{13}$ (public file):** block-diagonal across tracer bins. Within each anisotropic bin,
$D_M$ and $D_H$ are anti-correlated ($|\rho|\sim 0.35$–$0.50$). **Cross-bin blocks are zero** in this
compressed likelihood. That is a property of the public product, not a shortcut invented here.

## Projection onto isotropic $\alpha$ (7 bins)

Pipelines in this repo work on isotropic $\alpha(z)$ from figure6. We attach the measurement covariance

$$
C_\alpha = J\, C_{13}\, J^{\mathsf{T}},
$$

with $J$ the Jacobian of $\alpha_{\mathrm{iso}}$ w.r.t. the 13 distances (BGS: $\alpha\propto D_V$;
anisotropic: $\alpha\propto D_M^{2/3}D_H^{1/3}$). Because $C_{13}$ is block-diagonal across bins,
$C_\alpha$ is **diagonal across the 7 $\alpha$ bins**; each variance correctly folds the within-bin
$D_M$–$D_H$ correlation.

Loader: `scripts/desi_dr2_data.py`  
Public copy: `results/desi_cov/desi_cov_alpha_iso_7x7.txt`  
Also: `results/desi_cov/desi_gaussian_bao_ALL_GCcomb_cov.txt` (13×13).

## Residual likelihood

For OU / QNM / rank-1 residual models,

$$
C_{\mathrm{tot}} = C_\alpha + C_{\mathrm{signal}}(\theta,\sigma_X,\ldots).
$$

$C_{\mathrm{signal}}$ can (and does) introduce **off-diagonal residual correlations** between bins;
that is the model, not the measurement noise.

## What is still not public here

- **Multipole-level** $\xi_\ell(s)$ full covariances (RascalC) are not in the small figure5 text tables
  (those tables only list per-bin $\xi$ errors). T2 multipole residual $\chi^2$ therefore uses the
  published diagonal $\xi$ errors — honest limitation, not a mock.
- A single joint non-block-diagonal $\alpha$ covariance across all redshifts is **not** in the public
  compressed Gaussian BAO product used above.

## Scripts updated to use $C_\alpha$

| Script | Role |
|--------|------|
| `desi_dr2_data.py` | load + project |
| `desi_dr2_real_bao_test.py` | real $\alpha$ residual scan |
| `profile_sigma_x_desi.py` | profile $\sigma_X$ |
| `joint_w0wa_sigma_desi.py` | joint $\{w_0,w_a,\sigma_X\}$ |
| `ou_bao_likelihood.py` | OU/QNM comparison |
| `ou_bao_stochastic_test.py` | alternate OU/QNM |
| `tachyonic_rank1_mle.py` | rank-1 growth exclusion |
| `eos_efectiva.py` | CPL + nested extension |

## Reproduce

```bash
python -c "from desi_dr2_data import load_alpha_dv; d=load_alpha_dv(); print(d['cov_mode'], d['cov'].shape)"
python scripts/profile_sigma_x_desi.py
python scripts/desi_dr2_real_bao_test.py
python scripts/tachyonic_rank1_mle.py
```
