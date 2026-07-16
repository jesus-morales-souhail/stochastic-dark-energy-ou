# Repository structure

```
papers/       Scientific notes and numerical summary (English)
scripts/      Analysis code
  desqueezing/   Open-system desqueezing and timescale mapping
  gpe/           Gross–Pitaevskii / Bogoliubov numerics
figures/      Figures
results/      Numerical outputs (BAO fits, forecasts, desqueezing tables)
notes/        Technical notes (desqueezing synthesis, mapping tables)
```

## Entry points

| Path | Role |
|------|------|
| `scripts/ou_bao_likelihood.py` | OU / QNM residual kernel on DESI DR2 BAO |
| `scripts/eos_efectiva.py` | CPL background and nested (σ, θ) extension |
| `scripts/ou_bao_stochastic_test.py` | Alternate OU / QNM BAO pipeline |
| `scripts/cross_correlation_DESI.py` | Galaxy × residual cross-correlation |
| `scripts/desqueezing/desqueezing_relax_time.py` | Lindblad desqueezing half-life scan |
| `scripts/desqueezing/cosmological_mapping.py` | Map γ ↔ θ H(z) |
| `scripts/gpe/gpe_sim.py` | Negative-mass GPE illustration |
| `papers/resume.txt` | Compact numerical summary |
| `results/eos_cpl_desi_dr2/` | CPL fit outputs |

See `README.md` for installation and overview.
