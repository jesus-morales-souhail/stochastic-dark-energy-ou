# Repository structure

```
papers/       Scientific notes and summary of results (English)
scripts/      Analysis and simulation code
  desqueezing/   Open-system desqueezing and cosmological timescale mapping
  gpe/           Gross–Pitaevskii / Bogoliubov numerics
figures/      Figures
results/      Numerical outputs associated with the desqueezing analysis
notes/        Additional technical notes (desqueezing synthesis, mapping tables)
```

## Entry points

| Path | Role |
|------|------|
| `scripts/ou_bao_stochastic_test.py` | OU / QNM BAO likelihood pipeline |
| `scripts/ou_bao_likelihood.py` | Core BAO likelihood utilities |
| `scripts/cross_correlation_DESI.py` | Galaxy × residual cross-correlation pipeline |
| `scripts/desqueezing/desqueezing_relax_time.py` | Lindblad desqueezing half-life scan |
| `scripts/desqueezing/cosmological_mapping.py` | Map γ ↔ θ H(z) using published anchors |
| `scripts/gpe/gpe_sim.py` | Negative-mass GPE illustration |

See `README.md` for installation and scientific overview.
