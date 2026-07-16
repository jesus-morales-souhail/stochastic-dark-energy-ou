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
| `scripts/ou_bao_likelihood.py` | OU / QNM BAO residual kernel on DESI DR2 (main null + σ_X limit) |
| `scripts/eos_efectiva.py` | Clean CPL (+ optional nested σ,θ) on same BAO alphas; multi-start Nelder–Mead |
| `scripts/ou_bao_stochastic_test.py` | Alternate OU / QNM BAO likelihood pipeline |
| `scripts/cross_correlation_DESI.py` | Galaxy × residual cross-correlation (preliminary) |
| `scripts/desqueezing/desqueezing_relax_time.py` | Lindblad desqueezing half-life scan |
| `scripts/desqueezing/cosmological_mapping.py` | Map γ ↔ θ H(z) using published anchors |
| `scripts/gpe/gpe_sim.py` | Negative-mass GPE illustration |
| `results/REAL_DATA_SOLID.md` | Claim-safe summary of what real data support |

See `README.md` for installation and scientific overview.
