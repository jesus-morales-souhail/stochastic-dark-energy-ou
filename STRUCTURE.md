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
| `papers/anisotropic-slip-option0.md` | Option 0: η / anisotropic stress vs σ_X (literature + scaling; no Boltzmann) |
| `papers/no-go-superoscillation-tesseract.md` | No-go: Born tax + undeclared tesseract optical power |
| `papers/self-shielding-vs-untestability.md` | SDiff “shielding” vs amplitude/operator honesty |
| `papers/pattern-undeclared-physical-power.md` | Method hygiene: ln4 / B4 style claims |
| `scripts/superoscillation_energy_cost_demo.py` | 1D demo: E_core collapses as superoscillation grows |
| `scripts/car_drone_pupil_newton_einstein.py` | Car+drone@120 km/h, pupil diffraction: Newton vs SR + Hubble |
| `papers/car-drone-pupil-newton-einstein.md` | Write-up of the car/drone/pupil experiment (no tesseract) |
| `papers/scale-operator-experiment-map.md` | Closed wrong-scale results + correct next experiments |
| `papers/data-pack-option0-internet.md` | Internet-sourced numbers + arXiv links for Option 0 |
| `papers/optics-ou-analogies-and-limits.md` | Quantum lens OK; tesseract ABCD no; f↔1/θ analogy only |
| `papers/self-shielding-triple-barrier.md` | 0% physical / 100% math: Born, band-limit, OU shields |
| `papers/wavefront-shaping-vs-ou-vacuum.md` | Lab T-matrix optics ≠ BAO OU; no vacuum SLM |
| `results/eos_cpl_desi_dr2/` | CPL fit outputs |

See `README.md` for installation and overview.

### Separated (not this repo)

| Path | Content |
|------|--------|
| `~/Proyectos/04_Optica_medios_complejos/` | Lab T-matrix wavefront shaping experiment (not DE) |

