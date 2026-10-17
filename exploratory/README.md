# exploratory/ — desqueezing programme

This folder holds the vacuum-relaxation ("desqueezing") programme: QuTiP models of how a squeezed mode relaxes, the mapping from those half-lives to cosmological numbers, the SDiff discrimination grid, the Euclid design forecasts, and the Option 0 data pack.

**It is exploratory.** None of it is cited as a DESI result, and it is outside the preprint's claim set C1–C7 (see [`../manuscript/CLAIMS.md`](../manuscript/CLAIMS.md), N6, and [`../theory/EXPLORATORY_BOUNDARY.md`](../theory/EXPLORATORY_BOUNDARY.md)). The only input taken from the claim set is the DESI ceiling on \(\sigma_X\), read from `results/profile_sigma_x/profile_sigma_x.json` through `scripts/desi_bound.py`.

## Layout

```
scripts/   code (run from anywhere; outputs go to exploratory/results/)
results/   outputs, one folder per script
notes/     working notes and the Option 0 data pack
```

| Script | Output folder | What it does |
|:--|:--|:--|
| `scripts/desqueezing_relax_time.py` | `results/desqueezing_relax_time/` | QuTiP scan of the half-life \(t_{1/2}(\gamma, n_{\rm th})\) |
| `scripts/desqueezing_heavy_scan.py` | `results/desqueezing_heavy/` | Wider \((\gamma, n_{\rm th}, r)\) grid for the half-life of \(|\langle a^2\rangle|\). **Known issue:** the script still compares against the old formula \(\ln 2/[\gamma(1+2n_{\rm th})]\); the versioned outputs record the corrected law \(t_{1/2}=\ln 2/\gamma\) (see `analisis_cmb/REORGANIZACION_CHECK.md`) |
| `scripts/termico_vs_squeezed_recomendado.py` | `results/termico_vs_squeezed_recomendado/` | Thermal vs squeezed baselines |
| `scripts/cosmological_mapping.py` | `results/cosmological_mapping/` | Gyr / Sorkin / Euclid / DESI mapping tables |
| `scripts/sdiff_discrimination.py` | `results/sdiff_discrimination/` | Residual vs \(\theta\), minimum \(A_0\) for detection |
| `scripts/euclid_protocol_forecasts.py` | `results/euclid_protocol/` | Euclid design forecasts over \((\theta, A_0)\) |

| Note | Content |
|:--|:--|
| [`notes/desqueezing-relaxation-vacuum-fluctuations-note.md`](notes/desqueezing-relaxation-vacuum-fluctuations-note.md) | Main desqueezing note |
| [`notes/mapping_tables.md`](notes/mapping_tables.md) | Hand-edited mapping tables (the raw script output is `results/cosmological_mapping/mapping_tables_script.md`) |
| [`notes/data-pack-option0-internet.md`](notes/data-pack-option0-internet.md) | Literature data pack for the Option 0 slip discussion |

## Running

```bash
pip install -r requirements.txt      # includes qutip
python exploratory/scripts/desqueezing_relax_time.py
python exploratory/scripts/cosmological_mapping.py
python exploratory/scripts/sdiff_discrimination.py
python exploratory/scripts/euclid_protocol_forecasts.py
python exploratory/scripts/termico_vs_squeezed_recomendado.py
python exploratory/scripts/desqueezing_heavy_scan.py      # slowest
```

The `.npz` series these scripts write are not versioned (see `.gitignore`).
