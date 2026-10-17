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
| `scripts/desqueezing_heavy_scan.py` | `results/desqueezing_heavy/` | Wider \((\gamma, n_{\rm th}, r)\) grid checking \(t_{1/2}(|\langle a^2\rangle|)=\ln 2/\gamma\) (\(n_{\rm th}\) sets the stationary state, not the rate) |
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

## Correction log

**2026-10-17 — `desqueezing_heavy_scan.py`: reference half-life corrected** (same correction as block B1b in `desqueezing_relax_time.py`). The script compared the QuTiP half-life of \(|\langle a^2\rangle|\) with \(\ln 2/[\gamma(1+2n_{\rm th})]\) and printed "Law … holds" unconditionally; it also wrote the absolute local path of the CSV. Now the reference is \(t_{1/2}=\ln 2/\gamma\), "Law holds" is printed only if the median relative error is below 5 % ("Law does NOT hold" otherwise), and paths are relative to the repository. The simulated half-lives themselves are unchanged (the `t_half` column is identical).

| Output | Before | After |
|:--|:--|:--|
| `desqueezing_heavy_grid.csv`, `t_half_theory` | \(\ln 2/[\gamma(1+2n_{\rm th})]\) | \(\ln 2/\gamma\) |
| `desqueezing_heavy_grid.csv`, `rel_err` median / max | 0.600 / 2.003 | \(1.31\times 10^{-4}\) / \(8.98\times 10^{-4}\) |
| `desqueezing_heavy_summary.txt` | median/max \(1.31\times 10^{-4}\) / \(8.98\times 10^{-4}\) and "t_1/2 = ln2/gamma", written by hand: it did not match the CSV or the script, which gave 0.60 / 2.0 | same numbers, now produced by the script: "Law t_1/2 = ln2/γ holds across grid (median rel err < 5%)" |
| `desqueezing_heavy_summary.json` | `law_measured`, `note`, `*_vs_ln2_over_gamma` (hand-written) | `law` = "t_1/2 = ln2 / gamma", `law_tolerance_median_rel_err` = 0.05, `law_holds` = true, `csv` relative to the repository |
| `desqueezing_heavy.png` | dashed reference \(\ln 2/[\gamma(1+2n_{\rm th})]\) | dashed reference \(\ln 2/\gamma\) |
