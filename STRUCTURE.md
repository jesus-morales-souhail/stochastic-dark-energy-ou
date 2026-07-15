# Directory structure (canonical)

All SDE / DESI / desqueezing work for this project lives under:

`/home/ashpokemon/proyectos/stochastic-dark-energy-ou`

| Directory | Role | Git? |
|-----------|------|------|
| `papers/` | English scientific markdown + resume | Yes (core) |
| `scripts/` | Runnable analysis code | Yes (core) |
| `scripts/desqueezing/` | QuTiP + mapping | Yes |
| `scripts/gpe/` | GPE / Bogoliubov | Yes |
| `scripts/legacy/` | Older variants | Optional |
| `figures/` | PNG | Yes |
| `results/` | CSV / NPZ / txt outputs | Yes (selected) |
| `notes/` | Desqueezing synthesis | Yes |
| `outreach/` | Spanish accessible text | Optional |
| `literature/` | Reference PDFs | Optional (size) |
| `drafts/` | DOCX history + appendix draft | Local / optional |
| `data/` | Local DESI caches | **Local only** (gitignored bulk) |
| `local_archive/` | Everything moved from scattered PC locations | **Local only** (gitignored) |

## Sources consolidated (moved, not deleted)

- Windows Desktop `Nueva carpeta` → `literature/` + `drafts/docx/`
- Desktop `ou_bao_likelihood.py` → `local_archive/from_windows_desktop/`
- WSL `simulaciones_cuanticas/` → `notes/`, `scripts/desqueezing/`, `results/`
- Local `paper_Bogoliubov/`, `paper_OU/` → papers/scripts + archive remainder
- Reorg tools → `local_archive/reorg_tools/`

## Not part of this repo

- `~/proyectos/enzimas_microplasticos`
- `~/proyectos/especulacion_borradores` (geology, etc.)
- `~/proyectos/proyecto_unificacion` (PETase+DESI bridge — separate)
- Game / GTA folders on Desktop
