# Comprobaciones de la reorganización (2026-10-17)

Plan: [`REORGANIZACION.md`](REORGANIZACION.md). Movimientos en el commit «reorganización 1/2» (81 `git mv`, 22 `git rm`); rutas, cabeceras y READMEs en «reorganización 2/2». `figures/` conserva `profile_sigma_x.png` y `joint_w0wa_sigma.png` (decisión del autor), así que se quitaron 22 de las 24 copias, no 24.

Historial conservado (`git log --follow --oneline <fichero> | wc -l`): `exploratory/scripts/desqueezing_relax_time.py` 5 commits (hasta `33eb67b`), `theory/unimodular-gravity-vacuum-smoothness.md` 26 (hasta `c7aa79d` «Add files via upload»), `archive/scripts/ou_bao_stochastic_test.py` 17 (hasta `ab1d4bf`).

## (i) Enlaces relativos en Markdown

Script: `python analisis_cmb/check_links.py` (todos los `.md` versionados o nuevos; ignora http/mailto/anclas, bloques de código y código en línea).

Ficheros .md revisados: 64; enlaces relativos: 49; rotos: 0 (más 1 en snapshots históricos)

Rotos en snapshots históricos (no se editan):

| Fichero | Línea | Enlace |
|:--|--:|:--|
| `analisis_cmb/bloque_B/antes/results/amplification_routes/VERDICT.md` | 7 | `../../papers/amplification-gap.md` |

| Fichero | Enlaces relativos |
|:--|--:|
| `README.md` | 12 |
| `analisis_cmb/README.md` | 2 |
| `analisis_cmb/REORGANIZACION.md` | 1 |
| `analisis_cmb/REORGANIZACION_CHECK.md` | 1 |
| `analisis_cmb/bloque_B/antes/results/amplification_routes/VERDICT.md` | 1 |
| `exploratory/README.md` | 5 |
| `manuscript/CLAIMS.md` | 4 |
| `manuscript/PREPRINT.md` | 4 |
| `manuscript/README.md` | 3 |
| `papers/HONEST_HEADLINES.md` | 2 |
| `papers/amplification-gap.md` | 1 |
| `results/amplification_routes/VERDICT.md` | 1 |
| `theory/EXPLORATORY_BOUNDARY.md` | 2 |
| `theory/README.md` | 10 |

El enlace roto del snapshot ya lo estaba en `d71c4d2`: `analisis_cmb/bloque_B/antes/results/amplification_routes/VERDICT.md` es una copia de `results/amplification_routes/VERDICT.md` y su `../../papers/` resuelve a `analisis_cmb/bloque_B/antes/papers/`, que nunca existió. No se edita (snapshot histórico).

## (ii) Scripts del claim set: los resultados no cambian

Cada script se ejecutó desde la raíz del repo; después, `git status --porcelain` sobre `results/` y `figures/` no muestra ningún cambio. Los scripts reescribieron sus salidas (p. ej. `results/profile_sigma_x/*` con hora 04:39) y salieron idénticas byte a byte.

| Script | Código de salida | Tiempo |
|:--|:--:|--:|
| `scripts/verify_bao_pipelines.py` | 0 | 138 s |
| `scripts/profile_sigma_x_desi.py` | 0 | 1 s |
| `scripts/gap_two_targets.py` | 0 | 0 s |
| `scripts/slip_bridge.py` | 0 | 0 s |
| `scripts/amplifier_audit.py` | 0 | 1 s |
| `scripts/amplification/route1_local_causal_set_seed.py` | 0 | 0 s |

## (iii) Scripts de `exploratory/scripts/`: escriben dentro de `exploratory/results/`

Ejecutados con el directorio de trabajo en `/tmp`, para comprobar que no dependen del cwd. Ninguno tiene opción rápida, así que se ejecutaron completos. Después: ningún fichero nuevo fuera de `exploratory/results/`, nada en `figures/` (las copias a `figures/` se quitaron de los scripts) y no existe `scripts/results/`.

| Script | Código de salida | Tiempo | Salida | Diferencia con lo versionado |
|:--|:--:|--:|:--|:--|
| `termico_vs_squeezed_recomendado.py` | 0 | 10 s | `exploratory/results/termico_vs_squeezed_recomendado/` | resumen y PNG distintos (ver nota 2); no se versionan |
| `cosmological_mapping.py` | 0 | 0 s | `exploratory/results/cosmological_mapping/` | ninguna |
| `euclid_protocol_forecasts.py` | 0 | 1 s | `exploratory/results/euclid_protocol/` | ninguna |
| `sdiff_discrimination.py` | 0 | 1 s | `exploratory/results/sdiff_discrimination/` | ninguna |
| `desqueezing_relax_time.py` | 0 | 36 s | `exploratory/results/desqueezing_relax_time/` | solo las 3 rutas del resumen (`results/…` → `exploratory/results/…`); se versiona |
| `desqueezing_heavy_scan.py` | 0 | 135 s | `exploratory/results/desqueezing_heavy/` | resumen, JSON y PNG distintos (ver nota 1, ya resuelta) |

Notas (problemas anteriores a la reorganización):

1. ~~`desqueezing_heavy_scan.py` no reproducía su salida versionada~~ — **resuelto** el 2026-10-17 en el commit «exploratory: fórmula de t½ corregida en heavy_scan, salida regenerada» (antes/después en `exploratory/README.md`, sección *Correction log*).
2. **`termico_vs_squeezed_recomendado_resumen.txt` versionado** es un resumen corto en inglés escrito a mano; el script escribe uno largo en español. Los números coinciden (Var(0) con error 0.0497 %, termalización a Var(q) = 0.6 y ⟨n⟩ = 0.1). Se restauró la versión versionada.

`scripts/gpe/gpe_sim.py` (nuevo `savefig`): ejecutado desde `/tmp`, código 0, escribe en `results/gpe/` y nada en el cwd; los dos PNG salen idénticos a los versionados.

Ficheros ignorados (`.npz`/`.npy`, no versionados) que quedaron en las carpetas viejas de `results/` tras el `git mv`: las cadenas Euclid (3.2 MB) se movieron a `archive/results/euclid_mcmc/` y `archive/results/euclid_joint_mcmc/`, donde siguen ignoradas por las reglas nuevas del `.gitignore`. Las dos series antiguas de desqueezing se borraron porque (iii) las regeneró en `exploratory/results/`.

## (iv) Todo lo que queda en `scripts/` se importa

`python -c "import <módulo>"` con `PYTHONPATH` = carpeta del script + `scripts/`, desde `/tmp`. Cinco scripts no tienen `if __name__ == "__main__"` y al importarlos se ejecutarían: `descargarcovarianza.py` descarga de Zenodo, y los otros cuatro hacen el cálculo y escriben ficheros. Esos cinco se comprueban con `compile()` más el import de todas sus dependencias. Tras la comprobación, `git status` no cambia.

| Script | Comprobación | Resultado |
|:--|:--|:--|
| `scripts/amplification/route1_local_causal_set_seed.py` | import | OK |
| `scripts/amplification/route2_late_horizon_exit.py` | import | OK |
| `scripts/amplification/route3_nonlinear_avalanche.py` | import | OK |
| `scripts/amplification/run_all_routes.py` | import | OK |
| `scripts/amplifier_audit.py` | import | OK |
| `scripts/descargarcovarianza.py` | compile + import de sus dependencias (astropy.io, json, numpy, os, pathlib, requests, tarfile, zipfile) | OK |
| `scripts/desi_bound.py` | import | OK |
| `scripts/desi_dr2_data.py` | import | OK |
| `scripts/desi_dr2_real_bao_test.py` | import | OK |
| `scripts/eos_efectiva.py` | import | OK |
| `scripts/extraer_cov_cobaya.py` | compile + import de sus dependencias (numpy, os, pathlib) | OK |
| `scripts/gap_two_targets.py` | import | OK |
| `scripts/gpe/Bogoliubov_DESI_extension.py` | compile + import de sus dependencias (numpy) | OK |
| `scripts/gpe/gpe_sim.py` | compile + import de sus dependencias (matplotlib.animation, matplotlib.pyplot, numpy, pathlib) | OK |
| `scripts/joint_w0wa_sigma_desi.py` | import | OK |
| `scripts/ou_bao_likelihood.py` | compile + import de sus dependencias (desi_dr2_data, matplotlib.patches, matplotlib.pyplot, numpy, os, pathlib, scipy.integrate, scipy.linalg, scipy.optimize, scipy.stats, sys, warnings) | OK |
| `scripts/profile_sigma_x_desi.py` | import | OK |
| `scripts/slip_bridge.py` | import | OK |
| `scripts/tachyonic_rank1_mle.py` | import | OK |
| `scripts/verify_bao_pipelines.py` | import | OK |
