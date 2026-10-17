# Reorganización del repo — paso 1: inventario y destinos propuestos

Estado: **solo inventario**. No se ha movido ni borrado nada.

**Ejecutado el 2026-10-17** en los commits «reorganización 1/2» y «2/2», con las decisiones del autor sobre los ⚠: `figures/` conserva `profile_sigma_x.png` y `joint_w0wa_sigma.png` (se quitaron 22 copias, no 24); `cross_correlation_DESI.py` → `archive/scripts/`, con las rutas corregidas; `resume.txt` sigue en `papers/`; `gpe_*.png` → `results/gpe/`. Comprobaciones en [`REORGANIZACION_CHECK.md`](REORGANIZACION_CHECK.md). Las rutas de las tablas de abajo son las de antes de mover. Base: `git ls-files` en `d71c4d2` (472 ficheros versionados).

Leídos antes de clasificar: `README.md`, `manuscript/README.md`, `manuscript/CLAIMS.md` (N5, N6, N9) y `papers/EXPLORATORY_BOUNDARY.md`.

## Resumen

| Destino | Ficheros |
|:--|--:|
| Sin cambios (raíz, `manuscript/`, `papers/`, `scripts/`, `results/`, `analisis_cmb/`) | 367 |
| `theory/` | 8 |
| `exploratory/` (`scripts/`, `results/`, `notes/`) | 32 |
| `archive/` | 39 |
| `results/gpe/` (figuras que hoy están en `figures/`) | 2 |
| Eliminar de `figures/` (copias idénticas; el original sigue en `results/`) | 24 |

`figures/` queda **vacía**: `manuscript/PREPRINT.md` no cita ninguna figura (`grep -n 'png\|figures/\|Figure'` solo encuentra la frase genérica de la l.277, «…slip bridge, figures, and technical notes…»). Habrá que decidir si se borra la carpeta o si el preprint pasa a citar figuras (por ejemplo `results/profile_sigma_x/profile_sigma_x.png`).

## Puntos que necesitan decisión (⚠)

1. **`scripts/cross_correlation_DESI.py` → `archive/`**, pero el PREPRINT lo cita (§6, l.265, como test exploratorio y preliminar, de acuerdo con N5), y también `papers/stochastic-dark-energy-desi-dr2.md`. Al moverlo hay que cambiar esas dos rutas a `archive/scripts/…`; la frase del PREPRINT ya lo presenta como fuera del claim set, así que el contenido no cambia.
2. **`papers/resume.txt`** no estaba en ninguna lista. Es el resumen numérico de C1–C7 (lo enlaza el README como «Numbers snapshot») → propongo `papers/`.
3. **`figures/gpe_negative_mass_{evolution,instability}.png`**: no son copias ni figuras huérfanas. Las genera `scripts/gpe/gpe_sim.py`, que escribe en el directorio actual. Propongo `results/gpe/` y cambiar el `savefig` a esa ruta. La alternativa es `archive/figures/`.
4. **«Las 14 figuras» de Euclid en `figures/`**: en realidad hay **17** (`euclid_mcmc_*` ×10, `joint_mcmc_*` ×7). 14 son copias idénticas (cmp) de `results/euclid_mcmc/` y `results/euclid_joint_mcmc/`, que van enteros a `archive/`, así que se eliminan de `figures/` y no se duplican. Las otras 3 (`joint_mcmc_{E1_amplitude,E2_mean_reversion,null_F_E0}_BAO.png`) **difieren** de las de `results/euclid_joint_mcmc/` (cmp: byte 127). Son otra versión, de 470619a, y van a `archive/figures/`.
5. **Snapshots y logs de `analisis_cmb/bloque_A|B/`**: son registro histórico y se quedan como están. Sus rutas antiguas **no** se corrigen (se excluyen de la columna de referencias).

## `ou_bao_likelihood.py` vs `ou_bao_stochastic_test.py`

`diff` (591 vs 551 líneas; el resto es idéntico):

| | `ou_bao_likelihood.py` | `ou_bao_stochastic_test.py` |
|:--|:--|:--|
| Datos | DESI **DR2** (arXiv:2503.14738) mediante `desi_dr2_data.load_alpha_dv(prefer_file=True, use_full_cov=True)` | DESI **DR1** (arXiv:2404.03000, Tabla 3) con `z_eff`, `alpha`, `sigma` y `S_z` escritos a mano en el código |
| Covarianza | `C_MEAS` = Gaussian BAO oficial proyectada a α (completa) | diagonal `sigma²` |
| Import | `sys.path` → `scripts/`; `from desi_dr2_data import …` | ninguno del repo |
| Figura | `plots/test_desi_QNM.png` (ignorado por `.gitignore`) | `figures/test_desi_QNM.png` |
| Extra | define `generate_exclusion_plot()` (no se llama nunca) | — |
| Texto | «Consistency check» | «Falsification criterion F6», «Version 2.1» |

El **PREPRINT cita `scripts/ou_bao_likelihood.py`** (tabla l.125 y reproducción l.308). Por tanto `ou_bao_stochastic_test.py` → `archive/scripts/`.

## Figuras de `figures/` que ningún script genera (→ `archive/figures/`)

| Figura | Quién la creó (git log --follow) |
|:--|:--|
| `Figure_1.png`, `Figure_2.png` | `500f05b` 2026-07-13 «Add files via upload» (en la raíz); movidas a `figures/` en `33eb67b` |
| `falsos_positivos_detectores.png` | `826dd3e` 2026-07-06 «Add files via upload»; movida en `33eb67b` |
| `exclusion_plot.png` | `f247fcb` 2026-07-05 «Add files via upload»; movida en `33eb67b`. El código que la dibujaría (`generate_exclusion_plot` en `ou_bao_likelihood.py`) no se llama nunca y escribe en `plots/` |
| `test_desi_QNM.png` | añadida en `33eb67b` 2026-07-15 («Restructure repository…»); hoy solo la escribe en `figures/` `ou_bao_stochastic_test.py` (DR1), que va a `archive/` |
| `joint_mcmc_{E1_amplitude,E2_mean_reversion,null_F_E0}_BAO.png` | `470619a` 2026-07-15 «Add joint Euclid BAO + SN mock MCMC…»; el script que las generaba es hoy un stub (`euclid_joint_bao_sne_mcmc.py`) |

Las `gpe_*` (subidas en `500f05b`) sí tienen generador: ver ⚠3.

## Tabla de ficheros

Columna de referencias: `git grep -l -F` con la ruta completa, la ruta sin el primer directorio y el nombre del fichero (este último solo si no es genérico y aparece como mucho 3 veces en el repo). Para `.py` se busca además `import <módulo>` / `from <módulo> import`. Se excluyen `analisis_cmb/bloque_*/logs/` y `analisis_cmb/bloque_*/antes/` (registro histórico). Solo se rellena para los ficheros que cambian de sitio; `=` = se queda donde está, `—` = no aplica. Las referencias a **carpetas** enteras (p. ej. `results/euclid_mcmc/`) están en el apartado (b).

Hay algunos falsos positivos por nombre de fichero, y se dejan a la vista: `figures/test_desi_QNM.png` ← `ou_bao_likelihood.py` (escribe `plots/test_desi_QNM.png`, no esta), `figures/exclusion_plot.png` ← `ou_bao_likelihood.py` (ídem).

| Ruta actual | Destino | Motivo | Referencias que habrá que corregir |
|:--|:--|:--|:--|
| `.gitignore` | = | raíz | — |
| `LICENSE` | = | raíz | — |
| `README.md` | = | raíz | — |
| `analisis_cmb/AUDITORIA.md` | = | análisis CMB/Planck | — |
| `analisis_cmb/NOTAS_ANALISIS_CMB.md` | = | análisis CMB/Planck | — |
| `analisis_cmb/README.md` | = | análisis CMB/Planck | — |
| `analisis_cmb/RESULTADO.md` | = | análisis CMB/Planck | — |
| `analisis_cmb/bloque_A/COMPARACION.md` | = | análisis CMB/Planck | — |
| `analisis_cmb/bloque_A/antes/figures/joint_w0wa_sigma.png` | = | snapshot histórico (bloque A/B) | — |
| `analisis_cmb/bloque_A/antes/results/eos_cpl_desi_dr2/eos_cpl_summary.json` | = | snapshot histórico (bloque A/B) | — |
| `analisis_cmb/bloque_A/antes/results/eos_cpl_desi_dr2/eos_cpl_summary.txt` | = | snapshot histórico (bloque A/B) | — |
| `analisis_cmb/bloque_A/antes/results/joint_w0wa_sigma/joint_w0wa_sigma.json` | = | snapshot histórico (bloque A/B) | — |
| `analisis_cmb/bloque_A/antes/results/joint_w0wa_sigma/joint_w0wa_sigma.png` | = | snapshot histórico (bloque A/B) | — |
| `analisis_cmb/bloque_A/antes/results/joint_w0wa_sigma/joint_w0wa_sigma.txt` | = | snapshot histórico (bloque A/B) | — |
| `analisis_cmb/bloque_A/logs/eos_efectiva.log` | = | log histórico (bloque A/B) | — |
| `analisis_cmb/bloque_A/logs/joint_w0wa_sigma_desi.log` | = | log histórico (bloque A/B) | — |
| `analisis_cmb/bloque_A/logs/verify_bao_pipelines.log` | = | log histórico (bloque A/B) | — |
| `analisis_cmb/bloque_B/COMPARACION_B1.md` | = | análisis CMB/Planck | — |
| `analisis_cmb/bloque_B/COMPARACION_B2.md` | = | análisis CMB/Planck | — |
| `analisis_cmb/bloque_B/antes/analisis_cmb/verificacion_numeros.md` | = | snapshot histórico (bloque A/B) | — |
| `analisis_cmb/bloque_B/antes/figures/desqueezing_relax_time.png` | = | snapshot histórico (bloque A/B) | — |
| `analisis_cmb/bloque_B/antes/figures/euclid_A0_theta_plane.png` | = | snapshot histórico (bloque A/B) | — |
| `analisis_cmb/bloque_B/antes/figures/euclid_theta_shape_detectability.png` | = | snapshot histórico (bloque A/B) | — |
| `analisis_cmb/bloque_B/antes/figures/profile_sigma_x.png` | = | snapshot histórico (bloque A/B) | — |
| `analisis_cmb/bloque_B/antes/figures/sdiff_A0_min_for_detection.png` | = | snapshot histórico (bloque A/B) | — |
| `analisis_cmb/bloque_B/antes/figures/sdiff_residual_vs_theta.png` | = | snapshot histórico (bloque A/B) | — |
| `analisis_cmb/bloque_B/antes/figures/sdiff_t_half_vs_theta.png` | = | snapshot histórico (bloque A/B) | — |
| `analisis_cmb/bloque_B/antes/notes/mapping_tables.md` | = | snapshot histórico (bloque A/B) | — |
| `analisis_cmb/bloque_B/antes/results/amplification_routes/VERDICT.md` | = | snapshot histórico (bloque A/B) | — |
| `analisis_cmb/bloque_B/antes/results/amplification_routes/route1_N_eff_grid.csv` | = | snapshot histórico (bloque A/B) | — |
| `analisis_cmb/bloque_B/antes/results/amplification_routes/route1_N_eff_required.csv` | = | snapshot histórico (bloque A/B) | — |
| `analisis_cmb/bloque_B/antes/results/amplification_routes/route2_horizon_exit_scan.csv` | = | snapshot histórico (bloque A/B) | — |
| `analisis_cmb/bloque_B/antes/results/amplification_routes/route3_avalanche_scan.csv` | = | snapshot histórico (bloque A/B) | — |
| `analisis_cmb/bloque_B/antes/results/cosmological_mapping/mapping_summary.txt` | = | snapshot histórico (bloque A/B) | — |
| `analisis_cmb/bloque_B/antes/results/cosmological_mapping/table1_theta.csv` | = | snapshot histórico (bloque A/B) | — |
| `analisis_cmb/bloque_B/antes/results/cosmological_mapping/table2_A0_required.csv` | = | snapshot histórico (bloque A/B) | — |
| `analisis_cmb/bloque_B/antes/results/cosmological_mapping/table3_page_scrambling.csv` | = | snapshot histórico (bloque A/B) | — |
| `analisis_cmb/bloque_B/antes/results/cosmological_mapping/table4_DESI_S_z.csv` | = | snapshot histórico (bloque A/B) | — |
| `analisis_cmb/bloque_B/antes/results/desi_dr2_real_bao/desi_dr2_real_bao.json` | = | snapshot histórico (bloque A/B) | — |
| `analisis_cmb/bloque_B/antes/results/desi_dr2_real_bao/desi_dr2_real_bao.txt` | = | snapshot histórico (bloque A/B) | — |
| `analisis_cmb/bloque_B/antes/results/desqueezing_relax_time/desqueezing_relax_time.png` | = | snapshot histórico (bloque A/B) | — |
| `analisis_cmb/bloque_B/antes/results/desqueezing_relax_time/desqueezing_relax_time_resumen.txt` | = | snapshot histórico (bloque A/B) | — |
| `analisis_cmb/bloque_B/antes/results/desqueezing_relax_time/t_relax_series.npz` | = | snapshot histórico (bloque A/B) | — |
| `analisis_cmb/bloque_B/antes/results/desqueezing_relax_time/t_relax_tabla.csv` | = | snapshot histórico (bloque A/B) | — |
| `analisis_cmb/bloque_B/antes/results/euclid_protocol/euclid_A0_theta_plane.png` | = | snapshot histórico (bloque A/B) | — |
| `analisis_cmb/bloque_B/antes/results/euclid_protocol/euclid_forecast_grid.csv` | = | snapshot histórico (bloque A/B) | — |
| `analisis_cmb/bloque_B/antes/results/euclid_protocol/theta_shape_detectability.png` | = | snapshot histórico (bloque A/B) | — |
| `analisis_cmb/bloque_B/antes/results/ou_bao_desi_dr2_run.txt` | = | snapshot histórico (bloque A/B) | — |
| `analisis_cmb/bloque_B/antes/results/profile_sigma_x/profile_sigma_x.json` | = | snapshot histórico (bloque A/B) | — |
| `analisis_cmb/bloque_B/antes/results/profile_sigma_x/profile_sigma_x.png` | = | snapshot histórico (bloque A/B) | — |
| `analisis_cmb/bloque_B/antes/results/profile_sigma_x/profile_sigma_x.txt` | = | snapshot histórico (bloque A/B) | — |
| `analisis_cmb/bloque_B/antes/results/sdiff_discrimination/A0_min_for_detection.png` | = | snapshot histórico (bloque A/B) | — |
| `analisis_cmb/bloque_B/antes/results/sdiff_discrimination/residual_vs_theta.png` | = | snapshot histórico (bloque A/B) | — |
| `analisis_cmb/bloque_B/antes/results/sdiff_discrimination/sdiff_discrimination_grid.csv` | = | snapshot histórico (bloque A/B) | — |
| `analisis_cmb/bloque_B/antes/results/sdiff_discrimination/summary.txt` | = | snapshot histórico (bloque A/B) | — |
| `analisis_cmb/bloque_B/antes/results/sdiff_discrimination/t_half_vs_theta.png` | = | snapshot histórico (bloque A/B) | — |
| `analisis_cmb/bloque_B/logs/amplifier_audit.log` | = | log histórico (bloque A/B) | — |
| `analisis_cmb/bloque_B/logs/cosmological_mapping.log` | = | log histórico (bloque A/B) | — |
| `analisis_cmb/bloque_B/logs/cosmological_mapping_B1b.log` | = | log histórico (bloque A/B) | — |
| `analisis_cmb/bloque_B/logs/desi_dr2_real_bao_test.log` | = | log histórico (bloque A/B) | — |
| `analisis_cmb/bloque_B/logs/desqueezing_relax_time.log` | = | log histórico (bloque A/B) | — |
| `analisis_cmb/bloque_B/logs/desqueezing_relax_time_B1b.log` | = | log histórico (bloque A/B) | — |
| `analisis_cmb/bloque_B/logs/estadisticos_sigma_2.633e-02.log` | = | log histórico (bloque A/B) | — |
| `analisis_cmb/bloque_B/logs/euclid_protocol_forecasts.log` | = | log histórico (bloque A/B) | — |
| `analisis_cmb/bloque_B/logs/gap_two_targets.log` | = | log histórico (bloque A/B) | — |
| `analisis_cmb/bloque_B/logs/profile_sigma_x_desi.log` | = | log histórico (bloque A/B) | — |
| `analisis_cmb/bloque_B/logs/realizaciones_sigma_2.633e-02.log` | = | log histórico (bloque A/B) | — |
| `analisis_cmb/bloque_B/logs/realizaciones_sigma_2.633e-02_refs_json.log` | = | log histórico (bloque A/B) | — |
| `analisis_cmb/bloque_B/logs/route1_local_causal_set_seed.log` | = | log histórico (bloque A/B) | — |
| `analisis_cmb/bloque_B/logs/sdiff_discrimination.log` | = | log histórico (bloque A/B) | — |
| `analisis_cmb/bloque_B/logs/slip_bridge.log` | = | log histórico (bloque A/B) | — |
| `analisis_cmb/bloque_B/logs/verificacion_numeros.log` | = | log histórico (bloque A/B) | — |
| `analisis_cmb/diagnosticos/banda_delta.py` | = | análisis CMB/Planck | — |
| `analisis_cmb/diagnosticos/conv_fina.py` | = | análisis CMB/Planck | — |
| `analisis_cmb/diagnosticos/conv_tabla.py` | = | análisis CMB/Planck | — |
| `analisis_cmb/diagnosticos/descomposicion_dchi2.py` | = | análisis CMB/Planck | — |
| `analisis_cmb/diagnosticos/obsoleto/realizaciones_H0fijo_mead2020.py` | `archive/analisis_cmb/diagnosticos/obsoleto/realizaciones_H0fijo_mead2020.py` | diagnóstico H0-fijo/mead2020 obsoleto, ya marcado como tal | `analisis_cmb/NOTAS_ANALISIS_CMB.md` |
| `analisis_cmb/diagnosticos/obsoleto/tabla_H0fijo.py` | `archive/analisis_cmb/diagnosticos/obsoleto/tabla_H0fijo.py` | diagnóstico H0-fijo/mead2020 obsoleto, ya marcado como tal | ninguna |
| `analisis_cmb/estadisticos.py` | = | análisis CMB/Planck | — |
| `analisis_cmb/everpresent/.gitignore` | = | everpresent vs DESI | — |
| `analisis_cmb/everpresent/EVERPRESENT_vs_DESI.md` | = | everpresent vs DESI | — |
| `analisis_cmb/everpresent/LECTURA_EVERPRESENT.md` | = | everpresent vs DESI | — |
| `analisis_cmb/everpresent/N_COUNT.md` | = | everpresent vs DESI | — |
| `analisis_cmb/everpresent/ahmed_sim.py` | = | everpresent vs DESI | — |
| `analisis_cmb/everpresent/bao_sim.py` | = | everpresent vs DESI | — |
| `analisis_cmb/everpresent/desi_dr2_data.py` | = | everpresent vs DESI | — |
| `analisis_cmb/everpresent/n_count.py` | = | everpresent vs DESI | — |
| `analisis_cmb/everpresent/n_count_salida.txt` | = | everpresent vs DESI | — |
| `analisis_cmb/everpresent/profile_desi_everpresent.py` | = | everpresent vs DESI | — |
| `analisis_cmb/everpresent/rd_planck.py` | = | everpresent vs DESI | — |
| `analisis_cmb/everpresent/rd_real.py` | = | everpresent vs DESI | — |
| `analisis_cmb/everpresent/rdcalc.py` | = | everpresent vs DESI | — |
| `analisis_cmb/everpresent/resultados/ahmed_sim.json` | = | everpresent vs DESI | — |
| `analisis_cmb/everpresent/resultados/ahmed_sim.txt` | = | everpresent vs DESI | — |
| `analisis_cmb/everpresent/resultados/bao_m1_alpha0.005.json` | = | everpresent vs DESI | — |
| `analisis_cmb/everpresent/resultados/bao_m1_alpha0.005.txt` | = | everpresent vs DESI | — |
| `analisis_cmb/everpresent/resultados/bao_m1_alpha0.01.json` | = | everpresent vs DESI | — |
| `analisis_cmb/everpresent/resultados/bao_m1_alpha0.01.npz` | = | everpresent vs DESI | — |
| `analisis_cmb/everpresent/resultados/bao_m1_alpha0.01.txt` | = | everpresent vs DESI | — |
| `analisis_cmb/everpresent/resultados/bao_m2.json` | = | everpresent vs DESI | — |
| `analisis_cmb/everpresent/resultados/bao_m2.txt` | = | everpresent vs DESI | — |
| `analisis_cmb/everpresent/resultados/bao_m2_alphagrid.npz` | = | everpresent vs DESI | — |
| `analisis_cmb/everpresent/resultados/bao_m2_bestfit.npz` | = | everpresent vs DESI | — |
| `analisis_cmb/everpresent/resultados/das_model1.json` | = | everpresent vs DESI | — |
| `analisis_cmb/everpresent/resultados/das_model1.txt` | = | everpresent vs DESI | — |
| `analisis_cmb/everpresent/resultados/profile_desi_everpresent.json` | = | everpresent vs DESI | — |
| `analisis_cmb/everpresent/resultados/profile_desi_everpresent.txt` | = | everpresent vs DESI | — |
| `analisis_cmb/everpresent/resultados/rd_planck.json` | = | everpresent vs DESI | — |
| `analisis_cmb/everpresent/resultados/rd_real.json` | = | everpresent vs DESI | — |
| `analisis_cmb/everpresent/resultados/rd_real_m1_alpha0.005.npz` | = | everpresent vs DESI | — |
| `analisis_cmb/everpresent/resultados/rd_real_m1_alpha0.01.npz` | = | everpresent vs DESI | — |
| `analisis_cmb/everpresent/resultados/rd_real_m2.npz` | = | everpresent vs DESI | — |
| `analisis_cmb/everpresent/resultados/robustez.json` | = | everpresent vs DESI | — |
| `analisis_cmb/everpresent/resultados/robustez.txt` | = | everpresent vs DESI | — |
| `analisis_cmb/everpresent/robustez.py` | = | everpresent vs DESI | — |
| `analisis_cmb/malla.py` | = | análisis CMB/Planck | — |
| `analisis_cmb/oscilador/README.md` | = | oscilador | — |
| `analisis_cmb/oscilador/oscilador.py` | = | oscilador | — |
| `analisis_cmb/oscilador/resultados/oscilador.csv` | = | oscilador | — |
| `analisis_cmb/oscilador/resultados/oscilador.json` | = | oscilador | — |
| `analisis_cmb/oscilador/resultados/tabla.md` | = | oscilador | — |
| `analisis_cmb/ou_cmb.py` | = | análisis CMB/Planck | — |
| `analisis_cmb/planck_datos.py` | = | análisis CMB/Planck | — |
| `analisis_cmb/realizacion_unica.py` | = | análisis CMB/Planck | — |
| `analisis_cmb/realizaciones.py` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/diagnosticos/conv_fina.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/diagnosticos/conv_tabla.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/log_sigma_2.506e-02.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/log_takahashi_2.506e-02.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/log_takahashi_malla.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/realizacion_unica/convergencia.json` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/realizacion_unica/convergencia.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/realizacion_unica/realizacion_theta1.png` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/sigma_2.506e-02/banda_delta.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/sigma_2.506e-02/dchi2_TT_theta0.001_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/sigma_2.506e-02/dchi2_TT_theta0.001_Delta0.10.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/sigma_2.506e-02/dchi2_TT_theta1_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/sigma_2.506e-02/dchi2_TT_theta1_Delta0.10.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/sigma_2.506e-02/dchi2_TT_theta20_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/sigma_2.506e-02/dchi2_TT_theta20_Delta0.10.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/sigma_2.506e-02/dchi2_pp_theta0.001_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/sigma_2.506e-02/dchi2_pp_theta0.001_Delta0.10.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/sigma_2.506e-02/dchi2_pp_theta1_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/sigma_2.506e-02/dchi2_pp_theta1_Delta0.10.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/sigma_2.506e-02/dchi2_pp_theta20_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/sigma_2.506e-02/dchi2_pp_theta20_Delta0.10.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/sigma_2.506e-02/dchi2_thetastar_theta0.001_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/sigma_2.506e-02/dchi2_thetastar_theta0.001_Delta0.10.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/sigma_2.506e-02/dchi2_thetastar_theta1_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/sigma_2.506e-02/dchi2_thetastar_theta1_Delta0.10.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/sigma_2.506e-02/dchi2_thetastar_theta20_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/sigma_2.506e-02/dchi2_thetastar_theta20_Delta0.10.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/sigma_2.506e-02/observables_theta0.001_Delta0.05.npz` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/sigma_2.506e-02/observables_theta0.001_Delta0.10.npz` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/sigma_2.506e-02/observables_theta1_Delta0.05.npz` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/sigma_2.506e-02/observables_theta1_Delta0.10.npz` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/sigma_2.506e-02/observables_theta20_Delta0.05.npz` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/sigma_2.506e-02/observables_theta20_Delta0.10.npz` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/sigma_2.506e-02/resumen.json` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/sigma_2.506e-02/tabla.md` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/malla.json` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/malla.md` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_1.000e-04/dchi2_TT_theta0.001_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_1.000e-04/dchi2_TT_theta1_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_1.000e-04/dchi2_pp_theta0.001_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_1.000e-04/dchi2_pp_theta1_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_1.000e-04/fallos_theta0.001_Delta0.05.json` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_1.000e-04/fallos_theta1_Delta0.05.json` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_1.000e-04/meta.json` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_1.000e-04/observables_theta0.001_Delta0.05.npz` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_1.000e-04/observables_theta1_Delta0.05.npz` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_1.000e-04/resumen.json` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_1.000e-04/tabla.md` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_1.165e-03/dchi2_TT_theta0.001_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_1.165e-03/dchi2_TT_theta1_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_1.165e-03/dchi2_pp_theta0.001_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_1.165e-03/dchi2_pp_theta1_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_1.165e-03/fallos_theta0.001_Delta0.05.json` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_1.165e-03/fallos_theta1_Delta0.05.json` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_1.165e-03/meta.json` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_1.165e-03/observables_theta0.001_Delta0.05.npz` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_1.165e-03/observables_theta1_Delta0.05.npz` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_1.165e-03/resumen.json` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_1.165e-03/tabla.md` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_1.357e-02/dchi2_TT_theta0.001_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_1.357e-02/dchi2_TT_theta1_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_1.357e-02/dchi2_pp_theta0.001_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_1.357e-02/dchi2_pp_theta1_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_1.357e-02/fallos_theta0.001_Delta0.05.json` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_1.357e-02/fallos_theta1_Delta0.05.json` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_1.357e-02/meta.json` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_1.357e-02/observables_theta0.001_Delta0.05.npz` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_1.357e-02/observables_theta1_Delta0.05.npz` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_1.357e-02/resumen.json` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_1.357e-02/tabla.md` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_1.847e-04/dchi2_TT_theta0.001_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_1.847e-04/dchi2_TT_theta1_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_1.847e-04/dchi2_pp_theta0.001_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_1.847e-04/dchi2_pp_theta1_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_1.847e-04/fallos_theta0.001_Delta0.05.json` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_1.847e-04/fallos_theta1_Delta0.05.json` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_1.847e-04/meta.json` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_1.847e-04/observables_theta0.001_Delta0.05.npz` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_1.847e-04/observables_theta1_Delta0.05.npz` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_1.847e-04/resumen.json` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_1.847e-04/tabla.md` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_2.152e-03/dchi2_TT_theta0.001_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_2.152e-03/dchi2_TT_theta1_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_2.152e-03/dchi2_pp_theta0.001_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_2.152e-03/dchi2_pp_theta1_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_2.152e-03/fallos_theta0.001_Delta0.05.json` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_2.152e-03/fallos_theta1_Delta0.05.json` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_2.152e-03/meta.json` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_2.152e-03/observables_theta0.001_Delta0.05.npz` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_2.152e-03/observables_theta1_Delta0.05.npz` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_2.152e-03/resumen.json` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_2.152e-03/tabla.md` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_2.506e-02/dchi2_TT_theta0.001_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_2.506e-02/dchi2_TT_theta1_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_2.506e-02/dchi2_pp_theta0.001_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_2.506e-02/dchi2_pp_theta1_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_2.506e-02/fallos_theta0.001_Delta0.05.json` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_2.506e-02/fallos_theta1_Delta0.05.json` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_2.506e-02/meta.json` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_2.506e-02/observables_theta0.001_Delta0.05.npz` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_2.506e-02/observables_theta1_Delta0.05.npz` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_2.506e-02/resumen.json` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_2.506e-02/tabla.md` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_3.413e-04/dchi2_TT_theta0.001_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_3.413e-04/dchi2_TT_theta1_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_3.413e-04/dchi2_pp_theta0.001_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_3.413e-04/dchi2_pp_theta1_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_3.413e-04/fallos_theta0.001_Delta0.05.json` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_3.413e-04/fallos_theta1_Delta0.05.json` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_3.413e-04/meta.json` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_3.413e-04/observables_theta0.001_Delta0.05.npz` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_3.413e-04/observables_theta1_Delta0.05.npz` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_3.413e-04/resumen.json` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_3.413e-04/tabla.md` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_3.975e-03/dchi2_TT_theta0.001_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_3.975e-03/dchi2_TT_theta1_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_3.975e-03/dchi2_pp_theta0.001_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_3.975e-03/dchi2_pp_theta1_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_3.975e-03/fallos_theta0.001_Delta0.05.json` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_3.975e-03/fallos_theta1_Delta0.05.json` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_3.975e-03/meta.json` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_3.975e-03/observables_theta0.001_Delta0.05.npz` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_3.975e-03/observables_theta1_Delta0.05.npz` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_3.975e-03/resumen.json` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_3.975e-03/tabla.md` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_6.305e-04/dchi2_TT_theta0.001_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_6.305e-04/dchi2_TT_theta1_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_6.305e-04/dchi2_pp_theta0.001_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_6.305e-04/dchi2_pp_theta1_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_6.305e-04/fallos_theta0.001_Delta0.05.json` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_6.305e-04/fallos_theta1_Delta0.05.json` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_6.305e-04/meta.json` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_6.305e-04/observables_theta0.001_Delta0.05.npz` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_6.305e-04/observables_theta1_Delta0.05.npz` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_6.305e-04/resumen.json` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_6.305e-04/tabla.md` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_7.344e-03/dchi2_TT_theta0.001_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_7.344e-03/dchi2_TT_theta1_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_7.344e-03/dchi2_pp_theta0.001_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_7.344e-03/dchi2_pp_theta1_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_7.344e-03/fallos_theta0.001_Delta0.05.json` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_7.344e-03/fallos_theta1_Delta0.05.json` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_7.344e-03/meta.json` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_7.344e-03/observables_theta0.001_Delta0.05.npz` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_7.344e-03/observables_theta1_Delta0.05.npz` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_7.344e-03/resumen.json` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/malla/sigma_7.344e-03/tabla.md` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/sigma_2.506e-02/dchi2_TT_theta0.001_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/sigma_2.506e-02/dchi2_TT_theta0.001_Delta0.10.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/sigma_2.506e-02/dchi2_TT_theta1_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/sigma_2.506e-02/dchi2_TT_theta1_Delta0.10.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/sigma_2.506e-02/dchi2_TT_theta20_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/sigma_2.506e-02/dchi2_TT_theta20_Delta0.10.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/sigma_2.506e-02/dchi2_pp_theta0.001_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/sigma_2.506e-02/dchi2_pp_theta0.001_Delta0.10.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/sigma_2.506e-02/dchi2_pp_theta1_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/sigma_2.506e-02/dchi2_pp_theta1_Delta0.10.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/sigma_2.506e-02/dchi2_pp_theta20_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/sigma_2.506e-02/dchi2_pp_theta20_Delta0.10.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/sigma_2.506e-02/descomposicion_dchi2.json` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/sigma_2.506e-02/descomposicion_dchi2.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/sigma_2.506e-02/fallos_theta0.001_Delta0.05.json` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/sigma_2.506e-02/fallos_theta0.001_Delta0.10.json` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/sigma_2.506e-02/fallos_theta1_Delta0.05.json` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/sigma_2.506e-02/fallos_theta1_Delta0.10.json` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/sigma_2.506e-02/fallos_theta20_Delta0.05.json` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/sigma_2.506e-02/fallos_theta20_Delta0.10.json` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/sigma_2.506e-02/meta.json` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/sigma_2.506e-02/observables_theta0.001_Delta0.05.npz` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/sigma_2.506e-02/observables_theta0.001_Delta0.10.npz` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/sigma_2.506e-02/observables_theta1_Delta0.05.npz` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/sigma_2.506e-02/observables_theta1_Delta0.10.npz` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/sigma_2.506e-02/observables_theta20_Delta0.05.npz` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/sigma_2.506e-02/observables_theta20_Delta0.10.npz` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/sigma_2.506e-02/resumen.json` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/sigma_2.506e-02/tabla.md` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/sigma_2.633e-02/dchi2_TT_theta0.001_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/sigma_2.633e-02/dchi2_TT_theta1_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/sigma_2.633e-02/dchi2_pp_theta0.001_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/sigma_2.633e-02/dchi2_pp_theta1_Delta0.05.txt` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/sigma_2.633e-02/fallos_theta0.001_Delta0.05.json` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/sigma_2.633e-02/fallos_theta1_Delta0.05.json` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/sigma_2.633e-02/meta.json` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/sigma_2.633e-02/observables_theta0.001_Delta0.05.npz` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/sigma_2.633e-02/observables_theta1_Delta0.05.npz` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/sigma_2.633e-02/resumen.json` | = | análisis CMB/Planck | — |
| `analisis_cmb/resultados/takahashi/sigma_2.633e-02/tabla.md` | = | análisis CMB/Planck | — |
| `analisis_cmb/verificacion_numeros.md` | = | análisis CMB/Planck | — |
| `analisis_cmb/verificacion_numeros.py` | = | análisis CMB/Planck | — |
| `figures/Figure_1.png` | `archive/figures/Figure_1.png` | no citada en PREPRINT; ningún script la genera; subida a mano (500f05b «Add files via upload», 2026-07-13) | ninguna |
| `figures/Figure_2.png` | `archive/figures/Figure_2.png` | no citada en PREPRINT; ningún script la genera; subida a mano (500f05b «Add files via upload», 2026-07-13) | ninguna |
| `figures/desqueezing_heavy.png` | ✗ eliminar de figures/ | no citada en PREPRINT; copia idéntica (cmp) de `results/desqueezing_heavy/desqueezing_heavy.png` → original queda en `exploratory/results/desqueezing_heavy/desqueezing_heavy.png` | `scripts/desqueezing/desqueezing_heavy_scan.py` |
| `figures/desqueezing_relax_time.png` | ✗ eliminar de figures/ | no citada en PREPRINT; copia idéntica (cmp) de `results/desqueezing_relax_time/desqueezing_relax_time.png` → original queda en `exploratory/results/desqueezing_relax_time/desqueezing_relax_time.png` | `analisis_cmb/bloque_B/COMPARACION_B1.md` |
| `figures/euclid_A0_theta_plane.png` | ✗ eliminar de figures/ | no citada en PREPRINT; copia idéntica (cmp) de `results/euclid_protocol/euclid_A0_theta_plane.png` → original queda en `exploratory/results/euclid_protocol/euclid_A0_theta_plane.png` | `papers/euclid-protocol-vacuum-relaxation.md` |
| `figures/euclid_mcmc_E1_amplitude.png` | ✗ eliminar de figures/ | no citada en PREPRINT; copia idéntica (cmp) de `results/euclid_mcmc/euclid_mcmc_E1_amplitude.png` → original queda en `archive/results/euclid_mcmc/euclid_mcmc_E1_amplitude.png` | ninguna |
| `figures/euclid_mcmc_E1_amplitude_fixedBG.png` | ✗ eliminar de figures/ | no citada en PREPRINT; copia idéntica (cmp) de `results/euclid_mcmc/euclid_mcmc_E1_amplitude_fixedBG.png` → original queda en `archive/results/euclid_mcmc/euclid_mcmc_E1_amplitude_fixedBG.png` | ninguna |
| `figures/euclid_mcmc_E1_amplitude_freeBG.png` | ✗ eliminar de figures/ | no citada en PREPRINT; copia idéntica (cmp) de `results/euclid_mcmc/euclid_mcmc_E1_amplitude_freeBG.png` → original queda en `archive/results/euclid_mcmc/euclid_mcmc_E1_amplitude_freeBG.png` | ninguna |
| `figures/euclid_mcmc_E2_mean_reversion.png` | ✗ eliminar de figures/ | no citada en PREPRINT; copia idéntica (cmp) de `results/euclid_mcmc/euclid_mcmc_E2_mean_reversion.png` → original queda en `archive/results/euclid_mcmc/euclid_mcmc_E2_mean_reversion.png` | ninguna |
| `figures/euclid_mcmc_E2_mean_reversion_fixedBG.png` | ✗ eliminar de figures/ | no citada en PREPRINT; copia idéntica (cmp) de `results/euclid_mcmc/euclid_mcmc_E2_mean_reversion_fixedBG.png` → original queda en `archive/results/euclid_mcmc/euclid_mcmc_E2_mean_reversion_fixedBG.png` | ninguna |
| `figures/euclid_mcmc_E2_mean_reversion_freeBG.png` | ✗ eliminar de figures/ | no citada en PREPRINT; copia idéntica (cmp) de `results/euclid_mcmc/euclid_mcmc_E2_mean_reversion_freeBG.png` → original queda en `archive/results/euclid_mcmc/euclid_mcmc_E2_mean_reversion_freeBG.png` | ninguna |
| `figures/euclid_mcmc_null_F_E0.png` | ✗ eliminar de figures/ | no citada en PREPRINT; copia idéntica (cmp) de `results/euclid_mcmc/euclid_mcmc_null_F_E0.png` → original queda en `archive/results/euclid_mcmc/euclid_mcmc_null_F_E0.png` | ninguna |
| `figures/euclid_mcmc_null_F_E0_fixedBG.png` | ✗ eliminar de figures/ | no citada en PREPRINT; copia idéntica (cmp) de `results/euclid_mcmc/euclid_mcmc_null_F_E0_fixedBG.png` → original queda en `archive/results/euclid_mcmc/euclid_mcmc_null_F_E0_fixedBG.png` | ninguna |
| `figures/euclid_mcmc_null_F_E0_freeBG.png` | ✗ eliminar de figures/ | no citada en PREPRINT; copia idéntica (cmp) de `results/euclid_mcmc/euclid_mcmc_null_F_E0_freeBG.png` → original queda en `archive/results/euclid_mcmc/euclid_mcmc_null_F_E0_freeBG.png` | ninguna |
| `figures/euclid_mcmc_sigmaX_panels.png` | ✗ eliminar de figures/ | no citada en PREPRINT; copia idéntica (cmp) de `results/euclid_mcmc/euclid_mcmc_sigmaX_panels.png` → original queda en `archive/results/euclid_mcmc/euclid_mcmc_sigmaX_panels.png` | ninguna |
| `figures/euclid_theta_shape_detectability.png` | ✗ eliminar de figures/ | no citada en PREPRINT; copia idéntica (cmp) de `results/euclid_protocol/theta_shape_detectability.png` → original queda en `exploratory/results/euclid_protocol/theta_shape_detectability.png` | `papers/euclid-protocol-vacuum-relaxation.md`; `scripts/desqueezing/euclid_protocol_forecasts.py` |
| `figures/exclusion_plot.png` | `archive/figures/exclusion_plot.png` | no citada en PREPRINT; subida a mano (f247fcb, 2026-07-05); `generate_exclusion_plot()` de ou_bao_likelihood.py nunca se llama y escribe en plots/ | `scripts/ou_bao_likelihood.py` |
| `figures/falsos_positivos_detectores.png` | `archive/figures/falsos_positivos_detectores.png` | no citada en PREPRINT; ningún script la genera; subida a mano (826dd3e «Add files via upload», 2026-07-06) | ninguna |
| `figures/gpe_negative_mass_evolution.png` | `results/gpe/gpe_negative_mass_evolution.png` | no citada en PREPRINT ni copia; la genera scripts/gpe/gpe_sim.py (escribe en el cwd) → ⚠ decidir | `scripts/gpe/gpe_sim.py` |
| `figures/gpe_negative_mass_instability.png` | `results/gpe/gpe_negative_mass_instability.png` | no citada en PREPRINT ni copia; la genera scripts/gpe/gpe_sim.py (escribe en el cwd) → ⚠ decidir | `scripts/gpe/gpe_sim.py` |
| `figures/joint_mcmc_E1_amplitude_BAO.png` | `archive/figures/joint_mcmc_E1_amplitude_BAO.png` | no citada en PREPRINT; versión distinta (cmp difiere) de la de results/euclid_joint_mcmc/; creada en 470619a con el MCMC Euclid ya retirado | ninguna |
| `figures/joint_mcmc_E1_amplitude_BAO_SN.png` | ✗ eliminar de figures/ | no citada en PREPRINT; copia idéntica (cmp) de `results/euclid_joint_mcmc/joint_mcmc_E1_amplitude_BAO_SN.png` → original queda en `archive/results/euclid_joint_mcmc/joint_mcmc_E1_amplitude_BAO_SN.png` | ninguna |
| `figures/joint_mcmc_E2_mean_reversion_BAO.png` | `archive/figures/joint_mcmc_E2_mean_reversion_BAO.png` | no citada en PREPRINT; versión distinta (cmp difiere) de la de results/euclid_joint_mcmc/; creada en 470619a con el MCMC Euclid ya retirado | ninguna |
| `figures/joint_mcmc_E2_mean_reversion_BAO_SN.png` | ✗ eliminar de figures/ | no citada en PREPRINT; copia idéntica (cmp) de `results/euclid_joint_mcmc/joint_mcmc_E2_mean_reversion_BAO_SN.png` → original queda en `archive/results/euclid_joint_mcmc/joint_mcmc_E2_mean_reversion_BAO_SN.png` | ninguna |
| `figures/joint_mcmc_null_F_E0_BAO.png` | `archive/figures/joint_mcmc_null_F_E0_BAO.png` | no citada en PREPRINT; versión distinta (cmp difiere) de la de results/euclid_joint_mcmc/; creada en 470619a con el MCMC Euclid ya retirado | ninguna |
| `figures/joint_mcmc_null_F_E0_BAO_SN.png` | ✗ eliminar de figures/ | no citada en PREPRINT; copia idéntica (cmp) de `results/euclid_joint_mcmc/joint_mcmc_null_F_E0_BAO_SN.png` → original queda en `archive/results/euclid_joint_mcmc/joint_mcmc_null_F_E0_BAO_SN.png` | ninguna |
| `figures/joint_mcmc_w0_width_compare.png` | ✗ eliminar de figures/ | no citada en PREPRINT; copia idéntica (cmp) de `results/euclid_joint_mcmc/joint_mcmc_w0_width_compare.png` → original queda en `archive/results/euclid_joint_mcmc/joint_mcmc_w0_width_compare.png` | ninguna |
| `figures/joint_w0wa_sigma.png` | ✗ eliminar de figures/ | no citada en PREPRINT; copia idéntica (cmp) de `results/joint_w0wa_sigma/joint_w0wa_sigma.png` → original queda en `results/joint_w0wa_sigma/joint_w0wa_sigma.png` | `analisis_cmb/AUDITORIA.md`; `analisis_cmb/bloque_A/COMPARACION.md` |
| `figures/profile_sigma_x.png` | ✗ eliminar de figures/ | no citada en PREPRINT; copia idéntica (cmp) de `results/profile_sigma_x/profile_sigma_x.png` → original queda en `results/profile_sigma_x/profile_sigma_x.png` | `analisis_cmb/AUDITORIA.md`; `analisis_cmb/bloque_B/COMPARACION_B1.md` |
| `figures/sdiff_A0_min_for_detection.png` | ✗ eliminar de figures/ | no citada en PREPRINT; copia idéntica (cmp) de `results/sdiff_discrimination/A0_min_for_detection.png` → original queda en `exploratory/results/sdiff_discrimination/A0_min_for_detection.png` | `papers/fundamental-vs-emergent-vacuum-relaxation.md`; `scripts/desqueezing/sdiff_discrimination.py` |
| `figures/sdiff_residual_vs_theta.png` | ✗ eliminar de figures/ | no citada en PREPRINT; copia idéntica (cmp) de `results/sdiff_discrimination/residual_vs_theta.png` → original queda en `exploratory/results/sdiff_discrimination/residual_vs_theta.png` | `papers/fundamental-vs-emergent-vacuum-relaxation.md`; `scripts/desqueezing/sdiff_discrimination.py` |
| `figures/sdiff_t_half_vs_theta.png` | ✗ eliminar de figures/ | no citada en PREPRINT; copia idéntica (cmp) de `results/sdiff_discrimination/t_half_vs_theta.png` → original queda en `exploratory/results/sdiff_discrimination/t_half_vs_theta.png` | `papers/fundamental-vs-emergent-vacuum-relaxation.md`; `scripts/desqueezing/sdiff_discrimination.py` |
| `figures/termico_vs_squeezed_recomendado.png` | ✗ eliminar de figures/ | no citada en PREPRINT; copia idéntica (cmp) de `results/termico_vs_squeezed_recomendado/termico_vs_squeezed_recomendado.png` → original queda en `exploratory/results/termico_vs_squeezed_recomendado/termico_vs_squeezed_recomendado.png` | `scripts/desqueezing/termico_vs_squeezed_recomendado.py` |
| `figures/test_desi_QNM.png` | `archive/figures/test_desi_QNM.png` | no citada en PREPRINT; solo la escribe en figures/ ou_bao_stochastic_test.py (DR1, va a archive/); ou_bao_likelihood.py escribe plots/test_desi_QNM.png (DR2, ignorado) | `results/ou_bao_desi_dr2_run.txt`; `scripts/ou_bao_likelihood.py`; `scripts/ou_bao_stochastic_test.py` |
| `manuscript/CLAIMS.md` | = | preprint / claims | — |
| `manuscript/COVER_LETTER.md` | = | preprint / claims | — |
| `manuscript/PREPRINT.md` | = | preprint / claims | — |
| `manuscript/README.md` | = | preprint / claims | — |
| `notes/desqueezing-relaxation-vacuum-fluctuations-note.md` | `exploratory/notes/desqueezing-relaxation-vacuum-fluctuations-note.md` | nota de trabajo del programa desqueezing | `analisis_cmb/AUDITORIA.md`; `analisis_cmb/bloque_B/COMPARACION_B2.md`; `papers/amplification-gap.md`; `papers/anisotropic-slip-option0.md`; `papers/fundamental-vs-emergent-vacuum-relaxation.md` |
| `notes/mapping_tables.md` | `exploratory/notes/mapping_tables.md` | nota de trabajo del programa desqueezing | `analisis_cmb/AUDITORIA.md`; `analisis_cmb/bloque_B/COMPARACION_B1.md`; `notes/desqueezing-relaxation-vacuum-fluctuations-note.md`; `scripts/desqueezing/cosmological_mapping.py` |
| `papers/BOUND_CORRECTION_SIGMA_X.md` | = | sostiene C1–C7 | — |
| `papers/EXPLORATORY_BOUNDARY.md` | `theory/EXPLORATORY_BOUNDARY.md` | teoría/programa, fuera del claim set C1–C7 | `README.md` |
| `papers/HONEST_HEADLINES.md` | = | sostiene C1–C7 | — |
| `papers/THEORY_PROGRAMME_POINTER.md` | `theory/THEORY_PROGRAMME_POINTER.md` | teoría/programa, fuera del claim set C1–C7 | ninguna |
| `papers/amplification-gap.md` | = | sostiene C1–C7 | — |
| `papers/anisotropic-slip-option0.md` | = | sostiene C1–C7 | — |
| `papers/covariance-desi-dr2-full.md` | = | sostiene C1–C7 | — |
| `papers/data-pack-option0-internet.md` | `exploratory/notes/data-pack-option0-internet.md` | data pack Option 0 (desqueezing/internet), exploratorio | `analisis_cmb/AUDITORIA.md`; `analisis_cmb/bloque_B/COMPARACION_B2.md`; `analisis_cmb/verificacion_numeros.md`; `analisis_cmb/verificacion_numeros.py`; `papers/resume.txt`; `scripts/slip_bridge.py` |
| `papers/euclid-protocol-vacuum-relaxation.md` | `theory/euclid-protocol-vacuum-relaxation.md` | teoría/programa, fuera del claim set C1–C7 | `analisis_cmb/AUDITORIA.md`; `analisis_cmb/bloque_B/COMPARACION_B2.md`; `analisis_cmb/verificacion_numeros.md`; `analisis_cmb/verificacion_numeros.py`; `papers/fundamental-vs-emergent-vacuum-relaxation.md` |
| `papers/fundamental-vs-emergent-vacuum-relaxation.md` | `theory/fundamental-vs-emergent-vacuum-relaxation.md` | teoría/programa, fuera del claim set C1–C7 | `analisis_cmb/AUDITORIA.md`; `analisis_cmb/bloque_B/COMPARACION_B2.md`; `papers/amplification-gap.md`; `papers/euclid-protocol-vacuum-relaxation.md`; `papers/sdiff-fundamental-vs-emergent.md` |
| `papers/neutrino-de-degeneracy-and-sigma-x.md` | `theory/neutrino-de-degeneracy-and-sigma-x.md` | teoría/programa, fuera del claim set C1–C7 | `analisis_cmb/AUDITORIA.md`; `analisis_cmb/bloque_A/COMPARACION.md`; `manuscript/CLAIMS.md` |
| `papers/principle-of-vacuum-smoothness.md` | `theory/principle-of-vacuum-smoothness.md` | teoría/programa, fuera del claim set C1–C7 | `README.md`; `analisis_cmb/AUDITORIA.md`; `analisis_cmb/bloque_B/COMPARACION_B2.md`; `papers/quantum-fluid-instabilities-desi-dr2.md` |
| `papers/quantum-fluid-instabilities-desi-dr2.md` | = | sostiene C1–C7 | — |
| `papers/resume.txt` | = | ⚠ no listado; snapshot numérico C1–C7 | — |
| `papers/scope-and-mixups.md` | = | sostiene C1–C7 | — |
| `papers/sdiff-fundamental-vs-emergent.md` | `theory/sdiff-fundamental-vs-emergent.md` | teoría/programa, fuera del claim set C1–C7 | `analisis_cmb/AUDITORIA.md`; `analisis_cmb/bloque_B/COMPARACION_B2.md`; `papers/anisotropic-slip-option0.md`; `papers/euclid-protocol-vacuum-relaxation.md`; `papers/fundamental-vs-emergent-vacuum-relaxation.md`; `papers/unimodular-gravity-vacuum-smoothness.md`; `scripts/desqueezing/cosmological_mapping.py` |
| `papers/sensitivity_kernel_table.md` | = | sostiene C1–C7 | — |
| `papers/stochastic-dark-energy-desi-dr2.md` | = | sostiene C1–C7 | — |
| `papers/unimodular-gravity-vacuum-smoothness.md` | `theory/unimodular-gravity-vacuum-smoothness.md` | teoría/programa, fuera del claim set C1–C7 | `README.md`; `analisis_cmb/AUDITORIA.md`; `analisis_cmb/NOTAS_ANALISIS_CMB.md`; `analisis_cmb/RESULTADO.md`; `analisis_cmb/bloque_B/COMPARACION_B2.md`; `analisis_cmb/verificacion_numeros.md` |
| `requirements.txt` | = | raíz | — |
| `results/amplification_routes/VERDICT.md` | = | salida del pipeline DESI | — |
| `results/amplification_routes/route1_N_eff_grid.csv` | = | salida del pipeline DESI | — |
| `results/amplification_routes/route1_N_eff_required.csv` | = | salida del pipeline DESI | — |
| `results/amplification_routes/route2_horizon_exit_scan.csv` | = | salida del pipeline DESI | — |
| `results/amplification_routes/route3_avalanche_scan.csv` | = | salida del pipeline DESI | — |
| `results/cosmological_mapping/mapping_summary.txt` | `exploratory/results/cosmological_mapping/mapping_summary.txt` | salida del programa desqueezing (exploratorio) | `analisis_cmb/AUDITORIA.md`; `analisis_cmb/bloque_B/COMPARACION_B1.md`; `scripts/desqueezing/cosmological_mapping.py` |
| `results/cosmological_mapping/mapping_tables_script.md` | `exploratory/results/cosmological_mapping/mapping_tables_script.md` | salida del programa desqueezing (exploratorio) | `analisis_cmb/bloque_B/COMPARACION_B1.md`; `scripts/desqueezing/cosmological_mapping.py` |
| `results/cosmological_mapping/table1_theta.csv` | `exploratory/results/cosmological_mapping/table1_theta.csv` | salida del programa desqueezing (exploratorio) | `scripts/desqueezing/cosmological_mapping.py` |
| `results/cosmological_mapping/table2_A0_required.csv` | `exploratory/results/cosmological_mapping/table2_A0_required.csv` | salida del programa desqueezing (exploratorio) | `analisis_cmb/AUDITORIA.md`; `analisis_cmb/bloque_B/COMPARACION_B1.md`; `notes/mapping_tables.md`; `scripts/desqueezing/cosmological_mapping.py` |
| `results/cosmological_mapping/table3_page_scrambling.csv` | `exploratory/results/cosmological_mapping/table3_page_scrambling.csv` | salida del programa desqueezing (exploratorio) | `scripts/desqueezing/cosmological_mapping.py` |
| `results/cosmological_mapping/table4_DESI_S_z.csv` | `exploratory/results/cosmological_mapping/table4_DESI_S_z.csv` | salida del programa desqueezing (exploratorio) | `scripts/desqueezing/cosmological_mapping.py` |
| `results/desi_cov/desi_cov_alpha_iso_7x7.txt` | = | salida del pipeline DESI | — |
| `results/desi_cov/desi_gaussian_bao_ALL_GCcomb_cov.txt` | = | salida del pipeline DESI | — |
| `results/desi_dr2_real_bao/desi_dr2_real_bao.json` | = | salida del pipeline DESI | — |
| `results/desi_dr2_real_bao/desi_dr2_real_bao.txt` | = | salida del pipeline DESI | — |
| `results/desqueezing_heavy/desqueezing_heavy.png` | `exploratory/results/desqueezing_heavy/desqueezing_heavy.png` | salida del programa desqueezing (exploratorio) | `scripts/desqueezing/desqueezing_heavy_scan.py` |
| `results/desqueezing_heavy/desqueezing_heavy_grid.csv` | `exploratory/results/desqueezing_heavy/desqueezing_heavy_grid.csv` | salida del programa desqueezing (exploratorio) | `scripts/desqueezing/desqueezing_heavy_scan.py` |
| `results/desqueezing_heavy/desqueezing_heavy_summary.json` | `exploratory/results/desqueezing_heavy/desqueezing_heavy_summary.json` | salida del programa desqueezing (exploratorio) | `scripts/desqueezing/desqueezing_heavy_scan.py` |
| `results/desqueezing_heavy/desqueezing_heavy_summary.txt` | `exploratory/results/desqueezing_heavy/desqueezing_heavy_summary.txt` | salida del programa desqueezing (exploratorio) | `scripts/desqueezing/desqueezing_heavy_scan.py` |
| `results/desqueezing_relax_time/desqueezing_relax_time.png` | `exploratory/results/desqueezing_relax_time/desqueezing_relax_time.png` | salida del programa desqueezing (exploratorio) | `analisis_cmb/bloque_B/COMPARACION_B1.md`; `results/desqueezing_relax_time/desqueezing_relax_time_resumen.txt` |
| `results/desqueezing_relax_time/desqueezing_relax_time_resumen.txt` | `exploratory/results/desqueezing_relax_time/desqueezing_relax_time_resumen.txt` | salida del programa desqueezing (exploratorio) | `analisis_cmb/bloque_B/COMPARACION_B1.md`; `scripts/desqueezing/desqueezing_relax_time.py` |
| `results/desqueezing_relax_time/t_relax_tabla.csv` | `exploratory/results/desqueezing_relax_time/t_relax_tabla.csv` | salida del programa desqueezing (exploratorio) | `analisis_cmb/bloque_B/COMPARACION_B1.md`; `results/desqueezing_relax_time/desqueezing_relax_time_resumen.txt`; `scripts/desqueezing/desqueezing_relax_time.py` |
| `results/eos_cpl_desi_dr2/eos_cpl_summary.json` | = | salida del pipeline DESI | — |
| `results/eos_cpl_desi_dr2/eos_cpl_summary.txt` | = | salida del pipeline DESI | — |
| `results/euclid_joint_mcmc/RETIRED.md` | `archive/results/euclid_joint_mcmc/RETIRED.md` | MCMC Euclid sintético retirado (N9) | `papers/euclid-protocol-vacuum-relaxation.md` |
| `results/euclid_joint_mcmc/joint_mcmc_E1_amplitude_BAO.png` | `archive/results/euclid_joint_mcmc/joint_mcmc_E1_amplitude_BAO.png` | MCMC Euclid sintético retirado (N9) | ninguna |
| `results/euclid_joint_mcmc/joint_mcmc_E1_amplitude_BAO_SN.png` | `archive/results/euclid_joint_mcmc/joint_mcmc_E1_amplitude_BAO_SN.png` | MCMC Euclid sintético retirado (N9) | ninguna |
| `results/euclid_joint_mcmc/joint_mcmc_E2_mean_reversion_BAO.png` | `archive/results/euclid_joint_mcmc/joint_mcmc_E2_mean_reversion_BAO.png` | MCMC Euclid sintético retirado (N9) | ninguna |
| `results/euclid_joint_mcmc/joint_mcmc_E2_mean_reversion_BAO_SN.png` | `archive/results/euclid_joint_mcmc/joint_mcmc_E2_mean_reversion_BAO_SN.png` | MCMC Euclid sintético retirado (N9) | ninguna |
| `results/euclid_joint_mcmc/joint_mcmc_null_F_E0_BAO.png` | `archive/results/euclid_joint_mcmc/joint_mcmc_null_F_E0_BAO.png` | MCMC Euclid sintético retirado (N9) | ninguna |
| `results/euclid_joint_mcmc/joint_mcmc_null_F_E0_BAO_SN.png` | `archive/results/euclid_joint_mcmc/joint_mcmc_null_F_E0_BAO_SN.png` | MCMC Euclid sintético retirado (N9) | ninguna |
| `results/euclid_joint_mcmc/joint_mcmc_w0_width_compare.png` | `archive/results/euclid_joint_mcmc/joint_mcmc_w0_width_compare.png` | MCMC Euclid sintético retirado (N9) | ninguna |
| `results/euclid_joint_mcmc/joint_summaries.json` | `archive/results/euclid_joint_mcmc/joint_summaries.json` | MCMC Euclid sintético retirado (N9) | ninguna |
| `results/euclid_joint_mcmc/joint_summary.txt` | `archive/results/euclid_joint_mcmc/joint_summary.txt` | MCMC Euclid sintético retirado (N9) | ninguna |
| `results/euclid_mcmc/RETIRED.md` | `archive/results/euclid_mcmc/RETIRED.md` | MCMC Euclid sintético retirado (N9) | `papers/euclid-protocol-vacuum-relaxation.md` |
| `results/euclid_mcmc/euclid_mcmc_E1_amplitude.png` | `archive/results/euclid_mcmc/euclid_mcmc_E1_amplitude.png` | MCMC Euclid sintético retirado (N9) | ninguna |
| `results/euclid_mcmc/euclid_mcmc_E1_amplitude_fixedBG.png` | `archive/results/euclid_mcmc/euclid_mcmc_E1_amplitude_fixedBG.png` | MCMC Euclid sintético retirado (N9) | ninguna |
| `results/euclid_mcmc/euclid_mcmc_E1_amplitude_freeBG.png` | `archive/results/euclid_mcmc/euclid_mcmc_E1_amplitude_freeBG.png` | MCMC Euclid sintético retirado (N9) | ninguna |
| `results/euclid_mcmc/euclid_mcmc_E2_mean_reversion.png` | `archive/results/euclid_mcmc/euclid_mcmc_E2_mean_reversion.png` | MCMC Euclid sintético retirado (N9) | ninguna |
| `results/euclid_mcmc/euclid_mcmc_E2_mean_reversion_fixedBG.png` | `archive/results/euclid_mcmc/euclid_mcmc_E2_mean_reversion_fixedBG.png` | MCMC Euclid sintético retirado (N9) | ninguna |
| `results/euclid_mcmc/euclid_mcmc_E2_mean_reversion_freeBG.png` | `archive/results/euclid_mcmc/euclid_mcmc_E2_mean_reversion_freeBG.png` | MCMC Euclid sintético retirado (N9) | ninguna |
| `results/euclid_mcmc/euclid_mcmc_null_F_E0.png` | `archive/results/euclid_mcmc/euclid_mcmc_null_F_E0.png` | MCMC Euclid sintético retirado (N9) | ninguna |
| `results/euclid_mcmc/euclid_mcmc_null_F_E0_fixedBG.png` | `archive/results/euclid_mcmc/euclid_mcmc_null_F_E0_fixedBG.png` | MCMC Euclid sintético retirado (N9) | ninguna |
| `results/euclid_mcmc/euclid_mcmc_null_F_E0_freeBG.png` | `archive/results/euclid_mcmc/euclid_mcmc_null_F_E0_freeBG.png` | MCMC Euclid sintético retirado (N9) | ninguna |
| `results/euclid_mcmc/euclid_mcmc_sigmaX_panels.png` | `archive/results/euclid_mcmc/euclid_mcmc_sigmaX_panels.png` | MCMC Euclid sintético retirado (N9) | ninguna |
| `results/euclid_mcmc/mcmc_summaries.json` | `archive/results/euclid_mcmc/mcmc_summaries.json` | MCMC Euclid sintético retirado (N9) | ninguna |
| `results/euclid_mcmc/mcmc_summary.txt` | `archive/results/euclid_mcmc/mcmc_summary.txt` | MCMC Euclid sintético retirado (N9) | `analisis_cmb/AUDITORIA.md`; `analisis_cmb/bloque_B/COMPARACION_B1.md` |
| `results/euclid_protocol/euclid_A0_theta_plane.png` | `exploratory/results/euclid_protocol/euclid_A0_theta_plane.png` | salida del programa desqueezing (exploratorio) | `analisis_cmb/bloque_B/COMPARACION_B1.md` |
| `results/euclid_protocol/euclid_forecast_grid.csv` | `exploratory/results/euclid_protocol/euclid_forecast_grid.csv` | salida del programa desqueezing (exploratorio) | `analisis_cmb/AUDITORIA.md`; `analisis_cmb/bloque_B/COMPARACION_B1.md`; `papers/euclid-protocol-vacuum-relaxation.md`; `scripts/desqueezing/euclid_protocol_forecasts.py` |
| `results/euclid_protocol/theta_shape_detectability.png` | `exploratory/results/euclid_protocol/theta_shape_detectability.png` | salida del programa desqueezing (exploratorio) | `papers/euclid-protocol-vacuum-relaxation.md`; `scripts/desqueezing/euclid_protocol_forecasts.py` |
| `results/joint_w0wa_sigma/joint_w0wa_sigma.json` | = | salida del pipeline DESI | — |
| `results/joint_w0wa_sigma/joint_w0wa_sigma.png` | = | salida del pipeline DESI | — |
| `results/joint_w0wa_sigma/joint_w0wa_sigma.txt` | = | salida del pipeline DESI | — |
| `results/ou_bao_desi_dr2_run.txt` | `archive/results/ou_bao_desi_dr2_run.txt` | log de una ejecución antigua de ou_bao | `analisis_cmb/AUDITORIA.md`; `analisis_cmb/NOTAS_ANALISIS_CMB.md`; `analisis_cmb/bloque_B/COMPARACION_B1.md`; `analisis_cmb/bloque_B/COMPARACION_B2.md`; `papers/euclid-protocol-vacuum-relaxation.md`; `papers/resume.txt`; `scripts/desqueezing/cosmological_mapping.py` |
| `results/profile_sigma_x/profile_sigma_x.json` | = | salida del pipeline DESI | — |
| `results/profile_sigma_x/profile_sigma_x.png` | = | salida del pipeline DESI | — |
| `results/profile_sigma_x/profile_sigma_x.txt` | = | salida del pipeline DESI | — |
| `results/sdiff_discrimination/A0_min_for_detection.png` | `exploratory/results/sdiff_discrimination/A0_min_for_detection.png` | salida del programa desqueezing (exploratorio) | `papers/fundamental-vs-emergent-vacuum-relaxation.md`; `scripts/desqueezing/sdiff_discrimination.py` |
| `results/sdiff_discrimination/residual_vs_theta.png` | `exploratory/results/sdiff_discrimination/residual_vs_theta.png` | salida del programa desqueezing (exploratorio) | `papers/fundamental-vs-emergent-vacuum-relaxation.md`; `scripts/desqueezing/sdiff_discrimination.py` |
| `results/sdiff_discrimination/sdiff_discrimination_grid.csv` | `exploratory/results/sdiff_discrimination/sdiff_discrimination_grid.csv` | salida del programa desqueezing (exploratorio) | `analisis_cmb/bloque_B/COMPARACION_B1.md`; `papers/fundamental-vs-emergent-vacuum-relaxation.md`; `results/sdiff_discrimination/summary.txt`; `scripts/desqueezing/sdiff_discrimination.py` |
| `results/sdiff_discrimination/summary.txt` | `exploratory/results/sdiff_discrimination/summary.txt` | salida del programa desqueezing (exploratorio) | `analisis_cmb/AUDITORIA.md`; `analisis_cmb/bloque_B/COMPARACION_B1.md`; `papers/fundamental-vs-emergent-vacuum-relaxation.md` |
| `results/sdiff_discrimination/t_half_vs_theta.png` | `exploratory/results/sdiff_discrimination/t_half_vs_theta.png` | salida del programa desqueezing (exploratorio) | `papers/fundamental-vs-emergent-vacuum-relaxation.md`; `scripts/desqueezing/sdiff_discrimination.py` |
| `results/tachyonic_rank1/SUMMARY.json` | = | salida del pipeline DESI | — |
| `results/tachyonic_rank1/SUMMARY.txt` | = | salida del pipeline DESI | — |
| `results/tachyonic_rank1/scan_fixed_sigma0.txt` | = | salida del pipeline DESI | — |
| `results/termico_vs_squeezed_recomendado/termico_vs_squeezed_recomendado.png` | `exploratory/results/termico_vs_squeezed_recomendado/termico_vs_squeezed_recomendado.png` | salida del programa desqueezing (exploratorio) | `scripts/desqueezing/termico_vs_squeezed_recomendado.py` |
| `results/termico_vs_squeezed_recomendado/termico_vs_squeezed_recomendado_resumen.txt` | `exploratory/results/termico_vs_squeezed_recomendado/termico_vs_squeezed_recomendado_resumen.txt` | salida del programa desqueezing (exploratorio) | `scripts/desqueezing/termico_vs_squeezed_recomendado.py` |
| `scripts/amplification/route1_local_causal_set_seed.py` | = | pipeline DESI | — |
| `scripts/amplification/route2_late_horizon_exit.py` | = | pipeline DESI | — |
| `scripts/amplification/route3_nonlinear_avalanche.py` | = | pipeline DESI | — |
| `scripts/amplification/run_all_routes.py` | = | pipeline DESI | — |
| `scripts/amplifier_audit.py` | = | pipeline DESI | — |
| `scripts/cross_correlation_DESI.py` | `archive/scripts/cross_correlation_DESI.py` | correlación cruzada preliminar (N5) — ⚠ la cita PREPRINT §6 (línea 265) | `analisis_cmb/AUDITORIA.md`; `manuscript/PREPRINT.md`; `papers/stochastic-dark-energy-desi-dr2.md` |
| `scripts/descargarcovarianza.py` | = | pipeline DESI | — |
| `scripts/desi_bound.py` | = | pipeline DESI | — |
| `scripts/desi_dr2_data.py` | = | pipeline DESI | — |
| `scripts/desi_dr2_real_bao_test.py` | = | pipeline DESI | — |
| `scripts/desqueezing/cosmological_mapping.py` | `exploratory/scripts/cosmological_mapping.py` | programa desqueezing (exploratorio) | `analisis_cmb/AUDITORIA.md`; `analisis_cmb/bloque_B/COMPARACION_B1.md`; `notes/desqueezing-relaxation-vacuum-fluctuations-note.md` |
| `scripts/desqueezing/desqueezing_heavy_scan.py` | `exploratory/scripts/desqueezing_heavy_scan.py` | programa desqueezing (exploratorio) | ninguna |
| `scripts/desqueezing/desqueezing_relax_time.py` | `exploratory/scripts/desqueezing_relax_time.py` | programa desqueezing (exploratorio) | `analisis_cmb/AUDITORIA.md`; `analisis_cmb/bloque_B/COMPARACION_B1.md`; `notes/desqueezing-relaxation-vacuum-fluctuations-note.md`; `scripts/desqueezing/cosmological_mapping.py` |
| `scripts/desqueezing/euclid_protocol_forecasts.py` | `exploratory/scripts/euclid_protocol_forecasts.py` | programa desqueezing (exploratorio) | `analisis_cmb/AUDITORIA.md`; `analisis_cmb/bloque_B/COMPARACION_B1.md`; `papers/euclid-protocol-vacuum-relaxation.md` |
| `scripts/desqueezing/sdiff_discrimination.py` | `exploratory/scripts/sdiff_discrimination.py` | programa desqueezing (exploratorio) | `analisis_cmb/AUDITORIA.md`; `analisis_cmb/bloque_B/COMPARACION_B1.md`; `papers/fundamental-vs-emergent-vacuum-relaxation.md` |
| `scripts/desqueezing/termico_vs_squeezed_recomendado.py` | `exploratory/scripts/termico_vs_squeezed_recomendado.py` | programa desqueezing (exploratorio) | `notes/desqueezing-relaxation-vacuum-fluctuations-note.md`; `results/termico_vs_squeezed_recomendado/termico_vs_squeezed_recomendado_resumen.txt` |
| `scripts/eos_efectiva.py` | = | pipeline DESI | — |
| `scripts/euclid_joint_bao_sne_mcmc.py` | `archive/scripts/euclid_joint_bao_sne_mcmc.py` | stub DEPRECATED que solo reenvía a desi_dr2_real_bao_test.py | `papers/euclid-protocol-vacuum-relaxation.md` |
| `scripts/euclid_mock_mcmc.py` | `archive/scripts/euclid_mock_mcmc.py` | stub DEPRECATED que solo reenvía a desi_dr2_real_bao_test.py | `papers/euclid-protocol-vacuum-relaxation.md` |
| `scripts/extraer_cov_cobaya.py` | = | pipeline DESI | — |
| `scripts/gap_two_targets.py` | = | pipeline DESI | — |
| `scripts/gpe/Bogoliubov_DESI_extension.py` | = | pipeline DESI | — |
| `scripts/gpe/gpe_sim.py` | = | pipeline DESI | — |
| `scripts/joint_w0wa_sigma_desi.py` | = | pipeline DESI | — |
| `scripts/ou_bao_likelihood.py` | = | pipeline DESI | — |
| `scripts/ou_bao_stochastic_test.py` | `archive/scripts/ou_bao_stochastic_test.py` | versión DR1 hard-coded de ou_bao_likelihood.py (ver diff); el PREPRINT cita la otra | `analisis_cmb/AUDITORIA.md`; `papers/covariance-desi-dr2-full.md`; `papers/euclid-protocol-vacuum-relaxation.md`; `scripts/desqueezing/cosmological_mapping.py`; `scripts/ou_bao_likelihood.py` |
| `scripts/profile_sigma_x_desi.py` | = | pipeline DESI | — |
| `scripts/slip_bridge.py` | = | pipeline DESI | — |
| `scripts/tachyonic_rank1_mle.py` | = | pipeline DESI | — |
| `scripts/test_desi_dr2_correlations.py` | `archive/scripts/test_desi_dr2_correlations.py` | test de correlaciones DR2 suelto, escribe en plots/ (ignorado) | ninguna |
| `scripts/verify_bao_pipelines.py` | = | pipeline DESI | — |

## (a) Scripts con `ROOT` o rutas relativas que se romperían al moverse

| Script (destino) | Qué se rompe | Corrección |
|:--|:--|:--|
| `scripts/desqueezing/cosmological_mapping.py` → `exploratory/scripts/` | l.18 `sys.path.insert(parents[1])` (= `scripts/`, para `from desi_bound import …`); l.21 `ROOT` = tres `dirname` (raíz); l.22 `RESULTS = ROOT/results/cosmological_mapping` | `sys.path` → `ROOT/"scripts"`; `ROOT` no cambia (sigue a 3 `dirname`); `RESULTS` → `ROOT/exploratory/results/cosmological_mapping` |
| `scripts/desqueezing/desqueezing_relax_time.py` → ídem | l.33 `sys.path` `parents[1]` (`from desi_bound import`); l.38 `ROOT` (3 niveles); l.39 `RESULTS = ROOT/results/desqueezing_relax_time` | igual que arriba |
| `scripts/desqueezing/euclid_protocol_forecasts.py` → ídem | l.12 `sys.path` `parents[1]` (`desi_bound`); l.15 `ROOT`; l.16 `OUT = ROOT/results/euclid_protocol`; l.17 `FIG = ROOT/figures` (+ `savefig` l.75, l.90) | `OUT` → `exploratory/results/…`; quitar las copias a `figures/` |
| `scripts/desqueezing/sdiff_discrimination.py` → ídem | l.12 `sys.path` `parents[1]` (`desi_bound`); l.15 `ROOT`; l.16 `OUT`; l.17 `FIG` (+ `savefig` l.85, l.100, l.114) | ídem |
| `scripts/desqueezing/desqueezing_heavy_scan.py` → ídem | l.19 `ROOT = parents[2]`; l.20 `OUT = ROOT/results/desqueezing_heavy`; l.22 `FIG = ROOT/figures` (+ `savefig` l.153) | `OUT` → `exploratory/results/…`; quitar la copia a `figures/` |
| `scripts/desqueezing/termico_vs_squeezed_recomendado.py` → ídem | l.13 `ROOT` = **dos** `dirname` → `scripts/` (**ya roto hoy**: escribe en `scripts/results/…`, el mismo bug que B1 corrigió en los otros dos); l.14 `RESULTS` | `ROOT` = 3 niveles; `RESULTS` → `exploratory/results/termico_vs_squeezed_recomendado` |
| `scripts/euclid_mock_mcmc.py`, `scripts/euclid_joint_bao_sne_mcmc.py` → `archive/scripts/` | l.18–19 `REAL = Path(__file__).parent / "desi_dr2_real_bao_test.py"` (deja de existir junto a ellos) | ninguna; quedan archivados como stubs (o se apunta `REAL` a `../../scripts/`) |
| `scripts/ou_bao_stochastic_test.py` → `archive/scripts/` | l.550–551 `savefig('figures/test_desi_QNM.png')` relativo al cwd | ninguna (archivado) |
| `scripts/cross_correlation_DESI.py` → `archive/scripts/` | l.52–74, 619–620 `data/…`, l.668–669 `results/…` relativos al cwd | ninguna (archivado); ya dependía del cwd |
| `scripts/test_desi_dr2_correlations.py` → `archive/scripts/` | l.27 `../cobaya_packages_dr2/…`; l.200–201 `plots/` | ninguna (archivado) |
| `analisis_cmb/diagnosticos/obsoleto/*.py` → `archive/analisis_cmb/diagnosticos/obsoleto/` | `tabla_H0fijo.py` l.10 `import planck_datos` (necesita `analisis_cmb/` en `sys.path`); `HERE / "resultados"` en los dos (**ya roto hoy**: apunta a `obsoleto/resultados/`) | ninguna (archivado) |
| `scripts/gpe/gpe_sim.py` (se queda) | l.128, l.142 `savefig("gpe_….png")` en el cwd | `savefig` → `ROOT/results/gpe/` si se acepta ⚠3 |
| Scripts que **se quedan** y escriben copias en `figures/` | `profile_sigma_x_desi.py` l.148, `joint_w0wa_sigma_desi.py` l.27/311 | quitar la segunda `savefig` (o dejarla si el preprint pasa a citar figuras) |
| Scripts que **se quedan** y citan rutas que se mueven (solo en comentarios o docstrings) | `scripts/slip_bridge.py` l.33 (`papers/data-pack-option0-internet.md`); `analisis_cmb/verificacion_numeros.py` l.45, l.70 (data-pack y euclid-protocol, por nombre); `scripts/ou_bao_likelihood.py` l.4, l.38 (`ou_bao_stochastic_test.py`) | actualizar el texto |

Nota: en todos los scripts de `exploratory/scripts/` el número de niveles hasta la raíz es el mismo que hoy (`exploratory/scripts/x.py` frente a `scripts/desqueezing/x.py`), así que `ROOT` con tres `dirname` / `parents[2]` sigue valiendo. Lo que hay que cambiar es el `sys.path` (`parents[1]` ya no es `scripts/`) y las subrutas `results/…` y `figures/`.

`.gitignore`: las reglas `results/euclid_mcmc/*.np[zy]`, `results/euclid_joint_mcmc/*.npz` y `results/**/*.npz` dejan de cubrir `archive/results/…` y `exploratory/results/…` → añadir `archive/**/*.npz`, `exploratory/**/*.npz` (y `*.npy`).

## (b) Documentos con enlaces o rutas relativas que habrá que corregir

Enlaces Markdown `[..](ruta)` (se rompen de verdad en GitHub):

| Documento | Enlace | Corrección |
|:--|:--|:--|
| `README.md` | `papers/EXPLORATORY_BOUNDARY.md` | `theory/EXPLORATORY_BOUNDARY.md` |
| `manuscript/CLAIMS.md` (N8) | `../papers/neutrino-de-degeneracy-and-sigma-x.md` | `../theory/…` |
| `papers/EXPLORATORY_BOUNDARY.md` (→ `theory/`) | `scope-and-mixups.md` | `../papers/scope-and-mixups.md` (`../manuscript/CLAIMS.md` sigue valiendo) |

Rutas entre comillas invertidas (no son enlaces, pero quedan desactualizadas):

| Documento | Rutas que cambian |
|:--|:--|
| `README.md` | l.42 `papers/unimodular-gravity-vacuum-smoothness.md`, `papers/principle-of-vacuum-smoothness.md` → `theory/`; bloque *Layout* (l.46–53): añadir `theory/`, `exploratory/`, `archive/`, `analisis_cmb/`; quitar `notes/` y `figures/` |
| `manuscript/PREPRINT.md` | l.265 `scripts/cross_correlation_DESI.py` → `archive/scripts/…` (⚠1) |
| `manuscript/CLAIMS.md` | N9: `results/euclid_mcmc/`, `results/euclid_joint_mcmc/` → `archive/results/…` |
| `papers/stochastic-dark-energy-desi-dr2.md` | `scripts/cross_correlation_DESI.py` |
| `papers/amplification-gap.md` | `notes/desqueezing-…-note.md`; `fundamental-vs-emergent-vacuum-relaxation.md` |
| `papers/anisotropic-slip-option0.md` | `notes/desqueezing-…-note.md`; `sdiff-fundamental-vs-emergent.md` |
| `papers/quantum-fluid-instabilities-desi-dr2.md` | `principle-of-vacuum-smoothness.md` |
| `papers/covariance-desi-dr2-full.md` | `ou_bao_stochastic_test.py` |
| `papers/resume.txt` | `papers/data-pack-option0-internet.md`; `results/ou_bao_desi_dr2_run.txt` |
| `theory/euclid-protocol-vacuum-relaxation.md` (hoy en `papers/`) | `scripts/desqueezing/euclid_protocol_forecasts.py` (l.211, 214, 258); `results/euclid_protocol/…`; `results/euclid_mcmc/`, `results/euclid_joint_mcmc/` (l.253); `RETIRED.md`; `euclid_mock_mcmc.py`, `euclid_joint_bao_sne_mcmc.py`; `ou_bao_stochastic_test.py`; `ou_bao_desi_dr2_run.txt`; `figures/euclid_*.png`; `fundamental-vs-emergent-vacuum-relaxation.md` |
| `theory/fundamental-vs-emergent-vacuum-relaxation.md` | `scripts/desqueezing/sdiff_discrimination.py` (l.170); `results/sdiff_discrimination/…`; `figures/sdiff_*.png`; `notes/` (l.192); `euclid-protocol-…`, `sdiff-fundamental-…` (ahora en la misma carpeta, siguen valiendo) |
| `theory/sdiff-fundamental-vs-emergent.md`, `theory/unimodular-gravity-vacuum-smoothness.md` | nombres de otras notas de `theory/` (siguen valiendo); revisar las menciones a `papers/` |
| `exploratory/notes/desqueezing-relaxation-vacuum-fluctuations-note.md` | l.461–465 `scripts/desqueezing/…`, `results/desqueezing_relax_time/`, `results/cosmological_mapping/`; `mapping_tables.md` (misma carpeta, sigue valiendo) |
| `exploratory/notes/mapping_tables.md` | `results/cosmological_mapping/table2_A0_required.csv` |
| `archive/results/euclid_mcmc/RETIRED.md` | l.14 `results/euclid_protocol/` → `exploratory/results/euclid_protocol/` |
| `exploratory/results/sdiff_discrimination/summary.txt` | l.41 `results/sdiff_discrimination/` y `figures/sdiff_*.png` (se regenera con el script) |
| `exploratory/results/termico_vs_squeezed_recomendado/…_resumen.txt`, `exploratory/results/desqueezing_relax_time/…_resumen.txt` | `scripts/desqueezing/…`, `results/…` (se regeneran con el script) |
| `analisis_cmb/NOTAS_ANALISIS_CMB.md` | l.253 `diagnosticos/obsoleto/…` → `../archive/analisis_cmb/diagnosticos/obsoleto/`; `results/ou_bao_desi_dr2_run.txt`; `unimodular-gravity-…` |
| `analisis_cmb/realizaciones.py` | l.18 docstring «kept in diagnosticos/obsoleto/» |
| `analisis_cmb/RESULTADO.md`, `analisis_cmb/verificacion_numeros.md` | `papers/unimodular-gravity-…`, data-pack, euclid-protocol |
| `analisis_cmb/AUDITORIA.md`, `analisis_cmb/bloque_A/COMPARACION.md`, `analisis_cmb/bloque_B/COMPARACION_B1.md`, `COMPARACION_B2.md` | muchas rutas antiguas (`notes/`, `scripts/desqueezing/`, `figures/`, `results/euclid_*`). **Propuesta: no reescribirlas.** Son actas fechadas de lo que se hizo con el layout de entonces; basta con una línea al principio de cada una: «Rutas anteriores a la reorganización; ver `analisis_cmb/REORGANIZACION.md`». |

## (c) README nuevos que hacen falta (en inglés)

| README | Contenido mínimo |
|:--|:--|
| `README.md` (raíz, reescribir) | Mapa de carpetas nuevo (`manuscript/`, `papers/`, `theory/`, `scripts/`, `results/`, `analisis_cmb/`, `exploratory/`, `archive/`); qué está en el claim set C1–C7 y qué no; comandos de reproducción; enlaces corregidos; quitar `notes/` y `figures/` del layout |
| `theory/README.md` | Qué son estas notas (programa teórico: unimodular/SDiff, suavidad del vacío, ν–DE, protocolo Euclid), que **no** forman parte de C1–C7, y puntero al repo `measurable-stochastic-vacuum` |
| `exploratory/README.md` | Programa desqueezing (QuTiP, mapping cosmológico, discriminación SDiff, forecasts Euclid, data pack Option 0); que es exploratorio y no se cita como resultado DESI (N6 / EXPLORATORY_BOUNDARY); estructura `scripts/ results/ notes/`; cómo ejecutar (requiere `qutip`) |
| `archive/README.md` | Qué hay y por qué está retirado: MCMC Euclid sintéticos (N9), stubs, `ou_bao_stochastic_test.py` (DR1 hard-coded, sustituido por `ou_bao_likelihood.py`), correlación cruzada preliminar (N5), figuras sin generador (con commit de origen), diagnósticos CMB obsoletos. Aviso: no se mantiene y puede no ejecutar |
| `analisis_cmb/README.md` (ya existe, en español → traducir) | Análisis CMB/Planck del OU (realizaciones, malla Takahashi, everpresent, oscilador), auditoría y bloques A/B; aclarar si entra o no en el preprint |
