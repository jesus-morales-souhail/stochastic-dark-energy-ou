# CORRECCION_FILA_SN: impacto de la errata 129 (fila SN de la respuesta lineal A) en todo lo que depende de `build_A`

Exploración (no toca PREPRINT, CLAIMS ni las cifras de PROPUESTA_V3). 2026-11-13.

**Decisión del autor (12-11-2026):** la errata de la fila SN de A (fila 129 de `manuscript/ERRATAS_V3.md`) queda **confirmada** por el criterio
fijado antes (A = jacobiano exacto del modelo directo; gana la A del reproductor, `REPRODUCCION.md` §'P1: la fila SN de A', commit dc3a98d).
Ahora se mide su impacto en todo lo que depende de `build_A`. Todavía NO se cambian cifras de PROPUESTA_V3, PREPRINT ni CLAIMS.

## Predicciones (de Claude, escritas antes de ejecutar, copiadas literales)

- P7. Las 8 cotas OU de cabecera (θ = 0.3/1/3/10, sin y con ρ_DE > 0) cambian < 1e-3 relativo con la A corregida.
- P8. Familias (todas las familias y ℓ calculados) y muestras de SN (Pantheon+, DES-SN5YR, Dovekie, Union3): cada cota cambia < 1e-2 relativo y cada máximo de ln Z − ln Z(0) cambia < 0.05.
- P9. INYECCION, mismas semillas: la cobertura de la marginal −1.92 cambia como mucho en 1 universo por escenario; el 43 % de falsas exclusiones del escenario de tensión cambia como mucho 2 puntos.
- P10. En el modelo directo exacto, la interpolación lineal en u de D̃_M en los redshifts de las 1580 SN (exacto_sn.py:147) se aparta del cálculo directo (integrando hasta cada z_SN) menos de 1e-3 mag en μ en todas, en p0 y X = 0.

Reglas seguidas: `build_A` y los json existentes no se han modificado; todo lo nuevo está en `respuesta_v3A.py` (build_A_v3, make_linear_v3),
`interpolacion_mu_check.py/json`, `correccion_fila_sn_runner.py`, `correccion_fila_sn_tabla.py`, `CORRECCION_FILA_SN_v3A.json` y los json
`*_v3A.json`; mismo número de hilos BLAS que los originales (cada script se importa primero en un subproceso, de modo que su propio
ajuste de hilos se aplica como en la corrida original: positividad_lin con el valor por defecto, familias/muestras_sn/inyeccion con 1 por proceso).

## 1. Inventario: scripts que construyen o usan la A lineal, y qué producen

| script | cómo obtiene A | fila SN | salida | recalculado aquí |
|---|---|---|---|---|
| `marginal_guiada.py` | `build_A` (l. 113–149) | antigua (l. 147: interpola ∂ln D̃_M por nodo) | `marginal_guiada.json` (lineal y exacta guiada) | no: la lineal de θ = 1 está cubierta por positividad_lin; la exacta guiada no se recalcula (§5) |
| `positividad_lin.py` | `MG.build_A` | antigua | `positividad_lin.json` (cotas lineales con y sin ρ_DE > 0, 4 θ) | **sí** → `positividad_lin_v3A.json` |
| `positividad.py` (v2, exacta guiada) | `MG.build_A` | antigua | `positividad.json` | no (§5) |
| `familias.py` | `init_linear` → `MG.build_A` (l. 310–311); `_linear` (l. 270); `guided_cell` (l. 354) | antigua | `familias.json` → linear (20 celdas), guided | **sí, la parte lineal** → `familias_linear_v3A.json`; la guiada no (§5) |
| `muestras_sn.py` | `make_linear` (l. 160–206) construye A por muestra, misma fórmula | antigua | `muestras_sn.json` → linear (32 celdas), drho01, ratios | **sí, la parte lineal** → `muestras_sn_linear_v3A.json` |
| `inyeccion.py` | `make_full` → `MS.make_linear`; `make_binned` comprime esa A (l. 90–120) | antigua | `inyeccion.json/txt`, `inyeccion_check.json` | **sí, todo** → `inyeccion_v3A.json/txt` (el `inyeccion_check.json` solo se escribe en modo `--check`) |
| `cmb_final.py`, `cruces_w.py` | `MG.build_A` (muestreo de Matheron del posterior lineal) | antigua | `cmb_final*.json/npz`, `cruces_w.json` | no (§5) |
| `modos.py` | `MG.build_A` (A canónica, Fisher por tramos; fase 1) | antigua | `modos_pred.json` | no: fase 1 parada por diseño (MODOS.md); fuera del encargo |
| `exacto_sn.py` `linear_bounds` (l. 272–298) y `maximos_globales.py` (l. 90–115) | A propia con la **misma** interpolación de la derivada por nodo (l. 291–298 y 115) | antigua | `exacto_sn.json` → linear; `maximos_globales.json` → exacto_sn | **no (fuera del encargo): afectados; cifras de SN solas y BAO + SN de PROPUESTA_V3 (A34) pendientes de recalcular** |
| `cmb_fondo.py`, `exacto.py` | respuesta sin filas SN / sin A | — | `cmb_fondo.json`, `exacto.json`, `maximos_globales.json` → cmb_fondo | no afectados |
| `inyeccion_percentiles.py`, `fisher_tramos_pantheon.py` | leen `inyeccion.json` / usan solo Bm | — | — | no afectados por A (Fisher usa Bm); los percentiles cambiarían con `inyeccion_v3A.json` |
| `comparacion_reproduccion.py`, `jacobiano_sn_check.py`, `sensibilidad_fila_sn.py`, `marginal_guiada_check.py`, `marginal_guiada_cf.py` | diagnósticos | — | — | no |

## 2. `build_A_v3` (`respuesta_v3A.py`, comprobación en `respuesta_v3A_check.json`)

Definición fijada: ∂y_SN/∂X_j = (5/ln 10)·[(1 − t) ∂D̃_M(u_k) + t ∂D̃_M(u_{k+1})] / D̃_M(u_HD); las demás filas y el convenio B idénticos a `build_A`
(filas 0–14 iguales a 1.8e-15). Comprobación: coincide con el `A.npy` del reproductor a 8.2e-09 (todas las filas; filas SN 7.9e-10),
exigido < 10⁻⁸; con el jacobiano numérico de `jacobiano_sn_check.json` a 5.5e-07 en los 44 nodos convergidos, exigido < 2 × 10⁻⁶. Filas suman cero a
1.3e-15. Diferencia máxima con la fila SN antigua: 4.09e-03. (Al repetir el jacobiano con otro número de hilos, 44 de 46 nodos
convergen en vez de 45: el criterio está en el límite del redondeo; la comparación usa h/10 en los convergidos.)

## 3. P10: interpolación frente a integración directa (`interpolacion_mu_check.json`)

En p₀ y X = 0, μ de las 1580 SN por interpolación lineal en u de D̃_M (como `exacto_sn.py:147`) frente a la integral directa de (1 + z)/E(u)
hasta cada u_HD (quad, tolerancia relativa 10⁻¹²): máximo 5.23e-05 mag (en z_HD = 0.01096, la SN de menor z, a mitad de intervalo), mediana
9.7e-07 mag, p95 2.7e-05 mag. El trapecio de la malla en los nodos BAO difiere de la integral exacta en 3 × 10⁻⁷ relativo. **P10: sí** (< 10⁻³ mag): el
problema estaba solo en la derivada, no en el modelo directo.

## 4. Reproducción con la A antigua y tabla antiguo / nuevo

Antes de recalcular, cada script se corrió con la A antigua en el mismo entorno y reprodujo su json original **bit a bit** (diferencia máxima 0.0 en
positividad_lin, familias (lineal), muestras_sn (lineal) e inyeccion; `v3A_runs_repro/*_oldA_repro.json`). Después, con `build_A_v3` /
`make_linear_v3` (`*_v3A.json`). "cambio" es relativo para cotas y cocientes, absoluto para ln Z y fracciones.

### 4.1 positividad_lin (cotas de cabecera)

| cifra | antiguo | nuevo | cambio |
|---|---|---|---|
| theta_0.3 lineal −1.92 | 0.4916 | 0.4916 | -3.31e-05 |
| theta_0.3 lineal −1.35 | 0.4266 | 0.4266 | -1.19e-05 |
| theta_0.3 con ρ_DE > 0, −1.92 | 0.4606 | 0.4605 | -8.38e-05 |
| theta_0.3 σ del máximo | 0.2000 | 0.2000 | +0.00e+00 |
| theta_0.3 máx ln Z − ln Z(0) | 1.1651 | 1.1653 | +1.59e-04 (abs) |
| theta_0.3 ρ_DE ≤ 0 en el techo lineal | 0.5957 | 0.5959 | +1.62e-04 (abs) |
| theta_1 lineal −1.92 | 0.2719 | 0.2719 | -4.17e-05 |
| theta_1 lineal −1.35 | 0.2345 | 0.2345 | -1.79e-05 |
| theta_1 con ρ_DE > 0, −1.92 | 0.2625 | 0.2625 | -4.59e-05 |
| theta_1 σ del máximo | 0.1000 | 0.1000 | +0.00e+00 |
| theta_1 máx ln Z − ln Z(0) | 1.1394 | 1.1395 | +1.34e-04 (abs) |
| theta_1 ρ_DE ≤ 0 en el techo lineal | 0.4106 | 0.4106 | -3.79e-05 (abs) |
| theta_3 lineal −1.92 | 0.1727 | 0.1727 | -6.64e-05 |
| theta_3 lineal −1.35 | 0.1517 | 0.1517 | -3.62e-05 |
| theta_3 con ρ_DE > 0, −1.92 | 0.1741 | 0.1741 | -1.47e-04 |
| theta_3 σ del máximo | 0.0700 | 0.0700 | +0.00e+00 |
| theta_3 máx ln Z − ln Z(0) | 0.9624 | 0.9626 | +1.68e-04 (abs) |
| theta_3 ρ_DE ≤ 0 en el techo lineal | 0.0989 | 0.0990 | +1.56e-04 (abs) |
| theta_10 lineal −1.92 | 0.1151 | 0.1151 | -1.34e-04 |
| theta_10 lineal −1.35 | 0.1002 | 0.1002 | -8.15e-05 |
| theta_10 con ρ_DE > 0, −1.92 | 0.1159 | 0.1159 | -1.46e-04 |
| theta_10 σ del máximo | 0.0300 | 0.0300 | +0.00e+00 |
| theta_10 máx ln Z − ln Z(0) | 0.5898 | 0.5900 | +2.04e-04 (abs) |
| theta_10 ρ_DE ≤ 0 en el techo lineal | 0.0031 | 0.0031 | -2.98e-06 (abs) |
| P_prior(ρ_DE>0 | σ = cota pos θ = 1) | 0.7320 | 0.7320 | +0.00e+00 (abs) |

**P7: sí** (cambio relativo máximo de las 8 cotas 1.5e-04 < 10⁻³).

### 4.2 familias (parte lineal, 20 celdas)

| cifra | antiguo | nuevo | cambio |
|---|---|---|---|
| OU|l_0.1 192 | 0.1179 | 0.1179 | -1.45e-04 |
| OU|l_0.1 135 | 0.1003 | 0.1003 | -1.10e-04 |
| OU|l_0.1 máx ln Z | 0.5904 | 0.5906 | +2.14e-04 (abs) |
| OU|l_0.3 192 | 0.1675 | 0.1675 | -6.35e-05 |
| OU|l_0.3 135 | 0.1453 | 0.1453 | -4.13e-05 |
| OU|l_0.3 máx ln Z | 0.9547 | 0.9548 | +1.53e-04 (abs) |
| OU|l_1 192 | 0.2772 | 0.2772 | -3.50e-05 |
| OU|l_1 135 | 0.2417 | 0.2417 | -1.65e-05 |
| OU|l_1 máx ln Z | 1.1402 | 1.1403 | +1.47e-04 (abs) |
| OU|l_3 192 | 0.4676 | 0.4676 | -2.79e-05 |
| OU|l_3 135 | 0.4081 | 0.4081 | -1.00e-05 |
| OU|l_3 máx ln Z | 1.1677 | 1.1678 | +1.53e-04 (abs) |
| OU|l_10 192 | 0.8476 | 0.8476 | -2.32e-05 |
| OU|l_10 135 | 0.7398 | 0.7398 | -5.08e-06 |
| OU|l_10 máx ln Z | 1.1653 | 1.1655 | +1.47e-04 (abs) |
| Matern32|l_0.1 192 | 0.1259 | 0.1259 | +9.60e-05 |
| Matern32|l_0.1 135 | 0.1082 | 0.1082 | +1.02e-04 |
| Matern32|l_0.1 máx ln Z | 0.7017 | 0.7018 | +1.33e-04 (abs) |
| Matern32|l_0.3 192 | 0.2475 | 0.2475 | +6.34e-05 |
| Matern32|l_0.3 135 | 0.2137 | 0.2137 | +5.90e-05 |
| Matern32|l_0.3 máx ln Z | 1.4025 | 1.4025 | -2.52e-05 (abs) |
| Matern32|l_1 192 | 0.7370 | 0.7370 | +2.56e-05 |
| Matern32|l_1 135 | 0.6074 | 0.6074 | +1.70e-05 |
| Matern32|l_1 máx ln Z | 1.5029 | 1.5028 | -4.14e-05 (abs) |
| Matern32|l_3 192 | — | — | — |
| Matern32|l_3 135 | — | — | — |
| Matern32|l_3 máx ln Z | 1.2707 | 1.2706 | -5.24e-05 (abs) |
| Matern32|l_10 192 | — | — | — |
| Matern32|l_10 135 | — | — | — |
| Matern32|l_10 máx ln Z | 0.0000 | 0.0000 | +0.00e+00 (abs) |
| Gauss|l_0.1 192 | 0.1545 | 0.1545 | +6.35e-05 |
| Gauss|l_0.1 135 | 0.1325 | 0.1325 | +6.14e-05 |
| Gauss|l_0.1 máx ln Z | 0.8914 | 0.8914 | -7.30e-06 (abs) |
| Gauss|l_0.3 192 | 0.3348 | 0.3348 | -2.38e-06 |
| Gauss|l_0.3 135 | 0.2690 | 0.2690 | +8.99e-06 |
| Gauss|l_0.3 máx ln Z | 1.3793 | 1.3793 | -9.30e-06 (abs) |
| Gauss|l_1 192 | — | — | — |
| Gauss|l_1 135 | — | — | — |
| Gauss|l_1 máx ln Z | 1.3311 | 1.3310 | -9.38e-05 (abs) |
| Gauss|l_3 192 | — | — | — |
| Gauss|l_3 135 | — | — | — |
| Gauss|l_3 máx ln Z | 0.0000 | 0.0000 | +0.00e+00 (abs) |
| Gauss|l_10 192 | — | — | — |
| Gauss|l_10 135 | — | — | — |
| Gauss|l_10 máx ln Z | 0.0000 | 0.0000 | +0.00e+00 (abs) |
| Steps|l_0.1 192 | 0.0969 | 0.0969 | -2.15e-05 |
| Steps|l_0.1 135 | 0.0806 | 0.0806 | +2.22e-05 |
| Steps|l_0.1 máx ln Z | 0.3750 | 0.3752 | +1.73e-04 (abs) |
| Steps|l_0.3 192 | 0.1548 | 0.1548 | -3.76e-05 |
| Steps|l_0.3 135 | 0.1339 | 0.1339 | -1.57e-05 |
| Steps|l_0.3 máx ln Z | 1.0098 | 1.0100 | +1.62e-04 (abs) |
| Steps|l_1 192 | 0.2681 | 0.2681 | -2.66e-05 |
| Steps|l_1 135 | 0.2339 | 0.2339 | -8.75e-06 |
| Steps|l_1 máx ln Z | 1.2383 | 1.2384 | +1.18e-04 (abs) |
| Steps|l_3 192 | 0.4629 | 0.4629 | -2.51e-05 |
| Steps|l_3 135 | 0.4041 | 0.4041 | -7.51e-06 |
| Steps|l_3 máx ln Z | 1.1624 | 1.1626 | +1.62e-04 (abs) |
| Steps|l_10 192 | 0.8453 | 0.8453 | -2.24e-05 |
| Steps|l_10 135 | 0.7379 | 0.7379 | -4.31e-06 |
| Steps|l_10 máx ln Z | 1.1617 | 1.1619 | +1.50e-04 (abs) |

### 4.3 muestras_sn (parte lineal, 32 celdas)

| cifra | antiguo | nuevo | cambio |
|---|---|---|---|
| pantheon|OU|theta_0.3 192 | 0.4924 | 0.4924 | -2.55e-05 |
| pantheon|OU|theta_0.3 135 | 0.4297 | 0.4297 | -7.24e-06 |
| pantheon|OU|theta_0.3 máx ln Z | 1.1680 | 1.1681 | +1.45e-04 (abs) |
| pantheon|OU|theta_1 192 | 0.2772 | 0.2772 | -3.50e-05 |
| pantheon|OU|theta_1 135 | 0.2417 | 0.2417 | -1.65e-05 |
| pantheon|OU|theta_1 máx ln Z | 1.1402 | 1.1403 | +1.47e-04 (abs) |
| pantheon|OU|theta_3 192 | 0.1743 | 0.1742 | -5.74e-05 |
| pantheon|OU|theta_3 135 | 0.1513 | 0.1513 | -3.53e-05 |
| pantheon|OU|theta_3 máx ln Z | 0.9812 | 0.9813 | +1.43e-04 (abs) |
| pantheon|OU|theta_10 192 | 0.1179 | 0.1179 | -1.45e-04 |
| pantheon|OU|theta_10 135 | 0.1003 | 0.1003 | -1.10e-04 |
| pantheon|OU|theta_10 máx ln Z | 0.5904 | 0.5906 | +2.14e-04 (abs) |
| pantheon|Matern32|l_0.3 192 | 0.2475 | 0.2475 | +6.34e-05 |
| pantheon|Matern32|l_0.3 135 | 0.2137 | 0.2137 | +5.90e-05 |
| pantheon|Matern32|l_0.3 máx ln Z | 1.4025 | 1.4025 | -2.52e-05 (abs) |
| pantheon|Matern32|l_1 192 | 0.7370 | 0.7370 | +2.56e-05 |
| pantheon|Matern32|l_1 135 | 0.6074 | 0.6074 | +1.70e-05 |
| pantheon|Matern32|l_1 máx ln Z | 1.5029 | 1.5028 | -4.14e-05 (abs) |
| pantheon|Gauss|l_0.3 192 | 0.3348 | 0.3348 | -2.38e-06 |
| pantheon|Gauss|l_0.3 135 | 0.2690 | 0.2690 | +8.99e-06 |
| pantheon|Gauss|l_0.3 máx ln Z | 1.3793 | 1.3793 | -9.30e-06 (abs) |
| pantheon|Gauss|l_1 192 | — | — | — |
| pantheon|Gauss|l_1 135 | — | — | — |
| pantheon|Gauss|l_1 máx ln Z | 1.3311 | 1.3310 | -9.38e-05 (abs) |
| des5yr|OU|theta_0.3 192 | 0.5158 | 0.5159 | +8.23e-05 |
| des5yr|OU|theta_0.3 135 | 0.4476 | 0.4476 | +7.71e-05 |
| des5yr|OU|theta_0.3 máx ln Z | 4.8570 | 4.8568 | -2.12e-04 (abs) |
| des5yr|OU|theta_1 192 | 0.2896 | 0.2896 | +8.69e-05 |
| des5yr|OU|theta_1 135 | 0.2511 | 0.2511 | +8.22e-05 |
| des5yr|OU|theta_1 máx ln Z | 4.8569 | 4.8567 | -2.17e-04 (abs) |
| des5yr|OU|theta_3 192 | 0.1819 | 0.1820 | +9.75e-05 |
| des5yr|OU|theta_3 135 | 0.1574 | 0.1574 | +9.32e-05 |
| des5yr|OU|theta_3 máx ln Z | 4.7646 | 4.7644 | -2.10e-04 (abs) |
| des5yr|OU|theta_10 192 | 0.1249 | 0.1249 | +1.25e-04 |
| des5yr|OU|theta_10 135 | 0.1076 | 0.1076 | +1.19e-04 |
| des5yr|OU|theta_10 máx ln Z | 4.5929 | 4.5928 | -1.80e-04 (abs) |
| des5yr|Matern32|l_0.3 192 | 0.2542 | 0.2543 | +2.57e-05 |
| des5yr|Matern32|l_0.3 135 | 0.2167 | 0.2167 | +2.10e-05 |
| des5yr|Matern32|l_0.3 máx ln Z | 4.7648 | 4.7645 | -2.23e-04 (abs) |
| des5yr|Matern32|l_1 192 | 0.6952 | 0.6952 | +8.08e-08 |
| des5yr|Matern32|l_1 135 | 0.5655 | 0.5655 | -3.71e-06 |
| des5yr|Matern32|l_1 máx ln Z | 5.0112 | 5.0110 | -2.06e-04 (abs) |
| des5yr|Gauss|l_0.3 192 | 0.3198 | 0.3198 | -2.30e-05 |
| des5yr|Gauss|l_0.3 135 | 0.2522 | 0.2522 | -1.80e-05 |
| des5yr|Gauss|l_0.3 máx ln Z | 4.8238 | 4.8236 | -1.73e-04 (abs) |
| des5yr|Gauss|l_1 192 | — | — | — |
| des5yr|Gauss|l_1 135 | — | — | — |
| des5yr|Gauss|l_1 máx ln Z | 5.3525 | 5.3522 | -2.26e-04 (abs) |
| dovekie|OU|theta_0.3 192 | 0.4905 | 0.4905 | +6.63e-05 |
| dovekie|OU|theta_0.3 135 | 0.4226 | 0.4226 | +5.69e-05 |
| dovekie|OU|theta_0.3 máx ln Z | 1.8401 | 1.8398 | -2.86e-04 (abs) |
| dovekie|OU|theta_1 192 | 0.2753 | 0.2754 | +7.05e-05 |
| dovekie|OU|theta_1 135 | 0.2369 | 0.2369 | +6.12e-05 |
| dovekie|OU|theta_1 máx ln Z | 1.8319 | 1.8316 | -2.89e-04 (abs) |
| dovekie|OU|theta_3 192 | 0.1714 | 0.1714 | +8.46e-05 |
| dovekie|OU|theta_3 135 | 0.1468 | 0.1468 | +7.53e-05 |
| dovekie|OU|theta_3 máx ln Z | 1.7687 | 1.7684 | -2.94e-04 (abs) |
| dovekie|OU|theta_10 192 | 0.1137 | 0.1137 | +1.29e-04 |
| dovekie|OU|theta_10 135 | 0.0962 | 0.0962 | +1.18e-04 |
| dovekie|OU|theta_10 máx ln Z | 1.7412 | 1.7409 | -3.11e-04 (abs) |
| dovekie|Matern32|l_0.3 192 | 0.2417 | 0.2417 | -1.91e-05 |
| dovekie|Matern32|l_0.3 135 | 0.2046 | 0.2046 | -3.00e-05 |
| dovekie|Matern32|l_0.3 máx ln Z | 1.7007 | 1.7005 | -2.46e-04 (abs) |
| dovekie|Matern32|l_1 192 | 0.6446 | 0.6446 | -4.87e-05 |
| dovekie|Matern32|l_1 135 | 0.5162 | 0.5162 | -5.06e-05 |
| dovekie|Matern32|l_1 máx ln Z | 1.9498 | 1.9496 | -2.03e-04 (abs) |
| dovekie|Gauss|l_0.3 192 | 0.2837 | 0.2837 | -5.08e-05 |
| dovekie|Gauss|l_0.3 135 | 0.2201 | 0.2201 | -3.99e-05 |
| dovekie|Gauss|l_0.3 máx ln Z | 1.7845 | 1.7843 | -1.62e-04 (abs) |
| dovekie|Gauss|l_1 192 | — | — | — |
| dovekie|Gauss|l_1 135 | — | — | — |
| dovekie|Gauss|l_1 máx ln Z | 2.2199 | 2.2197 | -2.18e-04 (abs) |
| union3|OU|theta_0.3 192 | 0.6872 | 0.6872 | +2.97e-05 |
| union3|OU|theta_0.3 135 | 0.6025 | 0.6025 | +2.78e-05 |
| union3|OU|theta_0.3 máx ln Z | 3.5929 | 3.5930 | +1.02e-04 (abs) |
| union3|OU|theta_1 192 | 0.3866 | 0.3866 | +3.10e-05 |
| union3|OU|theta_1 135 | 0.3393 | 0.3393 | +2.94e-05 |
| union3|OU|theta_1 máx ln Z | 3.6089 | 3.6090 | +1.00e-04 (abs) |
| union3|OU|theta_3 192 | 0.2449 | 0.2449 | +3.39e-05 |
| union3|OU|theta_3 135 | 0.2157 | 0.2157 | +3.27e-05 |
| union3|OU|theta_3 máx ln Z | 3.5496 | 3.5497 | +1.03e-04 (abs) |
| union3|OU|theta_10 192 | 0.1720 | 0.1720 | +4.15e-05 |
| union3|OU|theta_10 135 | 0.1517 | 0.1517 | +4.14e-05 |
| union3|OU|theta_10 máx ln Z | 3.3176 | 3.3177 | +1.13e-04 (abs) |
| union3|Matern32|l_0.3 192 | 0.3435 | 0.3435 | +2.12e-05 |
| union3|Matern32|l_0.3 135 | 0.2903 | 0.2903 | +2.32e-05 |
| union3|Matern32|l_0.3 máx ln Z | 3.4024 | 3.4025 | +1.05e-04 (abs) |
| union3|Matern32|l_1 192 | 0.8283 | 0.8284 | +5.21e-05 |
| union3|Matern32|l_1 135 | 0.6615 | 0.6615 | +5.05e-05 |
| union3|Matern32|l_1 máx ln Z | 3.4090 | 3.4091 | +5.07e-05 (abs) |
| union3|Gauss|l_0.3 192 | 0.3232 | 0.3232 | +3.07e-05 |
| union3|Gauss|l_0.3 135 | 0.2626 | 0.2626 | +2.77e-05 |
| union3|Gauss|l_0.3 máx ln Z | 3.3025 | 3.3025 | +4.89e-05 (abs) |
| union3|Gauss|l_1 192 | — | — | — |
| union3|Gauss|l_1 135 | — | — | — |
| union3|Gauss|l_1 máx ln Z | 3.5484 | 3.5485 | +2.94e-05 (abs) |

**P8: sí** (cotas: cambio relativo máximo 1.4e-04 en familias y 1.4e-04 en muestras, < 10⁻²; máximos de ln Z: 2.1e-04 y
3.1e-04, < 0.05). Las celdas sin cota siguen sin cota; las clases de máximo no cambian.

### 4.4 inyeccion (todos los escenarios, mismas semillas)

| cifra | antiguo | nuevo | cambio |
|---|---|---|---|
| A marginal −1.92 mediana | 0.1095 | 0.1095 | +8.99e-05 |
| A 5 % mediana | 0.0948 | 0.0948 | +1.79e-04 |
| A σ_X = 0 excluido | 0.0000 | 0.0000 | +0.00e+00 (abs) |
| A marginal /A | 1.0000 | 1.0000 | +0.00e+00 |
| A Δln Z máx mediana | 0.0000 | 0.0000 | +0.00e+00 (abs) |
| B_0.1 marginal −1.92 mediana | 0.2217 | 0.2218 | +8.48e-05 |
| B_0.1 5 % mediana | 0.1949 | 0.1950 | +6.43e-04 |
| B_0.1 cobertura marginal | 0.9550 | 0.9550 | +0.00e+00 (abs) |
| B_0.1 cobertura 5 % | 0.9250 | 0.9250 | +0.00e+00 (abs) |
| B_0.1 σ_X = 0 excluido | 0.5850 | 0.5850 | +0.00e+00 (abs) |
| B_0.1 marginal /A | 2.0253 | 2.0252 | -5.16e-06 |
| B_0.1 Δln Z máx mediana | 2.8483 | 2.8479 | -3.44e-04 (abs) |
| B_0.2 marginal −1.92 mediana | 0.3756 | 0.3756 | +1.79e-04 |
| B_0.2 5 % mediana | 0.3478 | 0.3479 | +3.69e-04 |
| B_0.2 cobertura marginal | 0.9650 | 0.9650 | +0.00e+00 (abs) |
| B_0.2 cobertura 5 % | 0.8650 | 0.8650 | +0.00e+00 (abs) |
| B_0.2 σ_X = 0 excluido | 0.9450 | 0.9450 | +0.00e+00 (abs) |
| B_0.2 marginal /A | 3.4303 | 3.4306 | +8.92e-05 |
| B_0.2 Δln Z máx mediana | 17.6616 | 17.6612 | -3.41e-04 (abs) |
| B_0.3 marginal −1.92 mediana | 0.5222 | 0.5224 | +3.45e-04 |
| B_0.3 5 % mediana | 0.5066 | 0.5067 | +3.39e-04 |
| B_0.3 cobertura marginal | 0.9750 | 0.9750 | +0.00e+00 (abs) |
| B_0.3 cobertura 5 % | 0.8500 | 0.8500 | +0.00e+00 (abs) |
| B_0.3 σ_X = 0 excluido | 0.9950 | 0.9950 | +0.00e+00 (abs) |
| B_0.3 marginal /A | 4.7698 | 4.7710 | +2.55e-04 |
| B_0.3 Δln Z máx mediana | 46.6872 | 46.6868 | -3.70e-04 (abs) |
| C marginal −1.92 mediana | 0.1646 | 0.1647 | +2.36e-04 |
| C 5 % mediana | 0.1624 | 0.1624 | +3.52e-04 |
| C σ_X = 0 excluido | 0.4300 | 0.4300 | +0.00e+00 (abs) |
| C marginal /A | 1.5037 | 1.5040 | +1.46e-04 |
| C Δln Z máx mediana | 1.4647 | 1.4648 | +1.01e-04 (abs) |
| D marginal −1.92 mediana | 0.1111 | 0.1111 | +5.83e-05 |
| D 5 % mediana | 0.0984 | 0.0984 | +1.57e-04 |
| D σ_X = 0 excluido | 0.0000 | 0.0000 | +0.00e+00 (abs) |
| D marginal /A | 1.0151 | 1.0150 | -3.16e-05 |
| D Δln Z máx mediana | 0.0000 | 0.0000 | +0.00e+00 (abs) |
| datos reales marg_linear | 0.2771 | 0.2771 | -3.37e-05 |
| datos reales 5pct_linear | 0.1559 | 0.1559 | +2.30e-05 |
| datos reales dlnL_max_linear | 1.1413 | 1.1414 | +1.24e-04 |
| comprobación exacta: cociente mín/máx | 0.8847 / 1.1018 | 0.8847 / 1.1018 | — |
| δ del escenario C | 0.0512 | 0.0512 | +0.00e+00 |

**P9: sí** (coberturas de la marginal idénticas en los tres escenarios B: 0 universos de cambio; falsas exclusiones del escenario C 43 % → 43 %, 0 puntos).

## 5. Qué no se recalculó y por qué (nivel: asumido)

- **Marginal exacta guiada** (`marginal_guiada.json`, `positividad.json`, parte guiada de `familias.json`): A entra solo en la propuesta de muestreo (Matheron) y en el
  χ² lineal del peso; los pesos exacto/lineal corrigen la propuesta, así que el estimador es consistente con cualquier A razonable; cambia su varianza (ESS), no su
  valor esperado. Comprobado en el código: en `marginal_guiada.cell` A aparece en `x = x0 + σ² K Aᵀ S⁻¹ R` y en `chi_lin`; no en `chi_ex`.
- **CMB_FINAL y cruces_w**: las historias son muestras del posterior lineal; cambiar A cambia las 120 historias y exige volver a correr CAMB. En el código A entra
  igual que en la marginal guiada (muestreo), más los pesos de información; no entra en el CMB. Esperable: cambios de orden 10⁻⁴ en las historias.
- **EXACTO / EXACTO_SN (regla del 5 %)**: no usan A. **Pero** `exacto_sn.linear_bounds` y `maximos_globales.py` sí construyen una A propia con la misma
  interpolación de la derivada por nodo (l. 291–298 y 115): están afectados y no se han recalculado (fuera del encargo); alimentan las cifras de SN solas,
  BAO + SN y BAO + θ* + ω_m + SN de PROPUESTA_V3 (A34, A26). El tamaño esperable del cambio, por lo visto aquí, es ≤ 2 × 10⁻⁴ relativo.
- **modos.py**: fase 1 parada por diseño; usa A pero no produce cifras de la propuesta.

## 6. Cifras de PROPUESTA_V3 que cambiarían a la precisión con la que se citan (no se cambian aquí)

Con 3 decimales en cotas y 1 decimal en porcentajes, **ninguna cifra de PROPUESTA_V3 cambia**: las 8 cotas de cabecera (0.492 / 0.272 / 0.173 / 0.115 y
0.461 / 0.263 / 0.174 / 0.116), las de familias y muestras citadas, las coberturas (95.5 / 96.5 / 97.5 %), el 43 %, el ×1.50 del escenario C, el 41.1 % de
ρ_DE ≤ 0 en el techo lineal y el 0.732 de P_prior quedan iguales. Cifras que sí cambiarían en sus notas de origen, a su precisión: MUESTRAS_SN.md Union3
OU θ = 0.3, −1.35: 0.602 → 0.603; INYECCION.md "marginal /A" de B 0.2 y B 0.3 a tres decimales (3.430 → 3.431; 4.770 → 4.771). Pendiente de medir: las cifras
de `exacto_sn.linear_bounds` / `maximos_globales` (SN solas 0.246 / 0.147 / 0.092 y los cocientes ×1.85 / ×1.98 / ×2.22 de A34 [corregido 13-11-2026: es A35 y cita ×1.85 / ×1.93 / ×2.21; error de transcripción de esta nota, no de PROPUESTA_V3]) cuando se recalculen con la
fila SN corregida.

## Predicción, punto por punto

P7 sí (1.5e-04). P8 sí (1.4e-04 en cotas, 3.1e-04 en ln Z). P9 sí (0 universos, 0 puntos). P10 sí (5.2e-05 mag).

## 7. Cierre: `exacto_sn.linear_bounds` y `maximos_globales.py` (13-11-2026)

**Predicción (de Claude, escrita antes de ejecutar, copiada literal):**

- P11. Con la A corregida, las cotas de SN solas (0.246 / 0.147 / 0.092) y los cocientes de A34 (×1.85 / ×1.98 / ×2.22) no cambian a la precisión con la que se citan, ni las cifras de A26; los máximos globales siguen en el mismo nodo y su ln Z − ln Z(0) cambia < 1e-3.

**Cómo se hizo.** `exacto_sn.linear_bounds` (l. 272–363) y la parte `exacto_sn` de `maximos_globales.py` (l. 90–178) construyen su propia A en la malla
uniforme de `cmb_fondo` (fiducial `CF.OM_FID`), no en la de `build_A`, así que no pueden llamar a `build_A_v3` tal cual: el módulo nuevo
`exacto_sn_lineal_v3A.py` copia las dos rutinas sin cambios y aplica en su malla la **misma definición** de la fila SN de `build_A_v3`
(derivada del logaritmo de la D̃_M interpolada; modo `v3A`), o la fila antigua (modo `old`, reproducción). Hilos BLAS como los originales:
1 para `exacto_sn` (su `setdefault`), 16 para `maximos_globales` (importa numpy antes que los módulos que fijan hilos). Salidas:
`exacto_sn_lineal_v3A_{exacto_sn,maximos}_{old,v3A}.json` (más `_old_omp1` y `_old_omp16` para la prueba de hilos) y la tabla en
`CORRECCION_FILA_SN_cierre_v3A.json` (`correccion_fila_sn_cierre.py`).

**Paso 1, reproducción con la fila antigua.** `exacto_sn.linear_bounds`: **bit a bit** (diferencia máxima 0.0 en las 30 cotas, χ² y parámetros).
`maximos_globales` (parte `exacto_sn`): ln Z − ln Z(0) del máximo global reproducido a ≤ 5.7 × 10⁻⁹ y, en SN solas, exactamente (0.0); las cotas
de BAO + SN y BAO + θ* + ω_m + SN a ≤ 1.5 × 10⁻¹⁰; **pero las raíces de las cotas de SN solas difieren de las originales en hasta 9.4 × 10⁻⁸
(1.2 × 10⁻⁷ relativo) con 16 hilos y hasta 5.2 × 10⁻⁷ con 1 hilo, por encima del ≤ 10⁻⁸ fijado.** Con ln Z idéntica en los dos recuentos de hilos
y `exacto_sn.linear_bounds` bit a bit con la misma A, la causa es el orden de suma dentro de `brentq` (evaluaciones de ln Z fuera de la malla
que difieren en ~10⁻⁹) amplificado por el perfil casi plano de SN solas cerca de la raíz; no es un cambio de método ni de datos. Se deja
constancia como desviación del criterio y se sigue con esa salvedad:

| entrada | \|16 hilos − original\| cota 192 | \|1 hilo − original\| cota 192 | \|16 h − orig.\| ln Z máx | \|1 h − orig.\| ln Z máx |
|---|---|---|---|---|
| SN theta_0.1 | 9.4e-08 | 5.2e-07 | 0.0e+00 | 0.0e+00 |
| SN theta_0.3 | 5.4e-08 | 3.0e-07 | 0.0e+00 | 0.0e+00 |
| SN theta_1 | 2.9e-08 | 1.6e-07 | 0.0e+00 | 0.0e+00 |
| SN theta_3 | 1.6e-08 | 8.7e-08 | 0.0e+00 | 0.0e+00 |
| SN theta_10 | 9.0e-09 | 5.0e-08 | 0.0e+00 | 0.0e+00 |
| BAO+SN theta_0.1 | 2.3e-12 | 8.3e-11 | 1.7e-10 | 3.1e-10 |
| BAO+SN theta_0.3 | 1.0e-12 | 4.6e-11 | 1.7e-10 | 3.1e-10 |
| BAO+SN theta_1 | 1.7e-13 | 2.2e-11 | 1.6e-10 | 3.3e-10 |
| BAO+SN theta_3 | 9.0e-13 | 1.3e-11 | 1.3e-10 | 3.0e-10 |
| BAO+SN theta_10 | 1.5e-12 | 8.6e-12 | 1.1e-10 | 3.8e-10 |
| BAO+CMB+SN theta_0.1 | 1.3e-10 | 9.5e-11 | 3.9e-10 | 4.9e-10 |
| BAO+CMB+SN theta_0.3 | 1.3e-10 | 9.6e-11 | 8.9e-10 | 8.7e-10 |
| BAO+CMB+SN theta_1 | 1.1e-10 | 8.4e-11 | 1.6e-09 | 1.4e-09 |
| BAO+CMB+SN theta_3 | 1.5e-10 | 1.2e-10 | 4.3e-09 | 3.6e-09 |
| BAO+CMB+SN theta_10 | 3.1e-11 | 2.0e-11 | 5.7e-09 | 4.6e-09 |

**Paso 2, tabla antiguo / nuevo** (cambio relativo en cotas; para los máximos globales, nodo y ln Z − ln Z(0)):

| fuente | cifra | antiguo | nuevo | cambio |
|---|---|---|---|---|
| exacto_sn.linear_bounds | SN theta_0.1 192 | 0.772432 | 0.772469 | +4.76e-05 |
| exacto_sn.linear_bounds | SN theta_0.1 135 | 0.539704 | 0.539782 | +1.46e-04 |
| exacto_sn.linear_bounds | SN theta_0.3 192 | 0.446683 | 0.446705 | +4.87e-05 |
| exacto_sn.linear_bounds | SN theta_0.3 135 | 0.312096 | 0.312142 | +1.47e-04 |
| exacto_sn.linear_bounds | SN theta_1 192 | 0.246444 | 0.246457 | +5.17e-05 |
| exacto_sn.linear_bounds | SN theta_1 135 | 0.172245 | 0.172271 | +1.51e-04 |
| exacto_sn.linear_bounds | SN theta_3 192 | 0.146536 | 0.146544 | +5.61e-05 |
| exacto_sn.linear_bounds | SN theta_3 135 | 0.102775 | 0.102791 | +1.59e-04 |
| exacto_sn.linear_bounds | SN theta_10 192 | 0.092046 | 0.092049 | +2.65e-05 |
| exacto_sn.linear_bounds | SN theta_10 135 | 0.065846 | 0.065856 | +1.48e-04 |
| exacto_sn.linear_bounds | BAO+SN theta_0.1 192 | 0.832690 | 0.832705 | +1.76e-05 |
| exacto_sn.linear_bounds | BAO+SN theta_0.1 135 | 0.724570 | 0.724600 | +4.23e-05 |
| exacto_sn.linear_bounds | BAO+SN theta_0.3 192 | 0.447218 | 0.447231 | +2.95e-05 |
| exacto_sn.linear_bounds | BAO+SN theta_0.3 135 | 0.384041 | 0.384062 | +5.55e-05 |
| exacto_sn.linear_bounds | BAO+SN theta_1 192 | 0.232668 | 0.232673 | +2.53e-05 |
| exacto_sn.linear_bounds | BAO+SN theta_1 135 | 0.197165 | 0.197175 | +5.05e-05 |
| exacto_sn.linear_bounds | BAO+SN theta_3 192 | 0.144732 | 0.144732 | +1.16e-06 |
| exacto_sn.linear_bounds | BAO+SN theta_3 135 | 0.123323 | 0.123327 | +2.66e-05 |
| exacto_sn.linear_bounds | BAO+SN theta_10 192 | 0.104270 | 0.104262 | -7.19e-05 |
| exacto_sn.linear_bounds | BAO+SN theta_10 135 | 0.089561 | 0.089557 | -3.95e-05 |
| exacto_sn.linear_bounds | BAO+CMB+SN theta_0.1 192 | None | None | — |
| exacto_sn.linear_bounds | BAO+CMB+SN theta_0.1 135 | 0.917273 | 0.917275 | +2.71e-06 |
| exacto_sn.linear_bounds | BAO+CMB+SN theta_0.3 192 | 0.566548 | 0.566541 | -1.29e-05 |
| exacto_sn.linear_bounds | BAO+CMB+SN theta_0.3 135 | 0.504526 | 0.504530 | +8.85e-06 |
| exacto_sn.linear_bounds | BAO+CMB+SN theta_1 192 | 0.298185 | 0.298179 | -1.92e-05 |
| exacto_sn.linear_bounds | BAO+CMB+SN theta_1 135 | 0.262616 | 0.262616 | +3.36e-06 |
| exacto_sn.linear_bounds | BAO+CMB+SN theta_3 192 | 0.177840 | 0.177831 | -5.21e-05 |
| exacto_sn.linear_bounds | BAO+CMB+SN theta_3 135 | 0.154421 | 0.154417 | -2.84e-05 |
| exacto_sn.linear_bounds | BAO+CMB+SN theta_10 192 | 0.118970 | 0.118954 | -1.41e-04 |
| exacto_sn.linear_bounds | BAO+CMB+SN theta_10 135 | 0.100913 | 0.100902 | -1.05e-04 |
| maximos_globales exacto_sn | SN theta_0.1 192_global | 0.772432 | 0.772469 | +4.75e-05 |
| maximos_globales exacto_sn | SN theta_0.1 135_global | 0.539704 | 0.539783 | +1.45e-04 |
| maximos_globales exacto_sn | SN theta_0.1 192_clipped | 0.772432 | 0.772469 | +4.75e-05 |
| maximos_globales exacto_sn | SN theta_0.1 máximo global: nodo σ / ln Z − ln Z(0) | σ = 1e-06, +0.000000 | σ = 1e-06, +0.000000 | Δln Z +0.00e+00, mismo nodo: sí |
| maximos_globales exacto_sn | SN theta_0.3 192_global | 0.446683 | 0.446705 | +4.85e-05 |
| maximos_globales exacto_sn | SN theta_0.3 135_global | 0.312096 | 0.312142 | +1.47e-04 |
| maximos_globales exacto_sn | SN theta_0.3 192_clipped | 0.446683 | 0.446705 | +4.85e-05 |
| maximos_globales exacto_sn | SN theta_0.3 máximo global: nodo σ / ln Z − ln Z(0) | σ = 1e-06, +0.000000 | σ = 1e-06, +0.000000 | Δln Z +0.00e+00, mismo nodo: sí |
| maximos_globales exacto_sn | SN theta_1 192_global | 0.246444 | 0.246457 | +5.16e-05 |
| maximos_globales exacto_sn | SN theta_1 135_global | 0.172245 | 0.172271 | +1.51e-04 |
| maximos_globales exacto_sn | SN theta_1 192_clipped | 0.246444 | 0.246457 | +5.16e-05 |
| maximos_globales exacto_sn | SN theta_1 máximo global: nodo σ / ln Z − ln Z(0) | σ = 1e-06, +0.000000 | σ = 1e-06, +0.000000 | Δln Z +0.00e+00, mismo nodo: sí |
| maximos_globales exacto_sn | SN theta_3 192_global | 0.146536 | 0.146544 | +5.60e-05 |
| maximos_globales exacto_sn | SN theta_3 135_global | 0.102775 | 0.102791 | +1.59e-04 |
| maximos_globales exacto_sn | SN theta_3 192_clipped | 0.146536 | 0.146544 | +5.60e-05 |
| maximos_globales exacto_sn | SN theta_3 máximo global: nodo σ / ln Z − ln Z(0) | σ = 1e-06, +0.000000 | σ = 1e-06, +0.000000 | Δln Z +0.00e+00, mismo nodo: sí |
| maximos_globales exacto_sn | SN theta_10 192_global | 0.092046 | 0.092049 | +2.64e-05 |
| maximos_globales exacto_sn | SN theta_10 135_global | 0.065846 | 0.065856 | +1.48e-04 |
| maximos_globales exacto_sn | SN theta_10 192_clipped | 0.092046 | 0.092049 | +2.64e-05 |
| maximos_globales exacto_sn | SN theta_10 máximo global: nodo σ / ln Z − ln Z(0) | σ = 1e-06, +0.000000 | σ = 1e-06, +0.000000 | Δln Z +0.00e+00, mismo nodo: sí |
| maximos_globales exacto_sn | BAO+SN theta_0.1 192_global | 0.703044 | 0.703062 | +2.52e-05 |
| maximos_globales exacto_sn | BAO+SN theta_0.1 135_global | 0.591214 | 0.591242 | +4.78e-05 |
| maximos_globales exacto_sn | BAO+SN theta_0.1 192_clipped | 0.832690 | 0.832705 | +1.76e-05 |
| maximos_globales exacto_sn | BAO+SN theta_0.1 máximo global: nodo σ / ln Z − ln Z(0) | σ = 0.199272, +1.014742 | σ = 0.199272, +1.014805 | Δln Z +6.26e-05, mismo nodo: sí |
| maximos_globales exacto_sn | BAO+SN theta_0.3 192_global | 0.408588 | 0.408598 | +2.34e-05 |
| maximos_globales exacto_sn | BAO+SN theta_0.3 135_global | 0.343956 | 0.343972 | +4.60e-05 |
| maximos_globales exacto_sn | BAO+SN theta_0.3 192_clipped | 0.447218 | 0.447231 | +2.95e-05 |
| maximos_globales exacto_sn | BAO+SN theta_0.3 máximo global: nodo σ / ln Z − ln Z(0) | σ = 0.116393, +1.013019 | σ = 0.116393, +1.013086 | Δln Z +6.69e-05, mismo nodo: sí |
| maximos_globales exacto_sn | BAO+SN theta_1 192_global | 0.229720 | 0.229723 | +1.59e-05 |
| maximos_globales exacto_sn | BAO+SN theta_1 135_global | 0.194078 | 0.194085 | +3.86e-05 |
| maximos_globales exacto_sn | BAO+SN theta_1 192_clipped | 0.232668 | 0.232673 | +2.53e-05 |
| maximos_globales exacto_sn | BAO+SN theta_1 máximo global: nodo σ / ln Z − ln Z(0) | σ = 0.0679844, +1.002552 | σ = 0.0679844, +1.002639 | Δln Z +8.64e-05, mismo nodo: sí |
| maximos_globales exacto_sn | BAO+SN theta_3 192_global | 0.144732 | 0.144732 | +1.16e-06 |
| maximos_globales exacto_sn | BAO+SN theta_3 135_global | 0.123323 | 0.123327 | +2.66e-05 |
| maximos_globales exacto_sn | BAO+SN theta_3 192_clipped | 0.144732 | 0.144732 | +1.16e-06 |
| maximos_globales exacto_sn | BAO+SN theta_3 máximo global: nodo σ / ln Z − ln Z(0) | σ = 0.0397184, +0.942693 | σ = 0.0397184, +0.942809 | Δln Z +1.16e-04, mismo nodo: sí |
| maximos_globales exacto_sn | BAO+SN theta_10 192_global | 0.104270 | 0.104262 | -7.19e-05 |
| maximos_globales exacto_sn | BAO+SN theta_10 135_global | 0.089561 | 0.089557 | -3.95e-05 |
| maximos_globales exacto_sn | BAO+SN theta_10 192_clipped | 0.104270 | 0.104262 | -7.19e-05 |
| maximos_globales exacto_sn | BAO+SN theta_10 máximo global: nodo σ / ln Z − ln Z(0) | σ = 0.0315511, +0.691489 | σ = 0.0315511, +0.691755 | Δln Z +2.66e-04, mismo nodo: sí |
| maximos_globales exacto_sn | BAO+CMB+SN theta_0.1 192_global | 0.864302 | 0.864279 | -2.65e-05 |
| maximos_globales exacto_sn | BAO+CMB+SN theta_0.1 135_global | 0.754154 | 0.754148 | -8.40e-06 |
| maximos_globales exacto_sn | BAO+CMB+SN theta_0.1 192_clipped | None | None | — |
| maximos_globales exacto_sn | BAO+CMB+SN theta_0.1 máximo global: nodo σ / ln Z − ln Z(0) | σ = 0.341164, +1.113129 | σ = 0.341164, +1.113298 | Δln Z +1.69e-04, mismo nodo: sí |
| maximos_globales exacto_sn | BAO+CMB+SN theta_0.3 192_global | 0.502171 | 0.502158 | -2.55e-05 |
| maximos_globales exacto_sn | BAO+CMB+SN theta_0.3 135_global | 0.438101 | 0.438098 | -6.70e-06 |
| maximos_globales exacto_sn | BAO+CMB+SN theta_0.3 192_clipped | 0.566548 | 0.566541 | -1.29e-05 |
| maximos_globales exacto_sn | BAO+CMB+SN theta_0.3 máximo global: nodo σ / ln Z − ln Z(0) | σ = 0.184538, +1.113878 | σ = 0.184538, +1.114029 | Δln Z +1.51e-04, mismo nodo: sí |
| maximos_globales exacto_sn | BAO+CMB+SN theta_1 192_global | 0.282443 | 0.282433 | -3.51e-05 |
| maximos_globales exacto_sn | BAO+CMB+SN theta_1 135_global | 0.246176 | 0.246173 | -1.60e-05 |
| maximos_globales exacto_sn | BAO+CMB+SN theta_1 192_clipped | 0.298185 | 0.298179 | -1.92e-05 |
| maximos_globales exacto_sn | BAO+CMB+SN theta_1 máximo global: nodo σ / ln Z − ln Z(0) | σ = 0.107787, +1.086738 | σ = 0.107787, +1.086892 | Δln Z +1.54e-04, mismo nodo: sí |
| maximos_globales exacto_sn | BAO+CMB+SN theta_3 192_global | 0.177139 | 0.177129 | -5.70e-05 |
| maximos_globales exacto_sn | BAO+CMB+SN theta_3 135_global | 0.153678 | 0.153673 | -3.43e-05 |
| maximos_globales exacto_sn | BAO+CMB+SN theta_3 192_clipped | 0.177840 | 0.177831 | -5.21e-05 |
| maximos_globales exacto_sn | BAO+CMB+SN theta_3 máximo global: nodo σ / ln Z − ln Z(0) | σ = 0.0583028, +0.929204 | σ = 0.0583028, +0.929355 | Δln Z +1.50e-04, mismo nodo: sí |
| maximos_globales exacto_sn | BAO+CMB+SN theta_10 192_global | 0.118970 | 0.118954 | -1.41e-04 |
| maximos_globales exacto_sn | BAO+CMB+SN theta_10 135_global | 0.100913 | 0.100902 | -1.05e-04 |
| maximos_globales exacto_sn | BAO+CMB+SN theta_10 192_clipped | 0.118970 | 0.118954 | -1.41e-04 |
| maximos_globales exacto_sn | BAO+CMB+SN theta_10 máximo global: nodo σ / ln Z − ln Z(0) | σ = 0.0315511, +0.552450 | σ = 0.0315511, +0.552675 | Δln Z +2.24e-04, mismo nodo: sí |

**Cocientes de A35 (SN solas frente a BAO solo, `maximos_globales` → cmb_fondo.BAO|Lambda_absorbs_X0, sin filas SN, sin cambio):**

| θ | SN solas antiguo → nuevo | BAO solo | cociente antiguo → nuevo | citado antiguo → nuevo |
|---|---|---|---|---|
| 1 | 0.246444 → 0.246457 | 0.456421 | ×1.8520 → ×1.8519 | 0.246 / ×1.85 → 0.246 / ×1.85 |
| 3 | 0.146536 → 0.146544 | 0.282900 | ×1.9306 → ×1.9305 | 0.147 / ×1.93 → 0.147 / ×1.93 |
| 10 | 0.092046 → 0.092049 | 0.203306 | ×2.2087 → ×2.2087 | 0.092 / ×2.21 → 0.092 / ×2.21 |

**Cifra de A26** (`inyeccion.json` → real_check.marg_linear frente a la exacta 0.287): 0.277119 → 0.277109; citada 0.2771 → 0.2771, −3 % → −3 %.

**Nota sobre el enunciado de P11.** La predicción copia de §6 los cocientes "×1.85 / ×1.98 / ×2.22" y la etiqueta "A34"; PROPUESTA_V3 cita en A35 ×1.85 / ×1.93 / ×2.21
(`maximos_globales.json`: 1.8520 / 1.9306 / 2.2087). El error era de transcripción de esta nota (§6, corregido entre corchetes) y P11 se evalúa sobre lo que cita PROPUESTA_V3.

**P11: sí.** SN solas 0.246 / 0.147 / 0.092 y cocientes ×1.85 / ×1.93 / ×2.21 sin cambio a su precisión; A26 sin cambio; los 15 máximos globales siguen en el
mismo nodo y su ln Z − ln Z(0) cambia como mucho 2.7e-04 (< 10⁻³); cambio relativo máximo de cualquier cota 1.6e-04. Ninguna cifra citada a
tres decimales cambia en estas dos fuentes.

## 8. Regla para lo que sigue

**Todo cálculo nuevo debe usar `respuesta_v3A.build_A_v3` (o, en las rutinas con malla propia, la misma definición de la fila SN).** Los scripts
antiguos (`marginal_guiada.build_A`, `muestras_sn.make_linear`, `exacto_sn.linear_bounds`, `maximos_globales.py`) se conservan sin tocar por
trazabilidad de los json históricos; sus cifras quedan marcadas por la errata 129 y su impacto está medido aquí.
