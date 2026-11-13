# PREDICTIVO_BINS: ¿cómo de anómalo es cada bin de DESI DR2 bajo ΛCDM? (13-11-2026)

## Contexto

JACKKNIFE_BINS.md (commit 1105cdf) muestra que la preferencia por σ_X > 0 sin SN vive en LRG2 (máximo lineal +1.306 → 0 al quitarlo) y que
BGS empuja en sentido contrario (quitarlo sube el máximo a +1.94). Aquí se pregunta otra cosa: con σ_X = 0 (ΛCDM con priors de Planck
comprimido), ¿cómo de anómalo es cada bin cuando se predice con todos los demás datos? Comprobación predictiva dejando un bin fuera, con
corrección por mirar en 7 bins mediante simulaciones ΛCDM.

## Predicciones (de Claude, escritas antes de ejecutar, copiadas literales)

- P15. Sin SN (BAO + θ* + ω_m), el p-valor local de LRG2 está entre 0.003 y 0.05.
- P16. Sin SN, el p-valor global (peor bin de los 7, calibrado con simulaciones ΛCDM) es > 0.05.
- P17. Sin SN, el p-valor local de BGS es > 0.05.

## Método (`predictivo_bins.py` → `predictivo_bins.json`, log `predictivo_bins_run.log`)

**Montaje.** El de `jackknife_bins.py` (importado, sin modificar): malla G0 de exacto, datos = 13 números de DESI DR2 con su covarianza oficial
13 × 13 + 100 θ* y ω_m de Planck [+ Pantheon+, 1580 SN, STAT+SYS]; combinaciones "sin SN" (parámetros h, ω_m) y "Pantheon+" (h, ω_m, M). Los
7 bins y sus filas son los de JACKKNIFE_BINS.md (BGS 1 número; los otros seis 2 cada uno; total 13). No hace falta A (σ_X = 0).

**Paso 0, ruido del bin excluido.** Para cada bin y combinación se registra max |C[filas del bin, filas usadas en el ajuste]| (otros bins BAO,
θ*, ω_m y, con SN, las 1580 SN): **0.0 exacto en los 14 casos** (`step0` del json). La covarianza completa es diagonal por bloques (13 × 13
oficial, que a su vez es diagonal por bloques entre bins; σ_θ*²; σ_ω²; C_SN), así que la fórmula de V_k no necesita término adicional.

**Paso 1, observado (reajuste completo, no lineal).** Para cada bin k: ajuste ΛCDM a todas las filas menos las del bin k con el optimizador y
priors del fiducial original (sin SN `exacto.fit`, con Pantheon+ `exacto_sn.fit_combo`: los mismos ajustes de `jackknife_bins.py`, cuyos
ajustadores se reproducen aquí bit a bit con el vector de datos como argumento, `check_fitters_vs_jackknife_bins` = 0.0) → p0_(−k).
**Definición de Fisher, la misma en los 7 casos y en el ajuste completo:** F = Jᵀ C⁻¹ J, con J el jacobiano del modelo exacto por diferencias
finitas centradas (ε = 1e-4 en h, 1e-5 en ω_m, 1e-3 en M) sobre las filas conservadas, evaluado en el mínimo p0_(−k), y C la covarianza (diagonal
por bloques) de esas filas; sin término de segundas derivadas (Gauss-Newton). Convención de χ²: χ² = rᵀ C⁻¹ r sobre las filas conservadas, con M
como parámetro explícito en la combinación con SN (sin marginalizar analíticamente; en el mínimo es equivalente). Σ_(−k) = F⁻¹. Predicción del
bin: r_k = datos_k − modelo_k(p0_(−k)); V_k = C_kk + J_k Σ_(−k) J_kᵀ; d²_k = r_kᵀ V_k⁻¹ r_k; p local = P(χ²_{n_k} ≥ d²_k) con n_k = 1 (BGS) o 2.
**Tirón marginal** de cada componente: r_i / √V_ii con su signo. No tiene en cuenta la correlación entre las dos componentes del bin
(la correlación D_M–D_H de la medida, −0.35 a −0.49, y la que induce Σ_(−k)); el valor conjunto es d². Los dos tirones de un bin no son dos
pruebas independientes.

**Paso 2, corrección por mirar en 7 bins.** N = 20000 conjuntos de datos ΛCDM por combinación: modelo en el p0 del ajuste con todos los datos
más ruido gaussiano con la covarianza completa de la combinación (13 × 13 oficial por Cholesky, σ_θ*, σ_ω y, con SN, Pantheon+ con su Cholesky
STAT+SYS; semilla 20261113). En cada uno, los 7 reajustes dejando un bin fuera son **linealizados** (Gauss-Newton alrededor del p0 completo, J en
p0): r_k = n_k − J_k Σ_(−k) J_(−k)ᵀ C_(−k)⁻¹ n_(−k), todo lineal en el ruido, con V_k de la misma fórmula; d²_k, p local y el p local mínimo de
los 7. p global = (N_extremos + 1) / (N + 1), N_extremos = simulaciones con p mínimo ≤ el p mínimo observado del paso 1; incertidumbre
Monte Carlo binomial √(p(1 − p)/N). Comprobaciones internas: Cov(r_k) = L_k C L_kᵀ coincide con V_k a 2e-16; la media de d² simulada por bin
es ≈ n_k (1.00, 1.99, 2.02, 2.00, 1.98, 1.99, 2.00 sin SN); la fracción de simulaciones con p local ≤ el observado reproduce el p local analítico
(sin SN: BGS 0.862, LRG1 0.183, LRG2 0.029, LRG3+ELG1 0.375, ELG2 0.693, QSO 0.823, Lyα 0.999); el bin con el p mínimo se reparte casi uniformemente entre los 7 en las simulaciones
(sin SN: BGS 2864, LRG1 2834, LRG2 2874, LRG3+ELG1 2810, ELG2 2874, QSO 2772, Lyα 2972).
**Validación numérica del reajuste linealizado:** en las primeras 20 simulaciones se repiten los 7 reajustes completos no lineales y se compara
d²: mayor diferencia absoluta **0.018** (sin SN, simulación 6, bin BGS) y **0.021** (Pantheon+, simulación 10, bin LRG3+ELG1), ninguna
supera 0.1. Eso respalda el método linealizado **en esta comprobación** (estas covarianzas, estos tamaños de residuo), no en general. Sobre los
datos reales, el d² linealizado de LRG2 es 7.121 frente a 7.131 del reajuste completo (sin SN).

## Tabla (paso 1 observado; tirones marginales con signo, d² y p local)

| bin | z_eff | n_k | sin SN: tirón marginal | d² | p local | Pantheon+: tirón marginal | d² | p local |
|---|---|---|---|---|---|---|---|---|
| BGS | 0.295 | 1 | D_V -0.17 | 0.03 | 0.863 | D_V -0.30 | 0.09 | 0.762 |
| LRG1 | 0.51 | 2 | D_M +1.46, D_H -1.61 | 3.38 | 0.184 | D_M +1.36, D_H -1.65 | 3.32 | 0.190 |
| LRG2 | 0.706 | 2 | D_M -1.11, D_H -1.89 | 7.13 | 0.028 | D_M -1.21, D_H -1.92 | 7.70 | 0.021 |
| LRG3+ELG1 | 0.934 | 2 | D_M -1.39, D_H +0.59 | 1.98 | 0.372 | D_M -1.52, D_H +0.56 | 2.32 | 0.314 |
| ELG2 | 1.321 | 2 | D_M -0.82, D_H +0.53 | 0.73 | 0.695 | D_M -0.88, D_H +0.53 | 0.82 | 0.664 |
| QSO | 1.484 | 2 | D_M +0.60, D_H -0.12 | 0.40 | 0.820 | D_M +0.57, D_H -0.12 | 0.36 | 0.834 |
| Lyα | 2.33 | 2 | D_H -0.02, D_M +0.03 | 0.00 | 1.000 | D_H -0.01, D_M -0.01 | 0.00 | 1.000 |

p global (peor bin de los 7, calibrado con simulaciones ΛCDM; el p mínimo observado es el de LRG2 en las dos combinaciones):

| combinación | p mínimo observado (bin) | N | N_extremos | p global | incertidumbre MC |
|---|---|---|---|---|---|
| sin SN | 0.0283 (LRG2) | 20000 | 3545 | 0.1773 | ±0.0027 |
| Pantheon+ | 0.0213 (LRG2) | 20000 | 2787 | 0.1394 | ±0.0024 |

Cuantiles del p mínimo en las simulaciones sin SN (1 %, 5 %, 50 %): 0.0015, 0.0076, 0.0950: con 7 bins, la mediana del peor p local es ≈ 0.1.

## Veredicto

- **P15: sí.** Sin SN, p local de LRG2 = 0.0283 (d² = 7.13 con 2 g.d.l.), dentro de [0.003, 0.05].
- **P16: sí.** Sin SN, p global = 0.177 ± 0.003 (N_extremos = 3545 de 20000) > 0.05. Con Pantheon+, 0.139 ± 0.002.
- **P17: sí.** Sin SN, p local de BGS = 0.863 (tirón −0.17) > 0.05. BGS no es anómalo bajo ΛCDM: en JACKKNIFE_BINS su papel era
  frenar el máximo de σ_X (quitarlo lo sube), no desviarse de la predicción.

## Paso 3, descriptivo (sin test): signos de los vecinos de LRG2

Residuos observados (datos − predicción dejando el bin fuera, unidades de D/r_d), sin SN: LRG1 D_M +0.255, D_H -0.695; LRG2 D_M -0.210,
D_H -0.635; LRG3+ELG1 D_M -0.242, D_H +0.119. Con Pantheon+: LRG1 +0.238, -0.714; LRG2 -0.227, -0.645; LRG3+ELG1 -0.263, +0.114.
LRG2 tiene residuo negativo en D_M y en D_H. **LRG1 coincide en signo con LRG2 en D_H (negativo) y no en D_M (positivo); LRG3+ELG1 coincide
en D_M (negativo) y no en D_H (positivo)**, en las dos combinaciones. Ningún vecino repite el patrón completo de LRG2.

## Cautelas

1. Un p local bajo en un bin no demuestra error de medida ni física nueva: es una desviación de la predicción de ΛCDM con los demás datos, y
   un bin de 2 números con p ≈ 0.03 aparece una vez de cada 35 bajo ΛCDM.
2. P15 y P16 se cumplen: la desviación de LRG2 es llamativa mirando ese bin (p local 0.028), pero no excepcional habiendo mirado siete
   (p global 0.18: en una de cada seis realizaciones ΛCDM el peor de los 7 bins se desvía al menos tanto).
3. El ajuste es ΛCDM con priors de Planck comprimido (100 θ* y ω_m como dos filas gaussianas): cambiar el modelo de fondo o quitar esas dos
   filas cambia la predicción de cada bin y sus p.
4. Las simulaciones suponen covarianzas correctas y gaussianas (la 13 × 13 oficial, los dos errores de Planck y la STAT+SYS de Pantheon+); no
   incluyen sistemáticos fuera de esas matrices ni no gaussianidad.

## Ficheros

`predictivo_bins.py`, `predictivo_bins.json` (paso 0, comprobación de los ajustadores, por combinación: p0 completo, 7 ajustes sin bin con Σ,
r, V, tirones, d², p; versión linealizada sobre los datos; simulaciones con N, N_extremos, p global, error MC, comprobaciones; validación con las
140 comparaciones; paso 3; veredicto), `predictivo_bins_run.log`. Nada existente se modifica.
