# INYECCION: universos falsos, ¿qué regla da cotas correctas y cuál confunde tensión con temblor?

Exploración (no toca PREPRINT ni CLAIMS). Script `inyeccion.py`; salidas `inyeccion.{json,txt}`.

Pregunta: con universos falsos donde sabemos la verdad, ¿qué regla da cotas correctas y cuál confunde tensión con temblor?

La hipótesis sobre las dos reglas se formuló después de ver MUESTRAS_SN; este test es
su primera prueba.

## Predicciones (escritas antes de calcular nada)

1. A (ΛCDM puro): la caída de 1.92 excluye σ_X = 0 en menos del 5 % de los universos;
   la regla del 5 % y la marginal dan cotas medianas dentro de un ±15 % de las reales
   con Pantheon+.
2. B (temblor OU θ = 1 metido, σ_true = 0.1, 0.2 y 0.3): la marginal −1.92 cubre el
   σ_true en el 93–99 % de los universos; la regla del 5 % lo cubre en menos del 90 %.
3. C (tensión tipo DES, sin temblor): la regla del 5 % se afloja ×1.3 o más frente a A;
   la marginal menos de ×1.15; y la caída de 1.92 excluye σ_X = 0 en más del 50 % de
   los universos (una "detección" falsa).
4. C con menos información (SN agrupadas en 22 tramos como Union3, sin tensión): las
   dos reglas se aflojan, y la marginal más de ×1.15.

## Montaje (del prompt)

- Universos falsos: vector de datos de BAO + θ* + ω_m + Pantheon+ generado desde el
  mejor ajuste ΛCDM, con ruido tomado de las covarianzas reales. 200 universos por
  escenario y la misma semilla de ruido entre escenarios.
- B: se suma la respuesta de un camino OU θ = 1 con σ_true (un camino distinto por
  universo), con la respuesta exacta del convenio B.
- C: las SN se generan con Ω_m = Ω_m(mejor ajuste) + δ, eligiendo δ para que la Δχ² de
  tensión media sea 8.7 (como DES-SN5YR); BAO y CMB sin cambio. Da δ.
- C-tramos: las mismas SN de A comprimidas en 22 tramos de z con su covarianza.
- Para cada universo, OU θ = 1: regla del 5 % y marginal −1.92 en la versión lineal
  (χ²_min de cada camino como forma cuadrática con (h, ω_m, M) marginalizados; mismos
  2000 caminos para todos los universos), y Δln L máximo.
- Comprobaciones: en los datos reales, la versión lineal reproduce la cota exacta del
  5 % (0.169) y la marginal (0.287) a menos de un 10 %; y en 5 universos de cada
  escenario, la versión exacta (como familias.py) reproduce la lineal a menos de un 10 %.

Nota añadida el 2026-11-12 tras la revisión independiente (`REVISION_INYECCION.md`), sin cambiar el texto del prompt: el 0.169 es la
cota exacta del 5 % de FAMILIAS / MUESTRAS_SN (misma celda OU θ = 1 con Pantheon+, nodos 0.1/0.15/0.2; `muestras_sn.json`, 0.169 con
remuestreo 0.164–0.173); EXACTO_SN da 0.173 para la misma celda con otra malla de σ, y con ese valor el cociente de la comprobación previa
sería 0.901: la condición de parada se cumplió al borde. Y la segunda comprobación salió con 3 de 30 universos fuera del 10 % (sección
"Comprobación exacta").

## Decisiones fijadas antes de los escenarios

- "La marginal" de las predicciones es la marginal lineal analítica (la misma `lnZ_lin` de MUESTRAS_SN y FAMILIAS: C + σ² A K Aᵀ, priors
  planos en (h, ω_m, M)), perfilada en una malla densa de 91 σ (0 y 0.02 … 1.5) con spline cúbico y caída de 1.92 desde el máximo global.
  La marginal Monte Carlo sobre los mismos 2000 caminos (log de la media de exp(σ b_i − σ² c_i/2)) se da al lado como comprobación con su
  ESS, no se usa para juzgar las predicciones.
- Regla del 5 % lineal en forma cerrada: χ²_i(σ) − χ²_ΛCDM = −2σ b_i + σ² c_i ≤ 0 ⇔ σ ≤ 2 b_i/c_i, así que la cota es el percentil 95 de
  máx(2 b_i/c_i, 0) sobre los 2000 caminos (sin malla ni descarte por H² ≤ 0; la comprobación exacta sí descarta).
- "Cotas reales con Pantheon+" de la predicción 1: las exactas 0.169 (5 %, FAMILIAS / MUESTRAS_SN; EXACTO_SN da 0.173 con otra malla) y
  0.287 (marginal, MARGINAL_GUIADA); se dan también
  las lineales con los datos reales.
- Tensión Δχ² = χ²_min(conjunto, 3 parámetros) − χ²_min(15 filas, h y ω_m) − χ²_min(SN solas, ω_m y M con h fijo), lineal; sin tensión su
  media es 1, así que δ se elige con Δχ² sin ruido = 7.7.

## Comprobación previa con los datos reales (`inyeccion_check.json`; condición de parada no saltó)

| cota | lineal (aquí) | exacta | cociente |
|---|---|---|---|
| 5 % | 0.156 | 0.169 (FAMILIAS / MUESTRAS_SN; 0.173 en EXACTO_SN, cociente 0.901) | 0.92 |
| marginal −1.92 | 0.277 (spline 0.2771, brentq 0.2771; MUESTRAS_SN 0.277) | 0.287 | 0.97 |

Δln L máximo lineal +1.14 en σ_X = 0.104 (MUESTRAS_SN: +1.14 en 0.108); χ²_ΛCDM lineal = exacto = 1406.33; tensión lineal 2.66 (no lineal,
`muestras_sn_lcdm.json`: 2.68). La marginal Monte Carlo con 2000 caminos da 0.261 con ESS 3.4 en la cota: inservible por sí sola, como ya
se vio en MARGINAL_GUIADA.

## Montaje tal como se corrió (`inyeccion.py`, 221 s: 14 s lineal, 207 s la comprobación exacta con 16 procesos)

- Vector fiducial: B0 = modelo ΛCDM en el mejor ajuste P0 de los datos reales (h = 0.682, ω_m = 0.1407, Ω_m = 0.3037; `make_linear` de
  `muestras_sn.py`), 1595 filas; ruido n_u = L z_u con L el Cholesky de la covarianza completa (BAO 13 × 13, θ*, ω_m, Pantheon+
  1580 × 1580), z_u fijo por universo (semilla 20261118) en todos los escenarios.
- B: camino OU θ = 1 unitario por universo (otra semilla), escalado por σ_true y restado X(0); su respuesta por el modelo exacto
  (`exact_model`) menos B0; ninguno dio H² ≤ 0.
- C: δ = 0.0512 (Ω_m de las SN 0.3037 → 0.3549, cerca del 0.352 de DES-SN5YR); Δχ² sin ruido 7.70, media prevista 8.70, media obtenida
  8.09 (mediana 7.31; la media tiene error ≈ 0.4 con 200 universos).
- D: 22 tramos de z con bordes en las medias geométricas de los nodos de Union3 (primero 0.01, último 2.26); todos con SN (cuentas 579,
  69, 96, 147, 136, 117, 104, 57, 51, 38, 35, 39, 22, 25, 29, 8, 2, 5, 2, 5, 11, 3); compresión GLS T = (EᵀWE)⁻¹EᵀW de datos, B0, Bm y A,
  covarianza (EᵀWE)⁻¹.
- Comprobación exacta: 5 universos por escenario, 6 σ entre 0.6 y 1.6 veces la cota lineal, 2000 ajustes no lineales por σ (caminos
  σ(x − x(0)), descarte por H² ≤ 0 en el mejor ajuste ΛCDM del universo), P interpolada en ln σ hasta el 5 %.

## Tabla (200 universos por escenario; mediana y percentiles 16–84; cobertura = σ_true ≤ cota; "σ_X = 0 excluido" = Δln L máximo > 1.92)

| escenario | σ_true | 5 % | marginal −1.92 | MC 2000 caminos | cobertura 5 % / marginal | σ_X = 0 excluido | Δln L máx. mediana | Δχ² tensión media | 5 % /A | marginal /A |
|---|---|---|---|---|---|---|---|---|---|---|
| A (ΛCDM) | — | 0.095 (0.069–0.128) | 0.109 (0.082–0.165) | 0.111 | — | 0.000 | +0.00 | 0.93 | 1 | 1 |
| B | 0.1 | 0.195 (0.130–0.334) | 0.222 (0.161–0.320) | 0.215 | 0.925 / 0.955 | 0.585 | +2.85 | 11.4 | 2.06 | 2.03 |
| B | 0.2 | 0.348 (0.220–0.651) | 0.376 (0.282–0.513) | 0.292 | 0.865 / 0.965 | 0.945 | +17.7 | 45.8 | 3.67 | 3.43 |
| B | 0.3 | 0.507 (0.305–0.986) | 0.522 (0.404–0.693) | 0.347 | 0.850 / 0.975 | 0.995 | +46.7 | 110 | 5.34 | 4.77 |
| C (tensión) | — | 0.162 (0.115–0.204) | 0.165 (0.124–0.222) | 0.164 | — | 0.430 | +1.46 | 8.09 | 1.71 | 1.50 |
| D (22 tramos) | — | 0.098 (0.070–0.136) | 0.111 (0.084–0.166) | 0.112 | — | 0.000 | +0.00 | 0.95 | 1.04 | 1.02 |

Cocientes emparejados por universo (mismo ruido): C/A 1.80 (1.12–2.35) en el 5 % y 1.58 (1.07–1.94) en la marginal; D/A 1.05 (0.96–1.12) y
1.04 (0.95–1.11). Ninguna marginal quedó sin cruce; cotas del 5 % por encima de 1 (donde la exacta diría "sin cota"): 3.5 % de los
universos en B 0.2 y 15.5 % en B 0.3. La marginal Monte Carlo con 2000 caminos coincide con la analítica en A, C y D (ESS alta) y se
queda corta en B 0.2–0.3 (0.29 y 0.35 frente a 0.38 y 0.52: colapso de pesos, como en MARGINAL_GUIADA).

Los datos reales frente al escenario A: la cota lineal del 5 % real (0.156) está en el percentil 97 de la distribución de A (la exacta
0.169, en el 99.5; máximo de A 0.178), la marginal real (0.277) en el 97.5 (percentil 95 de A: 0.229), el Δln L máximo real (+1.14) en el
97 (percentil 95 de A: +0.82) y la tensión real (2.66) en el 91.5. Pantheon+ con BAO + θ* + ω_m no es una realización típica de ΛCDM puro
para estos estadísticos; qué la hace atípica (tensión, temblor u otra cosa) no se decide aquí.

## Comprobación exacta (`exact_check` del json)

La comprobación prevista en el montaje ("la versión exacta reproduce la lineal a menos de un 10 %") **no se cumple en 3 de los 30
universos**: 30 de 30 con cruce del 5 %, cociente exacta/lineal 0.885–1.102, 27 de 30 dentro del 10 %. Fuera: B 0.3 u0 (0.885) y u2
(0.898), con cotas 1.12 y 1.36 (> 1, donde la lineal ya no vale y la exacta descarta el 45–54 % de los caminos por H² ≤ 0), y C u1
(1.102, cota 0.31 frente a 0.28). Con cotas < 1 la lineal reproduce la exacta a ±8 % salvo ese caso al 10.2 %. La condición de parada
del prompt era solo la de los datos reales (cumplida, al borde: cociente 0.92 con el 0.169 de FAMILIAS / MUESTRAS_SN y 0.901 con el
0.173 de EXACTO_SN); esta comprobación no detuvo nada y se deja constancia.

## Predicción, punto por punto

1. **A: mitad.** La caída de 1.92 excluye σ_X = 0 en el 0 % de los universos (< 5 %, cumplido). Las cotas medianas no están dentro del
   ±15 % de las reales: 0.095 frente a 0.169 (−44 %) en el 5 % y 0.109 frente a 0.287 (−62 %) en la marginal (frente a las lineales
   reales 0.156 y 0.277: −39 % y −60 %). Los universos ΛCDM dan cotas mucho más apretadas que los datos reales (véase arriba).
2. **B: cumplida para la marginal, a medias para el 5 %.** La marginal −1.92 cubre σ_true en el 95.5 %, 96.5 % y 97.5 % (dentro de
   93–99 % las tres). El 5 % lo cubre en el 92.5 % (σ_true = 0.1, no es < 90 %), 86.5 % y 85.0 % (sí).
3. **C: un tercio.** El 5 % se afloja ×1.71 frente a A (≥ 1.3, cumplido). La marginal se afloja ×1.50 (predicho < 1.15: no). La caída de
   1.92 excluye σ_X = 0 en el 43 % de los universos (predicho > 50 %: no, por poco).
4. **D: no cumplida.** Ni el 5 % (×1.04) ni la marginal (×1.02) se aflojan. Causa: la compresión GLS en 22 tramos no quita información
   para un modelo suave: el error Fisher de Ω_m con SN solas es 0.01745 con las 1580 SN y 0.01754 con los 22 tramos (publicado 0.018), y
   el conjunto da σ(ω_m) = 0.000657 en los dos casos. Lo que hace a Union3 menos informativa es la muestra (su error publicado de Ω_m,
   0.027), no los tramos; el escenario no reproduce eso. Un escenario con menos información (menos SN o errores mayores) sería otro test.

## Respuesta

Con la verdad conocida, la marginal −1.92 da cotas correctas como cota superior: su cobertura nominal unilateral es ≈ 97.5 % (límite
gaussiano) y la empírica es 95.5–97.5 % del σ_true. La regla del 5 % no es un intervalo de confianza y no tiene cobertura nominal; su
cobertura empírica es 85–92.5 % (92.5 % con σ_true = 0.1, 86.5 % y 85.0 % con 0.2 y 0.3), así que como cota superior deja fuera el valor
verdadero más a menudo que la marginal cuando hay temblor; es una comparación empírica, no la violación de un objetivo. Pero ninguna de las dos distingue tensión de temblor: con una tensión tipo
DES-SN5YR y sin temblor, las dos cotas se aflojan (5 % ×1.71, marginal ×1.50) y la caída de 1.92 "detecta" σ_X > 0 en el 43 % de los
universos. La hipótesis de MUESTRAS_SN ("la marginal sigue la información y el 5 % mezcla información y tensión") no se sostiene en su
primera prueba: la marginal también se mueve con la tensión, solo algo menos que el 5 %. Y el escenario de "menos información" no se
construyó bien (los tramos no quitan información), así que la parte de información queda sin probar.

## Nota sobre "Δln L máximo" (comprobado en el código, nada recalculado)

Es el mismo estadístico que la curva de la marginal: en `inyeccion.py`, `prof = self.lnz_profile(R)` (línea 175) es lnZ_lin(σ) − lnZ_lin(0)
en la malla GRID (σ = 0 es el primer nodo; la resta está en la línea 161), `"dlnL_max": prof.max(axis=0)` (línea 177) es su máximo,
`"excl0": prof.max(axis=0) > DROP` (línea 178, DROP = 1.92), y la cota −1.92 sale de `drop_bound(prof[:, j])` (línea 181) sobre ese mismo
`prof`, con umbral `thr = prof[k0] − DROP` desde el máximo global (línea 207). Por tanto "Δln L máx > 1.92" equivale exactamente a que el
intervalo de la caída de 1.92 desde el máximo global de la marginal lineal excluya σ_X = 0 (máximo buscado en los nodos de la malla, no en
el spline; la malla llega a σ = 1.5). La versión Monte Carlo (`mc_dlnL_max`, `mc_excl0`) hace lo mismo sobre su propia curva.

## Addendum (13-11-2026): fila SN de la respuesta lineal (errata 129)

Todo lo de esta nota usa la A de `muestras_sn.make_linear` antigua (errata 129, confirmada). Recalculado con la fila corregida y las mismas semillas
(`inyeccion_v3A.json`): las coberturas de la marginal (95.5 / 96.5 / 97.5 %) y de la regla del 5 % no cambian en ningún universo, el 43 % de falsas
exclusiones del escenario C se mantiene, las medianas de las cotas cambian ≤ 6.5 × 10⁻⁴ relativo; solo cambia, en la tercera decimal, "marginal /A"
de B 0.2 (3.430 → 3.431) y B 0.3 (4.770 → 4.771), y la comprobación con datos reales pasa de 0.27712 a 0.27709 (sigue siendo 0.2771 y −3 %). Las cifras
originales no se reescriben; detalle en `CORRECCION_FILA_SN.md` (§4.4). Todo cálculo nuevo debe usar `respuesta_v3A`.
