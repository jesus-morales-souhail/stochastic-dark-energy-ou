# ¿Cuánto bajan los techos de EXACTO.md al añadir las supernovas de Pantheon+? (exploración)

**Predicción, escrita antes de correr `exacto_sn.py`:**
- BAO + θ* + ω_m + SN, cota exacta al 5 %: entre 0.06 y 0.20 en θ = 1; entre 0.05 y
  0.15 en θ = 3; entre 0.04 y 0.13 en θ = 10.
- θ = 0.3: aparece cota, entre 0.2 y 0.6. θ = 0.1: sin cota por debajo de 0.5.
- SN solas acotan más que BAO solo (0.47 / 0.28 / 0.20) en θ = 1, 3 y 10.
- En los techos nuevos, menos del 10 % de los supervivientes tiene ρ_DE ≤ 0 en algún punto.
- Con el corte z_HD > 0.023 la cota se afloja menos de un 30 % en θ = 1 y más de un
  30 % en θ = 10.
- La lineal marginalizada (−1.92) queda a menos de un factor 1.5 de la exacta.

Exploración: no cambia `results/`, `figures/` ni ningún paper. Script `exacto_sn.py`; salidas `exacto_sn.{json,txt}` y
`exacto_sn_chi2.npz` (χ²_min de cada realización en las cuatro combinaciones, marca de descarte y de ρ_DE ≤ 0). Tiempo:
10 min con 16 procesos (las celdas, 4 min; el resto, los techos y las cotas lineales).

**Nota de malla (añadida el 2026-11-12, sin recalcular nada aquí).** La cota del 5 % se obtiene interpolando P(χ²_min ≤ χ²_obs)
linealmente en ln σ_X entre las celdas, y en BAO + θ* + ω_m + SN las celdas que encierran la cota eran 0.1 / 0.2 (θ = 1 y 3) y
0.05 / 0.1 (θ = 10), con P cóncava hacia arriba entre ellas: la recta corta tarde. Con las celdas 0.07 y 0.15 en la malla
(`CONVERGENCIA_SN.md`, mismo montaje, N = 20000) las cotas son **0.170 / 0.102 / 0.073** en θ = 1 / 3 / 10 (−2 %, −10 %, −10 %), y
las celdas OU de `FAMILIAS.md` (caminos idénticos, N = 2000, con 0.07 y 0.15) dan 0.169 y 0.070 en θ = 1 y 10. Las cifras
históricas de esta nota (0.173 / 0.113 / 0.081) se dejan como se calcularon; cada cita lleva la marca [con celdas 0.07 y 0.15: …].

## Montaje

El de `exacto.py` sin cambios: convenio B (ρ_DE = Ω_Λ + X(z) − X(0)), misma malla en u = ln(1+z) con los z de DESI como nodos
exactos (3506 nodos), mismas semillas por celda, mismo descarte por H² ≤ 0 en el mejor ajuste ΛCDM de las 15 filas, mismo
criterio P(χ²_min ≤ χ²_obs) = 5 % interpolando en ln σ_X. Las SN entraron sin tocar nada de eso; lo añadido:
- **Datos SN** como `sn/perfil_sn.py`: Pantheon+ (`sn/datos/Pantheon+SH0ES.dat`, covarianza STAT+SYS), m_b_corr con
  z_HD > 0.01 y sin calibradoras: 1580 SN; con z_HD > 0.023, 1371. M marginalizada analíticamente con prior plano,
  χ²_SN = rᵀC⁻¹r − (1ᵀC⁻¹r)²/(1ᵀC⁻¹1) (las constantes independientes del modelo se omiten: se cancelan en P).
- **Modelo SN:** m = 5 log₁₀[(1 + z_HEL) D_M(z_HD)] + M, con D_M de la misma historia exacta (la H² de `exacto.py`),
  interpolada linealmente en u en la malla (error < 10⁻⁶ relativo).
- **Celdas:** las de `exacto.py` (σ_X ∈ {0.05, 0.1, 0.2, 0.3, 0.5, 0.7}) más σ_X = 0.01, 0.02 y 0.03 (semillas nuevas
  [SEED, 100 + índice]); θ ∈ {0.1, 0.3, 1, 3, 10}; 2000 realizaciones por celda.
- **Combinaciones**, cada una con su χ²_obs (ajuste con X = 0 a esos mismos datos): SN solas (h degenerado con M: se ajusta
  solo Ω_m, h = 0.68 fijo para Ω_r), χ²_obs = 1386.40; BAO + SN (h, ω_m), 1399.83
  (H0 = 68.76, ω_m = 0.14377); BAO + θ* + ω_m + SN (h, ω_m),
  1406.33 (H0 = 68.12, ω_m = 0.14092); variante
  con z_HD > 0.023, 1228.77. En los ajustes por mínimos cuadrados M es un tercer parámetro, lo
  que en el mínimo equivale a la marginalización analítica.
- Un detalle de ejecución, no de montaje: con 16 procesos y productos 1580 × 1580, cada proceso usa un solo hilo de BLAS
  (sin eso la máquina se satura).

**Comprobaciones por dos caminos.**
1. Con X = 0, SN solas dan Ω_m = **0.3324** (`SN.md`: 0.3326; la diferencia, 2 × 10⁻⁴, es la radiación y la
   malla).
2. Malla: las celdas de θ = 10 repetidas con paso 10⁻³ en u por debajo de z = 0.3 (3637 nodos; mismas semillas, nodos
   distintos, así que los caminos no son los mismos). Máximo de |ΔP| sobre el error de Monte Carlo (√2 × binomial) en las
   nueve celdas: 2.1 (SN), 1.5 (BAO + SN), 1.3 (BAO + θ* + ω_m + SN), 0.9 (z > 0.023):
   dentro del ruido; las cotas con malla fina van entre corchetes en la tabla y difieren en ≤ 3 %.
3. Cota lineal marginalizada (−1.92 y −1.35) con el núcleo de respuesta completa del convenio B, como `cmb_fondo.py`, con las
   filas de SN añadidas (δm_i = (5/ln 10) δlnD_M(z_HD,i), interpoladas entre nodos; M como plantilla lineal de perturbación
   con prior plano), para las tres combinaciones y los cinco θ: tabla más abajo.

## Resultados

**Cota exacta al 5 %** (entre corchetes, la de θ = 10 con la malla fina):

| θ | SN solas | BAO + SN | BAO + θ* + ω_m + SN | BAO + θ* + ω_m + SN, z_HD > 0.023 | EXACTO.md (sin SN) |
|---|---|---|---|---|---|
| 0.1 | 0.641 | 0.598 | 0.476 | 0.456 | sin cota (≤ 0.7) |
| 0.3 | 0.376 | 0.353 | 0.307 | 0.284 | sin cota (≤ 0.7) |
| 1 | 0.199 | 0.192 | 0.173 | 0.165 | 0.427 |
| 3 | 0.131 | 0.139 | 0.113 | 0.100 | 0.286 |
| 10 | 0.082 [0.082] | 0.089 [0.086] | 0.081 [0.080] | 0.075 [0.076] | 0.192 |

Nota de malla: la columna BAO + θ* + ω_m + SN en θ = 1 / 3 / 10 (0.173 / 0.113 / 0.081) es 0.170 / 0.102 / 0.073 con las celdas
0.07 y 0.15 (`CONVERGENCIA_SN.md`); las demás columnas no se han repetido con esa malla.

- **Con BAO + θ* + ω_m + SN los techos de `EXACTO.md` bajan a menos de la mitad:** 0.43 → 0.17 (θ = 1), 0.29 → 0.11 (θ = 3),
  0.19 → 0.081 (θ = 10) [con celdas 0.07 y 0.15: 0.170 / 0.102 / 0.073, CONVERGENCIA_SN.md; ver nota de malla]. Y aparecen cotas donde no las había: 0.31 en θ = 0.3 y 0.48 en θ = 0.1.
- **SN solas acotan más que BAO solo** (0.47 / 0.28 / 0.20 en θ = 1 / 3 / 10): 0.20 / 0.13 / 0.082, un factor 2.2–2.4. Son
  1580 distancias frente a 13, y su propia M absorbe la parte constante igual que (h, ω_m) la absorbían en BAO.
- BAO + SN casi no mejora a SN solas (0.19 / 0.14 / 0.089), y el CMB comprimido añade un 10–20 % más (0.17 / 0.11 / 0.081 [con celdas 0.07 y 0.15: 0.170 / 0.102 / 0.073, CONVERGENCIA_SN.md; ver nota de malla]).
- **El corte z_HD > 0.023 no afloja la cota: la estrecha un poco** (0.173 → 0.165 en θ = 1, −5 %; 0.081 → 0.075 en θ = 10,
  −7 %; 0.113 → 0.100 en θ = 3, −12 %; las de partida son 0.170 / 0.102 / 0.073 con celdas 0.07 y 0.15, `CONVERGENCIA_SN.md`, y el
  corte z_HD no se ha repetido con esa malla). Se pierden 209 SN a z < 0.023, pero cambia también χ²_obs (1406 → 1229),
  y la cota es la fracción de realizaciones que igualan a ΛCDM, no la varianza: las SN más cercanas, donde el residuo
  X(z) − X(0) es más pequeño, aportaban poco y algo de holgura.

**P(χ²_min ≤ χ²_obs) por celda** (filas θ; la última fila de cada tabla es θ = 10 con la malla fina):

SN solas:

| θ \ σ_X | 0.01 | 0.02 | 0.03 | 0.05 | 0.1 | 0.2 | 0.3 | 0.5 | 0.7 |
|---|---|---|---|---|---|---|---|---|---|
| 0.1 | 0.474 | 0.490 | 0.463 | 0.426 | 0.354 | 0.256 | 0.159 | 0.066 | 0.044 |
| 0.3 | 0.493 | 0.444 | 0.421 | 0.382 | 0.268 | 0.139 | 0.074 | 0.019 | 0.013 |
| 1 | 0.465 | 0.403 | 0.373 | 0.289 | 0.152 | 0.050 | 0.021 | 0.004 | 0.002 |
| 3 | 0.423 | 0.357 | 0.277 | 0.186 | 0.073 | 0.014 | 0.004 | 0.001 | 0.000 |
| 10 | 0.428 | 0.286 | 0.215 | 0.113 | 0.025 | 0.003 | 0.002 | 0.000 | 0.000 |
| 10, malla fina | 0.413 | 0.316 | 0.225 | 0.117 | 0.022 | 0.003 | 0.000 | 0.000 | 0.000 |

BAO + SN:

| θ \ σ_X | 0.01 | 0.02 | 0.03 | 0.05 | 0.1 | 0.2 | 0.3 | 0.5 | 0.7 |
|---|---|---|---|---|---|---|---|---|---|
| 0.1 | 0.471 | 0.483 | 0.485 | 0.460 | 0.427 | 0.328 | 0.189 | 0.068 | 0.035 |
| 0.3 | 0.502 | 0.473 | 0.473 | 0.448 | 0.355 | 0.154 | 0.069 | 0.011 | 0.004 |
| 1 | 0.498 | 0.454 | 0.440 | 0.374 | 0.202 | 0.041 | 0.007 | 0.001 | 0.001 |
| 3 | 0.476 | 0.410 | 0.358 | 0.226 | 0.086 | 0.010 | 0.001 | 0.000 | 0.000 |
| 10 | 0.460 | 0.334 | 0.270 | 0.131 | 0.034 | 0.000 | 0.000 | 0.000 | 0.000 |
| 10, malla fina | 0.437 | 0.351 | 0.279 | 0.136 | 0.027 | 0.000 | 0.000 | 0.000 | 0.000 |

BAO + θ* + ω_m + SN:

| θ \ σ_X | 0.01 | 0.02 | 0.03 | 0.05 | 0.1 | 0.2 | 0.3 | 0.5 | 0.7 |
|---|---|---|---|---|---|---|---|---|---|
| 0.1 | 0.470 | 0.462 | 0.474 | 0.439 | 0.367 | 0.231 | 0.131 | 0.042 | 0.025 |
| 0.3 | 0.485 | 0.470 | 0.450 | 0.403 | 0.270 | 0.112 | 0.052 | 0.011 | 0.004 |
| 1 | 0.482 | 0.438 | 0.418 | 0.303 | 0.144 | 0.025 | 0.007 | 0.001 | 0.000 |
| 3 | 0.474 | 0.401 | 0.325 | 0.182 | 0.059 | 0.005 | 0.001 | 0.000 | 0.000 |
| 10 | 0.455 | 0.317 | 0.243 | 0.107 | 0.024 | 0.000 | 0.000 | 0.000 | 0.000 |
| 10, malla fina | 0.438 | 0.337 | 0.255 | 0.117 | 0.018 | 0.000 | 0.000 | 0.000 | 0.000 |

BAO + θ* + ω_m + SN, z_HD > 0.023:

| θ \ σ_X | 0.01 | 0.02 | 0.03 | 0.05 | 0.1 | 0.2 | 0.3 | 0.5 | 0.7 |
|---|---|---|---|---|---|---|---|---|---|
| 0.1 | 0.473 | 0.460 | 0.465 | 0.418 | 0.324 | 0.201 | 0.118 | 0.035 | 0.021 |
| 0.3 | 0.473 | 0.464 | 0.434 | 0.379 | 0.238 | 0.099 | 0.043 | 0.007 | 0.003 |
| 1 | 0.482 | 0.423 | 0.383 | 0.267 | 0.124 | 0.022 | 0.004 | 0.002 | 0.000 |
| 3 | 0.464 | 0.361 | 0.287 | 0.157 | 0.050 | 0.005 | 0.001 | 0.000 | 0.000 |
| 10 | 0.433 | 0.294 | 0.211 | 0.099 | 0.017 | 0.000 | 0.000 | 0.000 | 0.000 |
| 10, malla fina | 0.423 | 0.304 | 0.222 | 0.105 | 0.013 | 0.000 | 0.000 | 0.000 | 0.000 |

Descartes por H² ≤ 0 (idénticos a `EXACTO.md`, mismas semillas): ninguno con θ ≤ 0.3 ni σ_X ≤ 0.2; 0.7 % en (θ = 1, 0.7);
2.5 % y 14 % en (θ = 3, 0.5 y 0.7); 0.9 %, 18 % y 43 % en (θ = 10, 0.3, 0.5 y 0.7). Ninguna de las cotas nuevas llega a una
celda con descartes.

**ρ_DE ≤ 0 entre los supervivientes en los techos nuevos** (caminos regenerados en el techo con la semilla de la celda más
próxima, sin suavizar, sobre toda la malla hasta z*; ρ_DE = Ω_Λ + X − X(0) en el mejor ajuste ΛCDM de las 15 filas):

| θ | combinación | techo σ_X | semilla de la celda | supervivientes | ρ_DE ≤ 0 entre ellos | H² ≤ 0 |
|---|---|---|---|---|---|---|
| 0.1 | BAO + θ* + ω_m + SN | 0.476 | 0.5 | 95 (4.8 %) | 18.9 % | 0.0 % |
| 0.1 | BAO + θ* + ω_m + SN, z_HD > 0.023 | 0.456 | 0.5 | 91 (4.5 %) | 15.4 % | 0.0 % |
| 0.3 | BAO + θ* + ω_m + SN | 0.307 | 0.3 | 94 (4.7 %) | 16.0 % | 0.0 % |
| 0.3 | BAO + θ* + ω_m + SN, z_HD > 0.023 | 0.284 | 0.3 | 97 (4.9 %) | 14.4 % | 0.0 % |
| 1 | BAO + θ* + ω_m + SN | 0.173 | 0.2 | 82 (4.1 %) | 2.4 % | 0.0 % |
| 1 | BAO + θ* + ω_m + SN, z_HD > 0.023 | 0.165 | 0.2 | 81 (4.0 %) | 1.2 % | 0.0 % |
| 3 | BAO + θ* + ω_m + SN | 0.113 | 0.1 | 89 (4.5 %) | 0.0 % | 0.0 % |
| 3 | BAO + θ* + ω_m + SN, z_HD > 0.023 | 0.100 | 0.1 | 100 (5.0 %) | 0.0 % | 0.0 % |
| 10 | BAO + θ* + ω_m + SN | 0.081 | 0.1 | 98 (4.9 %) | 0.0 % | 0.0 % |
| 10 | BAO + θ* + ω_m + SN, z_HD > 0.023 | 0.075 | 0.1 | 96 (4.8 %) | 0.0 % | 0.0 % |

Nota de malla: los techos 0.173 / 0.113 / 0.081 (θ = 1 / 3 / 10) son 0.170 / 0.102 / 0.073 con las celdas 0.07 y 0.15
(`CONVERGENCIA_SN.md`); la regeneración en el techo no se ha repetido en ellos.

En los techos de θ = 1, 3 y 10 casi ningún superviviente tiene ρ_DE ≤ 0 (2.4 %, 0, 0): a σ_X ≤ 0.17 el residuo no llega a
anular Ω_Λ ≈ 0.69. En θ = 0.3 y 0.1, con techos de 0.3–0.5, lo tienen el 15–19 %. Comparado con `CMB_TECHOS.md` (67 / 28 /
4.4 % en los techos sin SN, tras suavizar y sobre la tabla de CAMB), el problema de representar estos caminos en CAMB
desaparece en los techos nuevos de θ ≥ 1.

**Cota lineal marginalizada frente a la exacta:**

| θ | SN solas: exacta / lineal −1.92 / −1.35 / cociente | BAO + SN: exacta / lineal −1.92 / −1.35 / cociente | BAO + θ* + ω_m + SN: exacta / lineal −1.92 / −1.35 / cociente |
|---|---|---|---|
| 0.1 | 0.641 / 0.772 / 0.540 / 1.21 | 0.598 / 0.833 / 0.725 / 1.39 | 0.476 / sin cota / 0.917 / — |
| 0.3 | 0.376 / 0.447 / 0.312 / 1.19 | 0.353 / 0.447 / 0.384 / 1.27 | 0.307 / 0.567 / 0.505 / 1.84 |
| 1 | 0.199 / 0.246 / 0.172 / 1.24 | 0.192 / 0.233 / 0.197 / 1.21 | 0.173 / 0.298 / 0.263 / 1.72 |
| 3 | 0.131 / 0.147 / 0.103 / 1.11 | 0.139 / 0.145 / 0.123 / 1.04 | 0.113 / 0.178 / 0.154 / 1.57 |
| 10 | 0.082 / 0.092 / 0.066 / 1.12 | 0.089 / 0.104 / 0.090 / 1.17 | 0.081 / 0.119 / 0.101 / 1.48 |

Nota de malla: la exacta de BAO + θ* + ω_m + SN (0.173 / 0.113 / 0.081) es 0.170 / 0.102 / 0.073 con las celdas 0.07 y 0.15
(`CONVERGENCIA_SN.md`); los cocientes con esa malla serían 1.75 / 1.74 / 1.63 con las mismas lineales.


**Cotas lineales con el máximo global (corrección documental, 2026-11-12).** Las cifras de la tabla anterior marcadas aquí como
sustituidas están referidas a un máximo recortado en σ_X ≤ 0.05; estas son las correctas:

| θ | SN solas: −1.92 / −1.35 | BAO + SN: −1.92 / −1.35 (sustituida → global) | BAO + θ* + ω_m + SN: −1.92 / −1.35 (sustituida → global) | exacta (5 %) / lineal global: SN, BAO + SN, BAO + CMB + SN |
|---|---|---|---|---|
| 0.1 | 0.772 / 0.540 (igual) | 0.703 / 0.591 (sustituye a 0.833 / 0.725 (-15.6 %)) | 0.864 / 0.754 (sustituye a sin cota / 0.917) | 1.21, 1.18, 1.81 |
| 0.3 | 0.447 / 0.312 (igual) | 0.409 / 0.344 (sustituye a 0.447 / 0.384 (-8.6 %)) | 0.502 / 0.438 (sustituye a 0.567 / 0.505 (-11.4 %)) | 1.19, 1.16, 1.63 |
| 1 | 0.246 / 0.172 (igual) | 0.230 / 0.194 (sustituye a 0.233 / 0.197 (-1.3 %)) | 0.282 / 0.246 (sustituye a 0.298 / 0.263 (-5.3 %)) | 1.24, 1.20, 1.63 |
| 3 | 0.147 / 0.103 (igual) | 0.145 / 0.123 (igual) | 0.177 / 0.154 (igual) | 1.11, 1.04, 1.57 |
| 10 | 0.092 / 0.066 (igual) | 0.104 / 0.090 (igual) | 0.119 / 0.101 (igual) | 1.12, 1.17, 1.48 |

**Por qué cambian:** `perfil_fondo.bound` mide la caída de 1.92 (o 1.35) desde el máximo de la verosimilitud sobre su malla de σ_X, que llega solo a 0.05; en las combinaciones con tensión el máximo está más allá de 0.05, y la cota medida desde ese máximo recortado sale más floja de lo debido. Descubierto en `MARGINAL_GUIADA.md`; recálculo en `maximos_globales.py` → `maximos_globales.json` (los valores históricos se reproducen antes con el máximo recortado). Sin CMB, o donde el máximo está en σ_X → 0 o dentro de la malla, nada cambia. Las cotas del 5 % no dependen de ningún máximo y no cambian.

**Qué cambia en la lectura y en la predicción:** el cociente lineal/exacta (5 %) pasa a ×1.04–1.24 sin CMB (antes ×1.04–1.39) y a ×1.48–1.81 con el CMB comprimido (antes ×1.48–1.84; en θ = 0.1 ahora hay cota lineal, 0.864, y el cociente es 1.82). El veredicto de la
predicción "la lineal marginalizada a menos de un factor 1.5 de la exacta" no cambia: cumplida sin CMB, no cumplida con él.

- SN solas y BAO + SN: la lineal queda a ×1.04–1.39 de la exacta (siempre más floja). [Con el máximo global: ×1.04–1.20; ver la tabla corregida.]
- BAO + θ* + ω_m + SN: ×1.48–1.84 (y sin cota lineal en θ = 0.1, exacta 0.48). [Con el máximo global: ×1.48–1.82, con cota lineal en θ = 0.1; ver la tabla corregida.] Es lo mismo que en `EXACTO.md` con el CMB
  comprimido: el perfil lineal tiene el máximo en σ_X > 0 (absorbe en promedio la tensión DESI–Planck) y la cota se mide desde
  ahí; la exacta cuenta realizaciones concretas y no tiene ese sesgo.

## Predicción: qué se cumplió

- **BAO + θ* + ω_m + SN, cota exacta entre 0.06 y 0.20 (θ = 1), 0.05 y 0.15 (θ = 3), 0.04 y 0.13 (θ = 10): cumplida** (0.173,
  0.113, 0.081; con celdas 0.07 y 0.15: 0.170 / 0.102 / 0.073, `CONVERGENCIA_SN.md`, también dentro de los rangos).
- **θ = 0.3, aparece cota entre 0.2 y 0.6: cumplida** (0.307). **θ = 0.1, sin cota por debajo de 0.5: no se cumplió, por
  poco** (0.476; P = 4.2 % en σ_X = 0.5).
- **SN solas acotan más que BAO solo en θ = 1, 3 y 10: cumplida** (0.199 / 0.131 / 0.082 frente a 0.47 / 0.28 / 0.20).
- **Menos del 10 % de supervivientes con ρ_DE ≤ 0 en los techos nuevos: cumplida en θ = 1, 3 y 10** (2.4 %, 0, 0); en los
  techos de θ = 0.3 y 0.1, que no existían, es el 16 y el 19 %.
- **Con z_HD > 0.023 la cota se afloja menos de un 30 % en θ = 1 y más de un 30 % en θ = 10: no se cumplió en el sentido.**
  No se afloja en ningún θ: se estrecha un 5 % (θ = 1), 12 % (θ = 3) y 7 % (θ = 10).
- **La lineal marginalizada (−1.92) a menos de un factor 1.5 de la exacta: cumplida para SN solas y BAO + SN** (1.04–1.39);
  **no para BAO + θ* + ω_m + SN** (1.48 en θ = 10, 1.57–1.84 en θ ≤ 3, sin cota en θ = 0.1). [Con las lineales de máximo global: sin CMB ×1.04–1.20 y con CMB ×1.48–1.82; mismo veredicto.]

**Número para la v3.0 (convenio B, cota exacta, 5 %):** DESI DR2 + ω_m y θ* de Planck + Pantheon+ dejan pasar un temblor OU
de σ_X < 0.17 (θ = 1), 0.11 (θ = 3), 0.081 (θ = 10) [con celdas 0.07 y 0.15: 0.170 / 0.102 / 0.073 (`CONVERGENCIA_SN.md`), es decir
0.17 / 0.10 / 0.073], 0.31 (θ = 0.3) y 0.48 (θ = 0.1); las supernovas son lo que más acota, por
delante de los BAO y del CMB comprimido, y casi toda la cota viene de ellas.
