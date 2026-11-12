# ¿Cuánto error tienen las cotas de EXACTO_SN.md y cuánto dependen de la regla del 5 %? (exploración)

**Predicción, escrita antes de correr `precision_sn.py`:**
- Error de Monte Carlo (bootstrap) de las cotas del 5 %: entre el 4 y el 12 %.
- Con las celdas nuevas, ninguna cota se mueve más de un 15 %.
- La cota por verosimilitud marginal exacta (−1.92) queda entre 0.9 y 1.5 veces la
  del 5 % en SN solas y BAO + SN, y entre 1.2 y 2.0 veces con el CMB comprimido.
- En la combinación con CMB el máximo de L está en σ_X > 0.
- La cota de −1.35 queda entre un 10 y un 20 % por debajo de la de −1.92.
- Tamaño efectivo de muestra en la cota: más de 100 caminos de 2000.

Exploración: no cambia `results/`, `figures/` ni ningún paper; **no cambia la regla del 5 % ni las cotas de `EXACTO_SN.md`**:
esto va al lado. Script `precision_sn.py`; salidas `precision_sn.{json,txt}` y `precision_sn_chi2.npz` (χ²_min de cada
realización, celda vieja y nueva, y combinación). Tiempo: 2 min (las 15 celdas nuevas, 91 s).


> Corrección documental (2026-11-12): las cotas lineales analíticas de `EXACTO_SN.md` citadas aquí (0.298 / 0.178 / 0.119 en θ = 1 / 3 / 10 y 0.567 en θ = 0.3) estaban
> referidas a un máximo recortado en σ_X ≤ 0.05; con el máximo global son 0.282 / 0.177 / 0.119 y 0.502 (`maximos_globales.json`). Las cifras de esta nota no se recalculan.

## Montaje

El de `exacto_sn.py` sin cambios, mismas semillas por celda. Añadido:
- **Celdas nuevas:** σ_X = 0.07, 0.15 y 0.4 para los cinco θ (semillas [SEED, 200 + índice]); malla de σ_X con 12 valores
  (0.01 … 0.7). En ellas se calcula también el χ²_min de BAO + θ* + ω_m sin SN, pero la columna sin SN de la tabla usa solo las
  celdas de `exacto_chi2.npz`, como se pidió.
- **Cota del 5 %** sobre la malla completa, y su error por bootstrap sobre las realizaciones (1000 remuestreos por
  celda, independientes entre celdas; percentiles 16 y 84 de la cota). Se da también la cota con las celdas viejas, que es la
  de `EXACTO_SN.md` (y la de `EXACTO.md` sin SN).
- **Verosimilitud marginal exacta por Monte Carlo:** L(σ_X, θ) = media sobre los 2000 caminos de exp(−(χ²_min − χ²_obs)/2),
  con los descartes por H² ≤ 0 a peso 0; ln L en la malla contando σ_X = 0 (L = 1); cota donde ln L cae 1.92 (y 1.35) desde
  su máximo sobre la malla, interpolando en ln σ_X; error con el mismo bootstrap; tamaño efectivo de muestra
  ESS = (Σw)²/Σw² y peso máximo relativo en la celda más próxima a la cota; dónde está el máximo de L.
- **Lineal marginalizada:** la de `EXACTO_SN.md` (tres combinaciones con SN) y la de `CMB_FONDO.md` (variante Λ, sin SN).

## Resultados

| combinación | θ | 5 %, celdas de EXACTO_SN | 5 %, todas las celdas (16–84 %) | marginal −1.92 (16–84 %) | marginal −1.35 | lineal −1.92 / −1.35 | máx ln L: σ_X (valor) | ESS en la cota marginal (celda) | aciertos en la celda del 5 % |
|---|---|---|---|---|---|---|---|---|---|
| SN solas | 0.1 | 0.641 | 0.641 (+0.026 −0.037) | sin cota | 0.545 | 0.772 / 0.540 | 0 (+0.00) | — | 89 (0.7) |
| SN solas | 0.3 | 0.376 | 0.368 (+0.011 −0.011) | 0.426 (+0.011 −0.014) | 0.309 | 0.447 / 0.312 | 0 (+0.00) | 80 (0.4) | 80 (0.4) |
| SN solas | 1 | 0.199 | 0.199 (+0.011 −0.008) | 0.240 (+0.016 −0.015) | 0.169 | 0.246 / 0.172 | 0 (+0.00) | 94 (0.2) | 99 (0.2) |
| SN solas | 3 | 0.131 | 0.133 (+0.006 −0.006) | 0.164 (+0.007 −0.008) | 0.116 | 0.147 / 0.103 | 0 (+0.00) | 55 (0.15) | 81 (0.15) |
| SN solas | 10 | 0.082 | 0.080 (+0.002 −0.003) | 0.103 (+0.004 −0.009) | 0.067 | 0.092 / 0.066 | 0 (+0.00) | 25 (0.1) | 131 (0.07) |
| BAO + SN | 0.1 | 0.598 | 0.598 (+0.021 −0.021) | 0.668 (+0.008 −0.070) | 0.505 | 0.833 / 0.725 | 0.2 (+1.13) | 18 (0.7) | 69 (0.7) |
| BAO + SN | 0.3 | 0.353 | 0.342 (+0.008 −0.009) | 0.408 (+0.018 −0.033) | 0.332 | 0.447 / 0.384 | 0.1 (+1.08) | 13 (0.4) | 137 (0.3) |
| BAO + SN | 1 | 0.192 | 0.190 (+0.004 −0.003) | 0.220 (+0.010 −0.019) | 0.177 | 0.233 / 0.197 | 0.05 (+0.87) | 24 (0.2) | 81 (0.2) |
| BAO + SN | 3 | 0.139 | 0.129 (+0.003 −0.003) | 0.138 (+0.006 −0.007) | 0.122 | 0.145 / 0.123 | 0.05 (+0.87) | 12 (0.15) | 59 (0.15) |
| BAO + SN | 10 | 0.089 | 0.086 (+0.003 −0.003) | 0.107 (+0.003 −0.005) | 0.100 | 0.104 / 0.090 | 0.03 (+0.64) | 10 (0.1) | 68 (0.1) |
| BAO + θ* + ω_m + SN | 0.1 | 0.476 | 0.476 (+0.011 −0.009) | sin cota | 0.640 | sin cota / 0.917 | 0.3 (+1.23) | — | 83 (0.5) |
| BAO + θ* + ω_m + SN | 0.3 | 0.307 | 0.306 (+0.012 −0.011) | 0.418 (+0.015 −0.037) | 0.365 | 0.567 / 0.505 | 0.2 (+1.44) | 9 (0.4) | 104 (0.3) |
| BAO + θ* + ω_m + SN | 1 | 0.173 | 0.167 (+0.004 −0.005) | 0.230 (+0.016 −0.032) | 0.191 | 0.298 / 0.263 | 0.1 (+0.95) | 17 (0.2) | 129 (0.15) |
| BAO + θ* + ω_m + SN | 3 | 0.113 | 0.111 (+0.004 −0.005) | 0.149 (+0.008 −0.014) | 0.130 | 0.178 / 0.154 | 0.05 (+0.96) | 10 (0.15) | 119 (0.1) |
| BAO + θ* + ω_m + SN | 10 | 0.081 | 0.075 (+0.003 −0.004) | 0.105 (+0.002 −0.004) | 0.101 | 0.119 / 0.101 | 0.05 (+0.53) | 8 (0.1) | 113 (0.07) |
| BAO + θ* + ω_m (sin SN) | 0.1 | sin cota | sin cota | sin cota | sin cota | sin cota / sin cota | 0.7 (+1.27) | — | — |
| BAO + θ* + ω_m (sin SN) | 0.3 | sin cota | sin cota | sin cota | sin cota | sin cota / sin cota | 0.3 (+1.38) | — | — |
| BAO + θ* + ω_m (sin SN) | 1 | 0.427 | 0.435 (+0.011 −0.013) | 0.520 (+0.013 −0.029) | 0.447 | 0.624 / 0.547 | 0.15 (+1.36) | 20 (0.5) | 123 (0.4) |
| BAO + θ* + ω_m (sin SN) | 3 | 0.286 | 0.286 (+0.008 −0.007) | 0.324 (+0.008 −0.020) | 0.281 | 0.378 / 0.329 | 0.15 (+1.15) | 21 (0.3) | 86 (0.3) |
| BAO + θ* + ω_m (sin SN) | 10 | 0.192 | 0.190 (+0.006 −0.005) | 0.265 (+0.024 −0.028) | 0.221 | 0.269 / 0.232 | 0.1 (+0.88) | 6 (0.3) | 85 (0.2) |

(Las cotas "5 %, celdas de EXACTO_SN" reproducen `EXACTO_SN.md` y `EXACTO.md`; "aciertos" = realizaciones con χ²_min ≤ χ²_obs
en la celda más próxima a la cota del 5 %, el tamaño de muestra de esa regla.)

**P(χ²_min ≤ χ²_obs) y ln L por celda** (12 valores de σ_X: 0.01, 0.02, 0.03, 0.05, 0.07, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.7):

SN solas (P / ln L):

| θ \ σ_X | 0.01 | 0.02 | 0.03 | 0.05 | 0.07 | 0.1 | 0.15 | 0.2 | 0.3 | 0.4 | 0.5 | 0.7 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.1 | 0.474 / -0.00 | 0.490 / -0.00 | 0.463 / -0.01 | 0.426 / -0.03 | 0.407 / -0.05 | 0.354 / -0.12 | 0.296 / -0.23 | 0.256 / -0.34 | 0.159 / -0.68 | 0.100 / -0.99 | 0.066 / -1.28 | 0.044 / -1.54 |
| 0.3 | 0.493 / -0.00 | 0.444 / -0.02 | 0.421 / -0.04 | 0.382 / -0.09 | 0.307 / -0.17 | 0.268 / -0.27 | 0.188 / -0.55 | 0.139 / -0.78 | 0.074 / -1.31 | 0.040 / -1.72 | 0.019 / -2.43 | 0.013 / -2.65 |
| 1 | 0.465 / -0.01 | 0.403 / -0.06 | 0.373 / -0.09 | 0.289 / -0.26 | 0.214 / -0.44 | 0.152 / -0.66 | 0.082 / -1.17 | 0.050 / -1.62 | 0.021 / -2.30 | 0.007 / -3.01 | 0.004 / -3.56 | 0.002 / -4.33 |
| 3 | 0.423 / -0.04 | 0.357 / -0.13 | 0.277 / -0.26 | 0.186 / -0.54 | 0.121 / -0.89 | 0.073 / -1.16 | 0.041 / -1.67 | 0.014 / -2.44 | 0.004 / -1.88 | 0.002 / -3.80 | 0.001 / -5.21 | 0.000 / -8.24 |
| 10 | 0.428 / -0.04 | 0.286 / -0.26 | 0.215 / -0.44 | 0.113 / -0.92 | 0.066 / -1.41 | 0.025 / -1.81 | 0.005 / -3.38 | 0.003 / -4.10 | 0.002 / -5.05 | 0.001 / -5.69 | 0.000 / -13.77 | 0.000 / -17.49 |

BAO + SN (P / ln L):

| θ \ σ_X | 0.01 | 0.02 | 0.03 | 0.05 | 0.07 | 0.1 | 0.15 | 0.2 | 0.3 | 0.4 | 0.5 | 0.7 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.1 | 0.471 / +0.01 | 0.483 / +0.07 | 0.485 / +0.14 | 0.460 / +0.33 | 0.478 / +0.56 | 0.427 / +0.77 | 0.368 / +0.99 | 0.328 / +1.13 | 0.189 / +0.92 | 0.112 / +0.41 | 0.068 / -0.20 | 0.035 / -0.89 |
| 0.3 | 0.502 / +0.05 | 0.473 / +0.18 | 0.473 / +0.34 | 0.448 / +0.65 | 0.410 / +0.85 | 0.355 / +1.08 | 0.227 / +0.91 | 0.154 / +0.88 | 0.069 / +0.01 | 0.028 / -0.77 | 0.011 / -1.59 | 0.004 / -2.92 |
| 1 | 0.498 / +0.16 | 0.454 / +0.36 | 0.440 / +0.71 | 0.374 / +0.87 | 0.281 / +0.82 | 0.202 / +0.86 | 0.097 / -0.07 | 0.041 / -0.78 | 0.007 / -1.94 | 0.001 / -2.16 | 0.001 / -4.69 | 0.001 / -3.66 |
| 3 | 0.476 / +0.28 | 0.410 / +0.60 | 0.358 / +0.77 | 0.226 / +0.87 | 0.139 / +0.55 | 0.086 / +0.44 | 0.029 / -1.45 | 0.010 / -2.99 | 0.001 / -4.18 | 0.000 / -12.13 | 0.000 / -17.94 | 0.000 / -46.56 |
| 10 | 0.460 / +0.31 | 0.334 / +0.52 | 0.270 / +0.64 | 0.131 / +0.50 | 0.071 / +0.34 | 0.034 / -0.70 | 0.002 / -4.25 | 0.000 / -8.73 | 0.000 / -12.83 | 0.000 / -10.82 | 0.000 / -27.15 | 0.000 / -53.63 |

BAO + θ* + ω_m + SN (P / ln L):

| θ \ σ_X | 0.01 | 0.02 | 0.03 | 0.05 | 0.07 | 0.1 | 0.15 | 0.2 | 0.3 | 0.4 | 0.5 | 0.7 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.1 | 0.470 / +0.01 | 0.462 / +0.05 | 0.474 / +0.12 | 0.439 / +0.25 | 0.434 / +0.42 | 0.367 / +0.60 | 0.300 / +0.64 | 0.231 / +1.10 | 0.131 / +1.23 | 0.080 / +0.87 | 0.042 / +0.02 | 0.025 / -0.17 |
| 0.3 | 0.485 / +0.04 | 0.470 / +0.16 | 0.450 / +0.28 | 0.403 / +0.54 | 0.337 / +0.60 | 0.270 / +1.02 | 0.157 / +1.18 | 0.112 / +1.44 | 0.052 / +0.62 | 0.023 / -0.16 | 0.011 / -1.78 | 0.004 / -0.36 |
| 1 | 0.482 / +0.15 | 0.438 / +0.33 | 0.418 / +0.61 | 0.303 / +0.83 | 0.217 / +0.86 | 0.144 / +0.95 | 0.065 / +0.75 | 0.025 / -0.64 | 0.007 / -1.61 | 0.001 / -6.40 | 0.001 / -4.08 | 0.000 / -14.82 |
| 3 | 0.474 / +0.28 | 0.401 / +0.53 | 0.325 / +0.70 | 0.182 / +0.96 | 0.101 / +0.85 | 0.059 / +0.68 | 0.021 / -0.98 | 0.005 / -2.06 | 0.001 / -3.75 | 0.000 / -12.14 | 0.000 / -16.04 | 0.000 / -55.54 |
| 10 | 0.455 / +0.31 | 0.317 / +0.44 | 0.243 / +0.47 | 0.107 / +0.53 | 0.057 / +0.14 | 0.024 / -0.75 | 0.001 / -5.64 | 0.000 / -7.94 | 0.000 / -14.19 | 0.000 / -16.74 | 0.000 / -25.65 | 0.000 / -94.95 |

BAO + θ* + ω_m (sin SN) (P / ln L):

| θ \ σ_X | 0.05 | 0.07 | 0.1 | 0.15 | 0.2 | 0.3 | 0.4 | 0.5 | 0.7 |
|---|---|---|---|---|---|---|---|---|---|
| 0.1 | 0.448 / +0.04 | 0.422 / +0.06 | 0.401 / +0.15 | 0.368 / +0.29 | 0.337 / +0.61 | 0.279 / +1.01 | 0.247 / +1.07 | 0.182 / +1.20 | 0.140 / +1.27 |
| 0.3 | 0.427 / +0.11 | 0.366 / +0.15 | 0.366 / +0.44 | 0.287 / +0.76 | 0.264 / +1.19 | 0.199 / +1.38 | 0.146 / +1.31 | 0.110 / +0.95 | 0.053 / +0.25 |
| 1 | 0.362 / +0.28 | 0.341 / +0.57 | 0.292 / +0.80 | 0.217 / +1.36 | 0.168 / +1.33 | 0.093 / +0.80 | 0.061 / +0.36 | 0.030 / -0.35 | 0.012 / -2.17 |
| 3 | 0.340 / +0.58 | 0.285 / +0.79 | 0.242 / +1.14 | 0.163 / +1.15 | 0.103 / +0.91 | 0.043 / -0.41 | 0.020 / -1.74 | 0.010 / -2.53 | 0.003 / -2.73 |
| 10 | 0.323 / +0.57 | 0.247 / +0.76 | 0.176 / +0.88 | 0.085 / -0.00 | 0.043 / -0.15 | 0.018 / -1.44 | 0.004 / -4.04 | 0.004 / -2.29 | 0.003 / -5.74 |


## Lectura

- **Error de Monte Carlo de las cotas del 5 %** (semianchura 16–84 % sobre la cota): entre el 1.9 % y el 4.9 %
  (mediana 3.2 %). Es menor que lo que predije porque en la celda de la cota hay ~100 aciertos de 2000, y
  además la interpolación en ln σ_X entre dos celdas promedia dos binomiales.
- **Celdas nuevas:** la cota del 5 % se mueve como mucho un 7 % (BAO + θ* + ω_m + SN θ = 10: -7 %, BAO + SN θ = 3: -7 %, BAO + SN θ = 10: -4 %);
  mediana del desplazamiento 1.6 %. Las celdas σ_X = 0.07 y 0.15 caen justo donde están los techos de θ = 3 y 10 y
  afinan la interpolación.
- **Marginal exacta −1.92 frente al 5 %:** sin CMB, ×1.06–1.28 (SN solas y BAO + SN, donde hay cota; en SN solas con θ = 0.1 no
  hay cota marginal hasta 0.7); con CMB, ×1.13–1.40. La marginal es siempre más floja que el 5 %: promedia los caminos que
  ajustan mejor que ΛCDM con peso exponencial, y la regla del 5 % solo los cuenta.
- **Máximo de L:** en SN solas está en σ_X = 0 en los cinco θ. **En todas las combinaciones con BAO está en σ_X > 0**, con el CMB
  (0.05–0.3, ln L = +0.5 a +1.4) y también sin él (BAO + SN: 0.03–0.2, +0.6 a +1.1): el residuo absorbe en promedio la tensión
  entre los BAO y las SN y la de DESI con Planck. Un Δln L de +1 es Δχ² ≈ 2, sin significación.
- **−1.35 frente a −1.92:** con CMB, -13 %, -17 %, -13 %, -5 %, -14 %, -13 %, -16 %; BAO + SN,
  -24 %, -19 %, -19 %, -11 %, -6 %; SN solas, -28 %, -30 %, -29 %, -35 %.
  Cuanto más plano es ln L alrededor del máximo (combinaciones con tensión), menos separa el umbral.
- **Tamaño efectivo de muestra en la cota marginal:** SN solas θ = 0.3: 80, SN solas θ = 1: 94, SN solas θ = 3: 55, SN solas θ = 10: 25, BAO + SN θ = 0.1: 18, BAO + SN θ = 0.3: 13, BAO + SN θ = 1: 24, BAO + SN θ = 3: 12, BAO + SN θ = 10: 10, BAO + θ* + ω_m + SN θ = 0.3: 9, BAO + θ* + ω_m + SN θ = 1: 17, BAO + θ* + ω_m + SN θ = 3: 10, BAO + θ* + ω_m + SN θ = 10: 8, BAO + θ* + ω_m (sin SN) θ = 1: 20, BAO + θ* + ω_m (sin SN) θ = 3: 21, BAO + θ* + ω_m (sin SN) θ = 10: 6.
  Nunca llega a 100: la media de exp(−Δχ²/2) la dominan los pocos caminos que ajustan mejor que ΛCDM (peso máximo relativo
  del 5 % al 38 %). La cota marginal por Monte Carlo es, con 2000 caminos, menos precisa que la del 5 %
  (bootstrap: 3.0–10.5 %), y en las combinaciones con tensión depende de 6–20 caminos.
- La **lineal marginalizada** queda por encima de la marginal exacta en casi todo (×1.0–1.4 sin CMB, ×1.1–1.4 con él), y por
  encima del 5 % en todo: el orden es 5 % < marginal exacta < lineal.

## Predicción: qué se cumplió

- **Error de Monte Carlo del 5 % entre el 4 y el 12 %: no se cumplió.** Es menor: 1.9–4.9 %.
- **Con las celdas nuevas ninguna cota se mueve más de un 15 %: cumplida** (máximo 7 %).
- **Marginal −1.92 entre 0.9 y 1.5 veces el 5 % sin CMB: cumplida** (1.06–1.28, salvo SN solas con θ = 0.1, sin cota marginal).
  **Entre 1.2 y 2.0 con el CMB: cumplida en BAO + θ* + ω_m + SN** (1.37, 1.38, 1.35, 1.40) **y no del todo sin SN**
  (1.20, 1.13, 1.39: θ = 3 queda en 1.13).
- **Con el CMB el máximo de L está en σ_X > 0: cumplida** en las ocho combinaciones con CMB que tienen cota; también lo está en
  BAO + SN, que no predije.
- **−1.35 entre un 10 y un 20 % por debajo de −1.92: cumplida en parte.** Con CMB, entre el 4 y el 17 %; BAO + SN, entre el 7 y el
  24 %; SN solas, 27–35 %, fuera.
- **ESS en la cota mayor que 100: no se cumplió.** Máximo 94 (SN solas, θ = 1); con CMB, 6–21.

**Lo que cambia para la v3.0:** nada en las cotas; sí en cómo citarlas: la regla del 5 % tiene un error de Monte Carlo del
2–5 % y la marginal exacta la afloja entre un 6 % y un 40 % con un ESS de 6–94, así que la cota de cabecera sigue siendo la del 5 % de
`EXACTO_SN.md` con su error de bootstrap al lado, y la marginal se cita como alternativa más floja.
