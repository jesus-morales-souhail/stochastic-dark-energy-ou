# ¿Converge la cota por verosimilitud marginal al aumentar el número de caminos? (exploración)

**Predicción, escrita antes de correr `convergencia_sn.py`:**
- El ESS crece aproximadamente en proporción a N (×10 con N × 10, dentro de un factor 2).
- Con N = 20000 la cota marginal −1.92 queda a menos de un 10 % de la de N = 2000
  (0.230 / 0.149 / 0.105 en θ = 1 / 3 / 10).
- La marginal sigue entre 1.25 y 1.45 veces la del 5 %.
- El máximo de ln L sigue en σ_X > 0 y por debajo de +2.
- La cota del 5 % con N = 20000 queda a menos de un 5 % de la de N = 2000.

Exploración: no cambia `results/`, `figures/` ni ningún paper; **no cambia las cotas de `EXACTO_SN.md` ni de
`PRECISION_SN.md`**. Script `convergencia_sn.py`; salidas `convergencia_sn.{json,txt}` y `convergencia_sn_chi2.npz`
(χ²_min de los 20000 caminos de cada celda). Tiempo: 10 min las 44 celdas (16 procesos) y 2 min el análisis.

## Montaje

El de `precision_sn.py` sin cambios: solo BAO + θ* + ω_m + SN (χ²_obs = 1406.327), los cinco θ, N = 20000 caminos por
celda con semillas nuevas ([SEED, 300 + índice], distintas de las de N = 2000), generados por bloques de 2000 dentro de cada
proceso; celdas σ_X = 0.03, 0.05, 0.07, 0.1, 0.15, 0.2, 0.3, 0.4 (más 0.5 y 0.7 para θ = 0.1 y 0.3), 44 celdas. Para cada N
∈ {2000, 5000, 10000, 20000} (los primeros N caminos): cota del 5 %, marginal −1.92 y −1.35 (ln L sobre la malla contando
σ_X = 0, interpolación en ln σ_X; si cruza antes del primer nodo, "< 0.03"), máximo de ln L y ESS = (Σw)²/Σw² en la celda
más próxima a cada cota. Bootstrap de las tres cotas con N = 20000 (1000 remuestreos por celda). Convergencia: cambio de
la marginal −1.92 entre 10000 y 20000 menor del 5 %.

## Resultados

**Evolución con N:**

| θ | N | 5 % | marginal −1.92 | marginal −1.35 | máx ln L (σ_X) | ESS en la celda de la cota marginal | ESS en la celda de la cota del 5 % |
|---|---|---|---|---|---|---|---|
| 0.1 | 2000 | 0.482 | sin cota | sin cota | +1.12 (0.2) | — | 16 (0.5; 92 aciertos) |
| 0.1 | 5000 | 0.488 | sin cota | sin cota | +1.10 (0.4) | — | 12 (0.5; 230 aciertos) |
| 0.1 | 10000 | 0.501 | sin cota | sin cota | +1.01 (0.2) | — | 29 (0.5; 501 aciertos) |
| 0.1 | 20000 | 0.498 | sin cota | sin cota | +1.15 (0.4) | — | 53 (0.5; 991 aciertos) |
| 0.3 | 2000 | 0.281 | 0.472 | 0.441 | +1.35 (0.3) | 6 (0.5, peso máx 27 %) | 2 (0.3; 80 aciertos) |
| 0.3 | 5000 | 0.287 | 0.454 | 0.430 | +1.20 (0.15) | 10 (0.5, peso máx 21 %) | 5 (0.3; 216 aciertos) |
| 0.3 | 10000 | 0.290 | 0.470 | 0.441 | +1.20 (0.15) | 18 (0.5, peso máx 11 %) | 14 (0.3; 447 aciertos) |
| 0.3 | 20000 | 0.292 | 0.496 | 0.446 | +1.17 (0.15) | 9 (0.5, peso máx 23 %) | 22 (0.3; 916 aciertos) |
| 1 | 2000 | 0.177 | 0.357 | 0.339 | +1.24 (0.3) | 1 (0.4, peso máx 91 %) | 5 (0.2; 65 aciertos) |
| 1 | 5000 | 0.176 | 0.335 | 0.317 | +1.26 (0.1) | 2 (0.3, peso máx 81 %) | 14 (0.2; 159 aciertos) |
| 1 | 10000 | 0.172 | 0.322 | 0.303 | +1.18 (0.1) | 2 (0.3, peso máx 72 %) | 55 (0.15; 661 aciertos) |
| 1 | 20000 | 0.170 | 0.325 | 0.306 | +1.14 (0.1) | 5 (0.3, peso máx 34 %) | 86 (0.15; 1302 aciertos) |
| 3 | 2000 | 0.103 | 0.149 | 0.129 | +1.02 (0.05) | 5 (0.15, peso máx 39 %) | 13 (0.1; 106 aciertos) |
| 3 | 5000 | 0.102 | 0.138 | 0.122 | +1.25 (0.07) | 9 (0.15, peso máx 24 %) | 23 (0.1; 258 aciertos) |
| 3 | 10000 | 0.102 | 0.149 | 0.129 | +1.07 (0.07) | 24 (0.15, peso máx 10 %) | 45 (0.1; 518 aciertos) |
| 3 | 20000 | 0.102 | 0.173 | 0.136 | +1.07 (0.07) | 31 (0.15, peso máx 10 %) | 98 (0.1; 1040 aciertos) |
| 10 | 2000 | 0.070 | 0.100 | 0.083 | +0.56 (0.03) | 4 (0.1, peso máx 36 %) | 22 (0.07; 98 aciertos) |
| 10 | 5000 | 0.072 | 0.096 | 0.084 | +0.61 (0.03) | 11 (0.1, peso máx 18 %) | 56 (0.07; 264 aciertos) |
| 10 | 10000 | 0.071 | 0.107 | 0.088 | +0.58 (0.03) | 18 (0.1, peso máx 19 %) | 88 (0.07; 520 aciertos) |
| 10 | 20000 | 0.073 | 0.113 | 0.098 | +0.60 (0.03) | 28 (0.1, peso máx 14 %) | 138 (0.07; 1091 aciertos) |

**ESS en celdas fijas** (la de la cota marginal con N = 20000 y sus vecinas):

| θ | celda σ_X | ESS, N = 2000 | 5000 | 10000 | 20000 | cociente 20000/2000 |
|---|---|---|---|---|---|---|
| 0.1 | 0.4 | 13.6 | 42.0 | 84.2 | 74.7 | 5.5 |
| 0.1 | 0.5 | 15.9 | 12.0 | 29.4 | 52.6 | 3.3 |
| 0.1 | 0.7 | 2.9 | 4.5 | 9.1 | 22.6 | 7.9 |
| 0.3 | 0.4 | 2.6 | 8.2 | 9.7 | 17.5 | 6.8 |
| 0.3 | 0.5 | 6.2 | 9.8 | 18.2 | 9.3 | 1.5 |
| 0.3 | 0.7 | 2.9 | 2.5 | 4.8 | 6.0 | 2.1 |
| 1 | 0.2 | 5.4 | 13.6 | 26.9 | 49.3 | 9.1 |
| 1 | 0.3 | 1.1 | 1.5 | 1.9 | 4.9 | 4.5 |
| 1 | 0.4 | 1.2 | 2.2 | 4.7 | 4.3 | 3.6 |
| 3 | 0.1 | 13.2 | 22.7 | 45.1 | 97.6 | 7.4 |
| 3 | 0.15 | 4.6 | 8.8 | 23.9 | 30.7 | 6.7 |
| 3 | 0.2 | 4.6 | 6.2 | 9.0 | 7.2 | 1.6 |
| 10 | 0.07 | 22.0 | 55.7 | 87.9 | 137.9 | 6.3 |
| 10 | 0.1 | 3.9 | 11.0 | 17.6 | 27.7 | 7.1 |
| 10 | 0.15 | 4.0 | 7.5 | 3.6 | 4.7 | 1.2 |

**Bootstrap con N = 20000 y convergencia:**

| θ | 5 % (16–84 %) | marginal −1.92 (16–84 %) | marginal −1.35 (16–84 %) | cambio 10000 → 20000 de la −1.92 | convergida (< 5 %) |
|---|---|---|---|---|---|
| 0.1 | 0.498 (0.493–0.506) | sin cota | sin cota | — | no |
| 0.3 | 0.292 (0.289–0.295) | 0.496 (0.468–0.516) | 0.446 (0.425–0.461) | +5.6 % | no |
| 1 | 0.170 (0.169–0.172) | 0.325 (0.307–0.332) | 0.306 (0.258–0.315) | +0.8 % | sí |
| 3 | 0.102 (0.101–0.104) | 0.173 (0.153–0.200) | 0.136 (0.129–0.142) | +16.0 % | no |
| 10 | 0.073 (0.072–0.074) | 0.113 (0.107–0.118) | 0.098 (0.092–0.102) | +5.8 % | no |

## Lectura

- **La cota del 5 % es estable:** con N = 20000 vale 0.170 (θ = 1), 0.102 (θ = 3), 0.073 (θ = 10), frente a 0.177, 0.103, 0.070 con los
  primeros 2000 de esta corrida (-3.7 %, -0.8 %, +5.0 %) y 0.167, 0.111, 0.075 en `PRECISION_SN.md`
  (+2.0 %, -7.6 %, -2.8 %); en θ = 0.3 y 0.1, 0.292 y 0.498 (-4.6 % y +4.7 %). El
  bootstrap con 20000 la deja en ±1 %.
- **La marginal por Monte Carlo no ha convergido y sube con N:** 0.230 → 0.325 (θ = 1, +41 %), 0.149 → 0.173 (θ = 3, +16 %), 0.105 → 0.113 (θ = 10, +7 %) frente a `PRECISION_SN.md`;
  entre 10000 y 20000 cambia +5.6 % (θ = 0.3), +0.8 % (θ = 1), +16.0 % (θ = 3), +5.8 % (θ = 10). Solo θ = 1 pasa el criterio del 5 %, y lo pasa con un ESS de 5 en la
  celda de la cota: no es convergencia, es que dos muestras de pocos caminos dieron lo mismo. El bootstrap con 20000 da
  ±4–14 % en la marginal, y en θ = 3 el intervalo 16–84 % (0.153–0.200) no contiene el valor de N = 2000 (0.149).
- **Por qué sube:** L(σ_X) = ⟨exp(−Δχ²/2)⟩ la dominan los caminos raros que ajustan mucho mejor que ΛCDM (peso máximo relativo
  del 10 al 40 % en la celda de la cota con N = 20000). Al aumentar N aparecen caminos más raros y mejores, L crece en σ_X
  grande y la cota se afloja. Es el comportamiento típico de una media de exponenciales con cola pesada: la estimación con
  N finito está sesgada hacia abajo y converge despacio.
- **ESS:** en las celdas de las cotas crece casi en proporción a N: de 2000 a 20000 (×10) sube ×5 (θ = 1, celda 0.3), ×7
  (θ = 3, celda 0.15) y ×7 (θ = 10, celda 0.1), dentro del factor 2 de ×10; en las celdas vecinas de σ_X mayor es errático
  (×1.2–3.6, con ESS de 4–7 incluso con 20000 caminos), porque ahí un solo camino se lleva hasta el 40 % del peso. El ESS
  absoluto en la cota sigue siendo pequeño: 5, 31 y 28 con N = 20000.
- **Marginal frente al 5 %** con N = 20000: ×1.70, 1.91, 1.69, 1.55 (θ = 0.3, 1, 3, 10), más que el ×1.34–1.40 de N = 2000; y
  el máximo de ln L sigue en σ_X > 0 (0.15, 0.1, 0.07, 0.03) con valores +1.17, +1.14, +1.07, +0.60, estables con N.

## Predicción: qué se cumplió

- **ESS ∝ N dentro de un factor 2: cumplida en las celdas de las cotas** (×5 / ×7 / ×7 en θ = 1 / 3 / 10, en el borde en
  θ = 1); no en las celdas de σ_X mayor (×1.2–3.6).
- **Marginal −1.92 con N = 20000 a menos del 10 % de la de N = 2000 (0.230 / 0.149 / 0.105): no se cumplió** en θ = 1 y 3
  (+41 %, +16 %); **sí** en θ = 10 (+8 %).
- **Marginal entre 1.25 y 1.45 veces el 5 %: no se cumplió** (×1.91 / 1.69 / 1.55).
- **Máximo de ln L en σ_X > 0 y por debajo de +2: cumplida** (+1.17 en 0.15, +1.14 en 0.1, +1.07 en 0.07, +0.60 en 0.03).
- **5 % con N = 20000 a menos del 5 % de la de N = 2000: cumplida** frente a los primeros 2000 de esta corrida
  (-3.7 %, -0.8 %, +5.0 %) y frente a `PRECISION_SN.md` en θ = 1 y 10 (+2.0 %, -2.8 %); **no** en θ = 3
  (-7.6 %, con semillas y malla distintas).

**Conclusión.** La cota por verosimilitud marginal estimada por Monte Carlo con pesos exp(−Δχ²/2) no ha convergido con 20000
caminos (sube un 8–41 % respecto a 2000 y sigue moviéndose entre 10000 y 20000 en θ = 0.3, 3 y 10), porque la media la
dominan unos pocos caminos de cola; no sirve como cota de cabecera. La regla del 5 % es estable al 1–8 % entre 2000 y 20000 (θ = 3: 0.111 → 0.102)
caminos y es la que se queda: 0.17 / 0.10 / 0.073 en θ = 1 / 3 / 10 con BAO + θ* + ω_m + SN.
