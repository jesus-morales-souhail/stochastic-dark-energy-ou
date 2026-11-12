# POSITIVIDAD: ¿cuánto se estrechan las cotas de la marginal −1.92 si se exige ρ_DE > 0 en todo z?

Exploración (no toca PREPRINT ni CLAIMS). Script `positividad.py`; salida `positividad.json`.

Pregunta: ¿cuánto se estrechan las cotas de la marginal −1.92 si se exige ρ_DE > 0 en
todo z (prior físico)?

Exigir ρ_DE > 0 es una elección de prior, no una ley: la relatividad general admite
energía oscura total negativa mientras H² > 0. Que H² > 0 no basta para que una
historia con ρ_DE < 0 sea físicamente aceptable: también importan la conservación, la
estabilidad de las perturbaciones y la teoría que la genere. Se dan las dos cotas.

## Predicción (escrita antes de calcular nada)

1. Con positividad, la cota de θ = 1 baja entre un 20 y un 40 % (de 0.287).
2. La de θ = 3 baja menos de un 10 %; la de θ = 10, menos de un 2 %.
3. P_prior(ρ_DE > 0 en todo z | σ = cota nueva de θ = 1, Ω_Λ del mejor ajuste) es
   mayor que el 70 %.

## Método (del prompt)

BAO + θ* + ω_m + Pantheon+, OU θ = 0.3, 1, 3, 10, las celdas de σ_X de marginal_guiada.py. Los priors de (h, ω_m, M) son planos,
así que no se promedia sobre ellos: la truncación define un prior del camino para cada Ω_Λ, con su propia normalización
P(σ, Ω_Λ) = P_prior(ρ_DE > 0 | σ, Ω_Λ).
- Z_pos(σ) = Z_lin(σ) · E_q[ w · 1(ρ_DE > 0) / P(σ, Ω_Λ) ], con q el posterior lineal (Matheron, como en marginal_guiada.py), w el
  cociente exacto/lineal, y la positividad y P evaluadas con el Ω_Λ propio de cada muestra (sus h y ω_m).
- P(σ, Ω_Λ) con 20000 caminos del prior por punto, en una malla de 7 valores de Ω_Λ (mejor ajuste ± 3σ de su error) e interpolada;
  misma malla en z, mismo X(0), mismo borde y misma definición de positividad que el numerador. Se guarda el número de caminos
  positivos por punto y su error Monte Carlo. Se comprueba que el Ω_Λ de todas las muestras del posterior cae dentro de la malla; no se
  extrapola.
- Cota: caída de 1.92 desde el máximo global de ln Z_pos (malla que contenga el máximo). Por celda: ESS de la parte truncada y número de
  muestras con ρ_DE > 0.
- Sensibilidad de la cota al error de P: cota rehecha con P remuestreada por bootstrap (200 veces), rango 16–84 % de la cota; no se
  dan más decimales de los que ese rango permite.
- Condición de parada: si el ESS de la parte truncada baja del 1 % de N en la celda de la cota, o alguna muestra del posterior cae
  fuera de la malla de Ω_Λ, se para y se informa.

## Parada (2026-11-12): la segunda condición de parada saltó antes de calcular ninguna cota

`positividad.py` corrió el muestreo completo (42 celdas, 110 s; `positividad.json`, lista de condiciones en `positividad_stop.json`)
y se detuvo antes de las cotas: **en las 42 celdas hay muestras del posterior cuyo Ω_Λ cae fuera de la malla de 7 valores
(mejor ajuste de la celda ± 3σ de su error)**, 717 de 210 000 en total, entre 7 y 40 por celda (0.14–0.80 % de las 5000). No se ha
extrapolado P ni calculado cota alguna; las predicciones quedan sin evaluar.

Por qué salta: el posterior lineal de (h, ω_m) es gaussiano exacto (δp ~ N(δp̂, Σ)), así que Ω_Λ(h, ω_m) es gaussiano salvo la
curvatura de h⁻², y una malla de ± 3σ deja fuera el 0.27 % de las muestras por construcción: 13.5 esperadas de 5000 por celda, observadas
7–40 (más en las celdas de σ_X alto, donde la cola de h⁻² es más larga: |z| máximo hasta 5.6). El error de Ω_Λ por el método delta coincide
con la desviación típica de las muestras a menos de un 3 %: no es un fallo del error, es la malla. Con 5000 muestras gaussianas una
malla de ± 3σ no puede contener a todas; la condición "ninguna muestra fuera de la malla" era imposible de cumplir con la malla fijada.

Diagnósticos guardados por celda (sin cota): número de muestras con ρ_DE > 0, ESS de los pesos truncados v = w·1(ρ_DE > 0)/P, P en los
nodos con su recuento y error Monte Carlo, y la comprobación de la interpolación (P interpolada frente a P directa con los mismos 20000
caminos en las 200 primeras muestras: diferencia máxima 0.025, en θ = 10 con σ_X ≥ 0.2; ≤ 0.01 en el resto). En las dos celdas que
encierran la cota sin positividad de MARGINAL_GUIADA:

| θ | fuera de la malla por celda | fracción | \|z\| máx | cota sin positividad (MARGINAL_GUIADA) | celdas que la encierran | ESS truncado (fracción de N) | muestras con ρ_DE > 0 | P en el centro de la malla |
|---|---|---|---|---|---|---|---|---|
| θ = 0.3 | 9–24 | 0.18–0.48 % | 4.9 | 0.500 | 0.5 y 0.7 | 131 (2.6 %) y 39 (0.8 %) | 1987 y 1173 | 0.534 y 0.369 |
| θ = 1 | 7–24 | 0.14–0.48 % | 5.1 | 0.287 | 0.2 y 0.3 | 1165 (23.3 %) y 162 (3.2 %) | 4348 y 2572 | 0.901 y 0.610 |
| θ = 3 | 13–26 | 0.26–0.52 % | 4.9 | 0.177 | 0.15 y 0.2 | 840 (16.8 %) y 236 (4.7 %) | 4922 y 4055 | 0.964 y 0.782 |
| θ = 10 | 10–40 | 0.20–0.80 % | 5.6 | 0.116 | 0.1 y 0.15 | 187 (3.7 %) y 76 (1.5 %) | 5000 y 4952 | 1.000 y 0.917 |

Aviso para el rediseño: en θ = 0.3 la cota (0.500) cae entre las celdas 0.5 y 0.7, y en 0.7 el ESS truncado ya está por debajo del 1 % de N
(0.8 %): la primera condición de parada también saltaría ahí con la malla de σ_X actual. En θ = 10 la celda 0.15 queda en el 1.5 %.

Qué cambiaría si se vuelve a lanzar (decisión del autor, no hecha): (a) definir la malla de Ω_Λ para que contenga todas las muestras
(por ejemplo [mín, máx] de las muestras de la celda, o ± 6σ con 13 nodos), o (b) evaluar P directamente en el Ω_Λ de cada muestra con
los mismos 20000 caminos del prior (P(σ, Ω_Λ) = fracción con mín(x − x(0)) > −Ω_Λ, que es exacta para cualquier Ω_Λ sin malla ni
interpolación; la malla solo serviría para la tabla), y (c) añadir una celda de σ_X entre 0.5 y 0.7 para θ = 0.3 o aceptar que esa cota
no se da. El resto del método (muestreo, pesos, positividad, bootstrap de P) está implementado y corrido; no se cambió nada después de
ver los resultados.
