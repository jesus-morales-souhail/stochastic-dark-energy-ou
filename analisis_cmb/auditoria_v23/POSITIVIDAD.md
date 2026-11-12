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

## Versión 2 (diseño declarado antes de correr y antes de ver ninguna cota)

La versión 1 paró porque la malla de Ω_Λ (± 3σ) deja fuera el 0.27 % de un posterior gaussiano por construcción (fallo de diseño del
prompt). Las tres predicciones no cambian. Cambios de diseño, declarados antes de ver ninguna cota:

- Malla de Ω_Λ desde el mínimo hasta el máximo de Ω_Λ de TODAS las muestras del posterior de todas las celdas, más un 5 % de margen a
  cada lado (5 % del recorrido máx − mín), con 15 nodos. Comprobación: ninguna muestra fuera.
- Interpolación de P(σ, Ω_Λ) comprobada contra el cálculo directo con los mismos 20000 caminos del prior en 5 valores de Ω_Λ entre
  nodos, en las celdas de la cota; se da la diferencia máxima.
- Celda nueva σ_X = 0.6 en θ = 0.3.
- Si en un θ el ESS de la parte truncada es menor que el 1 % de N en la celda de la cota, ese θ se declara "no evaluable con este
  método" y no se da su cota; los demás θ siguen.

Todo lo demás, igual que la versión 1 (Z_pos con la normalización por muestra, bootstrap de P con 200 repeticiones, tabla, "Predicción"
punto por punto). Condición de parada: si alguna muestra queda fuera de la malla nueva, se para y se informa. La versión 1 del script
queda en el historial (commit 25d7a65); `positividad.json` y `positividad_stop.json` se sobrescriben.

## Resultados de la versión 2 (`positividad.py`, 109 s; `positividad.json`)

**Malla de Ω_Λ:** Ω_Λ de las muestras del posterior de las 43 celdas va de 0.4436 a 0.7598; con el 5 % de margen por lado la
malla de 15 nodos cubre 0.4278–0.7756; **ninguna de las 215000 muestras queda fuera** (condición de parada no saltó).
Ω_Λ del mejor ajuste ΛCDM: 0.6962. Interpolación de P frente al cálculo directo en 5 valores entre nodos, en las celdas de
las cotas: diferencia máxima 0.0014 (columna de la tabla).

**Cotas −1.92** (sin positividad: la de MARGINAL_GUIADA, de referencia, y la recalculada aquí con las muestras nuevas; con positividad;
cociente frente a la referencia / frente a la propia; P en la cota, interpolada en ln σ entre las dos celdas, con el error Monte Carlo mayor
de las dos; rango 16–84 % de la cota por el bootstrap de P; celdas que encierran la cota con su ESS truncado y sus muestras positivas):

| θ | sin positividad, referencia | sin positividad, aquí | con positividad | cociente ref. / propia | P en la cota | cota, rango 16–84 (bootstrap de P) | celdas | ESS truncado (fracción de N) | muestras positivas | interp. máx. dif. |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.3 | 0.500 | 0.507 | no evaluable | — | — | — | 0.4 y 0.5 | 32 (0.65 %) y 329 (6.59 %) | 2720 y 1964 | 0.0006 |
| 1 | 0.287 | 0.289 | no evaluable | — | — | — | 0.2 y 0.3 | 887 (17.75 %) y 44 (0.89 %) | 4343 y 2501 | 0.0008 |
| 3 | 0.177 | 0.212 | 0.177 | 1.000 / 0.836 | 0.855 ± 0.003 | 0.1773–0.1774 | 0.15 y 0.2 | 563 (11.26 %) y 342 (6.84 %) | 4923 y 4071 | 0.0014 |
| 10 | 0.116 | 0.113 | 0.113 | 0.973 / 1.005 | 0.975 ± 0.002 | 0.1131–0.1131 | 0.1 y 0.15 | 725 (14.50 %) y 203 (4.06 %) | 5000 y 4955 | 0.0010 |

- **θ = 0.3 y θ = 1: no evaluables con este método.** El ESS de la parte truncada en una de las dos celdas que encierran la cota es
  0.65 % de N (σ_X = 0.4, θ = 0.3) y 0.89 % (σ_X = 0.3, θ = 1), por debajo del 1 % fijado. Por la regla no se da su cota (el json
  conserva el número crudo en `bound_pos_192`, sin valor por esa regla). La celda nueva σ_X = 0.6 de θ = 0.3 no cambia eso: la cota cae
  entre 0.4 y 0.5, no entre 0.5 y 0.7.
- **θ = 3:** la cota con positividad, 0.177, coincide con la referencia sin positividad (0.177): cociente 1.000. La "sin positividad,
  aquí" (0.212) no sirve de comparación emparejada: en la celda σ_X = 0.2 el ESS de los pesos sin truncar es 2 (una muestra domina), y
  ln Z_exact relativo queda en -0.48 frente a -1.46 de MARGINAL_GUIADA; esa muestra tiene ρ_DE ≤ 0 y la truncación la quita, con lo que
  ln Z_pos en 0.2 (-1.46) recupera el valor de la referencia. Es decir: en θ = 3 la positividad no estrecha la cota; lo que hace es quitar
  una muestra atípica del cálculo sin truncar.
- **θ = 10:** 0.113 con positividad frente a 0.113 sin ella con las mismas muestras (cociente 1.005) y 0.116 de referencia
  (cociente 0.973). En la cota P = 0.98: casi todos los caminos del prior son positivos y la truncación no actúa.
- **Error de P frente a error de muestreo:** el bootstrap de P mueve la cota en menos de 0.0005 (cuarto decimal), pero la diferencia entre la
  cota sin positividad calculada aquí y la de MARGINAL_GUIADA (semillas distintas) es de 3 % en θ = 10 y de 20 % en θ = 3: el error que manda es
  el de los pesos (ESS de 200–700 en las celdas de la cota, y de 2 en la celda 0.2 de θ = 3), no el de P. Las cotas con positividad se dan con
  tres cifras por el rango del bootstrap, pero su error real es el de muestreo, del orden del 3 % o más.

## Predicción, punto por punto

1. **Sin evaluar.** θ = 1 no es evaluable con este método (ESS truncado 0.89 % de N en la celda 0.3, por debajo del 1 %).
2. **θ = 3: cumplida frente a la referencia** (cociente 1.000, baja 0 %; la comparación emparejada no es válida por el ESS = 2 de la celda 0.2).
   **θ = 10: cumplida en la comparación emparejada** (+0.5 %, es decir no baja) **y no frente a la referencia** (−2.7 %, por encima del 2 %
   predicho); la diferencia entre las dos comparaciones es el error de muestreo de la cota sin positividad, no la positividad.
3. **Sin evaluar** (requiere la cota nueva de θ = 1).

## Respuesta

Con este método, la positividad solo se puede evaluar en θ = 3 y θ = 10, y en ninguno de los dos estrecha la cota de la marginal −1.92:
0.177 frente a 0.177 (θ = 3) y 0.113 frente a 0.113 con las mismas muestras (θ = 10), con P_prior(ρ_DE > 0) en la cota de
0.86 y 0.98. En θ = 1 y θ = 0.3, donde el prior de positividad sí muerde (P = 0.60 en σ_X = 0.3 para θ = 1, 0.54–0.67 en 0.4–0.5 para
θ = 0.3), el peso truncado se concentra en muy pocas muestras (ESS < 1 % de N) y la importancia guiada por el posterior lineal no da una
cota fiable: haría falta muestrear el posterior truncado directamente. Las dos cotas (con y sin positividad) quedan dadas donde se pueden
dar; el error que manda en todas es el de muestreo, no el de P.

## Versión 3 (lineal; `positividad_lin.py`, salida `positividad_lin.json`; declarado antes de calcular)

- Hallazgo de la versión 2: la marginal exacta guiada es inestable cuando hay muestras con ρ_DE ≤ 0 (θ = 3, celda 0.2: ESS 2; cota sin
  positividad 0.177 → 0.212 solo con otras semillas). La cifra de cabecera pasa a ser la marginal lineal analítica: es estable y es la
  que validó INYECCION (cobertura 95.5–97.5 %). La exacta queda como comprobación, con esta inestabilidad dicha. Que una historia sea
  positiva no garantiza que la versión lineal sea buena para ella: se comprueba aparte.
- Método nuevo, declarado antes de ver cotas: Z_pos,lin(σ) = Z_lin(σ) · E_q[ 1(ρ_DE > 0) / P(σ, Ω_Λ) ], con q el posterior lineal
  (Matheron, sin pesos exacto/lineal), la positividad de cada muestra evaluada con su propio Ω_Λ (sus h y ω_m), y P(σ, Ω_Λ) la
  probabilidad de positividad del prior del camino a ese Ω_Λ (malla de Ω_Λ de la versión 2, interpolada). No se usa ninguna fracción
  del prior agrupada sobre Ω_Λ: los priors de (h, ω_m) son planos y no se promedia sobre ellos. 5000 muestras de posterior y 20000
  caminos de prior por nodo de Ω_Λ.
- Las tres predicciones originales no cambian; se evalúan contra la marginal lineal sin positividad (θ = 1: 0.272).
- Celdas: las de la versión 2. Cota: caída de 1.92 desde el máximo global de ln Z_pos,lin. Rango de la cota por bootstrap (200) sobre las
  muestras del posterior y sobre P (remuestreo conjunto en cada réplica; también se dan los dos por separado). Comprobación de la
  linealidad en las positivas: en la celda de la cota de θ = 1 y 3, mediana y percentiles 16–84 del cociente exacto/lineal solo entre
  las muestras con ρ_DE > 0. Se da también la fracción del posterior lineal con ρ_DE ≤ 0 en el techo sin positividad, para compararla
  con el 44.9 % de CMB_FINAL en θ = 1 (también del posterior lineal, en el techo exacto 0.287).
- Condición de parada: si en alguna celda de la cota hay menos de 100 muestras positivas del posterior, se para y se informa.

## Resultados de la versión 3 (`positividad_lin.py`, 105 s; `positividad_lin.json`)

Mismas muestras que la versión 2 (mismas semillas): malla de Ω_Λ 0.4278–0.7756, ninguna de las 215000 muestras fuera. La cota lineal sin
positividad reproduce la de MARGINAL_GUIADA en los cuatro θ (0.492, 0.272, 0.173, 0.115: es analítica). **Condición de parada: no saltó**
(mínimo de muestras positivas en una celda de la cota: 1964).

| θ | lineal sin positividad | lineal con positividad | cociente | fracción del posterior lineal con ρ_DE > 0 en la cota | P en la cota | rango 16–84 %, bootstrap conjunto | solo muestras | solo P |
|---|---|---|---|---|---|---|---|---|
| 0.3 | 0.492 | 0.461 | 0.937 (-6.3 %) | 0.448 | 0.586 ± 0.004 | 0.4590–0.4616 | 0.4592–0.4615 | 0.4601–0.4610 |
| 1 | 0.272 | 0.263 | 0.965 (-3.5 %) | 0.621 | 0.700 ± 0.003 | 0.2619–0.2631 | 0.2620–0.2631 | 0.2623–0.2627 |
| 3 | 0.173 | 0.174 | 1.008 (+0.8 %) | 0.896 | 0.867 ± 0.003 | 0.1739–0.1742 | 0.1740–0.1742 | 0.1740–0.1742 |
| 10 | 0.115 | 0.116 | 1.007 (+0.7 %) | 0.997 | 0.970 ± 0.002 | 0.1159–0.1160 | 0.1159–0.1159 | 0.1159–0.1159 |

Las cotas con positividad se dan con tres cifras (el rango del bootstrap está en la cuarta). Por qué sube en θ = 3 y 10 y baja en θ = 1 y
0.3: el factor E_q[1(ρ_DE > 0)/P] compara la fracción positiva del posterior con la del prior al mismo Ω_Λ; donde los datos hacen el
posterior más positivo que el prior (θ = 3, celda 0.2: 81.4 % frente a P = 0.778) el cociente pasa de 1 y la cota se afloja; donde lo hacen menos
positivo (θ = 1, celda 0.3: 50.0 % frente a 0.605) la cota se estrecha. En ningún caso el efecto pasa del 7 %.

**Fracción del posterior lineal con ρ_DE ≤ 0 en el techo** (interpolada en ln σ entre celdas), sin positividad: en el techo lineal y en el
exacto de MARGINAL_GUIADA, frente al recuento de CMB_FINAL (posterior lineal muestreado directamente en el techo exacto, otras semillas):

| θ | techo lineal | ρ_DE ≤ 0 | techo exacto | ρ_DE ≤ 0 (interp.) | CMB_FINAL |
|---|---|---|---|---|---|
| 0.3 | 0.492 | 59.6 % | 0.500 | 60.7 % | — |
| 1 | 0.272 | 41.1 % | 0.287 | 46.1 % | 44.9 % |
| 3 | 0.173 | 9.9 % | 0.177 | 11.4 % | 8.3 % |
| 10 | 0.115 | 0.3 % | 0.116 | 0.3 % | 0.0 % |

θ = 1: 46.1 % interpolado frente al 44.9 % de CMB_FINAL en 0.287 (entre las celdas 0.2, 13.1 %, y 0.3, 50.0 %): compatible. En θ = 3 la
interpolación en ln σ entre 1.5 % (0.15) y 18.6 % (0.2) sobreestima frente al 8.3 % medido en 0.177.

**Linealidad entre las positivas** (cociente exacto/lineal w = e^{−(χ²_ex − χ²_lin)/2} solo en las muestras con ρ_DE > 0, en las dos celdas
que encierran la cota de θ = 1 y θ = 3):

| θ | celda σ_X | positivas (fracción) | P | w mediana (16–84) | ESS de w entre positivas |
|---|---|---|---|---|---|
| 1 | 0.2 | 4343 (86.9 %) | 0.895 | 0.70 (0.204–1.55) | 898 |
| 1 | 0.3 | 2501 (50.0 %) | 0.605 | 0.45 (0.039–1.57) | 46 |
| 3 | 0.15 | 4923 (98.5 %) | 0.963 | 0.55 (0.097–1.55) | 567 |
| 3 | 0.2 | 4071 (81.4 %) | 0.778 | 0.33 (0.013–1.47) | 349 |

La mediana de w está en 0.3–0.7 (el χ² exacto queda 0.7–2.2 por encima del lineal en la muestra típica) y el rango 16–84 cubre un factor
8–100: ser positiva no hace a una historia "lineal". Lo que respalda la cifra lineal no es esa comprobación, sino la cobertura
95.5–97.5 % de INYECCION; la exacta guiada, con ESS de 46 entre las positivas de la celda 0.3 de θ = 1, queda como comprobación inestable.

## Predicción, punto por punto (contra la marginal lineal sin positividad)

1. **No cumplida.** La cota de θ = 1 baja un 3.5 % (0.272 → 0.263), no entre un 20 y un 40 %.
2. **Cumplida en magnitud, con el signo al revés.** θ = 3 cambia un +0.8 % (< 10 %) y θ = 10 un +0.7 % (< 2 %): las dos suben en vez de bajar.
3. **Cumplida por poco.** P_prior(ρ_DE > 0 en todo z | σ = 0.263, Ω_Λ = 0.6962) = 0.732 ± 0.003 (14640 de 20000), por encima del 70 %.

## Respuesta (versión 3)

Exigir ρ_DE > 0 en todo z casi no mueve las cotas de la marginal lineal −1.92: θ = 0.3, 0.492 → 0.461 (−6.3 %); θ = 1, 0.272 →
0.263 (−3.5 %); θ = 3, 0.173 → 0.174 (+0.8 %); θ = 10, 0.115 → 0.116 (+0.7 %), con errores de bootstrap en la cuarta cifra.
Eso ocurre aunque en el techo de θ = 1 el 41 % del posterior lineal tiene ρ_DE ≤ 0: la truncación quita esas historias pero la
normalización por el prior (P = 0.70 en la cota) devuelve casi todo, porque el prior a ese σ tampoco es positivo. El prior de positividad
cambia qué historias se admiten, no cuánto σ_X permiten los datos.
