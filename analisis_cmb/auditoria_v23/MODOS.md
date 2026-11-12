# MODOS: ¿los pocos modos que los datos miden bien explican cuánto acota cada familia en FAMILIAS?

Exploración (no toca PREPRINT ni CLAIMS). Script `modos.py`; salidas `modos_pred.json` (fase 1) y `modos.json` (fase 2).

Claude conocía las cotas de FAMILIAS al escribir esta predicción; por eso se fija la
fórmula y se calibra con una sola celda.

## Predicción (escrita antes de calcular)

1. Con σ(a_i) < 0.1 (en unidades de ρ_crit,0) hay entre 2 y 4 modos.
2. El primer modo vive a z bajo: el máximo de |e_1| está entre z = 0.05 y 0.4.
3. Con k = 2, la fórmula predice las cotas del 5 % de FAMILIAS para OU, Matérn y
   gaussiano dentro de un factor 1.3 en al menos el 80 % de las celdas con cota
   (sin contar la de calibración).
4. Con k = 2, lo mismo frente a las cotas lineales −1.92, dentro de un factor 1.2.
5. Todas las celdas "sin cota" de Matérn y gaussiano en FAMILIAS salen con
   σ_pred > 1.
6. Con el núcleo triangular y prior gaussiano la fórmula funciona como en 3; con los
   escalones de verdad falla en ℓ ≥ 3 (predice cota donde no la hay).
7. Pasar de 20 a 40 tramos cambia los tres primeros λ_i menos de un 10 %.
8. Marginalizar (h, ω_m, M) reduce la información del primer modo más de un 30 %.

## Fase 1 (sin abrir `familias.json` ni `FAMILIAS.md`)

### Montaje

- Respuesta lineal canónica A de `marginal_guiada.build_A` (BAO + θ* + ω_m + SN, 1595 filas, convenio B, malla de `exacto.py`,
  3506 nodos; la suma de cada fila es 0 a 2.7e-15). Fisher sobre las amplitudes de X por tramos: F = A_bᵀ P A_b,
  A_b = A M (M, indicador nodo → tramo), P = C⁻¹ − C⁻¹B(BᵀC⁻¹B)⁻¹BᵀC⁻¹ con B las derivadas de `build_B` (h, ω_m, M marginalizados);
  P = C⁻¹ para la versión sin marginalizar.
- Tramos: nb uniformes en u = ln(1+z) hasta z = 3 (u = 1.3863) más uno hasta z* (u = 6.995); nb = 20 (21 amplitudes,
  anchura 0.0693) y nb = 40 (41, anchura 0.0347); el último tramo tiene anchura 5.608.
- Normalización: D = diag(anchuras en u), F_norm = D^(−1/2) F D^(−1/2), K_norm = D^(1/2) K_b D^(1/2), con K_b la covarianza de los
  promedios de X en los tramos (media sobre los nodos de la malla, σ_X = 1). Comprobación tr(F K_b) = tr(F_norm K_norm): diferencia
  relativa máxima 1.2e-13 (nb = 20) y 1.6e-13 (nb = 40) sobre las 20 parejas familia × ℓ.
- Dirección nula: λ < 10⁻¹⁰ λ_max. La constante D^(1/2)·1 está en el espacio nulo (solape 1.000000; F_norm sobre
  ella es 2.3e-15 relativo) y se quita.
- Fórmula fijada σ_pred = c/√S_k, con S_k = Σ_(i≤k) λ_i v_i, v_i = e_iᵀ K_norm e_i. Dos constantes, calibradas en (OU, ℓ = 1) con
  k = 2 y nb = 20 y congeladas para todo k, celda y familia: S_2(OU, ℓ = 1) = 1048.3, **c_5% = 5.472** (contra 0.169),
  **c_lin = 8.807** (contra 0.272). (Por información, las constantes que calibrarían cada k por separado están en el json,
  `c_per_k_info_only`; no se usan.)

### Direcciones nulas: más de una

| nb | nulas (λ < 10⁻¹⁰ λ_max) | de ellas, la constante | las otras |
|---|---|---|---|
| 20 | 3 | 1 | 2: los tramos 18–20 (z ≥ 2.48, sin ningún nodo de BAO ni SN) solo entran por la fila de θ* → rango 1 de 3 |
| 40 | 8 | 1 | 7: los tramos 35–40 (z ≥ 2.36) solo por θ* (rango 1 de 6) y los tramos 31–33 (z 1.93–2.25) solo por D_M y D_H de Lyα, una SN en z = 2.26 y θ* (rango 3 de 10 columnas) |

Comprobado con el rango de las columnas de A_b más allá del último nodo: las únicas filas no nulas son la 13 (θ*) para z ≥ 2.48 y
las 11, 12, 13 y 1594 para z ≥ 2.25. Son degeneraciones del rebinado (tramos donde no hay datos, que solo entran por integrales),
no una propiedad nueva de los datos; pero **aparece más de una dirección nula** en las dos mallas, condición de parada.

### Espectro y modos (marginalizado)

| nb | λ₁…λ₆ | σ(a_i) = λ_i^(−1/2), i = 1…6 | modos con σ(a_i) < 0.1 |
|---|---|---|---|
| 20 | 88275, 13191, 7291.9, 4507, 3645.7, 3603.5 | 0.00337, 0.00871, 0.0117, 0.0149, 0.0166, 0.0167 | 11 de 18 |
| 40 | 1.6895e+05, 28564, 10453, 9054, 7737.9, 6733.4 | 0.00243, 0.00592, 0.00978, 0.0105, 0.0114, 0.0122 | 17 de 33 |

σ(a_i) es la desviación de la amplitud normalizada a_i = e_iᵀ D^(1/2) a (coeficiente en L²(du), unidades ρ_crit,0 × √u), que es lo
que define la normalización pedida.

Modos con nb = 20 (amplitud por tramo a = D^(−1/2) e, signo arbitrario, aquí positivo en el máximo):

- **Modo 1** (λ = 8.83e+04): +3.4 en el primer tramo (z 0–0.072) frente a −1.15, −1.0, −0.64, −0.26, −0.16 en z 0.072–0.52 y
  ≈ 0 después; un cambio de signo; el 0.1 % del peso e² por encima de z = 0.5. Es el contraste entre X en z → 0 y X en
  0.07 ≲ z ≲ 0.5: la pendiente de Δρ en el origen, que miden las 587 SN de z < 0.072.
- **Modo 2** (λ = 1.32e+04): +3.24 en z 0.072–0.149 (y +0.57 en el primer tramo) frente a −1.54, −0.79, −0.59 en z 0.23–0.62;
  0.4 % del peso por encima de z = 0.5.
- **Modo 3** (λ = 7.29e+03): +2.16 en z 0.41–0.52 y +0.89 en 0.32–0.41 frente a −2.1 y −1.65 en z 0.15–0.32 y −0.83 en el
  primer tramo; 6.6 % del peso por encima de z = 0.5.

Con nb = 40 el modo 1 es +5.07 en z 0–0.035 frente a −0.6 … −0.9 en z 0.035–0.3; el modo 2, +4.73 en z 0.035–0.072; el modo 3,
+2.96 y +3.08 en z 0.072–0.149. Los modos siguen al primer tramo: cuando el tramo se estrecha, el modo se estrecha con él.

### 20 frente a 40 tramos: los λ no convergen

λ(40)/λ(20) para i = 1…6: 1.91, 2.17, 1.43, 2.01, 2.12, 1.87; diferencia máxima en los tres primeros 117 %.
Solapes modo a modo en L²(du): 0.67, 0.41, 0.24, 0.22, 0.051, 0.7; grupos de λ a menos del 20 % entre sí: [[0], [1], [2, 3, 4, 5]], con ángulos
principales [48.1]°, [66.1]°, [84.8, 44, 42.1, 28.8]°. Causa: en el convenio B el modelo depende de X(z) − X(0), y X(0) es una evaluación
puntual, no una integral: la información sobre el primer tramo no crece con su anchura, así que F_norm = D^(−1/2) F D^(−1/2) va como
1/anchura en esa dirección y λ se duplica al partir los tramos. **Los tres primeros λ difieren más del 20 %**: condición de parada.

### Marginalizar (h, ω_m, M)

λ₁ sin marginalizar 1.674e+06 → marginalizado 8.828e+04: reducción del 94.7 % (nb = 20; 95.6 % con 40). La información
a lo largo del primer modo sin marginalizar cae un 95.0 %; el solape entre los dos primeros modos es 0.971; la traza
total se reduce al 7.2 %. Sin marginalizar hay 11 modos con σ(a_i) < 0.1 (nb = 20); el espacio nulo tiene las mismas
3 direcciones (son del rebinado, no de los parámetros).

### S_k y σ_pred (nb = 20, constantes congeladas; fase 1, sin comparar con nada)

| familia | ℓ | S_1, S_2, S_3, S_all | S_3/S_all | σ_pred 5 % (k = 1 / **2** / 3 / all) | σ_pred lineal −1.92 (k = 1 / **2** / 3 / all) |
|---|---|---|---|---|---|
| OU | 0.1 | 3357, 4088, 4651, 5328 | 0.87 | 0.094 / **0.086** / 0.080 / 0.075 | 0.152 / **0.138** / 0.129 / 0.121 |
| OU | 0.3 | 2098, 2605, 3083, 3653 | 0.84 | 0.119 / **0.107** / 0.099 / 0.091 | 0.192 / **0.173** / 0.159 / 0.146 |
| OU | 1 | 838.1, 1048, 1273, 1563 | 0.81 | 0.189 / **0.169** / 0.153 / 0.138 | 0.304 / **0.272** / 0.247 / 0.223 |
| OU | 3 | 306.4, 384, 470.6, 588.9 | 0.80 | 0.313 / **0.279** / 0.252 / 0.225 | 0.503 / **0.449** / 0.406 / 0.363 |
| OU | 10 | 95.09, 119.2, 146.7, 185.4 | 0.79 | 0.561 / **0.501** / 0.452 / 0.402 | 0.903 / **0.807** / 0.727 / 0.647 |
| Matern32 | 0.1 | 3746, 4611, 5298, 6059 | 0.87 | 0.089 / **0.081** / 0.075 / 0.070 | 0.144 / **0.130** / 0.121 / 0.113 |
| Matern32 | 0.3 | 1586, 2014, 2490, 3045 | 0.82 | 0.137 / **0.122** / 0.110 / 0.099 | 0.221 / **0.196** / 0.176 / 0.160 |
| Matern32 | 1 | 289, 370, 484.7, 658.5 | 0.74 | 0.322 / **0.284** / 0.249 / 0.213 | 0.518 / **0.458** / 0.400 / 0.343 |
| Matern32 | 3 | 41.35, 52.84, 70.95, 104.4 | 0.68 | 0.851 / **0.753** / 0.650 / 0.536 | 1.370 / **1.212** / 1.046 / 0.862 |
| Matern32 | 10 | 4.11, 5.233, 7.125, 11.43 | 0.62 | 2.699 / **2.392** / 2.050 / 1.619 | 4.344 / **3.850** / 3.299 / 2.605 |
| Gauss | 0.1 | 3700, 4613, 5366, 6136 | 0.87 | 0.090 / **0.081** / 0.075 / 0.070 | 0.145 / **0.130** / 0.120 / 0.112 |
| Gauss | 0.3 | 1116, 1434, 1846, 2363 | 0.78 | 0.164 / **0.144** / 0.127 / 0.113 | 0.264 / **0.233** / 0.205 / 0.181 |
| Gauss | 1 | 136.3, 174.3, 235.1, 342 | 0.69 | 0.469 / **0.414** / 0.357 / 0.296 | 0.754 / **0.667** / 0.574 / 0.476 |
| Gauss | 3 | 15.78, 20.11, 27.42, 43.52 | 0.63 | 1.377 / **1.220** / 1.045 / 0.829 | 2.217 / **1.964** / 1.682 / 1.335 |
| Gauss | 10 | 1.434, 1.82, 2.501, 4.296 | 0.58 | 4.570 / **4.056** / 3.460 / 2.640 | 7.355 / **6.528** / 5.569 / 4.249 |
| Triangular | 0.1 | 4169, 4981, 5509, 6234 | 0.88 | 0.085 / **0.078** / 0.074 / 0.069 | 0.136 / **0.125** / 0.119 / 0.112 |
| Triangular | 0.3 | 2903, 3638, 4320, 5035 | 0.86 | 0.102 / **0.091** / 0.083 / 0.077 | 0.163 / **0.146** / 0.134 / 0.124 |
| Triangular | 1 | 964.3, 1210, 1489, 1871 | 0.80 | 0.176 / **0.157** / 0.142 / 0.126 | 0.284 / **0.253** / 0.228 / 0.204 |
| Triangular | 3 | 321.6, 403.4, 496.9, 627.9 | 0.79 | 0.305 / **0.272** / 0.245 / 0.218 | 0.491 / **0.438** / 0.395 / 0.351 |
| Triangular | 10 | 96.5, 121, 149.1, 189.5 | 0.79 | 0.557 / **0.497** / 0.448 / 0.397 | 0.896 / **0.801** / 0.721 / 0.640 |

Los tres primeros modos llevan el 58–88 % de tr(F K_b); la fracción baja con ℓ y con la suavidad (gaussiano ℓ = 10: 58 %).
Las mismas tablas con nb = 40 están en el json; como los λ no convergen, tampoco S_k, aunque S_all = tr(F K_b) sí es independiente
del rebinado (es la traza en la malla, salvo el promedio por tramos).

### Condiciones de parada

Las dos se dan: los tres primeros λ difieren ×1.4–2.2 entre 20 y 40 tramos, y hay 3 y 8 direcciones nulas (constante + 2, + 7).
**La fase 2 no se ha corrido.** `modos_pred.json` queda guardado para que, si se decide seguir (con otra normalización, otros tramos,
o tal cual), la comparación se haga contra lo ya predicho.

sha256 de `modos_pred.json`: `f44c8b1da3acdb6548d9a1324306dca9c80bcc57c18d0d793f37bea7888363f6`

### Predicción, estado tras la fase 1

1. **No** (con la normalización pedida): 11 modos con σ(a_i) < 0.1 (nb = 20) y 17 (nb = 40), no 2–4.
2. **No**: el máximo de |e₁| está en el primer tramo, z 0–0.072 (nb = 20) y z 0–0.035 (nb = 40), por debajo de 0.05.
3–6. Pendientes de la fase 2.
7. **No**: los tres primeros λ cambian ×1.91, ×2.17 y ×1.43 (no < 10 %).
8. **Sí**: la información del primer modo cae un 95 % (> 30 %).

## Parada

Cierre del test sin fase 2 (decidido tras la fase 1; nada se ha recalculado).

### Las dos condiciones que saltaron

**1. Los tres primeros λ difieren más del 20 % entre 20 y 40 tramos** (umbral de la condición: 20 %; medido, hasta 117 %).

| | 20 + 1 tramos | 40 + 1 tramos | λ(40)/λ(20) | solape modo a modo en L²(du) |
|---|---|---|---|---|
| λ_1 | 88275 | 1.6895e+05 | 1.91 | 0.67 |
| λ_2 | 13191 | 28564 | 2.17 | 0.41 |
| λ_3 | 7291.9 | 10453 | 1.43 | 0.24 |
| λ_4 | 4507 | 9054 | 2.01 | 0.22 |
| λ_5 | 3645.7 | 7737.9 | 2.12 | 0.05 |
| λ_6 | 3603.5 | 6733.4 | 1.87 | 0.70 |

Grupos de λ a menos del 20 % entre sí (en cualquiera de las dos mallas) y ángulos principales entre los subespacios que forman:
grupo [0]: 48°; grupo [1]: 66°; grupo [2, 3, 4, 5]: 85°, 44°, 42°, 29°. Ni los modos individuales (solapes 0.05–0.70) ni los subespacios (ángulos 29–85°) coinciden entre mallas.

**2. Más de una dirección nula** (λ < 10⁻¹⁰ λ_max).

| nb | nulas | constante | otras | de dónde salen |
|---|---|---|---|---|
| 20 | 3 | 1 (solape 1.0000) | 2 | tramos 18–20 (z ≥ 2.48): sin nodo de BAO ni SN, solo entran por la fila de θ*; rango 1 de 3 columnas |
| 40 | 8 | 1 (solape 1.0000) | 7 | tramos 35–40 (z ≥ 2.36): solo θ*, rango 1 de 6; tramos 31–33 (z 1.93–2.25): solo D_M y D_H de Lyα, una SN en z = 2.26 y θ*, rango 3 de 10 columnas (31–40) |

Las direcciones nulas extra son del rebinado (tramos sin datos que solo entran por integrales), no de los datos; el espacio nulo sin
marginalizar tiene las mismas 3 y 8 direcciones.

### La causa

En el convenio B el modelo depende de X(z) − X(0). X(0) es un valor puntual, no una integral: la información sobre el primer tramo
no crece con su anchura, y la normalización F_norm = D^(−1/2) F D^(−1/2) solo es correcta para respuestas integrales (operadores con
núcleo en L², donde los autovalores convergen al refinar la malla). Con un valor puntual, λ va como 1/anchura en esa dirección y los
modos se pegan al primer tramo (modo 1: +3.4 en z 0–0.072 con 20 tramos, +5.07 en z 0–0.035 con 40). Fue un fallo de diseño del
prompt, no un resultado sobre los datos.

### Estado de las ocho predicciones

| # | estado |
|---|---|
| 1 | no cumplida por el fallo de diseño: 11 y 17 modos con σ(a_i) < 0.1 en la normalización pedida; no evaluable en sentido físico |
| 2 | no cumplida por el fallo de diseño: el máximo de \|e₁\| está en el primer tramo (z 0–0.072 y 0–0.035), que se mueve con la malla; no evaluable en sentido físico |
| 3 | sin evaluar (fase 2 no corrida) |
| 4 | sin evaluar |
| 5 | sin evaluar |
| 6 | sin evaluar |
| 7 | no cumplida por el fallo de diseño: λ₁…λ₃ cambian ×1.91, ×2.17, ×1.43; no evaluable en sentido físico |
| 8 | cumplida: λ₁ cae un 94.7 % al marginalizar (h, ω_m, M) (1.674e+06 → 8.828e+04); 95.6 % con 40 tramos |

Que la predicción 8 se cumpla numéricamente no valida el diseño ni es un resultado físico nuevo: solo confirma que el fondo absorbe
casi toda la información cerca de hoy (el primer modo es el contraste entre X(0) y X en 0.07 ≲ z ≲ 0.5, y (h, ω_m, M) se lo llevan
en un 95 %).

`familias.json` y `FAMILIAS.md` no se abrieron en ninguna fase; `modos_pred.json` (sha256 arriba) queda como registro de lo predicho.

Un rediseño (variable Y = X − X(0), con Y(0) = 0, tramos hasta z = 2.33 más uno a z*) sería un test nuevo, MODOS_2, con
predicciones nuevas.
