# ¿La verosimilitud marginal exacta coincide con la lineal analítica? Muestreo de importancia guiado (exploración)

**Predicción, escrita antes de correr `marginal_guiada.py`:**
- ESS mayor que el 10 % de N en θ = 1, 3 y 10, en las celdas de las cotas.
- La corrección ln(Z_exacta / Z_lineal) es menor que 0.5 en valor absoluto en las
  celdas de las cotas lineales.
- La cota marginal exacta (−1.92) queda a menos de un 15 % de la lineal analítica
  (0.298 / 0.178 / 0.119 en θ = 1 / 3 / 10).
- La cota exacta queda entre 1.4 y 1.9 veces la del 5 % (0.170 / 0.102 / 0.073).

Exploración: no cambia `results/`, `figures/` ni ningún paper ni las cotas de las notas anteriores. Script `marginal_guiada.py`;
salidas `marginal_guiada.{json,txt}`. Tiempo: 1.5 min (un solo proceso; todo vectorizado sobre las muestras).

## Montaje

Combinación BAO + θ* + ω_m + SN (1595 filas: 13 BAO, θ*, ω_m, 1580 SN), θ ∈ {0.3, 1, 3, 10}, celdas σ_X de
`convergencia_sn.py` (0.03 … 0.4; 0.5 y 0.7 en θ = 0.3) más σ_X = 0 y 0.01 como comprobación, N = 5000 muestras por celda.

1. **Modelo lineal del convenio B** alrededor del mejor ajuste ΛCDM de estos datos, p₀ = (h, ω_m, M) = (0.68118, 0.14092,
   23.7979), χ² = 1406.327: d = b₀ + B δp + A x + n. x es el camino en la malla de `exacto.py` (3506 nodos, los z de
   DESI exactos), la misma en la que se evalúa el modelo exacto; A es la respuesta lineal completa en esa malla y en p₀
   (δy = y δlnD en los BAO y en θ*; (5/ln 10) δlnD_M interpolada entre nodos en las SN; cero en ω_m), con las filas sumando
   cero (Λ absorbe X(0); máximo |Σ| = 3e-15); comprobada frente a diferencias finitas del modelo exacto con un camino
   de σ_X = 0.02: 2 × 10⁻⁴ relativo (5 × 10⁻⁵ en los BAO). B son las derivadas del modelo exacto respecto a (h, ω_m, M) por
   diferencias finitas; prior x ∼ N(0, σ_X² K_θ) con K el núcleo OU en la malla; priores planos en δp, iguales en todas las
   celdas; n ∼ N(0, C), C = diag(oficial 13 × 13, σ_θ², σ_ω², C_SN STAT+SYS).
2. **Muestreo del posterior lineal:** δp de su marginal gaussiana N(δp̂, (BᵀS⁻¹B)⁻¹); x por la regla de Matheron,
   x = x₀ + σ_X² K Aᵀ S⁻¹ (d − b₀ − B δp − A x₀ − n₀), con x₀ del prior (AR(1) exacta), n₀ del ruido, S = C + σ_X² A K Aᵀ.
3. **Peso:** w = exp(−(χ²_exacto − χ²_lineal)/2); χ²_exacto es el modelo no lineal de `exacto_sn.py` en (p₀ + δp, x), sin
   reajustar; χ²_lineal el modelo lineal en el mismo punto. Los dos usan los mismos datos, la misma C y el mismo tratamiento
   de (h, ω_m, M): M es un parámetro muestreado en los dos, sin marginalización analítica en ninguno. H² ≤ 0 pesa 0
   (no ocurrió en ninguna celda).
4. ln Z_exacta(σ_X) = ln Z_lineal(σ_X) + ln⟨w⟩, con ln Z_lineal la marginal analítica (priores planos en δp). Cotas: caída de
   1.92 y 1.35 desde el máximo sobre la malla contando σ_X = 0, interpolando en ln σ_X; ESS = (Σw)²/Σw²; error por bootstrap
   de los pesos (1000 remuestreos).

**Comprobaciones.**
- Con X = 0 y δp = 0, χ²_exacto = 1406.3271 y χ²_lineal = 1406.3271: diferencia 0.0e+00.
- En σ_X = 0.01 los pesos son ≈ 1: ln w entre -0.05 y +0.03 (16–84 %), ESS ≥ 4979 de 5000; en σ_X = 0
   (solo δp), ln⟨w⟩ = +0.005: la linealización en (h, ω_m, M) es exacta a ese nivel.
- ln Z_lineal frente a la lineal analítica de `EXACTO_SN.md`: 0.492 frente a 0.567 (θ = 0.3, -13 %), 0.272 frente a 0.298 (θ = 1, -9 %), 0.173 frente a 0.178 (θ = 3, -3 %), 0.115 frente a 0.119 (θ = 10, -3 %).
  Coincide al 3 % en θ ≥ 3 y difiere un 9–13 % en θ ≤ 1: allí A difiere de la de `cmb_fondo.py` (fiducial p₀ en vez del
  de la tabla S, Ω_m = 0.304 frente a 0.315, y malla con los nodos exactos de DESI), y la cota lineal es sensible a eso en
  los θ pequeños, donde el modo común pesa más.
- Las muestras ajustan los datos: ⟨χ²_lin⟩ coincide con χ²_lin(E[m]) + tr(C⁻¹ Cov_post(m)) a ≤ 0.2 en las 42 celdas
  (tablas de abajo); n_eff crece de 3.0 (σ_X = 0) a 10–20 en las celdas grandes.

## Resultados

| θ | 5 % (N = 20000) | lineal analítica de EXACTO_SN (−1.92 / −1.35) | lineal con este A (−1.92 / −1.35) | exacta guiada −1.92 (16–84 %) | exacta −1.35 | ESS en la celda de la cota (% de N) | ln(Z_ex/Z_lin) en la cota | máx ln Z_ex: σ_X (valor) | MC directo −1.92, N = 20000 (ESS) |
|---|---|---|---|---|---|---|---|---|---|
| 0.3 | 0.292 | 0.567 / 0.505 | 0.492 / 0.427 | **0.500** (0.497–0.503) | 0.440 | 948 (19.0 %) | +0.04 | 0.2 (+1.14) | 0.496 (9) |
| 1 | 0.170 | 0.298 / 0.263 | 0.272 / 0.234 | **0.287** (0.282–0.293) | 0.244 | 187 (3.7 %) | +0.21 | 0.1 (+1.12) | 0.325 (5) |
| 3 | 0.102 | 0.178 / 0.154 | 0.173 / 0.152 | **0.177** (0.175–0.179) | 0.154 | 213 (4.3 %) | +0.14 | 0.07 (+0.95) | 0.173 (31) |
| 10 | 0.073 | 0.119 / 0.101 | 0.115 / 0.100 | **0.116** (0.115–0.118) | 0.101 | 229 (4.6 %) | +0.01 | 0.03 (+0.56) | 0.113 (28) |

(5 % y MC directo: `CONVERGENCIA_SN.md`, N = 20000; lineal analítica: `EXACTO_SN.md`.)

**Por celda:**

θ = 0.3:

| σ_X | ln Z_lin − ln Z_lin(0) | ln⟨w⟩ | ln Z_exacta (rel.) | ESS | peso máx | ⟨χ²_lin⟩ | esperado | n_eff | ln w (16 / 50 / 84 %) |
|---|---|---|---|---|---|---|---|---|---|
| 0 | +0.000 | +0.005 | +0.000 | 4999 | 0.0 % | 1409.3 | 1409.3 | 3.0 | -0.00 / +0.00 / +0.01 |
| 0.01 | +0.043 | +0.003 | +0.042 | 4998 | 0.0 % | 1409.1 | 1409.1 | 3.1 | -0.01 / +0.00 / +0.01 |
| 0.03 | +0.293 | -0.006 | +0.282 | 4971 | 0.0 % | 1407.8 | 1407.8 | 3.5 | -0.05 / -0.00 / +0.03 |
| 0.05 | +0.546 | -0.010 | +0.531 | 4879 | 0.1 % | 1406.5 | 1406.5 | 4.0 | -0.11 / -0.01 / +0.07 |
| 0.07 | +0.728 | -0.012 | +0.711 | 4728 | 0.1 % | 1405.3 | 1405.3 | 4.5 | -0.17 / -0.01 / +0.10 |
| 0.1 | +0.919 | -0.008 | +0.905 | 4212 | 0.3 % | 1403.5 | 1403.6 | 5.2 | -0.28 / -0.03 / +0.15 |
| 0.15 | +1.118 | -0.020 | +1.093 | 3962 | 0.3 % | 1401.1 | 1401.1 | 6.2 | -0.45 / -0.07 / +0.21 |
| 0.2 | +1.165 | -0.024 | +1.136 | 3414 | 0.2 % | 1399.5 | 1399.4 | 7.1 | -0.65 / -0.12 / +0.28 |
| 0.3 | +0.811 | -0.010 | +0.796 | 2290 | 0.4 % | 1397.8 | 1397.8 | 8.6 | -1.09 / -0.25 / +0.37 |
| 0.4 | +0.074 | +0.141 | +0.210 | 332 | 4.2 % | 1397.4 | 1397.3 | 9.9 | -1.62 / -0.40 / +0.46 |
| 0.5 | -0.823 | +0.046 | -0.782 | 948 | 0.9 % | 1397.1 | 1397.1 | 10.9 | -2.43 / -0.59 / +0.48 |
| 0.7 | -2.698 | +0.066 | -2.637 | 364 | 1.8 % | 1396.7 | 1396.7 | 12.7 | -5.01 / -1.32 / +0.33 |

θ = 1:

| σ_X | ln Z_lin − ln Z_lin(0) | ln⟨w⟩ | ln Z_exacta (rel.) | ESS | peso máx | ⟨χ²_lin⟩ | esperado | n_eff | ln w (16 / 50 / 84 %) |
|---|---|---|---|---|---|---|---|---|---|
| 0 | +0.000 | +0.005 | +0.000 | 4999 | 0.0 % | 1409.3 | 1409.3 | 3.0 | -0.00 / +0.00 / +0.01 |
| 0.01 | +0.129 | +0.001 | +0.126 | 4994 | 0.0 % | 1408.6 | 1408.7 | 3.2 | -0.02 / +0.00 / +0.02 |
| 0.03 | +0.603 | -0.015 | +0.584 | 4858 | 0.1 % | 1406.1 | 1406.3 | 4.1 | -0.12 / -0.01 / +0.07 |
| 0.05 | +0.873 | -0.014 | +0.853 | 4525 | 0.1 % | 1404.3 | 1404.3 | 4.9 | -0.25 / -0.03 / +0.14 |
| 0.07 | +1.031 | -0.012 | +1.015 | 3791 | 0.5 % | 1402.5 | 1402.5 | 5.7 | -0.40 / -0.05 / +0.19 |
| 0.1 | +1.139 | -0.017 | +1.118 | 3677 | 0.3 % | 1400.2 | 1400.3 | 6.8 | -0.57 / -0.09 / +0.27 |
| 0.15 | +0.950 | +0.025 | +0.971 | 1900 | 1.0 % | 1398.4 | 1398.3 | 8.2 | -0.97 / -0.20 / +0.36 |
| 0.2 | +0.400 | +0.074 | +0.469 | 972 | 1.5 % | 1397.5 | 1397.6 | 9.4 | -1.52 / -0.33 / +0.45 |
| 0.3 | -1.158 | +0.211 | -0.951 | 187 | 5.1 % | 1397.0 | 1397.1 | 11.3 | -2.86 / -0.77 / +0.41 |
| 0.4 | -2.818 | +0.041 | -2.782 | 126 | 7.0 % | 1396.6 | 1396.7 | 12.9 | -5.38 / -1.51 / +0.26 |

θ = 3:

| σ_X | ln Z_lin − ln Z_lin(0) | ln⟨w⟩ | ln Z_exacta (rel.) | ESS | peso máx | ⟨χ²_lin⟩ | esperado | n_eff | ln w (16 / 50 / 84 %) |
|---|---|---|---|---|---|---|---|---|---|
| 0 | +0.000 | +0.005 | +0.000 | 4999 | 0.0 % | 1409.3 | 1409.3 | 3.0 | -0.00 / +0.00 / +0.01 |
| 0.01 | +0.256 | -0.005 | +0.247 | 4988 | 0.0 % | 1408.0 | 1408.0 | 3.4 | -0.03 / -0.00 / +0.02 |
| 0.03 | +0.794 | -0.020 | +0.769 | 4584 | 0.1 % | 1404.9 | 1404.8 | 4.8 | -0.24 / -0.03 / +0.12 |
| 0.05 | +0.961 | -0.021 | +0.936 | 3570 | 0.5 % | 1402.4 | 1402.3 | 6.1 | -0.48 / -0.07 / +0.21 |
| 0.07 | +0.962 | -0.010 | +0.948 | 2737 | 0.6 % | 1400.3 | 1400.3 | 7.1 | -0.77 / -0.13 / +0.29 |
| 0.1 | +0.673 | -0.010 | +0.658 | 2094 | 0.6 % | 1398.6 | 1398.7 | 8.5 | -1.24 / -0.25 / +0.39 |
| 0.15 | -0.337 | +0.049 | -0.293 | 586 | 2.7 % | 1397.5 | 1397.6 | 10.4 | -2.34 / -0.55 / +0.46 |
| 0.2 | -1.604 | +0.148 | -1.461 | 213 | 4.1 % | 1397.0 | 1397.1 | 11.9 | -4.24 / -1.10 / +0.36 |
| 0.3 | -4.148 | +0.293 | -3.860 | 28 | 13.0 % | 1396.1 | 1396.1 | 14.2 | -10.76 / -3.05 / -0.28 |
| 0.4 | -6.431 | -0.451 | -6.886 | 27 | 16.1 % | 1395.0 | 1395.0 | 16.0 | -25.49 / -6.67 / -1.60 |

θ = 10:

| σ_X | ln Z_lin − ln Z_lin(0) | ln⟨w⟩ | ln Z_exacta (rel.) | ESS | peso máx | ⟨χ²_lin⟩ | esperado | n_eff | ln w (16 / 50 / 84 %) |
|---|---|---|---|---|---|---|---|---|---|
| 0 | +0.000 | +0.005 | +0.000 | 4998 | 0.0 % | 1409.3 | 1409.3 | 3.0 | -0.00 / +0.00 / +0.01 |
| 0.01 | +0.276 | -0.006 | +0.266 | 4979 | 0.0 % | 1407.8 | 1407.8 | 3.6 | -0.05 / -0.00 / +0.03 |
| 0.03 | +0.590 | -0.029 | +0.557 | 4225 | 0.1 % | 1404.4 | 1404.4 | 5.6 | -0.39 / -0.06 / +0.17 |
| 0.05 | +0.455 | +0.000 | +0.450 | 1355 | 1.4 % | 1401.9 | 1401.9 | 7.3 | -0.94 / -0.17 / +0.28 |
| 0.07 | +0.081 | +0.052 | +0.127 | 508 | 2.7 % | 1400.2 | 1400.2 | 8.8 | -1.65 / -0.34 / +0.37 |
| 0.1 | -0.751 | +0.011 | -0.745 | 229 | 5.3 % | 1398.6 | 1398.6 | 10.6 | -3.35 / -0.81 / +0.33 |
| 0.15 | -2.418 | +0.007 | -2.416 | 136 | 5.1 % | 1397.0 | 1396.9 | 12.9 | -8.54 / -2.15 / +0.05 |
| 0.2 | -4.114 | -0.288 | -4.407 | 42 | 9.7 % | 1395.6 | 1395.6 | 14.8 | -19.54 / -4.95 / -0.99 |
| 0.3 | -7.256 | -0.404 | -7.664 | 5 | 25.0 % | 1393.5 | 1393.4 | 17.7 | -59.16 / -15.52 / -5.04 |
| 0.4 | -10.081 | -2.103 | -12.189 | 2 | 70.2 % | 1391.8 | 1391.9 | 19.9 | -140.16 / -35.13 / -14.16 |

## La identidad y sus normalizaciones

Lo que se calcula en cada celda (σ_X, θ) es la evidencia del modelo exacto con el prior del camino y priores planos en δp:

Z_exacta(σ_X) = ∫ dδp ∫ dx p_ex(d | p₀ + δp, x) π_σ(x) π(δp),

con p_ex(d | ·) = (2π)^{−n/2} |C|^{−1/2} exp(−χ²_exacto/2), π_σ(x) = N(x; 0, σ_X² K_θ) (con su |2πσ_X²K|^{−1/2}) y π(δp) plana.
El posterior lineal es q(δp, x) = p_lin(d | δp, x) π_σ(x) π(δp) / Z_lineal, con p_lin la misma gaussiana (misma C, mismo
factor (2π)^{−n/2} |C|^{−1/2}) evaluada en el modelo lineal, y

Z_lineal(σ_X) = ∫ dδp ∫ dx p_lin(d | δp, x) π_σ(x) π(δp) = (2π)^{−n/2} |C|^{−1/2} |2πσ_X²K|^{−1/2} ∫∫ exp(−χ²_lin/2 − xᵀ(σ_X²K)⁻¹x/2)
             = const · |S|^{−1/2} |BᵀS⁻¹B|^{−1/2} exp(−½[rᵀS⁻¹r − bᵀF⁻¹b]),

con r = d − b₀, S = C + σ_X² A K Aᵀ, F = BᵀS⁻¹B, b = BᵀS⁻¹r, y "const" = (2π)^{(3−n)/2} × el volumen del prior plano de δp,
independiente de σ_X y de θ (la integral gaussiana en x cancela |2πσ_X²K|^{−1/2} contra |2πσ_X²K|^{1/2} y deja el |S|^{−1/2}
por el lema de Woodbury; la de δp deja |F|^{−1/2}). Esa expresión, relativa a σ_X = 0 (S = C), es el ln Z_lineal de la tabla.
Multiplicando y dividiendo por p_lin dentro de la integral de Z_exacta:

Z_exacta = ∫∫ [p_ex / p_lin] p_lin π_σ π dδp dx = Z_lineal ∫∫ [p_ex / p_lin] q dδp dx = Z_lineal · E_q[w],   w = p_ex / p_lin = exp(−(χ²_exacto − χ²_lineal)/2).

No falta ningún determinante ni normalización: el (2π)^{−n/2} |C|^{−1/2} de las dos verosimilitudes se cancela en el
cociente w, el prior π_σ(x) con su |2πσ_X²K|^{−1/2} y el prior plano de δp son los mismos en el numerador y en el denominador
(están dentro de q, no de w), y Z_lineal los lleva todos en su expresión analítica. Lo único que exige la identidad es que
las muestras sigan exactamente q, que es lo que hace Matheron (comprobado abajo), y que w use los dos χ² con la misma C,
los mismos datos y los mismos (h, ω_m, M), como se hace.

## Matheron en una celda

Celda (θ = 1, σ_X = 0.3), la de la cota exacta; 5000 muestras nuevas (semilla aparte), comparadas en 79 nodos (uno de cada
50, los siete de DESI y z*) con la media y la covarianza analíticas del posterior lineal de x con δp marginalizado,
μ = σ_X² K Aᵀ S⁻¹ (r − B δp̂) y Σ = σ_X² K − σ_X⁴ K Aᵀ S⁻¹ A K + σ_X⁴ K Aᵀ S⁻¹ B Σ_p Bᵀ S⁻¹ A K (`marginal_guiada_check.py` →
`marginal_guiada_check.json`):

| magnitud | error | esperado por Monte Carlo (N = 5000) |
|---|---|---|
| media, rms sobre nodos, en unidades de la desviación posterior | 0.016 σ | 0.014 σ |
| media, máximo | 0.033 σ | — |
| varianza (diagonal), rms del error relativo | 2.1 % | 2.0 % |
| varianza (diagonal), máximo del error relativo | 5.7 % | — |
| covarianza completa, error relativo de Frobenius | 3.1 % | — |
| correlaciones, máximo error absoluto | 0.037 | — |

La media posterior va de -0.108 a +0.233 y la desviación de 0.241 a 0.300: la media está dentro del
error de Monte Carlo (0.016σ frente a 0.014σ esperado) y la covarianza también (2.1 % frente a 2.0 %). El error relativo
máximo frente a la propia media, 41 % en nodos con |μ| > 0.01, no dice nada: ahí la media es mucho menor que la
desviación y el cociente lo dominan nodos con μ ≈ 0.

## Comprobación con la otra guía, y de dónde sale el 9–13 %

**La misma marginal exacta con el A de `cmb_fondo.py` / `exacto_sn.py` como guía** (malla uniforme, fiducial de la tabla S, BAO en el
nodo más próximo, escalado con el ajuste sin SN; el camino, en esa malla, interpolado a la de `exacto.py` para el χ² exacto;
`marginal_guiada_cf.py` → `marginal_guiada_cf.json`). La identidad vale para cualquier propuesta lineal, así que la cota exacta
debe salir igual; solo cambian ln Z_lin, la corrección y el ESS:

| θ | lineal con el A de cmb_fondo (ln Z_lin, máx. global) | exacta guiada por ese A, −1.92 (16–84 %) | ESS en la cota | corrección | exacta guiada por el A canónico (16–84 %) | ESS | corrección |
|---|---|---|---|---|---|---|---|
| 0.3 | 0.502 | **0.506** (0.501–0.510) | 339 (6.8 %) | +0.100 | **0.500** (0.497–0.503) | 948 | +0.041 |
| 1 | 0.278 | **0.279** (0.276–0.283) | 384 (7.7 %) | +0.065 | **0.287** (0.282–0.293) | 187 | +0.206 |

Las dos guías dan la misma cota dentro del error de bootstrap: la diferencia es 0.8σ (θ = 0.3) y 1.3σ (θ = 1, 2.8 %),
con la σ combinada de los dos bootstraps. La guía cambia el reparto entre ln Z_lin y la corrección (+0.21 con el A canónico
frente a +0.07 con el otro en θ = 1), no el producto. Ambas cotas usan la misma malla de σ_X; el error de interpolar en esa
malla frente a la raíz exacta es −1.9 % (θ = 1) y −0.1 % (θ = 0.3) en la lineal, comprobado con `brentq` en una malla fina.

**De dónde sale el 9–13 % con `EXACTO_SN.md`, corregido.** No es el A: con el máximo global de ln Z_lin y raíz exacta, el A
canónico da 0.277 / 0.492 / 0.174 / 0.118 (θ = 1 / 0.3 / 3 / 10) y el de `cmb_fondo.py` 0.283 / 0.502 / 0.177 / 0.119: a un 2 %. Lo que
difiere es el punto de referencia: `perfil_fondo.bound` mide la caída de 1.92 desde el máximo de la verosimilitud **sobre su
malla de σ_X, que llega solo a 0.05**, y en las combinaciones con tensión el máximo está más allá: en σ_X = 0.11 (θ = 1,
ln Z = +1.09 frente a +0.83 en σ_X ≤ 0.05) y 0.19 (θ = 0.3). Midiendo la caída desde ese máximo recortado, con el A de `cmb_fondo.py`,
salen 0.2982 / 0.5665 / 0.1778 / 0.1190: los valores de `EXACTO_SN.md` (0.2982 / 0.5665 / 0.1778 / 0.1190) a cuatro
decimales (`marginal_guiada_clip.json`). **Las cotas lineales marginalizadas de `EXACTO_SN.md` con el CMB comprimido, y las de
`CMB_FONDO.md` en las combinaciones cuyo σ_at_max salía en el borde 0.05, están referidas a un máximo recortado y son un
6–13 % más flojas de lo debido en θ ≤ 1**; en θ ≥ 3 el máximo está en σ_X ≤ 0.057 y el efecto es < 1 %. Sin CMB (SN solas,
BAO solo) el máximo está en σ_X = 0 y no hay sesgo. El párrafo anterior de esta nota atribuía la diferencia al fiducial y a
la construcción de las filas: el fiducial explica un 2 % (0.272 → 0.278 en θ = 1 con Ω_m = 0.3153 en la misma malla,
`marginal_guiada_check.json`) y el resto es el máximo recortado. Ninguna cifra calculada cambia; cuál de las dos lineales se
cita, y si se corrigen esas columnas, es decisión aparte.

## Lectura

- **La marginal exacta coincide con la lineal analítica.** Con el mismo A, la corrección ln(Z_ex/Z_lin) en las celdas de
  las cotas es +0.04, +0.21, +0.14, +0.01 (θ = 0.3, 1, 3, 10) y la cota exacta −1.92 queda a ×1.02, 1.06, 1.03, 1.01 de la lineal; frente a la
  lineal de `EXACTO_SN.md`, -12 %, -4 %, -0 %, -2 %. La no linealidad del modelo pesa poco hasta σ_X ≈ 0.2–0.3; a σ_X ≥ 0.3 (θ = 3, 10) los
  pesos se hunden (ln w mediana −3 a −35) porque la respuesta exacta de D_M y θ* a un camino grande es mucho menor que
  la lineal, pero eso está ya más allá de las cotas.
- **El muestreo guiado funciona donde el directo no:** ESS de 948, 187, 213, 229 en las celdas de las cotas (frente a 9, 5, 31, 28 con 20000
  caminos del prior en `CONVERGENCIA_SN.md`), bootstrap de ±1–2 %, y la cota es estable. El MC directo daba 0.496, 0.325, 0.173, 0.113: en
  θ = 1 sobreestimaba un 13 % (su ESS era 5); en los demás coincide al 3 %.
- **Marginal frente al 5 %:** ×1.71, 1.69, 1.73, 1.59: la verosimilitud marginal es sistemáticamente ×1.6–1.75 más floja que la regla del 5 %,
  con los dos métodos ya de acuerdo. Son dos preguntas distintas: el 5 % cuenta caminos que igualan a ΛCDM; la marginal
  promedia la verosimilitud, y el máximo de ln Z_ex en σ_X > 0 (0.2 (+1.14), 0.1 (+1.12), 0.07 (+0.95), 0.03 (+0.56)) desplaza el punto de referencia.
- Dónde el ESS cae por debajo del 1 % de N: solo en las celdas σ_X ≥ 0.3 de θ = 10 (ESS 2–5) y σ_X ≥ 0.3 de θ = 3 (27–28),
  todas más allá de la cota; no se ha forzado nada ahí y no entran en ninguna cota.

## Predicción: qué se cumplió

- **ESS mayor que el 10 % de N en las celdas de las cotas (θ = 1, 3, 10): no se cumplió.** Es el 3.7 %, 4.3 %, 4.6 % (187, 213, 229 de 5000);
  en θ = 0.3, el 19 %. Muy por encima del 1 %, pero por debajo del 10 %.
- **|ln(Z_ex/Z_lin)| < 0.5 en las celdas de las cotas lineales: cumplida** (+0.04, +0.21, +0.14, +0.01).
- **Cota exacta a menos del 15 % de la lineal analítica (0.298 / 0.178 / 0.119): cumplida** (-4 %, -0 %, -2 %).
- **Cota exacta entre 1.4 y 1.9 veces la del 5 % (0.170 / 0.102 / 0.073): cumplida** (×1.69, 1.73, 1.59).

**Conclusión.** La verosimilitud marginal exacta de la combinación BAO + θ* + ω_m + SN es la lineal analítica con una corrección
de +0.01 a +0.21 en ln Z en las celdas de las cotas: 0.287 / 0.177 / 0.116 (θ = 1 / 3 / 10, ±1–2 %) frente a 0.298 / 0.178 /
0.119. El MC directo del prior con 20000 caminos no servía (ESS 5–31); el guiado con 5000 muestras sí (ESS 187–229). Y las dos
cotas de cabecera quedan como estaban: la del 5 % (0.17 / 0.10 / 0.073) y, si se quiere la marginal, la lineal analítica,
que ahora está respaldada por el cálculo exacto.
