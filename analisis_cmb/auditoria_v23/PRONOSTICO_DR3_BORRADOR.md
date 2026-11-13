# PRONÓSTICO DR3 para el residuo de LRG2 — BORRADOR (13-11-2026)

**BORRADOR — pendiente de revisión independiente antes de registrar (tag).** No se crea tag ni se publica nada.

## Contexto

PREDICTIVO_BINS.md (commit 2d98ea0): LRG2 tiene p local 0.028 sin SN (d² = 7.13; con Pantheon+ 0.021) y la preferencia por σ_X > 0 depende de
él (JACKKNIFE_BINS.md, 1105cdf), pero el p global con 7 bins es 0.177. DESI DR3 (esperado 2027) contiene los datos de DR2 más ~2 años de
observación (Kim, Mota y Tamosiunas, arXiv:2607.28918 v3: pesos de tiempo 3/9 para DR2 y 5/9 para DR3). Aquí se diseña un pronóstico
registrable para el residuo de LRG2 en DR3 bajo dos hipótesis, usando explícitamente que los datos están anidados.

## Predicciones (de Claude, escritas antes de ejecutar, copiadas literales)

- P18. Con ρ² = σ3²/σ2² = 0.6 y falsa alarma del 5 %, la potencia del test para distinguir H0 de H1 en LRG2 está entre 0.5 y 0.75.
- P19. Con ρ² = 0.6, la mediana de d²_DR3 de LRG2 está entre 3 y 6 bajo H0 y entre 9 y 15 bajo H1.

(Se evalúan en la combinación sin SN, que es la del contexto; la de Pantheon+ se da al lado.)

## Supuestos del pronóstico (son SUPUESTOS, no propiedades garantizadas del análisis futuro)

- **Anidamiento.** Para los 7 bins de BAO, el ruido de DR3 es n3 = ρ²·n2 + e, con e independiente de n2 y Cov(e) = ρ²(1 − ρ²)·C. Así
  Cov(n3) = ρ²·C y Cov(n3, n2) = ρ²·C, con C la covarianza oficial de DR2 (diagonal por bloques entre bins).
- **Datos externos sin cambio** entre DR2 y DR3: θ*, ω_m y, en la combinación con SN, Pantheon+.
- **ρ² = 0.6** (caso principal); sensibilidad con ρ² = 0.5 y 0.72. El mismo ρ² en todos los bins.
- **H0 (ruido):** el universo es ΛCDM; LRG2 no tiene desvío propio.
- **H1 (persistente):** ΛCDM más un desvío fijo μ (2 componentes: D_M/r_d, D_H/r_d) solo en LRG2, desconocido, con prior plano no acotado.
  μ es físico: no se reduce entre DR2 y DR3; solo el ruido se transmite con el factor ρ².
- **Residuos:** r = datos_LRG2 − predicción ΛCDM ajustada sin LRG2, con el mismo procedimiento que PREDICTIVO_BINS (optimizador, priors,
  Fisher F = JᵀC⁻¹J). En DR3 se reajusta con los otros seis bins de DR3, que también cambian, usando la covarianza de DR3 (ρ²·C en BAO, la
  original en los externos).
- Modelo de fondo ΛCDM plano con (h, ω_m) [+ M], priors de Planck comprimido (100 θ*, ω_m) como dos filas gaussianas.

## Cantidades de DR2 que fijan el pronóstico (`pronostico_dr3_lrg2.json` → `DR2`)

Sin SN: p2 (ajuste sin LRG2) = (h, ω_m) = (0.67974, 0.14129); r2 = (-0.2096, -0.6353) en D_M/r_d, D_H/r_d; Q = J Σ2 Jᵀ = [[0.00325, 0.00184], [0.00184, 0.00105]];
C_LRG2 = [[0.03238, -0.02374], [-0.02374, 0.11147]]; V2 = [[0.03562, -0.02191], [-0.02191, 0.11252]]; d²_DR2 = 7.131 (p local 0.0283), idéntico al de PREDICTIVO_BINS (diferencia 0.0).
Con Pantheon+: r2 = (-0.2270, -0.6450); d²_DR2 = 7.701 (p local 0.0213), idéntico a PREDICTIVO_BINS. Los ajustadores con covarianza como argumento
reproducen los de `predictivo_bins.py` bit a bit (diferencia 0.0 en p0 y χ²).

## Fórmulas de referencia (aproximación gaussiana; NO son las distribuciones definitivas ni se usan para calcular Λ)

Con Q = J Σ2 Jᵀ (Σ2 la inversa de la Fisher del ajuste DR2 sin LRG2), V2 = C22 + Q, V3 = ρ²·C22 + Q, Cov(r3, r2) = ρ²·C22 + Q:

- **H0:** r3 | r2 ~ N(K r2, V3 − K V2 Kᵀ), K = (ρ²·C22 + Q)·V2⁻¹. Para ρ² = 0.6, sin SN: media (-0.1453, -0.3922), covarianza [[0.00822, -0.00545], [-0.00545, 0.02690]].
- **H1 sin reajuste**, con r2 = μ + n2 y r3 = μ + ρ²·n2 + e:
  - condicionada a μ conocido: r3 | r2, μ ~ N(ρ²·r2 + (1 − ρ²)·μ, ρ²(1 − ρ²)·C22); equivalentemente r3 − r2 | r2, μ ~ N((1 − ρ²)(μ − r2), ρ²(1 − ρ²)·C22);
  - **marginalizada** sobre μ con prior plano (μ | r2 ~ N(r2, C22)): r3 | r2 ~ N(r2, (1 − ρ²)·C22). Comprobación de la integración (numérica, en el
    json): media ρ²r2 + (1 − ρ²)r2 = r2 (diferencia 0.0); covarianza ρ²(1 − ρ²)C22 + (1 − ρ²)²C22 = (1 − ρ²)C22 (diferencia 0.0). Para ρ² = 0.6, sin SN:
    media (-0.2096, -0.6353), covarianza [[0.01295, -0.00950], [-0.00950, 0.04459]].
  - **La referencia que se compara con las simulaciones es la marginalizada**, porque en las simulaciones μ se saca de su posterior.
- **H1 con reajuste:** r3 − r2 ≈ (ρ² − 1)·n2 + e − J2·(p3 − p2), con p2 y p3 los parámetros ajustados sin LRG2 en DR2 y DR3. Esta es la expresión
  que corresponde al procedimiento, y la que las simulaciones realizan.

## Comprobaciones del posterior de H1 (antes de simular; Fisher conjunta de (p, μ) en el máximo conjunto (p2, μ̂ = r2), todas las filas)

- **sinSN.** (1) C_LRG2 invertible: autovalores 0.02579, 0.118, det 3.045e-03. (2) Bloque μ de la Fisher (C_LRG2⁻¹) invertible: autovalores 8.471, 38.77. (3) Fisher conjunta de (p, μ) definida positiva, inversa finita: autovalores 8.452, 33.53, 1.676e+05, 1.354e+07; |F Σ − I| máx 1.2e-14. (4) Autovalores de la covarianza posterior conjunta Σ = F⁻¹: 7.386e-08, 5.967e-06, 0.02982, 0.1183; su bloque p coincide con Σ2 (dif. máx 4.2e-21) y su bloque μ con V2 (dif. máx 1.4e-17). (5) Posterior de μ propio al integrar sobre p: la Fisher de p con los otros seis bins y los externos es definida positiva, autovalores 1.447e+05, 1.354e+07. Todo ok: True.
- **Pantheon+.** (1) C_LRG2 invertible: autovalores 0.02579, 0.118, det 3.045e-03. (2) Bloque μ de la Fisher (C_LRG2⁻¹) invertible: autovalores 8.471, 38.77. (3) Fisher conjunta de (p, μ) definida positiva, inversa finita: autovalores 8.453, 33.74, 6.18e+04, 1.981e+05, 1.353e+07; |F Σ − I| máx 3.9e-14. (4) Autovalores de la covarianza posterior conjunta Σ = F⁻¹: 7.39e-08, 5.049e-06, 1.618e-05, 0.02964, 0.1183; su bloque p coincide con Σ2 (dif. máx 4.2e-21) y su bloque μ con V2 (dif. máx 1.0e-17). (5) Posterior de μ propio al integrar sobre p: la Fisher de p con los otros seis bins y los externos es definida positiva, autovalores 6.013e+04, 1.768e+05, 1.353e+07. Todo ok: True.

Ninguna comprobación falló; no se añadió ningún límite al prior de μ.

## Generación de DR3 (ecuaciones del procedimiento simulado)

- **H0.** p_true ~ N(p_full, Σ_full): aproximación gaussiana del posterior ΛCDM con todos los datos DR2 (p_full el ajuste completo, Σ_full la
  inversa de su Fisher). Ruido DR2 implícito en los 7 bins: n2 = datos_DR2 − modelo(p_true). Datos DR3 = modelo(p_true) + ρ²·n2 + e,
  e ~ N(0, ρ²(1 − ρ²)·C). Externos sin cambio.
- **H1.** (p_true, μ) del posterior conjunto con μ libre y prior plano no acotado: p_true ~ N(p2, Σ2) y μ | p_true ~ N(datos_LRG2 −
  modelo_LRG2(p_true), C_LRG2) (al marginalizar μ, las filas de LRG2 no informan sobre p). Con E₂μ el vector de los 7 bins con μ en las dos
  componentes de LRG2 y cero en los demás: ruido inferido n2 = datos_DR2 − modelo(p_true) − E₂μ; **datos DR3 = modelo(p_true) + E₂μ + ρ²·n2 + e**.
  Así μ permanece fijo entre DR2 y DR3 y solo el ruido se transmite con ρ². No se usa ρ²·(n2 + E₂μ), que reduciría el desvío persistente.
- **En cada universo DR3:** reajuste ΛCDM sin LRG2 con la covarianza de DR3 (ρ²·C en BAO, la original en los externos) → p3; r3 = residuo de LRG2
  frente a la predicción de ese reajuste; Σ3 = inversa de la Fisher de ese reajuste (con ρ²·C); V3,k = ρ²·C_LRG2 + J3 Σ3 J3ᵀ; d²_DR3 = r3ᵀ V3,k⁻¹ r3 con
  la V3,k de ese mismo universo.
- **Versión linealizada (la de los 100000 universos):** modelo, reajuste, J y Fisher linealizados alrededor del ajuste DR2 sin LRG2 (J evaluada en p2;
  con ello Σ3 y V3,k son las mismas en todos los universos de un caso). 100000 universos por hipótesis, combinación y ρ², en dos mitades
  independientes con semillas distintas: A (densidades) y B (gaussianidad, umbrales, potencia). Semillas en el json: [20261114, i_comb, i_ρ², i_hip, mitad]
  con mitad 0 = A, 1 = B, 2 = validación.

## Validación del linealizado frente al cálculo completo (12 casos × 20 universos; generación con el modelo exacto, reajuste con el optimizador original, J y Fisher recalculadas en el nuevo mínimo)

| combinación | ρ² | hipótesis | máx \|Δd²\| | universo | máx \|Δr3\| | ≤ 0.1 |
|---|---|---|---|---|---|---|
| sinSN | 0.6 | H0 | 0.0101 | 4 | 6.2e-04 | sí |
| sinSN | 0.6 | H1 | 0.0196 | 19 | 2.5e-04 | sí |
| sinSN | 0.5 | H0 | 0.0160 | 5 | 7.1e-04 | sí |
| sinSN | 0.5 | H1 | 0.0336 | 15 | 4.7e-04 | sí |
| sinSN | 0.72 | H0 | 0.0074 | 18 | 3.5e-04 | sí |
| sinSN | 0.72 | H1 | 0.0164 | 9 | 1.5e-04 | sí |
| Pantheon+ | 0.6 | H0 | 0.0205 | 7 | 5.6e-04 | sí |
| Pantheon+ | 0.6 | H1 | 0.0253 | 2 | 3.2e-04 | sí |
| Pantheon+ | 0.5 | H0 | 0.0232 | 3 | 7.1e-04 | sí |
| Pantheon+ | 0.5 | H1 | 0.0441 | 0 | 7.2e-04 | sí |
| Pantheon+ | 0.72 | H0 | 0.0153 | 18 | 3.7e-04 | sí |
| Pantheon+ | 0.72 | H1 | 0.0185 | 7 | 3.1e-04 | sí |

Mayor diferencia de todas: 0.0441 (Pantheon+, ρ² = 0.5, H1, universo 0). Ninguna supera 0.1. Esto respalda el linealizado **en esta
comprobación** (estos residuos y covarianzas), no en general.

## Comprobación de gaussianidad (mitad B, con media y covarianza estimadas en A)

| combinación | ρ² | hipótesis | Mahalanobis pct 5 / 50 / 95 (χ²₂: 0.103 / 1.386 / 5.991) | fracción > 5.991 (ref. 0.05) | gaussiana |
|---|---|---|---|---|---|
| sinSN | 0.6 | H0 | 0.104 / 1.380 / 5.973 | 0.0495 | sí |
| sinSN | 0.6 | H1 | 0.104 / 1.397 / 6.046 | 0.0516 | sí |
| sinSN | 0.5 | H0 | 0.103 / 1.397 / 6.052 | 0.0514 | sí |
| sinSN | 0.5 | H1 | 0.102 / 1.380 / 5.932 | 0.0485 | sí |
| sinSN | 0.72 | H0 | 0.101 / 1.390 / 6.020 | 0.0509 | sí |
| sinSN | 0.72 | H1 | 0.098 / 1.371 / 5.959 | 0.0493 | sí |
| Pantheon+ | 0.6 | H0 | 0.104 / 1.385 / 5.995 | 0.0501 | sí |
| Pantheon+ | 0.6 | H1 | 0.102 / 1.383 / 5.930 | 0.0484 | sí |
| Pantheon+ | 0.5 | H0 | 0.106 / 1.396 / 5.963 | 0.0493 | sí |
| Pantheon+ | 0.5 | H1 | 0.099 / 1.385 / 5.965 | 0.0496 | sí |
| Pantheon+ | 0.72 | H0 | 0.102 / 1.381 / 6.008 | 0.0504 | sí |
| Pantheon+ | 0.72 | H1 | 0.103 / 1.372 / 5.953 | 0.0492 | sí |

En los 12 casos la distancia de Mahalanobis sigue χ²₂ (desviaciones de los percentiles ≤ 4 %, fracción sobre el cuantil 95 % entre 0.048 y 0.052):
las densidades son gaussianas bidimensionales y no hace falta la estimación por núcleos. (Es lo esperado: en la versión linealizada r3 es lineal en
variables gaussianas.)

## Comparación con las fórmulas de referencia (media y covarianza simuladas en A frente a la referencia; H1: la marginalizada)

| combinación | ρ² | hipótesis | media simulada − referencia (D_M / D_H) | en unidades de σ | varianza diagonal, simulada/referencia − 1 |
|---|---|---|---|---|---|
| sinSN | 0.6 | H0 | +0.0139 / +0.0069 | +0.151 / +0.042 | +3.7 % / +1.1 % |
| sinSN | 0.6 | H1 | +0.0003 / +0.0023 | +0.003 / +0.011 | +7.6 % / +0.1 % |
| sinSN | 0.5 | H0 | +0.0185 / +0.0100 | +0.196 / +0.059 | +1.4 % / -0.1 % |
| sinSN | 0.5 | H1 | -0.0003 / -0.0012 | -0.003 / -0.005 | +8.7 % / +1.2 % |
| sinSN | 0.72 | H0 | +0.0095 / +0.0050 | +0.112 / +0.033 | +5.0 % / +0.5 % |
| sinSN | 0.72 | H1 | -0.0001 / -0.0004 | -0.001 / -0.002 | +7.5 % / +0.9 % |
| Pantheon+ | 0.6 | H0 | +0.0138 / +0.0072 | +0.150 / +0.044 | +3.4 % / -0.2 % |
| Pantheon+ | 0.6 | H1 | -0.0016 / +0.0019 | -0.013 / +0.009 | +6.6 % / +2.3 % |
| Pantheon+ | 0.5 | H0 | +0.0173 / +0.0111 | +0.184 / +0.066 | +0.7 % / -0.4 % |
| Pantheon+ | 0.5 | H1 | -0.0005 / +0.0000 | -0.003 / +0.000 | +6.8 % / +1.3 % |
| Pantheon+ | 0.72 | H0 | +0.0097 / +0.0048 | +0.116 / +0.032 | +4.2 % / -0.2 % |
| Pantheon+ | 0.72 | H1 | +0.0002 / -0.0011 | +0.002 / -0.006 | +7.1 % / +0.7 % |

Lectura: bajo H1 la media simulada coincide con r2 (≤ 0.013 σ) y la varianza de D_H con (1 − ρ²)C22 (≤ 2 %); la varianza de D_M es un 7–9 % mayor
que la referencia sin reajuste, que es el término −J2·(p3 − p2) de la fórmula con reajuste (Q es relativamente mayor en D_M). Bajo H0 la media
simulada difiere de K r2 en 0.11–0.20 σ en D_M (0.01–0.02 en D/r_d) y ≤ 0.07 σ en D_H, y la varianza en ≤ 5 %: la referencia condiciona en r2 con
Q de DR2 y prior plano, mientras la simulación saca p_true del posterior completo y reajusta con Σ3; ninguna de las dos se usa para Λ.

## Resultados: d²_DR3 y Λ

| combinación | ρ² | hipótesis | d²_DR3 pct 5 / 50 / 95 (mitad B) | media |
|---|---|---|---|---|
| sinSN | 0.6 | H0 | 1.40 / 4.76 / 10.18 | 5.14 |
| sinSN | 0.6 | H1 | 4.92 / 12.29 / 23.13 | 12.93 |
| sinSN | 0.5 | H0 | 0.91 / 4.15 / 10.04 | 4.64 |
| sinSN | 0.5 | H1 | 5.21 / 14.93 / 29.57 | 15.86 |
| sinSN | 0.72 | H0 | 2.23 / 5.45 / 10.14 | 5.73 |
| sinSN | 0.72 | H1 | 4.82 / 10.14 / 17.53 | 10.54 |
| Pantheon+ | 0.6 | H0 | 1.61 / 5.10 / 10.68 | 5.49 |
| Pantheon+ | 0.6 | H1 | 5.54 / 13.27 / 24.32 | 13.89 |
| Pantheon+ | 0.5 | H0 | 1.04 / 4.47 / 10.39 | 4.94 |
| Pantheon+ | 0.5 | H1 | 5.93 / 16.04 / 31.31 | 16.99 |
| Pantheon+ | 0.72 | H0 | 2.52 / 5.91 / 10.73 | 6.17 |
| Pantheon+ | 0.72 | H1 | 5.38 / 10.94 / 18.44 | 11.29 |

Λ(r3) = ln p(r3 | datos DR2, H1) − ln p(r3 | datos DR2, H0), densidades gaussianas estimadas en A, evaluada en B:

| combinación | ρ² | hipótesis | Λ pct 5 / 50 / 95 (mitad B) | P(Λ > 3) / P(Λ < −3) |
|---|---|---|---|---|
| sinSN | 0.6 | H0 | -3.28 / -1.65 / +1.01 | 0.005 / 0.102 |
| sinSN | 0.6 | H1 | -1.57 / +2.05 / +7.41 | 0.369 / 0.004 |
| sinSN | 0.5 | H0 | -3.79 / -2.17 / +0.76 | 0.005 / 0.249 |
| sinSN | 0.5 | H1 | -1.65 / +3.21 / +10.54 | 0.523 / 0.008 |
| sinSN | 0.72 | H0 | -2.77 / -1.17 / +1.19 | 0.004 / 0.025 |
| sinSN | 0.72 | H1 | -1.48 / +1.19 / +4.92 | 0.196 / 0.001 |
| Pantheon+ | 0.6 | H0 | -3.49 / -1.79 / +0.99 | 0.005 / 0.146 |
| Pantheon+ | 0.6 | H1 | -1.57 / +2.28 / +7.83 | 0.404 / 0.005 |
| Pantheon+ | 0.5 | H0 | -4.01 / -2.30 / +0.67 | 0.005 / 0.296 |
| Pantheon+ | 0.5 | H1 | -1.57 / +3.51 / +11.19 | 0.553 / 0.009 |
| Pantheon+ | 0.72 | H0 | -2.92 / -1.22 / +1.20 | 0.005 / 0.040 |
| Pantheon+ | 0.72 | H1 | -1.48 / +1.30 / +5.06 | 0.213 / 0.001 |

## Regla de decisión (calibrada con la mitad B)

| combinación | ρ² | c0 (pct 95 de Λ, H0) | c1 (pct 5 de Λ, H1) | c1 < c0 | potencia P(Λ > c0 \| H1) | P(Λ < c1 \| H0) | no concluyente H0 / H1 |
|---|---|---|---|---|---|---|---|
| sinSN | 0.6 | +1.011 | -1.575 | sí | 0.655 | 0.523 | 0.427 / 0.295 |
| sinSN | 0.5 | +0.765 | -1.648 | sí | 0.767 | 0.640 | 0.310 / 0.183 |
| sinSN | 0.72 | +1.187 | -1.484 | sí | 0.500 | 0.394 | 0.556 / 0.450 |
| Pantheon+ | 0.6 | +0.988 | -1.570 | sí | 0.684 | 0.562 | 0.388 / 0.266 |
| Pantheon+ | 0.5 | +0.671 | -1.567 | sí | 0.794 | 0.683 | 0.267 / 0.156 |
| Pantheon+ | 0.72 | +1.200 | -1.483 | sí | 0.521 | 0.417 | 0.533 / 0.429 |

En los 6 casos c1 < c0, así que la regla de tres resultados es aplicable: **Λ > c0 favorece H1** (falsa alarma del 5 % bajo H0); **Λ < c1 favorece H0**
(error del 5 % bajo H1); en medio, **no concluyente**. |Λ| > 3 se usa como indicador adicional de evidencia fuerte, sin sustituir a c0 y c1: bajo H0,
P(Λ > 3) ≈ 0.005 en todos los casos; bajo H1, P(Λ > 3) = 0.37 (ρ² = 0.6, sin SN). Con ρ² = 0.72 la potencia queda en 0.500: DR3 separa las
dos hipótesis solo a medias si el ruido de DR2 pesa tanto.

## Veredicto de P18–P19

- **P18: sí.** ρ² = 0.6, falsa alarma 5 %: potencia 0.655 sin SN (con Pantheon+ 0.684), dentro de [0.5, 0.75].
- **P19: sí.** ρ² = 0.6, sin SN: mediana de d²_DR3 4.76 bajo H0 (en [3, 6]) y 12.29 bajo H1 (en [9, 15]). Con Pantheon+: 5.10 y 13.27.

Nota sobre d² bajo H0: su mediana (4.8) es mayor que la de una χ²₂ (1.39) porque DR3 hereda el 60 % del ruido de DR2 (ρ²·n2) y la V3,k del reajuste
no incluye esa correlación con DR2: d²_DR3 no es un p-valor de DR3 por sí solo; lo que se compara es Λ, que sí condiciona en DR2.

## Qué rompería el modelo y debe declararse si ocurre (y cómo se comprobaría al salir DR3)

1. **Reprocesado de DR3 que cambie la medida de DR2** (nueva reducción, nuevos templates o ajustes de BAO aplicados también a los datos de DR2): el
   anidamiento n3 = ρ²n2 + e deja de valer. Comprobación: comparar el vector de 13 números y la covarianza de DR2 publicados con DR3 (si DESI los
   republica) con los de DR2 usados aquí (sha256 en `data/PROCEDENCIA_DESI_DR2.md`); si cambian, recalcular r2 y declarar el cambio.
2. **Cambio de z_eff o de bins** (otro agrupamiento, LRG2 con otro intervalo de z, bins combinados): r3 no es comparable con r2. Comprobación: leer
   z_eff y la definición de muestras de DR3; si LRG2 no existe con z_eff = 0.706 ± 0.01, el pronóstico no se evalúa.
3. **Cambio de cosmología fiducial de DESI** (la que fija las distancias de referencia y las compresiones): los números D/r_d cambian de forma
   sistemática. Comprobación: comparar la fiducial declarada en el paper de DR3 con la de DR2; si cambia, recalcular con la nueva o declarar no evaluable.
4. **ρ² distinto por bin** (tiempos de observación distintos por muestra): el reajuste con los otros seis bins usa otra covarianza. Comprobación:
   ρ²_k = σ²_DR3,k / σ²_DR2,k por bin a partir de las covarianzas publicadas; repetir las simulaciones con los ρ²_k medidos antes de evaluar Λ.
5. **Sistemáticos nuevos** (términos de sistemáticos añadidos a la covarianza o correcciones nuevas): Cov(n3) ≠ ρ²C. Comprobación: comparar la
   covarianza 13 × 13 de DR3 con ρ²·C_DR2 bin a bin; si la diferencia relativa de las diagonales supera lo que explica un ρ² por bin, declarar.
6. **Cambio en los datos externos** (otra compresión de Planck, otra muestra de SN o nueva versión de Pantheon+): el reajuste cambia. Comprobación:
   el pronóstico se evalúa con exactamente los externos de este json (θ*, ω_m, Pantheon+ 1580 SN con sus sha256 en `datos/FUENTES.json`); si se
   usan otros, es otro pronóstico.

En cualquiera de los seis casos la evaluación se declara "modelo roto" antes de mirar Λ, y la ruptura se documenta con los números que la muestran.

## Cautelas

- Las densidades y umbrales dependen de los supuestos de arriba (anidamiento exacto con un ρ² común, externos fijos, gaussianidad, ΛCDM con
  Planck comprimido) y de la linealización, validada solo en estas condiciones.
- H1 es un desvío constante en un solo bin: otras hipótesis (desvío que cambia con la muestra, física que afecte a varios bins) no están
  representadas, y Λ solo compara estas dos.
- Este es un borrador: nada se registra hasta la revisión independiente.

## Ficheros

`pronostico_dr3_lrg2.py`, `pronostico_dr3_lrg2.json` (DR2, comprobaciones del posterior, por combinación y ρ²: Σ3 y V3,k linealizadas, referencias,
semillas, medias y covarianzas simuladas, gaussianidad, comparación con referencias, percentiles de d² y Λ, regla de decisión, validación con las
20 comparaciones por caso), `pronostico_dr3_lrg2_run.log`. Nada existente se modifica.
