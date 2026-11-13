# PRONÓSTICO DR3 para el residuo de LRG2 — REGISTRO, versión 3 (13-11-2026)

**PENDIENTE DE COMPROBACIÓN DEL REVISOR — todavía sin tag.**

Parte de `PRONOSTICO_DR3_REGISTRO.md` (v2, commit 3cdb79b) con las correcciones pedidas por la segunda revisión `REVISION2_PRONOSTICO_DR3.md`
(commit c7df2d8, veredicto "no registrable hasta corregir"). Cada cambio está marcado **(cambio por la revisión 2)**; los de la primera revisión siguen
marcados **(cambio por la revisión)**. No cambian P18, P19, H0, H1, el anidamiento como idea ni la regla de decisión (c0 = percentil 95 de Λ | H0, c1 =
percentil 5 de Λ | H1, tres salidas). El borrador, v2 y sus json quedan como historia. Código de esta versión: `pronostico_dr3_lrg2_v3.py` (sha256
`4625081d0b9ad3363c6f1d670b6c03c31e97ae3b7486a88a08e3291a9cd462d2`, también en `pronostico_dr3_lrg2_v3.sha256`); resultados `pronostico_dr3_lrg2_v3.json` (sha256 `35936e8bfe48396d2a402ecb7961acdeb69c1c7180cb6daa387ca203bdc5b0ec`).
**(cambio por la revisión 2)** El campo `git_head_at_run` de los json es, inevitablemente, el del commit anterior a esta versión (c7df2d8); lo que
fija el registro es el tag que se ponga al commit que contiene estos ficheros, no ese campo.

## Predicciones (de Claude, escritas antes de ejecutar; copiadas literales)

- P18. Con ρ² = σ3²/σ2² = 0.6 y falsa alarma del 5 %, la potencia del test para distinguir H0 de H1 en LRG2 está entre 0.5 y 0.75.
- P19. Con ρ² = 0.6, la mediana de d²_DR3 de LRG2 está entre 3 y 6 bajo H0 y entre 9 y 15 bajo H1.
- **(nueva, antes de ejecutar v3)** P23. Con los umbrales de alta precisión, la potencia del caso principal (ρ² = 0.6, sin SN) queda en 0.656 ± 0.010.
  Se evalúa solo con: sin SN, ρ² = 0.6, umbrales de alta precisión (c0 y c1 nuevos) y la misma regla de decisión; no se compara una potencia con
  umbral antiguo con otra con umbral nuevo.

## Contexto, hipótesis y qué no es este pronóstico (sin cambio)

PREDICTIVO_BINS.md (commit 2d98ea0): LRG2 tiene p local 0.028 sin SN (d² = 7.13; con Pantheon+ 0.021) y la preferencia por σ_X > 0 depende de él
(JACKKNIFE_BINS.md, 1105cdf), pero el p global con 7 bins es 0.177. DESI DR3 (esperado 2027) contiene los datos de DR2 más ~2 años de observación.
**H0 (ruido):** el universo es ΛCDM; LRG2 no tiene desvío propio. **H1 (persistente):** ΛCDM más un desvío fijo μ (2 componentes: D_M/r_d, D_H/r_d) solo en
LRG2, desconocido, con prior plano no acotado; μ es físico, no se reduce entre DR2 y DR3; solo el ruido se transmite. **(cambio por la revisión)** H1 es un
desvío fijo en un bin, no el modelo OU: un resultado a favor de H1 no es evidencia de energía oscura estocástica ni de σ_X > 0; uno a favor de H0 no excluye
σ_X > 0 dentro de las cotas actuales. **Este pronóstico no contiene ninguna prueba del modelo OU.**

**(cambio por la revisión 2) Salida principal única.** La salida del pronóstico es la de la combinación **sin SN** (BAO + θ* + ω_m; la de P18, P19 y P23).
La combinación con Pantheon+ se calcula y se informa al lado como **secundaria**, y no puede cambiar el veredicto principal. El script devuelve
`result_principal` (sin SN) y `result_secondary_Pantheon+` por separado.

## Supuestos del pronóstico (son SUPUESTOS, no propiedades garantizadas del análisis futuro)

- **(cambio por la revisión 2) Anidamiento matricial.** Para los 7 bins de BAO, con C2 la covarianza oficial de DR2 y C3 la que publique DR3 (ambas
  diagonales por bloques entre bins): n3 = R·n2 + e, con R = C3·C2⁻¹ (bloque por bin), e independiente de n2 y Cov(e) = C3 − C3·C2⁻¹·C3. Deducción:
  Cov(n3, n2) = R·Cov(n2) = C3·C2⁻¹·C2 = C3; Cov(n3) = R·C2·Rᵀ + Cov(e) = C3·C2⁻¹·C2·C2⁻¹·C3 + C3 − C3·C2⁻¹·C3 = C3·C2⁻¹·C3 + C3 − C3·C2⁻¹·C3 = C3.
  Es la combinación por inversa de la varianza de DR2 con la parte nueva independiente (Cov(θ̂_sub, θ̂_total) = Var(θ̂_total)). Cuando C3 = ρ²·C2 se reduce
  al caso de v2: R = ρ²·I y Cov(e) = ρ²·C2 − ρ⁴·C2 = ρ²(1 − ρ²)·C2. Para el registro de ahora se usa C3 = ρ²·C2 con ρ² = 0.6 (principal), 0.5 y 0.72
  (sensibilidad) y un barrido 0.40–0.80; en la evaluación se usa la C3 publicada, bloque a bloque, sin reducirla a un escalar.
- **(cambio por la revisión 2) Comprobaciones obligatorias en cada bin**, con la tolerancia de autovalores: en una matriz que debe ser semidefinida
  positiva, un autovalor λ < −1e-10·λ_max cuenta como negativo de verdad (fallo); si −1e-10·λ_max ≤ λ < 0, es redondeo, se registra como tal y se
  trunca a 0 al factorizar. (i) C2 y C3 simétricas (|C − Cᵀ| máx registrado) y definidas positivas (Cholesky sin fallo); (ii) C2 − C3 semidefinida positiva
  (DR3 igual o más preciso que DR2 en todas las direcciones); (iii) Cov(e) = C3 − C3·C2⁻¹·C3 semidefinida positiva. Si alguna falla de verdad en algún
  bin: **modelo roto**. Los autovalores se registran en el json (`nesting`).
- **Datos externos sin cambio** entre DR2 y DR3: θ*, ω_m y, en la combinación con SN, Pantheon+.
- **(cambio por la revisión)** ρ² = 3/5 sale de los pesos de tiempo de observación de Kim, Mota y Tamosiunas (arXiv:2607.28918 v3, §4.4: 1/9, 3/9 y 5/9
  para DR1, DR2 y DR3) y **supone σ² ∝ 1/tiempo**, que no es exacto para BAO. En la evaluación eso ya no se supone: entra la C3 publicada.
- **(cambio por la revisión) Prior plano en p explícito:** las posteriores gaussianas N(p_full, Σ_full) y N(p2, Σ2) suponen prior plano en (h, ω_m[, M])
  e ignoran las cotas del optimizador.
- **Residuos:** r = datos_LRG2 − predicción ΛCDM ajustada sin LRG2, mismo procedimiento que PREDICTIVO_BINS (optimizador, priors, Fisher F = JᵀC⁻¹J).
  En DR3 se reajusta con los otros seis bins homólogos de DR3 (regla e) usando C3 en BAO y la covarianza original en los externos.
- Modelo de fondo ΛCDM plano con (h, ω_m) [+ M], priors de Planck comprimido (100 θ*, ω_m) como dos filas gaussianas.
- **(cambio por la revisión 2) Álgebra numérica:** en todo el script, factorizaciones de Cholesky y resolución de sistemas lineales (`cho_factor`,
  `cho_solve`, `solve_triangular`); no se calculan inversas explícitas. Las muestras de p ~ N(p̂, F⁻¹) se obtienen como p̂ + L⁻ᵀz con F = LLᵀ.

## Cantidades de DR2 que fijan el pronóstico (reproducidas en v3 a ≤ 4e-18 respecto a v2)

Sin SN: p2 = (h, ω_m) = (0.67974, 0.14129); r2 = (-0.2096, -0.6353) en D_M/r_d, D_H/r_d; Q = J Σ2 Jᵀ = [[0.00325, 0.00184], [0.00184, 0.00105]];
C_LRG2 = [[0.03238, -0.02374], [-0.02374, 0.11147]]; V2 = [[0.03562, -0.02191], [-0.02191, 0.11252]]; d²_DR2 = 7.131 (p local 0.0283). Con Pantheon+: r2 = (-0.2270, -0.6450); d²_DR2 = 7.701.

## (cambio por la revisión 2) Comprobaciones de matrices del registro (C3 = 0.6·C2; las de 0.5 y 0.72 son proporcionales y están en el json)

| bin | autovalores C2 | autovalores C3 = 0.6·C2 | autovalores C2 − C3 | autovalores Cov(e) = C3 − C3C2⁻¹C3 | asimetría máx | Cholesky C2, C3 | negativos reales / de redondeo | ρ²_eff = tr C3 / tr C2 |
|---|---|---|---|---|---|---|---|---|
| BGS | 0.00579 | 0.00347 | 0.00232 | 0.00139 | 0e+00 | sí | 0 / 0 | 0.600 |
| LRG1 | 0.02179, 0.19049 | 0.01307, 0.11429 | 0.00872, 0.07619 | 0.00523, 0.04572 | 0e+00 | sí | 0 / 0 | 0.600 |
| LRG2 | 0.02579, 0.11805 | 0.01548, 0.07083 | 0.01032, 0.04722 | 0.00619, 0.02833 | 0e+00 | sí | 0 / 0 | 0.600 |
| LRG3+ELG1 | 0.01994, 0.04665 | 0.01197, 0.02799 | 0.00798, 0.01866 | 0.00479, 0.01120 | 0e+00 | sí | 0 / 0 | 0.600 |
| ELG2 | 0.03792, 0.11784 | 0.02275, 0.07070 | 0.01517, 0.04714 | 0.00910, 0.02828 | 0e+00 | sí | 0 / 0 | 0.600 |
| QSO | 0.17495, 0.67641 | 0.10497, 0.40585 | 0.06998, 0.27056 | 0.04199, 0.16234 | 0e+00 | sí | 0 / 0 | 0.600 |
| Lya | 0.00826, 0.28464 | 0.00496, 0.17078 | 0.00330, 0.11385 | 0.00198, 0.06831 | 0e+00 | sí | 0 / 0 | 0.600 |

Ningún autovalor negativo, ni real ni de redondeo; Cholesky de C2 y C3 sin fallo en los 7 bins; ρ²_eff = 0.600 en todos. Las de ρ² = 0.5 y 0.72 pasan igual.

## Fórmulas de referencia (sin cambio; no se usan para Λ)

Con Q = J Σ2 Jᵀ, V2 = C22 + Q, V3 = ρ²·C22 + Q y Cov(r3, r2) **sin reajuste** = ρ²·C22 + Q: H0 sin reajuste r3 | r2 ~ N(K r2, V3 − K V2 Kᵀ),
K = (ρ²·C22 + Q)·V2⁻¹ (ρ² = 0.6, sin SN: media (-0.1453, -0.3922)); H1 sin reajuste, marginalizada en μ con prior plano: r3 | r2 ~ N(r2, (1 − ρ²)·C22)
(la igualdad ρ²r2 + (1 − ρ²)r2 = r2 es una identidad algebraica, no una comprobación); H1 con reajuste: r3 − r2 ≈ (ρ² − 1)·n2 + e − J2·(p3 − p2). La prueba
real es la distribución exacta de abajo.

## Densidades exactas del modelo linealizado (ahora con la C3 matricial) y comprobación frente a v2

En el procedimiento linealizado r3 es afín en variables gaussianas estándar z = (z_p, z_e[, z_μ]): r3 = a + M z, con a = r3(0) y columnas de M = r3(e_k) − a
(no linealidad residual ≤ 4e-15). **Bajo cada hipótesis r3 | datos DR2 ~ N(a, MMᵀ) exactamente** y **Λ se calcula con estas densidades.** **(cambio por la
revisión 2)** La generación usa R = C3·C2⁻¹ y un factor de Cov(e) por descomposición espectral (autovalores de redondeo a 0); con C3 = ρ²·C2 las densidades
**coinciden con las de v2**: máx |Δa| ≤ 9e-16 y máx |ΔΣ| ≤ 5e-16 en los 6 casos (tabla de abajo), y V3,k a ≤ 1.4e-17. Las simulaciones de 100000 universos de
v2 (medias a ≤ 0.0072 σ con A+B, varianzas a ≤ 1.64 %) siguen valiendo como comprobación de esas densidades **(cambio por la revisión 2: redondeos exactos)**.

**(cambio por la revisión) Gaussianidad.** En el linealizado es cierta por construcción; lo que respalda el r3 real es la validación con el modelo completo
(máx |Δr3| 7e-4 en D/r_d = 0.0084 σ de la componente más estrecha, √0.00702 = 0.084).

## Generación de DR3 (ecuaciones del procedimiento)

- **H0.** p_true ~ N(p_full, Σ_full) (prior plano en p). n2 = datos_DR2 − modelo(p_true) en los 13 números. **Datos DR3 = modelo(p_true) + R·n2 + e**, e ~ N(0, C3 − C3C2⁻¹C3).
  Externos sin cambio.
- **H1.** p_true ~ N(p2, Σ2) y μ | p_true ~ N(datos_LRG2 − modelo_LRG2(p_true), C_LRG2). Con E₂μ el vector con μ en LRG2 y cero en los demás: n2 = datos_DR2 −
  modelo(p_true) − E₂μ; **datos DR3 = modelo(p_true) + E₂μ + R·n2 + e**. μ queda fijo; solo el ruido se transmite con R.
- **Reajuste:** ΛCDM sin LRG2 con C3 en BAO y la covarianza original en los externos → p3; r3 = residuo de LRG2; V3,k = C3,LRG2 + J3 Σ3 J3ᵀ; d²_DR3 = r3ᵀ V3,k⁻¹ r3.
  **(cambio por la revisión)** En el linealizado J, Σ3 y V3,k son fijas por caso (J en p2); en la evaluación real el reajuste es completo (optimizador original, J y
  Fisher en el nuevo mínimo).
- Comprobaciones del posterior de H1 (autovalores) y validación del linealizado (12 casos × 20 universos, máx |Δd²| 0.0441 en Pantheon+, ρ² = 0.5, H1, universo 0): sin cambio
  respecto al borrador y a v2 (`pronostico_dr3_lrg2.json`).

| combinación | ρ² | hipótesis | máx \|Δd²\| | universo | máx \|Δr3\| |
|---|---|---|---|---|---|
| sinSN | 0.6 | H0 | 0.0101 | 4 | 6.2e-04 |
| sinSN | 0.6 | H1 | 0.0196 | 19 | 2.5e-04 |
| sinSN | 0.5 | H0 | 0.0160 | 5 | 7.1e-04 |
| sinSN | 0.5 | H1 | 0.0336 | 15 | 4.7e-04 |
| sinSN | 0.72 | H0 | 0.0074 | 18 | 3.5e-04 |
| sinSN | 0.72 | H1 | 0.0164 | 9 | 1.5e-04 |
| Pantheon+ | 0.6 | H0 | 0.0205 | 7 | 5.6e-04 |
| Pantheon+ | 0.6 | H1 | 0.0253 | 2 | 3.2e-04 |
| Pantheon+ | 0.5 | H0 | 0.0232 | 3 | 7.1e-04 |
| Pantheon+ | 0.5 | H1 | 0.0441 | 0 | 7.2e-04 |
| Pantheon+ | 0.72 | H0 | 0.0153 | 18 | 3.7e-04 |
| Pantheon+ | 0.72 | H1 | 0.0185 | 7 | 3.1e-04 |

## (cambio por la revisión 2) Umbrales de alta precisión

Como Λ depende solo de r3 y las dos densidades son gaussianas exactas, c0 y c1 se calculan con **10⁷ muestras de N(a, Σ)** bajo H0 y bajo H1 (semillas
[20261114, i_comb, índice de caso, i_hip, 7]: índice 0, 1, 2 para ρ² = 0.6, 0.5, 0.72; 100 + k en el barrido; 500 en la evaluación), y el error Monte Carlo es la
desviación típica de 10 bloques de 10⁶ dividida por √10. Objetivo ±0.002: conseguido (≤ 0.0024 en c0 y c1, ≤ 0.0004 en potencia). El mismo método se usa en
el caso principal, en la sensibilidad y en el barrido de ρ². **De dónde venía la diferencia 1.081 / 1.024 de v2:** v2 calibraba cada caso con una mitad B
distinta de 50000 universos (semillas distintas para el caso principal y para el punto 0.60 del barrido), y el ruido MC de ±0.02 en c0 producía dos valores
para el mismo ρ²; con 10⁷ muestras exactas el valor es único: c0 = 1.0274 ± 0.0013 (principal) y 1.0285 ± 0.0009 (punto 0.60 del barrido, otra semilla).

| combinación | ρ² | máx \|Δa\| / \|ΔΣ\| frente a v2 | v2 (mitad B de 50000): c0 / c1 / potencia ± MC | v3 (10⁷ muestras exactas): c0 / c1 / potencia ± MC | desviación máx. de v2 en sus σ_MC |
|---|---|---|---|---|---|
| sinSN | 0.6 | 3e-17 / 3e-16 | +1.024 ± 0.017 / -1.604 ± 0.015 / 0.6557 ± 0.0031 | +1.0274 ± 0.0013 / -1.5944 ± 0.0012 / 0.6575 ± 0.0002 | 0.9 |
| sinSN | 0.5 | 8e-16 / 5e-16 | +0.777 ± 0.018 / -1.638 ± 0.019 / 0.7669 ± 0.0027 | +0.7559 ± 0.0009 / -1.6115 ± 0.0020 / 0.7675 ± 0.0001 | 1.4 |
| sinSN | 0.72 | 3e-17 / 2e-16 | +1.184 ± 0.016 / -1.476 ± 0.015 / 0.4998 ± 0.0048 | +1.1971 ± 0.0012 / -1.4707 ± 0.0008 / 0.4997 ± 0.0003 | 2.2 |
| Pantheon+ | 0.6 | 6e-17 / 3e-16 | +0.993 ± 0.021 / -1.579 ± 0.020 / 0.6850 ± 0.0037 | +0.9999 ± 0.0014 / -1.5861 ± 0.0006 / 0.6801 ± 0.0002 | 1.3 |
| Pantheon+ | 0.5 | 9e-16 / 5e-16 | +0.667 ± 0.020 / -1.560 ± 0.024 / 0.7939 ± 0.0027 | +0.7061 ± 0.0014 / -1.5726 ± 0.0014 / 0.7878 ± 0.0002 | 2.3 |
| Pantheon+ | 0.72 | 8e-17 / 2e-16 | +1.195 ± 0.015 / -1.478 ± 0.013 / 0.5211 ± 0.0033 | +1.1968 ± 0.0006 / -1.4830 ± 0.0010 / 0.5208 ± 0.0002 | 0.4 |

Los umbrales de v2 están a ≤ 2.3 σ_MC (de v2) de los de alta precisión: los errores declarados en v2 eran correctos; los de v3 son 10–15 veces menores.

## Resultados (densidades exactas, 10⁷ muestras; ± error MC)

| combinación | ρ² | hipótesis | d²_DR3 pct 5 / 50 / 95 | mediana ± MC | E[d²] exacto |
|---|---|---|---|---|---|
| sinSN | 0.6 | H0 | 1.42 / 4.76 / 10.18 | 4.764 ± 0.001 | 5.146 |
| sinSN | 0.6 | H1 | 4.94 / 12.32 / 23.12 | 12.320 ± 0.001 | 12.956 |
| sinSN | 0.5 | H0 | 0.91 / 4.16 / 10.00 | 4.163 ± 0.001 | 4.639 |
| sinSN | 0.5 | H1 | 5.26 / 14.91 / 29.66 | 14.907 ± 0.003 | 15.852 |
| sinSN | 0.72 | H0 | 2.25 / 5.48 / 10.17 | 5.480 ± 0.001 | 5.748 |
| sinSN | 0.72 | H1 | 4.83 / 10.16 / 17.50 | 10.164 ± 0.001 | 10.535 |
| Pantheon+ | 0.6 | H0 | 1.61 / 5.12 / 10.70 | 5.115 ± 0.001 | 5.500 |
| Pantheon+ | 0.6 | H1 | 5.52 / 13.24 / 24.35 | 13.241 ± 0.002 | 13.871 |
| Pantheon+ | 0.5 | H0 | 1.04 / 4.46 / 10.47 | 4.460 ± 0.001 | 4.938 |
| Pantheon+ | 0.5 | H1 | 5.91 / 16.00 / 31.18 | 15.998 ± 0.003 | 16.936 |
| Pantheon+ | 0.72 | H0 | 2.52 / 5.90 / 10.74 | 5.899 ± 0.001 | 6.168 |
| Pantheon+ | 0.72 | H1 | 5.37 / 10.94 / 18.51 | 10.943 ± 0.001 | 11.308 |

| combinación | ρ² | hipótesis | Λ pct 5 / 50 / 95 | P(Λ > 3) / P(Λ < −3) ± MC |
|---|---|---|---|---|
| sinSN | 0.6 | H0 | -3.354 / -1.682 / +1.027 | 0.00513 ± 0.00002 / 0.11558 ± 0.00012 |
| sinSN | 0.6 | H1 | -1.594 / +2.095 / +7.494 | 0.37613 ± 0.00008 / 0.00416 ± 0.00002 |
| sinSN | 0.5 | H0 | -3.786 / -2.160 / +0.756 | 0.00490 ± 0.00002 / 0.24464 ± 0.00019 |
| sinSN | 0.5 | H1 | -1.612 / +3.210 / +10.584 | 0.52279 ± 0.00017 / 0.00785 ± 0.00003 |
| sinSN | 0.72 | H0 | -2.761 / -1.145 / +1.197 | 0.00391 ± 0.00002 / 0.02400 ± 0.00007 |
| sinSN | 0.72 | H1 | -1.471 / +1.196 / +4.863 | 0.19348 ± 0.00013 / 0.00097 ± 0.00001 |
| Pantheon+ | 0.6 | H0 | -3.543 / -1.790 / +1.000 | 0.00527 ± 0.00001 / 0.15353 ± 0.00012 |
| Pantheon+ | 0.6 | H1 | -1.586 / +2.271 / +7.826 | 0.40193 ± 0.00020 / 0.00518 ± 0.00002 |
| Pantheon+ | 0.5 | H0 | -4.006 / -2.296 / +0.706 | 0.00489 ± 0.00003 / 0.29434 ± 0.00013 |
| Pantheon+ | 0.5 | H1 | -1.573 / +3.470 / +11.057 | 0.54958 ± 0.00015 / 0.00879 ± 0.00003 |
| Pantheon+ | 0.72 | H0 | -2.910 / -1.221 / +1.197 | 0.00423 ± 0.00003 / 0.03939 ± 0.00008 |
| Pantheon+ | 0.72 | H1 | -1.483 / +1.300 / +5.082 | 0.21374 ± 0.00006 / 0.00150 ± 0.00001 |

## Regla de decisión (sin cambio de regla; umbrales de alta precisión)

| combinación | ρ² | c0 | c1 | c1 < c0 | potencia P(Λ > c0 \| H1) | P(Λ < c1 \| H0) | no concluyente H0 / H1 |
|---|---|---|---|---|---|---|---|
| sinSN | 0.6 | 1.0274 ± 0.0013 | -1.5944 ± 0.0012 | sí | 0.6575 ± 0.0002 | 0.5258 ± 0.0003 | 0.4242 ± 0.0003 / 0.2925 ± 0.0002 |
| sinSN | 0.5 | 0.7559 ± 0.0009 | -1.6115 ± 0.0020 | sí | 0.7675 ± 0.0001 | 0.6459 ± 0.0006 | 0.3041 ± 0.0006 / 0.1825 ± 0.0001 |
| sinSN | 0.72 | 1.1971 ± 0.0012 | -1.4707 ± 0.0008 | sí | 0.4997 ± 0.0003 | 0.3900 ± 0.0002 | 0.5600 ± 0.0002 / 0.4503 ± 0.0003 |
| Pantheon+ | 0.6 | 0.9999 ± 0.0014 | -1.5861 ± 0.0006 | sí | 0.6801 ± 0.0002 | 0.5576 ± 0.0002 | 0.3924 ± 0.0002 / 0.2699 ± 0.0002 |
| Pantheon+ | 0.5 | 0.7061 ± 0.0014 | -1.5726 ± 0.0014 | sí | 0.7878 ± 0.0002 | 0.6806 ± 0.0004 | 0.2694 ± 0.0004 / 0.1622 ± 0.0002 |
| Pantheon+ | 0.72 | 1.1968 ± 0.0006 | -1.4830 ± 0.0010 | sí | 0.5208 ± 0.0002 | 0.4149 ± 0.0003 | 0.5351 ± 0.0003 / 0.4292 ± 0.0002 |

**Regla.** En los 6 casos c1 < c0, así que la regla de tres resultados es aplicable **bajo los supuestos declarados** (anidamiento matricial con la C3
publicada y sus comprobaciones, externos fijos, linealización validada, ΛCDM con Planck comprimido): **Λ > c0 favorece H1** (falsa alarma del 5 % bajo H0, bajo
esos supuestos); **Λ < c1 favorece H0** (error del 5 % bajo H1, bajo esos supuestos); en medio, **no concluyente**. **(cambio por la revisión)** |Λ| > 3 es solo un
indicador adicional de **evidencia fuerte de una hipótesis frente a la otra, dentro de este par y bajo estos supuestos**; es asimétrico: P(Λ < −3 | H0) =
0.1156 ± 0.0001 frente a P(Λ > 3 | H1) = 0.3761 ± 0.0001 (ρ² = 0.6, sin SN). Con ρ² = 0.72 la potencia es 0.4997 ± 0.0003: DR3 separa las
dos hipótesis solo a medias si el ruido de DR2 pesa tanto. **(cambio por la revisión 2)** Con la C3 publicada, calibrar con el anidamiento matricial hace que la
falsa alarma del 5 % sea exacta dentro del modelo; se eliminan las tolerancias del 20 % entre componentes y de 0.1 en la correlación de v2, que dejaban la
falsa alarma real entre 2 % y 8 % (medido por la revisión 2).

## Potencia frente a ρ² común, 0.40–0.80 (10⁷ muestras exactas por hipótesis y punto; sin SN principal, Pantheon+ secundaria)

| ρ² común | sin SN (principal): c0 / c1 | potencia | P(Λ < c1 \| H0) | no concl. H0 / H1 | mediana d² H0 / H1 | Pantheon+ (secundaria): c0 / c1 | potencia | mediana d² H0 / H1 |
|---|---|---|---|---|---|---|---|---|
| 0.40 | +0.3615 ± 0.0009 / -1.5143 ± 0.0019 | 0.8542 ± 0.0001 | 0.7688 | 0.181 / 0.096 | 3.56 / 18.78 | +0.2808 / -1.4200 | 0.8707 ± 0.0001 | 3.80 / 20.12 |
| 0.42 | +0.4536 ± 0.0014 / -1.5470 ± 0.0019 | 0.8386 ± 0.0002 | 0.7441 | 0.206 / 0.111 | 3.68 / 17.86 | +0.3777 / -1.4687 | 0.8561 ± 0.0001 | 3.93 / 19.14 |
| 0.44 | +0.5369 ± 0.0011 / -1.5679 ± 0.0013 | 0.8224 ± 0.0001 | 0.7201 | 0.230 / 0.128 | 3.80 / 17.03 | +0.4680 / -1.5009 | 0.8407 ± 0.0002 | 4.07 / 18.25 |
| 0.46 | +0.6143 ± 0.0013 / -1.5869 ± 0.0019 | 0.8051 ± 0.0002 | 0.6957 | 0.254 / 0.145 | 3.92 / 16.26 | +0.5518 / -1.5290 | 0.8240 ± 0.0002 | 4.20 / 17.42 |
| 0.48 | +0.6879 ± 0.0017 / -1.6028 ± 0.0013 | 0.7867 ± 0.0002 | 0.6706 | 0.279 / 0.163 | 4.04 / 15.56 | +0.6312 / -1.5529 | 0.8064 ± 0.0001 | 4.33 / 16.68 |
| 0.50 | +0.7551 ± 0.0023 / -1.6112 ± 0.0015 | 0.7675 ± 0.0003 | 0.6459 | 0.304 / 0.182 | 4.16 / 14.91 | +0.7060 / -1.5706 | 0.7879 ± 0.0002 | 4.46 / 16.00 |
| 0.52 | +0.8175 ± 0.0016 / -1.6144 ± 0.0015 | 0.7473 ± 0.0002 | 0.6219 | 0.328 / 0.203 | 4.28 / 14.31 | +0.7736 / -1.5798 | 0.7686 ± 0.0003 | 4.59 / 15.36 |
| 0.54 | +0.8782 ± 0.0014 / -1.6157 ± 0.0007 | 0.7264 ± 0.0002 | 0.5973 | 0.353 / 0.224 | 4.40 / 13.76 | +0.8364 / -1.5880 | 0.7480 ± 0.0002 | 4.72 / 14.77 |
| 0.56 | +0.9340 ± 0.0015 / -1.6108 ± 0.0012 | 0.7043 ± 0.0003 | 0.5736 | 0.376 / 0.246 | 4.52 / 13.25 | +0.8955 / -1.5875 | 0.7266 ± 0.0001 | 4.86 / 14.23 |
| 0.58 | +0.9827 ± 0.0007 / -1.6025 ± 0.0009 | 0.6814 ± 0.0001 | 0.5500 | 0.400 / 0.269 | 4.64 / 12.77 | +0.9500 / -1.5919 | 0.7038 ± 0.0002 | 4.99 / 13.72 |
| 0.60 | +1.0285 ± 0.0009 / -1.5930 ± 0.0009 | 0.6577 ± 0.0002 | 0.5264 | 0.424 / 0.292 | 4.76 / 12.33 | +1.0003 / -1.5847 | 0.6802 ± 0.0002 | 5.12 / 13.24 |
| 0.62 | +1.0667 ± 0.0014 / -1.5794 ± 0.0012 | 0.6335 ± 0.0003 | 0.5034 | 0.447 / 0.316 | 4.88 / 11.91 | +1.0458 / -1.5758 | 0.6555 ± 0.0002 | 5.25 / 12.80 |
| 0.64 | +1.1043 ± 0.0010 / -1.5647 ± 0.0006 | 0.6078 ± 0.0002 | 0.4797 | 0.470 / 0.342 | 5.00 / 11.52 | +1.0842 / -1.5623 | 0.6307 ± 0.0003 | 5.38 / 12.38 |
| 0.66 | +1.1328 ± 0.0008 / -1.5451 ± 0.0010 | 0.5822 ± 0.0003 | 0.4571 | 0.493 / 0.368 | 5.12 / 11.15 | +1.1184 / -1.5476 | 0.6046 ± 0.0002 | 5.51 / 11.99 |
| 0.68 | +1.1585 ± 0.0011 / -1.5235 ± 0.0006 | 0.5558 ± 0.0002 | 0.4344 | 0.516 / 0.394 | 5.24 / 10.80 | +1.1495 / -1.5314 | 0.5774 ± 0.0003 | 5.64 / 11.62 |
| 0.70 | +1.1780 ± 0.0013 / -1.4982 ± 0.0013 | 0.5285 ± 0.0002 | 0.4122 | 0.538 / 0.422 | 5.36 / 10.47 | +1.1755 / -1.5093 | 0.5495 ± 0.0003 | 5.77 / 11.27 |
| 0.72 | +1.1964 ± 0.0009 / -1.4695 ± 0.0009 | 0.5000 ± 0.0002 | 0.3905 | 0.559 / 0.450 | 5.48 / 10.17 | +1.1929 / -1.4868 | 0.5209 ± 0.0003 | 5.90 / 10.94 |
| 0.74 | +1.2065 ± 0.0009 / -1.4392 ± 0.0005 | 0.4718 ± 0.0002 | 0.3684 | 0.582 / 0.478 | 5.60 / 9.87 | +1.2078 / -1.4573 | 0.4917 ± 0.0002 | 6.03 / 10.63 |
| 0.76 | +1.2110 ± 0.0016 / -1.4063 ± 0.0010 | 0.4427 ± 0.0004 | 0.3466 | 0.603 / 0.507 | 5.72 / 9.60 | +1.2160 / -1.4258 | 0.4621 ± 0.0003 | 6.16 / 10.33 |
| 0.78 | +1.2099 ± 0.0012 / -1.3680 ± 0.0006 | 0.4134 ± 0.0003 | 0.3252 | 0.625 / 0.537 | 5.83 / 9.33 | +1.2181 / -1.3910 | 0.4317 ± 0.0003 | 6.29 / 10.05 |
| 0.80 | +1.2023 ± 0.0008 / -1.3262 ± 0.0005 | 0.3838 ± 0.0003 | 0.3037 | 0.646 / 0.566 | 5.95 / 9.09 | +1.2117 / -1.3512 | 0.4013 ± 0.0002 | 6.42 / 9.78 |

## Veredicto de P18, P19 y P23 (caso principal: sin SN, ρ² = 0.6, umbrales de alta precisión)

- **P18: sí.** Potencia 0.6575 ± 0.0002, dentro de [0.5, 0.75]. (Pantheon+, secundaria: 0.6801 ± 0.0002.)
- **P19: sí.** Mediana de d²_DR3 4.764 ± 0.001 bajo H0 (en [3, 6]) y 12.320 ± 0.001 bajo H1 (en [9, 15]).
- **P23: sí.** Con c0 = 1.0274 ± 0.0013 y c1 = -1.5944 ± 0.0012, la potencia es 0.6575 ± 0.0002, dentro de 0.656 ± 0.010.

Nota sobre d² bajo H0: su mediana (4.8) es mayor que la de una χ²₂ (1.39) porque DR3 hereda el 60 % del ruido de DR2 y V3,k no incluye esa correlación:
d²_DR3 no es un p-valor de DR3 por sí solo; lo que decide es Λ, que condiciona en DR2.

## Modo "evaluación DR3" (escrito y probado con datos falsos; todavía sin usar)

- **(cambio por la revisión 2) Identificación de los bins en DR3 (regla f, igual en el código y en el texto):** cada fila de DR2 (tipo D_M/D_H/D_V, trazador,
  z_eff de la etiqueta) se empareja con la fila de DR3 del mismo tipo y trazador cuyo z_eff está en z_eff(DR2) ± 0.01, y el emparejamiento debe ser único;
  LRG2 = trazador LRG con z_eff en 0.706 ± 0.01. Si falta un homólogo o hay ambigüedad: **no evaluable**. Las filas de DR3 sin homólogo (bins nuevos) se
  ignoran (regla e). Si algún z_eff difiere (dentro de ±0.01) se registra el desplazamiento en la calibración y se declara antes de evaluar; la predicción
  se hace en el z_eff de DR2.
- **(cambio por la revisión 2) Regla b con tolerancia:** coeficientes de correlación entre bins distintos de DR3 con |r| ≥ 1e-6 → **modelo roto**.
- **(cambio por la revisión 2) Reglas c y d:** se aplican las comprobaciones de matrices (simetría, Cholesky, C2 − C3 y Cov(e) semidefinidas positivas con la
  tolerancia de autovalores) y el rango ρ²_LRG2,eff = tr(C3,LRG2)/tr(C2,LRG2) ∈ [0.40, 0.80]; fuera → **modelo roto**. No hay tolerancias de 20 % ni de 0.1:
  las cubre el anidamiento matricial.
- `calibrar COV.json TAG` lee **solo** el fichero de covarianza, aplica las reglas, construye las densidades exactas con esa C3, calcula c0 y c1 con 10⁷
  muestras (semillas 500) y escribe `..._calibracion_TAG.json` con los sha256 de la covarianza y del script más su propio sha256 en `.sha256`.
  `evaluar COV.json VALORES.json TAG` comprueba esos sha256, lee los valores, reajusta ΛCDM sin LRG2 con el optimizador original y C3 (cálculo completo),
  calcula r3, d²_DR3 y Λ, da `result_principal` (sin SN) y la secundaria, y se niega si ya existe el resultado de ese TAG.
- **(cambio por la revisión 2) Qué son esos controles.** El script comprueba su propio sha256 contra `pronostico_dr3_lrg2_v3.sha256` y el de la calibración
  contra su `.sha256`: son **controles de coherencia** (protegen de errores y accidentes), **no el compromiso**. Como mostró la revisión 2, se pueden saltar
  reescribiendo ficheros o usando otra etiqueta; lo que impide decidir después de ver los datos es el **orden público de los commits** del procedimiento de
  abajo.

**Prueba con datos falsos** (marcados "DATOS FALSOS" en su cabecera; `pronostico_dr3_lrg2_v3_falso_run.log`):

1. C3 = 0.6·C2, valores = DR2: estado "ok", ρ²_eff = 0.6 en los 7 bins, correlación entre bins 0.0; c0 / c1 (sin SN) = +1.0260 / -1.5921 (± 0.0012),
   frente a +1.0274 / -1.5944 del registro con otra semilla: **lo mismo que v2 con umbrales de alta precisión** dentro del error MC. Evaluación: r3 = (-0.182, -0.620),
   d² = 10.38, Λ = +1.127 → principal "favorece H1" (secundaria Pantheon+: "favorece H1"), idénticos a los de v2 con los mismos datos falsos. La segunda evaluación
   con el mismo TAG se niega.
2. C3 = 0.6·C2 con Var(D_H de LRG2) = 2 × la de DR2 (una dirección menos precisa que DR2): C2 − C3 tiene autovalor -0.1122 (negativo real) y Cov(e)
   -0.2635; ρ²_LRG2,eff = 1.685 → **"modelo roto"**, sin calibrar ni calcular Λ.
3. C3 = 0.6·C2 con correlación 1e-3 entre LRG1 y LRG2: |r| máx 1e-03 ≥ 1e-6 → **"modelo roto"**.
4. (extra) correlación 1e-9 entre bins: por debajo de la tolerancia → "ok", misma calibración que 1.

## (cambio por la revisión 2) PROCEDIMIENTO PÚBLICO — este es el compromiso

Al salir DESI DR3:

1. Descargar **solo** el fichero de covarianza de las medidas BAO de DR3 y guardarlo con su sha256 y su URL/fecha.
2. Calibrar con el script congelado (`pronostico_dr3_lrg2_v3.py calibrar`, commit de este registro) y hacer **commit y push públicos** de la calibración
   y de su `.sha256`.
3. Solo después, descargar los valores y evaluar **una sola vez** (`evaluar`).
4. Publicar el resultado con una de las cinco salidas: **favorece H1 / favorece H0 / no concluyente / modelo roto / no evaluado**, con estos umbrales y sin
   añadir umbrales después.

El orden lo prueban las fechas de los commits públicos. Si DESI publica valores y covarianza en un único fichero, se declara en el registro antes de abrirlo
y se mantiene el orden de lectura dentro del script (`calibrar` solo lee la clave `cov`). Fecha límite: si las medidas BAO de DR3 no son públicas antes del
31-12-2028, el pronóstico se declara **no evaluado** (no fallido).

## Rupturas (sin cambio salvo las reglas marcadas arriba)

1. Reprocesado de DR3 que cambie la medida de DR2: comparar con los sha256 de `data/PROCEDENCIA_DESI_DR2.md`; si cambian, recalcular r2 y declarar.
2. Cambio de bins o de z_eff: regla f (homólogo único por tipo, trazador y z_eff ± 0.01); si no, no evaluable.
3. Cambio de cosmología fiducial de DESI: comparar con la de DR2; recalcular con la nueva o declarar no evaluable.
4. Covarianza de DR3 no anidada (reglas b, c, d): modelo roto.
5. Sistemáticos nuevos: los cubre la comprobación C2 − C3 ⪰ 0 y Cov(e) ⪰ 0; si fallan, modelo roto.
6. Cambio en los datos externos: con otros externos es otro pronóstico.
- a. Formato: evaluable con (D_M/r_d, D_H/r_d) y covarianza; con (α_iso, α_AP) y fiducial declarada, conversión determinista (α_⊥ = α_iso·α_AP^{−1/3} → D_M,
  α_∥ = α_iso·α_AP^{2/3} → D_H, covarianza por el jacobiano) registrada antes de evaluar; solo D_V o LRG2 no separado: no evaluable.
- e. "Los otros seis bins" son los homólogos de DR2; los nuevos se ignoran. f. Erratas propias descubiertas después: se declaran y se recalcula con registro antes de mirar Λ.

## Cautelas

- Densidades y umbrales dependen de los supuestos y de la linealización, validada solo en estas condiciones.
- H1 es un desvío constante en un solo bin; Λ solo compara estas dos hipótesis. Ni H1 ni Λ prueban el modelo OU ni σ_X > 0.
- Pendiente de comprobación del revisor; sin tag.

## Ficheros nuevos de esta versión

`pronostico_dr3_lrg2_v3.py`, `pronostico_dr3_lrg2_v3.sha256`, `pronostico_dr3_lrg2_v3.json`, `pronostico_dr3_lrg2_v3_run.log`; prueba con datos falsos:
`pronostico_dr3_lrg2_v3_falso*_dr3_cov.json` (4), `pronostico_dr3_lrg2_v3_falso_dr3_valores.json`, `pronostico_dr3_lrg2_v3_calibracion_*.json` (+ `.sha256`),
`pronostico_dr3_lrg2_v3_resultado_*.json`, `pronostico_dr3_lrg2_v3_falso_run.log`. Nada existente se modifica.
