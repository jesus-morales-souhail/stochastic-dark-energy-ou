# PRONÓSTICO DR3 para el residuo de LRG2 — REGISTRO (versión 2, 13-11-2026)

**PENDIENTE DE COMPROBACIÓN DEL REVISOR — todavía sin tag.**

Parte de `PRONOSTICO_DR3_BORRADOR.md` (commit c742989) con los cambios pedidos por la revisión independiente `REVISION_PRONOSTICO_DR3.md`
(commit 99f879d, veredicto "registrable con estos cambios"). Cada cambio está marcado **(cambio por la revisión)**. No cambian P18, P19, ningún
supuesto ni ningún umbral. El borrador, su script y su json quedan como historia; el código de esta versión es `pronostico_dr3_lrg2_v2.py`
(sha256 `add9489055326e38b8f58f4d5f66687aa5058c09ad77b825e5d105c20feded1a`) y sus resultados `pronostico_dr3_lrg2_v2.json` (sha256 `4c831e0e69e4cb6eb8801d10c0342f73abb55d1e7e280a78b674c7c0a2c19099`).

## Contexto

PREDICTIVO_BINS.md (commit 2d98ea0): LRG2 tiene p local 0.028 sin SN (d² = 7.13; con Pantheon+ 0.021) y la preferencia por σ_X > 0 depende de
él (JACKKNIFE_BINS.md, 1105cdf), pero el p global con 7 bins es 0.177. DESI DR3 (esperado 2027) contiene los datos de DR2 más ~2 años de
observación. Aquí se registra un pronóstico para el residuo de LRG2 en DR3 bajo dos hipótesis, usando explícitamente que los datos están anidados.

**(cambio por la revisión) Qué es y qué no es este pronóstico.** H1 es un desvío fijo en un bin, no el modelo OU. Un resultado a favor de H1 no es
evidencia de energía oscura estocástica ni de σ_X > 0; uno a favor de H0 no excluye σ_X > 0 dentro de las cotas actuales (JACKKNIFE_BINS,
PROPUESTA_V3). **Este pronóstico no contiene ninguna prueba del modelo OU.** Λ compara dos hipótesis concretas sobre dos números de un bin.

## Predicciones (de Claude, escritas antes de ejecutar el borrador, copiadas literales; sin cambio)

- P18. Con ρ² = σ3²/σ2² = 0.6 y falsa alarma del 5 %, la potencia del test para distinguir H0 de H1 en LRG2 está entre 0.5 y 0.75.
- P19. Con ρ² = 0.6, la mediana de d²_DR3 de LRG2 está entre 3 y 6 bajo H0 y entre 9 y 15 bajo H1.

(Se evalúan en la combinación sin SN; la de Pantheon+ se da al lado.)

## Supuestos del pronóstico (son SUPUESTOS, no propiedades garantizadas del análisis futuro)

- **Anidamiento.** Para los 7 bins de BAO, el ruido de DR3 es n3 = ρ²·n2 + e, con e independiente de n2 y Cov(e) = ρ²(1 − ρ²)·C. Así
  Cov(n3) = ρ²·C y Cov(n3, n2) = ρ²·C, con C la covarianza oficial de DR2 (diagonal por bloques entre bins).
- **Datos externos sin cambio** entre DR2 y DR3: θ*, ω_m y, en la combinación con SN, Pantheon+.
- **ρ² = 0.6** (caso principal); sensibilidad con ρ² = 0.5 y 0.72, y barrido continuo 0.40–0.80. **(cambio por la revisión) Supuesto
  declarado:** ρ² = 3/5 sale de los pesos de tiempo de observación de Kim, Mota y Tamosiunas (arXiv:2607.28918 v3, §4.4: 1/9, 3/9 y 5/9 para DR1,
  DR2 y DR3) y **supone σ² ∝ 1/tiempo**, que no es exacto para BAO (área frente a profundidad, varianza cósmica, cambios de pipeline). El script
  admite ρ² distinto por bin (vector de 7) y el modo de evaluación lo mide de la covarianza publicada (regla c).
- **H0 (ruido):** el universo es ΛCDM; LRG2 no tiene desvío propio.
- **H1 (persistente):** ΛCDM más un desvío fijo μ (2 componentes: D_M/r_d, D_H/r_d) solo en LRG2, desconocido, con prior plano no acotado.
  μ es físico: no se reduce entre DR2 y DR3; solo el ruido se transmite con el factor ρ².
- **(cambio por la revisión) Prior plano en p explícito:** las dos posteriores gaussianas N(p_full, Σ_full) y N(p2, Σ2) suponen prior plano en
  (h, ω_m[, M]) e ignoran las cotas del optimizador (muy lejos de los mínimos).
- **Residuos:** r = datos_LRG2 − predicción ΛCDM ajustada sin LRG2, con el mismo procedimiento que PREDICTIVO_BINS (optimizador, priors,
  Fisher F = JᵀC⁻¹J). En DR3 se reajusta con los otros seis bins homólogos de DR3 (regla e), que también cambian, usando la covarianza de DR3
  (ρ²_k·C_kk en BAO, la original en los externos).
- Modelo de fondo ΛCDM plano con (h, ω_m) [+ M], priors de Planck comprimido (100 θ*, ω_m) como dos filas gaussianas.

## Cantidades de DR2 que fijan el pronóstico (reproducidas en v2 a 0.0 respecto al borrador)

Sin SN: p2 (ajuste sin LRG2) = (h, ω_m) = (0.67974, 0.14129); r2 = (-0.2096, -0.6353) en D_M/r_d, D_H/r_d; Q = J Σ2 Jᵀ = [[0.00325, 0.00184], [0.00184, 0.00105]];
C_LRG2 = [[0.03238, -0.02374], [-0.02374, 0.11147]]; V2 = [[0.03562, -0.02191], [-0.02191, 0.11252]]; d²_DR2 = 7.131 (p local 0.0283), idéntico a PREDICTIVO_BINS.
Con Pantheon+: r2 = (-0.2270, -0.6450); d²_DR2 = 7.701 (p local 0.0213).

## Fórmulas de referencia (aproximación gaussiana; NO son las distribuciones del procedimiento ni se usan para calcular Λ)

Con Q = J Σ2 Jᵀ (Σ2 la inversa de la Fisher del ajuste DR2 sin LRG2), V2 = C22 + Q, V3 = ρ²·C22 + Q y **(cambio por la revisión)** Cov(r3, r2)
**sin reajuste** = ρ²·C22 + Q (supone p3 − p2 = 0 y Σ3 = Σ2):

- **H0, sin reajuste:** r3 | r2 ~ N(K r2, V3 − K V2 Kᵀ), K = (ρ²·C22 + Q)·V2⁻¹. Para ρ² = 0.6, sin SN: media (-0.1453, -0.3922), covarianza [[0.00822, -0.00545], [-0.00545, 0.02690]].
- **H1 sin reajuste**, con r2 = μ + n2 y r3 = μ + ρ²·n2 + e:
  - condicionada a μ conocido: r3 | r2, μ ~ N(ρ²·r2 + (1 − ρ²)·μ, ρ²(1 − ρ²)·C22); equivalentemente r3 − r2 | r2, μ ~ N((1 − ρ²)(μ − r2), ρ²(1 − ρ²)·C22);
  - **marginalizada** sobre μ con prior plano (μ | r2 ~ N(r2, C22)): r3 | r2 ~ N(r2, (1 − ρ²)·C22). **(cambio por la revisión)** La igualdad
    ρ²r2 + (1 − ρ²)r2 = r2 y ρ²(1 − ρ²)C22 + (1 − ρ²)²C22 = (1 − ρ²)C22 es una **identidad algebraica**, no una comprobación; la prueba real de las
    simulaciones es la comparación con la distribución exacta del modelo linealizado (sección siguiente). Para ρ² = 0.6, sin SN: media
    (-0.2096, -0.6353), covarianza [[0.01295, -0.00950], [-0.00950, 0.04459]].
  - La referencia que se compara con las simulaciones es la marginalizada, porque en las simulaciones μ se saca de su posterior.
- **H1, con reajuste:** r3 − r2 ≈ (ρ² − 1)·n2 + e − J2·(p3 − p2), con p2 y p3 los parámetros ajustados sin LRG2 en DR2 y DR3. Esta es la expresión
  que corresponde al procedimiento, y la que realiza la distribución exacta de abajo.

## (cambio por la revisión) Densidades exactas del modelo linealizado y comparación con las simulaciones

En el procedimiento linealizado r3 es afín en variables gaussianas estándar z = (z_p, z_e[, z_μ]): r3 = a + M z, con a = r3(z = 0) y las columnas
de M = r3(e_k) − a (no linealidad residual ≤ 8e-15). **Bajo cada hipótesis r3 | datos DR2 ~ N(a, MMᵀ) exactamente**, y **Λ se calcula con estas
densidades** (no con medias y covarianzas estimadas). Las simulaciones (100000 universos por caso: mitad A, semillas [20261114, i_comb, i_ρ², i_hip, 0],
y mitad B, [.., 1], las mismas del borrador) quedan como comprobación:

| combinación | ρ² | hipótesis | a (media exacta) | Σ = MMᵀ (covarianza exacta) | media simulada (A+B) − a, en σ (1 σ_MC = 0.0032) | var. simulada (A) / exacta − 1 (1 σ_MC = 0.63 %) | media A del borrador − a, en σ |
|---|---|---|---|---|---|---|---|
| sinSN | 0.6 | H0 | (-0.1312, -0.3843) | [[0.00847, -0.00530], [-0.00530, 0.02698]] | -0.0042 / +0.0007 | +0.59 % / +0.79 % | -0.0012 / -0.0063 |
| sinSN | 0.6 | H1 | (-0.2096, -0.6353) | [[0.01389, -0.00897], [-0.00897, 0.04489]] | +0.0057 / +0.0054 | +0.30 % / -0.52 % | +0.0027 / +0.0107 |
| sinSN | 0.5 | H0 | (-0.1108, -0.3210) | [[0.00890, -0.00548], [-0.00548, 0.02813]] | -0.0029 / +0.0031 | +0.25 % / -0.19 % | +0.0014 / -0.0026 |
| sinSN | 0.5 | H1 | (-0.2096, -0.6353) | [[0.01742, -0.01118], [-0.01118, 0.05613]] | -0.0009 / -0.0036 | +1.02 % / +0.49 % | -0.0026 / -0.0053 |
| sinSN | 0.72 | H0 | (-0.1553, -0.4599) | [[0.00705, -0.00449], [-0.00449, 0.02264]] | +0.0072 / -0.0009 | +0.49 % / +0.03 % | +0.0022 / -0.0020 |
| sinSN | 0.72 | H1 | (-0.2096, -0.6353) | [[0.00969, -0.00630], [-0.00630, 0.03141]] | -0.0002 / -0.0023 | +0.63 % / +0.27 % | -0.0007 / -0.0022 |
| Pantheon+ | 0.6 | H0 | (-0.1422, -0.3904) | [[0.00843, -0.00532], [-0.00532, 0.02697]] | -0.0010 / +0.0011 | +0.55 % / -0.50 % | +0.0017 / -0.0033 |
| Pantheon+ | 0.6 | H1 | (-0.2270, -0.6450) | [[0.01382, -0.00901], [-0.00901, 0.04486]] | -0.0061 / +0.0022 | -0.12 % / +1.64 % | -0.0133 / +0.0091 |
| Pantheon+ | 0.5 | H0 | (-0.1201, -0.3262) | [[0.00886, -0.00550], [-0.00550, 0.02811]] | -0.0045 / +0.0029 | -0.36 % / -0.51 % | -0.0074 / +0.0054 |
| Pantheon+ | 0.5 | H1 | (-0.2270, -0.6450) | [[0.01733, -0.01123], [-0.01123, 0.05610]] | -0.0038 / -0.0022 | -0.22 % / +0.64 % | -0.0035 / +0.0002 |
| Pantheon+ | 0.72 | H0 | (-0.1682, -0.4671) | [[0.00702, -0.00451], [-0.00451, 0.02263]] | +0.0040 / -0.0011 | +0.07 % / -0.62 % | +0.0084 / -0.0018 |
| Pantheon+ | 0.72 | H1 | (-0.2270, -0.6450) | [[0.00964, -0.00633], [-0.00633, 0.03139]] | +0.0001 / -0.0004 | +0.73 % / +0.09 % | +0.0023 / -0.0062 |

Las simulaciones coinciden con la distribución exacta al nivel Monte Carlo (medias a ≤ 0.007 σ con A+B, varianzas a ≤ 1.6 %, E[d²] simulada frente a
exacta a ≤ 0.06). La media de A del borrador se desviaba hasta 0.013 σ (2.9 σ_MC) en los H1, lo que explicaba los fallos de los tests de ángulo de la
revisión; con las densidades exactas ese efecto desaparece. Frente a las fórmulas de referencia: bajo H1 la exacta tiene la media de la marginalizada
(r2) y una varianza de D_M un 7 % mayor (el término −J2(p3 − p2)); bajo H0 la exacta difiere de la de "sin reajuste" en 0.15 σ en D_M (condicionamiento
distinto). Esas diferencias son de las referencias, no de las simulaciones.

**(cambio por la revisión) Gaussianidad.** En el linealizado la gaussianidad de r3 es cierta por construcción (afín en gaussianas); la "comprobación de
gaussianidad" del borrador solo medía el ruido Monte Carlo de A. Lo que respalda la gaussianidad del r3 real es la validación con el modelo completo
(abajo): máx |Δr3| 7e-04 en D/r_d, es decir ≤ 0.009 σ de la componente más estrecha.

## Generación de DR3 (ecuaciones del procedimiento simulado)

- **H0.** p_true ~ N(p_full, Σ_full): aproximación gaussiana del posterior ΛCDM con todos los datos DR2 (p_full el ajuste completo, Σ_full la
  inversa de su Fisher; prior plano en p). Ruido DR2 implícito en los 7 bins: n2 = datos_DR2 − modelo(p_true). Datos DR3 = modelo(p_true) + ρ²·n2 + e,
  e ~ N(0, ρ²(1 − ρ²)·C). Externos sin cambio.
- **H1.** (p_true, μ) del posterior conjunto con μ libre y prior plano no acotado: p_true ~ N(p2, Σ2) y μ | p_true ~ N(datos_LRG2 −
  modelo_LRG2(p_true), C_LRG2). Con E₂μ el vector de los 7 bins con μ en las dos componentes de LRG2 y cero en los demás: ruido inferido
  n2 = datos_DR2 − modelo(p_true) − E₂μ; **datos DR3 = modelo(p_true) + E₂μ + ρ²·n2 + e**. Así μ permanece fijo entre DR2 y DR3 y solo el ruido se
  transmite con ρ². No se usa ρ²·(n2 + E₂μ), que reduciría el desvío persistente. Con ρ² por bin, ρ²·n2 y Cov(e) se aplican bloque a bloque.
- **En cada universo DR3:** reajuste ΛCDM sin LRG2 con la covarianza de DR3 (ρ²·C en BAO, la original en los externos) → p3; r3 = residuo de LRG2
  frente a la predicción de ese reajuste; V3,k = ρ²·C_LRG2 + J3 Σ3 J3ᵀ; d²_DR3 = r3ᵀ V3,k⁻¹ r3. **(cambio por la revisión) En los 100000 universos
  linealizados J, Σ3 y por tanto V3,k son fijas por caso** (J evaluada en p2; Σ3 la inversa de la Fisher con ρ²·C en p2); solo en la validación
  completa (20 universos por caso) se recalculan p3, J3, Σ3 y V3,k en el nuevo mínimo de cada universo.
- Comprobaciones del posterior de H1 (C_LRG2 invertible; bloque μ de la Fisher invertible; Fisher conjunta de (p, μ) definida positiva con inversa
  finita; autovalores registrados; posterior de μ propio porque la Fisher de p con los otros seis bins y los externos es definida positiva): sin
  cambio respecto al borrador, donde están los autovalores (`pronostico_dr3_lrg2.json` → `H1_posterior_checks`).

## Validación del linealizado frente al cálculo completo (del borrador, 12 casos × 20 universos; sin cambio)

| combinación | ρ² | hipótesis | máx \|Δd²\| | universo | máx \|Δr3\| (D/r_d) |
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

Mayor diferencia: 0.0441 (Pantheon+, ρ² = 0.5, H1, universo 0); ninguna supera 0.1. Respalda el linealizado en esta comprobación, no en general.

## Comparación con el borrador y errores Monte Carlo (cambio por la revisión)

El estimador del borrador (densidades gaussianas estimadas en A, evaluadas en B) se reproduce aquí a ≤ 1e-13 con las mismas semillas. Con las
densidades exactas y la misma mitad B:

| combinación | ρ² | borrador (densidades estimadas en A): c0 / c1 / potencia | reproducción de ese estimador aquí | densidades exactas (misma mitad B): c0 / c1 / potencia | error MC (20 reejecuciones de B) |
|---|---|---|---|---|---|
| sinSN | 0.6 | +1.011 / -1.575 / 0.6555 | +1.011 / -1.575 / 0.6555 (Δ ≤ 7e-16) | +1.024 / -1.604 / 0.6557 | ± 0.017 / ± 0.015 / ± 0.0031 |
| sinSN | 0.5 | +0.765 / -1.648 / 0.7670 | +0.765 / -1.648 / 0.7670 (Δ ≤ 8e-14) | +0.777 / -1.638 / 0.7669 | ± 0.018 / ± 0.019 / ± 0.0027 |
| sinSN | 0.72 | +1.187 / -1.484 / 0.5001 | +1.187 / -1.484 / 0.5001 (Δ ≤ 1e-15) | +1.184 / -1.476 / 0.4998 | ± 0.016 / ± 0.015 / ± 0.0048 |
| Pantheon+ | 0.6 | +0.988 / -1.570 / 0.6844 | +0.988 / -1.570 / 0.6844 (Δ ≤ 3e-15) | +0.993 / -1.579 / 0.6850 | ± 0.021 / ± 0.020 / ± 0.0037 |
| Pantheon+ | 0.5 | +0.671 / -1.567 / 0.7940 | +0.671 / -1.567 / 0.7940 (Δ ≤ 3e-14) | +0.667 / -1.560 / 0.7939 | ± 0.020 / ± 0.024 / ± 0.0027 |
| Pantheon+ | 0.72 | +1.200 / -1.483 / 0.5212 | +1.200 / -1.483 / 0.5212 (Δ ≤ 9e-16) | +1.195 / -1.478 / 0.5211 | ± 0.015 / ± 0.013 / ± 0.0033 |

Errores Monte Carlo: desviación típica (ddof = 1) sobre 20 reejecuciones de la mitad B con semillas nuevas [20261114, i_comb, i_ρ², i_hip, 1, 2000 + k],
k = 0..19 (registradas en el json), con las densidades exactas. c1 < c0 en las 20 reejecuciones de los 6 casos.

## Resultados: d²_DR3 y Λ (densidades exactas; mitad B; ± error MC)

| combinación | ρ² | hipótesis | d²_DR3 pct 5 / 50 / 95 (B) | mediana ± MC | E[d²] exacto |
|---|---|---|---|---|---|
| sinSN | 0.6 | H0 | 1.40 / 4.76 / 10.18 | 4.76 ± 0.01 | 5.146 |
| sinSN | 0.6 | H1 | 4.92 / 12.29 / 23.13 | 12.29 ± 0.03 | 12.956 |
| sinSN | 0.5 | H0 | 0.91 / 4.15 / 10.04 | 4.15 ± 0.01 | 4.639 |
| sinSN | 0.5 | H1 | 5.21 / 14.93 / 29.57 | 14.93 ± 0.03 | 15.852 |
| sinSN | 0.72 | H0 | 2.23 / 5.45 / 10.14 | 5.45 ± 0.01 | 5.748 |
| sinSN | 0.72 | H1 | 4.82 / 10.14 / 17.53 | 10.14 ± 0.03 | 10.535 |
| Pantheon+ | 0.6 | H0 | 1.61 / 5.10 / 10.68 | 5.10 ± 0.02 | 5.500 |
| Pantheon+ | 0.6 | H1 | 5.54 / 13.27 / 24.32 | 13.27 ± 0.05 | 13.871 |
| Pantheon+ | 0.5 | H0 | 1.04 / 4.47 / 10.39 | 4.47 ± 0.02 | 4.938 |
| Pantheon+ | 0.5 | H1 | 5.93 / 16.04 / 31.31 | 16.04 ± 0.04 | 16.936 |
| Pantheon+ | 0.72 | H0 | 2.52 / 5.91 / 10.73 | 5.91 ± 0.01 | 6.168 |
| Pantheon+ | 0.72 | H1 | 5.38 / 10.94 / 18.44 | 10.94 ± 0.02 | 11.308 |

| combinación | ρ² | hipótesis | Λ pct 5 / 50 / 95 (B) | P(Λ > 3) / P(Λ < −3), ± MC |
|---|---|---|---|---|
| sinSN | 0.6 | H0 | -3.36 / -1.68 / +1.02 | 0.005 ± 0.000 / 0.118 ± 0.001 |
| sinSN | 0.6 | H1 | -1.60 / +2.08 / +7.50 | 0.374 ± 0.002 / 0.004 ± 0.000 |
| sinSN | 0.5 | H0 | -3.78 / -2.17 / +0.78 | 0.005 ± 0.000 / 0.247 ± 0.002 |
| sinSN | 0.5 | H1 | -1.64 / +3.22 / +10.54 | 0.524 ± 0.002 / 0.008 ± 0.000 |
| sinSN | 0.72 | H0 | -2.77 / -1.16 / +1.18 | 0.004 ± 0.000 / 0.025 ± 0.001 |
| sinSN | 0.72 | H1 | -1.48 / +1.18 / +4.88 | 0.194 ± 0.002 / 0.001 ± 0.000 |
| Pantheon+ | 0.6 | H0 | -3.54 / -1.80 / +0.99 | 0.005 ± 0.000 / 0.154 ± 0.002 |
| Pantheon+ | 0.6 | H1 | -1.58 / +2.29 / +7.81 | 0.405 ± 0.003 / 0.005 ± 0.000 |
| Pantheon+ | 0.5 | H0 | -4.01 / -2.29 / +0.67 | 0.004 ± 0.000 / 0.293 ± 0.002 |
| Pantheon+ | 0.5 | H1 | -1.56 / +3.49 / +11.13 | 0.551 ± 0.002 / 0.009 ± 0.001 |
| Pantheon+ | 0.72 | H0 | -2.91 / -1.22 / +1.19 | 0.005 ± 0.000 / 0.038 ± 0.001 |
| Pantheon+ | 0.72 | H1 | -1.48 / +1.30 / +5.05 | 0.212 ± 0.002 / 0.001 ± 0.000 |

## Regla de decisión (calibrada con la mitad B, semillas del borrador; cambio por la revisión: con errores MC)

| combinación | ρ² | c0 | c1 | c1 < c0 (y en las 20 reejecuciones) | potencia P(Λ > c0 \| H1) | P(Λ < c1 \| H0) | no concluyente H0 / H1 |
|---|---|---|---|---|---|---|---|
| sinSN | 0.6 | 1.024 ± 0.017 | -1.604 ± 0.015 | sí | 0.656 ± 0.003 | 0.524 ± 0.005 | 0.426 ± 0.005 / 0.294 ± 0.003 |
| sinSN | 0.5 | 0.777 ± 0.018 | -1.638 ± 0.019 | sí | 0.767 ± 0.003 | 0.640 ± 0.005 | 0.310 ± 0.005 / 0.183 ± 0.003 |
| sinSN | 0.72 | 1.184 ± 0.016 | -1.476 ± 0.015 | sí | 0.500 ± 0.005 | 0.394 ± 0.005 | 0.556 ± 0.005 / 0.450 ± 0.005 |
| Pantheon+ | 0.6 | 0.993 ± 0.021 | -1.579 ± 0.020 | sí | 0.685 ± 0.004 | 0.561 ± 0.006 | 0.389 ± 0.006 / 0.265 ± 0.004 |
| Pantheon+ | 0.5 | 0.667 ± 0.020 | -1.560 ± 0.024 | sí | 0.794 ± 0.003 | 0.683 ± 0.005 | 0.267 ± 0.005 / 0.156 ± 0.003 |
| Pantheon+ | 0.72 | 1.195 ± 0.015 | -1.478 ± 0.013 | sí | 0.521 ± 0.003 | 0.416 ± 0.004 | 0.534 ± 0.004 / 0.429 ± 0.003 |

**Regla.** En los 6 casos c1 < c0, así que la regla de tres resultados es aplicable **bajo los supuestos declarados arriba** (anidamiento exacto con
ρ² común o por bin medido, externos fijos, linealización validada, ΛCDM con Planck comprimido): **Λ > c0 favorece H1** (falsa alarma del 5 % bajo H0,
bajo esos supuestos); **Λ < c1 favorece H0** (error del 5 % bajo H1, bajo esos supuestos); en medio, **no concluyente**. **(cambio por la revisión)**
|Λ| > 3 se usa como indicador adicional de **evidencia fuerte de una hipótesis frente a la otra, dentro de este par y bajo estos supuestos**, sin
sustituir a c0 y c1; es asimétrico: P(Λ < −3 | H0) = 0.118 ± 0.001 frente a P(Λ > 3 | H1) = 0.374 ± 0.002 (ρ² = 0.6, sin SN).
**(cambio por la revisión)** Con ρ² = 0.72 la potencia es 0.500 ± 0.005 (densidades exactas: 0.4998; con las densidades estimadas del borrador,
0.501 ± 0.004 como mide la revisión): tres cifras no son significativas, y DR3 separa las dos hipótesis solo a medias si el ruido de DR2 pesa tanto.

## (cambio por la revisión) Potencia frente a ρ² común, 0.40–0.80 (densidades exactas; c0 y c1 de 50000 universos B por hipótesis con semillas [20261114, i_comb, 100 + k, i_hip, 1])

| ρ² común | sin SN: c0 / c1 | potencia | P(Λ < c1 \| H0) | no concl. H0 / H1 | mediana d² H0 / H1 | Pantheon+: c0 / c1 | potencia | mediana d² H0 / H1 |
|---|---|---|---|---|---|---|---|---|
| 0.40 | +0.363 / -1.509 | 0.855 | 0.771 | 0.178 / 0.095 | 3.55 / 18.90 | +0.254 / -1.436 | 0.872 | 3.80 / 20.02 |
| 0.42 | +0.423 / -1.573 | 0.839 | 0.738 | 0.212 / 0.111 | 3.68 / 17.91 | +0.363 / -1.474 | 0.855 | 3.93 / 19.07 |
| 0.44 | +0.543 / -1.552 | 0.823 | 0.721 | 0.229 / 0.127 | 3.81 / 17.06 | +0.432 / -1.545 | 0.843 | 4.07 / 18.20 |
| 0.46 | +0.606 / -1.609 | 0.806 | 0.690 | 0.260 / 0.144 | 3.93 / 16.25 | +0.543 / -1.533 | 0.823 | 4.21 / 17.41 |
| 0.48 | +0.679 / -1.608 | 0.787 | 0.668 | 0.282 / 0.163 | 4.04 / 15.53 | +0.629 / -1.545 | 0.808 | 4.35 / 16.72 |
| 0.50 | +0.751 / -1.620 | 0.766 | 0.646 | 0.304 / 0.184 | 4.14 / 14.90 | +0.691 / -1.576 | 0.787 | 4.43 / 15.93 |
| 0.52 | +0.833 / -1.625 | 0.745 | 0.620 | 0.330 / 0.205 | 4.28 / 14.30 | +0.787 / -1.572 | 0.767 | 4.60 / 15.41 |
| 0.54 | +0.924 / -1.619 | 0.721 | 0.595 | 0.355 / 0.229 | 4.41 / 13.75 | +0.838 / -1.618 | 0.745 | 4.73 / 14.77 |
| 0.56 | +0.945 / -1.616 | 0.702 | 0.572 | 0.378 / 0.248 | 4.52 / 13.25 | +0.931 / -1.580 | 0.722 | 4.85 / 14.21 |
| 0.58 | +0.972 / -1.601 | 0.684 | 0.550 | 0.400 / 0.266 | 4.64 / 12.79 | +0.958 / -1.602 | 0.703 | 4.99 / 13.73 |
| 0.60 | +1.081 / -1.591 | 0.652 | 0.529 | 0.421 / 0.298 | 4.75 / 12.38 | +1.021 / -1.577 | 0.677 | 5.14 / 13.24 |
| 0.62 | +1.060 / -1.580 | 0.637 | 0.502 | 0.448 / 0.313 | 4.89 / 11.95 | +1.053 / -1.574 | 0.656 | 5.26 / 12.80 |
| 0.64 | +1.129 / -1.550 | 0.604 | 0.486 | 0.464 / 0.346 | 4.99 / 11.53 | +1.064 / -1.559 | 0.640 | 5.38 / 12.43 |
| 0.66 | +1.156 / -1.540 | 0.584 | 0.460 | 0.490 / 0.366 | 5.11 / 11.18 | +1.141 / -1.529 | 0.606 | 5.53 / 12.04 |
| 0.68 | +1.172 / -1.510 | 0.554 | 0.436 | 0.514 / 0.396 | 5.26 / 10.81 | +1.160 / -1.517 | 0.575 | 5.63 / 11.64 |
| 0.70 | +1.177 / -1.499 | 0.525 | 0.413 | 0.537 / 0.425 | 5.36 / 10.44 | +1.157 / -1.513 | 0.558 | 5.77 / 11.32 |
| 0.72 | +1.185 / -1.470 | 0.499 | 0.389 | 0.561 / 0.451 | 5.49 / 10.14 | +1.199 / -1.475 | 0.521 | 5.91 / 10.95 |
| 0.74 | +1.211 / -1.452 | 0.471 | 0.363 | 0.587 / 0.479 | 5.60 / 9.87 | +1.257 / -1.446 | 0.479 | 6.02 / 10.60 |
| 0.76 | +1.204 / -1.389 | 0.444 | 0.353 | 0.597 / 0.506 | 5.72 / 9.58 | +1.222 / -1.430 | 0.462 | 6.17 / 10.35 |
| 0.78 | +1.202 / -1.376 | 0.412 | 0.321 | 0.629 / 0.538 | 5.83 / 9.31 | +1.204 / -1.372 | 0.436 | 6.29 / 10.07 |
| 0.80 | +1.199 / -1.316 | 0.383 | 0.306 | 0.644 / 0.567 | 5.96 / 9.04 | +1.219 / -1.362 | 0.397 | 6.42 / 9.76 |

Las potencias tienen un error MC ≈ ±0.004 y c0, c1 ≈ ±0.02 (como en los casos principales).

## Veredicto de P18–P19 (sin cambio de enunciado; cifras con error MC)

- **P18: sí.** ρ² = 0.6, falsa alarma 5 %: potencia 0.656 ± 0.003 sin SN (con Pantheon+ 0.685 ± 0.004), dentro de [0.5, 0.75].
- **P19: sí.** ρ² = 0.6, sin SN: mediana de d²_DR3 4.76 ± 0.01 bajo H0 (en [3, 6]) y 12.29 ± 0.03 bajo H1 (en [9, 15]). Con Pantheon+:
  5.10 y 13.27.

Nota sobre d² bajo H0: su mediana (4.8) es mayor que la de una χ²₂ (1.39) porque DR3 hereda el 60 % del ruido de DR2 y la V3,k del reajuste no incluye
esa correlación con DR2: d²_DR3 no es un p-valor de DR3 por sí solo; lo que se compara es Λ, que sí condiciona en DR2.

## Modo "evaluación DR3" (escrito y probado con datos falsos; todavía sin usar)

`pronostico_dr3_lrg2_v2.py calibrar COV.json TAG` lee **solo** la covarianza 13 × 13 de DR3 (etiquetas en el orden de DR2), calcula ρ²_k, aplica las
reglas b y c, recalibra c0 y c1 para los dos combos con esos ρ²_k, las densidades exactas y las semillas fijadas [20261114, i_comb, 500, i_hip, 1], y
escribe `..._calibracion_TAG.json` con el sha256 de la covarianza, del script y del commit, más su propio sha256 en `..._calibracion_TAG.sha256`.
Solo después, `evaluar COV.json VALORES.json TAG` comprueba los tres sha256, lee los valores de DR3, reajusta ΛCDM sin LRG2 con el optimizador
original y la covarianza de DR3 (cálculo completo, no linealizado), calcula r3, d²_DR3 y **Λ una sola vez** y se niega a sobrescribir un resultado.
Prueba con **datos falsos** (DR2 reescalado; `pronostico_dr3_lrg2_v2_falso_run.log`): covarianza 0.6·C_DR2 y valores de DR2 sin cambio →
estado "ok", ρ²_k = 0.6 en los 7 bins, c0 / c1 = +1.013 / -1.599 (sin SN) y +0.978 / -1.573 (Pantheon+); evaluación: r3 = (-0.182, -0.620),
d² = 10.38, Λ = +1.127 → "favorece H1" (sin SN; Pantheon+ Λ = +1.169, "favorece H1"), que es lo esperado al repetir los valores de DR2 con menos ruido; la segunda
evaluación con el mismo TAG se niega ("Λ se calcula una sola vez"); con covarianza 0.3·C_DR2 (ρ²_LRG2 = 0.30 fuera de [0.40, 0.80]) la calibración
declara "modelo roto" y la evaluación devuelve "modelo roto" sin calcular Λ. Esos ficheros llevan "DATOS FALSOS" en su cabecera.

## Qué rompería el modelo y debe declararse si ocurre (y cómo se comprobaría al salir DR3)

1. **Reprocesado de DR3 que cambie la medida de DR2**: comparar el vector de 13 números y la covarianza de DR2 usados aquí (sha256 en
   `data/PROCEDENCIA_DESI_DR2.md`) con lo que DESI republique; si cambian, recalcular r2 y declarar.
2. **Cambio de z_eff o de bins**: si LRG2 no existe con z_eff = 0.706 ± 0.01, no evaluable.
3. **Cambio de cosmología fiducial de DESI**: comparar con la de DR2; si cambia, recalcular con la nueva o declarar no evaluable.
4. **ρ² distinto por bin**: se mide (regla c) y se recalibra con el script congelado.
5. **Sistemáticos nuevos**: comparar la covarianza 13 × 13 de DR3 con ρ²_k·C_kk bin a bin; si la diferencia supera lo que explica un ρ² por bin
   (regla c), declarar.
6. **Cambio en los datos externos**: el pronóstico se evalúa con exactamente los externos de este json (θ*, ω_m, Pantheon+ 1580 SN con sus sha256 en
   `datos/FUENTES.json`); con otros, es otro pronóstico.

**(cambio por la revisión) Rupturas añadidas, con su regla:**

- **a. Formato de DR3.** Evaluable si LRG2 viene como (D_M/r_d, D_H/r_d) con covarianza. Si viene como (α_iso, α_AP) con cosmología fiducial declarada,
  conversión determinista a (D_M/r_d, D_H/r_d) con esa fiducial (con α_iso = (α_∥α_⊥²)^{1/3} y α_AP = α_∥/α_⊥: D_M/r_d = α_⊥·(D_M/r_d)_fid = α_iso·α_AP^{−1/3}·(D_M/r_d)_fid, D_H/r_d = α_∥·(D_H/r_d)_fid = α_iso·α_AP^{2/3}·(D_H/r_d)_fid, y
  la covarianza por el jacobiano), registrada antes de evaluar. Si solo viene D_V o LRG2 no está separado: **no evaluable**.
- **b. Correlaciones entre bins** distintas de cero en la covarianza de DR3: **modelo roto**, no evaluable con este pronóstico (el script lo comprueba:
  `max_abs_offblock_cov3`).
- **c. ρ²_k operativo.** Cociente de diagonales σ²_DR3/σ²_DR2 por componente; ρ²_k del bin = media de sus componentes. **Modelo roto** si las componentes
  de un bin difieren más del 20 % (relativo a su media), si el coeficiente de correlación D_M–D_H de LRG2 cambia más de 0.1, o si ρ²_LRG2 queda fuera de
  [0.40, 0.80]. Si no hay ruptura: se recalibran c0 y c1 con los ρ²_k medidos usando el script congelado sin ningún otro cambio, antes de leer los
  valores de DR3 (modo `calibrar`).
- **d. Código y densidades congelados.** Se citan el hash del commit de `pronostico_dr3_lrg2_v2.py` y de su json (los sha256 de arriba; el commit es
  el de este registro) y Λ se calcula una sola vez con ese código (modo `evaluar`, que comprueba los sha256 y no sobrescribe).
- **e. "Los otros seis bins"** son los seis homólogos de DR2 (mismos trazadores y rangos de z: BGS, LRG1, LRG3+ELG1, ELG2, QSO, Lyα); los bins nuevos de
  DR3 se ignoran en el reajuste.
- **f. Erratas propias** descubiertas después en los datos de entrada (DR2 o externos): se declaran y se recalcula con registro (ERRATAS_V3.md), antes
  de mirar Λ.
- **g. Fecha límite.** Si las medidas BAO de DESI DR3 no son públicas antes del 31-12-2028, el pronóstico se declara **no evaluado** (no fallido).

En cualquiera de estos casos la evaluación se declara "modelo roto" o "no evaluable" antes de mirar Λ, y se documenta con los números que lo muestran.

**(cambio por la revisión) Compromiso.** Se informará el resultado con una de estas cinco salidas: **favorece H1 / favorece H0 / no concluyente /
modelo roto / no evaluado**, con estos umbrales (c0 y c1 calibrados como se describe, |Λ| > 3 solo como indicador adicional), sin añadir umbrales después.

## Cautelas

- Las densidades y umbrales dependen de los supuestos de arriba y de la linealización, validada solo en estas condiciones.
- H1 es un desvío constante en un solo bin: otras hipótesis (desvío que cambia con la muestra, física que afecte a varios bins) no están
  representadas, y Λ solo compara estas dos. Ni H1 ni Λ prueban el modelo OU ni σ_X > 0.
- Pendiente de comprobación del revisor; sin tag.

## Ficheros nuevos de esta versión

`pronostico_dr3_lrg2_v2.py`, `pronostico_dr3_lrg2_v2.json`, `pronostico_dr3_lrg2_v2_run.log`; prueba del modo evaluación con datos falsos:
`pronostico_dr3_lrg2_v2_falso_dr3_cov.json`, `..._falso_roto_dr3_cov.json`, `..._falso_dr3_valores.json`, `..._calibracion_falso.json` (+ `.sha256`),
`..._resultado_falso.json`, `..._calibracion_falso_roto.json` (+ `.sha256`), `..._resultado_falso_roto.json`, `pronostico_dr3_lrg2_v2_falso_run.log`.
Nada existente se modifica.
