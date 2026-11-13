# REVISIÓN INDEPENDIENTE del pronóstico DR3 de LRG2 (13-11-2026)

Revisor: Claude (sesión independiente; no participó en el borrador). Objeto: `PRONOSTICO_DR3_BORRADOR.md`, `pronostico_dr3_lrg2.py`,
`pronostico_dr3_lrg2.json`, `pronostico_dr3_lrg2_run.log` (commit c742989), con `PREDICTIVO_BINS.md` y `JACKKNIFE_BINS.md` como contexto.
Ningún fichero existente se modifica. Todo lo nuevo lleva prefijo `revision_pronostico_`:

- `revision_pronostico_reproduccion.py` → `.json`: reproducción con las semillas del json (reimplementación de la generación a partir de las
  ecuaciones del borrador, sin llamar a `simulate` del original), media y covarianza **exactas** de r3 en el modelo linealizado (r3 es afín en
  variables gaussianas), decisión con densidades exactas y error Monte Carlo (20 reejecuciones con otras semillas).
- `revision_pronostico_gaussianidad.py` → `.json`: tests bidimensionales en la mitad B blanqueada con media y covarianza de A (Rayleigh,
  Kuiper, KS del radio frente a χ con 2 g.d.l., radio por 8 sectores, independencia radio–ángulo, asimetría, curtosis, correlación) y, al fallar
  un test, c0, c1 y potencia con KDE (ancho de banda por validación cruzada solo dentro de A; B solo calibra y evalúa).
- `revision_pronostico_gauss_exacta.py` → `.json`: diagnóstico de los fallos de ángulo (blanqueo con el centro exacto y 20 mitades B nuevas).

## 1. Matemáticas (deducciones propias)

**Anidamiento.** Con n3 = ρ²·n2 + e, e ⊥ n2, Cov(n2) = C, Cov(e) = ρ²(1 − ρ²)·C:
Cov(n3) = ρ⁴C + ρ²(1 − ρ²)C = ρ²C; Cov(n3, n2) = ρ²·Cov(n2) = ρ²C. Correcto. Equivale a la combinación por inversa de la varianza de DR2
con dos años nuevos independientes (1/σ3² = 1/σ2² + 1/σ_nuevo²: el peso de DR2 es σ3²/σ2² = ρ² y la parte nueva tiene varianza (1 − ρ²)σ3² =
ρ²(1 − ρ²)σ2²), y también a la propiedad del estimador eficiente anidado Cov(θ̂_sub, θ̂_total) = Var(θ̂_total). La correlación entre n2 y n3
es ρ (no ρ²); el borrador usa ρ² solo como cociente de varianzas, de forma consistente.

**Generación.** H0: d3 = m(p_true) + ρ²(d2 − m(p_true)) + e = (1 − ρ²)m(p_true) + ρ²d2 + e. H1: con n2 = d2 − m(p_true) − E₂μ,
d3 = m(p_true) + E₂μ + ρ²n2 + e = (1 − ρ²)(m(p_true) + E₂μ) + ρ²d2 + e: μ entra con peso 1 (ρ² + (1 − ρ²)) y el ruido con ρ². Correcto, y así
está en el código (líneas 207–215 del script: `n2[:, LRG2] -= mu` antes del factor ρ² y `d3[:, LRG2] += mu` después). La posterior de μ bajo
H1 con prior plano, marginal en p, es N(r2, C22 + J2Σ2J2ᵀ) = N(r2, V2), que coincide con la comprobación (4) del borrador.

**Fórmulas de referencia.** H1 sin reajuste: r3 | r2, μ ~ N(ρ²r2 + (1 − ρ²)μ, ρ²(1 − ρ²)C22); con μ | r2 ~ N(r2, C22), la marginal es
N(r2, ρ²(1 − ρ²)C22 + (1 − ρ²)²C22) = N(r2, (1 − ρ²)C22). Correctas dentro de su modelo "sin reajuste y sin Q". La que se compara con las
simulaciones (la marginalizada) es la adecuada porque μ se saca de su posterior. **Pero** ninguna de las dos referencias del borrador es la
distribución del procedimiento simulado: la del procedimiento se obtiene exactamente porque, linealizado, r3 = a + Mz con z gaussiana estándar
(comprobado: no linealidad residual 7e-15). La he calculado (`exact_linear` en mi json) y:

| caso (ρ² = 0.6) | hip. | media simulada (A) − exacta, en σ (D_M / D_H); 1 σ_MC = 0.0045 | var. simulada/exacta − 1 (1 σ_MC = 0.006) | referencia del borrador − exacta: media en σ / var. |
|---|---|---|---|---|
| sinSN | H0 | −0.001 / −0.006 | +0.6 % / +0.8 % | −0.152 / −0.048 σ; −3.0 % / −0.3 % |
| sinSN | H1 | +0.003 / +0.011 | +0.3 % / −0.5 % | 0 / 0; −6.7 % / −0.7 % |
| Pantheon+ | H0 | +0.002 / −0.003 | +0.6 % / −0.5 % | −0.149 / −0.047 σ; −2.8 % / −0.3 % |
| Pantheon+ | H1 | −0.013 / +0.009 | −0.1 % / +1.6 % | 0 / 0; −6.3 % / −0.6 % |

Las simulaciones coinciden con la distribución exacta al nivel Monte Carlo (la peor desviación, −0.013 σ en D_M de Pantheon+ H1, es 2.9 σ_MC;
40 mitades nuevas dan media agrupada − exacta = −0.0006 σ ± 0.0007: la exacta es la correcta y la mitad A original es una fluctuación).
Las lecturas del borrador sobre sus referencias (−J2(p3 − p2) explica el +7 % en D_M bajo H1; el condicionamiento distinto explica los 0.15 σ
bajo H0) quedan confirmadas: son diferencias de la referencia, no de la simulación. E[d²] exacto (ρ² = 0.6, sinSN): 5.146 (H0), 12.956 (H1);
medias simuladas en B: 5.142 y 12.930.

## 2. Código frente a especificación

- Generación: coincide con lo escrito; E₂μ fuera del factor ρ² (véase arriba). Mi reimplementación desde las ecuaciones reproduce las 50000 × 2
  r3 de cada mitad a 1e-16.
- Reajuste DR3: `fisher(combo, KEEP, p2, rho2 * COV13)` usa ρ²C en los 11 números BAO conservados y σ_θ*, σ_ω y C_SN originales en los externos
  (`cinv_apply`); el vector de externos es el de DR2 (`const_rest`). En la validación, `fit_c(..., d3, cov3)` con `d3 = D.copy()` y solo
  `d3[:13]` sustituido. Correcto.
- V3,k: **en los 100000 universos linealizados es una sola matriz por caso** (Σ3 y J en p2), como declara el borrador en "Versión linealizada";
  solo en la validación (20 universos) se recalcula con el Σ3 de cada universo. El párrafo "En cada universo DR3 … con la V3,k de ese mismo
  universo" describe el procedimiento completo, no el que produce los números; véase (b).
- Mitades: semillas [20261114, i_comb, i_ρ², i_hip, 0/1/2] → `default_rng` con SeedSequence distintas: independientes. A se usa para media y
  covarianza (densidades y blanqueo); B para gaussianidad, percentiles de d², c0, c1, potencia. Única fuga: la elección gaussiana/KDE de la
  densidad (`gaussian_ok`) se decide con B (sin efecto: siempre gaussiana).
- Tiempo de ejecución 1.4 s: plausible (un ajuste completo tarda 4–5 ms, la Fisher 1–15 ms; 480 ajustes + 12 Fisher + álgebra vectorizada).

## 3. Reproducción (semillas del json)

ρ² = 0.6, sinSN: c0 = 1.0114702263158493 (Δ 1e-14), c1 = −1.5745851795775228 (Δ 8e-15), potencia 0.65546 (Δ 0), medianas de d² 4.7621 / 12.2913
(Δ 3e-14), medias y covarianzas de A (Δ ≤ 1e-15), V3,k (Δ 0). También idénticos ρ² = 0.5 y 0.72 (sinSN) y Pantheon+ 0.6. **Reproducible.**

Error Monte Carlo (20 reejecuciones completas con otras semillas, misma N): c0 = 1.030 ± 0.024; c1 = −1.591 ± 0.015; potencia 0.658 ± 0.004;
mediana d² H0 4.765 ± 0.016; H1 12.317 ± 0.026. ρ² = 0.5: potencia 0.768 ± 0.002; ρ² = 0.72: 0.501 ± 0.004; Pantheon+ 0.6: 0.680 ± 0.003.
Con las densidades exactas (sin estimar) en la misma B: c0 = 1.024, c1 = −1.604, potencia 0.6557 (sinSN); 0.993, −1.579, 0.685 (Pantheon+).
P18 (0.655 en [0.5, 0.75]) y P19 (4.76 en [3, 6]; 12.29 en [9, 15]) están lejos de los bordes respecto a estos errores.

## 4. Gaussianidad bidimensional (mitad B blanqueada con A; ρ² = 0.6; `revision_pronostico_gaussianidad.json`)

| caso | Rayleigh p | Kuiper p | KS radio p | Kruskal 8 sectores p | indep. 8×5 p | asimetría (/SE) | exceso curtosis (/SE) | corr (/SE) |
|---|---|---|---|---|---|---|---|---|
| sinSN H0 | 0.013 | **0.005** | 0.55 | 0.68 | 0.09 | −0.000, +0.017 (1.5) | +0.050 (2.3), −0.005 | +0.003 (0.6) |
| sinSN H1 | 0.13 | 0.07 | 0.12 | 0.08 | 0.29 | +0.015 (1.4), +0.010 | −0.020, −0.022 (−1.0) | −0.015 (−3.3) |
| Pantheon+ H0 | 0.54 | 0.76 | 0.99 | 0.28 | 0.71 | +0.002, −0.002 | −0.010, +0.031 (1.4) | −0.000 |
| Pantheon+ H1 | **0.001** | **< 0.001** | 0.12 | 0.67 | 0.98 | −0.006, +0.013 (1.2) | +0.029 (1.3), −0.002 | −0.008 (−1.7) |

Radio ~ χ₂, independencia radio–ángulo, asimetría, curtosis: todo compatible. Fallan los tests de **ángulo** en sinSN H0 y Pantheon+ H1. Diagnóstico
(`revision_pronostico_gauss_exacta.json`): es el desplazamiento del centro, no una no gaussianidad. La media de A difiere de la exacta en
0.007 σ (sinSN H0) y 0.013–0.014 σ (los dos H1), y con n = 50000 un desplazamiento así hace no uniforme el ángulo. Blanqueando B con la media y
covarianza exactas, todos los tests pasan (Rayleigh p = 0.10 / 0.25 / 0.90 / 0.36; Kuiper 0.05 / 0.40 / 0.89 / 0.06); en 20 mitades B nuevas, con el
centro de A las p de Rayleigh tienen mediana 0.003–0.006 en los H1 y con el centro exacto 0.38–0.40. La correlación −0.015 (sinSN H1) está a 2.3 σ del
error combinado de A y B (0.0063). Además r3 es exactamente afín en variables gaussianas (no linealidad 7e-15): la gaussianidad es cierta por construcción.

Efecto medido con KDE (ancho de banda por validación cruzada de verosimilitud dentro de A: 0.21 / 0.19 de la covarianza; B solo calibra y evalúa):

| combinación | c0 KDE / gauss | c1 KDE / gauss | potencia KDE / gauss |
|---|---|---|---|
| sinSN | 0.971 / 1.011 | −1.527 / −1.575 | 0.654 / 0.6555 |
| Pantheon+ | 0.943 / 0.988 | −1.516 / −1.570 | 0.685 / 0.6844 |

Δc0 ≈ −0.04 (1.7 σ_MC, sesgo de suavizado del KDE), potencia sin cambio. Ninguna conclusión cambia.

## (a) Errores que cambian conclusiones

**Ninguno encontrado.** Las matemáticas del anidamiento, de la generación y de las referencias son correctas; el código hace lo que dice el
borrador; la reproducción es bit a bit; la distribución simulada coincide con la exacta del modelo linealizado al nivel Monte Carlo; c0, c1,
potencia y medianas son estables frente a densidades exactas, KDE y semillas nuevas; P18 y P19 se mantienen con margen.

## (b) Errores menores (no cambian ninguna cifra de P18–P19)

1. **"Comprobación de la integración (numérica, en el json)"** (sección Fórmulas de referencia): el código evalúa la identidad algebraica
   ρ²r2 + (1 − ρ²)r2 − r2 y ρ²(1 − ρ²)C + (1 − ρ²)²C − (1 − ρ²)C, que dan 0.0 por construcción. No comprueba nada; debe decirse que es una
   identidad o sustituirse por la comparación con la distribución exacta del procedimiento (sección 1 de esta nota), que sí es una prueba y
   que las simulaciones superan.
2. **Sin errores Monte Carlo** en c0, c1, potencia, medianas y en los P(Λ > 3): c0 ± 0.024, c1 ± 0.015–0.019, potencia ± 0.004, medianas ± 0.02–0.04.
   "Potencia 0.500" con ρ² = 0.72 es 0.501 ± 0.004 (tres cifras no significativas). Para un registro con fecha deben constar.
3. **La "comprobación de gaussianidad" no puede fallar**: en el linealizado r3 es afín en gaussianas, así que solo mide el ruido Monte Carlo de A
   (el paréntesis "es lo esperado" lo insinúa; debe decirse claro). Con tolerancias del 10 % en percentiles no detectaría tampoco una no gaussianidad
   moderada. Lo que sí sostiene la gaussianidad del r3 real es la validación con el modelo completo (máx |Δr3| 7e-4, es decir ≤ 0.008 σ), y eso es
   lo que hay que citar.
4. **Media de A con error apreciable para tests finos**: 0.013–0.014 σ en los dos H1 (2.9 σ_MC). Sin efecto en Λ (densidades exactas: potencia
   0.6557 frente a 0.6555), pero las densidades registradas podrían ser las exactas del modelo linealizado (media a y covarianza MMᵀ, sin error
   Monte Carlo), dejando las simulaciones como comprobación.
5. **"En cada universo DR3 … Σ3 … V3,k de ese mismo universo"** describe el procedimiento completo; el que produce los números usa una V3,k
   fija por caso. Hay que decirlo en ese mismo párrafo, no dos párrafos después.
6. **Prior plano en p implícito** (y cotas del optimizador ignoradas) en las dos posteriores gaussianas N(p_full, Σ_full) y N(p2, Σ2). Debe constar.
7. **Cov(r3, r2) = ρ²·C22 + Q** no es la covarianza del procedimiento (supone p3 − p2 = 0 y Σ3 = Σ2); está etiquetada como referencia, pero
   conviene llamarla explícitamente "sin reajuste" también en la línea de H0.
8. **Elección gaussiana/KDE decidida con B** (`gaussian_ok`): fuga mínima de B hacia la densidad; sin efecto aquí.
9. **Cita de los pesos**: comprobado en arXiv:2607.28918 v3 §4.4: pesos 1/9, 3/9, 5/9 para DR1, DR2, DR3 en proporción al tiempo de observación
   (≈ 1, 3, 5 años). El borrador omite el 1/9 de DR1 (sin consecuencia) y no dice que ρ² = 3/5 supone σ² ∝ 1/tiempo, que es otro supuesto (los
   errores de BAO no escalan exactamente con el tiempo: área frente a profundidad, varianza cósmica).

## 5. Texto: frases que dicen más de lo que muestran los cálculos

- "**|Λ| > 3 se usa como indicador adicional de evidencia fuerte**": Λ es un cociente de densidades predictivas de dos números entre dos hipótesis
  concretas, una de ellas (H1) construida a partir del propio residuo de DR2 en LRG2, y es asimétrico (P(Λ < −3 | H0) = 0.10 frente a
  P(Λ > 3 | H1) = 0.37). Debe decir "evidencia fuerte de una frente a la otra dentro de este par y bajo los supuestos", no "evidencia fuerte".
- "**Λ > c0 favorece H1 (falsa alarma del 5 % bajo H0)**" y "potencia": los dos números son bajo los supuestos (anidamiento exacto con ρ² común,
  externos fijos, linealización, ΛCDM con Planck comprimido). Las cautelas lo dicen; la regla de decisión debería repetirlo en la misma frase.
- Falta una frase explícita de que **ni H1 ni Λ prueban el modelo OU** (σ_X > 0): H1 es un desvío fijo en un bin, no una fluctuación; un resultado
  a favor de H1 no es evidencia de energía oscura estocástica, y uno a favor de H0 no excluye σ_X > 0 al nivel de las cotas de JACKKNIFE/PROPUESTA.
  El borrador no afirma lo contrario, pero tampoco lo cierra, y en el registro público hace falta.
- "**Comprobación de la integración (numérica…)**" y "**gaussianidad confirmada**" (también en el mensaje del commit): véase (b) 1 y 3.
- "Potencia 0.500": véase (b) 2.
- No he encontrado frases que digan que DR3 decidirá si la energía oscura fluctúa ni que los límites sean universales.

## 6. Rupturas: qué falta en la lista

1. **Formato de publicación de DR3**: si DESI publica LRG2 solo como D_V (o en α_iso/α_AP, o sin LRG2 separado), r3 no se puede calcular como
   está definido. Declarar "no evaluable" o fijar ahora la conversión permitida.
2. **Correlaciones entre bins en la covarianza de DR3**: DR2 es exactamente diagonal por bloques (144 ceros); el anidamiento, el reajuste y
   V3,k lo suponen. La regla 5 habla de la diagonal; debe incluir "si aparecen bloques fuera de la diagonal, modelo roto o nueva simulación".
3. **Definición operativa de ρ²_k y rango**: cómo se mide (cociente de diagonales por componente; qué hacer si D_M y D_H dan ρ² distintos o la
   correlación D_M–D_H de LRG2 cambia, es decir C3,LRG2 no ∝ C22) y qué pasa si ρ²_LRG2 queda fuera de [0.5, 0.72]: la regla 4 dice "repetir",
   pero debe decir que se repite **con el script registrado sin otro cambio** (commit y semillas) antes de mirar DR3, y con qué ρ² se calibran c0 y c1.
4. **Congelar el código y las densidades**: el registro debe citar el hash del commit de `pronostico_dr3_lrg2.py` y del json, y que Λ se calcula
   una sola vez con ese código. Ahora el borrador no lo dice.
5. **"Los otros seis bins"**: si DR3 añade bins o tracers (más bins de Lyα, ELG separados), hay que fijar ahora si el reajuste usa solo los seis
   homólogos de DR2 o todos los de DR3. Afecta a Σ3 y por tanto a V3,k y a Λ.
6. **Cambios en la cadena propia**: erratas de transcripción de los inputs de DR2 o de los externos descubiertas después (hay un registro de
   erratas en el repositorio); la regla 1 cubre solo lo que DESI republique. Declarar y recalcular con registro.
7. **Fecha límite y qué pasa si DR3 no sale** (o sale después de otros datos que cambien los externos).
8. Sugerencia, no ruptura: comprometer el informe de los tres resultados (H1 / H0 / no concluyente) con los mismos umbrales, sin añadir
   umbrales después.

## (c) Sugerencias

- Registrar las densidades exactas del modelo linealizado (a, MMᵀ por caso) junto a las simuladas, y la comparación de la sección 1 como
  prueba real de las simulaciones.
- Añadir los errores Monte Carlo y una tabla de sensibilidad de potencia a ρ² en continuo (0.4–0.8), que con el modelo afín es inmediata.
- Explicitar el supuesto σ² ∝ tiempo y, si DESI publica tiempos efectivos por tracer, usar ρ² por bin desde el principio.
- Decir en una frase qué sería evidencia a favor del modelo OU y qué no (nada de este pronóstico).

## Veredicto

**Registrable con estos cambios:** (1) errores Monte Carlo en c0, c1, potencia y medianas; (2) corregir la "comprobación de la integración"
y la descripción de la gaussianidad (identidad / cierta por construcción; citar la validación completa como respaldo), y aclarar que V3,k es
fija en los universos linealizados; (3) matizar "evidencia fuerte" y añadir que H1/Λ no prueban el modelo OU ni σ_X > 0; (4) completar las
rupturas con el formato de publicación, las correlaciones entre bins, la definición operativa de ρ²_k con su rango y el script congelado,
los "otros seis bins" y la fecha límite. Ninguna cifra cambia; P18 y P19 se mantienen con margen (potencia 0.655 ± 0.004; medianas 4.76 ± 0.02 y
12.29 ± 0.03).
