# REVISIÓN 3 del pronóstico DR3 de LRG2: comprobación de v3 (6f1628b) y v4 (1439356) (13-11-2026)

Revisor: el mismo de las revisiones 1 (99f879d) y 2 (c7df2d8). Objeto: `PRONOSTICO_DR3_REGISTRO_v3.md`, `PRONOSTICO_DR3_REGISTRO_v4.md`,
`pronostico_dr3_lrg2_v3.py`, `pronostico_dr3_lrg2_v4.py`, sus json, logs y ficheros de datos falsos. Ningún fichero existente se modifica. Nuevos, con
prefijo `revision3_pronostico_`: `revision3_pronostico_comprobacion.py` → `.json` (anidamiento, densidades y umbrales de v3, modelo en z y unidades de v4,
conversión de μ, denominador, tablas) y `revision3_pronostico_romper.py` → `.json` (ataques al modo evaluación de v4 con el script cargado desde su ruta y las
salidas redirigidas a un directorio fuera del repositorio).

## 1. Los dos errores de la revisión 2: resueltos

- **Anidamiento matricial.** Deducción correcta: con R = C3C2⁻¹ y Cov(e) = C3 − C3C2⁻¹C3, Cov(n3, n2) = RC2 = C3 y Cov(n3) = RC2Rᵀ + Cov(e) = C3 (verificado en 166
  pares aleatorios a 5e-14). Cov(e) ⪰ 0 ⟺ C2 − C3 ⪰ 0 (las dos comprobaciones son equivalentes; en los 166 pares coinciden en signo). `nesting` de v3: con
  C3 = 0.6·C2, R − 0.6I = 2e-16 y FeFeᵀ − 0.24·C2 = 4e-17; con una C3 no proporcional (LRG2 0.54/0.66, corr −0.45) RC2 − C3 = 3e-17, FeFeᵀ − Cov(e) = 4e-17 y
  un Monte Carlo de 4·10⁵ ruidos reproduce Cov(n3) = Cov(n3, n2) = C3 al 0.5 % (error MC). Tolerancia de autovalores como se describe. Un caso informativo: si
  DR3 tuviese la misma varianza que DR2 en D_H de LRG2 (ρ² = 1 en una componente), C2 − C3 tiene un autovalor real negativo (−0.005) por el término cruzado y
  el modelo se declara roto; es coherente con "igual o más preciso en todas las direcciones", pero conviene saber que "igual" en una componente ya rompe.
- **Salida principal única sin SN**: `result_principal` y `result_secondary_Pantheon+` en el código, y así en el texto y en las pruebas falsas.

## 2. Los menores de la revisión 2: tratados

Umbrales con 10⁷ muestras exactas y error de 10 bloques; la diferencia 1.081/1.024 explicada (semillas distintas de v2) y resuelta (1.0274 ± 0.0013 frente a
1.0285 ± 0.0009 con otra semilla). Mi muestreo independiente (8·10⁶ por caso, 2·10⁷ en el principal, semillas distintas) frente a v3, en σ combinada:

| combinación | ρ² | c0 mío / v3 | c1 mío / v3 | potencia mía / v3 | desviación máx. (σ comb.) |
|---|---|---|---|---|---|
| sinSN | 0.6 | 1.0269 ± 0.0008 / 1.0274 ± 0.0013 | −1.5918 ± 0.0007 / −1.5944 ± 0.0012 | 0.6581 ± 0.0002 / 0.6575 ± 0.0002 | 2.0 |
| sinSN | 0.5 | 0.7579 ± 0.0022 / 0.7559 ± 0.0009 | −1.6101 ± 0.0020 / −1.6115 ± 0.0020 | 0.7674 ± 0.0003 / 0.7675 ± 0.0001 | 0.8 |
| sinSN | 0.72 | 1.1942 ± 0.0011 / 1.1971 ± 0.0012 | −1.4696 ± 0.0006 / −1.4707 ± 0.0008 | 0.5006 ± 0.0004 / 0.4997 ± 0.0003 | 1.8 |
| Pantheon+ | 0.6 | 0.9995 ± 0.0012 / 0.9999 ± 0.0014 | −1.5838 ± 0.0009 / −1.5861 ± 0.0006 | 0.6803 ± 0.0003 / 0.6801 ± 0.0002 | 2.2 |
| Pantheon+ | 0.5 | 0.7032 ± 0.0021 / 0.7061 ± 0.0014 | −1.5692 ± 0.0012 / −1.5726 ± 0.0014 | 0.7881 ± 0.0002 / 0.7878 ± 0.0002 | 1.8 |
| Pantheon+ | 0.72 | 1.1943 ± 0.0014 / 1.1968 ± 0.0006 | −1.4856 ± 0.0009 / −1.4830 ± 0.0010 | 0.5211 ± 0.0003 / 0.5208 ± 0.0002 | 2.0 |

Las 18 comparaciones quedan a ≤ 2.2 σ combinadas: los errores MC de v3 son correctos. Regla b con tolerancia 1e-6 (probada en v3: 1e-3 roto, 1e-9 ok). Texto:
controles del script como coherencia y orden público de commits como compromiso, nota sobre `git_head_at_run`, redondeos exactos (0.0072 σ, 1.64 %). z_eff ± 0.01
igual en texto y código **salvo en el borde** (véase (a) 1). Las densidades de v3 coinciden con mi `affine_map` de la revisión 1 a ≤ 4e-15 (12 casos).

## 3. Conversión de μ en v4 (foco)

**Relación exacta.** Bajo el supuesto de H1 fraccional, dato_L = m_L(p, z)(1 + μ_f) + ruido en los dos releases, con μ_f el mismo en z_DR2 y z_DR3. El desvío
absoluto de DR2 es μ2 = m_L(p_true, z_DR2) ⊙ μ_f y el de DR3 μ3 = m_L(p_true, z_DR3) ⊙ μ_f, luego, componente a componente (D_M, D_H):

  μ3,i = μ2,i · ρ_i(p_true),   ρ_i(p) = m_i(p, z_DR3) / m_i(p, z_DR2).

v4 usa ρ_i(p2) (`ratio` en `case_ops`), de modo que r3 siga siendo afín en gaussianas. **Término omitido:** Δ_i = μ2,i · [ρ_i(p_true) − ρ_i(p2)].

**Medida sobre el posterior de (p, μ)** (10⁵ muestras: p ~ N(p2, Σ2), μ2 | p ~ N(datos_L − m_L(p, z2), C_LRG2), ρ con el modelo exacto, sin SN, ρ² = 0.6):

| Δz de LRG2 | ρ(p2) (D_M, D_H) | máx \|ρ(p_true) − ρ(p2)\| | \|Δ\| mediana (abs; D_M, D_H) | \|Δ\| máx (abs) | \|Δ\|/σ mediana | \|Δ\|/σ máx |
|---|---|---|---|---|---|---|
| −0.01 | 0.98852, 1.00606 | 1.6e-4 | 2.1e-6, 1.2e-5 | 5.5e-5, 1.6e-4 | 1.8e-5, 5.8e-5 | 4.8e-4, 7.5e-4 |
| +0.01 | 1.01141, 0.99398 | 1.5e-4 | 2.2e-6, 1.2e-5 | 5.6e-5, 1.6e-4 | 1.8e-5, 5.7e-5 | 4.7e-4, 7.5e-4 |

(σ = desviación típica de r3 en unidades absolutas, 0.09–0.21 en D/r_d.) El "~1e-5" del texto es la mediana; el máximo en la cola del posterior es 1.6e-4
(≤ 7.5e-4 σ). En el caso principal registrado (z_DR3 = z_DR2) el término es idénticamente cero.

**Efecto en c0, c1 y potencia.** c0 no depende de μ. Con 4·10⁶ universos pareados (mismas z_p, z_e, z_μ) bajo H1, linealizado frente a expresión completa (ρ(p_true)
por universo; la expansión de primer orden en p que uso difiere del ρ exacto en ≤ 5.7e-7, 100 veces menos que el propio término), con las densidades y el c0 de v4
en cada Δz:

| Δz | c1 linealizado | c1 completo | Δc1 (en σ_MC) | potencia linealizada | potencia completa | Δpotencia (en σ_MC) | máx \|ΔΛ\| por universo |
|---|---|---|---|---|---|---|---|
| −0.01 | −1.5927 | −1.5927 | −2e-5 (−0.01) | 0.6571 | 0.6571 | +5e-6 (+0.02) | 6.2e-3 |
| +0.01 | −1.5939 | −1.5939 | +1e-5 (+0.00) | 0.6579 | 0.6579 | −3e-6 (−0.01) | 6.0e-3 |

Sustituir la linealización por la expresión completa no cambia c1 ni la potencia de forma apreciable (≪ error MC). La conversión está bien derivada e
implementada. Nota informativa: con el supuesto alternativo (μ absoluto, no fraccional) μ3 diferiría en (ρ − 1)·μ2 ≈ (0.002, 0.004) en D/r_d (0.02–0.03 σ) a
Δz = 0.01; el supuesto fraccional está declarado y es razonable.

## 4. Unidades

Generación en unidades absolutas (modelo, ruido heredado R·n2, e, reajuste), y solo el residuo final de LRG2 se divide por m(p2, z_DR3) (`s3`); V3,k se reescala
con s3 ⊗ s3; las densidades (a, Σ) y los umbrales son fraccionales; en `evaluar`, r3 = r_abs / m(p3, z_DR3) y V = V_abs / (m3 ⊗ m3). Comprobado: con Δz = 0,
a_frac · m(p2) − a_v3 ≤ 9e-15 y Σ_frac ⊙ (m ⊗ m) − Σ_v3 ≤ 3e-16, V3k ≤ 8e-16 (las dos combinaciones); m(p2) de LRG2 = (17.5603, 20.0907) como dice el texto.
Λ, c0, c1 y d² son invariantes bajo el reescalado (coinciden con v3 a ≤ 8e-14). La única mezcla es el denominador m(p3) en la evaluación frente a m(p2) en la
calibración: medida en (b) 3.

## 5. Redshifts condicionantes y de predicción

En `gen_r3` de v4: m_true2 = m(p2, z_DR2) + J(z_DR2)·dp y n2 = datos_DR2 − m_true2 [− E₂μ2] (condicionamiento en z_DR2); m_true3 = m(p2, z_DR3) + J(z_DR3)·dp y
d3 = m_true3 [+ E₂μ3] + R·n2 + e (predicción en z_DR3); reajuste con G3 de la Fisher en z_DR3 y C3, residuo frente a m(p2, z_DR3) + J(z_DR3)·dp3. Correcto. Modelo en z
arbitrario: en los nodos de DR2 reproduce la malla a 2e-16; fuera de los nodos, D_M interpolada frente a la integral directa de 1/E con el mismo E: 2–4e-7 relativo
(el texto dice < 1e-6: cierto); D_H cerrado exacto. `fit_cz` en z_DR2 reproduce `fit_c` a 6e-10 (p) y 7e-12 (χ²). Las 75 filas numéricas de los dos registros
coinciden con sus json; las cifras del texto (prueba z0712: c0/c1/potencia, r3, d², Λ; v3 falso; P24) también.

## 6. Modo evaluación de v4: orden de lectura y ataques (`revision3_pronostico_romper.json`)

`calibrar` lee solo `labels` (z_eff) y `cov` (probado con un fichero que además traía `mean`: `keys_read = [labels, cov]`), empareja por tipo, trazador y z_eff,
aplica las reglas, calibra con 10⁷ muestras y escribe la calibración con sha256. `evaluar` exige calibración, sha256, misma covarianza y script, y se niega a
repetir el TAG. Salidas probadas: LRG duplicado dentro de la tolerancia → no evaluable; trazador renombrado → no evaluable; LRG2 en 0.715 → ok con desplazamiento
+0.009 registrado; LRG1 en 0.519 → ok con desplazamiento registrado en LRG1; evaluar sin calibrar, repetir, otra covarianza, sha registrado incorrecto → se niega.

Vías que sí pasan:

| intento | resultado |
|---|---|
| **LRG2 en z = 0.716 o 0.696 (exactamente ± 0.01)** | **"no evaluable"**: 0.716 − 0.706 = 0.010000000000000009 > 0.01 en coma flotante; lo mismo 0.52 − 0.51. El texto dice ± 0.01 inclusive y P24 se calculó en ± 0.01 |
| sin fichero `pronostico_dr3_lrg2_v4.sha256` en el directorio | `self_check` devuelve `match: None` y `evaluar` continúa (`assert … is not False`) |
| `calibrar` con un sha registrado incorrecto | calibra y solo anota `match: False` (luego `evaluar` sí se niega) |
| valores etiquetados z0.706 con calibración hecha con covarianza z0.712 | **se ejecuta**: solo se comparan los índices del emparejamiento, no las etiquetas |
| valores con las filas en otro orden, o con un bin nuevo delante | `AssertionError` ("las etiquetas … no coinciden"): rechazo seguro pero no es una salida declarada; exige el mismo orden en los dos ficheros |
| reescribir la calibración y su `.sha256` | se ejecuta con el c0 alterado (agujero conocido; el compromiso son los commits públicos) |

## (a) Errores que impiden registrar

1. **Frontera de la regla f en coma flotante.** `map_labels` usa `abs(z3 − z2) <= 0.01` y, con z a tres decimales, 0.716 − 0.706 y 0.706 − 0.696 valen
   0.010000000000000009: un z_eff publicado exactamente en el borde (0.696 o 0.716, perfectamente posible con tres decimales) se declara **"no evaluable"** contra
   el texto registrado ("0.706 ± 0.01", inclusive) y contra P24 (calculada en ± 0.01). Efecto medido: z = 0.715 → ok; 0.716 → no evaluable; 0.696 → no evaluable;
   LRG1 0.520 → no evaluable. Corrección de una línea (comparar z redondeados a 3 decimales o añadir 1e-9 a la tolerancia), en v4 y en la regla de v3.

## (b) Errores menores

1. **Autocontrol del script**: `assert sc["match"] is not False` deja pasar `None` (fichero `.sha256` ausente o renombrado), y `calibrar` no se niega con
   `match: False`. Poner `is True` en `evaluar` y negarse también en `calibrar`. (Mismo código en v3.)
2. **Etiquetas del fichero de valores**: solo se comparan los índices del emparejamiento; un fichero de valores con z0.706 pasa con una calibración a z0.712
   (probado). Exigir etiquetas idénticas a las de la covarianza (o emparejar por etiqueta, lo que también cubriría el reordenamiento y los bins nuevos, que hoy
   terminan en `AssertionError` y no en una salida declarada). Documentar que los dos ficheros deben tener las mismas etiquetas en el mismo orden.
3. **Denominador m(p3, z3) en la evaluación frente a m(p2, z3) en las densidades.** Medido con Δz = 0.006 en 4·10⁵ universos: m3/m2 − 1 entre −0.4 % y +0.3 %;
   desplazamiento de r3 típico 1.5e-3 σ (como dice el texto), pero 1e-2 σ al 99 % y 3.6e-2 σ máximo; ΔΛ típico 0.004–0.009, 0.02–0.06 al 99 %, 0.09–0.19 máximo.
   En la prueba falsa z0712, Λ = 2.484 con m(p3) frente a 2.469 con m(p2). Efecto en la falsa alarma ≈ f(c0)·E|ΔΛ| ~ 3e-4: no cambia nada en la práctica, pero la
   calibración ya guarda `m2_bao3_LRG2`, y usarlo como denominador en `evaluar` haría Λ exactamente el mismo que el de v3 en unidades absolutas con el modelo en
   z_DR3, sin esta diferencia. El texto debería decir "≤ 0.04 σ (típico 1.5e-3 σ)".
4. Texto "el término de segundo orden dp·μ es ~1e-5": es la mediana; el máximo en el posterior es 1.6e-4 en D/r_d (≤ 7.5e-4 σ). Decir "≤ 2e-4 (≤ 1e-3 σ)".
5. Informativo (sección 1): una C3 con la misma varianza que DR2 en una componente de un bin rompe C2 − C3 ⪰ 0 por el término cruzado; conviene anotarlo junto a la
   regla d para que no sorprenda.

## Veredicto

**No registrable hasta corregir:** la frontera ± 0.01 de `map_labels` (coma flotante; (a) 1), en v4 y en v3. Recomiendo corregir a la vez (b) 1 y (b) 2, que son
dos líneas, y valorar (b) 3 (denominador de la calibración en `evaluar`). La conversión de μ, las unidades, la separación de redshifts, el anidamiento matricial, los
umbrales de alta precisión y las tablas están bien; ninguna cifra de P18, P19, P23 o P24 cambia.
