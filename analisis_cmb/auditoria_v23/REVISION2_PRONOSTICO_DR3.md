# REVISIÓN 2 del pronóstico DR3 de LRG2: comprobación de la versión 2 (13-11-2026)

Revisor: el mismo de `REVISION_PRONOSTICO_DR3.md` (commit 99f879d). Objeto: `PRONOSTICO_DR3_REGISTRO.md`, `pronostico_dr3_lrg2_v2.py`,
`pronostico_dr3_lrg2_v2.json`, sus logs y los ficheros de la prueba con datos falsos (commit 3cdb79b). Ningún fichero existente se modifica.
Ficheros nuevos, prefijo `revision2_pronostico_`:

- `revision2_pronostico_comprobacion.py` → `.json`: densidades exactas de v2 frente a mi `affine_map` de la primera revisión (12 casos y 42 puntos
  del barrido); reproducción de las decisiones registradas con mi simulador; c0, c1 y potencia de alta precisión muestreando 4·10⁶ universos
  directamente de N(a, Σ); cotejo numérico de las 69 filas de las tablas del registro con el json; efecto de la tolerancia de la regla c
  (anidamiento matricial exacto); camino ρ² por bin.
- `revision2_pronostico_romper.py` → `.json`: intentos de romper el modo evaluación, con el script cargado desde su ruta original pero con sus
  salidas redirigidas a un directorio de trabajo fuera del repositorio (nada se escribe ni se borra en `auditoria_v23`).

## 1. Cambios pedidos en la primera revisión: aplicados y correctos

Todos están en el texto y, donde toca, en el código: errores MC en c0, c1, potencia, P(Λ < c1 | H0), no concluyentes, medianas y P(|Λ| > 3)
(20 reejecuciones de B con semillas registradas); la "comprobación de la integración" llamada identidad algebraica; gaussianidad cierta por construcción
con la validación completa como respaldo (máx |Δr3| 7e-4 = 0.0084 σ de la componente más estrecha, √0.00702 = 0.084: "≤ 0.009 σ" correcto);
V3,k fija en los universos linealizados dicho en el mismo párrafo; prior plano en p explícito; Cov(r3, r2) etiquetada "sin reajuste"; σ² ∝ 1/tiempo
declarado como supuesto con la cita correcta (§4.4: 1/9, 3/9, 5/9); "evidencia fuerte de una hipótesis frente a la otra, dentro de este par y bajo
estos supuestos" con la asimetría 0.118 / 0.374; "ni H1 ni Λ prueban el modelo OU ni σ_X > 0" (dos veces); rupturas a–g con regla cada una; fecha
límite 31-12-2028 con salida "no evaluado"; compromiso de cinco salidas. La fórmula de conversión de la regla a es correcta (α_⊥ = α_iso α_AP^{−1/3},
α_∥ = α_iso α_AP^{2/3}).

## 2. Texto frente a código

- **Tablas**: las 69 filas numéricas del registro (decisión, comparación con el borrador, exactas frente a simuladas, d², Λ, barrido 21 × 2) coinciden
  con `pronostico_dr3_lrg2_v2.json` a medio dígito. P18 (0.65574 ± 0.0031), P19 (4.7621 ± 0.014; 12.2913 ± 0.031) y las cifras citadas de la prueba
  falsa (c0/c1 +1.013/−1.599 y +0.978/−1.573; r3 = (−0.182, −0.620); d² 10.38; Λ +1.127 y +1.169) coinciden. Los sha256 del script (add948…) y del json
  (4c831e…) coinciden con los ficheros. "Medias ≤ 0.007 σ" es 0.0072 y "varianzas ≤ 1.6 %" es 1.64 % (redondeo).
- **Densidades exactas**: `exact_density` de v2 hace lo que dice (a = r3(0), columnas r3(e_k) − a, Σ = MMᵀ, no linealidad ≤ 8e-15) y Λ usa
  `multivariate_normal` con esas (a, Σ) tanto en el registro como en `evaluar`.
- **c0, c1, potencia**: percentiles 95 (H0) y 5 (H1) de Λ sobre la mitad B con las semillas del borrador; potencia = fracción de Λ > c0 bajo H1. Como dice.
- **Barrido**: cada ρ² con su propia mitad B (semillas 100 + k): por eso el punto 0.60 del barrido (c0 1.081, potencia 0.652) no coincide con el caso
  principal (1.024, 0.656); es ruido MC, véase (b) 1.
- **Modo evaluación**: `calibrar` abre solo el fichero de covarianza; aplica las reglas b y c; recalibra con las semillas fijadas 500; escribe la
  calibración con sha256 de covarianza, script y commit y su propio sha256. `evaluar` exige la calibración, su sha256, la misma covarianza y el mismo
  script; reajusta con el optimizador original y C3 (cálculo completo); calcula Λ y se niega si ya existe el resultado de ese TAG. Como dice, con las
  matizaciones de la sección 4.
- Las etiquetas de DR2 llevan z_eff (`DM_LRG_z0.706`), así que `read_dr3_cov` exige la cadena exacta: en el código la regla 2 es "z_eff idéntico", en
  el texto "0.706 ± 0.01". Véase (b) 5.

## 3. Densidades exactas frente a mi cálculo

En los 12 casos principales y los 42 puntos del barrido, |a_v2 − a_mía| ≤ 4e-15 y |Σ_v2 − Σ_mía| ≤ 5e-17. V3,k idéntica. Con esas densidades y las
semillas del borrador, mi simulador reproduce c0, c1, potencia y medianas registradas a ≤ 3e-14. El camino "ρ² por bin" de v2 coincide con una
generación propia con R = diag(ρ²_fila) y Cov(e) = C3 − C3C⁻¹C3 (diferencia 0 / 3e-18).

**Valores de alta precisión** (4·10⁶ universos por hipótesis muestreados de N(a, Σ); error ≈ ±0.002 en c0, c1 y ±0.0003 en potencia), frente a los
registrados ± MC:

| combinación | ρ² | c0 registrado / alta precisión | c1 | potencia | desviación máx. en σ_MC |
|---|---|---|---|---|---|
| sinSN | 0.6 | 1.024 ± 0.017 / 1.028 | −1.604 ± 0.015 / −1.590 | 0.6557 ± 0.0031 / 0.6580 | 0.9 |
| sinSN | 0.5 | 0.777 ± 0.018 / 0.758 | −1.638 ± 0.019 / −1.607 | 0.7669 ± 0.0027 / 0.7675 | 1.7 |
| sinSN | 0.72 | 1.184 ± 0.016 / 1.198 | −1.476 ± 0.015 / −1.469 | 0.4998 ± 0.0048 / 0.5003 | 0.8 |
| Pantheon+ | 0.6 | 0.993 ± 0.021 / 1.001 | −1.579 ± 0.020 / −1.582 | 0.6850 ± 0.0037 / 0.6805 | 1.2 |
| Pantheon+ | 0.5 | 0.667 ± 0.020 / 0.706 | −1.560 ± 0.024 / −1.567 | 0.7939 ± 0.0027 / 0.7882 | 2.1 |
| Pantheon+ | 0.72 | 1.195 ± 0.015 / 1.196 | −1.478 ± 0.013 / −1.483 | 0.5211 ± 0.0033 / 0.5210 | 0.4 |

Los 24 valores registrados están a ≤ 2.1 σ_MC del valor de alta precisión: los errores MC declarados son correctos. Barrido: desviación típica de
(registrado − alta precisión) c0 0.022 / 0.021, c1 0.014 / 0.016, potencia 0.003 / 0.0045 (sin SN / Pantheon+), como dice el texto ("≈ ±0.02, ≈ ±0.004").

**Punto 5 del encargo (0.500 ± 0.005 con ρ² = 0.72).** Correcto y bien explicado: con las densidades exactas la mitad B del borrador da 0.4998 y el error MC
(20 reejecuciones de B con densidades fijas) 0.0048; mi 0.501 ± 0.004 era la media de 20 reejecuciones completas con densidades reestimadas en cada una.
El valor de alta precisión es 0.5003. Las tres cifras son coherentes.

## 4. Intentos de romper el modo evaluación (`revision2_pronostico_romper.json`)

| intento | resultado |
|---|---|
| evaluar sin calibrar | se niega ("primero calibrar") |
| evaluar dos veces con el mismo TAG (mismos u otros valores) | se niega ("Λ se calcula una sola vez") |
| evaluar con otra covarianza que la calibrada | se niega |
| cambiar c0 en la calibración sin tocar el `.sha256` | se niega |
| cambiar c0 en la calibración **y reescribir el `.sha256`** | **se ejecuta** y el veredicto cambia (sin SN: "favorece H1" → "no concluyente" con c0 = 5) |
| mismos datos con otro TAG (calibrar + evaluar de nuevo) | **se ejecuta**: el "una sola vez" es por etiqueta |
| **script modificado** (SEED_EVAL 500 → 501): calibrar y evaluar con él | **se ejecuta** sin aviso; registra otro sha256 (288d33…) pero nadie lo compara con el registrado (add948…); el script original sí se niega a usar esa calibración |
| etiquetas permutadas | se niega ("no evaluable sin conversión declarada") |
| un 1e-15 fuera de los bloques | "modelo roto" (umbral `!= 0.0`) |
| componentes de LRG2 con 0.9/1.1 (spread exactamente 0.20) | pasa (el límite es "> 0.20") |
| borrar el resultado y reevaluar | no probado (prohibido borrar); el script no puede impedirlo |

Conclusión: los controles del script son de consistencia interna (protegen de errores y de accidentes), no un compromiso. Lo que impide decidir
después de ver los datos es el commit + push públicos de la calibración y de su `.sha256` **antes** de descargar los valores, y el commit público del
resultado; el registro debe decirlo así (véase (b) 4). Por otro lado, con el script registrado y la semilla 500 el procedimiento es determinista, de modo
que un cambio de c0 solo es posible cambiando ficheros que quedan en el historial público. También: c0 cambia de 1.013 (semilla 500) a 1.059 (semilla 501):
la semilla fija pesa ±0.02–0.05 en el umbral, véase (b) 2.

## (a) Errores que impiden registrar

1. **Regla c: dentro de su propia tolerancia, la falsa alarma registrada (5 %) no es la del procedimiento.** `calibrar` construye las densidades con
   Cov(n3) = ρ²_k·C_kk (ρ²_k = media de las dos componentes), pero la regla admite componentes que difieran hasta un 20 % y una correlación D_M–D_H de
   LRG2 que cambie hasta 0.1. En ese régimen la covarianza real de DR3 no es ρ²_k·C_kk, y `evaluar` reajusta con la C3 real. Medido con el anidamiento
   matricial exacto (n3 = C3C⁻¹n2 + e, Cov(e) = C3 − C3C⁻¹C3, que es Cov(n3) = Cov(n3, n2) = C3 y se reduce al actual cuando C3 ∝ C; sin SN, 2·10⁶
   universos, calibración del script con semilla 500):

   | C3 de LRG2 (ρ²_DM / ρ²_DH / corr; DR2: −0.395) | estado regla c | falsa alarma real P(Λ > c0 \| H0) | potencia real | P(Λ < c1 \| H1) real |
   |---|---|---|---|---|
   | 0.54 / 0.66 / −0.49 | ok | **0.030** | 0.649 | **0.059** |
   | 0.66 / 0.54 / −0.30 | ok | **0.081** | 0.672 | 0.036 |
   | 0.60 / 0.60 / −0.49 (solo correlación) | ok | **0.022** | 0.644 | **0.067** |
   | 0.54 / 0.66 / −0.395 (solo spread) | ok | **0.065** | 0.667 | 0.043 |
   | 0.585 / 0.615 / −0.425 (spread 5 %, Δcorr 0.03) | ok | 0.042 | 0.656 | 0.053 |
   | 0.615 / 0.585 / −0.365 (spread 5 %, Δcorr 0.03) | ok | 0.060 | 0.664 | 0.045 |

   Con la tolerancia registrada, el "5 %" puede ser cualquier cosa entre 2 % y 8 %, y el "error del 5 % bajo H1" entre 4 % y 7 %. Corrección mínima
   (una de las dos, a elección del autor, antes de registrar): **(i)** calibrar con la covarianza publicada C3 directamente mediante el anidamiento
   matricial (cambio local en `case_ops`/`gen_r3`: R = C3C⁻¹ por bloques y Cov(e) = C3 − C3C⁻¹C3, exigiendo C − C3 semidefinida positiva), con lo que la
   regla c queda solo como control de rango de ρ²_LRG2 y el 5 % vuelve a ser exacto en el modelo; o **(ii)** estrechar la tolerancia (p. ej. spread ≤ 5 % y
   |Δcorr| ≤ 0.03, que da 4.2–6.0 %) y declarar en el registro ese intervalo real de falsa alarma. En cualquier caso, ninguna cifra de P18/P19 cambia.
2. **Combinación principal sin declarar para la salida única.** El compromiso promete "una de estas cinco salidas", pero `evaluar` devuelve un veredicto
   por combinación y pueden discrepar (en mis pruebas: sin SN "no concluyente" y Pantheon+ "favorece H1"). Hay que registrar que la salida es la de
   **sin SN** (la combinación de P18 y P19) y que Pantheon+ se informa al lado como secundaria, sin elegir después.

## (b) Errores menores

1. El punto 0.60 del barrido (c0 1.081, potencia 0.652) y el caso principal (1.024, 0.656) difieren solo por la mitad B; el valor de alta precisión es
   c0 1.028, potencia 0.658. Conviene decirlo junto a la tabla del barrido para que nadie lea dos valores de c0 para el mismo ρ².
2. **Calibración con 50000 universos y semilla fija**: c0 y c1 llevan ±0.02 de ruido MC (y cambian 0.05 al cambiar la semilla), que se traslada a la
   decisión cuando Λ cae cerca del umbral (en la prueba falsa, Λ = 1.127 frente a c0 = 1.013 con semilla 500 y 1.059 con 501). Como Λ es función solo de
   r3 y las dos densidades son gaussianas exactas, `calibrar` puede muestrear 4·10⁶ universos directamente de N(a, Σ) (un segundo) y fijar c0, c1 a ±0.002,
   con lo que la semilla deja de importar. No bloquea (la semilla está congelada), pero es gratis y elimina (b) 1 y la dependencia de semilla.
3. **Regla b con umbral `!= 0.0`**: un 1e-15 fuera de los bloques (redondeo de un fichero de texto) declara "modelo roto" (probado). DR2 es exactamente
   0.0, pero la regla sería más robusta con una tolerancia de correlación (p. ej. |ρ_ij| ≤ 1e-6) declarada ahora.
4. **El texto sobrevalora los controles del script**: "Λ se calcula una sola vez con ese código (modo evaluar, que comprueba los sha256 y no sobrescribe)".
   Como muestra la sección 4, los sha256 y el "una sola vez" se saltan reescribiendo el `.sha256`, con otra etiqueta o con un script modificado. Debe decir
   que el compromiso es el commit + push de la calibración antes de descargar los valores y el commit del resultado, y, si se quiere, que `evaluar`
   compare el sha256 del script con la constante registrada y exija árbol de git limpio y calibración ya confirmada.
5. Regla 2 en el texto ("z_eff = 0.706 ± 0.01") frente al código (etiqueta idéntica, `DM_LRG_z0.706`): si DR3 publica z_eff = 0.71 el script dice
   "no evaluable" mientras el texto lo admitiría. Elegir uno; el del código es más estricto y más simple.
6. Redondeos "≤ 0.007 σ" (0.0072) y "≤ 1.6 %" (1.64 %).
7. `git_head_at_run` del json y de las calibraciones falsas es 99f879d (el commit anterior al registro), como es inevitable; el registro lo resuelve citando
   los sha256. Sin consecuencia.

## Veredicto

**No registrable hasta corregir:** (1) la regla c (calibrar con la C3 publicada mediante el anidamiento matricial, o estrechar la tolerancia y declarar
el intervalo real de falsa alarma); (2) la combinación principal (sin SN) de la salida única. Son dos cambios pequeños que no tocan P18, P19, las
densidades exactas ni los umbrales de los casos registrados. Recomiendo además, aunque no bloquea, (b) 2 (calibración con 4·10⁶ universos exactos) y
(b) 4 (describir el compromiso como commit público, no como control del script).
