# PRONÓSTICO DR3 para el residuo de LRG2 — REGISTRO, versión 4 (13-11-2026)

**PENDIENTE DE COMPROBACIÓN DEL REVISOR — todavía sin tag.**

Parte de `PRONOSTICO_DR3_REGISTRO_v3.md` (commit 6f1628b) con una corrección limitada al tratamiento del redshift, marcada **(cambio por z_eff)**.
Todo lo demás de v3 sigue vigente sin cambio: P18, P19, P23, H0, H1, el anidamiento matricial con sus comprobaciones y tolerancias, la regla de decisión
(c0 = percentil 95 de Λ | H0, c1 = percentil 5 de Λ | H1, tres salidas), los umbrales de alta precisión, la salida principal sin SN, las reglas b–g, el
procedimiento público y las cautelas. Código: `pronostico_dr3_lrg2_v4.py` (importa v3 sin modificarlo; sha256 `c38e623539a6e7cc64aba9baacc30c14366ca797a335ea519e3b9fb059847e0d`, también en
`pronostico_dr3_lrg2_v4.sha256`); resultados `pronostico_dr3_lrg2_v4.json` (sha256 `78fc549ec9197f04eaf5b320c75970e0410e6b6000305531d4fdcf41dcaa5c8a`). El campo `git_head_at_run` es el del
commit anterior (6f1628b); el registro lo fija el tag.

## Predicción nueva (de Claude, escrita antes de ejecutar; copiada literal)

- P24. Desplazar el z_eff de LRG2 en ±0.01 (recalculando el modelo y el jacobiano en el nuevo z) cambia c0 y c1 menos de 0.05 y la potencia del caso
  principal menos de 0.02.

## (cambio por z_eff) El problema que corrige esta versión

La regla f de v3 admitía un z_eff de DR3 hasta ±0.01 del de DR2, pero la predicción de DR3 se hacía en el z_eff de DR2. Un desplazamiento de 0.01 en z
mueve D_M/r_d del orden del 1 %, comparable al error del bin, y habría entrado en el residuo como si fuera señal. Ahora el desplazamiento lo absorbe la
predicción: el modelo ΛCDM y su jacobiano se evalúan en el z_eff que publique DR3 para cada bin.

## (cambio por z_eff) Redshifts condicionantes (DR2) y redshifts de predicción (DR3)

- **Los datos condicionantes de DR2 conservan sus z_eff.** El ajuste completo p_full, el ajuste sin LRG2 p2, sus Fisher, el ruido de DR2
  n2 = datos_DR2 − modelo(p_true, z_DR2) y el residuo observado de DR2 en LRG2 (r2, en z = 0.706) no cambian. No se sustituyen los redshifts de los
  datos DR2 por los de DR3.
- **Las predicciones de DR3 se evalúan en z_DR3.** Modelo en un z arbitrario con la misma historia de la malla: D_M por interpolación lineal en
  u = ln(1+z) (error < 1e-6 relativo, como para las SN), D_H = 1/E(z) en forma cerrada, D_V = (z·D_M²·D_H)^(1/3); en los nodos de DR2 reproduce el modelo
  de la malla a 2e-16 relativo. Jacobiano de las filas BAO por diferencias centradas del mismo modelo en z_DR3 (ε = 1e-4, 1e-5, 1e-3); filas de θ*, ω_m y
  SN sin cambio. El ajustador con z como argumento reproduce el de v3 en z = z_DR2 a 6e-10 en los parámetros y 7e-12 en χ².
- **Generación simulada (linealizada en p2):** n2 = datos_DR2 − modelo(p_true, z_DR2) [− E₂μ2 bajo H1]; **datos DR3 = modelo(p_true, z_DR3)
  [+ E₂μ3 bajo H1] + R·n2 + e**, con R = C3·C2⁻¹ y Cov(e) = C3 − C3·C2⁻¹·C3 como en v3. El reajuste sin LRG2 en DR3 (Fisher con C3), V3,k y las densidades
  exactas se calculan con el modelo y el jacobiano en z_DR3.
- **(supuesto nuevo) H1 con μ fraccional.** El desvío persistente μ es **fraccional** y el mismo en z_eff DR2 y z_eff DR3 dentro de ±0.01:
  datos_LRG2 = modelo(z)·(1 + μ_f) + ruido en los dos releases. Con μ2 = modelo(z_DR2)·μ_f el desvío absoluto de DR2 (su posterior con prior plano es
  N(datos_LRG2 − modelo_LRG2(p_true, z_DR2), C_LRG2), como en v3), el de DR3 es μ3 = μ2 · m(p2, z_DR3)_LRG2 / m(p2, z_DR2)_LRG2, con el cociente
  evaluado en p2 (linealización; el término de segundo orden dp·μ es ~1e-5). Así r3 sigue siendo afín en gaussianas y las densidades exactas valen. El
  residuo de DR2 sigue en z = 0.706.
- **Residuos fraccionales.** r3 = (dato − modelo(z_eff DR3)) / modelo(z_eff DR3): en las densidades linealizadas el denominador es m(p2, z_DR3); en la
  evaluación real, m(p3, z_DR3) (diferencia de segundo orden, ~1e-3 σ). V3,k se reescala con el mismo factor. Λ, c0, c1, la potencia y d² son invariantes
  bajo este reescalado, de modo que con z_DR3 = z_DR2 v4 reproduce v3 exactamente salvo redondeo.
- **Regla de no evaluable:** si |z_eff DR3 − z_eff DR2| > 0.01 en LRG2 (trazador LRG sin homólogo en 0.706 ± 0.01) → **no evaluable** (la regla f de v3,
  igual en el código). En los otros seis bins el desplazamiento (dentro de ±0.01) solo cambia dónde se evalúa la predicción de DR3; el script lo registra
  por fila (`z_shifts_per_row`).

## (cambio por z_eff) Comprobación v4 = v3 con desplazamiento 0

Con los z_eff de DR2 en los seis casos (sin SN y Pantheon+; ρ² = 0.6, 0.5, 0.72), mismas semillas que v3:

| combinación | ρ² | máx \|Δa\| (unidades absolutas) | máx \|ΔΣ\| | Δc0 / Δc1 / Δpotencia | Δmediana d² H0 / H1 |
|---|---|---|---|---|---|
| sinSN | 0.6 | 7.2e-15 | 2.3e-16 | -2.4e-14 / 4.4e-14 / 0.0e+00 | -3.1e-14 / 4.8e-13 |
| sinSN | 0.5 | 6.8e-15 | 2.8e-16 | -6.6e-14 / -2.3e-14 / 0.0e+00 | -7.5e-14 / 5.8e-13 |
| sinSN | 0.72 | 7.5e-15 | 1.6e-16 | -4.3e-14 / -1.0e-14 / 0.0e+00 | -3.5e-14 / 3.1e-13 |
| Pantheon+ | 0.6 | 8.9e-15 | 2.0e-16 | 8.0e-14 / 1.3e-15 / 0.0e+00 | 9.9e-14 / -3.9e-13 |
| Pantheon+ | 0.5 | 8.4e-15 | 2.7e-16 | 6.6e-14 / 1.6e-14 / 0.0e+00 | 1.4e-13 / -5.6e-13 |
| Pantheon+ | 0.72 | 9.4e-15 | 1.4e-16 | -1.3e-15 / -2.4e-14 / 0.0e+00 | 1.0e-13 / -3.1e-13 |

v4 reproduce v3 a nivel de redondeo (≤ 1e-14 en a y Σ absolutas, ≤ 8e-14 en c0 y c1, 0 en potencia). Las densidades registradas en v4 están en unidades
fraccionales (para el caso principal, modelo de LRG2 en p2: D_M/r_d = 17.5603, D_H/r_d = 20.0907); a y Σ absolutas se recuperan multiplicando por esos valores.

## (cambio por z_eff) P24: z_eff de LRG2 desplazado (caso principal: sin SN, ρ² = 0.6; resto de bins sin desplazar; 10⁷ muestras exactas, semillas [20261114, 0, 600 + i, i_hip, 7])

| Δz_eff de LRG2 (z) | c0 | c1 | potencia | Δc0 | Δc1 | Δpotencia |
|---|---|---|---|---|---|---|
| 0 (0.706) | 1.0274 ± 0.0013 | -1.5944 ± 0.0012 | 0.6575 ± 0.0002 | — | — | — |
| -0.01 (0.696) | 1.0275 ± 0.0017 | -1.5923 ± 0.0007 | 0.6574 ± 0.0003 | +0.0001 | +0.0021 | -0.0002 |
| -0.005 (0.701) | 1.0259 ± 0.0013 | -1.5937 ± 0.0010 | 0.6580 ± 0.0003 | -0.0015 | +0.0006 | +0.0005 |
| +0.005 (0.711) | 1.0275 ± 0.0011 | -1.5941 ± 0.0010 | 0.6578 ± 0.0002 | +0.0001 | +0.0003 | +0.0002 |
| +0.01 (0.716) | 1.0278 ± 0.0017 | -1.5925 ± 0.0008 | 0.6581 ± 0.0003 | +0.0004 | +0.0019 | +0.0006 |

- **P24: sí.** Cambio máximo de c0 0.0015, de c1 0.0021 (< 0.05) y de la potencia 0.0006 (< 0.02); todos del orden del error MC
  (±0.001–0.002), porque el desplazamiento lo absorbe la predicción y las densidades fraccionales apenas cambian.

## (cambio por z_eff) Modo evaluación: orden de lectura

1. **z_eff y covarianza de DR3** → `calibrar COV.json TAG` lee **solo** las claves `labels` (de donde salen los z_eff: `TIPO_TRAZADOR_zZ`) y `cov`; empareja
   cada fila de DR2 con la de DR3 del mismo tipo y trazador con z_eff ± 0.01 (LRG2 fuera de 0.706 ± 0.01 → no evaluable), registra los desplazamientos,
   aplica las reglas de v3 (simetría, Cholesky, C2 − C3 ⪰ 0, Cov(e) ⪰ 0, correlación entre bins < 1e-6, ρ²_LRG2,eff ∈ [0.40, 0.80]), construye las
   densidades exactas fraccionales con el modelo en z_DR3 y C3, calcula c0 y c1 con 10⁷ muestras (semillas 500) y escribe la calibración con su sha256.
   **Si z_eff viene en el mismo fichero que los valores, el script lee solo las columnas de z_eff y covarianza en este paso** (nunca la clave `mean`), y se
   declara en el registro antes de abrirlo.
2. **Solo después**, `evaluar COV.json VALORES.json TAG` lee los valores, reajusta ΛCDM sin LRG2 con el optimizador original, C3 y el modelo en z_DR3
   (cálculo completo), calcula r3 fraccional, d²_DR3 y Λ una sola vez, y da la salida principal (sin SN) y la secundaria (Pantheon+).
   Los sha256 son controles de coherencia; el compromiso es el orden público de los commits.

**Prueba con datos falsos** (marcados "DATOS FALSOS"; `pronostico_dr3_lrg2_v4_falso_run.log`):

1. z_eff de LRG2 = 0.712, C3 = 0.6·C2, valores de DR2: emparejamiento único con desplazamiento +0.006 registrado en las dos filas de LRG2 (0 en las demás),
   estado "ok"; c0 / c1 (sin SN) = +1.0260 / -1.5918, potencia 0.6582 (a ±0.002 de los de z = 0.706). Evaluación: modelo de LRG2 evaluado en 0.712, r3 fraccional =
   (-0.0171, -0.0274), d² = 13.05, Λ = +2.484 → principal "favorece H1" (secundaria "favorece H1"). (Con los valores de DR2 repetidos en
   z = 0.712 el residuo crece porque el modelo en 0.712 es mayor: es un dato falso, no una predicción.)
2. z_eff de LRG2 = 0.720: las filas de LRG2 de DR2 no tienen homólogo en 0.706 ± 0.01 (las de DR3 en 0.720 se ignoran) → **"no evaluable"** en la
   calibración y en la evaluación, sin calcular Λ.

## PROCEDIMIENTO PÚBLICO — este es el compromiso (orden de lectura actualizado; cambio por z_eff)

Al salir DESI DR3:

1. Descargar **solo** los z_eff de los bins y el fichero de covarianza de las medidas BAO de DR3 (si vienen en el mismo fichero que los valores, declararlo
   en el registro antes de abrirlo; el script solo lee z_eff y covarianza en este paso), guardándolos con su sha256 y su URL/fecha.
2. Calibrar con el script congelado (`pronostico_dr3_lrg2_v4.py calibrar`, commit de este registro) y hacer **commit y push públicos** de la calibración y
   de su `.sha256`.
3. Solo después, descargar los valores y evaluar **una sola vez** (`evaluar`).
4. Publicar el resultado con una de las cinco salidas: **favorece H1 / favorece H0 / no concluyente / modelo roto / no evaluado**, con estos umbrales y sin
   añadir umbrales después ("no evaluable" por formato o por z_eff cuenta como "no evaluado").

El orden lo prueban las fechas de los commits públicos. Fecha límite: 31-12-2028 (después, "no evaluado", no fallido).

## Lo que no cambia respecto a v3

Supuestos (anidamiento matricial con C3 publicada y sus comprobaciones con tolerancia de autovalores, externos fijos, prior plano en p, álgebra por
Cholesky), fórmulas de referencia, validación del linealizado, umbrales de alta precisión de los casos registrados (sin SN, ρ² = 0.6: c0 = 1.0274 ± 0.0013,
c1 = -1.5944 ± 0.0012, potencia 0.6575 ± 0.0002; P18, P19, P23 sí), barrido de ρ², rupturas a–g y cautelas: véase `PRONOSTICO_DR3_REGISTRO_v3.md`. Ni H1 ni Λ prueban el
modelo OU ni σ_X > 0.

## Ficheros nuevos de esta versión

`pronostico_dr3_lrg2_v4.py`, `pronostico_dr3_lrg2_v4.sha256`, `pronostico_dr3_lrg2_v4.json`, `pronostico_dr3_lrg2_v4_run.log`; prueba con datos falsos:
`pronostico_dr3_lrg2_v4_falso_z0712_dr3_{cov,valores}.json`, `..._falso_z0720_dr3_{cov,valores}.json`, `pronostico_dr3_lrg2_v4_calibracion_falso_z07{12,20}.json`
(+ `.sha256`), `pronostico_dr3_lrg2_v4_resultado_falso_z07{12,20}.json`, `pronostico_dr3_lrg2_v4_falso_run.log`. Nada existente se modifica.
