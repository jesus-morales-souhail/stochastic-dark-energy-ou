# PRONÓSTICO DR3 para el residuo de LRG2 — REGISTRO, versión 5 (13-11-2026)

**PENDIENTE DE COMPROBACIÓN DEL REVISOR — todavía sin tag.**

**Esta versión sustituye formalmente a v3 (`PRONOSTICO_DR3_REGISTRO_v3.md`, commit 6f1628b) y a v4 (`PRONOSTICO_DR3_REGISTRO_v4.md`, commit 1439356) para el
registro; v3 y v4 quedan como historia.** Parte del texto de v4 con las correcciones pedidas por la tercera revisión `REVISION3_PRONOSTICO_DR3.md` (commit 53bcb70,
veredicto "no registrable hasta corregir"), cada una marcada **(cambio por la revisión 3)**. No cambian P18, P19, P23, P24, H0, H1, el anidamiento matricial, los umbrales
ni la regla de decisión; no se rehace ninguna simulación científica (las densidades, umbrales y tablas de v3 y v4 siguen siendo los registrados). Código:
`pronostico_dr3_lrg2_v5.py` (importa v4 sin modificarlo; sha256 `734a4c8de4ea85df6a911abe9cd67825a5c5ec1ba029c0a43993dde4f48b1839`, también en `pronostico_dr3_lrg2_v5.sha256`); `pronostico_dr3_lrg2_v5.json`
(sha256 `913c8936551081ba82bc0af3d58b9389266694bc80d68ae36dddeffcfa59d34f`). `git_head_at_run` = 53bcb70 (commit anterior; el registro lo fija el tag).

## Predicciones (sin cambio; copiadas literales de sus versiones)

- P18. Con ρ² = σ3²/σ2² = 0.6 y falsa alarma del 5 %, la potencia del test para distinguir H0 de H1 en LRG2 está entre 0.5 y 0.75. (sí: 0.6575 ± 0.0002, v3)
- P19. Con ρ² = 0.6, la mediana de d²_DR3 de LRG2 está entre 3 y 6 bajo H0 y entre 9 y 15 bajo H1. (sí: 4.764 y 12.320, v3)
- P23. Con los umbrales de alta precisión, la potencia del caso principal (ρ² = 0.6, sin SN) queda en 0.656 ± 0.010. (sí: 0.6575 ± 0.0002, v3)
- P24. Desplazar el z_eff de LRG2 en ±0.01 (recalculando el modelo y el jacobiano en el nuevo z) cambia c0 y c1 menos de 0.05 y la potencia del caso principal
  menos de 0.02. (sí: máx 0.0015, 0.0021 y 0.0006, v4)

## Lo que no cambia (véanse v3 y v4)

Contexto, H0 y H1, "este pronóstico no contiene ninguna prueba del modelo OU", salida principal única sin SN (Pantheon+ secundaria), anidamiento matricial
n3 = C3·C2⁻¹·n2 + e con sus comprobaciones por bin y la tolerancia de autovalores (λ < −1e-10·λ_max = fallo real), álgebra por Cholesky, umbrales de alta precisión
(10⁷ muestras exactas; sin SN, ρ² = 0.6: c0 = 1.0274 ± 0.0013, c1 = −1.5944 ± 0.0012, potencia 0.6575 ± 0.0002), barrido de ρ², redshifts condicionantes (DR2) y
de predicción (DR3) con el modelo y el jacobiano en z_DR3, H1 con μ fraccional, residuos fraccionales, reglas a–g, procedimiento público y cautelas.

## (cambio por la revisión 3) Frontera inclusiva de la regla f con tolerancia explícita

Un bin de DR3 es homólogo de uno de DR2 si tiene el mismo tipo (D_M, D_H o D_V) y el mismo trazador y **|z_DR3 − z_DR2| ≤ 0.01 + 1e-9**. La tolerancia 1e-9 es
numérica: en coma flotante 0.716 − 0.706 = 0.010000000000000009, y sin ella un z_eff publicado exactamente en el borde (0.696, 0.716, 0.500, 0.520) se declaraba
"no evaluable" contra el texto ("0.706 ± 0.01", inclusive) y contra P24 (calculada en ±0.01). Pruebas de frontera (`pronostico_dr3_lrg2_v5.json` → `boundary_tests`):

| bin | z_eff de DR3 | Δz | resultado |
|---|---|---|---|
| LRG2 | 0.696 | -0.0100 | evaluable |
| LRG2 | 0.6965 | -0.0095 | evaluable |
| LRG2 | 0.706 | +0.0000 | evaluable |
| LRG2 | 0.7155 | +0.0095 | evaluable |
| LRG2 | 0.716 | +0.0100 | evaluable |
| LRG2 | 0.6959 | -0.0101 | no evaluable |
| LRG2 | 0.7161 | +0.0101 | no evaluable |
| LRG2 | 0.72 | +0.0140 | no evaluable |
| LRG1 | 0.5 | -0.0100 | evaluable |
| LRG1 | 0.52 | +0.0100 | evaluable |
| LRG1 | 0.4999 | -0.0101 | no evaluable |
| LRG1 | 0.5201 | +0.0101 | no evaluable |

## (cambio por la revisión 3) Controles de coherencia que fallan cerrados

`calibrar` y `evaluar` comprueban el sha256 del script contra `pronostico_dr3_lrg2_v5.sha256` y `evaluar` además el de la calibración contra su `.sha256`, la
covarianza y el script con que se calibró. Si falta cualquiera de esos ficheros o el hash no coincide (antes `match: None` dejaba pasar y `calibrar` solo anotaba
`match: False`), **los dos modos se niegan con una salida declarada "control de coherencia fallido"** (escrita en el json de calibración o de resultado), no con una
excepción. Siguen siendo **controles de coherencia, no el compromiso**: se pueden saltar reescribiendo ficheros, y lo que impide decidir después de ver los datos es
el orden público de los commits (procedimiento público de v4). La única negativa que sigue siendo una excepción es repetir un TAG ya evaluado (el fichero de
resultado existe y no se sobrescribe).

## (cambio por la revisión 3) Emparejamiento por identidad en la evaluación

`calibrar` guarda las etiquetas de DR3 emparejadas con cada fila de DR2 (`labels_calibrated`). `evaluar` empareja **cada fila de valores con la calibración por
trazador, cantidad (D_M/r_d o D_H/r_d) y z_eff idéntico (|Δz| ≤ 1e-9)**, no por índice: los valores deben estar en el z_eff que se calibró (mismo release).
Salidas declaradas: valores en un z distinto del calibrado o fila ausente → **"no evaluable"** (antes se ejecutaba por índice); filas reordenadas → evaluable
(mismo resultado); bins nuevos en los valores → se ignoran (regla e) y se listan; filas ambiguas → "no evaluable". Nunca `AssertionError`. Los dos ficheros ya
no tienen que compartir el orden.

## (cambio por la revisión 3) Denominador congelado

En `evaluar` el residuo fraccional y su covarianza usan como denominador **m(p2, z_DR3) guardado en la calibración** (`m2_bao3_LRG2`), el mismo de las densidades,
en lugar de m(p3, z_DR3) (v4). Así Λ es exactamente el de v3 en unidades absolutas con el modelo en z_DR3, sin la mezcla de denominadores (la revisión 3 midió
para esa mezcla un desplazamiento de r3 de 1.5e-3 σ típico y ≤ 0.04 σ, y ΔΛ hasta 0.09–0.19 en la cola). Cambio de Λ medido en las pruebas falsas (denominador
congelado frente al de v4), sin SN / Pantheon+: C3 = 0.6·C2 en z = 0.706: +1.1269 frente a +1.1372 (Δ = -0.0102) / +1.1692 frente a +1.1819 (Δ = -0.0127); z_LRG2 = 0.712:
+2.4693 frente a +2.4843 (Δ = -0.0151) / +2.5955 frente a +2.6140 (Δ = -0.0185). El Λ de v5 con z = 0.706 (1.1269) coincide con el de v3 (1.1269) como debe. El script informa también
del Λ con el denominador de v4 (`Lambda_v4_denominator_m_p3`).

## (cambio por la revisión 3) Término omitido en la conversión de μ

Bajo H1 fraccional, μ3,i = μ2,i · m_i(p_true, z3)/m_i(p_true, z2) componente a componente (D_M, D_H); v4 usa el cociente en p2 para que r3 siga siendo afín en
gaussianas. El término omitido Δ_i = μ2,i·[ρ_i(p_true) − ρ_i(p2)], medido por la revisión 3 sobre el posterior de (p, μ) con Δz = ±0.01, tiene **mediana 2e-6 (D_M) y
1.2e-5 (D_H) y máximo 1.6e-4 en D/r_d (≤ 7.5e-4 σ)**; en el caso principal registrado (z_DR3 = z_DR2) es idénticamente cero. Sustituir la linealización por la
expresión completa cambia **Δc1 = −2e-5 y Δpotencia = +5e-6** (≪ error MC). La frase de v4 "el término de segundo orden dp·μ es ~1e-5" se refería a la mediana.

## (cambio por la revisión 3) Advertencia sobre la regla d

Una C3 con la **misma** varianza que DR2 en una componente de un bin (por ejemplo ρ² = 1 en D_H de LRG2 con 0.6 en D_M) ya rompe C2 − C3 ⪰ 0 por el término
cruzado (autovalor real negativo, −0.005 en el caso medido por la revisión) y el modelo se declara roto. DR3 debe ser **más preciso que DR2 en todas las
direcciones**, no "parecido": "igual" en una componente ya rompe.

## Comprobación v5 = v4 y pruebas con datos falsos

Con los z_eff de DR2 (ρ² = 0.6, mismas semillas que el registro de v4), la calibración de v5 reproduce la de v4 exactamente: Δc0 = Δc1 = Δpotencia = 0 y
densidades idénticas (máx |Δa| = |ΔΣ| = 0) en las dos combinaciones (v5 llama al mismo código de v4 para densidades y umbrales; solo cambian el emparejamiento,
los controles y el denominador de la evaluación).

| prueba (DATOS FALSOS) | calibrar | evaluar (salida principal) |
|---|---|---|
| C3 = 0.6·C2, z = 0.706, valores DR2 (prueba de v3) | ok | favorece H1; Λ = +1.1269 |
| Var(D_H de LRG2) = 2 × DR2 (v3) | modelo roto (LRG2: C2_minus_C3, Cov_e; rho2_LRG2,eff = 1.685 fuera de [0.4, 0.8]) | modelo roto (['LRG2: C2_minus_C3, Cov_e', 'rho2_LRG2,eff = 1.685 fuera de [0.4, 0.8]']) |
| correlación 1e-3 entre bins (v3) | modelo roto (correlación entre bins |r| máx 1.00e-03 ≥ 1e-06) | modelo roto (['correlación entre bins |r| máx 1.00e-03 ≥ 1e-06']) |
| correlación 1e-9 (v3) | ok | favorece H1; Λ = +1.1269 |
| z_LRG2 = 0.712, C3 = 0.6·C2 (v4) | ok | favorece H1; Λ = +2.4693 |
| z_LRG2 = 0.720 (v4) | no evaluable (bins de DR2 sin homólogo único en DR3 (mismo tipo y trazador, |Δz| ≤ 0.01 + 1e-9)) | no evaluable (bins de DR2 sin homólogo único en DR3 (mismo tipo y trazador, |Δz| ≤ 0.01 + 1e-9) |
| fichero .sha256 del script ausente | control de coherencia fallido (falta el fichero no_existe.sha256) | control de coherencia fallido (falta el fichero no_existe.sha256) |
| hash del script alterado | control de coherencia fallido (el sha256 del script no coincide con el registrado) | control de coherencia fallido (el sha256 del script no coincide con el registrado) |
| valores en z = 0.706 con calibración en z = 0.712 | ok | no evaluable (filas de valores sin pareja en la calibración (trazador, cantidad y z_eff idénti) |
| filas de valores en orden inverso (calibración z = 0.712) | ok | favorece H1; Λ = +2.4693 |
| un bin nuevo delante en los valores (DV_ELG2b_z1.100) | ok | favorece H1; Λ = +2.4693; filas ignoradas ['DV_ELG2b_z1.100'] |

Además, la segunda evaluación con un TAG ya evaluado se niega ("Λ se calcula una sola vez", excepción). Ficheros: `pronostico_dr3_lrg2_v5_falso_run.log`,
`pronostico_dr3_lrg2_v5_calibracion_*.json` (+ `.sha256`), `pronostico_dr3_lrg2_v5_resultado_*.json`; datos falsos nuevos
`pronostico_dr3_lrg2_v5_falso_reordenado_dr3_valores.json`, `..._falso_binnuevo_dr3_valores.json`, `..._falso_sha256_alterado.sha256` (los de v3 y v4 se reutilizan).

## Procedimiento público (sin cambio respecto a v4; con el script congelado de esta versión)

1. Descargar solo los z_eff y la covarianza de las medidas BAO de DR3 (si vienen con los valores en un fichero, declararlo antes de abrirlo; `calibrar` solo lee `labels` y `cov`).
2. Calibrar con `pronostico_dr3_lrg2_v5.py calibrar` (commit de este registro) y hacer commit y push públicos de la calibración y de su `.sha256`.
3. Solo después, descargar los valores y evaluar una sola vez (`evaluar`): los valores se emparejan por identidad con la calibración.
4. Publicar una de las cinco salidas: **favorece H1 / favorece H0 / no concluyente / modelo roto / no evaluado** ("no evaluable" y "control de coherencia fallido"
   cuentan como "no evaluado" si no se resuelven con registro público antes de evaluar). Fecha límite 31-12-2028.

## Cautelas

Las de v3 y v4. Ni H1 ni Λ prueban el modelo OU ni σ_X > 0. Pendiente de comprobación del revisor; sin tag.
