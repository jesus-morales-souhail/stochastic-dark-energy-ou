# EBOSS_LRG: ¿es el desvío de DESI LRG2 compatible con SDSS/eBOSS DR16 LRG? (13-11-2026)

## Contexto

La preferencia por σ_X > 0 depende del bin LRG2 de DESI DR2 (z_eff = 0.706; JACKKNIFE_BINS.md), con p local 0.028 bajo ΛCDM y p global 0.177 con 7 bins
(PREDICTIVO_BINS.md). Comprobación con otro instrumento ya existente: SDSS/eBOSS DR16 LRG en z_eff = 0.698 (Bautista et al. 2021, arXiv:2007.08993;
Gil-Marín et al. 2020, arXiv:2007.08994). Esta prueba **no puede confirmar** LRG2: eBOSS tiene errores mayores (1.8 % y 2.8 % frente a 1.0 % y 1.7 % de DESI) y
comparte parte del volumen con DESI. Puede decir si es compatible o si lo contradice. Lo esperado bajo cada hipótesis se calculó y se publicó (commit
`b8d48a5`) **antes** de leer los valores observados de eBOSS.

## Predicciones (de Claude, escritas antes de ejecutar, copiadas literales) y declaración

- P20. Bajo H1 (desvío de LRG2 real, con el tamaño y la incertidumbre estimados en DESI), la probabilidad de que eBOSS dé p local < 0.05 frente a ΛCDM es < 0.4.
  (Se evalúa con ρ_x = 0; los demás ρ_x son sensibilidad, no nuevas predicciones.)
- P21. El p local observado de eBOSS frente a ΛCDM es > 0.05.
- P22. Con ρ_x = 0, DESI LRG2 y eBOSS son compatibles entre sí (p > 0.05 para su diferencia).

**Declaración:** P20 no depende de los valores observados de eBOSS. P21 y P22 NO son predicciones ciegas: quien las escribió (Claude) había visto los valores
publicados de eBOSS (los de consenso BAO+RSD en el mismo mensaje; los BAO-only fueron citados después por otra IA), aunque no los había comparado
numéricamente con ΛCDM.

## Lo que se encontró de la ejecución cancelada

Una ejecución anterior de este encargo se canceló a medias. Al retomar: `git status` mostraba tres ficheros sin seguimiento (`eboss_lrg_.py`,
`eboss_lrg_potencia.json`, `eboss_lrg_potencia_run.log`) y **ningún commit** nuevo (`git log` terminaba en 1439356, el de v4 del pronóstico). Se renombraron
con el sufijo `_cancelado` y no se usan ni se incluyen en ningún commit. Esa ejecución **sí llegó a leer e imprimir los valores observados de eBOSS**: en su
exploración de los datos hizo `cat` del fichero `sdss_DR16_LRG_BAO_DMDH.dat` (dos veces: al localizarlo y al compararlo byte a byte con las copias remotas), de
modo que los valores aparecieron en la salida de esas órdenes antes de que existiera ningún commit con lo esperado. Su script de potencia, en cambio, solo
leía la covarianza, y su json no contiene los valores. Consecuencia: **lo esperado bajo H0 y H1 (paso 3) no se vio alterado** (es el mismo cálculo con la misma
semilla, y no usa los valores), pero la ejecución retomada no fue ciega a los valores en el sentido estricto: Claude los había visto en la salida de la
ejecución cancelada antes del commit del paso 4. Se declara aquí; P21 y P22 ya no eran ciegas por la declaración anterior.

## Procedencia (detalle en `data/PROCEDENCIA_EBOSS_DR16_LRG.md`)

Ficheros `sdss_DR16_LRG_BAO_DMDH.dat` (sha256 `b3317e75…2d91be`) y `sdss_DR16_LRG_BAO_DMDH_covtot.txt` (sha256 `1a45e106…f6160`) del clon local de
`CobayaSampler/bao_data` (commit bb0c1c9, ficheros de 2022-02-15), idénticos a los de `raw.githubusercontent.com/CobayaSampler/bao_data/master/`; la covarianza
es idéntica a la de la publicación de SDSS (`svn.sdss.org/.../DR16cosmo/tags/v1_0_1/likelihoods/BAO-only/`), cuyo README describe los `*_DMDH` como aproximación
gaussiana a (D_M/r_d, D_H/r_d); la URL del vector de datos en SDSS dio 404. **Producto BAO-only** (no BAO+RSD ni el consenso final): consenso BAO-only de los
análisis en espacio de configuración (Bautista et al.) y de Fourier (Gil-Marín et al.), según la literatura (Alam et al. 2021, Tabla 3); ni el README ni el clon
lo dicen fichero a fichero. Muestra: 174 816 LRG de eBOSS + 202 642 galaxias CMASS de BOSS, 0.6 < z < 1, z_eff = 0.698 (resumen de Bautista et al.). Valores y
covarianza del mismo producto. En el paso 3 se leyeron solo los metadatos (z_eff, tipo de fila: columnas 0 y 2) y la covarianza: σ = (0.328, 0.533), coeficiente
-0.334. Valores (leídos en el paso 5): D_M/r_d = 17.858, D_H/r_d = 19.326. No existe covarianza cruzada publicada entre DESI DR2 LRG2 y eBOSS DR16 LRG.

## Modelo (aproximación gaussiana) y deducciones

- **Desvío fraccional** de cada medida respecto a la predicción ΛCDM con p_ref (ajuste sin DESI LRG2 de PREDICTIVO_BINS, sin SN: (h, ω_m) =
  (0.67974, 0.14129); con Pantheon+ por separado), por componente (D_M/r_d, D_H/r_d), con el modelo evaluado en el z_eff de cada medida (0.706 para DESI,
  0.698 para eBOSS; `pronostico_dr3_lrg2_v4.model_bao`): δ = (datos − modelo(p_ref))/modelo(p_ref). Modelo en p_ref: DESI (17.5603, 20.0907), eBOSS (17.3992, 20.1879).
  δ_D = (-0.0119, -0.0316) (DESI LRG2, el r de PREDICTIVO_BINS dividido por el modelo; reproducido a 1e-9).
- **Descomposición:** δ_D = μ + J_D q + n_D; δ_E = μ + J_E q + n_E. μ = desvío físico (0 bajo H0; bajo H1 el mismo desvío fraccional en z = 0.706 y 0.698:
  **supuesto físico del contraste, no se deduce de que los z estén próximos**); q = p_verdadero − p_ref ~ N(0, Σ) con Σ del ajuste sin DESI LRG2
  (√diag = (0.00254, 0.00071) en h, ω_m); n_D, n_E errores de medida con covarianzas fraccionales C_D = C_LRG2/(m_D m_Dᵀ), C_E = C_eBOSS/(m_E m_Eᵀ)
  (√diag: DESI (0.0102, 0.0166); eBOSS (0.0189, 0.0264)).
- **J es el jacobiano del RESIDUO fraccional**, no de la distancia predicha: J_D = −∂ ln modelo_D/∂p y J_E = −∂ ln modelo_E/∂p en p_ref, por diferencias
  centradas (ε = 1e-4 en h, 1e-5 en ω_m; la columna de M es 0 con Pantheon+), con la misma definición y convención para DESI y eBOSS:
  J_D = [[1.060e+00, -8.240e-01], [6.647e-01, 1.252e-01]], J_E = [[1.063e+00, -8.327e-01], [6.712e-01, 1.097e-01]] (filas D_M, D_H; columnas h, ω_m). Las dos difieren en < 1 %, como corresponde a Δz = 0.008.
- **Covarianza cruzada de los errores:** C_x = ρ_x · C_D^{1/2} · C_E^{1/2} (raíces simétricas), ρ_x = 0, 0.3, 0.5. **Es una familia de correlaciones hipotéticas para
  sensibilidad, no una medida de cuánto volumen comparten DESI y eBOSS.** No hay covarianza cruzada publicada que usar además.
- **H0:** V_D = C_D + J_D Σ J_Dᵀ; V_E = C_E + J_E Σ J_Eᵀ; V_DE = C_x + J_D Σ J_Eᵀ (covarianza conjunta de (δ_D, δ_E) gaussiana con media 0); condicionando,
  **δ_E | δ_D, H0 ~ N(V_DEᵀ V_D⁻¹ δ_D, V_E − V_DEᵀ V_D⁻¹ V_DE)**. Esto condiciona también los parámetros cosmológicos con LRG2 (el término J_D Σ J_Eᵀ): bajo H0 el
  desvío de DESI se reparte entre ruido y un q que eBOSS hereda.
- **H1 (μ libre, prior plano):** δ_E − δ_D = (J_E − J_D) q + (n_E − n_D), que no contiene μ. **Al marginalizar el desplazamiento común μ con prior plano, observar δ_D no
  aporta información sobre q** (μ absorbe δ_D: para cualquier (q, n_D) hay un μ que lo explica con la misma verosimilitud), así que δ_E − δ_D es independiente de
  δ_D y **δ_E | δ_D, H1 ~ N(δ_D, C_Δ)**, con C_Δ = Var[(J_E − J_D) q + n_E − n_D] = (J_E − J_D) Σ (J_E − J_D)ᵀ + C_E + C_D − Cov(n_E, n_D) − Cov(n_D, n_E) =
  C_D + C_E − C_x − C_xᵀ + (J_E − J_D) Σ (J_E − J_D)ᵀ.
- **p local de eBOSS frente a ΛCDM:** d²_E = δ_Eᵀ V_E⁻¹ δ_E con 2 g.d.l. (eBOSS por sí solo, sin condicionar en DESI). Tirones marginales = δ_E,i / √V_E,ii, **con las
  desviaciones típicas de V_E, no con las de la distribución condicionada tras DESI**.

## Comprobaciones numéricas (antes de las simulaciones grandes; `eboss_lrg_potencia.json` → `checks`, `check_small_sample`)

| combinación | ρ_x | autovalores de [[C_D, C_x], [C_xᵀ, C_E]] | autovalores de C_Δ | autovalores de Cov(δ_E \| δ_D, H0) | muestra de 2000: máx \|Δmedia\|/σ, máx \|Δvar\| rel. (1 σ_MC = 0.022 / 0.032) | ok |
|---|---|---|---|---|---|---|
| sinSN | 0 | 8.17e-05, 2.88e-04, 2.99e-04, 7.64e-04 | 3.70e-04, 1.06e-03 | 3.00e-04, 7.64e-04 | 0.047 / 0.040 | sí |
| sinSN | 0.3 | 7.19e-05, 2.58e-04, 2.98e-04, 8.05e-04 | 2.78e-04, 7.76e-04 | 2.64e-04, 6.95e-04 | 0.029 / 0.047 | sí |
| sinSN | 0.5 | 5.63e-05, 1.98e-04, 3.14e-04, 8.65e-04 | 2.17e-04, 5.85e-04 | 2.16e-04, 5.73e-04 | 0.025 / 0.018 | sí |
| Pantheon+ | 0 | 8.15e-05, 2.88e-04, 2.99e-04, 7.63e-04 | 3.70e-04, 1.06e-03 | 2.99e-04, 7.63e-04 | 0.029 / 0.047 | sí |
| Pantheon+ | 0.3 | 7.18e-05, 2.58e-04, 2.98e-04, 8.04e-04 | 2.78e-04, 7.75e-04 | 2.64e-04, 6.94e-04 | 0.015 / 0.027 | sí |
| Pantheon+ | 0.5 | 5.62e-05, 1.98e-04, 3.13e-04, 8.64e-04 | 2.16e-04, 5.84e-04 | 2.16e-04, 5.72e-04 | 0.042 / 0.035 | sí |

Las seis covarianzas conjuntas, las seis C_Δ y las doce condicionadas son definidas positivas (Cholesky sin fallo; V_E también); la muestra pequeña reproduce
media y covarianza analíticas dentro de 2.2 σ_MC en todos los casos.

## Lo esperado bajo H0 y H1 (paso 3; commit `b8d48a5`, antes de leer los valores; 100000 δ_E por distribución, semilla 20261116)

| combinación | ρ_x | δ_E \| δ_D, H0: media ± √diag | δ_E \| δ_D, H1: media ± √diag | P(p local < 0.05 \| H0) | potencia P(p local < 0.05 \| H1) | mediana de d²_E: H0 / H1 |
|---|---|---|---|---|---|---|
| sinSN | 0 | (-0.0028, -0.0014) ± (0.0191, 0.0264) | (-0.0119, -0.0316) ± (0.0215, 0.0312) | 0.052 | 0.331 | 1.41 / 3.98 |
| sinSN | 0.3 | (-0.0086, -0.0160) ± (0.0180, 0.0252) | (-0.0119, -0.0316) ± (0.0186, 0.0266) | 0.096 | 0.270 | 1.88 / 3.57 |
| sinSN | 0.5 | (-0.0124, -0.0257) ± (0.0163, 0.0229) | (-0.0119, -0.0316) ± (0.0164, 0.0231) | 0.162 | 0.223 | 2.76 / 3.34 |
| Pantheon+ | 0 | (-0.0028, -0.0014) ± (0.0191, 0.0264) | (-0.0129, -0.0321) ± (0.0214, 0.0312) | 0.053 | 0.345 | 1.40 / 4.14 |
| Pantheon+ | 0.3 | (-0.0091, -0.0162) ± (0.0180, 0.0252) | (-0.0129, -0.0321) ± (0.0186, 0.0266) | 0.099 | 0.291 | 1.93 / 3.77 |
| Pantheon+ | 0.5 | (-0.0133, -0.0261) ± (0.0163, 0.0228) | (-0.0129, -0.0321) ± (0.0163, 0.0231) | 0.178 | 0.241 | 2.92 / 3.52 |

Lectura: bajo H1 eBOSS debería reproducir el desvío de DESI (media = δ_D = (−1.2 %, −3.2 %)) pero con sus errores mayores, de modo que solo en el 33 % de los
casos daría p < 0.05 frente a ΛCDM: **la prueba tiene poca potencia**. Con ρ_x > 0, eBOSS repite parte del desvío de DESI incluso bajo H0 (la media condicionada
bajo H0 pasa de (−0.3 %, −0.1 %) con ρ_x = 0 a (−1.2 %, −2.6 %) con ρ_x = 0.5): la falsa alarma P(p < 0.05 | H0) sube de 0.052 a 0.162 y la diferencia entre
hipótesis se reduce (0.223 frente a 0.162 con ρ_x = 0.5). **P20: sí** (potencia 0.331 < 0.4 con ρ_x = 0, sin SN; Pantheon+ 0.345).

## Resultados (paso 6; valores leídos después del commit `b8d48a5`)

δ_E observado = (+0.0264, -0.0427) (sin SN; Pantheon+ (+0.0254, -0.0432)): eBOSS mide D_M/r_d un 2.6 % **por encima** de ΛCDM y D_H/r_d un 4.3 % por debajo.
**Tirones marginales** (con √V_E): D_M +1.38, D_H -1.62 (sin SN; Pantheon+ +1.33, -1.63). **Signos (+, −): no coinciden con los de DESI LRG2 (−, −)**: el signo de D_H
coincide y el de D_M es el contrario. **p local de eBOSS frente a ΛCDM: 0.179** (d²_E = 3.44 con 2 g.d.l.; Pantheon+ 0.183). **P21: sí** (> 0.05; no ciega, véase la declaración).

Posición del d²_E observado en las distribuciones del paso 3:

| combinación | ρ_x | d²_E observado | p local | percentil bajo H0 | percentil bajo H1 |
|---|---|---|---|---|---|
| sinSN | 0 | 3.442 | 0.1789 | 81.7 % | 44.6 % |
| sinSN | 0.3 | 3.442 | 0.1789 | 72.9 % | 48.5 % |
| sinSN | 0.5 | 3.442 | 0.1789 | 59.8 % | 51.4 % |
| Pantheon+ | 0 | 3.395 | 0.1831 | 81.2 % | 42.4 % |
| Pantheon+ | 0.3 | 3.395 | 0.1831 | 71.3 % | 45.7 % |
| Pantheon+ | 0.5 | 3.395 | 0.1831 | 56.9 % | 48.4 % |

El d²_E observado está en el percentil 82 de lo esperado bajo H0 y en el 45 de lo esperado bajo H1 (ρ_x = 0): ni raro bajo H0 ni raro bajo H1.

## Compatibilidad DESI LRG2–eBOSS (paso 7) y cociente de verosimilitudes (paso 8)

Δ = δ_D − δ_E, d² = Δᵀ C_Δ⁻¹ Δ con 2 g.d.l. Bajo H1 esta diferencia es la misma que la predictiva (δ_E | δ_D, H1 ~ N(δ_D, C_Δ)): la prueba de compatibilidad y
la posición de δ_E en la predictiva de H1 son la misma cantidad. Λ_eBOSS = ln p(δ_E | δ_D, H1) − ln p(δ_E | δ_D, H0), descriptivo, sin regla de decisión nueva.

| combinación | ρ_x | Δ = δ_D − δ_E | d² | p (2 g.d.l.) | Λ_eBOSS |
|---|---|---|---|---|---|
| sinSN | 0 | (-0.0383, +0.0111) | 3.269 | 0.1951 | -0.096 |
| sinSN | 0.3 | (-0.0383, +0.0111) | 4.355 | 0.1133 | -0.282 |
| sinSN | 0.5 | (-0.0383, +0.0111) | 5.594 | 0.0610 | +0.012 |
| Pantheon+ | 0 | (-0.0383, +0.0111) | 3.268 | 0.1951 | -0.131 |
| Pantheon+ | 0.3 | (-0.0383, +0.0111) | 4.353 | 0.1134 | -0.318 |
| Pantheon+ | 0.5 | (-0.0383, +0.0111) | 5.592 | 0.0611 | +0.006 |

**P22: sí** (p = 0.195 > 0.05 con ρ_x = 0, sin SN; Pantheon+ 0.195; no ciega). Con ρ_x = 0.3 y 0.5 la compatibilidad baja a 0.113 y 0.061 porque una correlación
positiva hace esperar diferencias menores. Λ_eBOSS está entre −0.3 y +0.01 en los seis casos: **eBOSS no distingue H0 de H1** (|Λ| < 0.35, mucho menos que 1).

## Tabla de lectura

| caso | criterio | ¿se da? |
|---|---|---|
| eBOSS compatible con ΛCDM y con DESI | p local > 0.05 y p de compatibilidad > 0.05 | **sí** con ρ_x = 0 y 0.3 (p = 0.179 y 0.195 / 0.113); con ρ_x = 0.5 la compatibilidad queda en 0.061, en el límite |
| eBOSS compatible con ΛCDM pero no con DESI | p local > 0.05 y p de compatibilidad < 0.05 | no (con ρ_x = 0.5 se acerca: 0.061) |
| eBOSS discrepante con ΛCDM y compatible con DESI | p local < 0.05 y p de compatibilidad > 0.05 | no |

Lectura corta: eBOSS DR16 LRG no contradice el desvío de DESI LRG2 ni lo confirma. Su propio residuo frente a ΛCDM (p = 0.18) tiene el signo de D_H de DESI
pero el contrario en D_M, está donde se espera tanto bajo H0 (percentil 82) como bajo H1 (percentil 45), y Λ_eBOSS ≈ 0. Con la potencia de 0.33, este resultado
era el más probable bajo cualquiera de las dos hipótesis.

## Cautelas

1. **Solape de volumen y ρ_x hipotético.** DESI y eBOSS observan en parte el mismo cielo y las mismas estructuras; ρ_x = 0, 0.3, 0.5 es una familia de sensibilidad,
   no una medida de ese solape. Con ρ_x > 0, coincidir con DESI es menos informativo (eBOSS repite parte del desvío incluso bajo H0) y la compatibilidad se exige
   más estricta (p de 0.195 a 0.061).
2. **Producto BAO-only frente al consenso.** Se usa el BAO-only (consenso de configuración y Fourier) con su propia covarianza; el consenso BAO+RSD (17.65, 19.77)
   daría otro δ_E. La identificación del fichero como consenso de los dos espacios viene de la literatura, no de un README fichero a fichero.
3. **Muestra mixta eBOSS + BOSS CMASS (z > 0.6)**, no un sondeo independiente de BOSS; su z_eff y selección no son los de DESI LRG2 (0.6 < z < 0.8).
4. **Diferencia de z y μ común como supuesto.** Δz = 0.008 entre las dos medidas; que el desvío fraccional sea el mismo en 0.706 y 0.698 es un supuesto físico de H1,
   no una consecuencia de la cercanía. El modelo sí se evalúa en el z de cada medida.
5. **Dos preguntas distintas.** El p local de eBOSS frente a ΛCDM (eBOSS solo) y la comparación con lo esperado tras DESI (condicionada) no son lo mismo; la tabla de
   percentiles mezcla ambas a propósito para leerlas juntas.
6. **P21 y P22 son pruebas distintas y no son ciegas** (declaración arriba); P20 sí es independiente de los valores.
7. **Que eBOSS no confirme no refuta H1:** con potencia 0.33, dos de cada tres universos H1 darían exactamente esto. **Que coincida no demuestra física nueva ni es
   confirmación independiente si ρ_x > 0.**
8. Todo bajo ΛCDM plano con Planck comprimido y la aproximación gaussiana (linealización en p_ref). Ni este contraste ni LRG2 prueban el modelo OU ni σ_X > 0.

## Ficheros

`eboss_lrg_analisis.py` (modos `potencia` y `resultado`), `eboss_lrg_potencia.json` + `eboss_lrg_potencia_run.log` (commit b8d48a5, antes de los valores),
`eboss_lrg_resultado.json` + `eboss_lrg_resultado_run.log`, `data/PROCEDENCIA_EBOSS_DR16_LRG.md`. De la ejecución cancelada quedan sin seguimiento
`eboss_lrg__cancelado.py`, `eboss_lrg_potencia_cancelado.json` y `eboss_lrg_potencia_run_cancelado.log` (no se usan ni se commitean). Nada existente se modifica.
