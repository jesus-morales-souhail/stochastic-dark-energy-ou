# ¿Acota el CMB lo que DESI solo no puede? Núcleo completo con fondo de Planck (exploración v3)

**Predicción, escrita antes de correr `cmb_fondo.py`:** con BAO + CMB y el fondo marginalizado, la cota en θ = 1 queda
entre 3 y 8 × 10⁻²; en θ → 0 aparece una cota (sin CMB no la había) entre 1.5 y 4 × 10⁻²; en θ = 10 entre 2 y 5 × 10⁻²;
la mayor parte de la ganancia viene de θ*.

Exploración: no cambia `results/`, `figures/` ni ningún paper. Script `cmb_fondo.py`; salidas `cmb_fondo.{json,txt}`.
Tiempo: 2 min.

## Montaje

**Datos (15 filas).** Los 13 números de DESI DR2 con la covarianza oficial, 100θ* = 1.04110 ± 0.00031 y ω_m = 0.14301
± 0.0011 (Planck 2018, tabla 2), con la calibración de θ* de `escalon.py` (B) (÷1.000278): todo importado de `escalon.py`.
**Fondo:** ΛCDM plano (h, ω_m) con radiación, r_d(ω_m) y r_s(z*; ω_m) de `escalon.py`. ΛCDM con las 15 filas: H0 = 68.19,
ω_m = 0.14072, χ² = 17.24 (idéntico a `ESCALON.md`); con BAO solo, χ² = 10.28; BAO + θ*, 10.31; BAO + ω_m, 10.62.

**Residuo.** X(z) OU estacionario en u = ln(1+z), respuesta lineal completa de `verificacion.py`: δE² = X(z) − X(0)(1+z)³
(planitud con Ω_Λ fijo, el convenio de la tabla S), δlnD_H, δlnD_M integrada, δlnD_V = ⅔δlnD_M + ⅓δlnD_H, y la fila de
θ*: δlnθ* = −δlnD_M(z*) integrando X de 0 a z*. Malla uniforme en u con paso 2 × 10⁻³ hasta z* (3499 nodos; trapecio;
los z de DESI en el nodo más próximo, |δu| ≤ 8 × 10⁻⁴). La contribución de z > 5 a Var(δlnθ*) es 10⁻⁷ (muy por debajo del
1 %: el peso de X(z') cae como (1+z)⁻⁴·⁵). La fila de ω_m no responde a X (ω_m es un parámetro del universo temprano; el
residuo vive en la energía oscura). Covarianza σ_X² W K Wᵀ, K = e^{−θ|Δu|}, calculada en O(N) por θ con la recursión del
núcleo exponencial (frente a la matriz densa: 10⁻¹⁴). Σ_a W_ia reproduce −S_tipo de `perfil_13.s_type` a 10⁻³ en las 13
filas; en la fila de θ*, Σ_a W_a = −1.264 y el coeficiente de X(0) es −1.361: un X constante de 10⁻² mueve θ* un 1.3 %,
40σ de Planck.

**Variante "Λ absorbe X(0)":** δE² = X(z) − X(0) (la parte constante de un residuo de energía oscura es Λ; ω_m es la
densidad de materia). Σ_a W_ia = 0 en todas las filas: un X constante no hace nada. Coeficiente de X(0): 0.04–0.41 en los
BAO y −0.10 en θ*.

**Tratamientos del fondo:** perfilado (mínimos cuadrados sobre (h, ω_m) en cada (σ_X, θ)) y marginalizado (linealizado
alrededor del mejor ajuste ΛCDM de cada conjunto, prior plano). Cotas con ΔlnL = −1.92 y −1.35, θ ∈ {0.001, 0.1, 0.3, 1, 3,
10} y θ libre; conjuntos BAO, BAO + θ*, BAO + ω_m, BAO + CMB.

**Comprobación:** BAO solo, marginalizado, θ = 1: **0.4758** (VERIFICACION.md, núcleo completo, 13 números: 0.475); θ = 3:
0.287 (0.287); θ = 10: 0.203 (0.203). Perfilado, θ libre: 0.0232 (0.0232). La reparametrización (h, ω_m) con r_d(ω_m) y la
radiación no cambian nada.

## Resultados

**Fondo marginalizado, ΔlnL = −1.92** (σ_X al 95 %):

| datos | θ = 0.001 | 0.1 | 0.3 | 1 | 3 | 10 | θ libre |
|---|---|---|---|---|---|---|---|
| BAO | sin cota | sin cota | 0.87 | 0.476 | 0.287 | 0.203 | sin cota |
| BAO + θ* | sin cota | sin cota | 0.65 | 0.374 | 0.242 | 0.182 | sin cota |
| BAO + ω_m | **0.042** | 0.043 | 0.044 | **0.049** | 0.061 | **0.064** | 0.065 |
| BAO + CMB | **0.045** | 0.046 | 0.047 | **0.052** | 0.062 | **0.063** | **0.062** |

**Fondo marginalizado, ΔlnL = −1.35:**

| datos | θ = 0.001 | 0.1 | 0.3 | 1 | 3 | 10 | θ libre |
|---|---|---|---|---|---|---|---|
| BAO | sin cota | sin cota | 0.72 | 0.394 | 0.237 | 0.169 | sin cota |
| BAO + θ* | sin cota | 0.86 | 0.53 | 0.308 | 0.201 | 0.153 | sin cota |
| BAO + ω_m | 0.023 | 0.023 | 0.023 | 0.024 | 0.025 | 0.028 | 0.029 |
| BAO + CMB | 0.025 | 0.025 | 0.025 | 0.026 | 0.029 | 0.032 | 0.031 |

**Fondo perfilado, ΔlnL = −1.92** (con aviso, ver abajo):

| datos | θ = 0.001 | 0.1 | 0.3 | 1 | 3 | 10 | θ libre |
|---|---|---|---|---|---|---|---|
| BAO | 0.0232 | 0.0224 | 0.0213 | 0.0194 | 0.0181 | 0.0183 | 0.0232 |
| BAO + θ* | 0.0016 | 0.0016 | 0.0016 | 0.0015 | 0.0015 | 0.0015 | 0.0016 |
| BAO + ω_m | 0.0273 | 0.0261 | 0.0244 | 0.0218 | 0.0199 | 0.0205 | 0.0273 |
| BAO + CMB | 0.0410 | 0.0363 | 0.0325 | 0.0287 | 0.0279 | 0.0320 | 0.0410 |

(−1.35: BAO 0.0128–0.0102; BAO + θ* 0.0009–0.0008; BAO + ω_m 0.0150–0.0115; BAO + CMB 0.0230–0.0182.)

**Variante "Λ absorbe X(0)", BAO + CMB:**

| tratamiento, umbral | θ = 0.001 | 0.1 | 0.3 | 1 | 3 | 10 | θ libre |
|---|---|---|---|---|---|---|---|
| marginalizado, −1.92 | sin cota | sin cota | sin cota | 0.62 | 0.378 | 0.269 | sin cota |
| marginalizado, −1.35 | sin cota | sin cota | sin cota | 0.55 | 0.329 | 0.232 | sin cota |
| perfilado, −1.92 | sin cota | 0.12 | 0.072 | 0.044 | 0.032 | 0.025 | sin cota |
| perfilado, −1.35 | 0.54 | 0.055 | 0.033 | 0.020 | 0.015 | 0.012 | 0.54 |

(BAO solo en esta variante, marginalizado, −1.92: sin cota hasta θ = 0.1, 0.89 / 0.482 / 0.290 / 0.204 en θ = 0.3 / 1 / 3 / 10:
igual que el convenio principal, porque con el fondo marginalizado la diferencia entre los dos convenios, el modo X(0), está en
el espacio de plantillas.)


**Fondo marginalizado con el máximo global (corrección documental, 2026-11-12).** Donde una celda lleva "sust.", la cifra de las
tablas anteriores está referida a un máximo recortado en σ_X ≤ 0.05 y queda sustituida por esta. Las filas BAO + ω_m y BAO + CMB del
convenio principal no cambian (máximo en σ_X → 0 o en 0.004, dentro de la malla), así que las predicciones de esta nota no se ven afectadas.

ΔlnL = −1.92:

| datos | θ = 0.001 | 0.1 | 0.3 | 1 | 3 | 10 | θ libre |
|---|---|---|---|---|---|---|---|
| BAO | sin cota | sin cota | 0.807 (sust. 0.872 (-7.4 %)) | 0.452 (sust. 0.476 (-4.9 %)) | 0.281 (sust. 0.287 (-2.0 %)) | 0.202 | sin cota |
| BAO + θ* | sin cota | sin cota | 0.641 (sust. 0.647 (-0.9 %)) | 0.372 (sust. 0.374 (-0.6 %)) | 0.242 | 0.182 | sin cota |
| BAO + ω_m | 0.042 | 0.043 | 0.044 | 0.049 | 0.061 | 0.064 | 0.065 |
| BAO + CMB | 0.045 | 0.046 | 0.047 | 0.052 | 0.062 | 0.063 | 0.062 |

ΔlnL = −1.35:

| datos | θ = 0.001 | 0.1 | 0.3 | 1 | 3 | 10 | θ libre |
|---|---|---|---|---|---|---|---|
| BAO | sin cota | sin cota | 0.660 (sust. 0.724 (-8.8 %)) | 0.371 (sust. 0.394 (-5.8 %)) | 0.231 (sust. 0.237 (-2.4 %)) | 0.168 (sust. 0.169 (-0.6 %)) | sin cota |
| BAO + θ* | sin cota | 0.857 (sust. 0.862 (-0.6 %)) | 0.525 (sust. 0.531 (-1.1 %)) | 0.306 (sust. 0.308 (-0.8 %)) | 0.200 | 0.153 | sin cota |
| BAO + ω_m | 0.023 | 0.023 | 0.023 | 0.024 | 0.025 | 0.028 | 0.029 |
| BAO + CMB | 0.025 | 0.025 | 0.025 | 0.026 | 0.029 | 0.032 | 0.031 |

Variante "Λ absorbe X(0)", marginalizado (−1.92 y −1.35):

| datos | θ = 0.001 | 0.1 | 0.3 | 1 | 3 | 10 | θ libre |
|---|---|---|---|---|---|---|---|
| BAO + CMB, −1.92 | sin cota | sin cota | 0.867 (sust. sin cota) | 0.495 (sust. 0.624 (-20.6 %)) | 0.327 (sust. 0.378 (-13.5 %)) | 0.254 (sust. 0.269 (-5.6 %)) | sin cota |
| BAO + CMB, −1.35 | sin cota | sin cota | 0.736 (sust. sin cota) | 0.421 (sust. 0.547 (-23.0 %)) | 0.279 (sust. 0.329 (-15.4 %)) | 0.217 (sust. 0.232 (-6.5 %)) | sin cota |
| BAO, −1.92 | sin cota | sin cota | 0.818 (sust. 0.890 (-8.1 %)) | 0.456 (sust. 0.482 (-5.3 %)) | 0.283 (sust. 0.289 (-2.3 %)) | 0.203 (sust. 0.204 (-0.5 %)) | sin cota |
| BAO, −1.35 | sin cota | sin cota | 0.669 (sust. 0.740 (-9.5 %)) | 0.374 (sust. 0.399 (-6.3 %)) | 0.233 (sust. 0.240 (-2.7 %)) | 0.169 (sust. 0.170 (-0.6 %)) | sin cota |

**Por qué cambian:** `perfil_fondo.bound` mide la caída de 1.92 (o 1.35) desde el máximo de la verosimilitud sobre su malla de σ_X, que llega solo a 0.05; en las combinaciones con tensión el máximo está más allá de 0.05, y la cota medida desde ese máximo recortado sale más floja de lo debido. Descubierto en `MARGINAL_GUIADA.md`; recálculo en `maximos_globales.py` → `maximos_globales.json` (los valores históricos se reproducen antes con el máximo recortado). Sin CMB, o donde el máximo está en σ_X → 0 o dentro de la malla, nada cambia. Las cotas del 5 % no dependen de ningún máximo y no cambian.

## Lectura

1. **El CMB sí acota lo que DESI solo no puede, y lo hace ω_m, no θ*.** Con el fondo marginalizado, BAO + ω_m da
   0.042 (θ → 0), 0.049 (θ = 1), 0.064 (θ = 10); añadir θ* lo deja en 0.045, 0.052, 0.063. BAO + θ* solo apenas mueve
   la cota de BAO (0.476 → 0.374 en θ = 1; sin cota en θ → 0). El mecanismo: en el convenio de la tabla S, el modo
   constante de X es un desplazamiento de la densidad de materia a baja z (δΩ_m = −X(0) a h fijo); los BAO fijan
   Ω_m y h·r_d, y el prior de Planck fija ω_m y, con él, r_d y h, así que Ω_m = ω_m/h² queda clavada y el modo X(0)
   deja de ser absorbible por el fondo. θ* no ayuda porque su respuesta a X(0) (−1.36) es la misma que la de (h, ω_m)
   libres: 14 números, 2 parámetros, sin tensión (χ² = 10.31), y el fondo se lo come.
2. **La cota resultante es plana en θ (0.045–0.064) y es, otra vez, el modo común.** Con el fondo marginalizado y ω_m
   fijado, lo que acota el CMB es X(0); la parte local del residuo (coeficiente 0.12 → 0.014) sigue sin verse, por eso la
   cota apenas depende de θ y no mejora con θ grande (empeora: 0.045 → 0.063). Sigue sin ser una cota sobre la varianza
   estocástica del residuo entre bins; es una cota sobre la densidad de energía oscura de hoy a través de la de materia.
3. **El tratamiento perfilado con la fila de θ* es un artefacto de Occam.** BAO + θ* perfilado da 1.6 × 10⁻³ en todos
   los θ: con coeficiente 1.36 sobre un dato al 0.03 %, el ln det de la fila de θ* crece como ln[1 + (σ_X·1.41/3.1 × 10⁻⁴)²]
   y vale −1.92 en σ_X = 1.6 × 10⁻³ sin que el ajuste pierda nada (el fondo lo absorbe). Con BAO + CMB el mismo término
   compite con la tensión DESI–Planck (7 unidades de χ², `ESCALON.md`): el perfil tiene un máximo en σ_X ≈ 3 × 10⁻³
   (ΔlnL = +0.21) y la "cota" 0.041 es donde el castigo de inflar el error de θ* por un factor 190 supera la ganancia.
   Ninguno de los dos números es una cota sobre el residuo; con un dato tan preciso en el vector, solo el tratamiento
   marginalizado (ln det F compensa las direcciones del fondo) da algo interpretable.
4. **Con BAO + CMB marginalizado, el perfil prefiere σ_X ≈ 4 × 10⁻³ con ΔlnL = +2.0 (Δχ² ≈ 4) sobre σ_X = 0**, en todos
   los θ. Es la tensión DESI–Planck en ω_m/θ* absorbida por el modo X(0), lo mismo que absorbe el escalón de `ESCALON.md`
   (−6.6 con un parámetro) y la cuesta w₀wₐ. Las cotas de la tabla se miden desde ese máximo.
5. **Variante "Λ absorbe X(0)".** Si la parte constante del residuo es Λ (lo natural cuando ω_m es la densidad de materia
   medida por el CMB), un X constante no es observable y la cota con el fondo marginalizado es 0.62 (θ = 1), 0.38 (θ = 3),
   0.27 (θ = 10), sin cota en θ ≤ 0.3 [sustituidas por las de máximo global: 0.495, 0.327, 0.254, y 0.867 en θ = 0.3]: el CMB no añade nada (BAO solo: 0.48, 0.29, 0.20 [máximo global: 0.456, 0.283, 0.203]), porque lo que acota es la parte
   local, y esa no la ve ni DESI ni Planck comprimido. En esta lectura, los datos de hoy no acotan un residuo estocástico
   de energía oscura por debajo de 0.2–0.6.

## Predicción: qué se cumplió

- θ = 1, BAO + CMB marginalizado, entre 3 y 8 × 10⁻²: **cumplida** (0.052).
- θ → 0: aparece una cota (sin CMB no la había): **cumplida**; entre 1.5 y 4 × 10⁻²: **no, por poco** (0.045).
- θ = 10, entre 2 y 5 × 10⁻²: **no se cumplió** (0.063; la cota no mejora con θ, empeora).
- La mayor parte de la ganancia viene de θ*: **no se cumplió.** Viene de ω_m (BAO + ω_m: 0.042–0.064; BAO + θ*: sin
  cota en θ → 0, 0.37 en θ = 1). θ* solo aporta cuando ω_m ya está fijado, y entonces lo que aporta es la tensión.

**Cota de cabecera honesta para la v3.0:** con la respuesta lineal completa y el fondo (h, ω_m) marginalizado, DESI DR2
(13 números) más ω_m y θ* de Planck 2018 dan **σ_X < 6 × 10⁻² (95 %, ΔlnL = −1.92, θ libre; 3 × 10⁻² con −1.35; 4.5–6.3 ×
10⁻² según θ)** en la lectura en la que la parte constante del residuo desplaza la densidad de materia tardía (el convenio
de la tabla S), siendo ω_m de Planck lo que la fija; en la lectura en la que Λ absorbe la parte constante, estos datos no
acotan el residuo por debajo de 0.27–0.6 (θ = 10–1) y no hay cota en θ ≤ 0.3 [con el máximo global: 0.254–0.495 y 0.867 en θ = 0.3; sustituye a lo anterior]. La cota de 2.6 × 10⁻² de BAO solo con el
fondo fijo no sobrevive como cota de cabecera: es el término de Occam del modo de Ω_m y desaparece al liberar el fondo.
