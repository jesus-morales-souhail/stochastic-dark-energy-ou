# Auditoría v2.3: verificación del taquiónico, del núcleo OU y del umbral (exploración)

**Predicción, escrita antes de correr `verificacion.py`:**
1. Los dos puntos del taquiónico se confirman, y la cota resultante es una cota de amplitud, no una exclusión del modo.
2. El núcleo completo coincide en θ → 0 y afloja las cotas de θ ≥ 1 con fondo fijo entre ×1.5 y ×1.9.
3. ΔlnL = −1.92 cubre cerca del 97.5 %.

Exploración: no cambia `results/`, `figures/` ni ningún paper. Script `verificacion.py`; salidas `verificacion.{json,txt}` y
`diagnostico_nucleo.json`. Tiempo: 3.5 min con 16 procesos. Datos: los 7 α_DV de DESI DR2 con la covarianza oficial
proyectada (la de C2) y los 13 números con la 13 × 13; maquinaria de `fondo/perfil_fondo.py` y `fondo/perfil_13.py`
(tratamientos (a) fondo fijo, (b) perfilado, (c) marginalizado; malla de σ_X, búsqueda en θ, raíz por `brentq`).

## 1. Modelo taquiónico (C4)

**Qué hace `scripts/tachyonic_rank1_mle.py`.** C = C_meas + u uᵀ con u_i = σ₀ S_i exp(+t_i/t_c), t_i = tiempo de
retroceso del bin i (3.5 Gyr en z = 0.295, 11.0 Gyr en z = 2.33), σ₀ = 10⁻⁶¹ fijo, barrido en t_c ∈ [10⁻⁴, 10⁴] Gyr
(500 puntos). Línea 97: `i_near = argmin(|ΔlnL + 11.35|)` entre los puntos "activos" (amplitud > 10⁻⁶), y ese punto es
el que el preprint llama "representative active t_c".

**(i) Confirmado.** El −12.68 (Δχ² = 25.36) es el punto del barrido más cercano a −11.35, un objetivo escrito en el
código. No es un mínimo ni un máximo de nada: el barrido pasa de ΔlnL = −0.002 a −30 en un factor 1.25 de t_c, y el
punto elegido tiene una amplitud del modo de **2.7 × 10³ en unidades de α** (un residuo 2700 veces el valor medido).

| t_c (Gyr) | ΔlnL | amplitud máxima del modo (α) |
|---|---|---|
| 0.0663 | −30.05 | 9.6 × 10¹⁰ |
| 0.0688 | −24.04 | 2.4 × 10⁸ |
| 0.0714 | −18.26 | 7.3 × 10⁵ |
| **0.0741** | **−12.68** | **2.7 × 10³** |
| 0.0769 | −7.30 | 13 |
| 0.0798 | −2.13 | 0.072 |
| 0.0828 | −0.002 | 4.9 × 10⁻⁴ |

El mejor punto activo tiene ΔlnL = −1 × 10⁻⁷ (t_c = 0.086 Gyr, amplitud 4 × 10⁻⁶): el modelo no está "fuertemente
desfavorecido"; lo que ocurre es que σ₀ e^{t_max/t_c} se hace enorme para t_c < 0.081 Gyr (= 11.0 Gyr / ln(0.03/10⁻⁶¹)),
y cualquier residuo de amplitud 10³ es incompatible con cualquier dato. Es una cota de amplitud disfrazada.

**(ii) Confirmado.** Con t = tiempo de retroceso, exp(+t/t_c) crece hacia el pasado: en el t_c "representativo", el
cociente de amplitudes entre z = 2.33 y z = 0.295 es e^{7.5/0.074} = 9 × 10⁴³. En tiempo cósmico el modo **decrece**:
es lo contrario de un modo que crece por inestabilidad.

**Lo que sí es una cota.** Modo que crece en tiempo cósmico, u_i = A₀ S_i exp(−t_i/t_c), con A₀ la amplitud hoy (en α)
y t_c el tiempo de crecimiento; mismo S, misma covarianza, mismo criterio de perfil. El máximo de la verosimilitud está
en A₀ → 0 con ΔlnL = 0 en todos los t_c (ninguna preferencia, ninguna exclusión). Cotas al 95 % sobre A₀:

| t_c (Gyr) | amplitud en z = 2.33 / A₀ | A₀ < (ΔlnL = −1.92) | A₀ < (−1.35) |
|---|---|---|---|
| 0.1 | 2 × 10⁻⁴⁸ | sin cota (hasta A₀ = 1) | sin cota |
| 0.3 | 1 × 10⁻¹⁶ | sin cota | sin cota |
| 1 | 1.7 × 10⁻⁵ | sin cota | sin cota |
| 3 | 0.026 | 0.32 | 0.18 |
| 10 | 0.33 | 0.061 | 0.033 |
| 100 | 0.90 | 0.0288 | 0.0159 |
| ∞ (modo constante) | 1 | **0.02644** | 0.0146 |
| t_c perfilado (≥ 0.03 Gyr) | — | sin cota | sin cota |
| t_c perfilado (≥ 1 Gyr) | — | sin cota (ΔlnL = −0.25 en A₀ = 1) | sin cota |

- Con t_c ≤ 1 Gyr el modo es invisible: en z ≥ 0.295 (t ≥ 3.5 Gyr) su amplitud es < 10⁻⁵ A₀ y DESI no dice nada de
  A₀. Perfilando t_c, la cota desaparece. Solo hay cota para tiempos de crecimiento largos: A₀ < 0.06 (t_c = 10 Gyr),
  A₀ < 0.029 (100 Gyr), A₀ < 0.0264 (t_c → ∞).
- **Comparación con el OU en θ → 0:** el límite t_c → ∞ es exactamente el modo constante S_i X, es decir, el núcleo OU
  en θ → 0, y da la misma cota: 0.02644 = √(e^{3.84} − 1) σ_v (el término de Occam de `FONDO.md`; C2 = 0.02633 en
  θ = 10⁻³). Un modo coherente que crece despacio no se distingue del modo de fondo, y uno que crece deprisa no se ve.
  En ningún caso hay una exclusión del "crecimiento coherente": hay una cota de amplitud, y solo para t_c ≳ 3 Gyr.

## 2. Núcleo del OU: respuesta lineal completa

**Respuesta.** Plano, H0 fijo, fiducial Ω_Λ = 0.6847 (el de la tabla S): δE²(z) = X(z) − X(0)(1+z)³ (la planitud
convierte X(0) en −δΩ_m), δlnD_H = −δE²/(2E²), δlnD_M(z) = −(1/D̂_M) ∫₀ᶻ δE²/(2E³) dz', δlnD_V = ⅔ δlnD_M + ⅓ δlnD_H.
Escrita como δlnD_i = Σ_a W_ia X(z_a) sobre una malla de Simpson (897 nodos: z = 0 y 128 subintervalos por tramo entre
los z de DESI), la covarianza OU se propaga como C_ij = σ_X² (W K Wᵀ)_ij, K_ab = e^{−θ|x_a − x_b|}, x = ln(1+z).
El núcleo actual es C_ij = σ_X² S_i S_j e^{−θ|x_i − x_j|}, es decir, W_ia = −S_i δ(a, i): toda la respuesta del bin i
atribuida al valor de X en z_i.

**Comprobaciones.** Σ_a W_ia = ∂lnD_i/∂Ω_Λ|_plano = −S_i: con los 7 α coincide con la tabla S a 8.6 × 10⁻⁴ (la tabla
va con 3 decimales), con los 13 números a 4 × 10⁻⁸ frente a `s_type`. En θ = 10⁻³ las dos covarianzas coinciden al
0.24 % (7 α) y 0.10 % (13). **El núcleo completo coincide con el actual en θ → 0**, como se predijo.

**Descomposición de la respuesta** (7 α, `diagnostico_nucleo.json`): W_i = [coef. de X(0)] + [coef. local de X(z_i),
la parte D_H] + [parte integrada]:

| z | −S_i (actual) | coef. de X(0) | coef. local X(z_i) | parte integrada |
|---|---|---|---|---|
| 0.295 | 0.284 | 0.694 | −0.122 | −0.287 |
| 0.51 | 0.462 | 0.815 | −0.094 | −0.260 |
| 0.706 | 0.595 | 0.907 | −0.074 | −0.238 |
| 0.934 | 0.719 | 0.992 | −0.056 | −0.217 |
| 1.321 | 0.870 | 1.095 | −0.036 | −0.190 |
| 1.484 | 0.917 | 1.127 | −0.030 | −0.181 |
| 2.33 | 1.070 | 1.232 | −0.014 | −0.149 |

La respuesta de cada bin está dominada por X(0), el valor de hoy, común a todos los bins, con coeficiente 0.7–1.2; la
sensibilidad al valor local X(z_i) es 0.12 en z = 0.3 y 0.014 en z = 2.33 (en D_H sola, 0.37 → 0.04). El núcleo actual
pone todo el S_i (0.28–1.07) sobre el valor local, y por eso un residuo con θ ≥ 1 parece dispersión independiente entre
bins cuando en realidad es, casi todo, un modo común. Esa dirección X(0) está entera en el espacio de plantillas del
fondo plano (Ω_m, h·r_d): fracción de norma dentro del espacio, 1.000 (como S, 0.99996); la parte local, 0.9993.
A θ = 1 y σ_X = 0.015, el núcleo actual tiene cuatro modos relevantes (autovalores 11.3, 2.0, 0.69, 0.37, en unidades del
ruido; ln det = 4.69) y el completo es casi rango 1 (17.9, 0.11, 0.013, 0.005; ln det = 3.07): más varianza total, menos
castigo de Occam.

**Cotas al 95 % (ΔlnL = −1.92), núcleo actual / completo (cociente):**

7 α_DV:

| fondo | θ libre | θ = 0.001 | θ = 0.1 | θ = 1 | θ = 3 | θ = 10 |
|---|---|---|---|---|---|---|
| (a) fijo | 0.0263 / 0.0265 (1.00) | 0.0263 / 0.0265 (1.00) | 0.0210 / 0.0255 (1.21) | 0.0144 / 0.0215 (**1.49**) | 0.0129 / 0.0192 (**1.49**) | 0.0129 / 0.0186 (**1.44**) |
| (b) perfilado | 0.0263 / 0.0265 (1.00) | 0.0263 / 0.0265 (1.00) | 0.0213 / 0.0255 (1.20) | 0.0150 / 0.0216 (1.44) | 0.0134 / 0.0193 (1.44) | 0.0133 / 0.0186 (1.40) |
| (c) marginalizado | sin cota / sin cota | sin cota / sin cota | 0.140 / sin cota | 0.0461 / **0.497** (10.8) | 0.0295 / **0.319** (10.8) | 0.0220 / **0.244** (11.1) |

13 números:

| fondo | θ libre | θ = 0.001 | θ = 0.1 | θ = 1 | θ = 3 | θ = 10 |
|---|---|---|---|---|---|---|
| (a) fijo | 0.0231 / 0.0231 (1.00) | 0.0231 / 0.0231 (1.00) | 0.0183 / 0.0223 (1.22) | 0.0128 / 0.0192 (**1.50**) | 0.0116 / 0.0177 (**1.53**) | 0.0116 / 0.0180 (**1.55**) |
| (b) perfilado | 0.0231 / 0.0232 (1.00) | 0.0231 / 0.0232 (1.00) | 0.0185 / 0.0224 (1.21) | 0.0133 / 0.0193 (1.46) | 0.0120 / 0.0181 (1.50) | 0.0120 / 0.0182 (1.53) |
| (c) marginalizado | sin cota / sin cota | sin cota / sin cota | 0.112 / sin cota | 0.0377 / **0.475** (12.6) | 0.0249 / **0.287** (11.5) | 0.0193 / **0.203** (10.6) |

- Fondo fijo o perfilado, θ ≥ 1: el núcleo completo afloja **×1.40–1.55** (7 α: 1.44–1.49; 13: 1.46–1.55). La cota pasa
  de 1.3–1.5 × 10⁻² a 1.9–2.2 × 10⁻², y sigue siendo, como en θ → 0, un término de Occam sobre un modo casi común.
- Fondo marginalizado, θ ≥ 1: **×10.6–12.6**. La cota sobre un residuo que varía entre bins pasa de 2.2–4.6 × 10⁻² a
  **0.20–0.50**: con la respuesta correcta, DESI DR2 no acota un residuo estocástico de energía oscura al nivel del
  por ciento una vez que el fondo es libre, porque la parte que no absorbe el fondo (local, coeficiente ≤ 0.12 en D_V)
  es cinco a diez veces menor de lo que supone S_i. En θ = 0.1 no hay cota con el completo (había 0.14 / 0.11).
- La predicción (×1.5–1.9 con fondo fijo) se cumple en el borde inferior (1.44–1.55) y se quedó muy corta con el fondo
  marginalizado, que no predije: ahí el efecto es un orden de magnitud.

## 3. Umbral y cobertura

Cotas principales (7 α, núcleo actual) con ΔlnL = −1.92 (el del pipeline, 95 % bilateral de una χ²₁) y −1.35
(95 % unilateral con el parámetro en el borde):

| fondo | θ libre | θ = 0.001 | θ = 1 | θ = 10 |
|---|---|---|---|---|
| (a) fijo, −1.92 | 0.02633 | 0.02633 | 0.01444 | 0.01290 |
| (a) fijo, −1.35 | **0.01458** | 0.01458 | **0.01042** | 0.00994 |
| (b) perfilado, −1.92 | 0.02633 | 0.02633 | 0.01500 | 0.01327 |
| (b) perfilado, −1.35 | 0.01458 | 0.01458 | 0.01069 | 0.01021 |
| (c) marginalizado, −1.92 | sin cota | sin cota | 0.0461 | 0.0220 |
| (c) marginalizado, −1.35 | sin cota | sin cota | 0.0349 | 0.0171 |

En θ → 0 el cociente es √(e^{2.70} − 1)/√(e^{3.84} − 1) = 0.554 (ambas son el término de Occam 6.75 σ_v o 3.74 σ_v).

**Cobertura por inyección** (2000 universos α = 1 + S·X + n, n con la covarianza 7 × 7, semilla 20261105; X un sorteo OU
en (σ_t, θ_t) con σ_t = la cota correspondiente; cobertura = fracción de universos cuya cota, con el mismo tratamiento,
queda por encima de σ_t; fondo fijo):

| θ_t | σ_t | umbral | θ en el ajuste | cobertura |
|---|---|---|---|---|
| 0.001 | 0.02633 | −1.92 | fijo en θ_t | 99.45 ± 0.17 % |
| 0.001 | 0.02633 | −1.92 | libre | 99.15 ± 0.21 % |
| 0.001 | 0.01458 | −1.35 | fijo | 98.15 ± 0.30 % |
| 0.001 | 0.01458 | −1.35 | libre | 97.80 ± 0.33 % |
| 1 | 0.01444 | −1.92 | fijo | **97.60 ± 0.34 %** |
| 1 | 0.02633 | −1.92 | libre | 99.10 ± 0.21 % |
| 1 | 0.01042 | −1.35 | fijo | **95.30 ± 0.47 %** |
| 1 | 0.01458 | −1.35 | libre | 97.85 ± 0.32 % |

- Con θ_t = 1 y θ fijo, los dos umbrales hacen lo que dice la teoría asintótica: −1.92 cubre el 97.6 % (es el 97.5 %
  unilateral de un intervalo bilateral del 95 %) y −1.35 el 95.3 %.
- En θ → 0 los dos sobrecubren (99.5 % y 98.2 %): la cota ahí es un término de Occam sobre un modo que el ruido no
  excita, y el estadístico no es una χ²₁ en el borde. Con θ libre, igual (99.1 % y 97.8 %).
- Límites esperados sin señal (mediana, 16–84 %): θ libre, −1.92: 0.035 (0.027–0.064); −1.35: 0.021 (0.015–0.037);
  θ = 1, −1.92: 0.021 (0.015–0.033); −1.35: 0.015 (0.011–0.025). El 0.0263 observado está en el percentil ~15 de la
  distribución de solo ruido: DESI quedó algo mejor que lo típico (χ² = 4.5 para 7 bins).
- Predicción (3), "−1.92 cubre cerca del 97.5 %": cumplida donde el estadístico es regular (θ = 1 fijo, 97.6 %) y
  sobrecubre en el resto (97.8–99.5 %). El umbral −1.92 es conservador en todas las celdas; −1.35 es exacto en θ = 1 y
  conservador en θ → 0.

## Predicción: qué se cumplió

1. **Cumplida.** Los dos puntos del taquiónico se confirman en el código (línea 97; exp(+t/t_c) con t de retroceso). Con
   el signo correcto no hay exclusión: cota de amplitud A₀ < 0.026–0.06 para t_c ≥ 10 Gyr, ninguna para t_c ≤ 1 Gyr ni
   con t_c perfilado, e igual al OU en θ → 0 para t_c → ∞.
2. **Cumplida en parte.** Coincide en θ → 0 (ratio 1.00). Con fondo fijo afloja ×1.44–1.55, en el borde inferior de lo
   predicho (×1.5–1.9). Con fondo marginalizado, no predicho, afloja ×10.6–12.6, y es el resultado que importa.
3. **Cumplida** en θ = 1 con θ fijo (97.6 %); en las demás celdas −1.92 sobrecubre (97.8–99.5 %).

## Frases del PREPRINT y de CLAIMS que cada punto obliga a cambiar (no cambiadas aquí)

**Punto 1 (taquiónico):**
- PREPRINT l. 25 (abstract): "a coherent tachyonic growth model … is strongly disfavoured: ΔlnL ≈ −12.7 … (Δχ² ≈ +25)".
- PREPRINT l. 68 (guía de lectura): "Exclusions: coherent tachyonic growth (§4.1)".
- PREPRINT §4.1, l. 214 ("σ_X(t) = σ₀ e^{t/t_c}", el signo), l. 217 (la ecuación con −12.7, +25 y "representative active
  t_c"), l. 219 ("Any finite collapse / growth time t_c is disfavoured"), l. 222–224 ("Claim (model kill)" y la cita),
  l. 226 ("This exclusion is logically independent…").
- PREPRINT l. 263 (tabla §5.2): "Coherent / tachyonic growth | exponential | Excluded (Section 4)".
- PREPRINT l. 384 (conclusión 2): "Model kill: coherent tachyonic growth with rank-1 covariance is excluded (Δχ² ≈ +25)".
- CLAIMS C4 (l. 20): "Coherent tachyonic growth (rank-1 cov) is excluded: ΔlnL ≈ −12.7".
- COVER_LETTER l. 9 ("a controlled exclusion of a coherent tachyonic growth model") y l. 11 ("a rank-1 coherent growing
  mode is strongly disfavoured (Δχ² ≈ +25)"); README l. 36 ("A coherent tachyonic growing mode is excluded when it is
  active"); manuscript/README.md l. 32 ("Tachyonic coherent growth exclusion").
- `results/tachyonic_rank1/` y `scripts/tachyonic_rank1_mle.py` (signo de t, σ₀ fijo, línea 97) tendrían que
  regenerarse con el modo que crece en tiempo cósmico y A₀ libre.

**Punto 2 (núcleo):**
- PREPRINT l. 97 ("enters the BAO residual covariance through the sensitivity kernel S_i … δα_i = S_i X(z_i)") y la
  ecuación l. 101 (C^OU_ij = S_i S_j σ_X² e^{−θΔx}): el núcleo debe ser W K Wᵀ, o decirse que S_i X(z_i) es la
  aproximación de modo común y solo vale en θ → 0.
- PREPRINT l. 22 (abstract): "A residual that varies between redshifts (θ ≥ 1) is bounded at σ_X < 1.5 × 10⁻² with the
  background fixed or profiled, and < 4.6 × 10⁻² with it marginalised" → 1.9–2.2 × 10⁻² y 0.20–0.50.
- PREPRINT §3.2: la tabla l. 169–173 (filas θ = 0.1, 1, 10) y el texto l. 175 ("The bound on a stochastic residual that
  varies between redshifts (θ ≳ 1) is 1.3–1.5 × 10⁻² … and 2.2–4.6 × 10⁻² marginalised"; los factores 1.11–1.25 de los
  13 números; "(i) … 11–25 % above the marginalised limits, which it reproduces without using the kernel" deja de
  cuadrar, porque la cota por bondad de ajuste sí usa S_i X(z_i) en los universos; "(ii) A < 0.013–0.028" del armónico
  usa el mismo núcleo local).
- Todos los números derivados "at the marginalised θ = 1 ceiling (4.611 × 10⁻², ×1.751)": PREPRINT l. 197 (§3.5), l. 201
  y 203 (PICOS, RAPIDO a las cotas de θ = 3, 10), l. 254 (§5.1, G_DESI = 3.92 × 10⁵⁹), l. 272 (§5.3, N_eff ≈ 470), l. 314
  (§6.2, slip 8.67 × 10⁻², ≈ 1.7×); CLAIMS C7 (último inciso), C8 ("fast residuals (θ = 3, 10)" a sus cotas) y C9
  (A < 0.013–0.028).
- PREPRINT §7, l. 322 ("Using the full 13-number vector tightens the limits by a factor 1.11–1.25") y l. 323 (la lectura
  de que los límites de θ ≥ 1 sí son cotas sobre un residuo).
- PREPRINT §9, l. 382 (conclusión 1, los mismos números) y l. 390 ("separating θ = 1 from θ = 10⁻³ needs σ_X ≥ 2.2 × 10⁻²,
  above the θ = 1 limit of 1.44 × 10⁻²; with the background marginalised a narrow window opens (2.8–4.6 × 10⁻²)";
  también el 3.8 % de fσ₈ y el 2.2 % en D_L se calcularon con δΩ_Λ = S-modo, que es el caso θ → 0 y no cambia).
- CLAIMS C2 ("For a residual varying between redshifts (θ ≥ 1): < 1.5 × 10⁻² …, < 4.6 × 10⁻²") y N14 ("the θ ≥ 1 limits of
  §3.2 are [bounds on stochastic variance]"); README l. 36 ("for a residual that varies between redshifts (θ ≥ 1) the
  limits are 1.3–1.5 × 10⁻² … and 2.2–4.6 × 10⁻²") y `papers/resume.txt`.
- `fondo/FONDO.md`, `TRECE.md`, `cuesta/GOF.md`, `cuesta/CUESTA.md`, `fondo/FISHER_DISENO.md`, `RESOLVER_THETA.md`,
  `sn/SN.md` usan el mismo δα_i = S_i X(z_i) y heredan el cambio en θ ≥ 1 (no en θ → 0).

**Punto 3 (umbral):**
- PREPRINT l. 18 (abstract, "the profile likelihood over θ (ΔlnL = −1.92 …) gives the 95% confidence upper limit") y
  l. 166 (§3.2); CLAIMS C2 ("95% CL, profile over θ"); README l. 30 ("Profile 95% CL"). O bien se dice que −1.92 es el
  intervalo bilateral y cubre ≥ 97.5 % (97.6–99.5 % por inyección), o bien se pasa a −1.35 (1.46 × 10⁻² con θ libre,
  1.04 × 10⁻² en θ = 1, fondo fijo; 95.3 % en θ = 1) y cambian todos los números derivados de C2. Los números no son
  erróneos con −1.92; su etiqueta ("95 %") es conservadora.
