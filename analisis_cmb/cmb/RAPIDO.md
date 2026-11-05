# El temblor rápido (θ = 3, 10) en Planck y en los picos (exploración)

**Predicción, escrita antes de correr `rapido.py`:**
- p95 del Δχ² en TT ℓ < 30 + lente < 0.5 para θ = 10 a σ_X = 2.2 × 10⁻², con |w + 1| de orden 0.3–0.5.

Exploración: no cambia ningún resultado ni ningún paper. Script `rapido.py` (venv de `analisis_cmb`);
salidas `rapido.{json,txt}`. Hueco de §3.5: allí solo se corrieron θ = 0.001 y θ = 1 (θ = 20 se descartó
por estar dominado por el suavizado).

## Montaje

**El de §3.5 / `realizaciones.py` y de `picos.py`:** camino OU en la malla fina, suavizado gaussiano Δ, tabla
w(a) de 2000 puntos en `DarkEnergyPPF`, H0 resuelto para que 100θ* = 1.041198 (ΛCDM), las 60 primeras
realizaciones de la semilla 20261003 (las mismas de `PICOS.md`).
- **Casos:** θ ∈ {3, 10} × σ_X ∈ {1.290 × 10⁻², 2.204 × 10⁻²} (cotas de θ = 10 con fondo fijo y
  marginalizado, `FONDO.md` / `perfil_fondo.json`).
- **Planck:** Δχ² = χ²(curva) − χ²(ΛCDM) frente a TT 2 ≤ ℓ < 30 y frente a las bandas φφ (L = 8–400), con
  los errores de Planck 2018 (`estadisticos.chi2`). Además el ℓ (o la banda) del máximo de |ΔD_ℓ/σ_ℓ|, con
  Δ = curva − ΛCDM y σ el error de Planck.
- **Picos:** la distancia de `PICOS.md` (TT/TE/EE con lente, ℓ ∈ [30, 2500], covariancia 3 × 3 por ℓ,
  f_sky = 0.6, ruido blanco de 33 µK·arcmin), y la misma sin lente.
- **w_eff:** w = −1 − X′_s/(3(Ω_Λ + X_s − X_s(1))) sobre la tabla; se da el p95 de |w + 1| por realización
  (mediana y p95 sobre realizaciones).

**¿Resuelve Δ = 0.05 un temblor con θ = 10?** No. Un suavizado gaussiano de desviación Δ conserva una
fracción e^{(θΔ)²} erfc(θΔ) de la varianza de un OU estacionario:

| θ | Δ = 0.05 | Δ = 0.025 | Δ = 0.0125 |
|---|---|---|---|
| 1 | 0.946 | 0.972 | 0.986 |
| 3 | **0.851** | 0.921 | 0.959 |
| 10 | 0.616 | 0.770 | **0.873** |

En la malla (60 realizaciones) sale 0.83 (θ = 3, Δ = 0.05) y 0.60 / 0.76 / 0.87 (θ = 10). Con Δ = 0.05, el
suavizado se lleva el 38 % de la varianza del temblor de θ = 10, así que **θ = 10 se corre con Δ = 0.0125**
(87 %, el mismo θΔ que θ = 3 con 0.05; la memoria 1/θ = 0.1 queda en 8Δ). CAMB reproduce nuestra ρ_DE a
3.6 × 10⁻⁴ (relativo, máximo) con ese Δ (1.2–2 × 10⁻⁵ con Δ = 0.05). Como sensibilidad, θ = 10 se corre
también con Δ = 0.05 y 0.025, y θ = 3 con 0.025.

**Un artefacto del halofit, encontrado al hacer la comprobación de Δ.** Con Takahashi (el halofit de §3.5),
el Δχ² de la lente a θ = 10 crecía sin converger al bajar Δ: p95 = 0.077 / 0.32 / 0.65 con Δ = 0.05 / 0.025 /
0.0125 (σ_X = 2.2 × 10⁻²). En la realización 0 a Δ = 0.0125, C_L^φφ sube un 1.4 % en todo L = 8–400 (0.2 %
en L = 8, 1.8 % en L = 400), **mientras el potencial de Weyl, P(k, z) y σ₈(z) no cambian** (< 5 × 10⁻⁴ en 26 z
entre 0.5 y 1, y en z = 0, 2, 3, 5), ni las distancias (χ(z = 2) a 3 × 10⁻⁵). Con `AccuracyBoost = 2` y con
una tabla de 8000 puntos sale lo mismo. **Con lente lineal, el efecto desaparece** (−0.024 en los dos Δ) y
HMcode (mead2020) falla directamente ("Integration timed out"). Es la dependencia explícita en w(a) de la
fórmula de Takahashi, evaluada con un w que oscila ±0.4 en Δln a = 0.0125: un artefacto del ajuste, no
física. Por eso **cada caso se corre con Takahashi y con lente lineal**; la lectura física es la lineal.
A θ = 3 con Δ = 0.05 el artefacto ya resta un 30 % al Δχ² de φφ (0.047 frente a 0.067).

## Resultados

**Casos principales (60 realizaciones; mediana / p95):**

| θ, Δ | σ_X | lente | Δχ² TT ℓ < 30 | Δχ² φφ | **TT + φφ, p95** | picos con lente | picos sin lente, p95 | p95 \|w + 1\| (med / p95) | σ(H0) |
|---|---|---|---|---|---|---|---|---|---|
| 3, 0.05 | 1.290e-2 | Takahashi | −0.0000 / 0.038 | −0.001 / 0.027 | **0.061** | 8.1e-5 / 1.1e-3 | 1.6e-5 | 0.061 / 0.071 | 0.32 |
| 3, 0.05 | 1.290e-2 | lineal | −0.0000 / 0.038 | −0.003 / 0.039 | **0.070** | 5.0e-5 / 8.4e-4 | 1.6e-5 | 0.061 / 0.071 | 0.32 |
| 3, 0.05 | 2.204e-2 | Takahashi | 0.0001 / 0.065 | −0.002 / 0.047 | **0.105** | 3.1e-4 / 3.6e-3 | 4.6e-5 | 0.104 / 0.123 | 0.54 |
| 3, 0.05 | 2.204e-2 | lineal | 0.0001 / 0.065 | −0.005 / 0.067 | **0.122** | 1.5e-4 / 2.5e-3 | 4.6e-5 | 0.104 / 0.123 | 0.54 |
| 10, 0.0125 | 1.290e-2 | Takahashi | −0.0006 / 0.033 | −0.003 / 0.54 | 0.54 (artefacto) | 1.2e-3 / 4.3e-2 | 8.9e-6 | 0.228 / 0.251 | 0.37 |
| 10, 0.0125 | 1.290e-2 | **lineal** | −0.0007 / 0.033 | −0.001 / 0.032 | **0.066** | 3.0e-5 / 3.4e-4 | 8.9e-6 | 0.228 / 0.251 | 0.37 |
| 10, 0.0125 | 2.204e-2 | Takahashi | −0.0008 / 0.056 | 0.010 / 0.65 | 0.63 (artefacto) | 2.0e-3 / 6.1e-2 | 3.2e-5 | 0.390 / 0.434 | 0.63 |
| 10, 0.0125 | 2.204e-2 | **lineal** | −0.0008 / 0.056 | −0.001 / 0.055 | **0.115** | 8.6e-5 / 1.0e-3 | 3.2e-5 | 0.390 / 0.434 | 0.63 |

Para comparar, §3.5 / `PICOS.md` con θ = 1 y σ_X = 2.633 × 10⁻² (Takahashi, Δ = 0.05): TT p95 = 0.067,
φφ p95 = 0.054, picos p95 = 3.5 × 10⁻³, σ(H0) = 0.48.

**Sensibilidad a Δ (p95; σ_X = 2.204 × 10⁻²):**

| θ | Δ | varianza conservada | p95 \|w + 1\| | TT | φφ Takahashi | φφ lineal | TT + φφ lineal | picos Takahashi | picos lineal |
|---|---|---|---|---|---|---|---|---|---|
| 10 | 0.05 | 0.60 | 0.150 | 0.057 | 0.077 | 0.053 | 0.111 | 4.7e-3 | 1.0e-3 |
| 10 | 0.025 | 0.76 | 0.248 | 0.054 | 0.32 | 0.054 | 0.115 | 3.4e-2 | 1.0e-3 |
| 10 | 0.0125 | 0.87 | 0.390 | 0.056 | 0.65 | 0.055 | 0.115 | 6.1e-2 | 1.0e-3 |
| 3 | 0.05 | 0.83 | 0.104 | 0.065 | 0.047 | 0.067 | 0.122 | 3.6e-3 | 2.5e-3 |
| 3 | 0.025 | 0.91 | 0.157 | 0.067 | 0.16 | 0.068 | 0.123 | 1.9e-2 | 2.5e-3 |

Con lente lineal, nada depende de Δ (TT, φφ, picos: ≤ 4 % entre Δ = 0.05 y 0.0125), aunque el suavizado
cambie la varianza del 60 al 87 % y |w + 1| de 0.15 a 0.39. Con Takahashi, φφ y los picos crecen con 1/Δ.
σ(H0) tampoco depende de Δ (0.50–0.63).

**Dónde está el máximo de |ΔC_ℓ/σ_ℓ|** (frente a los errores de Planck; mediana sobre realizaciones):

| caso (lineal) | máx \|Δ/σ\| en TT | ℓ del máximo (p16 / p50 / p84) | máx \|Δ/σ\| en φφ | banda del máximo, L_av | fracción con el máximo en φφ |
|---|---|---|---|---|---|
| θ = 3, 1.290e-2 | 0.004 σ | 2 / 2 / 6 | 0.007 σ | 64 / 64 / 106 | 0.82 |
| θ = 3, 2.204e-2 | 0.006 σ | 2 / 2 / 6 | 0.012 σ | 64 / 64 / 106 | 0.82 |
| θ = 10, 1.290e-2 | 0.004 σ | 2 / 2 / 6 | 0.006 σ | 28 / 64 / 64 | 0.70 |
| θ = 10, 2.204e-2 | 0.006 σ | 2 / 2 / 6 | 0.009 σ | 28 / 64 / 64 | 0.70 |

- El máximo está **en la lente** (bandas L ≈ 30–100, el 70–82 % de las veces), no en el ISW: en TT el
  máximo cae en ℓ = 2, la cola del ISW tardío, pero es más pequeño. Ambos son ≤ 0.012σ de Planck: el Δχ²
  de 0.03–0.07 contra los datos viene casi todo del término cruzado con los residuos del dato (`RESULTADO.md`,
  sección 3), no de la señal propia.
- En los picos (ℓ ≥ 30), con lente lineal el máximo está en ℓ = 30 (borde, ISW) para θ = 10, y en 830
  (valle de EE) o 1613 (flanco de TT) para θ = 3, como en `PICOS.md`; sin lente, siempre en ℓ = 30.
- Las 60 realizaciones dejan 100θ* a ≤ 5 × 10⁻⁸ del objetivo; ninguna da ρ_DE ≤ 0.

## Por qué

1. **El CMB integra: un temblor con memoria 0.1 en ln a se promedia.** Con lente lineal, pasar de θ = 3 a
   θ = 10 baja la distancia en los picos un factor 2.5–3 (2.5e-3 → 1.0e-3) y el Δχ² de φφ un 15–20 %, a
   σ_X igual: cuanto más rápido el temblor, menos ve el CMB. En TT ℓ < 30 ocurre lo mismo (0.065 → 0.056).
2. **Por eso el resultado no depende de Δ aunque Δ = 0.05 no resuelva el temblor.** El suavizado se lleva
   el 38 % de la varianza de θ = 10, pero esa varianza está en escalas Δln a < 0.05 que ninguna observable
   del CMB resuelve. Lo que sí depende de Δ es w_eff: |w + 1| crece como Δ^{−1/2} (0.15 → 0.25 → 0.39), porque
   un OU no es derivable y w es su derivada. Es lo que dice `NOTAS_ANALISIS_CMB.md`: Δ es un parámetro del
   modelo, no un detalle numérico; aquí se ve que el CMB no lo nota, pero w sí.
3. **El artefacto de Takahashi es consecuencia directa de lo anterior:** su fórmula usa w(a) en cada z, y un
   w que oscila ±0.4 cambia el ajuste no lineal aunque la densidad y el crecimiento no cambien. Cualquier
   cálculo con halofit y un temblor rápido hereda este sesgo; a θ = 1 (|w + 1| ≈ 0.05) es despreciable, a
   θ = 3 con Δ = 0.05 ya es del 30 % en φφ.
4. **La huella que queda es la de siempre:** σ(H0) de 0.3–0.6 km/s/Mpc por la renormalización en a = 1, y un
   residuo de lente en L ≈ 30–100 de ≤ 0.01σ.

## Conclusión

Un temblor con θ = 3 o θ = 10 a las cotas de DESI (1.3 y 2.2 × 10⁻²) es invisible para Planck: p95 del Δχ²
en TT ℓ < 30 + lente de 0.066–0.12 (lente lineal) y en los picos ≤ 2.5 × 10⁻³. Es menor que lo que dejaba
θ = 1 en §3.5. El único número que supera 0.5 es el de Takahashi a Δ = 0.0125, y es un artefacto de la
fórmula del halofit, no del universo.

## Predicción

- **p95 Δχ² < 0.5 en TT ℓ < 30 + lente para θ = 10 a 2.2 × 10⁻²: cumplida con la física (0.115 con lente
  lineal, estable en Δ), rota por el número bruto del montaje de §3.5 (0.63 con Takahashi a Δ = 0.0125),
  que es el artefacto descrito arriba.** Con Takahashi y Δ = 0.05, el montaje literal de §3.5, también se
  cumple (0.105).
- **|w + 1| de orden 0.3–0.5: cumplida en el Δ elegido** (p95 = 0.39 mediana, 0.43 p95, con Δ = 0.0125),
  pero el número es de Δ, no de θ: con Δ = 0.05 habría sido 0.15.
