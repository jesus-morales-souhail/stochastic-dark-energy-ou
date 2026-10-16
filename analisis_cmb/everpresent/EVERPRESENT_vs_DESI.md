# Everpresent Λ frente a DESI DR2 BAO

Referencias: [A] Ahmed, Dodelson, Greene, Sorkin, astro-ph/0209274; [Z] Zwane, Afshordi, Sorkin,
arXiv:1703.06265 (lectura en `LECTURA_EVERPRESENT.md`). Scripts y salidas, todos en esta carpeta:

- `ahmed_sim.py` → `resultados/ahmed_sim.json`, `resultados/ahmed_sim.txt`
- `profile_desi_everpresent.py` (copia de `scripts/profile_sigma_x_desi.py`) con
  `desi_dr2_data.py` (copia de `scripts/desi_dr2_data.py`, sin escrituras a disco, más el
  kernel de [Z]) → `resultados/profile_desi_everpresent.json`, `resultados/profile_desi_everpresent.txt`

Datos DESI: los 7 α isotrópicos de DR2, residuo α − 1, covarianza `full_gaussian_bao_projected_13x13`,
criterio del 95 %: ΔlnL = −1.92 desde el máximo del perfil (1 grado de libertad). Los límites se dan
en la malla de 48 puntos del perfil vigente ("malla") y refinados por raíz entre nodos ("refinado").

---

## 1. Reproducción de Ahmed+2004 (Sec. III): no reproducible con las convenciones declaradas

Implementación literal de [A] ec. (5): S_{i+1} = S_i + α ξ √δN_i, ρ_Λ = S/N, con N = V el
4-volumen del cono de luz pasado (ec. 3, integrado sin cancelaciones por momentos
J_k = ∫ a³(η − η')^k dt), y la Friedmann de la ec. (6) con radiación y materia de Planck 2018.
Arranca en a = 1e-32 con S_0 = 0, y un run muere si ρ_tot ≤ 0. Se usan α = 0.02, 200 realizaciones
y la semilla 20261016. Si se escribe todo en unidades de H0, las potencias de ℓ_P se cancelan, y las
dos convenciones solo difieren en el factor K del término de energía oscura,
u = K α Ŝ/V̂: K = 1/3 en (R), 8π/3 en (NR).

Comprobaciones: V sin Λ sale igual que por cuadratura directa (0.3434609 frente a 0.3434610 (c/H0)⁴),
y el resultado no cambia con el paso (h = 0.01 y h = 0.003 en ln a).

| convención | completan | Ω_Λ(a=1) p16 / p50 / p84 | σ(Ω_Λ hoy) | fichero |
|---|---|---|---|---|
| (R) 8πG = ħ = c = 1, la declarada en [A] | 200/200 | −0.0352 / +0.0000 / +0.0351 | 0.0352 | `resultados/ahmed_sim.json` → `conventions.R` |
| (NR) G = ħ = c = 1 | 0/200 (mueren en a ≈ 1.05–2.5 × 10⁻³² , p16–p84) | — | — | `conventions.NR` |
| (R) desde a = 1/1100 | 200/200 | | 0.0331 | `conventions.R_from_a0.000909` |
| (NR) desde a = 1/1100 | 0/200 | | — | `conventions.NR_from_a0.000909` |
| [A] Sec. III | "about a 1 in 3 chance" | "amplitude of order 0.5" | | |

Barrido diagnóstico K = f × (1/3), con α = 0.02 (`resultados/ahmed_sim.json` → `f_scan`):

| f | 1 | 2 | 3 | 4 | 5 | 6 | 6.5 | **7** | 7.5 | 8 | 12 | 8π |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| completan /200 | 200 | 200 | 200 | 200 | 198 | 164 | 110 | **63** | 25 | 5 | 0 | 0 |
| σ(Ω_Λ hoy) | 0.035 | 0.071 | 0.107 | 0.144 | 0.179 | 0.204 | 0.233 | **0.273** | 0.322 | 0.531 | — | — |

Los dos resultados de Sec. III ("1 de cada 3 completa", amplitud de Ω_Λ ~0.5) se reproducen solo con
un factor efectivo f ≈ 7 sobre la convención declarada, equivalente a un V unas 50 veces menor que el
del cono pasado de su ec. (3). Ese factor no se usa como "convención de Ahmed"; queda registrado.
La integral de [Z] ec. (2.6) (interior ∫_0^{t'}) va en sentido contrario: da un V 230 veces mayor
(28.41 frente a 0.1235 (c/H0)⁴, `n_count_salida.txt`), es decir una fluctuación menor, así que no
explica el factor.

---

## 2. Forma del proceso de Ahmed, independiente de la normalización

Este es el resultado que sí vale de la parte 1. Con la convención R (régimen lineal, sin muertes),
el proceso de [A] escala con su amplitud, así que los cocientes no dependen ni de α ni de la
normalización del 4-volumen. Fichero: `resultados/ahmed_sim.json` → `process`.

| caso | n | σ(Ω_Λ hoy) | σ(ΔΩ_Λ), 0.3<z<2.3 | r | θ ajustado |
|---|---|---|---|---|---|
| R, α = 0.02 | 200 | 0.0352 | 0.0466 | **1.323** | **2.720** |
| R, α = 0.01 | 200 | 0.0176 | 0.0233 | 1.324 | 2.701 |
| f = 7, α = 0.02 (no lineal, solo los runs que completan) | 63 | 0.2735 | 0.3907 | 1.429 | 2.933 |

Definiciones:
- ΔΩ_Λ(z) = Ω_Λ(z) − Ω_Λ(0).
- σ(ΔΩ_Λ) es la desviación típica conjunta sobre realizaciones y sobre 50 puntos uniformes en
  ln(1+z) entre z = 0.3 y z = 2.3.
- r = σ(ΔΩ_Λ en rango DESI) / σ(Ω_Λ hoy).

Si en vez de la diferencia con hoy se toma el nivel Ω_Λ(z), el cociente es
σ(Ω_Λ(z))/σ(Ω_Λ hoy) = 1.059 (`process.R_alpha0.02.r_level`).

Matriz de correlación de Ω_Λ entre los 7 z_eff de DESI (R, α = 0.02;
`process.R_alpha0.02.corr_7x7`):

| z_eff | 0.295 | 0.510 | 0.706 | 0.934 | 1.321 | 1.484 | 2.330 |
|---|---|---|---|---|---|---|---|
| 0.295 | 1 | 0.658 | 0.357 | 0.267 | 0.122 | 0.153 | 0.050 |
| 0.510 | 0.658 | 1 | 0.659 | 0.521 | 0.321 | 0.343 | 0.140 |
| 0.706 | 0.357 | 0.659 | 1 | 0.740 | 0.438 | 0.411 | 0.235 |
| 0.934 | 0.267 | 0.521 | 0.740 | 1 | 0.576 | 0.535 | 0.290 |
| 1.321 | 0.122 | 0.321 | 0.438 | 0.576 | 1 | 0.855 | 0.388 |
| 1.484 | 0.153 | 0.343 | 0.411 | 0.535 | 0.855 | 1 | 0.439 |
| 2.330 | 0.050 | 0.140 | 0.235 | 0.290 | 0.388 | 0.439 | 1 |

Ajuste por mínimos cuadrados de exp(−θ|Δln(1+z)|) a los 21 elementos fuera de la diagonal:
θ* = 2.720, con residuo rms 0.050.

Límite DESI DR2 del kernel OU con θ fijo a θ* = 2.720 (`resultados/profile_desi_everpresent.json`
→ `T2_ou_theta_fit`): σ_X < **1.294 × 10⁻²** (refinado; 1.256 × 10⁻² en la malla). El máximo del
perfil está en σ_X → 0.

**Inferencia.** Si Everpresent Λ exige σ(Ω_Λ hoy) ≈ 0.7 ([Z] Sec. 3: "the typical value of |Ω_de|
is 0.7"), la variación esperada en el rango de DESI es

  σ(ΔΩ_Λ)_DESI = r × σ(Ω_Λ hoy) = 1.323 × 0.7 = 0.926,

y el factor frente al límite DESI con la misma forma de correlación es

  F = r × 0.7 / σ_X,95(θ*) = 0.926 / 1.294 × 10⁻² = **72**  (74 con el límite de la malla).

Con el nivel en lugar de la diferencia: 1.059 × 0.7 / 1.294 × 10⁻² = 57. En el régimen no lineal
(f = 7) el cociente sube a r = 1.43, lo que daría F = 77.

---

## 3. Cota de DESI DR2 sobre el Modelo 2 de Zwane

### 3.1 Prueba de la copia (`resultados/profile_desi_everpresent.json` → `T0`, `T1`)

| kernel OU | límite σ_X 95 % (malla) | refinado | vigente |
|---|---|---|---|
| θ perfilado en (10⁻³, 20) | **2.506316 × 10⁻²** | 2.633 × 10⁻² | 2.506316 × 10⁻² (`results/profile_sigma_x/profile_sigma_x.json`) |
| θ = 1 fijo (información) | 1.256 × 10⁻² | 1.444 × 10⁻² | — |

La copia reproduce exactamente el valor vigente. Nota: el 2.506 × 10⁻² vigente es el último nodo
de la malla que pasa el umbral; la raíz entre nodos está en 2.633 × 10⁻².

### 3.2 Kernel de [Z] ec. (2.9)

C_ij += S_i S_j α² exp[−μ (λ_i − λ_j)²/2], con λ = ln a = −ln(1+z) y los mismos S_z que el
kernel OU. Se ignora el tanh de (2.10), que es irrelevante para la cota a α pequeño; abajo se
comprueba su efecto en α = 0.8824. Ω̂_de se trata como residuo sobre el fiducial ΛCDM, igual que
X en el repo (en [Z] la media de Ω_de es 0, no Ω_Λ0).

| modo | mejor α | α 95 % (malla) | α 95 % refinado | ΔlnL(α = 0.8824) − lnL(α = 0) | Δχ² | √Δχ² | fichero |
|---|---|---|---|---|---|---|---|
| μ = 0.9804 fijo | → 0 (10⁻⁶, borde) | 1.688 × 10⁻² | **1.849 × 10⁻²** | −11.58 | 23.15 | **4.8σ** | `Z1_zwane_mu_fixed` |
| μ libre en [0.1, 10] | → 0 (10⁻⁶, borde) | 2.028 × 10⁻² | **2.407 × 10⁻²** | −8.57 (μ = 0.1) | 17.15 | **4.1σ** | `Z2_zwane_mu_free` |
| μ = 0.9804, con tanh exacto | — | — | — | −10.78 | 21.56 | 4.6σ | `Z1_zwane_mu_fixed.alpha_zwane_with_tanh` |

La fila "con tanh" usa la verosimilitud gaussiana con la covarianza exacta de tanh(Ω̂) (Gauss–Hermite);
la σ efectiva baja a 0.591.

**μ libre: el perfil prefiere μ → 0.1, el borde inferior.** Para α ≳ 10⁻² el μ que maximiza la
verosimilitud es siempre 0.1 (`Z2_zwane_mu_free.mu_along_profile`). Para α ≲ 3 × 10⁻³ el μ es
indiferente (el kernel no pesa frente a la covarianza de medida) y el optimizador cae en 10.
Correlación larga significa que el ruido se convierte en un desplazamiento casi constante en el
rango de DESI, que la normalización absorbe. Es la única salida del modelo, y aun así α = 0.8824
queda a 4.1σ y el límite es 2.4 × 10⁻².

---

## Conclusiones

(i) Everpresent Λ en la forma de Ahmed+2004 (Modelo 1) exige una variación de Ω_Λ en 0.3<z<2.3 de
≈ r × 0.7 = 1.323 × 0.7 = 0.93 (`resultados/ahmed_sim.json`, `process.R_alpha0.02.r`); DESI DR2 con
la covarianza completa limita esa variación, con la misma correlación (θ* = 2.72), a
< 1.29 × 10⁻² (`resultados/profile_desi_everpresent.json`, `T2_ou_theta_fit.upper95_refined`):
factor 72, independiente de la normalización del 4-volumen.

(ii) El mejor ajuste de Zwane+2018 (α = 0.8824, μ = 0.9804) está fuera de la cota de DESI DR2:
ΔlnL = −11.6 frente a α = 0, es decir √Δχ² = 4.8σ con μ fijo (4.6σ con el tanh incluido; 4.1σ
dejando μ libre, que se va al borde μ = 0.1). Está 48 veces por encima del límite α < 1.85 × 10⁻² con
μ fijo y 37 veces por encima de α < 2.41 × 10⁻² con μ libre (`resultados/profile_desi_everpresent.json`).

(iii) El modelo del repo es una constante cosmológica Ω_Λ0 ≈ 0.685 más un ruido estacionario
pequeño X de media cero, mientras que Everpresent Λ no tiene constante: Λ es solo ruido de media
cero cuya amplitud tiene que ser del orden de Ω_Λ0 para explicar la aceleración actual.
