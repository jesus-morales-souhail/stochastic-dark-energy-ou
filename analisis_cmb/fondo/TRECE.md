# El perfil de σ_X con los 13 números de DESI (exploración)

**Predicción, escrita antes de correr `perfil_13.py`:**
- Las cotas con θ ≥ 1 mejoran entre ×1.2 y ×1.4: con el fondo fijo, θ = 1 entre 1.0 y
  1.25 × 10⁻²; marginalizado, θ = 1 entre 3.3 y 3.8 × 10⁻².
- La de θ → 0 sigue siendo 6.75 σ_v, con un σ_v menor (cota entre 1.9 y 2.3 × 10⁻²), y sigue sin
  cota al marginalizar.

Exploración: no cambia `results/`, `figures/` ni ningún paper. Script `perfil_13.py`; salidas
`perfil_13.{json,txt}`.

## Montaje

**Datos:** el vector de 13 puntos de `desi_dr2_data.load_gaussian_bao_vector` (D_V/r_d de BGS;
D_M/r_d y D_H/r_d de los otros seis bins), con su covarianza 13 × 13.
- Se pasa a unidades de α dividiendo por la predicción del mejor ajuste ΛCDM de DESI
  (Ω_m = 0.29740, h·r_d = 101.538 Mpc, el mismo de `FONDO.md`), y la covarianza se reescala igual.

**Modelo:**
- α_th = [(h·r_d)_f/(h·r_d)] · D̂_tipo(z; Ω_m)/D̂_tipo(z; Ω_m,f), con D̂_M = ∫dz/E, D̂_H = 1/E y
  D̂_V = (z D̂_M² D̂_H)^{1/3}.
- **S por tipo:** S = −∂ln D_tipo/∂Ω_Λ, plano, H0 fijo, Ω_Λ = 0.6847; es el convenio de signo de la
  tabla S_z.
- **Núcleo:** K_ij = S_iS_j σ_X² e^{−θ|Δx|}, con x = ln(1+z) de cada fila. Las filas D_M y D_H del
  mismo bin tienen Δx = 0.
- Los mismos tres tratamientos del fondo, la misma malla, la misma búsqueda de θ y el mismo
  criterio que `perfil_fondo.py`.

**Comprobaciones:**
- El mismo código, alimentado con los 7 α_DV, la tabla S_z y su covarianza proyectada, reproduce
  las 15 cotas de `FONDO.md` (`perfil_fondo.json`) a ≤ 2 × 10⁻¹⁴ relativo, y los dos casos sin cota.
- S_DV calculada reproduce la tabla S_z a 8.6 × 10⁻⁴, que es su redondeo.
- Con 13 números y σ_X = 0, el ajuste da Ω_m = 0.2975 y σ(Ω_m) = 0.0086 con h·r_d libre:
  exactamente lo que publica DESI DR2 con solo DESI (0.2975 ± 0.0086, arXiv:2503.14738).

## Resultados

**Cota de σ_X al 95 %, 7 α_DV frente a 13 números (cociente 7/13):**

| tratamiento | θ libre | θ = 0.001 | θ = 0.1 | θ = 1 | θ = 10 |
|---|---|---|---|---|---|
| fijo, 7 | 2.633e-2 | 2.633e-2 | 2.101e-2 | 1.444e-2 | 1.290e-2 |
| **fijo, 13** | **2.305e-2** | 2.305e-2 | 1.833e-2 | **1.280e-2** | 1.160e-2 |
| cociente | ×1.14 | ×1.14 | ×1.15 | ×1.13 | ×1.11 |
| perfilado, 7 | 2.633e-2 | 2.633e-2 | 2.129e-2 | 1.500e-2 | 1.327e-2 |
| **perfilado, 13** | **2.305e-2** | 2.305e-2 | 1.849e-2 | **1.325e-2** | 1.195e-2 |
| cociente | ×1.14 | ×1.14 | ×1.15 | ×1.13 | ×1.11 |
| marginalizado, 7 | sin cota | sin cota | 0.140 | 4.611e-2 | 2.204e-2 |
| **marginalizado, 13** | **sin cota** | sin cota | 0.112 | **3.774e-2** | 1.925e-2 |
| cociente | — | — | ×1.25 | ×1.22 | ×1.14 |

En todos los casos, θ en la cota con θ libre es 10⁻³ (el borde).

**El modo θ → 0 y el fondo:**

| | 7 α_DV | 13 números |
|---|---|---|
| σ_v = (SᵀC⁻¹S)^{−1/2} | 0.00392 | **0.00343** |
| X̂/σ_v (residuo de los datos a lo largo del modo) | +0.020 | +0.031 |
| 6.75 σ_v (supremo con θ → 0) | 0.02644 | **0.02314** |
| σ(Ω_m) con h·r_d libre | 0.0124 | **0.0086** (= DESI DR2) |
| σ(Ω_m) con h·r_d fijo | 0.0038 | 0.0033 |
| Ω_m de mejor ajuste con σ_X = 0 | 0.2978 | 0.2975 |

## Por qué

- **El significado de la cota no cambia.** El límite con θ → 0 sigue siendo el factor de Occam del
  modo de Ω_m: 2.305 × 10⁻² en el borde θ = 10⁻³ y 2.314 × 10⁻² = 6.75 σ_v con θ → 0. El residuo de
  los datos a lo largo del modo sigue siendo nulo (0.03σ), y marginalizando el fondo sigue sin
  haber cota.
- **La mejora es modesta: ×1.11–1.15 con el fondo fijo o perfilado.**
  - σ_v solo baja un 12 %.
  - Separar D_M y D_H añade sobre todo información sobre la forma de H(z), y eso ayuda mucho a
    σ(Ω_m) con h·r_d libre (0.0124 → 0.0086).
  - Pero con el fondo fijo, la dirección de S ya estaba bien medida por D_V.
- **Al marginalizar, la mejora es algo mayor** (×1.14–1.25): ahí sí cuenta que el fondo esté mejor
  determinado.

## Predicción

- **Cotas con θ ≥ 1 mejoran entre ×1.2 y ×1.4: no se cumplió en general.**
  - Con el fondo fijo o perfilado la mejora es ×1.11–1.13, y la cota fija con θ = 1 queda en
    1.28 × 10⁻², algo por encima del rango 1.0–1.25.
  - **Marginalizado con θ = 1 sí se cumplió:** ×1.22, con 3.77 × 10⁻² dentro de 3.3–3.8.
- **θ → 0 sigue siendo 6.75 σ_v con σ_v menor: cumplida.**
  - Cota entre 1.9 y 2.3 × 10⁻²: en el borde. Da 2.305 × 10⁻² en θ = 10⁻³ (2.314 con θ → 0), justo
    por encima de 2.3.
  - Sigue sin cota al marginalizar: cumplida.
