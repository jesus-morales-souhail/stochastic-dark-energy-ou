# Cota σ_X con supernovas (Pantheon+), exploración

**Predicción, escrita antes de correr `perfil_sn.py`:**
- SN no acota el modo θ → 0, porque es degenerado con M_B.
- SN acota σ_X entre 1 y 3 × 10⁻² para θ ≳ 1.

Exploración: no cambia `results/` ni ningún paper. Script `perfil_sn.py`; salidas
`perfil_sn.{json,txt}`.

## Datos

Ficheros de PantheonPlusSH0ES/DataRelease, `Pantheon+_Data/4_DISTANCES_AND_COVAR/`, en `datos/`. El
último commit del repositorio es `c447f0f` (2022-12-21), y los ficheros se descargaron el
2026-10-31.

| fichero | tamaño | sha256 |
|---|---|---|
| `Pantheon+SH0ES.dat` | 0.58 MB | `1cb0fc379ef066afdc2ffd1857681cc478024570d8a3eba284fb645775198cf8` |
| `Pantheon+SH0ES_STAT+SYS.cov` | 33.3 MB | `abf806d966485e64afdb359c87bffc0ecc00d05eff0a31ced66f247385df0fdc` |

Los dos pesan menos de 50 MB, así que están en git y no en `.gitignore`.

**Uso de los datos:**
- `m_b_corr` con la covarianza STAT+SYS.
- Corte z_HD > 0.01 e `IS_CALIBRATOR` = 0: quedan 1580 de 1701, con z entre 0.0102 y 2.26.
- Los ficheros descargados se tratan como entrada no fiable: solo se leen como números y el script
  se corre con `python3 -I`.

**Comprobación:** ΛCDM plano con solo estas SN da Ω_m = **0.3326**. Pantheon+ publica
0.334 ± 0.018 con SN Ia solas (Brout et al. 2022, arXiv:2202.04077, resumen).

## Modelo (el mismo que en DESI)

- **Residuo:** X(x), con x = ln(1+z), es δΩ_Λ en z, y δln D_L(z_i) = S_i X(z_i).
- **S_i:** S_i = −∂ln D_L/∂Ω_Λ, plano, con H0 fijo y Ω_Λ = 0.6847. Es la definición que reproduce
  a 0.1 % la tabla S_z de `scripts/desi_dr2_data.py`, que en el código solo se describe como
  "∂ln D_V/∂Ω_Λ (flat)". El signo no entra en S_iS_j.
- **Núcleo:** en magnitudes es (5/ln 10)² S_iS_j σ_X² e^{−θ|Δx|}, y se suma a la covarianza.
- **Fondo:** m = 5 log₁₀[(1 + z_HEL) ∫₀^{z_HD} dz/E] + M, ΛCDM plano.
  - M se marginaliza analíticamente con prior plano.
  - Ω_m se perfila en [0.1, 0.6] para cada (σ_X, θ).
- **Perfil, como en `scripts/profile_sigma_x_desi.py`:**
  - malla de σ_X: `geomspace(1e-6, 5e-2, 48)`;
  - θ ∈ (10⁻³, 20), con tres búsquedas acotadas más una malla de 40 puntos;
  - ΔlnL = −1.92 desde el máximo, con raíz por `brentq`.
  - Ninguna cota cae dentro de esa malla (≤ 5 × 10⁻²), así que se extiende hasta σ_X = 1. Es una
    extensión, no la malla de DESI.

## Resultados

| caso | σ_X al 95 % | θ en la cota | máximo del perfil |
|---|---|---|---|
| **θ libre en (10⁻³, 20)** | **0.117** (fuera de la malla de DESI) | 10⁻³ (borde inferior); Ω_m = 0.341 | σ_X ≈ 0.010, ΔlnL = +0.03 frente a ΛCDM |
| θ = 0.001 | 0.120 | — | σ_X → 0 |
| θ = 0.1 | 0.119 | — | σ_X → 0 |
| θ = 1 | 0.096 | — | σ_X → 0 |
| θ = 10 | 0.063 | — | σ_X ≈ 0.010, ΔlnL = +0.02 |
| *DESI DR2 (C2)* | *0.02633* | *10⁻³ (borde)* | *σ_X → 0* |

**Diagnósticos del modo θ → 0:**

| diagnóstico | resultado |
|---|---|
| coseno entre S y ∂ln D_L/∂Ω_m (en el mejor ajuste) | **0.99997** |
| θ = 0.001 con Ω_m **fijo** en 0.3326 | σ_X < 0.120 (igual que perfilado) |
| θ = 0.001 con Ω_m **marginalizado** (lineal, prior plano) | **sin cota**: ΔlnL = −0.07 en σ_X = 1 |
| θ = 1 con Ω_m marginalizado | σ_X < 0.177 |
| θ = 10 con Ω_m marginalizado | σ_X < 0.073 |

## Por qué

**El modo θ → 0 es la dirección de Ω_m, no la de M_B.**
- Con H0 fijo y universo plano, −∂/∂Ω_Λ = ∂/∂Ω_m, así que el residuo totalmente correlacionado
  S_iX es, a primer orden, un cambio de Ω_m: coseno 0.99997.
- M_B absorbe solo la parte constante, y S_i no es constante (va de −0.008 a −0.87).
- Como Ω_m está libre, los datos no pueden medir ese modo.

**La "cota" de 0.12 a θ → 0 sale de perfilar, no de información.**
- Con Ω_m en su mejor ajuste, el residuo a lo largo del modo es ≈ 0. Añadir varianza σ² en esa
  dirección solo cambia ln det C, y ΔlnL = −½ ln(1 + σ²/σ_v²), donde σ_v ≈ 0.018 es el error de
  Ω_m de las SN.
- −1.92 se alcanza en σ ≈ √(e^{3.84} − 1) σ_v ≈ 6.7 σ_v ≈ 0.12, que es justo lo que sale. Fijar
  Ω_m da lo mismo por la misma razón.
- Marginalizando Ω_m, los términos de ln det se compensan y la cota desaparece. Es un factor de
  Occam del perfil, no una medida.

**Con θ ≳ 1, el modo se descorrelaciona entre redshifts** y deja de parecerse a un cambio de Ω_m.
Ahí las SN sí lo acotan: 0.096 (θ = 1) y 0.063 (θ = 10) perfilando Ω_m; 0.18 y 0.07
marginalizándolo. Aun así, quedan 2.4–4.5 veces por encima de la cota de DESI.

**Comparación con DESI.** En DESI la cota (2.633 × 10⁻²) también está en el borde θ = 10⁻³, pero
allí el fondo está fijo (α = 1: ni Ω_m ni h·r_d libres), así que el modo totalmente correlacionado
se mide directamente frente al fiducial. Las SN, con Ω_m libre, no pueden hacerlo.

**Pregunta abierta, no comprobada aquí:** si en el perfil de DESI se dejaran libres Ω_m y h·r_d,
¿se volvería también degenerado el modo θ → 0 y pasaría la cota a depender del tratamiento de esos
parámetros (perfilados o marginalizados)? No se ha corrido, y no cambia nada de C2 por ahora.

## Predicción

- **"SN no acota el modo θ → 0 (degenerado con M_B)": cumplida en el fondo, no en el mecanismo.**
  - El modo no tiene información: marginalizando Ω_m no hay cota.
  - Perfilando sale un número (0.12), pero es el factor de Occam del perfil.
  - El modo es degenerado con **Ω_m** (coseno 0.99997), no con M_B, que solo absorbe la parte
    constante.
- **"SN acota σ_X entre 1 y 3 × 10⁻² para θ ≳ 1": no se cumplió.** Las cotas son 0.096 (θ = 1) y
  0.063 (θ = 10) perfilando Ω_m, y 0.18 y 0.07 marginalizándolo: entre 2 y 6 veces por encima del
  rango predicho.
