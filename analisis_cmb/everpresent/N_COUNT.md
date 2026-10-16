# Conteo N y amplitud ingenua de Everpresent Λ

Script: `n_count.py` (salida en `n_count_salida.txt`). ΛCDM plano Planck 2018: H0 = 67.36,
ω_b = 0.02237, ω_c = 0.1200, radiación = fotones (T0 = 2.7255 K) + 3.046 ν sin masa
(masa de ν ignorada). ℓ_P = √(ħG/c³) = 1.616 × 10⁻³⁵ m (sin reducir).

Fórmulas:

- 4-volumen del pasado causal de (aquí, hoy), [A] ec. (3) con t = t0:
  V_4 = (4π/3) c ∫_0^{t0} dt a(t)³ r_c(t)³,  r_c(t) = c ∫_t^{t0} dt'/a(t')
  (radio comóvil del cono de luz pasado en t; r_c(0) = horizonte de partículas hoy,
  3.185 c/H0 = 14175 Mpc). Resultado: V_4 = 0.12347 (c/H0)⁴ = 4.392 × 10¹⁰³ m⁴.
  N_4 = V_4/ℓ_P⁴.
- N_área = A_H/(4ℓ_P²), A_H = 4π(c/H0)².
- δΛ/Λ ~ 1/√N;  X = Ω_Λ δΛ/Λ con Ω_Λ = 0.686.

> **DESCARTADA (2026-10-16).** La tabla siguiente lee 1/√N como δΛ/Λ_obs. Era un error del
> enunciado: en [A] y [Z], 1/√V es la fluctuación *absoluta* δρ_Λ = α/√V en unidades de
> Planck, no una fracción de Λ_obs. La cantidad correcta y su conversión a X están en la
> "Nota de normalización" de más abajo y en `EVERPRESENT_vs_DESI.md`. Se deja la tabla como
> registro de por qué no se usa.

| conteo | N | δΛ/Λ | X equivalente | X / 2.506e-2 |
|---|---|---|---|---|
| 4-volumen V_4/ℓ_P⁴ | 6.44 × 10²⁴² | 3.94 × 10⁻¹²² | 2.71 × 10⁻¹²² | 1.08 × 10⁻¹²⁰ |
| área A_H/(4ℓ_P²) | 2.27 × 10¹²² | 6.64 × 10⁻⁶² | 4.56 × 10⁻⁶² | 1.82 × 10⁻⁶⁰ |

Nota de normalización. En [A] ec. (1)–(2) y [Z] ec. (2.3), 1/√V es la fluctuación *absoluta*
de ρ_Λ en unidades de Planck, es decir, relativa a la densidad de Planck, no a Λ_obs. Por
eso las filas de arriba (1/√N leído como δΛ/Λ_obs) no son la amplitud de Everpresent Λ.
Con la normalización de los papers (coeficiente 1), δρ_Λ/ρ_Λ0 = (c/H0)²/(3Ω_Λ√V_4), que no
depende de ℓ_P:

| convención | δρ_Λ/ρ_Λ0 | X | X / 2.506e-2 |
|---|---|---|---|
| unidades reducidas, 8πG = ħ = c = 1 (las que declaran [A] y [Z]) | 1.38 | 0.95 | 38 |
| G = ħ = c = 1 (sin reducir; ×8π) | 34.7 | 23.8 | 951 |

Hay que multiplicar por el coeficiente de amplitud (α en [A] ec. 5). Con α = 0.02 sale
X ≈ 0.019 (reducidas) o ≈ 0.48 (sin reducir); [A] Sec. III da para α = 0.02 una amplitud de
Ω_Λ "of order 0.5".

Variante de [Z] ec. (2.6) tal como está impresa (integral interior ∫_0^{t'}):
V = 28.41 (c/H0)⁴. No se usa en las tablas.
