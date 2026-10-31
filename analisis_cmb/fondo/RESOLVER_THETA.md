# σ_X mínimo para resolver θ (exploración)

**Predicción, escrita antes de correr `resolver_theta.py`:**
- DESI a σ_X = 1.5e-2 da un ΔlnL esperado < 1 entre θ = 0.001 y θ = 1, y necesita σ_X > 3e-2 para
  llegar a 1.92.
- Un cartografiado tipo Euclid necesita σ_X > 1e-2.

Exploración: no cambia ningún resultado ni ningún paper. Script `resolver_theta.py`; salidas
`resolver_theta.{json,txt}`.

## Montaje

**Separación esperada entre dos memorias del residuo.** Se comparan dos hipótesis,
H₁ = N(0, C₀ + K(σ_X, θ₁)) y H₂ = N(0, C₀ + K(σ_X, θ₂)), con K_ij = S_iS_j σ_X² e^{−θ|Δx|} y
x = ln(1+z).
- La divergencia de Kullback–Leibler KL(1‖2) es el ΔlnL esperado si la verdad es H₁.
- Se da en las dos direcciones, su **media** (la "simétrica" de la tabla principal: el ΔlnL
  esperado promediando sobre cuál es la verdadera, que es lo que se compara con 1.92) y su suma
  (Jeffreys).
- Pares (θ₁, θ₂): (0.001, 1), (0.001, 10) y (1, 10).
- σ_X en `geomspace(1e-3, 1e-1, 241)`, y el cruce con 1.92 se calcula por `brentq`.

**Cartografiados:**
- **DESI DR2:** los 7 α_DV, la tabla S_z y la covarianza oficial proyectada C₀, igual que en
  `scripts/profile_sigma_x_desi.py`.
- **Tipo Euclid:** 13 tramos equiespaciados en ln(1+z) entre z = 0.1 y 2.0, errores del 0.5 % en α
  sin correlación, y S(z) = −∂ln D_V/∂Ω_Λ (plano, H0 fijo, Ω_Λ = 0.6847). Esta definición
  reproduce la tabla de DESI a 8.6 × 10⁻⁴, que es su redondeo.

**Fondo:**
- **Fijo**, como en el perfil de C2.
- **Marginalizado:** verosimilitud restringida a lo ortogonal a las plantillas T = [S, 1], es decir,
  ∂α/∂Ω_m y la escala h·r_d. Es exacto para la escala y de primer orden para Ω_m (coseno 0.99996,
  `FONDO.md`).

## Resultados

**σ_X al que la divergencia media llega a 1.92 (ΔlnL esperado, 95 %):**

| par θ₁ – θ₂ | DESI, fondo fijo | DESI, fondo marginalizado | tipo Euclid, fijo | tipo Euclid, marginalizado |
|---|---|---|---|---|
| 0.001 – 1 | **2.61 × 10⁻²** | **3.23 × 10⁻²** | **1.49 × 10⁻²** | **1.73 × 10⁻²** |
| 0.001 – 10 | 1.40 × 10⁻² | 1.65 × 10⁻² | 7.5 × 10⁻³ | 8.6 × 10⁻³ |
| 1 – 10 | 2.86 × 10⁻² | 3.41 × 10⁻² | 1.14 × 10⁻² | 1.29 × 10⁻² |

**Divergencia media a σ_X = 1.5 × 10⁻² y 2.633 × 10⁻²:**

| par | DESI fijo | DESI marginalizado | tipo Euclid fijo | tipo Euclid marginalizado |
|---|---|---|---|---|
| 0.001 – 1 | **0.48** / 1.96 | 0.25 / 1.16 | 1.95 / 7.2 | 1.33 / 5.3 |
| 0.001 – 10 | 2.24 / 7.5 | 1.51 / 5.8 | 9.0 / 29.6 | 7.6 / 26.6 |
| 1 – 10 | 0.79 / 1.74 | 0.54 / 1.38 | 2.95 / 6.0 | 2.50 / 5.4 |

**Por direcciones, comparado con la cota de `FONDO.md` en el θ verdadero.** Esto es lo que decide si
hay ventana. Cuando la verdad es θ_v, la separación esperada es KL(θ_v‖otra), y σ_X tiene que estar
por debajo de la cota actual en θ_v. Para DESI:

| verdad | alternativa | fondo | σ_X para KL = 1.92 | cota actual en la verdad (FONDO.md) | ¿hay ventana? |
|---|---|---|---|---|---|
| θ = 1 | 0.001 | fijo | 2.24 × 10⁻² | 1.44 × 10⁻² | **no** |
| θ = 1 | 0.001 | marginalizado | 2.78 × 10⁻² | 4.61 × 10⁻² | **sí**: 2.8–4.6 × 10⁻² |
| θ = 0.001 | 1 | fijo | 4.12 × 10⁻² | 2.63 × 10⁻² | **no** |
| θ = 0.001 | 1 | marginalizado | 5.13 × 10⁻² | sin cota | sí, por encima de 5.1 × 10⁻² (sin límite superior) |
| θ = 10 | 0.001 | fijo | 1.29 × 10⁻² | 1.29 × 10⁻² | en el borde |
| θ = 10 | 0.001 | marginalizado | 1.45 × 10⁻² | 2.20 × 10⁻² | sí: 1.5–2.2 × 10⁻² |
| θ = 10 | 1 | fijo | 2.43 × 10⁻² | 1.29 × 10⁻² | no |
| θ = 10 | 1 | marginalizado | 2.66 × 10⁻² | 2.20 × 10⁻² | no |

## Por qué

- **Distinguir θ cuesta más que acotar σ_X.** La cota de C2 pregunta si hay varianza extra; resolver
  θ pregunta cómo se reparte esa varianza entre 7 bins correlacionados. Con 7 puntos y errores del
  0.4–1 %, la forma de la correlación solo se nota cuando σ_X es comparable o mayor que los errores.
- **Con el fondo fijo no hay ventana para θ = 1:** distinguirlo de θ = 0.001 exige 2.2 × 10⁻², y los
  datos ya excluyen θ = 1 por encima de 1.44 × 10⁻².
- **Marginalizar el fondo** debilita la separación (necesita σ_X un 20–25 % mayor), pero afloja mucho
  más la cota con θ ≥ 1: de 1.44 a 4.61 × 10⁻². Por eso se abre una ventana estrecha.
- **El extremo θ = 10 frente a 0.001 es más fácil.** Un residuo casi blanco entre bins es muy
  distinto de uno rígido, pero también queda justo en el borde de lo permitido.
- **Un cartografiado tipo Euclid** (13 bins, 0.5 %) baja el umbral para 0.001 frente a 1 a
  1.5–1.7 × 10⁻²: unas 1.7–1.9 veces menos que DESI. Ese umbral está a la altura de las cotas
  actuales con θ = 1, no muy por debajo.

## Conclusión

**Con DESI DR2, θ no es resoluble dentro de las cotas si el fondo está fijo, que es como se calcula
C2.** Distinguir θ = 1 de θ = 0.001 exige σ_X ≈ 2.2–2.6 × 10⁻², por encima de la cota de 1.44 × 10⁻²
para θ = 1. Solo con el fondo marginalizado se abre una ventana estrecha (2.8–4.6 × 10⁻² si la verdad
es θ = 1). Un cartografiado tipo Euclid lo resolvería solo si σ_X está cerca de las cotas actuales
(≥ 1.5 × 10⁻²).

## Predicción

- **DESI a σ_X = 1.5e-2 da ΔlnL esperado < 1 entre θ = 0.001 y 1: cumplida.** Da 0.48 con el fondo
  fijo y 0.25 marginalizado.
- **DESI necesita σ_X > 3e-2 para 1.92: a medias.** Se cumple con el fondo marginalizado
  (3.23 × 10⁻²). Con el fondo fijo basta con 2.61 × 10⁻², algo menos de lo predicho.
- **Tipo Euclid necesita > 1e-2 (para θ = 0.001 frente a 1): cumplida.** Necesita 1.49 × 10⁻²
  (fijo) y 1.73 × 10⁻² (marginalizado). Para el par extremo 0.001 frente a 10 basta con
  7.5–8.6 × 10⁻³, por debajo de 1e-2.
