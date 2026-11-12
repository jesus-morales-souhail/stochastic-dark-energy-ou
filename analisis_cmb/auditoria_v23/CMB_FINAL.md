# CMB_FINAL: ¿descarta el CMB completo las historias que permiten los techos de la marginal −1.92?

Exploración (no toca PREPRINT ni CLAIMS). Script `cmb_final.py`; salida `cmb_final.json`.

Pregunta: las historias que permiten BAO + θ* + ω_m + Pantheon+ en los techos de la
marginal −1.92 (OU θ = 1: 0.287; θ = 3: 0.177; θ = 10: 0.116), ¿las descarta el CMB
completo que tenemos (TT ℓ < 30, bandas φφ de Planck y distancia en los picos)?

## Predicción (escrita antes de calcular nada)

1. Historias con ρ_DE ≤ 0 en algún punto: menos del 5 % en θ = 1 y ninguna en θ = 3 y 10.
2. Más del 95 % de las que entran en CAMB son fieles en H² (|Δρ_DE|/ρ_total < 1e-3)
   en todo el rango 1e-6 ≤ a ≤ 1, con el borde arreglado.
3. Δχ² total (TT ℓ < 30 + φφ + picos): mediana < 0.5, p95 < 2 y ninguna > 4 en los tres techos.
4. El máximo de |Δ/σ| cae más veces en φφ que en TT.
5. Lectura: el CMB completo que tenemos no mueve las cotas de la marginal.

## Montaje (del prompt)

- Historias: 40 por techo, sacadas del posterior lineal de marginal_guiada.py (Matheron)
  con σ_X en el techo, es decir historias compatibles con los datos, no del prior.
  Cuenta antes cuántas tienen ρ_DE ≤ 0 (no entran en CAMB; se cuentan aparte).
- CAMB como en cmb_techos.py (DarkEnergyPPF, lente lineal, H0 resuelto para el θ* de
  cada historia), con el borde arreglado: por encima de z* el camino se mantiene
  constante (X(z) = X(z*)), así que ρ_DE no cambia hacia atrás y w = −1 por debajo de
  la tabla. Declara en la nota que este arreglo es nuevo y por qué.
- Referencia de cada historia: ΛCDM con el mismo ω_b y ω_m y H0 resuelto para su θ*.

## Montaje tal como se corrió (`cmb_final.py prep` 10 s, `run` 90 s con 8 procesos; `cmb_final_prep.json`, `cmb_final.json`)

- **Techos:** los `exact_192` de `marginal_guiada.json` (θ = 1: 0.2875; θ = 3: 0.1773; θ = 10: 0.1162), con BAO + θ* + ω_m + Pantheon+
  (1595 filas).
- **Historias del posterior lineal:** 5000 sorteos por techo de (δp, x) con la regla de Matheron de `marginal_guiada.py` (mismo código:
  δp de su marginal gaussiana, x = x₀ + σ² K Aᵀ S⁻¹(d − b₀ − B δp − A x₀ − n₀), σ_X fijado en el techo; semilla [20261119, celda]).
  ρ_DE/ρ_crit,0 = Ω_Λ(h, ω_m) + x(z) − x(0) con (h, ω_m) = p₀ + δp del sorteo. Se cuentan, sobre los 5000, las que tienen H² ≤ 0 y las
  que tienen ρ_DE ≤ 0 en algún nodo de la malla (0 ≤ z ≤ z*; más allá el camino es constante, así que la malla basta). Para
  información, el peso exacto/lineal de MARGINAL_GUIADA (w = e^{−(χ²_ex − χ²_lin)/2}) de cada historia corrida.
- **Suavizado** como `cmb_techos.py` (gaussiana de Δ = 0.05 en ln a para θ = 1 y 3, 0.0125 para θ = 10; prolongación AR(1) hacia u < 0 con
  normales nuevas, que es el condicional del prior para un OU de Markov), sobre las primeras 160 historias con ρ_DE > 0; ninguna pasa a
  ρ_DE ≤ 0 al suavizar; se corren las 40 primeras (por índice del sorteo) de cada techo.
- **Borde arreglado (nuevo en este test, declarado aquí):** por encima de z* el camino se mantiene constante, X(z) = X(z*), así que ρ_DE no
  cambia hacia atrás y w = −1 por debajo de la tabla. En la práctica: (i) antes de suavizar, el camino se prolonga hacia u > u* con el
  valor constante X(u*) durante 0.4 e-folds (más que los 6Δ = 0.3 del núcleo), para que el suavizado no tenga borde en z* y w llegue
  a −1 de forma continua; (ii) la tabla w(a) que recibe CAMB cubre toda la malla suavizada (a de 1 a e^{−(u*+0.4)} = 6.1 × 10⁻⁴) y se
  prolonga hasta a = 10⁻⁶ con 59 nodos en w = −1. Por qué: `cmb_techos.py` cortaba la tabla en a = 10⁻³ y CAMB extrapolaba w(10⁻³)
  constante hacia el pasado; con w(10⁻³) > 1/3 (7 de 38 en θ = 10), ρ_DE ∝ a^{−3(1+w)} crecía hasta dominar antes de la recombinación,
  un artefacto de la tabla y no del camino, que obligó a separar "fieles" de "no fieles". Con el camino constante más allá de z* la
  historia tiene w = −1 ahí por construcción y CAMB no tiene nada que extrapolar. Comprobado: w en el suelo de la malla es −1 con error
  < 10⁻¹³ en las 120 historias.
- **CAMB y referencia** como `cmb_techos.py` (funciones reutilizadas: `make_params`, `spectra`, `joint_per_l`, `_stats`): DarkEnergyPPF con
  la tabla, lente lineal, H0 = 100 h y ω_m del sorteo, ω_b = 0.02237, τ, n_s, A_s de Planck 2018; referencia ΛCDM con los mismos ω_b y ω_m y
  H0 resuelto por bisección para el 100 θ* de la historia (`ou_cmb.solve_H0_for_thetastar`).
- **Fidelidad por dos caminos** en 10⁻⁶ ≤ a ≤ 1 (tabla más prolongación): ρ_DE que integra CAMB frente a la del camino; fiel si
  máx_a |Δρ_DE|/ρ_total < 10⁻³.
- **Estadísticos**, historia menos referencia: Δχ² de TT (2 ≤ ℓ < 30) y de las bandas φφ frente a Planck 2018 (`estadisticos.chi2`),
  distancia conjunta TT/TE/EE en los picos (ℓ ∈ [30, 2500], `picos.distance`); total = suma; dónde cae el máximo de |Δ/σ|.

## Resultados

**Historias con ρ_DE ≤ 0 (sobre 5000 sorteos del posterior en cada techo; H² ≤ 0: 0 en los tres):**

| techo | σ_X | ρ_DE ≤ 0 en algún punto | ρ_DE mín / Ω_Λ (p5, mediana) | ESS del peso exacto/lineal | peso de las 40 corridas (16–50–84) | χ²_ex − χ²_lin de las 40 (16–50–84) |
|---|---|---|---|---|---|---|
| θ = 1 | 0.287 | 2246 (44.9 %) | −0.63, 0.05 | 420 | 0.016 – 0.38 – 2.8 | −2.1 – +1.9 – +9.4 |
| θ = 3 | 0.177 | 414 (8.3 %) | −0.05, 0.30 | 182 | 0.033 – 0.35 – 1.4 | −0.6 – +2.1 – +6.9 |
| θ = 10 | 0.116 | 0 (0 %) | 0.32, 0.51 | 517 | 0.009 – 0.32 – 1.0 | 0.0 – +2.3 – +10.7 |

**CAMB (40 historias por techo, 0 fallos, 120 de 120 fieles):**

| techo | fieles | error en H² mediana / máx | \|w + 1\| máx mediana / p95 / máx | H0 mediana | H0_ref − H0 mediana (mín, máx) |
|---|---|---|---|---|---|
| θ = 1 | 40 / 40 | 1.2 × 10⁻⁵ / 2.7 × 10⁻⁵ | 1.94 / 3.90 / 5.69 | 66.3 | +0.21 (−2.5, +3.3) |
| θ = 3 | 40 / 40 | 1.2 × 10⁻⁵ / 1.8 × 10⁻⁵ | 1.72 / 2.89 / 3.69 | 66.4 | +0.96 (−2.8, +2.9) |
| θ = 10 | 40 / 40 | 2.1 × 10⁻⁴ / 4.0 × 10⁻⁴ | 3.96 / 5.36 / 8.22 | 66.3 | +1.00 (−2.1, +3.2) |

**Δχ² total y por bloque (historia − referencia), 40 por techo:**

| techo | total mediana / p95 / máx / mín | TT ℓ < 30 mediana / p95 / máx | φφ mediana / p95 / máx / mín | picos mediana / p95 / máx | total > 1 | total > 4 |
|---|---|---|---|---|---|---|
| θ = 1 | −0.15 / 0.18 / 0.76 / −0.64 | −0.02 / 0.12 / 0.20 | −0.22 / 0.14 / 0.65 / −0.74 | 0.04 / 0.10 / 0.16 | 0 | 0 |
| θ = 3 | −0.23 / 0.11 / 0.24 / −1.15 | 0.02 / 0.10 / 0.15 | −0.29 / 0.05 / 0.17 / −1.24 | 0.02 / 0.11 / 0.12 | 0 | 0 |
| θ = 10 | −0.23 / 0.03 / 0.13 / −0.61 | 0.01 / 0.11 / 0.19 | −0.21 / −0.03 / 0.08 / −0.58 | 0.01 / 0.03 / 0.04 | 0 | 0 |

**Dónde cae el máximo de |Δ/σ|:** φφ en el 97.5 % / 95 % / 87.5 % (θ = 1 / 3 / 10), TT ℓ < 30 en el 2.5 % / 5 % / 12.5 % (siempre en ℓ = 2),
picos nunca (|Δ/σ| por ℓ ≤ 0.014; el ℓ del máximo está en 1616–1620). Máximos de |Δ/σ|, mediana / p95: TT 0.05 / 0.08, 0.05 / 0.08,
0.03 / 0.05; φφ 0.18 / 0.29, 0.14 / 0.30, 0.09 / 0.17.

## Predicción, punto por punto

1. **No cumplida.** Con ρ_DE ≤ 0 en algún punto: 44.9 % en θ = 1 (predicho < 5 %), 8.3 % en θ = 3 (predicho ninguna) y 0 % en θ = 10
   (cumplido). En el techo de θ = 1 la mediana del mínimo de ρ_DE es el 5 % de Ω_Λ: las historias que el posterior admite ahí rozan o
   cruzan ρ_DE = 0. Esas historias no entran en CAMB (PPF integra ρ a partir de w) y quedan fuera de todo lo que sigue.
2. **Cumplida.** 120 de 120 fieles en H² en 10⁻⁶ ≤ a ≤ 1 (error máximo 4.0 × 10⁻⁴, en θ = 10, donde |w + 1| llega a 8). Con el borde
   arreglado no hay grupo "no fiel": el artefacto de CMB_TECHOS era de la tabla.
3. **Cumplida.** Mediana −0.15 / −0.23 / −0.23 (< 0.5), p95 0.18 / 0.11 / 0.03 (< 2), máximo 0.76 / 0.24 / 0.13 (ninguna > 4).
4. **Cumplida.** El máximo de |Δ/σ| cae en φφ en el 87.5–97.5 % de las historias y en TT en el resto; nunca en los picos.
5. **Lectura: cumplida con una salvedad.** Entre las historias que CAMB puede integrar, el CMB completo que tenemos añade Δχ² ≤ 0.8 a
   historias que ya pasan BAO + θ* + ω_m + SN en el techo, con |Δ/σ| ≤ 0.3 en el peor punto: no mueve las cotas de la marginal. La
   salvedad es el punto 1: en θ = 1 casi la mitad de las historias del techo (y una de doce en θ = 3) tienen ρ_DE ≤ 0 en algún punto y
   este test no dice nada de ellas; si ρ_DE < 0 se considera físicamente inadmisible, la cota de la marginal en θ = 1 cambiaría por un
   prior, no por el CMB.

## Respuesta

No. Las historias que permiten BAO + θ* + ω_m + Pantheon+ en los techos de la marginal −1.92 pasan el CMB completo que tenemos sin
pena: Δχ² total mediana negativa en los tres techos, p95 ≤ 0.18, máximo 0.76, ninguna por encima de 1; lo poco que hay está en las
bandas de lente φφ (|Δ/σ| ≤ 0.3), no en TT ℓ < 30 ni en los picos (|Δ/σ| ≤ 0.014). Con el borde arreglado las 120 historias son fieles
en H² y el artefacto de CMB_TECHOS desaparece. Lo que sí aparece es que el techo de θ = 1 (0.287) está poblado en un 45 % por historias
con ρ_DE ≤ 0 en algún punto (8 % en θ = 3, 0 % en θ = 10), que el CMB con PPF no puede juzgar; para ellas el test no se aplica.

## Pesos exacto/lineal de las 40 historias (añadido después; `cmb_final_pesos.json`, sin recalcular CAMB)

Los 40 pesos por historia no están en `cmb_final.json` (solo sus percentiles en `prep`): se leen de `cmb_final_theta_*.npz`, la salida de
`prep` que guarda w, χ²_lin y χ²_ex de cada historia corrida, en el mismo orden `idx_run` que `runs` del json (comprobado). Los Δχ²
totales por historia sí están en `runs`. Cuantiles ponderados: el valor donde la suma acumulada de pesos cruza 0.5 y 0.95; el máximo
ponderado coincide con el máximo (todos los pesos son > 0).

| techo | w, percentiles 16–84 | w máximo (peso del total) | ESS de las 40 | Δχ² total sin ponderar: mediana / p95 / máx | ponderado: mediana / p95 / máx |
|---|---|---|---|---|---|
| θ = 1 | 0.016 – 2.81 | 8.9 (19 % del total) | 11.4 | -0.15 / +0.18 / +0.76 | -0.29 / +0.10 / +0.76 |
| θ = 3 | 0.033 – 1.35 | 16.4 (35 % del total) | 5.8 | -0.23 / +0.11 / +0.24 | -0.24 / -0.01 / +0.24 |
| θ = 10 | 0.009 – 1.00 | 3.2 (16 % del total) | 16.1 | -0.23 / +0.03 / +0.13 | -0.22 / +0.13 / +0.13 |

El ESS de las 40 es de 6–16 historias: la media ponderada la dominan unas pocas (en θ = 3 una historia lleva el 35 % del peso). Ponderar
mueve las medianas y los p95 en menos de 0.3 y no cambia el máximo ni la lectura. Y con 40 historias sin ningún Δχ² > 4, lo que se
puede decir es que la probabilidad de un caso así queda por debajo de ~7.5 % (regla de 3/N, N = 40), no que sea cero; ponderando por w el
tamaño efectivo es menor y esa cota es más floja.
