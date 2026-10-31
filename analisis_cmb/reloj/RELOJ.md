# Cota σ_X frente al reloj de estacionariedad (exploración)

**Predicción, escrita antes de correr `perfil_reloj.py`:** las cuatro cotas son iguales a ±10 % de
2.633 × 10⁻².

Exploración: no cambia ningún resultado ni ningún paper.

## Pregunta

El núcleo OU de `scripts/desi_dr2_data.py::add_ou_kernel` es
C_ij += S_i S_j σ_X² exp(−θ |x_i − x_j|), con x = ln(1+z). Es decir, el proceso es estacionario en
ln a. ¿Cambia la cota de σ_X si la estacionariedad se define con otro reloj?

## Montaje

`perfil_reloj.py` importa `scripts/profile_sigma_x_desi.py` y reutiliza sus funciones sin llamar a su
`main()`, así que no escribe en `results/` ni en `figures/`. Mantiene:
- los mismos datos: DESI DR2, 7 α_DV, con la covarianza oficial 13 × 13 proyectada;
- la misma malla de σ_X: 48 puntos en [10⁻⁶, 5 × 10⁻²];
- la misma búsqueda de θ: (10⁻³, 20), igual para todos los relojes;
- el mismo criterio: ΔlnL = −1.92 desde el máximo del perfil, con raíz por `brentq`.

Solo cambia la coordenada x. Los relojes (b) y (c) se calculan con ΛCDM-Planck: H0 = 67.36,
Ω_m = 0.3153, Ω_r = ω_r/h² con el ω_r de `everpresent/resultados/rd_planck.json`, plano.

| reloj | x |
|---|---|
| (a) | ln(1+z) = −ln a (el actual) |
| (b) | tiempo cósmico t, en unidades de 1/H0 |
| (c) | tiempo conforme η, en unidades de 1/H0 |
| (d) | z |

Salidas: `perfil_reloj.{json,txt}`.

## Resultados

| reloj | σ_X al 95 % | frente a 2.633e-2 | θ en la cota (unidades del reloj) | θ · extensión | máx. ΔlnL frente a ΛCDM |
|---|---|---|---|---|---|
| (a) ln a | 2.6328 × 10⁻² | **+0.000 %** (idéntico bit a bit a `results/profile_sigma_x/profile_sigma_x.json`) | 10⁻³ (borde) | 9.4 × 10⁻⁴ | −2.0 × 10⁻⁸ |
| (b) t | 2.6391 × 10⁻² | +0.24 % | 10⁻³ (borde) | 5.2 × 10⁻⁴ | −2.0 × 10⁻⁸ |
| (c) η | 2.6323 × 10⁻² | −0.02 % | 10⁻³ (borde) | 1.0 × 10⁻³ | −2.0 × 10⁻⁸ |
| (d) z | 2.6163 × 10⁻² | −0.63 % | 10⁻³ (borde) | 2.0 × 10⁻³ | −1.9 × 10⁻⁸ |
| θ → 0 (θ = 10⁻¹², cualquier reloj) | 2.6441 × 10⁻² | +0.43 % | — | 0 | — |

Separaciones entre bins consecutivos, en orden de z (z_eff = 0.295, 0.510, 0.706, 0.934, 1.321,
1.484, 2.330):

| reloj | dx₁₂ | dx₂₃ | dx₃₄ | dx₄₅ | dx₅₆ | dx₆₇ | extensión total |
|---|---|---|---|---|---|---|---|
| (a) ln a | 0.1536 | 0.1220 | 0.1254 | 0.1824 | 0.0679 | 0.2931 | 0.9445 |
| (b) t [1/H0] | 0.1234 | 0.0865 | 0.0782 | 0.0951 | 0.0302 | 0.1030 | 0.5165 |
| (c) η [1/H0] | 0.1725 | 0.1388 | 0.1419 | 0.2011 | 0.0725 | 0.2943 | 1.0212 |
| (d) z | 0.2150 | 0.1960 | 0.2280 | 0.3870 | 0.1630 | 0.8460 | 2.0350 |

Las separaciones sí cambian mucho de un reloj a otro. La extensión total varía un factor 4, de 0.52
a 2.04, y la forma también: en t, el último intervalo (1.48 → 2.33) pesa poco; en z es el mayor con
diferencia.

## Por qué las diferencias son las que son

- **El máximo del perfil está en σ_X → 0 en los cuatro casos** (ΔlnL máximo ≈ −2 × 10⁻⁸, en el
  primer nodo, σ_X = 10⁻⁶). Los datos no piden residuo, y la referencia del criterio es la misma
  para todos.
- **En la cota, el θ óptimo es el borde inferior de la búsqueda (10⁻³) en los cuatro relojes.** Cerca
  de la cota, el perfil prefiere el proceso totalmente correlacionado. Cuando θ·dx → 0, el núcleo
  tiende a σ_X² S_i S_j, de rango 1, que no depende del reloj. Ese límite da 2.6441 × 10⁻².
- **Las diferencias residuales, de −0.6 % a +0.2 %, vienen de que el borde θ_mín = 10⁻³ es fijo en
  unidades de cada reloj.** θ_mín · extensión es distinto en cada uno (de 5.2 × 10⁻⁴ en t a
  2.0 × 10⁻³ en z), y eso deja a cada reloj a una distancia distinta del límite de rango 1. El
  orden lo confirma: cuanto menor es θ·extensión, más cerca de 2.6441 × 10⁻² queda la cota.
  t (0.52) da 2.6391; ln a (0.94) y η (1.02) dan 2.6328 y 2.6323; z (2.04) da 2.6163.

Observación lateral, sin cambiar nada: la cota publicada, 2.633 × 10⁻², está definida con θ ≥ 10⁻³.
El supremo con θ → 0 es 2.644 × 10⁻², un 0.43 % mayor. El redondeo a 2.6 × 10⁻² no cambia.

## Conclusión

Las cuatro cotas están dentro del ±0.7 % de 2.633 × 10⁻² (y del ±1.1 % del límite θ → 0): **la cota no
depende del reloj; la fija el límite totalmente correlacionado (θ → 0).**

**Predicción: cumplida.** Las cuatro están a ±10 %; de hecho, a ±0.7 %.
