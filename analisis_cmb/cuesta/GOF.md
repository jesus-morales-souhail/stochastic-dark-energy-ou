# Cota por bondad de ajuste de ΛCDM (exploración)

**Predicción, escrita antes de correr `gof.py`:**
- θ = 1 entre 1.5 y 2.0 × 10⁻²; θ = 10 entre 1.2 y 1.6 × 10⁻². Ambas cerca de la de fondo fijo y 2–3 veces por
  debajo de la marginalizada.
- θ = 0.1 cerca de la marginalizada.

Exploración: no cambia `results/`, `figures/` ni ningún paper. Script `gof.py`; salidas `gof.{json,txt}` y
`gof_chi2.npz` (χ²_min de cada universo de cada celda).

## Montaje

- **El mismo de `cuesta.py`**, reutilizado por importación: 13 números de DESI DR2, universo base ΛCDM con
  el mejor ajuste de DESI (Ω_m = 0.29740, h·r_d = 101.538 Mpc), residuo OU S_i X(z_i) con S por tipo, ruido
  con la covarianza oficial 13 × 13, **los mismos 2000 vectores de ruido en todas las celdas**, realización
  OU con semilla propia por celda.
- **Malla:** θ ∈ {0.1, 0.3, 1, 3, 10} × σ_X ∈ geomspace(3 × 10⁻³, 0.15, 16); 2000 universos por celda.
- **Ajuste:** solo ΛCDM plano (Ω_m, h·r_d) por mínimo χ² (`cuesta.fit_lcdm`); se guarda χ²_min.
- **Cota por bondad de ajuste:** el σ_X al que P(χ²_min ≤ 10.27 | σ_X, θ) = 0.05, interpolando linealmente
  en ln σ_X. 10.27 es el χ² real de ΛCDM con los 13 números (11 grados de libertad, `CUESTA.md`). En el
  control (σ_X = 0), P = 0.493 y ⟨χ²⟩ = 10.97.
- **Las otras tres cotas** (fondo fijo, perfilado, marginalizado) se recalculan en los cinco θ con la
  maquinaria de `perfil_13.py` sobre los mismos 13 números; en θ = 0.1, 1 y 10 son las de la fila
  "13 números" de `TRECE.md`.

## Resultados

**Cota de σ_X al 95 % por θ (13 números):**

| θ | **bondad de ajuste** | fondo fijo (a) | perfilado (b) | marginalizado (c) | gof / marg | gof / fijo |
|---|---|---|---|---|---|---|
| 0.1 | **0.127** | 0.0183 | 0.0185 | 0.112 | ×1.14 | ×6.9 |
| 0.3 | **0.0821** | 0.0155 | 0.0158 | 0.0656 | ×1.25 | ×5.3 |
| 1 | **0.0439** | 0.0128 | 0.0133 | 0.0377 | ×1.16 | ×3.4 |
| 3 | **0.0281** | 0.0116 | 0.0120 | 0.0249 | ×1.13 | ×2.4 |
| 10 | **0.0215** | 0.0116 | 0.0119 | 0.0193 | ×1.11 | ×1.9 |

**Qué ve el χ² de ΛCDM en cada cota** (P = P(χ²_min ≤ 10.27), ⟨χ²⟩ = media de χ²_min, 11 gdl):

| θ | en la cota fija: P, ⟨χ²⟩ | en la marginalizada: P, ⟨χ²⟩ | en la de bondad de ajuste: P, ⟨χ²⟩ |
|---|---|---|---|
| 0.1 | 0.45, 11.5 | 0.070, 33.8 | 0.05, 43 |
| 0.3 | 0.40, 12.2 | 0.082, 31.8 | 0.05, 44 |
| 1 | 0.35, 13.3 | 0.075, 31.6 | 0.05, 39 |
| 3 | 0.26, 15.3 | 0.070, 30.4 | 0.05, 36 |
| 10 | 0.19, 17.7 | 0.064, 29.9 | 0.05, 35 |

La tabla completa de P y ⟨χ²⟩ en las 80 celdas está en `gof.txt`.

## Por qué

1. **La cota por bondad de ajuste va con la marginalizada, no con la fija: un 11–25 % más floja en todos
   los θ.** Las dos miden lo mismo: cuánto residuo cabe en la parte de los datos que (Ω_m, h·r_d) no puede
   absorber, dado que el χ² de ΛCDM es el esperado (10.3 para 11 gdl). En la cota marginalizada, P(χ²_min ≤
   10.27) ya es 0.06–0.08 y ⟨χ²⟩ ≈ 30–34; la de bondad de ajuste solo pide un poco más (P = 0.05, ⟨χ²⟩ ≈
   35–44). Es algo más floja porque usa un solo número, el χ² total, y tira la forma: la verosimilitud
   marginalizada sabe en qué dirección (el núcleo S_iS_j e^{−θ|Δx|}) aparecería el exceso, y el χ² total
   reparte el mismo exceso entre 11 grados de libertad sin mirar dónde está.
2. **La cota fija es otra cosa.** En ella, el χ² de ΛCDM casi no se mueve (P = 0.45 → 0.19, ⟨χ²⟩ = 11.5 →
   17.7): un residuo de esa amplitud es invisible para la bondad de ajuste. Lo que la fija no es el exceso
   de varianza sino la consistencia con el fondo de referencia a lo largo del modo S, el factor de Occam de
   `FONDO.md`. Por eso la bondad de ajuste queda ×7 por encima de ella en θ = 0.1 y solo ×1.9 en θ = 10: a
   θ pequeño casi todo el residuo es el modo de Ω_m, que la cota fija cuenta y el χ² de ΛCDM absorbe en Ω_m;
   a θ grande el residuo es temblor entre bins, que sí infla el χ², y las tres cotas se acercan.
3. **La cota por bondad de ajuste de un temblor cae como θ^{−1/2}, aproximadamente** (0.127 → 0.0215 entre
   θ = 0.1 y 10, factor 5.9 por factor 100): el exceso de χ² que produce un OU de amplitud σ_X crece con el
   número de bins que descorrelaciona.

## Predicción

- **θ = 1 entre 1.5 y 2.0 × 10⁻², cerca de la fija: no se cumplió.** Da 4.4 × 10⁻², ×3.4 por encima de la
  fija y ×1.16 por encima de la marginalizada.
- **θ = 10 entre 1.2 y 1.6 × 10⁻², cerca de la fija: no se cumplió.** Da 2.15 × 10⁻², ×1.9 por encima de
  la fija y ×1.11 por encima de la marginalizada.
- **θ = 0.1 cerca de la marginalizada: cumplida.** 0.127 frente a 0.112 (×1.14).
- Me equivoqué en qué cota se parece a la bondad de ajuste: creí que la fija medía exceso de varianza y
  resulta que mide consistencia con el fondo; la que mide exceso de varianza es la marginalizada, y la
  bondad de ajuste es su versión sin forma.

**¿Es demasiado floja la cota marginalizada de §3.2?** No: es un 11–25 % más estricta que lo que el χ² de
ΛCDM de DESI permite por sí solo, porque conoce la forma del residuo, y el único dato que DESI tiene contra
un residuo que el fondo no absorbe es precisamente ese χ²; la que no mide estocasticidad es la fija, que
está ×2–7 por debajo porque cuenta la coincidencia con el mejor ajuste de DESI.
