# Cota exacta (no lineal) del temblor en el convenio B (exploración v3)

**Decisión registrada:** convenio principal B (δE² = X(z) − X(0): la parte constante del residuo la absorbe Λ; ω_m es la
materia y no responde a X), confirmado por una segunda opinión independiente; el convenio A (tabla S) se presentará como
análisis de sensibilidad.

**Predicción, escrita antes de correr `exacto.py`:** con BAO + CMB la cota exacta en θ = 1 está entre 0.2 y 0.5 y a menos
de un factor 1.5 de la lineal; en θ ≤ 0.3 no hay cota por debajo de 0.7.

Exploración: no cambia `results/`, `figures/` ni ningún paper. Script `exacto.py`; salidas `exacto.{json,txt}` y
`exacto_chi2.npz` (χ²_min de cada realización, 15 filas y 13 filas, y la marca de descarte). Tiempo: 21 s con 16 procesos.

## Montaje

- **Historia:** ρ_DE(z)/ρ_crit,0 = Ω_Λ + X(z) − X(0), H² = 100²[ω_m(1+z)³ + ω_r(1+z)⁴ + (h² − ω_m − ω_r) + h²(X(z) − X(0))]:
  la de `escalon.py` (B) con el escalón sustituido por el camino OU. Nada linealizado.
- **Realizaciones:** OU estacionario en u = ln(1+z) sobre la malla de `cmb_fondo.py` (paso 2 × 10⁻³ hasta z*) con los
  siete z de DESI insertados como nodos exactos (3506 nodos), AR(1) con decaimiento por paso e^{−θΔu}; varianza y
  correlación comprobadas (0.090 frente a σ² = 0.09; e^{−θΔu} a 4 %, en 2000 caminos). 2000 por celda, semilla por celda.
- **Distancias exactas:** trapecio acumulado en u (exacto en los nodos de los datos), D_M y D_H en los 13 números,
  D_M(z*) con la corrección de z*(ω_m), θ* = 100 r_s(ω_m)/D_M(z*)/CAL, r_d(ω_m): todo de `escalon.py`. Con X = 0 reproduce
  `model15` de `cmb_fondo.py` a 3 × 10⁻⁷ y los ajustes ΛCDM: 15 filas χ² = 17.238, H0 = 68.19, ω_m = 0.14072
  (`ESCALON.md`: 17.24, 68.19, 0.14072); 13 filas χ² = 10.284, H0 = 68.67 (`CMB_FONDO.md`: 10.28, 68.67).
- **Descarte:** realización con H² ≤ 0 en algún nodo, evaluada en el mejor ajuste ΛCDM de las 15 filas; en el ajuste, los
  (h, ω_m) que dan H² ≤ 0 se penalizan.
- **Ajuste por realización:** (h, ω_m) por mínimos cuadrados a BAO + θ* + ω_m de Planck (15 filas, covarianza oficial 13 × 13
  más los dos errores de Planck) y a BAO solo (13 filas); χ²_min de cada uno.
- **Cota (i):** σ_X al que P(χ²_min ≤ χ²_obs) = 5 %, χ²_obs = 17.24 (15 filas) y 10.28 (BAO solo), interpolando
  linealmente en (ln σ_X, P) entre las seis celdas. Variante: los descartes contados como fallo (no cambia nada: los cruces
  ocurren antes de que haya descartes).
- **Cota (ii):** la lineal del convenio B de `CMB_FONDO.md`, marginalizada (−1.92 y −1.35) en las mismas celdas.

## Resultados

**BAO + θ* + ω_m (χ²_obs = 17.24):**

| θ | exacta, 5 % | lineal marg. −1.92 | lineal marg. −1.35 | exacta/lineal(−1.92) | P(χ²_min ≤ 17.24) en σ_X = 0.05 / 0.1 / 0.2 / 0.3 / 0.5 / 0.7 | descartes H² ≤ 0 (σ_X = 0.3 / 0.5 / 0.7) |
|---|---|---|---|---|---|---|
| 0.1 | sin cota (≤ 0.7) | sin cota | sin cota | — | 0.45 / 0.40 / 0.34 / 0.28 / 0.18 / 0.14 | 0 / 0 / 0 |
| 0.3 | sin cota (≤ 0.7; 5.4 % en 0.7) | sin cota | sin cota | — | 0.43 / 0.37 / 0.26 / 0.20 / 0.11 / 0.054 | 0 / 0 / 0 |
| 1 | **0.43** | 0.62 | 0.55 | 0.68 | 0.36 / 0.29 / 0.17 / 0.093 / 0.031 / 0.012 | 0 / 0 / 0.7 % |
| 3 | **0.29** | 0.38 | 0.33 | 0.76 | 0.34 / 0.24 / 0.10 / 0.043 / 0.010 / 0.003 | 0 / 2.5 % / 14 % |
| 10 | **0.19** | 0.27 | 0.23 | 0.71 | 0.32 / 0.18 / 0.043 / 0.018 / 0.004 / 0.003 | 0.9 % / 18 % / 43 % |

**BAO solo (χ²_obs = 10.28):**

| θ | exacta, 5 % | lineal marg. −1.92 | lineal marg. −1.35 | exacta/lineal(−1.92) | P(χ²_min ≤ 10.28) en σ_X = 0.05 / … / 0.7 |
|---|---|---|---|---|---|
| 0.1 | sin cota (≤ 0.7) | sin cota | sin cota | — | 0.48 / 0.46 / 0.43 / 0.36 / 0.26 / 0.19 |
| 0.3 | sin cota (≤ 0.7; 6.5 % en 0.7) | 0.89 | 0.74 | — | 0.47 / 0.45 / 0.34 / 0.25 / 0.14 / 0.065 |
| 1 | **0.47** | 0.48 | 0.40 | 0.97 | 0.45 / 0.37 / 0.21 / 0.11 / 0.041 / 0.012 |
| 3 | **0.28** | 0.29 | 0.24 | 0.97 | 0.38 / 0.28 / 0.12 / 0.037 / 0.012 / 0.006 |
| 10 | **0.20** | 0.20 | 0.17 | 0.97 | 0.36 / 0.21 / 0.046 / 0.015 / 0.004 / 0.001 |

Descartes: ninguno con θ ≤ 0.3 ni con σ_X ≤ 0.2; empiezan en θ = 10 a σ_X = 0.3 (0.9 %) y llegan al 43 % en (θ = 10, σ_X = 0.7)
y al 14 % en (θ = 3, 0.7). Todas las cotas quedan por debajo de la primera celda con descartes, salvo θ = 10, donde 0.19 está
justo bajo la celda 0.3 (0.9 %); contar los descartes como fallo no mueve ningún número.

**Lectura.**
- Con BAO solo la cota exacta y la lineal marginalizada coinciden al 3 %: la linealización y la marginalización gaussiana
  del fondo son buenas aquí, y la cota es la misma 0.47 / 0.28 / 0.20 (θ = 1 / 3 / 10) de `VERIFICACION.md` y
  `CMB_FONDO.md`.
- Con el CMB comprimido la exacta es más estrecha que la lineal por ×1.3–1.5 (0.43 frente a 0.62 en θ = 1). La lineal se
  aflojaba con el CMB porque su perfil tenía el máximo en σ_X > 0 (la tensión DESI–Planck absorbida en promedio, `CMB_FONDO.md`
  punto 4) y la cota se medía desde ahí; la exacta pregunta cuántas realizaciones concretas ajustan tan bien como ΛCDM,
  y con 15 filas y la tensión de por medio son menos: el CMB hace la cota más estrecha, no más ancha (0.43 frente a 0.47 de
  BAO solo), aunque poco.
- En θ ≤ 0.3 no hay cota por debajo de 0.7 en ninguno de los dos caminos; en θ = 0.3 y σ_X = 0.7 la fracción es 5.4 %
  (BAO + CMB) y 6.5 % (BAO), en el borde. Un residuo con memoria larga en ln(1+z) es, casi todo, Λ más una pendiente que el
  fondo absorbe.
- La cota exacta es monótona en θ (0.43 → 0.19) porque lo que se ve es la parte que varía entre bins, y esa crece con θ a
  σ_X fija; en ningún caso baja de 0.19, muy por encima del 10⁻² de la cabecera actual.

## Predicción: qué se cumplió

- Cota exacta en θ = 1 entre 0.2 y 0.5: **cumplida** (0.43).
- A menos de un factor 1.5 de la lineal: **cumplida, en el borde** (0.62/0.43 = 1.46).
- En θ ≤ 0.3 ninguna cota por debajo de 0.7: **cumplida** (P = 5.4 % y 14 % en σ_X = 0.7 para θ = 0.3 y 0.1).

Número para la v3.0 (convenio B, cota exacta): DESI DR2 + ω_m y θ* de Planck dejan pasar al 5 % un temblor OU de
σ_X < 0.43 (θ = 1), 0.29 (θ = 3), 0.19 (θ = 10), y no acotan θ ≤ 0.3 por debajo de 0.7.
