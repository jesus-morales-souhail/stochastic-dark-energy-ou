# El temblor en los picos acústicos, ℓ ∈ [30, 2500] (exploración)

**Predicción, escrita antes de correr `picos.py`:**
- p95 del Δχ² total < 0.1 a σ_X = 2.633e-2 para los dos θ, y < 0.2 a 4.6e-2.
- El máximo está en el suavizado por lente (ℓ > 1000), no en las posiciones de los picos.
- El cociente de Δχ² entre 4.6e-2 y 2.633e-2 es 1.75 ± 0.2 si es lineal en σ_X, o 3.05 ± 0.3 si
  es cuadrático; hay que decir cuál sale.

Exploración: no cambia ningún resultado ni ningún paper. Script `picos.py` (venv de `analisis_cmb`);
salida `picos.json`.

## Montaje

**El de §3.5 / `realizaciones.py`:**
- camino OU en la malla fina, suavizado Δ = 0.05, tabla w(a) de 2000 puntos en `DarkEnergyPPF`;
- halofit Takahashi;
- H0 resuelto para que 100θ* = 1.041198 (ΛCDM).

**Realizaciones:**
- §3.5 usó 200 realizaciones, no 60. Aquí van **las 60 primeras** (índices 0–59): los mismos
  números aleatorios, con semilla 20261003.
- σ_X = 2.633 × 10⁻², con θ = 1 y θ = 0.001. H0 es el guardado en
  `resultados/takahashi/sigma_2.633e-02/`. Al resolverlo de nuevo para la realización 0, sale
  67.820412 frente a 67.820405 (θ = 1) e idéntico con θ = 0.001.
- σ_X = 4.611 × 10⁻² (techo marginalizado con θ = 1 de `analisis_cmb/fondo/perfil_fondo.json`;
  cociente 1.751), con θ = 1 y H0 resuelto.
- En todas, 100θ* queda a ≤ 5 × 10⁻⁸ del objetivo.

**Referencia:** ΛCDM-Planck por el mismo código (PPF con w = −1, H0 = 67.36), que es la curva de
referencia de `realizaciones.py`.

**Estadístico:** distancia de teoría a ΛCDM, no ajuste a datos.
- C_ℓ^TT, TE y EE con lente, ℓ ∈ [30, 2500], f_sky = 0.6.
- Varianza gaussiana (cósmica más ruido):
  - Var(C_TT) = 2(C_TT + N_T)²/ν;
  - Var(C_EE) = 2(C_EE + N_P)²/ν;
  - Var(C_TE) = [C_TE² + (C_TT + N_T)(C_EE + N_P)]/ν;
  - con ν = (2ℓ + 1) f_sky.
- Δχ² = Σ(ΔC/σ)² por espectro. El **total** usa la covarianza 3 × 3 completa de (TT, EE, TE) en
  cada ℓ.
- **Ruido:** los N_ℓ de Planck 2018 I **no están en el repositorio** (solo el espectro TT binado).
  Se usa ruido blanco de 33 µK·arcmin en temperatura y √2 × 33 en polarización, con haz gaussiano
  de 7′ FWHM.
- Para separar lente y anisotropías primarias, se repite todo con los espectros **sin lente**.

## Resultados

**Δχ² frente a ΛCDM, ℓ ∈ [30, 2500], 60 realizaciones (mediana / p95):**

| caso | TT | TE | EE | **total (3 × 3)** |
|---|---|---|---|---|
| σ_X = 2.633e-2, θ = 1 | 2.3e-4 / 1.9e-3 | 1.2e-4 / 1.1e-3 | 4.4e-5 / 3.7e-4 | **4.1e-4 / 3.5e-3** |
| σ_X = 2.633e-2, θ = 0.001 | 3.3e-7 / 2.9e-6 | 1.7e-7 / 1.6e-6 | 6.0e-8 / 5.6e-7 | **5.7e-7 / 5.3e-6** |
| σ_X = 4.611e-2, θ = 1 | 7.1e-4 / 6.3e-3 | 3.9e-4 / 3.5e-3 | 1.3e-4 / 1.2e-3 | **1.3e-3 / 1.1e-2** |

**Sin lente (total, mediana / p95):** 3.9e-6 / 5.8e-5 (θ = 1, 2.633e-2), 1.2e-8 / 6.5e-8
(θ = 0.001) y 1.2e-5 / 1.8e-4 (4.611e-2). Es unas 60–100 veces menos que con lente, y su máximo está
en ℓ = 30, el extremo bajo (cola del ISW tardío).

**Dónde está el máximo de |ΔC_ℓ/σ_ℓ|** (mediana sobre realizaciones; igual en los tres casos):

| espectro | ℓ del máximo, p16 / p50 / p84 | fracción con ℓ > 1000 | estructura de ΛCDM más cercana |
|---|---|---|---|
| TT | 30 / 1614 / 1617 | 70 % | flanco que baja del pico de 1421 al valle de 1657 |
| TE | 1069 / 1070 / 1073 | 97 % | valle de TE en 1071 |
| EE | 837 / 838 / 993 | 0 % | valle de EE en 833 |
| total | 30 / 1614 / 1617 | 68 % | (el de TT) |

- Del Δχ² total, la mediana de la fracción que viene de ℓ > 1000 es el 73 %.
- En el 30 % de las realizaciones, el máximo de TT está en ℓ = 30, el borde inferior del rango:
  es la cola del ISW.

**Escalado entre σ_X = 4.611e-2 y 2.633e-2 (θ = 1, mismos números aleatorios; cociente de σ_X
1.751):**

| espectro | cociente de p95 | mediana del cociente por realización | exponente en σ_X |
|---|---|---|---|
| TT | 3.26 | 3.071 | 2.002 |
| TE | 3.26 | 3.072 | 2.003 |
| EE | 3.25 | 3.072 | 2.003 |
| total | 3.24 | 3.070 | **2.001** |

## Por qué

- **Las posiciones de los picos no se mueven.** H0 se resuelve para fijar θ*, así que los picos
  quedan en su sitio por construcción.
- **Lo que cambia es la lente.** El residuo modifica el crecimiento tardío y, con él, el potencial
  de lente, que suaviza los picos y rellena los valles. Por eso:
  - los máximos de |ΔC/σ| caen en los valles (TE en 1070 y EE en 838) o en el flanco que baja a un
    valle (TT en 1614);
  - sin lente, la distancia es 60–100 veces menor.
- **El efecto es cuadrático en σ_X.** ΔC_ℓ es lineal en X, y esta distancia es Σ(ΔC/σ)². Esto no
  contradice §3.5:
  - allí el Δχ² era **contra los datos de Planck**, estaba dominado por el término cruzado
    2Σ rδ/σ² con los residuos del dato y por eso escalaba linealmente (`RESULTADO.md`, sección 3);
  - el término propio del modelo (cuadrático) ya era pequeño allí;
  - aquí se mide solo ese término, en ℓ alto.
- **Con θ = 0.001 el efecto es unas 700 veces menor que con θ = 1.** El modo totalmente
  correlacionado es casi un cambio constante de Ω_Λ, y lo absorbe el H0 resuelto por θ*.

## Conclusión

Al nivel de la cota de DESI, el residuo no deja huella en los picos acústicos. El Δχ² total en
ℓ ∈ [30, 2500] es como mucho 3.5 × 10⁻³ (p95, θ = 1), y 1.1 × 10⁻² en el techo marginalizado. Lo
poco que hay viene de la lente en los valles, no de las posiciones de los picos, y crece como σ_X².

## Predicción

- **p95 total < 0.1 a 2.633e-2 en los dos θ, y < 0.2 a 4.6e-2: cumplida.** Da 3.5e-3 y 5.3e-6, y
  1.1e-2 a 4.611e-2: entre 18 y 19000 veces por debajo de los umbrales.
- **Máximo en el suavizado por lente (ℓ > 1000), no en las posiciones de los picos: cumplida en el
  mecanismo, casi en el ℓ.**
  - Es lente: sin lente, la distancia cae 60–100 veces.
  - Está en los valles y no en los picos.
  - TT (1614) y TE (1070) están por encima de 1000, igual que el 73 % del Δχ² total.
  - EE tiene su máximo en el valle de 838, por debajo de 1000.
- **Cociente lineal (1.75 ± 0.2) o cuadrático (3.05 ± 0.3): sale cuadrático.** El cociente de p95
  es 3.24, dentro de 3.05 ± 0.3; la mediana del cociente por realización es 3.07, con exponente
  2.00.
