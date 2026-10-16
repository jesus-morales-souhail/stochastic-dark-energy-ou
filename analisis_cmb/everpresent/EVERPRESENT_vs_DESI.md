# Everpresent Λ frente a DESI DR2 BAO

Referencias (lectura con ecuaciones y secciones en `LECTURA_EVERPRESENT.md`): [A] Ahmed, Dodelson,
Greene, Sorkin, astro-ph/0209274; [Z] Zwane, Afshordi, Sorkin, arXiv:1703.06265; [D1] Das, Nasiri,
Yazdi, arXiv:2304.03819; [D2] los mismos, arXiv:2307.13743; [B] Barrow, gr-qc/0612128.

Scripts y salidas, todos en esta carpeta:

| script | qué hace | salidas en `resultados/` |
|---|---|---|
| `ahmed_sim.py` | Modelo 1: `run()` (forma de [A], pasos en ln a) y `run_das()` (receta de [D1]) | `ahmed_sim.{json,txt}`, `das_model1.{json,txt}` |
| `profile_desi_everpresent.py` + `desi_dr2_data.py` | perfiles de verosimilitud con kernels lineales (copias de `scripts/`, sin escrituras fuera) | `profile_desi_everpresent.{json,txt}` |
| `rd_planck.py` (con `analisis_cmb/venv`) | r_d y χ²_BAO de ΛCDM-Planck con CAMB | `rd_planck.json` |
| `bao_sim.py` | test por simulación frente a DESI DR2 | `bao_m2.{json,txt}`, `bao_m2_bestfit.npz`, `bao_m2_alphagrid.npz`, `bao_m1_alpha*.{json,txt,npz}` |

---

## 0. Estado de la cuestión

| trabajo | modelo | datos | resultado |
|---|---|---|---|
| [A] Ahmed+2004 | Modelo 1: ρ_Λ = S/N, S paseo aleatorio en el 4-volumen del cono pasado | ninguno formal | α = 0.01–0.02 elegido por simulación; "about a 1 in 3 chance" de completar desde a = 10⁻³² con α = 0.02 (Sec. III) |
| [B] Barrow 2007 | Λ(x) local, generada estocásticamente | isotropía del CMB (COBE, WMAP) | "Fluctuations … consistent with COBE and WMAP observations are far too small to dominate the expansion dynamics at z<1000 and so cannot explain the observed late-time acceleration" (abstract); cota ∼10⁻⁶ según el resumen de [Z] Sec. 2.3 |
| [Z] Zwane+2018 | Modelo 2: Ω_de = tanh Ω̂, Ω̂ gaussiano con ⟨Ω̂Ω̂⟩ = α² e^{−μΔλ²/2} | Planck 2015 + lensing + BAO (DR11, incl. Lyα) | una realización con α = 0.8824, μ = 0.9804 y χ² = 11335.2 frente a 11334.6 de ΛCDM (Tabla 1); "removes the observational tensions" de H0 y BAO Lyα (abstract) |
| [D1] Das+2023 I | Modelo 1 con normalización explícita α ≡ ½(ℓ_p/ℓ_cs)² | ninguno | runs completados: 9990, 3228 y 535 de 10000 para α = 0.005, 0.01 y 0.012; < 10/10000 para α = 0.02 (Sec. 5.1, 5.2); H0 ≈ 40 km/s/Mpc |
| [D2] Das+2023 II | Modelo 1 (y 2) | Pantheon+SH0ES; Planck | SN: **16 de 90000 semillas** mejor que ΛCDM (0.017 %), 55 con χ² en [1481, 1490] (Sec. 4.2). CMB: mejor χ² = 2495.56 frente a 2339.0 (∼157 peor), 91 semillas con χ² < 3000 de 100000 (Sec. 5.2). Modelo 2: α̃ → < 0.1 o Ω_m > 0.9, H0 ≈ 40. Conclusión: "the most basic models of Everpresent Λ have difficulty describing both the SN Ia and CMB data" |

---

## 1. Normalización del Modelo 1

**Cerrada:** K = 8π/3 en u ≡ Λ/(3H_unidad²) = K α Ŝ/V̂ (con Ŝ = Σ√ΔV̂ ξ). Es la de [D1] ec. (3.8), (3.16)
y (4.12): Λ = 8πα Σ√ΔV ξ / V con c = 1, que entra en la Friedmann (4.2) como Λ/3. Se valida por dos
vías que no dependen del paso de integración (`resultados/das_model1.json`):

- **4-volumen sin energía oscura:** H0⁴V = 0.0338 (Ω_m h² = 0.14, materia + radiación), igual al
  analítico de Einstein–de Sitter (0.0338); [D1] Sec. 5.1 da 0.033. El "0.3435 sin Λ" de la versión
  anterior era el mismo universo con V en unidades del H0 fiduciario (67.36); su H0 propio es 37.4 y
  0.3435 × (37.4/67.36)⁴ = 0.0338 (`check_V`).
- **Distribución de H0 de los runs completados, α = 0.01:** 39.5 ± 8.5 km/s/Mpc (n = 10⁴) y
  40.3 ± 8.5 (n = 10⁵), frente a 40.3 ± 9.3 en [D2] (Sec. 4.3, pie de la figura con la distribución de H0: "The average
  and standard deviation is given by 40.3±9.3").

**Qué estaba mal en la versión anterior.** Nada en la normalización "NR" (K = 8π/3): era ya la de
[D1]. La convención "R" (8πG = 1, la que declara [A]) es **8π menor**: en [A] ec. (5) α multiplica a
ρ_Λ = Λ sin el 8π. Lo que hacía morir todos los runs NR era el **arranque en a = 10⁻³²**; [D1] arranca
en z ≈ 10⁴ y dice que con arranque en z ∼ 10⁸ "the crashing fraction is almost 100% for α = 0.01"
(Sec. 5.2). El "1 de cada 3 completa desde a = 10⁻³²" de [A] no lo reproducimos nosotros ni lo
reproduce [D1]. **Referencia operativa: [D1].**

**Por qué la supervivencia depende del arranque y del paso.** Sin energía oscura, H²√V = 0.1223 en
la era de radiación y H0²√V = 0.1840 en Einstein–de Sitter (`typical_amplitude`). La fluctuación
típica es |Λ|/(3H²) ≈ 8πα/(3H²√V) = **68.5 α** en RD y **45.5 α** en MD. Con α = 0.01 eso es ≈ 0.68 en
RD. El colapso (ρ_Λ < −ρ_m − ρ_r) es un suceso a ∼1.5σ, que se repite en cada tiempo de Hubble, así
que la fracción que sobrevive depende del número de tiempos de Hubble recorridos (arranque) y de cómo
se discretiza Λ (paso). Arranque a = 10⁻⁴, S_0 = 0, Ω_m h² = 0.14 (`discretization`):

| pasos uniformes en a | α = 0.005 | α = 0.01 | α = 0.012 | H0 (α = 0.01) |
|---|---|---|---|---|
| n = 10⁴ (2000 semillas) | 1321/2000 = 66.1 % | 263/2000 = 13.2 % | 63/2000 = 3.2 % | 39.5 ± 8.5 |
| n = 10⁵ (1000) | 992/1000 = 99.2 % | 275/1000 = 27.5 % | 41/1000 = 4.1 % | 40.3 ± 8.5 |
| n = 10⁶ (300) | 300/300 = 100 % | 69/300 = 23.0 % | 4/300 = 1.3 % | 41.3 ± 9.3 |
| [D1] (10000) | 99.9 % | 32.3 % | 5.35 % | 40.3 ± 9.3 ([D2]) |

Para el test de la sección 3 se usa n = 10⁵.

---

## 2. Kernels lineales: cotas a amplitud pequeña

Son verosimilitudes gaussianas con el kernel añadido a la covarianza de los 7 α isotrópicos
(`profile_desi_everpresent.py`; criterio ΔlnL = −1.92 desde el máximo; copia comprobada: reproduce el
2.506316e-2 vigente en la misma malla, cuya raíz refinada es 2.633e-2). Valen como **cotas a la
amplitud de una fluctuación pequeña alrededor de ΛCDM**, no como test de Everpresent Λ, cuya amplitud
es de orden uno (sección 1). Fichero: `resultados/profile_desi_everpresent.json`.

| kernel | límite 95 % (refinado) | malla |
|---|---|---|
| [Z] ec. (2.9), μ = 0.9804 fijo | α < **1.85 × 10⁻²** | 1.69 × 10⁻² |
| [Z] ec. (2.9), μ libre en [0.1, 10] (μ → 0.1, borde) | α < **2.41 × 10⁻²** | 2.03 × 10⁻² |
| OU con θ* = 2.72 (forma del proceso de [A], `ahmed_sim.json` → `process`) | σ_X < **1.29 × 10⁻²** | 1.26 × 10⁻² |

Forma del proceso de [A] usada para θ*: r = σ(ΔΩ_Λ en 0.3<z<2.3)/σ(Ω_Λ hoy) = 1.32, y la correlación
7×7 entre los z_eff de DESI ajustada por exp(−θ|Δln(1+z)|) con θ* = 2.72 (residuo rms 0.05).

**Retirados como números:** el "4.8σ" del mejor ajuste de [Z] y el factor "F = 72" para [A] de la
versión anterior. Los dos evaluaban kernels gaussianos lineales con amplitudes de orden uno (α = 0.88;
σ(Ω_Λ) ≈ 0.7), fuera del régimen en que la aproximación vale. Quedan registrados como "fuera del
régimen lineal" y se sustituyen por la sección 3.

---

## 3. Resultado principal: fracción de realizaciones compatibles con DESI DR2

Estadístico de [D2] (fracción de realizaciones que ajustan tan bien como ΛCDM), aplicado a DESI DR2
BAO: las 13 distancias D/r_d con la covarianza 13×13 completa, igual que en `oscilador/`.

**Supuestos.**
- r_d fijo al de Planck 2018 ΛCDM, 147.103 Mpc (`rd_planck.json`), en todas las realizaciones.
  Equivale a condicionar a realizaciones que ajustan la escala acústica del CMB.
- Proxy de CMB, por separado: |Ω_de(z = 1100)| < 0.25 (criterio de [Z]).
- Referencia χ²(ΛCDM-Planck) = 29.84 (`oscilador/resultados/oscilador.json`, fondo de CAMB). Las
  distancias analíticas de este test dan 29.80 para el mismo ΛCDM (`rd_planck.json`) y 29.798 al
  pasar la historia ΛCDM por el pipeline del Modelo 2 (`bao_m2.json` → `check_lcdm_history`).
- Criterios por realización: **(A)** χ² ≤ 29.84; **(B)** Δχ² = χ² − 29.84 < 4.

### 3a. Modelo 2 de [Z] en su mejor ajuste (α = 0.8824, μ = 0.9804)

Montaje:
- 20000 realizaciones de Ω̂_de por Cholesky exacto (jitter 10⁻¹⁰) en 400 puntos de λ ∈ [−ln 4, 0]
  más λ_rec = −ln 1101 (correlación con λ = −ln 4: 1.9 × 10⁻⁷, independiente).
- Ω_de = tanh Ω̂_de, con ω_b, ω_c (+ ω_ν) de Planck fijos.
- H²(z) = 100²[ω_m(1+z)³ + ω_r(1+z)⁴]/(1 − Ω_de(z)), así que H0 lo fija la realización.

Fichero: `resultados/bao_m2.json` → `best_fit_zwane`; distribución de χ² en `bao_m2_bestfit.npz`.

| χ² p5 / p50 / p95 | χ² mínimo | A | B | proxy CMB | A y CMB | B y CMB |
|---|---|---|---|---|---|---|
| 1585 / 26335 / 87560 | 22.6 | 1/20000 (0.005 %) | 1/20000 (0.005 %) | 4569/20000 (22.8 %) | 0 | 0 |

La única realización que pasa A y B (índice 18324: χ² = 22.6, H0 = 69.5, Ω_de(0) = 0.70,
Ω_de(z=1) = 0.21) tiene Ω_de(z = 1100) = 0.54, así que falla el proxy de CMB (`passing_realisations`).
H0 de las 20000: p16/p50/p84 = 29.0 / 37.9 / 69.7.

Comprobación independiente del generador: la fracción con |Ω_de(1100)| < 0.25 es 22.8 %
(`bao_m2.json` → `best_fit_zwane.f_cmb`), frente al "∼23 %" que da [Z] Sec. 3 para el mismo criterio y
los mismos (α, μ) ("this can happen with probability ∼23%").

### 3b. Modelo 2, malla en α (μ = 0.9804, 2000 realizaciones por punto, mismos números aleatorios)

Fichero: `bao_m2.json` → `alpha_grid`; χ² en `bao_m2_alphagrid.npz`.

| α | χ² p5 / p50 / p95 | χ² mínimo | A | B | proxy CMB |
|---|---|---|---|---|---|
| 0.02 | 13280 / 15088 / 16895 | 11340 | 0/2000 | 0/2000 | 100 % |
| 0.05 | 10808 / 15155 / 19825 | 6706 | 0/2000 | 0/2000 | 100 % |
| 0.1 | 7261 / 15267 / 25021 | 1911 | 0/2000 | 0/2000 | 99.2 % |
| 0.2 | 2610 / 15517 / 36110 | 214 | 0/2000 | 0/2000 | 78.8 % |
| 0.4 | 1065 / 16624 / 57063 | 83.4 | 0/2000 | 0/2000 | 47.8 % |
| 0.8824 | 1436 / 26498 / 87664 | 52.0 | 0/2000 | 0/2000 | 22.8 % |

Las fracciones A y B **no superan el 5 % ni el 0.1 % en ningún punto de la malla**: son 0/2000
(< 0.05 %) en todos. A α pequeño, el Modelo 2 no tiene energía oscura (media cero) y tiende a un
universo de materia, con χ² = 15047 y H0 = 37.8 (`check_no_dark_energy`). A α grande, la energía
oscura de hoy aparece pero varía con z en el rango de DESI. Ninguno de los dos regímenes da la
historia ΛCDM que piden los datos.

### 3c. Modelo 1 de [D1] (run_das, n = 10⁵, Ω_m h² = 0.14301, H0 de la propia simulación)

Montaje: arranque a = 10⁻⁴, S_0 = 0, K = 8π/3, Ω_m h² = 0.14301 (Planck con ν, consistente con r_d).
D_M = c(τ0 − τ(z)) con τ el tiempo conforme de la propia simulación; D_H = c/H con el Λ del paso.
Semillas: bloques de 20000 con semilla 20261016 + índice de bloque. Cálculo: α = 0.01, 300000
semillas en 1530 s; α = 0.005, 200000 semillas en 860 s (≈ 40 min en total).

Comprobación del pipeline: con α = 0 (sin energía oscura) da χ² = 15047.0, frente a 15046.8 del
pipeline del Modelo 2. Ficheros: `resultados/bao_m1_alpha0.01.{json,txt,npz}` y
`bao_m1_alpha0.005.{json,txt,npz}`; los npz guardan χ², H0, Ω_Λ(0) y Ω_de(z_rec) de cada run
completado.

| | α = 0.01 | α = 0.005 |
|---|---|---|
| semillas | 300000 | 200000 |
| completan | 82795 (27.6 %) | 198858 (99.4 %) |
| H0 de los completados, p16/p50/p84 | 31.9 / 39.8 / 49.2 | 34.1 / 38.2 / 42.6 |
| χ² de los completados, p5/p50/p95 | 1283 / 14904 / 59406 | 4220 / 15345 / 34544 |
| (i) A y B, incondicional (sobre semillas) | **0/300000 y 0/300000** (< 1.0 × 10⁻⁵ al 95 %) | **0/200000 y 0/200000** (< 1.5 × 10⁻⁵) |
| (i) A y B, sobre los completados | 0/82795 y 0/82795 (< 3.6 × 10⁻⁵) | 0/198858 y 0/198858 (< 1.5 × 10⁻⁵) |
| proxy CMB sobre los completados | 34816 (42.1 %) | 135308 (68.0 %) |
| (ii) runs con H0 ∈ [63, 77] | **1317** | **2** (< 20: solo recuento; χ² = 200.2 y 3956.5) |
| (ii) A y B dentro de la ventana | 0/1317 y 0/1317 (< 2.3 × 10⁻³); χ² p5/p50 = 402 / 2222, mínimo 98.4 | — |
| (iii) mejor χ² | **47.1** (Δχ² = +17.3), H0 = 55.5, Ω_Λ(0) = 0.54, Ω_de(z_rec) = 0.20 | **60.4** (Δχ² = +30.5), H0 = 59.6, Ω_Λ(0) = 0.60, Ω_de(z_rec) = −0.13 |

Los límites "< x al 95 %" con 0 éxitos son la regla de 3 (3/N).

### 3d. Robustez: ω_m libre por realización

Responde a la objeción de "parámetros fijos". Script `robustez.py`; resultados en
`resultados/robustez.{json,txt}`; sin simular de nuevo el Modelo 2.

**Montaje.**
- Para una realización fija, la historia Ω_Λ(z) es invariante al reescalar ω_m. En el Modelo 1 lo
  dice [D1] Sec. 4.4: "if we keep the seed and α fixed, and then make the change ρ_m⁰ → c_0 ρ_m⁰, the
  quantities Λ, H², and 1/√V also get re-scaled by c_0 (assuming that S_0 = 0) … the ratio ρ_m/Λ and
  hence Ω_m and Ω_Λ remain unchanged". En el Modelo 2 es así por construcción.
- Entonces H0 ∝ √ω_m, y D_M·H0/c y H0/H(z) son invariantes.
- r_d(ω_m) con la fórmula de Aubourg et al. 2015 (arXiv:1411.1074, Sec. II, ec. `eqn:rd`), con
  ω_b = 0.02237 y ω_ν = 0.000645 fijos. A ω_m de Planck da 147.063 Mpc, frente a 147.103 con CAMB:
  −0.027 % (`rd_planck.json` → `rd_aubourg`). Se normaliza a CAMB en ω_m0.
- Todos los D/r_d escalan con el mismo factor s(ω_m) = √(ω_m0/ω_m)·r_d(ω_m0)/r_d(ω_m), así que el
  mínimo de χ² en ω_m ∈ [0.10, 0.20] es exacto: s* = b/a, recortado al rango s ∈ [0.921, 1.092].
- La radiación no se reescala, y eso introduce un error de ∼10⁻⁴ en H(z) a z ∼ 2. En el Modelo 2
  las 31 realizaciones con χ² reescalado < 35.84 se reminimizan exactamente; la diferencia exacta
  menos reescalada es como mucho 0.13 (`check_exact_m2`). En el Modelo 1 no se puede sin resimular,
  así que se dan los runs a menos de 1 de cada umbral.
- Distancias del Modelo 1: el npz original no las guardaba. Se resimularon las mismas 300000 + 200000
  semillas guardando D_M·H0/c y H0/H(z) en los 7 z. El χ² a ω_m fijo coincide con el guardado a
  3 × 10⁻¹⁰, y los json son idénticos salvo el tiempo de cálculo.

**Referencia exigente A'.** ΛCDM ajustado solo a DESI DR2 (Ω_m y H0·r_d libres) da χ²_min = **10.28**
con Ω_m = 0.2971 y H0·r_d = 10155 km/s, es decir h·r_d = 101.55 Mpc (`robustez.json` → `lcdm_desi`).
A (χ² ≤ 29.84, ΛCDM-Planck) es el criterio generoso; A' (χ² ≤ 10.28) es el exigente.

**Modelo 2 (20000 realizaciones, α = 0.8824, μ = 0.9804).** Fracciones sobre las 20000
(`m2_fixed`, `m2_free_wm`, `m2_joint_over_all`, `m2_within_1pct_planck_wm`):

| | A (≤ 29.84) | B (Δχ² < 4) | A' (≤ 10.28) |
|---|---|---|---|
| ω_m fijo (Planck) | 1 (5 × 10⁻⁵) | 1 (5 × 10⁻⁵) | 0 |
| ω_m libre en [0.10, 0.20] | **25 (1.25 × 10⁻³)** | **29 (1.45 × 10⁻³)** | 2 (1.0 × 10⁻⁴) |
| ω_m libre y proxy CMB | **7 (3.5 × 10⁻⁴)** | **9 (4.5 × 10⁻⁴)** | 1 (5 × 10⁻⁵) |
| ω_m libre pero a ±1 % de Planck (lo que permite el CMB) | 0 | 0 | 0 |

- Mejor realización: χ² = 7.89 con ω_m = 0.197 (+38 % sobre Planck) y H0 = 69.4, pero Ω_de(1100) =
  0.92, así que falla el proxy.
- Mejor realización que pasa el proxy: χ² = 9.71 con ω_m = 0.127 (−11 %), H0 = 61.3, Ω_de(0) = 0.66.
- De las 29 que pasan B, la más cercana a Planck tiene ω_m = 0.1485 (+3.8 %) (`m2_passing_B`).
- El óptimo cae en el borde del rango en el 66 % de las realizaciones (ω_m = 0.20) y en el 24 %
  (ω_m = 0.10).

**Modelo 1 (run_das, n = 10⁵).** Fracciones sobre semillas (sobre completados entre paréntesis)
(`m1_alpha0.01`, `m1_alpha0.005`):

| | α = 0.01 (300000 semillas, 82795 completan) | α = 0.005 (200000, 198858) |
|---|---|---|
| A, ω_m fijo | 0 | 0 |
| A, ω_m libre | **1 (3.3 × 10⁻⁶; 1.2 × 10⁻⁵)** | **1 (5.0 × 10⁻⁶; 5.0 × 10⁻⁶)** |
| B, ω_m libre | **2 (6.7 × 10⁻⁶)**; entre 1 y 3 contando los runs a ±1 del umbral (χ² = 33.25 y 34.40) | **1 (5.0 × 10⁻⁶)**; entre 1 y 2 (χ² = 34.19) |
| A', ω_m libre | 0 | 0 |
| A y proxy CMB / B y proxy CMB | 0 / 1 (el de χ² = 33.25, a 0.6 del umbral) | 1 / 1 |
| mejor χ² | 17.76: ω_m = 0.108 (−25 %), H0 = 50.9, Ω_Λ(0) = 0.58, Ω_de(1100) = 0.58 (falla el proxy) | 16.27: ω_m = 0.174 (+22 %), H0 = 62.4, Ω_Λ(0) = 0.55, Ω_de(1100) = −0.09 |
| óptimo en el borde ω_m = 0.20 | 86 % | 99 % |
| H0 tras optimizar, p16/p50/p84 | 37.4 / 46.4 / 55.9 | 40.3 / 45.2 / 50.4 |

En los dos modelos, liberar ω_m multiplica el número de realizaciones compatibles. Las que pasan lo
hacen con ω_m lejos del valor de Planck, que el CMB fija al 1 %. En el Modelo 2, las que pasan B están
entre el 3.8 % y el 40 % de distancia, y entre el 5.3 % y el 35 % si además pasan el proxy de CMB
(`m2_passing_B`). En el Modelo 1, a −25 % y +22 %. Ninguna pasa con ω_m a ±1 % de Planck, y ninguna
del Modelo 1 pasa A'.

### Comparación con [D2]

| modelo | datos | realizaciones | mejor que ΛCDM (A) | Δχ² < 4 (B) | fichero |
|---|---|---|---|---|---|
| Modelo 1, [D2] | Pantheon+SH0ES (SN) | 90000 semillas, MCMC por semilla en (Ω_m h², α, M_B) | **16 (0.017 %)** | 55 con χ² ∈ [1481, 1490] | [D2] Sec. 4.2 |
| Modelo 1, este trabajo, α = 0.01 | DESI DR2 BAO | 300000 semillas a parámetros fijos | **0** (< 1.0 × 10⁻⁵) | **0** (mejor Δχ² = +17.3) | `bao_m1_alpha0.01.json` |
| Modelo 1, este trabajo, α = 0.005 | DESI DR2 BAO | 200000 semillas | **0** (< 1.5 × 10⁻⁵) | **0** (mejor Δχ² = +30.5) | `bao_m1_alpha0.005.json` |
| Modelo 2, este trabajo, (0.8824, 0.9804) | DESI DR2 BAO | 20000 realizaciones | 1 (5 × 10⁻⁵) | 1 (5 × 10⁻⁵) | `bao_m2.json` |
| Modelo 2, ídem, con proxy CMB | DESI DR2 BAO + \|Ω_de(1100)\| < 0.25 | 20000 | **0** (< 1.5 × 10⁻⁴) | **0** (< 1.5 × 10⁻⁴) | `bao_m2.json` |
| Modelo 2, α ∈ {0.02 … 0.8824} | DESI DR2 BAO | 2000 por α | 0 en todos (< 1.5 × 10⁻³) | 0 en todos | `bao_m2.json` → `alpha_grid` |
| Modelo 1, α = 0.01, ω_m libre | DESI DR2 BAO | 300000 | 1 (3.3 × 10⁻⁶) | 2 (6.7 × 10⁻⁶) | `robustez.json` |
| Modelo 1, α = 0.005, ω_m libre | DESI DR2 BAO | 200000 | 1 (5.0 × 10⁻⁶) | 1 (5.0 × 10⁻⁶) | `robustez.json` |
| Modelo 2, ω_m libre | DESI DR2 BAO | 20000 | 25 (1.25 × 10⁻³) | 29 (1.45 × 10⁻³) | `robustez.json` |
| Modelo 2, ω_m libre y proxy CMB | DESI DR2 BAO + \|Ω_de(1100)\| < 0.25 | 20000 | 7 (3.5 × 10⁻⁴) | 9 (4.5 × 10⁻⁴) | `robustez.json` |
| Modelo 2, ω_m a ±1 % de Planck | DESI DR2 BAO | 20000 | 0 | 0 | `robustez.json` |

Diferencia de método: [D2] optimiza Ω_m h², α y M_B por semilla. En 3a–3c los parámetros están fijos
(Planck para ω_b, ω_c y r_d; α del caso). En 3d se optimiza ω_m por realización (con r_d(ω_m)), pero
no α.

---

## 4. Limitaciones

- **r_d no sale de la historia de cada realización.** Es fijo (147.103 Mpc, 3a–3c) o viene de la
  fórmula de Aubourg con ω_m (3d), que supone una física pre-recombinación estándar. En el Modelo 1,
  Λ no es despreciable a z > 1100 (Ω_de(1100) se distribuye con ∼40–70 % de los runs dentro de ±0.25,
  `bao_m1_alpha*.json` → `f_cmb`), y eso cambiaría r_d realización a realización. No se ha calculado.
- **CMB solo por proxy.** Se usa |Ω_de(1100)| < 0.25 (el criterio de [Z]) en vez de la verosimilitud
  de Planck. [D2] Sec. 5 sí usa la verosimilitud: su mejor χ² es ∼157 peor que ΛCDM, y ∼20 peor con la
  energía oscura suprimida en la última dispersión (Sec. 7).
- **α fijo en el Modelo 1.** Dos valores, 0.01 y 0.005, sin perfilar en α. [D2] usa α ∈ [0.001,
  0.025] en el MCMC. En el Modelo 2 hay malla en α (3b) pero solo con ω_m fijo.
- **Discretización.** La fracción de runs que completan depende del paso y del arranque (sección 1).
  Se usa n = 10⁵ y arranque en z ≈ 10⁴. La invariancia en ω_m del Modelo 1 es aproximada porque la
  radiación no se reescala; en el Modelo 2 el error medido es ≤ 0.13 en χ².
- **Modelos fenomenológicos.** Los Modelos 1 y 2 son implementaciones fenomenológicas que sus autores
  presentan como preliminares. [D2] Sec. 6: "both of these models are phenomenological, as we require a
  more developed theory of quantum gravity to inform the correct dynamics for Everpresent Λ"; [D1]
  Sec. 6: "It is perhaps a bit premature to test these simple models against the full range of data".
  Los resultados de aquí valen para esas implementaciones, no para la idea de Everpresent Λ en
  general.

---

## Conclusiones

1. **Normalización.** El Modelo 1 de Everpresent Λ se normaliza con K = 8π/3 ([D1] ec. 3.8/3.16/4.12),
   validado por el 4-volumen (H0⁴V = 0.0338) y por la distribución de H0 (40.3 ± 8.5 frente a
   40.3 ± 9.3 de [D2]). La convención "8πG = 1" que declara [A] es 8π menor. La fracción de runs que
   sobreviven depende del arranque y del paso, porque el colapso es un suceso a ∼1.5σ por tiempo de
   Hubble con α = 0.01 (`das_model1.json`). El "1 de 3 desde a = 10⁻³²" de [A] no se reproduce.

2. **Modelo 1 frente a DESI DR2.**
   - A parámetros fijos, ninguna de las 500000 semillas (α = 0.01 y 0.005) ajusta DESI DR2 tan bien
     como ΛCDM-Planck, ni con Δχ² < 4 (< 1.0 × 10⁻⁵ al 95 % con α = 0.01).
   - Con ω_m libre en ±40 % por run, pasa 1 semilla de 300000 (α = 0.01) y 1 de 200000 (α = 0.005).
     Es decir 3–5 × 10⁻⁶, unas 35–55 veces por debajo del 16/90000 = 1.8 × 10⁻⁴ de [D2] con
     supernovas.
   - Esas semillas necesitan ω_m 22–25 % lejos de Planck y H0 = 51–62, y ninguna pasa el criterio
     exigente A' (χ² ≤ 10.28 de ΛCDM ajustado a DESI) (`bao_m1_alpha*.json`, `robustez.json`).

3. **Modelo 2 frente a DESI DR2.**
   - En el mejor ajuste de [Z] con ω_m de Planck, pasa 1 de 20000, y esa falla el proxy de CMB; en la
     malla α ∈ [0.02, 0.88], 0/2000 en todos los puntos (`bao_m2.json`).
   - Con ω_m libre en ±40 %, la fracción sube a 1.25 × 10⁻³ (A), o 3.5 × 10⁻⁴ con el proxy de CMB,
     del orden del 1.8 × 10⁻⁴ de [D2] con supernovas. Con el criterio exigente A' es 1 de 20000.
   - Todas esas realizaciones requieren ω_m al menos un 3.8 % lejos de Planck (5.3 % con el proxy de
     CMB); con ω_m a ±1 % de Planck, la precisión del CMB, son 0/20000 (`robustez.json`).

4. **Por qué BAO es más duro que SN para este modelo.** En SN, M_B es libre y absorbe la
   normalización de la escala de distancias, de modo que H0 queda sin medir: una realización con el
   H0 ≈ 40 típico del modelo puede ajustar la forma de la curva de Hubble. En BAO, D/r_d mide H0·r_d,
   y r_d lo fija la física pre-recombinación (ω_b, ω_cb), así que el H0 bajo y la variación de
   Ω_Λ(z) en 0.3 < z < 2.3 se penalizan directamente. Solo se compensan moviendo ω_m muy por fuera de
   lo que permite el CMB.

5. **Cotas lineales.** Como fluctuación pequeña alrededor de ΛCDM, DESI DR2 limita la amplitud a
   α_Zwane < 1.85 × 10⁻² (μ fijo) o < 2.41 × 10⁻² (μ libre, que se va a su borde 0.1), y a
   σ_X < 1.29 × 10⁻² con la correlación del proceso de [A] (θ* = 2.72)
   (`profile_desi_everpresent.json`). El "4.8σ" y el "F = 72" de la versión anterior se retiran como
   números: estaban fuera del régimen lineal.

6. **Diferencia con el modelo del repo.** El modelo del repo es una constante Ω_Λ0 ≈ 0.685 más un ruido
   estacionario pequeño X de media cero. Everpresent Λ no tiene constante: Λ es solo ruido de media
   cero, cuya amplitud tiene que ser del orden de Ω_Λ0 para dar la aceleración actual. DESI DR2 admite
   el primero con σ_X < 2.5 × 10⁻². El segundo, en las dos implementaciones publicadas, solo deja
   pasar realizaciones raras, y todas necesitan un ω_m incompatible con el CMB. Valen las limitaciones
   de la sección 4.
