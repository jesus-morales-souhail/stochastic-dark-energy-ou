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

Validación del pipeline frente al paper de DESI DR2 (arXiv:2503.14738, fuente TeX):

| | Ω_m | h·r_d [Mpc] | χ²/dof |
|---|---|---|---|
| DESI DR2, ΛCDM con solo DESI: Sec. 6, ec. `constraints_lcdm_desi`, y tabla `tab:cosmo_constraints`, fila "DESI" bajo ΛCDM (Ω_m) | 0.2975 ± 0.0086 | 101.54 ± 0.73 | 10.2/(13 − 2) (Sec. 3, "DESI Distance Measurements") |
| este pipeline (mejor ajuste, `robustez.json` → `lcdm_desi`) | 0.2971 | 101.55 | 10.28/11 |

DESI da medias posteriores y nosotros el mejor ajuste; las diferencias son 0.05σ en Ω_m y 0.01σ en
h·r_d.

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

### 3e. r_d por realización

Hasta aquí r_d era fijo (3a–3c) o salía de una fórmula con ω_m (3d), suponiendo que no hay energía
oscura antes de la recombinación. Aquí se calcula con la historia de cada realización. Scripts
`rdcalc.py` y `rd_real.py`; resultados en `resultados/rd_real.json` y
`rd_real_m2.npz` / `rd_real_m1_alpha*.npz`.

**Montaje.**
- r_d = ∫_{z_d}^∞ c_s/H dz, con c_s = c/√(3(1 + R)), R = 3ω_b/(4ω_γ(1+z)), ω_b = 0.02237.
- z_d = 1059.93 (CAMB, `rd_planck.json`), escalado con ω_m por la fórmula de Eisenstein–Hu; no se
  recalcula con la historia de cada realización.
- Validación: con Ω_de ≡ 0, la integral sin normalizar da 146.936 Mpc, −0.11 % frente a CAMB (147.103)
  y −0.087 % frente a Aubourg (147.063) (`rdcalc.py`). Todo r_d se normaliza a CAMB con ese único
  factor (1.00114), y la cuadratura de `rd_real.py` reproduce CAMB a 4 × 10⁻⁷ con Ω_de ≡ 0
  (`model2.check_rd_no_de`).
- **Modelo 2:** el proceso se extiende a λ ∈ [−ln 10⁵, −ln 1000] con paso 0.1, más los 10 puntos de
  la ventana: 57 puntos nuevos. Se muestrean condicionados a Ω̂(λ_rec); los 401 puntos anteriores
  conservan sus valores, y la correlación con el resto es < 3.2 × 10⁻⁷. La integral llega a z = 10⁵,
  y más allá se usa el integrando estándar (2.65 Mpc de ∼147).
- **Modelo 1:** r_d se acumula dentro de cada run para 25 valores de ω_m (escalado de [D1] Sec. 4.4,
  radiación sin reescalar), desde z_d hasta z = 9999, donde arranca la simulación. Para z > 9999 se
  usa el integrando estándar (Λ = 0, coherente con S_0 = 0 al arranque): 24.4 Mpc de 146.9, el 16.6 %.
  Comprobación con α = 0: la r_d acumulada más la cola coincide con la integral estándar a
  1.5 × 10⁻⁵.
- **Proxy de CMB en ventana:** |Ω_de(z)| < 0.25 en todo z ∈ [1000, 5000]. En el Modelo 2, en 10 puntos
  log-espaciados; en el Modelo 1, en los puntos de la malla de 25 en [1000, 10⁴] que caen en
  [1000, 5000].
- **Criterios:** (a) r_d propio, ω_m fijo; (b) r_d propio, ω_m libre en [0.10, 0.20]; (c) r_d propio,
  proxy en ventana y ω_m a ±1.6 % de Planck (±2σ). **(c) es el criterio final.** Con 0 éxitos, el
  límite al 95 % es 3/N.

**Modelo 2 (20000 realizaciones, α = 0.8824, μ = 0.9804)** (`rd_real.json` → `model2`):

- r_d propio, p16/p50/p84 = 96 / 143 / 179 Mpc. Con α = 0.88 la energía oscura no es despreciable
  antes de la recombinación.
- Proxy: el de un solo punto (z = 1100) lo pasa el 22.8 %; el de ventana, el **2.2 % (444/20000)**.

| | A (≤ 29.84) | B (Δχ² < 4) | A' (≤ 10.28) | mejor χ² |
|---|---|---|---|---|
| (a) ω_m fijo | 1 (5 × 10⁻⁵) | 2 (1.0 × 10⁻⁴) | 0 | 10.6 |
| (a) + ventana | 0 | 0 | 0 | 254 |
| (b) ω_m libre | 37 (1.85 × 10⁻³) | 46 (2.3 × 10⁻³) | 4 (2 × 10⁻⁴) | 6.8 (ω_m = 0.186, falla la ventana) |
| (b) + ventana | 1 (5 × 10⁻⁵) | 1 (5 × 10⁻⁵) | 0 | 16.5 |
| (c) ventana y ω_m ±1.6 % | **0 (< 1.5 × 10⁻⁴)** | **0** | **0** | **230.7** (ω_m = 0.1453, r_d = 156.5) |
| (c) sin ventana (referencia) | 6 (3 × 10⁻⁴) | 6 (3 × 10⁻⁴) | 1 | 9.0 |

**Modelo 1 (run_das, n = 10⁵; 300000 + 200000 semillas)** (`rd_real.json` → `model1_alpha*`).
Fracciones sobre semillas. La resimulación con las mismas semillas da json idénticos salvo el tiempo
de cálculo.

| | α = 0.01 | α = 0.005 |
|---|---|---|
| r_d propio a ω_m0, p16/p50/p84 | 132 / 149 / 165 Mpc | 139 / 148 / 156 Mpc |
| pasan el proxy en ventana (17 puntos en [1000, 5000]) | **1 de 82795 completados** | **2033 de 198858 (1.0 %)** |
| (a) ω_m fijo: A / B / A' | 1 / 1 / 0 (χ² = 18.1) | 0 / 0 / 0 (mejor 51.9) |
| (a) + ventana | 0 / 0 / 0 | 0 / 0 / 0 (mejor 93.9) |
| (b) ω_m libre: A / B / A' | 2 / 2 / 0 (mejor 17.76, ω_m = 0.142, H0 = 58.5, falla la ventana) | 1 / 1 / 0 (mejor 16.27, ω_m = 0.134, H0 = 54.8, falla la ventana) |
| (b) + ventana | 0 / 0 / 0 | 0 / 0 / 0 (mejor 88.8) |
| **(c) ventana y ω_m ±1.6 %** | **0 / 0 / 0 (< 1.0 × 10⁻⁵)** | **0 / 0 / 0 (< 1.5 × 10⁻⁵); mejor χ² = 89.3, H0 = 51.6, Ω_Λ(0) = 0.47** |
| (c) sin ventana (referencia) | 1 / 1 / 0 (χ² = 17.76) | 0 / 0 / 0 (mejor 44.0) |

Notas:
- Con α = 0.01, en el 49 % de los runs completados algún ω_m de la malla da H² < 0 al reescalar Λ con
  ω_m sin reescalar la radiación (`frac_rd_nan_any_wm`); esos ω_m quedan excluidos de la minimización.
- El único run con α = 0.01 que pasa la ventana tiene Ω_Λ(0) ≈ 0.02 y H0 ≈ 38: un universo casi sin
  energía oscura, con χ² = 11889.
- Con el r_d propio, el mejor run con α = 0.01 pasa A con ω_m dentro de ±1.6 % (χ² = 17.76), cosa que
  con r_d de fórmula requería ω_m = 0.108. Pero ese run falla el proxy (Ω_de(1100) = 0.58).

**Resultado final (criterio c):** con r_d calculado con la historia de cada realización, energía
oscura pequeña en toda la ventana de recombinación y ω_m a ±2σ de Planck, **ninguna** realización de
ninguno de los dos modelos ajusta DESI DR2 tan bien como ΛCDM-Planck, ni con Δχ² < 4:

- Modelo 2: 0/20000, < 1.5 × 10⁻⁴.
- Modelo 1: 0/300000 (α = 0.01), < 1.0 × 10⁻⁵; y 0/200000 (α = 0.005), < 1.5 × 10⁻⁵.

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
| **Modelo 2, criterio final (c)** | DESI DR2 BAO + ventana CMB, r_d propio, ω_m ±1.6 % | 20000 | **0 (< 1.5 × 10⁻⁴)** | **0** | `rd_real.json` |
| **Modelo 1, α = 0.01, criterio final (c)** | ídem | 300000 | **0 (< 1.0 × 10⁻⁵)** | **0** | `rd_real.json` |
| **Modelo 1, α = 0.005, criterio final (c)** | ídem | 200000 | **0 (< 1.5 × 10⁻⁵)** | **0** | `rd_real.json` |

Diferencia de método: [D2] optimiza Ω_m h², α y M_B por semilla. En 3a–3c los parámetros están fijos
(Planck para ω_b, ω_c y r_d; α del caso). En 3d se optimiza ω_m por realización (con r_d(ω_m)), pero
no α.

---

## 4. Limitaciones

- **r_d en el Modelo 1 por encima de z = 10⁴.** La simulación arranca en z ≈ 10⁴; para z > 9999 se
  usa el integrando estándar sin energía oscura, que es el 16.6 % de r_d. Si el modelo tuviera Λ de
  orden uno también a z > 10⁴, ese tramo cambiaría. En el Modelo 2 el tramo análogo es z > 10⁵, el
  1.8 %.
- **z_d no se recalcula por realización.** Es el de CAMB escalado con ω_m (Eisenstein–Hu); la energía
  oscura cerca de z_d lo movería ligeramente.
- **CMB solo por proxy.** Se usa |Ω_de| < 0.25 (en z = 1100 en 3a–3d, en la ventana z ∈ [1000, 5000]
  en 3e) en vez de la verosimilitud de Planck. [D2] Sec. 5 sí usa la verosimilitud: su mejor χ² es ∼157 peor que ΛCDM, y ∼20 peor con la
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
   - Criterio final (c): r_d calculado con la historia de cada run, |Ω_de| < 0.25 en todo
     z ∈ [1000, 5000] y ω_m a ±1.6 % de Planck. Pasan **0 de 300000** semillas con α = 0.01
     (< 1.0 × 10⁻⁵ al 95 %) y **0 de 200000** con α = 0.005 (< 1.5 × 10⁻⁵), con los tres criterios
     A, B y A'.
   - El mejor run en ventana tiene χ² = 89.3 con H0 = 51.6 (α = 0.005).
   - Sin el proxy, o con ω_m libre en ±40 %, pasan 1–2 runs, todos con H0 = 51–62 y energía oscura de
     orden uno en la recombinación (`rd_real.json`).

3. **Modelo 2 frente a DESI DR2.**
   - Criterio final (c): **0 de 20000** (< 1.5 × 10⁻⁴), y el mejor χ² en ventana con ω_m a ±1.6 % es
     230.7 (`rd_real.json`).
   - El proxy en ventana lo pasa solo el 2.2 % de las realizaciones, frente al 22.8 % del proxy en un
     solo punto.
   - Con ω_m libre y sin proxy, el r_d propio deja pasar 37/20000 (1.85 × 10⁻³) en A; con ω_m a
     ±1.6 % y sin proxy, 6/20000. Todas tienen energía oscura de orden uno en la recombinación.

4. **BAO frente a SN.** DESI DR2 BAO solo, con ω_m libre, deja pasar una fracción del orden de la de
   SN en [D2]: 1.25–1.85 × 10⁻³ en el Modelo 2, frente a 1.8 × 10⁻⁴. En SN, M_B absorbe la
   normalización y H0 queda libre. En BAO, D/r_d mide H0·r_d con r_d fijado por la física
   pre-recombinación, pero esa restricción se puede esquivar moviendo ω_m o dejando que la energía
   oscura cambie r_d. **Lo que excluye el modelo es BAO junto con lo que fija el CMB:** ω_m a ±2σ de
   Planck y energía oscura pequeña alrededor de la recombinación. Con esas dos condiciones, 0
   realizaciones en los dos modelos.

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
