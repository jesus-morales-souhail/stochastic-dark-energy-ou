# Sustituir el proxy del CMB por CAMB + Planck 2018: exploración

Exploración, sin resultados. No cambia ningún número de `EVERPRESENT_vs_DESI.md` ni de ningún paper.
Pregunta: ¿se puede sustituir el proxy |Ω_de| < 0.25 (en z = 1100, o en la ventana z ∈ [1000, 5000])
por un cálculo con CAMB, solo para las realizaciones que sobreviven a DESI DR2 BAO con algún criterio
de 3a, 3d o 3e?

Ficheros de esta carpeta (todos de exploración):

| script | qué hace | salida |
|---|---|---|
| `supervivientes.py` | lista de supervivientes, signo de Ω_de en los puntos guardados | `supervivientes.{json,txt}` |
| `resim_m1.py` | resimula un run del Modelo 1 desde (α, bloque, índice) y lo compara con el npz | `resim_m1_alpha*.txt`, `historias_m1/*.npz` |
| `prueba_camb.py` (venv de `analisis_cmb`) | pasa dos supervivientes del Modelo 2 a CAMB con `set_w_a_table` | `prueba_camb.{json,txt}` |
| `tiempos.py` (venv, `OMP_NUM_THREADS=2`) | tiempos con la configuración de `realizaciones.py` | `tiempos.{json,txt}` |

Resumen:
- Hay **65 supervivientes distintas del Modelo 2 y 6 del Modelo 1** (contando las 3 de los ±1 del
  umbral de 3d).
- El Modelo 1 se puede resimular run a run, y la resimulación es idéntica bit a bit.
- El Modelo 2 **no tiene historia completa**: z ∈ (3, 1000) no se muestreó nunca, y hay que sortearlo
  de nuevo.
- `set_w_a_table` funciona si ρ_DE no cambia de signo. Si cambia, CAMB **corre sin avisar con un fondo
  equivocado**, o falla en Recfast.
- Solo 1 de las 71 supervivientes tiene Ω_de > 0 en todos los puntos guardados, y su tramo sin muestrear
  es positivo en el 0.2 % de los sorteos. La opción "solo ρ_DE > 0" deja la muestra vacía.
- Zwane et al. sí trataron historias con cambio de signo: su mejor ajuste cruza cero varias veces. No
  dicen cómo lo hicieron.

---

## 1. Supervivientes

Criterios (`EVERPRESENT_vs_DESI.md`, sección 3):
- A: χ² ≤ 29.84.
- B: Δχ² < 4.
- 3a: ω_m y r_d fijos.
- 3d: ω_m libre, r_d de fórmula.
- 3e(a): r_d propio, ω_m fijo.
- 3e(b): r_d propio, ω_m libre.
- 3e(c)\*: r_d propio, ω_m a ±1.6 %, **sin** ventana. Es la fila "(c) sin ventana" de 3e; con
  ventana no pasa ninguna.

Los recuentos de `supervivientes.json` reproducen los del documento:

| | A | B |
|---|---|---|
| 3a | 1 | 1 |
| 3d | 25 | 29 |
| 3e(a) | 1 | 2 |
| 3e(b) | 37 | 46 |
| 3e(c)\* | 6 | 6 |

La unión son **65 realizaciones** distintas.

### Modelo 2 (α = 0.8824, μ = 0.9804; 20000 realizaciones, semilla 20261016)

El "índice" es la columna j de `rng.standard_normal((401, 20000))` con semilla 20261016
(`bao_sim.main_m2`). Entre paréntesis va el ω_m óptimo. Ω_de en los puntos guardados: 400 puntos en
z ∈ [0, 3], z = 1100 y 57 puntos en z ∈ [10³, 10⁵] (`rd_real.py`, semilla 20261017). "Ventana": pasa
|Ω_de| < 0.25 en z ∈ [1000, 5000].

| índice | pasa (χ², ω_m) | Ω_de(0) | Ω_de(1100) | mín Ω_de, z ∈ [0, 3] | mín Ω_de, z ∈ [10³, 10⁵] | ventana |
|---|---|---|---|---|---|---|
| 286 | 3e(b) B 31.5 (0.106); 3d B 31.5 (0.168) | +0.58 | +0.16 | -0.14 | -0.54 | no |
| 336 | 3d B 31.0 (0.151) | +0.68 | +0.18 | +0.09 | -0.99 | no |
| 1431 | 3e(b) A 24.9 (0.119) | +0.52 | -0.45 | -0.19 | -0.77 | no |
| 1540 | 3e(b) A 15.9 (0.135); 3e(c)\* A 27.9 (0.141) | +0.51 | -0.96 | -0.47 | -0.97 | no |
| 1860 | 3d A 11.2 (0.125) | +0.68 | -0.44 | +0.11 | -0.64 | no |
| 1864 | 3e(b) A 23.1 (0.200) | +0.51 | -0.94 | -0.64 | -0.96 | no |
| 2114 | 3d A 22.4 (0.105) | +0.70 | +0.02 | +0.13 | -0.79 | no |
| 2180 | 3e(b) A 10.9 (0.100); 3d A 16.8 (0.200) | +0.57 | +0.55 | -0.12 | -0.96 | no |
| 2510 | 3e(b) A 9.8 (0.176) | +0.39 | -0.87 | -0.83 | -0.91 | no |
| 3249 | 3e(b) A 22.9 (0.175) | +0.44 | -0.81 | -0.80 | -0.93 | no |
| 4369 | 3e(b) B 33.1 (0.162) | +0.79 | -0.02 | +0.31 | -0.72 | no |
| 5329 | 3e(b) B 30.4 (0.175) | +0.41 | -0.45 | -0.71 | -0.98 | no |
| 5331 | 3e(b) B 33.6 (0.100) | +0.86 | +0.40 | +0.62 | -0.65 | no |
| **5492** | 3d A 20.7 (0.112) | +0.69 | +0.08 | **+0.08** | **+0.07** | no |
| 5508 | 3d A 27.0 (0.200) | +0.57 | -0.77 | -0.20 | -0.88 | no |
| 5598 | 3d A 12.0 (0.119) | +0.64 | -0.43 | +0.14 | -0.93 | no |
| 6078 | 3e(b) A 25.6 (0.133); 3d A 25.5 (0.128) | +0.66 | +0.11 | -0.00 | -0.81 | no |
| 6126 | 3e(b) A 16.2 (0.109) | +0.61 | -0.35 | -0.15 | -0.68 | no |
| 6233 | 3d A 26.7 (0.163) | +0.64 | +0.22 | -0.10 | -0.84 | no |
| 6342 | 3e(b) B 33.8 (0.200) | +0.62 | +0.45 | -0.20 | -0.87 | no |
| 6455 | 3e(a) B 30.4; 3e(b) A 24.5 (0.147); 3e(c)\* A 25.5 (0.145) | +0.54 | -0.14 | -0.33 | -0.93 | no |
| 6630 | 3e(b) A 26.6 (0.114); 3d A 26.6 (0.152) | +0.72 | +0.29 | +0.03 | -0.91 | no |
| 6701 | 3d A 29.8 (0.199) | +0.58 | +0.86 | -0.29 | -0.67 | no |
| 7339 | 3e(b) B 33.8 (0.109) | +0.52 | -0.69 | -0.28 | -0.92 | no |
| 7525 | 3d A 20.5 (0.162) | +0.64 | +0.11 | +0.05 | -0.45 | no |
| 7907 | 3e(b) A 18.1 (0.118) | +0.42 | -0.82 | -0.51 | -0.96 | no |
| 8163 | 3e(b) A 19.2 (0.138); 3e(c)\* A 24.7 (0.141) | +0.79 | +0.85 | +0.37 | -0.83 | no |
| 8529 | 3e(b) A 21.9 (0.135); 3d A 21.9 (0.190) | +0.56 | -0.32 | -0.11 | -0.88 | no |
| 9421 | 3e(b) A 26.0 (0.161) | +0.45 | -0.76 | -0.59 | -0.83 | no |
| 9444 | 3e(b) B 29.9 (0.200) | +0.42 | +0.02 | -0.69 | -0.97 | no |
| 9525 | 3d B 33.3 (0.186) | +0.65 | -0.78 | -0.00 | -0.93 | no |
| 9526 | 3e(b) A 7.9 (0.174) | +0.37 | -0.85 | -0.86 | -0.97 | no |
| 10569 | 3e(b) A 6.8 (0.186) | +0.50 | -0.69 | -0.35 | -0.70 | no |
| 10788 | 3d A 20.3 (0.151) | +0.68 | -0.91 | -0.06 | -0.97 | no |
| 10885 | 3e(b) A 28.5 (0.111) | +0.53 | -0.81 | -0.31 | -0.90 | no |
| 11072 | 3e(b) A 11.3 (0.152); 3e(c)\* A 23.9 (0.145) | +0.53 | -0.87 | -0.26 | -0.87 | no |
| **11555** | 3e(b) A 16.5 (0.172) | +0.56 | -0.20 | **-0.23** | -0.71 | **sí** |
| 11561 | 3e(b) A 29.7 (0.200) | +0.42 | -0.77 | -0.77 | -0.87 | no |
| 12423 | 3e(b) A 25.4 (0.198) | +0.43 | -0.88 | -0.90 | -0.89 | no |
| 12909 | 3e(b) A 22.6 (0.149); 3e(c)\* A 27.2 (0.145) | +0.50 | -0.26 | -0.51 | -0.79 | no |
| 13052 | 3d A 9.7 (0.127) | +0.66 | +0.19 | +0.11 | -0.93 | no |
| 13826 | 3d A 13.0 (0.152) | +0.67 | +0.28 | +0.04 | -0.95 | no |
| 14006 | 3e(b) A 22.5 (0.166); 3d A 22.5 (0.200) | +0.57 | -0.89 | -0.19 | -0.91 | no |
| 14159 | 3e(b) A 13.3 (0.160) | +0.51 | -0.33 | -0.35 | -0.91 | no |
| 15180 | 3e(b) A 22.6 (0.107) | +0.44 | -0.91 | -0.45 | -0.91 | no |
| 15208 | 3d A 7.9 (0.197) | +0.59 | +0.92 | -0.13 | -0.70 | no |
| 15266 | 3e(b) A 19.0 (0.185) | +0.43 | -0.67 | -0.81 | -0.96 | no |
| 15478 | 3d A 25.6 (0.160) | +0.65 | -0.59 | -0.13 | -0.65 | no |
| 15606 | 3d A 24.5 (0.117) | +0.67 | +0.93 | +0.15 | -0.62 | no |
| 16309 | 3e(b) A 14.6 (0.175) | +0.47 | -0.18 | -0.45 | -0.95 | no |
| 16616 | 3e(a) A 10.6; 3e(b) A 9.0 (0.141); 3e(c)\* A 9.0 (0.141) | +0.48 | -0.77 | -0.48 | -0.98 | no |
| 16813 | 3d A 19.4 (0.178) | +0.67 | +0.42 | -0.04 | -0.47 | no |
| 17070 | 3e(b) A 12.1 (0.186) | +0.52 | -0.75 | -0.58 | -0.77 | no |
| 17132 | 3e(b) A 24.5 (0.132) | +0.60 | -0.84 | -0.10 | -0.87 | no |
| 17148 | 3e(b) B 33.1 (0.100) | +0.78 | +0.39 | +0.32 | -0.39 | no |
| 17158 | 3e(b) A 28.7 (0.117) | +0.81 | +0.83 | +0.46 | -0.68 | no |
| 17299 | 3e(b) A 13.1 (0.100); 3d A 13.1 (0.193) | +0.65 | -0.04 | -0.09 | -0.79 | no |
| 17316 | 3e(b) A 29.8 (0.128); 3d B 29.9 (0.177) | +0.62 | +0.69 | -0.17 | -0.76 | no |
| 17497 | 3d A 25.7 (0.191) | +0.65 | -0.89 | -0.21 | -0.95 | no |
| 17590 | 3e(b) A 24.0 (0.160) | +0.44 | -0.50 | -0.42 | -0.88 | no |
| 17853 | 3e(b) A 10.5 (0.123); 3d A 10.5 (0.195) | +0.63 | -0.28 | -0.18 | -0.52 | no |
| 17886 | 3d A 10.7 (0.182) | +0.57 | +0.34 | -0.10 | -0.58 | no |
| 17907 | 3e(b) B 32.5 (0.187) | +0.78 | -0.26 | +0.17 | -0.55 | no |
| 18324 | 3a A 22.6; 3e(b) A 11.6 (0.132); 3d A 11.6 (0.148) | +0.70 | +0.54 | -0.01 | -0.96 | no |
| 19827 | 3e(b) A 28.2 (0.100) | +0.82 | -0.03 | +0.49 | -0.83 | no |

- Las 65 tienen Ω_de(0) > 0.
- 64 tienen Ω_de < 0 en algún punto guardado. La excepción es 5492.
- Las 9 que en 3d pasan además el proxy de un punto (286, 336, 2114, 5492, 6078, 6233, 7525, 13052,
  17299) cambian de signo en z ∈ [10³, 10⁵], salvo 5492.

### Modelo 1 (run_das, n = 10⁵)

El "bloque" c corresponde a la semilla 20261016 + c. El "índice" es la columna dentro del bloque de
20000. "Fila" es la posición en `bao_m1_alpha*.npz`. H0 y Ω_Λ(0) van a ω_m de Planck. ±1: en 3d, a
menos de 1 del umbral B (el χ² del Modelo 1 en 3d es solo reescalado).

| α | fila | bloque (semilla) | índice | pasa (χ², ω_m) | H0 | Ω_Λ(0) | Ω_de(1100) | Ω_de mín / máx, z ∈ [10³, 10⁴] |
|---|---|---|---|---|---|---|---|---|
| 0.01 | 5367 | 0 (20261016) | 19735 | 3d A 17.76 (0.108); 3e(a) A 18.10; 3e(b) A 17.76 (0.142); 3e(c)\* A 17.76 (0.142) | 58.7 | 0.58 | +0.58 | -1.82 / +0.65 |
| 0.01 | 29757 | 5 (20261021) | 7692 | 3e(b) A 27.32 (0.163) | 71.0 | 0.72 | +0.54 | -2.68 / +0.72 |
| 0.01 | 63491 | 11 (20261027) | 9957 | 3d B 33.25 (0.102) | 63.6 | 0.65 | -0.21 | -1.08 / +0.61 |
| 0.01 | 64524 | 11 (20261027) | 13734 | 3d ±1 34.40 (0.137) | 52.2 | 0.48 | +0.06 | -3.61 / +0.43 |
| 0.005 | 51223 | 2 (20261018) | 11542 | 3d A 16.27 (0.174); 3e(b) A 16.27 (0.134) | 56.6 | 0.55 | -0.09 | -1.66 / +0.17 |
| 0.005 | 196680 | 9 (20261025) | 17809 | 3d ±1 34.19 (0.197) | 57.0 | 0.56 | +0.19 | -0.12 / +0.39 |

Las 6 tienen Λ < 0 en algún z ∈ [10³, 10⁴]. En el Modelo 1, Ω_de = Λ/(3H²) baja de −1 cuando Λ < 0
domina sobre la materia. Ninguna mantiene ρ_DE > 0.

### ¿Es recuperable la historia completa?

**Modelo 2. No de los npz:**
- `bao_m2_bestfit.npz` guarda χ², H0 y Ω_de(z_rec).
- `rd_real_m2.npz` guarda r_d y χ²; no guarda historias.

**Pero los puntos que se usaron se regeneran exactamente** desde las semillas:
- Los 401 puntos: semilla 20261016 y Cholesky con jitter 10⁻¹⁰. `robustez.json` →
  `m2_check_regenerated_vs_saved_max_abs` = 3.3 × 10⁻¹⁰ en χ².
- Los 57 puntos de z ∈ [10³, 10⁵]: semilla 20261017.

**El problema es lo que nunca se muestreó:**
- **z ∈ (3, 1000):** BAO solo necesita z ≤ 2.33 y r_d solo z ≥ z_d, así que ese tramo no existe en
  ninguna realización. Para el CMB sí hace falta: distancia a la última dispersión, ISW y lensing.
- **z > 10⁵.**

El proceso gaussiano permite **sortearlos condicionados a los puntos guardados**
(`prueba_camb.gap_sampler`). La media condicional empalma con los puntos guardados a 10⁻³ (`prueba_camb.json`
→ `gap_edge_check`). Pero es aleatoriedad nueva, no la de la semilla original. Para cada superviviente,
el resultado del CMB es una **distribución sobre sorteos del hueco**, no un número.

Lo correcto sería dar, por superviviente, P(CMB aceptable | puntos guardados) con N sorteos. La correlación
del núcleo cae a e⁻¹ en Δλ ≈ 1.4, y el hueco mide 5.5 e-folds, así que el centro del hueco queda casi
libre.

**Modelo 1. Hay que resimular, y se puede:**
- El npz guarda Λ solo en los 7 z de DESI, en z_rec y en 25 puntos de z ∈ [10³, 10⁴].
- `run_das` saca en cada paso `rng.standard_normal(20000)` para todo el bloque, y las columnas
  evolucionan por separado: el número de subpasos RK4 depende solo de la malla en a.
- `resim_m1.py` corre `run_das` con `n_real = 1` y un generador que, en cada paso, saca los 20000
  normales del bloque y devuelve la columna del run. Con `z_rec = 9999` guarda u(a) entero.

Resultado, en los 6 runs: χ², H0, Ω_Λ(0), Ω_de(z_rec), `win_ode` y `rd_part` **idénticos a los del
npz (diferencia 0.00)**. Cada run tarda 22–25 s. Las historias quedan en `historias_m1/` (5.8 MB).

La historia del Modelo 1 es completa: Λ = 0 para z > 9999, porque S_0 = 0 al arranque.

---

## 2. Cómo lo hicieron Zwane, Afshordi y Sorkin (arXiv:1703.06265)

**Fuente.** TeX descargado de arXiv el 2026-10-31 (`https://arxiv.org/e-print/1703.06265`, v4,
tarball del 2018-05-24). El `fcc_3.tex` del tarball es **idéntico byte a byte** a
`../fuentes/zwane2018_1703.06265.tex`. Las figuras se leyeron del mismo tarball. Numeración de la
versión compilada (estilo JHEP): secciones 1–8, ecuaciones (sección.n).

**Código: CosmoMC + CAMB modificado** (Sec. 3, primer párrafo):
> "we use CosmoMc together with CAMB, some parts of which we had to rewrite in order to accommodate a
> stochastic dark energy."

**Parámetros:**
- Muestrean (α, μ, seed) como parámetros más de la cadena MCMC. Sobre la semilla: "the likelihood
  surface will not be a smooth function of seed" (Sec. 3, segundo párrafo).
- Generan Ω̂ por Fourier, con espectro α² e^{−ω²/4μ}√(π/μ) (Sec. 2.2, ec. 2.11). No dicen qué
  frecuencias ni qué rango de λ usan.
- La Fig. 3 usa "1000 ω's sampled".

**Fondo y perturbaciones: los dos** (Sec. 2.3, "Inhomogeneities"):
> "we shall simply adopt the model provided by the CAMB software package, which assumes that dark
> energy is a perfect fluid with perturbations that propagate with the speed of light […] equivalent to
> a quintessence dark energy […] Then, fixing the history of the homogeneous dark energy […] and
> assuming also energy-momentum conservation and adiabatic initial conditions, uniquely determines the
> behavior of linear perturbations."

Es decir, un fluido con c_s² = 1 y perturbaciones de la energía oscura activadas. En la Sec. 8
proponen como prueba de robustez "turn off the dark-energy perturbations", que no hacen.

**Cómo entra la historia en CAMB.** Hay un apéndice **comentado** en el TeX, que no aparece en el
PDF: "CAMB Cook Book", líneas 1029–1080. Dice:
- Ω_DE(a) → ρ_DE(a) por la ecuación de Friedmann → w_DE(a) = p/ρ.
- "w_DE(a) is then feed to CAMB as w_ppf and a_ppf in equation-ppf.f90". Es decir, la tabla w(a) del
  módulo PPF, lo mismo que `set_w_a_table` hoy.
- Además modificaron `halofit_ppf.f90` (funciones `omega_m`, `omega_v` con w) y `derivsv` en
  `equations_ppf.f90` (`grhov_t=grhov*a**(-1-3*w_de)`).

Observaciones:
- El texto habla de "perfect fluid", pero el cook book usa PPF.
- La línea de `derivsv` es la ley de w constante aplicada con la w local, no ρ_0 exp(−3∫(1+w) d ln a).
- Esto se deduce del TeX comentado; no es una descripción publicada.

**H0 y la normalización:**
- En la Tabla 1, θ_MC es un parámetro de la cadena (100θ_MC = 1.041 ± 0.00063), y H0 = 68.48 ± 0.67 y
  Ω_Λ = 0.6781 ± 0.0081 salen derivados con barras de error.
- Eso indica (es una inferencia) que Ω_de(0) no se fijaba a tanh Ω̂(0), sino por planitud a partir de
  θ_MC. Es decir: se usaba la forma de la historia vía w(a), y CAMB normalizaba ρ_DE hoy.

**Λ negativa: el texto no dice cómo la trataron.** Pero sus historias la tienen:
- **Fig. 1** (`CMBomg2017f.pdf`): la historia del Modelo 2 que ajusta Planck + BAO (línea roja
  discontinua) cruza cero cerca de z ≈ 1 y varias veces más hasta z = 10⁴. Leído a ojo de la figura,
  es negativa aproximadamente en z ≈ 1–12, ≈ 60–200 y ≈ 1000–10⁴, con un mínimo ≈ −0.65 cerca de
  z ≈ 5000.
- **Fig. 6** (`firstBossR3.pdf`, eje "ρ_de/ρ_c0"): ρ_de del mejor ajuste es negativa para z ≳ 1.85.
- **Sec. 5**: "there is a 78% chance that Ω_de is negative somewhere in the range 1.5 < z < 3".
- Una tabla w(a) no puede representar un cruce de ρ_DE por cero: w = p/ρ diverge ahí.
- Que su mejor ajuste cruce cero varias veces implica que el código "reescrito" trataba ρ_DE con
  signo, o que hicieron algo con w en los cruces. El TeX no lo documenta.
- La nota al pie de la Sec. 2.3 dice que la compresión tanh "rules out any such 'phantom' behavior".
  No es así para w: en la tabla de la superviviente 5492 de aquí, w_mín = −1.33 (`prueba_camb.json`).

**Likelihood: Planck 2015** (Tablas 1 y 2):
- `plik` (alto ℓ; el pie de la Tabla 1 habla de "temperature and polarization spectra at high
  multipoles"), `lowTEB`, lensing y BKPLANCK.
- Más BAO: 6dF, MGS, DR11 CMASS, DR11 LOWZ, DR11 Lyα auto y cross.
- χ² total: 11335.2 frente a 11334.6 de ΛCDM, con α = 0.8824 y μ = 0.9804.
- El proxy |Ω_de| < 0.25 en z ≈ 1000 es una descripción a posteriori de los buenos ajustes (Sec. 3,
  último párrafo: "∼23 %"), no el criterio que usaron.

---

## 3. Prueba técnica con CAMB (venv de `analisis_cmb`, CAMB 2.0.4)

`prueba_camb.py` hace lo siguiente:
- Ω̂ en los 458 puntos guardados, más el hueco z ∈ (3.09, 979): 275 puntos con paso 0.02 en λ,
  sorteados condicionados (semilla nueva 20261031).
- Spline cúbico en λ hasta una malla uniforme de 4000 puntos en ln a ∈ [ln 10⁻⁵, 0].
- ρ_DE ∝ Ω_de/(1 − Ω_de)(ω_m a⁻³ + ω_r a⁻⁴) y w = −1 − ⅓ d ln|ρ_DE|/d ln a.
- CAMB con ω_b = 0.02237, ω_c = ω_m − ω_b − ω_ν, m_ν = 0.06, plano, H0 = H(z = 0) de la realización
  y `DarkEnergyPPF.set_w_a_table`.

### Caso 1: 5492, Ω_de > 0 en todos los puntos guardados (3d A, ω_m = 0.1125, H0 = 59.82)

- **Hueco:** solo 1 de 500 sorteos (0.2 %) es positivo en todo z ∈ (3, 1000). Se usa ese (el sorteo
  195). En los demás, ρ_DE cambia de signo.
- **Funciona:**
  - Sin errores ni avisos de Python.
  - H(z) de CAMB / H(z) de la realización = 1 − (1.6 × 10⁻⁵ … 1.8 × 10⁻³). El máximo está en
    z ≈ 1000; probablemente viene del neutrino masivo de CAMB frente a ω_r sin masa en `bao_sim`
    (no comprobado).
  - Ω_DE(0) de CAMB = 0.68565, frente a 0.68563 de la realización.
  - TT ℓ = 2–29 y C_L^φφ finitos.
- **Avisos de Fortran** por stdout: `Integrate_Romberg failed to converge` (error 3–4 × 10⁻⁷ con
  tolerancia 10⁻⁷).
- **100θ\* = 0.6845** a su propio H0, frente a 1.0412 de ΛCDM-Planck. Es un valor de prueba, no un
  resultado: depende del único sorteo positivo del hueco y de la extrapolación de CAMB para z > 10⁵.
  Aun así indica que θ\* es, con diferencia, la observable que más discrimina.
- **Con la configuración de `realizaciones.py`**, que resuelve H0 por θ\*, **falla**:
  - Con el intervalo [30, 120]: "Universe stops expanding before today (recollapse not supported)".
    A H0 bajo, Ω_DE(0) = 1 − Ω_m < 0 y la forma de la tabla hace colapsar el universo.
  - Con [40, 150]: 100θ\* va de 0.78 a 0.57 y **decrece** con H0 (en ΛCDM crece), así que nunca
    alcanza 1.0412: "theta_* not bracketed" (`tiempos.json` → `thetastar_scan`). Al subir H0, la
    planitud sube ρ_DE hoy, y con la forma fija de la tabla sube también la energía oscura temprana,
    lo que encoge r_s.

### Caso 2: 11555, ρ_DE cambia de signo (3e(b) A, ω_m = 0.1719, H0 = 62.53; la única que pasa la ventana)

Ω_de ∈ [−0.71, +0.95] con 6 cambios de signo en z < 10⁵ (incluido el sorteo 0 del hueco). Ninguno de
los 500 sorteos del hueco es positivo.

| entrada a CAMB | qué pasa |
|---|---|
| w de ln\|ρ_DE\|: finita, \|w\| hasta 178 junto a los cruces | **Corre sin error ni aviso, y el fondo está mal.** H_CAMB/H_real = 1.21 (z = 2), 1.15 (z = 3), 1.18 (z = 1000–1100), 1.27 (z = 10⁴); 1.0000 donde Ω_de > 0. CAMB integra ρ ∝ exp(−3∫(1+w)), que no puede cambiar de signo: calcula la historia con \|ρ_DE\|. TT y φφ finitos. 100θ\* = 1.231. La corrida completa tarda 78 s (16 hilos), por la rigidez de w |
| w = +∞ en los dos puntos de cada cruce | `CAMBError: Recfast: Rosenbrock integration failed` |
| w = NaN en los mismos puntos | ídem |

El caso más peligroso es el primero: un pipeline ingenuo daría espectros finitos para una cosmología
que no es la de la realización.

### Opciones si hace falta tratar ρ_DE con signo

| opción | qué supone | coste estimado | qué cubre |
|---|---|---|---|
| **(a)** CAMB con fondo modificado en Fortran | Nueva clase de energía oscura que lea una tabla de ρ_DE(a) con signo y de ρ + p (finitos en el cruce), en vez de w. En PPF, dρ/d ln a = −3(ρ + p). Hay que compilar CAMB desde la fuente (gfortran está instalado), envolver la clase en Python y validarla frente a `set_w_a_table` en casos con ρ > 0 (5492). Hay que revisar si el módulo PPF divide por (1 + w) o por ρ_DE en algún sitio | 2–4 días de trabajo + validación. Por realización, el mismo coste que la sección 4 | Todas las supervivientes. Es lo que Zwane et al. tuvieron que hacer de algún modo ("some parts … we had to rewrite") |
| **(b)** CLASS | CLASS 3.3.1 (comprobado en el sdist de `classy`) solo admite `fluid_equation_of_state` = CLP o EDE. Integra ρ_fld con w igual que CAMB, y las perturbaciones del fluido usan δ = δρ/ρ, singular en ρ = 0 salvo con `use_ppf`. Hay que modificar `background.c` (tabla de ρ con signo) y comprobar el PPF de CLASS. Además hay que instalar classy y montar las mismas observables que en `planck_datos.py` | Parecido a (a) o algo mayor, porque el pipeline actual es de CAMB | Todas, con el mismo trabajo de validación |
| **(c)** Solo supervivientes con ρ_DE > 0 en todo z | Ninguna modificación | Trivial | **0–1 de 71.** Con ρ_DE > 0 en todos los puntos guardados: 1 de 65 del Modelo 2 (5492) y 0 de 6 del Modelo 1. El hueco de 5492 es positivo en el 0.2 % de los sorteos. Además, un proceso estacionario en λ cruza cero infinitas veces en z → ∞, así que "en todo z" exige un corte (aquí z = 10⁵). La muestra queda vacía y sesgada: justo la historia más parecida a ΛCDM |
| **(d1)** Solo fondo: θ\* propio (y distance priors) | Calcular r_s(z\*) y D_M(z\*) con el integrador de `rdcalc.py`/`rd_real.py`, que ya admite Ω_de con signo, y comparar con 100θ\* = 1.04110 ± 0.00031 (Planck 2018). Como paso siguiente, los distance priors (R, ℓ_A, ω_b) de Planck 2018 | Horas. Segundos de cálculo para todo, incluidos N sorteos del hueco | Todas. No da TT de bajo ℓ ni lensing, y los distance priors suponen física estándar antes de la recombinación, lo que aquí no se cumple. Mejora el proxy en mucho, sin sustituir la likelihood |
| **(d2)** Quintaesencia con potencial V(φ) reconstruido | Es lo que Zwane et al. dicen que su fluido "equivale" a. ρ < 0 se puede tener con V < 0 y los cruces son suaves. Pero φ̇² = ρ + p ≥ 0 excluye las fases fantasma, y 5492 ya tiene w_mín = −1.33. También exige Fortran (la clase `Quintessence` de CAMB es abstracta en Python) | Como (a), más la reconstrucción | Solo historias sin fase fantasma |

Mi recomendación, si se sigue: **(d1) primero**. Es barata, cubre todas las supervivientes y todos los
sorteos del hueco, y con 0.6845 frente a 1.0412 en el único caso probado es probable que decida casi
todo. Después, **(a)** solo para las que pasen θ\*.

---

## 4. Tiempo por realización y total

Configuración de `realizaciones.py`:
- `OMP_NUM_THREADS = 2` por proceso, 8 procesos.
- Takahashi, lmax 2500 y `lens_potential_accuracy = 1`.
- H0 por bisección en `get_background` con tolerancia 10⁻⁵, que en un intervalo de 60 son unas 23
  iteraciones más 2 evaluaciones de los extremos.
- TT ℓ < 30 y C_L^φφ en bandas.

Medido con `tiempos.py` sobre 5492, en esta máquina (16 núcleos):

| paso | Everpresent (tabla w de 5492) | ΛCDM, mismas llamadas |
|---|---|---|
| `get_background` + θ\* | **3.0 s** (2.96–3.9; una vez 0.43). No depende del tamaño de la tabla: 500, 1000, 2000 y 4000 puntos dan 3.0 s. Coincide con los avisos de Romberg; la causa no está comprobada | 0.003 s |
| bisección completa de H0 | ≈ 25 × 3.0 ≈ **75 s** (en 5492 no se puede hacer: no hay raíz) | 0.07 s |
| corrida completa (TT + lensing) | **6.5–9.6 s** | 1.1 s |

Por realización:
- **≈ 85 s** con H0 resuelto por θ\*.
- **≈ 10–13 s** con H0 fijado por la propia realización. Esa es la elección coherente con el Modelo 2,
  donde H0 no es libre: H0 = 100√((ω_m + ω_r)/(1 − Ω_de(0))).
- En la configuración de `realizaciones.py`, θ\* no entra en el χ², porque se iguala por construcción.
  Para Everpresent habría que meterlo como observable (`planck_datos.py` ya carga su error) en vez de
  usarlo para fijar H0.

Totales para las 71 supervivientes (65 + 6), 8 procesos:

| esquema | Modelo 2 (65) | Modelo 1 (6) | total |
|---|---|---|---|
| H0 de la realización, 1 sorteo del hueco | 65 × 12 s / 8 ≈ 1.6 min | 6 × (24 s de resimulación + 12 s) / 8 ≈ 0.5 min | **≈ 2 min** |
| H0 por θ\* (si hubiera raíz), 1 sorteo | 65 × 85 s / 8 ≈ 12 min | ≈ 1.5 min | **≈ 13 min** |
| H0 de la realización, N = 100 sorteos del hueco | 6500 × 12 s / 8 ≈ 2.7 h | sin hueco: ≈ 0.5 min | **≈ 2.7 h** |
| (d1) solo θ\*, N = 100 sorteos | segundos (numpy, como `rd_real.py`) | segundos | **< 1 min** |

Con `set_w_a_table` tal como está, de estas 71 **solo se puede calcular correctamente 5492, y solo con
1 de cada 500 sorteos del hueco**. Las otras 70 necesitan (a) o (b) antes de que estos tiempos apliquen.
La corrida de 11555 con w rígida tardó 78 s con 16 hilos. Si (a) conserva PPF, las historias con
cruces podrían costar bastante más que 10 s.

---

## Limitaciones de esta exploración

- **El hueco z ∈ (3, 1000) del Modelo 2 se rellena con un sorteo nuevo.** Cualquier número de CAMB
  sobre una superviviente del Modelo 2 depende de ese sorteo. Lo mismo vale para z > 10⁵, donde CAMB
  extrapola la tabla con el w del extremo.
- **El 0.6845 de θ\* es de un solo caso y un solo sorteo.** No es un resultado.
- **Las conclusiones sobre el método de Zwane et al. salen de su TeX y sus figuras.** El tratamiento de
  ρ_DE < 0 y la normalización por planitud son inferencias, no descripciones suyas.
- **Los costes de (a) y (b) son estimaciones**, sin haber abierto el Fortran de CAMB ni el C de CLASS
  más allá de lo citado.
