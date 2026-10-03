# Notas del análisis CMB (OU a varianza máxima DESI)

## Fuente de verdad del ajuste BAO

`results/profile_sigma_x/profile_sigma_x.json` (generado por `scripts/profile_sigma_x_desi.py`,
`cov_mode: full_gaussian_bao_projected_13x13`).

- Mejor ajuste: σ_X → 0 (`sigma_X_at_max = 1e-06`, `max_dlogL_vs_lcdm ≈ -2e-8`).
- Límite 95 % CL (perfil, θ libre): `sigma_X_95CL_upper_profile = 0.025063155506929828`
  (redondeado 2.5e-2). Es el último punto de la malla con ΔlnL ≥ −1.92 (en ese punto
  ΔlnL = −1.8715; el siguiente, 3.155e-2, da −2.0995).

## Ficheros antiguos / no vigentes — NO usar sus números

| Fichero | Valor que NO hay que usar | Por qué |
|---|---|---|
| `results/ou_bao_desi_dr2_run.txt` | θ = 0.001, σ_X = 5e-5, logL(ΛCDM) = 27.013, AIC −50.03 / −48.03 | Ejecución antigua OU+QNM; su logL no coincide con la covarianza completa (24.8996) → casi seguro covarianza diagonal |
| `papers/resume.txt` (sección 1) | θ = 0.001, σ_X = 5e-5, AIC = −50.03 / −48.03 | Copia los parámetros de la ejecución antigua de arriba (el límite 2.5e-2 que cita sí es el vigente) |
| `results/desi_dr2_real_bao/desi_dr2_real_bao.txt` (y `.json`) | σ_X 95 % "approx" = 2.175e-2 | Usa covarianza completa pero es una estimación aproximada, no el perfil oficial |
| `../measurable-stochastic-vacuum/results/r1_lambda_profile/lambda_profile.json` | `lambda_working_from_1p5e-4` | Basado en el antiguo σ_X < 1.5e-4 (objetivo de trabajo, no cota real) |
| `../measurable-stochastic-vacuum/results/r1_lambda_fullcov/lambda_fullcov_profile.json` | `lambda_working_1p5e-4`, `dlogL_at_working_1p5e-4` | Ídem, basado en 1.5e-4 |

Cualquier aparición de σ_X < 1.5e-4 es la cifra antigua (ver `papers/BOUND_CORRECTION_SIGMA_X.md`).

## Elección de (σ_X, θ) para las realizaciones

σ_X = 2.506e-2 (punto exacto de la malla) en los tres casos.

- θ = 0.001 — es el `theta_best` del perfil a esa σ_X, **pero es el borde inferior de la
  malla de θ** del script (1e-3 … 20): DESI no determina θ ahí. Es el "preferido por DESI"
  sólo en ese sentido.
- θ = 1 — exploración (tiempo de correlación de un e-fold).
- θ = 20 — exploración (borde superior de la malla, ruido rápido).

## Definición de X en el código del ajuste BAO

- `scripts/desi_dr2_data.py::add_ou_kernel`: C_ij += S_i S_j σ_X² exp(−θ|x_i − x_j|), con
  x = ln(1+z) (`scripts/profile_sigma_x_desi.py`). σ_X² es la varianza estacionaria de X
  (no σ²/2θ). No hay discretización de SDE: es un proceso gaussiano con covarianza analítica.
- S(z) = ∂ ln D_V/∂Ω_Λ en la dirección plana Ω_m = 1 − Ω_Λ (H0 = 67.4, Ω_Λ = 0.685), valores
  fijos `_FALLBACK_S_Z`; recalculado aquí con la misma receta: coincide en |S| a 3 decimales
  (el código usa el signo opuesto, irrelevante para la covarianza).
- ⇒ en el código X = δΩ_Λ evaluado localmente en cada z_i (con Ω_m compensando).
- Los papers dan definiciones distintas: X = δΩ_Λ (`papers/stochastic-dark-energy-desi-dr2.md` l.70),
  X = δρ_DE/ρ_crit,0 sumado a E² (mismo paper, l.562) y X = δρ_Λ/ρ̄_Λ
  (`papers/unimodular-gravity-vacuum-smoothness.md` l.68).

## Inconsistencias a corregir en los papers

Tres definiciones distintas de X en el repo:

1. `papers/stochastic-dark-energy-desi-dr2.md` l.70 — X ≡ δΩ_Λ(x), x = ln a.
2. `papers/stochastic-dark-energy-desi-dr2.md` l.562 — H² = H0²[Ω_m(1+z)³ + Ω_Λ + X(z)], i.e.
   X = δρ_DE/ρ_crit,0 ⇒ δ ln ρ_DE = X/Ω_Λ0.
3. `papers/unimodular-gravity-vacuum-smoothness.md` l.68 — X ≡ δρ_Λ/ρ̄_Λ = δ ln ρ_DE.

Además, el kernel BAO del código (δα_i = S(z_i)·X(z_i), S = ∂ln D_V/∂Ω_Λ con Ω_m = 1 − Ω_Λ
compensando punto a punto) no define una cosmología única: cada bin ve un Ω_Λ (y un Ω_m)
distinto. Lo que se construye en `analisis_cmb/` es la primera versión física (una sola
historia ρ_DE(a), H(a)) del modelo.

Ajustes CPL del repo (no reajustados):
- `results/eos_cpl_desi_dr2/eos_cpl_summary.json`: w0 = −0.98990, wa = −0.01590, hecho con
  **errores diagonales** ("7 bins, diagonal errors"), no con la covarianza completa del perfil OU.
- `results/joint_w0wa_sigma/joint_w0wa_sigma.json`: w0 = −1.00165, wa = +0.01725. **Se contradice**:
  `"data": "DESI DR2 BAO 7 bins, full meas. cov, S(z) fixed"` pero `"cov_mode": "diagonal_figure6_errors"`.

## Decisiones de implementación (2026-10-03)

- **Lectura (a):** δ ln ρ_DE = X/Ω_Λ0, que conserva la amplitud del kernel que DESI acotó, así que la
  cota de Planck queda en las mismas unidades de σ_X. La lectura (b) (δ ln ρ_DE = X) difiere sólo
  en un factor 1/Ω_Λ0 ≈ 1.46 de amplitud: la cota bajo (b) = cota bajo (a) × 1.46. No se corre aparte.
- **Ω_Λ0** = 0.684727 (CAMB, Planck 2018 base ΛCDM con mnu = 0.06 eV, el valor por defecto
  de `set_cosmology`; la masa de neutrinos no venía especificada).
- **Normalización:** H0, ω_b, ω_c fijos a Planck 2018; ρ_DE(a) = ρ_Λ0·[1 + (X_s(a) − X_s(1))/Ω_Λ0],
  así que ρ_DE(1) es la de ΛCDM y el universo es plano. Se pierde el modo constante X(1), que es degenerado con
  Ω_Λ/H0 (el θ→0 que DESI prefiere y que w0-wa absorbe); lo que queda es la parte
  potencialmente distinguible. (CAMB sólo recibe w(a) y fija Ω_DE hoy por cierre, así que es
  consistente con esta normalización.)
- **Parámetros fijos** (H0=67.36, Ωb h²=0.02237, Ωc h²=0.1200, τ=0.0544, ns=0.9649,
  As=2.1e-9) en todas las curvas: es una simplificación; el paso siguiente sería dejarlos libres.
- **Fuera del rango z de DESI** el proceso corre estacionario hasta a = 1e-3 (no se congela).
- **Trayectorias:** GP estacionario en x = ln(1+z) con Cov = σ_X² exp(−θ|Δx|), semilla 20261003.
  Se muestrea con la recursión AR(1), que es exactamente el factor de Cholesky de esa covarianza
  (proceso markoviano). Comprobado: max|Cholesky − AR(1)|/σ_X = 2.5e-14 (N=500, θ=1). Se usa la
  recursión porque el Cholesky directo está mal condicionado a θ = 0.001.
- **Suavizado:** kernel gaussiano con desviación típica Δ en ln a; en los bordes de la malla
  (a = 1e-3) el kernel se trunca y se renormaliza; en a = 1 el kernel es simétrico gracias al relleno (ver abajo).
- **w(a)** = −1 − (1/3) d ln ρ_DE/d ln a con `np.gradient` sobre la curva suavizada →
  `DarkEnergyPPF.set_w_a_table`. Comprobación: ρ_DE que integra CAMB frente a la nuestra:
  max rel 1.2e-5 (N=2000) y 1.9e-4 (N=500).
- **θ = 20:** la longitud de correlación (1/θ = 0.05) es igual a Δ ⇒ resultado dominado por el
  suavizado, no por los datos (marcar en la tabla).
- **mnu = 0.06 eV** (un neutrino masivo), la que asume Planck 2018 base ΛCDM (Planck 2018 VI, pie
  de la tabla 1: "Ω_m includes the contribution from one neutrino with a mass of 0.06 eV").

## Escala de resolución Δ

El valor de X en un punto (y por tanto X_s(1), que fija el cero de toda la historia ρ_DE(a))
depende de la escala de suavizado como ~σ_X·√(θΔ): la trayectoria OU no es derivable y no
tiene valor puntual "bien definido" por debajo de su escala de rugosidad. No hay forma de
evitarlo con un proceso OU. Por tanto **Δ es un parámetro declarado del modelo físico**, no un
detalle numérico; en la teoría original correspondería al volumen de la celda causal de Sorkin.

Las 200 realizaciones se corren con Δ = 0.05 y Δ = 0.10 para cada θ; la tabla final da, por
observable, la diferencia entre las bandas 68 % de ambos Δ en unidades de σ_Planck. Si en θ_*
supera 0.5σ, es el segundo resultado principal: el modelo no predice θ_* sin fijar Δ.

## Relleno más allá de a = 1

El proceso estacionario se genera también en ln a ∈ (0, 0.4] (= 4·Δ_max, x negativo), sólo para
que el kernel de suavizado sea simétrico en a = 1. La derivada de w también se toma con a = 1
como punto interior. A CAMB sólo se pasa a ≤ 1. El borde a = 1e-3 sigue truncado (ρ_DE/ρ_tot ~ 1e-9 allí).

Efecto en la realización única θ=1 (N=2000; `resultados/realizacion_unica/convergencia.json`):
100θ_*(Δ=0.10) − 100θ_*(Δ=0.05) pasó de +2.5 σ_Planck (sin relleno) a −0.515 σ_Planck.

## Convergencia de malla

### Sin relleno (superado)
500 vs 2000, θ=1: TT 1.6e-4 / 7.4e-4, θ_* 1.08e-2 / 5.10e-2 (Δ = 0.05 / 0.10). Falla θ_*.

### Con relleno, 500 vs 2000 (criterio pedido) — FALLA en θ_*
`resultados/realizacion_unica/convergencia.txt` (θ=1, semilla 20261003):

| Δ | TT max/σ | φφ max/σ (bandpowers) | θ_*/σ |
|---|---|---|---|
| 0.05 | 2.7e-4 ✔ | 5.0e-4 ✔ | 1.38e-2 ✘ |
| 0.10 | 7.4e-4 ✔ | 1.0e-3 ✔ | 5.48e-2 ✘ |

Causa: la misma malla servía para (i) integrar el suavizado sobre una trayectoria rugosa y (ii)
tabular w(a) para CAMB. El error de (i) en X_s(1) es ruido de cuadratura de media cero: con 400
trayectorias, X_s(1)[500] − X_s(1)[2000] tiene media −1.3e-5 / −6.3e-6 y desviación típica 3.2e-4 /
2.4e-4 (Δ = 0.05 / 0.10). Como θ_* responde a ese desplazamiento de ρ_DE en toda la historia,
es la única observable afectada.

### Criterio vigente: separar malla del camino y tabla de CAMB (aprobado 2026-10-03)

El criterio original (500 vs 2000 puntos con la misma malla para camino y tabla) **se sustituye**
por los dos diagnósticos siguientes, `conv_fina` (resolución del camino) y `conv_tabla`
(resolución de la tabla de CAMB). Motivo: la prueba original mezclaba dos errores distintos y
el dominante (ruido de cuadratura al suavizar una trayectoria rugosa muestreada a 500 puntos)
no es un error de la tabla que recibe CAMB. Producción: camino con paso h_tab/16, tabla de 2000 puntos.
Camino generado en malla uniforme con paso h_tab/k (h_tab = paso de la tabla de 2000 puntos),
suavizado con `scipy.ndimage.gaussian_filter1d` renormalizado, y la tabla para CAMB submuestreada
de él. 10 semillas (20261003…20261012), θ=1, máximo sobre semillas, en unidades de σ_Planck:

| Prueba | Δ | TT | φφ | θ_* |
|---|---|---|---|---|
| camino k=4 vs k=16 | 0.05 | 1.35e-4 | 1.17e-4 | 8.80e-3 |
| camino k=16 vs k=64 | 0.05 | 4.84e-5 | 3.96e-5 | 3.48e-3 |
| camino k=4 vs k=16 | 0.10 | 1.09e-4 | 9.90e-5 | 7.21e-3 |
| camino k=16 vs k=64 | 0.10 | 3.12e-5 | 2.88e-5 | 2.22e-3 |
| tabla 2000 vs 7997 (camino k=64) | 0.05 | 1.54e-5 | 1.06e-4 | 1.11e-3 |
| tabla 2000 vs 7997 (camino k=64) | 0.10 | 9.25e-6 | 1.37e-5 | 5.91e-4 |

Fuentes: `resultados/diagnosticos/conv_fina.txt`, `resultados/diagnosticos/conv_tabla.txt`
(scripts en `diagnosticos/`). Con camino k=16 y tabla 2000 todo queda < 1 % de σ_Planck
(peor caso θ_*: 3.5e-3 por el camino y 1.1e-3 por la tabla).

## Ficheros de Planck y valores de referencia

- TT: `datos_planck/COM_PowerSpect_CMB-TT-full_R3.01.txt` (PLA, sha256 ccf31136…5522).
- Lensing: `COM_Lensing_4096_R3.00.tgz` no trae bandpowers (sólo mapas klm, nlkk y máscara); se
  borró. Se usa `smicadx12_Dec5_ftl_mv2_ndclpp_p_teb_consext8_bandpowers.dat` (+ ventanas
  `..._window/window{1..9}.dat`) de `cobaya-install planck_2018_lensing.native`
  (planck_supp_data_and_covmats v2.1): 9 bins, L = 8–400, unidades [L(L+1)]²C_L/2π. Comprobación:
  ΛCDM binado con las ventanas queda a ≤ 1.4σ de cada bandpower. Se usa la columna `Error`
  como σ (diagonal) y **no** se aplica la corrección lineal de fiducial (se cancela casi entera
  en diferencias OU − CPL; no se ha verificado cuantitativamente).
- θ_*: 100θ_* = 1.04110 ± 0.00031 verificado en Planck 2018 VI (arXiv:1807.06209v4), tabla 2,
  columna TT,TE,EE+lowE+lensing (fuente LaTeX `ms.tex`, fila `100\theta_\ast`). Mismo valor en
  la tabla 1 (Plik, columna [1]).
- θ_MC ≠ θ_*: en la tabla 1 de Planck, 100θ_MC = 1.04092 ± 0.00031 frente a 100θ_* = 1.04110; en CAMB
  ΛCDM con los parámetros fijados: `cosmomc_theta` → 1.040982, `thetastar` → 1.041198
  (diferencia 0.00022 ≈ 0.7σ). Se usa `thetastar`.

## Curvas de referencia w0-wa

- **CPL del repo** (primaria): `results/eos_cpl_desi_dr2/eos_cpl_summary.json`, w0 = −0.98990,
  wa = −0.01590. ⚠ Ese ajuste usa **errores diagonales** ("7 bins, diagonal errors"), no la
  covarianza completa del perfil OU. El otro CPL del repo, `results/joint_w0wa_sigma/joint_w0wa_sigma.json`
  (w0 = −1.0017, wa = +0.0172), dice "full meas. cov" en `data` pero `cov_mode:
  diagonal_figure6_errors`: contradictorio, no se usa.
- **DESI DR2 publicado**: arXiv:2503.14738, ec. `eq:w0wa_DESI_CMB` (main.tex l.883):
  w0 = −0.42 ± 0.21, wa = −1.75 ± 0.58, combinación **DESI+CMB** sin SNe, donde CMB = Planck PR4
  (simall + Commander ℓ<30 + CamSpec ℓ≥30) + lensing Planck+ACT DR6. Son medias marginalizadas de un
  ajuste con los demás parámetros libres; aquí se evalúan con los parámetros de Planck 2018 fijos,
  así que esta curva no es la cosmología de DESI, sólo su w(a).
- **ΛCDM** como control.

## Estadístico y cota (decisión 2026-10-03)

La media de las realizaciones es ciega a σ_X (X tiene media cero), así que:
- Referencia = **ΛCDM** (con σ_X→0 preferido, el CPL del repo es ΛCDM a efectos prácticos y
  además está ajustado con errores diagonales). Se da también la comparación con los dos CPL.
- Por realización r y observable: Δχ²_r = χ²(OU_r) − χ²(ΛCDM) sobre los bins de Planck con los
  errores de Planck (TT: error superior si modelo > dato, inferior si no).
- Se guardan los 200 Δχ²_r por (θ, Δ, σ_X) y los percentiles 50 y 95.
- Cota: el mayor σ_X con percentil95(Δχ²_r) < 4, es decir, Planck no excluye a más de 2σ al
  menos el 95 % de los universos que el modelo genera con esa varianza.
- Números aleatorios comunes (semilla 20261003) en todo θ, Δ y σ_X: X_r ∝ σ_X exactamente.
- Banda del 68 % frente a ΛCDM: "hueco" = distancia de ΛCDM a la banda [p16, p84] (0 si ΛCDM
  queda dentro), y "borde" = borde de la banda más alejado de ΛCDM, ambos en σ_Planck, máximo sobre bins.
- Lensing: no lineal de CAMB por defecto con `lens_potential_accuracy=1` (NonLinear_lens, HMCode mead2020).

## Ejecución σ_X = 2.506e-2 (200 realizaciones × 6 (θ, Δ)) — problemas abiertos

Fuentes: `resultados/sigma_2.506e-02/{resumen.json, tabla.md, banda_delta.txt, observables_*.npz, dchi2_*.txt}`.

- **ρ_DE ≤ 0 (lectura (a) rota): 0 realizaciones** en los seis casos.
- **Fallos de CAMB** (`CAMBError: HMCode INTEGRATE, Integration timed out`, lensing no lineal
  mead2020): θ=1/Δ=0.05: 1; θ=20/Δ=0.05: 46; resto: 0. Sin espectros ni θ_* para esas filas.
- **NaN silenciosos**: CAMB devuelve TT y φφ enteros en NaN sin lanzar error (θ_* sí válido):
  θ=20/Δ=0.10: 1 (índice 38); θ=20/Δ=0.05: 5 (70, 89, 120, 183, 193). En `resumen.json` los
  percentiles de TT/φφ los excluyeron sin contarlos y las columnas "vs CPL" salen NaN en θ=20:
  **pendiente de corregir en `realizaciones.py`** (contarlos como fallo por observable).
- θ=20/Δ=0.05 queda con 149 realizaciones válidas en TT/φφ y 154 en θ_*: **no es una muestra
  representativa** (los fallos no son aleatorios a priori).
- **θ_* con parámetros fijos**: la señal en θ_* (p95 Δχ² ≈ 67–70 con θ=1) se calcula con H0, ω_b, ω_c
  fijos. θ_* es casi degenerado con H0 en la dirección geométrica; dejando H0 libre buena parte del
  desplazamiento podría absorberse. Indicio: el CPL del repo (w0 = −0.990) ya da Δχ²(θ_*) = 3.8 y
  DESI DR2+CMB da 1.4e3 sólo por fijar los parámetros.

## Cambio de método: H0 resuelto por θ_* (decisión 2026-10-03)

- Con H0 fijo, el desplazamiento de θ_* infla el resultado: θ_* es casi degenerado con H0 (dirección
  geométrica). Prueba: CPL repo (w0 = −0.990) daba Δχ²(θ_*) = 3.8 y DESI DR2+CMB, 1.4e3.
- **Ahora**, para cada curva (cada realización OU y las de referencia ΛCDM, CPL repo, CPL joint,
  DESI DR2+CMB) se resuelve H0 por bisección sobre `camb.get_background` (sin espectros; ω_b,
  ω_c, mnu fijos, plano, Ω_DE = 1 − Ω_m − Ω_r) tal que 100θ_* = 1.041198 (ΛCDM de referencia,
  calculado con los mismos ajustes). Tolerancia en H0: 1e-5. Luego TT y φφ con ese H0.
  `ou_cmb.solve_H0_for_thetastar`; θ_* de fondo frente al de `get_results`: diferencia 4e-8.
- **La cota ya no puede venir de θ_*** (queda fijado por construcción). Observables: TT 2≤l<30,
  φφ (L 8–400) y la **distribución de H0** resuelto.
- Comprobación externa: DESI DR2+CMB (w0 = −0.42, wa = −1.75) da H0 = 63.71, en línea con el H0 bajo
  que DESI obtiene en w0wa.
- θ_* con H0 = 67.36 se sigue calculando (por `get_background`) para todas las realizaciones,
  incluidas las que fallen en espectros, sólo como información.
- **H0 de referencia**:
  - Planck 2018: 67.36 ± 0.54, arXiv:1807.06209 tabla 2, columna TT,TE,EE+lowE+lensing (también
    ec. `BAO_H01`). Verificado en la fuente LaTeX.
  - DESI DR2 BAO + BBN (ΛCDM): 68.51 ± 0.58, arXiv:2503.14738, sección `sec:cosmologcal_constraints`
    (ecuación sin etiqueta justo antes de `eq:H0_BAO+BBN+thetastar`, y tabla de resultados ΛCDM,
    fila DESI+BBN). Las conclusiones del mismo paper dicen 68.50 ± 0.58 (incoherencia menor del
    paper). **Sólo comparación, no cota**: el BAO ya entró en el ajuste de σ_X, y el valor asume ΛCDM.
- **Lensing no lineal**: halofit `takahashi` en todas las curvas, ΛCDM incluida (mead2020 fallaba en
  47/200 con θ=20/Δ=0.05). Regla: si fallan > 5 de 200 en alguna combinación, pasar a lineal en todas.
  No verificado: cómo trata halofit una tabla w(a) en PPF (puede usar un w efectivo); sólo afecta a la
  corrección no lineal de φφ.
- **Fallos**: NaN silenciosos en espectros = fallo de esa observable, contados junto a los errores
  de CAMB (`fallos_*.json`). Versión anterior del script (H0 fijo, mead2020):
  `diagnosticos/obsoleto/realizaciones_H0fijo_mead2020.py`; sus resultados en `resultados/sigma_2.506e-02/`
  son **obsoletos**.

## Resultado σ_X = 2.506e-2, método H0 resuelto + takahashi (2026-10-03)

Fuente: `resultados/takahashi/sigma_2.506e-02/{tabla.md, resumen.json}`.
- Fallos: 0 en las seis combinaciones (ρ≤0, CAMB, NaN silenciosos, H0). Se queda takahashi.
- TT 2≤l<30: Δχ² p95 ≤ 0.064 en todo (θ, Δ); φφ: p95 ≤ 0.099. Fracción con Δχ² > 4: 0.
- H0 resuelto, θ=1: banda 68 % [66.94, 67.85] (Δ=0.05); fuera de Planck ±2σ: 2.5 %.
- Fuera de DESI+BBN ±2σ: 28–30 % ya con θ=0.001 (≈ΛCDM), porque el propio H0 de Planck ΛCDM
  (67.36) está a 1.98σ de 68.51 ± 0.58: esa columna mide la tensión Planck–DESI en ΛCDM, no el OU.
- Sensibilidad a Δ: todas < 0.5σ_Planck (máximo H0 con θ=20: 0.36σ).

## Malla de σ_X (decisión 2026-10-03)

- 10 puntos logarítmicos entre 2.506e-2 y 1e-4, sólo θ = 0.001 y θ = 1, sólo Δ = 0.05.
- **Descartados como no informativos**: θ = 20 (dominado por el suavizado, 1/θ = Δ) y Δ = 0.10
  (sensibilidad a Δ < 0.5σ_Planck en todas las observables a σ_X = 2.506e-2).
- Salida: `resultados/takahashi/malla/sigma_*/` y `resultados/takahashi/malla/{malla.json, malla.md}` (`malla.py`).
- σ(H0) ≡ (p84 − p16)/2. Leyes de escala por mínimos cuadrados en log10–log10.
- SH0ES para la conclusión (c): H0 = 73.04 ± 1.04, arXiv:2112.04510 (resumen). Planck VI cita el
  valor anterior 74.03 ± 1.42 (ms.tex l.3122); se usa el más reciente.

## Resultado de la malla (2026-10-03)

Fuente: `resultados/takahashi/malla/{malla.json, malla.md}`. 0 fallos en los 20 casos.
- p95 < 4 en toda la malla; máx p95 = 0.0636 (TT) y 0.0514 (φφ), θ = 1, σ_X = 2.506e-2.
- p95(Δχ²) ∝ σ_X^n con **n ≈ 1** (0.998–1.007), no n ≈ 2 como se esperaba. σ(H0) ∝ σ_X^m con m ≈ 1.
- n ≈ 1 **verificado** (θ=1, σ_X = 2.506e-2; `diagnosticos/descomposicion_dchi2.py` →
  `resultados/takahashi/sigma_2.506e-02/descomposicion_dchi2.{txt,json}`): con r = ΛCDM − dato y δ = OU − ΛCDM,
  Δχ² = 2·Σ r_b δ_b/σ_b² + Σ δ_b²/σ_b² exactamente (resto < 1e-14; ningún bin cambia de lado del error
  asimétrico de TT). p95 lineal / cuadrático: TT 0.0632 / 0.00063; φφ 0.0502 / 0.0074. El p95 medido está
  dominado por el término cruzado con los residuos de Planck (∝ σ_X); la contribución propia del modelo
  (cuadrática, ∝ σ_X²) es 1–2 órdenes de magnitud menor.
- Resumen final: `RESULTADO.md`.

## Oscilador determinista (`oscilador/`, 2026-10-03, sin commit)

- X(ln a) = A sin(2π ln a / P + φ), 60 curvas; mismo pipeline; derivada analítica, sin Δ.
- BAO: **las funciones de distancia del repo no se pueden reutilizar sin modificarlas** (sólo CPL, H0 y Ω_m
  fijos, sin r_d). Se reutiliza el cargador `load_gaussian_bao_vector` (13 D/r_d, cov 13×13) y las
  distancias salen de CAMB. El cargador reescribe `results/desi_cov/desi_gaussian_bao_ALL_GCcomb_cov.txt`
  con contenido idéntico (git limpio).
- ΛCDM con parámetros Planck: χ²_BAO = 29.84 (13 puntos). "Permitida" = Δχ²_BAO < 4 frente a esa ΛCDM.
- Resultado (`oscilador/resultados/tabla.md`): 0 fallos; 58/60 permitidas (fuera: A=0.025, P=1, φ=0,
  Δχ²_BAO = 4.23; A=0.025, P=2, φ=3π/2, 5.67); 0/60 distinguibles en CMB (máx Δχ²_TT = 0.070,
  máx Δχ²_φφ = 0.631).
