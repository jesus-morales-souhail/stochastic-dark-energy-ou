# Lectura: Everpresent Λ

Fuentes LaTeX descargadas de arXiv el 2026-10-16 (copia local en `fuentes/`, no versionada):

- **[A]** Ahmed, Dodelson, Greene, Sorkin, "Everpresent Λ", arXiv:astro-ph/0209274
  (fichero `v6.tex`, fechado 2002-09-13). Secciones: I Introduction, II Causet Theory,
  III Simulations, IV Complications, V Conclusion. Ecuaciones numeradas (1)–(8).
- **[Z]** Zwane, Afshordi, Sorkin, "Cosmological Tests of Everpresent Λ", arXiv:1703.06265
  (fichero `fcc_3.tex`, fechado 2018-05-23). Ecuaciones numeradas por sección (2.1), (2.2)…
  Secciones: 1 Introduction, 2 Everpresent Λ (2.1 Model 1, 2.2 Model 2, 2.3 Inhomogeneities),
  3 CMB, 4 Local Measurements of the Hubble Constant, 5 BAO Measurements, 6 Ultramassive
  Black Holes, 7 Lithium-7, 8 Conclusion.

Las citas literales van entre comillas. Lo marcado "(derivado)" no está escrito en el paper:
se sigue directamente de la ecuación citada.

---

## [A] Ahmed, Dodelson, Greene, Sorkin

### (a) Variable estocástica, tiempo, tipo de proceso, varianza, correlación

- **Variable.** Λ, implementada como ρ_Λ. Unidades ħ = c = m_p = 1 con m_p la masa de Planck
  *reducida*, (8πG)^{-1/2} (Sec. I, nota 1), de modo que Λ y ρ_Λ coinciden. La ec. (6) es la
  ecuación de Friedmann H² = (ρ_matter + ρ_radiation + ρ_Λ)/3; la ec. (5) define el proceso
  sobre ρ_Λ.
- **"Tiempo".** El número de elementos N, identificado con el 4-volumen V: "the ever growing
  number N of causet elements plays the role of a kind of parameter time" (Sec. II). En la
  simulación se avanza en tiempo cósmico t_i y el incremento es
  δN_i = V(t_{i+1}) − V(t_i) (ec. 4).
- **Tipo de proceso.** "let the fluctuations in Λ be driven by those of an unadorned random
  walk" (Sec. II). Λ(N) = S(N)/N, con S(N) la suma de los N primeros números aleatorios
  (Sec. II, último párrafo). Ec. (5):

  ρ_{Λ,i+1} = S_{i+1}/N_{i+1} = (S_i + α ξ_{i+1} √δN_i) / (N_i + δN_i),
  ξ de media 0 y desviación 1, S_0 = 0 en un instante t_0 muy temprano.

  Es decir: S es un paseo aleatorio en N; Λ es ese paseo dividido por N.
- **Varianza.** Escrita solo como orden de magnitud: ΔΛ ∼ 1/ΔV ∼ 1/√V (ec. 1), y
  ρ_Λ ∼ V^{-1/2} ∼ H² ∼ ρ_critical (ec. 2). (Derivado de la ec. 5: Var ρ_Λ(N) = α² (N − N_0)/N²
  ≈ α²/N.)
- **Correlación entre épocas.** No dan fórmula. Dicen que hace falta: "we need to know, not only
  their typical magnitude at each moment of cosmic time, but also how the fluctuations at one
  moment correlate with those at other moments" (Sec. II). (Derivado de la ec. 5, con
  N_0 → 0: Cov[ρ_Λ(N_1), ρ_Λ(N_2)] = α² min(N_1,N_2)/(N_1 N_2) = α²/max(N_1,N_2); proceso no
  estacionario.)

### (b) Parámetro de amplitud

- **Nombre.** α, "an unknown dimensionless parameter which governs the dynamics of the theory"
  (tras ec. 5).
- **Valor "natural".** No dan un valor natural de α. Lo más cercano: Λ como "action per element"
  "would be expected to be of order unity in fundamental units", y tras la cancelación de signos
  queda un residuo 1/√N (Sec. II). En Sec. III: "If α≈1, there will inevitably be times during
  which the total effective energy density … goes negative, thereby invalidating the equation."
- **Valor usado / "ajuste".** No ajustan datos. Eligen α por simulación: α = 0.01 (Fig. 2,
  histograma de 6000 realizaciones en Fig. 4) y α = 0.02 (Fig. 3). Criterio: que la simulación no
  llegue a ρ_tot < 0 y que Ω_Λ hoy sea suficiente: "this happens reasonably often for α in the
  range 0.01−0.02" (Sec. III). Con α = 0.02 desde a = 10^{-32}, "only about a 1 in 3 chance of
  completing the run". Con α = 0.02, Ω_Λ "fluctuates about zero with an amplitude of order 0.5".
  Combinación de datos: ninguna formal; la comparación con observaciones es cualitativa
  (Ω_Λ ≳ 0.5 hoy y casi constante hasta z = 1; BBN mencionado cualitativamente).

### (c) Tensión de H0 y BAO Lyα

- H0: no lo dicen.
- BAO Lyα: no lo dicen. (Lo único observacional además de la aceleración es BBN, Sec. III:
  un ρ_Λ negativo reduciría el ⁴He predicho.)

### (d) Conteo de N

- **4-volumen.** "the four-volume of a region of spacetime reflects the number of discrete
  elements" (Sec. II); "the equality N = V between volume and number of elements" (Sec. II).
  V se identifica "(up to an unknown factor of order unity) with the volume of the past light
  cone of any representative point on the hypersurface" (Sec. II). Fórmula, ec. (3):

  V(t) = (4π/3) ∫_0^t dt' a(t')³ [ ∫_{t'}^{t} dt''/a(t'') ]³

  en unidades ħ = c = m_p(reducida) = 1, y N = V (ec. 4).
- Área: no lo usan.

---

## [Z] Zwane, Afshordi, Sorkin

Dos modelos. Las simulaciones y el MCMC usan el Modelo 2: "Since Model 2 is somewhat easier to
simulate than Model 1, we shall adopt it in our Monte Carlo Markov Chain (MCMC) exploration"
(Sec. 2.2). Del Modelo 1 dicen que lo simularon "only enough to confirm that it behaves
qualitatively like Model 2" (Sec. 8).

### (a) Variable estocástica, tiempo, tipo de proceso, varianza, correlación

**Modelo 1 (Sec. 2.1, atribuido a [A]).**
- Variable: ρ_Λ. Ec. (2.8): ρ_Λ(t+1) = [S(t) + β ξ √(N(t+1) − N(t))] / V(t). Ojo: aquí el
  denominador es V(t), no N como en [A] ec. 5. La ecuación de Friedmann usada es (2.5):
  (ȧ/a)² = ρ/3 − Λ(t)/3 (con signo menos delante de Λ, tal cual está escrita).
- Tiempo: iteraciones en tiempo cósmico t, con N(t) = V(t)/ℓ_p⁴ (2.7).
- Proceso: "ξ is a gaussian random number which drives the random walk" (tras 2.8).
- Varianza y correlación: no las escriben.

**Modelo 2 (Sec. 2.2).**
- Variable: Ω̂_de, "the auxiliary classical gaussian random process with mean zero", y la
  densidad física es Ω_de(λ) = tanh[Ω̂_de(λ)] (2.10), para que |Ω_de| < 1. "We conveniently
  build in the 'everpresent' nature of Λ by expressing the fluctuations directly in terms of
  Ω_Λ rather than in terms of the more basic ρ_Λ as in Model 1."
- Tiempo: λ = log(a) (tras 2.9). Lo describen como "a random gaussian function of logarithmic
  time".
- Tipo de proceso: gaussiano, media cero, estacionario en λ, definido por su función de
  correlación. No le dan nombre. Se simula en Fourier (2.11), muestreando
  Ω̃_de(ω) como gaussianas independientes con ⟨Ω̃Ω̃⟩ = α² e^{−ω²/4μ} √(π/μ) δ(ω_1 − ω_2).
- Correlación (2.9):

  ⟨Ω̂_de(λ_1) Ω̂_de(λ_2)⟩ = α² exp[ −μ (λ_1 − λ_2)² / 2 ]

  (gaussiana en Δλ, no exponencial).
- Varianza: α² según (2.9). El texto dice "with zero mean, variance α, and a correlation time of
  Hubble Time/√μ" (Sec. 2.2) — escrito así, "variance α", aunque la ecuación da α².

### (b) Parámetro de amplitud

- **Modelo 1:** β. "Ahmed et al. found that one needs roughly 0.01 < β < 0.02 to generate a
  universe consistent with astrophysical observations" (Sec. 2.1). El criterio que citan para
  β pequeño es la edad: "a value of β which is too small results in a universe which is too
  young compared to the oldest known galaxies, which are about 12 Gyr" (el criterio de [A]
  era que Ω_Λ hoy no fuese suficiente). No ajustan β.
- **Modelo 2:** α (amplitud) y μ (coherencia).
  - Valor natural de α: no lo dicen. De μ: la coherencia se "presumed to be comparable to the
    Hubble time, just as one sees in Model 1" (Sec. 2.2). (Una frase "On general grounds one
    would expect μ to be of order unity" existe en la fuente pero comentada; no sale en el
    paper.)
  - Valor de ajuste: α = 0.8824, μ = 0.9804 (Tabla 1; Sec. 3: "for α = 0.8824 and μ = 0.9804
    there exists a history of dark energy that is a good fit to the CMB + BAO"). Sin barras de
    error. Se exploran (α, μ, seed) con CosmoMC + CAMB modificado.
  - Datos (pie de Tabla 1): Planck 2015 baseline + lensing reconstruction + BAO (6DF, MGS,
    DR11CMASS, DR11LOWZ, DR11LyaAuto, DR11LyaCross). χ²: 11335.2 (Everpresent) vs 11334.6
    (ΛCDM). En la misma tabla: H0 = 68.48 ± 0.67, Ω_Λ = 0.6781 ± 0.0081 (Everpresent) frente a
    H0 = 67.79 ± 0.54, Ω_Λ = 0.6920 ± 0.0072 (ΛCDM).
  - "Since the typical value of |Ω_de| is 0.7 for this model, there is nothing unusual about the
    current density of dark energy." Probabilidad de |Ω_de| < 0.25 cerca de z ≃ 1000: "∼23%"
    (Sec. 3). La mejor realización está a 1.8σ en χ²_red, ΛCDM a ∼4.3σ (Sec. 3, Fig. 3).

### (c) Tensión de H0 y BAO Lyα

- **Abstract:** "it removes the observational tensions that ΛCDM experiences in relation to low
  redshift measurements of the Hubble constant, and to the Baryonic Acoustic Oscillations (BAO)
  in Lyman-α forest at z∼2−3."
- **H0 (Sec. 4):** citan Riess et al. H0 = 73.24 ± 1.74 frente a 67.8 ± 0.9 de Planck, "3.4σ".
  "Everpresent Λ fares better than ΛCDM, because with a fluctuating dark energy, local estimates
  of H0 are less constrained by high redshift data like CMB or BAO. As seen in Figure 4, higher
  values of H0 are consistent with the Planck+BAO data; and the SNe Ia measurements fall within
  2σ of Model 2." Pie de Fig. 4: "the tension between Planck+BAO and SNe Ia is relieved in the
  Everpresent Λ model". Intro (Sec. 1): "this tension is removed in the Everpresent Λ model,
  thanks to its stochastic nature."
- **BAO Lyα (Sec. 1):** BOSS "suggest that Ω_de(z) is negative at z=2.34 in ∼2.5σ tension with
  standard ΛCDM … This is not a problem for Everpresent Λ, since it allows for Ω_de(z) to change
  sign roughly once per Hubble time."
- **BAO Lyα (Sec. 5):** "the dark energy history in Everpresent Λ that best fits the Planck 2015
  CMB data falls well within the BOSS bounds, which on the other hand are in 2.5σ tension with
  ΛCDM. Indeed, if we fix the cosmological parameters to their best fit values in Table 1, there
  is a 78% chance that Ω_de is negative somewhere in the range 1.5 < z < 3." Tabla 2: χ² de
  DR11LyaAuto 4.265 → 2.636 y DR11LyaCross 4.748 → 3.945 (ΛCDM → Everpresent). También
  (Sec. 5): los análisis BAO usan un modelo fiducial sin energía oscura antes de la
  recombinación, "it seems possible that most analyses are introducing an implicit bias against
  an early dark energy … the recent BAO analyses are compatible with Everpresent Lambda or even
  favour it."

### (d) Conteo de N

- **4-volumen.** (2.1): N ∼ V ± √V "In Planck units (8πG = ħ = c = 1)". (2.6)–(2.7):

  V(t) = (4π/3) ∫_0^t dt' a(t')³ ( ∫_0^{t'} dt''/a(t'') )³,   N(t) = V(t)/ℓ_p⁴

  El límite de la integral interior es ∫_0^{t'} (horizonte de partículas en t'), mientras que en
  [A] ec. (3) es ∫_{t'}^{t} (radio del cono de luz pasado en t'). El texto dice "the volume of
  the past lightcone of a point".
- Área: no la usan. (En la fuente hay una línea comentada `N(t) = V(t)/ℓ_p^2` con la nota
  "RDS i assume this was supposed to be l^4 not l^2"; no sale en el paper.)

---

## Comparación con el OU del repo

Repo: `scripts/desi_dr2_data.py::add_ou_kernel` añade a la covarianza de medida de los 7 α
isotrópicos de DESI DR2 el término S_z,i S_z,j σ_X² exp(−θ|x_i − x_j|), con x = ln(1+z_eff);
X se lee como δΩ_Λ; θ se perfila libre en (10⁻³, 20)
(`scripts/profile_sigma_x_desi.py`); cota σ_X < 2.506 × 10⁻² (95 %).

| | [A] | [Z] Modelo 2 | Repo |
|---|---|---|---|
| Variable | ρ_Λ = Λ (unidades 8πG = 1), densidad absoluta (ec. 5) | Ω̂_de, con Ω_de = tanh Ω̂_de (2.9–2.10) | X = δΩ_Λ, residuo sobre ΛCDM |
| Media | 0 ("fluctuations about zero", Sec. II) | 0 (Ω̂_de de media cero) | 0 para X; Ω_Λ total = Ω_Λ0 + X con Ω_Λ0 ≈ 0.685 |
| Variable temporal | N = V (4-volumen del cono pasado); pasos en t | λ = ln a | x = ln(1+z) = −ln a |
| Tipo de proceso | S paseo aleatorio en N, Λ = S/N; no estacionario | gaussiano estacionario, correlación gaussiana (squared-exponential) en Δλ | gaussiano estacionario, Ornstein–Uhlenbeck (exponencial en Δx) |
| Amplitud | α = 0.01–0.02 elegida por simulación, sin ajuste | α = 0.8824 (σ de Ω̂_de) en el mejor ajuste Planck15+lensing+BAO; sin barras | σ_X < 2.506 × 10⁻² (95 %, DESI DR2) |
| Correlación | no la escriben; (derivado) α²/max(N_1,N_2) | α² exp[−μ(Δλ)²/2], μ = 0.9804; "correlation time Hubble Time/√μ" | σ_X² exp(−θ|Δx|), θ libre |
| Acotación | ninguna; si ρ_tot < 0 se corta la simulación | tanh comprime a (−1, 1) | ninguna |

Cifra para la fila "Amplitud": 0.8824 / 2.506 × 10⁻² ≈ 35 (α de [Z] frente a la cota σ_X; α es
la σ de Ω̂_de antes del tanh, tanh(0.8824) ≈ 0.71).

---

## [D1] Das, Nasiri, Yazdi, "Aspects of Everpresent Λ (I)", arXiv:2304.03819

Fuente `EverLam1.tex` (fichero de arXiv fechado 2024-02-22). Ecuaciones numeradas por sección.
Secciones: 2 Background (2.1 Causal sets, 2.2 Unimodular gravity), 3 Everpresent Λ (3.1 Model 1),
4 Methodology (4.1 FLRW, 4.2 Numerical Recipe, 4.3 H²<0 Crashed Runs, 4.4 ρ_m as a Scaling
Parameter), 5 Results (5.1 Statistical Distributions, 5.2 Frequency of Completed Runs,
5.3 Initial Conditions, 5.4 Auto-Correlation Time), 6 Discussion.

- **N y V.** Correspondencia número–volumen: ⟨N⟩ = ρ_cs V, ρ_cs ≡ 1/ℓ_cs⁴, con N poissoniano
  (2.1) y δV = √N/ρ_cs (2.2). V es el 4-volumen del cono de luz pasado de un punto (Sec. 3.1, tras
  3.11), (4.4): V(τ) = (4π/3) ∫_{τ_i}^{τ} dτ' a(τ')⁴ (τ − τ')³, equivalente en tiempo cósmico a [A]
  ec. (3) (nota al pie de 4.4).
- **Incertidumbre.** (3.4)–(3.6): δΛ/(8πG) · δV ≥ ħ/2 → δΛ = 4π (ℓ_p/ℓ_cs)² /√V, con ℓ_p = √(Għ).
- **α.** (3.8) δΛ = (8πα/ℓ_cs²)(1/√N), con (3.9) **α ≡ ½ (ℓ_p/ℓ_cs)²**.
- **S y Λ.** S_Λ = (ℓ_cs⁴/8πℓ_p²) Σ_x Λ|_x (3.11); S_Λ|_pasado de x = Λ(x)V(x)/(8πG) (3.12);
  Λ(x) ∼ 𝒩(0, 8πα/√V) (3.13). Recursión (3.14)–(3.15):
  Λ(x_n) = [Λ(x_{n−1}) V(x_{n−1}) + 8πα √ΔV_n ξ_n] / V(x_n), y sumada (3.16):
  Λ(x_n) = [8πG S_0 + 8πα Σ_{i=1}^{n} √ΔV_i ξ_i] / V(x_n).
  En la receta numérica (Sec. 4.2): S_i = S_{i−1} + α √(V_i − V_{i−1}) ξ_i / G (4.11),
  Λ_i = 8πG S_i / V_i (4.12). Es decir, α entra multiplicando a 8π en la fluctuación de Λ
  (geométrica, c = 1), sin ℓ_p: Λ = 8πα Σ√ΔV ξ / V.
- **Friedmann.** (4.2): H² = (8πG/3)(ρ_m + ρ_r) + Λ/3. Solo se usa esta ecuación (Sec. 4.1).
  En la receta, dτ/da = √(3 / [8πG a⁴(ρ_m + ρ_r) + a⁴ Λ_{i−1}]) (4.6)–(4.7), con Λ_{i−1} fijo en el
  paso. Pasos uniformes en a: "the value of a is i/n" (Sec. 4.2).
- **Criterio de muerte.** H² < 0 ("H² < 0 crash", Sec. 4.3): "we will disregard universes that reach
  H² < 0 before the present, and terminate their evolutions at this point."
- **Condición inicial.** Parámetros: S_0, V_0, n, i_initial (Sec. 4.2). Sec. 5.2: "we start the
  simulations from z ≈ 10⁴, and except in the next subsection, we calculate the volume up to this time
  using only the matter and the radiation component of the universe"; con arranque en z ∼ 10⁸ y paso
  logarítmico "the crashing fraction is almost 100% for α = 0.01". Sec. 5.3: V_0 solo entra en el
  denominador de (4.12), no en el numerador fluctuante; los resultados "do not appear to be sensitive"
  a V_0; S_0 > 0 añade un sesgo positivo persistente. El estudio de S_0 usa α = 0.015, V_0 = 1,
  i_initial = 1, n = 10000, S_0 de 0 a 1. El α = 0.02 de Sec. 5.2 es "if we take S_0 = 0".
- **Ω_m h².** Fijado a 0.14 en Sec. 5.1 ("For our universe Ω_m⁰h² ≈ 0.14"). Sec. 4.4: con S_0 = 0 y
  semilla fija, Λ, H² y 1/√V escalan con ρ_m⁰, así que Ω_m y Ω_Λ no dependen de Ω_m⁰h².
- **Fracciones de runs completados.** Pie de Fig. 3 (Sec. 5.1), 10000 semillas cada una:
  **9990** (α = 0.005), **3228** (α = 0.01), **535** (α = 0.012). Sec. 5.2: α = 0.02, "less than 10
  out of 10,000 for any given Δτ" (pasos constantes en tiempo conforme). Sec. 5.1 no da el número de
  pasos; [D2] Sec. 4.1 dice "when dividing the history of the Universe into 10,000 timesteps". Sec. 5.2:
  "as we decrease the size of the time step, the fraction reaching the present day increases (in fact
  it starts to plateau for the smallest values of Δτ)".
- **Otros números.** H0 medio cerca de 40 km/s/Mpc en los tres α; sin energía oscura H0 ≈ 37.7
  (Sec. 5.1). "H0⁴V = 0.033" para un universo sin energía oscura (Sec. 5.1). Autocorrelación ajustada
  ≈ exp(−3.2 (Δλ)^0.48) (Sec. 5.4, Fig. 5).

## [D2] Das, Nasiri, Yazdi, "Aspects of Everpresent Λ (II)", arXiv:2307.13743

Fuente `EverLam2.tex` (fichero de arXiv fechado 2024-10-07). Secciones: 2 Overview (2.1 Model 1,
2.2 Model 2), 3 Cosmological Observations, 4 Fitting Model 1 to Supernova Data (4.1 MCMC, 4.2
Diagnosing Good Seeds, 4.3 Mock Data), 5 Fitting Model 1 to CMB Data, 6 Constraining Dark Energy
Using CMB, 7 Discussion and Conclusion.

- **Modelo 1** (Sec. 2.1): V_t (2.1) igual que [A] ec. (3); Λ_t = [Λ_{t−dt}V_{t−dt} + 8πα ξ_t √(V_t −
  V_{t−dt})]/V_t (2.2); ρ_Λ = Λ/(8πG) (2.3); Friedmann (2.4): ȧ²/a² = (8πG/3)(ρ_t + ρ_Λ) =
  (8πG ρ_cr⁰/3)(Ω_m⁰a⁻³ + Ω_r⁰a⁻⁴ + ρ_Λ/ρ_cr⁰), con ρ_cr⁰ = 3H0²/(8πG) y "H0 … taken to be
  approximately 70 km/s/Mpc", una constante de normalización, no el H0 de la simulación. c = 1.
  Condición inicial: "we set it during the radiation-dominated era by specifying an initial non-zero
  four-volume V_0 (calculated using a radiation-dominated universe without dark energy)".
  Mismo criterio de muerte H² < 0.
- **Modelo 2** (Sec. 2.2): Ω_Λ = tanh Ω̂_Λ, ⟨Ω̂Ω̂⟩ = α̃² e^{−μ(λ_1−λ_2)²/2}, λ = log a.
- **Supernovas** (Sec. 4.1–4.2): Pantheon+SH0ES; 90000 semillas, un MCMC por semilla en
  (Ω_m⁰h², α, M_B), α ∈ [0.001, 0.025]. ΛCDM: χ² ≈ 1485.3. 55 semillas con χ² en [1481, 1490];
  **16 mejores que ΛCDM: "16 out of 90000, or about 0.017%"** (Sec. 4.2). Con mocks del propio
  Modelo 1: f_D entre 0.026 y 0.205, media 0.076 (Sec. 4.3).
- **CMB** (Sec. 5): 100000 semillas; tras cortes (H0 > 60 para algún α, Ω_m^min ≤ 0.45), 3106; 91
  con χ² < 3000; mejor χ² = 2495.56 frente a 2339.0 de ΛCDM (∼157 peor), con H0 entre 40 y 51.
  Con supresión de la energía oscura en la última dispersión, ∼20 peor (Sec. 7). Modelo 2 con 3000
  semillas: α̃ → < 0.1 o Ω_m⁰ > 0.9, H0 ≈ 40. Sobre Zwane: "their equation of state parameter w(z) was
  restricted to a narrow negative range … Therefore, our results should not be compared to those of
  this previous work" (Sec. 5.2).
- **Conclusión** (Sec. 7): "the most basic models of Everpresent Λ have difficulty describing both the
  SN Ia and CMB data."

## [B] Barrow 2007, "A Strong Constraint on Ever-Present Lambda", gr-qc/0612128, PRD 75, 067301

Solo el abstract (API de arXiv; fuente no descargada): "strongly constrained by the isotropy of the
microwave background. Fluctuations generated by stochastic lambda generation which are consistent with
COBE and WMAP observations are far too small to dominate the expansion dynamics at z<1000 and so cannot
explain the observed late-time acceleration of the universe." [Z] Sec. 2.3 lo resume como una cota
∼10⁻⁶ en las fluctuaciones de Λ, bajo la hipótesis de Λ(x) local con parches del CMB no correlacionados.
