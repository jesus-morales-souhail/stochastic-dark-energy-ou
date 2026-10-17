# Bloque C: definición única de X, factor S_iS_j, normalización del ruido, factor 1/Ω_Λ0; m_φ retirado (2026-10-17)

**Decisión del autor:** definición estándar **[X-std]**, tal como la propone `AUDITORIA.md` §1: *X ≡ δρ_DE/ρ_crit,0 (= δΩ_Λ en unidades de hoy); δ ln ρ_DE = X/Ω_Λ0; σ_X² ≡ Var X = σ²/(2θ).* Se valida por separado, numérica y analíticamente, la normalización del ruido ⟨ξ(t)ξ(t')⟩ = σ²H ρ_crit,0² δ_D(t−t').

**Fuentes:** `AUDITORIA.md` tabla 1 (§1) y plan (C); `verificacion_numeros.md` filas 2, 6–8, 10–15 (las de "fórmula deducida" que B2 dejó aparte). Constantes: σ_X = 2.632828e-2 (`results/profile_sigma_x/profile_sigma_x.json`), Ω_Λ0 = 0.685, Ω_m = 0.315 (`scripts/slip_bridge.py`), σ₀ = L_P/L_H = 1.1776e-61 (`scripts/gap_two_targets.py`).

**Reglas:** cambios sobre el texto existente, línea a línea. Cada cambio se aplicó con un script que exige que el texto de origen aparezca una sola vez y registra la línea, el antes y el después: lo de abajo es ese registro, no está escrito de memoria. **Los números de línea son los del fichero antes del bloque C** (commit `7327c69`).

**Rutas de la tabla 1 traducidas** (la tabla es anterior a la reorganización): `papers/unimodular-gravity-vacuum-smoothness.md` → `theory/`; `papers/fundamental-vs-emergent-vacuum-relaxation.md` → `theory/`; `papers/principle-of-vacuum-smoothness.md` → `theory/`; `notes/desqueezing-relaxation-vacuum-fluctuations-note.md` → `exploratory/notes/`. Los demás ficheros no se movieron.

## Resumen de números que cambian

| Magnitud | Antes | Después | Origen |
|:--|:--|:--|:--|
| rms δρ_Λ/ρ_Λ en la cota (σ_X/Ω_Λ0) | 2.6e-2 (se usaba σ_X directamente) | 3.84e-2 | 2.633e-2 / 0.685 = 3.8435e-2 |
| Órdenes de magnitud de supresión (stochastic l.372) | "more than four" | ∼1.6 respecto a ρ_crit (1.58); 1.42 respecto a ρ_Λ | log₁₀(1/σ_X), log₁₀(Ω_Λ0/σ_X) |
| Gap en densidad, techo DESI (amplification-gap l.54) | σ_tgt/σ₀ ∼ 10⁵⁹ | 3.3e59 (3.264e59) | σ_X/(Ω_Λ0 σ₀) |
| H.5: cota sobre ζ V/(M_Pl V') | 1.0e-4 | 7.0e-5 | 2.7e-6 / 3.84e-2 |
| λ (potencial exponencial) | ≳ 9.8e3 | ≳ 1.4e4 | 1 / 7.0e-5 |
| M_cut | ≳ 9.8e3 M_Pl | ≳ 1.4e4 M_Pl | M_Pl / ζ |
| slip máx. ε = 1, z = 0.5 / 1.0 / 1.5 | 3.39e-2 / 1.43e-2 / 7.33e-3 | 4.95e-2 / 2.09e-2 / 1.07e-2 | slip_bridge.py regenerado (×1/Ω_Λ0 = ×1.460) |
| ratio frente a Maus (0.17) | 5.0 / 12 / 23 | 3.4 / 8.1 / 16 | 0.17 / slip |
| ratio frente al suelo de forecast (0.05) | 1.5 / 3.5 / 6.8 | 1.0 / 2.4 / 4.7 | 0.05 / slip |
| C7 | slip ≲ 3e-2 | slip ≲ 5e-2, en el mejor caso a la altura del suelo 0.05 | ídem |
| m_φ (stochastic l.476, l.621–625; unimodular l.130) | ≲ 1e-5 eV, ≲ 9.45e-5 eV; V'/V M_Pl ≲ 1e-2 | retirado | verificacion_numeros filas 6–8 (no se reproduce) |

## Avisos

1. **Dos números de la instrucción no corresponden a la cota vigente; se usó el valor con 2.633e-2.**
   - El gap en densidad que se pidió, **3.1e59**, es σ_X/(Ω_Λ0 σ₀) con el nodo antiguo 2.506e-2 (3.107e59). Con 2.633e-2 sale **3.26e59** (`verificacion_numeros.md` fila 2), y en el texto se escribe 3.3×10⁵⁹.
   - El **1.44** de "∼1.6 órdenes (1.44 respecto a ρ_Λ)" es log₁₀(Ω_Λ0/2.5e-2). Con 2.633e-2 sale **1.42** (fila 13).
   - El 3.84e-2 de la instrucción sí corresponde a 2.633e-2.
   - Si el autor prefiere los valores pedidos, son dos cambios de una línea.
2. **Propagación no listada en la instrucción: H.3–H.5 y los casos 1 y 2 del apéndice H** (7 líneas). Poner 1/Ω_Λ0 en H.2 (l.833), como se pidió, sin llevarlo a H.3 y a la cota H.5 dejaba el apéndice incoherente. La cota pasa de 1.0e-4 a 7.0e-5, que es el valor "con σ_X/Ω_Λ0 en H.2" de `verificacion_numeros.md` filas 10–11. La conclusión de l.906 ("|ζ| ≲ 10⁻⁴") sigue valiendo en orden de magnitud y no se toca.
3. **El slip cambia (×1.460) y eso afecta a la lectura de C7.** `slip_bridge.py` usaba π_T = εσ_X ρ_X con ρ_X/ρ_m = Ω_Λ0/(Ω_m(1+z)³), o sea ρ_X = Ω_Λ0 ρ_crit,0 (lectura D3), no ρ_crit,0. Con [X-std], π_T = εσ_X ρ_crit,0 y el slip sube un factor 1/0.685.
   - En el mejor caso (ε = 1, δ_m = 1, z = 0.5) queda en 4.95e-2: a 1.01× del suelo de forecast 0.05. Antes estaba 1.5× por debajo.
   - Se ajustaron el resumen y la tabla del PREPRINT, C7, las tablas de `anisotropic-slip-option0.md` y `resume.txt`.
   - En anisotropic-slip l.158 "**still above** by ∼0.9–1.5×" pasa a "**comparable**: ∼0.6–1.0×".
   - No se tocó la frase interpretativa de anisotropic-slip l.160 ("nor automatically at Euclid"), que sigue siendo cierta salvo en ese mejor caso. Conviene que el autor la revise.

## Cambios

### a. Definición única de X ([X-std]), σ_X² ≡ σ²/(2θ) y factor S_iS_j (tabla 1)

**`manuscript/PREPRINT.md` l.82** — tabla 1 (PREPRINT l.82): Ω_X no está definido; definición → [X-std]

```diff
- Let $X(x)\equiv\delta\Omega_X(x)$ with $x=\ln a$. The stationary OU process is
+ Let $X(x)\equiv\delta\rho_{\mathrm{DE}}(x)/\rho_{\mathrm{crit},0}$ ($=\delta\Omega_\Lambda$ in today's units) with $x=\ln a$, so that $\delta\ln\rho_{\mathrm{DE}}=X/\Omega_{\Lambda0}$. The stationary OU process is
```

**`manuscript/PREPRINT.md` l.86** — tabla 1 (PREPRINT l.86–94): σ_X² ≡ Var X = σ²/(2θ)

```diff
- \mathrm{d}X=-\theta X\,\mathrm{d}x+\sigma\,\mathrm{d}W_x,\quad \mathrm{Var}(X)=\frac{\sigma^{2}}{2\theta}\quad(\theta>0).
+ \mathrm{d}X=-\theta X\,\mathrm{d}x+\sigma\,\mathrm{d}W_x,\quad \sigma_X^2\equiv\mathrm{Var}(X)=\frac{\sigma^{2}}{2\theta}\quad(\theta>0).
```

**`manuscript/PREPRINT.md` l.90** — tabla 1 (PREPRINT l.86–94): σ_X enlazado a σ y definición de S_i para el factor S_iS_j

```diff
- I work with an effective residual amplitude $\sigma_X$ entering the BAO residual covariance. For bins $i,j$ with logarithmic separation $\Delta x_{ij}$,
+ I work with the effective residual amplitude $\sigma_X$ (the stationary standard deviation of $X$, not the diffusion $\sigma$), which enters the BAO residual covariance through the sensitivity kernel $S_i\equiv S(z_i)=\partial\ln D_V(z_i)/\partial\Omega_\Lambda$ (flat), $\delta\alpha_i=S_iX(z_i)$. For bins $i,j$ with logarithmic separation $\Delta x_{ij}$,
```

**`manuscript/PREPRINT.md` l.94** — tabla 1 (PREPRINT l.86–94): falta S_iS_j (scripts/desi_dr2_data.py::add_ou_kernel)

```diff
- C^{\mathrm{OU}}_{ij}=\sigma_X^2 e^{-\theta\Delta x_{ij}},\quad C^{\mathrm{QNM}}_{ij}=\sigma_X^2 e^{-\theta\Delta x_{ij}}\cos(\omega_R\Delta x_{ij}).
+ C^{\mathrm{OU}}_{ij}=S_iS_j\,\sigma_X^2 e^{-\theta\Delta x_{ij}},\quad C^{\mathrm{QNM}}_{ij}=S_iS_j\,\sigma_X^2 e^{-\theta\Delta x_{ij}}\cos(\omega_R\Delta x_{ij}).
```

**`papers/stochastic-dark-energy-desi-dr2.md` l.70** — tabla 1 (l.70): definición → [X-std]

```diff
- Define $X(x) \equiv \delta\Omega_\Lambda(x)$, where $x = \ln a$ is the logarithmic scale factor. I model $X$ as an Ornstein–Uhlenbeck process:
+ Define $X(x) \equiv \delta\rho_{\rm DE}(x)/\rho_{{\rm crit},0}$ ($= \delta\Omega_\Lambda$ in today's units; $\delta\ln\rho_{\rm DE} = X/\Omega_{\Lambda0}$), where $x = \ln a$ is the logarithmic scale factor. I model $X$ as an Ornstein–Uhlenbeck process:
```

**`papers/stochastic-dark-energy-desi-dr2.md` l.81** — tabla 1 (l.82): añadir σ_X² ≡ σ²/(2θ); σ_X es la desviación típica estacionaria, no la difusión σ

```diff
- $$
- \mathrm{Var}(X) = \frac{\sigma^2}{2\theta}.
+ $$
+ \sigma_X^2 \equiv \mathrm{Var}(X) = \frac{\sigma^2}{2\theta} \qquad (\sigma_X:\ \text{stationary standard deviation, not the diffusion } \sigma).
```

**`papers/stochastic-dark-energy-desi-dr2.md` l.176** — tabla 1 (l.176): D3 → [X-std]; y el cambio de variable que da un solo H (fila 14)

```diff
- Converting to the logarithmic scale factor $x = \ln a$ (using $dt = dx/H$ and defining $X = \delta\rho_\Lambda / \rho_{\Lambda,0}$), this becomes exactly Axiom A3. The mean-reversion rate $\theta$ is the effective damping of fluctuations by the Hubble friction term $3H(1+w_\Lambda)$; the diffusion amplitude $\sigma$ encodes the strength of the microscopic noise source from Axiom A2.
+ Converting to the logarithmic scale factor $x = \ln a$ (using $dt = dx/H$, so that $dW_x = \sqrt{H}\,dW_t$, with $X = \delta\rho_{\rm DE} / \rho_{{\rm crit},0}$), this becomes exactly Axiom A3. The mean-reversion rate $\theta$ is the effective damping of fluctuations by the Hubble friction term $3H(1+w_\Lambda)$; the diffusion amplitude $\sigma$ encodes the strength of the microscopic noise source from Axiom A2.
```

**`papers/stochastic-dark-energy-desi-dr2.md` l.181** — tabla 1 (l.182): σ_X² ≡ σ²/(2θ)

```diff
- $$
- \mathrm{Var}(X) = \frac{\sigma^2}{2\theta},
+ $$
+ \sigma_X^2 \equiv \mathrm{Var}(X) = \frac{\sigma^2}{2\theta},
```

**`papers/stochastic-dark-energy-desi-dr2.md` l.272** — tabla 1 (l.264/272): σ_X² ≡ σ²/(2θ) en la covarianza (el factor S(z_i)S(z_j) ya está en l.264)

```diff
- \mathrm{Cov}[X(x_i), X(x_j)] = \frac{\sigma^2}{2\theta} \exp\left[-\theta \lvert x_i - x_j \rvert\right],
+ \mathrm{Cov}[X(x_i), X(x_j)] = \sigma_X^2 \exp\left[-\theta \lvert x_i - x_j \rvert\right], \qquad \sigma_X^2 \equiv \frac{\sigma^2}{2\theta},
```

**`papers/stochastic-dark-energy-desi-dr2.md` l.566** — tabla 1 (l.562, D2): X(z) ≡ δρ_DE/ρ_crit,0 y, en lectura (a), X(z) − X(0) para que E(0) = 1

```diff
- one obtains, to first order in $X$:
+ with $X(z) \equiv \delta\rho_{\rm DE}(z)/\rho_{{\rm crit},0}$ (in reading (a), $X(z) - X(0)$, so that $E(0) = 1$), one obtains, to first order in $X$:
```

**`papers/stochastic-dark-energy-desi-dr2.md` l.653** — tabla 1 (l.653): D3 y 'fractional' → [X-std]

```diff
- Define the fractional fluctuation field $X(x, \mathbf{r}) \equiv \delta\rho_\Lambda(x, \mathbf{r}) / \rho_{\Lambda,0}$, where $x = \ln a$ and $\mathbf{r}$ is the comoving position vector. The two-point correlation function is:
+ Define the fluctuation field $X(x, \mathbf{r}) \equiv \delta\rho_{\rm DE}(x, \mathbf{r}) / \rho_{{\rm crit},0}$ (as in Axiom A3), where $x = \ln a$ and $\mathbf{r}$ is the comoving position vector. The two-point correlation function is:
```

**`theory/unimodular-gravity-vacuum-smoothness.md` l.64** — tabla 1 (unimodular l.68): 'fractional' y D3 → [X-std]

```diff
- I define the fractional fluctuation of the dark-energy density as:
+ I define the dark-energy density fluctuation in units of today's critical density as:
```

**`theory/unimodular-gravity-vacuum-smoothness.md` l.68** — tabla 1 (unimodular l.68): X ≡ δρ_DE/ρ_crit,0 = Ω_Λ0·δρ_Λ/ρ̄_Λ

```diff
- X(x) \equiv \frac{\delta\rho_\Lambda(x)}{\bar{\rho}_\Lambda},
+ X(x) \equiv \frac{\delta\rho_{\rm DE}(x)}{\rho_{{\rm crit},0}} = \Omega_{\Lambda 0}\,\frac{\delta\rho_\Lambda(x)}{\bar{\rho}_\Lambda},
```

**`theory/unimodular-gravity-vacuum-smoothness.md` l.72** — tabla 1 (unimodular l.72): σ_X² ≡ ⟨X²⟩ = σ²/2θ (⇔ rms δρ_Λ/ρ̄_Λ = σ_X/Ω_Λ0 < 3.84e-2)

```diff
- with variance $\sigma_X^2 \equiv \langle X^2 \rangle$. The DESI DR2 BAO data, combined with the sensitivity kernel $S(z) = \partial \ln D_V / \partial \Omega_\Lambda$, yield [1]:
+ with variance $\sigma_X^2 \equiv \langle X^2 \rangle$ ($= \sigma^2/2\theta$ for the OU process; equivalently, rms $\delta\rho_\Lambda/\bar{\rho}_\Lambda = \sigma_X/\Omega_{\Lambda 0}$, below $3.84\times 10^{-2}$ for the bound that follows). The DESI DR2 BAO data, combined with the sensitivity kernel $S(z) = \partial \ln D_V / \partial \Omega_\Lambda$, yield [1]:
```

**`theory/fundamental-vs-emergent-vacuum-relaxation.md` l.32** — tabla 1 (fundamental-vs-emergent l.32): quitar 'fractional' (sugiere D3)

```diff
- I take fractional fluctuations $X\equiv\delta\Omega_\Lambda$ in $x=\ln a$:
+ I take fluctuations $X\equiv\delta\Omega_\Lambda=\delta\rho_{\rm DE}/\rho_{{\rm crit},0}$ in $x=\ln a$:
```

**`theory/principle-of-vacuum-smoothness.md` l.15** — tabla 1 (principle l.15–19): σ_X = desviación típica de δρ_DE/ρ_crit,0 (⇔ rms δρ_Λ/ρ_Λ < 3.84e-2)

```diff
- Any physical mechanism that could generate fluctuations in the dark energy density $\rho_\Lambda$ — whether from spacetime discreteness, horizon quasi-normal modes, scalar field dynamics, or quantum gravitational effects — must be naturally suppressed below what current surveys can detect:
+ Any physical mechanism that could generate fluctuations in the dark energy density $\rho_\Lambda$ — whether from spacetime discreteness, horizon quasi-normal modes, scalar field dynamics, or quantum gravitational effects — must be naturally suppressed below what current surveys can detect. With $\sigma_X$ the standard deviation of $X\equiv\delta\rho_{\rm DE}/\rho_{{\rm crit},0}$ (equivalently, rms $\delta\rho_\Lambda/\rho_\Lambda=\sigma_X/\Omega_{\Lambda 0}<3.84\times 10^{-2}$):
```

**`papers/quantum-fluid-instabilities-desi-dr2.md` l.103** — tabla 1 (quantum-fluid l.103): D3 frente a la verosimilitud D1 (u_i = σ_0 S_i e^{t/tc}); → δρ_DE/ρ_crit,0

```diff
- \sigma_{X}^{2}(t)\equiv\frac{\langle(\delta\rho_{\Lambda})^{2}\rangle}{\rho_{\Lambda}^{2}}=\frac{1}{\rho_{\Lambda}^{2}}\int\frac{\mathrm{d}^{3}q}{(2\pi)^{3}}P_{0}(q)\,e^{2\Gamma_{q}t}.
+ \sigma_{X}^{2}(t)\equiv\frac{\langle(\delta\rho_{\mathrm{DE}})^{2}\rangle}{\rho_{\mathrm{crit},0}^{2}}=\frac{1}{\rho_{\mathrm{crit},0}^{2}}\int\frac{\mathrm{d}^{3}q}{(2\pi)^{3}}P_{0}(q)\,e^{2\Gamma_{q}t}.
```

**`exploratory/notes/desqueezing-relaxation-vacuum-fluctuations-note.md` l.47** — tabla 1 (nota desqueezing l.47): D3 → [X-std] y quitar 'or equivalent' (no son equivalentes, ×1.46)

```diff
- where $X\equiv\delta\Omega_\Lambda/\bar{\Omega}_\Lambda$ (or an equivalent fractional density contrast in the notation of the companion papers). This is not merely a non-detection: it is a structural constraint on how “noisy” the late-time vacuum is allowed to be at BAO scales.
+ where $X\equiv\delta\rho_{\rm DE}/\rho_{{\rm crit},0}$ ($=\delta\Omega_\Lambda$ in today's units; $\delta\ln\rho_{\rm DE}=X/\Omega_{\Lambda0}$) and $\sigma_X^2\equiv\mathrm{Var}\,X=\sigma^2/(2\theta)$. This is not merely a non-detection: it is a structural constraint on how “noisy” the late-time vacuum is allowed to be at BAO scales.
```

**`exploratory/notes/desqueezing-relaxation-vacuum-fluctuations-note.md` l.165** — tabla 1 (nota l.165): 'fractional fluctuations' → X ≡ δρ_DE/ρ_crit,0

```diff
- - **A3 (OU closure):** fractional fluctuations $X(x)$ in logarithmic scale factor $x=\ln a$ obey
+ - **A3 (OU closure):** fluctuations $X(x)\equiv\delta\rho_{\rm DE}/\rho_{{\rm crit},0}$ in logarithmic scale factor $x=\ln a$ obey
```

**`exploratory/notes/desqueezing-relaxation-vacuum-fluctuations-note.md` l.173** — tabla 1 (nota l.173): añadir σ_X² ≡ σ²/(2θ)

```diff
- with stationary variance $\mathrm{Var}(X)=\sigma^2/(2\theta)$ when $\theta>0$.
+ with stationary variance $\sigma_X^2\equiv\mathrm{Var}(X)=\sigma^2/(2\theta)$ when $\theta>0$ ($\sigma_X$ is the stationary standard deviation, not the diffusion $\sigma$).
```

**`exploratory/notes/desqueezing-relaxation-vacuum-fluctuations-note.md` l.209** — tabla 1 (nota l.209): σ_X² en la covarianza (S(z_i)S(z_j) ya está en l.200)

```diff
- \mathrm{Cov}[X(x_i),X(x_j)]=\frac{\sigma^2}{2\theta} e^{-\theta\lvert x_i-x_j\rvert},
+ \mathrm{Cov}[X(x_i),X(x_j)]=\sigma_X^2 e^{-\theta\lvert x_i-x_j\rvert},\qquad \sigma_X^2\equiv\frac{\sigma^2}{2\theta},
```

### a. Scripts: solo nombres, comentarios y docstrings

**`scripts/desi_dr2_data.py` l.225** — tabla 1 (desi_dr2_data l.224–233): es la fuente de verdad y no definía X

```diff
-     """C_total = C_meas + OU residual signal."""
+     """C_total = C_meas + OU residual signal: C_ij += S_i S_j σ_X² exp(−θ |x_i − x_j|).
+ 
+     X(x) ≡ δρ_DE(x)/ρ_crit,0 (= δΩ_Λ in today's units); δ ln ρ_DE = X/Ω_Λ0;
+     σ_X² ≡ Var X = σ²/(2θ): sigma_X is the stationary standard deviation, not the diffusion σ.
+     S_z = ∂ln D_V/∂Ω_Λ (flat), so δα_i = S_i X(z_i).
+     """
```

**`scripts/ou_bao_likelihood.py` l.87** — tabla 1 (ou_bao_likelihood l.86): el nombre invitaba a confundirlo con la difusión σ

```diff
- SIGMA_OU   = 2.31e-2    # sigma_X (calibrated from DESI DR1 anchor bin)
- SIGMA_X2   = SIGMA_OU**2
+ SIGMA_X_DR1 = 2.31e-2   # stationary std of X (sigma_X, not the diffusion sigma); obsolete DR1 calibration, used only as a starting value
+ SIGMA_X2   = SIGMA_X_DR1**2
```

**`scripts/eos_efectiva.py` l.135** — tabla 1 (eos_efectiva l.129/150–156): es la amplitud de un toy de δw, no σ_X; se renombra (las claves 'sigma' de la salida no cambian, para no tocar results/)

```diff
- # w_eff(z) = w_CPL(z) + δw(z),  δw ~ σ²/(2θ) (1 - e^{-2θ t}) / Ω_DE  (toy)
+ # w_eff(z) = w_CPL(z) + δw(z),  δw ~ σ_w²/(2θ) (1 - e^{-2θ t}) / Ω_DE  (toy)
+ # sigma_w is the amplitude of this toy δw correction; it is NOT σ_X (the BAO residual amplitude).
+ # The output key "sigma" (results/eos_cpl_desi_dr2/) is sigma_w; kept for compatibility.
```

**`scripts/eos_efectiva.py` l.147, 152, 160, 162, 172, 173, 176, 177, 183, 186, 194, 212, 214, 228, 229, 280** — renombre `sigma` → `sigma_w` en las funciones del toy δw (16 líneas; las claves `"sigma"` de la salida no cambian).

```diff
- def Omega_DE_eff(z: float, w0: float, wa: float, sigma: float, theta: float) -> float:
+ def Omega_DE_eff(z: float, w0: float, wa: float, sigma_w: float, theta: float) -> float:
-         if sigma <= 0.0:
+         if sigma_w <= 0.0:
-                 corr = (sigma**2) * tz / max(omega_bg, 1e-12)
+                 corr = (sigma_w**2) * tz / max(omega_bg, 1e-12)
-                 corr = (sigma**2 / (2.0 * theta * max(omega_bg, 1e-12))) * (
+                 corr = (sigma_w**2 / (2.0 * theta * max(omega_bg, 1e-12))) * (
- def E2_eff(z: float, w0: float, wa: float, sigma: float, theta: float) -> float:
+ def E2_eff(z: float, w0: float, wa: float, sigma_w: float, theta: float) -> float:
-     return Om * (1.0 + z) ** 3 + Omega_DE_eff(z, w0, wa, sigma, theta)
+     return Om * (1.0 + z) ** 3 + Omega_DE_eff(z, w0, wa, sigma_w, theta)
- def H_eff(z: float, w0: float, wa: float, sigma: float, theta: float) -> float:
+ def H_eff(z: float, w0: float, wa: float, sigma_w: float, theta: float) -> float:
-     e2 = E2_eff(z, w0, wa, sigma, theta)
+     e2 = E2_eff(z, w0, wa, sigma_w, theta)
- def alpha_pred_eff(w0: float, wa: float, sigma: float, theta: float) -> np.ndarray:
+ def alpha_pred_eff(w0: float, wa: float, sigma_w: float, theta: float) -> np.ndarray:
-             h = H_eff(zp, w0, wa, sigma, theta)
+             h = H_eff(zp, w0, wa, sigma_w, theta)
-         h = H_eff(z, w0, wa, sigma, theta)
+         h = H_eff(z, w0, wa, sigma_w, theta)
- def chi2_eff(w0: float, wa: float, sigma: float, theta: float) -> float:
+ def chi2_eff(w0: float, wa: float, sigma_w: float, theta: float) -> float:
-         pred = alpha_pred_eff(w0, wa, sigma, theta)
+         pred = alpha_pred_eff(w0, wa, sigma_w, theta)
- def logL_eff(w0: float, wa: float, sigma: float, theta: float) -> float:
+ def logL_eff(w0: float, wa: float, sigma_w: float, theta: float) -> float:
-     c2 = chi2_eff(w0, wa, sigma, theta)
+     c2 = chi2_eff(w0, wa, sigma_w, theta)
-     """Joint {w0, wa, sigma, theta}; sigma>=0, theta>0."""
+     """Joint {w0, wa, sigma_w, theta}; sigma_w>=0, theta>0 (output key "sigma")."""
```

### b. Normalización del ruido (fórmula deducida, fila 14)

**`papers/stochastic-dark-energy-desi-dr2.md` l.172** — fórmula deducida (verificacion_numeros fila 14), validada por el autor: ⟨ξξ⟩ = σ²H ρ_crit,0² δ_D

```diff
- \langle \xi(t) \rangle = 0, \qquad \langle \xi(t) \xi(t') \rangle = \frac{\sigma^2 H^2 \rho_{\Lambda,0}^2}{2\theta} \cdot \delta_D(t - t').
+ \langle \xi(t) \rangle = 0, \qquad \langle \xi(t) \xi(t') \rangle = \sigma^2 H \rho_{{\rm crit},0}^2 \, \delta_D(t - t').
```

### b. Órdenes de magnitud de supresión (fórmula deducida, fila 13)

**`papers/stochastic-dark-energy-desi-dr2.md` l.372** — fórmula deducida (verificacion_numeros fila 13), validada: log₁₀(1/σ_X) = 1.58; log₁₀(Ω_Λ0/σ_X) = 1.42 con 2.633e-2 (el autor escribió 1.44, que es el valor con 2.5e-2)

```diff
- The upper limit $\sigma_X < 2.6\times 10^{-2}$ implies that the dark energy density is not only fine-tuned to a specific value, but that this value is realized with extraordinary homogeneity across the observable universe. Whatever cancellation mechanism is responsible for the smallness of $\Lambda$ (whether anthropic selection, a symmetry, or a dynamical attractor), it must operate with such precision that it suppresses local fluctuations by more than four orders of magnitude relative to the background density.
+ The upper limit $\sigma_X < 2.6\times 10^{-2}$ implies that the dark energy density is not only fine-tuned to a specific value, but that this value is realized with extraordinary homogeneity across the observable universe. Whatever cancellation mechanism is responsible for the smallness of $\Lambda$ (whether anthropic selection, a symmetry, or a dynamical attractor), it must operate with such precision that it suppresses local fluctuations by $\sim 1.6$ orders of magnitude relative to the critical density ($\log_{10}(1/\sigma_X) = 1.58$), or $1.42$ relative to $\rho_\Lambda$ ($\log_{10}(\Omega_{\Lambda0}/\sigma_X)$).
```

### b. Factor 1/Ω_Λ0: rms δρ_Λ/ρ_Λ = σ_X/Ω_Λ0 (§7 y apéndice H)

**`papers/stochastic-dark-energy-desi-dr2.md` l.400** — tabla 1 (l.404): 'fractional' sugiere D3; con [X-std] la fluctuación relativa de ρ_Λ es σ_X/Ω_Λ0

```diff
- where $\beta$ is the effective coupling constant. Assuming $\rho_\Lambda \approx V(\phi)$, the fractional fluctuation in $\rho_\Lambda$ is:
+ where $\beta$ is the effective coupling constant. Assuming $\rho_\Lambda \approx V(\phi)$, the rms relative fluctuation of $\rho_\Lambda$ is, with $\delta\ln\rho_{\rm DE} = X/\Omega_{\Lambda0}$:
```

**`papers/stochastic-dark-energy-desi-dr2.md` l.404** — tabla 1 (l.404): rms(δρ_Λ)/ρ_Λ = σ_X/Ω_Λ0 ≈ (V'/V) rms(δφ)

```diff
- \sigma_X \equiv \frac{\delta \rho_\Lambda}{\rho_\Lambda} \approx \frac{V'(\phi)}{V(\phi)} \delta \phi.
+ \frac{\mathrm{rms}(\delta \rho_\Lambda)}{\rho_\Lambda} = \frac{\sigma_X}{\Omega_{\Lambda0}} \approx \frac{V'(\phi)}{V(\phi)}\, \mathrm{rms}(\delta \phi).
```

**`papers/stochastic-dark-energy-desi-dr2.md` l.412** — propagar l.404 a l.412

```diff
- \frac{\Delta \alpha}{\alpha} \approx \beta \frac{V(\phi)}{V'(\phi)} \sigma_X.
+ \frac{\Delta \alpha}{\alpha} \approx \beta \frac{V(\phi)}{V'(\phi)} \frac{\sigma_X}{\Omega_{\Lambda0}}.
```

**`papers/stochastic-dark-energy-desi-dr2.md` l.416** — propagar a l.416: σ_X/Ω_Λ0 = 2.633e-2/0.685 = 3.84e-2

```diff
- My upper bound $\sigma_X < 2.6\times 10^{-2}$ therefore implies:
+ My upper bound $\sigma_X < 2.6\times 10^{-2}$, i.e. $\sigma_X/\Omega_{\Lambda0} < 3.84\times 10^{-2}$, therefore implies:
```

**`papers/stochastic-dark-energy-desi-dr2.md` l.420** — propagar a l.420

```diff
- \left| \frac{\Delta \alpha}{\alpha}\right| \lesssim \lvert \beta \rvert \left| \frac{V}{V'}\right| (2.6\times 10^{-2}).
+ \left| \frac{\Delta \alpha}{\alpha}\right| \lesssim \lvert \beta \rvert \left| \frac{V}{V'}\right| (3.84\times 10^{-2}).
```

**`papers/stochastic-dark-energy-desi-dr2.md` l.821** — tabla 1 (l.825): el texto principal define X ≡ δρ_DE/ρ_crit,0, no δρ_Λ/ρ_Λ

```diff
- where $V'(\phi) = dV/d\phi$. The dimensionless stochastic amplitude defined in the main text is:
+ where $V'(\phi) = dV/d\phi$. With the definition of the main text, $X \equiv \delta\rho_{\rm DE}/\rho_{{\rm crit},0}$, the rms relative fluctuation of $\rho_\Lambda$ is:
```

**`papers/stochastic-dark-energy-desi-dr2.md` l.825** — tabla 1 (l.825): rms(δρ_Λ)/ρ_Λ = σ_X/Ω_Λ0

```diff
- \sigma_X \equiv \frac{\delta\rho_\Lambda}{\rho_\Lambda}.
+ \frac{\mathrm{rms}(\delta\rho_\Lambda)}{\rho_\Lambda} = \frac{\sigma_X}{\Omega_{\Lambda 0}}.
```

**`papers/stochastic-dark-energy-desi-dr2.md` l.833** — tabla 1 (l.833): en H.2, δφ ≃ (σ_X/Ω_Λ0)V/V'

```diff
- \sigma_X \simeq \frac{V'(\phi)}{V(\phi)} \delta\phi \quad \Longrightarrow \quad \delta\phi \simeq \sigma_X \frac{V(\phi)}{V'(\phi)}. \tag{H.2}
+ \frac{\sigma_X}{\Omega_{\Lambda 0}} \simeq \frac{V'(\phi)}{V(\phi)} \delta\phi \quad \Longrightarrow \quad \delta\phi \simeq \frac{\sigma_X}{\Omega_{\Lambda 0}} \frac{V(\phi)}{V'(\phi)}. \tag{H.2}
```

### b. Propagación necesaria de H.2 a H.3–H.5 (no estaba en la lista)

**`papers/stochastic-dark-energy-desi-dr2.md` l.841** — propagación necesaria de H.2 a H.3 (no estaba en la lista; sin ella el apéndice es incoherente)

```diff
- \boxed{ \frac{\Delta\alpha}{\alpha} \simeq \zeta \left( \frac{V(\phi)}{M_{\rm Pl} V'(\phi)} \right) \sigma_X }. \tag{H.3}
+ \boxed{ \frac{\Delta\alpha}{\alpha} \simeq \zeta \left( \frac{V(\phi)}{M_{\rm Pl} V'(\phi)} \right) \frac{\sigma_X}{\Omega_{\Lambda 0}} }. \tag{H.3}
```

**`papers/stochastic-dark-energy-desi-dr2.md` l.865** — propagación: la cota que entra en H.5 es σ_X/Ω_Λ0 < 3.84e-2

```diff
- Combining this with the upper bound on the stochastic amplitude (Sec. 4.2, $\sigma_X < 2.6\times 10^{-2}$ at 95% CL) and using Eq. (H.3), I obtain a combined constraint on the product of the coupling and the potential slope:
+ Combining this with the upper bound on the stochastic amplitude (Sec. 4.2, $\sigma_X < 2.6\times 10^{-2}$ at 95% CL, i.e. $\sigma_X/\Omega_{\Lambda 0} < 3.84\times 10^{-2}$) and using Eq. (H.3), I obtain a combined constraint on the product of the coupling and the potential slope:
```

**`papers/stochastic-dark-energy-desi-dr2.md` l.869** — propagación: 2.7e-6/3.84e-2 = 7.0e-5 (verificacion_numeros fila 10, 'con σ_X/Ω_Λ0 en H.2')

```diff
- \left| \zeta \left( \frac{V(\phi)}{M_{\rm Pl} V'(\phi)} \right) \right| \lesssim \frac{2.7 \times 10^{-6}}{2.633\times 10^{-2}} \approx 1.0 \times 10^{-4}. \tag{H.5}
+ \left| \zeta \left( \frac{V(\phi)}{M_{\rm Pl} V'(\phi)} \right) \right| \lesssim \frac{2.7 \times 10^{-6}}{3.84\times 10^{-2}} \approx 7.0 \times 10^{-5}. \tag{H.5}
```

**`papers/stochastic-dark-energy-desi-dr2.md` l.883** — propagación de H.5 (caso 1)

```diff
- \left| \frac{V(\phi)}{M_{\rm Pl} V'(\phi)} \right| \lesssim 1.0 \times 10^{-4}.
+ \left| \frac{V(\phi)}{M_{\rm Pl} V'(\phi)} \right| \lesssim 7.0 \times 10^{-5}.
```

**`papers/stochastic-dark-energy-desi-dr2.md` l.887** — propagación: λ ≳ 1/7.0e-5 = 1.4e4 (fila 11, 'con D3')

```diff
- For an exponential potential $V(\phi) = V_0 \exp(-\lambda \phi / M_{\rm Pl})$, the above ratio is $1/\lambda$, hence $\lambda \gtrsim 9.8\times 10^{3}$. However, a scalar field with $\lambda > \sqrt{24\pi} \approx 8.7$ (or even $\lambda > 3$ in the slow-roll limit) cannot produce acceleration at the present epoch; instead it enters a tracking solution where the field energy density scales as matter or radiation, and the equation of state deviates significantly from $-1$ [2]. This scenario is therefore incompatible with DESI DR2 observations that favour a background close to $\Lambda$CDM.
+ For an exponential potential $V(\phi) = V_0 \exp(-\lambda \phi / M_{\rm Pl})$, the above ratio is $1/\lambda$, hence $\lambda \gtrsim 1.4\times 10^{4}$. However, a scalar field with $\lambda > \sqrt{24\pi} \approx 8.7$ (or even $\lambda > 3$ in the slow-roll limit) cannot produce acceleration at the present epoch; instead it enters a tracking solution where the field energy density scales as matter or radiation, and the equation of state deviates significantly from $-1$ [2]. This scenario is therefore incompatible with DESI DR2 observations that favour a background close to $\Lambda$CDM.
```

**`papers/stochastic-dark-energy-desi-dr2.md` l.895** — propagación de H.5 (caso 2)

```diff
- \lvert \zeta \rvert \lesssim 1.0 \times 10^{-4}.
+ \lvert \zeta \rvert \lesssim 7.0 \times 10^{-5}.
```

**`papers/stochastic-dark-energy-desi-dr2.md` l.899** — propagación: M_cut = M_Pl/|ζ| ≳ 1.4e4 M_Pl

```diff
- This implies a suppression scale $M_{\rm cut} \equiv M_{\rm Pl} / \lvert \zeta \rvert \gtrsim 9.8\times 10^{3} M_{\rm Pl}$, placing the new-physics scale well above the Planck scale. Such a suppression is technically unnatural and requires a protective symmetry mechanism, such as an axionic shift symmetry, to suppress quantum corrections from charged-particle loops [3].
+ This implies a suppression scale $M_{\rm cut} \equiv M_{\rm Pl} / \lvert \zeta \rvert \gtrsim 1.4\times 10^{4} M_{\rm Pl}$, placing the new-physics scale well above the Planck scale. Such a suppression is technically unnatural and requires a protective symmetry mechanism, such as an axionic shift symmetry, to suppress quantum corrections from charged-particle loops [3].
```

### b. Gap en densidad (fórmula deducida, fila 2)

**`papers/amplification-gap.md` l.54** — fórmula deducida (verificacion_numeros fila 2), validada: |δρ|max = σ_X ρ_crit,0 = (σ_X/Ω_Λ0)ρ_Λ; σ_X/(Ω_Λ0 σ_0) = 3.26e59 con 2.633e-2 (el autor escribió 3.1e59, que es el valor con 2.506e-2)

```diff
- - Density form of the same ratio: $\lvert\delta\rho\rvert_{\max}/(\sigma_0\rho_\Lambda)=\sigma_{\mathrm{tgt}}/\sigma_0$ — so the DESI density gap is $\sim 10^{59}$, not $10^{56}$.
+ - Density form: with $X\equiv\delta\rho_{\rm DE}/\rho_{{\rm crit},0}$, $\lvert\delta\rho\rvert_{\max}=\sigma_X\rho_{{\rm crit},0}=(\sigma_X/\Omega_{\Lambda0})\rho_\Lambda$, so $\lvert\delta\rho\rvert_{\max}/(\sigma_0\rho_\Lambda)=\sigma_X/(\Omega_{\Lambda0}\sigma_0)\approx 3.3\times 10^{59}$ for the DESI ceiling — still $\sim 10^{59}$, not $10^{56}$.
```

**`papers/amplification-gap.md` l.56** — tabla 1 (amplification-gap l.56): 'fractional residual amplitude' (D3 implícita)

```diff
- **What DESI actually measured** is a **null residual bound**, not a new $\rho_\Lambda$. Mean vacuum density $\rho_\Lambda\sim\Omega_{\Lambda 0}\rho_{\mathrm{crit}}$ comes from standard cosmology (e.g. Planck-class $H_0,\Omega_\Lambda$); this corpus constrains the **fractional residual amplitude** $\sigma_X$ on top of that background.
+ **What DESI actually measured** is a **null residual bound**, not a new $\rho_\Lambda$. Mean vacuum density $\rho_\Lambda\sim\Omega_{\Lambda 0}\rho_{\mathrm{crit}}$ comes from standard cosmology (e.g. Planck-class $H_0,\Omega_\Lambda$); this corpus constrains the **residual amplitude** $\sigma_X$ (in units of $\rho_{{\rm crit},0}$) on top of that background.
```

### b. slip_bridge.py: π_T = εσ_X ρ_crit,0

**`scripts/slip_bridge.py` l.20** — [X-std]: X ≡ δρ_DE/ρ_crit,0, so the shear is ε σ_X ρ_crit,0 (fórmula deducida b, validada por el autor)

```diff
-   π_T = eps · σ_X · ρ_X
-   |γ − 1| = 2 · eps · σ_X · (ρ_X/ρ_m) / |δ_m|
+   π_T = eps · σ_X · ρ_crit,0          (X ≡ δρ_DE/ρ_crit,0, so rms δρ_DE = σ_X ρ_crit,0)
+   |γ − 1| = 2 · eps · σ_X · (ρ_crit,0/ρ_m) / |δ_m|
```

**`scripts/slip_bridge.py` l.46** — Ω_Λ0 deja de entrar en el slip con [X-std]

```diff
- OM_L0 = 0.685
+ OM_L0 = 0.685  # not used by the slip map: π_T is σ_X ρ_crit,0, not σ_X ρ_Λ (block C)
```

**`scripts/slip_bridge.py` l.65** — [X-std]: X ≡ δρ_DE/ρ_crit,0, so the shear is ε σ_X ρ_crit,0 (fórmula deducida b, validada por el autor)

```diff
- def rho_ratio_x_over_m(z: float) -> float:
-     """ρ_DE / ρ_m at redshift z (DE ~ cosmological constant)."""
-     return OM_L0 / (OM_M0 * (1.0 + z) ** 3)
+ def rho_crit0_over_m(z: float) -> float:
+     """ρ_crit,0 / ρ_m(z) = 1 / [Ω_m0 (1+z)^3]; σ_X is in units of ρ_crit,0 ([X-std])."""
+     return 1.0 / (OM_M0 * (1.0 + z) ** 3)
```

**`scripts/slip_bridge.py` l.77** — renombre de la función

```diff
-     return 2.0 * eps * sigma_x * rho_ratio_x_over_m(z) / abs(delta_m)
+     return 2.0 * eps * sigma_x * rho_crit0_over_m(z) / abs(delta_m)
```

**`scripts/slip_bridge.py` l.82** — renombre de la función

```diff
-     denom = 2.0 * sigma_x * rho_ratio_x_over_m(z)
+     denom = 2.0 * sigma_x * rho_crit0_over_m(z)
```

**`scripts/slip_bridge.py` l.103** — renombre de la función y etiqueta

```diff
-         r = rho_ratio_x_over_m(z)
-         print(f"--- z = {z}   (ρ_X/ρ_m = {r:.3f}) ---")
+         r = rho_crit0_over_m(z)
+         print(f"--- z = {z}   (ρ_crit,0/ρ_m = {r:.3f}) ---")
```

### b. Slip regenerado: propagación a PREPRINT, CLAIMS y las tablas de B2

**`manuscript/PREPRINT.md` l.28** — propagación del slip con π_T = εσ_X ρ_crit,0 (slip_bridge regenerado)

```diff
- Gravitational slip $\gamma=\Phi/\Psi$ is the correct operator for anisotropic leakage past volume-preserving (SDiff) structure, but inherits the same amplitude starvation: even full anisotropy of the BAO-bounded residual gives $\lvert\gamma-1\rvert\sim 3\times 10^{-2}$, of order $1.5$–$7\times$ below current and forecast floors.
+ Gravitational slip $\gamma=\Phi/\Psi$ is the correct operator for anisotropic leakage past volume-preserving (SDiff) structure, but inherits the same amplitude starvation: even full anisotropy of the BAO-bounded residual gives at most $\lvert\gamma-1\rvert\sim 5\times 10^{-2}$ ($z=0.5$, $\delta_m=1$): $\sim 3$–$16\times$ below the current measurement and $1.0$–$4.7\times$ below the forecast floor, which it only reaches in that best case.
```

**`manuscript/PREPRINT.md` l.238** — tabla 1 (anisotropic-slip l.137/141 y slip_bridge): π_T = εσ_X ρ_crit,0

```diff
- \lvert\gamma-1\rvert=\left| \frac{\Phi-\Psi}{\Psi}\right| =2\varepsilon\sigma_X\frac{\rho_X}{\rho_m\lvert\delta_m\rvert},
+ \lvert\gamma-1\rvert=\left| \frac{\Phi-\Psi}{\Psi}\right| =2\varepsilon\sigma_X\frac{\rho_{\mathrm{crit},0}}{\rho_m\lvert\delta_m\rvert},
```

**`manuscript/PREPRINT.md` l.242** — π_T = εσ_X ρ_crit,0 explícito

```diff
- with $\varepsilon\in[0,1]$ the anisotropic fraction of the residual. With $\sigma_X=2.633\times 10^{-2}$, $\varepsilon=1$, $\delta_m=1$:
+ with $\varepsilon\in[0,1]$ the anisotropic fraction of the residual, $\pi_T=\varepsilon\sigma_X\rho_{\mathrm{crit},0}$ (since $X\equiv\delta\rho_{\mathrm{DE}}/\rho_{\mathrm{crit},0}$). With $\sigma_X=2.633\times 10^{-2}$, $\varepsilon=1$, $\delta_m=1$:
```

**`manuscript/PREPRINT.md` l.246** — slip_bridge regenerado: 3.393e-2 → 4.953e-2; 0.17/|γ−1| = 3.43; 0.05/|γ−1| = 1.01

```diff
- | 0.5 | $3.39\times 10^{-2}$ | $\sim 5.0\times$ short | $\sim 1.5\times$ short |
+ | 0.5 | $4.95\times 10^{-2}$ | $\sim 3.4\times$ short | $\sim 1.0\times$ (at the floor) |
```

**`manuscript/PREPRINT.md` l.247** — slip_bridge regenerado: 1.431e-2 → 2.090e-2; 8.13×; 2.39×

```diff
- | 1.0 | $1.43\times 10^{-2}$ | $\sim 12\times$ short | $\sim 3.5\times$ short |
+ | 1.0 | $2.09\times 10^{-2}$ | $\sim 8.1\times$ short | $\sim 2.4\times$ short |
```

**`manuscript/PREPRINT.md` l.248** — slip_bridge regenerado: 7.328e-3 → 1.070e-2; 15.9×; 4.67×

```diff
- | 1.5 | $7.33\times 10^{-3}$ | $\sim 23\times$ short | $\sim 6.8\times$ short |
+ | 1.5 | $1.07\times 10^{-2}$ | $\sim 16\times$ short | $\sim 4.7\times$ short |
```

**`manuscript/CLAIMS.md` l.22** — C7 con el slip regenerado (π_T = εσ_X ρ_crit,0)

```diff
- | **C7** | Even $\varepsilon=1$ anisotropic fraction of the BAO-bounded residual gives $\|\gamma-1\|\lesssim 3\times 10^{-2}$ ($z=0.5$, $\delta_m=1$) | `scripts/slip_bridge.py` |
+ | **C7** | Even $\varepsilon=1$ anisotropic fraction of the BAO-bounded residual gives $\|\gamma-1\|\lesssim 5\times 10^{-2}$ ($z=0.5$, $\delta_m=1$), at best level with the forecast floor $\sim 0.05$ | `scripts/slip_bridge.py` |
```

**`papers/anisotropic-slip-option0.md` l.118** — slip regenerado: escala de referencia 3e-2 → 5e-2

```diff
- With $\sigma(\gamma)\sim 0.11$, current sensitivity to $\lvert \eta-1 \rvert$ (or $\lvert \gamma-1 \rvert$) is still $\mathcal{O}(0.1)$, far above the $\sim 3\times 10^{-2}$ scale of the $f=1$, $\sigma_X=2.6\times 10^{-2}$ estimate below.
+ With $\sigma(\gamma)\sim 0.11$, current sensitivity to $\lvert \eta-1 \rvert$ (or $\lvert \gamma-1 \rvert$) is still $\mathcal{O}(0.1)$, far above the $\sim 5\times 10^{-2}$ scale of the $f=1$, $\sigma_X=2.6\times 10^{-2}$ estimate below.
```

**`papers/anisotropic-slip-option0.md` l.137** — tabla 1 (anisotropic-slip l.137): π_T = εσ_X ρ_crit,0; quitar 'fractional'

```diff
- Let $\varepsilon\in[0,1]$ (also written $f $ below) be the fraction of the DE-sector residual stress that is anisotropic (shear), and let $\sigma_X$ be the fractional residual amplitude constrained by BAO. With $\pi_T=\varepsilon \sigma_X \rho_X$ and $\delta_m\sim 1$:
+ Let $\varepsilon\in[0,1]$ (also written $f $ below) be the fraction of the DE-sector residual stress that is anisotropic (shear), and let $\sigma_X$ be the residual amplitude constrained by BAO, in units of today's critical density ($X\equiv\delta\rho_{\rm DE}/\rho_{{\rm crit},0}$). With $\pi_T=\varepsilon \sigma_X \rho_{{\rm crit},0}$ and $\delta_m\sim 1$:
```

**`papers/anisotropic-slip-option0.md` l.141** — tabla 1 (anisotropic-slip l.141)

```diff
- \lvert \gamma-1 \rvert = 2 \varepsilon \sigma_X \frac{\rho_X}{\rho_m \lvert \delta_m \rvert}.
+ \lvert \gamma-1 \rvert = 2 \varepsilon \sigma_X \frac{\rho_{{\rm crit},0}}{\rho_m \lvert \delta_m \rvert}.
```

**`papers/anisotropic-slip-option0.md` l.147** — columna ρ_X/ρ_m → ρ_crit,0/ρ_m = 1/[Ω_m(1+z)³]

```diff
- | $z$ | $\rho_X/\rho_m$ | max $\|\gamma-1\|$ | gap vs Maus $\sim 0.17$ | gap vs Sakr floor $\sim 0.05$ |
+ | $z$ | $\rho_{{\rm crit},0}/\rho_m$ | max $\|\gamma-1\|$ | gap vs Maus $\sim 0.17$ | gap vs Sakr floor $\sim 0.05$ |
```

**`papers/anisotropic-slip-option0.md` l.149** — slip_bridge regenerado (z = 0.5)

```diff
- | 0.5 | 0.644 | $3.39\times 10^{-2}$ | $\sim 5.0\times$ | $\sim 1.5\times$ |
+ | 0.5 | 0.941 | $4.95\times 10^{-2}$ | $\sim 3.4\times$ | $\sim 1.0\times$ |
```

**`papers/anisotropic-slip-option0.md` l.150** — slip_bridge regenerado (z = 1.0)

```diff
- | 1.0 | 0.272 | $1.43\times 10^{-2}$ | $\sim 12\times$ | $\sim 3.5\times$ |
+ | 1.0 | 0.397 | $2.09\times 10^{-2}$ | $\sim 8.1\times$ | $\sim 2.4\times$ |
```

**`papers/anisotropic-slip-option0.md` l.151** — slip_bridge regenerado (z = 1.5)

```diff
- | 1.5 | 0.139 | $7.33\times 10^{-3}$ | $\sim 23\times$ | $\sim 6.8\times$ |
+ | 1.5 | 0.203 | $1.07\times 10^{-2}$ | $\sim 16\times$ | $\sim 4.7\times$ |
```

**`papers/anisotropic-slip-option0.md` l.155** — slip regenerado: referencia 3e-2 → 5e-2

```diff
- | Instrument / status | Rough sensitivity to $\eta-1$ | vs $\sim 3\times 10^{-2}$ (best case) |
+ | Instrument / status | Rough sensitivity to $\eta-1$ | vs $\sim 5\times 10^{-2}$ (best case) |
```

**`papers/anisotropic-slip-option0.md` l.157** — O(0.1) frente a 4.95e-2: ∼2×

```diff
- | Current DESI-era slip (Maus $\gamma=1.17\pm 0.11$) | $\mathcal{O}(0.1)$ | **far above** prediction |
+ | Current DESI-era slip (Maus $\gamma=1.17\pm 0.11$) | $\mathcal{O}(0.1)$ | **above** prediction ($\sim 2\times$) |
```

**`papers/anisotropic-slip-option0.md` l.158** — 0.03/4.95e-2 = 0.61; 0.05/4.95e-2 = 1.01: el mejor caso queda a la altura del suelo de Euclid

```diff
- | Euclid-like forecasts (Sakr et al.) | $\mathcal{O}(0.03\text{–}0.05)$ | **still above** by $\sim 0.9$–$1.5\times$ |
+ | Euclid-like forecasts (Sakr et al.) | $\mathcal{O}(0.03\text{–}0.05)$ | **comparable**: $\sim 0.6$–$1.0\times$ |
```

**`papers/anisotropic-slip-option0.md` l.166** — slip regenerado (ε = 1)

```diff
- | 1 | $3.39\times 10^{-2}$ |
+ | 1 | $4.95\times 10^{-2}$ |
```

**`papers/anisotropic-slip-option0.md` l.167** — slip regenerado (ε = 0.1)

```diff
- | 0.1 | $3.39\times 10^{-3}$ |
+ | 0.1 | $4.95\times 10^{-3}$ |
```

**`papers/anisotropic-slip-option0.md` l.168** — slip regenerado (ε = 0.01)

```diff
- | 0.01 | $3.39\times 10^{-4}$ |
+ | 0.01 | $4.95\times 10^{-4}$ |
```

**`papers/anisotropic-slip-option0.md` l.177** — slip regenerado: 3e-2 → 5e-2

```diff
- 2. Published growth/MG analyses with DESI-like data are GR-compatible at the $\mathcal{O}(0.1)$ level, which is still too coarse to see $\lvert \eta-1 \rvert\sim 3\times 10^{-2}$.
+ 2. Published growth/MG analyses with DESI-like data are GR-compatible at the $\mathcal{O}(0.1)$ level, which is still too coarse to see $\lvert \eta-1 \rvert\sim 5\times 10^{-2}$.
```

**`papers/resume.txt` l.16** — slip regenerado: 3e-2 → 5e-2; 1.0–4.7× por debajo del suelo de forecast

```diff
-   Act IV bridge: scripts/slip_bridge.py  (|gamma-1|~3e-2 best case; ~1.5-7x starved)
+   Act IV bridge: scripts/slip_bridge.py  (|gamma-1|~5e-2 best case; ~1-5x below forecast floor)
```

### c. m_φ retirado

**`papers/stochastic-dark-energy-desi-dr2.md` l.476** — m_φ retirado (verificacion_numeros filas 6–7: el valor no se reproduce)

```diff
- The upper limit translates into a constraint on the mass of a possible scalar field mediating dark energy: $m_\phi \lesssim 10^{-5}\ \text{eV}$ [valor no reproducible; pendiente de decisión del autor]. This is consistent with ultralight boson scenarios and with astrophysical bounds on the variation of fundamental constants ($\Delta\alpha/\alpha \lesssim 10^{-5}$ from quasar spectroscopy).
+ A scalar-field reading would translate $\sigma_X$ into a bound on $m_\phi$; it is not derived here.
```

**`papers/stochastic-dark-energy-desi-dr2.md` l.621** — m_φ retirado (verificacion_numeros fila 6)

```diff
- Assuming that dark energy is a canonical scalar field $\phi$ with $\rho_\Lambda \approx V(\phi)$ and $\phi \sim M_{\rm Pl}$, and using $\sigma_X \sim m_\phi^2 M_{\rm Pl} H_0 / (2\pi \rho_\Lambda)$, the limit $\sigma_X < 2.5\times 10^{-2}$ implies:
- 
- 
- $$
- m_\phi \lesssim 9.45 \times 10^{-5} \text{eV}. \quad \text{[valor no reproducible; pendiente de decisión del autor]}
- $$
+ A scalar-field reading would translate $\sigma_X$ into a bound on $m_\phi$; it is not derived here.
```

**`papers/stochastic-dark-energy-desi-dr2.md` l.629** — consecuencia de retirar m_φ: la frase siguiente hablaba de esa traducción

```diff
- This is consistent with ultralight boson scenarios. However, this translation depends on the assumed normalization of $\phi$ and the relation between $\sigma_X$ and $\delta\phi$. It is not a direct measurement from DESI data.
+ Such a translation would depend on the assumed normalization of $\phi$ and on the relation between $\sigma_X$ and $\delta\phi$; it would not be a direct measurement from DESI data.
```

**`theory/unimodular-gravity-vacuum-smoothness.md` l.126** — m_φ y |V'/V|M_Pl retirados (verificacion_numeros fila 8: ninguno se reproduce)

```diff
- If $m_\phi \ll H_0$ (where $H_0 \sim 10^{-33}\ \text{eV}$), the friction term $3H \delta\dot{\phi}$ forces the mode into the “frozen” regime and effectively erases inhomogeneities. The limits I derived in [1] give:
- 
- 
- $$
- m_\phi \lesssim 10^{-5}\ \text{eV}, \qquad \left| \frac{V'}{V}\right| M_{\mathrm{Pl}} \lesssim 10^{-2}. \quad \text{[valor no reproducible; pendiente de decisión del autor]}
- $$
+ If $m_\phi \ll H_0$ (where $H_0 \sim 10^{-33}\ \text{eV}$), the friction term $3H \delta\dot{\phi}$ forces the mode into the “frozen” regime and effectively erases inhomogeneities. A scalar-field reading would translate $\sigma_X$ into a bound on $m_\phi$; it is not derived here.
```


## Filas de la tabla 1 que no se cambian

| Fichero | Motivo |
|:--|:--|
| `papers/scope-and-mixups.md` l.12, 18, 23, 54 | La tabla 1 lo da por correcto ("OK; opcional"). Ya usa X ≡ δΩ_Λ, S_iS_j σ_X² y σ_X = σ/√(2θ). |
| `papers/stochastic-dark-energy-desi-dr2.md` l.160 | "X ≡ δΩ_Λ" es coherente con [X-std]. |
| `papers/stochastic-dark-energy-desi-dr2.md` l.164 | La ecuación de continuidad en δρ_Λ no lleva normalización; el cambio de la fila 14 está en l.172. |
| `papers/stochastic-dark-energy-desi-dr2.md` l.264 | Ya tenía S(z_i)S(z_j); la definición de σ_X² se añade en l.272. |
| `archive/scripts/ou_bao_stochastic_test.py` l.83 | `SIGMA_OU`: el script está archivado (reorganización 1/2) y la instrucción solo pedía `ou_bao_likelihood.py`. |
| `exploratory/notes/data-pack-option0-internet.md` l.87, 91–94 | Estimación de orden de magnitud \|η−1\| ∼ f·(Ω_DE/Ω_m)·σ_X ∼ 5e-2. No está en la tabla 1 ni en la instrucción. Con [X-std] el factor sería 1/Ω_m (∼8e-2 con f = 1), pero el texto es explícitamente "rough". Queda para el autor. |
| `analisis_cmb/verificacion_numeros.md` filas 13–18 | Registro de B1. Sus valores "recalculados" del slip (3.393e-2, …) usan ρ_X = Ω_Λ0 ρ_crit; los válidos tras C están en el resumen de arriba y en `logs/slip_bridge_despues.log`. |

## d. Comprobación

**`git grep` en `manuscript/ papers/ theory/ exploratory/notes/`** (después de los cambios):

- `fractional`: **0 apariciones**.
- `δΩ_Λ/Ω` (unicode o LaTeX `\delta\Omega_\Lambda/\Omega`): **0 apariciones**.
- `δρ_Λ/ρ_Λ`: la forma unicode no aparece. Con LaTeX (`\delta\rho_\Lambda` seguido de `/` o `}{` y `\rho_\Lambda`), quedan estas apariciones, todas del tipo "rms δρ_Λ/ρ_Λ = σ_X/Ω_Λ0" o la conversión X = Ω_Λ0 δρ_Λ/ρ̄_Λ:

| Fichero:línea (después de C) | Texto | Por qué se queda |
|:--|:--|:--|
| `papers/stochastic-dark-energy-desi-dr2.md:404` | rms(δρ_Λ)/ρ_Λ = σ_X/Ω_Λ0 ≈ (V'/V) rms(δφ) | Es la relación pedida (factor 1/Ω_Λ0). |
| `papers/stochastic-dark-energy-desi-dr2.md:820` | rms(δρ_Λ)/ρ_Λ = σ_X/Ω_Λ0 | Ídem (apéndice H, antes l.825). |
| `theory/principle-of-vacuum-smoothness.md:15` | rms δρ_Λ/ρ_Λ = σ_X/Ω_Λ0 < 3.84e-2 | Ídem. |
| `theory/unimodular-gravity-vacuum-smoothness.md:68` | X ≡ δρ_DE/ρ_crit,0 = Ω_Λ0 · δρ_Λ/ρ̄_Λ | Conversión de [X-std] a la fluctuación relativa de ρ_Λ (la propone la tabla 1). |
| `theory/unimodular-gravity-vacuum-smoothness.md:72` | rms δρ_Λ/ρ̄_Λ = σ_X/Ω_Λ0 | Relación pedida. |

Otras apariciones de δρ_Λ que no son una definición de X ni de σ_X (no se cambian): `papers/quantum-fluid-instabilities-desi-dr2.md:52, 97, 100` (la perturbación δρ_Λ = ψ₀*δψ + c.c. y su semilla δρ_Λ⁽⁰⁾ en el modelo de fluido cuántico); `papers/stochastic-dark-energy-desi-dr2.md:164, 633` (la ecuación de continuidad en δρ_Λ); `:808, 812` (δρ_Λ ≃ V'δφ, apéndice H).

Otras comprobaciones:
- **Factor S_iS_j:** ninguna covarianza OU/QNM de BAO queda sin S_iS_j (`PREPRINT.md:94`, `stochastic…:264`, `scope-and-mixups.md:23`, `euclid-protocol…:54`, nota desqueezing `:201`). La de `quantum-fluid…:127` es la de rango 1, con S(z_i)S(z_j).
- **m_φ:** solo queda la frase nueva (`stochastic…:476` y `:621`) y, en unimodular l.122/126, la ecuación del modo δφ y la condición m_φ ≪ H_0, que no son cotas. No quedan "9.45", "valor no reproducible" ni "pendiente de decisión".
- **Scripts del claim set:** se volvieron a ejecutar desde la raíz del repo. Después, `git status --porcelain results/ figures/` no muestra nada. Los logs están en `bloque_C/logs/`, con las rutas locales quitadas.

| Script | Código de salida | Tiempo | Nota |
|:--|:--:|--:|:--|
| `scripts/verify_bao_pipelines.py` | 0 | 134 s | |
| `scripts/profile_sigma_x_desi.py` | 0 | 1 s | |
| `scripts/gap_two_targets.py` | 0 | 0 s | |
| `scripts/slip_bridge.py` | 0 | 0 s | solo imprime; números nuevos en `logs/slip_bridge_despues.log` (antes: `logs/slip_bridge_antes.log`) |
| `scripts/amplifier_audit.py` | 0 | 0 s | |
| `scripts/amplification/route1_local_causal_set_seed.py` | 0 | 0 s | |
| `scripts/eos_efectiva.py` | 0 | 134 s | renombre `sigma` → `sigma_w`; `results/eos_cpl_desi_dr2/` idéntico |
| `scripts/tachyonic_rank1_mle.py` | 0 | 0 s | |
| `scripts/joint_w0wa_sigma_desi.py` | 0 | 2 s | |
| `scripts/ou_bao_likelihood.py` | 0 | 1 s | renombre `SIGMA_OU` → `SIGMA_X_DR1`; escribe en `plots/` (ignorado) |
