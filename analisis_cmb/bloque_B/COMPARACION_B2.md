# Bloque B2: números de los papers alineados con la cota 2.633e-2 (2026-10-17)

**Reglas:**
- Las ediciones se hacen línea a línea, sobre el texto existente.
- La cota se escribe 2.6×10⁻², o 2.633×10⁻² donde el texto da tres cifras o el número entra en una cuenta: tablas de slip, H.5, BOUND_CORRECTION.
- Cada número sale de `analisis_cmb/verificacion_numeros.md` (columna "fuente": "verificacion fila N") o de los json/csv regenerados en B1:
  - `results/profile_sigma_x/profile_sigma_x.json`;
  - `results/amplification_routes/route1_N_eff_required.csv` y `route3_avalanche_scan.csv`;
  - `results/joint_w0wa_sigma/joint_w0wa_sigma.json`.

**Ficheros:** 20 ficheros y 139 cambios, en manuscript/, papers/, README.md y notes/.

**Excluido, como se pidió:**
- **m_φ:**
  - stochastic l.476 y l.625, y unimodular l.130: solo se añade "[valor no reproducible; pendiente de decisión del autor]", sin cambiar el número.
  - stochastic l.621, la derivación de m_φ, que contiene "σ_X < 2.5×10⁻²", **no se toca**: cambiar la cota ahí sin cambiar el m_φ resultante dejaría la cuenta incoherente.
- **Fórmula deducida (C):**
  - stochastic l.372 "more than four orders of magnitude": solo se cambia el número de la cota en la misma línea; la frase de órdenes de magnitud queda para C.
  - stochastic l.164/172 (normalización del ruido).
  - amplification-gap l.54/56 (gap en densidad con 1/Ω_Λ0).
  - Todas las filas de definición de X (punto 1 de AUDITORIA), incluida stochastic l.404/412.
- **Euclid y Everpresent (D):** filas de las tablas 4 y 5 de AUDITORIA.
  - Las filas de `euclid-protocol-vacuum-relaxation.md` que están en la tabla 2a (l.81, l.148) sí se aplican: son números de σ_X y θ_min, no fechas ni bins de Euclid.

**Filas de 2a hechas antes:**
- En el bloque A: PREPRINT l.76, stochastic l.430/932, l.440, l.926.
- En B1: los ficheros R y los scripts.

## Cambios línea a línea

### README.md

| línea | antes | después | fuente |
|---|---|---|---|
| 32 | \sigma_X < 2.5\times 10^{-2}. | \sigma_X < 2.6\times 10^{-2}. | cota |

### manuscript/PREPRINT.md

| línea | antes | después | fuente |
|---|---|---|---|
| 18 | maximum-likelihood estimation drives the stochastic amplitude to the numerical floor. The working 95% confidence upper limit is | the best fit drives the stochastic amplitude to zero, and the profile likelihood over $\theta$ ($\Delta\ln\mathcal{L}=-1.92$, full measurement covariance) gives the 95% confidence upper limit | 2a l.74; profile_sigma_x.json criterion |
| 20 | $$\sigma_X < 2.5\times 10^{-2}$$ | $$\sigma_X < 2.6\times 10^{-2}$$ | cota |
| 26 | below the DESI residual ceiling measured here ($2.5\times 10^{-2}$) | below the DESI residual ceiling measured here ($2.6\times 10^{-2}$) | cota |
| 28 | gives $\lvert\gamma-1\rvert\sim 10^{-4}$, of order $10^{2}$–$10^{3}$ below current and forecast floors. | gives $\lvert\gamma-1\rvert\sim 3\times 10^{-2}$, of order $1.5$–$7\times$ below current and forecast floors. | verificacion_numeros.md filas 15–17 |
| 144 | \| OU (H0) \| $0.001$ (boundary) \| $5\times 10^{-5}$ \| | \| OU (H0; old run, diagonal cov) \| $0.001$ (boundary) \| $5\times 10^{-5}$ \| | 2a l.77; results/ou_bao_desi_dr2_run.txt |
| 145 | \| QNM (H1) \| $0.001$ \| $5\times 10^{-5}$ \| | \| QNM (H1; old run, diagonal cov) \| $0.001$ \| $5\times 10^{-5}$ \| | 2a l.77 |
| 151 | public pipeline): $\sigma_X < 2.5\times 10^{-2}$ (95% CL). | public pipeline): $\sigma_X < 2.6\times 10^{-2}$ (95% CL). | cota |
| 191 | \| DESI residual bound (this work, measured) \| $\sigma_X<2.5\times 10^{-2}$ (95% CL)\| | \| DESI residual bound (this work, measured) \| $\sigma_X<2.6\times 10^{-2}$ (95% CL)\| | cota |
| 194 | \| $G_{\mathrm{DESI}}=2.5\times 10^{-2}/\sigma_0$ \| $1.27\times 10^{59}\approx 10^{59}$ \| | \| $G_{\mathrm{DESI}}=2.6\times 10^{-2}/\sigma_0$ \| $2.24\times 10^{59}\approx 10^{59}$ \| | verificacion_numeros.md fila 1 |
| 196 | is incorrect by a factor $\sim 15$. Both gaps remain unbridgeable by soft linear amplifiers ($r\sim 64$–$66$). | is incorrect by a factor $\sim 2\times 10^{3}$. Both gaps remain unbridgeable by soft linear amplifiers ($r\sim 64$–$68$). | verificacion_numeros.md fila 1 y nota |
| 214 | $N_{\mathrm{eff}}\sim 4.44\times 10^{7}$ (DESI) | $N_{\mathrm{eff}}\sim 1.44\times 10^{3}$ (DESI) | verificacion_numeros.md fila 3; route1_N_eff_required.csv l.3 |
| 216 | BAO-safe if $\sigma\lesssim 4\times 10^{-5}$. No free gain. | All 288 jobs BAO-safe (max p95 $=2.1\times 10^{-3}$). No free gain. | route3_avalanche_scan.csv (288/288, l.289); VERDICT.md l.55 |
| 222 | 3. **Amplification gap:** $\sigma_0\to 10^{-5}$ is $\sim 10^{56}$ (Euclid target); $\sigma_0\to 2.5\times 10^{-2}$ is $\sim 10^{59}$ (DESI ceiling); linear free amplifiers fail both. | 3. **Amplification gap:** $\sigma_0\to 10^{-5}$ is $\sim 10^{56}$ (Euclid target); $\sigma_0\to 2.6\times 10^{-2}$ is $\sim 10^{59}$ (DESI ceiling); linear free amplifiers fail both. | cota |
| 242 | With $\sigma_X=2.5\times 10^{-2}$, $\varepsilon=1$, $\delta_m=1$: | With $\sigma_X=2.633\times 10^{-2}$, $\varepsilon=1$, $\delta_m=1$: | cota (tres cifras: entrada de la tabla) |
| 246 | \| 0.5 \| $1.93\times 10^{-4}$ \| $\sim 880\times$ short \| $\sim 260\times$ short \| | \| 0.5 \| $3.39\times 10^{-2}$ \| $\sim 5.0\times$ short \| $\sim 1.5\times$ short \| | verificacion_numeros.md filas 15–17 |
| 247 | \| 1.0 \| $8.16\times 10^{-5}$ \| $\sim 2\times 10^{3}$ short \| $\sim 610\times$ short \| | \| 1.0 \| $1.43\times 10^{-2}$ \| $\sim 12\times$ short \| $\sim 3.5\times$ short \| | verificacion_numeros.md filas 15–17 |
| 248 | \| 1.5 \| $4.18\times 10^{-5}$ \| $\sim 4\times 10^{3}$ short \| $\sim 1.2\times 10^{3}$ short \| | \| 1.5 \| $7.33\times 10^{-3}$ \| $\sim 23\times$ short \| $\sim 6.8\times$ short \| | verificacion_numeros.md filas 15–17 |
| 262 | baseline residual analysis uses published BAO $\alpha$ and diagonal errors. Full DESI covariance matrices may change the numerical factor | baseline residual analysis uses published BAO $\alpha$ with the official Gaussian BAO covariance projected to $\alpha$ (block-diagonal across bins). Cross-bin covariances, if released, may change the numerical factor | 2a l.83; profile_sigma_x.json cov_mode |
| 319 | do not favour stationary OU/QNM stochastic dark energy; $\sigma_X<2.5\times 10^{-2}$ (95% CL) | do not favour stationary OU/QNM stochastic dark energy; $\sigma_X<2.6\times 10^{-2}$ (95% CL) | cota |
| 321 | 3. **Amplification gap:** $\sigma_0\to 10^{-5}$ is $\sim 10^{56}$ (Euclid target); $\sigma_0\to 2.5\times 10^{-2}$ is $\sim 10^{59}$ (DESI ceiling); linear free amplifiers fail both. | 3. **Amplification gap:** $\sigma_0\to 10^{-5}$ is $\sim 10^{56}$ (Euclid target); $\sigma_0\to 2.6\times 10^{-2}$ is $\sim 10^{59}$ (DESI ceiling); linear free amplifiers fail both. | cota |

### manuscript/CLAIMS.md

| línea | antes | después | fuente |
|---|---|---|---|
| 17 | \| **C2** \| $\sigma_X < 2.5\times 10^{-2}$ | \| **C2** \| $\sigma_X < 2.6\times 10^{-2}$ | cota |
| 20 | gap $\to 2.5\times 10^{-2}$ is $\sim 10^{59}$ | gap $\to 2.6\times 10^{-2}$ is $\sim 10^{59}$ | cota |
| 22 | gives $\\|\gamma-1\\|\sim 10^{-4}$ \| | gives $\\|\gamma-1\\|\lesssim 3\times 10^{-2}$ ($z=0.5$, $\delta_m=1$) \| | verificacion_numeros.md fila 15 |

### manuscript/COVER_LETTER.md

| línea | antes | después | fuente |
|---|---|---|---|
| 11 | ($\sigma_X < 2.5\times 10^{-2}$ (95% CL) | ($\sigma_X < 2.6\times 10^{-2}$ (95% CL) | cota |

### manuscript/README.md

| línea | antes | después | fuente |
|---|---|---|---|
| 31 | - $\sigma_X < 2.5\times 10^{-2}$ (95% CL) | - $\sigma_X < 2.6\times 10^{-2}$ (95% CL) | cota |

### papers/stochastic-dark-energy-desi-dr2.md

| línea | antes | después | fuente |
|---|---|---|---|
| 14 | I perform a Maximum Likelihood Estimation (MLE) to constrain the amplitude $\sigma_X$ of these fluctuations, assuming the CPL background is fixed to the best-fit values. | I constrain the amplitude $\sigma_X$ of these fluctuations with a profile likelihood over $\theta$ (official Gaussian BAO covariance projected to $\alpha$), with the background fixed at $\alpha=1$ ($\Lambda$CDM fiducial). | 2a l.100; profile_sigma_x.json notes |
| 16 | I find that the MLE drives the stochastic amplitude to the numerical floor: $\sigma_X \to 0$ and $\omega_R \to 0$. | I find that the best fit drives the stochastic amplitude to zero: $\sigma_X \to 0$ and $\omega_R \to 0$. | 2a l.100 |
| 16 | This result places a phenomenological upper limit on the amplitude of such fluctuations: ⏎  ⏎  ⏎ $$ ⏎ \sigma_X < 2.5\times 10^{-2} ⏎ $$ | The profile ($\Delta\ln\mathcal{L}=-1.92$) places a phenomenological upper limit on the amplitude of such fluctuations: ⏎  ⏎  ⏎ $$ ⏎ \sigma_X < 2.6\times 10^{-2} ⏎ $$ | cota; 2a l.100 |
| 224 | \sigma_X < 2.5\times 10^{-2}. | \sigma_X < 2.6\times 10^{-2}. | cota |
| 278 | I compute the likelihood using the public DESI DR2 BAO measurements (7 bins from arXiv:2503.14738, Table 1) and the standard diagonal covariance from measurement errors. | I compute the likelihood using the public DESI DR2 BAO measurements (7 bins from arXiv:2503.14738, Table 1); the runs in the table below use the diagonal covariance from measurement errors (old runs), while the bound of Sec. 7.1 uses the official covariance projected to $\alpha$. | 2a l.101 |
| 287 | \| H0: OU free MLE \| 0.0010 \| | \| H0: OU free MLE (old run, diagonal cov) \| 0.0010 \| | 2a l.102 |
| 288 | \| H1: QNM free MLE \| 0.0010 \| | \| H1: QNM free MLE (old run, diagonal cov) \| 0.0010 \| | 2a l.102 |
| 368 | The Maximum Likelihood Estimation (MLE) drives the stochastic noise parameter to the numerical floor, yielding a conservative $95\%$ confidence level upper limit of $\sigma_X < 2.5\times 10^{-2}$. This limit applies to the specific OU/QNM kernel tested here, assuming a fixed CPL background. | The best fit drives the stochastic noise parameter to zero, and the profile likelihood over $\theta$ ($\Delta\ln\mathcal{L}=-1.92$) gives the $95\%$ confidence level upper limit $\sigma_X < 2.6\times 10^{-2}$. This limit applies to the specific OU kernel tested here, with the background fixed at $\alpha=1$ ($\Lambda$CDM fiducial). | cota; 2a l.100 |
| 370 | The upper limit $\sigma_X < 2.5\times 10^{-2}$ is a conservative estimate based on the numerical floor of the MLE ($\sigma_X \approx 5 \times 10^{-5}$), multiplied by a factor of 3 to account for the flatness of the likelihood in $\sigma_X$ with only 7 bins. A formal 95% CL profile likelihood scan over $\sigma_X$ is deferred to Euclid DR1; with the current dataset, the likelihood is dominated by the optimizer's numerical floor. | The upper limit $\sigma_X < 2.6\times 10^{-2}$ comes from a formal profile likelihood scan over $\sigma_X$ (`scripts/profile_sigma_x_desi.py`): root of $\Delta\ln\mathcal{L}=-1.92$ at $\sigma_X=2.633\times 10^{-2}$, between grid nodes; with only 7 bins the profile is flat at small $\sigma_X$ and the best fit is $\sigma_X\to 0$. | 2a l.103; profile_sigma_x.json |
| 372 | The upper limit $\sigma_X < 2.5\times 10^{-2}$ implies that the dark energy density | The upper limit $\sigma_X < 2.6\times 10^{-2}$ implies that the dark energy density | cota (el resto de l.372 va en C) |
| 376 | any such model must have an amplitude $\sigma_X < 2.5\times 10^{-2}$ to be compatible | any such model must have an amplitude $\sigma_X < 2.6\times 10^{-2}$ to be compatible | cota |
| 416 | My upper bound $\sigma_X < 2.5\times 10^{-2}$ therefore implies: | My upper bound $\sigma_X < 2.6\times 10^{-2}$ therefore implies: | cota |
| 420 | \right\| (2.5\times 10^{-2}). | \right\| (2.6\times 10^{-2}). | cota |
| 470 | \sigma_X < 2.5\times 10^{-2} | \sigma_X < 2.6\times 10^{-2} | cota |
| 476 | $m_\phi \lesssim 10^{-5}\ \text{eV}$. This is consistent | $m_\phi \lesssim 10^{-5}\ \text{eV}$ [valor no reproducible; pendiente de decisión del autor]. This is consistent | m_φ: verificacion_numeros.md filas 6–8 |
| 478 | More broadly, $\sigma_X < 2.5\times 10^{-2}$ implies | More broadly, $\sigma_X < 2.6\times 10^{-2}$ implies | cota |
| 617 | of the phenomenological limit $\sigma_X < 2.5\times 10^{-2}$. These | of the phenomenological limit $\sigma_X < 2.6\times 10^{-2}$. These | cota |
| 625 | m_\phi \lesssim 9.45 \times 10^{-5} \text{eV}. | m_\phi \lesssim 9.45 \times 10^{-5} \text{eV}. \quad \text{[valor no reproducible; pendiente de decisión del autor]} | m_φ: verificacion_numeros.md fila 6 |
| 865 | (Sec. 4.2, $\sigma_X < 2.5\times 10^{-2}$ at 95% CL) | (Sec. 4.2, $\sigma_X < 2.6\times 10^{-2}$ at 95% CL) | cota |
| 869 | \lesssim \frac{2.7 \times 10^{-6}}{2.5\times 10^{-2}} \approx 1.8 \times 10^{-2}. \tag{H.5} | \lesssim \frac{2.7 \times 10^{-6}}{2.633\times 10^{-2}} \approx 1.0 \times 10^{-4}. \tag{H.5} | verificacion_numeros.md fila 10 |
| 883 | \left\| \frac{V(\phi)}{M_{\rm Pl} V'(\phi)} \right\| \lesssim 1.8 \times 10^{-2}. | \left\| \frac{V(\phi)}{M_{\rm Pl} V'(\phi)} \right\| \lesssim 1.0 \times 10^{-4}. | verificacion_numeros.md fila 10 |
| 887 | hence $\lambda \gtrsim 55.5$. | hence $\lambda \gtrsim 9.8\times 10^{3}$. | verificacion_numeros.md fila 11 |
| 895 | \lvert \zeta \rvert \lesssim 1.8 \times 10^{-2}. | \lvert \zeta \rvert \lesssim 1.0 \times 10^{-4}. | verificacion_numeros.md fila 12 |
| 899 | \gtrsim 55 M_{\rm Pl}$ | \gtrsim 9.8\times 10^{3} M_{\rm Pl}$ | verificacion_numeros.md fila 12 (l.899, no marcada en 2a) |
| 906 | ($\lvert \zeta \rvert \lesssim 10^{-2}$) | ($\lvert \zeta \rvert \lesssim 10^{-4}$) | verificacion_numeros.md fila 12 |
| 908 | interpretation of the bound $\sigma_X < 2.5\times 10^{-2}$ as a robust | interpretation of the bound $\sigma_X < 2.6\times 10^{-2}$ as a robust | cota |
| 924 | $\sigma_X \approx 9.2 \times 10^{-3}$ — nearly 200 times larger than the reported limit, | $\sigma_X \approx 9.2 \times 10^{-3}$ — about 0.35 of the 95% CL limit $2.6\times 10^{-2}$ (≈180× the $5\times 10^{-5}$ floor of the old diagonal run), | 2a l.108 |

### papers/resume.txt

| línea | antes | después | fuente |
|---|---|---|---|
| 13 | R1: N_eff ~ 4.4e7 (DESI) | R1: N_eff ~ 1.4e3 (DESI) | verificacion_numeros.md fila 3 |
| 15 | BAO-safe if sigma ≲ 4e-5 | all 288 jobs BAO-safe (max p95 2.1e-3) | route3_avalanche_scan.csv |
| 16 | (\|gamma-1\|~1e-4 best case; O(1e2-1e3)x starved) | (\|gamma-1\|~3e-2 best case; ~1.5-7x starved) | verificacion_numeros.md filas 15–17 |
| 39 |   OU model (H0): ⏎     θ     = 0.001  (≈ boundary, effectively undamped) ⏎     σ_X   = 5 × 10⁻⁵ ⏎  |   OU model (H0): ⏎     θ     = 0.001  (≈ boundary, effectively undamped) ⏎     σ_X   = 5 × 10⁻⁵   (old run, diagonal cov: results/ou_bao_desi_dr2_run.txt) ⏎  | 2a l.98 |
| 46 |     θ     = 0.001 ⏎     σ_X   = 5 × 10⁻⁵ ⏎  |     θ     = 0.001 ⏎     σ_X   = 5 × 10⁻⁵   (old run, diagonal cov) ⏎  | 2a l.98 |
| 55 | - 95% CL upper limit: σ_X < 2.5 × 10⁻² | - 95% CL upper limit: σ_X < 2.6 × 10⁻² (profile over θ, background fixed at α = 1) | cota; 2a l.99 |
| 83 | - When marginalizing over the background cosmology, the stochastic amplitude is consistent with zero. | - When the background (w0, wa) is left free (results/joint_w0wa_sigma), the stochastic amplitude is consistent with zero. | 2a l.99; joint_w0wa_sigma.json |
| 84 | - This constitutes a **strict upper limit**, not a detection. | - The 95% CL bound is a profile over θ with the background fixed at α = 1: an upper limit, not a detection. | 2a l.99 |

### papers/HONEST_HEADLINES.md

| línea | antes | después | fuente |
|---|---|---|---|
| 18 | glued to DESI $\sigma_X=2.5\times 10^{-2}$ | glued to DESI $\sigma_X=2.6\times 10^{-2}$ | cota |
| 25 | 1. Residual bound $\sigma_X < 2.5\times 10^{-2}$ | 1. Residual bound $\sigma_X < 2.6\times 10^{-2}$ | cota |
| 27 | soft $r\sim 64$ (Euclid target) / $\sim 66$ (DESI ceiling) | soft $r\sim 64$ (Euclid target) / $\sim 68$ (DESI ceiling) | verificacion_numeros.md fila 1 |
| 36 | no fluctuations above $2.5\times 10^{-2}$ | no fluctuations above $2.6\times 10^{-2}$ | cota |
| 42 | with DESI DR2 ($\sigma_X<2.5\times 10^{-2}$ (95% CL)) | with DESI DR2 ($\sigma_X<2.6\times 10^{-2}$ (95% CL)) | cota |

### papers/BOUND_CORRECTION_SIGMA_X.md

| línea | antes | después | fuente |
|---|---|---|---|
| 10 | \| Number \| $\sigma_X < 1.5\times 10^{-4}$ \| $\sigma_X < 2.5\times 10^{-2}$ \| | \| Number \| $\sigma_X < 1.5\times 10^{-4}$ \| $\sigma_X < 2.6\times 10^{-2}$ \| | cota |
| 11 | \| Exact grid \| — \| $2.506\times 10^{-2}$ \| | \| Exact \| — \| $2.633\times 10^{-2}$ (root of the profile; last grid node $2.506\times 10^{-2}$, previous version) \| | profile_sigma_x.json sigma_X_95CL_upper_profile y _grid |
| 29 | $\mathrm{ceiling}/\sigma_{\mathrm{free}}\approx 2.5\times 10^{-2}/8.5\times 10^{-5}\approx 294$ | $\mathrm{ceiling}/\sigma_{\mathrm{free}}\approx 2.633\times 10^{-2}/8.5\times 10^{-5}\approx 310$ | cota / σ_free (misma línea) |

### papers/amplification-gap.md

| línea | antes | después | fuente |
|---|---|---|---|
| 37 | \| DESI DR2 BAO residual bound \| $\sigma_X < 2.5\times 10^{-2}$ (95% CL) \| | \| DESI DR2 BAO residual bound \| $\sigma_X < 2.6\times 10^{-2}$ (95% CL) \| | cota |
| 43 | The answers differ by about **one decade**. | The answers differ by about **three decades** ($10^{56}$ vs $10^{59}$). | verificacion_numeros.md fila 1 |
| 48 | \| **$G_{\mathrm{DESI}}$** \| $2.5\times 10^{-2}/\sigma_0$ \| $\mathbf{1.27\times 10^{59}\approx 10^{59}}$ \| $r\approx 65.7$ \| | \| **$G_{\mathrm{DESI}}$** \| $2.6\times 10^{-2}/\sigma_0$ \| $\mathbf{2.24\times 10^{59}\approx 10^{59}}$ \| $r\approx 68.3$ \| | verificacion_numeros.md fila 1 |
| 53 | for the **measured DESI ceiling** $\sigma_X=2.5\times 10^{-2}$. | for the **measured DESI ceiling** $\sigma_X=2.6\times 10^{-2}$. | cota |
| 114 | must still satisfy $\sigma_X < 2.5\times 10^{-2}$ (95% CL) \| | must still satisfy $\sigma_X < 2.6\times 10^{-2}$ (95% CL) \| | cota |
| 191 | \| DESI 95% ceiling \| $2.5\times 10^{-2}$ \| $\mathbf{4.44\times 10^{7}}$ \| smaller by $\sim 10^{114}$ \| | \| DESI 95% ceiling \| $2.6\times 10^{-2}$ \| $\mathbf{1.44\times 10^{3}}$ \| smaller by $\sim 10^{119}$ \| | verificacion_numeros.md fila 3 |
| 199 | \| DESI $N_{\mathrm{eff}}\sim 4.4\times 10^{7}$ \| $\sim 10^{-31} \mathrm{m}$ \| $\sim 10^{-12} \mathrm{m}$ \| $\sim\mathrm{mm}$ \| | \| DESI $N_{\mathrm{eff}}\sim 1.4\times 10^{3}$ \| $\sim 5\times 10^{-34} \mathrm{m}$ \| $\sim 3\times 10^{-14} \mathrm{m}$ \| $\sim 0.3\,\mathrm{mm}$ \| | verificacion_numeros.md fila 4 |
| 227 | \| BAO-safe ($\mathrm{p95}<2.5\times 10^{-2}$) \| **216 / 288** \| | \| BAO-safe ($\mathrm{p95}<2.6\times 10^{-2}$) \| **288 / 288** \| | route3_avalanche_scan.csv |
| 228 | \| BAO-unsafe ($\mathrm{p95}\ge 2.5\times 10^{-2}$) \| **72 / 288** — **all** have $\sigma\gtrsim 1.2\times 10^{-4}$ \| | \| BAO-unsafe ($\mathrm{p95}\ge 2.6\times 10^{-2}$) \| **0 / 288** (max p95 $=2.1\times 10^{-3}$) \| | route3_avalanche_scan.csv l.289 |
| 229 | \| Euclid window $10^{-5}\le\mathrm{p95}<2.5\times 10^{-2}$ \| 72 jobs | \| Euclid window $10^{-5}\le\mathrm{p95}<2.6\times 10^{-2}$ \| 144 jobs | route3_avalanche_scan.csv (recuento) |
| 231 | **BAO-safe rule of thumb in this scan:** keep drive noise $\sigma\lesssim 4\times 10^{-5}$ (p95 stays below DESI ceiling for all $ (a,b)$ tested). ⏎ **BAO-unsafe:** $\sigma\gtrsim 1.2\times 10^{-4}$ → p95 exceeds DESI for **every** $ (a,b)$. | **BAO-safe in this scan:** the whole grid ($\sigma\le 10^{-3}$) keeps p95 below the DESI ceiling for all $ (a,b)$ tested. | route3_avalanche_scan.csv; 2a l.90 (borrar l.232) |

### papers/anisotropic-slip-option0.md

| línea | antes | después | fuente |
|---|---|---|---|
| 14 | $\sigma_X < 2.5\times 10^{-2}$, | $\sigma_X < 2.6\times 10^{-2}$, | cota |
| 25 | but with $\sigma_X\sim 10^{-4}$ it inherits | but with $\sigma_X\lesssim 2.6\times 10^{-2}$ it inherits | cota |
| 118 | far above the $\sim 10^{-4}$ scale of the $f=1$, $\sigma_X\sim 10^{-4}$ estimate below. | far above the $\sim 3\times 10^{-2}$ scale of the $f=1$, $\sigma_X=2.6\times 10^{-2}$ estimate below. | verificacion_numeros.md fila 15 |
| 145 | **Machine check (best case $\varepsilon=1$, $\sigma_X=2.5\times 10^{-2}$, $\delta_m=1$):** | **Machine check (best case $\varepsilon=1$, $\sigma_X=2.633\times 10^{-2}$, $\delta_m=1$):** | cota (tres cifras: entrada de la tabla) |
| 149 | \| 0.5 \| 0.644 \| $1.93\times 10^{-4}$ \| $\sim 880\times$ \| $\sim 260\times$ \| | \| 0.5 \| 0.644 \| $3.39\times 10^{-2}$ \| $\sim 5.0\times$ \| $\sim 1.5\times$ \| | verificacion_numeros.md filas 15–17 |
| 150 | \| 1.0 \| 0.272 \| $8.16\times 10^{-5}$ \| $\sim 2\times 10^{3}$ \| $\sim 610\times$ \| | \| 1.0 \| 0.272 \| $1.43\times 10^{-2}$ \| $\sim 12\times$ \| $\sim 3.5\times$ \| | verificacion_numeros.md filas 15–17 |
| 151 | \| 1.5 \| 0.139 \| $4.18\times 10^{-5}$ \| $\sim 4\times 10^{3}$ \| $\sim 1.2\times 10^{3}$ \| | \| 1.5 \| 0.139 \| $7.33\times 10^{-3}$ \| $\sim 23\times$ \| $\sim 6.8\times$ \| | verificacion_numeros.md filas 15–17 |
| 155 | \| vs $\sim 10^{-4}$ (best case) \| | \| vs $\sim 3\times 10^{-2}$ (best case) \| | verificacion_numeros.md fila 15 |
| 158 | **still above** by $\sim 10^{2}$–$10^{3}$ \| | **still above** by $\sim 0.9$–$1.5\times$ \| | verificacion_numeros.md fila 15 (0.03/0.0339, 0.05/0.0339) |
| 164 | \| $\varepsilon$ \| $\\|\gamma-1\\|$ at $z=0.5$, $\sigma_X=2.5\times 10^{-2}$ \| | \| $\varepsilon$ \| $\\|\gamma-1\\|$ at $z=0.5$, $\sigma_X=2.633\times 10^{-2}$ \| | cota (tres cifras) |
| 166 | \| 1 \| $1.93\times 10^{-4}$ \| ⏎ \| 0.1 \| $1.93\times 10^{-5}$ \| ⏎ \| 0.01 \| $1.93\times 10^{-6}$ \| | \| 1 \| $3.39\times 10^{-2}$ \| ⏎ \| 0.1 \| $3.39\times 10^{-3}$ \| ⏎ \| 0.01 \| $3.39\times 10^{-4}$ \| | verificacion_numeros.md filas 15, 18 |
| 177 | too coarse to see $\lvert \eta-1 \rvert\sim 10^{-4}$. | too coarse to see $\lvert \eta-1 \rvert\sim 3\times 10^{-2}$. | verificacion_numeros.md fila 15 |
| 193 | (use $2.5\times 10^{-2}$ and $5\times 10^{-5}$ as bookends from `resume.txt`) | (use $2.6\times 10^{-2}$ and $\sigma_X\to 0$, the profile best fit, as bookends; `results/profile_sigma_x/profile_sigma_x.json`) | 2a l.93 |

### papers/data-pack-option0-internet.md

| línea | antes | después | fuente |
|---|---|---|---|
| 13 | \| $< 2.5\times 10^{-2}$ \| | \| $< 2.6\times 10^{-2}$ \| | cota |
| 83 | Using $\sigma_X \sim 10^{-4}$, | Using $\sigma_X = 2.6\times 10^{-2}$, | cota |
| 87 | \cdot\sigma_X\right)\sim 2\times 10^{-4}. | \cdot\sigma_X\right)\sim 5\times 10^{-2}. | verificacion_numeros.md fila 19 |
| 91 | \| vs $2\times 10^{-4}$ \| | \| vs $5\times 10^{-2}$ \| | verificacion_numeros.md fila 19 |
| 93 | \| ~500× larger than prediction \| | \| ~2× larger than prediction \| | verificacion_numeros.md fila 19 |
| 94 | \| $\sim 5\% = 0.05$ \| ~250× larger \| | \| $\sim 5\% = 0.05$ \| ~1× (comparable) \| | verificacion_numeros.md fila 19 |
| 95 | \| $\gtrsim 30\%$ \| still ≫ \| | \| $\gtrsim 30\%$ \| ~6× \| | 0.30/0.053 (misma tabla; l.95, fuera de 83–94) |

### papers/euclid-protocol-vacuum-relaxation.md

| línea | antes | después | fuente |
|---|---|---|---|
| 11 | extending the DESI DR2 working limit $\sigma_X < 2.5\times 10^{-2}$ | extending the DESI DR2 working limit $\sigma_X < 2.6\times 10^{-2}$ | cota |
| 81 | \| Spans Euclid targets and the DESI limit \| | \| Spans Euclid targets (the DESI limit $2.6\times 10^{-2}$ lies above the upper end) \| | 2a l.95 |
| 121 | \| **E1** \| $10^{-5}\lesssim \sigma_X \lesssim 2.5\times 10^{-2}$ | \| **E1** \| $10^{-5}\lesssim \sigma_X \lesssim 2.6\times 10^{-2}$ | cota |
| 123 | \| **E3** \| $\sigma_X \gtrsim 2.5\times 10^{-2}$ | \| **E3** \| $\sigma_X \gtrsim 2.6\times 10^{-2}$ | cota |
| 148 | \| $5\times 10^{-5}$ \| $2.5\times 10^{-2}$ \| $\sim 0.4$ \| | \| $5\times 10^{-5}$ \| $2.6\times 10^{-2}$ \| $\sim 2\times 10^{-3}$ \| | verificacion_numeros.md fila 9 |
| 247 | this pipeline): $\sigma_X < 2.5\times 10^{-2}$ (95% CL | this pipeline): $\sigma_X < 2.6\times 10^{-2}$ (95% CL | cota |

### papers/fundamental-vs-emergent-vacuum-relaxation.md

| línea | antes | después | fuente |
|---|---|---|---|
| 14 | null working limit $\sigma_X < 2.5\times 10^{-2}$ | null working limit $\sigma_X < 2.6\times 10^{-2}$ | cota: profile_sigma_x.json sigma_X_95CL_upper_profile |
| 86 | \| DESI working limit \| $\sigma_X<2.5\times 10^{-2}$ (95% CL) \| | \| DESI working limit \| $\sigma_X<2.6\times 10^{-2}$ (95% CL) \| | cota: profile_sigma_x.json sigma_X_95CL_upper_profile |
| 140 | \| **E1** \| Emergent + $A_0\in[10^{-5},2.5\times 10^{-2}]$, small $\theta$ \| Compatible if residual $<2.5\times 10^{-2}$ \| | \| **E1** \| Emergent + $A_0\in[10^{-5},2.6\times 10^{-2}]$, small $\theta$ \| Compatible if residual $<2.6\times 10^{-2}$ \| | cota: profile_sigma_x.json sigma_X_95CL_upper_profile |
| 142 | \| **E3** \| Small $\theta$ + $A_0>2.5\times 10^{-2}$ \| | \| **E3** \| Small $\theta$ + $A_0>2.6\times 10^{-2}$ \| | cota: profile_sigma_x.json sigma_X_95CL_upper_profile |

### papers/principle-of-vacuum-smoothness.md

| línea | antes | después | fuente |
|---|---|---|---|
| 19 | \sigma_X < 2.5\times 10^{-2} | \sigma_X < 2.6\times 10^{-2} | cota: profile_sigma_x.json sigma_X_95CL_upper_profile |

### papers/quantum-fluid-instabilities-desi-dr2.md

| línea | antes | después | fuente |
|---|---|---|---|
| 12 | 2.5\times 10^{-2}$ (exact grid $2.506\times 10^{-2}$) | 2.6\times 10^{-2}$ (root of the profile $2.633\times 10^{-2}$) | profile_sigma_x.json |
| 18 | is $\sigma_X < 2.5\times 10^{-2}$ (7 bins, full meas. cov.) [2] | is $\sigma_X < 2.6\times 10^{-2}$ (7 bins, full meas. cov.) [2] | cota: profile_sigma_x.json sigma_X_95CL_upper_profile |
| 158 | ceiling $\sigma_{X}<2.5\times 10^{-2}$ is | ceiling $\sigma_{X}<2.6\times 10^{-2}$ is | cota: profile_sigma_x.json sigma_X_95CL_upper_profile |
| 190 | Profile 95% CL $\sigma_X<2.5\times 10^{-2}$; | Profile 95% CL $\sigma_X<2.6\times 10^{-2}$; | cota: profile_sigma_x.json sigma_X_95CL_upper_profile |

### papers/scope-and-mixups.md

| línea | antes | después | fuente |
|---|---|---|---|
| 31 | \sigma_X < 2.5\times 10^{-2}\quad | \sigma_X < 2.6\times 10^{-2}\quad | cota: profile_sigma_x.json sigma_X_95CL_upper_profile |
| 49 | \| \(\sigma_X < 2.5\times 10^{-2}\) \| | \| \(\sigma_X < 2.6\times 10^{-2}\) \| | cota: profile_sigma_x.json sigma_X_95CL_upper_profile |
| 58 | the ceiling \(2.5\times 10^{-2}\) is | the ceiling \(2.6\times 10^{-2}\) is | cota: profile_sigma_x.json sigma_X_95CL_upper_profile |
| 58 | “baryons at amplitude 0.025” | “baryons at amplitude 0.026” | cota (decimal) |

### papers/sdiff-fundamental-vs-emergent.md

| línea | antes | después | fuente |
|---|---|---|---|
| 19 | current data ($\sigma_X < 2.5\times 10^{-2}$ from DESI DR2) | current data ($\sigma_X < 2.6\times 10^{-2}$ from DESI DR2) | cota: profile_sigma_x.json sigma_X_95CL_upper_profile |
| 72 | smoothness of the vacuum ($\sigma_X < 2.5\times 10^{-2}$) but | smoothness of the vacuum ($\sigma_X < 2.6\times 10^{-2}$) but | cota: profile_sigma_x.json sigma_X_95CL_upper_profile |
| 95 | smoothness of the vacuum ($\sigma_X < 2.5\times 10^{-2}$) without | smoothness of the vacuum ($\sigma_X < 2.6\times 10^{-2}$) without | cota: profile_sigma_x.json sigma_X_95CL_upper_profile |

### papers/unimodular-gravity-vacuum-smoothness.md

| línea | antes | después | fuente |
|---|---|---|---|
| 14 | establish a strict upper limit on stochastic fluctuations of the dark-energy density: $\sigma_X < 2.5\times 10^{-2}$ (95% CL) [1] | give a 95% CL profile upper limit (over $\theta$) on stochastic fluctuations of the dark-energy density: $\sigma_X < 2.6\times 10^{-2}$ [1] | cota; 2a l.109 |
| 32 | Maximum Likelihood Estimation drives the amplitude to zero: | The best fit drives the amplitude to zero: | 2a l.109 |
| 40 | and establishes a conservative upper limit: | and the profile likelihood over $\theta$ gives the 95% CL upper limit: | 2a l.109 |
| 44 | \sigma_X < 2.5\times 10^{-2} | \sigma_X < 2.6\times 10^{-2} | cota |
| 76 | \sigma_X < 2.5\times 10^{-2} | \sigma_X < 2.6\times 10^{-2} | cota |
| 82 | level of the residual bound (profile 95% CL $\sigma_X<2.5\times 10^{-2}$; | level of the residual bound (profile 95% CL $\sigma_X<2.6\times 10^{-2}$; | cota: profile_sigma_x.json sigma_X_95CL_upper_profile |
| 130 | m_\phi \lesssim 10^{-5}\ \text{eV}, \qquad \left\| \frac{V'}{V}\right\| M_{\mathrm{Pl}} \lesssim 10^{-2}. | m_\phi \lesssim 10^{-5}\ \text{eV}, \qquad \left\| \frac{V'}{V}\right\| M_{\mathrm{Pl}} \lesssim 10^{-2}. \quad \text{[valor no reproducible; pendiente de decisión del autor]} | m_φ: verificacion_numeros.md filas 6–8 |
| 296 | I have traced the empirical limit $\sigma_X < 2.5\times 10^{-2}$ | I have traced the empirical limit $\sigma_X < 2.6\times 10^{-2}$ | cota: profile_sigma_x.json sigma_X_95CL_upper_profile |
| 300 | the observational constraint that $\sigma_X$ is below $10^{-4}$; | the observational constraint that $\sigma_X$ is below $2.6\times 10^{-2}$ (95% CL); | cota; 2a l.111 |

### notes/desqueezing-relaxation-vacuum-fluctuations-note.md

| línea | antes | después | fuente |
|---|---|---|---|
| 31 | the DESI bound $\sigma_X<2.5\times 10^{-2}$ (95% CL) | the DESI bound $\sigma_X<2.6\times 10^{-2}$ (95% CL) | cota: profile_sigma_x.json sigma_X_95CL_upper_profile |
| 43 | \sigma_X < 2.5\times 10^{-2}, | \sigma_X < 2.6\times 10^{-2}, | cota |
| 233 | \| $\sigma_X$ (95% CL) \| $<2.5\times 10^{-2}$ \| | \| $\sigma_X$ (95% CL) \| $<2.6\times 10^{-2}$ \| | cota: profile_sigma_x.json sigma_X_95CL_upper_profile |
| 234 | \| MLE floor \| $\theta\sim 10^{-3}$, $\sigma_X\sim 5\times 10^{-5}$ \| numerical boundary \| | \| MLE floor (old run, diagonal cov) \| $\theta\sim 10^{-3}$, $\sigma_X\sim 5\times 10^{-5}$ \| numerical boundary \| | 2a l.112 |
| 327 | \ll 10^{-5}\ll 2.5\times 10^{-2}. | \ll 10^{-5}\ll 2.6\times 10^{-2}. | cota: profile_sigma_x.json sigma_X_95CL_upper_profile |
| 335 | The limit $\sigma_X<2.5\times 10^{-2}$ (95% CL) | The limit $\sigma_X<2.6\times 10^{-2}$ (95% CL) | cota: profile_sigma_x.json sigma_X_95CL_upper_profile |
| 393 | compare to DESI limit $2.5\times 10^{-2}$ \| | compare to DESI limit $2.6\times 10^{-2}$ \| | cota: profile_sigma_x.json sigma_X_95CL_upper_profile |
| 474 | \sigma_X^{\rm DESI}<2.5\times 10^{-2}, | \sigma_X^{\rm DESI}<2.6\times 10^{-2}, | cota: profile_sigma_x.json sigma_X_95CL_upper_profile |

## Afirmaciones que cambian de sentido

Solo se han cambiado números y palabras que describen números; la conclusión en prosa se deja tal cual y decide el autor en D. Todas vienen del slip: con σ_X = 2.633e-2, |γ−1| (ε = 1, δ_m = 1) es 3.39e-2 / 1.43e-2 / 7.33e-3 a z = 0.5 / 1 / 1.5 (`verificacion_numeros.md` filas 15–17). Queda:
- 5.0× / 12× / 23× por debajo de Maus (0.17);
- 1.5× / 3.5× / 6.8× por debajo del suelo de Sakr (0.05);
- 0.9–1.5× frente a Euclid-like 0.03–0.05 a z = 0.5.

Antes la distancia era 10²–10³.

| fichero:línea | frase actual | por qué cambia de sentido | redacción propuesta |
|---|---|---|---|
| manuscript/PREPRINT.md:28 | "…but inherits the same amplitude starvation: even full anisotropy of the BAO-bounded residual gives \|γ−1\|∼3×10⁻², of order 1.5–7× below current and forecast floors." | "Amplitude starvation" describía 10²–10³; ahora a z = 0.5 el mejor caso está 1.5× por debajo del suelo de forecast | "…but with the BAO-bounded residual it sits near forecast reach: full anisotropy gives \|γ−1\|∼3×10⁻² at z = 0.5, ≈5× below Maus et al. and ≈1.5× below the Sakr et al. constant-η floor, falling as ρ_X/ρ_m at higher z." |
| manuscript/PREPRINT.md:232 | "### 6.2 Option 0: anisotropic gap is real in principle, starved in amplitude" | ídem | "### 6.2 Option 0: anisotropic gap is real in principle, marginal in amplitude" |
| manuscript/PREPRINT.md:254 | "> Slip is the correct *kind* of operator…, but with the BAO-bounded residual amplitude it is not an observational shortcut. It inherits the amplification problem." | Con ε = 1 y δ_m = 1 no hace falta amplificación para quedar a 1.5× del suelo de forecast a z = 0.5 | "> Slip is the correct *kind* of operator…; with the BAO-bounded residual amplitude and ε = 1 it sits at the edge of forecast sensitivity (≈1.5× below at z = 0.5). It is a marginal channel, not a demonstrated shortcut: ε < 1, δ_m > 1 or z ≳ 1 push it below reach." |
| manuscript/PREPRINT.md:322 | "**Slip honesty:** … an amplitude-starved channel in practice." | ídem | "… a marginal channel in practice (best case ≈1.5× below forecast floors at z = 0.5)." |
| manuscript/README.md:34 | "- Gravitational-slip amplitude starvation (Option 0)" | ídem | "- Gravitational slip: marginal amplitude (Option 0)" |
| manuscript/CLAIMS.md:22 (C7) | "Even ε = 1 anisotropic fraction of the BAO-bounded residual gives \|γ−1\|≲3×10⁻² (z = 0.5, δ_m = 1)" | "Even" presenta el número como insuficiente; ahora está a 1.5× del suelo de forecast | "With ε = 1, the BAO-bounded residual gives \|γ−1\|≲3×10⁻² (z = 0.5, δ_m = 1): ≈5× below Maus et al., ≈1.5× below the Sakr et al. constant-η floor" |
| papers/resume.txt:16 | "(\|gamma-1\|~3e-2 best case; ~1.5-7x starved)" | "starved" para 1.5–7× | "(\|gamma-1\|~3e-2 best case; ~1.5-7x below floors at z=0.5)" |
| papers/anisotropic-slip-option0.md:25 | "…with σ_X≲2.6×10⁻² it inherits the same amplitude problem as the isotropic residual. It is an additional channel that turns on with amplification, not free amplification." | Con la cota vigente no hace falta amplificación para acercarse a la sensibilidad de forecast | "…with σ_X≲2.6×10⁻² the best case (ε = 1) reaches \|γ−1\|∼3×10⁻² at z = 0.5, at the edge of Euclid-like forecasts; its visibility depends on ε and δ_m rather than on amplification." |
| papers/anisotropic-slip-option0.md:41 | "…quantified as amplitude-starved via `slip_bridge.py`…" | ídem | "…quantified as marginal in amplitude via `slip_bridge.py` (best case ≈1.5× below forecast floors at z = 0.5)…" |
| papers/anisotropic-slip-option0.md:118 | "current sensitivity … is still O(0.1), far above the ∼3×10⁻² scale…" | O(0.1) frente a 3×10⁻² es un factor 3–5, no "far above" | "…is still O(0.1), a factor ∼3–5 above the ∼3×10⁻² scale…" |
| papers/anisotropic-slip-option0.md:157 | "\| Current DESI-era slip (Maus γ = 1.17 ± 0.11) \| O(0.1) \| **far above** prediction \|" | ídem | "\| … \| O(0.1) \| **3–5× above** prediction \|" |
| papers/anisotropic-slip-option0.md:158 | "\| Euclid-like forecasts (Sakr et al.) \| O(0.03–0.05) \| **still above** by ∼0.9–1.5× \|" | A 0.03 el forecast queda por debajo de la predicción (0.9×) | "\| … \| O(0.03–0.05) \| **comparable** (0.9–1.5×) \|" |
| papers/anisotropic-slip-option0.md:160 | "**Interpretation.** Even the most optimistic anisotropic fraction does not make the gap visible at current DESI precision, nor automatically at Euclid, unless there is amplification…" | La parte de Euclid deja de valer para ε ≈ 1, z ≈ 0.5 | "**Interpretation.** The most optimistic anisotropic fraction does not make the gap visible at current DESI precision (≈5× short at z = 0.5); at Euclid-like forecast precision it is marginal (0.9–1.5×) for ε ≈ 1, δ_m ≈ 1, and out of reach for ε ≲ 0.3 or z ≳ 1." |
| papers/anisotropic-slip-option0.md:162 | "…The anisotropic channel does not evade that problem; it inherits it." | El canal ya no necesita el factor 10⁵⁶ para quedar cerca de forecast; el problema de amplificación sigue siendo el de la semilla σ₀, no el de σ_X | "…The anisotropic channel does not evade the seed problem (σ₀ → σ_X); given σ_X at the DESI bound, it is marginal rather than starved." |
| papers/anisotropic-slip-option0.md:181 | "- That Euclid will “easily” see the gap without amplification." | Sigue siendo cierto ("easily"), pero ya no por mucho | Opcional: "- That Euclid will see the gap without amplification for ε < 1 or z ≳ 1 (best case ε = 1, z = 0.5 is marginal)." |
| papers/data-pack-option0-internet.md:96 | "**Conclusion (unchanged):** … the amplitude still requires amplification before a detection story makes sense." | La tabla de E da ahora ∼2× (Maus) y ∼1× (Sakr): no hace falta amplificación para el orden de magnitud | "**Conclusion (revised):** the correct operator exists and is measured or forecast; with σ_X at the DESI bound, the order-of-magnitude \|η−1\| is comparable to the Sakr et al. constant-η forecast floor." |
| papers/amplification-gap.md:83 | "Every closed channel that looked like a shortcut — anisotropic slip, … — **inherits this same gap**." | Vale para la semilla σ₀ → σ_X, pero no para slip una vez que σ_X está en la cota | "…— anisotropic slip, … — inherits the seed gap (σ₀ → σ_X); slip at the DESI bound itself is marginal, not starved (see anisotropic-slip-option0.md)." |
| papers/stochastic-dark-energy-desi-dr2.md:924 | "…the null result disappears. θ = 0.146, σ_X ≈ 9.2×10⁻³ — about 0.35 of the 95% CL limit 2.6×10⁻²…" | Antes decía "nearly 200 times larger than the reported limit", lo que implicaba exceder la cota; ahora la amplitud queda dentro. "The null result disappears" se refiere al best fit, pero l.928/930 ("a non-negligible amplitude reappears", "reintroduces an apparent signal") se leen como tensión con la cota | "…the best fit moves away from zero: θ = 0.146, σ_X ≈ 9.2×10⁻³, still below the 95% CL limit 2.6×10⁻²…", y en l.928/930 "a non-zero best-fit amplitude appears, within the bound". |

Sin cambio de sentido, aunque el número cambie mucho:
- stochastic l.887/899/906: la cota H.5 se endurece (1.8e-2 → 1.0e-4; λ ≳ 9.8e3), así que la conclusión de l.887 se refuerza.
- amplification-gap l.43 ("one decade" → "three decades") era un número.
- amplification-gap l.231 (regla BAO-safe) describe el escaneo.
- HONEST_HEADLINES l.27 (r 66 → 68).

## Comprobación: grep de los patrones de B1 punto 6

```
git grep -n -E "1\.5e-4|1\.5e-04|1\.5 ?×10⁻⁴|1\.5\\times ?10^{-4}|4\.44e7|4\.4e7|4\.4 ?×10|4\.44\\times|4\.4\\times|1\.27e59|1\.27 ?×10⁵⁹|1\.27\\times ?10^{59}|2\.506|2\.5e-2|2\.5e-02|2\.5 ?×10⁻²|2\.5\\times ?10^{-2}|0\.025\b" -- manuscript papers README.md notes
```

| fichero:línea | texto | por qué se queda |
|---|---|---|
| manuscript/CLAIMS.md:17 | "Former 1.5×10⁻⁴ was a working target **not** equal to this profile — corrected." | nota histórica; en la lista de obsoletos de AUDITORIA 2b ("CLAIMS 17") |
| papers/BOUND_CORRECTION_SIGMA_X.md:10 | columna "antes": σ_X < 1.5×10⁻⁴ | documento de corrección; en 2b ("BOUND_CORRECTION") |
| papers/BOUND_CORRECTION_SIGMA_X.md:11 | "…last grid node 2.506×10⁻², previous version" | nodo de rejilla marcado como versión anterior |
| papers/BOUND_CORRECTION_SIGMA_X.md:24 | "The former 1.5×10⁻⁴ was carried as a **working target**…" | nota histórica; en 2b ("BOUND_CORRECTION") |
| papers/scope-and-mixups.md:35 | "The old 1.5×10⁻⁴ was a working target, not this profile." | nota histórica; en 2b ("scope-and-mixups 35") |
| papers/stochastic-dark-energy-desi-dr2.md:621 | "…the limit σ_X < 2.5×10⁻² implies:" (derivación de m_φ) | excluido: fila m_φ, pendiente de decisión del autor (ver arriba) |

Fuera de estos patrones, siguen en el texto valores que pertenecen a C o D:
- Las definiciones de X y el factor 1/Ω_Λ0 (C).
- Las fechas y los bins de Euclid (D).
- La redacción de Everpresent y del Axioma A2 (D).
