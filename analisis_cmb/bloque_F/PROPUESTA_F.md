# Bloque F, paso 1: propuesta de texto (sin aplicar)

Fondo de referencia de los α, significado de C2, SN, picos acústicos, sirenas. **No se ha aplicado nada**: este fichero es lo único nuevo.

Fuentes de los números: `analisis_cmb/fondo/FONDO.md` y `perfil_fondo.json`, `analisis_cmb/sn/SN.md` y `perfil_sn.json`, `analisis_cmb/cmb/PICOS.md` y `picos.json`, `results/eos_cpl_desi_dr2/` y `scripts/eos_efectiva.py`, y el TeX de arXiv:2503.14738 (descargado el 2026-10-31; `main.tex`, l. 800 y l. 876) y el resumen de arXiv:2501.07477.

Las líneas son las de `manuscript/PREPRINT.md` y `manuscript/CLAIMS.md` en el commit 8bfe628. Los textos actuales se copian literalmente de los ficheros.

## Puntos a decidir

- **Cambio 4:** fundir las limitaciones 2 y 3 renumera la lista (8 → 7). Se da una alternativa sin renumerar.
- **Cambio 7:** las cifras de sirenas (10–25 %, 0.2–2 %) necesitan una referencia o hay que quitarlas.
- **Cambio 5:** solo el slip, el hueco y el Δχ² contra datos escalan linealmente; N_eff va como σ_X⁻² y la distancia en los picos como σ_X² (medida).
- **Cambio 8:** fecha de la v2.1.

## Índice

- 1. §2.1: la referencia de los α y el ajuste CPL — `manuscript/PREPRINT.md`, l. 82
- 2. §3.2: significado de la cota, tabla tratamiento × θ y SN — `manuscript/PREPRINT.md`, l. 165 (se inserta detrás)
- 3a. Resumen (l. 18–22): los dos números con su significado — `manuscript/PREPRINT.md`, l. 22
- 3b. CLAIMS C2 — `manuscript/CLAIMS.md`, l. 18
- 3c. Conclusión 1 — `manuscript/PREPRINT.md`, l. 369
- 4. Limitaciones 2 y 3 fundidas — `manuscript/PREPRINT.md`, l. 309–310
- 5a. Escalado al techo marginalizado: Δχ² de §3.5 contra datos (lineal) — `manuscript/PREPRINT.md`, l. 187 (se añade al final)
- 5b. Escalado: el hueco (lineal), §5.1 — `manuscript/PREPRINT.md`, l. 240 (se añade al final)
- 5c. Escalado: N_eff (∝ σ_X⁻²), §5.3 — `manuscript/PREPRINT.md`, l. 258
- 5d. Escalado: slip (lineal), §6.2 — `manuscript/PREPRINT.md`, l. 300 (se añade al final de la cita)
- 5e. Escalado: slip en C7 — `manuscript/CLAIMS.md`, l. 23
- 5f. Escalado: distancia a ΛCDM en los picos (cuadrática, medida) — `manuscript/PREPRINT.md`, l. —
- 6. §3.5: los picos acústicos — `manuscript/PREPRINT.md`, l. 191 (se inserta delante)
- 7. Trabajo futuro: sirenas estándar y si Euclid resuelve θ — `manuscript/PREPRINT.md`, l. 377
- 8. Historial de versiones: v2.1 — `manuscript/PREPRINT.md`, l. 384 (se inserta detrás)

## 1. §2.1: la referencia de los α y el ajuste CPL

**Fichero:** `manuscript/PREPRINT.md`  **Línea:** 82

**Texto actual:**

```text
**Important separation:** DESI analyses that combine BAO with CMB+SN often prefer dynamical $w_0$–$w_a$. The reference for $\alpha=1$ in my BAO-only residual analysis is the fiducial $\Lambda$CDM cosmology of the published $\alpha$ products ($w_0=-1$, $w_a=0$). As a check, a CPL fit to the same data lands at $w_0\approx -1.001$, $w_a\approx +0.014$ (`results/eos_cpl_desi_dr2/`, projected covariance, figure6 α), indistinguishable from that reference. The stochastic question is about *residuals*, not a re-fit of the full multi-probe DE equation of state.
```

**Texto propuesto:**

```text
**Important separation:** DESI analyses that combine BAO with CMB+SN often prefer dynamical $w_0$–$w_a$. The reference for $\alpha=1$ in my BAO-only residual analysis is the one of the published $\alpha$ products: the DESI DR2 best-fit flat $\Lambda$CDM ($\Omega_m = 0.2975$, $h r_\mathrm{d} = 101.54$ Mpc; arXiv:2503.14738, Sec. VI, eq. `constraints_lcdm_desi`), which I recover from the data themselves to $3\times 10^{-4}$ in α ([`analisis_cmb/fondo/FONDO.md`](../analisis_cmb/fondo/FONDO.md)); it is not the Planck cosmology. The CPL check of `results/eos_cpl_desi_dr2/` ($w_0\approx -1.001$, $w_a\approx +0.014$) computed its α with a Planck background ($\Omega_m = 0.315$, $H_0 = 67.4$) and with $\Omega_m$ and $h r_\mathrm{d}$ fixed, against data that are relative to the DESI best fit: at $w_0=-1$, $w_a=0$ it reproduces the reference by construction, so only a departure of shape could show up, and its closeness to $\Lambda$CDM is not an independent result. For comparison, DESI DR2 BAO alone with all parameters free gives $w_0=-0.48^{+0.35}_{-0.17}$, $w_a<-1.34$ (68 %), prior-dependent and only 1.7σ from $\Lambda$CDM (arXiv:2503.14738, Sec. VII A, eq. `w0wa_DESI`). The stochastic question is about *residuals*, not a re-fit of the full multi-probe DE equation of state.
```

**Nota:** Comprobado en el código: `scripts/eos_efectiva.py` usa `Om = 0.315`, `H0 = 67.4` y α_pred = D_V^CPL/D_V^ΛCDM con el mismo Ω_m, de modo que H0 y r_d se cancelan y α_pred(−1, 0) = 1. Valores de DESI copiados de su TeX (`main.tex`, l. 800 y l. 876; la Sec. VII A dice: "the posteriors in this poorly constrained case are cut off by the priors, so the marginalized means and limits quoted above are prior-dependent", y "a preference of just 1.7σ").

## 2. §3.2: significado de la cota, tabla tratamiento × θ y SN

**Fichero:** `manuscript/PREPRINT.md`  **Línea:** 165 (se inserta detrás)

**Texto actual:**

```text
From the residual analysis (profile / likelihood ratio under the OU kernel as implemented in the public pipeline): $\sigma_X < 2.6\times 10^{-2}$ (95% CL). The profile is shown in Fig. 1 ([`figures/profile_sigma_x.png`](../figures/profile_sigma_x.png)). At the limit the profile sits at the lower edge of the θ search (θ = 10⁻³, fully correlated bins); the θ → 0 supremum is 2.644 × 10⁻² (+0.4 %). Redefining the stationarity clock as cosmic time, conformal time or z instead of ln a moves the limit by less than 0.7 % ([`analisis_cmb/reloj/RELOJ.md`](../analisis_cmb/reloj/RELOJ.md)): it is set by the rank-1, fully correlated limit of the kernel and does not depend on the clock.
```

**Texto propuesto:**

```text
From the residual analysis (profile / likelihood ratio under the OU kernel as implemented in the public pipeline): $\sigma_X < 2.6\times 10^{-2}$ (95% CL). The profile is shown in Fig. 1 ([`figures/profile_sigma_x.png`](../figures/profile_sigma_x.png)). At the limit the profile sits at the lower edge of the θ search (θ = 10⁻³, fully correlated bins); the θ → 0 supremum is 2.644 × 10⁻² (+0.4 %). Redefining the stationarity clock as cosmic time, conformal time or z instead of ln a moves the limit by less than 0.7 % ([`analisis_cmb/reloj/RELOJ.md`](../analisis_cmb/reloj/RELOJ.md)): it is set by the rank-1, fully correlated limit of the kernel and does not depend on the clock.

**What the limit measures** ([`analisis_cmb/fondo/FONDO.md`](../analisis_cmb/fondo/FONDO.md)). The fully correlated mode (θ → 0) is, to first order, a change of $\Omega_m$: $S_i$ and $\partial\ln\alpha_i/\partial\Omega_m$ have cosine 0.99996. Because the α are referred to the DESI best fit, the data place this mode at $0.02\sigma$, and the limit is the $\ln\det$ (Occam) term of the profile, $\sqrt{e^{3.84}-1}\,\sigma_v = 6.75\,\sigma_v$ with $\sigma_v = (S^{\mathrm T}C^{-1}S)^{-1/2} = 0.0039$, i.e. $2.644\times 10^{-2}$ at θ → 0: consistency with the reference background, not a measurement of a residual. 95 % limits on $\sigma_X$ by treatment of the background ($\Omega_m$, $h r_\mathrm{d}$):

| background | θ free | θ = 0.001 | θ = 0.1 | θ = 1 | θ = 10 |
|:--|:--|:--|:--|:--|:--|
| fixed (DESI best fit) | $2.63\times 10^{-2}$ | $2.63\times 10^{-2}$ | $2.10\times 10^{-2}$ | $1.44\times 10^{-2}$ | $1.29\times 10^{-2}$ |
| profiled | $2.63\times 10^{-2}$ | $2.63\times 10^{-2}$ | $2.13\times 10^{-2}$ | $1.50\times 10^{-2}$ | $1.33\times 10^{-2}$ |
| marginalised (flat prior) | none | none | 0.140 | $4.6\times 10^{-2}$ | $2.2\times 10^{-2}$ |

The θ → 0 part of the limit disappears when the background is marginalised. The bound on a stochastic residual that varies between redshifts (θ ≳ 1) is $1.3$–$1.5\times 10^{-2}$ with the background fixed or profiled, and $2.2$–$4.6\times 10^{-2}$ marginalised. Pantheon+ supernovae with the same kernel bound $\sigma_X$ only at 0.06–0.12, 2.4–4.5 times weaker than DESI, with the same θ → 0 mode degenerate with $\Omega_m$ ([`analisis_cmb/sn/SN.md`](../analisis_cmb/sn/SN.md)).
```

**Nota:** La tabla y los números salen de `analisis_cmb/fondo/perfil_fondo.json` y `analisis_cmb/sn/perfil_sn.json`. 2.633 × 10⁻² es el valor en el borde θ = 10⁻³; 2.644 × 10⁻² (= 6.75 σ_v) es el supremo con θ → 0, ya citado en esta misma sección.

## 3a. Resumen (l. 18–22): los dos números con su significado

**Fichero:** `manuscript/PREPRINT.md`  **Línea:** 22

**Texto actual:**

```text
under the stated phenomenological kernel and public BAO summary statistics (official Gaussian BAO covariance projected to isotropic α in the baseline pipeline). Nested extensions that free $\{w_0,w_a,\sigma_X,\theta\}$ are not preferred over a smooth background ($\Delta\mathrm{AIC}\approx +4$).
```

**Texto propuesto:**

```text
under the stated phenomenological kernel and public BAO summary statistics (official Gaussian BAO covariance projected to isotropic α in the baseline pipeline). This number is reached in the fully correlated limit of the kernel, which is a shift of $\Omega_m$: it measures consistency with the DESI best-fit background and disappears if $\Omega_m$ and $h r_\mathrm{d}$ are marginalised. A residual that varies between redshifts (θ ≥ 1) is bounded at $\sigma_X < 1.5\times 10^{-2}$ with the background fixed or profiled, and $< 4.6\times 10^{-2}$ with it marginalised. Nested extensions that free $\{w_0,w_a,\sigma_X,\theta\}$ are not preferred over a smooth background ($\Delta\mathrm{AIC}\approx +4$).
```

**Nota:** Las líneas 18 y 20 (la ecuación $\sigma_X < 2.6\times 10^{-2}$) no cambian.

## 3b. CLAIMS C2

**Fichero:** `manuscript/CLAIMS.md`  **Línea:** 18

**Texto actual:**

```text
| **C2** | $\sigma_X < 2.6\times 10^{-2}$ (95% CL, profile over $\theta$, 7 isotropic $\alpha$ bins, full meas. cov). Best fit $\sigma_X\to 0$ (null preferred). Former $1.5\times 10^{-4}$ was a working target **not** equal to this profile — corrected. | `profile_sigma_x_desi.py` / `results/profile_sigma_x/`; robustness: `analisis_cmb/reloj/RELOJ.md` (stationarity clock) |
```

**Texto propuesto:**

```text
| **C2** | $\sigma_X < 2.6\times 10^{-2}$ (95% CL, profile over $\theta$, 7 isotropic $\alpha$ bins, full meas. cov), reached at θ → 0: a shift of $\Omega_m$ relative to the DESI best-fit background (ln det term, 6.75 σ_v; no bound with $\Omega_m$, $h r_\mathrm{d}$ marginalised). For a residual varying between redshifts (θ ≥ 1): $< 1.5\times 10^{-2}$ (background fixed or profiled), $< 4.6\times 10^{-2}$ (marginalised). Best fit $\sigma_X\to 0$ (null preferred). Former $1.5\times 10^{-4}$ was a working target **not** equal to this profile — corrected. | `profile_sigma_x_desi.py` / `results/profile_sigma_x/`; meaning: `analisis_cmb/fondo/FONDO.md`, `analisis_cmb/sn/SN.md`; robustness: `analisis_cmb/reloj/RELOJ.md` (stationarity clock) |
```

## 3c. Conclusión 1

**Fichero:** `manuscript/PREPRINT.md`  **Línea:** 369

**Texto actual:**

```text
1. **Null residual:** DESI DR2 BAO summary statistics do not favour stationary OU/QNM stochastic dark energy; $\sigma_X<2.6\times 10^{-2}$ (95% CL) under the stated kernel.
```

**Texto propuesto:**

```text
1. **Null residual:** DESI DR2 BAO summary statistics do not favour stationary OU/QNM stochastic dark energy; $\sigma_X<2.6\times 10^{-2}$ (95% CL) under the stated kernel. That limit is set by the fully correlated mode, i.e. consistency with the DESI best-fit background; for a residual varying between redshifts (θ ≥ 1), $\sigma_X<1.5\times 10^{-2}$ with the background fixed or profiled and $<4.6\times 10^{-2}$ marginalised (§3.2).
```

## 4. Limitaciones 2 y 3 fundidas

**Fichero:** `manuscript/PREPRINT.md`  **Línea:** 309–310

**Texto actual:**

```text
2. **Kernel dependence:** bounds are for the OU/QNM residual kernels stated; other residual models need their own covariance and MLE. Within the OU family the bound is insensitive to the stationarity clock (ln a, t, η, z within ±0.7 %) because it is reached in the fully correlated limit; DESI DR2 does not resolve θ and hence cannot select the clock.
3. **Background degeneracy:** with only seven BAO bins, flexible $w(z)$ can absorb residual variance (Section 3.3). Multi-probe fits are a separate analysis. The CMB check of §3.5 is a consistency test with parameters fixed except H0, not a multi-probe fit.
```

**Texto propuesto:**

```text
2. **Kernel dependence and the fully correlated limit:** bounds are for the OU/QNM residual kernels stated; other residual models need their own covariance and MLE. The headline limit is reached in the fully correlated limit (θ → 0), which is a background parameter, not a residual: a shift of $\Omega_m$ relative to the DESI best fit. That is why it does not depend on the stationarity clock (ln a, t, η, z within ±0.7 %), why DESI DR2 does not resolve θ, and why it vanishes when $\Omega_m$ and $h r_\mathrm{d}$ are marginalised (§3.2); with only seven BAO bins, flexible $w(z)$ can likewise absorb residual variance (Section 3.3). Multi-probe fits are a separate analysis. The CMB check of §3.5 is a consistency test with parameters fixed except H0, not a multi-probe fit.
```

**Nota:** Al fundir 2 y 3, las limitaciones 4–8 pasan a 3–7 ("Everpresent Λ test" pasa de 8 a 7). Ningún texto vivo cita las limitaciones por número; solo lo hacen registros de bloques anteriores (`COMPARACION_3f.md`, `COMPARACION_RELOJ_VALVULA.md`), que no se editan. Alternativa sin renumerar: dejar la 3 como "Background degeneracy" con solo la frase de w(z) y multi-sonda, y poner en la 2 lo del límite totalmente correlacionado.

## 5a. Escalado al techo marginalizado: Δχ² de §3.5 contra datos (lineal)

**Fichero:** `manuscript/PREPRINT.md`  **Línea:** 187 (se añade al final)

**Texto actual:**

```text
On a grid of ten values from 2.506×10⁻² down to 10⁻⁴ the 95th percentile stays below 4 everywhere and scales linearly with σ_X, so Planck does not bound σ_X below the DESI ceiling. The only CMB imprint is a scatter in H0, σ(H0) = 0.48 km/s/Mpc at θ = 1 and 0.020 at θ = 0.001, i.e. 12 and 280 times smaller than the Planck–SH0ES difference (5.68 km/s/Mpc): the residual does not explain the H0 tension.
```

**Texto propuesto:**

```text
On a grid of ten values from 2.506×10⁻² down to 10⁻⁴ the 95th percentile stays below 4 everywhere and scales linearly with σ_X, so Planck does not bound σ_X below the DESI ceiling. The only CMB imprint is a scatter in H0, σ(H0) = 0.48 km/s/Mpc at θ = 1 and 0.020 at θ = 0.001, i.e. 12 and 280 times smaller than the Planck–SH0ES difference (5.68 km/s/Mpc): the residual does not explain the H0 tension. At the marginalised θ = 1 ceiling of §3.2 ($4.611\times 10^{-2}$, ×1.751), this linear scaling gives a p95 of $0.067\to 0.117$ (TT, θ = 1): still invisible.
```

**Nota:** Extrapolado con el exponente 1 medido en la malla (que llega hasta 2.633 × 10⁻²); no recalculado.

## 5b. Escalado: el hueco (lineal), §5.1

**Fichero:** `manuscript/PREPRINT.md`  **Línea:** 240 (se añade al final)

**Texto actual:**

```text
The headline “$\sim 10^{56}$” always means the **Euclid-scale** target. Using the DESI ceiling in a density ratio and labelling it $10^{56}$ is incorrect by a factor $\sim 2\times 10^{3}$. Both gaps remain unbridgeable by soft linear amplifiers ($r\sim 64$–$68$).
```

**Texto propuesto:**

```text
The headline “$\sim 10^{56}$” always means the **Euclid-scale** target. Using the DESI ceiling in a density ratio and labelling it $10^{56}$ is incorrect by a factor $\sim 2\times 10^{3}$. Both gaps remain unbridgeable by soft linear amplifiers ($r\sim 64$–$68$). Against the marginalised θ = 1 ceiling ($4.611\times 10^{-2}$) the DESI gap becomes $G_{\mathrm{DESI}}=3.92\times 10^{59}$ (×1.751); the conclusion is unchanged.
```

**Nota:** 2.24 × 10⁵⁹ × 1.751 = 3.92 × 10⁵⁹ (equivalente a 4.611 × 10⁻²/1.18 × 10⁻⁶¹).

## 5c. Escalado: N_eff (∝ σ_X⁻²), §5.3

**Fichero:** `manuscript/PREPRINT.md`  **Línea:** 258

**Texto actual:**

```text
| **1. Local $N_{\mathrm{eff}}$** | $\sigma_{0,\mathrm{eff}}=1/\sqrt{N_{\mathrm{eff}}}$, $N_{\mathrm{eff}}\ll N_{\mathrm{BH}}$ | $N_{\mathrm{eff}}\sim 1.44\times 10^{3}$ (DESI), $1.00\times 10^{10}$ (Euclid). Only amplitude card—redefines the seed. |
```

**Texto propuesto:**

```text
| **1. Local $N_{\mathrm{eff}}$** | $\sigma_{0,\mathrm{eff}}=1/\sqrt{N_{\mathrm{eff}}}$, $N_{\mathrm{eff}}\ll N_{\mathrm{BH}}$ | $N_{\mathrm{eff}}\sim 1.44\times 10^{3}$ (DESI), $1.00\times 10^{10}$ (Euclid). Only amplitude card—redefines the seed. $N_{\mathrm{eff}}\propto\sigma_X^{-2}$: at the marginalised θ = 1 ceiling ($4.611\times 10^{-2}$), $\approx 470$ (DESI). |
```

**Nota:** 1.44 × 10³ × (2.633/4.611)² = 470. No es lineal en σ_X: σ_{0,eff} = 1/√N_eff.

## 5d. Escalado: slip (lineal), §6.2

**Fichero:** `manuscript/PREPRINT.md`  **Línea:** 300 (se añade al final de la cita)

**Texto actual:**

```text
> Slip is the correct *kind* of operator for anisotropic leakage past SDiff. With the BAO-bounded residual amplitude it is marginal, not starved: at ε = 1 it reaches $\lvert\gamma-1\rvert\approx 5\times 10^{-2}$ at $z=0.5$, level with the forecast floor ($\approx 1.0\times$) and $\approx 3.4\times$ below the current measurement. It is not a demonstrated shortcut. The σ_X → slip map is an order-of-magnitude estimate: it holds for δ_m ≈ 1 (quasi-linear scales) and does not extrapolate to δ_m ≪ 1; ε < 1 or z ≳ 1 push it below reach.
```

**Texto propuesto:**

```text
> Slip is the correct *kind* of operator for anisotropic leakage past SDiff. With the BAO-bounded residual amplitude it is marginal, not starved: at ε = 1 it reaches $\lvert\gamma-1\rvert\approx 5\times 10^{-2}$ at $z=0.5$, level with the forecast floor ($\approx 1.0\times$) and $\approx 3.4\times$ below the current measurement. It is not a demonstrated shortcut. The σ_X → slip map is an order-of-magnitude estimate: it holds for δ_m ≈ 1 (quasi-linear scales) and does not extrapolate to δ_m ≪ 1; ε < 1 or z ≳ 1 push it below reach. At the marginalised θ = 1 ceiling ($4.611\times 10^{-2}$, ×1.751) the slip at $z=0.5$ becomes $8.67\times 10^{-2}$, $\approx 1.7\times$ the Sakr et al. (Euclid-like) forecast floor: within reach.
```

**Nota:** 4.95 × 10⁻² × 1.751 = 8.67 × 10⁻²; 8.67/5 = 1.73. Sakr et al. (arXiv:2501.07477) son previsiones para "a Euclid-like photometric or a combined one with a DESI-like spectroscopic survey" (resumen; comprobado). El slip es lineal en σ_X (π_T = εσ_X ρ_crit,0).

## 5e. Escalado: slip en C7

**Fichero:** `manuscript/CLAIMS.md`  **Línea:** 23

**Texto actual:**

```text
| **C7** | With ε = 1, the BAO-bounded residual gives $\|\gamma-1\|\lesssim 5\times 10^{-2}$ ($z=0.5$, $\delta_m=1$): ≈3.4× below Maus et al. and ≈1.0× the Sakr et al. constant-η forecast floor (order-of-magnitude map; holds for δ_m ≈ 1, quasi-linear scales, and does not extrapolate to δ_m ≪ 1; ε < 1 or z ≳ 1 push it below reach) | `scripts/slip_bridge.py` |
```

**Texto propuesto:**

```text
| **C7** | With ε = 1, the BAO-bounded residual gives $\|\gamma-1\|\lesssim 5\times 10^{-2}$ ($z=0.5$, $\delta_m=1$): ≈3.4× below Maus et al. and ≈1.0× the Sakr et al. constant-η forecast floor (order-of-magnitude map; holds for δ_m ≈ 1, quasi-linear scales, and does not extrapolate to δ_m ≪ 1; ε < 1 or z ≳ 1 push it below reach). At the marginalised θ = 1 ceiling ($4.611\times 10^{-2}$): $8.67\times 10^{-2}$, ≈ 1.7× the floor | `scripts/slip_bridge.py` |
```

## 5f. Escalado: distancia a ΛCDM en los picos (cuadrática, medida)

**Fichero:** `manuscript/PREPRINT.md`  **Línea:** —

**Texto actual:**

```text
(ver el cambio 6)
```

**Texto propuesto:**

```text
(incluido en la frase del cambio 6: p95 3.5 × 10⁻³ → 1.1 × 10⁻², medido; ∝ σ_X²)
```

## 6. §3.5: los picos acústicos

**Fichero:** `manuscript/PREPRINT.md`  **Línea:** 191 (se inserta delante)

**Texto actual:**

```text
Figure and data: `analisis_cmb/RESULTADO.md`, `analisis_cmb/oscilador/`.
```

**Texto propuesto:**

```text
At high ℓ the residual is equally invisible: for the first 60 of these realizations, the Gaussian distance to $\Lambda$CDM in TT, TE and EE over ℓ ∈ [30, 2500] (cosmic variance plus Planck-like noise, $f_\mathrm{sky}=0.6$) has p95 $\Delta\chi^2 = 3.5\times 10^{-3}$ (θ = 1) and $5.3\times 10^{-6}$ (θ = 0.001); it comes from lensing in the acoustic troughs, not from peak positions (which θ_* fixes), and grows as $\sigma_X^2$: $1.1\times 10^{-2}$ at the marginalised θ = 1 ceiling ([`analisis_cmb/cmb/PICOS.md`](../analisis_cmb/cmb/PICOS.md)).

Figure and data: `analisis_cmb/RESULTADO.md`, `analisis_cmb/oscilador/`.
```

**Nota:** Ruido: blanco, 33 µK·arcmin y haz de 7′ (los N_ℓ de Planck 2018 I no están en el repositorio); por eso "Planck-like" y no "Planck".

## 7. Trabajo futuro: sirenas estándar y si Euclid resuelve θ

**Fichero:** `manuscript/PREPRINT.md`  **Línea:** 377

**Texto actual:**

```text
Future data (Euclid BAO multi-bin products; full DESI covariances) can tighten $\sigma_X$ with the same residual likelihood. The CMB gives no additional leverage at this amplitude (Section 3.5); BAO remains the route. Future theory work must pay the price of Routes 1 or 3 (or an inflation-style seed redefinition)—not invent linear gain on white noise.
```

**Texto propuesto:**

```text
Future data (Euclid BAO multi-bin products; full DESI covariances) can tighten $\sigma_X$ with the same residual likelihood. The CMB gives no additional leverage at this amplitude (Section 3.5); BAO remains the route. Future theory work must pay the price of Routes 1 or 3 (or an inflation-style seed redefinition)—not invent linear gain on white noise. Two questions follow from §3.2. Standard sirens measure $D_L$ without the $M_B$ calibration of supernovae; at present per-event precision (10–25 %) their effect on the bound would be at the 0.2–2 % level, and, like supernovae, they would not by themselves break the degeneracy of the fully correlated mode with $\Omega_m$. And a residual bound proper requires resolving θ: whether Euclid's BAO bins can separate a residual that varies between redshifts from a background shift is the measurement that would turn the headline limit into a bound on stochasticity.
```

**Nota:** **Fuente pendiente:** "10–25 %" y "0.2–2 %" son números del encargo, no calculados ni citados en el repositorio. Hace falta una referencia (por ejemplo, el H0 de GW170817 y de las sirenas oscuras de GWTC-3) o quitar las cifras. La segunda mitad de la frase sobre las sirenas viene de `SN.md`: el modo θ → 0 es degenerado con Ω_m, no con M_B, así que quitar M_B no basta.

## 8. Historial de versiones: v2.1

**Fichero:** `manuscript/PREPRINT.md`  **Línea:** 384 (se inserta detrás)

**Texto actual:**

```text
- **v2 (October 2026):** bound = root of the profile, 2.6×10⁻² (2.633×10⁻²); all derived numbers regenerated; CPL fits with the official DESI covariance; CMB consistency and harmonic residual (§3.5) and Everpresent Λ (§4.2) added; slip map corrected (π_T = εσ_X ρ_crit,0).
```

**Texto propuesto:**

```text
- **v2 (October 2026):** bound = root of the profile, 2.6×10⁻² (2.633×10⁻²); all derived numbers regenerated; CPL fits with the official DESI covariance; CMB consistency and harmonic residual (§3.5) and Everpresent Λ (§4.2) added; slip map corrected (π_T = εσ_X ρ_crit,0).
- **v2.1 (November 2026):** meaning of the σ_X limit stated: the θ → 0 part is the background mode relative to the DESI best fit (Occam term, no bound with the background marginalised), θ ≥ 1 limits by background treatment (§3.2, Pantheon+ cross-check); reference of the α named and the CPL check qualified (§2.1); high-ℓ CMB check (§3.5); derived numbers at the marginalised θ = 1 ceiling (§3.5, §5, §6.2); stationarity clock (§3.2), tanh map of Model 2 and convergence of Model 1 survival (§4.2), θ* per realization (§4.2). Before/after: block [F](../analisis_cmb/bloque_F/COMPARACION_F.md).
```

**Nota:** "November 2026" a confirmar. El enlace a `COMPARACION_F.md` existirá cuando se aplique el bloque (paso 2); hasta entonces, el comprobador de enlaces lo daría por roto.
