# Bloque F: antes / después

Aplica `PROPUESTA_F.md` con las decisiones del paso 2:

- Cambio 4: opción A (fundir las limitaciones 2 y 3 y renumerar).
- Cambio 7: sirenas con fuentes y frase de resolución de θ.
- Cambio 5: tal como estaba.
- Cambio 8: "November 2026".
- N14 nuevo, barrido del número de cabecera y entradas en `analisis_cmb/README.md`.

No se tocan `results/` ni `figures/`.

**Diferencias con lo pedido:**
- La cita de D_L de GW170817 es arXiv:1710.05832 (PRL 119, 161101), no 1710.05835: esa es la medida de H0 en Nature 551, 85, y su resumen no da ese D_L.
- El H0 = 68^{+8}_{−6} de GWTC-3 lleva su método (catálogo de galaxias, población BBH fija).
- Se añaden a Referencias GW170817, GWTC-3 y Pantheon+ (Brout et al. 2022, citado en §3.2).
- La v2.1 menciona además las sirenas y la resolución de θ (§9).

**Barrido del número de cabecera:**
- Tocados: `README.md`, `manuscript/README.md`, `manuscript/COVER_LETTER.md` y `papers/resume.txt` (la línea de la cota).
- No tocados: las líneas de `papers/resume.txt` con 2.633e-2 del chequeo del CMB, que no son la cota de cabecera.
- Ninguna nota viva citaba "limitation 8" por número; solo registros de bloques anteriores, que no se editan.

## Propuesta 1

Fichero: `manuscript/PREPRINT.md`

**Antes:**

```text
**Important separation:** DESI analyses that combine BAO with CMB+SN often prefer dynamical $w_0$–$w_a$. The reference for $\alpha=1$ in my BAO-only residual analysis is the fiducial $\Lambda$CDM cosmology of the published $\alpha$ products ($w_0=-1$, $w_a=0$). As a check, a CPL fit to the same data lands at $w_0\approx -1.001$, $w_a\approx +0.014$ (`results/eos_cpl_desi_dr2/`, projected covariance, figure6 α), indistinguishable from that reference. The stochastic question is about *residuals*, not a re-fit of the full multi-probe DE equation of state.
```

**Después:**

```text
**Important separation:** DESI analyses that combine BAO with CMB+SN often prefer dynamical $w_0$–$w_a$. The reference for $\alpha=1$ in my BAO-only residual analysis is the one of the published $\alpha$ products: the DESI DR2 best-fit flat $\Lambda$CDM ($\Omega_m = 0.2975$, $h r_\mathrm{d} = 101.54$ Mpc; arXiv:2503.14738, Sec. VI, eq. `constraints_lcdm_desi`), which I recover from the data themselves to $3\times 10^{-4}$ in α ([`analisis_cmb/fondo/FONDO.md`](../analisis_cmb/fondo/FONDO.md)); it is not the Planck cosmology. The CPL check of `results/eos_cpl_desi_dr2/` ($w_0\approx -1.001$, $w_a\approx +0.014$) computed its α with a Planck background ($\Omega_m = 0.315$, $H_0 = 67.4$) and with $\Omega_m$ and $h r_\mathrm{d}$ fixed, against data that are relative to the DESI best fit: at $w_0=-1$, $w_a=0$ it reproduces the reference by construction, so only a departure of shape could show up, and its closeness to $\Lambda$CDM is not an independent result. For comparison, DESI DR2 BAO alone with all parameters free gives $w_0=-0.48^{+0.35}_{-0.17}$, $w_a<-1.34$ (68 %), prior-dependent and only 1.7σ from $\Lambda$CDM (arXiv:2503.14738, Sec. VII A, eq. `w0wa_DESI`). The stochastic question is about *residuals*, not a re-fit of the full multi-probe DE equation of state.
```

## Propuesta 2

Fichero: `manuscript/PREPRINT.md`

**Antes:**

```text
From the residual analysis (profile / likelihood ratio under the OU kernel as implemented in the public pipeline): $\sigma_X < 2.6\times 10^{-2}$ (95% CL). The profile is shown in Fig. 1 ([`figures/profile_sigma_x.png`](../figures/profile_sigma_x.png)). At the limit the profile sits at the lower edge of the θ search (θ = 10⁻³, fully correlated bins); the θ → 0 supremum is 2.644 × 10⁻² (+0.4 %). Redefining the stationarity clock as cosmic time, conformal time or z instead of ln a moves the limit by less than 0.7 % ([`analisis_cmb/reloj/RELOJ.md`](../analisis_cmb/reloj/RELOJ.md)): it is set by the rank-1, fully correlated limit of the kernel and does not depend on the clock.
```

**Después:**

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

## Propuesta 3a

Fichero: `manuscript/PREPRINT.md`

**Antes:**

```text
under the stated phenomenological kernel and public BAO summary statistics (official Gaussian BAO covariance projected to isotropic α in the baseline pipeline). Nested extensions that free $\{w_0,w_a,\sigma_X,\theta\}$ are not preferred over a smooth background ($\Delta\mathrm{AIC}\approx +4$).
```

**Después:**

```text
under the stated phenomenological kernel and public BAO summary statistics (official Gaussian BAO covariance projected to isotropic α in the baseline pipeline). This number is reached in the fully correlated limit of the kernel, which is a shift of $\Omega_m$: it measures consistency with the DESI best-fit background and disappears if $\Omega_m$ and $h r_\mathrm{d}$ are marginalised. A residual that varies between redshifts (θ ≥ 1) is bounded at $\sigma_X < 1.5\times 10^{-2}$ with the background fixed or profiled, and $< 4.6\times 10^{-2}$ with it marginalised. Nested extensions that free $\{w_0,w_a,\sigma_X,\theta\}$ are not preferred over a smooth background ($\Delta\mathrm{AIC}\approx +4$).
```

## Propuesta 3b

Fichero: `manuscript/CLAIMS.md`

**Antes:**

```text
| **C2** | $\sigma_X < 2.6\times 10^{-2}$ (95% CL, profile over $\theta$, 7 isotropic $\alpha$ bins, full meas. cov). Best fit $\sigma_X\to 0$ (null preferred). Former $1.5\times 10^{-4}$ was a working target **not** equal to this profile — corrected. | `profile_sigma_x_desi.py` / `results/profile_sigma_x/`; robustness: `analisis_cmb/reloj/RELOJ.md` (stationarity clock) |
```

**Después:**

```text
| **C2** | $\sigma_X < 2.6\times 10^{-2}$ (95% CL, profile over $\theta$, 7 isotropic $\alpha$ bins, full meas. cov), reached at θ → 0: a shift of $\Omega_m$ relative to the DESI best-fit background (ln det term, 6.75 σ_v; no bound with $\Omega_m$, $h r_\mathrm{d}$ marginalised). For a residual varying between redshifts (θ ≥ 1): $< 1.5\times 10^{-2}$ (background fixed or profiled), $< 4.6\times 10^{-2}$ (marginalised). Best fit $\sigma_X\to 0$ (null preferred). Former $1.5\times 10^{-4}$ was a working target **not** equal to this profile — corrected. | `profile_sigma_x_desi.py` / `results/profile_sigma_x/`; meaning: `analisis_cmb/fondo/FONDO.md`, `analisis_cmb/sn/SN.md`; robustness: `analisis_cmb/reloj/RELOJ.md` (stationarity clock) |
```

## Propuesta 3c

Fichero: `manuscript/PREPRINT.md`

**Antes:**

```text
1. **Null residual:** DESI DR2 BAO summary statistics do not favour stationary OU/QNM stochastic dark energy; $\sigma_X<2.6\times 10^{-2}$ (95% CL) under the stated kernel.
```

**Después:**

```text
1. **Null residual:** DESI DR2 BAO summary statistics do not favour stationary OU/QNM stochastic dark energy; $\sigma_X<2.6\times 10^{-2}$ (95% CL) under the stated kernel. That limit is set by the fully correlated mode, i.e. consistency with the DESI best-fit background; for a residual varying between redshifts (θ ≥ 1), $\sigma_X<1.5\times 10^{-2}$ with the background fixed or profiled and $<4.6\times 10^{-2}$ marginalised (§3.2).
```

## Propuesta 4

Fichero: `manuscript/PREPRINT.md`

**Antes:**

```text
2. **Kernel dependence:** bounds are for the OU/QNM residual kernels stated; other residual models need their own covariance and MLE. Within the OU family the bound is insensitive to the stationarity clock (ln a, t, η, z within ±0.7 %) because it is reached in the fully correlated limit; DESI DR2 does not resolve θ and hence cannot select the clock.
3. **Background degeneracy:** with only seven BAO bins, flexible $w(z)$ can absorb residual variance (Section 3.3). Multi-probe fits are a separate analysis. The CMB check of §3.5 is a consistency test with parameters fixed except H0, not a multi-probe fit.
```

**Después:**

```text
2. **Kernel dependence and the fully correlated limit:** bounds are for the OU/QNM residual kernels stated; other residual models need their own covariance and MLE. The headline limit is reached in the fully correlated limit (θ → 0), which is a background parameter, not a residual: a shift of $\Omega_m$ relative to the DESI best fit. That is why it does not depend on the stationarity clock (ln a, t, η, z within ±0.7 %), why DESI DR2 does not resolve θ, and why it vanishes when $\Omega_m$ and $h r_\mathrm{d}$ are marginalised (§3.2); with only seven BAO bins, flexible $w(z)$ can likewise absorb residual variance (Section 3.3). Multi-probe fits are a separate analysis. The CMB check of §3.5 is a consistency test with parameters fixed except H0, not a multi-probe fit.
```

## Propuesta 5a

Fichero: `manuscript/PREPRINT.md`

**Antes:**

```text
On a grid of ten values from 2.506×10⁻² down to 10⁻⁴ the 95th percentile stays below 4 everywhere and scales linearly with σ_X, so Planck does not bound σ_X below the DESI ceiling. The only CMB imprint is a scatter in H0, σ(H0) = 0.48 km/s/Mpc at θ = 1 and 0.020 at θ = 0.001, i.e. 12 and 280 times smaller than the Planck–SH0ES difference (5.68 km/s/Mpc): the residual does not explain the H0 tension.
```

**Después:**

```text
On a grid of ten values from 2.506×10⁻² down to 10⁻⁴ the 95th percentile stays below 4 everywhere and scales linearly with σ_X, so Planck does not bound σ_X below the DESI ceiling. The only CMB imprint is a scatter in H0, σ(H0) = 0.48 km/s/Mpc at θ = 1 and 0.020 at θ = 0.001, i.e. 12 and 280 times smaller than the Planck–SH0ES difference (5.68 km/s/Mpc): the residual does not explain the H0 tension. At the marginalised θ = 1 ceiling of §3.2 ($4.611\times 10^{-2}$, ×1.751), this linear scaling gives a p95 of $0.067\to 0.117$ (TT, θ = 1): still invisible.
```

## Propuesta 5b

Fichero: `manuscript/PREPRINT.md`

**Antes:**

```text
The headline “$\sim 10^{56}$” always means the **Euclid-scale** target. Using the DESI ceiling in a density ratio and labelling it $10^{56}$ is incorrect by a factor $\sim 2\times 10^{3}$. Both gaps remain unbridgeable by soft linear amplifiers ($r\sim 64$–$68$).
```

**Después:**

```text
The headline “$\sim 10^{56}$” always means the **Euclid-scale** target. Using the DESI ceiling in a density ratio and labelling it $10^{56}$ is incorrect by a factor $\sim 2\times 10^{3}$. Both gaps remain unbridgeable by soft linear amplifiers ($r\sim 64$–$68$). Against the marginalised θ = 1 ceiling ($4.611\times 10^{-2}$) the DESI gap becomes $G_{\mathrm{DESI}}=3.92\times 10^{59}$ (×1.751); the conclusion is unchanged.
```

## Propuesta 5c

Fichero: `manuscript/PREPRINT.md`

**Antes:**

```text
| **1. Local $N_{\mathrm{eff}}$** | $\sigma_{0,\mathrm{eff}}=1/\sqrt{N_{\mathrm{eff}}}$, $N_{\mathrm{eff}}\ll N_{\mathrm{BH}}$ | $N_{\mathrm{eff}}\sim 1.44\times 10^{3}$ (DESI), $1.00\times 10^{10}$ (Euclid). Only amplitude card—redefines the seed. |
```

**Después:**

```text
| **1. Local $N_{\mathrm{eff}}$** | $\sigma_{0,\mathrm{eff}}=1/\sqrt{N_{\mathrm{eff}}}$, $N_{\mathrm{eff}}\ll N_{\mathrm{BH}}$ | $N_{\mathrm{eff}}\sim 1.44\times 10^{3}$ (DESI), $1.00\times 10^{10}$ (Euclid). Only amplitude card—redefines the seed. $N_{\mathrm{eff}}\propto\sigma_X^{-2}$: at the marginalised θ = 1 ceiling ($4.611\times 10^{-2}$), $\approx 470$ (DESI). |
```

## Propuesta 5d

Fichero: `manuscript/PREPRINT.md`

**Antes:**

```text
> Slip is the correct *kind* of operator for anisotropic leakage past SDiff. With the BAO-bounded residual amplitude it is marginal, not starved: at ε = 1 it reaches $\lvert\gamma-1\rvert\approx 5\times 10^{-2}$ at $z=0.5$, level with the forecast floor ($\approx 1.0\times$) and $\approx 3.4\times$ below the current measurement. It is not a demonstrated shortcut. The σ_X → slip map is an order-of-magnitude estimate: it holds for δ_m ≈ 1 (quasi-linear scales) and does not extrapolate to δ_m ≪ 1; ε < 1 or z ≳ 1 push it below reach.
```

**Después:**

```text
> Slip is the correct *kind* of operator for anisotropic leakage past SDiff. With the BAO-bounded residual amplitude it is marginal, not starved: at ε = 1 it reaches $\lvert\gamma-1\rvert\approx 5\times 10^{-2}$ at $z=0.5$, level with the forecast floor ($\approx 1.0\times$) and $\approx 3.4\times$ below the current measurement. It is not a demonstrated shortcut. The σ_X → slip map is an order-of-magnitude estimate: it holds for δ_m ≈ 1 (quasi-linear scales) and does not extrapolate to δ_m ≪ 1; ε < 1 or z ≳ 1 push it below reach. At the marginalised θ = 1 ceiling ($4.611\times 10^{-2}$, ×1.751) the slip at $z=0.5$ becomes $8.67\times 10^{-2}$, $\approx 1.7\times$ the Sakr et al. (Euclid-like) forecast floor: within reach.
```

## Propuesta 5e

Fichero: `manuscript/CLAIMS.md`

**Antes:**

```text
| **C7** | With ε = 1, the BAO-bounded residual gives $\|\gamma-1\|\lesssim 5\times 10^{-2}$ ($z=0.5$, $\delta_m=1$): ≈3.4× below Maus et al. and ≈1.0× the Sakr et al. constant-η forecast floor (order-of-magnitude map; holds for δ_m ≈ 1, quasi-linear scales, and does not extrapolate to δ_m ≪ 1; ε < 1 or z ≳ 1 push it below reach) | `scripts/slip_bridge.py` |
```

**Después:**

```text
| **C7** | With ε = 1, the BAO-bounded residual gives $\|\gamma-1\|\lesssim 5\times 10^{-2}$ ($z=0.5$, $\delta_m=1$): ≈3.4× below Maus et al. and ≈1.0× the Sakr et al. constant-η forecast floor (order-of-magnitude map; holds for δ_m ≈ 1, quasi-linear scales, and does not extrapolate to δ_m ≪ 1; ε < 1 or z ≳ 1 push it below reach). At the marginalised θ = 1 ceiling ($4.611\times 10^{-2}$): $8.67\times 10^{-2}$, ≈ 1.7× the floor | `scripts/slip_bridge.py` |
```

## Propuesta 6

Fichero: `manuscript/PREPRINT.md`

**Antes:**

```text
Figure and data: `analisis_cmb/RESULTADO.md`, `analisis_cmb/oscilador/`.
```

**Después:**

```text
At high ℓ the residual is equally invisible: for the first 60 of these realizations, the Gaussian distance to $\Lambda$CDM in TT, TE and EE over ℓ ∈ [30, 2500] (cosmic variance plus Planck-like noise, $f_\mathrm{sky}=0.6$) has p95 $\Delta\chi^2 = 3.5\times 10^{-3}$ (θ = 1) and $5.3\times 10^{-6}$ (θ = 0.001); it comes from lensing in the acoustic troughs, not from peak positions (which θ_* fixes), and grows as $\sigma_X^2$: $1.1\times 10^{-2}$ at the marginalised θ = 1 ceiling ([`analisis_cmb/cmb/PICOS.md`](../analisis_cmb/cmb/PICOS.md)).

Figure and data: `analisis_cmb/RESULTADO.md`, `analisis_cmb/oscilador/`.
```

## Propuesta 8

Fichero: `manuscript/PREPRINT.md`

**Antes:**

```text
- **v2 (October 2026):** bound = root of the profile, 2.6×10⁻² (2.633×10⁻²); all derived numbers regenerated; CPL fits with the official DESI covariance; CMB consistency and harmonic residual (§3.5) and Everpresent Λ (§4.2) added; slip map corrected (π_T = εσ_X ρ_crit,0).
```

**Después:**

```text
- **v2 (October 2026):** bound = root of the profile, 2.6×10⁻² (2.633×10⁻²); all derived numbers regenerated; CPL fits with the official DESI covariance; CMB consistency and harmonic residual (§3.5) and Everpresent Λ (§4.2) added; slip map corrected (π_T = εσ_X ρ_crit,0).
- **v2.1 (November 2026):** meaning of the σ_X limit stated: the θ → 0 part is the background mode relative to the DESI best fit (Occam term, no bound with the background marginalised), θ ≥ 1 limits by background treatment (§3.2, Pantheon+ cross-check); reference of the α named and the CPL check qualified (§2.1); high-ℓ CMB check (§3.5); derived numbers at the marginalised θ = 1 ceiling (§3.5, §5, §6.2); stationarity clock (§3.2), tanh map of Model 2 and convergence of Model 1 survival (§4.2), θ* per realization (§4.2); outlook on standard sirens and on resolving θ (§9). Before/after: block [F](../analisis_cmb/bloque_F/COMPARACION_F.md).
```

Nota: v2.1 añade "outlook on standard sirens and on resolving θ (§9)", por el cambio 7.

## Cambio 4: renumeración 4 → 3

Fichero: `manuscript/PREPRINT.md`

**Antes:**

```text
4. **Cross-correlation
```

**Después:**

```text
3. **Cross-correlation
```

## Cambio 4: renumeración 5 → 4

Fichero: `manuscript/PREPRINT.md`

**Antes:**

```text
5. **Sorkin seed
```

**Después:**

```text
4. **Sorkin seed
```

## Cambio 4: renumeración 6 → 5

Fichero: `manuscript/PREPRINT.md`

**Antes:**

```text
6. **Independent research
```

**Después:**

```text
5. **Independent research
```

## Cambio 4: renumeración 7 → 6

Fichero: `manuscript/PREPRINT.md`

**Antes:**

```text
7. **CMB check
```

**Después:**

```text
6. **CMB check
```

## Cambio 4: renumeración 8 → 7

Fichero: `manuscript/PREPRINT.md`

**Antes:**

```text
8. **Everpresent Λ test
```

**Después:**

```text
7. **Everpresent Λ test
```

## Cambio 7 (sirenas con fuentes, resolución de θ)

Fichero: `manuscript/PREPRINT.md`

**Antes:**

```text
Future data (Euclid BAO multi-bin products; full DESI covariances) can tighten $\sigma_X$ with the same residual likelihood. The CMB gives no additional leverage at this amplitude (Section 3.5); BAO remains the route. Future theory work must pay the price of Routes 1 or 3 (or an inflation-style seed redefinition)—not invent linear gain on white noise.
```

**Después:**

```text
Future data (Euclid BAO multi-bin products; full DESI covariances) can tighten $\sigma_X$ with the same residual likelihood. The CMB gives no additional leverage at this amplitude (Section 3.5); BAO remains the route. Future theory work must pay the price of Routes 1 or 3 (or an inflation-style seed redefinition)—not invent linear gain on white noise. Two questions follow from §3.2. Standard sirens measure $D_L$ without the $M_B$ calibration of supernovae, but their present precision is far from the effect at the limit: GW170817 gave $D_L = 40^{+8}_{-14}$ Mpc (Abbott et al. 2017, PRL 119, 161101, arXiv:1710.05832), and the GWTC-3 dark sirens $H_0 = 68^{+8}_{-6}$ km/s/Mpc with the galaxy-catalogue method and a fixed BBH population (Abbott et al. 2023, ApJ 949, 76, arXiv:2111.03604), i.e. ≈ 10–35 % per measurement, whereas the fully correlated residual at the limit changes $D_L$ by 0.2 % at $z = 0.1$ to 2.2 % at $z = 2$ ($S(z)\times 2.633\times 10^{-2}$; [`analisis_cmb/sn/SN.md`](../analisis_cmb/sn/SN.md)); like supernovae, they would not by themselves break the degeneracy of that mode with $\Omega_m$. Turning the headline limit into a bound on stochasticity requires resolving θ, and that needs a residual near the present limits: with DESI DR2 and the background fixed, separating θ = 1 from θ = 10⁻³ needs $\sigma_X \geq 2.2\times 10^{-2}$, above the θ = 1 limit of $1.44\times 10^{-2}$; with the background marginalised a narrow window opens ($2.8$–$4.6\times 10^{-2}$ if the truth is θ = 1); a Euclid-like survey (13 bins at 0.5 %) needs $\sigma_X \geq 1.5\times 10^{-2}$, i.e. it resolves θ only if the residual sits at the present ceiling ([`analisis_cmb/fondo/RESOLVER_THETA.md`](../analisis_cmb/fondo/RESOLVER_THETA.md)).
```

Nota: Fuentes comprobadas en arXiv y por DOI. El D_L = 40^{+8}_{−14} Mpc es de arXiv:1710.05832 (PRL 119, 161101), no de 1710.05835 (Nature 551, 85, la medida de H0). El H0 = 68^{+8}_{−6} de GWTC-3 es el del método de catálogo con población BBH fija; el resumen avisa que depende de la distribución de masas. Los porcentajes de D_L salen de SN.md (S de perfil_sn.py × 2.633e-2).

## Referencias: sirenas y Pantheon+

Fichero: `manuscript/PREPRINT.md`

**Antes:**

```text
9. Morales Souhail, J., this repository and companion notes (October 2026).
```

**Después:**

```text
9. Morales Souhail, J., this repository and companion notes (October 2026).
10. Abbott, B. P., et al. (LIGO Scientific and Virgo Collaborations), "GW170817: Observation of Gravitational Waves from a Binary Neutron Star Inspiral," Phys. Rev. Lett. 119, 161101 (2017), [arXiv:1710.05832](https://arxiv.org/abs/1710.05832).
11. Abbott, R., et al. (LIGO, Virgo and KAGRA Collaborations), "Constraints on the Cosmic Expansion History from GWTC-3," Astrophys. J. 949, 76 (2023), [arXiv:2111.03604](https://arxiv.org/abs/2111.03604).
12. Brout, D., et al., "The Pantheon+ Analysis: Cosmological Constraints," Astrophys. J. 938, 110 (2022), [arXiv:2202.04077](https://arxiv.org/abs/2202.04077).
```

Nota: Revista, volumen y página comprobados por DOI (10.1103/PhysRevLett.119.161101, 10.3847/1538-4357/ac74bb, 10.3847/1538-4357/ac8e04).

## N14

Fichero: `manuscript/CLAIMS.md`

**Antes:**

```text
| **N13** | No full CMB MCMC: parameters other than H0 are fixed to Planck 2018 in the CMB check |
```

**Después:**

```text
| **N13** | No full CMB MCMC: parameters other than H0 are fixed to Planck 2018 in the CMB check |
| **N14** | The $2.6\times 10^{-2}$ profile limit is not a bound on stochastic variance (it is the fully correlated mode, a background shift); the θ ≥ 1 limits of §3.2 are. |
```

## Barrido: README.md

Fichero: `README.md`

**Antes:**

```text
Best fit still goes to \(\sigma_X\to 0\).
```

**Después:**

```text
This number is reached in the fully correlated limit (\(\theta\to 0\)), a shift of \(\Omega_m\) relative to the DESI best-fit background, not a bound on stochastic variance; for a residual that varies between redshifts (\(\theta\ge 1\)) the limits are \(1.3\)–\(1.5\times 10^{-2}\) (background fixed or profiled) and \(2.2\)–\(4.6\times 10^{-2}\) (marginalised) — see [`manuscript/PREPRINT.md`](manuscript/PREPRINT.md) §3.2. Best fit still goes to \(\sigma_X\to 0\).
```

## Barrido: manuscript/README.md

Fichero: `manuscript/README.md`

**Antes:**

```text
- $\sigma_X < 2.6\times 10^{-2}$ (95% CL) under the stated kernel
```

**Después:**

```text
- $\sigma_X < 2.6\times 10^{-2}$ (95% CL) under the stated kernel — set by the fully correlated mode, a background shift; the θ ≥ 1 limits on a residual that varies with redshift are in §3.2
```

## Barrido: COVER_LETTER.md

Fichero: `manuscript/COVER_LETTER.md`

**Antes:**

```text
($\sigma_X < 2.6\times 10^{-2}$ (95% CL) under a stated OU kernel)
```

**Después:**

```text
($\sigma_X < 2.6\times 10^{-2}$ (95% CL) under a stated OU kernel, set by the fully correlated, background-like mode; §3.2 gives the θ ≥ 1 limits on a residual that varies with redshift)
```

## Barrido: papers/resume.txt

Fichero: `papers/resume.txt`

**Antes:**

```text
- 95% CL upper limit: σ_X < 2.6 × 10⁻² (profile over θ, background fixed at α = 1)
```

**Después:**

```text
- 95% CL upper limit: σ_X < 2.6 × 10⁻² (profile over θ, background fixed at α = 1)
  This is the fully correlated mode (theta -> 0), a shift of Omega_m relative to the DESI best fit, not a bound on
  stochastic variance; theta >= 1: 1.3-1.5e-2 (background fixed or profiled), 2.2-4.6e-2 (marginalised).
  PREPRINT §3.2; analisis_cmb/fondo/FONDO.md.
```

## analisis_cmb/README.md: fondo/, sn/, cmb/

Fichero: `analisis_cmb/README.md`

**Antes:**

```text
| `reloj/` | σ_X profile with four stationarity clocks (ln a, cosmic time, conformal time, z): `perfil_reloj.py`, `perfil_reloj.{json,txt}`, `RELOJ.md`; before/after of the related text in `COMPARACION_RELOJ_VALVULA.md` |
```

**Después:**

```text
| `reloj/` | σ_X profile with four stationarity clocks (ln a, cosmic time, conformal time, z): `perfil_reloj.py`, `perfil_reloj.{json,txt}`, `RELOJ.md`; before/after of the related text in `COMPARACION_RELOJ_VALVULA.md` |
| `fondo/` | The DESI σ_X profile with the background free (`FONDO.md`, `perfil_fondo.py`): the α reference is the DESI DR2 best fit, the θ → 0 limit is the Ω_m direction (ln det term), θ ≥ 1 limits by background treatment; σ_X needed to resolve θ with DESI and a Euclid-like survey (`RESOLVER_THETA.md`, `resolver_theta.py`) |
| `sn/` | σ_X profile with Pantheon+ (`SN.md`, `perfil_sn.py`, data in `sn/datos/` with sha256); effect of the limit on D_L for standard sirens |
| `cmb/` | High-ℓ CMB check (`PICOS.md`, `picos.py`): TT, TE, EE distance to ΛCDM over ℓ ∈ [30, 2500], lensing in the troughs, ∝ σ_X² |
```

## SN.md: efecto en D_L a la cota

Fichero: `analisis_cmb/sn/SN.md`

**Antes:**

```text
## Predicción
```

**Después:**

```text
## Efecto en D_L a la cota (para las sirenas estándar)

El residuo totalmente correlacionado a la cota cambia ln D_L en S(z) × 2.633 × 10⁻², con la S de
`perfil_sn.py` (−∂ln D_L/∂Ω_Λ, plano, H0 fijo, Ω_Λ = 0.6847):

| z | S(z) | δln D_L a 2.633 × 10⁻² | a 4.611 × 10⁻² (techo marginalizado con θ = 1) |
|---|---|---|---|
| 0.1 | −0.074 | 0.20 % | 0.34 % |
| 0.5 | −0.341 | 0.90 % | 1.57 % |
| 1.0 | −0.578 | 1.52 % | 2.66 % |
| 2.0 | −0.828 | 2.18 % | 3.82 % |

Frente a la precisión actual de las sirenas:
- GW170817: D_L = 40^{+8}_{−14} Mpc (+20/−35 %; Abbott et al. 2017, PRL 119, 161101,
  arXiv:1710.05832).
- Sirenas oscuras de GWTC-3: H0 = 68^{+8}_{−6} km/s/Mpc (+12/−9 %), con el método de catálogo de
  galaxias y población BBH fija; el resumen avisa que depende de la distribución de masas (Abbott et
  al. 2023, ApJ 949, 76, arXiv:2111.03604).

Como en las SN, una sirena mide la forma de D_L(z), y el modo θ → 0 es la dirección de Ω_m: quitar
M_B no rompe esa degeneración.

## Predicción
```
