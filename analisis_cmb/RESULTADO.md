# RESULTADO — OU a varianza máxima DESI frente a Planck 2018

## 1. Qué se construyó

- σ_X = 2.506e-2: límite 95 % del perfil DESI DR2 BAO (`../results/profile_sigma_x/profile_sigma_x.json`). NOTAS: "Fuente de verdad".
- X: proceso gaussiano estacionario con la covarianza del kernel BAO, σ_X² exp(−θ|Δln(1+z)|) (`../scripts/desi_dr2_data.py::add_ou_kernel`). NOTAS: "Definición de X".
- Lectura (a): δ ln ρ_DE = X/Ω_Λ0; la lectura (b) multiplica la amplitud por 1.46. NOTAS: "Decisiones de implementación".
- ρ_DE(a) = ρ_Λ0[1 + (X_s(a) − X_s(1))/Ω_Λ0], w = −1 − (1/3) d ln ρ_DE/d ln a → CAMB `DarkEnergyPPF.set_w_a_table`. NOTAS: "Decisiones de implementación".
- Δ (suavizado gaussiano en ln a) es un parámetro del modelo: X puntual depende de Δ como σ_X√(θΔ). NOTAS: "Escala de resolución Δ".
- Relleno estacionario hasta ln a = 0.4; camino con paso h_tab/16, tabla de 2000 puntos. NOTAS: "Relleno más allá de a = 1", "Convergencia de malla".
- H0 de cada curva resuelto por bisección para 100θ_* = 1.041198 (ΛCDM); ω_b, ω_c, mnu fijos, plano. NOTAS: "Cambio de método".
- Lensing no lineal halofit takahashi en todas las curvas. NOTAS: "Cambio de método".
- 200 realizaciones, semilla 20261003, números aleatorios comunes en θ, Δ y σ_X. NOTAS: "Estadístico y cota".

## 2. σ_X = 2.506e-2, Δ = 0.05

Fuente: `resultados/takahashi/sigma_2.506e-02/resumen.json`. Δχ² = χ²(OU_r) − χ²(ΛCDM), datos Planck (`COM_PowerSpect_CMB-TT-full_R3.01.txt`; `smicadx12_Dec5_ftl_mv2_ndclpp_p_teb_consext8_bandpowers.dat`). Fallos (ρ≤0, CAMB, NaN, H0): 0.

| θ | Δχ² TT p50 / p95 | Δχ² φφ p50 / p95 | frac Δχ²>4 | hueco banda68–ΛCDM | H0 p16 / p50 / p84 | H0 p2.5 / p97.5 | fuera Planck ±2σ |
|---|---|---|---|---|---|---|---|
| 0.001 | −1.4e-4 / 2.5e-3 | −1.2e-4 / 2.2e-3 | 0 | 0 | 67.34 / 67.36 / 67.38 | 67.32 / 67.40 | 0 |
| 1 | −2.4e-3 / 0.064 | −2.7e-3 / 0.051 | 0 | 0 | 66.94 / 67.36 / 67.85 | 66.55 / 68.28 | 0.025 |

Planck H0 = 67.36 ± 0.54 (arXiv:1807.06209, tabla 2, TT,TE,EE+lowE+lensing).

## 3. Malla de σ_X (10 puntos, 2.506e-2 → 1e-4)

Fuente: `resultados/takahashi/malla/malla.json` (`malla.md`). Fallos: 0 en los 20 casos.

- p95(Δχ²) < 4 en toda la malla. Máximo p95: TT 0.0636, φφ 0.0514 (θ = 1, σ_X = 2.506e-2). Planck no acota σ_X por debajo del límite de DESI.
- p95(Δχ²) ∝ σ_X^n: θ = 1: n = 0.998 (TT), 1.007 (φφ); θ = 0.001: n = 1.002 (TT), 0.998 (φφ). rms del residuo log10 ≤ 0.0037. El valor esperado era n ≈ 2; sale n ≈ 1.
- El p95 medido está dominado por el término cruzado con los residuos de Planck, 2Σ r_b δ_b/σ_b² (r = ΛCDM − dato); la contribución propia del modelo (cuadrática, Σ δ_b²/σ_b²) es p95 = 6.3e-4 (TT) y 7.4e-3 (φφ) a θ = 1, σ_X = 2.506e-2, frente a 0.063 y 0.050 del término cruzado (`resultados/takahashi/sigma_2.506e-02/descomposicion_dchi2.json`).
- σ(H0) ≡ (p84 − p16)/2 ∝ σ_X^m: θ = 1: m = 1.001; θ = 0.001: m = 0.999. rms del residuo log10 ≤ 0.0028.

## 4. Conclusiones

(a) Planck (TT 2≤l<30 y C_L^φφ L 8–400) no distingue el OU de ΛCDM a ninguna σ_X ≤ 2.506e-2: p95(Δχ²) ≤ 0.064 (`malla.json`).

(b) La única huella del OU en el CMB es un desplazamiento de H0, con dispersión σ(H0) = 0.454 × (σ_X/2.5e-2) km/s/Mpc para θ = 1 y 0.019 × (σ_X/2.5e-2) km/s/Mpc para θ = 0.001 (prefactores del ajuste σ(H0) ∝ σ_X^m evaluado en σ_X = 2.5e-2, `malla.json`; valores medidos en σ_X = 2.506e-2: 0.4558 y 0.0190, `malla.md`).

(c) Esa dispersión (prefactor del ajuste) es 12.5 veces (θ = 1) y 300 veces (θ = 0.001) menor que la diferencia Planck–SH0ES, 73.04 − 67.36 = 5.68 km/s/Mpc (SH0ES: arXiv:2112.04510; Planck: arXiv:1807.06209, tabla 2), así que el modelo no la explica.

## 5. Inconsistencias a corregir en los papers

Tres definiciones distintas de X en el repo:

1. `papers/stochastic-dark-energy-desi-dr2.md` l.70 — X ≡ δΩ_Λ(x), x = ln a.
2. `papers/stochastic-dark-energy-desi-dr2.md` l.562 — H² = H0²[Ω_m(1+z)³ + Ω_Λ + X(z)], i.e. X = δρ_DE/ρ_crit,0 ⇒ δ ln ρ_DE = X/Ω_Λ0.
3. `papers/unimodular-gravity-vacuum-smoothness.md` l.68 — X ≡ δρ_Λ/ρ̄_Λ = δ ln ρ_DE.

Además, el kernel BAO del código (δα_i = S(z_i)·X(z_i), S = ∂ln D_V/∂Ω_Λ con Ω_m = 1 − Ω_Λ compensando punto a punto) no define una cosmología única: cada bin ve un Ω_Λ (y un Ω_m) distinto. Lo que se construye en `analisis_cmb/` es la primera versión física (una sola historia ρ_DE(a), H(a)) del modelo.

Ajustes CPL del repo (no reajustados):
- `results/eos_cpl_desi_dr2/eos_cpl_summary.json`: w0 = −0.98990, wa = −0.01590, hecho con errores diagonales ("7 bins, diagonal errors"), no con la covarianza completa del perfil OU.
- `results/joint_w0wa_sigma/joint_w0wa_sigma.json`: w0 = −1.00165, wa = +0.01725. Se contradice: `"data": "DESI DR2 BAO 7 bins, full meas. cov, S(z) fixed"` pero `"cov_mode": "diagonal_figure6_errors"`.

## 6. Limitaciones

- Parámetros cosmológicos fijos a Planck 2018 salvo H0; lo riguroso sería un MCMC completo, no corrido porque Δχ² ≤ 0.064 en el peor caso (`malla.json`).
- Δ es un parámetro introducido en este análisis, no de la teoría original (NOTAS: "Escala de resolución Δ").
- Lectura (a) elegida entre dos posibles; la (b) multiplica la amplitud por 1.46 (NOTAS: "Decisiones de implementación").
- No verificado cómo trata halofit una w(a) tabulada en PPF (NOTAS: "Cambio de método").
