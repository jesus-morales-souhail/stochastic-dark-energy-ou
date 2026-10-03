# σ_X = 3.4129e-04 — 200 realizaciones, semilla 20261003, lensing no lineal: takahashi

H0 de cada curva resuelto para 100θ_* = 1.041198 (ΛCDM de referencia), ω_b, ω_c, mnu fijos, plano.
Δχ² = χ²(curva) − χ²(ΛCDM) frente a datos Planck. Fuente: `resumen.json`, `dchi2_*.txt`, `fallos_*.json` de esta carpeta.

## Curvas de referencia

| curva | H0 resuelto | 100θ_* con H0 = 67.36 | Δχ² TT 2≤l<30 | Δχ² φφ L 8–400 | fallos |
|---|---|---|---|---|---|
| LCDM | 67.360 | 1.04120 | 0 | 0 | 0 |
| CPL_repo | 67.189 | 1.04171 | 0.0129 | 0.00877 | 0 |
| CPL_joint | 67.270 | 1.04147 | 0.00892 | 0.0176 | 0 |
| DESI_DR2_CMB | 63.711 | 1.05287 | 0.161 | 0.114 | 0 |

Ficheros Planck: TT `COM_PowerSpect_CMB-TT-full_R3.01.txt`; φφ `smicadx12_Dec5_ftl_mv2_ndclpp_p_teb_consext8_bandpowers.dat`.

## OU: TT y φφ (Δχ² frente a ΛCDM)

| obs. | θ | Δ | ρ≤0 | fallos (CAMB+NaN) | válidas | Δχ² p50 | Δχ² p95 | frac >4 | hueco banda68 [σ] | borde banda68 [σ] | p95 vs CPL repo | p95 vs CPL joint | nota |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| TT l<30 | 0.001 | 0.05 | 0 | 0+0 | 200 | -2.72e-06 | 3.39e-05 | 0.000 | 0 | 4.89e-06 | -0.0128 | -0.00889 | θ de DESI a esta σ_X (borde de malla) |
| TT l<30 | 1 | 0.05 | 0 | 0+0 | 200 | -3.53e-05 | 0.000879 | 0.000 | 0 | 0.000126 | -0.012 | -0.00804 |  |
| φφ | 0.001 | 0.05 | 0 | 0+0 | 200 | -1.55e-06 | 3e-05 | 0.000 | 0 | 1.55e-05 | -0.00874 | -0.0176 | θ de DESI a esta σ_X (borde de malla) |
| φφ | 1 | 0.05 | 0 | 0+0 | 200 | -4.27e-05 | 0.000667 | 0.000 | 0 | 0.000345 | -0.0081 | -0.0169 |  |

## OU: H0 resuelto (Planck 67.36 ± 0.54; DESI+BBN 68.51 ± 0.58, sólo comparación)

| θ | Δ | fallos H0 | p2.5 | p16 | p50 | p84 | p97.5 | frac fuera Planck ±2σ | frac fuera DESI+BBN ±2σ | 100θ_* con H0=67.36 [p16, p50, p84] |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.001 | 0.05 | 0 | 67.36 | 67.36 | 67.36 | 67.36 | 67.36 | 0.000 | 0.000 | [1.04120, 1.04120, 1.04120] |
| 1 | 0.05 | 0 | 67.35 | 67.35 | 67.36 | 67.37 | 67.37 | 0.000 | 0.050 | [1.04118, 1.04120, 1.04122] |

## Sensibilidad a Δ (0.05 frente a 0.10)

| θ | obs. | Δ banda68 [σ_Planck] | Δ p95 de Δχ² |
|---|---|---|---|

σ de H0 en la sensibilidad a Δ: error de Planck (0.54).
