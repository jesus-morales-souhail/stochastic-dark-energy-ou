# oscilador — X(ln a) = A sin(2π ln a / P + φ) frente a DESI DR2 BAO y Planck 2018

Mismo pipeline que `../` (lectura (a), cero en a = 1, H0 resuelto por θ_*, halofit takahashi, mismos
ficheros de Planck), sin realizaciones: 60 curvas, A ∈ {0.005, 0.01, 0.025}, P ∈ {0.1, 0.2, 0.5, 1, 2},
φ ∈ {0, π/2, π, 3π/2}. X es suave: w(a) con derivada analítica, sin suavizado Δ.

BAO: datos y covarianza 13×13 con `scripts/desi_dr2_data.py::load_gaussian_bao_vector` (sin modificar).
Las funciones de distancia del repo (`joint_w0wa_sigma_desi.DM/DH/DV`, `eos_efectiva.DM_cpl/...`) sólo
aceptan CPL con H0 = 67.4, Ω_m = 0.315 fijos y sin r_d, así que las distancias D_M, D_H, D_V y r_d salen
de CAMB con la misma cosmología que la parte CMB (H0 resuelto). Validación: DESI DR2+CMB (w0 = −0.42,
wa = −1.75) da Δχ²_BAO = −22.6 frente a ΛCDM.

```bash
../venv/bin/python oscilador.py      # → resultados/{oscilador.json, oscilador.csv, tabla.md}
```
