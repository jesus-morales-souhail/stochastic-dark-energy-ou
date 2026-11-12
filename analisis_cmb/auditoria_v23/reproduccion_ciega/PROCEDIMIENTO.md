# PROCEDIMIENTO — pasos seguidos, versiones y resultados intermedios

Fecha: 2026-11-12. Carpeta de trabajo: esta (`reproduccion_ciega/`). No se salió de ella ni se consultó nada externo.

## 0. Entorno

| componente | versión |
|---|---|
| Python | 3.13.13 (conda-forge, GCC 14.3.0) |
| numpy | 2.5.0 |
| scipy | 1.18.0 |
| BLAS/LAPACK | scipy-openblas 0.3.33 |
| sistema | Linux 6.6.87.2-microsoft-standard-WSL2 (Ubuntu en WSL2), x86_64 |

Código propio (todo en `mi_codigo/`):

- `modelo.py`: constantes, tablas (splines), malla en u, lectura de BAO/Pantheon+, modelo exacto convenio B, χ², M analítico.
- `paso1_datos_fiducial.py`: comprobaciones de datos y malla, ajuste del fiducial p₀. Salida `salidas/paso1.json`, `paso1.log`, `p0.npy`.
- `paso2_lineal_lnZ.py`: B, A analítica, comprobación de A por diferencias finitas, A K Aᵀ, ln Z en la malla, regla 1.92. Salida `salidas/paso2.json`, `paso2.log`, `A.npy`, `B.npy`, `AKAt.npy`.
- `paso2b_fd_nodos_especiales.py`: diferencias finitas en el nodo 0, los 7 nodos BAO y el último. Salida `salidas/paso2b.json`.
- `paso3_tablas_vs_formulas.py`: tarea aparte (sección 13). Salida `salidas/paso3.json`.

Orden de ejecución: `cd mi_codigo && python3 paso1_datos_fiducial.py && python3 paso2_lineal_lnZ.py && python3 paso2b_fd_nodos_especiales.py && python3 paso3_tablas_vs_formulas.py`. Tiempo total ≈ 5 s.

## 1. Manifiesto
`sha256sum -c MANIFIESTO_SHA256.txt`: los 7 ficheros OK (ESPECIFICACION.md, 2 de DESI, 2 de Pantheon+, constantes.json, rd_rs_zstar_vs_omega_m.csv).

## 2. Lectura de la especificación
Leí ESPECIFICACION.md completa y `constantes.json`. No encontré discrepancias entre el texto y el fichero: σ_X_grid, extensión,
pasos de diferencias finitas, límites y arranque del fiducial, Δu = 0.002, `grid_u_max` = 6.994764490536918 y `grid_n_nodes` = 3506
coinciden. Las 13 filas BAO de `constantes.json` (z y tipo) coinciden con `mean.txt` (lo comprueba un `assert` en `modelo.py`).

## 3. Datos y malla (resultados intermedios)

| comprobación | resultado | esperado |
|---|---|---|
| filas en Pantheon+SH0ES.dat | 1701 | 1701 |
| SN seleccionadas (zHD > 0.01, IS_CALIBRATOR == 0) | **1580** | 1580 |
| zHD de las seleccionadas | [0.01016, 2.26137] | dentro de la malla |
| covarianza SN: cabecera / nº de valores | 1701 / 1701² | — |
| asimetría máx. de C_SN restringida | 0.0 | — |
| C₁₃ simetrizada: ½(C + Cᵀ) | hecho | — |
| u* = ln(1 + z*_ref) | 6.994764490536918 | = grid_u_max |
| n lineal = ceil(u*/Δu) + 1 | 3499 | — |
| nodos tras añadir los 7 ln(1+z_BAO) y `unique` | **3506** | 3506 |
| Δu_k mínimo / máximo | 3.09·10⁻⁵ / 1.9996·10⁻³ | no uniforme cerca de BAO |
| índices de los nodos BAO (0.295, 0.51, 0.706, 0.934, 1.321, 1.484, 2.33) | 130, 208, 270, 333, 426, 461, 608 | nodos exactos (|Δu| < 10⁻¹²) |
| longitud del vector de datos | 1595 | 13 + 2 + 1580 |
| Cholesky de C (1595×1595) | sin error (definida positiva) | — |

## 4. Fiducial p₀ (sección 5)

- Arranque: h = 0.682, ω_m = 0.1407, M₀ analítico = 23.797477150 (fórmula 1ᵀC_SN⁻¹r_SN / 1ᵀC_SN⁻¹1 con Cholesky de C_SN). χ² en el arranque = 1406.452732.
- `scipy.optimize.least_squares(method='trf', bounds, x_scale=[0.01, 0.001, 0.01], ftol=1e-10, xtol=1e-10, max_nfev=200)` sobre el residuo blanqueado L⁻¹(d − y).
- Resultado: status 2 (`ftol` satisfecho), nfev = 4. **No hubo parada por máximo de evaluaciones**, no hubo que repetir.
- p₀ = (h, ω_m, M) = (0.6811798064334976, 0.14092328770632676, 23.797874740711453) → **(0.681180, 0.140923, 23.797875)**.
- χ²(p₀) = 1406.3271474273733 → **1406.3271**. Por bloques: BAO 13.4976, θ* 0.2271, ω_m 3.6156, SN 1388.9869.
- Sensibilidad (ver DUDAS D2): con gtol = 10⁻¹⁴ desde p₀: mismo punto (dif. 0.0). Desde (0.70, 0.15, M₀): dif. < 3·10⁻⁹, χ² igual.
- Modelo en p₀: r_d = 147.65992 Mpc, r_s(z*) = 144.98557 Mpc, z*(ω_m0) = 1089.72656, D̃_M(u_final) = 3.16294704, D̃_M(z*) = 3.16293915, y_θ = 1.0412477 (dato 1.0411), y_BAO = [7.9630, 13.3527, 22.5804, 17.5282, 20.0723, 21.8060, 17.5191, 27.8828, 14.0564, 30.0762, 12.8805, 8.6322, 38.9925].
- Comprobación de consistencia posterior: en σ_X = 0 el desplazamiento lineal δp = F⁻¹b vale (−3.3·10⁻¹⁰, 1.0·10⁻¹⁰, 1.6·10⁻¹⁰): p₀ es efectivamente el mínimo.

## 5. Modelo lineal (sección 6)

- B (1595×3) por diferencias centradas con pasos 10⁻⁴, 10⁻⁵, 10⁻³. Columna M: exactamente 1 en SN (desviación máx. 0.0) y 0 en el resto. Fila ω_m: (0, 1.000000000001, 0).
- A (1595×3506) analítica, convenio B. Construcción: g_m = ½(1+z_m)E_m⁻³; ∂D̃_M(u_k)/∂X_j = −w^{(k)}_j g_j para 1 ≤ j ≤ k, y columna 0 = −Σ_{j≥1}; ∂D̃_H(u_k)/∂X_j = −½E_k⁻³(δ_kj − δ_0j). Filas BAO/θ*/SN según 6.2. Fila ω_m = 0.
- Suma de filas: máx |Σ_j A_ij| = 1.33·10⁻¹⁵ (**< 10⁻¹²: OK**).
- Diferencias finitas (±10⁻⁴, X(u₀) = 0; para j = 0 se perturban todos los demás nodos en ∓10⁻⁴), criterio de fila no nula |A_ij| > 10⁻⁶·máx|A_·j| (DUDAS D4):
  - 20 nodos al azar de toda la malla (semilla 20261112): **15/20 cumplen el 1 %** (máx. dif. rel. 2.5·10⁻⁶ hasta z ≈ 19, 4·10⁻⁴ en z ≈ 41, 0.9 % en z ≈ 93, 0.5 % en z ≈ 119–128); **5/20 no cumplen**: nodos con z = 201, 438, 478, 942, 1004, donde la única fila no nula es θ* con valor analítico 1.5·10⁻¹¹, 9·10⁻¹³, 7·10⁻¹³, 5·10⁻¹⁴, 4·10⁻¹⁴ y la diferencia finita es ruido de redondeo (dif. absoluta ≤ 1.5·10⁻¹²). Tabla completa en RESULTADO_REPRODUCCION.md §5.3.
  - 20 nodos al azar con z ≤ 2.33: **20/20 cumplen**, máx. dif. rel. 2.50·10⁻⁶ en todos (siempre fila θ*; en filas BAO y SN < 10⁻⁷).
  - Nodo 0: 1594 filas no nulas, máx. dif. rel. 2.5·10⁻⁶ (OK). Nodos BAO 130, 208, 270, 333, 426, 461, 608: OK (2.5·10⁻⁶). Nodo final 3505: analítico 1.5·10⁻¹⁴, FD = 0 (redondeo, mismo motivo).
  - El 2.5·10⁻⁶ sistemático en la fila θ* es D̃_M(u_final)/D̃_M(z*) − 1 (DUDAS D3).

## 6. ln Z(σ_X) (secciones 7 y 8)

- K_ij = exp(−|u_i − u_j|) (θ = 1), A K Aᵀ calculado una vez y simetrizado.
- Para cada σ_X: S = C + σ_X² AKAᵀ, Cholesky de S (`scipy.linalg.cholesky` + `cho_solve`), S⁻¹r y S⁻¹B, F = BᵀS⁻¹B, b = BᵀS⁻¹r, ln det S y ln det F de las diagonales de los factores. ln Z = −½[rᵀS⁻¹r − bᵀF⁻¹b + ln det S + ln det F].
- Cholesky de S tuvo éxito en los 10 nodos (S definida positiva). ln Z finito en todos. Valores: ver tabla en RESULTADO_REPRODUCCION.md §2.
- Máximo en σ_X = 0.1 (nodo 5 de 10), no en el último: **no se extendió la malla**.

## 7. Regla de la caída de 1.92 (sección 9)

v_k₀ = 1.139519, t = −0.780481. k = 7 (σ = 0.2, v = 0.400261 ≥ t > v = −1.157912 en σ = 0.3). Interpolación en ln σ:
σ_cota = exp[ln 0.2 + (0.400261 + 0.780481)/(0.400261 + 1.157912)·(ln 0.3 − ln 0.2)] = **0.2719** (0.27193636).

## 8. Tarea aparte (sección 13), realizada después de cerrar la cota

- r_d tabla vs EH98 ec. 26: dif. rel. máx. −3.17 % en ω_m = 0.26 (tabla 124.583, fórmula 128.663); −1.21 % en 0.08, −1.78 % en 0.14. Diagnóstico con la expresión cerrada completa de EH98: máx. −2.75 % en ω_m = 0.08.
- z* tabla vs HS96: dif. rel. −0.1833 % **constante** en los 361 nodos (cociente 0.9981670226, dispersión 9·10⁻¹³); máx. nominal en ω_m = 0.084. z*_tabla(ω_m,Planck) = 1089.906728619 = z*_ref; HS96 sin reescalar daría 1091.908.
- La cota no se modificó.

## 9. Incidencias
- Dos errores de sintaxis triviales en `paso3_tablas_vs_formulas.py` (nombres de clave con punto) corregidos antes de obtener resultados; no afectan a nada anterior.
- Ninguna comprobación numérica falló salvo las 5 (+1) columnas de diferencias finitas a z > 200 explicadas arriba, que se entregan tal cual.
