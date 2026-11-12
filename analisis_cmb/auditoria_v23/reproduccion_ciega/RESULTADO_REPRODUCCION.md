# RESULTADO_REPRODUCCION — reproducción a ciegas de la cota −1.92 sobre σ_X (OU, θ = 1, convenio B)

Reproducción independiente hecha el 2026-11-12 siguiendo únicamente `ESPECIFICACION.md` (versión 2) y los ficheros de `datos/` y
`tablas/` de esta carpeta (manifiesto SHA-256 verificado: 7/7 OK). Código propio en `mi_codigo/`; pasos y comprobaciones
intermedias en `PROCEDIMIENTO.md`; decisiones no fijadas por la especificación en `DUDAS.md`.

## 1. Punto fiducial p₀ y χ²(p₀)  (sección 5)

| parámetro | valor (seis cifras) | valor completo |
|---|---|---|
| h | 0.681180 | 0.6811798064334976 |
| ω_m | 0.140923 | 0.14092328770632676 |
| M | 23.797875 | 23.797874740711453 |

χ²(p₀) = **1406.3271** (1595 filas, X = 0).

Descomposición por bloques (solo informativa): BAO 13.4976, θ* 0.2271, ω_m 3.6156, SN 1388.9869.

El optimizador (trust-region reflective, scipy `least_squares`) terminó por `ftol` en el primer intento (status 2,
4 evaluaciones de la función); no hubo que repetir por máximo de evaluaciones.

## 2. Tabla completa de ln Z(σ_X) − ln Z(0) en la malla  (secciones 7 y 8)

Malla usada: la de `constantes.json` (10 nodos). El máximo NO cayó en el último nodo, así que **no hubo que extender la malla**.

| σ_X | ln Z(σ_X) − ln Z(0) | ln Z(σ_X) (sin constante) | rᵀS⁻¹r − bᵀF⁻¹b | ln det S | ln det F |
|---|---|---|---|---|---|
| 0.0   | +0.000000 | 2299.169651 | 1406.3271 | -6044.3173 | 39.650898 |
| 0.01  | +0.129491 | 2299.299142 | 1405.8851 | -6042.4026 | 37.919238 |
| 0.03  | +0.603123 | 2299.772775 | 1403.8072 | -6039.0802 | 35.727401 |
| 0.05  | +0.872562 | 2300.042213 | 1401.7659 | -6036.5177 | 34.667370 |
| 0.07  | +1.030975 | 2300.200626 | 1399.8856 | -6034.2821 | 33.995235 |
| 0.1   | +1.139519 | 2300.309170 | 1397.3673 | -6031.2833 | 33.297705 |
| 0.15  | +0.950253 | 2300.119904 | 1394.1348 | -6026.8865 | 32.511903 |
| 0.2   | +0.400261 | 2299.569912 | 1391.9057 | -6023.0014 | 31.955863 |
| 0.3   | -1.157912 | 2298.011739 | 1389.0687 | -6016.2617 | 31.169528 |
| 0.4   | -2.818158 | 2296.351493 | 1387.1645 | -6010.4715 | 30.604012 |

ln Z es finito en todos los nodos; ln Z(0) − ln Z(0) = 0 exactamente.

## 3. Máximo global

ln Z máximo en **σ_X = 0.1** (índice k₀ = 5), con ln Z − ln Z(0) = **+1.139519**.

## 4. Cota por la caída de 1.92  (sección 9)

Umbral t = v_k₀ − 1.92 = -0.780481. Primer k ≥ k₀ con v_k ≥ t > v_{k+1}: k = 7 (σ_k = 0.2, σ_{k+1} = 0.3; v_k = +0.400261,
v_{k+1} = -1.157912). Interpolando linealmente en ln σ_X:

**σ_X (cota, caída 1.92) = 0.2719**  (valor completo 0.2719363570360128).

## 5. Comprobaciones de la sección 10

1. SN seleccionadas: **1580** de 1701 (zHD > 0.01 e IS_CALIBRATOR == 0). Nodos de la malla: **3506** (esperados 3506); nodos de la parte lineal 3499, u* = 6.994764490536918, idéntico al `grid_u_max` del fichero.
2. χ²(p₀) y p₀: apartado 1.
3. A por diferencias finitas: ver tablas abajo. Suma de filas de A: máx |Σ_j A_ij| = 1.33e-15 (< 10⁻¹²: sí).
4. ln Z finito en todos los nodos: sí. ln Z(0) − ln Z(0) = 0: sí.
5. Versiones: Python 3.13.13, numpy 2.5.0, scipy 1.18.0 (BLAS/LAPACK: scipy-openblas 0.3.33), Linux-6.6.87.2-microsoft-standard-WSL2-x86_64-with-glibc2.39.

### 5.3 Comprobación de A por diferencias finitas (±10⁻⁴ en X_j con X(u₀) = 0)

Criterio: se comparan las filas donde |A_ij| > 10⁻⁶ · máx_i |A_ij| (ver `DUDAS.md`, D4). 20 nodos al azar de **toda** la malla
(semilla 20261112):

| nodo j | z_j | filas no nulas | máx. dif. rel. | A analítica (fila peor) | A dif. finitas (fila peor) | máx. dif. abs. en la columna | < 1 % |
|---|---|---|---|---|---|---|---|
| 318 | 0.8774 | 32 | 2.5e-06 | 0.0001387 | 0.0001387 | 3.5e-10 | sí |
| 348 | 0.9895 | 29 | 2.5e-06 | 0.0001206 | 0.0001206 | 3e-10 | sí |
| 350 | 0.9975 | 29 | 2.5e-06 | 0.0001194 | 0.0001194 | 3e-10 | sí |
| 377 | 1.108 | 24 | 2.49e-06 | 0.000104 | 0.000104 | 2.6e-10 | sí |
| 426 | 1.321 | 18 | 2.51e-06 | 4.006e-05 | 4.006e-05 | 4.9e-10 | sí |
| 448 | 1.425 | 10 | 2.5e-06 | 7.074e-05 | 7.074e-05 | 1.8e-10 | sí |
| 1093 | 7.772 | 1 | 1.95e-06 | 9.747e-07 | 9.747e-07 | 1.9e-12 | sí |
| 1122 | 8.296 | 1 | 2.35e-06 | 7.961e-07 | 7.961e-07 | 1.9e-12 | sí |
| 1152 | 8.871 | 1 | 4.39e-06 | 6.456e-07 | 6.456e-07 | 2.8e-12 | sí |
| 1200 | 9.865 | 1 | 1.73e-06 | 4.616e-07 | 4.616e-07 | 8e-13 | sí |
| 1497 | 18.68 | 1 | 9.64e-06 | 5.765e-08 | 5.765e-08 | 5.6e-13 | sí |
| 1879 | 41.24 | 1 | 0.000405 | 3.941e-09 | 3.939e-09 | 1.6e-12 | sí |
| 2277 | 92.62 | 1 | 0.00866 | 2.377e-10 | 2.398e-10 | 2.1e-12 | sí |
| 2400 | 118.7 | 1 | 0.00529 | 9.939e-11 | 9.992e-11 | 5.3e-13 | sí |
| 2436 | 127.7 | 1 | 0.00463 | 7.696e-11 | 7.661e-11 | 3.6e-13 | sí |
| 2662 | 201.2 | 1 | 0.0591 | 1.534e-11 | 1.443e-11 | 9.1e-13 | NO |
| 3050 | 438.2 | 1 | 0.205 | 9.217e-13 | 1.11e-12 | 1.9e-13 | NO |
| 3093 | 477.6 | 1 | 2.31 | 6.717e-13 | 2.22e-12 | 1.5e-12 | NO |
| 3432 | 941.7 | 1 | 1 | 5.279e-14 | 0 | 5.3e-14 | NO |
| 3464 | 1004 | 1 | 1 | 4.13e-14 | 0 | 4.1e-14 | NO |

Resultado: 15 de 20 cumplen el 1 %. Los 5 que no cumplen son nodos con z > 200, donde la única fila no nula es la de θ* y el
valor analítico es de 10⁻¹¹ a 10⁻¹⁴; la diferencia finita con paso 10⁻⁴ tiene ahí un ruido de redondeo ≈ 10⁻¹² (y_θ ≈ 1.04, ε_máquina
· y_θ / 2·10⁻⁴ ≈ 5·10⁻¹³), de modo que la comparación relativa no es significativa. La diferencia absoluta en esas columnas es
≤ 1.5·10⁻¹² en todos los casos. **No es una discrepancia de la fórmula**: en todos los nodos con señal medible (z ≲ 120) la
coincidencia es mejor del 1 %.

20 nodos al azar entre los que tienen z ≤ 2.33 (donde las filas BAO y SN son no nulas):

| nodo j | z_j | filas no nulas | máx. dif. rel. | A analítica (fila peor) | A dif. finitas (fila peor) | máx. dif. abs. en la columna | < 1 % |
|---|---|---|---|---|---|---|---|
| 30 | 0.06183 | 1032 | 2.5e-06 | 0.0003203 | 0.0003203 | 8e-10 | sí |
| 52 | 0.1096 | 957 | 2.5e-06 | 0.0003118 | 0.0003118 | 7.8e-10 | sí |
| 82 | 0.1782 | 834 | 2.49e-06 | 0.0002976 | 0.0002976 | 7.4e-10 | sí |
| 146 | 0.3364 | 411 | 2.5e-06 | 0.0002596 | 0.0002596 | 6.5e-10 | sí |
| 167 | 0.3937 | 319 | 2.49e-06 | 0.0002451 | 0.0002451 | 6.1e-10 | sí |
| 198 | 0.4828 | 227 | 2.49e-06 | 0.0002228 | 0.0002228 | 5.6e-10 | sí |
| 200 | 0.4887 | 223 | 2.49e-06 | 0.0002213 | 0.0002213 | 5.5e-10 | sí |
| 210 | 0.5158 | 202 | 2.5e-06 | 0.0002147 | 0.0002147 | 5.4e-10 | sí |
| 211 | 0.5188 | 199 | 2.49e-06 | 0.000214 | 0.000214 | 5.3e-10 | sí |
| 232 | 0.5839 | 146 | 2.49e-06 | 0.0001985 | 0.0001985 | 4.9e-10 | sí |
| 234 | 0.5903 | 139 | 2.49e-06 | 0.000197 | 0.000197 | 4.9e-10 | sí |
| 257 | 0.6651 | 99 | 2.5e-06 | 0.0001802 | 0.0001802 | 4.5e-10 | sí |
| 298 | 0.8038 | 35 | 2.5e-06 | 0.0001521 | 0.0001521 | 3.8e-10 | sí |
| 304 | 0.8256 | 35 | 2.5e-06 | 0.000148 | 0.000148 | 3.7e-10 | sí |
| 310 | 0.8476 | 33 | 2.49e-06 | 0.000144 | 0.000144 | 3.6e-10 | sí |
| 375 | 1.1 | 24 | 2.5e-06 | 0.0001051 | 0.0001051 | 2.6e-10 | sí |
| 411 | 1.257 | 21 | 2.5e-06 | 8.665e-05 | 8.665e-05 | 2.2e-10 | sí |
| 507 | 1.723 | 5 | 2.5e-06 | 5.016e-05 | 5.016e-05 | 1.3e-10 | sí |
| 516 | 1.773 | 5 | 2.49e-06 | 4.748e-05 | 4.748e-05 | 1.2e-10 | sí |
| 524 | 1.817 | 4 | 2.49e-06 | 4.52e-05 | 4.52e-05 | 1.1e-10 | sí |

Resultado: 20 de 20 cumplen el 1 %. La diferencia relativa máxima es sistemáticamente ≈ 2.5·10⁻⁶ y siempre está en la fila θ*
(fila 13); es exactamente el cociente D̃_M(u_final)/D̃_M(z*) − 1 = 2.49e-06, es decir, la aproximación que la propia
especificación prescribe en la fila θ* (despreciar la dependencia del término de corrección en z*). En las filas BAO y SN la
coincidencia es mejor que 10⁻⁷.

Nodos especiales (nodo 0, los 7 nodos de BAO y el último nodo):

| nodo j | z_j | filas no nulas | máx. dif. rel. | A analítica (fila peor) | A dif. finitas (fila peor) | máx. dif. abs. en la columna | < 1 % |
|---|---|---|---|---|---|---|---|
| 0 | 0 | 1594 | 2.5e-06 | -0.1013 | -0.1013 | 2.5e-07 | sí |
| 130 | 0.295 | 514 | 2.5e-06 | 0.000135 | 0.000135 | 2.2e-09 | sí |
| 208 | 0.51 | 207 | 2.48e-06 | 0.0001081 | 0.0001081 | 1.3e-08 | sí |
| 270 | 0.706 | 78 | 2.51e-06 | 8.575e-05 | 8.575e-05 | 5.9e-09 | sí |
| 333 | 0.934 | 33 | 2.51e-06 | 6.463e-05 | 6.463e-05 | 2.3e-09 | sí |
| 426 | 1.321 | 18 | 2.51e-06 | 4.006e-05 | 4.006e-05 | 4.9e-10 | sí |
| 461 | 1.484 | 11 | 2.51e-06 | 3.299e-05 | 3.299e-05 | 2.6e-10 | sí |
| 608 | 2.33 | 3 | 2.51e-06 | 1.331e-05 | 1.331e-05 | 3.3e-11 | sí |
| 3505 | 1090 | 1 | 1 | 1.505e-14 | 0 | 1.5e-14 | NO |

Nodo 0 (columna "menos la suma de las demás"): 1594 filas no nulas, coincidencia 2.5·10⁻⁶. El último nodo (z = z*_ref) solo
afecta a θ* con valor analítico 1.5·10⁻¹⁴, indetectable por diferencias finitas (mismo motivo que arriba).

## 6. Tarea aparte (sección 13): tablas frente a fórmulas de ajuste

Fórmulas escritas por mí (en `mi_codigo/paso3_tablas_vs_formulas.py`), con ω_b = 0.02237 y la ω_m de la tabla (incluye neutrinos):

- r_d: Eisenstein y Hu (1998), ec. 26: r_d = 44.5 ln(9.83/ω_m) / √(1 + 10 ω_b^{3/4}) Mpc.
- z*: Hu y Sugiyama (1996): z* = 1048 [1 + 0.00124 ω_b^{−0.738}] [1 + g₁ ω_m^{g₂}], g₁ = 0.0783 ω_b^{−0.238}/(1 + 39.5 ω_b^{0.763}), g₂ = 0.560/(1 + 21.1 ω_b^{1.81}).

Diferencia relativa (tabla − fórmula)/fórmula sobre ω_m ∈ [0.08, 0.26] (361 nodos):

| tabla | fórmula | máx. dif. rel. (en valor absoluto) | en ω_m | tabla allí | fórmula allí | dif. rel. en ω_m = 0.08 / 0.14 / 0.26 |
|---|---|---|---|---|---|---|
| r_d | EH98 ec. 26 | **3.17 %** (tabla por debajo) | 0.26 | 124.583 Mpc | 128.663 Mpc | -1.21 % / -1.78 % / -3.17 % |
| z* | HS96 | **0.183 %** (tabla por debajo) | 0.084 | 1084.240 | 1086.231 | -0.1833 % / -0.1833 % / -0.1833 % |

Observaciones:

- **r_d**: la diferencia crece monótonamente con ω_m, de −1.2 % en ω_m = 0.08 a −3.2 % en ω_m = 0.26 (media −2.2 %); en el
  fiducial (ω_m ≈ 0.141) es −1.8 %. Esto es del orden de la precisión conocida de la ec. 26 de EH98 frente a códigos de Boltzmann
  (que la sobreestima en ~2 % para cosmologías tipo Planck), y la tendencia con ω_m es la esperable de un ajuste sin neutrinos
  masivos. No lo considero un error de la tabla, pero lo comunico porque supera ligeramente el "orden del 1 %" del enunciado.
  Diagnóstico adicional con la expresión cerrada completa de EH98 (ecs. 2–6, con z_d, z_eq y R): diferencia máxima
  2.75 % (en ω_m = 0.08), −2.5 % en ω_m = 0.14. La tabla en ω_m = ω_m,Planck da r_d = 147.10 Mpc.
- **z***: la diferencia relativa es **constante en toda la malla a 10⁻¹²**: z*_tabla / z*_HS96 = 0.99816702 en los 361 nodos.
  Es decir, la tabla z*(ω_m) es exactamente la fórmula de Hu–Sugiyama reescalada por un factor fijo, elegido de modo que
  z*(ω_m,Planck) = z*_ref = 1089.9067 (lo comprobé: el spline de la tabla en ω_m,Planck devuelve 1089.906728619 frente a
  1091.908 de HS96 sin reescalar). Es una observación sobre cómo se construyó la tabla, no una discrepancia de 0.18 %.

**La cota del apartado 4 no se ha tocado tras esta comparación.**
