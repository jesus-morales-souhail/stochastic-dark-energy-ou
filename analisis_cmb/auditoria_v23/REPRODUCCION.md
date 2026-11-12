# REPRODUCCION: preparación de la reproducción a ciegas (2026-11-12)

Qué se preparó: una carpeta fuera de este repositorio, `~/Proyectos/01_Fisica_y_Cosmologia/reproduccion_ciega/`, para que una
persona sin acceso a nuestro código calcule la cota de la marginal lineal −1.92 sobre σ_X del OU con θ = 1 con BAO + θ* + ω_m +
Pantheon+. La carpeta no está en git ni se añadirá; este fichero solo deja constancia de qué contiene y del sha256 de la
especificación, sin el resultado esperado ni ninguna cifra nuestra.

Contenido de la carpeta:

- `ESPECIFICACION.md`: especificación autosuficiente (datos y columnas, covarianzas, modelo exacto en convenio B, distancias, θ* y
  su calibración, selección de supernovas y marginalización de M, priors, punto fiducial, malla en u, núcleo OU y cómo entra X(0),
  respuesta lineal con fórmulas y comprobación por diferencias finitas, verosimilitud marginal completa, malla de σ_X, máximo global,
  regla de la caída de 1.92, entregables y comprobaciones). Sin cifras de resultados.
- `datos/pantheon_plus/`: copias de `Pantheon+SH0ES.dat` y `Pantheon+SH0ES_STAT+SYS.cov` (release público de Pantheon+SH0ES).
- `datos/desi_dr2/`: copias de `desi_gaussian_bao_ALL_GCcomb_mean.txt` y `_cov.txt` (compresión gaussiana oficial de DESI DR2).
- `tablas/rd_rs_zstar_vs_omega_m.csv`: r_d(ω_m), r_s(z*; ω_m) y z*(ω_m) en 361 nodos, exportados de nuestras tablas; y
  `tablas/constantes.json` (c, Ω_r h², z* de referencia, los dos números de Planck con sus errores, la calibración de θ*, la malla en u, la malla de σ_X, los ajustes del fiducial y los pasos de diferencias finitas).
- `MANIFIESTO_SHA256.txt`: sha256 de todos los ficheros anteriores, comprobado contra los ficheros presentes al crearlo.

No contiene ningún script de este proyecto, ninguna ruta a `analisis_cmb/` ni a `auditoria_v23/`, ni ninguna cota, ln Z o máximo
nuestros (comprobado con grep en los ficheros de texto; en los de datos no se busca).

sha256 de `ESPECIFICACION.md` (versión 2, 2026-11-12: fiducial detallado, constantes y mallas solo de constantes.json, sección 'Qué comprueba esta reproducción', tarea aparte de comparación de tablas): `e9db3e095bb2d8fc2535968a4a938a7e4d0b66418625ca476831b61c2a2866db`

Regla fijada: la especificación se escribió antes de que el segundo investigador haga el cálculo y sin darle acceso a nuestros
scripts; cualquier duda sobre el método se responde por escrito y sin cifras. La reproducción no ha empezado.

## Reproducción ciega (12-11-2026)

**Cómo se hizo.** Una sesión nueva de Claude Code, sin contexto del proyecto y sin acceso a nuestros scripts, reprodujo la cota de
cabecera (OU, θ = 1, convenio B, BAO + θ* + ω_m + Pantheon+) solo con `ESPECIFICACION.md` versión 2 (sha256
`e9db3e095bb2d8fc2535968a4a938a7e4d0b66418625ca476831b61c2a2866db`), los datos y las tablas de la carpeta ciega (manifiesto 7/7 OK). Sus
ficheros (`RESULTADO_REPRODUCCION.md`, `PROCEDIMIENTO.md`, `DUDAS.md`, su código `mi_codigo/*.py` y sus json de salida) están copiados en
`auditoria_v23/reproduccion_ciega/` con `SHA256SUMS` comprobado contra los originales (12/12 OK). Sus json guardan la malla de ln Z, la
cota y el fiducial en doble precisión completa (no hubo que usar cifras redondeadas). No se ha modificado nada en la carpeta ciega.

**Predicciones (escritas por Claude antes de ejecutar, copiadas literales):**

- P1. La malla de ln Z(σ_X) − ln Z(0) del reproductor coincide con la nuestra en los 10 nodos con |diferencia| < 1e-4.
- P2. La cota del reproductor, SIN redondear, coincide con 0.27194769579562567 con diferencia relativa < 1e-4.
- P3. Nuestro r_d tabulado coincide con CAMB a mejor del 0.3 % en el punto de Planck y a mejor del 1 % en todo ω_m ∈ [0.12, 0.16].
- P4. Nuestra z* tabulada coincide con la z_star de CAMB a mejor del 0.5 % en todo ω_m ∈ [0.12, 0.16].

**Qué se compara (`comparacion_reproduccion.py` → `comparacion_reproduccion.json`).** La misma magnitud, ln Z(σ_X) − ln Z(0) de la marginal
lineal analítica con priors planos en (h, ω_m, M), en los mismos 10 nodos y el mismo orden; la nuestra es `positividad_lin.json` →
thetas.theta_1.lnZ_lin_rel (idéntica a `marginal_guiada.json` → thetas.theta_1.lnZ_lin_rel salvo 10⁻¹²). Las constantes son las mismas
(`constantes.json` se exportó de nuestro código).

| σ_X | nuestro | reproductor | diferencia | < 1e-4 |
|---|---|---|---|---|
| 0 | +0.000000 | +0.000000 | +0.00e+00 | sí |
| 0.01 | +0.129500 | +0.129491 | -9.10e-06 | sí |
| 0.03 | +0.603128 | +0.603123 | -5.15e-06 | sí |
| 0.05 | +0.872527 | +0.872562 | +3.54e-05 | sí |
| 0.07 | +1.030898 | +1.030975 | +7.73e-05 | sí |
| 0.1 | +1.139385 | +1.139519 | +1.34e-04 | no |
| 0.15 | +0.950065 | +0.950253 | +1.89e-04 | no |
| 0.2 | +0.400098 | +0.400261 | +1.63e-04 | no |
| 0.3 | -1.157825 | -1.157912 | -8.64e-05 | sí |
| 0.4 | -2.817695 | -2.818158 | -4.64e-04 | no |

Cota: nuestra 0.2719476958, reproductor 0.2719363570, diferencia relativa -4.17e-05. Fiducial: nuestro h = 0.6811798073,
ω_m = 0.1409232875, M = 23.797874740, χ² = 1406.3271474; reproductor h = 0.6811798064, ω_m = 0.1409232877, M = 23.797874741, χ² = 1406.3271474. Máximo de ln Z en σ_X = 0.1 en los dos.

**Qué valida:** la implementación estadística de la cota a partir de la especificación (selección y orden de datos, modelo exacto,
fiducial, respuesta lineal, marginal analítica, malla y regla de la caída de 1.92): la cota coincide a 4 × 10⁻⁵ y el fiducial a 10⁻⁹.

**Qué NO valida:** las tablas de r_d, r_s y z* y los datos (compartidos con el reproductor); el modelo físico; la validez de la
aproximación lineal (eso lo cubre INYECCION); θ = 0.3, 3 y 10 no se reprodujeron; y la especificación la escribimos nosotros.

**P1 no se cumple del todo: 6 de 10 nodos dentro de 10⁻⁴, diferencia máxima 4.6e-04 (en σ_X = 0.4).** Diagnóstico hecho
(sin cambiar nada): no es la doble inversión de la covarianza de Pantheon+ (nuestro `C = inv(inv(C))` difiere del fichero en
3.9e-16 y ln Z no cambia). Causa candidata, **no verificada numéricamente** (se para aquí, como se pidió): la fila SN de la
respuesta lineal. Nuestro `marginal_guiada.build_A` interpola entre los dos nodos la derivada ∂ln D̃_M de cada nodo
(MAG·[(1 − t) W_k + t W_{k+1}]); el reproductor (`paso2_lineal_lnZ.py`, línea 65) deriva el logaritmo de la D̃_M interpolada
([(1 − t) ∂D̃_M,k + t ∂D̃_M,k+1] / D̃_M(u_HD)), que es la derivada exacta del modelo lineal a trozos y por eso su comprobación por
diferencias finitas da 2.5 × 10⁻⁶. Las dos lecturas caben en la frase de la especificación ("la derivada se interpola con los mismos
pesos"); la diferencia es de orden Δu y explica el tamaño de 10⁻⁴ en ln Z. Efecto en la cota: 4 × 10⁻⁵ relativo (P2 cumplida).

**Fila θ* de A.** Nuestro código usa como denominador D̃_M(u_final) (`build_A`: `Wdm = −(T g)/dm` con `dm` la cumulativa hasta el último
nodo, y `A[13] = −Wdm[−1]·B0[13]`), igual que la lectura literal del reproductor (D3 de `DUDAS.md`); la alternativa D̃_M(z*) difiere en
2.5 × 10⁻⁶.

**Origen de las tablas** (`escalon/escalon.py`, commit 12a7aa8, 2026-11-11; `everpresent/rdcalc.py`, a65c4f8, 2026-10-16; valores de CAMB
en `everpresent/resultados/rd_planck.json`, generado por `rd_planck.py`, ba4a0b2 / 15ad950):
- r_d(ω_m): integral de c_s/H sin energía oscura desde z_d(ω_m), con z_d = z_d^CAMB · EH98(ω_m)/EH98(ω_m,Planck) (fórmula de la época de
  arrastre de Eisenstein y Hu como cociente) y la integral multiplicada por un factor único NORM = r_d^CAMB / integral(ω_m,Planck).
  **Es una integral física reescalada en un punto**: en el punto de Planck coincide con CAMB por construcción; lo que cuenta es el rango.
- r_s(z*; ω_m): la misma integral desde z*(ω_m), con el mismo NORM (no se renormaliza aparte: en el punto de Planck difiere de CAMB en
  9 × 10⁻⁶).
- z*(ω_m): fórmula de Hu y Sugiyama (1996) reescalada por un factor único para que z*(ω_m,Planck) = z*_CAMB; es lo que vio el reproductor
  (cociente constante 0.99816702 frente a HS96). Reescalada en un punto: por construcción coincide en Planck; cuenta el rango.

**Frente a CAMB** (`tablas_camb_check.py` → `tablas_camb_check.json`; CAMB 2.0.4, Python 3.13.13; parámetros fijos de `rd_planck.py`:
H0 = 67.36, ω_b = 0.02237, m_ν = 0.06 eV (ω_ν = 0.0006449), τ = 0.0544, n_s = 0.9649, A_s = 2.1 × 10⁻⁹, plano, w = −1 (PPF), T_CMB
2.7255 K, Y_He 0.2457 (BBN), N_eff por defecto, Recfast, precisión por defecto; solo cambia ω_c = ω_m − ω_b − ω_ν; tablas evaluadas con los
splines del código):

| ω_m | r_d tabla | r_d CAMB | dif. | z* tabla | z* CAMB | dif. | r_s(z*) tabla | r_s(z*) CAMB | dif. |
|---|---|---|---|---|---|---|---|---|---|
| 0.12 | 153.7020 | 153.6531 | +0.032 % | 1087.853 | 1087.819 | +0.0031 % | 150.8202 | 150.8408 | -0.014 % |
| 0.13 | 150.7027 | 150.6836 | +0.013 % | 1088.765 | 1088.747 | +0.0017 % | 147.9312 | 147.9413 | -0.007 % |
| 0.14 | 147.9084 | 147.9060 | +0.002 % | 1089.647 | 1089.643 | +0.0004 % | 145.2267 | 145.2279 | -0.001 % |
| 0.140923 | 147.6599 | 147.6583 | +0.001 % | 1089.727 | 1089.724 | +0.0002 % | 144.9856 | 144.9860 | -0.000 % |
| 0.143015 | 147.1027 | 147.1027 | -0.000 % | 1089.907 | 1089.907 | -0.0000 % | 144.4444 | 144.4431 | +0.001 % |
| 0.15 | 145.2967 | 145.2988 | -0.001 % | 1090.500 | 1090.508 | -0.0008 % | 142.6867 | 142.6801 | +0.005 % |
| 0.16 | 142.8483 | 142.8446 | +0.003 % | 1091.328 | 1091.346 | -0.0017 % | 140.2942 | 140.2808 | +0.010 % |

Máximos en ω_m ∈ [0.12, 0.16] (nodos cada 0.005 más los anteriores): r_d 0.032 %, r_s(z*) 0.014 %, z* 0.0031 %; en [0.13, 0.15]: 0.013 %,
0.007 %, 0.0017 %. En el punto de Planck r_d y z* coinciden a 10⁻¹² (por construcción). El −1.8 % a −3.2 % que vio el reproductor
frente a la ec. 26 de EH98 es la precisión de esa fórmula de ajuste, no de la tabla.

**Resultado de las predicciones:** P1 no (6 de 10 nodos; máximo 4.6 × 10⁻⁴); P2 sí (4.2 × 10⁻⁵); P3 sí (Planck 3 × 10⁻¹²; rango 0.032 %); P4 sí (0.003 %).
El paso de sensibilidad de las tablas no hizo falta (P3 y P4 cumplidas en todo el rango).

**Puntos a fijar en la especificación v3** (cada duda de `DUDAS.md` y la opción que tomamos nosotros):
- D1 spline: not-a-knot (`scipy.interpolate.CubicSpline` por defecto), que es lo que usa `escalon.py`.
- D2 optimizador: `gtol` por defecto de scipy (10⁻⁸) y jacobiano '2-point'; es lo que usa `exacto_sn.fit_combo`. Fijarlo.
- D3 fila θ*: denominador D̃_M(u_final), como nuestro `build_A`. Decirlo explícitamente.
- D4 comprobación de A: fijar el umbral relativo por columna (10⁻⁶ del máximo) y excluir de la comparación relativa los nodos con z > 200
  (solo la fila θ*, con valores ≤ 10⁻¹¹ donde la diferencia finita es ruido); exigir el 1 % en los nodos con z ≤ 2.33 y en los especiales.
- D5 nodos al azar: fijar semilla (p. ej. 20261112) y los dos conjuntos (20 de toda la malla, 20 con z ≤ 2.33) más los especiales.
- D6 interpolación SN: nodo k = el mayor con u_k ≤ u_HD (searchsorted por la derecha menos 1); nuestro código usa searchsorted por la
  izquierda menos 1, equivalente salvo empate exacto, que no ocurre.
- **Nueva, por P1:** fila SN de A: fijar una de las dos lecturas (derivada del logaritmo de la D̃_M interpolada, la del reproductor, o
  interpolación de la derivada por nodo, la nuestra); decisión del autor.
- D7 malla: `numpy.unique` sobre el concatenado, como nosotros.
- D8 fórmulas de la tarea aparte: decir que r_d se compare con la ec. 26 de EH98 y que la precisión esperable de esa fórmula es del
  2–3 %; y que nuestra z* es HS96 reescalada, así que esa comparación solo mide el factor.
- D9 covarianza de Pantheon+: fila mayor; nuestro `load_sn` hace `reshape(n, n)` del vector leído.
- D10 constante de ln Z: ninguna; como nosotros.

## P1: la fila SN de A (12-11-2026)

**Predicciones (de Claude, escritas antes de ejecutar, copiadas literales):**

- P5. Nuestra A coincide con el jacobiano numérico del modelo exacto en las filas SN con error relativo < 1e-5, y mejor que la A del reproductor.
- P6. Con la A que mejor coincida, la cota de θ = 1 cambia menos de 1e-4 relativo respecto a 0.27194769579562567.

**Criterio fijado antes de mirar:** A debe ser el jacobiano exacto, en X = 0 y p₀, del mismo modelo directo que usa la verosimilitud exacta;
gana la versión que coincida con ese jacobiano. **Criterio de convergencia:** diferencias centrales con h = 10⁻³ y h/10 = 10⁻⁴; convergido si
en los elementos no despreciables (|A_ij| > 10⁻⁶ · máximo de la columna) la diferencia relativa entre h y h/10 es < 10⁻⁶ (o la absoluta
< 10⁻¹² si el elemento es menor de 10⁻⁶); se usa h/10. Nada se ha cambiado: ni la cota de cabecera, ni `build_A`, ni el código del
reproductor; todo en scripts y json nuevos (`jacobiano_sn_check.py/json`, `sensibilidad_fila_sn.py/json`).

**Paso 1, cómo obtiene cada código la distancia de cada supernova:**
- Modelo directo exacto (el de la verosimilitud exacta): `exacto_sn.py` línea 147, `mu_sn`: `np.interp(sn["uhd"], g.U, dm_u)`, es decir D̃_M
  interpolada **linealmente en u = ln(1 + z)** entre los dos nodos que encierran u_HD, y μ = 5 log₁₀[(1 + z_hel) D̃_M] + M. `marginal_guiada.py`
  líneas 66–67 y 85 (`exact_model`) hacen exactamente lo mismo con pesos (1 − t, t). No integra hasta cada z_SN.
- Nuestra `build_A` (`marginal_guiada.py` líneas 130–131 y 147): calcula ∂ln D̃_M(u_k)/∂X_j en cada nodo (`Wdm`) y **interpola esa derivada**
  entre los dos nodos: A_SN = MAG·[(1 − t) Wdm_k + t Wdm_{k+1}].
- Reproductor (`paso2_lineal_lnZ.py` líneas 63–65; `modelo.py` 115–118 y 188): interpola ∂D̃_M (no su logaritmo) y divide por la D̃_M
  interpolada: A_SN = MAG·[(1 − t) ∂D̃_M,k + t ∂D̃_M,k+1] / D̃_M(u_HD), que es la derivada exacta del modelo lineal a trozos.

**Paso 2, misma magnitud:** las tres cosas están en ∂μ/∂X_j (magnitudes por unidad de X), sobre las mismas 1580 filas en el orden del
fichero, en el convenio B (X entra como X − X(0); el nodo 0 se perturba moviendo todos los demás en sentido contrario) y en el mismo p₀
(el nuestro; el del reproductor difiere en 10⁻⁹). No hizo falta ninguna conversión. La A del reproductor se tomó de su `A.npy` (solo lectura,
carpeta ciega) y además se reconstruyó su fórmula de la línea 65 con nuestros propios arrays: las dos coinciden a 7.9e-10; sus filas 0–14 coinciden
con las nuestras a 8.2e-09.

**Paso 3, jacobiano numérico** (`jacobiano_sn_check.json`): 46 nodos (los 20 con z ≤ 2.33 de su comprobación, los que rodean a las 10 SN de
menor z, 10 de z intermedia y 10 de mayor z, y el nodo 0). **Convergencia: 44 de 46 nodos cumplen el criterio; no converge el nodo
86 (z = 0.1876), con diferencia relativa máxima entre h y h/10 de 1.043e-06 (un 4 % por encima del umbral 10⁻⁶) en alguno de sus 789 elementos
no despreciables.** Por el criterio fijado se documenta y no se elige h por conveniencia: la comparación de abajo usa h/10 en todos los nodos,
y la conclusión se sostiene en los 44 nodos convergidos por sí solos (en el nodo 86 la discrepancia de nuestra A, 1.2e-02, es cuatro órdenes mayor
que la no convergencia).

| nodo | z | elementos no despreciables | convergido (máx. dif. rel. h frente a h/10) | nuestra A: error rel. máx. / mediana | A del reproductor: máx. / mediana |
|---|---|---|---|---|---|
| 0 | 0.0000 | 1580 | sí (3.5e-07) | 7.9e-04 / 2.1e-06 | 3.5e-09 / 3.2e-09 |
| 5 | 0.0100 | 1580 | sí (6.1e-07) | 2.3e-02 / 2.4e-05 | 1.1e-08 / 5.9e-09 |
| 6 | 0.0121 | 1580 | sí (6.1e-07) | 1.6e-01 / 2.4e-05 | 1.6e-08 / 5.9e-09 |
| 30 | 0.0618 | 1024 | sí (5.6e-07) | 3.3e-02 / 1.0e-05 | 3.7e-08 / 5.7e-09 |
| 52 | 0.1096 | 949 | sí (5.1e-07) | 1.6e-02 / 9.1e-06 | 1.2e-08 / 5.4e-09 |
| 82 | 0.1782 | 826 | sí (4.4e-07) | 8.2e-03 / 7.6e-06 | 1.0e-08 / 5.2e-09 |
| 84 | 0.1829 | 812 | sí (4.7e-07) | 1.1e-02 / 7.5e-06 | 2.9e-08 / 5.0e-09 |
| 85 | 0.1853 | 800 | sí (4.4e-07) | 8.6e-03 / 7.4e-06 | 1.3e-08 / 5.0e-09 |
| 86 | 0.1876 | 789 | NO (1.0e-06) | 1.2e-02 / 7.2e-06 | 1.9e-06 / 5.1e-09 |
| 87 | 0.1900 | 786 | sí (4.3e-07) | 9.7e-03 / 7.1e-06 | 1.0e-08 / 4.9e-09 |
| 146 | 0.3364 | 404 | sí (3.1e-07) | 5.0e-03 / 4.1e-06 | 2.3e-08 / 4.5e-09 |
| 167 | 0.3937 | 312 | sí (2.8e-07) | 2.6e-03 / 3.6e-06 | 1.2e-08 / 4.1e-09 |
| 198 | 0.4828 | 220 | sí (3.4e-07) | 5.3e-03 / 2.9e-06 | 1.5e-07 / 4.3e-09 |
| 200 | 0.4887 | 216 | sí (2.3e-07) | 2.2e-04 / 2.8e-06 | 1.3e-08 / 4.4e-09 |
| 210 | 0.5158 | 196 | sí (2.1e-07) | 3.2e-03 / 2.7e-06 | 1.4e-08 / 3.8e-09 |
| 211 | 0.5188 | 193 | sí (2.1e-07) | 2.6e-03 / 2.7e-06 | 1.7e-08 / 4.2e-09 |
| 232 | 0.5839 | 140 | NO (1.6e-06) | 4.5e-03 / 2.1e-06 | 1.8e-06 / 4.0e-09 |
| 234 | 0.5903 | 133 | sí (2.9e-07) | 4.3e-03 / 2.0e-06 | 1.1e-07 / 3.7e-09 |
| 257 | 0.6651 | 93 | sí (1.6e-07) | 3.7e-03 / 1.9e-06 | 5.7e-08 / 5.5e-09 |
| 298 | 0.8038 | 30 | sí (1.1e-07) | 2.1e-06 / 1.2e-06 | 2.1e-08 / 6.0e-09 |
| 304 | 0.8256 | 30 | sí (1.1e-07) | 2.1e-06 / 1.2e-06 | 1.2e-08 / 5.9e-09 |
| 310 | 0.8476 | 28 | sí (1.0e-07) | 2.1e-06 / 1.2e-06 | 2.1e-08 / 4.5e-09 |
| 375 | 1.0999 | 20 | sí (7.0e-08) | 1.5e-06 / 8.1e-07 | 2.6e-08 / 6.5e-09 |
| 411 | 1.2566 | 17 | sí (6.5e-08) | 1.4e-06 / 8.1e-07 | 4.0e-08 / 7.3e-09 |
| 433 | 1.3534 | 11 | sí (6.1e-08) | 9.8e-04 / 7.2e-07 | 6.2e-08 / 2.5e-08 |
| 434 | 1.3581 | 10 | sí (7.0e-08) | 9.7e-04 / 6.3e-07 | 3.3e-08 / 8.3e-09 |
| 440 | 1.3865 | 9 | sí (4.4e-08) | 6.6e-05 / 7.2e-07 | 6.0e-08 / 1.3e-08 |
| 441 | 1.3913 | 9 | sí (5.0e-08) | 1.4e-04 / 7.5e-07 | 4.3e-08 / 2.7e-08 |
| 446 | 1.4154 | 8 | sí (3.3e-08) | 3.0e-04 / 6.2e-07 | 5.0e-08 / 6.8e-09 |
| 447 | 1.4202 | 8 | sí (7.3e-08) | 1.8e-03 / 6.3e-07 | 9.0e-08 / 1.1e-08 |
| 472 | 1.5392 | 7 | sí (5.1e-08) | 2.3e-04 / 5.0e-07 | 2.6e-08 / 1.4e-08 |
| 473 | 1.5442 | 7 | sí (4.2e-08) | 5.5e-04 / 7.4e-07 | 4.3e-08 / 1.7e-08 |
| 474 | 1.5493 | 6 | sí (6.7e-08) | 1.3e-04 / 6.6e-07 | 4.0e-08 / 2.6e-08 |
| 486 | 1.6113 | 5 | sí (5.0e-08) | 2.3e-04 / 4.8e-07 | 3.6e-08 / 1.2e-08 |
| 487 | 1.6165 | 5 | sí (7.0e-08) | 5.4e-04 / 5.3e-07 | 4.7e-08 / 2.0e-08 |
| 502 | 1.6961 | 4 | sí (8.3e-08) | 2.3e-04 / 5.9e-07 | 7.0e-08 / 2.5e-08 |
| 503 | 1.7015 | 4 | sí (5.4e-08) | 1.6e-03 / 5.6e-07 | 7.6e-08 / 2.8e-08 |
| 507 | 1.7232 | 3 | sí (1.2e-08) | 7.4e-07 / 3.3e-07 | 1.0e-08 / 7.0e-09 |
| 516 | 1.7727 | 3 | sí (4.7e-08) | 7.7e-07 / 3.5e-07 | 2.7e-08 / 1.1e-08 |
| 521 | 1.8006 | 3 | sí (6.2e-08) | 1.6e-04 / 7.9e-07 | 4.4e-08 / 3.0e-08 |
| 522 | 1.8062 | 3 | sí (2.4e-07) | 1.6e-03 / 6.8e-07 | 2.2e-07 / 6.2e-08 |
| 524 | 1.8174 | 2 | sí (2.3e-08) | 7.1e-07 / 5.1e-07 | 3.6e-08 / 3.6e-08 |
| 540 | 1.9090 | 2 | sí (4.2e-08) | 2.9e-04 / 1.5e-04 | 6.0e-08 / 5.5e-08 |
| 541 | 1.9148 | 2 | sí (4.3e-08) | 9.4e-04 / 4.7e-04 | 5.6e-08 / 3.6e-08 |
| 597 | 2.2602 | 1 | sí (6.2e-08) | 1.9e-04 / 1.9e-04 | 5.3e-08 / 5.3e-08 |
| 598 | 2.2667 | 1 | sí (5.5e-07) | 1.2e-03 / 1.2e-03 | 5.5e-07 / 5.5e-07 |

Agregado: nuestra A, error relativo máximo 0.158 (nodo 6, z = 0.0121, las SN de menor z) y mediana de medianas 1.6e-06;
A del reproductor, máximo 1.9e-06 y mediana 6.7e-09 (su fórmula reconstruida con nuestros arrays: 1.9e-06 y 6.0e-09).
**El jacobiano exacto del modelo directo es la A del reproductor; la nuestra no lo es en las filas SN** (interpolar la derivada del logaritmo por
nodo no es la derivada del logaritmo del interpolado; el error es grande donde ln D̃_M cambia deprisa entre nodos, z ≲ 0.02).

**Paso 4:** no aplica (el modelo directo interpola, no integra hasta cada z_SN; no hay tercera A).

**Paso 5, sensibilidad** (`sensibilidad_fila_sn.json`; misma marginal, prior, malla y regla; Msig = A K Aᵀ sin simetrizar como en `positividad.py`;
corrido con 16 hilos BLAS, que es con lo que se generó `positividad_lin.json`: con otro número de hilos ln Z cambia en 10⁻⁸ y la cota en
≤ 6 × 10⁻⁹ por el orden de suma, lo que se comprobó antes de dar por reproducida la referencia):

| σ_X | (a) nuestra A | (b) nuestra A con la fila SN del reproductor | (c) A.npy del reproductor | (b) − (a) |
|---|---|---|---|---|
| 0 | +0.000000000 | +0.000000000 | +0.000000000 | +0.00e+00 |
| 0.01 | +0.129499973 | +0.129490871 | +0.129490871 | -9.10e-06 |
| 0.03 | +0.603128482 | +0.603123330 | +0.603123330 | -5.15e-06 |
| 0.05 | +0.872526675 | +0.872562045 | +0.872562044 | +3.54e-05 |
| 0.07 | +1.030897522 | +1.030974806 | +1.030974804 | +7.73e-05 |
| 0.1 | +1.139384714 | +1.139518957 | +1.139518955 | +1.34e-04 |
| 0.15 | +0.950064525 | +0.950253031 | +0.950253030 | +1.89e-04 |
| 0.2 | +0.400097894 | +0.400260997 | +0.400260999 | +1.63e-04 |
| 0.3 | -1.157825499 | -1.157911909 | -1.157911902 | -8.64e-05 |
| 0.4 | -2.817694762 | -2.818158499 | -2.818158489 | -4.64e-04 |

Cotas: (a) 0.27194769579562567 (= referencia, bit a bit); (b) 0.27193635653203496 (-4.17e-05 relativo);
(c) 0.27193635707638547 (-4.17e-05). Con (b) y (c) la malla de ln Z coincide con la del reproductor a 6.7e-09 y
2.8e-09: **la fila SN explica toda la diferencia de P1.**

**Conclusión por el criterio fijado:** gana la A del reproductor. P5: **no** (nuestro error máximo 0.158, no < 10⁻⁵, y peor que la del
reproductor, 1.9e-06). P6: **sí** (la cota cambia -4.17e-05 relativo, 0.27195 → 0.27194). Como la ganadora no es la nuestra, no se
corrige nada: queda anotada como errata candidata en `manuscript/ERRATAS_V3.md` (fila 129, categoría J) con el impacto medido, pendiente de
decisión del autor. Impacto no medido: las demás cotas lineales de la auditoría (θ = 0.3, 3, 10; familias; muestras de SN; positividad;
INYECCION) usan la misma `build_A`.

**Punto a fijar en la especificación v3 (definición exacta de la fila SN de A):** ∂y_SN/∂X_j = (5/ln 10) · [(1 − t) ∂D̃_M(u_k)/∂X_j +
t ∂D̃_M(u_{k+1})/∂X_j] / D̃_M(u_HD), con D̃_M(u_HD) = (1 − t) D̃_M(u_k) + t D̃_M(u_{k+1}), k el mayor nodo con u_k ≤ u_HD y t = (u_HD − u_k)/(u_{k+1} − u_k):
la derivada del logaritmo de la D̃_M interpolada (la del reproductor), que es el jacobiano exacto del modelo directo. La frase actual de la
sección 6.2 ("la derivada se interpola con los mismos pesos") se sustituye por esa fórmula.
