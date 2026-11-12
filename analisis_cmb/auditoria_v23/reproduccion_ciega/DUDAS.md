# DUDAS — puntos no fijados por ESPECIFICACION.md y lectura adoptada

En cada caso tomé la lectura más literal de la especificación y anoto el efecto que tiene (o que medí).

## D1. Tipo de spline cúbico para las tablas (sección 2)
La especificación admite "natural o not-a-knot". Usé **not-a-knot** (el valor por defecto de `scipy.interpolate.CubicSpline`).
Comprobé en p₀ (ω_m = 0.140923) que ambas condiciones dan r_d, r_s y z* idénticos a precisión de máquina (diferencia relativa 0.0),
como anticipa la especificación.

## D2. Tolerancia de gradiente (`gtol`) del optimizador (sección 5)
La especificación fija `ftol = 1e-10`, `xtol = 1e-10` y `max_nfev = 200`, pero no `gtol` ni el método de cálculo del jacobiano.
Lectura literal: dejé `gtol` en su valor por defecto de scipy (10⁻⁸) y el jacobiano por diferencias finitas hacia delante
(`jac='2-point'`, por defecto). Medí la sensibilidad: repitiendo desde p₀ con `gtol = 1e-14` el punto no cambia en ningún dígito
(diferencia 0.0 en los tres parámetros); desde otro arranque (h = 0.70, ω_m = 0.15) el punto cambia en < 3·10⁻⁹. p₀ a seis cifras
y χ² a cuatro decimales no dependen de esta elección.

## D3. Fila θ* en A: denominador de ∂ ln D̃_M (sección 6.2)
La fórmula dice ∂y_θ/∂X_j = −y_θ · ∂ ln D̃_M(u_final)/∂X_j "despreciando la dependencia del término de corrección en z*".
Lectura literal: ∂ ln D̃_M(u_final) = ∂D̃_M(u_final)/D̃_M(u_final), es decir, también el **denominador** es D̃_M(u_final) y no
D̃_M(z*) = D̃_M(u_final) + D̃_H(u_final)(z* − z*_ref). La alternativa (denominador D̃_M(z*)) difiere en un factor
D̃_M(u_final)/D̃_M(z*) − 1 = 2.5·10⁻⁶. La comprobación por diferencias finitas lo muestra con nitidez: en todos los nodos, la fila
θ* difiere del modelo exacto justo en 2.5·10⁻⁶ relativo, muy por debajo del 1 % exigido. Efecto sobre la cota: despreciable.

## D4. Qué se entiende por "filas donde no son cero" en la comprobación de A (sección 6.2)
Muchas entradas de A son estructuralmente no nulas pero numéricamente minúsculas (p. ej. ∂y_θ/∂X_j para nodos con z ≳ 100 vale
10⁻¹⁰…10⁻¹⁴), y la diferencia finita con paso 10⁻⁴ tiene un ruido de redondeo de ≈ 10⁻¹² en esa fila. Criterio adoptado: comparar
las filas con |A_ij| > 10⁻⁶ · máx_i |A_ij| (umbral relativo por columna). Con ese criterio, de 20 nodos aleatorios de toda la malla,
los 5 con z > 200 fallan el 1 % por puro redondeo (diferencia absoluta ≤ 1.5·10⁻¹²); los 15 restantes y los 20 nodos aleatorios con
z ≤ 2.33, el nodo 0, los 7 nodos BAO y el nodo final (salvo por el mismo motivo de redondeo) cumplen el 1 % con margen (≈ 2.5·10⁻⁶).
Lo entrego así, sin excluir los nodos que fallan, para que se vea el motivo.

## D5. Elección de los "al menos 20 nodos al azar"
Elegí 20 nodos con `numpy.random.default_rng(20261112)` de toda la malla y, como muchos caen a z alto donde solo la fila θ* es no
nula, añadí otros 20 nodos aleatorios entre los de z ≤ 2.33 y los nodos especiales (0, BAO, final). La especificación no fija la
semilla ni la distribución.

## D6. Interpolación de la fila SN cuando u_HD coincide exactamente con un nodo
Uso `searchsorted(side='right') − 1`, de modo que el nodo k es el mayor con u_k ≤ u_HD y el peso t ∈ [0, 1). Si u_HD cayese
exactamente en un nodo, t = 0 y el resultado es el valor del nodo (no hay ambigüedad numérica). zHD de las SN está en
[0.01016, 2.26137], siempre dentro de la malla.

## D7. Eliminación de duplicados en la malla (sección 3)
Usé `numpy.unique` (duplicados exactos en coma flotante). Ninguno de los ln(1 + z_BAO) coincidió exactamente con un nodo del
linspace; resultaron 3499 + 7 = 3506 nodos, como pide `grid_n_nodes`.

## D8. Fórmulas de ajuste de la sección 13
Para r_d usé la ec. 26 de Eisenstein y Hu (1998), que es la "fórmula de ajuste del horizonte de sonido en la época de arrastre"
con ω_m y ω_b como únicos argumentos; como diagnóstico adicional (no pedido) calculé también la expresión cerrada completa de EH98
(ecs. 2–6 con z_d, z_eq, R, T_CMB = 2.7255 K). Para z* usé la fórmula de Hu y Sugiyama (1996) en la forma estándar con g₁, g₂. En
ambos casos tomé ω_m de la tabla tal cual (con neutrinos) y ω_b = 0.02237.

## D9. Lectura de la covarianza de Pantheon+
"1701 × 1701 números, uno por línea, por filas": leí los 1701² valores en orden y los reordené fila-mayor. La matriz resultante es
exactamente simétrica (asimetría máxima 0.0), lo que confirma el orden (fila o columna serían indistinguibles de todos modos).

## D10. Constante de ln Z
Entrego ln Z(σ_X) − ln Z(0). En la columna "ln Z(σ_X) (sin constante)" de la tabla doy −½[rᵀS⁻¹r − bᵀF⁻¹b + ln det S + ln det F]
tal cual, sin añadir ningún término (ni ln det C, ni priors), como pide la sección 12.
