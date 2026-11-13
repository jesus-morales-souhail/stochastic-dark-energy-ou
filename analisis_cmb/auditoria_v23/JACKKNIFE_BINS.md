# JACKKNIFE_BINS: ¿la preferencia por σ_X > 0 vive en un solo bin de DESI DR2? (13-11-2026)

## Contexto

Kim, Mota y Tamosiunas, arXiv:2607.28918 (v3, 22-10-2026), *"An E-value Reanalysis of DESI's Dynamical Dark Energy Signal"*, encuentran que la
evidencia de DESI DR2 BAO contra ΛCDM (w0wa, fondo Planck fijo) se concentra en el bin LRG2: quitándolo, su e-valor cae de 33.97 a 0.49. Su
modelo y su análisis no son los nuestros (nosotros: residuo OU, fiducial reajustado, marginal sobre (h, ω_m[, M])). La pregunta aquí es si
nuestra preferencia por σ_X > 0 sin supernovas (el máximo de ln Z(σ_X) − ln Z(0); es ese máximo, no un factor de Bayes) también vive en un
solo bin.

## Predicciones (de Claude, escritas antes de ejecutar, copiadas literales)

- P12. Sin SN (BAO + θ* + ω_m), OU θ = 1: quitar el bin LRG2 deja el máximo de ln Z(σ_X) − ln Z(0) en menos de la mitad del valor con todos los bins.
- P13. Sin SN, θ = 1: quitar cualquier otro bin individual cambia ese máximo menos de ±0.5.
- P14. Con Pantheon+ (BAO + θ* + ω_m + Pantheon+), θ = 1: quitar cualquier bin individual cambia ese máximo menos de ±0.5.

Decisión posterior del autor (tras la parada del paso 1, ver método): P12–P14 se mantienen literales; P12 se evalúa contra +1.306, el máximo
lineal con build_A_v3 y todos los bins, calculado con el mismo procedimiento que cada caso sin bin.

## Método (`jackknife_bins.py` → `jackknife_bins.json`; paso 1 en `jackknife_bins_paso1.json`)

**Qué se comprueba y qué no.** El jackknife comprueba la **marginal lineal** (analítica, convenio B, priors planos en (h, ω_m) y en M con SN),
no la marginal exacta por Monte Carlo. El "+1.4 sin SN" del encargo venía de PRECISION_SN.md, que es la marginal **exacta por Monte Carlo**
(2000 caminos, ESS ≈ 20 en la celda del máximo): +1.36 en θ = 1 (σ_X = 0.15) y +1.38 en θ = 0.3. La lineal antigua de `maximos_globales.json`
(malla de cmb_fondo, fiducial fijo OM_FID) da +1.278 en θ = 1. La referencia usada aquí es la lineal con build_A_v3 y todos los bins:
**+1.306** (θ = 1, σ_X = 0.171). Como estaba previsto en el encargo, el cálculo **se paró en el paso 1** porque con todos los bins no salía
~+1.4; el autor decidió seguir con +1.306 como referencia.

**Montaje.** El de marginal_guiada / respuesta_v3A: malla G0 de exacto (3506 nodos en u = ln(1+z), z de DESI como nodos exactos); datos = los
13 números de DESI DR2 con su covarianza oficial 13 × 13 + 100 θ* y ω_m de Planck [+ Pantheon+, 1580 SN, STAT+SYS]; prior x ~ N(0, σ_X² K_θ) con
K el núcleo OU en la malla. La respuesta A se construye con la **definición de `respuesta_v3A.build_A_v3`** (mismo código, con p0 y B0
inyectados; filas SN jacobiano exacto; filas suman cero) en el fiducial p0 **reajustado a las filas conservadas** con el optimizador original
de cada combinación: sin SN, `exacto.fit` (least_squares trf, cotas [0.4, 1.0] × [0.05, 0.30], x_scale [0.01, 0.001], ftol = xtol = 1e-10,
max_nfev 300, arranque (0.682, 0.1407), con la misma función `exacto.model`); con Pantheon+, `exacto_sn.fit_combo` (trf, cotas
[0.4, 1.0] × [0.05, 0.30] × [−50, 50], x_scale [0.01, 0.001, 0.01], ftol = xtol = 1e-10, max_nfev 200, M de arranque analítico), con el mismo
residuo blanqueado. Plantillas B por diferencias finitas del modelo exacto (ε = 1e-4, 1e-5, 1e-3), como `marginal_guiada.build_B`.
ln Z_lin(σ) = `marginal_guiada.lnZ_lin` con S = C + σ² A K Aᵀ en la malla de σ de `maximos_globales` (perfil_fondo.SIG hasta 0.05 y
geomspace(0.05, 1, 40): 87 nodos, de 1e-6 a 1); ln Z(0) en σ = 0 exactamente; máximo (valor y nodo) sobre la malla; cota por caída de 1.92 (y
1.35) desde el máximo global con `maximos_globales.bound_from` (brentq entre nodos). Hilos BLAS: los de numpy por defecto (16).

**Comprobaciones con todos los bins** (`checks` del json): el p0 con Pantheon+ reproduce `marginal_guiada.P0` **bit a bit** (diferencia 0.0,
χ² 1406.3271); el p0 sin SN reproduce `exacto_sn.P0` bit a bit (χ² 17.2378, el 17.24 de EXACTO.md); la A con todos los bins
es idéntica a `build_A_v3()` (diferencia 0.0) y el código con p0 inyectado reproduce `build_A_v3()` a 0.0; filas de A suman cero a
1.3e-15. Con todos los bins y Pantheon+, θ = 1: máximo +1.1403 en σ_X = 0.1078 (A8: 1.1402 en 0.1078) y cota
0.2772 (A30 cita 0.272 de `positividad_lin_v3A.json`, que interpola en ln σ entre celdas gruesas; `exacto_sn_lineal_v3A_maximos_v3A.json` da 0.2824 en
la malla de cmb_fondo con fiducial fijo: la diferencia es de malla y fiducial, no de datos).

**Los 7 bins y la covarianza.** Filas del vector de 13 (ficheros `desi_gaussian_bao_ALL_GCcomb_mean.txt` y `_cov.txt` del clon de
`CobayaSampler/bao_data`, `data/desi_dr2_local/dr2_data/dr2-bao-zenodo/cosmology_chains/bao_data/desi_bao_dr2/`, sha256 `9ac154ab…0640585` y
`252a1432…886e509`, iguales a los de `data/PROCEDENCIA_DESI_DR2.md` §2, comprobados hoy; leídos como números por `cuesta.load_gaussian_bao_vector`):
BGS D_V (fila 0, z = 0.295); LRG1 D_M, D_H (1–2, z = 0.510); LRG2 (3–4, z = 0.706); LRG3+ELG1 (5–6, z = 0.934); ELG2 (7–8, z = 1.321); QSO (9–10,
z = 1.484); Lyα D_H, D_M (11–12, z = 2.33). Total 13. **La covarianza oficial de DESI es diagonal por bloques entre bins:** las
144 entradas entre bins distintos son exactamente 0.0 (0 no nulas); dentro de cada bin la correlación D_M–D_H es
LRG1 -0.452, LRG2 -0.395, LRG3+ELG1 -0.347, ELG2 -0.398, QSO -0.494, Lya -0.431. Al retirar un bin se eliminan sus filas del vector de datos y se toma la **submatriz principal de S = C + σ² A K Aᵀ** (y de C, de B y
de A): eso es la distribución **marginal** de los datos conservados, **no la condicional** (complemento de Schur). Como C no correlaciona bins
distintos, en la parte de medida marginal y condicional coinciden; la diferencia está solo en el término σ² A K Aᵀ, que correlaciona todos los
bins a través del residuo, y ahí se usa la marginal.

## Tabla: máximo de ln Z(σ_X) − ln Z(0) (σ_X del máximo) / cota −1.92

| bin quitado | sin SN θ = 0.3 | sin SN θ = 1 | sin SN θ = 3 | sin SN θ = 10 | Pantheon+ θ = 0.3 | Pantheon+ θ = 1 | Pantheon+ θ = 3 | Pantheon+ θ = 10 |
|---|---|---|---|---|---|---|---|---|
| todos | +1.317 (0.316) / 0.844 | +1.306 (0.171) / 0.482 | +1.192 (0.116) / 0.318 | +0.809 (0.0856) / 0.248 | +1.168 (0.185) / 0.492 | +1.140 (0.108) / 0.277 | +0.981 (0.0583) / 0.174 | +0.591 (0.0316) / 0.118 |
| sin BGS | +1.912 (0.341) / 0.941 | +1.939 (0.199) / 0.535 | +1.893 (0.136) / 0.351 | +1.509 (0.0998) / 0.280 | +1.512 (0.215) / 0.543 | +1.511 (0.126) / 0.306 | +1.399 (0.0793) / 0.194 | +0.931 (0.05) / 0.139 |
| sin LRG1 | +1.334 (0.316) / 0.924 | +1.287 (0.185) / 0.533 | +1.122 (0.126) / 0.359 | +0.742 (0.0998) / 0.291 | +1.030 (0.215) / 0.549 | +0.961 (0.116) / 0.311 | +0.709 (0.063) / 0.198 | +0.379 (0.0251) / 0.131 |
| sin LRG2 | +0.005 (0.05) / 0.584 | +0.000 (0.005) / 0.335 | -0.000 (1e-06) / 0.218 | -0.000 (1e-06) / 0.159 | +0.168 (0.054) / 0.324 | +0.195 (0.0316) / 0.181 | +0.201 (0.0199) / 0.110 | +0.139 (0.0126) / 0.073 |
| sin LRG3+ELG1 | +1.096 (0.316) / 0.984 | +1.115 (0.185) / 0.562 | +1.099 (0.126) / 0.372 | +1.022 (0.0998) / 0.290 | +0.739 (0.171) / 0.509 | +0.745 (0.0924) / 0.286 | +0.680 (0.05) / 0.181 | +0.441 (0.0251) / 0.129 |
| sin ELG2 | +1.458 (0.341) / 0.950 | +1.476 (0.199) / 0.544 | +1.423 (0.136) / 0.363 | +1.045 (0.108) / 0.294 | +1.147 (0.199) / 0.528 | +1.124 (0.116) / 0.299 | +0.962 (0.068) / 0.192 | +0.560 (0.0316) / 0.130 |
| sin QSO | +1.341 (0.316) / 0.865 | +1.324 (0.185) / 0.495 | +1.203 (0.116) / 0.328 | +0.802 (0.0856) / 0.257 | +1.188 (0.199) / 0.498 | +1.159 (0.108) / 0.280 | +0.992 (0.0583) / 0.176 | +0.601 (0.0316) / 0.118 |
| sin Lyα | +1.611 (0.341) / 0.936 | +1.571 (0.199) / 0.538 | +1.376 (0.126) / 0.358 | +0.888 (0.0924) / 0.274 | +1.309 (0.215) / 0.530 | +1.254 (0.116) / 0.298 | +1.033 (0.063) / 0.186 | +0.600 (0.0316) / 0.122 |

χ²_ΛCDM del ajuste con X = 0 (sin SN / con Pantheon+): todos 17.24 / 1406.33; sin BGS 17.21 / 1406.24; sin LRG1 13.85 / 1403.01;
sin LRG2 10.11 / 1398.63; sin LRG3+ELG1 15.26 / 1404.01; sin ELG2 16.51 / 1405.51; sin QSO 16.84 / 1405.97;
sin Lyα 17.24 / 1406.33. Sin LRG2 y sin SN, el máximo en θ = 3 y 10 está en σ → 0 (clase "sigma->0"); en θ = 1 el máximo sobre la malla es +6.8e-6
(σ = 0.005) y en θ = 0.3 +0.005 (σ = 0.05): el perfil es plano hasta σ ≈ 0.05 y luego cae.

**Cambio del máximo (Δ) y de la cota −1.92 (relativo) al quitar cada bin, los 7 sin elegir a posteriori:**

Sin SN:

| θ | BGS | LRG1 | LRG2 | LRG3+ELG1 | ELG2 | QSO | Lyα | el que más mueve el máximo / la cota |
|---|---|---|---|---|---|---|---|---|
| 0.3 | +0.59 / +11 % | +0.02 / +9 % | -1.31 / -31 % | -0.22 / +17 % | +0.14 / +12 % | +0.02 / +2 % | +0.29 / +11 % | LRG2 / LRG2 |
| 1 | +0.63 / +11 % | -0.02 / +10 % | -1.31 / -30 % | -0.19 / +17 % | +0.17 / +13 % | +0.02 / +3 % | +0.26 / +12 % | LRG2 / LRG2 |
| 3 | +0.70 / +10 % | -0.07 / +13 % | -1.19 / -32 % | -0.09 / +17 % | +0.23 / +14 % | +0.01 / +3 % | +0.18 / +12 % | LRG2 / LRG2 |
| 10 | +0.70 / +13 % | -0.07 / +17 % | -0.81 / -36 % | +0.21 / +17 % | +0.24 / +18 % | -0.01 / +4 % | +0.08 / +10 % | LRG2 / LRG2 |

Con Pantheon+:

| θ | BGS | LRG1 | LRG2 | LRG3+ELG1 | ELG2 | QSO | Lyα | el que más mueve el máximo / la cota |
|---|---|---|---|---|---|---|---|---|
| 0.3 | +0.34 / +10 % | -0.14 / +12 % | -1.00 / -34 % | -0.43 / +3 % | -0.02 / +7 % | +0.02 / +1 % | +0.14 / +8 % | LRG2 / LRG2 |
| 1 | +0.37 / +10 % | -0.18 / +12 % | -0.95 / -35 % | -0.40 / +3 % | -0.02 / +8 % | +0.02 / +1 % | +0.11 / +8 % | LRG2 / LRG2 |
| 3 | +0.42 / +11 % | -0.27 / +14 % | -0.78 / -37 % | -0.30 / +4 % | -0.02 / +10 % | +0.01 / +1 % | +0.05 / +7 % | LRG2 / LRG2 |
| 10 | +0.34 / +18 % | -0.21 / +11 % | -0.45 / -38 % | -0.15 / +9 % | -0.03 / +10 % | +0.01 / +0 % | +0.01 / +3 % | LRG2 / LRG2 |

En las 8 filas (2 combinaciones × 4 θ) el bin que más mueve el máximo y el que más mueve la cota es **LRG2**. Los dos siguientes, sin SN
θ = 1: BGS (+0.63, quitarlo sube el máximo) y Lyα (+0.26); con Pantheon+ θ = 1: LRG3+ELG1 (-0.40) y BGS (+0.37).

## Veredicto

- **P12: sí.** Sin SN, θ = 1: con todos los bins +1.306; sin LRG2 +6.8e-06 (el máximo desaparece: ln Z es plana desde σ = 0 y
  decrece), muy por debajo de la mitad (0.653). La cota −1.92 pasa de 0.482 a 0.335 (-30 %).
- **P13: no.** Sin SN, θ = 1, los otros seis bins: BGS +0.63, LRG1 -0.02, LRG3+ELG1 -0.19, ELG2 +0.17, QSO +0.02,
  Lyα +0.26. Quitar BGS mueve el máximo más de 0.5 (lo sube: +1.94); los otros cinco quedan dentro de ±0.5.
- **P14: no.** Con Pantheon+, θ = 1: BGS +0.37, LRG1 -0.18, LRG2 -0.95, LRG3+ELG1 -0.40, ELG2 -0.02, QSO +0.02,
  Lyα +0.11. Quitar LRG2 mueve el máximo -0.95 (de +1.14 a +0.19); los otros seis quedan dentro de ±0.5. Predije que las SN
  diluirían el peso de LRG2 por debajo de 0.5 y no es así: con 1580 SN el máximo sigue dependiendo sobre todo de LRG2, aunque no desaparece.

Lectura corta: en nuestra marginal lineal, igual que en el e-valor de Kim, Mota y Tamosiunas para w0wa, la preferencia por σ_X > 0 sin SN
está en LRG2; con Pantheon+ queda un máximo de +0.2 sin LRG2. Los máximos son Δln Z ≈ 1 (Δχ² ≈ 2) y no excluyen σ_X = 0 en ningún caso.

## Tres cautelas

1. **Mirar en varios sitios.** Con 7 bins, que uno destaque por azar es esperable (efecto de mirar en varios sitios); esta prueba no corrige por
   ello, solo describe lo que pasa al quitar cada bin.
2. **Depender de un bin no es que el bin esté mal.** Que la preferencia dependa de LRG2 no demuestra que LRG2 esté mal medido: un bin que lleva
   la señal también se comporta así.
3. **Sobrevivir no es ser real.** Que un máximo sobreviva al quitar un bin (como con Pantheon+ al quitar cualquiera salvo LRG2) no demuestra que
   sea real: faltan la tasa de falsos positivos bajo ΛCDM de este estadístico jackknife, la sensibilidad a los supuestos (ℓ, familia, convenio,
   fiducial) y la recuperación de fluctuaciones inyectadas.

## Ficheros

`jackknife_bins.py` (script), `jackknife_bins.json` (todo: comprobaciones, covarianza, p0 y χ² por caso, perfiles ln Z − ln Z(0) en los 87
nodos, máximos, cotas, resumen, veredicto de P12–P14, referencias), `jackknife_bins_paso1.json` (solo todos los bins, lo que motivó la
parada), logs `jackknife_bins_run_paso1.log` y `jackknife_bins_run_todo.log`. Nada existente se modifica; todo cálculo usa build_A_v3.
