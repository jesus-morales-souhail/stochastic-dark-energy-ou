# Barrido: perfiles y cotas por caída desde un máximo, ¿dónde está el máximo? (exploración)

**Predicción, escrita antes de buscar:**
- Fuera de EXACTO_SN.md y CMB_FONDO.md no hay ningún perfil con el máximo en el
  borde superior de su malla que cambie una cota más de un 5 %.
- Los perfiles del paper (results/profile_sigma_x, FONDO, TRECE, SN) tienen el
  máximo en σ_X → 0 y no están afectados.
- Ampliar el rango de θ por debajo de 1e-3 cambia la cota de profile_sigma_x menos
  de un 2 %.

Exploración: no cambia `results/`, `figures/` ni ningún paper. Script `barrido_maximos.py` → `barrido_maximos.json` (238 entradas;
10 min). Recorre todos los json de `results/` y `analisis_cmb/` (sin `venv/` ni instantáneas `antes/`) que guardan un perfil o
una cota por caída desde un máximo, más los que guardan cotas de otro tipo (cruces de P, bondad de ajuste), que se listan como
"no aplica". Para cada entrada: fichero, malla (mínimo y máximo), dónde está el máximo, clase y, si es un fallo o el máximo o la
cota están en el borde de un rango explorado, la cota recalculada con la malla ampliada dos veces.

**Clases.** Interior; borde inferior físico (σ_X → 0: no es fallo); borde superior artificial (fallo). Para θ, μ y f: si la cota o
el máximo caen en el borde del rango explorado, se amplía el rango dos veces y se dice si la cota cambia (la segunda ampliación
tiene que dejar quieto el máximo). Recuento: 163 entradas con el máximo en σ_X → 0, 37 interiores, 35 en un borde
superior artificial (fallo), y el resto sin máximo que caracterizar (cruces de P, bondad de ajuste, ajustes).

## Resumen de lo que cambia

**1. Fallos: máximo recortado en el borde superior de la malla de `perfil_fondo.bound` (σ_X ≤ 0.05).** Ya corregidos
documentalmente en `EXACTO_SN.md`, `CMB_FONDO.md` y `EXACTO.md` (commit anterior, `maximos_globales.json`): BAO + CMB + SN lineal
(−5 % en θ = 1, −11 % en θ = 0.3), BAO + SN lineal (−1 a −16 %), `CMB_FONDO` marginalizado BAO (−0.5 a −7 %), BAO + θ* (≤ −1 %),
variante Λ de BAO + CMB (−6 a −21 %) y de BAO (−0.5 a −8 %). El barrido añade uno más del mismo tipo, fuera de esas notas:
`verificacion.json`, núcleo completo marginalizado con los 13 números (`VERIFICACION.md`, tabla "13 números", fila (c)): el
máximo está en σ_X = 0.108 / 0.073 / 0.058 (θ = 1 / 3 / 10), más allá de 0.05, y las cotas con el máximo global son
0.452 frente a 0.475 (-4.8 %, theta_1), 0.281 frente a 0.287 (-2.0 %, theta_3); con los 7 α el máximo está en 0.040 (dentro de la malla) y nada cambia.
El 0.475 de θ = 1 es el que `EXACTO.md` y `CMB_FONDO.md` usan como comprobación (0.4758 reproducido con el mismo máximo recortado):
la comprobación sigue siendo válida como reproducción, y el valor correcto de ambos es 0.452. Ninguna de estas cotas está en el
PREPRINT.

**2. Escalón: mejor ajuste o intervalo en el borde superior f = 3** (`escalon.json`, `ESCALON.md`), ampliado a [−2, 6] y [−2, 12]:

| celda | histórico (f̂, Δχ², 95 % sup.) | [−2, 6] | [−2, 12] | clase |
|---|---|---|---|---|
| A|zt_1.6|w_0.3 | 3.00, -3.15, 3.00 (borde) | 3.86, -3.28, 3.00 (borde) | 3.86, -3.28, 9.71 (interior) | borde superior artificial (fallo): el mejor ajuste estaba en f = 3 |
| A|zt_2|w_0.3 | 1.63, -0.96, 3.00 (borde) | 1.63, -0.96, 5.53 (interior) | 1.63, -0.96, 5.53 (interior) | interior; solo el intervalo superior tocaba el borde |
| B|zt_3|w_0.05 | 2.72, -6.62, 3.00 (borde) | 2.72, -6.62, 4.64 (interior) | 2.72, -6.62, 4.64 (interior) | interior; solo el intervalo superior tocaba el borde |
| B|zt_3|w_0.3 | 1.91, -6.64, 3.00 (borde) | 1.91, -6.64, 3.30 (interior) | 1.91, -6.64, 3.30 (interior) | interior; solo el intervalo superior tocaba el borde |

En (A) z_t = 1.6, w = 0.3 el mejor ajuste se desplaza a f = 3.86 y el Δχ² mejora de −3.15 a −3.28 (4 %); el intervalo superior
necesita la segunda ampliación (9.7). En las otras tres celdas el mejor ajuste no se mueve y el intervalo superior pasa de "borde"
a 5.5 / 4.6 / 3.3, estable en la segunda ampliación. El p global de `ESCALON.md` (0.26 en (A), 0.047 en (B)) no se ha
recalculado con el rango ampliado: el mejor Δχ² observado de (A) pasaría de −3.15 a −3.28 y los universos de ruido también
tendrían más rango; el cambio esperable es de centésimas en p.

**3. Bordes de rango de θ, μ, t_c (no son fallos de máximo, pero mueven cotas):**

| perfil | parámetro en el borde | cota histórica | ampliación 1 | ampliación 2 | ¿cambia? |
|---|---|---|---|---|---|
| `results/profile_sigma_x` (= FONDO (a), everpresent T0, desi_dr2_real_bao) | θ = 10⁻³ en la cota | 0.02633 | 0.02644 (θ ∈ (10⁻⁵, 200)) | 0.02644 (θ ∈ (10⁻⁷, 2000)) | +0.43 %, igual al supremo θ → 0 (0.02644) ya escrito en el PREPRINT (§3.2, "2.644 × 10⁻² (+0.4 %)") |
| `reloj/perfil_reloj.json`, reloj a_ln_a | θ = 10⁻³ en la cota | 0.02633 | 0.02644 | 0.02644 | +0.43 %: los cuatro relojes convergen al mismo 0.02644 |
| `reloj/perfil_reloj.json`, reloj b_cosmic_time | θ = 10⁻³ en la cota | 0.02639 | 0.02644 | 0.02644 | +0.19 %: los cuatro relojes convergen al mismo 0.02644 |
| `reloj/perfil_reloj.json`, reloj c_conformal_time | θ = 10⁻³ en la cota | 0.02632 | 0.02644 | 0.02644 | +0.45 %: los cuatro relojes convergen al mismo 0.02644 |
| `reloj/perfil_reloj.json`, reloj d_redshift | θ = 10⁻³ en la cota | 0.02616 | 0.02644 | 0.02644 | +1.06 %: los cuatro relojes convergen al mismo 0.02644 |
| `fondo/perfil_fondo.json` (b), (c) y `perfil_13.json` (a), (b), (c), θ libre | θ = 10⁻³ en la cota | ver tabla completa | ≤ +0.5 % | ≤ +0.5 % | no (dentro del 0.5 %; mismo supremo θ → 0) |
| `sn/perfil_sn.json`, θ libre | θ = 10⁻³ en la cota | 0.1167 | 0.1132 | 0.1132 | **−3.0 %**: el máximo sube (ln L 2290.007 → 2290.037) con θ < 10⁻³; estable en la segunda ampliación |
| `everpresent/profile_desi_everpresent.json`, Z2 (Zwane, μ libre) | μ = 0.1 en la cota (borde inferior de [0.1, 10]) | 0.02407 | 0.02612 ([0.01, 10]) | 0.02641 ([0.001, 10]) | **+9.7 %**; la segunda ampliación aún mueve +1.1 %, así que se amplió dos veces más: 0.02644 ([10⁻⁴, 10]), 0.02644 ([10⁻⁵, 10]); converge al límite μ → 0 de rango 1, 0.02644, que es el mismo supremo que el del OU en θ → 0 |
| `results/tachyonic_rank1/SUMMARY.json` | máximo en la meseta t_c → ∞ (borde superior 10⁴) | −12.7 "representativo" | — | — | borde físico (sin crecimiento = ΛCDM); rehecho en `VERIFICACION.md` |

**4. El resto está bien:** máximos en σ_X → 0 (perfiles del paper con el núcleo local: `profile_sigma_x`, FONDO (a)/(b)/(c), TRECE,
`cuesta` bounds_13, everpresent T1/T2/Z1, `ou_qnm_mle`, `joint_w0wa_sigma`, SN con θ fijo salvo θ = 10) o interiores
(SN θ libre y θ = 10 en σ_X ≈ 0.01, Fisher ecualizador con A_máx ≤ 0.0085, `precision_sn`, `convergencia_sn`, `marginal_guiada`,
`CMB_FONDO` BAO + CMB del convenio principal en 0.004, `verificacion` núcleo completo con 7 α en ≤ 0.040), y las cotas por cruce de
P o por bondad de ajuste (`exacto`, `exacto_sn`, `gof`, `malla`) no dependen de ningún máximo. Las entradas "sin cota (≤ 0.7)" de
las notas de la auditoría son bordes superiores explícitos de la malla de σ_X, no máximos.

**Cotas del PREPRINT.** Con el núcleo local y la malla de `profile_sigma_x_desi.py`, el máximo de todos los perfiles del PREPRINT
(§3.1, §3.2: 2.6 × 10⁻², la tabla por tratamiento y θ, TRECE) está en σ_X → 0 y las cotas no cambian con la malla de σ_X; con el
rango de θ ampliado cambian +0.4 % (el supremo θ → 0, ya escrito). **Una cifra citada en el PREPRINT sí cambia:** la cota de
Pantheon+ con θ libre, "σ_X only at 0.06–0.12" (§3.2, de `SN.md`), pasa de 0.1167 a 0.1132 (−3 %) al dejar θ bajar de 10⁻³
(el máximo de la verosimilitud de las SN está en θ < 10⁻³, y la cota se mide desde él); el 0.12 se convertiría en 0.11 y el
"2.4–4.5 veces más débil que DESI" en 2.4–4.3. No se corrige aquí: se para y se avisa.

## Tabla completa

**`results/profile_sigma_x/profile_sigma_x.json`**

| magnitud | malla | máximo en | clase | histórico | recalculado | nota |
|---|---|---|---|---|---|---|
| sigma_X 95 % (theta free) | sigma 1e-6..5e-2 (48), theta (1e-3, 20) | sigma = 1e-6 | borde inferior físico (σ_X → 0); θ en el borde inferior de su rango | 0.02633 | ext1: 0.02644; ext2: 0.02644 | theta at the bound = 1e-3 (lower edge); extending theta: 0.02633 -> 0.02644 (+0.42 %) -> 0.02644 (+0.43 %); theta->0 supremum 0.02644; max at sigma = 1e-6 where |

**`results/desi_dr2_real_bao/desi_dr2_real_bao.json`**

| magnitud | malla | máximo en | clase | histórico | recalculado | nota |
|---|---|---|---|---|---|---|
| sigma_X_95_profile_upper_approx (old grid node) | sigma 1e-6..5e-2, theta (1e-3, 20) | sigma = 1e-6 | borde inferior físico (σ_X → 0); θ en el borde inferior | 0.02175 | ext1: 0.02644; ext2: 0.02644 | same likelihood as profile_sigma_x; same theta extension |

**`analisis_cmb/everpresent/resultados/profile_desi_everpresent.json`**

| magnitud | malla | máximo en | clase | histórico | recalculado | nota |
|---|---|---|---|---|---|---|
| T0_ou_theta_profiled upper95_refined | sigma 1e-6..5e-2, theta (1e-3, 20) | sigma = 1e-6 | borde inferior físico (σ_X → 0); θ en el borde inferior | 0.02633 | ext1: 0.02644; ext2: 0.02644 | same likelihood as profile_sigma_x; same theta extension |
| Z2_zwane_mu_free upper95_refined (alpha) | alpha 1e-6..2 (80), mu [0.1, 10] | alpha = 1e-6 | borde inferior físico (σ_X → 0) (alpha); μ en la cota = 0.1 (borde inferior de su rango) | 0.02407 | ext1: 0.02612; ext2: 0.02641; ext3 [1e-4, 10]: 0.02644; ext4 [1e-5, 10]: 0.02644; mu->0: 0.02644 | base reproduces 0.02407; converges to the rank-1 (mu -> 0) limit |
| T1_ou_theta1 upper95_refined | amp 1e-6..0.05 or 2 (Z1) | amp = 1e-06 | borde inferior físico (σ_X → 0); nuisance interior or fixed | 0.01444 | — | nuisance at bound 1.0 |
| T2_ou_theta_fit upper95_refined | amp 1e-6..0.05 or 2 (Z1) | amp = 1e-06 | borde inferior físico (σ_X → 0); nuisance interior or fixed | 0.01294 | — | nuisance at bound 2.7197306671361585 |
| Z1_zwane_mu_fixed upper95_refined | amp 1e-6..0.05 or 2 (Z1) | amp = 1e-06 | borde inferior físico (σ_X → 0); nuisance interior or fixed | 0.01849 | — | nuisance at bound 0.9804 |

**`analisis_cmb/reloj/perfil_reloj.json`**

| magnitud | malla | máximo en | clase | histórico | recalculado | nota |
|---|---|---|---|---|---|---|
| clock a_ln_a: sigma_X 95 % (theta free) | sigma 1e-6..5e-2, theta (1e-3, 20) | sigma = 1e-6 | borde inferior físico (σ_X → 0); θ en el borde inferior | 0.02633 | ext1: 0.02644; ext2: 0.02644 | base reproduces 0.02633; ext changes +0.43 % |
| clock b_cosmic_time: sigma_X 95 % (theta free) | sigma 1e-6..5e-2, theta (1e-3, 20) | sigma = 1e-6 | borde inferior físico (σ_X → 0); θ en el borde inferior | 0.02639 | ext1: 0.02644; ext2: 0.02644 | base reproduces 0.02639; ext changes +0.19 % |
| clock c_conformal_time: sigma_X 95 % (theta free) | sigma 1e-6..5e-2, theta (1e-3, 20) | sigma = 1e-6 | borde inferior físico (σ_X → 0); θ en el borde inferior | 0.02632 | ext1: 0.02644; ext2: 0.02644 | base reproduces 0.02632; ext changes +0.45 % |
| clock d_redshift: sigma_X 95 % (theta free) | sigma 1e-6..5e-2, theta (1e-3, 20) | sigma = 1e-6 | borde inferior físico (σ_X → 0); θ en el borde inferior | 0.02616 | ext1: 0.02644; ext2: 0.02644 | base reproduces 0.02616; ext changes +1.06 % |

**`analisis_cmb/fondo/perfil_fondo.json`**

| magnitud | malla | máximo en | clase | histórico | recalculado | nota |
|---|---|---|---|---|---|---|
| a_fixed theta_free | sigma 1e-6..5e-2 (ext. to 1 for the crossing), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0); θ en el borde inferior | 0.02633 | base (1e-3, 20): 0.02633; ext1 (1e-5, 200): 0.02644; ext2 (1e-7, 2000): 0.02644 | theta free: base 0.026328275699446086 -> ext1 0.026439525495473317 -> ext2 0.02644065283048112 |
| a_fixed theta_0.001 | sigma 1e-6..5e-2 (ext. to 1 for the crossing), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.02633 | global_max: 0.02633 |  |
| a_fixed theta_0.1 | sigma 1e-6..5e-2 (ext. to 1 for the crossing), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.02101 | global_max: 0.02101 |  |
| a_fixed theta_0.3 | sigma 1e-6..5e-2 (ext. to 1 for the crossing), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | n/a | global_max: 0.01766 |  |
| a_fixed theta_1 | sigma 1e-6..5e-2 (ext. to 1 for the crossing), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.01444 | global_max: 0.01444 |  |
| a_fixed theta_3 | sigma 1e-6..5e-2 (ext. to 1 for the crossing), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | n/a | global_max: 0.01287 |  |
| a_fixed theta_10 | sigma 1e-6..5e-2 (ext. to 1 for the crossing), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.0129 | global_max: 0.0129 |  |
| b_profiled theta_free | sigma 1e-6..5e-2 (ext. to 1 for the crossing), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0); θ en el borde inferior | 0.02633 | base (1e-3, 20): 0.02633; ext1 (1e-5, 200): 0.02644; ext2 (1e-7, 2000): 0.02644 | theta free: base 0.02633115611426577 -> ext1 0.026439426999087835 -> ext2 0.02644052544708822 |
| b_profiled theta_0.001 | sigma 1e-6..5e-2 (ext. to 1 for the crossing), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.02633 | global_max: 0.02633 |  |
| b_profiled theta_0.1 | sigma 1e-6..5e-2 (ext. to 1 for the crossing), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.02129 | global_max: 0.02129 |  |
| b_profiled theta_0.3 | sigma 1e-6..5e-2 (ext. to 1 for the crossing), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | n/a | global_max: 0.01812 |  |
| b_profiled theta_1 | sigma 1e-6..5e-2 (ext. to 1 for the crossing), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.015 | global_max: 0.015 |  |
| b_profiled theta_3 | sigma 1e-6..5e-2 (ext. to 1 for the crossing), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | n/a | global_max: 0.0134 |  |
| b_profiled theta_10 | sigma 1e-6..5e-2 (ext. to 1 for the crossing), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.01327 | global_max: 0.01327 |  |
| c_marginalised theta_free | sigma 1e-6..5e-2 (ext. to 1 for the crossing), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0); θ en el borde inferior | sin cota | base (1e-3, 20): sin cota; ext1 (1e-5, 200): sin cota; ext2 (1e-7, 2000): sin cota | theta free: base None -> ext1 None -> ext2 None |
| c_marginalised theta_0.001 | sigma 1e-6..5e-2 (ext. to 1 for the crossing), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | sin cota | global_max: sin cota |  |
| c_marginalised theta_0.1 | sigma 1e-6..5e-2 (ext. to 1 for the crossing), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.1398 | global_max: 0.1398 |  |
| c_marginalised theta_0.3 | sigma 1e-6..5e-2 (ext. to 1 for the crossing), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | n/a | global_max: 0.08143 |  |
| c_marginalised theta_1 | sigma 1e-6..5e-2 (ext. to 1 for the crossing), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.04611 | global_max: 0.04611 |  |
| c_marginalised theta_3 | sigma 1e-6..5e-2 (ext. to 1 for the crossing), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | n/a | global_max: 0.02951 |  |
| c_marginalised theta_10 | sigma 1e-6..5e-2 (ext. to 1 for the crossing), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.02204 | global_max: 0.02204 |  |

**`analisis_cmb/fondo/perfil_13.json`**

| magnitud | malla | máximo en | clase | histórico | recalculado | nota |
|---|---|---|---|---|---|---|
| a_fixed theta_free | sigma 1e-6..5e-2 (ext. to 1 for the crossing), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0); θ en el borde inferior | 0.02305 | base (1e-3, 20): 0.02305; ext1 (1e-5, 200): 0.02315; ext2 (1e-7, 2000): 0.02315 | theta free: base 0.023046960120132635 -> ext1 0.02315133581892691 -> ext2 0.02315239543221778 |
| a_fixed theta_0.001 | sigma 1e-6..5e-2 (ext. to 1 for the crossing), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.02305 | global_max: 0.02305 |  |
| a_fixed theta_0.1 | sigma 1e-6..5e-2 (ext. to 1 for the crossing), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.01833 | global_max: 0.01833 |  |
| a_fixed theta_1 | sigma 1e-6..5e-2 (ext. to 1 for the crossing), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.0128 | global_max: 0.0128 |  |
| a_fixed theta_10 | sigma 1e-6..5e-2 (ext. to 1 for the crossing), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.0116 | global_max: 0.0116 |  |
| b_profiled theta_free | sigma 1e-6..5e-2 (ext. to 1 for the crossing), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0); θ en el borde inferior | 0.02305 | base (1e-3, 20): 0.02305; ext1 (1e-5, 200): 0.02315; ext2 (1e-7, 2000): 0.02315 | theta free: base 0.02304852053452651 -> ext1 0.023153317488662912 -> ext2 0.02315438300550139 |
| b_profiled theta_0.001 | sigma 1e-6..5e-2 (ext. to 1 for the crossing), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.02305 | global_max: 0.02305 |  |
| b_profiled theta_0.1 | sigma 1e-6..5e-2 (ext. to 1 for the crossing), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.01849 | global_max: 0.01849 |  |
| b_profiled theta_1 | sigma 1e-6..5e-2 (ext. to 1 for the crossing), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.01325 | global_max: 0.01325 |  |
| b_profiled theta_10 | sigma 1e-6..5e-2 (ext. to 1 for the crossing), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.01195 | global_max: 0.01195 |  |
| c_marginalised theta_free | sigma 1e-6..5e-2 (ext. to 1 for the crossing), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0); θ en el borde inferior | sin cota | base (1e-3, 20): sin cota; ext1 (1e-5, 200): sin cota; ext2 (1e-7, 2000): sin cota | theta free: base None -> ext1 None -> ext2 None |
| c_marginalised theta_0.001 | sigma 1e-6..5e-2 (ext. to 1 for the crossing), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | sin cota | global_max: sin cota |  |
| c_marginalised theta_0.1 | sigma 1e-6..5e-2 (ext. to 1 for the crossing), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.1121 | global_max: 0.1121 |  |
| c_marginalised theta_1 | sigma 1e-6..5e-2 (ext. to 1 for the crossing), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.03774 | global_max: 0.03774 |  |
| c_marginalised theta_10 | sigma 1e-6..5e-2 (ext. to 1 for the crossing), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.01925 | global_max: 0.01925 |  |

**`analisis_cmb/fondo/perfil_13.json (+ cuesta bounds_13)`**

| magnitud | malla | máximo en | clase | histórico | recalculado | nota |
|---|---|---|---|---|---|---|
| a_fixed theta_0.3 | sigma 1e-6..5e-2 (ext. to 1 for the crossing), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | n/a | global_max: 0.01548 |  |
| a_fixed theta_3 | sigma 1e-6..5e-2 (ext. to 1 for the crossing), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | n/a | global_max: 0.01155 |  |
| b_profiled theta_0.3 | sigma 1e-6..5e-2 (ext. to 1 for the crossing), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | n/a | global_max: 0.0158 |  |
| b_profiled theta_3 | sigma 1e-6..5e-2 (ext. to 1 for the crossing), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | n/a | global_max: 0.01202 |  |
| c_marginalised theta_0.3 | sigma 1e-6..5e-2 (ext. to 1 for the crossing), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | n/a | global_max: 0.06559 |  |
| c_marginalised theta_3 | sigma 1e-6..5e-2 (ext. to 1 for the crossing), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | n/a | global_max: 0.02494 |  |

**`analisis_cmb/auditoria_v23/verificacion.json`**

| magnitud | malla | máximo en | clase | histórico | recalculado | nota |
|---|---|---|---|---|---|---|
| 7alpha local a_fixed theta_0.001 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.02633 | global_max: 0.02633 |  |
| 7alpha local a_fixed theta_0.1 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.02101 | global_max: 0.02101 |  |
| 7alpha local a_fixed theta_1 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.01444 | global_max: 0.01444 |  |
| 7alpha local a_fixed theta_3 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.01287 | global_max: 0.01287 |  |
| 7alpha local a_fixed theta_10 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.0129 | global_max: 0.0129 |  |
| 7alpha local a_fixed theta_free | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0); θ en el borde inferior | 0.02633 | global_max: 0.02633 |  |
| 7alpha local b_profiled theta_0.001 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.02633 | global_max: 0.02633 |  |
| 7alpha local b_profiled theta_0.1 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.02129 | global_max: 0.02129 |  |
| 7alpha local b_profiled theta_1 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.015 | global_max: 0.015 |  |
| 7alpha local b_profiled theta_3 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.0134 | global_max: 0.0134 |  |
| 7alpha local b_profiled theta_10 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.01327 | global_max: 0.01327 |  |
| 7alpha local b_profiled theta_free | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0); θ en el borde inferior | 0.02633 | global_max: 0.02633 |  |
| 7alpha local c_marginalised theta_0.001 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | sin cota | global_max: sin cota |  |
| 7alpha local c_marginalised theta_0.1 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.1398 | global_max: 0.1398 |  |
| 7alpha local c_marginalised theta_1 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.04611 | global_max: 0.04611 |  |
| 7alpha local c_marginalised theta_3 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.02951 | global_max: 0.02951 |  |
| 7alpha local c_marginalised theta_10 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.02204 | global_max: 0.02204 |  |
| 7alpha local c_marginalised theta_free | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0); θ en el borde inferior | sin cota | global_max: sin cota |  |
| 7alpha full a_fixed theta_0.001 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.02645 | global_max: 0.02645 |  |
| 7alpha full a_fixed theta_0.1 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.02552 | global_max: 0.02552 |  |
| 7alpha full a_fixed theta_1 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.02154 | global_max: 0.02154 |  |
| 7alpha full a_fixed theta_3 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.01924 | global_max: 0.01924 |  |
| 7alpha full a_fixed theta_10 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.01855 | global_max: 0.01855 |  |
| 7alpha full a_fixed theta_free | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0); θ en el borde inferior | 0.02645 | global_max: 0.02645 |  |
| 7alpha full b_profiled theta_0.001 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.02645 | global_max: 0.02645 |  |
| 7alpha full b_profiled theta_0.1 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.02553 | global_max: 0.02553 |  |
| 7alpha full b_profiled theta_1 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.02161 | global_max: 0.02161 |  |
| 7alpha full b_profiled theta_3 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.01933 | global_max: 0.01933 |  |
| 7alpha full b_profiled theta_10 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.01857 | global_max: 0.01857 |  |
| 7alpha full b_profiled theta_free | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0); θ en el borde inferior | 0.02645 | global_max: 0.02645 |  |
| 7alpha full c_marginalised theta_0.001 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 0.251 | interior | sin cota | global_max: sin cota |  |
| 7alpha full c_marginalised theta_0.1 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 0.0998 | interior | sin cota | global_max: sin cota |  |
| 7alpha full c_marginalised theta_1 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 0.0397 | interior | 0.4967 | global_max: 0.4967 |  |
| 7alpha full c_marginalised theta_3 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 0.0316 | interior | 0.3185 | global_max: 0.3185 |  |
| 7alpha full c_marginalised theta_10 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 0.0316 | interior | 0.2441 | global_max: 0.2441 |  |
| 7alpha full c_marginalised theta_free | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 0.0397 | interior; θ en el borde inferior | sin cota | global_max: sin cota |  |
| 13 local a_fixed theta_0.001 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.02305 | global_max: 0.02305 |  |
| 13 local a_fixed theta_0.1 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.01833 | global_max: 0.01833 |  |
| 13 local a_fixed theta_1 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.0128 | global_max: 0.0128 |  |
| 13 local a_fixed theta_3 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.01155 | global_max: 0.01155 |  |
| 13 local a_fixed theta_10 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.0116 | global_max: 0.0116 |  |
| 13 local a_fixed theta_free | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0); θ en el borde inferior | 0.02305 | global_max: 0.02305 |  |
| 13 local b_profiled theta_0.001 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.02305 | global_max: 0.02305 |  |
| 13 local b_profiled theta_0.1 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.01849 | global_max: 0.01849 |  |
| 13 local b_profiled theta_1 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.01325 | global_max: 0.01325 |  |
| 13 local b_profiled theta_3 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.01202 | global_max: 0.01202 |  |
| 13 local b_profiled theta_10 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.01195 | global_max: 0.01195 |  |
| 13 local b_profiled theta_free | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0); θ en el borde inferior | 0.02305 | global_max: 0.02305 |  |
| 13 local c_marginalised theta_0.001 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | sin cota | global_max: sin cota |  |
| 13 local c_marginalised theta_0.1 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.1121 | global_max: 0.1121 |  |
| 13 local c_marginalised theta_1 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.03774 | global_max: 0.03774 |  |
| 13 local c_marginalised theta_3 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.02494 | global_max: 0.02494 |  |
| 13 local c_marginalised theta_10 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.01925 | global_max: 0.01925 |  |
| 13 local c_marginalised theta_free | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0); θ en el borde inferior | sin cota | global_max: sin cota |  |
| 13 full a_fixed theta_0.001 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.02314 | global_max: 0.02314 |  |
| 13 full a_fixed theta_0.1 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.02234 | global_max: 0.02234 |  |
| 13 full a_fixed theta_1 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.01915 | global_max: 0.01915 |  |
| 13 full a_fixed theta_3 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.0177 | global_max: 0.0177 |  |
| 13 full a_fixed theta_10 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.01802 | global_max: 0.01802 |  |
| 13 full a_fixed theta_free | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0); θ en el borde inferior | 0.02314 | global_max: 0.02314 |  |
| 13 full b_profiled theta_0.001 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.02315 | global_max: 0.02315 |  |
| 13 full b_profiled theta_0.1 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.02236 | global_max: 0.02236 |  |
| 13 full b_profiled theta_1 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.01934 | global_max: 0.01934 |  |
| 13 full b_profiled theta_3 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.01805 | global_max: 0.01805 |  |
| 13 full b_profiled theta_10 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.01824 | global_max: 0.01824 |  |
| 13 full b_profiled theta_free | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0); θ en el borde inferior | 0.02315 | global_max: 0.02315 |  |
| 13 full c_marginalised theta_0.001 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 0.584 | interior | sin cota | global_max: sin cota |  |
| 13 full c_marginalised theta_0.1 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 0.316 | interior | sin cota | global_max: sin cota |  |
| 13 full c_marginalised theta_1 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 0.108 | interior | 0.4751 | global_max: 0.4522 |  |
| 13 full c_marginalised theta_3 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 0.0734 | interior | 0.2866 | global_max: 0.2809 |  |
| 13 full c_marginalised theta_10 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 0.0583 | interior | 0.2032 | global_max: 0.2023 |  |
| 13 full c_marginalised theta_free | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 0.063 | interior; θ en el borde inferior | sin cota | global_max: sin cota |  |

**`analisis_cmb/sn/perfil_sn.json`**

| magnitud | malla | máximo en | clase | histórico | recalculado | nota |
|---|---|---|---|---|---|---|
| free_theta sigma_95 (Omega_m profiled) | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 0.0100 | interior (lnL +0.026); θ en la cota = 10⁻³ (borde inferior) | 0.1167 | ext1: 0.1132; ext2: 0.1132 | base reproduces 0.1167; the stored extension profile decreases monotonically beyond 0.05 (max interior) |
| fixed theta = 0.001: sigma_95 | sigma 1e-6..5e-2 (ext. to 1) | sigma = 1.99e-06 | borde inferior físico (σ_X → 0) | 0.1198 | — | stored extension profile decreases monotonically beyond 0.05 |
| fixed theta = 0.1: sigma_95 | sigma 1e-6..5e-2 (ext. to 1) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.1186 | — | stored extension profile decreases monotonically beyond 0.05 |
| fixed theta = 1.0: sigma_95 | sigma 1e-6..5e-2 (ext. to 1) | sigma = 1.26e-06 | borde inferior físico (σ_X → 0) | 0.09634 | — | stored extension profile decreases monotonically beyond 0.05 |
| fixed theta = 10.0: sigma_95 | sigma 1e-6..5e-2 (ext. to 1) | sigma = 0.00998 | interior | 0.06289 | — | stored extension profile decreases monotonically beyond 0.05 |

**`analisis_cmb/escalon/escalon.json`**

| magnitud | malla | máximo en | clase | histórico | recalculado | nota |
|---|---|---|---|---|---|---|
| (A) zt_1.6|w_0.3: f_best, dchi2, 95 % upper | f in [-2, 3] step 0.05 | f = 3.000 | borde superior artificial (f = 3): fallo | f_best: 3; dchi2: -3.152; hi95: 3 | base [-2, 3]: f_best: 3; dchi2: -3.152; hi95: 3; ext1 [-2, 6]: f_best: 3.857; dchi2: -3.281; hi95: 3; ext2 [-2, 12]: f_best: 3.857; dchi2: -3.281; hi95: 9.706 | global p-value of ESCALON.md not recomputed with the extended f range |
| (A) zt_2|w_0.3: f_best, dchi2, 95 % upper | f in [-2, 3] step 0.05 | f = 1.635 | interior; intervalo 95 % superior en el borde | f_best: 1.635; dchi2: -0.9582; hi95: 3 | base [-2, 3]: f_best: 1.635; dchi2: -0.9582; hi95: 3; ext1 [-2, 6]: f_best: 1.635; dchi2: -0.9582; hi95: 5.531; ext2 [-2, 12]: f_best: 1.635; dchi2: -0.9582; hi95: 5.531 | global p-value of ESCALON.md not recomputed with the extended f range |
| (B) zt_3|w_0.05: f_best, dchi2, 95 % upper | f in [-2, 3] step 0.05 | f = 2.715 | interior; intervalo 95 % superior en el borde | f_best: 2.715; dchi2: -6.619; hi95: 3 | base [-2, 3]: f_best: 2.715; dchi2: -6.619; hi95: 3; ext1 [-2, 6]: f_best: 2.715; dchi2: -6.619; hi95: 4.639; ext2 [-2, 12]: f_best: 2.715; dchi2: -6.619; hi95: 4.639 | global p-value of ESCALON.md not recomputed with the extended f range |
| (B) zt_3|w_0.3: f_best, dchi2, 95 % upper | f in [-2, 3] step 0.05 | f = 1.914 | interior; intervalo 95 % superior en el borde | f_best: 1.914; dchi2: -6.642; hi95: 3 | base [-2, 3]: f_best: 1.914; dchi2: -6.642; hi95: 3; ext1 [-2, 6]: f_best: 1.914; dchi2: -6.642; hi95: 3.301; ext2 [-2, 12]: f_best: 1.914; dchi2: -6.642; hi95: 3.301 | global p-value of ESCALON.md not recomputed with the extended f range |

**`analisis_cmb/auditoria_v23/cmb_fondo.json`**

| magnitud | malla | máximo en | clase | histórico | recalculado | nota |
|---|---|---|---|---|---|---|
| BAO profiled theta_0.001 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.02316 | global_max_192: 0.02316; global_max_135: 0.01278 | change +0.0 % |
| BAO profiled theta_0.1 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.02237 | global_max_192: 0.02237; global_max_135: 0.01254 | change +0.0 % |
| BAO profiled theta_0.3 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.02126 | global_max_192: 0.02126; global_max_135: 0.01215 | change +0.0 % |
| BAO profiled theta_1.0 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.01936 | global_max_192: 0.01936; global_max_135: 0.01132 | change +0.0 % |
| BAO profiled theta_3.0 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.01808 | global_max_192: 0.01808; global_max_135: 0.01053 | change +0.0 % |
| BAO profiled theta_10.0 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.01827 | global_max_192: 0.01827; global_max_135: 0.01023 | change +0.0 % |
| BAO profiled theta_free | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0); θ en el borde inferior | 0.02316 | global_max_192: 0.02316; global_max_135: 0.01278 | change +0.0 % |
| BAO marginalised theta_0.001 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 0.584 | borde superior artificial del máximo (0.05): fallo | sin cota | global_max_192: sin cota; global_max_135: sin cota |  |
| BAO marginalised theta_0.1 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 0.316 | borde superior artificial del máximo (0.05): fallo | sin cota | global_max_192: sin cota; global_max_135: sin cota |  |
| BAO marginalised theta_0.3 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 0.185 | borde superior artificial del máximo (0.05): fallo | 0.872 | global_max_192: 0.8072; global_max_135: 0.6602 | change -7.4 % |
| BAO marginalised theta_1.0 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 0.108 | borde superior artificial del máximo (0.05): fallo | 0.4758 | global_max_192: 0.4524; global_max_135: 0.3706 | change -4.9 % |
| BAO marginalised theta_3.0 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 0.0734 | borde superior artificial del máximo (0.05): fallo | 0.2868 | global_max_192: 0.2809; global_max_135: 0.2313 | change -2.0 % |
| BAO marginalised theta_10.0 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 0.0583 | borde superior artificial del máximo (0.05): fallo | 0.2033 | global_max_192: 0.2023; global_max_135: 0.1683 | change -0.5 % |
| BAO marginalised theta_free | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 0.063 | borde superior artificial del máximo (0.05): fallo | sin cota | global_max_192: sin cota; global_max_135: sin cota |  |
| BAO+theta profiled theta_0.001 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.001585 | global_max_192: 0.001585; global_max_135: 0.0008753 | change +0.0 % |
| BAO+theta profiled theta_0.1 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.00158 | global_max_192: 0.00158; global_max_135: 0.0008723 | change +0.0 % |
| BAO+theta profiled theta_0.3 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.001569 | global_max_192: 0.001569; global_max_135: 0.000867 | change +0.0 % |
| BAO+theta profiled theta_1.0 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.001544 | global_max_192: 0.001544; global_max_135: 0.0008536 | change +0.0 % |
| BAO+theta profiled theta_3.0 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.001511 | global_max_192: 0.001511; global_max_135: 0.0008355 | change +0.0 % |
| BAO+theta profiled theta_10.0 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.001486 | global_max_192: 0.001486; global_max_135: 0.0008209 | change +0.0 % |
| BAO+theta profiled theta_free | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0); θ en el borde inferior | 0.001585 | global_max_192: 0.001585; global_max_135: 0.0008753 | change +0.0 % |
| BAO+theta marginalised theta_0.001 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | sin cota | global_max_192: sin cota; global_max_135: sin cota |  |
| BAO+theta marginalised theta_0.1 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 0.116 | borde superior artificial del máximo (0.05): fallo | sin cota | global_max_192: sin cota; global_max_135: 0.8573 |  |
| BAO+theta marginalised theta_0.3 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 0.0998 | borde superior artificial del máximo (0.05): fallo | 0.647 | global_max_192: 0.6414; global_max_135: 0.5249 | change -0.9 % |
| BAO+theta marginalised theta_1.0 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 0.0734 | borde superior artificial del máximo (0.05): fallo | 0.3743 | global_max_192: 0.3719; global_max_135: 0.306 | change -0.6 % |
| BAO+theta marginalised theta_3.0 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 0.0583 | borde superior artificial del máximo (0.05): fallo | 0.2422 | global_max_192: 0.2416; global_max_135: 0.2003 | change -0.2 % |
| BAO+theta marginalised theta_10.0 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 0.05 | borde superior artificial del máximo (0.05): fallo | 0.1824 | global_max_192: 0.1824; global_max_135: 0.1526 | change +0.0 % |
| BAO+theta marginalised theta_free | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 0.05 | borde superior artificial del máximo (0.05): fallo | sin cota | global_max_192: sin cota; global_max_135: sin cota |  |
| BAO+omega_m profiled theta_0.001 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.02729 | global_max_192: 0.02729; global_max_135: 0.01498 | change +0.0 % |
| BAO+omega_m profiled theta_0.1 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.02607 | global_max_192: 0.02607; global_max_135: 0.01462 | change +0.0 % |
| BAO+omega_m profiled theta_0.3 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.02444 | global_max_192: 0.02444; global_max_135: 0.01404 | change +0.0 % |
| BAO+omega_m profiled theta_1.0 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.02176 | global_max_192: 0.02176; global_max_135: 0.01284 | change +0.0 % |
| BAO+omega_m profiled theta_3.0 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.01987 | global_max_192: 0.01987; global_max_135: 0.0117 | change +0.0 % |
| BAO+omega_m profiled theta_10.0 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.02053 | global_max_192: 0.02053; global_max_135: 0.01151 | change +0.0 % |
| BAO+omega_m profiled theta_free | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0); θ en el borde inferior | 0.02729 | global_max_192: 0.02729; global_max_135: 0.01498 | change +0.0 % |
| BAO+omega_m marginalised theta_0.001 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.04216 | global_max_192: 0.04216; global_max_135: 0.0233 | change +0.0 % |
| BAO+omega_m marginalised theta_0.1 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.0428 | global_max_192: 0.0428; global_max_135: 0.0233 | change +0.0 % |
| BAO+omega_m marginalised theta_0.3 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.04414 | global_max_192: 0.04414; global_max_135: 0.02332 | change +0.0 % |
| BAO+omega_m marginalised theta_1.0 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.04924 | global_max_192: 0.04924; global_max_135: 0.02361 | change +0.0 % |
| BAO+omega_m marginalised theta_3.0 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.06076 | global_max_192: 0.06076; global_max_135: 0.02509 | change -0.0 % |
| BAO+omega_m marginalised theta_10.0 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.0636 | global_max_192: 0.0636; global_max_135: 0.02848 | change -0.0 % |
| BAO+omega_m marginalised theta_free | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0); θ en el borde inferior | 0.06456 | global_max_192: 0.06456; global_max_135: 0.02853 | change +0.0 % |
| BAO+CMB profiled theta_0.001 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 0.00316 | interior | 0.04097 | global_max_192: 0.04097; global_max_135: 0.02296 | change +0.0 % |
| BAO+CMB profiled theta_0.1 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 0.00316 | interior | 0.03626 | global_max_192: 0.03626; global_max_135: 0.02184 | change +0.0 % |
| BAO+CMB profiled theta_0.3 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 0.00316 | interior | 0.0325 | global_max_192: 0.0325; global_max_135: 0.02051 | change +0.0 % |
| BAO+CMB profiled theta_1.0 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 0.00316 | interior | 0.0287 | global_max_192: 0.0287; global_max_135: 0.01876 | change +0.0 % |
| BAO+CMB profiled theta_3.0 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 0.00316 | interior | 0.02786 | global_max_192: 0.02786; global_max_135: 0.01817 | change +0.0 % |
| BAO+CMB profiled theta_10.0 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 0.00316 | interior | 0.03199 | global_max_192: 0.03199; global_max_135: 0.02 | change +0.0 % |
| BAO+CMB profiled theta_free | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 0.00316 | interior; θ en el borde inferior | 0.04097 | global_max_192: 0.04097; global_max_135: 0.02296 | change +0.0 % |
| BAO+CMB marginalised theta_0.001 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 0.00397 | interior | 0.04488 | global_max_192: 0.04488; global_max_135: 0.02513 | change +0.0 % |
| BAO+CMB marginalised theta_0.1 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 0.00397 | interior | 0.04559 | global_max_192: 0.04559; global_max_135: 0.02523 | change +0.0 % |
| BAO+CMB marginalised theta_0.3 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 0.00397 | interior | 0.04707 | global_max_192: 0.04707; global_max_135: 0.02546 | change +0.0 % |
| BAO+CMB marginalised theta_1.0 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 0.00397 | interior | 0.0524 | global_max_192: 0.0524; global_max_135: 0.0264 | change +0.0 % |
| BAO+CMB marginalised theta_3.0 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 0.00397 | interior | 0.06182 | global_max_192: 0.06182; global_max_135: 0.02894 | change -0.0 % |
| BAO+CMB marginalised theta_10.0 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 0.00397 | interior | 0.06264 | global_max_192: 0.06264; global_max_135: 0.03193 | change -0.0 % |
| BAO+CMB marginalised theta_free | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 0.00397 | interior; θ en el borde inferior | 0.06188 | global_max_192: 0.06188; global_max_135: 0.03079 | change -0.0 % |
| BAO+CMB|Lambda_absorbs_X0 profiled theta_0.001 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | sin cota | global_max_192: sin cota; global_max_135: 0.5371 |  |
| BAO+CMB|Lambda_absorbs_X0 profiled theta_0.1 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.12 | global_max_192: 0.12; global_max_135: 0.05468 | change +0.0 % |
| BAO+CMB|Lambda_absorbs_X0 profiled theta_0.3 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.07191 | global_max_192: 0.07191; global_max_135: 0.03269 | change -0.0 % |
| BAO+CMB|Lambda_absorbs_X0 profiled theta_1.0 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.04413 | global_max_192: 0.04413; global_max_135: 0.01997 | change +0.0 % |
| BAO+CMB|Lambda_absorbs_X0 profiled theta_3.0 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.0319 | global_max_192: 0.0319; global_max_135: 0.01457 | change +0.0 % |
| BAO+CMB|Lambda_absorbs_X0 profiled theta_10.0 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.02542 | global_max_192: 0.02542; global_max_135: 0.0123 | change +0.0 % |
| BAO+CMB|Lambda_absorbs_X0 profiled theta_free | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1e-06 | borde inferior físico (σ_X → 0); θ en el borde inferior | sin cota | global_max_192: sin cota; global_max_135: 0.5371 |  |
| BAO+CMB|Lambda_absorbs_X0 marginalised theta_0.001 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1 | borde superior artificial del máximo (0.05): fallo | sin cota | global_max_192: sin cota; global_max_135: sin cota |  |
| BAO+CMB|Lambda_absorbs_X0 marginalised theta_0.1 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 0.541 | borde superior artificial del máximo (0.05): fallo | sin cota | global_max_192: sin cota; global_max_135: sin cota |  |
| BAO+CMB|Lambda_absorbs_X0 marginalised theta_0.3 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 0.316 | borde superior artificial del máximo (0.05): fallo | sin cota | global_max_192: 0.8667; global_max_135: 0.7358 |  |
| BAO+CMB|Lambda_absorbs_X0 marginalised theta_1.0 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 0.185 | borde superior artificial del máximo (0.05): fallo | 0.6237 | global_max_192: 0.4955; global_max_135: 0.4211 | change -20.6 % |
| BAO+CMB|Lambda_absorbs_X0 marginalised theta_3.0 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 0.116 | borde superior artificial del máximo (0.05): fallo | 0.3782 | global_max_192: 0.3271; global_max_135: 0.2786 | change -13.5 % |
| BAO+CMB|Lambda_absorbs_X0 marginalised theta_10.0 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 0.0856 | borde superior artificial del máximo (0.05): fallo | 0.2694 | global_max_192: 0.2544; global_max_135: 0.2168 | change -5.6 % |
| BAO+CMB|Lambda_absorbs_X0 marginalised theta_free | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 0.293 | borde superior artificial del máximo (0.05): fallo | sin cota | global_max_192: sin cota; global_max_135: sin cota |  |
| BAO|Lambda_absorbs_X0 marginalised theta_0.001 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 1 | borde superior artificial del máximo (0.05): fallo | sin cota | global_max_192: sin cota; global_max_135: sin cota |  |
| BAO|Lambda_absorbs_X0 marginalised theta_0.1 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 0.341 | borde superior artificial del máximo (0.05): fallo | sin cota | global_max_192: sin cota; global_max_135: sin cota |  |
| BAO|Lambda_absorbs_X0 marginalised theta_0.3 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 0.199 | borde superior artificial del máximo (0.05): fallo | 0.8898 | global_max_192: 0.8179; global_max_135: 0.6694 | change -8.1 % |
| BAO|Lambda_absorbs_X0 marginalised theta_1.0 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 0.116 | borde superior artificial del máximo (0.05): fallo | 0.4822 | global_max_192: 0.4564; global_max_135: 0.3742 | change -5.3 % |
| BAO|Lambda_absorbs_X0 marginalised theta_3.0 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 0.0734 | borde superior artificial del máximo (0.05): fallo | 0.2895 | global_max_192: 0.2829; global_max_135: 0.2331 | change -2.3 % |
| BAO|Lambda_absorbs_X0 marginalised theta_10.0 | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 0.0583 | borde superior artificial del máximo (0.05): fallo | 0.2044 | global_max_192: 0.2033; global_max_135: 0.1692 | change -0.5 % |
| BAO|Lambda_absorbs_X0 marginalised theta_free | sigma 1e-6..5e-2 (ext. to 1), theta (1e-3, 20) | sigma = 0.063 | borde superior artificial del máximo (0.05): fallo | sin cota | global_max_192: sin cota; global_max_135: sin cota |  |

**`analisis_cmb/auditoria_v23/exacto_sn.json`**

| magnitud | malla | máximo en | clase | histórico | recalculado | nota |
|---|---|---|---|---|---|---|
| linear SN theta_0.1 | sigma 1e-6..5e-2 (ext. to 1) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.7724 | global_max_192: 0.7724; global_max_135: 0.5397 | change +0.0 % |
| linear SN theta_0.3 | sigma 1e-6..5e-2 (ext. to 1) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.4467 | global_max_192: 0.4467; global_max_135: 0.3121 | change +0.0 % |
| linear SN theta_1 | sigma 1e-6..5e-2 (ext. to 1) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.2464 | global_max_192: 0.2464; global_max_135: 0.1722 | change +0.0 % |
| linear SN theta_3 | sigma 1e-6..5e-2 (ext. to 1) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.1465 | global_max_192: 0.1465; global_max_135: 0.1028 | change +0.0 % |
| linear SN theta_10 | sigma 1e-6..5e-2 (ext. to 1) | sigma = 1e-06 | borde inferior físico (σ_X → 0) | 0.09205 | global_max_192: 0.09205; global_max_135: 0.06585 | change +0.0 % |
| linear BAO+SN theta_0.1 | sigma 1e-6..5e-2 (ext. to 1) | sigma = 0.199 | borde superior artificial del máximo (0.05): fallo | 0.8327 | global_max_192: 0.703; global_max_135: 0.5912 | change -15.6 % |
| linear BAO+SN theta_0.3 | sigma 1e-6..5e-2 (ext. to 1) | sigma = 0.116 | borde superior artificial del máximo (0.05): fallo | 0.4472 | global_max_192: 0.4086; global_max_135: 0.344 | change -8.6 % |
| linear BAO+SN theta_1 | sigma 1e-6..5e-2 (ext. to 1) | sigma = 0.068 | borde superior artificial del máximo (0.05): fallo | 0.2327 | global_max_192: 0.2297; global_max_135: 0.1941 | change -1.3 % |
| linear BAO+SN theta_3 | sigma 1e-6..5e-2 (ext. to 1) | sigma = 0.0397 | interior | 0.1447 | global_max_192: 0.1447; global_max_135: 0.1233 | change +0.0 % |
| linear BAO+SN theta_10 | sigma 1e-6..5e-2 (ext. to 1) | sigma = 0.0316 | interior | 0.1043 | global_max_192: 0.1043; global_max_135: 0.08956 | change +0.0 % |
| linear BAO+CMB+SN theta_0.1 | sigma 1e-6..5e-2 (ext. to 1) | sigma = 0.341 | borde superior artificial del máximo (0.05): fallo | sin cota | global_max_192: 0.8643; global_max_135: 0.7542 |  |
| linear BAO+CMB+SN theta_0.3 | sigma 1e-6..5e-2 (ext. to 1) | sigma = 0.185 | borde superior artificial del máximo (0.05): fallo | 0.5665 | global_max_192: 0.5022; global_max_135: 0.4381 | change -11.4 % |
| linear BAO+CMB+SN theta_1 | sigma 1e-6..5e-2 (ext. to 1) | sigma = 0.108 | borde superior artificial del máximo (0.05): fallo | 0.2982 | global_max_192: 0.2824; global_max_135: 0.2462 | change -5.3 % |
| linear BAO+CMB+SN theta_3 | sigma 1e-6..5e-2 (ext. to 1) | sigma = 0.0583 | borde superior artificial del máximo (0.05): fallo | 0.1778 | global_max_192: 0.1771; global_max_135: 0.1537 | change -0.4 % |
| linear BAO+CMB+SN theta_10 | sigma 1e-6..5e-2 (ext. to 1) | sigma = 0.0316 | interior | 0.119 | global_max_192: 0.119; global_max_135: 0.1009 | change +0.0 % |

**`analisis_cmb/fondo/fisher_diseno.json`**

| magnitud | malla | máximo en | clase | histórico | recalculado | nota |
|---|---|---|---|---|---|---|
| c_equaliser desi_real: A_95 (40 periods) | A in [0, 1] (4001) | A_at_max <= 0.0085 | borde inferior físico (σ_X → 0) / interior near 0 | A_95 max 0.085 | — | no A_95 at the edge 1 |
| c_equaliser desi_expected: A_95 (40 periods) | A in [0, 1] (4001) | A_at_max <= 0.0000 | borde inferior físico (σ_X → 0) / interior near 0 | A_95 max 0.066 | — | no A_95 at the edge 1 |
| c_equaliser survey_N20: A_95 (40 periods) | A in [0, 1] (4001) | A_at_max <= 0.0000 | borde inferior físico (σ_X → 0) / interior near 0 | A_95 max 0.087 | — | no A_95 at the edge 1 |

**`analisis_cmb/auditoria_v23/precision_sn.json`**

| magnitud | malla | máximo en | clase | histórico | recalculado | nota |
|---|---|---|---|---|---|---|
| marginal MC bounds (4 combos x 5 theta) | sigma 0.01..0.7 (12) | lnL max at sigma = 0 (SN) or 0.03-0.2 (others) | interior / borde inferior físico (σ_X → 0); 'sin cota (<= 0.7)' entries are explicit (upper edge = sin cota, not a max) | see table | — | none |

**`analisis_cmb/auditoria_v23/convergencia_sn.json`**

| magnitud | malla | máximo en | clase | histórico | recalculado | nota |
|---|---|---|---|---|---|---|
| marginal MC bounds vs N | sigma 0.03..0.4 (0.7 for theta <= 0.3) | lnL max at sigma = 0.03-0.4 (interior) | interior | see table | — | none |

**`analisis_cmb/auditoria_v23/marginal_guiada.json`**

| magnitud | malla | máximo en | clase | histórico | recalculado | nota |
|---|---|---|---|---|---|---|
| exact guided marginal bounds | sigma 0..0.4 (0.7 for theta = 0.3) | lnZ max at sigma = 0.03-0.2 (interior) | interior | see table | — | none |

**`analisis_cmb/auditoria_v23/exacto.json / exacto_sn.json`**

| magnitud | malla | máximo en | clase | histórico | recalculado | nota |
|---|---|---|---|---|---|---|
| 5 % bounds (P crossings) | sigma 0.05..0.7 / 0.01..0.7 | not a drop from a maximum | no aplica (cruce de P; 'sin cota (<= 0.7)' explicit) | see tables | — | none |

**`results/tachyonic_rank1/SUMMARY.json`**

| magnitud | malla | máximo en | clase | histórico | recalculado | nota |
|---|---|---|---|---|---|---|
| scan dlnL(t_c) | t_c 1e-4..1e4 (500) | plateau dlnL = 0 for t_c >= 0.09 up to the upper edge | borde físico (no growth = LCDM; amplitude -> 0); redone in VERIFICACION.md | -12.7 'representative' (not a bound) | — | none |

**`results/ou_qnm_mle/ou_qnm_mle.json`**

| magnitud | malla | máximo en | clase | histórico | recalculado | nota |
|---|---|---|---|---|---|---|
| MLE (sigma_X, theta, omega_R) | ln sigma_X^2 in [-20, -1], theta >= 1e-3, omega_R [0, 10] | sigma_X at the lower bound 4.54e-5 | borde inferior físico (σ_X → 0) (physical lower edge); theta, omega_R interior | no bound (MLE) | — | none |

**`results/joint_w0wa_sigma/joint_w0wa_sigma.json`**

| magnitud | malla | máximo en | clase | histórico | recalculado | nota |
|---|---|---|---|---|---|---|
| joint fits | sigma_X free | sigma_X = 0 at the best fit | borde inferior físico (σ_X → 0) | no bound | — | none |

**`analisis_cmb/resultados/takahashi/malla/malla.json`**

| magnitud | malla | máximo en | clase | histórico | recalculado | nota |
|---|---|---|---|---|---|---|
| bound = largest sigma with p95(dchi2) < 4 | sigma 1e-4..2.506e-2 (10) | p95 < 4 on the whole grid | sin cota on the grid (explicit); no es caída desde un máximo | none | — | none |

**`analisis_cmb/cuesta/gof.json`**

| magnitud | malla | máximo en | clase | histórico | recalculado | nota |
|---|---|---|---|---|---|---|
| GOF bound P(chi2 <= 10.27) = 5 % | sigma 0.003..0.15 (16) | not a drop from a maximum | no aplica (cruce de P) | see GOF.md | — | none |

## Predicción

- **Fuera de EXACTO_SN.md y CMB_FONDO.md no hay ningún perfil con el máximo en el borde superior de su malla que cambie una cota
  más de un 5 %: cumplida, por poco.** Los hay: `verificacion.json` (13 números, núcleo completo, marginalizado: −4.8 % en θ = 1,
  −2.0 % en θ = 3, −0.5 % en θ = 10) y las celdas del escalón con f̂ = 3 (Δχ² −3.15 → −3.28, 4 %; intervalos superiores que
  pasan de "borde" a 9.7 / 5.5 / 4.6 / 3.3). `EXACTO.md` cambia hasta −21 %, pero sus cotas lineales son las de `CMB_FONDO.md`.
- **Los perfiles del paper (profile_sigma_x, FONDO, TRECE, SN) tienen el máximo en σ_X → 0 y no están afectados: cumplida para
  profile_sigma_x, FONDO y TRECE; no para SN.** El máximo de la SN con θ libre está en σ_X = 0.010 (interior, ln L +0.03), no en
  σ_X → 0, y su cota cambia −3 % al ampliar el rango de θ.
- **Ampliar θ por debajo de 10⁻³ cambia la cota de profile_sigma_x menos de un 2 %: cumplida** (+0.43 %, estable en la segunda
  ampliación, igual al supremo θ → 0).
