# PROPUESTA_V3: esqueleto de la v3.0 (sin prosa)

Creado el 2026-11-12. No toca PREPRINT.md ni CLAIMS.md. Cada afirmación lleva un identificador **A<n>**, su nivel (1 observado, 2 deducido, 3 supuesto, 4 no se puede saber), el número copiado del json donde vive (ruta bajo `analisis_cmb/auditoria_v23/` salvo que se diga otra cosa, con su clave) y la nota y el commit que la respaldan. Las cifras de cabecera son la marginal lineal −1.92 (`positividad_lin.json`, `marginal_guiada.json`); la regla del 5 % solo aparece en la fe de erratas (Apéndice A).

## Título provisional y mensaje

**Título provisional:** How much stochastic dark-energy dynamics can hide under ΛCDM: bounds from DESI DR2 BAO, compressed Planck and Type Ia supernovae.

**Mensaje (una frase):** cuánta dinámica estocástica de la energía oscura cabe escondida bajo ΛCDM con DESI DR2 + Planck comprimido + supernovas, sin afirmar que exista.

## Secciones y afirmaciones

### 0 Nota de versión (qué cambia respecto a v2.2.1)

- **A1** (nivel 2, deducido) — La v2.2.1 contenía dos errores que afectan a resultados: la exclusión del taquiónico (C4) no se sostiene y las cotas de θ ≥ 1 usaban un núcleo válido solo para θ → 0 — número: 128 filas de errata ← `manuscript/ERRATAS_V3.md` (recuento por dos caminos) — nota ERRATAS_V3.md (726cb3d, fced750).
- **A2** (nivel 2, deducido) — El −12.7 de C4 es el punto del barrido más cercano a un objetivo escrito en el código, con amplitud del modo 2.7 × 10³ en α — número: dlnL = −12.6789, max_amp_alpha = 2743.83 ← `verificacion.json` → tachyonic.scan_near_minus_11.35 — nota VERIFICACION.md (334e8f4).
- **A3** (nivel 2, deducido) — El núcleo local δα_i = S_i X(z_i) subestima la cota de θ ≥ 1 frente a la respuesta completa en ×1.49 (7 α) y ×1.50 (13 números) con fondo fijo — número: 1.4911 y 1.4955 ← `verificacion.json` → full_kernel.7alpha.ratio_full_over_local.a_fixed.theta_1, full_kernel.13.ratio_full_over_local.a_fixed.theta_1 — nota VERIFICACION.md (334e8f4).
- **A4** (nivel 2, deducido) — La cifra 2.6 × 10⁻² de C2 es el término de Occam del modo de Ω_m (σ → 0, θ en el borde inferior), no una cota del temblor — número: 0.026328 ← `verificacion.json` → threshold.192.a_fixed.theta_free; clase 'sigma->0; theta at lower edge' ← `barrido_maximos.json` → rows[0].class — notas VERIFICACION.md (334e8f4), BARRIDO_MAXIMOS.md (838c1c9).
- **A5** (nivel 2, deducido) — La etiqueta '95 % CL' de la caída de 1.92 cubre ≥ 97.5 % como cota unilateral — número: 0.976 (θ = 1, fijo), 0.9945 (θ → 0, fijo) ← `verificacion.json` → coverage.cases.theta_1|192|fixed.coverage, coverage.cases.theta_0.001|192|fixed.coverage — nota VERIFICACION.md (334e8f4).
- **A6** (nivel 3, supuesto) — Las cifras de cabecera de la v3.0 son las de la marginal lineal analítica con caída de 1.92 (convenio B); la regla del 5 % y la marginal exacta guiada quedan fuera de la cabecera (fe de erratas y comprobación) — decisión declarada en POSITIVIDAD.md, sección Versión 3 (33783e1).

### 1 Introducción

- **A7** (nivel 3, supuesto) — Pregunta del paper: cuánta dinámica estocástica de la energía oscura cabe bajo ΛCDM con DESI DR2 + Planck comprimido + supernovas, sin afirmar que exista — formulación; sin número.
- **A8** (nivel 1, observado) — Con los datos reales la marginal lineal tiene su máximo en σ_X > 0 con Δln L = +1.14 (θ = 1), por debajo de 1.92: no excluye σ_X = 0 — número: lnZ_max_rel = 1.1402 en σ_X = 0.1078 ← `muestras_sn.json` → linear.pantheon|OU|theta_1.lnZ_max_rel, .sigma_at_max — nota MUESTRAS_SN.md (bbf5a96).
- **A9** (nivel 2, deducido) — Esa mejora de ajuste se lee como tensión entre datos absorbida, no como detección (regla fijada antes de ver MUESTRAS_SN) — regla ← MUESTRAS_SN.md, cabecera (bbf5a96); respaldo: INYECCION.md, escenario C (e35d4c7).

### 2 Modelo (OU en convenio B; convenio A como sensibilidad; límites: OU no derivable, PPF, conservación, sin fluctuaciones espaciales)

- **A10** (nivel 3, supuesto) — Residuo X(u) OU en u = ln(1+z) con ρ_DE = Ω_Λ + X(z) − X(0) (convenio B: Λ absorbe X(0)); el convenio A (X(0) libre) solo como sensibilidad — definición ← EXACTO.md (d2a7487), EXACTO_SN.md (7686a0d).
- **A11** (nivel 3, supuesto) — Supuestos declarados: el OU no se deriva de ninguna teoría; PPF para el CMB; conservación de ρ_DE con w efectivo; sin fluctuaciones espaciales — supuestos; sin número — notas CMB_FINAL.md (8b9621a), CMB_TECHOS.md (6587cde).
- **A12** (nivel 2, deducido) — Las historias del posterior en el techo de θ = 1 cruzan w = −1 con |w + 1| máximo mediano 1.94 (p95 3.90): el modelo implica cruce fantasma en casi todas — número: 1.9379 / 3.8989 ← `cmb_final.json` → cells.theta_1.faithful.abs_wp1_max.p50, .p95 — nota CMB_FINAL.md (8b9621a).
- **A13** (nivel 2, deducido) — En el techo exacto de θ = 1 el 44.9 % de las historias del posterior lineal tiene ρ_DE ≤ 0 en algún z (8.3 % en θ = 3, 0 % en θ = 10) — número: 0.4492 / 0.0828 / 0.0 ← `cmb_final.json` → prep.cells.theta_1.frac_rho_le_0_raw, …theta_3…, …theta_10… — nota CMB_FINAL.md (8b9621a).
- **A14** (nivel 3, supuesto) — Exigir ρ_DE > 0 es un prior, no una ley: la relatividad general admite ρ_DE < 0 con H² > 0, y H² > 0 no basta para la admisibilidad física — texto fijado ← POSITIVIDAD.md, cabecera (25d7a65).

### 3 Datos y procedencia

- **A15** (nivel 1, observado) — BAO DESI DR2: 7 α_DV de figure6 del paquete Zenodo 16644577 (v1.0, 2025-08-20, CC-BY-4.0; md5 del zip idéntico al publicado) y vector gaussiano de 13 números con covarianza 13 × 13 del clon CobayaSampler/bao_data (commit bb0c1c9), no incluido en el zip — números: md5 425623648b4b12515fa1943a5857fd46; sha256 de cada fichero ← `data/PROCEDENCIA_DESI_DR2.md` §1–2 (e4c766e).
- **A16** (nivel 1, observado) — Planck comprimido: 100 θ* = 1.0411 ± 0.00031 y ω_m = 0.14301 ± 0.0011 como dos filas gaussianas — número ← `cmb_fondo.json` → planck.theta100, planck.sigma, planck.omega_m, planck.sigma_omega_m — nota CMB_FONDO.md (9504bf7).
- **A17** (nivel 1, observado) — Pantheon+: 1580 SN (z_HD > 0.01, sin calibradores), covarianza STAT+SYS; Ω_m con SN solas 0.3324 frente a 0.334 ± 0.018 publicado — número ← `muestras_sn.json` → omega_check.pantheon.Omega_m_fit, .Omega_m_published; sha256 y licencia ← `datos/FUENTES.json` → pantheon_plus — nota MUESTRAS_SN.md (bbf5a96, 618d148).
- **A18** (nivel 1, observado) — DES-SN5YR (1829), Dovekie (1820) y Union3 (22 nodos) reproducen su Ω_m publicado a menos de 0.002 — número: diff −0.0002 / +0.0002 / −0.0002 ← `muestras_sn.json` → omega_check.des5yr.diff, .dovekie.diff, .union3.diff — nota MUESTRAS_SN.md (bbf5a96).
- **A19** (nivel 1, observado) — DES-SN5YR y Dovekie comparten 1718 SN: Dovekie es una prueba de calibración, no una muestra independiente — número: 1718 ← `muestras_sn.json` → shared_sne.des5yr_vs_dovekie — nota MUESTRAS_SN.md (618d148).
- **A20** (nivel 1, observado) — Licencias: Union3 MIT; DES-SN5YR CC-BY-4.0 (Zenodo 12720778); Dovekie y Pantheon+ sin licencia declarada, fuera del repositorio con script de descarga y sha256 — ← `datos/FUENTES.json` → *.licencia; `descargar_datos.py` (618d148).
- **A21** (nivel 1, observado) — Planck 2018 TT ℓ < 30 y bandas φφ: ficheros del Planck Legacy Archive con sha256 en `analisis_cmb/README.md` (ignorados por git) — ← `analisis_cmb/README.md` §'Downloading the data again'.

### 4 Método (respuesta exacta, marginal lineal con caída de 1.92, validación con universos falsos, positividad)

- **A22** (nivel 2, deducido) — La respuesta lineal completa A (filas que suman cero: Λ absorbe X(0)) reproduce el χ² del modelo exacto en p0 con diferencia 0 y sus filas suman cero a 10⁻¹⁵ — número ← `marginal_guiada.json` → checks.difference, checks.A_rows_sum_max_abs — nota MARGINAL_GUIADA.md (dcdfc82).
- **A23** (nivel 2, deducido) — Marginal lineal analítica: ln Z_lin(σ) con C + σ² A K Aᵀ y priors planos en (h, ω_m, M); cota por caída de 1.92 desde el máximo global, interpolada en ln σ — definición ← `marginal_guiada.py` lnZ_lin, `precision_sn.py` lnl_bound (dcdfc82, ed5cbf7).
- **A24** (nivel 1, observado) — Validación con universos falsos: la marginal −1.92 cubre σ_true en el 95.5 / 96.5 / 97.5 % (σ_true = 0.1 / 0.2 / 0.3, θ = 1) — número ← `inyeccion.json` → scenarios.B_0.1.coverage_marg, B_0.2, B_0.3 — nota INYECCION.md (e35d4c7).
- **A25** (nivel 1, observado) — La versión lineal reproduce la cota exacta del 5 % y la marginal exacta con los datos reales a −8 % y −3 % — número: 0.1559 / 0.169 y 0.2771 / 0.287 ← `inyeccion.json` → real_check.5pct_linear, real_check.marg_linear, real_check.ref — nota INYECCION.md (e35d4c7, 840f1da).
- **A26** (nivel 2, deducido) — La marginal exacta guiada (Matheron + pesos exacto/lineal) es inestable cuando hay muestras con ρ_DE ≤ 0: en θ = 3, celda 0.2, ESS = 2, y la cota sin positividad pasa de 0.177 a 0.212 solo con otras semillas — número: 0.2122 frente a 0.1773 ← `positividad.json` → thetas.theta_3.bound_exact_192, .reference_exact_192_marginal_guiada — nota POSITIVIDAD.md v2 (416996d).
- **A27** (nivel 2, deducido) — Positividad como prior: Z_pos,lin = Z_lin · E_q[1(ρ_DE > 0)/P(σ, Ω_Λ)], P del prior del camino al Ω_Λ propio de cada muestra (malla de 15 nodos, ninguna muestra fuera) — número: n_outside = 0 de 215000 ← `positividad_lin.json` → omega_l_grid.n_outside, .n_samples_total — nota POSITIVIDAD.md v3 (33783e1).
- **A28** (nivel 2, deducido) — Ser positiva no hace a una historia lineal: entre las positivas de la celda 0.3 de θ = 1 el cociente exacto/lineal tiene mediana 0.45 (16–84 %: 0.039–1.57) — número ← `positividad_lin.json` → thetas.theta_1.bound_cells[1].w_positive_p16_p50_p84 — nota POSITIVIDAD.md v3 (33783e1).
- **A29** (nivel 2, deducido) — Regla del 5 % (P(χ²_min ≤ χ²_obs) = 5 %): no es un intervalo de confianza y no tiene cobertura nominal; empírica 92.5 / 86.5 / 85.0 % — número ← `inyeccion.json` → scenarios.B_0.1.coverage_5pct, B_0.2, B_0.3 — nota INYECCION.md (840f1da). Solo aparece en la fe de erratas.

### 5 Resultados (cotas por θ, familias, muestras de SN, CMB, positividad)

- **A30** (nivel 1, observado) — Cotas de cabecera (marginal lineal −1.92, BAO + θ* + ω_m + Pantheon+, convenio B): σ_X < 0.492 / 0.272 / 0.173 / 0.115 para θ = 0.3 / 1 / 3 / 10 — número: 0.4916 / 0.2719 / 0.1727 / 0.1151 ← `positividad_lin.json` → thetas.theta_0.3.bound_lin_192, theta_1, theta_3, theta_10 (= `marginal_guiada.json` → thetas.*.lin_192) — notas POSITIVIDAD.md v3 (33783e1), MARGINAL_GUIADA.md (dcdfc82).
- **A31** (nivel 1, observado) — Con caída de 1.35 (unilateral 95 % en el límite gaussiano): 0.427 / 0.234 / 0.152 / 0.100 — número ← `positividad_lin.json` → thetas.*.bound_lin_135 — nota POSITIVIDAD.md v3 (33783e1).
- **A32** (nivel 2, deducido) — Comprobación con la marginal exacta guiada: 0.500 / 0.287 / 0.177 / 0.116 (remuestreo 0.282–0.293 en θ = 1), con la inestabilidad de A-método — número ← `marginal_guiada.json` → thetas.*.exact_192, .exact_192_boot_p16_p84 — nota MARGINAL_GUIADA.md (dcdfc82).
- **A33** (nivel 1, observado) — Con positividad (ρ_DE > 0 en todo z): 0.461 / 0.263 / 0.174 / 0.116; cambio −6.3 / −3.5 / +0.8 / +0.7 % — número ← `positividad_lin.json` → thetas.*.bound_pos_lin_192, .ratio_pos_over_lin — nota POSITIVIDAD.md v3 (33783e1).
- **A34** (nivel 1, observado) — P_prior(ρ_DE > 0 | σ = 0.263, Ω_Λ = 0.696) = 0.732 ± 0.003 — número ← `positividad_lin.json` → prediction_3.P_prior_pos, .err_mc — nota POSITIVIDAD.md v3 (33783e1).
- **A35** (nivel 1, observado) — Solo BAO (13 números) y BAO + θ* + ω_m, marginal lineal −1.92, θ = 1: 0.482 y 0.624 — número ← `exacto.json` → bounds.theta_1.BAO.linear_marginalised_192, bounds.theta_1.BAO+CMB.linear_marginalised_192 — nota EXACTO.md (d2a7487); pendiente: confirmar frente a `maximos_globales.json` (ebda524) el máximo global.
- **A36** (nivel 1, observado) — Solo SN (Pantheon+), marginal lineal −1.92: 0.447 / 0.246 / 0.147 / 0.092 — número ← `exacto_sn.json` → linear.SN.theta_0.3.192, theta_1, theta_3, theta_10 — nota EXACTO_SN.md (7686a0d).
- **A37** (nivel 2, deducido) — El modo θ → 0 es un desplazamiento de Ω_m que el fondo libre absorbe: con fondo marginalizado BAO solo no acota θ → 0 y BAO + ω_m lo acota en 0.042 — número: 0.0422 ← `cmb_fondo.json` → fits.BAO+omega_m.marginalised.theta_0.001.192 — nota CMB_FONDO.md (9504bf7); pendiente: máximo global ← `maximos_globales.json` → cmb_fondo (ebda524).
- **A38** (nivel 1, observado) — Familias con el mismo montaje (ℓ = 1): Matérn 3/2 0.737, escalones 0.268, gaussiano sin cota (marginal lineal −1.92); OU 0.277 — número ← `familias.json` → linear.Matern32|l_1.192, linear.Steps|l_1.192, linear.Gauss|l_1.192 (null), linear.OU|l_1.192 — nota FAMILIAS.md (0fb9ee3).
- **A39** (nivel 2, deducido) — No hay cifra universal de Δρ independiente del modelo: la rms de Δρ en la cota entre familias varía ×2.4–3.4 en z ≈ 0.5 y ×1.86 en z = 0.2 (5 %), ×1.18 en z = 0.1 (−1.92) — número ← `familias.json` → universality.5pct.l_0.3_1_3.ratio_vs_z, .ratio_at_z_star; universality.192.l_0.3_1_3.ratio_at_z_star — nota FAMILIAS.md (0fb9ee3).
- **A40** (nivel 1, observado) — Muestras de SN, marginal lineal −1.92, OU θ = 1: Pantheon+ 0.277, DES-SN5YR 0.290, Dovekie 0.275, Union3 0.387 — número ← `muestras_sn.json` → linear.pantheon|OU|theta_1.192, des5yr, dovekie, union3 — nota MUESTRAS_SN.md (bbf5a96).
- **A41** (nivel 2, deducido) — La dispersión de Δρ(0.1) en la cota entre familias es ×1.24 (Pantheon+), ×1.38 (DES-SN5YR), ×1.42 (Dovekie), ×1.88 (Union3): no es una cifra universal — número ← `muestras_sn.json` → drho01_spread.pantheon|192.ratio, des5yr|192, dovekie|192, union3|192 — nota MUESTRAS_SN.md (bbf5a96).
- **A42** (nivel 1, observado) — CMB completo (TT ℓ < 30, φφ, picos) sobre 40 historias del posterior por techo: Δχ² total mediana −0.15 / −0.23 / −0.23, p95 0.18 / 0.11 / 0.03, ninguna > 1 — número ← `cmb_final.json` → cells.theta_1.faithful.total.p50, .p95; theta_3; theta_10; .frac_total_gt_1 — nota CMB_FINAL.md (8b9621a).
- **A43** (nivel 1, observado) — 120 de 120 historias fieles en H² con el borde arreglado (error máximo 4.0 × 10⁻⁴); el máximo de |Δ/σ| cae en φφ en el 97.5 / 95 / 87.5 % — número ← `cmb_final.json` → cells.*.n_faithful, cells.theta_10.faithful.rho_err_H2.max, cells.*.faithful.where_max_abs_z.phiphi — nota CMB_FINAL.md (8b9621a).
- **A44** (nivel 2, deducido) — Con 40 historias sin ningún Δχ² > 4 la probabilidad de un caso así queda por debajo de ~7.5 % (regla de 3/N, orientativa con pesos desiguales; ESS 6–16) — número ← `cmb_final_pesos.json` → cells.*.ess, .rule_of_three_upper — nota CMB_FINAL.md (367975f).
- **A45** (nivel 1, observado) — Planck 2018 (TT ℓ < 30 + lente) no distingue el residuo de ΛCDM en 2.633 × 10⁻²: p95(Δχ²) = 0.067 (TT) — número ← `analisis_cmb/RESULTADO.md` tabla, fila σ_X = 2.633e-2, θ = 1 (f041018); ficheros de salida ← `analisis_cmb/resultados/log_takahashi_*.txt`; pendiente: no hay json con esa cifra en git.
- **A46** (nivel 1, observado) — Everpresent Λ: 0/20000 (modelo 2) y 0/300000 + 0/200000 (modelo 1) con el criterio final — número ← `analisis_cmb/everpresent/resultados/rd_real.json` → model2.c_rd_real_window_wm_2sigma.nA, model1_alpha0.01…, model1_alpha0.005… — notas EVERPRESENT_vs_DESI.md, PROCEDENCIA_RESULTADOS.md (a65c4f8, e4c766e). Sin cambio respecto a v2.2.1.
- **A47** (nivel 2, deducido) — Taquiónico: la cota es de amplitud del modo que crece en tiempo cósmico, no una exclusión del modo; el punto 'activo' mejor tiene Δln L ≈ 0 con amplitud 4 × 10⁻⁶ — número ← `verificacion.json` → tachyonic.scan_best_active.dlnL, .max_amp_alpha — nota VERIFICACION.md (334e8f4).

### 6 Tensión frente a temblor

- **A48** (nivel 1, observado) — Tensión Δχ² (conjunto − 15 filas − SN solas) de ΛCDM por muestra: Pantheon+ 2.68, Dovekie 3.22, Union3 4.15, DES-SN5YR 8.73 — número ← `muestras_sn_lcdm.json` → pantheon.tension.delta_chi2, dovekie, union3, des5yr — nota MUESTRAS_SN.md §2b (618d148).
- **A49** (nivel 1, observado) — El máximo de ln L en σ_X > 0 sigue a la calibración: DES-SN5YR +4.86, Dovekie +1.83 con 1718 SN en común; Pantheon+ +1.14; Union3 +3.61 (θ = 1) — número ← `muestras_sn.json` → linear.des5yr|OU|theta_1.lnZ_max_rel, dovekie, pantheon, union3 — nota MUESTRAS_SN.md §2a (618d148).
- **A50** (nivel 1, observado) — En universos falsos con tensión tipo DES-SN5YR y sin temblor, la marginal −1.92 se afloja ×1.50 y la caída de 1.92 excluye σ_X = 0 en el 43 % (detección falsa) — número ← `inyeccion.json` → scenarios.C.ratio_vs_A.marg, scenarios.C.frac_excl0 — nota INYECCION.md (e35d4c7).
- **A51** (nivel 2, deducido) — Ninguna de las dos reglas distingue tensión de temblor; la hipótesis 'la marginal sigue la información y el 5 % mezcla información y tensión' no se sostiene — número: 1.504 frente a 1.713 ← `inyeccion.json` → scenarios.C.ratio_vs_A.marg, .5pct — nota INYECCION.md (e35d4c7).
- **A52** (nivel 2, deducido) — El p-valor del χ²_min conjunto (≈ 0.999 en tres muestras) no mide tensión: lo domina la normalización de la covarianza de SN — número ← `muestras_sn_lcdm.json` → pantheon.p_value, des5yr.p_value, dovekie.p_value — nota MUESTRAS_SN.md §2e (73aa74a).
- **A53** (nivel 1, observado) — Pantheon+ real cae en el percentil 97 de los universos ΛCDM puros para la cota lineal del 5 %, 97.5 para la marginal y 91.5 para la tensión — número ← INYECCION.md §'Los datos reales frente al escenario A' (e35d4c7), calculado de `inyeccion.json` → per_universe.A.b5, .marg, .tension y real_check.

### 7 Lo que no se puede saber

- **A54** (nivel 4, no se puede saber) — θ no se mide: el perfil con θ libre cae al borde inferior (θ → 0), que es el modo de fondo; las cotas se dan por θ fijado — ← `barrido_maximos.json` → rows[*].class = 'sigma->0; theta at lower edge' (838c1c9); `marginal_guiada.json` → thetas.*.lnZ_lin_max_at_sigma (dcdfc82).
- **A55** (nivel 4, no se puede saber) — No hay cifra puntual de Δρ independiente del modelo (familias y muestras) — respaldo: A-universalidad (familias.json → universality) y A-dispersión (muestras_sn.json → drho01_spread).
- **A56** (nivel 4, no se puede saber) — Las historias con ρ_DE ≤ 0 (45 % del techo de θ = 1) no pueden juzgarse con el CMB en PPF: el test no dice nada de ellas — número ← `cmb_final.json` → prep.cells.theta_1.frac_rho_le_0_raw — nota CMB_FINAL.md (8b9621a).
- **A57** (nivel 4, no se puede saber) — La marginal exacta en θ = 1 y θ = 0.3 con positividad no es evaluable por importancia (ESS < 1 % de N): haría falta muestrear el posterior truncado — número: 0.89 % y 0.65 % ← `positividad.json` → thetas.theta_1.bound_cells[1].ess_v_frac_N, thetas.theta_0.3.bound_cells[0].ess_v_frac_N — nota POSITIVIDAD.md v2 (416996d).
- **A58** (nivel 4, no se puede saber) — Si la mejora de ajuste con σ_X > 0 es tensión entre conjuntos o dinámica real no se puede decidir con estos datos (misma firma en el escenario C y en el B de INYECCION) — ← `inyeccion.json` → scenarios.C.frac_excl0 (0.43), scenarios.B_0.1.frac_excl0 (0.585) — nota INYECCION.md (e35d4c7).
- **A59** (nivel 4, no se puede saber) — La compresión en 22 tramos no reproduce la pérdida de información de Union3 (Fisher de Ω_m 0.01745 → 0.01754): qué hace a Union3 menos informativa queda sin probar aquí — número ← INYECCION.md, predicción 4 (e35d4c7).

### 8 Discusión física (cruce de w = −1, constante que cambia de signo, sin mecanismo microfísico)

- **A60** (nivel 2, deducido) — Cruce de w = −1: casi todas las historias admitidas en el techo cruzan w = −1 (|w + 1| máx ≥ 1.02 en las 40 de θ = 1) — número ← `cmb_final.json` → cells.theta_1.faithful.abs_wp1_max.min — nota CMB_FINAL.md (8b9621a). Cita por verificar: Vikman 2005 (no cruce en un campo escalar único).
- **A61** (nivel 3, supuesto) — Constante que cambia de signo: las historias con ρ_DE < 0 a z alto que el posterior admite son del tipo 'Λ con signo cambiante'; se discute con la literatura (Akarsu et al.) sin afirmar que exista — respaldo: `cmb_final.json` → prep.cells.theta_1.rho_min_over_ol_p05_p50 (−0.633, 0.048) (8b9621a).
- **A62** (nivel 4, no se puede saber) — No hay mecanismo microfísico que produzca el OU: el modelo es fenomenológico y la amplitud acotada no se conecta con ninguna escala de la teoría — sin número.
- **A63** (nivel 3, supuesto) — Everpresent Λ no es ese mecanismo con los datos actuales: sus dos implementaciones publicadas no ajustan BAO + condiciones CMB — respaldo: A-Everpresent (rd_real.json).

### 9 Pronóstico y datos futuros

- **A64** (nivel 3, supuesto) — Qué datos moverían las cotas: covarianzas BAO completas entre tramos, SN con calibración consolidada (DES-SN5YR frente a Dovekie: la preferencia sigue a la calibración), más tramos en z < 0.3 — sin número; respaldo cualitativo MUESTRAS_SN.md (618d148).
- **A65** (nivel 4, no se puede saber) — No se da pronóstico numérico de cotas futuras: no hay montaje validado con universos falsos para Euclid/DESI DR3 en este repositorio (los mocks Euclid de v2.2.1 están retirados, CLAIMS N9).

### 10 Conclusiones

- **A66** (nivel 1, observado) — Cabe bajo ΛCDM un temblor OU de hasta σ_X ≈ 0.27 (θ = 1) en la densidad de energía oscura relativa a la crítica, 0.17 (θ = 3) y 0.12 (θ = 10), con cobertura ≥ 95 % validada por inyección — números ← A-cabecera y A-cobertura.
- **A67** (nivel 2, deducido) — Los datos no requieren el temblor: el máximo en σ_X > 0 (+1.1 en Pantheon+, +4.9 en DES-SN5YR) es tensión absorbida que sigue a la calibración — ← A-máximos por muestra.
- **A68** (nivel 2, deducido) — La positividad de ρ_DE cambia qué historias se admiten (hasta el 45 % del techo de θ = 1), no cuánto σ_X permiten los datos (≤ 6.3 %) — ← A-positividad.
- **A69** (nivel 2, deducido) — El CMB que tenemos (comprimido y completo) no mueve estas cotas; lo que las fija es ω_m de Planck con los BAO y las SN — ← A-CMB completo, A-ω_m.

### Apéndice A: fe de erratas v2.2.1 → v3.0

Lista completa en la sección 5 de este documento (128 filas de `manuscript/ERRATAS_V3.md`); la tabla del apéndice se copia de allí con la columna v3.0. Aquí también vivirán las cotas de la regla del 5 % (exacto.json, exacto_sn.json, convergencia_sn.json) como cifras históricas de la v2.3, no de cabecera.

### Apéndice B: tabla de afirmaciones

La tabla se genera de la lista A1–A69 de arriba (identificador, nivel, número, fichero → clave, nota, commit). Sin afirmaciones nuevas.

## Reglas que la propuesta debe cumplir y comprobar

- Cada regla estadística con su nombre completo y ninguna llamada 'estándar': 'marginal lineal analítica con caída de 1.92 desde el máximo global (cobertura unilateral ≥ 97.5 % en el límite gaussiano; 95.5–97.5 % por inyección)', 'caída de 1.35', 'regla del 5 % (P(χ²_min ≤ χ²_obs) = 5 %, sin cobertura nominal)', 'marginal exacta guiada por importancia (Matheron)'. Comprobación: grep de 'standard|estándar|95 % CL' en la v3.0 debe dar cero aciertos sin la aclaración.
- Una mejora de ajuste con σ_X > 0 se presenta como tensión entre datos absorbida, nunca como detección. Comprobación: toda frase con 'prefer', 'detect', 'evidence' va acompañada del Δln L máximo y de la lectura de INYECCION (escenario C: 43 % de falsas detecciones).
- No hay cifra puntual de Δρ independiente del modelo. Comprobación: ninguna frase da 'Δρ = x' sin familia, ℓ, z y regla; A-universalidad y A-dispersión lo respaldan.
- θ no se mide: las cotas se dan por θ fijado (0.3, 1, 3, 10) y nunca se escribe 'θ = …' como resultado. Comprobación: grep de 'theta =' y 'θ =' en frases de resultado.
- Ningún número sin fichero: cada cifra del texto tiene fichero json/txt seguido por git y clave, y una nota con commit. Comprobación: script que extrae los números del texto y los busca en los json citados (como las relecturas de las notas de auditoria_v23).
- Las cifras de cabecera son de la marginal lineal −1.92; la exacta guiada solo como comprobación con su inestabilidad declarada (ESS 2 en θ = 3, 0.2); la regla del 5 % solo en el Apéndice A.
- Datos: cada fichero de entrada con URL oficial, fecha, sha256 y licencia (`data/PROCEDENCIA_DESI_DR2.md`, `datos/FUENTES.json`); nada se redistribuye sin licencia declarada.
- Predicciones prerregistradas: cada test cita su predicción escrita antes de calcular y su resultado punto por punto (como las notas); nada se recalcula después de ver los resultados.

## Lo que desaparece o cambia de la v2.2.1 (tomado de `manuscript/ERRATAS_V3.md`, con la fila)

- **T** — C4 taquiónico: desaparece la exclusión (ΔlnL −12.7); sustituida por cota de amplitud del modo que crece en tiempo cósmico — filas 1–33 (33): 1 manuscript/PREPRINT.md:25; 2 manuscript/PREPRINT.md:68; 3 manuscript/PREPRINT.md:136; 4 manuscript/PREPRINT.md:212; 5 manuscript/PREPRINT.md:214; 6 manuscript/PREPRINT.md:217; 7 manuscript/PREPRINT.md:220; 8 manuscript/PREPRINT.md:222; 9 manuscript/PREPRINT.md:224; 10 manuscript/PREPRINT.md:226; 11 manuscript/PREPRINT.md:263; 12 manuscript/PREPRINT.md:357; 13 manuscript/PREPRINT.md:384; 14 manuscript/CLAIMS.md:20; 15 README.md:40; 16 manuscript/COVER_LETTER.md:9; 17 manuscript/COVER_LETTER.md:11; 18 manuscript/README.md:24; 19 papers/BOUND_CORRECTION_SIGMA_X.md:19; 20 papers/HONEST_HEADLINES.md:29; 21 papers/HONEST_HEADLINES.md:46; 22 papers/amplification-gap.md:74; 23 papers/amplification-gap.md:100; 24 papers/amplification-gap.md:136; 25 papers/covariance-desi-dr2-full.md:71; 26 papers/quantum-fluid-instabilities-desi-dr2.md:155; 27 papers/quantum-fluid-instabilities-desi-dr2.md:158; 28 papers/quantum-fluid-instabilities-desi-dr2.md:166; 29 papers/quantum-fluid-instabilities-desi-dr2.md:168; 30 papers/quantum-fluid-instabilities-desi-dr2.md:178; 31 theory/euclid-protocol-vacuum-relaxation.md:244; 32 theory/unimodular-gravity-vacuum-smoothness.md:50; 33 theory/unimodular-gravity-vacuum-smoothness.md:54
- **B** — cotas θ ≥ 1 con núcleo local (1.3–1.5 × 10⁻², 2.2–4.6 × 10⁻²): sustituidas por las del convenio B con respuesta completa — filas 34–56 (23): 34 manuscript/PREPRINT.md:23; 35 manuscript/PREPRINT.md:172; 36 manuscript/PREPRINT.md:173; 37 manuscript/PREPRINT.md:174; 38 manuscript/PREPRINT.md:176; 39 manuscript/PREPRINT.md:198; 40 manuscript/PREPRINT.md:202; 41 manuscript/PREPRINT.md:204; 42 manuscript/PREPRINT.md:255; 43 manuscript/PREPRINT.md:273; 44 manuscript/PREPRINT.md:323; 45 manuscript/PREPRINT.md:324; 46 manuscript/PREPRINT.md:383; 47 manuscript/PREPRINT.md:391; 48 manuscript/CLAIMS.md:18; 49 manuscript/CLAIMS.md:23; 50 manuscript/CLAIMS.md:24; 51 manuscript/CLAIMS.md:25; 52 manuscript/CLAIMS.md:47; 53 README.md:40; 54 manuscript/COVER_LETTER.md:11; 55 manuscript/README.md:23; 56 papers/resume.txt:58
- **S** — comparación DESI–Pantheon+ invertida: las SN acotan más que los BAO — filas 57–58 (2): 57 manuscript/PREPRINT.md:176; 58 manuscript/PREPRINT.md:391
- **L** — tabla del slip con 2.633 × 10⁻² como entrada: recalcular con el convenio B — filas 59–80 (22): 59 manuscript/PREPRINT.md:31; 60 manuscript/PREPRINT.md:307; 61 manuscript/PREPRINT.md:308; 62 manuscript/PREPRINT.md:309; 63 manuscript/PREPRINT.md:315; 64 manuscript/PREPRINT.md:388; 65 manuscript/CLAIMS.md:23; 66 manuscript/README.md:27; 67 papers/resume.txt:16; 68 papers/BOUND_CORRECTION_SIGMA_X.md:30; 69 papers/amplification-gap.md:87; 70 papers/anisotropic-slip-option0.md:25; 71 papers/anisotropic-slip-option0.md:41; 72 papers/anisotropic-slip-option0.md:118; 73 papers/anisotropic-slip-option0.md:149; 74 papers/anisotropic-slip-option0.md:150; 75 papers/anisotropic-slip-option0.md:151; 76 papers/anisotropic-slip-option0.md:157; 77 papers/anisotropic-slip-option0.md:158; 78 papers/anisotropic-slip-option0.md:160; 79 papers/anisotropic-slip-option0.md:166; 80 papers/anisotropic-slip-option0.md:183
- **K** — núcleo S: signo de la definición y fiducial — filas 81–88 (8): 81 manuscript/PREPRINT.md:97; 82 manuscript/PREPRINT.md:101; 83 papers/quantum-fluid-instabilities-desi-dr2.md:123; 84 papers/sensitivity_kernel_table.md:11; 85 papers/stochastic-dark-energy-desi-dr2.md:210; 86 papers/stochastic-dark-energy-desi-dr2.md:236; 87 papers/stochastic-dark-energy-desi-dr2.md:505; 88 papers/stochastic-dark-energy-desi-dr2.md:529
- **H** — apéndice H.5: cota sobre ζ V/(M_Pl V′) con 2.6 × 10⁻²: recalcular — filas 89–95 (7): 89 papers/stochastic-dark-energy-desi-dr2.md:869; 90 papers/stochastic-dark-energy-desi-dr2.md:873; 91 papers/stochastic-dark-energy-desi-dr2.md:879; 92 papers/stochastic-dark-energy-desi-dr2.md:887; 93 papers/stochastic-dark-energy-desi-dr2.md:891; 94 papers/stochastic-dark-energy-desi-dr2.md:899; 95 papers/stochastic-dark-energy-desi-dr2.md:907
- **M** — etiqueta '95 % CL' de la caída de 1.92: corregir a cobertura unilateral ≥ 97.5 % o pasar a −1.35 — filas 96–123 (28): 96 manuscript/PREPRINT.md:19; 97 manuscript/PREPRINT.md:166; 98 manuscript/PREPRINT.md:250; 99 manuscript/PREPRINT.md:383; 100 manuscript/CLAIMS.md:18; 101 README.md:34; 102 manuscript/COVER_LETTER.md:11; 103 manuscript/README.md:23; 104 papers/resume.txt:56; 105 papers/resume.txt:94; 106 papers/HONEST_HEADLINES.md:28; 107 papers/HONEST_HEADLINES.md:46; 108 papers/amplification-gap.md:37; 109 papers/amplification-gap.md:118; 110 papers/quantum-fluid-instabilities-desi-dr2.md:12; 111 papers/quantum-fluid-instabilities-desi-dr2.md:192; 112 papers/stochastic-dark-energy-desi-dr2.md:372; 113 papers/stochastic-dark-energy-desi-dr2.md:374; 114 papers/BOUND_CORRECTION_SIGMA_X.md:12; 115 theory/euclid-protocol-vacuum-relaxation.md:11; 116 theory/euclid-protocol-vacuum-relaxation.md:253; 117 theory/fundamental-vs-emergent-vacuum-relaxation.md:14; 118 theory/fundamental-vs-emergent-vacuum-relaxation.md:86; 119 theory/principle-of-vacuum-smoothness.md:23; 120 theory/unimodular-gravity-vacuum-smoothness.md:14; 121 theory/unimodular-gravity-vacuum-smoothness.md:40; 122 theory/unimodular-gravity-vacuum-smoothness.md:82; 123 theory/unimodular-gravity-vacuum-smoothness.md:303
- **Z** — referencia Zwane: completar volumen y artículo — filas 124–125 (2): 124 manuscript/PREPRINT.md:428; 125 papers/stochastic-dark-energy-desi-dr2.md:964
- **X** — remisiones fuera de git: notas del repo hermano con 2.5 × 10⁻² y test de correlación cruzada — filas 126–128 (3): 126 manuscript/PREPRINT.md:350; 127 README.md:15; 128 manuscript/PREPRINT.md:325

Total: 128 filas (= recuento por dos caminos de ERRATAS_V3.md).

## Citas que hay que verificar contra la fuente antes de escribir (sin verificar aún)

- DESI DR2 Results II (DESI Collaboration 2025, arXiv:2503.14738) y el paquete Zenodo 16644577 (figure6, 13 números gaussianos); DESI DR2 Lyα (arXiv:2503.14739).
- Pantheon+: Brout et al. 2022 (arXiv:2202.04077), Scolnic et al. 2022; valor publicado Ω_m = 0.334 ± 0.018 (abstract).
- DES-SN5YR: DES Collaboration 2024 (arXiv:2401.02929), Vincenzi et al. 2024, Sánchez et al. 2024; Ω_m = 0.352 ± 0.017; Zenodo 12720778 (CC-BY-4.0).
- Dovekie: Popovic et al. 2026 (arXiv:2511.07517); Ω_m = 0.330 ± 0.015; verosimilitud oficial con covarianza inversa.
- Union3 / UNITY1.5: Rubin et al. 2023 (arXiv:2311.12098), Tabla 5 (Ω_m = 0.356 +0.028/−0.026, χ² = 24.0); repositorio MIT.
- Planck 2018: resultados VI (parámetros), TT ℓ < 30 y lensing bandpowers (PLA); valores de 100 θ* y ω_m con sus errores usados como compresión.
- Everpresent Λ: Ahmed, Dodelson, Greene & Sorkin 2004 (astro-ph/0209274); Das et al. 2023 I y II (arXiv:2304.03819, 2307.13743); Zwane, Afshordi & Sorkin 2018, Class. Quantum Grav. 35, 194002 (ya comprobada en ERRATAS_V3, categoría Z).
- Vikman 2005 (cruce de w = −1 en campos escalares, astro-ph/0407107): enunciado exacto del teorema de no cruce.
- Constante cosmológica que cambia de signo: Akarsu et al. (ΛsCDM; arXiv:1912.08751, 2108.09239 y posteriores): qué modelo y qué z de cambio de signo.
- Reconstrucciones no paramétricas de la energía oscura con DESI DR2 (DESI DR2 extended dark-energy analysis, arXiv:2503.14743; reconstrucciones independientes publicadas en 2025–2026): qué dicen de cruces de w = −1 y de ρ_DE < 0.
- Matheron (regla de muestreo condicional) y la marginal gaussiana: referencia estándar de geoestadística para citar la regla usada en MARGINAL_GUIADA.
- CAMB (Lewis, Challinor & Lasenby 2000) y DarkEnergyPPF (Fang, Hu & Lewis 2008): versión usada y límites de la parametrización PPF con ρ_DE ≤ 0.

## Recuento de afirmaciones por nivel, por dos caminos

| nivel | camino 1: etiquetas `(nivel k, …)` en el texto | camino 2: suma de las listas por sección |
|---|---|---|
| 1 (observado) | 27 | 27 |
| 2 (deducido) | 26 | 26 |
| 3 (supuesto) | 8 | 8 |
| 4 (no se puede saber) | 8 | 8 |
| total | 69 (líneas `- **A<n>**`: 69) | 69 (por sección: 0 → 6, 1 → 3, 2 → 5, 3 → 7, 4 → 8, 5 → 18, 6 → 6, 7 → 6, 8 → 4, 9 → 2, 10 → 4) |

