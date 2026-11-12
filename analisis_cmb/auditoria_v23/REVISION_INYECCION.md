# REVISION_INYECCION: revisión independiente de `analisis_cmb/auditoria_v23/inyeccion.py` e `INYECCION.md`

Revisión independiente hecha el 2026-11-12 en una sesión nueva de Claude Code sin contexto del trabajo.

Revisión de solo lectura (2026-11-12), sin participación previa en el test. Ficheros leídos: `inyeccion.py`,
`INYECCION.md`, `inyeccion.{json,txt}`, `muestras_sn.py`, `exacto_sn.py`, `familias.py` (sample, cell_seed, kernel_matrix),
`marginal_guiada.py` (kernel, lnZ_lin), `exacto.py` (constantes y malla), `exacto_sn.json`, `marginal_guiada.json`,
`muestras_sn_lcdm.json`, MUESTRAS_SN.md y MARGINAL_GUIADA.md (solo para cruzar los números de referencia).
Las líneas citadas son de `inyeccion.py` salvo que se diga otra cosa. Nada se ha ejecutado salvo lecturas de los json ya
escritos (recuento de coberturas, cocientes emparejados y percentiles, que coinciden con la tabla de INYECCION.md).

**Veredicto global: no encuentro ningún error que cambie una conclusión de la tabla.** Los cinco puntos del encargo salen
"correcto"; hay tres asuntos "dudoso" que son de atribución o de interpretación, no de cálculo, y se listan al final.

---

## 1. Generación de los universos falsos — **correcto**

**1a. Ruido con el Cholesky de la covarianza completa.** `make_full` (l. 88–92) toma `C` de `MS.make_linear`
(`muestras_sn.py` l. 165–168: bloque BAO 13×13 = `X.COV13`, θ* y ω_m diagonales con `SIG_TH²`, `SIG_WM²`, bloque SN =
`inv(sn["Ci"])`, es decir la covarianza STAT+SYS de Pantheon+ con el corte z_HD > 0.01 y sin calibradores) y guarda
`Lc = np.linalg.cholesky(C)` (l. 92). El ruido es `noise = Lc @ Z` con `Z ~ N(0,1)` de forma (1595, 200) (l. 348–349), así
que `noise ~ N(0, C)`. Correcto. La doble inversión del bloque SN (`inv(inv(Cs))`) es numéricamente inocua: el χ²_ΛCDM lineal
y el exacto coinciden a 13 cifras (`real_check`), lo que confirma que `C` y el `L` de `exacto_sn` describen la misma matriz.

**1b. Misma semilla entre escenarios.** `Z` se genera una sola vez (l. 348) y el mismo `noise` entra en A (l. 360), en las
tres B (l. 372), en C (l. 373) y en D a través de `compress(binned, data["A"])` (l. 374). Correcto: el universo u tiene
idéntico vector de ruido en los seis escenarios, lo que justifica los cocientes emparejados C/A y D/A. Los tres generadores
usan semillas independientes y fijas: caminos de análisis `[SEED_I, 1]` (l. 314), ruido `[SEED_I, 2]` (l. 347), caminos
inyectados `[SEED_I, 3]` (l. 350).

**1c. Temblor de B con la respuesta exacta.** `Xinj = S.ou_paths(G0, 1.0, 1.0, N_UNI, rng_p)` son caminos OU unitarios
(σ = 1, θ = 1, estacionarios: `exacto_sn.py` l. 87–95) y `dXinj = Xinj − Xinj[:, [0]]` (l. 352) aplica el convenio B.
Para cada σ_true y universo, `MS.exact_model(sn, h0, wm0, M0, (st · dXinj[u])[None, :])` (l. 366) evalúa el modelo
**no lineal completo** (BAO, θ*, ω_m, μ_SN) en P0 con la historia perturbada, y se inyecta `out − B0` (l. 370). Es la
respuesta exacta, no la lineal `A`. El descarte por H² ≤ 0 se cuenta (`inj_bad`, l. 367–371) y fue 0 en los tres σ_true.
El mismo camino unitario por universo se escala en los tres σ_true, así que las B también están emparejadas entre sí.

**1d. Escenario D.** La compresión GLS `T = (EᵀWE)⁻¹EᵀW` con `W = sn["Ci"]` y covarianza `Cb = (EᵀWE)⁻¹` (l. 106–109)
se aplica de forma consistente a datos, B0, Bm y A (l. 114–116, 126). El ruido comprimido `T·noise_SN` tiene covarianza
`T C_SN Tᵀ = (EᵀWE)⁻¹ = Cb`, luego el setting binned es internamente coherente. 23 bordes → 22 tramos, todos con SN.

## 2. Regla del 5 % lineal frente a la de `exacto_sn.py` — **correcto, con las diferencias listadas**

Con `R = D − B0`, `W⊥ = W − W Bm (BmᵀWBm)⁻¹ BmᵀW` (`wperp`, l. 82–85) y `Y_i = A x_i` (l. 139):
`a = RᵀW⊥R` (l. 171), `b_i = Y_iᵀW⊥R` (l. 172), `c_i = Y_iᵀW⊥Y_i` (l. 141). Entonces
χ²_i(σ) = min_dp ‖R − σY_i − Bm dp‖²_W = a − 2σ b_i + σ² c_i, y χ²_i(σ) ≤ χ²_ΛCDM = a ⇔ σ ≤ s_i = 2b_i/c_i (c_i ≥ 0 por ser
W⊥ semidefinida positiva). Luego P(σ) = frac(s_i ≥ σ) es no creciente y P(σ) = 0.05 ⇔ σ = percentil 95 de s_i
(l. 173–174). El recorte `max(·, 0)` no afecta al percentil 95 mientras menos del 95 % de los s_i sean ≤ 0, que es siempre el
caso aquí. Esto es exactamente el criterio de `exacto_sn.py` (P(χ²_min ≤ χ²_obs) = 5 % con χ²_obs = χ²_min de ΛCDM,
l. 391 y 258–266 de ese fichero) en su límite lineal.

Dónde difiere (ninguna diferencia es un error; todas están en la dirección de "la exacta es la referencia"):

- **Descarte por H² ≤ 0**: la exacta quita del numerador y del denominador los caminos con H² ≤ 0 (`exacto_sn.py` l. 230–235,
  `_exact_task` l. 257–261); la lineal cuenta los 2000. INYECCION.md l. 47 lo declara. Solo pesa con cotas ≳ 1 (B 0.3 u0/u2,
  45–54 % descartados), donde el propio documento dice que la lineal ya no vale.
- **Caminos por celda**: `exacto_sn.py` genera 2000 caminos distintos por (σ, θ) (semilla por celda); la lineal escala los
  mismos 2000 unitarios para todo σ (y la comprobación exacta de aquí hace lo mismo, l. 254). Es equivalente en distribución;
  solo cambia la estructura del ruido Monte Carlo (la P(σ) lineal es monótona por construcción).
- **Interpolación**: la exacta interpola P en ln σ entre nodos (0.1, 0.15, 0.2 …); la lineal da el percentil sin malla. Parte
  del −8 % de los datos reales (0.156 frente a 0.169) es interpolación, no no-linealidad: con la malla de `exacto_sn.json`
  (0.1, 0.2) la misma celda da 0.173.
- **Percentil**: `np.quantile(…, 0.95)` con interpolación lineal entre estadísticos de orden frente al "primer σ con P ≤ 0.05"
  de `crossing`: difieren en menos de un estadístico de orden (1/2000). Despreciable.
- **Punto de descarte en la comprobación exacta**: `_exact_task` descarta H² ≤ 0 en el mejor ajuste ΛCDM **del universo**
  (l. 253–257), no en el P0 fijo de las 15 filas como `exacto_sn.py` (l. 228–231). Es lo razonable para universos falsos.
- **Ajustes no lineales con `max_nfev=200`** (l. 243): si alguno no converge, χ² sube, P baja y la exacta sale más apretada.
  No lo he verificado (no he ejecutado nada); coincide con la dirección de los cocientes 0.885–0.90 en B 0.3, pero esos casos
  ya están fuera del rango de validez y no cambian la tabla.

## 3. Marginal lineal: normalización entre σ y misma curva para cota y Δln L máximo — **correcto**

`MG.lnZ_lin` (`marginal_guiada.py` l. 155–161) devuelve −½[rᵀS⁻¹r − bᵀF⁻¹b + ln det S + ln det F] con S = C + σ²·AKAᵀ,
F = BmᵀS⁻¹Bm, b = BmᵀS⁻¹r. Es exactamente ln ∫ N(r | Bm·dp, S) d(dp) con prior plano en dp = (δh, δω_m, δM) salvo la constante
−(n−3)/2·ln 2π, que no depende de σ. Los dos términos que varían con σ (ln det S y ln det F) están incluidos, así que la curva
está bien normalizada entre σ. `lnz_profile` (l. 153–161) reproduce la misma expresión con factorizaciones compartidas por
universo (`self.fac`, l. 143–148: `ld = 2 Σ ln diag(L) + slogdet(F)`; con `lower=True` la diagonal de `Sf[0]` es la de L) y
resta el nodo σ = 0 (l. 161; `GRID[0] = 0.0`, l. 65). Comprobaciones internas que lo respaldan: (i) `lnz_single` → `MG.lnZ_lin`
y el brentq sobre esa función dan la misma cota real que el spline (0.27712 / 0.27712, `real_check`); (ii) la marginal Monte
Carlo sobre caminos (l. 183–188, media de exp(σ b_i − σ² c_i/2), que es la misma integral con dp perfilado por camino —válido
porque el hessiano en dp no depende de x) coincide con la analítica donde la ESS es alta (A 0.111/0.109, C 0.164/0.165,
D 0.112/0.111); un fallo en ln det S o ln det F las separaría.

El kernel `MG.kernel(1.0) = exp(−|u_i − u_j|)` sobre `G.U = S.G0.U` (`marginal_guiada.py` l. 44, 151–152) es la
autocovarianza unitaria del OU θ = 1 que genera `S.ou_paths(G0, 1.0, 1.0, …)`, así que σ²·AKAᵀ es la covarianza de σ·A·x. Misma
malla que `A` de `make_linear`. Prior plano en (h, ω_m, M): sí, `Bm` son las derivadas respecto a esas tres (`muestras_sn.py`
l. 172–177).

**Misma curva.** `prof = self.lnz_profile(R)` (l. 175); `dlnL_max = prof.max(axis=0)` y `excl0 = prof.max(axis=0) > DROP`
(l. 177–178); `marg = drop_bound(prof[:, j])` (l. 181) con umbral `prof[k0] − DROP` desde el máximo global en los nodos
(l. 205–207) y cruce por spline cúbico + brentq a la derecha del máximo (l. 208–213). Tal como dice la nota de INYECCION.md
(l. 130–137). He verificado en `per_universe` que `excl0 == (dlnL_max > 1.92)` en los 1200 universos. El máximo se busca en los
91 nodos y no en el spline: con una malla geométrica tan densa la diferencia es despreciable. Único matiz: si el máximo está en
σ = 0 (`k0 = 0`) el brentq puede empezar en `GRID[0] = 0` (l. 213), a diferencia de `muestras_sn._linear`, que devuelve −1 en
ese caso (l. 240 de ese fichero); aquí es más correcto, no menos.

## 4. Cobertura y colocación de la tensión en C — **correcto**

**Cobertura.** `coverage_5pct = mean(b5 ≥ σ_true)`, `coverage_marg = mean(~finite(marg) | marg ≥ σ_true)` (l. 303–305):
es "σ_true ≤ cota", con "sin cruce" contado como cota infinita (cubre). `marg_no_crossing` = 0 en los seis escenarios, así que
la segunda cláusula no actúa. Recomputado desde `per_universe`: 0.925/0.955, 0.865/0.965, 0.850/0.975, igual que la tabla.
Las cotas del 5 % > 1 (3.5 % y 15.5 % en B 0.2/0.3) cuentan como "cubre", que es lo lógico y lo que dice INYECCION.md l. 91–92.

**Escenario C.** `tension_noise_free` (l. 278–287) evalúa el modelo exacto con ω_m' = (Ω_m + δ)·h0² a h0 y M0 fijos, y después
`Bc[:15] = full["B0"][:15]` (l. 284): BAO, θ* y ω_m quedan en el fiducial y solo cambian las 1580 filas SN. `data["C"] = Bc +
noise` (l. 373). Sí, la tensión está solo en las SN. El estadístico Δχ² = χ²_joint(h, ω_m, M) − χ²_15(h, ω_m) − χ²_SN(ω_m, M | h)
(l. 194–197; `Bm[:15, :2]` y `Bm[15:, 1:]`) es una forma cuadrática PSD de rango 1 (el modelo separado, 4 parámetros, anida al
conjunto, 3; la columna h de Bm restringida a las SN está en el span de las de ω_m y M), luego sin tensión es χ²₁ con media
tr(QC) = (1595−3) − (15−2) − (1580−2) = 1, y con desplazamiento sin ruido λ su media es λ + 1. Por eso `f = λ(δ) − 7.7`
(l. 354) y la media prevista es 8.7. Medias obtenidas: A 0.93 ± 0.09, D 0.95 ± 0.09, C 8.09 ± 0.38 (σ teórica de la
no-central χ²₁(λ = 7.7) es 5.73, /√200 = 0.40): C está a 1.5σ del 8.7, compatible. El 8.7 de referencia es el Δχ² de DES-SN5YR
de `muestras_sn_lcdm.json` (8.73), con la misma definición (MUESTRAS_SN.md l. 202–210). Ω_m 0.3037 → 0.3549 frente al 0.352 de
DES: coherente.

## 5. Errores que cambien una conclusión de la tabla — **no encontrados**

Recomputado desde `inyeccion.json` (sin volver a correr nada): coberturas, fracciones σ_X = 0 excluido, medias y medianas de
tensión, cocientes emparejados C/A (1.80 [1.12–2.35] y 1.58 [1.07–1.94]) y D/A (1.05 [0.96–1.12] y 1.04 [0.95–1.11]), y los
percentiles de los datos reales dentro de A (0.156 → 97.0 %, 0.169 → 99.5 %, 0.277 → 97.5 %, +1.14 → 97.0 %, 2.66 → 91.5 %;
p95 de A: 0.153, 0.229, +0.82, 4.28). Todo coincide con INYECCION.md l. 83–98. Las cuatro predicciones están juzgadas con los
números correctos. El pipeline D es consistente (misma T para datos, B0, Bm, A; Cb = (EᵀWE)⁻¹) y la explicación de por qué no
pierde información (la proyección GLS conserva casi toda la Fisher de un modelo suave) es correcta: el escenario no reproduce la
pérdida de información de Union3, como dice el documento.

### Dudoso (no cambia la tabla)

- **Atribución del 0.169 (INYECCION.md l. 36 y 48, `REAL` l. 66).** `exacto_sn.json` da 0.173 para θ = 1 BAO+CMB+SN
  (interpolando entre 0.1 y 0.2); el 0.169 (0.164–0.173) es la reproducción con Pantheon+ de `muestras_sn.json`, con nodos
  0.1/0.15/0.2, y MARGINAL_GUIADA.md da 0.170. Las tres son la misma celda con distinta malla de σ. Con 0.173 el cociente de la
  comprobación previa sería 0.901, y `ok_within_10pct` (l. 334) pasaría por una milésima. No afecta a nada posterior, pero
  conviene citar MUESTRAS_SN y no EXACTO_SN para el 0.169, y saber que la condición de parada estaba al borde.
- **"Cobertura" de la regla del 5 %.** La regla del 5 % no es un intervalo de confianza y no tiene cobertura nominal; el
  "se queda corto" del documento (l. 123) es una comparación empírica, no la violación de un objetivo. La caída de 1.92 desde el
  máximo tiene cobertura nominal unilateral ≈ 97.5 % en el límite gaussiano, y los 95.5–97.5 % observados están en el blanco.
- **Comprobación exacta fuera del 10 % en 3 de 30** (B 0.3 u0/u2 con cota > 1; C u1 al 10.2 %). El documento lo deja
  constar (l. 102–105) y la condición de parada del prompt no la incluía; correcto como está, pero la frase "la versión exacta
  reproduce la lineal a menos de un 10 %" de INYECCION.md l. 37–38 (sección Montaje) queda matizada por la sección posterior.

### Menores (sin efecto)

- `cross_5pct` usa `P[k−1] ≥ 0.05 > P[k]` (l. 269) y `exacto_sn.crossing` usa `p[k] > 0.05 ≥ p[k+1]` (l. 263 de ese fichero):
  distinta convención en la igualdad; irrelevante con P continua.
- `ratio_vs_A` en la tabla (l. 394–397) es cociente de medianas entre escenarios, no mediana de cocientes emparejados; el
  documento da ambos (l. 90) y son compatibles (1.71 frente a 1.80, 1.50 frente a 1.58).
- El máximo de `prof` se busca en los nodos, no en el spline (l. 205); con 91 nodos geométricos hasta 1.5 no cambia nada.

---

Esta revisión no corrige nada: solo informa. No se ha modificado ningún fichero del repositorio ni se ha hecho commit.
