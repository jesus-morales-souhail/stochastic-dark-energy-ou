> Rutas a fecha del commit d71c4d2, antes de la reorganización del 2026-10-17.

# Bloque B1: cota = raíz del perfil (2.633e-2), constantes desde el json, resultados regenerados (2026-10-17)

**Decisión del autor:** σ_X < 2.6e-2, raíz refinada del perfil (anotada en AUDITORIA.md, tabla 2a última fila y principio del bloque B). No se han tocado papers/, manuscript/ ni README.

- Copia de todo lo regenerado, tal como estaba en `db8c938`: `analisis_cmb/bloque_B/antes/`.
- Salidas de consola: `analisis_cmb/bloque_B/logs/`.

## 1. Perfil (`scripts/profile_sigma_x_desi.py`)

Se añade la raíz de ΔlnL(rel. al máximo) = −1.92, con brentq sobre la función del perfil (máximo en θ) entre los dos nodos que la encierran, 2.506e-2 y 3.155e-2. Es el mismo procedimiento que `analisis_cmb/everpresent/profile_desi_everpresent.py` (`upper95_refined`). La malla `geomspace(1e-6, 5e-2, 48)` y los datos no cambian.

| magnitud | fichero, clave | antes | después |
|---|---|---|---|
| cota 95 % | profile_sigma_x.json `sigma_X_95CL_upper_profile` | 0.025063155506929828 (nodo) | **0.0263282756994459** (raíz) |
| nodo de malla | `sigma_X_95CL_upper_profile_grid` | — | 0.025063155506929828 |
| θ en la cota | `theta_best_at_upper` | — | 0.001 (borde inferior de la malla de θ) |
| claim redondeado | `sigma_X_95CL_upper_profile_claim` | 2.5e-2 | 2.6e-2 |
| criterio | `criterion` | "Delta logL >= -1.92 …" | "Delta logL = -1.92 …; root by brentq between grid nodes" |
| 48 filas del perfil, `ll_lcdm` | `profile`, `ll_lcdm` | — | idénticas |
| txt | profile_sigma_x.txt l.5–7 | "σ_X ≤ 2.506e-02"; "claim σ_X < 2.5e-2" | "σ_X < 2.6328e-02 (θ_best = 0.001)"; "last grid node 2.5063e-02"; "claim σ_X < 2.6e-2" |
| figura | results/profile_sigma_x/profile_sigma_x.png, figures/profile_sigma_x.png | axvline 2.5e-2 con etiqueta "claim 1.5×10⁻⁴" | nodo 2.506e-2 y raíz 2.633e-2 rotulados |

**Comprobación:** raíz = 0.0263282756994459, igual a `profile_desi_everpresent.json` → `T0_ou_theta_profiled.upper95_refined` = 0.0263282756994459. Diferencia relativa: 0.0, mejor que 1e-5.

## 2. Constantes escritas a mano → `desi_sigma_x_bound()`

Función común en `scripts/desi_bound.py` (nuevo). Lee `sigma_X_95CL_upper_profile` de `results/profile_sigma_x/profile_sigma_x.json`. Todos los scripts la importan sin problemas, así que ninguno lleva el valor fijo 2.633e-2.

| fichero:línea (antes) | nombre | valor antes | ahora |
|---|---|---|---|
| scripts/gap_two_targets.py:23 | `SIGMA_DESI` | 2.5e-2 | `desi_sigma_x_bound()` |
| scripts/gap_two_targets.py:7, 56, 57 | docstring "G_DESI = 2.5e-2/sigma_0 ~ 10^57"; dos prints | 2.5e-2 | "~10^59" y f-strings con la cota |
| scripts/amplifier_audit.py:20 | `SIGMA_X_DESI` | 2.5e-2 | `desi_sigma_x_bound()` |
| scripts/slip_bridge.py:45 | `SIGMA_X_95` | 2.5e-2 | `desi_sigma_x_bound()` |
| scripts/slip_bridge.py:141 | print "σ_X < 2.5×10⁻²" | 2.5e-2 | f-string |
| scripts/amplification/route1_local_causal_set_seed.py:42 | `TARGETS["desi_1p5e-4"]` | 2.5e-2 (etiqueta 1.5e-4) | `TARGETS["desi_2p6e-2"] = desi_sigma_x_bound()` |
| scripts/amplification/route1_local_causal_set_seed.py:16 | docstring | 2.5e-2 | "DESI 95% CL bound" |
| scripts/desqueezing/cosmological_mapping.py:38 | `SIGMA_X_LIMIT_95` | 2.5e-2 | `desi_sigma_x_bound()` |
| scripts/desqueezing/cosmological_mapping.py:209, 267, 269 | etiqueta "DESI_limit_2.5e-2"; dos prints | 2.5e-2 | f-strings con la cota |
| scripts/desqueezing/cosmological_mapping.py:39, 137 | `SIGMA_X_MLE_FLOOR` | 5e-5 | sin cambio de valor; comentario "old run, diagonal cov" |
| scripts/desqueezing/euclid_protocol_forecasts.py:21 | `SIGMA_DESI` | 2.5e-2 | `desi_sigma_x_bound()` |
| scripts/desqueezing/euclid_protocol_forecasts.py:36 | entrada de `A0s` | 2.5e-2 | `SIGMA_DESI`; comentario: 0.018 = calibración DR1 obsoleta |
| scripts/desqueezing/euclid_protocol_forecasts.py:66 | leyenda "σ=1.5×10⁻⁴" | 1.5e-4 | "σ=0.0263 (DESI 95%)" |
| scripts/desqueezing/sdiff_discrimination.py:21 | `SIGMA_X_DESI` | 2.5e-2 | `desi_sigma_x_bound()` |
| scripts/desqueezing/sdiff_discrimination.py:67 | leyenda "DESI 2.5e-2" | 2.5e-2 | f-string |
| scripts/desqueezing/desqueezing_relax_time.py:311 | texto "σ_X < 2.5e-2" | 2.5e-2 | f-string con `desi_sigma_x_bound()` |
| scripts/profile_sigma_x_desi.py:109, 128, 146 | claim; axvline "1.5×10⁻⁴"; línea del txt | 2.5e-2 / 1.5e-4 | sección 1 |
| scripts/desi_dr2_real_bao_test.py:99 | clave `sigma_X_95_profile_upper` | 0.02175 (aproximada) | clave `sigma_X_95_profile_upper_approx`, valor igual |
| analisis_cmb/realizaciones.py:51, 53 | `REFS["CPL_repo"]`, `REFS["CPL_joint"]` | (−0.98990, −0.01590), (−1.00165, +0.01725) | leídas de `results/eos_cpl_desi_dr2/eos_cpl_summary.json` (`pure_CPL`) y `results/joint_w0wa_sigma/joint_w0wa_sigma.json` (modelo `cpl`) |
| analisis_cmb/verificacion_numeros.py:22–23 | θ en la cota | buscaba el nodo en `profile` | `theta_best_at_upper` |

## 3. Ficheros de resultados: antes / después

**Orden de ejecución:**
1. gap_two_targets, amplifier_audit y slip_bridge. Solo imprimen: `logs/{gap_two_targets,amplifier_audit,slip_bridge}.log`.
2. route1.
3. cosmological_mapping, euclid_protocol_forecasts, sdiff_discrimination y desi_dr2_real_bao_test.
4. desqueezing_relax_time (punto 7b).

**Solo route1, no `run_all_routes.py`:** las rutas 2 y 3 no usan la cota. route2 barre σ ∈ {1e-6…1e-3}; route3 compara con umbrales fijos 1e-6…1e-3. Sus csv no cambian; el recuento BAO-safe de VERDICT se recalcula sobre el csv de route3.

| fichero | línea / clave | antes | después |
|---|---|---|---|
| results/amplification_routes/route1_N_eff_required.csv | l.3 | `desi_1p5e-4,0.025,1599.99` | `desi_2p6e-2,0.0263282756994459,1442.63` |
| results/amplification_routes/route1_N_eff_grid.csv | — | — | sin cambios (no depende de la cota) |
| results/amplification_routes/VERDICT.md (a mano) | l.22 | 2.5×10⁻², N_eff 1.6×10³ | 2.633×10⁻², **1.44×10³** (route1_N_eff_required.csv l.3) |
| ídem | l.26 | N_BH mayor en 10¹¹²–10¹¹⁴ | **10¹¹²–10¹¹⁹**: 10^112.0 (l.2), 10^114.0 (l.4), 10^118.8 (l.3) |
| ídem | l.55–56 | BAO-safe 216/288; unsafe 72/288, σ ≳ 1.2×10⁻⁴ | **288/288**; **0/288**. p95 máximo 2.08×10⁻³ en route3_avalanche_scan.csv l.289 (a = 10, b = 0.1, σ = 10⁻³) |
| ídem | l.58 | "σ ≲ 4×10⁻⁵" | "the whole grid (σ ≤ 10⁻³) is BAO-safe", con nota del valor anterior |
| results/cosmological_mapping/table2_A0_required.csv | l.2, 4, 6, 8, 10, 12, 14, 16 (filas DESI) | target `DESI_limit_1.5e-4`, 0.00015; A0 = 1.50e-4, 1.51e-4, 4.63e-4, 0.663, 2.88e-4, 1.91e-2, 1.65e-4, 3.00e-4 | target `DESI_limit_2.633e-02`, 0.026328; A0 = 2.635e-2, 2.651e-2, 8.134e-2, 116.4, 5.051e-2, 3.358, 2.890e-2, 5.266e-2 |
| results/cosmological_mapping/mapping_summary.txt | l.11 y formato | resumen corto escrito a mano: "sigma_X < 2.5e-2", "MLE floor … 5e-5" | salida completa del script: "sigma_X < 2.6e-02", "MLE floor … 5.0e-05 (old run, diagonal cov)" (+82 líneas de la salida real del script) |
| notes/mapping_tables.md | Tabla 2, filas DESI (l.22–36) | "DESI_limit_2.5e-2", A0 = 1.501e-4 … 3.000e-4 (calculados con 1.5e-4) | "DESI_limit_2.633e-02", A0 = 2.635e-2 … 5.266e-2 |
| notes/mapping_tables.md | cabecera y LaTeX | prosa pulida a mano, `$…$` | salida literal del script, `\(…\)`, y una línea "**Seed A₀=0.01.**" que imprime el script |
| results/euclid_protocol/euclid_forecast_grid.csv | l.68–78 | 11 filas con A0 = 0.00015 | 11 filas con A0 = 0.0263282756994459, por ejemplo θ = 0.001: σ_DESI_path 1.4986e-4 → 2.6304e-2 |
| results/euclid_protocol/euclid_A0_theta_plane.png, figures/ | figura | contorno 2.5e-2, leyenda "σ=1.5×10⁻⁴" | contorno y leyenda con 0.0263 |
| results/sdiff_discrimination/sdiff_discrimination_grid.csv | — | — | sin cambios (no depende de la cota) |
| results/sdiff_discrimination/*.png, figures/sdiff_*.png | líneas DESI | 2.5e-2 | 2.633e-2 |
| results/sdiff_discrimination/summary.txt (a mano) | l.5 | "DESI lim = 1.5e-04" | "DESI lim = 2.633e-02" (profile json) |
| ídem | l.16–20 | A0_min sin fuente | sin cambio de valor; cada uno con su fila: csv l.27 (θ = 0.001), 30 (0.1), 32 (0.693: 1.918e-5 frente a 1.919e-5 para θ = ln 2), 34 (1.2), 36 (5.0); A0_min = 10⁻¹⁰/sigma_res_DESI_path(amp_1e-5) |
| ídem | l.25–26, 29 | "[1e-5, 2.5e-2]", "< 2.5e-2", "A0 > 2.5e-2" | 2.633e-2 |
| results/desi_dr2_real_bao/desi_dr2_real_bao.json | l.39 | `sigma_X_95_profile_upper`: 0.021753 | `sigma_X_95_profile_upper_approx`: 0.021753 |
| results/ou_bao_desi_dr2_run.txt | l.1 | — | cabecera "cov diagonal; superseded por results/profile_sigma_x/ (2026-10-17)"; el resto no se toca |
| results/desqueezing_relax_time/desqueezing_relax_time_resumen.txt | todo | resumen en inglés escrito a mano ("t_1/2 = ln(2)/gamma … n_th does not change t_1/2") | salida completa del script; l.17 "σ_X < 2.633e-02" |
| results/desqueezing_relax_time/desqueezing_relax_time.png; figures/desqueezing_relax_time.png | figura | — | regenerada; la copia en figures/ se actualiza a mano (el script solo escribe en results/; antes eran idénticas) |
| results/desqueezing_relax_time/t_relax_tabla.csv | — | — | sin cambios |
| analisis_cmb/verificacion_numeros.md | todas las filas con σ_X | σ_X = 2.506e-2 | σ_X = 2.633e-2: r_DESI 68.33; N_eff 1.44e3; slip 3.39e-2 / 1.43e-2 / 7.33e-3; θ_min 2.02e-3; H.5 1.03e-4; λ ≳ 9.75e3 |

Observaciones sobre las salidas regeneradas:
- **mapping_summary.txt, notes/mapping_tables.md y desqueezing_relax_time_resumen.txt.** Las versiones versionadas no eran la salida del script, sino resúmenes reescritos a mano. Ahora son la salida literal.
  - Se pierde la prosa pulida de notes/mapping_tables.md.
  - La salida de desqueezing_relax_time incluye rutas absolutas de esta máquina (`/home/ashpokemon/...`).
  - Imprime además una "Teoría Markoviana" t½ = ln2/[γ(1+2n_th)] y un "Ratio sim/teoría 1.47 ± 0.58" que su propia tabla contradice: t½ = 0.28883 para todo n_th a γ = 2.4, y t½_th(n_th = 0) = 0.28881. El resumen anterior, escrito a mano, decía lo correcto: t½ no depende de n_th.
  - No he tocado ni el script ni su salida.

## 4. Corrección de ROOT

| script | antes | después | desde cuándo |
|---|---|---|---|
| scripts/desqueezing/cosmological_mapping.py l.21–23 | `ROOT = dirname(dirname(abspath(__file__)))` → `scripts/`; `DOCS = ROOT/docs` | `ROOT` sube tres niveles (raíz del repo); `DOCS = ROOT/notes` | **33eb67b** ("Restructure repository for review…"): el script entró ya en `scripts/desqueezing/` (como `cosmological_mapping_from_repo.py`) con `ROOT` de dos niveles. Se renombró en 99e4ce6 sin arreglarlo |
| scripts/desqueezing/desqueezing_relax_time.py l.38 | ídem → `scripts/results/desqueezing_relax_time` | `ROOT` sube tres niveles | **33eb67b**, mismo caso |

**Consecuencia del fallo:** desde 33eb67b, ejecutar cualquiera de los dos scripts escribía en `scripts/results/…`. `cosmological_mapping` además fallaba al abrir `scripts/docs/mapping_tables.md`. Los ficheros versionados en `results/` no salían de estas versiones de los scripts.

**En B1:** la primera ejecución de cosmological_mapping creó `scripts/results/cosmological_mapping/` y falló. Borré ese directorio, corregí las rutas y la volví a ejecutar. Solo cambian las rutas, ninguna fórmula ni dato.

## 5. CMB a 2.633e-2 frente a 2.506e-2

Configuración idéntica a `sigma_2.506e-02`: 200 realizaciones, semilla 20261003, takahashi, H0 resuelto por θ_*, Δ = 0.05, θ ∈ {0.001, 1}. Salida: `analisis_cmb/resultados/takahashi/sigma_2.633e-02/`; la carpeta `sigma_2.506e-02` no se ha tocado. 0 fallos (ρ ≤ 0, CAMB, NaN, H0).

| θ | magnitud | 2.506e-2 (nodo de rejilla, versión anterior) | 2.633e-2 | cociente |
|---|---|---|---|---|
| 0.001 | Δχ² TT p95 | 2.518e-3 | 2.645e-3 | 1.051 |
| 0.001 | Δχ² φφ p95 | 2.193e-3 | 2.304e-3 | 1.051 |
| 0.001 | H0 p2.5 / p16 / p50 / p84 / p97.5 | 67.32 / 67.34 / 67.36 / 67.38 / 67.40 | 67.32 / 67.34 / 67.36 / 67.38 / 67.40 | — |
| 0.001 | fuera de Planck ±2σ | 0 | 0 | — |
| 1 | Δχ² TT p95 | 0.06364 | 0.06682 | 1.050 |
| 1 | Δχ² φφ p95 | 0.05137 | 0.05411 | 1.053 |
| 1 | H0 p2.5 / p16 / p50 / p84 / p97.5 | 66.55 / 66.94 / 67.36 / 67.85 / 68.28 | 66.51 / 66.92 / 67.36 / 67.88 / 68.33 | — |
| 1 | fuera de Planck ±2σ | 0.025 | 0.035 | — |

- Δχ² > 4: 0 en todo.
- 2.633/2.506 = 1.051, coherente con el exponente 1 medido en la malla.
- Fuente: `resumen.json` de cada carpeta.

**Curvas de referencia (punto 7c).**
- `realizaciones.py` no tiene modo para recalcular solo las referencias, así que **se repitió el punto 2.633e-2 completo** con la misma configuración (`logs/realizaciones_sigma_2.633e-02_refs_json.log`).
- Los observables OU de las dos ejecuciones son **idénticos bit a bit** (`observables_*.npz` y `dchi2_*.txt`); solo cambian las curvas CPL.

| curva | (w0, wa) antes → después | H0 resuelto | Δχ² TT | Δχ² φφ |
|---|---|---|---|---|
| CPL_repo | (−0.98990, −0.01590) → (−1.00095, +0.01370) | 67.189 → 67.278 | 0.0129 → 0.00798 | 0.00877 → 0.0149 |
| CPL_joint | (−1.00165, +0.01725) → (−1.00095, +0.01370) | 67.270 → 67.278 | 0.00892 → 0.00798 | 0.0176 → 0.0149 |

p95 de OU frente a CPL_repo, a θ = 1: TT 0.0540 → 0.0588, φφ 0.0453 → 0.0392 (`tabla.md` de la carpeta).

**RESULTADO.md:**
- §2: filas de 2.633e-2 para los dos θ; las de 2.506e-2 quedan marcadas "(nodo de rejilla, versión anterior)".
- §3: línea del ×1.05.

## 6. Salida del punto 6

```
git grep -n -E "1\.5e-4|1\.5e-04|1\.5×10⁻⁴|1\.5\\times ?10^{-4}|4\.44e7|4\.4e7|4\.4×10|1\.27e59|1\.27×10⁵⁹|1\.27\\times ?10^{59}|2\.506e-2|2\.506e-02|2\.5e-2|2\.5e-02|0\.00015|0\.025\b|0\.02506" -- scripts/ results/
```

Quedan estas apariciones:

| fichero:línea | texto | por qué se queda |
|---|---|---|
| results/profile_sigma_x/profile_sigma_x.json:8 | `sigma_X_95CL_upper_profile_grid`: 0.025063… | clave `_grid` |
| results/profile_sigma_x/profile_sigma_x.json:278 | `"sigma_X": 0.025063…` | fila del perfil en ese nodo de malla; añadida a la lista 2b |
| results/profile_sigma_x/profile_sigma_x.json:304 | nota "Former 1.5e-4 … corrected to 2.5e-2 (grid node), then to the refined root 2.6e-2" | lista de obsoletos 2b ("profile json 302", ahora 304; actualizada) |
| scripts/profile_sigma_x_desi.py:124 | la misma nota | la escribe el script en el json; añadida a 2b |
| results/amplification_routes/VERDICT.md:58 | "(… the former 216/288, 72/288 and σ ≲ 4×10⁻⁵ were computed with the old 1.5×10⁻⁴ target)" | anotación del valor anterior; añadida a 2b |
| results/sdiff_discrimination/summary.txt:44 | "(was 1.5e-4 / 2.5e-2)" | anotación del valor anterior; añadida a 2b |
| results/euclid_mcmc/mcmc_summary.txt:8 | "targets 1e-5..2.5e-2 require multi-survey power beyond this mock" | mock retirado (cabecera RETIRED en l.1; ya en 2b, "Otros"); añadido a la lista de obsoletos de 2b |
| results/amplification_routes/route1_N_eff_grid.csv:8 | `…,0.025060996194449348,…` | **no es la cota**: es un punto de la malla de route1 (L/R_H = 3.99e-60, p = 2) cuyo σ_eff coincide numéricamente. El grep lo saca por `0.02506` |

En papers/, manuscript/ y notes/ quedan muchas más; son B2.

## 7. Ficheros cambiados

- **Scripts:**
  - `scripts/desi_bound.py` (nuevo);
  - `profile_sigma_x_desi.py`, `gap_two_targets.py`, `amplifier_audit.py`, `slip_bridge.py`;
  - `amplification/route1_local_causal_set_seed.py`;
  - `desqueezing/{cosmological_mapping,euclid_protocol_forecasts,sdiff_discrimination,desqueezing_relax_time}.py`;
  - `desi_dr2_real_bao_test.py`;
  - `analisis_cmb/realizaciones.py`, `analisis_cmb/verificacion_numeros.py`.
- **Resultados:** los de la sección 3, más `results/profile_sigma_x/*`, `figures/{profile_sigma_x,euclid_A0_theta_plane,sdiff_A0_min_for_detection,sdiff_residual_vs_theta,desqueezing_relax_time}.png` y `analisis_cmb/resultados/takahashi/sigma_2.633e-02/`.
- **Texto en analisis_cmb/:** `AUDITORIA.md` (decisión en 2a y en B; lista 2b), `RESULTADO.md` (§2, §3), `verificacion_numeros.md`.
- **Nuevos en `analisis_cmb/bloque_B/`:** `antes/`, `logs/` y este fichero.

---

## B1b (2026-10-17)

### 1. `scripts/desqueezing/desqueezing_relax_time.py`: fórmula de t½ y rutas

**La fórmula buena es t½ = ln2/γ.**
- En el baño térmico markoviano, d⟨a²⟩/dt = −(2iω + γ)⟨a²⟩: los términos en n_th se cancelan.
- n_th fija el estado estacionario, Var(q)_ss = n_th + 1/2, no la tasa.
- La propia simulación lo confirma: t½ = 0.28883 para todo n_th a γ = 2.4.

`gamma_eff` = γ(1+2n_th) se conserva como columna del csv con un comentario, pero ya no se usa como tasa. La figura (panel 3) pasa de t½ frente a 1/γ_eff a t½ frente a 1/γ.

Líneas del script:

| línea | antes | después |
|---|---|---|
| 110–113 | `# Para baño térmico Markoviano, \|⟨a²⟩\| ~ e^{-γ(1+2 n_th) t}`; `t_half_theory = ln2 / gamma_eff` | `# … d⟨a²⟩/dt = -(2iω + γ)⟨a²⟩ … n_th fija el estado estacionario …, no la tasa`; `t_half_theory = ln2 / gamma` |
| 256 | `inv = 1.0 / geff` | `inv = 1.0 / gamma` (por fila) |
| 260–263 | etiqueta "teoría ln2 / γ_eff", eje "1 / γ_eff [γ_eff = γ(1+2n_th)]", título "t_relax ∝ 1/γ_eff" | "teoría ln2 / γ", eje "1 / γ", título "t_relax ∝ 1/γ (independiente de n_th)" |
| 266, 271 | fit "a / γ_eff" | fit "a / γ" |
| 302–305 | "Teoría Markoviana: \|⟨a²⟩\|(t) ∝ exp[-γ(1+2n_th) t]"; "⇒ t_½ = ln(2) / [γ(1+2n_th)] = ln2 / γ_eff"; "Fit lineal … / γ_eff" | "… ∝ exp[-γ t] (n_th fija el estado estacionario, no la tasa)"; "⇒ t_½ = ln(2) / γ"; "Fit lineal … / γ" |
| 313 | "Aquí: γ_eff = γ(1+2n_th) juega el rol de tasa…" | "Aquí: γ juega el rol de tasa…" |
| 327–329, 338 | rutas absolutas (`fig_path`, `csv_path`, `npz_path`, `RESULTS`) | `os.path.relpath(…, ROOT)` |

**Resumen regenerado** (`results/desqueezing_relax_time/desqueezing_relax_time_resumen.txt`; log `logs/desqueezing_relax_time_B1b.log`):

| línea | antes (B1) | después (B1b) |
|---|---|---|
| 4 | Teoría Markoviana: \|⟨a²⟩\|(t) ∝ exp[-γ(1+2n_th) t] | Teoría Markoviana: \|⟨a²⟩\|(t) ∝ exp[-γ t] (n_th fija el estado estacionario, no la tasa) |
| 5 | ⇒ t_½ = ln(2) / [γ(1+2n_th)] = ln2 / γ_eff | ⇒ t_½ = ln(2) / γ |
| 7 | Fit lineal t_½ ≈ 0.7650 / γ_eff + 0.0560 | Fit lineal t_½ ≈ 0.6931 / γ + 0.0000 |
| 9 | Ratio sim/teoría (media ± std): 1.4701 ± 0.5798 | Ratio sim/teoría (media ± std): **1.0001 ± 0.0001** |
| 13 | Aquí: γ_eff = γ(1+2n_th) juega el rol de tasa de restauración de lisura | Aquí: γ juega el rol de tasa de restauración de lisura |
| 20–24 | t½_th = 0.96270 / 0.48135 / 0.24068 / 0.16045 / 0.12034 (γ = 0.6 … 4.8, n_th = 0.1) | t½_th = 1.15525 / 0.57762 / 0.28881 / 0.19254 / 0.14441 |
| 26–29 | t½_th = 0.26256 / 0.18051 / 0.14441 / 0.09627 (γ = 2.4, n_th = 0.05 … 1) | t½_th = 0.28881 en las cuatro |
| 31–33 | `figura: /home/ashpokemon/Proyectos/…/results/desqueezing_relax_time/desqueezing_relax_time.png` (y csv, npz) | `figura: results/desqueezing_relax_time/desqueezing_relax_time.png` (y csv, npz) |

Otros ficheros:
- `results/desqueezing_relax_time/t_relax_tabla.csv`: solo cambia la columna `t_half_theory`; las columnas de la simulación son idénticas.
- La figura se regenera y la copia en `figures/` se actualiza.

### 2. `notes/mapping_tables.md` restaurada y salida literal del script aparte

- **`cosmological_mapping.py`:** `DOCS = RESULTS` y el fichero pasa a llamarse `mapping_tables_script.md`. La salida literal va ahora a `results/cosmological_mapping/mapping_tables_script.md` (nuevo). Los csv y `mapping_summary.txt` no cambian al volver a ejecutarlo.
- **`notes/mapping_tables.md`:** vuelve a la versión de `bloque_B/antes/notes/` (prosa original). Solo se sustituyen los números de las 8 filas DESI de la Tabla 2, tomados de `results/cosmological_mapping/table2_A0_required.csv`, y se añade la línea de fuente.

Diff de la nota frente a la versión anterior a B1:

```diff
-| theta_MLE=0.001 | DESI_0.94 | DESI_limit_2.5e-2 | 1.501e-04 | 1.501e+57 |
+| theta_MLE=0.001 | DESI_0.94 | DESI_limit_2.633e-2 | 2.635e-02 | 2.635e+59 |
-| theta_MLE=0.001 | recomb | DESI_limit_2.5e-2 | 1.511e-04 | 1.511e+57 |
+| theta_MLE=0.001 | recomb | DESI_limit_2.633e-2 | 2.651e-02 | 2.651e+59 |
-| theta=1.2 (old) | DESI_0.94 | DESI_limit_2.5e-2 | 4.634e-04 | 4.634e+57 |
+| theta=1.2 (old) | DESI_0.94 | DESI_limit_2.633e-2 | 8.134e-02 | 8.134e+59 |
-| theta=1.2 (old) | recomb | DESI_limit_2.5e-2 | 6.630e-01 | 6.630e+60 |
+| theta=1.2 (old) | recomb | DESI_limit_2.633e-2 | 1.164e+02 | 1.164e+63 |
-| theta=ln2~0.693 (Hubble half-life) | DESI_0.94 | DESI_limit_2.5e-2 | 2.878e-04 | 2.878e+57 |
+| theta=ln2~0.693 (Hubble half-life) | DESI_0.94 | DESI_limit_2.633e-2 | 5.051e-02 | 5.051e+59 |
-| theta=ln2~0.693 (Hubble half-life) | recomb | DESI_limit_2.5e-2 | 1.913e-02 | 1.913e+59 |
+| theta=ln2~0.693 (Hubble half-life) | recomb | DESI_limit_2.633e-2 | 3.358e+00 | 3.358e+61 |
-| theta=0.099 (half over recomb path) | DESI_0.94 | DESI_limit_2.5e-2 | 1.646e-04 | 1.646e+57 |
+| theta=0.099 (half over recomb path) | DESI_0.94 | DESI_limit_2.633e-2 | 2.890e-02 | 2.890e+59 |
-| theta=0.099 (half over recomb path) | recomb | DESI_limit_2.5e-2 | 3.000e-04 | 3.000e+57 |
+| theta=0.099 (half over recomb path) | recomb | DESI_limit_2.633e-2 | 5.266e-02 | 5.266e+59 |
+
+DESI rows: `results/cosmological_mapping/table2_A0_required.csv` (DESI bound 2.633e-2 = `sigma_X_95CL_upper_profile` of `results/profile_sigma_x/profile_sigma_x.json`, 2026-10-17).
```
