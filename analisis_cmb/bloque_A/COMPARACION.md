# Bloque A: etiquetas de covarianza, CPL y joint reajustados (2026-10-17)

Bloque A de `analisis_cmb/AUDITORIA.md` → "Orden propuesto de corrección". No se ha tocado nada de B, C ni D.

- **Copia de los resultados anteriores:** `analisis_cmb/bloque_A/antes/`, es decir `results/eos_cpl_desi_dr2/`, `results/joint_w0wa_sigma/` y `figures/joint_w0wa_sigma.png`, tal como estaban en `859b33e`.
- **Salidas de consola:** `analisis_cmb/bloque_A/logs/`.

## 1. Estado de partida en esta máquina

`load_alpha_dv(prefer_file=True, use_full_cov=True)` (`scripts/desi_dr2_data.py` l.145–212) es la llamada que hacen los dos scripts. Ejecutada antes de cambiar nada da:

```
BAO_DR2 is_dir: True  FIG6 file: True
source : data/desi_dr2_local/dr2_data/dr2-bao-zenodo/figure6/DESI_DR2_alpha_DV_over_rs.txt
cov_mode : full_gaussian_bao_projected_13x13
gaussian_bao_source : data/.../desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_mean.txt
cov_fallback_reason : None
alpha : [1.004996 1.007276 0.988378 0.999921 0.999649 1.01097  1.002294]
fallback alpha: [1.003  0.9947 1.0016 0.996  1.002  0.9963 1.0008]
sigma (proj): [0.009629 0.007945 0.006937 0.005121 0.007505 0.015595 0.008224]
cov off-diag max: 0.0
```

- **Partida actual:** con el código de hoy, los dos scripts ya arrancaban en modo completo y con α figure6. Los datos locales (`data/desi_dr2_local/`, sin versionar) lo permiten.
- **Procedencia de los resultados archivados**, reproducida aquí bit a bit:
  - eos: α fallback y σ fallback diagonales. Da χ²_ΛCDM = 1.6647026678954178, igual que su json.
  - joint: α figure6 y σ figure6 diagonales (`use_full_cov=False`). Da lnL_ΛCDM = 24.895484644425657, igual que su json.
- **Proyección:** la covarianza proyectada a α es diagonal entre bins (el 13×13 público es diagonal por bloques). Lo que cambia frente a "diagonal" son las σ: proyectadas frente a figure6, con diferencias < 1.2 %.

## 2. Cambios en scripts: solo etiquetas

- **`scripts/eos_efectiva.py`:**
  - l.5 (docstring) y l.33 (comentario).
  - `data` del json (l.375) y la nota del json (l.404).
  - Línea "Data:" del txt (l.414) y la última línea de interpretación del txt (l.443, antes "diagonal approximation").
  - Todas salen ahora de `COV_LABEL`, que se calcula a partir del `cov_mode` del loader. El json y el txt añaden `cov_mode` y `alpha_source`.
- **`scripts/joint_w0wa_sigma_desi.py`:**
  - `data` del json (l.263): ahora sale de `COV_LABEL`, y el json añade `alpha_source`.
  - Título de la figura (l.300): ahora muestra el `cov_mode`.
- Ninguna fórmula ni dato cambia: el diff de los dos scripts solo toca cadenas, el diccionario `COV_LABELS` y dos claves nuevas del json.

## 3. Ejecución

| paso | comando | tiempo | salida |
|---|---|---|---|
| 1 | `MPLBACKEND=Agg python3 scripts/eos_efectiva.py` | 2 min 12 s | `logs/eos_efectiva.log` |
| 2 | `MPLBACKEND=Agg python3 scripts/joint_w0wa_sigma_desi.py` | 2 s | `logs/joint_w0wa_sigma_desi.log` |
| 3 | `python3 scripts/verify_bao_pipelines.py` | 2 min 13 s | `logs/verify_bao_pipelines.log`: **VERIFY: OK** |

`verify_bao_pipelines.py` vuelve a ejecutar `eos_efectiva.py`; el json resultante es idéntico byte a byte al del paso 1. También ejecuta `ou_bao_likelihood.py`, que solo escribe en `plots/`, una carpeta no versionada.

## 4. Comprobaciones

**(a) cov_mode y fuente.** Los dos json nuevos dicen:
- `"cov_mode": "full_gaussian_bao_projected_13x13"`;
- `"alpha_source": "data/desi_dr2_local/dr2_data/dr2-bao-zenodo/figure6/DESI_DR2_alpha_DV_over_rs.txt"`;
- `data` coherente con lo anterior: "…full measurement cov (DESI Gaussian BAO 13x13 projected to alpha)…".

Ficheros: `results/joint_w0wa_sigma/joint_w0wa_sigma.json`, `results/eos_cpl_desi_dr2/eos_cpl_summary.json`.

**(b) logL de ΛCDM del joint frente al perfil.**
- `results/joint_w0wa_sigma/joint_w0wa_sigma.json` → `models[lcdm].logL` = 24.899562520639115.
- `results/profile_sigma_x/profile_sigma_x.json` → `ll_lcdm` = 24.899562520639115.
- Diferencia: **0.0** exacto.

La predicción del joint para ΛCDM, `alpha_pred_cpl(−1, 0)`, es exactamente 1 en los 7 bins (diferencia 0.0 en cada uno). Por tanto el residuo es el mismo que usa el perfil, α − 1 (`profile_sigma_x_desi.py` l.38), con la misma covarianza y la misma `logL_gaussian`.

**(c) Cota del perfil intacta.**
- sha256 de `results/profile_sigma_x/profile_sigma_x.json`: `83558423a3690489…` antes y después.
- `git status results/profile_sigma_x/` está limpio.
- `results/desi_cov/` también está limpio: el loader lo reescribe con contenido idéntico.

## 5. Antes / después

Convenciones:
- **lnL del joint** = log-verosimilitud gaussiana completa, con log det C y n ln 2π (`logL_gaussian`).
- **lnL de eos** = −χ²/2, sin normalización (`eos_efectiva.py` l.120–124). Las dos no son comparables entre sí.
- **ΔAIC frente a ΛCDM** = AIC(modelo) − AIC(ΛCDM), con AIC = 2k − 2 lnL y k(ΛCDM) = 0.
- **El σ de eos** es el parámetro del toy δw (`eos_efectiva.py` l.129), no σ_X.

### Joint (`scripts/joint_w0wa_sigma_desi.py`)

Antes: `analisis_cmb/bloque_A/antes/results/joint_w0wa_sigma/joint_w0wa_sigma.json` (= commit 9e8b6bf/9e7c765). Después: `results/joint_w0wa_sigma/joint_w0wa_sigma.json`.

| modelo | magnitud | antes | después |
|---|---|---|---|
| — | cov_mode | diagonal_figure6_errors (con `data` = "full meas. cov") | full_gaussian_bao_projected_13x13 |
| — | vector α | figure6 (`DESI_DR2_alpha_DV_over_rs.txt`) | figure6 (el mismo) |
| ΛCDM | lnL | 24.895485 | 24.899563 |
| CPL | w0 | −1.00165 | **−1.00095** |
| CPL | wa | +0.01725 | **+0.01370** |
| CPL | lnL | 24.903781 | 24.907349 |
| CPL | ΔAIC frente a ΛCDM | +3.983 | +3.984 |
| OU, fondo fijo | θ | 1.890 | 1.898 (θ no está determinado con σ_X → 0) |
| OU, fondo fijo | σ_X | 3.0×10⁻¹⁰ | 3.8×10⁻¹⁰ |
| OU, fondo fijo | lnL | 24.895485 | 24.899563 |
| OU, fondo fijo | ΔAIC frente a ΛCDM | +4.000 | +4.000 |
| joint | w0 | −1.00165 | **−1.00095** |
| joint | wa | +0.01725 | **+0.01370** |
| joint | θ | 0.555 | 0.734 (no determinado) |
| joint | σ_X | 4.2×10⁻¹⁰ | **5.3×10⁻¹⁰** |
| joint | lnL | 24.903781 | 24.907349 |
| joint | ΔAIC frente a ΛCDM | +7.983 | **+7.984** |
| joint | ΔAIC frente a CPL | +4.000 | **+4.000** |
| joint | ΔlnL (joint − CPL) | 7×10⁻¹³ | 3×10⁻¹³ |

### CPL y extensión anidada (`scripts/eos_efectiva.py`)

Antes: `analisis_cmb/bloque_A/antes/results/eos_cpl_desi_dr2/eos_cpl_summary.json` (generado en f36cc79). Después: `results/eos_cpl_desi_dr2/eos_cpl_summary.json`.

| modelo | magnitud | antes | después |
|---|---|---|---|
| — | cov_mode | diagonal (σ fallback; `data` = "diagonal errors"; sin clave cov_mode) | full_gaussian_bao_projected_13x13 |
| — | vector α | fallback hard-coded [1.0030, 0.9947, 1.0016, 0.9960, 1.0020, 0.9963, 1.0008] | figure6 [1.004996, 1.007276, 0.988378, 0.999921, 0.999649, 1.010970, 1.002294] |
| ΛCDM | χ² / lnL | 1.6647 / −0.8324 | 4.4902 / −2.2451 |
| CPL | w0 | −0.98990 | **−1.00095** |
| CPL | wa | −0.01590 | **+0.01370** |
| CPL | χ² / lnL | 1.4116 / −0.7058 | 4.4746 / −2.2373 |
| CPL | ΔAIC frente a ΛCDM | +3.747 | +3.984 |
| anidado | w0 / wa | −0.98990 / −0.01590 | −1.00095 / +0.01370 |
| anidado | σ (toy δw) | 6.3×10⁻⁷ | 8.8×10⁻⁸ |
| anidado | θ | 1.734 | 4.071 (no determinado con σ → 0) |
| anidado | lnL | −0.7058 | −2.2373 |
| anidado | ΔAIC frente a ΛCDM | +7.747 | **+7.984** |
| anidado | ΔAIC frente a CPL | +4.000 | **+4.000** |

**Resultado:**
- Con la misma covarianza y el mismo vector α, los dos ajustes CPL dan el mismo punto, (w0, wa) = (−1.00095, +0.01370), con ΔAIC idénticos.
- La discrepancia anterior (−0.990/−0.016 frente a −1.002/+0.017) venía del vector α, como decía AUDITORIA §3.

**CLAIMS.md C3** ("Nested free background plus stochastic extension is not preferred, ΔAIC ≈ +4") **sigue en pie** con los números nuevos:
- ΔAIC(extensión − CPL) = **+4.000** en los dos ajustes: joint +3.9999999999995, eos +4.000000000000027.
- La amplitud estocástica va a 0: σ_X = 5.3×10⁻¹⁰ en el joint, σ_toy = 8.8×10⁻⁸ en eos.
- La ganancia en lnL es ≤ 3×10⁻¹³.
- El +4 es exactamente la penalización por los dos parámetros extra.
- Frente a ΛCDM, la extensión anidada da ΔAIC = +7.984 y el CPL solo, +3.984.

## 6. Líneas de texto corregidas (tabla del punto 3 de AUDITORIA.md)

Solo las líneas de la tabla, sobre el texto existente. Matices:
- `papers/resume.txt` l.25 ya decía "Gaussian BAO 13x13 projected" y no se toca; ahora concuerda con l.28–29 y l.33.
- `PREPRINT.md` l.22 no estaba en la lista de cambios y no se toca.
- En `NOTAS_ANALISIS_CMB.md` y `RESULTADO.md` se conservan los valores antiguos, porque son los que usó el análisis CMB como curva de referencia. Se añade la nota del α fallback y los valores reajustados.
- En stochastic l.922 se sustituye "[reproduces Sec. 4.2 exactly]". Con los números nuevos ya no reproduce el run antiguo de Sec. 4.2; ahora cita el modelo `ou_fixed_bg` del json.
- En l.932 se quita "with the real DESI covariance matrix", que ya está aplicada.
- En PREPRINT l.76 se quita "Full off-diagonal BAO covariances are a natural extension for tightening the limit (Section 7)", porque la covarianza ya es la oficial y no tiene bloques entre bins.

| fichero:línea | antes | después |
|---|---|---|
| papers/resume.txt:28 |   w0 ≈ −0.99 ⏎   wa ≈ −0.02 ⏎  |   w0 ≈ −1.001 ⏎   wa ≈ +0.014 ⏎  |
| papers/resume.txt:33 | - Baseline: ΛCDM (diagonal covariance only) | - Baseline: ΛCDM (same projected 13x13 covariance) |
| manuscript/PREPRINT.md:76 | The baseline pipeline adopts **diagonal** measurement errors in the residual likelihood. Full off-diagonal BAO covariances are a natural extension for tightening the limit (Section 7); they are not required | The baseline pipeline adopts the official DESI Gaussian BAO 13×13 covariance projected to isotropic α (block-diagonal across bins) in the residual likelihood. Off-diagonal covariances between bins are absent in the public product; they are not required |
| manuscript/PREPRINT.md:78 | uses a CPL background near $\Lambda$CDM ($ w_0\approx -0.99$, $w_a\approx -0.02$) as the smooth reference | uses a CPL background near $\Lambda$CDM ($ w_0\approx -1.001$, $w_a\approx +0.014$; `results/eos_cpl_desi_dr2/`, projected covariance, figure6 α) as the smooth reference |
| papers/stochastic-dark-energy-desi-dr2.md:430 | With only 7 BAO bins and using a diagonal covariance approximation, the model space | With only 7 BAO bins and the official Gaussian BAO covariance projected to α (block-diagonal across bins), the model space |
| papers/stochastic-dark-energy-desi-dr2.md:440 | yields best-fit near $w_0 \approx -0.99$, $w_a \approx -0.02$ (essentially $\Lambda$CDM) and $\sigma_X \approx 4.5 \times 10^{-5}$ (numerical floor). | yields best-fit near $w_0 \approx -1.001$, $w_a \approx +0.014$ (essentially $\Lambda$CDM) and $\sigma_X \approx 5 \times 10^{-10}$ ($\to 0$), $\Delta\text{AIC} = +4.00$ vs CPL (`results/joint_w0wa_sigma/`, projected covariance). |
| papers/stochastic-dark-energy-desi-dr2.md:922 | **(a)** Background fixed at $w_0 = -1$, $w_a = 0$ [reproduces Sec. 4.2 exactly]: $\sigma_X \to 5 \times 10^{-5}$, $\theta \to 0.0010$, $\ln L = 27.013$. | **(a)** Background fixed at $w_0 = -1$, $w_a = 0$ [`results/joint_w0wa_sigma/`, model `ou_fixed_bg`]: $\sigma_X \to 4 \times 10^{-10}$, $\theta = 1.90$ (unconstrained at $\sigma_X \to 0$), $\ln L = 24.8996$. |
| papers/stochastic-dark-energy-desi-dr2.md:926 | The optimum falls near $w_0 \approx -0.99$, $w_a \approx -0.02$ (nearly pure $\Lambda$CDM), with $\sigma_X \approx 4.5 \times 10^{-5}$ — the same numerical floor as (a) — and a log-likelihood improvement of only $(+0.13)$ over the 2-parameter fit, | The optimum falls near $w_0 \approx -1.001$, $w_a \approx +0.014$ (nearly pure $\Lambda$CDM), with $\sigma_X \approx 5 \times 10^{-10}$ — $\sigma_X \to 0$ as in (a) — and a log-likelihood improvement of only $(+3 \times 10^{-13})$ over the 2-parameter fit, |
| papers/stochastic-dark-energy-desi-dr2.md:932 | ($S(z)$ held fixed rather than recomputed at each $w_0, w_a$; diagonal DESI covariance) pending a full re-analysis with the real DESI covariance matrix and a proper MCMC sampler | ($S(z)$ held fixed rather than recomputed at each $w_0, w_a$; DESI Gaussian BAO covariance projected to α) pending a full re-analysis with a proper MCMC sampler |
| papers/neutrino-de-degeneracy-and-sigma-x.md:30 | The MLE I already have is $ (w_0,w_a)\approx(-0.99,-0.016)$ and $\sigma_X\to 0$. | The MLE I already have is $ (w_0,w_a)\approx(-1.001,+0.014)$ and $\sigma_X\to 0$ ($\approx 5\times10^{-10}$; `results/joint_w0wa_sigma/`, projected covariance). |
| papers/covariance-desi-dr2-full.md:68 | \| `joint_w0wa_sigma_desi.py` \| joint $\{w_0,w_a,\sigma_X\}$ \| | \| `joint_w0wa_sigma_desi.py` \| joint $\{w_0,w_a,\sigma_X\}$ (results regenerated with $C_\alpha$, 2026-10-17) \| |
| papers/covariance-desi-dr2-full.md:72 | \| `eos_efectiva.py` \| CPL + nested extension \| | \| `eos_efectiva.py` \| CPL + nested extension (results regenerated with $C_\alpha$ and figure6 α, 2026-10-17) \| |
| analisis_cmb/NOTAS_ANALISIS_CMB.md:62 | Ajustes CPL del repo (no reajustados): | Ajustes CPL del repo (no reajustados en este análisis; reajustados en `bloque_A/`, 2026-10-17: w0 = −1.00095, wa = +0.01370 en los dos, cov completa y α figure6): |
| analisis_cmb/NOTAS_ANALISIS_CMB.md:64 |   **errores diagonales** ("7 bins, diagonal errors"), no con la covarianza completa del perfil OU. ⏎ - `results/joint |   **errores diagonales** ("7 bins, diagonal errors") y vector α de fallback hard-coded (no figure6), no con la covarianza completa del perfil OU. ⏎ - `results/joint |
| analisis_cmb/NOTAS_ANALISIS_CMB.md:181 |   wa = −0.01590. ⚠ Ese ajuste usa **errores diagonales** ("7 bins, diagonal errors"), no la |   wa = −0.01590. ⚠ Ese ajuste usa **errores diagonales** ("7 bins, diagonal errors") y α de fallback hard-coded (no figure6), no la |
| analisis_cmb/RESULTADO.md:53 | Ajustes CPL del repo (no reajustados): | Ajustes CPL del repo (no reajustados en este análisis; reajustados en `bloque_A/`, 2026-10-17: w0 = −1.00095, wa = +0.01370 en los dos, cov completa y α figure6): |
| analisis_cmb/RESULTADO.md:54 | hecho con errores diagonales ("7 bins, diagonal errors"), no con la covarianza completa del perfil OU. | hecho con errores diagonales ("7 bins, diagonal errors") y vector α de fallback hard-coded (no figure6), no con la covarianza completa del perfil OU. |

## 7. Ficheros cambiados

- **Scripts:** `scripts/eos_efectiva.py` y `scripts/joint_w0wa_sigma_desi.py`, solo etiquetas.
- **Resultados regenerados:**
  - `results/eos_cpl_desi_dr2/eos_cpl_summary.{json,txt}`;
  - `results/joint_w0wa_sigma/joint_w0wa_sigma.{json,txt,png}`;
  - `figures/joint_w0wa_sigma.png`.
- **Texto:** los de la tabla 6.
- **Nuevos en `analisis_cmb/bloque_A/`:** `antes/`, `logs/` y este fichero.
