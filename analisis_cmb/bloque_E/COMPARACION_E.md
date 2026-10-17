# Bloque E: lectura de revisor del PREPRINT (2026-10-17)

Cambios sobre el texto existente, aplicados con el mismo motor que C: cada cambio exige que el texto de origen aparezca una sola vez y registra la línea, el antes y el después. Las líneas son las del fichero **antes de E** (commit `e285ec3`).

## Punto 2: `scripts/ou_bao_likelihood.py` y `results/ou_qnm_mle/`

**Cambios en el script** (`git diff`):
- La línea base ΛCDM usaba una covarianza **diagonal** (`np.diag(sigma**2)`), mientras que los ajustes OU/QNM usaban la oficial (`C_MEAS`). Ahora las tres verosimilitudes usan `C_MEAS`.
- El script escribe además `results/ou_qnm_mle/ou_qnm_mle.{json,txt}`, con θ, σ_X, ω_R, lnL, ΔlnL frente a ΛCDM, AIC, BIC, ΔAIC y ΔBIC de OU y QNM, los límites de búsqueda y si el ajuste cae en uno de ellos. La figura sigue en `plots/`.

```diff
diff --git a/scripts/ou_bao_likelihood.py b/scripts/ou_bao_likelihood.py
index 31743e1..2f67eaa 100644
--- a/scripts/ou_bao_likelihood.py
+++ b/scripts/ou_bao_likelihood.py
@@ -328,8 +328,8 @@ def lag_correlation_whitened(residuals, cov_matrix, lag_k):
 
 residuals = alpha - 1.0
 
-# ΛCDM baseline (diagonal covariance)
-C_LCDM    = np.diag(sigma**2)
+# ΛCDM baseline: same official measurement covariance as the OU/QNM fits
+C_LCDM    = np.array(C_MEAS, dtype=float, copy=True)
 logL_LCDM = log_likelihood_gaussian(residuals, C_LCDM)
 
 # OU calibrated (fixed parameters from paper Section 3.2)
@@ -472,6 +472,45 @@ print("  F7:  ω_R → 0 in Euclid fit                  → OU pure recovered, n
 print()
 print(SEP)
 
+# ============================================================
+# SECTION 9b: MACHINE-READABLE SUMMARY (results/ou_qnm_mle/)
+# ============================================================
+import json as _json
+_OUT = _ROOT / "results" / "ou_qnm_mle"
+_OUT.mkdir(parents=True, exist_ok=True)
+_SIGMA_X_MIN = float(np.exp(-10.0))   # search bound: log(sigma_X^2) >= -20
+
+def _summ(r):
+    d = {k: float(r[k]) for k in ("theta", "sigma_X", "omega_R", "logL", "AIC", "BIC")}
+    d["k"] = int(r["k"])
+    d["dlogL_vs_LCDM"] = float(r["logL"] - logL_LCDM)
+    d["dAIC_vs_LCDM"] = float(r["AIC"] - (-2.0 * logL_LCDM))
+    d["dBIC_vs_LCDM"] = float(r["BIC"] - (-2.0 * logL_LCDM))
+    d["theta_at_lower_bound"] = bool(r["theta"] < 1.01 * THETA_MIN)
+    d["sigma_X_at_lower_bound"] = bool(r["sigma_X"] < 1.01 * _SIGMA_X_MIN)
+    return d
+
+_summary = {
+    "data": _d["source"], "cov_mode": _d.get("cov_mode"), "n_bins": int(len(z_eff)),
+    "LCDM": {"logL": float(logL_LCDM), "AIC": float(-2.0 * logL_LCDM),
+             "BIC": float(-2.0 * logL_LCDM), "cov": "same measurement covariance as OU/QNM"},
+    "OU": _summ(result_OU), "QNM": _summ(result_QNM),
+    "search_bounds": {"theta_min": THETA_MIN, "sigma_X_min": _SIGMA_X_MIN,
+                      "omega_R": [0.0, OMEGA_R_MAX], "log_sigma_X2": [-20.0, -1.0]},
+}
+(_OUT / "ou_qnm_mle.json").write_text(_json.dumps(_summary, indent=2) + "\n", encoding="utf-8")
+_lines = [f"OU / QNM free MLE on {_summary['data']} (cov_mode: {_summary['cov_mode']}, {_summary['n_bins']} bins)",
+          f"LCDM: logL = {logL_LCDM:.6f}  (AIC = BIC = {-2.0 * logL_LCDM:.3f}, k = 0)",
+          f"search bounds: theta >= {THETA_MIN:g}, sigma_X >= {_SIGMA_X_MIN:.3e}, omega_R in [0, {OMEGA_R_MAX:g}]",
+          f"{'model':<6} {'theta':>9} {'sigma_X':>11} {'omega_R':>8} {'logL':>10} {'dlogL':>9} {'AIC':>9} {'BIC':>9} {'dAIC':>7} {'dBIC':>7}  bounds"]
+for _name in ("OU", "QNM"):
+    _r = _summary[_name]
+    _b = ",".join(n for n, f in (("theta", _r["theta_at_lower_bound"]), ("sigma_X", _r["sigma_X_at_lower_bound"])) if f) or "-"
+    _lines.append(f"{_name:<6} {_r['theta']:>9.4g} {_r['sigma_X']:>11.4e} {_r['omega_R']:>8.4f} {_r['logL']:>10.6f} "
+                  f"{_r['dlogL_vs_LCDM']:>+9.2e} {_r['AIC']:>9.3f} {_r['BIC']:>9.3f} {_r['dAIC_vs_LCDM']:>+7.3f} {_r['dBIC_vs_LCDM']:>+7.3f}  at lower bound: {_b}")
+(_OUT / "ou_qnm_mle.txt").write_text("\n".join(_lines) + "\n", encoding="utf-8")
+print("Summary written: results/ou_qnm_mle/ou_qnm_mle.{json,txt}")
+
 # ============================================================
 # SECTION 10: FIGURES (4 panels)
 # ============================================================
```

**Ejecución:** `python scripts/ou_bao_likelihood.py`, con código de salida 0. El log está en `bloque_E/ou_bao_likelihood.log`, sin rutas locales. Salida:

```
OU / QNM free MLE on data/desi_dr2_local/dr2_data/dr2-bao-zenodo/figure6/DESI_DR2_alpha_DV_over_rs.txt (cov_mode: full_gaussian_bao_projected_13x13, 7 bins)
LCDM: logL = 24.899563  (AIC = BIC = -49.799, k = 0)
search bounds: theta >= 0.001, sigma_X >= 4.540e-05, omega_R in [0, 10]
model      theta     sigma_X  omega_R       logL     dlogL       AIC       BIC    dAIC    dBIC  bounds
OU         266.1  4.5400e-05   0.0000  24.899523 -3.95e-05   -45.799   -45.907  +4.000  +3.892  at lower bound: sigma_X
QNM      0.04303  4.5400e-05   6.9313  24.899524 -3.82e-05   -43.799   -43.961  +6.000  +5.838  at lower bound: sigma_X
```

- **Comprobación:** lnL(ΛCDM) = **24.899563**, que coincide con 24.899562520639115 (`ll_lcdm` de `results/profile_sigma_x/profile_sigma_x.json`; bloque A). El 27.013 del run antiguo venía de la covarianza diagonal.
- **σ_X al borde:** en los dos modelos σ_X cae en el **borde inferior del rango de búsqueda**, 4.540×10⁻⁵ (ln σ_X² = −20). En el PREPRINT se dice así, y no "numerical floor". Con σ_X en el borde, θ y ω_R no quedan acotados: θ = 266 en OU; θ = 0.043 y ω_R = 6.93 en QNM.
- **Cuidado:** `verify_bao_pipelines.py` comprueba "ω_R = 0". Lo sigue encontrando en la línea "ω_R = 0 (fixed)" de OU, pero el ω_R libre de QNM ya no es 0, sino 6.93, sin acotar. La comprobación pasa por la línea de OU, no porque se cumpla en QNM.

## Cambios de texto (por fichero y línea)

### `manuscript/PREPRINT.md`

**l.9** — punto 4. Versión v2: Versión v2 en la cabecera

```diff
- | Date | July 2026 |
+ | Date | v2, October 2026 (v1: July 2026) |
```

**l.24** — punto 1. Everpresent como §4.2 y C10: Abstract: una frase de Everpresent tras la del tachyonic (C10, §4.2)

```diff
- Separately, a coherent tachyonic growth model with the correct rank-1 covariance induced by a single growing mode is strongly disfavoured: $\Delta\ln\mathcal{L}\approx -12.7$ relative to $\Lambda$CDM ($\Delta\chi^2\approx +25$).
- 
- At the DESI ceiling the OU residual is invisible to Planck 2018 (TT at ℓ < 30 and lensing; 95th percentile of Δχ² ≤ 0.067), so the CMB gives no leverage on σ_X, and its only imprint, a scatter σ(H0) ≈ 0.48 km/s/Mpc, is 12 times smaller than the H0 tension. A percent-level harmonic residual is likewise allowed by DESI and invisible in the CMB (§3.5).
+ Separately, a coherent tachyonic growth model with the correct rank-1 covariance induced by a single growing mode is strongly disfavoured: $\Delta\ln\mathcal{L}\approx -12.7$ relative to $\Lambda$CDM ($\Delta\chi^2\approx +25$). The two published implementations of Everpresent Λ, a zero-mean Λ with order-one fluctuations, do not fit the same BAO data together with the CMB: with r_d computed for each realization, small dark energy around recombination and ω_m within ±2σ of Planck, none fits as well as ΛCDM (0/20000 and 0/500000 realizations).
+ 
+ A CMB consistency check OU residual is invisible to Planck 2018 (TT at ℓ < 30 and lensing; 95th percentile of Δχ² ≤ 0.067), so the CMB gives no leverage on σ_X, and its only imprint, a scatter σ(H0) ≈ 0.48 km/s/Mpc, is 12 times smaller than the H0 tension. A percent-level harmonic residual is likewise allowed by DESI and invisible in the CMB (§3.5).
```

**l.26** — punto 5. Abstract (CMB): Abstract: párrafo del CMB a dos frases (los números detallados están en §3.5)

```diff
- A CMB consistency check with Planck 2018 (TT at ℓ < 30 and lensing bandpowers, via CAMB with the OU history mapped to w(a), H0 solved to keep θ_* fixed) does not distinguish the OU residual from ΛCDM at the DESI ceiling: the 95th percentile of Δχ² is ≤ 0.067 (TT) and ≤ 0.054 (lensing) at σ_X = 2.633×10⁻², and no realization reaches Δχ² > 4. Its only imprint is a scatter in the inferred H0 of 0.48 km/s/Mpc (θ = 1), about 12 times smaller than the Planck–SH0ES difference, so it neither tightens σ_X nor explains the H0 tension. A deterministic harmonic modulation X = A sin(2π ln a/P + φ) with A ≤ 0.025 and P ∈ [0.1, 2] is allowed by DESI DR2 BAO in 58 of 60 cases (Δχ²_BAO < 4) and is indistinguishable from ΛCDM in the CMB in all 60.
+ At the DESI ceiling the OU residual is invisible to Planck 2018 (TT at ℓ < 30 and lensing; 95th percentile of Δχ² ≤ 0.067), so the CMB gives no leverage on σ_X, and its only imprint, a scatter σ(H0) ≈ 0.48 km/s/Mpc, is 12 times smaller than the H0 tension. A percent-level harmonic residual is likewise allowed by DESI and invisible in the CMB (§3.5).
```

**l.30** — punto 3. Slip: δ_m: Slip: cautela sobre δ_m (abstract)

```diff
- Gravitational slip $\gamma=\Phi/\Psi$ is the correct operator for anisotropic leakage past volume-preserving (SDiff) structure, and with the BAO-bounded residual it is marginal, not starved: full anisotropy gives at most $\lvert\gamma-1\rvert\approx 5\times 10^{-2}$ at $z=0.5$, about $3.4\times$ below the current measurement (Maus et al.) and at the level of the Sakr et al. forecast floor ($\approx 1.0\times$), falling to $2.4\times$ and $4.7\times$ below that floor at $z=1$ and $1.5$. The σ_X → slip map is an order-of-magnitude estimate and the best case requires ε = 1, δ_m = 1 and z ≈ 0.5.
+ Gravitational slip $\gamma=\Phi/\Psi$ is the correct operator for anisotropic leakage past volume-preserving (SDiff) structure, and with the BAO-bounded residual it is marginal, not starved: full anisotropy gives at most $\lvert\gamma-1\rvert\approx 5\times 10^{-2}$ at $z=0.5$, about $3.4\times$ below the current measurement (Maus et al.) and at the level of the Sakr et al. forecast floor ($\approx 1.0\times$), falling to $2.4\times$ and $4.7\times$ below that floor at $z=1$ and $1.5$. The σ_X → slip map is an order-of-magnitude estimate: it holds for δ_m ≈ 1 (quasi-linear scales) and does not extrapolate to δ_m ≪ 1; ε < 1 or z ≳ 1 push it below reach.
```

**l.49** — punto 1. Everpresent como §4.2 y C10: §1.1: una sola frase de Everpresent; el resto pasa a §4.2

```diff
- These arguments come in two versions that should not be confused (Ahmed et al. 2004; Sorkin 2007; Zwane et al. 2018). The seed σ_0 ∼ 10⁻⁶¹ uses N ∼ 10¹²² from horizon-area (Bekenstein–Hawking) counting. With the four-volume counting of Everpresent Λ, δΛ tracks H² and the relative amplitude is of order one with no amplifier; the ∼10⁵⁶ gap is specific to the area-counting choice. Everpresent Λ (Ahmed et al. 2004; Zwane et al. 2018) is a different model from the one tested here: Λ has zero mean and fluctuates, with an amplitude fixed by causal-set counting of the past four-volume (δΛ ∼ V^{-1/2}), at the level of the critical density in every epoch and with no free amplitude. The residual X of this work is a small stationary perturbation around a non-zero mean Λ, with free σ_X and θ. The bound on σ_X does not, by itself, test Everpresent Λ.
+ Everpresent Λ (Ahmed et al. 2004; Zwane et al. 2018), in which the same counting applied to the four-volume gives order-one fluctuations, is a different model, tested directly in §4.2.
```

**l.67** — punto 1. Everpresent como §4.2 y C10: Guía de lectura §1.3

```diff
- | 4 | Exclusion of coherent tachyonic growth |
+ | 4 | Exclusions: coherent tachyonic growth (§4.1) and Everpresent Λ (§4.2) |
```

**l.138** — punto 6. Menores: §2.4: scripts de analisis_cmb/

```diff
- | `scripts/slip_bridge.py` | map $\sigma_X\to\lvert\gamma-1\rvert$ (order of magnitude) |
+ | `scripts/slip_bridge.py` | map $\sigma_X\to\lvert\gamma-1\rvert$ (order of magnitude) |
+ | `analisis_cmb/realizaciones.py`, `analisis_cmb/estadisticos.py` | CMB consistency: OU realizations through CAMB, $\Delta\chi^2$ against Planck 2018 (§3.5) |
+ | `analisis_cmb/oscilador/oscilador.py` | Harmonic residual against DESI BAO and Planck 2018 (§3.5) |
+ | `analisis_cmb/everpresent/bao_sim.py`, `rd_real.py`, `robustez.py` | Everpresent Λ realizations against DESI DR2 BAO (§4.2) |
```

**l.146** — punto 2. §3.1 con el run vigente: §3.1 con el run vigente (covarianza oficial; σ_X al borde inferior del rango de búsqueda)

```diff
- Under fixed smooth background and public BAO summary statistics, the MLE drives the stochastic amplitude to the numerical floor:
- 
- | Model | $\theta$ | $\sigma_X$ | $\omega_R$ | $\Delta\ln\mathcal{L}$ vs $\Lambda$CDM | AIC | BIC |
- |:------|:-----------|:-------------|:-------------|:-------------------------------------------|:----|:----|
- | OU (H0; old run, diagonal cov) | $0.001$ (boundary) | $5\times 10^{-5}$ | — | $0.00$ | $-50.03$ | $-50.14$ |
- | QNM (H1; old run, diagonal cov) | $0.001$ | $5\times 10^{-5}$ | $0.00$ | $0.00$ | $-48.03$ | $-48.19$ |
- 
- The data do not require a stationary stochastic residual. The preferred stochastic amplitude is consistent with zero within the resolution of the likelihood surface; reported floor values reflect optimizer bounds, not a detection.
+ Under a fixed smooth background and public BAO summary statistics, with the official DESI measurement covariance in all three likelihoods, the MLE drives the stochastic amplitude to the lower bound of the search range, $\sigma_X = 4.54\times 10^{-5}$ ($\ln\sigma_X^2 = -20$):
+ 
+ | Model | $\theta$ | $\sigma_X$ | $\omega_R$ | $\ln\mathcal{L}$ | $\Delta\ln\mathcal{L}$ vs $\Lambda$CDM | AIC | BIC |
+ |:------|:-----------|:-------------|:-------------|:------------------|:-------------------------------------------|:----|:----|
+ | $\Lambda$CDM | — | — | — | $24.8996$ | $0$ (ref) | $-49.80$ | $-49.80$ |
+ | OU (H0) | $266$ (unconstrained) | $4.54\times 10^{-5}$ (lower bound) | — | $24.8995$ | $-4.0\times 10^{-5}$ | $-45.80$ | $-45.91$ |
+ | QNM (H1) | $0.043$ (unconstrained) | $4.54\times 10^{-5}$ (lower bound) | $6.93$ (unconstrained) | $24.8995$ | $-3.8\times 10^{-5}$ | $-43.80$ | $-43.96$ |
+ 
+ Relative to $\Lambda$CDM, $\Delta$AIC $= +4.00$ (OU) and $+6.00$ (QNM), $\Delta$BIC $= +3.89$ and $+5.84$: the extra parameters bring no gain in likelihood. With $\sigma_X$ at the bound, $\theta$ and $\omega_R$ are not constrained. Source: [`results/ou_qnm_mle/ou_qnm_mle.json`](../results/ou_qnm_mle/ou_qnm_mle.json) (`scripts/ou_bao_likelihood.py`).
+ 
+ An earlier diagonal-covariance run (old run, [`archive/results/ou_bao_desi_dr2_run.txt`](../archive/results/ou_bao_desi_dr2_run.txt)) gave $\theta = 0.001$, $\sigma_X = 5\times 10^{-5}$ and AIC $-50.03$ / $-48.03$ (OU / QNM).
+ 
+ The data do not require a stationary stochastic residual. The preferred stochastic amplitude is consistent with zero within the resolution of the likelihood surface; values at the bound reflect the search range, not a detection.
```

**l.187** — punto 1. Everpresent como §4.2 y C10: §4 con subsecciones: §4.1 tachyonic

```diff
- ## 4. Results II — Coherent tachyonic growth excluded
- 
- A distinct hypothesis is that a dark-energy quantum fluid with negative effective mass undergoes Bogoliubov instability, yielding a coherent growing mode $\sigma_X(t)=\sigma_0 e^{t/t_c}$. The induced BAO residual covariance is **rank-1**, not OU-stationary. Using that covariance in an MLE on the same BAO summary statistics yields
+ ## 4. Results II — Model exclusions
+ 
+ ### 4.1 Coherent tachyonic growth
+ 
+ A distinct hypothesis is that a dark-energy quantum fluid with negative effective mass undergoes Bogoliubov instability, yielding a coherent growing mode $\sigma_X(t)=\sigma_0 e^{t/t_c}$. The induced BAO residual covariance is **rank-1**, not OU-stationary. Using that covariance in an MLE on the same BAO summary statistics yields
```

**l.201** — punto 1. Everpresent como §4.2 y C10: §4.2 nueva: Everpresent Λ (párrafo movido de §5.1 completo, tabla y condiciones; fuente EVERPRESENT_vs_DESI.md §3e)

```diff
- This exclusion is logically independent of the stationary OU bound: the covariance structure is different, and the correct structure must be used.
- 
+ This exclusion is logically independent of the stationary OU bound: the covariance structure is different, and the correct structure must be used.
+ 
+ ### 4.2 Everpresent Λ tested directly
+ 
+ The seed σ_0 ∼ 10⁻⁶¹ uses N ∼ 10¹²² from horizon-area (Bekenstein–Hawking) counting. With the four-volume counting of Everpresent Λ, δΛ tracks H² and the relative amplitude is of order one with no amplifier; the ∼10⁵⁶ gap is specific to the area-counting choice. Zwane et al. (2018) find that Everpresent Λ fits the CMB, BAO and SN data of that time about as well as ΛCDM and eases the low-redshift H0 and Lyα BAO tensions. The order-one alternative to the area-counting seed has been tested directly. In its two published implementations (Zwane et al. 2018, Model 2; Das et al. 2023, Model 1), computing r_d with each realization's own history, requiring small dark energy over z ∈ [1000, 5000] and ω_m within ±2σ of Planck, none of the realizations fits DESI DR2 BAO as well as ΛCDM, nor within Δχ² < 4: 0/20000 for Model 2 (< 1.5×10⁻⁴ at 95 %) and 0/300000 and 0/200000 for Model 1 with α = 0.01 and 0.005 (< 1.0×10⁻⁵, < 1.5×10⁻⁵). Without the CMB conditions, rare realizations pass, all with ω_m far from Planck or order-one dark energy at recombination. As small fluctuations around ΛCDM, DESI DR2 bounds the Zwane et al. amplitude to α < 1.85×10⁻² (μ fixed). These results hold for those phenomenological implementations, with the CMB only through a proxy (|Ω_de| < 0.25), not for the idea of Everpresent Λ in general.
+ 
+ | Implementation | Realizations | Parameters fixed to Planck (r_d, ω_m): pass A / B | Final criterion (c): pass A / B | 95 % upper limit, (c) |
+ |:--|--:|:--|:--|:--|
+ | Model 2 (Zwane et al. 2018; α = 0.8824, μ = 0.9804) | 20000 | 1 / 1 (fails the CMB proxy) | 0 / 0 | < 1.5×10⁻⁴ |
+ | Model 1 (Das et al. 2023), α = 0.01 | 300000 | 0 / 0 (best Δχ² = +17.3) | 0 / 0 | < 1.0×10⁻⁵ |
+ | Model 1, α = 0.005 | 200000 | 0 / 0 (best Δχ² = +30.5) | 0 / 0 | < 1.5×10⁻⁵ |
+ | Final criterion (c), all realizations | 520000 | — | **0 / 0** | — |
+ 
+ Conditions of criterion (c): r_d computed from each realization's own expansion history; small dark energy around recombination, $\lvert\Omega_{\rm de}\rvert < 0.25$ over $z\in[1000, 5000]$ (CMB proxy); $\omega_m$ within ±1.6 % (±2σ) of Planck 2018. Data: DESI DR2 BAO, 13 distances with full covariance. Pass A: χ² ≤ 29.84 (ΛCDM with Planck parameters); pass B: Δχ² < 4. Source: [`analisis_cmb/everpresent/EVERPRESENT_vs_DESI.md`](../analisis_cmb/everpresent/EVERPRESENT_vs_DESI.md) §3a, §3c and §3e.
+ 
```

**l.219** — punto 1. Everpresent como §4.2 y C10: §5.1: el párrafo de Everpresent se mueve a §4.2

```diff
- DESI constrains the **effective** residual amplitude after any amplification—not $\sigma_0$, and not a new value of $\rho_\Lambda$. Without a physical map $\sigma_0\to A_0$, “Euclid will see vacuum noise” is not a theorem of this analysis.
- 
- **Everpresent Λ.** The seed σ_0 ∼ 10⁻⁶¹ uses N ∼ 10¹²² from horizon-area (Bekenstein–Hawking) counting. With the four-volume counting of Everpresent Λ, δΛ tracks H² and the relative amplitude is of order one with no amplifier; the ∼10⁵⁶ gap is specific to the area-counting choice. Zwane et al. (2018) find that Everpresent Λ fits the CMB, BAO and SN data of that time about as well as ΛCDM and eases the low-redshift H0 and Lyα BAO tensions. The order-one alternative to the area-counting seed has been tested directly. In its two published implementations (Zwane et al. 2018, Model 2; Das et al. 2023, Model 1), computing r_d with each realization's own history, requiring small dark energy over z ∈ [1000, 5000] and ω_m within ±2σ of Planck, none of the realizations fits DESI DR2 BAO as well as ΛCDM, nor within Δχ² < 4: 0/20000 for Model 2 (< 1.5×10⁻⁴ at 95 %) and 0/300000 and 0/200000 for Model 1 with α = 0.01 and 0.005 (< 1.0×10⁻⁵, < 1.5×10⁻⁵). Without the CMB conditions, rare realizations pass, all with ω_m far from Planck or order-one dark energy at recombination. As small fluctuations around ΛCDM, DESI DR2 bounds the Zwane et al. amplitude to α < 1.85×10⁻² (μ fixed). These results hold for those phenomenological implementations, with the CMB only through a proxy (|Ω_de| < 0.25), not for the idea of Everpresent Λ in general.
+ DESI constrains the **effective** residual amplitude after any amplification—not $\sigma_0$, and not a new value of $\rho_\Lambda$. Without a physical map $\sigma_0\to A_0$, “Euclid will see vacuum noise” is not a theorem of this analysis.
```

**l.279** — punto 3. Slip: δ_m: Slip: cautela sobre δ_m (§6.2)

```diff
- > Slip is the correct *kind* of operator for anisotropic leakage past SDiff. With the BAO-bounded residual amplitude it is marginal, not starved: at ε = 1 it reaches $\lvert\gamma-1\rvert\approx 5\times 10^{-2}$ at $z=0.5$, level with the forecast floor ($\approx 1.0\times$) and $\approx 3.4\times$ below the current measurement. It is not a demonstrated shortcut: ε < 1, δ_m > 1 or z ≳ 1 push it below reach. The σ_X → slip map is an order-of-magnitude estimate and the best case requires ε = 1, δ_m = 1 and z ≈ 0.5.
+ > Slip is the correct *kind* of operator for anisotropic leakage past SDiff. With the BAO-bounded residual amplitude it is marginal, not starved: at ε = 1 it reaches $\lvert\gamma-1\rvert\approx 5\times 10^{-2}$ at $z=0.5$, level with the forecast floor ($\approx 1.0\times$) and $\approx 3.4\times$ below the current measurement. It is not a demonstrated shortcut. The σ_X → slip map is an order-of-magnitude estimate: it holds for δ_m ≈ 1 (quasi-linear scales) and does not extrapolate to δ_m ≪ 1; ε < 1 or z ≳ 1 push it below reach.
```

**l.291** — punto 1. Everpresent como §4.2 y C10: Limitación 5: solo el Sorkin seed + remisión a §4.2

```diff
- 5. **Sorkin seed:** used as motivational UV scale (horizon-area counting) for the amplification discussion, not as a fitted parameter. Everpresent Λ (Ahmed et al. 2004; Zwane et al. 2018) is a different model from the one tested here: Λ has zero mean and fluctuates, with an amplitude fixed by causal-set counting of the past four-volume (δΛ ∼ V^{-1/2}), at the level of the critical density in every epoch and with no free amplitude. The residual X of this work is a small stationary perturbation around a non-zero mean Λ, with free σ_X and θ. The bound on σ_X does not, by itself, test Everpresent Λ. Zwane et al. (2018) find that Everpresent Λ fits the CMB, BAO and SN data of that time about as well as ΛCDM and eases the low-redshift H0 and Lyα BAO tensions. Tested against DESI DR2 BAO with r_d computed for each realization, small dark energy around recombination and ω_m within ±2σ of Planck, none of the realizations of its two published implementations fits as well as ΛCDM (0/20000 for the model of Zwane et al.; 0/500000 for the model of Das et al., α = 0.01 and 0.005; `analisis_cmb/everpresent/EVERPRESENT_vs_DESI.md`). This applies to those phenomenological implementations, not to the idea of Everpresent Λ in general.
+ 5. **Sorkin seed:** used as motivational UV scale (horizon-area counting) for the amplification discussion, not as a fitted parameter. See §4.2 for Everpresent Λ.
```

**l.293** — punto 1. Everpresent como §4.2 y C10: Limitación 8 nueva: implementaciones fenomenológicas, CMB por proxy

```diff
- 7. **CMB check:** parameters other than H0 fixed to Planck 2018 (no MCMC); smoothing scale Δ introduced for the check; reading (a), δ ln ρ_DE = X/Ω_Λ0 (reading (b) would scale the amplitude by 1.46); halofit with tabulated w(a) in PPF not verified.
- 
+ 7. **CMB check:** parameters other than H0 fixed to Planck 2018 (no MCMC); smoothing scale Δ introduced for the check; reading (a), δ ln ρ_DE = X/Ω_Λ0 (reading (b) would scale the amplitude by 1.46); halofit with tabulated w(a) in PPF not verified.
+ 8. **Everpresent Λ test:** the results of §4.2 hold for the two published phenomenological implementations, with the CMB only through a proxy ($\lvert\Omega_{\rm de}\rvert < 0.25$ over $z\in[1000,5000]$), not for the idea of Everpresent Λ in general.
+ 
```

**l.339** — punto 6. Menores: §8.3 sin ./venv/bin/ y con nota sobre CAMB

```diff
- # CMB consistency (requires CAMB; see analisis_cmb/README.md)
- cd analisis_cmb && ./venv/bin/python realizaciones.py --halofit takahashi && ./venv/bin/python estadisticos.py resultados/takahashi/sigma_2.633e-02
+ # CMB consistency: requires CAMB (pip install camb) and the Planck files listed in analisis_cmb/README.md
+ cd analisis_cmb && python realizaciones.py --halofit takahashi && python estadisticos.py resultados/takahashi/sigma_2.633e-02
```

**l.347** — punto 6. Menores: Errata: espacio tras (95% CL)

```diff
- 1. **Null residual:** DESI DR2 BAO summary statistics do not favour stationary OU/QNM stochastic dark energy; $\sigma_X<2.6\times 10^{-2}$ (95% CL)under the stated kernel.
+ 1. **Null residual:** DESI DR2 BAO summary statistics do not favour stationary OU/QNM stochastic dark energy; $\sigma_X<2.6\times 10^{-2}$ (95% CL) under the stated kernel.
```

**l.349** — punto 1. Everpresent como §4.2 y C10: Conclusiones: punto nuevo tras el 2 (Everpresent) y renumeración

```diff
- 3. **CMB consistency:** at the DESI ceiling the OU residual is invisible to Planck 2018 (p95 Δχ² ≤ 0.067) and does not explain the H0 tension (σ(H0) ≈ 0.48 km/s/Mpc); a percent-level harmonic residual is equally allowed by DESI (58/60) and invisible in the CMB.
+ 3. **Everpresent Λ:** its two published implementations do not fit DESI DR2 BAO together with the CMB conditions (0/20000 and 0/500000 realizations; §4.2).
+ 4. **CMB consistency:** at the DESI ceiling the OU residual is invisible to Planck 2018 (p95 Δχ² ≤ 0.067) and does not explain the H0 tension (σ(H0) ≈ 0.48 km/s/Mpc); a percent-level harmonic residual is equally allowed by DESI (58/60) and invisible in the CMB.
```

**l.350** — punto 1. Everpresent como §4.2 y C10: Renumeración de conclusiones

```diff
- 4. **Amplification gap:** $\sigma_0\to 10^{-5}$ is $\sim 10^{56}$ (Euclid target); $\sigma_0\to 2.6\times 10^{-2}$ is $\sim 10^{59}$ (DESI ceiling); linear free amplifiers fail both.
- 5. **Slip honesty:** anisotropic leakage past SDiff is a real structural crack in principle and a marginal channel in practice — at best level with forecast floors ($\approx 1.0\times$ the Sakr et al. floor, $\approx 3.4\times$ below Maus et al., at $z=0.5$). The σ_X → slip map is an order-of-magnitude estimate: it holds for δ_m ≈ 1 (quasi-linear scales) and does not extrapolate to δ_m ≪ 1; ε < 1 or z ≳ 1 push it below reach.
+ 5. **Amplification gap:** $\sigma_0\to 10^{-5}$ is $\sim 10^{56}$ (Euclid target); $\sigma_0\to 2.6\times 10^{-2}$ is $\sim 10^{59}$ (DESI ceiling); linear free amplifiers fail both.
+ 6. **Slip honesty:** anisotropic leakage past SDiff is a real structural crack in principle and a marginal channel in practice — at best level with forecast floors ($\approx 1.0\times$ the Sakr et al. floor, $\approx 3.4\times$ below Maus et al., at $z=0.5$). The σ_X → slip map is an order-of-magnitude estimate: it holds for δ_m ≈ 1 (quasi-linear scales) and does not extrapolate to δ_m ≪ 1; ε < 1 or z ≳ 1 push it below reach.
```

**l.351** — punto 3. Slip: δ_m: Slip: cautela sobre δ_m (conclusiones)

```diff
- 5. **Slip honesty:** anisotropic leakage past SDiff is a real structural crack in principle and a marginal channel in practice — at best level with forecast floors ($\approx 1.0\times$ the Sakr et al. floor, $\approx 3.4\times$ below Maus et al., at $z=0.5$). The σ_X → slip map is an order-of-magnitude estimate and the best case requires ε = 1, δ_m = 1 and z ≈ 0.5.
- 6. **Scientific posture:** the strength of this programme is what it closes—limits, exclusions, and an honest map of open theory prices—not a promised detection.
+ 5. **Slip honesty:** anisotropic leakage past SDiff is a real structural crack in principle and a marginal channel in practice — at best level with forecast floors ($\approx 1.0\times$ the Sakr et al. floor, $\approx 3.4\times$ below Maus et al., at $z=0.5$). The σ_X → slip map is an order-of-magnitude estimate: it holds for δ_m ≈ 1 (quasi-linear scales) and does not extrapolate to δ_m ≪ 1; ε < 1 or z ≳ 1 push it below reach.
+ 7. **Scientific posture:** the strength of this programme is what it closes—limits, exclusions, and an honest map of open theory prices—not a promised detection.
```

**l.358** — punto 4. Versión v2: Historial de versiones antes de Acknowledgments

```diff
- ## Acknowledgments
- 
+ ## Version history
+ 
+ - **v1 (July 2026):** σ_X bound quoted as the grid node 2.5×10⁻²; derived numbers (amplification gaps, N_eff, slip, appendix bounds) computed with an earlier working target of 1.5×10⁻⁴; CPL fits with diagonal covariance.
+ - **v2 (October 2026):** bound = root of the profile, 2.6×10⁻² (2.633×10⁻²); all derived numbers regenerated; CPL fits with the official DESI covariance; CMB consistency and harmonic residual (§3.5) and Everpresent Λ (§4.2) added; slip map corrected (π_T = εσ_X ρ_crit,0).
+ - Audit and before/after records: [`analisis_cmb/AUDITORIA.md`](../analisis_cmb/AUDITORIA.md); blocks [A](../analisis_cmb/bloque_A/COMPARACION.md), [B1](../analisis_cmb/bloque_B/COMPARACION_B1.md), [B2](../analisis_cmb/bloque_B/COMPARACION_B2.md), [C](../analisis_cmb/bloque_C/COMPARACION_C.md), [D](../analisis_cmb/bloque_D/COMPARACION_D.md), [E](../analisis_cmb/bloque_E/COMPARACION_E.md).
+ 
+ ---
+ 
+ ## Acknowledgments
+ 
```

**l.382** — punto 4. Versión v2: Referencia 9 → October 2026

```diff
- 9. Morales Souhail, J., this repository and companion notes (July 2026).
+ 9. Morales Souhail, J., this repository and companion notes (October 2026).
```

### `manuscript/CLAIMS.md`

**l.23** — punto 3. Slip: δ_m: Slip: cautela sobre δ_m (C7)

```diff
- | **C7** | With ε = 1, the BAO-bounded residual gives $\|\gamma-1\|\lesssim 5\times 10^{-2}$ ($z=0.5$, $\delta_m=1$): ≈3.4× below Maus et al. and ≈1.0× the Sakr et al. constant-η forecast floor (order-of-magnitude map) | `scripts/slip_bridge.py` |
+ | **C7** | With ε = 1, the BAO-bounded residual gives $\|\gamma-1\|\lesssim 5\times 10^{-2}$ ($z=0.5$, $\delta_m=1$): ≈3.4× below Maus et al. and ≈1.0× the Sakr et al. constant-η forecast floor (order-of-magnitude map; holds for δ_m ≈ 1, quasi-linear scales, and does not extrapolate to δ_m ≪ 1; ε < 1 or z ≳ 1 push it below reach) | `scripts/slip_bridge.py` |
```

**l.25** — punto 1. Everpresent como §4.2 y C10: C10: Everpresent Λ como resultado propio (§4.2)

```diff
- | **C9** | A harmonic residual X = A sin(2π ln a/P + φ) with A ≤ 0.025, P ∈ [0.1, 2] is allowed by DESI DR2 BAO in 58/60 cases and is CMB-indistinguishable in 60/60 | `analisis_cmb/oscilador/resultados/tabla.md` |
- 
+ | **C9** | A harmonic residual X = A sin(2π ln a/P + φ) with A ≤ 0.025, P ∈ [0.1, 2] is allowed by DESI DR2 BAO in 58/60 cases and is CMB-indistinguishable in 60/60 | `analisis_cmb/oscilador/resultados/tabla.md` |
+ | **C10** | The two published implementations of Everpresent Λ do not fit DESI DR2 BAO together with the CMB: with r_d from each realization's own history, $\lvert\Omega_{\rm de}\rvert < 0.25$ over $z\in[1000,5000]$ and $\omega_m$ within ±2σ of Planck, 0/20000 (Model 2) and 0/500000 (Model 1, α = 0.01 and 0.005) realizations fit as well as ΛCDM | `analisis_cmb/everpresent/EVERPRESENT_vs_DESI.md` §3e; PREPRINT §4.2 |
+ 
```

**l.43** — punto 1. Everpresent como §4.2 y C10: N11 reescrita: la prueba directa es C10

```diff
- | **N11** | The $\sigma_X$ bound does not by itself test Everpresent Λ (Ahmed et al. 2004; Zwane et al. 2018), a zero-mean Λ with order-one fluctuations; that model is tested separately in `analisis_cmb/everpresent/` |
+ | **N11** | The $\sigma_X$ bound alone does not test Everpresent Λ; the direct test is C10 |
```

### `papers/anisotropic-slip-option0.md`

**l.25** — punto 3. Slip: δ_m: Slip: cautela sobre δ_m (preview)

```diff
- **Conclusion (preview).** The anisotropic “gap” is a real structural crack in SDiff protection, and $\eta$ is the right kind of observable. With $\sigma_X\lesssim 2.6\times 10^{-2}$ the best case (ε = 1, δ_m = 1) reaches $\lvert\gamma-1\rvert\approx 5\times 10^{-2}$ at $z=0.5$, level with Euclid-like forecasts: marginal, not starved. Its visibility depends on ε, δ_m and z rather than on amplification.
+ **Conclusion (preview).** The anisotropic “gap” is a real structural crack in SDiff protection, and $\eta$ is the right kind of observable. With $\sigma_X\lesssim 2.6\times 10^{-2}$ the best case (ε = 1, δ_m = 1) reaches $\lvert\gamma-1\rvert\approx 5\times 10^{-2}$ at $z=0.5$, level with Euclid-like forecasts: marginal, not starved. Its visibility depends on ε, δ_m and z rather than on amplification; the estimate holds for δ_m ≈ 1 (quasi-linear scales) and does not extrapolate to δ_m ≪ 1, and ε < 1 or z ≳ 1 push it below reach.
```

**l.160** — punto 3. Slip: δ_m: Slip: cautela sobre δ_m (variante ε ≲ 0.5 / z ≳ 1)

```diff
- **Interpretation.** The most optimistic anisotropic fraction does not make the gap visible at current DESI precision (≈3.4× short at z = 0.5); at Euclid-like forecast precision it is marginal (≈1.0× the Sakr et al. floor) for ε ≈ 1, δ_m ≈ 1, z ≈ 0.5, and out of reach for ε ≲ 0.5 or z ≳ 1. The map is order-of-magnitude.
+ **Interpretation.** The most optimistic anisotropic fraction does not make the gap visible at current DESI precision (≈3.4× short at z = 0.5); at Euclid-like forecast precision it is marginal (≈1.0× the Sakr et al. floor) for ε ≈ 1, δ_m ≈ 1, z ≈ 0.5, and out of reach for ε ≲ 0.5 or z ≳ 1. The map is order-of-magnitude: it holds for δ_m ≈ 1 (quasi-linear scales) and does not extrapolate to δ_m ≪ 1.
```

### `papers/resume.txt`

**l.16** — punto 3. Slip: δ_m: Slip: cautela sobre δ_m

```diff
-   Act IV bridge: scripts/slip_bridge.py  (|gamma-1|~5e-2 best case at z=0.5: ~1.0x forecast floor, ~3.4x below Maus; marginal, not starved)
+   Act IV bridge: scripts/slip_bridge.py  (|gamma-1|~5e-2 best case at z=0.5: ~1.0x forecast floor, ~3.4x below Maus; marginal, not starved; holds for delta_m ~ 1, not for delta_m << 1; eps < 1 or z >~ 1 push it below reach)
```

### `README.md`

**l.19** — punto 1. Everpresent como §4.2 y C10: Coherencia con C10

```diff
- 2. [`manuscript/CLAIMS.md`](manuscript/CLAIMS.md) — what I claim (C1–C9) and what I do not (N1–N13)  
+ 2. [`manuscript/CLAIMS.md`](manuscript/CLAIMS.md) — what I claim (C1–C10) and what I do not (N1–N13)  
```

### `manuscript/README.md`

**l.32** — punto 1. Everpresent como §4.2 y C10: Coherencia con §4.2

```diff
- - Tachyonic coherent growth exclusion
- - Amplification bottleneck $10^{56}$ + Routes 1–3 numerical verdicts
+ - Tachyonic coherent growth exclusion
+ - Everpresent Λ tested directly: both published implementations excluded with the CMB conditions (§4.2)
+ - Amplification bottleneck $10^{56}$ + Routes 1–3 numerical verdicts
```

**l.55** — punto 4. Versión v2: Estado v2

```diff
- **Preprint .** Independent research, July 2026.
+ **Preprint, v2 (October 2026; v1: July 2026).** Independent research. Version history at the end of `PREPRINT.md`.
```

### `manuscript/COVER_LETTER.md`

**l.11** — punto 4. Versión v2: Coherencia con v2: CMB, Everpresent y slip

```diff
- The main results are null and exclusionary: the data prefer smooth residual evolution ($\sigma_X < 2.6\times 10^{-2}$ (95% CL) under a stated OU kernel); a rank-1 coherent growing mode is strongly disfavoured ($\Delta\chi^2\approx +25$); and I show that linear free amplifiers of a Sorkin-type seed fail by many decades, with machine scans of three open cards documenting their theoretical price.
+ The main results are null and exclusionary: the data prefer smooth residual evolution ($\sigma_X < 2.6\times 10^{-2}$ (95% CL) under a stated OU kernel); a rank-1 coherent growing mode is strongly disfavoured ($\Delta\chi^2\approx +25$); the two published implementations of Everpresent Λ do not fit the same data together with the CMB; a Planck 2018 consistency check shows that the residual at the DESI ceiling is invisible in the CMB and does not explain the H0 tension; and I show that linear free amplifiers of a Sorkin-type seed fail by many decades, with machine scans of three open cards documenting their theoretical price. This is version 2 (October 2026): the bound is now the root of the profile likelihood, every derived number has been regenerated, and the changes from v1 are listed in the manuscript's version history.
```

## Punto 7: comprobación final

- **Numeración:**
  - PREPRINT: §1–§9 con §4.1 (tachyonic) y §4.2 (Everpresent), §3.1–§3.5, §5.1–§5.4, §6.1–§6.2. La guía de lectura §1.3, fila 4, cita las dos subsecciones.
  - Conclusiones 1–7, con la 3 nueva (Everpresent); limitaciones 1–8, con la 8 nueva.
  - CLAIMS: C1–C10 y N1–N13, consecutivos. README raíz: "C1–C10 / N1–N13".
- **Enlaces:** `python analisis_cmb/check_links.py`: 68 ficheros .md, 62 enlaces relativos, **0 rotos** (más el del snapshot histórico de siempre, `bloque_B/antes/…/VERDICT.md:7`).
- **results/:** solo aparece `results/ou_qnm_mle/` nuevo (`ou_qnm_mle.json`, `ou_qnm_mle.txt`); el resto de `results/` y `figures/` no cambia.
- **`git grep -n "old run" -- manuscript`:** una sola línea, `manuscript/PREPRINT.md:159` (la cita del run diagonal en §3.1).
