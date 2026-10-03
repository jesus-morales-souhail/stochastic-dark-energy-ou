# AUDITORÍA del repo stochastic-dark-energy-ou (2026-10-03)

Alcance: todos los ficheros versionados fuera de `analisis_cmb/` (README.md, manuscript/, papers/, notes/, results/, scripts/), más dos notas sobre `analisis_cmb/` al final del punto 3. **No se ha modificado ningún fichero.** Líneas comprobadas con `grep -n` en el estado de `main` = `f12527d`.

**Columna "tipo":**
- **T** = solo texto o etiqueta.
- **N** = cambia un número que aparece como resultado.
- **R** = fichero de resultados o figura que hay que **regenerar** con el script, no editar a mano.

**Comprobado directamente** (además de la lectura de líneas):
- La aritmética de los números derivados de 1.5e-4: G_DESI, N_eff, slip, H.5 y recuento route 3 en `results/amplification_routes/route3_avalanche_scan.csv` (p95 máx = 2.08e-3; 288/288 < 2.5e-2; 216/288 < 1.5e-4).
- Que `results/joint_w0wa_sigma/joint_w0wa_sigma.json` tiene logL idénticos bit a bit en `9e8b6bf` (cov diagonal) y `9e7c765` (etiqueta "full").
- Que `eos_cpl_summary.json` sale de `eos_efectiva.py` en `f36cc79`, que usa α hard-coded.
- Las referencias arXiv del punto 5.
- La página de Euclid del punto 4.

Lo demás procede de la lectura de líneas, sin recálculo. Se marca **(verificar)** donde la propuesta implica un número que no he recomputado.

Ω_Λ0 = 0.685 → factor entre "δΩ_Λ" y "δρ_Λ/ρ_Λ": 1/Ω_Λ0 ≈ 1.46. σ_X < 2.5e-2 equivale a rms(δρ_Λ/ρ̄_Λ) < 3.65e-2.

---

## 1. Definición de X

Referencia:
- `scripts/desi_dr2_data.py::add_ou_kernel` (l.224–233), con δα_i = S(z_i)·X(z_i) y S = ∂ln D_V/∂Ω_Λ en la dirección plana (Ω_m = 1 − Ω_Λ). Por tanto **X = δΩ_Λ**, y σ_X² = Var X (desviación típica estacionaria, = σ²/2θ con σ la difusión).
- Lectura (a) del análisis CMB: X = δρ_DE/ρ_crit,0, δ ln ρ_DE = X/Ω_Λ0.

Texto propuesto para unificar (en adelante **[X-std]**): *"X(x) ≡ δρ_DE(x)/ρ_crit,0 (= δΩ_Λ en unidades de hoy); δ ln ρ_DE = X/Ω_Λ0; σ_X² ≡ Var X = σ²/(2θ)."*

Definiciones: **D1** δΩ_Λ · **D2** δρ_DE/ρ_crit,0 sumado a E² (≡ D1 en unidades de hoy) · **D3** δρ_Λ/ρ̄_Λ = δ ln ρ_DE (lectura (b), ×1.46).

| fichero | línea | qué dice ahora | qué debería decir | por qué (fuente) | tipo |
|---|---|---|---|---|---|
| manuscript/PREPRINT.md | 82 | "X(x) ≡ δΩ_X(x)" | [X-std] | Ω_X no está definido (errata); kernel = D1 | T |
| manuscript/PREPRINT.md | 86–94 | Var(X)=σ²/2θ; "C^OU_ij = σ_X² e^{−θΔx}" (también QNM) | "C^OU_ij = S_i S_j σ_X² e^{−θΔx_ij}, σ_X² ≡ Var X" | Falta S_iS_j (add_ou_kernel l.232); σ_X no enlazado a σ | T |
| papers/stochastic-dark-energy-desi-dr2.md | 70 | "X(x) ≡ δΩ_Λ(x)" (D1) | [X-std] | Fija unidades y lectura (a) | T |
| papers/stochastic-dark-energy-desi-dr2.md | 82, 182, 264, 272 | Var/Cov con σ²/(2θ); σ_X sin definir | Añadir "σ_X² ≡ σ²/(2θ): σ_X es la desviación típica estacionaria, no la difusión σ" | Tablas (l.287) y código dan σ_X como desviación típica | T |
| papers/stochastic-dark-energy-desi-dr2.md | 164, 172 | ⟨ξξ⟩ = σ²H²ρ_Λ,0²/(2θ)·δ_D | ⟨ξ(t)ξ(t')⟩ = σ²H ρ_crit,0² δ_D(t−t') con X = δρ/ρ_crit,0 **(verificar)** | dW_x² = dx = H dt da un solo H; 1/(2θ) es la varianza estacionaria, no la del ruido | T |
| papers/stochastic-dark-energy-desi-dr2.md | 176 | "defining X = δρ_Λ/ρ_Λ,0" (D3) | "X = δρ_DE/ρ_crit,0" | Contradice l.70/160 por ×1.46 | T |
| papers/stochastic-dark-energy-desi-dr2.md | 404 | "σ_X ≡ δρ_Λ/ρ_Λ ≈ (V'/V)δφ" (D3; además confunde σ_X con la fluctuación) | "rms(δρ_Λ)/ρ_Λ = σ_X/Ω_Λ0 ≈ (V'/V) rms(δφ)"; propagar a l.412 y l.420 (3.65e-2) | Factor 1.46 en la cota que se deriva | N |
| papers/stochastic-dark-energy-desi-dr2.md | 562 | H² = H0²[Ω_m(1+z)³ + Ω_Λ + X(z)] (D2) | Añadir "X(z) ≡ δρ_DE/ρ_crit,0"; en lectura (a), X(z) − X(0) para que E(0)=1 | D2 = D1 en unidades de hoy; S(z) mueve también Ω_m, así que el χ(z,z') de C.4 no es exactamente el jacobiano del kernel | T |
| papers/stochastic-dark-energy-desi-dr2.md | 653 | "X(x,r) ≡ δρ_Λ/ρ_Λ,0" (D3) | "≡ δρ_DE(x,r)/ρ_crit,0" | Coherencia con l.70 | T |
| papers/stochastic-dark-energy-desi-dr2.md | 825 (+833) | "σ_X ≡ δρ_Λ/ρ_Λ", "defined in the main text" | "σ_X ≡ rms(δρ_DE)/ρ_crit,0; rms(δρ_Λ)/ρ_Λ = σ_X/Ω_Λ0"; en H.2, δφ ≃ (σ_X/Ω_Λ0)V/V' | D3, y el texto principal (l.70) no lo define así | N |
| papers/unimodular-gravity-vacuum-smoothness.md | 68, 72 | "X ≡ δρ_Λ/ρ̄_Λ"; σ_X² ≡ ⟨X²⟩ (D3) | "X ≡ δρ_DE/ρ_crit,0 = Ω_Λ0·δρ_Λ/ρ̄_Λ; σ_X² ≡ ⟨X²⟩ = σ²/2θ (⇔ rms δρ_Λ/ρ̄_Λ < 3.65e-2)" | ×1.46 | T |
| papers/fundamental-vs-emergent-vacuum-relaxation.md | 32 | "fractional fluctuations X ≡ δΩ_Λ" | "fluctuaciones X ≡ δΩ_Λ = δρ_DE/ρ_crit,0" (quitar "fractional") | "fractional" sugiere D3 | T |
| papers/scope-and-mixups.md | 12, 18, 23, 54 | X ≡ δΩ_Λ; S_iS_jσ_X²; "σ_X = σ/√(2θ)" (D1) | OK; opcional añadir "δ ln ρ_DE = X/Ω_Λ0" | Único sitio que separa σ de σ_X: modelo a seguir | T |
| papers/amplification-gap.md | 54, 56 | \|δρ\|max/(σ0ρ_Λ); "fractional residual amplitude σ_X" (D3 implícita) | "\|δρ\|max = σ_X ρ_crit,0 = (σ_X/Ω_Λ0)ρ_Λ"; el gap en densidad queda ≈1.85e59 **(verificar)** | ×1.46 | N |
| papers/anisotropic-slip-option0.md | 137, 141, 149–151 | π_T = εσ_X ρ_X, "fractional residual amplitude" (D3 implícita) | "π_T = εσ_X ρ_crit,0"; columna de la tabla ×1/Ω_Λ0 si ρ_X ≈ ρ_Λ0 | ×1.46; también scripts/slip_bridge.py l.45 y l.73 | N |
| papers/quantum-fluid-instabilities-desi-dr2.md | 103 vs 127/130 | σ_X² ≡ ⟨δρ_Λ²⟩/ρ_Λ² (D3), pero la verosimilitud usa u_i = σ_0 S_i e^{t/tc} (D1) | "σ_X² ≡ ⟨δρ_DE²⟩/ρ_crit,0²" | Incoherencia interna ×Ω_Λ0 (scripts/tachyonic_rank1_mle.py l.9/33 usa D1) | T |
| papers/principle-of-vacuum-smoothness.md | 15–19 | σ_X como cota sobre "fluctuations in ρ_Λ" (D3 implícita) | "σ_X = desviación típica de δρ_DE/ρ_crit,0 (⇔ rms δρ_Λ/ρ_Λ < 3.65e-2)" | ×1.46 | T |
| notes/desqueezing-relaxation-vacuum-fluctuations-note.md | 47 | "X ≡ δΩ_Λ/Ω̄_Λ (or an equivalent fractional density contrast…)" (D3) | [X-std], y quitar "or equivalent" | No son equivalentes (×1.46) | T |
| notes/desqueezing-relaxation-vacuum-fluctuations-note.md | 165, 173, 209 | "fractional fluctuations X(x)"; Var/Cov con σ²/(2θ) | "X(x) ≡ δρ_DE/ρ_crit,0"; añadir "σ_X² ≡ σ²/(2θ)" | La cota de l.43 es sobre σ_X | T |
| scripts/desi_dr2_data.py | 224–233 | docstring "OU residual signal", sin definir X | Docstring con [X-std] | Es la fuente de verdad y no lo documenta | T |
| scripts/ou_bao_likelihood.py; scripts/ou_bao_stochastic_test.py | 86; 83 | `SIGMA_OU = 2.31e-2 # sigma_X (calibrated…)` | Renombrar o comentar como desviación típica de la calibración DR1 obsoleta (ver punto 2) | El nombre invita a confundirlo con la difusión σ | T |
| scripts/eos_efectiva.py | 129, 150–156 | "sigma" de un toy de δw ~ σ²/(2θ)(…)/Ω_DE | Renombrar (p. ej. sigma_w) y decir que no es σ_X | Es otro parámetro; papers/neutrino-de-degeneracy-and-sigma-x.md l.30 los mezcla (punto 3) | T |
| papers/stochastic-dark-energy-desi-dr2.md; papers/euclid-protocol-vacuum-relaxation.md; scripts/cross_correlation_DESI.py; scripts/desqueezing/cosmological_mapping.py | 14, 36, 160, 214, 300–310, 689, 707; 43, 54; 215–281; 175, 280–290 | Usos coherentes con D1 o separación correcta σ/σ_X | OK, no cambiar | — | — |

## 2. Valores de σ_X

Categorías:
- **V** = 2.5e-2 vigente (`results/profile_sigma_x/profile_sigma_x.json`, profile sobre θ, cov gaussiana 13×13 proyectada a α).
- **O** = 1.5e-4, objetivo viejo. **O-cons** = número derivado de 1.5e-4 bajo una etiqueta de 2.5e-2.
- **A** = 2.175e-2 aproximado.
- **D** = 5e-5 / 4.5e-5 del run antiguo con covarianza diagonal.
- **C** = calibración antigua (0.018 / 0.023 / 2.31e-2).

**Nota:** la cov 7×7 "completa" es diagonal entre bins (`results/desi_cov/desi_cov_alpha_iso_7x7.txt`; el 13×13 público es diagonal por bloques). Sus σ difieren de las de figure6 en < 1.2 %. "Completa" frente a "diagonal" cambia poco los números, pero las etiquetas tienen que decir la verdad.

### 2a. Apariciones que hay que cambiar

| fichero | línea | qué dice ahora | qué debería decir | por qué (fuente) | tipo |
|---|---|---|---|---|---|
| manuscript/PREPRINT.md | 18–20 | "MLE drives… to the numerical floor. The working 95% … upper limit is σ_X<2.5e-2" (V mal descrita) | "profile likelihood over θ (ΔlnL ≥ −1.92), full meas. cov: σ_X < 2.5e-2" | La cota no sale de la MLE (profile_sigma_x_desi.py) | T |
| manuscript/PREPRINT.md | 28 | "full anisotropy… gives \|γ−1\|∼10^-4, of order 10^2–10^3 below…" (O-cons) | "\|γ−1\|∼3×10^-2 (∼5× bajo Maus, ∼1.6× bajo el suelo de forecast)" | 2·σ_X·0.644 con σ_X = 2.5e-2 → 0.032 (slip_bridge.py) | N |
| manuscript/PREPRINT.md | 76 | "The baseline pipeline adopts **diagonal** measurement errors…" | "official Gaussian BAO 13×13 projected to α (block-diagonal across bins)" | Contradice V y PREPRINT l.22 | T |
| manuscript/PREPRINT.md | 144–145 | Tabla MLE OU/QNM: σ_X = 5×10^-5 (D) | Etiquetar "(run antiguo, cov diagonal)" o sustituir por el profile (máximo en σ_X = 1e-6, borde del grid) | results/ou_bao_desi_dr2_run.txt es el run diagonal | N |
| manuscript/PREPRINT.md | 194 | "G_DESI = 2.5e-2/σ0 \| 1.27×10^59" (O-cons) | "2.1×10^59 ≈ 10^59" | 2.5e-2/1.18e-61 = 2.12e59; 1.27e57 era 1.5e-4/σ0 | N |
| manuscript/PREPRINT.md | 196 | "…incorrect by a factor ∼15… (r∼64–66)" (O-cons) | "factor ∼2.5×10^3; r∼64–68" **(verificar r)** | Mismo origen | N |
| manuscript/PREPRINT.md | 214 | "N_eff ∼ 4.44×10^7 (DESI)" (O-cons) | "N_eff ∼ 1.6×10^3 (DESI)" | 1/(2.5e-2)² = 1600; coincide con results/amplification_routes/VERDICT.md l.22 | N |
| manuscript/PREPRINT.md | 216 | "BAO-safe if σ ≲ 4×10^-5" (O-cons) | "Con umbral 2.5e-2, los 288 jobs son BAO-safe (p95 máx 2.1e-3)" | route3_avalanche_scan.csv | N |
| manuscript/PREPRINT.md | 242, 246–248 (+254) | "With σ_X=2.5e-2, ε=1": 1.93e-4 / 8.16e-5 / 4.18e-5; 880× / 2e3 / 4e3; 260× / 610× / 1.2e3 (V en la etiqueta + O-cons) | 3.2e-2 / 1.36e-2 / 6.95e-3; ∼5× / 12× / 24× vs Maus; ∼1.6× / 3.7× / 7.2× vs 0.05 **(verificar las dos últimas)**; revisar la conclusión de l.254 | 1.93e-4 = 2·1.5e-4·0.644 | N |
| manuscript/PREPRINT.md | 262 | "uses published BAO α and diagonal errors. Full DESI covariance matrices may change…" | "baseline usa ya la cov gaussiana oficial proyectada a α" | Contradice V | T |
| manuscript/CLAIMS.md | 22 | "C7 … ε=1 … gives \|γ−1\|∼10^-4" (O-cons) | "\|γ−1\| ≲ 3×10^-2 (z=0.5, ε=1, δ_m=1)" | Igual que PREPRINT l.28 | N |
| papers/HONEST_HEADLINES.md | 27 | "soft r∼64 (Euclid target) / ∼66 (DESI ceiling)" (O-cons) | "/ ∼68 (DESI ceiling)" **(verificar)** | Mismo origen que G_DESI | N |
| papers/amplification-gap.md | 43 | "The answers differ by about **one decade**" (O-cons) | "about three decades (10^56 vs 10^59)" | 10^57 → 10^59 | N |
| papers/amplification-gap.md | 48 | "G_DESI \| 1.27×10^59 \| r≈65.7" (O-cons) | "2.1×10^59 \| r≈68.0" **(verificar r)** | 2.5e-2/1.18e-61 | N |
| papers/amplification-gap.md | 191 | "DESI 95% ceiling \| 2.5e-2 \| 4.44×10^7 \| smaller by ∼10^114" | "1.6×10^3 \| ∼10^119" | 1/(2.5e-2)² | N |
| papers/amplification-gap.md | 199 | "DESI N_eff∼4.4×10^7 \| ∼10^-31 m \| ∼10^-12 m \| ∼mm" | "N_eff∼1.6×10^3" y longitudes recalculadas **(verificar)** | Mismo origen | N |
| papers/amplification-gap.md | 227–232 | "BAO-safe (p95<2.5e-2) 216/288"; "unsafe 72/288… σ≳1.2e-4"; "BAO-safe rule… σ≲4×10^-5" | "288/288 safe; 0/288 unsafe; todo el grid (σ ≤ 1e-3) es BAO-safe"; borrar l.232 | route3_avalanche_scan.csv: 288 con p95 < 2.5e-2 (216 eran con < 1.5e-4) | N |
| papers/anisotropic-slip-option0.md | 25, 118 | "with σ_X∼10^-4 it inherits…"; "the ∼10^-4 scale of the f=1, σ_X∼10^-4 estimate" (O) | "σ_X ≲ 2.5e-2 ⇒ \|γ−1\| ≲ 3e-2" | Cota vigente | N |
| papers/anisotropic-slip-option0.md | 145, 149–151, 155, 158, 164–168, 177 | Tabla slip con 1.93e-4…; "vs ∼10^-4"; "∼10^2–10^3"; ε=1/0.1/0.01 → 1.93e-4/1.93e-5/1.93e-6; "\|η−1\|∼10^-4" | 3.2e-2 / 1.36e-2 / 6.95e-3; "vs ∼3e-2"; "∼1.6–7×"; 3.2e-2 / 3.2e-3 / 3.2e-4; "∼10^-2" | Igual que PREPRINT l.246 | N |
| papers/anisotropic-slip-option0.md | 193 | "use 2.5e-2 and 5e-5 as bookends from resume.txt" | "2.5e-2 y σ_X→0 (best fit del profile)" | 5e-5 es el run diagonal | T |
| papers/data-pack-option0-internet.md | 83–94 | "Using σ_X∼10^-4… ∼2×10^-4"; "~500× larger"; "~250× larger" (O) | Recalcular con 2.5e-2 **(verificar)** | Cota vigente | N |
| papers/euclid-protocol-vacuum-relaxation.md | 81 | "σ_X log-uniform on [1e-6, 1e-2]… Spans Euclid targets and the DESI limit" | "[1e-6, 5e-2]" o quitar "and the DESI limit" | 1e-2 < 2.5e-2 | T |
| papers/euclid-protocol-vacuum-relaxation.md | 148 | "s=5e-5 \| A0=2.5e-2 \| min θ ∼0.4" (θ calculado con A0=1.5e-4) | "min θ ≈ 2×10^-3" **(verificar)** | O-cons | N |
| papers/resume.txt | 13, 15, 16 | "N_eff ~ 4.4e7 (DESI)"; "BAO-safe if sigma ≲ 4e-5"; "\|gamma-1\|~1e-4 best case; O(1e2-1e3)x starved" | "N_eff ~ 1.6e3"; "todo el grid BAO-safe"; "~3e-2; ~1.6–7x starved" | O-cons | N |
| papers/resume.txt | 41, 47 | "σ_X = 5 × 10⁻⁵" (D); además l.25 dice cov completa y l.33 "diagonal covariance only" | "(run antiguo, cov diagonal: results/ou_bao_desi_dr2_run.txt)" o sustituir por el profile | Run diagonal | N |
| papers/resume.txt | 55, 83–84 | "95% CL upper limit: σ_X < 2.5e-2"; "When marginalizing over the background… **strict upper limit**" | "(profile over θ, full cov, background fijo α=1); el null con background libre está en joint_w0wa_sigma (cov diagonal)" | El profile no marginaliza el background | T |
| papers/stochastic-dark-energy-desi-dr2.md | 14, 16–20, 368 | "MLE… assuming the CPL background is fixed to the best-fit values"; "MLE drives… upper limit"; "conservative 95% CL… fixed CPL background" | "profile likelihood over θ, full meas. cov, background fijo a α=1 (ΛCDM fiducial)" | profile_sigma_x_desi.py; el Apéndice G (l.920) ya lo admite | T |
| papers/stochastic-dark-energy-desi-dr2.md | 278, 430, 932 | "standard diagonal covariance from measurement errors"; "diagonal covariance approximation"; "diagonal DESI covariance" | "la cota usa la cov completa; §4.2 y el Apéndice G son runs diagonales" | Contradice V | T |
| papers/stochastic-dark-energy-desi-dr2.md | 287–288 | "H0/H1 free MLE \| 5×10^-5" (D) | "(cov diagonal, run antiguo)" o actualizar | Run diagonal | N |
| papers/stochastic-dark-energy-desi-dr2.md | 370 | "σ_X<2.5e-2 is a conservative estimate based on the numerical floor (5e-5), multiplied by a factor of 3… profile… deferred to Euclid DR1" | "Formal profile scan (scripts/profile_sigma_x_desi.py): 2.506e-2, ΔlnL ≥ −1.92" | Descripción falsa (5e-5×3 = 1.5e-4); el profile ya existe | T |
| papers/stochastic-dark-energy-desi-dr2.md | 372 | "…suppresses local fluctuations by more than four orders of magnitude" (O-cons) | "by ∼1.6 orders of magnitude" **(verificar)** | Con 2.5e-2 | N |
| papers/stochastic-dark-energy-desi-dr2.md | 440, 926 | "σ_X ≈ 4.5×10^-5 (numerical floor)" en el joint (D) | "σ_X → 4.2e-10 (results/joint_w0wa_sigma, cov diagonal)" | No coincide con el resultado archivado | N |
| papers/stochastic-dark-energy-desi-dr2.md | 476, 621–625 | "m_φ ≲ 10^-5 eV"; "σ_X<2.5e-2 implies m_φ ≲ 9.45×10^-5 eV" | Recalcular con 2.5e-2 (∼1.2×10^-3 eV si σ ∝ m²) **(verificar)** | Probable O-cons | N |
| papers/stochastic-dark-energy-desi-dr2.md | 869, 883, 887, 895, 906 | "2.7e-6/2.5e-2 ≈ 1.8×10^-2"; "≲1.8e-2"; "λ ≳ 55.5"; "\|ζ\|≲1.8e-2"; "\|ζ\|≲10^-2" | "≈1.1×10^-4"; "≲1.1e-4"; "λ ≳ 9.3e3" **(verificar)**; "\|ζ\|≲1.1e-4"; "≲10^-4" | Aritmética: 2.7e-6/2.5e-2 = 1.08e-4 (1.8e-2 = 2.7e-6/1.5e-4) | N |
| papers/stochastic-dark-energy-desi-dr2.md | 924 | "σ_X ≈ 9.2×10^-3 — nearly 200 times larger than the reported limit" | "≈180× el floor 5e-5 del run diagonal, pero por debajo del 95% CL 2.5e-2" | Compara con D, no con V | T |
| papers/unimodular-gravity-vacuum-smoothness.md | 14, 36–44 | "strict upper limit"; "MLE drives the amplitude to zero… conservative upper limit" | "profile 95% CL upper limit (over θ)" | V mal descrita | T |
| papers/unimodular-gravity-vacuum-smoothness.md | 130 | "m_φ ≲ 10^-5 eV, \|V'/V\| M_Pl ≲ 10^-2" | Recalcular con 2.5e-2 (igual que l.621 y H.5 del paper principal) **(verificar)** | O-cons | N |
| papers/unimodular-gravity-vacuum-smoothness.md | 300 | "the observational constraint that σ_X is below 10^-4" (O) | "below 2.5×10^-2 (95% CL)" | Cota vigente | N |
| notes/desqueezing-relaxation-vacuum-fluctuations-note.md | 234 | "MLE floor θ∼1e-3, σ_X∼5×10^-5" (D) | Añadir "(run antiguo, cov diagonal)" | Run diagonal | T |
| notes/mapping_tables.md | 24–38 (pares) | target "DESI_limit_2.5e-2" pero A0 = 1.501e-4 … 3.000e-4 (O-cons) | Regenerar desde scripts/desqueezing/cosmological_mapping.py (ya usa 2.5e-2) | Valores calculados con 1.5e-4 | R |
| results/ou_bao_desi_dr2_run.txt | 19–20 | "θ = 0.0010, σ_X = 0.00005" (D), sin marca | Cabecera "cov diagonal; superseded por results/profile_sigma_x/" | Run diagonal | R |
| results/desi_dr2_real_bao/desi_dr2_real_bao.json | 39 | "sigma_X_95_profile_upper": 0.02175… (A), la clave no dice approx | Clave "…_approx" o nota (el .txt l.8 sí dice "approx") | Aproximado | R |
| results/euclid_protocol/euclid_forecast_grid.csv | 68–78; 90–100 | A0 = 0.00015 (O); A0 = 0.018 (C) | Regenerar (el script usa 2.5e-2); etiquetar o retirar 0.018 | Sin regenerar | R |
| results/cosmological_mapping/table2_A0_required.csv | 2–16 (pares) | "DESI_limit_1.5e-4, 0.00015, …" (O) | Regenerar | El script ya dice 2.5e-2 | R |
| results/cosmological_mapping/mapping_summary.txt | 11 | "MLE floor: theta ~ 1e-3, sigma_X ~ 5e-5" (D) | "(run antiguo, cov diagonal)" | Run diagonal | R |
| results/sdiff_discrimination/summary.txt | 5 | "DESI lim = 1.5e-04" (O); l.25–29 ya dicen 2.5e-2 | Regenerar | Sin regenerar | R |
| results/amplification_routes/VERDICT.md | 26, 55–58 | "∼10^112–10^114"; "216/288 safe, 72 unsafe, σ≳1.2e-4"; "σ≲4×10^-5" | "∼10^112–10^119"; "288/288 safe, 0 unsafe" | O-cons; l.22 ya dice 1.6e3 (correcto) | R |
| results/amplification_routes/route1_N_eff_required.csv; scripts/amplification/route1_local_causal_set_seed.py | 3; 42 | Etiqueta "desi_1p5e-4" con valor 0.025 | Etiqueta "desi_2p5e-2" | Etiqueta O con valor V | T |
| scripts/gap_two_targets.py | 7 | Docstring "G_DESI = 2.5e-2/sigma_0 ~ 10^57" | "~ 10^59" | Los prints de l.55–57 ya son correctos | T |
| scripts/profile_sigma_x_desi.py (+ figures/profile_sigma_x.png) | 128 | axvline en 2.5e-2 con label "95% CL profile claim $1.5\times10^{-4}$" | label "$2.5\times10^{-2}$"; regenerar el PNG | Etiqueta O | T/R |
| scripts/desqueezing/euclid_protocol_forecasts.py | 36, 66 | A0s incluye 0.018 sin marca; label "$\sigma=1.5\times10^{-4}$" sobre SIGMA_DESI=2.5e-2 | Comentar que 0.018 es la calibración DR1 obsoleta; label "2.5e-2" | C y O | T |
| scripts/desqueezing/cosmological_mapping.py | 39, 137 | `SIGMA_X_MLE_FLOOR = 5.0e-5` | Comentario "(run antiguo, cov diagonal)" | D | T |
| scripts/ou_bao_likelihood.py; scripts/ou_bao_stochastic_test.py | 86, 583; 83 | `SIGMA_OU = 2.31e-2` "calibrated…"; scatter (0.001, 5e-5) "MLE best fit (DR2)" | "# SUPERSEDED DR1 calibration, NOT the 95% CL 2.5e-2"; label "(old diagonal-cov run)" | C (se confunde fácilmente con 2.5e-2) y D | T |
| scripts/cross_correlation_DESI.py | 32, 280, 650 | `SIGMA_X = 0.018`; "σ_X calibrado del paper" | "ilustrativo, calibración DR1 obsoleta" | C | T |

### 2b. Apariciones correctas (no cambiar)

Vigente 2.5e-2, bien descrita:
- README.md 30–35
- manuscript/CLAIMS.md 17, 20
- manuscript/COVER_LETTER.md 11
- manuscript/PREPRINT.md 22, 26, 151, 191, 222, 319, 321
- manuscript/README.md 31
- notes/desqueezing-… 31, 39–43, 233, 327, 335, 393, 474
- papers/BOUND_CORRECTION_SIGMA_X.md (todo)
- papers/HONEST_HEADLINES.md 18, 25, 36, 42
- papers/amplification-gap.md 37, 53–54, 114
- papers/anisotropic-slip-option0.md 14
- papers/data-pack-option0-internet.md 13
- papers/euclid-protocol-… 11, 121, 123, 247
- papers/fundamental-vs-emergent-… 14, 86, 140–142
- papers/principle-of-vacuum-smoothness.md 19
- papers/quantum-fluid-instabilities-desi-dr2.md 12, 18, 158, 190
- papers/scope-and-mixups.md 27–58
- papers/sdiff-fundamental-vs-emergent.md 19, 72, 95
- papers/stochastic-dark-energy-desi-dr2.md 224, 376, 416–420, 470, 478, 617, 865, 908
- papers/unimodular-… 76, 82, 296
- results/profile_sigma_x/* (fuente de verdad)
- results/cosmological_mapping/mapping_summary.txt 10
- results/sdiff_discrimination/summary.txt 25–29
- results/amplification_routes/VERDICT.md 22
- scripts: amplifier_audit.py 20; slip_bridge.py 45, 141; profile_sigma_x_desi.py 109, 113, 146; sdiff_discrimination.py; cosmological_mapping.py 38, 136, 209, 267, 269; desqueezing_relax_time.py 311; gap_two_targets.py 23, 56

Marcados ya como obsoletos:
- 1.5e-4: CLAIMS 17, BOUND_CORRECTION, scope-and-mixups 35, profile json 302
- 0.018: stochastic 220, sdiff grid "old_calib"
- 0.023: notes 235
- 2.31e-2: cosmological_mapping.py 43
- 2.175e-2: real_bao.txt 8

Otros, correctos: best fit ~0 (profile, joint, eos); mocks de Euclid en results/euclid_mcmc/ y euclid_joint_mcmc/ (marcados RETIRED en el .txt pero no en el .json); objetivo Euclid 1e-5; ventana 1e-4; σ_free ≈ 8.5e-5.

## 3. CPL (w0, wa de ajustes del repo)

Hechos comprobados:
- **(i) `results/eos_cpl_desi_dr2/eos_cpl_summary.{json,txt}`** (w0 = −0.98990, wa = −0.01590) se generó con `eos_efectiva.py` en `f36cc79` (2026-07-16). Esa versión tenía α hard-coded en l.27–29 (1.0030, 0.9947, … = `_FALLBACK_A` de desi_dr2_data.py l.38), **no** los α de figure6 de Zenodo (1.00500, 1.00728, …), y χ² diagonal. No se regeneró después.
- **(ii) `results/joint_w0wa_sigma/joint_w0wa_sigma.json`** (w0 = −1.00165, wa = +0.01725) se ejecutó con **cov diagonal** y α de figure6.
  - Su `"data": "…full meas. cov…"` está escrito a mano en joint_w0wa_sigma_desi.py l.263 (y en el título de la figura, l.300).
  - `cov_mode` (l.264) sí sale de `load_alpha_dv`, que arranca en "diagonal_figure6_errors" (desi_dr2_data.py l.168/182–183) y solo pasa a "full…" si se cumplen las condiciones de l.184–194.
  - Los logL del json son idénticos bit a bit a los de la ejecución diagonal `9e8b6bf`: el commit `9e7c765` solo cambió la etiqueta.

Los dos CPL discrepan (−0.99/−0.016 frente a −1.002/+0.017) por el vector α, no por la covarianza.

| fichero | línea | qué dice ahora | qué debería decir | por qué (fuente) | tipo |
|---|---|---|---|---|---|
| results/eos_cpl_desi_dr2/eos_cpl_summary.json, .txt | json 2, 5–6; txt 3, 7–8, 43 | "7 bins, diagonal errors"; −0.98990 / −0.01590 | Añadir "alpha = tabla fallback hard-coded (no figure6)", o regenerar con el script actual | (i) | R |
| scripts/eos_efectiva.py | 5, 33, 375, 404, 414 vs 443 | "full measurement covariance when available" frente a "diagonal approximation" hard-coded en l.443 | Usar `cov_mode` del loader en todas las etiquetas | Contradicción interna | T |
| results/joint_w0wa_sigma/joint_w0wa_sigma.json, .txt; figures/joint_w0wa_sigma.png | json 2–3; txt 4, 6; título | "full meas. cov" + "cov_mode": "diagonal_figure6_errors"; −1.0017 / +0.0172 | "data: …diagonal figure6 errors…", o regenerar con la cov proyectada | (ii) | R |
| scripts/joint_w0wa_sigma_desi.py | 263, 300 | Etiqueta "full meas. cov" fija | Tomar la etiqueta de `cov_mode` | (ii) | T |
| papers/resume.txt | 25 vs 28–29, 33 | "cov = Gaussian BAO 13x13 projected"; w0≈−0.99, wa≈−0.02; "ΛCDM (diagonal covariance only)" | "w0≈−0.99, wa≈−0.016 (eos_efectiva, α fallback, cov diagonal)", o citar el joint (−1.002/+0.017, diagonal) | (i); l.25 contradice l.33 | N |
| manuscript/PREPRINT.md | 76, 78 (y 22) | l.78 "CPL background near ΛCDM (w0≈−0.99, wa≈−0.02) as the smooth reference for α=1"; l.22 dice cov proyectada | Indicar fuente y cov; además α=1 es el fiducial ΛCDM (w0=−1, wa=0), no el ajuste CPL | (i) | T |
| papers/stochastic-dark-energy-desi-dr2.md | 440 (F5) | "joint… w0≈−0.99, wa≈−0.02… σ_X≈4.5e-5" | "w0=−1.002, wa=+0.017, σ_X→0 (4e-10), ΔAIC=+4, cov diagonal (results/joint_w0wa_sigma)" | Mezcla eos (w0, wa) con el run diagonal (σ_X); no es reproducible | N |
| papers/stochastic-dark-energy-desi-dr2.md | 922, 926 | (a) lnL=27.013; (c) w0≈−0.99, wa≈−0.02, σ_X≈4.5e-5, ΔlnL=+0.13 | Números del json archivado (lnL ΛCDM 24.8955; ΔlnL joint−cpl ≈ 7e-13) o regenerar | Run no archivado | N |
| papers/stochastic-dark-energy-desi-dr2.md | 430, 932 | "diagonal covariance approximation"; "S(z) fixed; diagonal DESI covariance" | Mantener (es cierto para los CPL archivados), o actualizar si se regenera | Choca con PREPRINT l.22 y resume l.25 | T |
| papers/neutrino-de-degeneracy-and-sigma-x.md | 30 | "{w0,wa,θ,σ_X}… MLE (w0,wa)≈(−0.99,−0.016) and σ_X→0" | "(w0,wa)≈(−1.002,+0.017), σ_X≈4e-10 (joint_w0wa_sigma, cov diagonal)" | Atribuye a un ajuste joint OU el (w0, wa) de eos, cuyo "σ" es otro parámetro | N |
| papers/covariance-desi-dr2-full.md | 68, 72 | joint y eos en "Scripts updated to use C_α" | Nota: "resultados archivados aún diagonales / α fallback; regenerar" | (i), (ii) | T |
| scripts/verify_bao_pipelines.py | 41–56 | Comprueba w0 ≈ −0.99 ± 0.05 en el json de eos | Sin cambio; la tolerancia no distingue los dos CPL | — | — |
| papers/stochastic-dark-energy-desi-dr2.md | 12, 924 | w0=−0.87±0.05, wa=−0.41±0.28; −0.785/−0.43 | Son externos (DESI), no del repo; no contrastados con los papers de DESI **(verificar)** | — | — |
| **analisis_cmb/NOTAS_ANALISIS_CMB.md, analisis_cmb/RESULTADO.md** | NOTAS "Curvas de referencia w0-wa", "Inconsistencias…"; RESULTADO §5 | "CPL repo… errores diagonales" | Añadir "y vector α de fallback hard-coded (no figure6)" | (i). Afecta solo a la curva de referencia secundaria; la cota se definió frente a ΛCDM | T |

## 4. Euclid

Fuente: cosmos.esa.int/web/euclid/dr1-timeline, "Update published on 15 June 2026" (página descargada y leída el 2026-10-03). **DR1-Foundation en noviembre de 2026**: datos brutos, imágenes calibradas, catálogos y espectros, ~1900 deg². **DR1 completo**, con los productos de alto nivel de clustering de galaxias y lensing débil, **a mediados de 2027**. DR1-Foundation no trae vector BAO.

Texto propuesto (**[E]**): *"Euclid DR1 completo con productos de cosmología (mediados de 2027; el DR1-Foundation de noviembre de 2026 no incluye LE3)"*.

| fichero | línea | qué dice ahora | qué debería decir | por qué (fuente) | tipo |
|---|---|---|---|---|---|
| papers/stochastic-dark-energy-desi-dr2.md | 446 | "## 9. Near-Term Observational Program: Euclid DR1" | "…Euclid DR1 completo (mediados de 2027)" | Página ESA | T |
| papers/stochastic-dark-energy-desi-dr2.md | 448 | "Euclid Data Release 1 (expected H2 2026) will provide >20 redshift bins… pipeline is fully ready" | [E]; quitar ">20 bins" o citar una fuente | Fecha errónea; el número de bins no tiene fuente en el repo | T |
| papers/stochastic-dark-energy-desi-dr2.md | 442 | "…ideally, the >20 bins of Euclid DR1" | "…el vector BAO del DR1 completo (mediados de 2027), con el número de bins que se publique" | Supuesto sin fuente | T |
| papers/stochastic-dark-energy-desi-dr2.md | 450 (+609) | "narrower redshift baseline (z∈[0.9,1.8])… ω_R,min≈16.2… >20 bins" | Marcar rango en z y bins como supuestos | papers/sensitivity_kernel_table.md l.37 da Δx≈0.5 y ω_R,min≈12.6 para "Euclid DR1": los dos ficheros no coinciden | N |
| papers/stochastic-dark-energy-desi-dr2.md | 452 | "Upon release of Euclid DR1, I will:" | "Cuando se publique el DR1 completo con productos LE3 (mediados de 2027)…" | DR1-Foundation no permite el paso 1 (l.453) | T |
| papers/stochastic-dark-energy-desi-dr2.md | 480, 611, 761 | "falsifiable prediction for Euclid DR1"; "Euclid DR1 geometrically cannot"; "**Euclid DR1.** If a residual…" | Basta con la aclaración [E] una vez en §9 | — | T |
| papers/resume.txt | 88 | "Euclid DR1 (expected H2 2026) will provide >20 BAO bins… σ_X > 10⁻⁵" | [E]; presentar ">20 bins" y "10⁻⁵" como supuesto/forecast | Fecha; 10⁻⁵ es el objetivo aspiracional (amplification-gap l.38) | T |
| results/ou_bao_desi_dr2_run.txt | 30, 50–51 | "Decisive test: >20 bins (Euclid DR1, expected H2 2026)"; "20+ bins" | Regenerar o anotar (los scripts actuales, ou_bao_likelihood.py l.415 y ou_bao_stochastic_test.py l.407, ya no ponen fecha) | Salida antigua | R |
| papers/unimodular-gravity-vacuum-smoothness.md | 280, 286, 308 (22, 290 sin fecha) | "§8… The Test of Euclid DR1"; "Euclid DR1 BAO (2026)"; "Euclid DR1 (expected 2026)" | [E] | Fecha | T |
| papers/sdiff-fundamental-vs-emergent.md | 83 (82, 87, 110 sin fecha) | "Euclid DR1 (2026)" | "Euclid DR1 completo (mediados de 2027)" | Fecha | T |
| papers/euclid-protocol-vacuum-relaxation.md | 100, 150, 159–160, 260, 269 | "≳20 bins when available"; "narrower Euclid DR1 path (Δx∼0.39)"; "Euclid-like 20/40"; "When real Euclid BAO summary statistics are public…"; ref. "[4] Euclid… DR1 BAO forecasts (as available…)" | Añadir [E]; marcar bins y Δx como supuestos de diseño; citar la página ESA | Fecha y supuestos | T |
| papers/fundamental-vs-emergent-vacuum-relaxation.md; results/sdiff_discrimination/summary.txt | 151; 34 | "Euclid null σ_X≪10^-5 with ≳20 bins" | "…con el vector BAO del DR1 completo (mediados de 2027)" | Supuesto de bins | T / R |
| notes/desqueezing-relaxation-vacuum-fluctuations-note.md | 374, 385 | "## 7. What Euclid DR1 can and cannot see"; "once N_bins≳20" | [E] | Fecha y supuesto de bins | T |
| papers/principle-of-vacuum-smoothness.md | 50 | "(for example Euclid DR1 or DESI Year 3+)" | Opcional [E] | — | T |

## 5. Everpresent Λ

**¿Se citan? No.** Cero coincidencias de Ahmed, Dodelson, Greene, Zwane, Afshordi, everpresent, 103523, astro-ph/0209274 y 1703.06265 en los ficheros versionados.

Referencias comprobadas en arXiv (2026-10-03):
- Ahmed, Dodelson, Greene & Sorkin, "Everpresent Λ", Phys. Rev. D 69, 103523 (2004), arXiv:astro-ph/0209274.
- Zwane, Afshordi & Sorkin, "Cosmological Tests of Everpresent Λ", Class. Quantum Grav., DOI 10.1088/1361-6382/aadc36, arXiv:1703.06265. El volumen y la página (35, 194002) que diste no aparecen en la ficha de arXiv: **verificar**.

**Hallazgo colateral:** la única referencia de Sorkin del repo está **mal emparejada**. "Sorkin, R. D., 'Is the Cosmological Constant a Nonlocal Quantum Residual?', arXiv:gr-qc/0503057 (2005)" aparece en stochastic-dark-energy-desi-dr2.md l.944 y quantum-fluid-instabilities-desi-dr2.md l.194; en fundamental-vs-emergent-… l.194 y notes/desqueezing-… l.449 sin título.
- gr-qc/0503057 es "Big extra dimensions make Λ too small" (Braz. J. Phys. 35, 280, 2005).
- El título citado corresponde a arXiv:0710.1675, "Is the cosmological 'constant' a nonlocal quantum residue of discreteness of the causal set type?" (AIP Conf. Proc. 957, 142, 2007).

Diferencias a enunciar, según los resúmenes de arXiv:
- **(D1) Objeto distinto.**
  - Everpresent Λ: Λ de media nula, que oscila entre valores positivos y negativos con magnitud ~ρ_crit en cada época. Amplitud fijada por el conteo de causal sets (δΛ ~ 1/√V, V = 4-volumen del pasado causal), sin amplitud libre.
  - Repo: residuo pequeño X = δΩ_Λ (σ_X < 2.5e-2) alrededor de una Λ media no nula; OU estacionario en ln a con σ_X y θ libres.
  - Por eso la cota σ_X **no excluye** Everpresent Λ.
- **(D2) Conclusión opuesta desde el mismo punto de partida.** Ahmed et al. usan la gravedad unimodular más causal sets para **predecir** fluctuaciones de Λ. El paper unimodular del repo concluye σ_X = 0.
- **(D3) Escala de la semilla.** La "semilla Sorkin ~10⁻⁶¹" del repo usa N ~ 10¹²² (conteo de área). Con el conteo de 4-volumen de Everpresent, δΛ sigue a H² y la amplitud relativa es O(1) sin amplificador. Incoherencia interna en el Axioma A2: stochastic-dark-energy-desi-dr2.md l.58 define N = V/L_P⁴ (4-volumen), pero l.66 usa N ~ 10¹²² (l.54: área).
- **(D4) Datos.** Zwane et al. encuentran, por MCMC, que Everpresent Λ ajusta las observaciones tan bien como ΛCDM y alivia las tensiones de H0 a bajo z y del BAO de Lyα a z~2–3. No ayuda con los agujeros negros ultramasivos ni con el litio.

| fichero | línea | qué dice ahora | qué debería decir / dónde citar | por qué (fuente) | tipo |
|---|---|---|---|---|---|
| manuscript/PREPRINT.md | 41 | §1.1 "…δΛ∼1/√N (Sorkin-type arguments in unimodular / causal-set settings). That seed is many orders…" | Citar Ahmed et al. 2004 y Zwane et al. 2018 tras esta frase; D1 + D3 en una frase | Primera mención de Sorkin | T |
| manuscript/PREPRINT.md | 26 | Abstract "a Sorkin–Bekenstein Poisson seed σ₀∼10⁻⁶¹ lies ∼10⁵⁶ below…" | Precisar "horizon-area counting; Everpresent Λ (4-volume) is O(1) and a different model" | D3 | T |
| manuscript/PREPRINT.md | 198 | §5.1, tras la tabla de semillas | D3 con la cita | D3 | T |
| manuscript/PREPRINT.md | 228 | §6.1 "Volume-preserving diffeomorphisms / unimodular structure project out…" | D2 con la cita de Ahmed et al. | D2 | T |
| manuscript/PREPRINT.md | 266 | §7 item 5 "Sorkin seed: used as motivational UV scale…" | Añadir D1 y D4 (la cota σ_X no excluye Everpresent Λ) | D1/D4 | T |
| manuscript/PREPRINT.md | 350 | Ref. 6 "Sorkin, causal-set / unimodular vacuum fluctuation arguments (see repository notes…)" | Referencias completas: Ahmed 2004, Zwane 2018, Sorkin 2007 (0710.1675) | Referencia vaga | T |
| manuscript/CLAIMS.md | 39 | Tras N10 | N11: "The σ_X bound does not test or exclude Everpresent Λ (Ahmed et al. 2004; Zwane et al. 2018)" | D1/D4 | T |
| papers/stochastic-dark-energy-desi-dr2.md | 34 | §1 "If spacetime is fundamentally discrete (as in causal-set theory)… δΛ∼1/√N" | Citar Ahmed et al. tras l.34; D1 | Origen de la idea | T |
| papers/stochastic-dark-energy-desi-dr2.md | 58, 66 | A2: "N = V/L_P⁴"; "With N∼10¹²², δΛ∼10⁻⁶¹" | Tras l.66: citas + D3; resolver la incoherencia 4-volumen/área | D3; l.58 frente a l.54/66 | N |
| papers/stochastic-dark-energy-desi-dr2.md | 384 | §7.3 Implications for Quantum Gravity | Tras l.384: D1 + D4 | D4 | T |
| papers/stochastic-dark-energy-desi-dr2.md | 633, 643, 944 | Apéndice E.2 "Sorkin mechanism"; nota [4]; Ref. [4] gr-qc/0503057 | Corregir [4] → arXiv:0710.1675 (2007) y añadir las dos referencias | Referencia mal emparejada | T |
| papers/quantum-fluid-instabilities-desi-dr2.md | 100, 112, 194 | "Sorkin-type Poisson fluctuations [4,5]"; σ₀∼10⁻⁶¹; Ref. [4] gr-qc/0503057 | Tras l.112: D3; corregir [4] | Ídem | T |
| papers/amplification-gap.md | 36–39 | Tabla "Sorkin / Bekenstein seed σ₀∼1.18×10⁻⁶¹ (holographic d=2)" | Tras l.39: D3 con la cita (el hueco 10⁵⁶ es propio de la elección d=2) | D3 | T |
| papers/amplification-gap.md | 202 | Verdict R1 "…justify a meso-scale causal-set / correlation volume" | Everpresent Λ es la redefinición publicada (conteo de 4-volumen); citar las dos y D4 | D3/D4 | T |
| papers/amplification-gap.md | 276–281 | §7 "Bare seed motivation σ₀∼1/√N, N∼10¹²²" | Nota de una línea: área frente a 4-volumen | D3 | T |
| papers/unimodular-gravity-vacuum-smoothness.md | 22, 246, 275 | Abstract "…its fluctuations are rigorously zero"; §6 "…the denominator V_4 is enormous…"; §7 "…σ_X = 0" | Tras l.246: D2 (es el mismo argumento de Ahmed et al., con la conclusión contraria). Tras l.275: la cancelación afecta a V(x)g_μν local, no al Λ global conjugado a V₄; citar Zwane et al. Una frase tras l.22. | D2/D4 | T |
| papers/scope-and-mixups.md | 52, 54 | Tabla "σ₀∼10⁻⁶¹ \| Motivational Sorkin / holographic seed"; "It is not the Sorkin seed unless…" | Fila nueva tras l.52: "δΛ/Λ ~ O(1) (Everpresent Λ, Ahmed et al. 2004) — not σ_X, not σ₀"; D1 tras l.54 | Confusión típica | T |
| papers/principle-of-vacuum-smoothness.md | 42 | "Mechanisms that predict residual fluctuations in Λ (… causal set theory or the Sorkin mechanism) need a strong suppression mechanism at late times" | Matizar con D1 y D4 (Zwane et al.: ajusta tan bien como ΛCDM) | Afirmación más fuerte de lo que dicen los datos citables | T |
| papers/sdiff-fundamental-vs-emergent.md | 108 | "Do the DESI limits on σ_X rule out… discrete causal sets with large variance?" | Citar las dos y responder con D1/D4: no directamente | D1/D4 | T |
| papers/fundamental-vs-emergent-vacuum-relaxation.md | 88–101, 194 | "Sorkin-only corollary… invisible to DESI and to Euclid"; ref. sin título | Tras l.101: "invisible" vale solo para el conteo de área (D3); corregir la ref. | D3 | T |
| notes/desqueezing-relaxation-vacuum-fluctuations-note.md | 161, 331, 449, 480 | A2 "δΛ∼1/√N"; "Euclid cannot see pure Sorkin noise"; ref.; "Sorkin's 10⁻⁶¹ never reaches Euclid" | Tras l.161: D3 con la cita; matizar "Sorkin's" en l.331 y l.480; corregir la ref. | El resultado de Sorkin/Ahmed et al. es O(1) | T |
| papers/euclid-protocol-vacuum-relaxation.md | 22–25 | "Microscopic seed \| Sorkin / Bekenstein–Hawking σ₀∼10⁻⁶¹ \| Always ≪ BAO noise" | Tras la tabla: D1 + D3 | D1/D3 | T |
| papers/data-pack-option0-internet.md | 16 | "σ₀∼10⁻⁶¹… Bekenstein–Hawking / Sorkin (motivational)" | Añadir las entradas arXiv de las dos | Lista de enlaces del repo | T |
| README.md; papers/HONEST_HEADLINES.md | README 35; HH tras 19 | "Sorkin seed ∼10⁻⁶¹ … ∼10⁵⁶" | README: "horizon-area seed; not Everpresent Λ". HH: fila "I will not say 'causal sets / Sorkin are excluded'" | D1/D3 | T |
| papers/anisotropic-slip-option0.md; results/amplification_routes/VERDICT.md; results/cosmological_mapping/mapping_summary.txt; notes/mapping_tables.md; scripts (varios) | varias | Etiquetas "Sorkin seed", "x Sorkin" | Opcional: "Sorkin (area count)" | D3 | T |

## 6. Dónde añadir los resultados de `analisis_cmb/` y del oscilador

- **R1** = `analisis_cmb/RESULTADO.md`.
- **R2** = `analisis_cmb/oscilador/resultados/tabla.md` (58/60 permitidas por DESI BAO, 0/60 distinguibles en CMB).

Sin redactar el texto.

| fichero | línea | qué dice ahora | dónde añadir / qué matiza o contradice | por qué (fuente) | tipo |
|---|---|---|---|---|---|
| manuscript/PREPRINT.md | 24 | Abstract "Separately, a coherent tachyonic growth model…" | Párrafo R1+R2 tras l.24, antes de "I then quantify…" (l.26) | Resultado de primer nivel | T |
| manuscript/PREPRINT.md | 59–68 | §1.3 Reading guide | Fila para la sección CMB nueva | Navegación | T |
| manuscript/PREPRINT.md | 78 | §2.1 "**Important separation:** DESI analyses that combine BAO with CMB+SN…" | Tras l.78: dataset secundario (Planck 2018 TT ℓ<30 + bandpowers de lensing; parámetros fijos salvo H0, resuelto por θ_*) | Método | T |
| manuscript/PREPRINT.md | 162 | §3.4 "Under the phenomenological OU/QNM residual kernels…" | §3.5 nueva "CMB consistency" tras l.162 (antes de `---` l.164), con R1 y R2 | Consecuencia directa de la cota | T |
| manuscript/PREPRINT.md | 256 | §6.2 "I do not implement a homemade Boltzmann hierarchy… should use community codes" | Matiza: la comprobación CMB usa CAMB (PPF, w(a) tabulada) | Coherencia | T |
| manuscript/PREPRINT.md | 264, 267 | §7 item 3 "Multi-probe fits are a separate analysis." | Matizar con R1; item 7 nuevo con las limitaciones de RESULTADO §6 | Limitaciones | T |
| manuscript/PREPRINT.md | 304–313 | §8.3 Reproducibility | Comandos de `analisis_cmb/` | Reproducibilidad | T |
| manuscript/PREPRINT.md | 320, 325 | §9 item 2; "Future data (Euclid BAO multi-bin products…) can tighten σ_X" | Item nuevo tras l.320 (R1+R2); en l.325, que el CMB no da palanca y Euclid BAO sigue siendo la vía | Conclusiones | T |
| manuscript/CLAIMS.md | 22, 39, 47 | Tras C7; tras N10; "**Out of primary claim:** full multi-probe…" | C8 (R1), C9 (R2); N11 "does not explain the H0 tension (σ(H0)≈0.45 ≪ 5.68 km/s/Mpc)", N12 "no full CMB MCMC"; en l.47, línea "Secondary check (in): Planck 2018 TT ℓ<30 + lensing via CAMB" | R1(c), RESULTADO §6 | T |
| papers/stochastic-dark-energy-desi-dr2.md | 26 | Abstract "I caution that this result is subject to degeneracies…" | R1+R2 tras l.26 | Abstract | T |
| papers/stochastic-dark-energy-desi-dr2.md | 90–98 | A4 "To preserve early-universe constraints (CMB, BBN), a smooth activation factor g(z)…" | **Matiza:** con σ_X ≤ 2.5e-2 el CMB no exige g(z) a esa amplitud (R1); nota tras l.98 | R1(a) | T |
| papers/stochastic-dark-energy-desi-dr2.md | 378 | §7.2 "this constraint applies specifically to the additive OU/QNM kernel…" | R1 tras l.378 | Implicaciones | T |
| papers/stochastic-dark-energy-desi-dr2.md | 384, 386 | §7.3 "A cosmic vacuum modulated by horizon quasi-normal modes would induce harmonic variations… macroscopically ruled out"; "…smooth, dissipative OU evolution as the only mathematically viable mechanism" | **Contradicho por R2:** 58/60 modulaciones armónicas con A ≤ 0.025, P ∈ [0.1, 2] están permitidas por DESI BAO. Reescribir l.384 y l.386 y poner R2 aquí | oscilador/resultados/tabla.md | T |
| papers/stochastic-dark-energy-desi-dr2.md | 440, 442 | §8 tabla F1–F5; "**Conclusion.**" | Fila F6 "CMB (Planck 2018)" tras l.440; mencionar R1 en l.442 | Estado del modelo | T |
| papers/stochastic-dark-energy-desi-dr2.md | 458 | Final de §9 Euclid | Una frase: el CMB no da palanca (R1); la prueba sigue siendo BAO | Programa observacional | T |
| papers/stochastic-dark-energy-desi-dr2.md | 478 | §10 "More broadly, σ_X < 2.5×10⁻² implies…" | R1 (incluida la tensión de H0) y R2, tras l.478 | Conclusión | T |
| papers/stochastic-dark-energy-desi-dr2.md | 934 | Fin del Apéndice I | Nuevo Apéndice J "CMB consistency (Planck 2018)" antes de References (l.936) | Detalle técnico | T |
| papers/HONEST_HEADLINES.md | 19, 29 | Tabla "Things I will not put in a subject line"; "5. Process…" | Filas tras l.19: "The residual explains the Hubble tension" (no) y "Planck tightens σ_X" (no); item 6 tras l.29 con R1+R2 | R1(a)(c) | T |
| README.md | 23, 35, 45–58 | "Read these first"; "Best fit still goes to σ_X→0…"; Layout; comandos | Item 7 `analisis_cmb/RESULTADO.md` tras l.23; párrafo R1+R2 tras l.35; `analisis_cmb/` en Layout y comandos | Resumen | T |
| papers/resume.txt | 17, 57, 85, 87–89 | Cabecera narrativa; fin de §1; §3 "Late-Time Vacuum Homogeneity"; "Future tests:" | Línea en la cabecera tras l.17; sección "1b. CMB CONSISTENCY" tras l.57; línea R1 tras l.85; en "Future tests", el CMB no añade palanca | Resumen numérico | T |
| papers/unimodular-gravity-vacuum-smoothness.md | 82, 288–290 | §2 "…the dark-energy density is homogeneous…"; §8 fila "ISW effect \| No late-time ISW anomaly \| Detection of a residual ISW signal…" | Tras l.82: el CMB no añade cota (R1). **Matiza** la fila ISW: TT ℓ<30 y lensing no distinguen el OU con σ_X ≤ 2.5e-2, así que no puede falsar nada con la precisión actual (nota tras l.290) | R1(a) | T |
| papers/euclid-protocol-vacuum-relaxation.md | 186–189, 247 | §8 fila "ISW / CMB lensing \| DE perturbations if present"; §12 "Working residual ceiling…" | **Matiza** tras l.189 (sin palanca con Planck 2018); R1 como comprobación de consistencia tras l.247 | R1 | T |
| papers/amplification-gap.md | 56–59 | §2 "**What DESI actually measured**…" | Tras l.59: el CMB no aprieta la amplitud por debajo del techo DESI (R1) | R1(a) | T |
| papers/principle-of-vacuum-smoothness.md | 13, 36 | "perfectly homogeneous"; "That robustness across different statistical frameworks…" | Tras l.36: R1 y R2; "perfectly homogeneous" vale solo hasta ≲2.5e-2, también para una X armónica (R2) | R1/R2 | T |
| (todo el repo) | — | No hay ninguna afirmación sobre la tensión de H0 ni sobre θ_* | Introducir R1(c) (CLAIMS N11, HONEST_HEADLINES, PREPRINT §3.5) contrastándolo con Zwane et al. (alivio de H0 en Everpresent Λ); explicar que H0 se re-resuelve fijando 100θ_* = 1.041198 | R1(c); D4 | T |

---

## Recuento

Calculado sobre las filas de las tablas 1–6 con tipo distinto de "—". Una fila puede cubrir varias líneas del mismo fichero; la lista 2b no cuenta.

| punto | filas con cambio | N (cambia número de resultado) | R (regenerar; incl. mixtas T/R) | T (solo texto) |
|---|---|---|---|---|
| 1. Definición de X | 22 | 4 | 0 | 18 |
| 2. Valores de σ_X | 54 | 27 | 9 | 18 |
| 3. CPL | 12 | 4 | 2 | 6 |
| 4. Euclid | 14 | 1 | 2 | 11 |
| 5. Everpresent Λ | 25 | 1 | 0 | 24 |
| 6. Resultados nuevos | 25 | 0 | 0 | 25 |
| **Total** | **152** | **37** | **13** | **102** |

**Ficheros afectados: 49**, todos versionados:
- papers/ 16
- results/ 12
- scripts/ 12
- manuscript/ 2 (PREPRINT, CLAIMS)
- notes/ 2
- figures/ 2 (profile_sigma_x.png, joint_w0wa_sigma.png)
- README.md 1
- analisis_cmb/ 2 (NOTAS, RESULTADO)

Incluye scripts/slip_bridge.py y papers/sensitivity_kernel_table.md, que aparecen en la columna "por qué".

**Cambios N, agrupados por origen:**
1. Aritmética hecha con 1.5e-4 bajo la etiqueta 2.5e-2: slip/γ, N_eff, G_DESI/r, recuentos BAO-safe, H.5/λ/ζ, m_φ, θ_min, "orders of magnitude". Es la mayoría del punto 2.
2. Factor 1.46 entre D1 y D3 (punto 1).
3. Valores CPL mezclados o no archivados (punto 3).
4. Axioma A2, área frente a 4-volumen (punto 5).
5. Rango y bins de Euclid sin fuente, que dos ficheros dan distinto (punto 4).

Ningún cambio N afecta a la cota σ_X < 2.5e-2 ni a los resultados de `analisis_cmb/` (RESULTADO.md, malla, oscilador).

**Pendiente de verificar antes de editar:**
- los números marcados **(verificar)**: r, m_φ, λ, longitudes de N_eff, θ_min, ratios del slip frente a Maus/Sakr, órdenes de magnitud de l.372, normalización del ruido en l.164;
- el volumen y la página de Zwane et al.;
- los valores DESI externos de stochastic l.12/924.
