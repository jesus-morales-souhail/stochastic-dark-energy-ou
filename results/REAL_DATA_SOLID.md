# What is solid from real data (no numerology)

Updated: 2026-07-16 after fixing `scripts/eos_efectiva.py`.

## Data used

| Probe | Source | Vector |
|-------|--------|--------|
| BAO alphas | DESI DR2 (arXiv:2503.14738), repo arrays | 7 bins, diagonal σ |
| Fiducial | Flat ΛCDM, Ωm=0.315, H0=67.4, w0=-1, wa=0 | Defines α=1 |

Limitations (honest): diagonal covariance only; N=7; BAO-only (not CMB+SNe joint).

---

## Solid results (reproduced today)

### 1. Pure CPL background on the same 7 BAO alphas

Script: `scripts/eos_efectiva.py` (rewritten)  
Output: `results/eos_cpl_desi_dr2/`

| Model | w0 | wa | χ² | logL | AIC |
|-------|----|----|-----|------|-----|
| ΛCDM (fixed) | −1 | 0 | 1.665 | −0.832 | — |
| **CPL free** | **−0.990** | **−0.016** | **1.412** | **−0.706** | 5.41 |
| Old bug (do not use) | −2 | −2 | — | −0.832 (fake position) | — |

- Parameters **re-evaluated**: logL at reported (w0,wa) matches the reported logL.
- Pure CPL at the box corner (w0=wa=−2) has logL≈−1663 (not ΛCDM). The old script’s clamp made that corner look like ΛCDM.

### 2. Nested effective EoS (CPL + σ,θ toy term)

| | Value |
|--|-------|
| Best (w0,wa) | same as pure CPL |
| σ | ∼6×10⁻⁷ (→0) |
| ΔlogL (eff−CPL) | ≈0 |
| ΔAIC | **+4** → prefer pure CPL / null |

**No evidence** for the stochastic correction on these BAO points.

### 3. Additive OU/QNM kernel on residuals (fixed ΛCDM background)

Script: `scripts/ou_bao_likelihood.py`  
Log: `results/ou_bao_desi_dr2_run.txt`

| Model | θ | σ_X | ω_R | ΔlogL vs ΛCDM | AIC |
|-------|---|-----|-----|---------------|-----|
| ΛCDM | — | 0 | 0 | 0 | ref |
| H0 OU free | floor | floor (5×10⁻⁵) | 0 | 0 | worse by +4 |
| H1 QNM free | floor | floor | **0** | 0 | worse by +6 |

**Working upper limit (phenomenological):** σ_X < 1.5×10⁻⁴ (95% CL, conservative; see paper).

### 4. What is *not* solid (do not quote as discovery)

| Item | Why |
|------|-----|
| ω_R ≈ 1.388 ≈ ln 4 | Free H1 peak on **DR1-like** vector with θ≈0; **DR2 → ω_R=0**; post-hoc ln4 |
| G_coupling = 5.04×10⁻⁴ | Hardcoded **prediction** in cross-corr script, not a measured coupling |
| Old `eos_efectiva` w0=wa=−2 | Optimizer/report bug + `w_eff≥−1` clamp |
| Cross-corr r≈0.17 | Preliminary; systematics (`WEIGHT_SYS`) not controlled |
| Euclid MCMC figures | Forecasts / mocks, not real Euclid data |
| Magnetar 8 s / “Richter 32” | Real astrophysics; **no** link to this BAO analysis |

---

## Claim-safe one-paragraph summary

Using the public DESI DR2 BAO alpha vector (7 bins, diagonal errors), a clean CPL fit prefers **w0≈−0.99, wa≈−0.02** (near ΛCDM). Adding a nested stochastic effective-EoS term does not improve the likelihood (σ→0, ΔAIC=+4). Separately, an additive OU/QNM kernel on BAO residuals is driven to the numerical floor (σ_X→0, ω_R→0); AIC/BIC prefer the null. The only robust quantitative product for external review remains the **phenomenological upper limit** σ_X < 1.5×10⁻⁴ (95% CL) under the stated kernel and assumptions—not a detection of ln(4), G_coup, or w0=−2.

---

## How to reproduce

```bash
cd stochastic-dark-energy-ou
python3 scripts/eos_efectiva.py          # CPL + nested eff EoS
python3 scripts/ou_bao_likelihood.py     # OU/QNM residual kernel
```
