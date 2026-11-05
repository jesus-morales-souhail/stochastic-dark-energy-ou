# Bloque G: antes / después

Commit de texto: cuesta (falsos positivos de w₀–wₐ), bondad de ajuste de ΛCDM, diseño por Fisher, ecualizador, temblor rápido en el CMB, trabajo futuro; v2.2.
Cada número se copió del `.md` que lo contiene (`analisis_cmb/cuesta/CUESTA.md`, `analisis_cmb/cuesta/GOF.md`, `analisis_cmb/fondo/FISHER_DISENO.md`, `analisis_cmb/cmb/RAPIDO.md`, `analisis_cmb/sn/SN.md`).
No se tocan `results/` ni `figures/`.

**Diferencias con lo pedido:**
- §3.3: "once the ΛCDM fit is required to be as good as the real one" se cambia por "once the ΛCDM fit is also required to remain acceptable (χ² ≤ 19.7 for 11 d.o.f.)", que es el criterio de `CUESTA.md` (el 95 % de χ²₁₁, no el χ² real de 10.27), y se añade el control con ese mismo criterio (3.75 %).
- §9: "against ≤ 15 % at z < 2.5" se cambia por "against ≤ 15 % anywhere in 0.3 < z < 2.3": `FISHER_DISENO.md` da +32–36 % en z = 0.05 (repunte en z → 0), así que el ≤ 15 % solo vale dentro del rango de DESI. Se añade "for θ = 1" a +114 % y +273 %, que son los números de θ = 1 (con θ = 0.3 son +132 % y +416 %).
- §9, TDCOSMO: los resúmenes no dan "5–10 % per lens". H0LiCOW XIII (arXiv:1907.04869) da 73.3 +1.7/−1.8 (2.4 %, seis lentes, perfiles de masa paramétricos); TDCOSMO IV (arXiv:2007.02941) da 74.5 +5.6/−6.1 (≈ 8 %, siete lentes) y 67.4 +4.1/−3.2 (≈ 5 %, con 33 lentes de SLACS) dejando la degeneración de la lámina de masa a la cinemática. Se escribe eso. Loeb 1998 (astro-ph/9802122) es la deriva de redshift, confirmado.
- §9: "the effect of the residual at the limit is 1–2 %" se cambia por los números de `SN.md`: 0.2 % (z = 0.1) a 2.2 % (z = 2) en distancias, los mismos de la frase de las sirenas.
- Se añaden las referencias 13–15 (Wong et al. 2020, Birrer et al. 2020, Loeb 1998).
- La cabecera del preprint y `manuscript/README.md` decían "v2, October 2026" (no se actualizaron en v2.1); se ponen en v2.2.
- `analisis_cmb/README.md`: tres filas nuevas (`cuesta/`, `fondo/FISHER_DISENO.md`, `cmb/RAPIDO.md`); el historial de versiones vive en `PREPRINT.md`, donde se añade la entrada v2.2.

## Propuesta 1: §3.2: bondad de ajuste y ecualizador

Fichero: `manuscript/PREPRINT.md`

**Antes:**

```text
the baseline pipeline keeps the 7 α ([`analisis_cmb/fondo/TRECE.md`](../analisis_cmb/fondo/TRECE.md)).
```

**Después:**

```text
the baseline pipeline keeps the 7 α ([`analisis_cmb/fondo/TRECE.md`](../analisis_cmb/fondo/TRECE.md)). Two independent checks follow. (i) The goodness of fit of $\Lambda$CDM to the 13 numbers ($\chi^2 = 10.27$ for 11 d.o.f.) limits the extra scatter a residual may add: requiring $P(\chi^2 \le 10.27) \ge 5\,\%$ gives $\sigma_X < 4.4\times 10^{-2}$ (θ = 1), $2.8\times 10^{-2}$ (θ = 3) and $2.2\times 10^{-2}$ (θ = 10), 11–25 % above the marginalised limits, which it reproduces without using the kernel; at the fixed-background limits the $\chi^2$ does not move (P from 0.49 to 0.45–0.19), confirming that those limits measure consistency with the reference, not scatter ([`analisis_cmb/cuesta/GOF.md`](../analisis_cmb/cuesta/GOF.md)). (ii) A harmonic residual $A\sin(2\pi\ln a/P + \phi)$ with phase and background marginalised is bounded at $A < 0.013$–$0.028$ for $P \le 1$ (median 0.017), the bound growing ∝ P beyond, where only the slope survives the background ([`analisis_cmb/fondo/FISHER_DISENO.md`](../analisis_cmb/fondo/FISHER_DISENO.md)).
```

## Propuesta 2: §3.3: falsos positivos de w0–wa

Fichero: `manuscript/PREPRINT.md`

**Antes:**

```text
shows the joint fit of $\{w_0, w_a, \sigma_X\}$.
```

**Después:**

```text
shows the joint fit of $\{w_0, w_a, \sigma_X\}$. A stochastic residual does not mimic the $w_0$–$w_a$ preference of DESI DR2 ($\Delta\chi^2 = -4.7$, 1.7σ with BAO alone): in 2000 simulated $\Lambda$CDM universes per case with an OU residual at the limits, the best-fit $(w_0, w_a)$ falls in the DESI quadrant ($w_0 > -1$, $w_a < 0$) half of the time, as without residual; the residual only inflates $\Delta\chi^2$ because the fit covariance ignores it, raising the fraction of universes with $\Delta\chi^2 \ge 4.7$ in that quadrant from 4.75 % (noise alone) to 6–18 % at the fixed-background limits (θ = 0.1–3), and to 5–9 % (against 3.75 %) once the $\Lambda$CDM fit is also required to remain acceptable ($\chi^2 \le 19.7$ for 11 d.o.f.) ([`analisis_cmb/cuesta/CUESTA.md`](../analisis_cmb/cuesta/CUESTA.md)).
```

## Propuesta 3: §3.5: temblor rápido

Fichero: `manuscript/PREPRINT.md`

**Antes:**

```text
$1.1\times 10^{-2}$ at the marginalised θ = 1 ceiling ([`analisis_cmb/cmb/PICOS.md`](../analisis_cmb/cmb/PICOS.md)).
```

**Después:**

```text
$1.1\times 10^{-2}$ at the marginalised θ = 1 ceiling ([`analisis_cmb/cmb/PICOS.md`](../analisis_cmb/cmb/PICOS.md)).

Fast residuals (θ = 3 and 10 at their limits, $1.3$ and $2.2\times 10^{-2}$) are equally invisible: p95 $\Delta\chi^2 = 0.066$–$0.12$ in TT ℓ < 30 + lensing and $\le 2.5\times 10^{-3}$ in the peaks with linear lensing, less than at θ = 1, because the CMB integrates over the memory 1/θ. The implied equation of state oscillates around −1 with an amplitude that depends on the smoothing Δ (p95 |w + 1| = 0.07 at θ = 3, Δ = 0.05; 0.43 at θ = 10, Δ = 0.0125), because an OU path is not differentiable, while no CMB observable depends on Δ (≤ 4 % between Δ = 0.05 and 0.0125). The Takahashi halofit, which uses w(a) explicitly, produces a spurious growth of $C_L^{\phi\phi}$ for such w (p95 $\Delta\chi^2$ up to 0.65 at θ = 10, Δ = 0.0125, with P(k, z) unchanged); the physical reading is linear lensing, and at θ = 1 the effect is negligible ([`analisis_cmb/cmb/RAPIDO.md`](../analisis_cmb/cmb/RAPIDO.md)).
```

## Propuesta 4: Limitación 6: halofit

Fichero: `manuscript/PREPRINT.md`

**Antes:**

```text
halofit with tabulated w(a) in PPF not verified.
```

**Después:**

```text
the Takahashi halofit misbehaves for a rapidly oscillating w(a) ([`analisis_cmb/cmb/RAPIDO.md`](../analisis_cmb/cmb/RAPIDO.md)); results at θ ≥ 3 use linear lensing.
```

## Propuesta 5: §9: diseño de cartografiados, retardos temporales y deriva de z

Fichero: `manuscript/PREPRINT.md`

**Antes:**

```text
while Euclid preparation VII (arXiv:1910.09273) projects the per-bin $f\sigma_8$ onto cosmological parameters without quoting them ([`analisis_cmb/fondo/CRECIMIENTO.md`](../analisis_cmb/fondo/CRECIMIENTO.md)).
```

**Después:**

```text
while Euclid preparation VII (arXiv:1910.09273) projects the per-bin $f\sigma_8$ onto cosmological parameters without quoting them ([`analisis_cmb/fondo/CRECIMIENTO.md`](../analisis_cmb/fondo/CRECIMIENTO.md)). Survey design: with the background marginalised, the Fisher information on $(\sigma_X, θ)$ is set by the survey volume, not by the number of bins (N = 7 already holds 90 % of the information on θ for θ ≤ 1; finer binning pays only for θ ≳ 3), and a new BAO point at z > 2.5 adds more than any within the DESI range (+114 % at z = 3 and +273 % at z = 4 for θ = 1, against ≤ 15 % anywhere in 0.3 < z < 2.3): the next measurement belongs to z = 2.5–5 surveys (DESI-II, Spec-S5) rather than to the Euclid range ([`analisis_cmb/fondo/FISHER_DISENO.md`](../analisis_cmb/fondo/FISHER_DISENO.md)). Strong-lens time delays (TDCOSMO: H0 to 2.4 % from six lenses with parametric mass profiles, Wong et al. 2020, arXiv:1907.04869; 8 % from seven lenses, or 5 % adding 33 SLACS lenses, once the mass-sheet degeneracy is constrained by kinematics alone, Birrer et al. 2020, arXiv:2007.02941) and redshift drift (Loeb 1998, arXiv:astro-ph/9802122; ELT-ANDES) would measure the expansion with neither ruler nor candle; the effect of the residual at the limit is the same 0.2 % (z = 0.1) to 2.2 % (z = 2) on distances as for the sirens.
```

## Propuesta 6: Referencias 13–15

Fichero: `manuscript/PREPRINT.md`

**Antes:**

```text
12. Brout, D., et al., "The Pantheon+ Analysis: Cosmological Constraints," Astrophys. J. 938, 110 (2022), [arXiv:2202.04077](https://arxiv.org/abs/2202.04077).
```

**Después:**

```text
12. Brout, D., et al., "The Pantheon+ Analysis: Cosmological Constraints," Astrophys. J. 938, 110 (2022), [arXiv:2202.04077](https://arxiv.org/abs/2202.04077).
13. Wong, K. C., et al. (H0LiCOW), "H0LiCOW XIII. A 2.4% measurement of $H_0$ from lensed quasars: 5.3σ tension between early- and late-Universe probes," Mon. Not. R. Astron. Soc. 498, 1420 (2020), [arXiv:1907.04869](https://arxiv.org/abs/1907.04869).
14. Birrer, S., et al. (TDCOSMO), "TDCOSMO IV: Hierarchical time-delay cosmography — joint inference of the Hubble constant and galaxy density profiles," Astron. Astrophys. 643, A165 (2020), [arXiv:2007.02941](https://arxiv.org/abs/2007.02941).
15. Loeb, A., "Direct Measurement of Cosmological Parameters from the Cosmic Deceleration of Extragalactic Objects," Astrophys. J. 499, L111 (1998), [arXiv:astro-ph/9802122](https://arxiv.org/abs/astro-ph/9802122).
```

## Propuesta 7: Historial de versiones: v2.2

Fichero: `manuscript/PREPRINT.md`

**Antes:**

```text
- Audit and before/after records: [`analisis_cmb/AUDITORIA.md`]
```

**Después:**

```text
- **v2.2 (November 2026):** goodness-of-fit and harmonic-residual checks of the limits (§3.2); a stochastic residual does not mimic the $w_0$–$w_a$ preference of DESI DR2 (§3.3); fast residuals θ = 3, 10 in the CMB and the Takahashi halofit caveat (§3.5, limitation 6); survey design by Fisher, strong-lens time delays and redshift drift (§9); evidence of C2, C8 and C9 extended. Before/after: block [G](../analisis_cmb/bloque_G/COMPARACION_G.md).
- Audit and before/after records: [`analisis_cmb/AUDITORIA.md`]
```

## Propuesta 8: Cabecera: fecha/versión

Fichero: `manuscript/PREPRINT.md`

**Antes:**

```text
| Date | v2, October 2026 (v1: July 2026) |
```

**Después:**

```text
| Date | v2.2, November 2026 (v2: October 2026; v1: July 2026) |
```

## Propuesta 9: manuscript/README.md: versión

Fichero: `manuscript/README.md`

**Antes:**

```text
**Preprint, v2 (October 2026; v1: July 2026).**
```

**Después:**

```text
**Preprint, v2.2 (November 2026; v2: October 2026; v1: July 2026).**
```

## Propuesta 10: CLAIMS C2: evidencia

Fichero: `manuscript/CLAIMS.md`

**Antes:**

```text
robustness: `analisis_cmb/reloj/RELOJ.md` (stationarity clock) |
```

**Después:**

```text
robustness: `analisis_cmb/reloj/RELOJ.md` (stationarity clock); checks: `analisis_cmb/cuesta/GOF.md` (goodness of fit of ΛCDM), `analisis_cmb/cuesta/CUESTA.md` (w₀–wₐ false positives) |
```

## Propuesta 11: CLAIMS C8: evidencia

Fichero: `manuscript/CLAIMS.md`

**Antes:**

```text
the CMB adds no bound on σ_X | `analisis_cmb/RESULTADO.md` |
```

**Después:**

```text
the CMB adds no bound on σ_X | `analisis_cmb/RESULTADO.md`; fast residuals (θ = 3, 10): `analisis_cmb/cmb/RAPIDO.md` |
```

## Propuesta 12: CLAIMS C9: ecualizador

Fichero: `manuscript/CLAIMS.md`

**Antes:**

```text
is allowed by DESI DR2 BAO in 58/60 cases and is CMB-indistinguishable in 60/60 | `analisis_cmb/oscilador/resultados/tabla.md` |
```

**Después:**

```text
is allowed by DESI DR2 BAO in 58/60 cases and is CMB-indistinguishable in 60/60; with the background marginalised, A < 0.013–0.028 for P ≤ 1 (`analisis_cmb/fondo/FISHER_DISENO.md`) | `analisis_cmb/oscilador/resultados/tabla.md`; `analisis_cmb/fondo/FISHER_DISENO.md` |
```

## Propuesta 13: README: FISHER_DISENO.md

Fichero: `analisis_cmb/README.md`

**Antes:**

```text
| `fondo/CRECIMIENTO.md` | Change of fσ₈(z) at the C2 limit (θ → 0 mode, δΩ_Λ = 2.633e-2 on Planck ΛCDM, growth equation solved; `crecimiento.py`): 3.8 % at z = 0.5, 3.5–5 times below the DESI DR1 full-shape errors per bin |
```

**Después:**

```text
| `fondo/CRECIMIENTO.md` | Change of fσ₈(z) at the C2 limit (θ → 0 mode, δΩ_Λ = 2.633e-2 on Planck ΛCDM, growth equation solved; `crecimiento.py`): 3.8 % at z = 0.5, 3.5–5 times below the DESI DR1 full-shape errors per bin |
| `fondo/FISHER_DISENO.md` | Survey design by Fisher information on (σ_X, θ) with the background marginalised (`fisher_diseno.py`): gain of a new BAO point versus z (optimum at z > 2.5), information versus number of bins at fixed volume (flat for θ ≤ 1), and the 95 % bound on a harmonic residual A sin(2π ln a/P + φ) versus P with phase and background marginalised (equaliser) |
```

## Propuesta 14: README: RAPIDO.md y cuesta/

Fichero: `analisis_cmb/README.md`

**Antes:**

```text
| `cmb/` | High-ℓ CMB check (`PICOS.md`, `picos.py`): TT, TE, EE distance to ΛCDM over ℓ ∈ [30, 2500], lensing in the troughs, ∝ σ_X² |
```

**Después:**

```text
| `cmb/` | High-ℓ CMB check (`PICOS.md`, `picos.py`): TT, TE, EE distance to ΛCDM over ℓ ∈ [30, 2500], lensing in the troughs, ∝ σ_X² |
| `cmb/RAPIDO.md` | Fast residuals θ = 3 and 10 at the θ = 10 limits (`rapido.py`): Planck TT ℓ < 30 + lensing and the high-ℓ distance, smoothing Δ check, amplitude of w_eff; the Takahashi halofit misbehaves for a rapidly oscillating w(a), so every case is also run with linear lensing |
| `cuesta/` | Does an OU residual mimic the DESI DR2 w₀–wₐ preference? 2000 simulated ΛCDM universes per (σ_X, θ) with the 13 numbers, fitted with ΛCDM and w₀wₐCDM (`CUESTA.md`, `cuesta.py`); bound on σ_X from the goodness of fit of ΛCDM, P(χ² ≤ 10.27) = 5 % (`GOF.md`, `gof.py`) |
```
