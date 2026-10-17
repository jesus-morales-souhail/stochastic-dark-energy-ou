# Bloque D, paso 1: propuestas de texto, sin aplicar (2026-10-17)

**Estado:** solo propuestas. **Ningún paper editado.** Base: commit `84fbe28` (después de la reorganización y de los bloques A, B, C).

**Cómo se ha hecho.**
- Cada punto lleva un ancla, un trozo único de la línea. El generador busca el ancla en el fichero actual y copia de ahí el número de línea y el texto, así que no están escritos a mano.
- Las líneas de AUDITORIA (tablas 4, 5 y 6) son anteriores a la reorganización y a los bloques B2 y C; aquí van las de hoy.
- "Insertar tras" significa un párrafo o línea nuevos después de la línea citada.
- En los textos citados, `]\(` es un enlace Markdown `](` escapado, para que se vea tal cual y no se interprete en este documento. Al aplicar, se quita la barra.
- Los números de R1 son los de σ_X = 2.633e-2 (`sigma_2.633e-02/resumen.json`); 2.506e-2 se usa solo para la malla. Los del slip son los del bloque C.

**Leído antes:**
- `analisis_cmb/everpresent/EVERPRESENT_vs_DESI.md` (entero);
- `analisis_cmb/RESULTADO.md`;
- `analisis_cmb/oscilador/resultados/tabla.md`;
- `analisis_cmb/AUDITORIA.md` §4, §5 y §6;
- `analisis_cmb/bloque_B/COMPARACION_B2.md`, "Afirmaciones que cambian de sentido".

**Referencias completas de Everpresent Λ** (para todas las citas de D2):

- Ahmed, M., Dodelson, S., Greene, P. B., & Sorkin, R. D., "Everpresent Λ," Phys. Rev. D 69, 103523 (2004), arXiv:astro-ph/0209274.
- Zwane, N., Afshordi, N., & Sorkin, R. D., "Cosmological tests of Everpresent Λ," Class. Quantum Grav. (2018), doi:10.1088/1361-6382/aadc36, arXiv:1703.06265.
- Sorkin, R. D., "Is the cosmological 'constant' a nonlocal quantum residue of discreteness of the causal set type?," AIP Conf. Proc. 957, 142 (2007), arXiv:0710.1675.
- Das, Nasiri & Yazdi, "Aspects of Everpresent Λ (I)," arXiv:2304.03819 (2023).
- Das, Nasiri & Yazdi, "Aspects of Everpresent Λ (II)," arXiv:2307.13743 (2023).
- Barrow, J. D., "A strong constraint on ever-present Lambda," Phys. Rev. D 75, 067301 (2007), arXiv:gr-qc/0612128.

Sobre las referencias:
- Las de Ahmed et al. y Sorkin 2007 están comprobadas en arXiv (AUDITORIA §5).
- Zwane et al.: AUDITORIA §5 dice que el volumen y la página (35, 194002) no aparecen en la ficha de arXiv. Van solo el DOI y el arXiv, pendientes de verificar.
- Das et al. y Barrow: las de `LECTURA_EVERPRESENT.md`. Allí Das et al. no tienen iniciales ni revista, así que se dejan solo los apellidos y el arXiv.

## D1. Euclid: fecha del DR1 con productos de cosmología y supuestos de diseño

Texto estándar **[E]** (inglés): *"the full Euclid DR1 release with cosmology products (expected mid-2027; the DR1-Foundation release of November 2026 does not include the LE3 high-level products)"*. Los bins y rangos en z de Euclid se marcan como *design assumption*: no hay binning BAO de Euclid publicado. Fuente de [E]: AUDITORIA.md §4 (página ESA cosmos.esa.int/web/euclid/dr1-timeline, "Update published on 15 June 2026", leída el 2026-10-03): DR1-Foundation en noviembre de 2026 sin vector BAO; DR1 completo con clustering y lensing a mediados de 2027.

### D1.1 — `papers/stochastic-dark-energy-desi-dr2.md` l.446 (reemplazar)

**Texto actual:**

> ## 9. Near-Term Observational Program: Euclid DR1

**Propuesto:**

> ## 9. Near-Term Observational Program: Euclid DR1 (full release, expected mid-2027)

**Fuente:** AUDITORIA.md §4 (página ESA cosmos.esa.int/web/euclid/dr1-timeline, "Update published on 15 June 2026", leída el 2026-10-03): DR1-Foundation en noviembre de 2026 sin vector BAO; DR1 completo con clustering y lensing a mediados de 2027.

### D1.2 — `papers/stochastic-dark-energy-desi-dr2.md` l.448 (reemplazar)

**Texto actual:**

> Euclid Data Release 1 (expected H2 2026) will provide $>20$ redshift bins, resolving the current degeneracies with substantially higher statistical power. The analysis pipeline is fully ready.

**Propuesto:**

> The decisive test will come from the full Euclid DR1 release with cosmology products (expected mid-2027; the DR1-Foundation release of November 2026 does not include the LE3 high-level products). Assuming its BAO vector has many more redshift bins than the 7 of DESI DR2 (a design assumption: no Euclid BAO binning has been published), it will resolve the current degeneracies with substantially higher statistical power. The analysis pipeline is ready to take that vector.

**Fuente:** AUDITORIA.md §4 (página ESA cosmos.esa.int/web/euclid/dr1-timeline, "Update published on 15 June 2026", leída el 2026-10-03): DR1-Foundation en noviembre de 2026 sin vector BAO; DR1 completo con clustering y lensing a mediados de 2027. AUDITORIA.md tabla 4: el número de bins y el rango en z no tienen fuente en el repo; `papers/sensitivity_kernel_table.md` l.37 da Δx≈0.5 y ω_R,min≈12.6 para "Euclid DR1", y `papers/stochastic-dark-energy-desi-dr2.md` l.609 da 0.388 y 16.2: no coinciden.

### D1.3 — `papers/stochastic-dark-energy-desi-dr2.md` l.442 (reemplazar)

**Texto actual:**

> **Conclusion.** The stochastic model is not favored by the current DESI DR2 BAO data when the background is allowed to vary, and I place a phenomenological upper limit on its amplitude. A definitive falsification would require a simultaneous fit of the stochastic parameters with the CPL background and, ideally, the $>20$ bins of Euclid DR1. This work provides the necessary pipeline and a reference limit for that future analysis.

**Propuesto:**

> **Conclusion.** The stochastic model is not favored by the current DESI DR2 BAO data when the background is allowed to vary, and I place a phenomenological upper limit on its amplitude. A definitive falsification would require a simultaneous fit of the stochastic parameters with the CPL background and, ideally, the BAO vector of the full Euclid DR1 release (expected mid-2027), with whatever number of bins it publishes. This work provides the necessary pipeline and a reference limit for that future analysis.

**Fuente:** AUDITORIA.md §4 (página ESA cosmos.esa.int/web/euclid/dr1-timeline, "Update published on 15 June 2026", leída el 2026-10-03): DR1-Foundation en noviembre de 2026 sin vector BAO; DR1 completo con clustering y lensing a mediados de 2027. AUDITORIA.md tabla 4: el número de bins y el rango en z no tienen fuente en el repo; `papers/sensitivity_kernel_table.md` l.37 da Δx≈0.5 y ω_R,min≈12.6 para "Euclid DR1", y `papers/stochastic-dark-energy-desi-dr2.md` l.609 da 0.388 y 16.2: no coinciden.

### D1.4 — `papers/stochastic-dark-energy-desi-dr2.md` l.450 (reemplazar)

**Texto actual:**

> Notably, Euclid DR1's narrower redshift baseline ($z \in [0.9, 1.8]$) yields a higher Rayleigh frequency limit ($\omega_{R,\min} \approx 16.2$) than DESI (6.66), meaning Euclid cannot geometrically resolve intermediate-frequency quasi-normal mode (QNM) oscillations that DESI could potentially detect (Appendix D). The decisive contribution of Euclid will be statistical power through $>20$ bins, not oscillatory frequency resolution.

**Propuesto:**

> Notably, if the Euclid BAO baseline is narrower than DESI's — I assume $z \in [0.9, 1.8]$ as a design case, not a published Euclid binning — the Rayleigh frequency limit is higher ($\omega_{R,\min} \approx 16.2$) than DESI's (6.66), so such a survey could not geometrically resolve intermediate-frequency QNM oscillations that DESI could potentially detect (Appendix D). Euclid's decisive contribution would then be statistical power from more bins, not oscillatory frequency resolution.

**Fuente:** AUDITORIA.md tabla 4: el número de bins y el rango en z no tienen fuente en el repo; `papers/sensitivity_kernel_table.md` l.37 da Δx≈0.5 y ω_R,min≈12.6 para "Euclid DR1", y `papers/stochastic-dark-energy-desi-dr2.md` l.609 da 0.388 y 16.2: no coinciden.

### D1.5 — `papers/stochastic-dark-energy-desi-dr2.md` l.609 (reemplazar)

**Texto actual:**

> | Euclid DR1 | [0.9, 1.8] | 0.388 | 16.2 | ✗ No (geometric limit) |

**Propuesto:**

> | Euclid DR1 (design assumption) | [0.9, 1.8] (assumed) | 0.388 | 16.2 | ✗ No (geometric limit) |

**Fuente:** AUDITORIA.md tabla 4: el número de bins y el rango en z no tienen fuente en el repo; `papers/sensitivity_kernel_table.md` l.37 da Δx≈0.5 y ω_R,min≈12.6 para "Euclid DR1", y `papers/stochastic-dark-energy-desi-dr2.md` l.609 da 0.388 y 16.2: no coinciden.

### D1.6 — `papers/stochastic-dark-energy-desi-dr2.md` l.452 (reemplazar)

**Texto actual:**

> Upon release of Euclid DR1, I will:

**Propuesto:**

> Upon release of the full Euclid DR1 with LE3 cosmology products (expected mid-2027; the November 2026 DR1-Foundation release has no BAO vector), I will:

**Fuente:** AUDITORIA.md §4 (página ESA cosmos.esa.int/web/euclid/dr1-timeline, "Update published on 15 June 2026", leída el 2026-10-03): DR1-Foundation en noviembre de 2026 sin vector BAO; DR1 completo con clustering y lensing a mediados de 2027.

### D1.7 — `papers/stochastic-dark-energy-desi-dr2.md` l.26 (reemplazar)

**Texto actual:**

> I caution that this result is subject to degeneracies between the stochastic component and the CPL parameters, which cannot be fully resolved with only 7 BAO bins. The primary contribution of this work is a validated analysis pipeline and a benchmark for future analyses with the $>20$ bins of Euclid DR1 (expected H2 2026).

**Propuesto:**

> I caution that this result is subject to degeneracies between the stochastic component and the CPL parameters, which cannot be fully resolved with only 7 BAO bins. The primary contribution of this work is a validated analysis pipeline and a benchmark for future analyses with the full Euclid DR1 release with cosmology products (expected mid-2027; the DR1-Foundation release of November 2026 does not include the LE3 high-level products).

**Fuente:** AUDITORIA.md §4 (página ESA cosmos.esa.int/web/euclid/dr1-timeline, "Update published on 15 June 2026", leída el 2026-10-03): DR1-Foundation en noviembre de 2026 sin vector BAO; DR1 completo con clustering y lensing a mediados de 2027.

### D1.8 — `papers/stochastic-dark-energy-desi-dr2.md` l.480 (matizar (sin cambio; basta [E] una vez en §9))

**Texto actual:**

> The correlation function derived in Appendix F provides a falsifiable prediction for Euclid DR1: if a stochastic component exists, the whitened BAO residuals should follow:

**Propuesto:**

> No change: the date is given once in §9 (D1.1–D1.2). Same for l.611 ("Euclid DR1 geometrically cannot", now conditional on the design case of D1.4) and l.756 ("**Euclid DR1.**").

**Fuente:** AUDITORIA.md tabla 4, fila l.480/611/761: "Basta con la aclaración [E] una vez en §9".

### D1.9 — `papers/resume.txt` l.88 (reemplazar)

**Texto actual:**

> - Euclid DR1 (expected H2 2026) will provide >20 BAO bins, offering significantly higher sensitivity to fluctuations at the level of σ_X > 10⁻⁵.

**Propuesto:**

> - The full Euclid DR1 release with cosmology products (expected mid-2027; the November 2026 DR1-Foundation release has no LE3 products) should provide many more BAO bins than DESI DR2 (design assumption: no Euclid binning published), with sensitivity aimed at σ_X ∼ 10⁻⁵ (forecast target, not a measurement).

**Fuente:** AUDITORIA.md §4 (página ESA cosmos.esa.int/web/euclid/dr1-timeline, "Update published on 15 June 2026", leída el 2026-10-03): DR1-Foundation en noviembre de 2026 sin vector BAO; DR1 completo con clustering y lensing a mediados de 2027. 10⁻⁵ es el objetivo aspiracional de `papers/amplification-gap.md` (tabla de objetivos).

### D1.10 — `archive/results/ou_bao_desi_dr2_run.txt` l.31 (matizar (salida archivada: no se edita))

**Texto actual:**

>   Decisive test: >20 bins (Euclid DR1, expected H2 2026).

**Propuesto:**

> No change: archived console output (`archive/README.md`). The current scripts no longer print a date.

**Fuente:** AUDITORIA.md tabla 4 (tipo R); el fichero está en `archive/` desde la reorganización.

### D1.11 — `theory/unimodular-gravity-vacuum-smoothness.md` l.275 (reemplazar)

**Texto actual:**

> ## 8. Observational Prospects: The Test of Euclid DR1

**Propuesto:**

> ## 8. Observational Prospects: The Test of the Full Euclid DR1 Release (expected mid-2027)

**Fuente:** AUDITORIA.md §4 (página ESA cosmos.esa.int/web/euclid/dr1-timeline, "Update published on 15 June 2026", leída el 2026-10-03): DR1-Foundation en noviembre de 2026 sin vector BAO; DR1 completo con clustering y lensing a mediados de 2027.

### D1.12 — `theory/unimodular-gravity-vacuum-smoothness.md` l.281 (reemplazar)

**Texto actual:**

> | Euclid DR1 BAO (2026) | $\sigma_X = 0$ exactly (within experimental uncertainty) | Any statistically significant deviation ($\sigma_X > 0$) |

**Propuesto:**

> | Euclid DR1 BAO (full release, expected mid-2027) | $\sigma_X = 0$ exactly (within experimental uncertainty) | Any statistically significant deviation ($\sigma_X > 0$) |

**Fuente:** AUDITORIA.md §4 (página ESA cosmos.esa.int/web/euclid/dr1-timeline, "Update published on 15 June 2026", leída el 2026-10-03): DR1-Foundation en noviembre de 2026 sin vector BAO; DR1 completo con clustering y lensing a mediados de 2027.

### D1.13 — `theory/unimodular-gravity-vacuum-smoothness.md` l.303 (reemplazar (solo la fecha))

**Texto actual:**

> Future surveys, such as Euclid DR1 (expected 2026), may further tighten the limit on $\sigma_X$ or even detect a residual signal. If the smoothness persists, I think the theoretical community must take seriously the possibility that the universe is geometrically unimodular at its foundation.

**Propuesto:**

> "Euclid DR1 (expected 2026)" → "the full Euclid DR1 release with cosmology products (expected mid-2027)"; the rest of the sentence unchanged.

**Fuente:** AUDITORIA.md §4 (página ESA cosmos.esa.int/web/euclid/dr1-timeline, "Update published on 15 June 2026", leída el 2026-10-03): DR1-Foundation en noviembre de 2026 sin vector BAO; DR1 completo con clustering y lensing a mediados de 2027.

### D1.14 — `theory/sdiff-fundamental-vs-emergent.md` l.83 (reemplazar)

**Texto actual:**

> | Euclid DR1 (2026) | No detection of stochastic signal | Possible detection of small $\sigma_X$ | High |

**Propuesto:**

> | Euclid DR1, full release (expected mid-2027) | No detection of stochastic signal | Possible detection of small $\sigma_X$ | High |

**Fuente:** AUDITORIA.md §4 (página ESA cosmos.esa.int/web/euclid/dr1-timeline, "Update published on 15 June 2026", leída el 2026-10-03): DR1-Foundation en noviembre de 2026 sin vector BAO; DR1 completo con clustering y lensing a mediados de 2027.

### D1.15 — `theory/euclid-protocol-vacuum-relaxation.md` l.100 (reemplazar)

**Texto actual:**

> For Euclid, I replace the DESI 7-bin vector by the Euclid BAO vector ($\gtrsim 20$ bins when available) and recompute $S(z)$.

**Propuesto:**

> For Euclid, I replace the DESI 7-bin vector by the BAO vector of the full Euclid DR1 release with cosmology products (expected mid-2027; the DR1-Foundation release of November 2026 does not include the LE3 high-level products) — I assume $\gtrsim 20$ bins as a design case; no Euclid binning has been published — and recompute $S(z)$.

**Fuente:** AUDITORIA.md §4 (página ESA cosmos.esa.int/web/euclid/dr1-timeline, "Update published on 15 June 2026", leída el 2026-10-03): DR1-Foundation en noviembre de 2026 sin vector BAO; DR1 completo con clustering y lensing a mediados de 2027. AUDITORIA.md tabla 4: el número de bins y el rango en z no tienen fuente en el repo; `papers/sensitivity_kernel_table.md` l.37 da Δx≈0.5 y ω_R,min≈12.6 para "Euclid DR1", y `papers/stochastic-dark-energy-desi-dr2.md` l.609 da 0.388 y 16.2: no coinciden.

### D1.16 — `theory/euclid-protocol-vacuum-relaxation.md` l.150 (matizar (añadir al inicio de la frase))

**Texto actual:**

> On a narrower Euclid DR1 path ($\Delta x\sim 0.39$), the same amplitude criterion requires larger $\theta$ by a factor $\sim 0.94/0.39\sim 2.4$. Thus $\theta\sim\mathcal{O}(1)$ is relevant when $A_0$ is a few times the noise $s$.

**Propuesto:**

> Prefix: "In the design case of a narrower Euclid path ($\Delta x\sim 0.39$, assumed, not published), …"; rest unchanged.

**Fuente:** AUDITORIA.md tabla 4: el número de bins y el rango en z no tienen fuente en el repo; `papers/sensitivity_kernel_table.md` l.37 da Δx≈0.5 y ω_R,min≈12.6 para "Euclid DR1", y `papers/stochastic-dark-energy-desi-dr2.md` l.609 da 0.388 y 16.2: no coinciden.

### D1.17 — `theory/euclid-protocol-vacuum-relaxation.md` l.159 (insertar tras (nota bajo la tabla de l.156–160))

**Línea tras la que se inserta:**

> | Euclid-like | 20 | $\sim 0.24$ | $\gtrsim 5$ |

**Propuesto:**

> *The Euclid-like rows (20 and 40 bins) are design assumptions, not a published Euclid binning.*

**Fuente:** AUDITORIA.md tabla 4: el número de bins y el rango en z no tienen fuente en el repo; `papers/sensitivity_kernel_table.md` l.37 da Δx≈0.5 y ω_R,min≈12.6 para "Euclid DR1", y `papers/stochastic-dark-energy-desi-dr2.md` l.609 da 0.388 y 16.2: no coinciden.

### D1.18 — `theory/euclid-protocol-vacuum-relaxation.md` l.260 (matizar (añadir al inicio))

**Texto actual:**

> When real Euclid BAO summary statistics are public, replace the retired synthetic MCMC with a loader analogous to `desi_dr2_data.py` and re-run. Until then, empirical residual claims in this repo rest on **DESI DR2 only**.

**Propuesto:**

> Prefix: "BAO summary statistics are expected with the full Euclid DR1 release (mid-2027; the November 2026 DR1-Foundation release does not include them)." Then the existing sentence.

**Fuente:** AUDITORIA.md §4 (página ESA cosmos.esa.int/web/euclid/dr1-timeline, "Update published on 15 June 2026", leída el 2026-10-03): DR1-Foundation en noviembre de 2026 sin vector BAO; DR1 completo con clustering y lensing a mediados de 2027.

### D1.19 — `theory/euclid-protocol-vacuum-relaxation.md` l.269 (reemplazar)

**Texto actual:**

> [4] Euclid Collaboration, DR1 BAO forecasts (as available at analysis time).

**Propuesto:**

> [4] ESA, Euclid DR1 timeline, https://www.cosmos.esa.int/web/euclid/dr1-timeline (update of 15 June 2026): DR1-Foundation November 2026; full DR1 with cosmology products mid-2027.

**Fuente:** AUDITORIA.md §4 (página ESA cosmos.esa.int/web/euclid/dr1-timeline, "Update published on 15 June 2026", leída el 2026-10-03): DR1-Foundation en noviembre de 2026 sin vector BAO; DR1 completo con clustering y lensing a mediados de 2027.

### D1.20 — `theory/fundamental-vs-emergent-vacuum-relaxation.md` l.151 (reemplazar (solo el inciso))

**Texto actual:**

> 2. **Euclid null** $\sigma_X\ll 10^{-5}$ with $\gtrsim 20$ bins: forces $A_0 e^{-\theta\Delta x}$ small; compatible with **F** and **E0**.

**Propuesto:**

> "with $\gtrsim 20$ bins" → "with the BAO vector of the full Euclid DR1 release (expected mid-2027; $\gtrsim 20$ bins is a design assumption)".

**Fuente:** AUDITORIA.md §4 (página ESA cosmos.esa.int/web/euclid/dr1-timeline, "Update published on 15 June 2026", leída el 2026-10-03): DR1-Foundation en noviembre de 2026 sin vector BAO; DR1 completo con clustering y lensing a mediados de 2027. AUDITORIA.md tabla 4: el número de bins y el rango en z no tienen fuente en el repo; `papers/sensitivity_kernel_table.md` l.37 da Δx≈0.5 y ω_R,min≈12.6 para "Euclid DR1", y `papers/stochastic-dark-energy-desi-dr2.md` l.609 da 0.388 y 16.2: no coinciden.

### D1.21 — `exploratory/results/sdiff_discrimination/summary.txt` l.34 (reemplazar)

**Texto actual:**

>   2) Euclid null << 1e-5 with >20 bins:

**Propuesto:**

>   2) Euclid null << 1e-5 (full DR1, expected mid-2027; >20 bins assumed):

**Fuente:** AUDITORIA.md §4 (página ESA cosmos.esa.int/web/euclid/dr1-timeline, "Update published on 15 June 2026", leída el 2026-10-03): DR1-Foundation en noviembre de 2026 sin vector BAO; DR1 completo con clustering y lensing a mediados de 2027. AUDITORIA.md tabla 4: el número de bins y el rango en z no tienen fuente en el repo; `papers/sensitivity_kernel_table.md` l.37 da Δx≈0.5 y ω_R,min≈12.6 para "Euclid DR1", y `papers/stochastic-dark-energy-desi-dr2.md` l.609 da 0.388 y 16.2: no coinciden.

### D1.22 — `exploratory/notes/desqueezing-relaxation-vacuum-fluctuations-note.md` l.374 (reemplazar)

**Texto actual:**

> ## 7. What Euclid DR1 can and cannot see (emergent scenario)

**Propuesto:**

> ## 7. What the full Euclid DR1 release (expected mid-2027) can and cannot see (emergent scenario)

**Fuente:** AUDITORIA.md §4 (página ESA cosmos.esa.int/web/euclid/dr1-timeline, "Update published on 15 June 2026", leída el 2026-10-03): DR1-Foundation en noviembre de 2026 sin vector BAO; DR1 completo con clustering y lensing a mediados de 2027.

### D1.23 — `exploratory/notes/desqueezing-relaxation-vacuum-fluctuations-note.md` l.385 (reemplazar (solo el inciso))

**Texto actual:**

> - Path-dependent damping if $\theta\sim\mathcal{O}(0.1$–$1)$ and the effective kick is not negligible: then the BAO covariance kernel shape (not only the amplitude) becomes informative once $N_{\rm bins}\gtrsim 20$.

**Propuesto:**

> "once $N_{\rm bins}\gtrsim 20$" → "once $N_{\rm bins}\gtrsim 20$ (design assumption for Euclid)".

**Fuente:** AUDITORIA.md tabla 4: el número de bins y el rango en z no tienen fuente en el repo; `papers/sensitivity_kernel_table.md` l.37 da Δx≈0.5 y ω_R,min≈12.6 para "Euclid DR1", y `papers/stochastic-dark-energy-desi-dr2.md` l.609 da 0.388 y 16.2: no coinciden.

### D1.24 — `theory/principle-of-vacuum-smoothness.md` l.50 (reemplazar (opcional, solo el inciso))

**Texto actual:**

> I propose this principle on the basis of current observational constraints. It is falsifiable: a statistically significant detection of stochastic fluctuations in future high-precision surveys (for example Euclid DR1 or DESI Year 3+) would force me to revise or abandon it.

**Propuesto:**

> "(for example Euclid DR1 or DESI Year 3+)" → "(for example the full Euclid DR1 release, expected mid-2027, or later DESI data)".

**Fuente:** AUDITORIA.md §4 (página ESA cosmos.esa.int/web/euclid/dr1-timeline, "Update published on 15 June 2026", leída el 2026-10-03): DR1-Foundation en noviembre de 2026 sin vector BAO; DR1 completo con clustering y lensing a mediados de 2027.

## D2. Everpresent Λ

Textos base, usados en las filas de la tabla 5 (AUDITORIA §5, diferencias D1–D4):

- **[EV-D1]** Everpresent Λ (Ahmed et al. 2004; Zwane et al. 2018) is a different model from the one tested here: Λ has zero mean and fluctuates, with an amplitude fixed by causal-set counting of the past four-volume (δΛ ∼ V^{-1/2}), at the level of the critical density in every epoch and with no free amplitude. The residual X of this work is a small stationary perturbation around a non-zero mean Λ, with free σ_X and θ. The bound on σ_X does not, by itself, test Everpresent Λ.
- **[EV-D3]** The seed σ_0 ∼ 10⁻⁶¹ uses N ∼ 10¹²² from horizon-area (Bekenstein–Hawking) counting. With the four-volume counting of Everpresent Λ, δΛ tracks H² and the relative amplitude is of order one with no amplifier; the ∼10⁵⁶ gap is specific to the area-counting choice.
- **[EV-D4]** Zwane et al. (2018) find that Everpresent Λ fits the CMB, BAO and SN data of that time about as well as ΛCDM and eases the low-redshift H0 and Lyα BAO tensions. Tested against DESI DR2 BAO with r_d computed for each realization, small dark energy around recombination and ω_m within ±2σ of Planck, none of the realizations of its two published implementations fits as well as ΛCDM (0/20000 for the model of Zwane et al.; 0/500000 for the model of Das et al., α = 0.01 and 0.005; `analisis_cmb/everpresent/EVERPRESENT_vs_DESI.md`). This applies to those phenomenological implementations, not to the idea of Everpresent Λ in general.
- **[EV-D2]** véase la propuesta única D2.A (Axioma A2) + D2.B (paper unimodular), marcada como **decisión del autor**.

Fuentes de los textos base: [EV-D1] y [EV-D3], AUDITORIA §5 (resúmenes de arXiv de Ahmed et al. 2004 y Zwane et al. 2018) y EVERPRESENT_vs_DESI.md §0 y conclusión 6. [EV-D4]: Zwane et al., abstract (vía AUDITORIA §5 D4) y `analisis_cmb/everpresent/EVERPRESENT_vs_DESI.md` §3e y Conclusiones 2–3, 6 (`resultados/rd_real.json`): criterio final (c) — r_d calculado con la historia de cada realización, |Ω_de| < 0.25 en z ∈ [1000, 5000] y ω_m a ±1.6 % de Planck —: Modelo 2 de Zwane et al. 0/20000 (< 1.5×10⁻⁴), Modelo 1 de Das et al. 0/300000 (α = 0.01) y 0/200000 (α = 0.005), 0/500000 en total (< 1.0×10⁻⁵ y < 1.5×10⁻⁵). Limitaciones §4: implementaciones fenomenológicas, CMB por proxy.

### D2.A — Axioma A2 (1/2), DECISIÓN DEL AUTOR — `papers/stochastic-dark-energy-desi-dr2.md` l.58 (reemplazar)

**Texto actual:**

> In unimodular gravity, $\Lambda$ can be interpreted as a constant of integration conjugate to the spacetime four-volume. If spacetime consists of $N$ discrete elements with $N = V / L_P^4$, Poisson fluctuations yield

**Propuesto:**

> In unimodular gravity, Λ can be interpreted as a constant of integration conjugate to the spacetime four-volume. If spacetime consists of N discrete elements, Poisson fluctuations in N yield

**Fuente:** AUDITORIA §5 D3: l.58 define N = V/L_P⁴ (4-volumen) y l.66 usa N ∼ 10¹²² (área, l.54). Incoherencia interna.

### D2.A — Axioma A2 (2/2), DECISIÓN DEL AUTOR — `papers/stochastic-dark-energy-desi-dr2.md` l.66 (reemplazar)

**Texto actual:**

> With $N \sim 10^{122}$, this gives $\delta\Lambda \sim 10^{-61}$ in Planck units. This motivates the existence of a small, non-zero stochastic component, though I do not claim a full derivation of $\rho_\Lambda$ from this argument.

**Propuesto:**

> Two countings give very different amplitudes, and I keep them apart. Counting the elements of the past four-volume, N = V/L_P⁴, gives δΛ ∼ V^{-1/2} ∼ H²: fluctuations of the order of the critical density in every epoch. This is Everpresent Λ (Ahmed et al. 2004; Sorkin 2007; Zwane et al. 2018), and it is not the model tested here. Counting horizon-area degrees of freedom, N ∼ 10¹²², gives the much smaller relative seed σ_0 ∼ 1/√N ∼ 10⁻⁶¹ that I use below only as a motivational UV scale. Whether either applies is an empirical question: the order-one version, in its two published implementations, does not fit DESI DR2 BAO together with the CMB (0/20000 and 0/500000 realizations; `analisis_cmb/everpresent/EVERPRESENT_vs_DESI.md`), and a small stationary version is bounded here by σ_X < 2.6×10⁻². I do not claim a full derivation of ρ_Λ from either counting.

**Fuente:** AUDITORIA §5 D3; Ahmed et al. 2004, Sorkin 2007, Zwane et al. 2018 (resúmenes de arXiv); `analisis_cmb/everpresent/EVERPRESENT_vs_DESI.md` §3e y Conclusiones 2–3, 6 (`resultados/rd_real.json`): criterio final (c) — r_d calculado con la historia de cada realización, |Ω_de| < 0.25 en z ∈ [1000, 5000] y ω_m a ±1.6 % de Planck —: Modelo 2 de Zwane et al. 0/20000 (< 1.5×10⁻⁴), Modelo 1 de Das et al. 0/300000 (α = 0.01) y 0/200000 (α = 0.005), 0/500000 en total (< 1.0×10⁻⁵ y < 1.5×10⁻⁵). Limitaciones §4: implementaciones fenomenológicas, CMB por proxy. Cota σ_X: `results/profile_sigma_x/profile_sigma_x.json`.

### D2.B — paper unimodular (1/3), DECISIÓN DEL AUTOR — `theory/unimodular-gravity-vacuum-smoothness.md` l.22 (reemplazar (solo la frase central))

**Texto actual:**

> which is a consistency condition, not a predictive equation. The value of $\Lambda$ remains a free global parameter, but its fluctuations are rigorously zero. I conclude that the observed smoothness of the vacuum is not an accident; it is a geometric necessity that singles out unimodular gravity as the only known structure compatible with the data. Future surveys, such as Euclid DR1, will provide the decisive test.

**Propuesto:**

> "…but its fluctuations are rigorously zero. I conclude that the observed smoothness of the vacuum is not an accident; it is a geometric necessity that singles out unimodular gravity as the only known structure compatible with the data." → "…but its local fluctuations of the form $V(x)g_{\mu\nu}$ are cancelled identically. The same framework with causal-set discreteness predicts order-one fluctuations of the global $\Lambda$ (Everpresent Λ, Ahmed et al. 2004); which case holds is decided by data (Section 7)."

**Fuente:** AUDITORIA §5 D2 (mismo punto de partida, conclusión opuesta); tabla 5, fila unimodular l.22.

### D2.B — paper unimodular (2/3), DECISIÓN DEL AUTOR — `theory/unimodular-gravity-vacuum-smoothness.md` l.241 (insertar tras)

**Línea tras la que se inserta:**

> The smallness of $\Lambda$($\sim 10^{-122}$ in Planck units) is not a fine-tuning problem in this framework. It is a consequence of the vastness of the universe: the denominator $V_4$ is enormous, so the average of any finite quantity over it tends to zero. The discrepancy between the quantum prediction and the observed value is not an inconsistency; it is the signature that the quantum fluctuations are local, while $\Lambda$ is global.

**Propuesto:**

> The same four-volume, built from discrete elements, is the origin of the Everpresent Λ fluctuation δΛ ∼ V_4^{-1/2} (Ahmed et al. 2004; Sorkin 2007); see the end of Section 7.

**Fuente:** AUDITORIA §5, tabla 5, fila unimodular l.246 (hoy l.241): "es el mismo argumento de Ahmed et al., con la conclusión contraria".

### D2.B — paper unimodular (3/3), DECISIÓN DEL AUTOR — `theory/unimodular-gravity-vacuum-smoothness.md` l.269 (insertar tras (cierre de §7, tras la ecuación σ_X = 0))

**Línea tras la que se inserta:**

> \sigma_X = 0.

**Propuesto:**

> **Relation to Everpresent Λ (author's decision).** Ahmed, Dodelson, Greene and Sorkin (2004) start from the same unimodular framework — Λ conjugate to the four-volume V_4 — and, adding causal-set discreteness, reach the opposite conclusion: Λ fluctuates, δΛ ∼ V_4^{-1/2}, at the level of the critical density in every epoch (Everpresent Λ; Zwane et al. 2018). The two statements concern different objects. The cancellation in this section acts on the local term V(x) g_μν of the trace-free equations; it says nothing about the global Λ conjugate to V_4, which is the object whose fluctuations Everpresent Λ predicts. Which picture describes our universe is decided empirically, not by this algebra. The order-one version, in its two published implementations, does not fit DESI DR2 BAO together with the CMB: once r_d is computed for each realization, dark energy is small around recombination and ω_m is within ±2σ of Planck, no realization fits as well as ΛCDM (0/20000 and 0/500000; `analisis_cmb/everpresent/EVERPRESENT_vs_DESI.md`). A small stationary version is bounded by σ_X < 2.6×10⁻². In this sense σ_X = 0 is the prediction of the local cancellation for local fluctuations, not a theorem about Λ.

**Fuente:** AUDITORIA §5 D2/D4 y tabla 5 (unimodular l.275, hoy l.269); `analisis_cmb/everpresent/EVERPRESENT_vs_DESI.md` §3e y Conclusiones 2–3, 6 (`resultados/rd_real.json`): criterio final (c) — r_d calculado con la historia de cada realización, |Ω_de| < 0.25 en z ∈ [1000, 5000] y ω_m a ±1.6 % de Planck —: Modelo 2 de Zwane et al. 0/20000 (< 1.5×10⁻⁴), Modelo 1 de Das et al. 0/300000 (α = 0.01) y 0/200000 (α = 0.005), 0/500000 en total (< 1.0×10⁻⁵ y < 1.5×10⁻⁵). Limitaciones §4: implementaciones fenomenológicas, CMB por proxy. Cota σ_X: `results/profile_sigma_x/profile_sigma_x.json`.

### D2.1 — `manuscript/PREPRINT.md` l.41 (insertar tras)

**Línea tras la que se inserta:**

> A finite information bound on the observable universe ($N\sim 10^{122}$ Planck units) has long motivated order-of-magnitude vacuum fluctuations $\delta\Lambda\sim 1/\sqrt{N}$ (Sorkin-type arguments in unimodular / causal-set settings). That seed is many orders of magnitude below BAO precision. The programme I follow is therefore:

**Propuesto:**

> These arguments come in two versions that should not be confused (Ahmed et al. 2004; Sorkin 2007; Zwane et al. 2018). The seed σ_0 ∼ 10⁻⁶¹ uses N ∼ 10¹²² from horizon-area (Bekenstein–Hawking) counting. With the four-volume counting of Everpresent Λ, δΛ tracks H² and the relative amplitude is of order one with no amplifier; the ∼10⁵⁶ gap is specific to the area-counting choice. Everpresent Λ (Ahmed et al. 2004; Zwane et al. 2018) is a different model from the one tested here: Λ has zero mean and fluctuates, with an amplitude fixed by causal-set counting of the past four-volume (δΛ ∼ V^{-1/2}), at the level of the critical density in every epoch and with no free amplitude. The residual X of this work is a small stationary perturbation around a non-zero mean Λ, with free σ_X and θ. The bound on σ_X does not, by itself, test Everpresent Λ.

**Fuente:** AUDITORIA §5 tabla 5 (PREPRINT l.41); [EV-D1], [EV-D3].

### D2.2 — `manuscript/PREPRINT.md` l.26 (reemplazar (solo el inciso))

**Texto actual:**

> I then quantify an amplification bottleneck: a Sorkin–Bekenstein Poisson seed $\sigma_0\sim 10^{-61}$ lies $\sim 10^{56}$ below Euclid-scale residual amplitudes ($\sim 10^{-5}$) and $\sim 10^{59}$ below the DESI residual ceiling measured here ($2.6\times 10^{-2}$). Audited linear amplifiers are excluded, short by many decades, or redefine the seed. Machine scans of three allowed non-linear cards show that late-time freeze-out alone yields unit gain, a soft double-well avalanche yields gain $\sim 2$, and only a local effective count $N_{\mathrm{eff}}\sim 4\times 10^{7}$ (DESI) / $10^{10}$ (Euclid) reaches target amplitudes—by redefining the microphysical seed, not by free amplification.

**Propuesto:**

> "a Sorkin–Bekenstein Poisson seed $\sigma_0\sim 10^{-61}$" → "a Sorkin–Bekenstein Poisson seed with horizon-area counting, $\sigma_0\sim 10^{-61}$ (four-volume counting gives order-one fluctuations instead — Everpresent Λ, a different model),".

**Fuente:** AUDITORIA §5 tabla 5 (PREPRINT l.26); [EV-D3].

### D2.3 — `manuscript/PREPRINT.md` l.198 (insertar tras (cierre de §5.1))

**Línea tras la que se inserta:**

> DESI constrains the **effective** residual amplitude after any amplification—not $\sigma_0$, and not a new value of $\rho_\Lambda$. Without a physical map $\sigma_0\to A_0$, “Euclid will see vacuum noise” is not a theorem of this analysis.

**Propuesto:**

> The seed σ_0 ∼ 10⁻⁶¹ uses N ∼ 10¹²² from horizon-area (Bekenstein–Hawking) counting. With the four-volume counting of Everpresent Λ, δΛ tracks H² and the relative amplitude is of order one with no amplifier; the ∼10⁵⁶ gap is specific to the area-counting choice. Zwane et al. (2018) find that Everpresent Λ fits the CMB, BAO and SN data of that time about as well as ΛCDM and eases the low-redshift H0 and Lyα BAO tensions. Tested against DESI DR2 BAO with r_d computed for each realization, small dark energy around recombination and ω_m within ±2σ of Planck, none of the realizations of its two published implementations fits as well as ΛCDM (0/20000 for the model of Zwane et al.; 0/500000 for the model of Das et al., α = 0.01 and 0.005; `analisis_cmb/everpresent/EVERPRESENT_vs_DESI.md`). This applies to those phenomenological implementations, not to the idea of Everpresent Λ in general.

**Fuente:** AUDITORIA §5 tabla 5 (PREPRINT l.198); [EV-D3], [EV-D4]; `analisis_cmb/everpresent/EVERPRESENT_vs_DESI.md` §3e y Conclusiones 2–3, 6 (`resultados/rd_real.json`): criterio final (c) — r_d calculado con la historia de cada realización, |Ω_de| < 0.25 en z ∈ [1000, 5000] y ω_m a ±1.6 % de Planck —: Modelo 2 de Zwane et al. 0/20000 (< 1.5×10⁻⁴), Modelo 1 de Das et al. 0/300000 (α = 0.01) y 0/200000 (α = 0.005), 0/500000 en total (< 1.0×10⁻⁵ y < 1.5×10⁻⁵). Limitaciones §4: implementaciones fenomenológicas, CMB por proxy.

### D2.4 — `manuscript/PREPRINT.md` l.230 (insertar tras)

**Línea tras la que se inserta:**

> Volume-preserving diffeomorphisms / unimodular structure project out isotropic vacuum stress of the form $T_{\mu\nu}=V(x)g_{\mu\nu}$. That is a candidate explanation of late-time residual smoothness, not a particle detection and not a substitute for the BAO likelihood.

**Propuesto:**

> The same unimodular starting point, combined with causal-set discreteness, predicts order-one fluctuations of the global Λ (Everpresent Λ; Ahmed et al. 2004). The projection above removes local terms of the form $V(x)g_{\mu\nu}$; it does not act on the global Λ conjugate to the four-volume.

**Fuente:** AUDITORIA §5 D2, tabla 5 (PREPRINT l.228, hoy l.230).

### D2.5 — `manuscript/PREPRINT.md` l.266 (reemplazar)

**Texto actual:**

> 5. **Sorkin seed:** used as motivational UV scale for the amplification discussion, not as a fitted parameter.

**Propuesto:**

> 5. **Sorkin seed:** used as motivational UV scale (horizon-area counting) for the amplification discussion, not as a fitted parameter. Everpresent Λ (Ahmed et al. 2004; Zwane et al. 2018) is a different model from the one tested here: Λ has zero mean and fluctuates, with an amplitude fixed by causal-set counting of the past four-volume (δΛ ∼ V^{-1/2}), at the level of the critical density in every epoch and with no free amplitude. The residual X of this work is a small stationary perturbation around a non-zero mean Λ, with free σ_X and θ. The bound on σ_X does not, by itself, test Everpresent Λ. Zwane et al. (2018) find that Everpresent Λ fits the CMB, BAO and SN data of that time about as well as ΛCDM and eases the low-redshift H0 and Lyα BAO tensions. Tested against DESI DR2 BAO with r_d computed for each realization, small dark energy around recombination and ω_m within ±2σ of Planck, none of the realizations of its two published implementations fits as well as ΛCDM (0/20000 for the model of Zwane et al.; 0/500000 for the model of Das et al., α = 0.01 and 0.005; `analisis_cmb/everpresent/EVERPRESENT_vs_DESI.md`). This applies to those phenomenological implementations, not to the idea of Everpresent Λ in general.

**Fuente:** AUDITORIA §5 tabla 5 (PREPRINT l.266); [EV-D1], [EV-D4].

### D2.6 — `manuscript/PREPRINT.md` l.350 (reemplazar)

**Texto actual:**

> 6. Sorkin, causal-set / unimodular vacuum fluctuation arguments (see repository notes for curated pointers).

**Propuesto:**

> 6. Ahmed, M., Dodelson, S., Greene, P. B., & Sorkin, R. D., "Everpresent Λ," Phys. Rev. D 69, 103523 (2004), arXiv:astro-ph/0209274.
> 6b. Zwane, N., Afshordi, N., & Sorkin, R. D., "Cosmological tests of Everpresent Λ," Class. Quantum Grav. (2018), arXiv:1703.06265.
> 6c. Sorkin, R. D., "Is the cosmological 'constant' a nonlocal quantum residue of discreteness of the causal set type?," AIP Conf. Proc. 957, 142 (2007), arXiv:0710.1675.
> (Renumber the following references.)

**Fuente:** AUDITORIA.md §5 (referencias comprobadas en arXiv el 2026-10-03) y `analisis_cmb/everpresent/LECTURA_EVERPRESENT.md` (Das et al. I y II, Barrow 2007).

### D2.7 — `manuscript/CLAIMS.md` l.39 (insertar tras)

**Línea tras la que se inserta:**

> | **N10** | Residual ceiling \(\sigma_X\) is not the baryon density \(\omega_b\); VLBI / photon sphere are not residual claims. See [`papers/scope-and-mixups.md`]\(../papers/scope-and-mixups.md) |

**Propuesto:**

> | **N11** | The $\sigma_X$ bound does not by itself test Everpresent Λ (Ahmed et al. 2004; Zwane et al. 2018), a zero-mean Λ with order-one fluctuations; that model is tested separately in `analisis_cmb/everpresent/` |
>
> (Note: D3 also proposes N11/N12 for the CMB; final numbering is the author's choice — see D3.6.)

**Fuente:** AUDITORIA §5 tabla 5 (CLAIMS l.39); [EV-D1].

### D2.8 — `papers/stochastic-dark-energy-desi-dr2.md` l.34 (insertar tras (cita))

**Línea tras la que se inserta:**

> This work proposes one such extension. The starting point is the Bekenstein–Hawking entropy bound, which implies that the observable universe has a finite number of effective degrees of freedom ($N \sim 10^{122}$ in Planck units). In unimodular gravity, the cosmological constant $\Lambda$ appears as an integration constant conjugate to spacetime four-volume. If spacetime is fundamentally discrete (as in causal-set theory), Poisson fluctuations in the number of elements $N$ induce residual fluctuations $\delta\Lambda \sim 1/\sqrt{N}$. This is not a derivation of the observed value of $\Lambda$, but it provides a plausible order-of-magnitude motivation for a small stochastic component.

**Propuesto:**

> Add after the sentence: "(For the causal-set version of this idea, see Ahmed et al. 2004 and Sorkin 2007; it predicts order-one fluctuations, Everpresent Λ, which is a different model from the small residual studied here — Axiom A2.)"

**Fuente:** AUDITORIA §5 tabla 5 (stochastic l.34); [EV-D1].

### D2.9 — `papers/stochastic-dark-energy-desi-dr2.md` l.380 (insertar tras (además de la reescritura D3.9))

**Línea tras la que se inserta:**

> ### 7.3 Implications for Quantum Gravity and Information-Theoretic Models

**Propuesto:**

> Everpresent Λ (Ahmed et al. 2004; Zwane et al. 2018) is a different model from the one tested here: Λ has zero mean and fluctuates, with an amplitude fixed by causal-set counting of the past four-volume (δΛ ∼ V^{-1/2}), at the level of the critical density in every epoch and with no free amplitude. The residual X of this work is a small stationary perturbation around a non-zero mean Λ, with free σ_X and θ. The bound on σ_X does not, by itself, test Everpresent Λ. Zwane et al. (2018) find that Everpresent Λ fits the CMB, BAO and SN data of that time about as well as ΛCDM and eases the low-redshift H0 and Lyα BAO tensions. Tested against DESI DR2 BAO with r_d computed for each realization, small dark energy around recombination and ω_m within ±2σ of Planck, none of the realizations of its two published implementations fits as well as ΛCDM (0/20000 for the model of Zwane et al.; 0/500000 for the model of Das et al., α = 0.01 and 0.005; `analisis_cmb/everpresent/EVERPRESENT_vs_DESI.md`). This applies to those phenomenological implementations, not to the idea of Everpresent Λ in general.

**Fuente:** AUDITORIA §5 tabla 5 (stochastic l.384); [EV-D1], [EV-D4].

### D2.10 — `papers/stochastic-dark-energy-desi-dr2.md` l.939 (reemplazar)

**Texto actual:**

> [4] Sorkin, R. D., "Is the Cosmological 'Constant' a Nonlocal Quantum Residual?," arXiv:gr-qc/0503057 (2005).

**Propuesto:**

> [4] Sorkin, R. D., "Is the cosmological 'constant' a nonlocal quantum residue of discreteness of the causal set type?," AIP Conf. Proc. 957, 142 (2007), arXiv:0710.1675.
> [4b] Ahmed, M., Dodelson, S., Greene, P. B., & Sorkin, R. D., "Everpresent Λ," Phys. Rev. D 69, 103523 (2004), arXiv:astro-ph/0209274.
> [4c] Zwane, N., Afshordi, N., & Sorkin, R. D., "Cosmological tests of Everpresent Λ," Class. Quantum Grav. (2018), arXiv:1703.06265.

**Fuente:** AUDITORIA §5, hallazgo colateral: gr-qc/0503057 es "Big extra dimensions make Λ too small" (Braz. J. Phys. 35, 280, 2005); el título citado es arXiv:0710.1675. AUDITORIA.md §5 (referencias comprobadas en arXiv el 2026-10-03) y `analisis_cmb/everpresent/LECTURA_EVERPRESENT.md` (Das et al. I y II, Barrow 2007).

### D2.11 — `papers/quantum-fluid-instabilities-desi-dr2.md` l.194 (reemplazar)

**Texto actual:**

> [4] Sorkin, R. D., "Is the Cosmological 'Constant' a Nonlocal Quantum Residual?," arXiv:gr-qc/0503057 (2005).

**Propuesto:**

> [4] Sorkin, R. D., "Is the cosmological 'constant' a nonlocal quantum residue of discreteness of the causal set type?," AIP Conf. Proc. 957, 142 (2007), arXiv:0710.1675.

**Fuente:** Ídem D2.10.

### D2.12 — `papers/quantum-fluid-instabilities-desi-dr2.md` l.112 (insertar tras)

**Línea tras la que se inserta:**

> with seed amplitude $\sigma_{0}\sim 10^{-61}$ (Bekenstein–Hawking / causal-set counting, same order used in the amplification-gap notes).

**Propuesto:**

> The seed σ_0 ∼ 10⁻⁶¹ uses N ∼ 10¹²² from horizon-area (Bekenstein–Hawking) counting. With the four-volume counting of Everpresent Λ, δΛ tracks H² and the relative amplitude is of order one with no amplifier; the ∼10⁵⁶ gap is specific to the area-counting choice.

**Fuente:** AUDITORIA §5 tabla 5 (quantum-fluid l.112); [EV-D3].

### D2.13 — `papers/amplification-gap.md` l.36 (insertar tras (bajo la tabla))

**Línea tras la que se inserta:**

> | Sorkin / Bekenstein seed | $\sigma_0 \sim 1.18\times 10^{-61}$ ($H_0=67.4$) | Motivational UV amplitude (holographic $d=2$) |

**Propuesto:**

> The seed σ_0 ∼ 10⁻⁶¹ uses N ∼ 10¹²² from horizon-area (Bekenstein–Hawking) counting. With the four-volume counting of Everpresent Λ, δΛ tracks H² and the relative amplitude is of order one with no amplifier; the ∼10⁵⁶ gap is specific to the area-counting choice. (Ahmed et al. 2004; Zwane et al. 2018.)

**Fuente:** AUDITORIA §5 tabla 5 (amplification-gap l.36–39); [EV-D3].

### D2.14 — `papers/amplification-gap.md` l.202 (insertar tras)

**Línea tras la que se inserta:**

> **Verdict R1:** The **only** card that reaches DESI/Euclid residual *amplitudes* without linear gain. **Price:** abandon global horizon counting; justify a meso-scale causal-set / correlation volume. That is a **new theory claim**, not a free amplification on $\sigma_0\sim 10^{-61}$.

**Propuesto:**

> The published version of this redefinition is Everpresent Λ (Ahmed et al. 2004; Zwane et al. 2018): four-volume counting, order-one fluctuations. Zwane et al. (2018) find that Everpresent Λ fits the CMB, BAO and SN data of that time about as well as ΛCDM and eases the low-redshift H0 and Lyα BAO tensions. Tested against DESI DR2 BAO with r_d computed for each realization, small dark energy around recombination and ω_m within ±2σ of Planck, none of the realizations of its two published implementations fits as well as ΛCDM (0/20000 for the model of Zwane et al.; 0/500000 for the model of Das et al., α = 0.01 and 0.005; `analisis_cmb/everpresent/EVERPRESENT_vs_DESI.md`). This applies to those phenomenological implementations, not to the idea of Everpresent Λ in general.

**Fuente:** AUDITORIA §5 tabla 5 (amplification-gap l.202); [EV-D4].

### D2.15 — `papers/amplification-gap.md` l.279 (insertar tras)

**Línea tras la que se inserta:**

> \sigma_0\sim\frac{1}{\sqrt{N}},\quad N\sim 10^{122}.

**Propuesto:**

> (Horizon-area counting. Four-volume counting gives order-one fluctuations — Everpresent Λ, Ahmed et al. 2004.)

**Fuente:** AUDITORIA §5 tabla 5 (amplification-gap l.276–281); [EV-D3].

### D2.16 — `papers/scope-and-mixups.md` l.52 (insertar tras)

**Línea tras la que se inserta:**

> | \(\sigma_0 \sim 10^{-61}\) | Motivational Sorkin / holographic seed |

**Propuesto:**

> | \(\delta\Lambda/\Lambda \sim \mathcal{O}(1)\) | Everpresent Λ (Ahmed et al. 2004): not \(\sigma_X\), not \(\sigma_0\) |

**Fuente:** AUDITORIA §5 tabla 5 (scope-and-mixups l.52).

### D2.17 — `papers/scope-and-mixups.md` l.54 (insertar tras)

**Línea tras la que se inserta:**

> \(\sigma_X\) in the likelihood is the residual amplitude. It is not \(\omega_b\). It is not the Sorkin seed unless someone builds a map. In the OU SDE, the diffusion coefficient \(\sigma\) and the reported \(\sigma_X\) only match after \(\sigma_X = \sigma/\sqrt{2\theta}\); the three-gate number \(\sigma_{\mathrm{free}}\) is already compared to \(\sigma_X\) as an *amplitude*, so I do not divide it by \(\sqrt{2\theta}\) again.

**Propuesto:**

> Everpresent Λ (Ahmed et al. 2004; Zwane et al. 2018) is a different model from the one tested here: Λ has zero mean and fluctuates, with an amplitude fixed by causal-set counting of the past four-volume (δΛ ∼ V^{-1/2}), at the level of the critical density in every epoch and with no free amplitude. The residual X of this work is a small stationary perturbation around a non-zero mean Λ, with free σ_X and θ. The bound on σ_X does not, by itself, test Everpresent Λ.

**Fuente:** AUDITORIA §5 tabla 5 (scope-and-mixups l.54); [EV-D1].

### D2.18 — `theory/principle-of-vacuum-smoothness.md` l.42 (matizar (añadir al final del punto))

**Texto actual:**

> - **For quantum gravity and emergent spacetime models:** Mechanisms that predict residual fluctuations in $\Lambda$ (for example causal set theory or the Sorkin mechanism) need a strong suppression mechanism at late times.

**Propuesto:**

> Append: "This applies to small residuals around a mean Λ. Everpresent Λ (Ahmed et al. 2004) is a zero-mean, order-one model; Zwane et al. (2018) found that it fits the data of that time about as well as ΛCDM, while its published implementations do not fit DESI DR2 BAO together with the CMB (`analisis_cmb/everpresent/`)."

**Fuente:** AUDITORIA §5 tabla 5 (principle l.42); [EV-D4]; `analisis_cmb/everpresent/EVERPRESENT_vs_DESI.md` §3e y Conclusiones 2–3, 6 (`resultados/rd_real.json`): criterio final (c) — r_d calculado con la historia de cada realización, |Ω_de| < 0.25 en z ∈ [1000, 5000] y ω_m a ±1.6 % de Planck —: Modelo 2 de Zwane et al. 0/20000 (< 1.5×10⁻⁴), Modelo 1 de Das et al. 0/300000 (α = 0.01) y 0/200000 (α = 0.005), 0/500000 en total (< 1.0×10⁻⁵ y < 1.5×10⁻⁵). Limitaciones §4: implementaciones fenomenológicas, CMB por proxy.

### D2.19 — `theory/sdiff-fundamental-vs-emergent.md` l.108 (insertar tras)

**Línea tras la que se inserta:**

> 3. **Structure of spacetime:** Do the DESI limits on $\sigma_X$ rule out models with significant local metric fluctuations (for example spacetime foam or discrete causal sets with large variance)?

**Propuesto:**

> Not directly. Everpresent Λ (Ahmed et al. 2004; Zwane et al. 2018) is a different model from the one tested here: Λ has zero mean and fluctuates, with an amplitude fixed by causal-set counting of the past four-volume (δΛ ∼ V^{-1/2}), at the level of the critical density in every epoch and with no free amplitude. The residual X of this work is a small stationary perturbation around a non-zero mean Λ, with free σ_X and θ. The bound on σ_X does not, by itself, test Everpresent Λ. Zwane et al. (2018) find that Everpresent Λ fits the CMB, BAO and SN data of that time about as well as ΛCDM and eases the low-redshift H0 and Lyα BAO tensions. Tested against DESI DR2 BAO with r_d computed for each realization, small dark energy around recombination and ω_m within ±2σ of Planck, none of the realizations of its two published implementations fits as well as ΛCDM (0/20000 for the model of Zwane et al.; 0/500000 for the model of Das et al., α = 0.01 and 0.005; `analisis_cmb/everpresent/EVERPRESENT_vs_DESI.md`). This applies to those phenomenological implementations, not to the idea of Everpresent Λ in general.

**Fuente:** AUDITORIA §5 tabla 5 (sdiff l.108); [EV-D1], [EV-D4].

### D2.20 — `theory/fundamental-vs-emergent-vacuum-relaxation.md` l.101 (insertar tras)

**Línea tras la que se inserta:**

> DESI therefore constrains the **effective** late-time amplitude after amplification or non-Sorkin coupling — not $\sigma_0$ itself.

**Propuesto:**

> "Invisible" holds for the horizon-area seed only. The seed σ_0 ∼ 10⁻⁶¹ uses N ∼ 10¹²² from horizon-area (Bekenstein–Hawking) counting. With the four-volume counting of Everpresent Λ, δΛ tracks H² and the relative amplitude is of order one with no amplifier; the ∼10⁵⁶ gap is specific to the area-counting choice.

**Fuente:** AUDITORIA §5 tabla 5 (fundamental-vs-emergent l.88–101); [EV-D3].

### D2.21 — `theory/fundamental-vs-emergent-vacuum-relaxation.md` l.194 (reemplazar)

**Texto actual:**

> [5] Sorkin, R. D., arXiv:gr-qc/0503057 (2005).

**Propuesto:**

> [5] Sorkin, R. D., "Is the cosmological 'constant' a nonlocal quantum residue of discreteness of the causal set type?," AIP Conf. Proc. 957, 142 (2007), arXiv:0710.1675.

**Fuente:** Ídem D2.10.

### D2.22 — `exploratory/notes/desqueezing-relaxation-vacuum-fluctuations-note.md` l.157 (insertar tras (al final del punto A2))

**Línea tras la que se inserta:**

> - **A2 (Sorkin seed):** discreteness induces Poisson fluctuations $\delta\Lambda\sim 1/\sqrt{N}$, i.e.

**Propuesto:**

> The seed σ_0 ∼ 10⁻⁶¹ uses N ∼ 10¹²² from horizon-area (Bekenstein–Hawking) counting. With the four-volume counting of Everpresent Λ, δΛ tracks H² and the relative amplitude is of order one with no amplifier; the ∼10⁵⁶ gap is specific to the area-counting choice. (Ahmed et al. 2004; Zwane et al. 2018.)

**Fuente:** AUDITORIA §5 tabla 5 (nota l.161, hoy l.157); [EV-D3].

### D2.23 — `exploratory/notes/desqueezing-relaxation-vacuum-fluctuations-note.md` l.426 (reemplazar (solo el inciso))

**Texto actual:**

> 4. With the repository seed $\sigma_0\sim 10^{-61}$, residuals stay unobservable: **Euclid cannot see pure Sorkin noise**.

**Propuesto:**

> "**Euclid cannot see pure Sorkin noise**" → "**Euclid cannot see a pure horizon-area Sorkin seed**". The same change on l.480 ("Sorkin’s $10^{-61}$" → "the horizon-area $10^{-61}$").

**Fuente:** AUDITORIA §5 tabla 5 (nota l.331/480, hoy l.426/480); [EV-D3].

### D2.24 — `exploratory/notes/desqueezing-relaxation-vacuum-fluctuations-note.md` l.449 (reemplazar)

**Texto actual:**

> [8] Sorkin, R. D., arXiv:gr-qc/0503057 (2005).

**Propuesto:**

> [8] Sorkin, R. D., "Is the cosmological 'constant' a nonlocal quantum residue of discreteness of the causal set type?," AIP Conf. Proc. 957, 142 (2007), arXiv:0710.1675.

**Fuente:** Ídem D2.10.

### D2.25 — `theory/euclid-protocol-vacuum-relaxation.md` l.22 (insertar tras (bajo la tabla))

**Línea tras la que se inserta:**

> | Microscopic seed | Sorkin / Bekenstein–Hawking $\sigma_0\sim 10^{-61}$ | Always $\ll$ BAO noise |

**Propuesto:**

> Everpresent Λ (Ahmed et al. 2004; Zwane et al. 2018) is a different model from the one tested here: Λ has zero mean and fluctuates, with an amplitude fixed by causal-set counting of the past four-volume (δΛ ∼ V^{-1/2}), at the level of the critical density in every epoch and with no free amplitude. The residual X of this work is a small stationary perturbation around a non-zero mean Λ, with free σ_X and θ. The bound on σ_X does not, by itself, test Everpresent Λ. The seed σ_0 ∼ 10⁻⁶¹ uses N ∼ 10¹²² from horizon-area (Bekenstein–Hawking) counting. With the four-volume counting of Everpresent Λ, δΛ tracks H² and the relative amplitude is of order one with no amplifier; the ∼10⁵⁶ gap is specific to the area-counting choice.

**Fuente:** AUDITORIA §5 tabla 5 (euclid-protocol l.22–25); [EV-D1], [EV-D3].

### D2.26 — `exploratory/notes/data-pack-option0-internet.md` l.16 (insertar tras)

**Línea tras la que se inserta:**

> | Motivational seed | $\sigma_0\sim 10^{-61}\sim 1/\sqrt{N}$, $N\sim 10^{122}$ | Bekenstein–Hawking / Sorkin (motivational, not a detection) |

**Propuesto:**

> | Everpresent Λ (four-volume counting, order-one) | Ahmed et al., arXiv:astro-ph/0209274; Zwane et al., arXiv:1703.06265 | Different model; tested in `analisis_cmb/everpresent/` |

**Fuente:** AUDITORIA §5 tabla 5 (data-pack l.16); AUDITORIA.md §5 (referencias comprobadas en arXiv el 2026-10-03) y `analisis_cmb/everpresent/LECTURA_EVERPRESENT.md` (Das et al. I y II, Barrow 2007).

### D2.27 — `README.md` l.35 (matizar (solo el inciso))

**Texto actual:**

> Best fit still goes to \(\sigma_X\to 0\). A coherent tachyonic growing mode is excluded when it is active (rank-1 covariance). Soft amplification of a pure Sorkin seed \(\sim 10^{-61}\) to Euclid-scale residuals needs \(\sim 10^{56}\) (to the DESI ceiling, \(\sim 10^{59}\)) and is not free — see the amplification note.

**Propuesto:**

> "a pure Sorkin seed \(\sim 10^{-61}\)" → "a pure horizon-area Sorkin seed \(\sim 10^{-61}\) (not Everpresent Λ, which is a different, order-one model)".

**Fuente:** AUDITORIA §5 tabla 5 (README l.35); [EV-D3].

### D2.28 — `papers/HONEST_HEADLINES.md` l.19 (insertar tras)

**Línea tras la que se inserta:**

> | Residual $\sigma_X$ sold as baryon density, or as an EHT photon-ring result | Different quantities; see `scope-and-mixups.md` |

**Propuesto:**

> | “Causal sets / Sorkin are excluded” | The σ_X bound tests a small residual, not Everpresent Λ; the published Everpresent implementations are tested separately (`analisis_cmb/everpresent/`) and only under their own assumptions. |

**Fuente:** AUDITORIA §5 tabla 5 (HONEST_HEADLINES tras l.19); [EV-D1], [EV-D4].

### D2.29 — Etiquetas "Sorkin seed" en anisotropic-slip, VERDICT.md, mapping_summary.txt, mapping_tables.md y scripts (opcional)

**Propuesto:**

> Optional: "Sorkin seed" → "Sorkin seed (area count)" where it labels σ_0 ∼ 10⁻⁶¹. Not itemized here; AUDITORIA marks it optional.

**Fuente:** AUDITORIA §5 tabla 5, última fila.

## D3. Resultados nuevos: R1 (CMB) y R2 (oscilador)

R1 = `analisis_cmb/RESULTADO.md` (números a σ_X = 2.633e-2; la malla, a 2.506e-2). R2 = `analisis_cmb/oscilador/resultados/tabla.md`. Fuentes comunes:

- R1: `analisis_cmb/resultados/takahashi/sigma_2.633e-02/resumen.json` y `tabla.md` (σ_X = 2.633e-2, Δ = 0.05, 200 realizaciones): θ = 1: Δχ²_TT p95 = 0.0668, Δχ²_φφ p95 = 0.0541; θ = 0.001: 0.00265 y 0.0023; fracción con Δχ² > 4 = 0; H0 p16/p50/p84 = 66.92/67.36/67.88 (θ = 1) → σ(H0) = 0.48 km/s/Mpc, y 67.34/67.36/67.38 (θ = 0.001) → 0.020; fuera de Planck ±2σ: 3.5 % (θ = 1), 0 (θ = 0.001).
- Malla: `analisis_cmb/resultados/takahashi/malla/malla.json` (`malla.md`; RESULTADO.md §3): 10 valores de σ_X de 2.506e-2 a 1e-4, θ ∈ {0.001, 1}; p95(Δχ²) < 4 en los 20 casos; p95 ∝ σ_X^n con n = 0.998–1.007; σ(H0) ∝ σ_X^m con m = 0.999–1.001.
- H0: RESULTADO.md §4(c): 73.04 − 67.36 = 5.68 km/s/Mpc (SH0ES arXiv:2112.04510; Planck arXiv:1807.06209 tabla 2). Con σ(H0) a 2.633e-2: 5.68/0.48 = 11.9 (θ = 1) y 5.68/0.020 = 284 (θ = 0.001).
- R2: `analisis_cmb/oscilador/resultados/tabla.md` (`oscilador.json`): 60 curvas X = A sin(2π ln a/P + φ), A ∈ {0.005, 0.01, 0.025}, P ∈ {0.1, 0.2, 0.5, 1, 2}, φ ∈ {0, π/2, π, 3π/2}; permitidas por DESI BAO (Δχ²_BAO < 4): 58/60; distinguibles en CMB: 0/60; máximo entre las permitidas Δχ²_TT = 0.070, Δχ²_φφ = 0.631; H0 resuelto entre 66.55 y 68.23.
- Everpresent: `analisis_cmb/everpresent/EVERPRESENT_vs_DESI.md` §3e y Conclusiones 2–3, 6 (`resultados/rd_real.json`): criterio final (c) — r_d calculado con la historia de cada realización, |Ω_de| < 0.25 en z ∈ [1000, 5000] y ω_m a ±1.6 % de Planck —: Modelo 2 de Zwane et al. 0/20000 (< 1.5×10⁻⁴), Modelo 1 de Das et al. 0/300000 (α = 0.01) y 0/200000 (α = 0.005), 0/500000 en total (< 1.0×10⁻⁵ y < 1.5×10⁻⁵). Limitaciones §4: implementaciones fenomenológicas, CMB por proxy.

### D3.1 (abstract) — `manuscript/PREPRINT.md` l.24 (insertar tras)

**Línea tras la que se inserta:**

> Separately, a coherent tachyonic growth model with the correct rank-1 covariance induced by a single growing mode is strongly disfavoured: $\Delta\ln\mathcal{L}\approx -12.7$ relative to $\Lambda$CDM ($\Delta\chi^2\approx +25$).

**Propuesto:**

> A CMB consistency check with Planck 2018 (TT at ℓ < 30 and lensing bandpowers, via CAMB with the OU history mapped to w(a), H0 solved to keep θ_* fixed) does not distinguish the OU residual from ΛCDM at the DESI ceiling: the 95th percentile of Δχ² is ≤ 0.067 (TT) and ≤ 0.054 (lensing) at σ_X = 2.633×10⁻², and no realization reaches Δχ² > 4. Its only imprint is a scatter in the inferred H0 of 0.48 km/s/Mpc (θ = 1), about 12 times smaller than the Planck–SH0ES difference, so it neither tightens σ_X nor explains the H0 tension. A deterministic harmonic modulation X = A sin(2π ln a/P + φ) with A ≤ 0.025 and P ∈ [0.1, 2] is allowed by DESI DR2 BAO in 58 of 60 cases (Δχ²_BAO < 4) and is indistinguishable from ΛCDM in the CMB in all 60.

**Fuente:** R1, R2, H0 (arriba).

### D3.2 (guía de lectura) — `manuscript/PREPRINT.md` l.62 (reemplazar)

**Texto actual:**

> | 3 | Stationary OU / QNM null and $\sigma_X$ limit |

**Propuesto:**

> | 3 | Stationary OU / QNM null, $\sigma_X$ limit and CMB consistency (§3.5) |

**Fuente:** Navegación (AUDITORIA §6, PREPRINT l.59–68).

### D3.3 (método) — `manuscript/PREPRINT.md` l.78 (insertar tras (además de D5.1))

**Línea tras la que se inserta:**

> **Important separation:** DESI analyses that combine BAO with CMB+SN often prefer dynamical $w_0$–$w_a$. My BAO-only residual analysis uses a CPL background near $\Lambda$CDM ($ w_0\approx -1.001$, $w_a\approx +0.014$; `results/eos_cpl_desi_dr2/`, projected covariance, figure6 α) as the smooth reference for $\alpha=1$ under the fiducial cosmology of the published $\alpha$ products. The stochastic question is about *residuals*, not a re-fit of the full multi-probe DE equation of state.

**Propuesto:**

> **Secondary dataset (Section 3.5 only).** Planck 2018 TT at 2 ≤ ℓ < 30 and the Planck lensing bandpowers (L = 8–400), with cosmological parameters fixed to Planck 2018 except H0, which is solved for each curve to keep θ_* at its ΛCDM value. It is used as a consistency check, not in the σ_X likelihood.

**Fuente:** RESULTADO.md §1; R1.

### D3.4 (§3.5 nueva) — `manuscript/PREPRINT.md` l.162 (insertar tras (antes del `---` de cierre de §3))

**Línea tras la que se inserta:**

> > Under the phenomenological OU/QNM residual kernels and public DESI DR2 BAO summary statistics, there is no statistical evidence for stationary stochastic dark-energy fluctuations. The data prefer smooth evolution; $\sigma_X$ is bounded.

**Propuesto:**

> ### 3.5 CMB consistency (Planck 2018)
>
> I map the OU residual onto a single expansion history, δ ln ρ_DE = X/Ω_Λ0 with X smoothed on Δ = 0.05 in ln a, convert it to w(a) and run CAMB (PPF) with ω_b, ω_c and Σm_ν fixed to Planck 2018 and H0 solved so that 100θ_* = 1.041198. For 200 realizations at σ_X = 2.633×10⁻²:
>
> | θ | Δχ²_TT (ℓ<30) p50 / p95 | Δχ²_φφ p50 / p95 | fraction Δχ² > 4 | H0 p16 / p50 / p84 [km/s/Mpc] |
> |:--|:--|:--|:--|:--|
> | 0.001 | −1.5×10⁻⁴ / 2.7×10⁻³ | −1.2×10⁻⁴ / 2.3×10⁻³ | 0 | 67.34 / 67.36 / 67.38 |
> | 1 | −2.5×10⁻³ / 0.067 | −2.8×10⁻³ / 0.054 | 0 | 66.92 / 67.36 / 67.88 |
>
> On a grid of ten values from 2.506×10⁻² down to 10⁻⁴ the 95th percentile stays below 4 everywhere and scales linearly with σ_X, so Planck does not bound σ_X below the DESI ceiling. The only CMB imprint is a scatter in H0, σ(H0) = 0.48 km/s/Mpc at θ = 1 and 0.020 at θ = 0.001, i.e. 12 and 280 times smaller than the Planck–SH0ES difference (5.68 km/s/Mpc): the residual does not explain the H0 tension.
>
> A deterministic harmonic residual, X = A sin(2π ln a/P + φ) with A ∈ {0.005, 0.01, 0.025}, P ∈ {0.1, …, 2} and four phases, is allowed by DESI DR2 BAO (13 distances, full covariance, Δχ²_BAO < 4) in 58 of 60 cases, and none of the 60 is distinguishable in the CMB (largest Δχ² among the allowed: 0.070 in TT, 0.63 in lensing).
>
> Figure and data: `analisis_cmb/RESULTADO.md`, `analisis_cmb/oscilador/`.

**Fuente:** R1, malla, H0, R2 (arriba).

### D3.5 (§6.2) — `manuscript/PREPRINT.md` l.256 (matizar (añadir al final))

**Texto actual:**

> I do not implement a homemade Boltzmann hierarchy. Future work, if amplitudes become reachable, should use community codes (e.g. hi_class / MGCAMB) with a pre-specified target.

**Propuesto:**

> Append: "(The CMB consistency check of §3.5 uses CAMB with the residual mapped to a tabulated w(a) in the PPF module; it is not a Boltzmann treatment of anisotropic stress.)"

**Fuente:** AUDITORIA §6 (PREPRINT l.256); RESULTADO.md §1.

### D3.6 (CLAIMS) — `manuscript/CLAIMS.md` l.22 (insertar tras)

**Línea tras la que se inserta:**

> | **C7** | Even $\varepsilon=1$ anisotropic fraction of the BAO-bounded residual gives $\|\gamma-1\|\lesssim 5\times 10^{-2}$ ($z=0.5$, $\delta_m=1$), at best level with the forecast floor $\sim 0.05$ | `scripts/slip_bridge.py` |

**Propuesto:**

> | **C8** | Planck 2018 (TT ℓ<30 + lensing) does not distinguish the OU residual from ΛCDM at the DESI ceiling: p95(Δχ²) ≤ 0.067 at σ_X = 2.633×10⁻², and < 4 on the whole grid down to 10⁻⁴; the CMB adds no bound on σ_X | `analisis_cmb/RESULTADO.md` |
> | **C9** | A harmonic residual X = A sin(2π ln a/P + φ) with A ≤ 0.025, P ∈ [0.1, 2] is allowed by DESI DR2 BAO in 58/60 cases and is CMB-indistinguishable in 60/60 | `analisis_cmb/oscilador/resultados/tabla.md` |

**Fuente:** R1, malla, R2.

### D3.7 (CLAIMS, no-claims) — `manuscript/CLAIMS.md` l.39 (insertar tras (tras la N11 de D2.7, si se acepta))

**Línea tras la que se inserta:**

> | **N10** | Residual ceiling \(\sigma_X\) is not the baryon density \(\omega_b\); VLBI / photon sphere are not residual claims. See [`papers/scope-and-mixups.md`]\(../papers/scope-and-mixups.md) |

**Propuesto:**

> | **N12** | The residual does not explain the H0 tension: its CMB imprint is σ(H0) ≈ 0.48 km/s/Mpc (θ = 1) ≪ 5.68 km/s/Mpc |
> | **N13** | No full CMB MCMC: parameters other than H0 are fixed to Planck 2018 in the CMB check |
>
> (If D2.7 is not accepted, these become N11 and N12, as the instruction names them.)

**Fuente:** H0 (arriba); RESULTADO.md §6.

### D3.8 (CLAIMS, dataset) — `manuscript/CLAIMS.md` l.47 (insertar tras)

**Línea tras la que se inserta:**

> - **Out of primary claim:** full multi-probe DESI+CMB+SN dynamical DE preference (I cite it as context only).

**Propuesto:**

> - **Secondary check (in, not in the σ_X likelihood):** Planck 2018 TT ℓ < 30 + lensing bandpowers via CAMB (PPF), parameters fixed except H0 (`analisis_cmb/`).

**Fuente:** RESULTADO.md §1.

### D3.9 (§7.3, contradicha por R2) — `papers/stochastic-dark-energy-desi-dr2.md` l.384 (reemplazar (l.382–386: los tres párrafos de §7.3 se sustituyen por este))

**Texto actual:**

> While the localized $\sim 2.3\sigma$ tension observed in the DESI high-redshift tracers allows a theoretical baseline for dynamical dark energy, the simultaneous flat constraints imposed by the lower-redshift LRG bins severely penalize any oscillatory behavior. A cosmic vacuum modulated by horizon quasi-normal modes would induce harmonic variations in the expansion rate that are macroscopically ruled out by the continuity between the LRG data and the Lyman-$\alpha$ forest measurements.

**Propuesto:**

> ### 7.3 Implications for Quantum Gravity and Information-Theoretic Models
>
> The residual analysis bounds stationary stochastic variance around a smooth background; it does not single out a mechanism. In particular, it does not exclude oscillatory modulations of the expansion history at the amplitudes probed: a deterministic harmonic residual X = A sin(2π ln a/P + φ) with A ≤ 0.025 and P ∈ [0.1, 2] is allowed by DESI DR2 BAO in 58 of 60 cases (Δχ²_BAO < 4; Appendix J), and is indistinguishable from ΛCDM in the CMB. A vacuum modulated by horizon quasi-normal modes is therefore not ruled out at this level by the continuity between the LRG and Lyman-α data, and the smooth, dissipative OU closure is a modelling choice, not the only viable mechanism. What the data do show is that any such residual, stochastic or harmonic, must stay at the percent level in X.

**Fuente:** R2 (`analisis_cmb/oscilador/resultados/tabla.md` (`oscilador.json`): 60 curvas X = A sin(2π ln a/P + φ), A ∈ {0.005, 0.01, 0.025}, P ∈ {0.1, 0.2, 0.5, 1, 2}, φ ∈ {0, π/2, π, 3π/2}; permitidas por DESI BAO (Δχ²_BAO < 4): 58/60; distinguibles en CMB: 0/60; máximo entre las permitidas Δχ²_TT = 0.070, Δχ²_φφ = 0.631; H0 resuelto entre 66.55 y 68.23.). AUDITORIA §6: "Contradicho por R2". La afirmación de l.382 sobre "Pantheon+ … rigorous Bayesian assessment" no tiene soporte en el repo y se retira junto con el resto.

### D3.10 (abstract largo) — `papers/stochastic-dark-energy-desi-dr2.md` l.26 (insertar tras)

**Línea tras la que se inserta:**

> I caution that this result is subject to degeneracies between the stochastic component and the CPL parameters, which cannot be fully resolved with only 7 BAO bins. The primary contribution of this work is a validated analysis pipeline and a benchmark for future analyses with the $>20$ bins of Euclid DR1 (expected H2 2026).

**Propuesto:**

> A CMB consistency check with Planck 2018 (TT at ℓ < 30 and lensing bandpowers, via CAMB with the OU history mapped to w(a), H0 solved to keep θ_* fixed) does not distinguish the OU residual from ΛCDM at the DESI ceiling: the 95th percentile of Δχ² is ≤ 0.067 (TT) and ≤ 0.054 (lensing) at σ_X = 2.633×10⁻², and no realization reaches Δχ² > 4. Its only imprint is a scatter in the inferred H0 of 0.48 km/s/Mpc (θ = 1), about 12 times smaller than the Planck–SH0ES difference, so it neither tightens σ_X nor explains the H0 tension. A deterministic harmonic modulation X = A sin(2π ln a/P + φ) with A ≤ 0.025 and P ∈ [0.1, 2] is allowed by DESI DR2 BAO in 58 of 60 cases (Δχ²_BAO < 4) and is indistinguishable from ΛCDM in the CMB in all 60. (Appendix J.)

**Fuente:** R1, R2, H0.

### D3.11 (A4) — `papers/stochastic-dark-energy-desi-dr2.md` l.98 (insertar tras)

**Línea tras la que se inserta:**

> where $x = \ln(1/(1+z))$. However, for the redshift range probed by DESI ($z \lesssim 2.3$) and Euclid ($z \lesssim 2.0$), $g(z)$ is essentially unconstrained by data. I set $z_* = 1.5$ (illustrative), with the understanding that this parameter is degenerate with $\sigma$ and cannot be independently calibrated. For all numerical results in Sections 4–6, $g(z) \approx 1$ in the observed range, meaning the effective amplitude $\sigma_X$ absorbs any early-time suppression.

**Propuesto:**

> Note: at the DESI ceiling the CMB does not require this activation factor: with g ≡ 1 the OU residual is indistinguishable from ΛCDM in Planck 2018 TT and lensing (p95 Δχ² ≤ 0.067; Appendix J). g(z) is kept only as a modelling option.

**Fuente:** R1; AUDITORIA §6 ("Matiza" A4).

### D3.12 (§7.2) — `papers/stochastic-dark-energy-desi-dr2.md` l.378 (insertar tras)

**Línea tras la que se inserta:**

> I emphasize that this constraint applies specifically to the additive OU/QNM kernel tested here, assuming a fixed CPL background. The results do not rule out stochastic models that are degenerate with the background evolution or that operate on scales not probed by BAO. Future analyses that jointly fit $\{w_0, w_a, \sigma_X\}$ will be needed to resolve this degeneracy.

**Propuesto:**

> The CMB does not tighten this bound: at σ_X = 2.633×10⁻² the residual is indistinguishable from ΛCDM in Planck 2018 TT and lensing, and its only imprint is an H0 scatter of 0.48 km/s/Mpc (Appendix J).

**Fuente:** R1.

### D3.13 (§8 tabla F) — `papers/stochastic-dark-energy-desi-dr2.md` l.440 (insertar tras)

**Línea tras la que se inserta:**

> | **F5** (Degeneracy) | $\sigma_X = 0$ when $w_0, w_a$ are free | **Robust when background is free.** Joint MLE on the same 7 BAO bins with free $\{w_0, w_a, \theta, \sigma_X\}$ yields best-fit near $w_0 \approx -1.001$, $w_a \approx +0.014$ (essentially $\Lambda$CDM) and $\sigma_X \approx 5 \times 10^{-10}$ ($\to 0$), $\Delta\text{AIC} = +4.00$ vs CPL (`results/joint_w0wa_sigma/`, projected covariance). The null result is therefore not an artifact of fixing the background. However, when the background is fixed to the external CPL values preferred by CMB+SNe, a larger amplitude appears. See Appendix I for details. |

**Propuesto:**

> | **F6** (CMB) | No CMB imprint beyond Planck 2018 errors at the DESI ceiling | **Consistent.** p95(Δχ²) ≤ 0.067 (TT ℓ<30) and ≤ 0.054 (lensing) at σ_X = 2.633×10⁻²; H0 scatter 0.48 km/s/Mpc (Appendix J). |

**Fuente:** R1.

### D3.14 (§8 conclusión) — `papers/stochastic-dark-energy-desi-dr2.md` l.442 (matizar (añadir al final; además de D1.3))

**Texto actual:**

> **Conclusion.** The stochastic model is not favored by the current DESI DR2 BAO data when the background is allowed to vary, and I place a phenomenological upper limit on its amplitude. A definitive falsification would require a simultaneous fit of the stochastic parameters with the CPL background and, ideally, the $>20$ bins of Euclid DR1. This work provides the necessary pipeline and a reference limit for that future analysis.

**Propuesto:**

> Append: "The Planck 2018 CMB does not add a constraint at this amplitude (F6, Appendix J)."

**Fuente:** R1.

### D3.15 (§9) — `papers/stochastic-dark-energy-desi-dr2.md` l.458 (insertar tras)

**Línea tras la que se inserta:**

> Both outcomes are scientifically valuable.

**Propuesto:**

> The CMB gives no extra leverage on σ_X at the DESI ceiling (Appendix J); the decisive test remains BAO.

**Fuente:** R1, malla.

### D3.16 (§10) — `papers/stochastic-dark-energy-desi-dr2.md` l.478 (insertar tras)

**Línea tras la que se inserta:**

> More broadly, $\sigma_X < 2.6\times 10^{-2}$ implies that any microscopic mechanism responsible for the observed value of $\Lambda$ must either strongly protect the coupling of the scalar field to the Standard Model (e.g., via shift symmetries), or operate at scales that leave an absolutely minimal cosmological imprint.

**Propuesto:**

> A CMB consistency check with Planck 2018 (TT at ℓ < 30 and lensing bandpowers, via CAMB with the OU history mapped to w(a), H0 solved to keep θ_* fixed) does not distinguish the OU residual from ΛCDM at the DESI ceiling: the 95th percentile of Δχ² is ≤ 0.067 (TT) and ≤ 0.054 (lensing) at σ_X = 2.633×10⁻², and no realization reaches Δχ² > 4. Its only imprint is a scatter in the inferred H0 of 0.48 km/s/Mpc (θ = 1), about 12 times smaller than the Planck–SH0ES difference, so it neither tightens σ_X nor explains the H0 tension. A deterministic harmonic modulation X = A sin(2π ln a/P + φ) with A ≤ 0.025 and P ∈ [0.1, 2] is allowed by DESI DR2 BAO in 58 of 60 cases (Δχ²_BAO < 4) and is indistinguishable from ΛCDM in the CMB in all 60.

**Fuente:** R1, R2, H0.

### D3.17 (Apéndice J) — `papers/stochastic-dark-energy-desi-dr2.md` l.927 (insertar tras (antes de `---` y de References))

**Línea tras la que se inserta:**

> This should be treated as a preliminary, linearized-fiducial check ($S(z)$ held fixed rather than recomputed at each $w_0, w_a$; DESI Gaussian BAO covariance projected to α) pending a full re-analysis with a proper MCMC sampler (e.g., emcee) to confirm the $w_0 \approx -0.99$ minimum is not a local artifact.

**Propuesto:**

> ## Appendix J — CMB Consistency (Planck 2018)
>
> **Construction.** X is a stationary Gaussian process with the BAO kernel covariance σ_X² exp(−θ|Δ ln(1+z)|), smoothed on Δ = 0.05 in ln a. Reading (a), δ ln ρ_DE = X/Ω_Λ0 (with [X-std]), gives ρ_DE(a) = ρ_Λ0[1 + (X_s(a) − X_s(1))/Ω_Λ0] and w = −1 − (1/3) d ln ρ_DE/d ln a, passed to CAMB through `DarkEnergyPPF.set_w_a_table`. ω_b, ω_c and Σm_ν are fixed to Planck 2018, the universe is flat, and H0 is solved for each curve so that 100θ_* = 1.041198 (the ΛCDM value). Non-linear lensing uses halofit (Takahashi). Data: Planck 2018 TT (2 ≤ ℓ < 30) and the lensing bandpowers (L = 8–400). 200 realizations, common random numbers across θ, Δ and σ_X.
>
> **Results at σ_X = 2.633×10⁻² (Δ = 0.05).** θ = 1: Δχ²_TT p50/p95 = −2.5×10⁻³ / 0.067, Δχ²_φφ = −2.8×10⁻³ / 0.054; θ = 0.001: −1.5×10⁻⁴ / 2.7×10⁻³ and −1.2×10⁻⁴ / 2.3×10⁻³. No realization has Δχ² > 4. H0 p16/p50/p84 = 66.92/67.36/67.88 (θ = 1) and 67.34/67.36/67.38 (θ = 0.001); 3.5 % and 0 % of the realizations fall outside Planck's ±2σ.
>
> **Grid (2.506×10⁻² → 10⁻⁴, ten values, θ ∈ {0.001, 1}).** p95(Δχ²) < 4 in all twenty cases; p95 ∝ σ_X^n with n = 0.998–1.007, and σ(H0) ∝ σ_X^m with m = 0.999–1.001. The 95th percentile is dominated by the cross term with the Planck residuals, 2Σ r_b δ_b/σ_b²; the model's own quadratic term is 6.3×10⁻⁴ (TT) and 7.4×10⁻³ (lensing) at θ = 1, σ_X = 2.506×10⁻².
>
> **Harmonic residual.** For X = A sin(2π ln a/P + φ) (A ∈ {0.005, 0.01, 0.025}, P ∈ {0.1, 0.2, 0.5, 1, 2}, φ ∈ {0, π/2, π, 3π/2}), 58 of 60 curves are allowed by DESI DR2 BAO (Δχ²_BAO < 4 with the 13 distances and full covariance; distances and r_d from CAMB), and none is distinguishable in the CMB; among the allowed, the largest Δχ² is 0.070 (TT) and 0.63 (lensing).
>
> **Limitations.** Cosmological parameters other than H0 are fixed (no full MCMC; not run because Δχ² ≤ 0.067 in the worst case). Δ is a parameter of this analysis, not of the theory. Reading (b), δ ln ρ_DE = X, would multiply the amplitude by 1/Ω_Λ0 ≈ 1.46. The treatment of a tabulated w(a) by halofit in PPF has not been verified.
>
> Code and outputs: `analisis_cmb/` (`realizaciones.py`, `estadisticos.py`, `RESULTADO.md`) and `analisis_cmb/oscilador/`.

**Fuente:** R1, malla, R2; RESULTADO.md §1, §3, §6; `analisis_cmb/resultados/takahashi/sigma_2.506e-02/descomposicion_dchi2.json` (término cruzado).

### D3.18 (Everpresent en el PREPRINT) — `manuscript/PREPRINT.md` l.198 (insertar tras (párrafo nuevo en §5.1; alternativa: §7, item 5 — ver D2.5))

**Línea tras la que se inserta:**

> DESI constrains the **effective** residual amplitude after any amplification—not $\sigma_0$, and not a new value of $\rho_\Lambda$. Without a physical map $\sigma_0\to A_0$, “Euclid will see vacuum noise” is not a theorem of this analysis.

**Propuesto:**

> **Everpresent Λ.** The order-one alternative to the area-counting seed has been tested directly. In its two published implementations (Zwane et al. 2018, Model 2; Das et al. 2023, Model 1), computing r_d with each realization's own history, requiring small dark energy over z ∈ [1000, 5000] and ω_m within ±2σ of Planck, none of the realizations fits DESI DR2 BAO as well as ΛCDM, nor within Δχ² < 4: 0/20000 for Model 2 (< 1.5×10⁻⁴ at 95 %) and 0/300000 and 0/200000 for Model 1 with α = 0.01 and 0.005 (< 1.0×10⁻⁵, < 1.5×10⁻⁵). Without the CMB conditions, rare realizations pass, all with ω_m far from Planck or order-one dark energy at recombination. As small fluctuations around ΛCDM, DESI DR2 bounds the Zwane et al. amplitude to α < 1.85×10⁻² (μ fixed). These results hold for those phenomenological implementations, with the CMB only through a proxy (|Ω_de| < 0.25), not for the idea of Everpresent Λ in general.

**Fuente:** `analisis_cmb/everpresent/EVERPRESENT_vs_DESI.md` §3e y Conclusiones 2–3, 6 (`resultados/rd_real.json`): criterio final (c) — r_d calculado con la historia de cada realización, |Ω_de| < 0.25 en z ∈ [1000, 5000] y ω_m a ±1.6 % de Planck —: Modelo 2 de Zwane et al. 0/20000 (< 1.5×10⁻⁴), Modelo 1 de Das et al. 0/300000 (α = 0.01) y 0/200000 (α = 0.005), 0/500000 en total (< 1.0×10⁻⁵ y < 1.5×10⁻⁵). Limitaciones §4: implementaciones fenomenológicas, CMB por proxy. EVERPRESENT_vs_DESI.md §2 (`profile_desi_everpresent.json`): como fluctuación pequeña alrededor de ΛCDM, α_Zwane < 1.85×10⁻² (μ fijo), < 2.41×10⁻² (μ libre), σ_X < 1.29×10⁻² con la correlación del proceso de Ahmed et al. (θ* = 2.72).

### D3.19 (§7 limitaciones) — `manuscript/PREPRINT.md` l.264 (matizar (añadir al final) + item nuevo)

**Texto actual:**

> 3. **Background degeneracy:** with only seven BAO bins, flexible $w(z)$ can absorb residual variance (Section 3.3). Multi-probe fits are a separate analysis.

**Propuesto:**

> Append to item 3: "The CMB check of §3.5 is a consistency test with parameters fixed except H0, not a multi-probe fit."
> New item 7: "7. **CMB check:** parameters other than H0 fixed to Planck 2018 (no MCMC); smoothing scale Δ introduced for the check; reading (a), δ ln ρ_DE = X/Ω_Λ0 (reading (b) would scale the amplitude by 1.46); halofit with tabulated w(a) in PPF not verified."

**Fuente:** RESULTADO.md §6.

### D3.20 (§8.3) — `manuscript/PREPRINT.md` l.312 (insertar tras (dentro del bloque de comandos))

**Línea tras la que se inserta:**

> python scripts/slip_bridge.py

**Propuesto:**

> # CMB consistency (requires CAMB; see analisis_cmb/README.md)
> cd analisis_cmb && ./venv/bin/python realizaciones.py --halofit takahashi && ./venv/bin/python estadisticos.py resultados/takahashi/sigma_2.633e-02

**Fuente:** `analisis_cmb/README.md` (comandos).

### D3.21 (§9) — `manuscript/PREPRINT.md` l.320 (insertar tras)

**Línea tras la que se inserta:**

> 2. **Model kill:** coherent tachyonic growth with rank-1 covariance is excluded ($\Delta\chi^2\approx +25$).

**Propuesto:**

> 3. **CMB consistency:** at the DESI ceiling the OU residual is invisible to Planck 2018 (p95 Δχ² ≤ 0.067) and does not explain the H0 tension (σ(H0) ≈ 0.48 km/s/Mpc); a percent-level harmonic residual is equally allowed by DESI (58/60) and invisible in the CMB. (Renumber the following items.)

**Fuente:** R1, R2, H0.

### D3.22 (§9, futuro) — `manuscript/PREPRINT.md` l.325 (matizar (añadir tras la primera frase))

**Texto actual:**

> Future data (Euclid BAO multi-bin products; full DESI covariances) can tighten $\sigma_X$ with the same residual likelihood. Future theory work must pay the price of Routes 1 or 3 (or an inflation-style seed redefinition)—not invent linear gain on white noise.

**Propuesto:**

> Insert after the first sentence: "The CMB gives no additional leverage at this amplitude (Section 3.5); BAO remains the route."

**Fuente:** R1, malla.

### D3.23 (HONEST_HEADLINES) — `papers/HONEST_HEADLINES.md` l.19 (insertar tras (además de D2.28))

**Línea tras la que se inserta:**

> | Residual $\sigma_X$ sold as baryon density, or as an EHT photon-ring result | Different quantities; see `scope-and-mixups.md` |

**Propuesto:**

> | “The residual explains the Hubble tension” | Its CMB imprint is σ(H0) ≈ 0.48 km/s/Mpc, ≈12× smaller than 5.68 km/s/Mpc. |
> | “Planck tightens σ_X” | It does not: p95(Δχ²) ≤ 0.067 at the DESI ceiling. |

**Fuente:** R1, H0.

### D3.24 (HONEST_HEADLINES) — `papers/HONEST_HEADLINES.md` l.29 (insertar tras)

**Línea tras la que se inserta:**

> 5. Process: tests, gap-label rule, no undeclared free parameters.

**Propuesto:**

> 6. CMB consistency (Planck 2018): the OU residual at the DESI ceiling and a percent-level harmonic residual are both invisible in TT ℓ<30 and lensing; no H0-tension explanation.

**Fuente:** R1, R2, H0.

### D3.25 (README) — `README.md` l.23 (insertar tras)

**Línea tras la que se inserta:**

> 6. [`papers/scope-and-mixups.md`]\(papers/scope-and-mixups.md) — residual \(\sigma_X\) is not \(\omega_b\), not EHT  

**Propuesto:**

> 7. [`analisis_cmb/RESULTADO.md`]\(analisis_cmb/RESULTADO.md) — CMB consistency check (Planck 2018)

**Fuente:** Navegación.

### D3.26 (README) — `README.md` l.35 (insertar tras)

**Línea tras la que se inserta:**

> Best fit still goes to \(\sigma_X\to 0\). A coherent tachyonic growing mode is excluded when it is active (rank-1 covariance). Soft amplification of a pure Sorkin seed \(\sim 10^{-61}\) to Euclid-scale residuals needs \(\sim 10^{56}\) (to the DESI ceiling, \(\sim 10^{59}\)) and is not free — see the amplification note.

**Propuesto:**

> Planck 2018 (TT ℓ<30 + lensing) does not distinguish this residual from ΛCDM at the DESI ceiling (p95 Δχ² ≤ 0.067) and adds no bound; its H0 imprint, 0.48 km/s/Mpc, does not explain the H0 tension. A percent-level harmonic residual is equally allowed (58/60) and CMB-invisible. See `analisis_cmb/RESULTADO.md`.

**Fuente:** R1, R2, H0. (El Layout del README ya incluye `analisis_cmb/` desde la reorganización.)

### D3.27 (resume.txt) — `papers/resume.txt` l.17 (insertar tras)

**Línea tras la que se inserta:**

>   External lab optics: separate optics project (NOT this paper)

**Propuesto:**

>   CMB check: Planck 2018 TT l<30 + lensing — OU at sigma_X=2.633e-2 invisible (p95 dchi2 <= 0.067); sigma(H0)=0.48 (theta=1); harmonic X: 58/60 BAO-allowed, 0/60 CMB-distinguishable

**Fuente:** R1, R2.

### D3.28 (resume.txt) — `papers/resume.txt` l.57 (insertar tras (sección nueva 1b))

**Línea tras la que se inserta:**

>   Nested CPL + (σ, θ) extensions are not preferred (ΔAIC ≈ +4).

**Propuesto:**

>
> 1b. CMB CONSISTENCY (Planck 2018, analisis_cmb/)
> - 200 OU realizations at sigma_X = 2.633e-2, H0 solved for theta_*: p95 dchi2 TT = 0.067, phiphi = 0.054 (theta = 1); 0 with dchi2 > 4.
> - Grid 2.506e-2 -> 1e-4: p95 < 4 everywhere, linear in sigma_X. Planck does not bound sigma_X below DESI.
> - sigma(H0) = 0.48 km/s/Mpc (theta = 1), 0.020 (theta = 0.001): 12x and 280x below 73.04 - 67.36 = 5.68.
> - Harmonic X = A sin(2pi ln a / P + phi), A <= 0.025: 58/60 allowed by DESI BAO, 0/60 distinguishable in the CMB.

**Fuente:** R1, malla, H0, R2.

### D3.29 (resume.txt) — `papers/resume.txt` l.85 (insertar tras)

**Línea tras la que se inserta:**

> - The result supports the working hypothesis of **Late-Time Vacuum Homogeneity** at currently observable scales.

**Propuesto:**

> - The Planck 2018 CMB adds no constraint at this amplitude (section 1b).

**Fuente:** R1.

### D3.30 (resume.txt) — `papers/resume.txt` l.89 (insertar tras)

**Línea tras la que se inserta:**

> - The current analysis pipeline is ready to be applied to upcoming datasets.

**Propuesto:**

> - The CMB gives no extra leverage on sigma_X at the DESI ceiling; BAO remains the test.

**Fuente:** R1, malla.

### D3.31 (unimodular §2) — `theory/unimodular-gravity-vacuum-smoothness.md` l.82 (insertar tras)

**Línea tras la que se inserta:**

> This is a direct empirical statement: the dark-energy density is homogeneous at the level of the residual bound (profile 95% CL $\sigma_X<2.6\times 10^{-2}$; null preferred) across the observable universe.

**Propuesto:**

> The Planck 2018 CMB does not add a bound at this level: at σ_X = 2.633×10⁻² the residual is indistinguishable from ΛCDM in TT (ℓ < 30) and lensing (`analisis_cmb/RESULTADO.md`).

**Fuente:** R1.

### D3.32 (unimodular §8, ISW) — `theory/unimodular-gravity-vacuum-smoothness.md` l.283 (insertar tras (nota bajo la tabla))

**Línea tras la que se inserta:**

> | ISW effect | No late-time ISW anomaly from dark-energy clustering | Detection of a residual ISW signal from dark-energy perturbations |

**Propuesto:**

> *Note on the ISW row: at the DESI ceiling, Planck 2018 TT at ℓ < 30 and lensing do not distinguish the residual from ΛCDM (p95 Δχ² ≤ 0.067), so with current precision this row cannot falsify anything.*

**Fuente:** R1; AUDITORIA §6 ("Matiza la fila ISW").

### D3.33 (euclid-protocol §8) — `theory/euclid-protocol-vacuum-relaxation.md` l.186 (insertar tras (nota bajo la tabla))

**Línea tras la que se inserta:**

> | ISW / CMB lensing | DE perturbations if present |

**Propuesto:**

> *With Planck 2018 this channel has no leverage at σ_X ≲ 2.6×10⁻²: p95 Δχ² ≤ 0.067 (TT ℓ<30) and ≤ 0.054 (lensing) (`analisis_cmb/RESULTADO.md`).*

**Fuente:** R1.

### D3.34 (euclid-protocol §12) — `theory/euclid-protocol-vacuum-relaxation.md` l.247 (insertar tras)

**Línea tras la que se inserta:**

> Working residual ceiling under the OU kernel (this pipeline): $\sigma_X < 2.6\times 10^{-2}$ (95% CL, free-$\theta$ profile, 7 bins). Profile and joint runs on the real $\alpha$ file are consistent with a null residual (amplitude driven to the floor when the background is free).

**Propuesto:**

> Consistency check: at that ceiling the residual is invisible to Planck 2018 (TT ℓ<30 + lensing) and its H0 imprint is 0.48 km/s/Mpc (`analisis_cmb/RESULTADO.md`).

**Fuente:** R1.

### D3.35 (amplification-gap §2) — `papers/amplification-gap.md` l.56 (insertar tras)

**Línea tras la que se inserta:**

> **What DESI actually measured** is a **null residual bound**, not a new $\rho_\Lambda$. Mean vacuum density $\rho_\Lambda\sim\Omega_{\Lambda 0}\rho_{\mathrm{crit}}$ comes from standard cosmology (e.g. Planck-class $H_0,\Omega_\Lambda$); this corpus constrains the **residual amplitude** $\sigma_X$ (in units of $\rho_{{\rm crit},0}$) on top of that background.

**Propuesto:**

> The CMB does not tighten this amplitude below the DESI ceiling (Planck 2018 TT ℓ<30 + lensing, `analisis_cmb/RESULTADO.md`).

**Fuente:** R1, malla.

### D3.36 (principle) — `theory/principle-of-vacuum-smoothness.md` l.36 (insertar tras)

**Línea tras la que se inserta:**

> That robustness across different statistical frameworks suggests that the smoothness of the dark energy sector is not an accident of one model assumption. It is a genuine feature that future theories should accommodate from the start.

**Propuesto:**

> Two further checks (`analisis_cmb/`): the Planck 2018 CMB does not distinguish the residual from ΛCDM at the DESI ceiling, and a deterministic harmonic residual of percent amplitude is equally allowed by DESI BAO (58/60) and CMB-invisible. "Perfectly homogeneous" (above) therefore means homogeneous to σ_X ≲ 2.6×10⁻² in X, for stochastic and harmonic residuals alike.

**Fuente:** R1, R2; AUDITORIA §6 ("perfectly homogeneous vale solo hasta …").

## D4. Slip: "marginal, not starved"

Números del bloque C: `analisis_cmb/bloque_C/logs/slip_bridge_despues.log` (bloque C): ε = 1, δ_m = 1, |γ−1| = 4.953e-2 / 2.090e-2 / 1.070e-2 a z = 0.5 / 1.0 / 1.5; 0.17/|γ−1| = 3.4 / 8.1 / 15.9; 0.05/|γ−1| = 1.01 / 2.39 / 4.67. Las 17 frases de COMPARACION_B2 ("Afirmaciones que cambian de sentido", sin la fila l.924, que no es de slip), reescritas. En las tres principales (PREPRINT resumen, cita de §6.2 y conclusión 4) va la cautela: *"The σ_X → slip map is an order-of-magnitude estimate and the best case requires ε = 1, δ_m = 1 and z ≈ 0.5."*

### D4.1 — PRINCIPAL (resumen) — `manuscript/PREPRINT.md` l.28 (reemplazar)

**Texto actual:**

> Gravitational slip $\gamma=\Phi/\Psi$ is the correct operator for anisotropic leakage past volume-preserving (SDiff) structure, but inherits the same amplitude starvation: even full anisotropy of the BAO-bounded residual gives at most $\lvert\gamma-1\rvert\sim 5\times 10^{-2}$ ($z=0.5$, $\delta_m=1$): $\sim 3$–$16\times$ below the current measurement and $1.0$–$4.7\times$ below the forecast floor, which it only reaches in that best case.

**Propuesto:**

> Gravitational slip $\gamma=\Phi/\Psi$ is the correct operator for anisotropic leakage past volume-preserving (SDiff) structure, and with the BAO-bounded residual it is marginal, not starved: full anisotropy gives at most $\lvert\gamma-1\rvert\approx 5\times 10^{-2}$ at $z=0.5$, about $3.4\times$ below the current measurement (Maus et al.) and at the level of the Sakr et al. forecast floor ($\approx 1.0\times$), falling to $2.4\times$ and $4.7\times$ below that floor at $z=1$ and $1.5$. The σ_X → slip map is an order-of-magnitude estimate and the best case requires ε = 1, δ_m = 1 and z ≈ 0.5.

**Fuente:** `analisis_cmb/bloque_C/logs/slip_bridge_despues.log` (bloque C): ε = 1, δ_m = 1, |γ−1| = 4.953e-2 / 2.090e-2 / 1.070e-2 a z = 0.5 / 1.0 / 1.5; 0.17/|γ−1| = 3.4 / 8.1 / 15.9; 0.05/|γ−1| = 1.01 / 2.39 / 4.67.

### D4.2 — `manuscript/PREPRINT.md` l.232 (reemplazar)

**Texto actual:**

> ### 6.2 Option 0: anisotropic gap is real in principle, starved in amplitude

**Propuesto:**

> ### 6.2 Option 0: anisotropic gap is real in principle, marginal in amplitude

**Fuente:** `analisis_cmb/bloque_C/logs/slip_bridge_despues.log` (bloque C): ε = 1, δ_m = 1, |γ−1| = 4.953e-2 / 2.090e-2 / 1.070e-2 a z = 0.5 / 1.0 / 1.5; 0.17/|γ−1| = 3.4 / 8.1 / 15.9; 0.05/|γ−1| = 1.01 / 2.39 / 4.67. COMPARACION_B2 (PREPRINT:232).

### D4.3 — PRINCIPAL (cita de §6.2) — `manuscript/PREPRINT.md` l.254 (reemplazar)

**Texto actual:**

> > Slip is the correct *kind* of operator for anisotropic leakage past SDiff, but with the BAO-bounded residual amplitude it is not an observational shortcut. It inherits the amplification problem.

**Propuesto:**

> > Slip is the correct *kind* of operator for anisotropic leakage past SDiff. With the BAO-bounded residual amplitude it is marginal, not starved: at ε = 1 it reaches $\lvert\gamma-1\rvert\approx 5\times 10^{-2}$ at $z=0.5$, level with the forecast floor ($\approx 1.0\times$) and $\approx 3.4\times$ below the current measurement. It is not a demonstrated shortcut: ε < 1, δ_m > 1 or z ≳ 1 push it below reach. The σ_X → slip map is an order-of-magnitude estimate and the best case requires ε = 1, δ_m = 1 and z ≈ 0.5.

**Fuente:** `analisis_cmb/bloque_C/logs/slip_bridge_despues.log` (bloque C): ε = 1, δ_m = 1, |γ−1| = 4.953e-2 / 2.090e-2 / 1.070e-2 a z = 0.5 / 1.0 / 1.5; 0.17/|γ−1| = 3.4 / 8.1 / 15.9; 0.05/|γ−1| = 1.01 / 2.39 / 4.67.

### D4.4 — PRINCIPAL (conclusión 4) — `manuscript/PREPRINT.md` l.322 (reemplazar)

**Texto actual:**

> 4. **Slip honesty:** anisotropic leakage past SDiff is a real structural crack in principle and an amplitude-starved channel in practice.

**Propuesto:**

> 4. **Slip honesty:** anisotropic leakage past SDiff is a real structural crack in principle and a marginal channel in practice — at best level with forecast floors ($\approx 1.0\times$ the Sakr et al. floor, $\approx 3.4\times$ below Maus et al., at $z=0.5$). The σ_X → slip map is an order-of-magnitude estimate and the best case requires ε = 1, δ_m = 1 and z ≈ 0.5.

**Fuente:** `analisis_cmb/bloque_C/logs/slip_bridge_despues.log` (bloque C): ε = 1, δ_m = 1, |γ−1| = 4.953e-2 / 2.090e-2 / 1.070e-2 a z = 0.5 / 1.0 / 1.5; 0.17/|γ−1| = 3.4 / 8.1 / 15.9; 0.05/|γ−1| = 1.01 / 2.39 / 4.67.

### D4.5 — `manuscript/README.md` l.34 (reemplazar)

**Texto actual:**

> - Gravitational-slip amplitude starvation (Option 0)

**Propuesto:**

> - Gravitational slip: marginal amplitude (Option 0)

**Fuente:** `analisis_cmb/bloque_C/logs/slip_bridge_despues.log` (bloque C): ε = 1, δ_m = 1, |γ−1| = 4.953e-2 / 2.090e-2 / 1.070e-2 a z = 0.5 / 1.0 / 1.5; 0.17/|γ−1| = 3.4 / 8.1 / 15.9; 0.05/|γ−1| = 1.01 / 2.39 / 4.67.

### D4.6 — `manuscript/CLAIMS.md` l.22 (reemplazar (solo la columna de la afirmación))

**Texto actual:**

> | **C7** | Even $\varepsilon=1$ anisotropic fraction of the BAO-bounded residual gives $\|\gamma-1\|\lesssim 5\times 10^{-2}$ ($z=0.5$, $\delta_m=1$), at best level with the forecast floor $\sim 0.05$ | `scripts/slip_bridge.py` |

**Propuesto:**

> With ε = 1, the BAO-bounded residual gives $\|\gamma-1\|\lesssim 5\times 10^{-2}$ ($z=0.5$, $\delta_m=1$): ≈3.4× below Maus et al. and ≈1.0× the Sakr et al. constant-η forecast floor (order-of-magnitude map)

**Fuente:** `analisis_cmb/bloque_C/logs/slip_bridge_despues.log` (bloque C): ε = 1, δ_m = 1, |γ−1| = 4.953e-2 / 2.090e-2 / 1.070e-2 a z = 0.5 / 1.0 / 1.5; 0.17/|γ−1| = 3.4 / 8.1 / 15.9; 0.05/|γ−1| = 1.01 / 2.39 / 4.67.

### D4.7 — `papers/resume.txt` l.16 (reemplazar)

**Texto actual:**

>   Act IV bridge: scripts/slip_bridge.py  (|gamma-1|~5e-2 best case; ~1-5x below forecast floor)

**Propuesto:**

>   Act IV bridge: scripts/slip_bridge.py  (|gamma-1|~5e-2 best case at z=0.5: ~1.0x forecast floor, ~3.4x below Maus; marginal, not starved)

**Fuente:** `analisis_cmb/bloque_C/logs/slip_bridge_despues.log` (bloque C): ε = 1, δ_m = 1, |γ−1| = 4.953e-2 / 2.090e-2 / 1.070e-2 a z = 0.5 / 1.0 / 1.5; 0.17/|γ−1| = 3.4 / 8.1 / 15.9; 0.05/|γ−1| = 1.01 / 2.39 / 4.67.

### D4.8 — `papers/anisotropic-slip-option0.md` l.25 (reemplazar)

**Texto actual:**

> **Conclusion (preview).** The anisotropic “gap” is a real structural crack in SDiff protection, and $\eta$ is the right kind of observable — but with $\sigma_X\lesssim 2.6\times 10^{-2}$ it inherits the same amplitude problem as the isotropic residual. It is an additional channel that turns on with amplification, not free amplification.

**Propuesto:**

> **Conclusion (preview).** The anisotropic “gap” is a real structural crack in SDiff protection, and $\eta$ is the right kind of observable. With $\sigma_X\lesssim 2.6\times 10^{-2}$ the best case (ε = 1, δ_m = 1) reaches $\lvert\gamma-1\rvert\approx 5\times 10^{-2}$ at $z=0.5$, level with Euclid-like forecasts: marginal, not starved. Its visibility depends on ε, δ_m and z rather than on amplification.

**Fuente:** `analisis_cmb/bloque_C/logs/slip_bridge_despues.log` (bloque C): ε = 1, δ_m = 1, |γ−1| = 4.953e-2 / 2.090e-2 / 1.070e-2 a z = 0.5 / 1.0 / 1.5; 0.17/|γ−1| = 3.4 / 8.1 / 15.9; 0.05/|γ−1| = 1.01 / 2.39 / 4.67.

### D4.9 — `papers/anisotropic-slip-option0.md` l.41 (reemplazar)

**Texto actual:**

> **Implication.** The “grieta” (anisotropic leakage past SDiff) is quantified as amplitude-starved via `slip_bridge.py`, not tested as a new MCMC. Option 0 remains literature + scaling.

**Propuesto:**

> **Implication.** The “grieta” (anisotropic leakage past SDiff) is quantified as marginal in amplitude via `slip_bridge.py` (best case ≈1.0× the forecast floor at z = 0.5), not tested as a new MCMC. Option 0 remains literature + scaling.

**Fuente:** `analisis_cmb/bloque_C/logs/slip_bridge_despues.log` (bloque C): ε = 1, δ_m = 1, |γ−1| = 4.953e-2 / 2.090e-2 / 1.070e-2 a z = 0.5 / 1.0 / 1.5; 0.17/|γ−1| = 3.4 / 8.1 / 15.9; 0.05/|γ−1| = 1.01 / 2.39 / 4.67.

### D4.10 — `papers/anisotropic-slip-option0.md` l.118 (reemplazar (solo el inciso))

**Texto actual:**

> With $\sigma(\gamma)\sim 0.11$, current sensitivity to $\lvert \eta-1 \rvert$ (or $\lvert \gamma-1 \rvert$) is still $\mathcal{O}(0.1)$, far above the $\sim 5\times 10^{-2}$ scale of the $f=1$, $\sigma_X=2.6\times 10^{-2}$ estimate below.

**Propuesto:**

> "far above the $\sim 5\times 10^{-2}$ scale" → "a factor $\sim 2$ above the $\sim 5\times 10^{-2}$ scale".

**Fuente:** `analisis_cmb/bloque_C/logs/slip_bridge_despues.log` (bloque C): ε = 1, δ_m = 1, |γ−1| = 4.953e-2 / 2.090e-2 / 1.070e-2 a z = 0.5 / 1.0 / 1.5; 0.17/|γ−1| = 3.4 / 8.1 / 15.9; 0.05/|γ−1| = 1.01 / 2.39 / 4.67. O(0.1)/4.95e-2 ≈ 2.

### D4.11 — `papers/anisotropic-slip-option0.md` l.157 (reemplazar)

**Texto actual:**

> | Current DESI-era slip (Maus $\gamma=1.17\pm 0.11$) | $\mathcal{O}(0.1)$ | **above** prediction ($\sim 2\times$) |

**Propuesto:**

> | Current DESI-era slip (Maus $\gamma=1.17\pm 0.11$) | $\mathcal{O}(0.1)$ | **∼2× above** prediction (0.17 is ≈3.4× the best case) |

**Fuente:** `analisis_cmb/bloque_C/logs/slip_bridge_despues.log` (bloque C): ε = 1, δ_m = 1, |γ−1| = 4.953e-2 / 2.090e-2 / 1.070e-2 a z = 0.5 / 1.0 / 1.5; 0.17/|γ−1| = 3.4 / 8.1 / 15.9; 0.05/|γ−1| = 1.01 / 2.39 / 4.67. Ya cambiado en C a "above (∼2×)"; aquí solo se añade el 3.4×.

### D4.12 — `papers/anisotropic-slip-option0.md` l.158 (sin cambio (ya dice "comparable: ∼0.6–1.0×" desde el bloque C))

**Línea tras la que se inserta:**

> | Euclid-like forecasts (Sakr et al.) | $\mathcal{O}(0.03\text{–}0.05)$ | **comparable**: $\sim 0.6$–$1.0\times$ |

**Propuesto:**

> No change.

**Fuente:** `analisis_cmb/bloque_C/logs/slip_bridge_despues.log` (bloque C): ε = 1, δ_m = 1, |γ−1| = 4.953e-2 / 2.090e-2 / 1.070e-2 a z = 0.5 / 1.0 / 1.5; 0.17/|γ−1| = 3.4 / 8.1 / 15.9; 0.05/|γ−1| = 1.01 / 2.39 / 4.67.

### D4.13 — `papers/anisotropic-slip-option0.md` l.160 (reemplazar)

**Texto actual:**

> **Interpretation.** Even the most optimistic anisotropic fraction does not make the gap visible at current DESI precision, nor automatically at Euclid, unless there is amplification of the bare residual amplitude (or a shear sector decoupled from the BAO-bounded isotropic $\sigma_X$).

**Propuesto:**

> **Interpretation.** The most optimistic anisotropic fraction does not make the gap visible at current DESI precision (≈3.4× short at z = 0.5); at Euclid-like forecast precision it is marginal (≈1.0× the Sakr et al. floor) for ε ≈ 1, δ_m ≈ 1, z ≈ 0.5, and out of reach for ε ≲ 0.5 or z ≳ 1. The map is order-of-magnitude.

**Fuente:** `analisis_cmb/bloque_C/logs/slip_bridge_despues.log` (bloque C): ε = 1, δ_m = 1, |γ−1| = 4.953e-2 / 2.090e-2 / 1.070e-2 a z = 0.5 / 1.0 / 1.5; 0.17/|γ−1| = 3.4 / 8.1 / 15.9; 0.05/|γ−1| = 1.01 / 2.39 / 4.67. Con ε = 0.5, |γ−1| = 2.5e-2, 2× bajo el suelo.

### D4.14 — `papers/anisotropic-slip-option0.md` l.162 (reemplazar (solo la última frase))

**Texto actual:**

> That word — amplification — is the same one that appears in amplification (`amplification-gap.md`) when mapping Sorkin-scale seeds ($\sim 10^{-61}$) up to BAO-relevant $\sigma_X$ (factors of order $A_0/\sigma_0\sim 10^{56}$). The anisotropic channel does not evade that problem; it inherits it.

**Propuesto:**

> "The anisotropic channel does not evade that problem; it inherits it." → "The anisotropic channel does not evade the seed problem (σ_0 → σ_X); given σ_X at the DESI bound, it is marginal rather than starved."

**Fuente:** `analisis_cmb/bloque_C/logs/slip_bridge_despues.log` (bloque C): ε = 1, δ_m = 1, |γ−1| = 4.953e-2 / 2.090e-2 / 1.070e-2 a z = 0.5 / 1.0 / 1.5; 0.17/|γ−1| = 3.4 / 8.1 / 15.9; 0.05/|γ−1| = 1.01 / 2.39 / 4.67.

### D4.15 — `papers/anisotropic-slip-option0.md` l.183 (reemplazar)

**Texto actual:**

> - That Euclid will “easily” see the gap without amplification.

**Propuesto:**

> - That Euclid will see the gap: the best case (ε = 1, δ_m = 1, z ≈ 0.5) is only level with the forecast floor, and anything less is below it.

**Fuente:** `analisis_cmb/bloque_C/logs/slip_bridge_despues.log` (bloque C): ε = 1, δ_m = 1, |γ−1| = 4.953e-2 / 2.090e-2 / 1.070e-2 a z = 0.5 / 1.0 / 1.5; 0.17/|γ−1| = 3.4 / 8.1 / 15.9; 0.05/|γ−1| = 1.01 / 2.39 / 4.67.

### D4.16 — `exploratory/notes/data-pack-option0-internet.md` l.97 (reemplazar)

**Texto actual:**

> **Conclusion (unchanged):** the external literature confirms Option 0. The correct operator exists and is measured or forecast; the amplitude still requires amplification before a detection story makes sense.

**Propuesto:**

> **Conclusion (revised):** the external literature confirms Option 0. The correct operator exists and is measured or forecast; with σ_X at the DESI bound, the best-case slip (≈5×10⁻², `scripts/slip_bridge.py`) is comparable to the Sakr et al. constant-η forecast floor — marginal, not starved, and order-of-magnitude only.

**Fuente:** `analisis_cmb/bloque_C/logs/slip_bridge_despues.log` (bloque C): ε = 1, δ_m = 1, |γ−1| = 4.953e-2 / 2.090e-2 / 1.070e-2 a z = 0.5 / 1.0 / 1.5; 0.17/|γ−1| = 3.4 / 8.1 / 15.9; 0.05/|γ−1| = 1.01 / 2.39 / 4.67. (El |η−1| ∼ f·(Ω_DE/Ω_m)·σ_X de l.87 de este fichero no se tocó en C; ver COMPARACION_C, "Filas que no se cambian".)

### D4.17 — `papers/amplification-gap.md` l.83 (reemplazar)

**Texto actual:**

> Every closed channel that looked like a shortcut — anisotropic slip, casual GRB phase stories, pupil/tesseract optics — **inherits this same gap**. That is why they stay amplitude-starved even when the **operator** is right (Option 0).

**Propuesto:**

> Every closed channel that looked like a shortcut — anisotropic slip, casual GRB phase stories, pupil/tesseract optics — inherits the seed gap (σ_0 → σ_X). Slip with σ_X already at the DESI bound is a separate question: there it is marginal, not starved (see `anisotropic-slip-option0.md`).

**Fuente:** `analisis_cmb/bloque_C/logs/slip_bridge_despues.log` (bloque C): ε = 1, δ_m = 1, |γ−1| = 4.953e-2 / 2.090e-2 / 1.070e-2 a z = 0.5 / 1.0 / 1.5; 0.17/|γ−1| = 3.4 / 8.1 / 15.9; 0.05/|γ−1| = 1.01 / 2.39 / 4.67.

## D5. PREPRINT: referencia de α = 1 y figuras de §3

### D5.1 — `manuscript/PREPRINT.md` l.78 (reemplazar)

**Texto actual:**

> **Important separation:** DESI analyses that combine BAO with CMB+SN often prefer dynamical $w_0$–$w_a$. My BAO-only residual analysis uses a CPL background near $\Lambda$CDM ($ w_0\approx -1.001$, $w_a\approx +0.014$; `results/eos_cpl_desi_dr2/`, projected covariance, figure6 α) as the smooth reference for $\alpha=1$ under the fiducial cosmology of the published $\alpha$ products. The stochastic question is about *residuals*, not a re-fit of the full multi-probe DE equation of state.

**Propuesto:**

> **Important separation:** DESI analyses that combine BAO with CMB+SN often prefer dynamical $w_0$–$w_a$. The reference for $\alpha=1$ in my BAO-only residual analysis is the fiducial $\Lambda$CDM cosmology of the published $\alpha$ products ($w_0=-1$, $w_a=0$). As a check, a CPL fit to the same data lands at $w_0\approx -1.001$, $w_a\approx +0.014$ (`results/eos_cpl_desi_dr2/`, projected covariance, figure6 α), indistinguishable from that reference. The stochastic question is about *residuals*, not a re-fit of the full multi-probe DE equation of state.

**Fuente:** `results/eos_cpl_desi_dr2/eos_cpl_summary.json` (w0 = −1.00095, wa = +0.01370 tras el bloque A); `results/profile_sigma_x/profile_sigma_x.json` (perfil con α = 1). La instrucción: la referencia es el fiducial ΛCDM; el CPL es una comprobación.

### D5.2 — `manuscript/PREPRINT.md` l.151 (matizar (añadir al final))

**Texto actual:**

> From the residual analysis (profile / likelihood ratio under the OU kernel as implemented in the public pipeline): $\sigma_X < 2.6\times 10^{-2}$ (95% CL).

**Propuesto:**

> Append: "The profile is shown in Fig. 1 ([`figures/profile_sigma_x.png`]\(../figures/profile_sigma_x.png))."

**Fuente:** `figures/profile_sigma_x.png` = copia de `results/profile_sigma_x/profile_sigma_x.png` que escribe `scripts/profile_sigma_x_desi.py` (decisión del autor en la reorganización).

### D5.3 — `manuscript/PREPRINT.md` l.158 (matizar (añadir al final))

**Texto actual:**

> When $\{w_0,w_a\}$ are freed together with residual parameters, the stochastic extension is not preferred ($\Delta\mathrm{AIC}\approx +4$). Marginalising over a flexible smooth background drives $\sigma_X$ toward zero: extra residual variance is absorbed by background freedom when the number of BAO bins is small ($N_{\mathrm{bin}}=7$).

**Propuesto:**

> Append: "Fig. 2 ([`figures/joint_w0wa_sigma.png`]\(../figures/joint_w0wa_sigma.png)) shows the joint fit of $\{w_0, w_a, \sigma_X\}$."

**Fuente:** `figures/joint_w0wa_sigma.png` = copia de `results/joint_w0wa_sigma/joint_w0wa_sigma.png` que escribe `scripts/joint_w0wa_sigma_desi.py`.


---

114 puntos. Ningún paper editado.
