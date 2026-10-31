# Texto: reloj de estacionariedad y válvula tanh, antes / después

Comprobación previa: la tanh está en el TeX de Zwane et al. 2018 (`everpresent/fuentes/zwane2018_1703.06265.tex`, línea 232): `\beq \Omega_{de}(\lambda) = \tanh[\hat\Omega_{de}(\lambda)] \eeq`, justo después del núcleo `\label{m2}` (ec. 2.9). Con `\numberwithin{equation}{section}`, es la **ec. (2.10)**.

Los enlaces del PREPRINT van como enlaces relativos a `manuscript/`, como el resto del documento.

## 1. PREPRINT §3.2

Fichero: `manuscript/PREPRINT.md`

**Antes:**

```text
The profile is shown in Fig. 1 ([`figures/profile_sigma_x.png`](../figures/profile_sigma_x.png)).
```

**Después:**

```text
The profile is shown in Fig. 1 ([`figures/profile_sigma_x.png`](../figures/profile_sigma_x.png)). At the limit the profile sits at the lower edge of the θ search (θ = 10⁻³, fully correlated bins); the θ → 0 supremum is 2.644 × 10⁻² (+0.4 %). Redefining the stationarity clock as cosmic time, conformal time or z instead of ln a moves the limit by less than 0.7 % ([`analisis_cmb/reloj/RELOJ.md`](../analisis_cmb/reloj/RELOJ.md)): it is set by the rank-1, fully correlated limit of the kernel and does not depend on the clock.
```

## 2. PREPRINT limitación 2

Fichero: `manuscript/PREPRINT.md`

**Antes:**

```text
2. **Kernel dependence:** bounds are for the OU/QNM residual kernels stated; other residual models need their own covariance and MLE.
```

**Después:**

```text
2. **Kernel dependence:** bounds are for the OU/QNM residual kernels stated; other residual models need their own covariance and MLE. Within the OU family the bound is insensitive to the stationarity clock (ln a, t, η, z within ±0.7 %) because it is reached in the fully correlated limit; DESI DR2 does not resolve θ and hence cannot select the clock.
```

## 3. PREPRINT §4.2

Fichero: `manuscript/PREPRINT.md`

**Antes:**

```text
The Model 2 exclusion therefore holds with the recombination condition; the full CMB likelihood with a sign-changing dark-energy density was not run.
```

**Después:**

```text
The Model 2 exclusion therefore holds with the recombination condition; the full CMB likelihood with a sign-changing dark-energy density was not run. Model 2 maps Ω_de = tanh Ω̂ (Zwane et al. 2018, eq. 2.10), which keeps H² > 0 by construction; |Ω̂| exceeds 1 somewhere in 87.5 % of the 20000 realizations. Removing the map, which is our modification and not a published model, 58.7 % of the realizations reach H² ≤ 0, two of the five θ*-survivors are discarded and the other three fit BAO worse (χ² from 9–28 to 63–1209), while the count with the recombination condition stays at 0: the exclusion does not rest on the map ([`analisis_cmb/everpresent/cmb_real/VALVULA.md`](../analisis_cmb/everpresent/cmb_real/VALVULA.md)).
```

## 4. EVERPRESENT_vs_DESI.md, sección 3g

Fichero: `analisis_cmb/everpresent/EVERPRESENT_vs_DESI.md`

**Antes:**

```text
### Comparación con [D2]
```

**Después:**

```text
### 3g. Válvula tanh

- [Z] ec. (2.10): Ω_de = tanh Ω̂ mantiene H² > 0 por construcción. |Ω̂| > 1 en algún punto guardado en
  el 87.5 % de las 20000 realizaciones.
- Sin la tanh (Ω_de = Ω̂), que es una modificación nuestra y no un modelo publicado, el 58.7 % llega a
  H² ≤ 0 y se descarta.
- De las 5 del criterio (c) que ajustan θ* (3f), 2 se descartan (8163, 11072). Las otras 3 ajustan
  peor BAO: χ²_c pasa de 9–28 a 63–1209.
- El criterio (c) con ventana sigue en 0 con y sin válvula; la exclusión no depende de la tanh.
- Detalle y comprobaciones: [`cmb_real/VALVULA.md`](cmb_real/VALVULA.md).

### Comparación con [D2]
```

## 4. EVERPRESENT_vs_DESI.md, Limitaciones

Fichero: `analisis_cmb/everpresent/EVERPRESENT_vs_DESI.md`

**Antes:**

```text
- **Modelos fenomenológicos.** Los Modelos 1 y 2
```

**Después:**

```text
- **Test sin válvula (3g).** Ω_de = Ω̂, sin la tanh de [Z] ec. (2.10), es una modificación nuestra para
  medir el papel de la tanh, no un modelo publicado.
- **Modelos fenomenológicos.** Los Modelos 1 y 2
```

## 4. CLAIMS C10, evidencia

Fichero: `manuscript/CLAIMS.md`

**Antes:**

```text
| `analisis_cmb/everpresent/EVERPRESENT_vs_DESI.md` §3e, §3f; PREPRINT §4.2 |
```

**Después:**

```text
| `analisis_cmb/everpresent/EVERPRESENT_vs_DESI.md` §3e, §3f; PREPRINT §4.2; robustness: `analisis_cmb/everpresent/cmb_real/VALVULA.md`, `analisis_cmb/reloj/RELOJ.md` |
```

## 5. analisis_cmb/README.md

Fichero: `analisis_cmb/README.md`

**Antes:**

```text
| `oscilador/` | Deterministic oscillating X(ln a) against DESI DR2 BAO and Planck 2018; `oscilador/README.md` |
```

**Después:**

```text
| `everpresent/cmb_real/valvula.py`, `valvula.json`, `VALVULA.md` | Role of the tanh map of Model 2 (Zwane et al. eq. 2.10): fraction of realizations where it acts, survivors and criterion (c) without it (our modification, not a published model) |
| `reloj/` | σ_X profile with four stationarity clocks (ln a, cosmic time, conformal time, z): `perfil_reloj.py`, `perfil_reloj.{json,txt}`, `RELOJ.md`; before/after of the related text in `COMPARACION_RELOJ_VALVULA.md` |
| `oscilador/` | Deterministic oscillating X(ln a) against DESI DR2 BAO and Planck 2018; `oscilador/README.md` |
```
