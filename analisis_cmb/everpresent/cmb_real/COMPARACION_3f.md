# 3f: cierre con la tasa de expansión en la recombinación, antes / después

Cambios de texto de este commit. Números de `cmb_real/theta_star.json` → `expansion_at_recombination`; los resultados anteriores de `theta_star.json` no cambian (comprobado: el json sin la clave nueva es idéntico).

Notas:
- N_eff = 2.99 ± 0.17 no está en la Tabla 2 de arXiv:1807.06209 (ΛCDM base), sino en el resumen, la ec. `equ:neff1` y la Tabla 4 (`tab:grid_1paramext`, 95 %), con BAO. Sin BAO, en la columna de θ*: 2.89^{+0.36}_{−0.38} (95 %).
- La sensibilidad que sale es ΔH/H ≈ 0.55 %, por debajo del 1 %, no "few-percent". En los textos en inglés va como "≈ 0.6 %".
- X = 11 %: es el factor en z = 1100 de la menos desviada del criterio (c), 12909, con Ω_de(1100) = −0.26 (1 − 1/√1.26 = 10.7 %). En z ∈ [800, 1500], ≥ 14 %. Sobre todos los sorteos a < 10σ de θ*, ≥ 6.2 %.
- Conclusión 2 de EVERPRESENT_vs_DESI.md: solo se añade "(≥ 21σ)" y que no necesita la condición de recombinación.

## (a) EVERPRESENT_vs_DESI.md, conclusión de 3f

Fichero: `analisis_cmb/everpresent/EVERPRESENT_vs_DESI.md`

**Antes:**

```text
**Conclusión de 3f:**
- θ* sustituye al proxy en el Modelo 1: el resultado es el mismo, 0.
- En el Modelo 2 no lo sustituye. θ* es un solo número, y el tramo z ∈ (3, 1000), que los datos de
  BAO no fijan, tiene libertad para ajustarlo.
- El 0 del criterio (c) para el Modelo 2 descansa en la condición de energía oscura pequeña alrededor
  de la recombinación, que θ* por sí solo no comprueba.
- Saber si esas historias, con Ω_de ≈ ±0.8 en z ≈ 1100, son compatibles con los espectros de Planck
  exige la verosimilitud completa con ρ_DE con signo: opción (a) de `cmb_real/EXPLORACION.md`.
```

**Después:**

```text
**Expansión en la recombinación** (`cmb_real/theta_star.json` → `expansion_at_recombination`;
resumen en `theta_star.txt`). Se compara H(z) con H_ΛCDM(z) en z ∈ [800, 1500] (ΛCDM con el mismo
ω_m, ω_r y el H0 de Planck) en cada sorteo que deja θ* a < 10σ: 174 sorteos de las 25 supervivientes.
- En todos, max |H/H_ΛCDM − 1| es al menos del **6.2 %**. El mínimo es el de la superviviente 6078,
  con Ω_de(1100) = +0.11; la mediana es 27 %.
- Las 5 del criterio (c) sin ventana:

| índice | ω_m | Ω_de(1100) | 1/√(1 − Ω_de(1100)) | H/H_ΛCDM en z ∈ [800, 1500] | mín. sobre sorteos de max \|H/H_ΛCDM − 1\| |
|---|---|---|---|---|---|
| 1540 | 0.1407 | −0.96 | 0.714 | 0.712–0.720 | 28.8 % |
| 8163 | 0.1407 | +0.85 | 2.616 | 1.988–3.072 | 207 % |
| 11072 | 0.1453 | −0.87 | 0.732 | 0.731–0.739 | 26.9 % |
| 12909 | 0.1453 | −0.26 | 0.893 | 0.856–0.909 | 14.4 % |
| 16616 | 0.1410 | −0.77 | 0.751 | 0.747–0.778 | 25.3 % |

**Sensibilidad de Planck a la tasa de expansión en la recombinación** (analogía de fondo, generosa;
no es una cota formal):
- Referencia: N_eff = 2.99 ± 0.17 (68 %), 2.99^{+0.34}_{−0.33} (95 %), con TT,TE,EE+lowE+lensing+BAO
  (arXiv:1807.06209, resumen, ec. `equ:neff1` y Tabla 4 `tab:grid_1paramext`; comprobado en el TeX).
  No está en la Tabla 2, que es ΛCDM base. Sin BAO, en la columna de θ*: 2.89^{+0.36}_{−0.38} (95 %).
- ΔN_eff = 0.34 (2σ) da Δρ_rad/ρ_rad = 0.2271 · 0.34/(1 + 0.2271 · 3.046) = 4.56 %.
- La radiación es el 24.2 % del total en z = 1090.
- Así que **ΔH/H = ½ · 0.242 · 0.0456 = 0.55 %**.
- N_eff también cambia el amortiguamiento de Silk y la fase de los picos, no solo H; por eso la
  analogía es generosa.

Las historias que ajustan θ* tienen H en la recombinación desviado ≥ 6.2 % (≥ 11 veces esa
sensibilidad). Las del criterio (c), ≥ 14 % en z ∈ [800, 1500] y ≥ 10.7 % en z = 1100 (≥ 19 veces).

**Conclusión de 3f:**
- **Modelo 1:** θ* sustituye al proxy y da el mismo 0. El más cercano está a 21σ.
- **Modelo 2:** θ* solo no lo sustituye. θ* es un solo número, y el tramo z ∈ (3, 1000), que BAO no
  fija, tiene libertad para ajustarlo. La fracción esperada que pasa (c) con θ* es 8.5 × 10⁻⁶ a 5σ y
  ≤ 2 × 10⁻⁵ a 10σ.
- Todas esas historias ajustan θ* con una tasa de expansión en la recombinación que se aleja de ΛCDM
  al menos un orden de magnitud más que la sensibilidad de Planck por la analogía con N_eff.
- El 0 del Modelo 2 se mantiene con la condición de recombinación (energía oscura pequeña alrededor
  de z ≈ 1100), y el argumento de la tasa de expansión la respalda.
- **La verosimilitud completa con ρ_DE con signo no se ha corrido.** Es la opción (a) de
  `cmb_real/EXPLORACION.md`, pendiente; la comparación con N_eff no la sustituye.
```

## (b) EVERPRESENT_vs_DESI.md, conclusión 2

Fichero: `analisis_cmb/everpresent/EVERPRESENT_vs_DESI.md`

**Antes:**

```text
   - Con θ* calculado con la historia de cada run (3f), ninguno de los 6 runs que pasan A o B en
     algún criterio se acerca al valor de Planck: el más cercano está a 21σ. Para el Modelo 1, θ* da
     el mismo 0 que el proxy de ventana.
```

**Después:**

```text
   - Con θ* calculado con la historia de cada run (3f), ninguno de los 6 runs que pasan A o B en
     algún criterio se acerca al valor de Planck: el más cercano está a 21σ (≥ 21σ). Para el
     Modelo 1, θ* da el mismo 0 que el proxy de ventana, sin necesidad de la condición de
     recombinación.
```

## (b) EVERPRESENT_vs_DESI.md, conclusión 3

Fichero: `analisis_cmb/everpresent/EVERPRESENT_vs_DESI.md`

**Antes:**

```text
     - Las cinco tienen |Ω_de(1100)| = 0.26–0.96.
     - Con θ* y la ventana a la vez, 0.
     - El 0 del Modelo 2 depende, por tanto, de exigir energía oscura pequeña en la recombinación,
       que θ* solo no comprueba.
```

**Después:**

```text
     - Las cinco tienen |Ω_de(1100)| = 0.26–0.96, es decir, H en z = 1100 desviado de ΛCDM entre un
       10.7 % y un 162 % (factor 1/√(1 − Ω_de) entre 0.71 y 2.62).
     - La sensibilidad de Planck a la tasa de expansión en la recombinación es ≈ 0.55 % (ΔN_eff =
       0.34 al 95 %; analogía de fondo, generosa, no una cota formal).
     - Con θ* y la ventana a la vez, 0.
     - El 0 del Modelo 2 se mantiene con la condición de recombinación, respaldada por ese argumento.
       θ* solo admite una fracción esperada ≤ 2 × 10⁻⁵.
     - La verosimilitud completa con ρ_DE con signo no se ha corrido (opción (a) de
       `cmb_real/EXPLORACION.md`, pendiente).
```

## (b) EVERPRESENT_vs_DESI.md, conclusión 4

Fichero: `analisis_cmb/everpresent/EVERPRESENT_vs_DESI.md`

**Antes:**

```text
   - En el Modelo 2 pasan 5/20000 con algún sorteo del tramo no muestreado (fracción esperada
     8.5 × 10⁻⁶): el ángulo acústico se puede ajustar con ese tramo, y la energía oscura en la
     recombinación sigue siendo de orden uno.
   - Que esas historias sean o no compatibles con los espectros de Planck requiere la verosimilitud
     completa.
```

**Después:**

```text
   - En el Modelo 2 pasan 5/20000 con algún sorteo del tramo no muestreado (fracción esperada
     8.5 × 10⁻⁶): el ángulo acústico se puede ajustar con ese tramo, y la energía oscura en la
     recombinación sigue siendo de orden uno.
   - Esas historias tienen H en la recombinación desviado ≥ 6.2 % de ΛCDM (≥ 10.7 % en z = 1100 en
     las del criterio (c)), frente a una sensibilidad de Planck de ≈ 0.55 % por la analogía con N_eff.
     Lo que excluye en el Modelo 2 sigue siendo BAO junto con la condición de recombinación.
   - La verosimilitud completa con ρ_DE con signo, que lo comprobaría formalmente, no se ha corrido.
```

## (c) PREPRINT §4.2

Fichero: `manuscript/PREPRINT.md`

**Antes:**

```text
Replacing the proxy by the acoustic angle θ* computed from each realization's own history (the full Planck likelihood is not available because CAMB and CLASS cannot follow a dark-energy density that changes sign, as 70 of the 71 BAO survivors do), Model 1 still gives 0, while for Model 2 the stretch z ∈ (3, 1000), unconstrained by BAO and sampled here (100 draws per realization), can tune θ*: 5/20000 realizations pass criterion (c) with θ* within 5σ of Planck for some draws (expected fraction 8.5×10⁻⁶), all with |Ω_de(z ≈ 1100)| = 0.26–0.96, so the Model 2 exclusion rests on requiring small dark energy at recombination, which only the full CMB likelihood can test.
```

**Después:**

```text
Replacing the proxy by the acoustic angle θ* computed from each realization's own history (the full Planck likelihood is not available because CAMB and CLASS cannot follow a dark-energy density that changes sign, as 70 of the 71 BAO survivors do), Model 1 still gives 0 (closest at 21σ). For Model 2 the stretch z ∈ (3, 1000), unconstrained by BAO and sampled here (100 draws per realization), can tune θ*: 5/20000 realizations pass criterion (c) with θ* within 5σ of Planck for some draws (expected fraction 8.5×10⁻⁶; ≤ 2×10⁻⁵ within 10σ), all with |Ω_de(z ≈ 1100)| ≥ 0.26, i.e. H at recombination off from ΛCDM by ≥ 11 % (≥ 6 % in every draw that brings any survivor within 10σ of θ*), whereas Planck's sensitivity to the expansion rate there is ≈ 0.6 % (ΔN_eff = 0.34 at 95 %; a generous background analogy, not a formal bound). The Model 2 exclusion therefore holds with the recombination condition; the full CMB likelihood with a sign-changing dark-energy density was not run.
```

## (c) PREPRINT abstract (l. 24)

Fichero: `manuscript/PREPRINT.md`

**Antes:**

```text
none fits as well as ΛCDM (0/20000 and 0/500000 realizations).
```

**Después:**

```text
none fits as well as ΛCDM (0/20000 and 0/500000 realizations). With the acoustic angle θ* computed for each realization in place of the recombination condition, Model 1 stays at 0 (≥ 21σ), while Model 2 admits an expected fraction ≤ 2×10⁻⁵, all with the expansion rate at recombination off by ≥ 11 % against a ≈ 0.6 % Planck sensitivity (N_eff analogy); the full CMB likelihood with a sign-changing dark-energy density was not run.
```

## (c) PREPRINT conclusión 3 (l. 371)

Fichero: `manuscript/PREPRINT.md`

**Antes:**

```text
3. **Everpresent Λ:** its two published implementations do not fit DESI DR2 BAO together with the CMB conditions (0/20000 and 0/500000 realizations; §4.2).
```

**Después:**

```text
3. **Everpresent Λ:** its two published implementations do not fit DESI DR2 BAO together with the CMB conditions (0/20000 and 0/500000 realizations; §4.2). Model 1: 0 also with θ* computed for each realization (≥ 21σ). Model 2: 0 with the recombination condition; θ* alone admits an expected fraction ≤ 2×10⁻⁵, all with H at recombination off by ≥ 11 % against a ≈ 0.6 % Planck sensitivity (N_eff analogy); the full likelihood with signed ρ_DE was not run.
```

## (c) CLAIMS C10

Fichero: `manuscript/CLAIMS.md`

**Antes:**

```text
. With θ* computed from each history in place of the proxy, Model 1 stays at 0 but 5/20000 Model 2 realizations pass for some draws of the unsampled stretch z ∈ (3, 1000) (expected fraction 8.5×10⁻⁶, all with order-one Ω_de at z ≈ 1100), so the Model 2 part rests on the small-dark-energy-at-recombination condition |
```

**Después:**

```text
. Model 1: 0 with θ* (≥ 21σ). Model 2: 0 with the recombination condition; θ* alone admits an expected fraction ≤ 2×10⁻⁵, all with $\lvert\Omega_{\rm de}(1100)\rvert \ge 0.26$, i.e. H at recombination off by ≥ 11 % against a ≈ 0.6 % Planck sensitivity (N_eff analogy, not a formal bound); the full likelihood with signed ρ_DE was not run |
```

## (d) PREPRINT limitación 8

Fichero: `manuscript/PREPRINT.md`

**Antes:**

```text
the stretch $z\in(3,1000)$ of Model 2 is sampled (N = 100 draws per realization).
```

**Después:**

```text
the stretch $z\in(3,1000)$ of Model 2 is sampled (N = 100 draws per realization). The comparison of H at recombination with Planck's N_eff sensitivity is a background analogy, not a formal bound; the full likelihood with signed ρ_DE (a modified CAMB) was not run.
```

## (d) analisis_cmb/README.md

Fichero: `analisis_cmb/README.md`

**Antes:**

```text
| `everpresent/` | Everpresent Λ (Ahmed–Dodelson–Greene–Sorkin and follow-ups) against DESI DR2 BAO; `EVERPRESENT_vs_DESI.md` |
```

**Después:**

```text
| `everpresent/` | Everpresent Λ (Ahmed–Dodelson–Greene–Sorkin and follow-ups) against DESI DR2 BAO; `EVERPRESENT_vs_DESI.md`. The CMB enters through a recombination proxy and through θ* computed for each realization, with the expansion rate at recombination (§3f, `everpresent/cmb_real/theta_star.py`); the full Planck likelihood with a sign-changing ρ_DE was not run (`everpresent/cmb_real/EXPLORACION.md`) |
```
