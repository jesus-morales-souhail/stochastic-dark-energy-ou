# Texto: supervivencia del Modelo 1 no convergida, antes / después

Números de producción leídos de `resultados/bao_m1_alpha0.01.json` (n_seeds = 300000, n_complete = 82795, frac_complete = 0.275983 → X = 27.6 %) y `bao_m1_alpha0.005.json` (n_seeds = 200000, n_complete = 198858, frac_complete = 0.994290).
En la frase del PREPRINT se añade "for α = 0.01", porque el 32.3 % de [D1] y el ≈ 16 % del continuo son de ese α.

## 1. EVERPRESENT_vs_DESI.md §1

Fichero: `analisis_cmb/everpresent/EVERPRESENT_vs_DESI.md`

**Antes:**

```text
Para el test de la sección 3 se usa n = 10⁵.
```

**Después:**

```text
Para el test de la sección 3 se usa n = 10⁵.

La fracción que completa no está convergida en ninguna fila. Con pasos uniformes en ln a baja de
forma monótona al refinar: 26.3, 19.4, 18.0 y 17.3 % con n = 10⁴, 10⁵, 3×10⁵ y 10⁶ (2000 semillas).
El límite continuo es ≈ 16 %: 16.2 ± 0.7 % por extrapolación √Δ y 15.7 ± 0.8 % con la corrección de
continuidad de Broadie–Glasserman–Kou. La fracción de producción (uniforme en a, n = 10⁵) y el 32.3 %
de [D1] quedan por encima por el primer paso largo con Λ congelada. H0 y Ω_Λ(0) de las supervivientes
no dependen del reloj, y los recuentos de la sección 3 tampoco
([`cmb_real/RELOJ_M1.md`](cmb_real/RELOJ_M1.md)).
```

## 2. EVERPRESENT_vs_DESI.md §3c

Fichero: `analisis_cmb/everpresent/EVERPRESENT_vs_DESI.md`

**Antes:**

```text
Los límites "< x al 95 %" con 0 éxitos son la regla de 3 (3/N).

### 3d.
```

**Después:**

```text
Los límites "< x al 95 %" con 0 éxitos son la regla de 3 (3/N).

Semillas y runs completados (`bao_m1_alpha*.json`):
- α = 0.01: 300000 semillas, 82795 completan (27.6 %).
- α = 0.005: 200000 semillas, 198858 completan (99.4 %).

Los 0 de (i) y sus límites (< 1.0 × 10⁻⁵ y < 1.5 × 10⁻⁵) son **sobre semillas**: las que llegan a
H² ≤ 0 cuentan como realizaciones que no reproducen el universo. La fracción que completa es
artefacto de la discretización (§1).

### 3d.
```

## 3. PREPRINT §4.2

Fichero: `manuscript/PREPRINT.md`

**Antes:**

```text
(< 1.0×10⁻⁵, < 1.5×10⁻⁵).
```

**Después:**

```text
(< 1.0×10⁻⁵, < 1.5×10⁻⁵). Counts are over seeds, with realizations that reach H² ≤ 0 counted as not reproducing the universe; the completed fraction itself is a discretization artefact (for α = 0.01, 27.6 % with Das et al.'s uniform-in-a steps in the production run, 32.3 % in their paper, ≈ 16 % in the continuum limit with steps in ln a, where a √Δ extrapolation and the Broadie–Glasserman–Kou continuity correction agree to 0.5 points), while the survivors' H0 and Ω_Λ(0) and the BAO counts do not depend on the stepping ([`analisis_cmb/everpresent/cmb_real/RELOJ_M1.md`](../analisis_cmb/everpresent/cmb_real/RELOJ_M1.md)).
```

## 3. CLAIMS C10, robustness

Fichero: `manuscript/CLAIMS.md`

**Antes:**

```text
robustness: `analisis_cmb/everpresent/cmb_real/VALVULA.md` |
```

**Después:**

```text
robustness: `analisis_cmb/everpresent/cmb_real/VALVULA.md`, `analisis_cmb/everpresent/cmb_real/RELOJ_M1.md` |
```
