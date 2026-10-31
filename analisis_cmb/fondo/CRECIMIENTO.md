# El ritmo de crecimiento fσ₈ a la cota de DESI (exploración)

**Predicción, escrita antes de correr `crecimiento.py`:**
- El cambio de fσ₈ es del 1–2 % a z = 0.5.
- Queda entre 3 y 5 veces por debajo del error por tramo de DESI DR1 full-shape.

Exploración: no cambia `results/` ni `figures/`. Script `crecimiento.py`; salidas `crecimiento.{json,txt}`.

## Montaje

- **Fondo:** ΛCDM plano de Planck 2018 (Ω_m = 0.3153, Ω_Λ = 0.6847, el Ω_Λ de la tabla S_z), sin radiación.
- **Temblor:** un cambio constante δΩ_Λ = ±2.633 × 10⁻², con universo plano y H0 fijo, así que Ω_m → Ω_m ∓ δ.
  Es el modo θ → 0, el que fija la cota de C2 (`FONDO.md`).
- **f(z) = Ω_m(z)^0.55**, como se pidió. La f exacta (d ln D/d ln a de la ecuación) coincide a ≤ 0.07 % en el
  fiducial, y el cambio relativo de fσ₈ con una u otra difiere en ≤ 0.04 puntos.
- **σ₈(z) = σ₈ D(z):** D sale de resolver la ecuación de crecimiento lineal en ln a,
  D'' + (2 + d ln H/d ln a) D' − (3/2) Ω_m(a) D = 0, con D = a en a = 10⁻³.
  - **Normalización de base:** D → a en la era de materia. El universo temprano y la amplitud primordial no
    cambian, porque el temblor es un cambio tardío de energía oscura.
  - **Variante:** σ₈(z = 0) fijo. Es una cota inferior del efecto, y solo queda el cambio de f y de la forma
    de D(z).

**Errores de DESI DR1 full-shape.** arXiv:2411.12022 (DESI 2024 VII) da solo parámetros cosmológicos; no tiene
una tabla de fσ₈ por tramo. Los valores por tramo están en el artículo compañero, DESI 2024 V, *Full-Shape Galaxy
Clustering from Galaxies and Quasars* (arXiv:2411.12021), **sección 7.1, tabla 9**: ajuste ShapeFit de base,
valores MAP, columna "ShapeFit only", en forma de fσ_s8/(fσ_s8)_fid. El error relativo se toma como la media del
error superior y el inferior.

| tramo | z_eff | fσ_s8/(fσ_s8)_fid | error relativo usado |
|---|---|---|---|
| BGS | 0.295 | 0.80 ± 0.20 | 20 % |
| LRG1 | 0.510 | 1.09 +0.12 −0.14 | 13 % |
| LRG2 | 0.706 | 1.05 ± 0.12 | 12 % |
| LRG3 | 0.919 | 0.96 +0.11 −0.10 | 10.5 % |
| ELG2 | 1.317 | 0.95 +0.11 −0.08 | 9.5 % |
| QSO | 1.491 | 1.16 ± 0.12 | 12 % |

La precisión combinada es del 4.7 % (sección 8 del mismo artículo).

## Resultados

**Cambio relativo de fσ₈ a la cota (f = Ω_m(z)^0.55, D normalizado en la era de materia):**

| z | +δΩ_Λ | −δΩ_Λ | por f | por D | con f exacta | con σ₈(0) fijo | tramo DESI | error | error / cambio |
|---|---|---|---|---|---|---|---|---|---|
| 0.3 | −4.77 % | +4.47 % | −3.46 % | −1.35 % | −4.81 % | −2.82 % | BGS | 20 % | 4.2 |
| **0.5** | **−3.76 %** | +3.45 % | −2.75 % | −1.03 % | −3.78 % | −1.79 % | LRG1 | 13 % | **3.5** |
| 0.7 | −2.95 % | +2.66 % | −2.17 % | −0.79 % | −2.95 % | −0.97 % | LRG2 | 12 % | 4.1 |
| 1.0 | −2.06 % | +1.82 % | −1.53 % | −0.54 % | −2.06 % | −0.06 % | LRG3 | 10.5 % | 5.1 |

"error / cambio" usa el mayor de los dos signos. Las columnas "por f" y "por D" son para +δΩ_Λ.

## Por qué

- **Domina f.** Con H0 fijo, quitar 0.026 a Ω_m es un −8.4 % en Ω_m. Eso baja Ω_m(z) un (1 − Ω_m(z)) × 8.4 %,
  y f baja 0.55 veces eso: −2.8 % a z = 0.5.
- **D se suma en el mismo sentido.** Más Λ frena antes el crecimiento, así que D(z) baja un 1 % a z = 0.5
  (con la amplitud temprana fija).
- **El efecto decrece con z,** porque Ω_m(z) → 1 y el temblor en Ω_Λ pesa menos.
- **Con σ₈(0) fijo, el efecto se reduce a la mitad** (1.8 % a z = 0.5) y desaparece a z ≈ 1. Ese caso no
  corresponde a un temblor tardío con un universo temprano fijo, pero da el rango: 1.8–3.8 % a z = 0.5.
- **Frente a DESI DR1 el efecto queda 3.5–5 veces por debajo del error de cada tramo.** Un error por tramo del
  2–3 % (lo que se espera de DESI Y5 y Euclid; esa cifra no se comprueba aquí) lo dejaría en 1.3–1.9 σ por
  tramo a z = 0.5: marginal.
- **Cuidado: es una huella del mismo modo.** El modo θ → 0 es la dirección de Ω_m (`FONDO.md`), así que en fσ₈
  también se confunde con un cambio de Ω_m. Solo distingue algo si Ω_m está fijado por el CMB, que es justo
  como está montado este cálculo.

## Predicción

- **Cambio del 1–2 % a z = 0.5: no se cumplió con la normalización de base.** Da 3.8 % (+δΩ_Λ) y 3.5 % (−δΩ_Λ),
  el doble de lo predicho. Solo con σ₈(0) fijo cae en el rango (1.8 %).
- **Entre 3 y 5 veces por debajo del error por tramo: cumplida.** Es 3.5 veces en LRG1 (z = 0.5) y 3.5–5.1 en
  los cuatro z pedidos.

## Previsiones por tramo (DESI Y5, Euclid)

**DESI 2016** (DESI Collaboration, arXiv:1611.00036), sección 2.4.3, "Summary of Forecasts": error relativo
previsto de fσ₈ por tramo, 14 000 deg², con k_max = 0.1 y 0.2 h Mpc⁻¹. Los centros están separados 0.1; el texto
no da el ancho, así que se lee como Δz = 0.1. La última columna es el mayor cambio de fσ₈ a la cota (de los dos
signos de δΩ_Λ), calculado con `crecimiento.py` en el centro de cada tramo.

| tabla | z | error, k_max = 0.1 | error, k_max = 0.2 | cambio a la cota | cambio / error (0.1; 0.2) |
|---|---|---|---|---|---|
| 2.5 (BGS) | 0.35 | 5.83 % | 2.60 % | 4.50 % | 0.77; 1.73 |
| 2.5 (BGS) | 0.45 | 6.35 % | 3.77 % | 3.99 % | 0.63; 1.06 |
| 2.3 | 0.65 | 3.31 % | 1.57 % | 3.13 % | 0.95; 2.00 |
| 2.3 | 0.75 | 2.10 % | 1.01 % | 2.78 % | 1.32; 2.75 |
| 2.3 | 0.85 | 2.12 % | 1.01 % | 2.46 % | 1.16; 2.44 |
| 2.3 | 0.95 | 2.09 % | 0.99 % | 2.18 % | 1.05; 2.21 |
| 2.3 | 1.05 | 2.23 % | 1.11 % | 1.94 % | 0.87; 1.75 |
| 2.3 | 1.15 | 2.25 % | 1.14 % | 1.73 % | 0.77; 1.52 |
| 2.3 | 1.25 | 2.25 % | 1.16 % | 1.55 % | 0.69; 1.33 |
| 2.3 | 1.35 | 2.90 % | 1.73 % | 1.39 % | 0.48; 0.80 |
| 2.3 | 1.45 | 3.06 % | 1.87 % | 1.25 % | 0.41; 0.67 |
| 2.3 | 1.55 | 3.53 % | 2.27 % | 1.12 % | 0.32; 0.49 |

La precisión agregada prevista es ≈ 0.74 % (k_max = 0.1) y ≈ 0.38 % (k_max = 0.2), en la sección 2.4.1.

**Euclid preparation VII** (arXiv:1910.09273): las fσ₈(z_i) de los cuatro tramos espectroscópicos (z = 1.0, 1.2,
1.4 y 1.65) son parámetros intermedios del Fisher. Se dejan libres en cada tramo y se proyectan sobre los
parámetros cosmológicos (sección 3.2), y el artículo no da sus errores. Por eso Euclid se cita sin número.

**Lectura:** en 0.75 ≤ z ≤ 1.25 el error previsto para DESI es del 1.0–2.3 % por tramo, y el cambio a la cota es
0.7–2.8 σ por tramo. Es marginal por tramo y detectable solo combinando tramos con k_max = 0.2.

## Antes / después de la frase (PREPRINT, trabajo futuro)

**Antes:**
> The growth rate is the velocity-side signature: at the DESI limit the residual shifts $f\sigma_8$ by ≈ 3.8 % at $z = 0.5$, 3.5 times below the DESI DR1 full-shape errors (arXiv:2411.12021, Table 9; cosmology in arXiv:2411.12022); it becomes marginal at the 2–3 % precision of DESI Y5 and Euclid (`analisis_cmb/fondo/CRECIMIENTO.md`).

**Después:**
> The growth rate is the velocity-side signature: at the DESI limit the residual shifts $f\sigma_8$ by ≈ 3.8 % at $z = 0.5$, 3.5 times below the DESI DR1 full-shape errors (arXiv:2411.12021, Table 9; cosmology in arXiv:2411.12022), with $\Omega_m$ fixed by the CMB (the constant mode is otherwise degenerate with $\Omega_m$, as for the distances); it becomes marginal at the per-bin precision forecast for DESI Y5, 1.0–2.3 % per $\Delta z = 0.1$ bin at $0.75 \le z \le 1.25$ ($k_{\max} = 0.2$–$0.1\,h$ Mpc$^{-1}$; arXiv:1611.00036, Table 2.3), where the shift is 0.7–2.8 σ per bin, while Euclid preparation VII (arXiv:1910.09273) projects the per-bin $f\sigma_8$ onto cosmological parameters without quoting them (`analisis_cmb/fondo/CRECIMIENTO.md`).

(En las citas, el enlace al propio fichero se deja como texto, porque desde aquí apuntaría mal.)

Cambios:
1. "at the 2–3 % precision of DESI Y5 and Euclid" pasa a ser el error por tramo de DESI 2016 (tabla 2.3) y la
   nota de que Euclid preparation VII no da errores por tramo.
2. Se añade la condición de Ω_m fijado por el CMB.
