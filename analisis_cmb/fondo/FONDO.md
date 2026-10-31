# La cota de DESI con el fondo libre (exploración)

**Predicción, escrita antes de correr `perfil_fondo.py`:**
- (c) marginalizado: θ → 0 sin cota (degenerado con Ω_m, como en SN).
- (b) perfilado: θ → 0 entre 0.04 y 0.07 (≈ 6.7 × σ(Ω_m)).
- Para θ = 10, la cota queda entre 1.5 y 4 × 10⁻², y cambia menos de un factor 2 entre fondo fijo
  y libre.

Exploración: no cambia `results/`, `figures/` ni ningún paper. Script `perfil_fondo.py`; salidas
`perfil_fondo.{json,txt}`.

## Montaje

**Igual que `scripts/profile_sigma_x_desi.py`**, cuyas funciones se reutilizan sin llamar a su
`main()`:
- los mismos 7 α_DV y la misma covarianza oficial proyectada;
- la misma S y el mismo núcleo;
- la misma malla de σ_X y la misma búsqueda de θ en (10⁻³, 20);
- ΔlnL = −1.92 con `brentq`.

**Hallazgo previo: el "fondo fijo" de C2 es el mejor ajuste de DESI, no Planck.**
- Los α de `figure6` son (D_V/r_d)/(D_V/r_d)^best, referidos al mejor ajuste ΛCDM de DESI DR2
  (`plot_figure6.ipynb` de la publicación de datos: "DESI DR2 best-fit ΛCDM").
- Ajustando los α frente a las D_V/r_d del vector de 13 puntos, ese fondo es Ω_m = 0.29740 y
  h·r_d = 101.538 Mpc. Reproduce los 7 α a ≤ 3 × 10⁻⁴, frente a errores de 0.5–1 %.
- El artículo de DESI DR2 da 0.2975 ± 0.0086 y 101.54 ± 0.73 Mpc con solo DESI (arXiv:2503.14738,
  Sec. 6, ec. `constraints_lcdm_desi`; comprobado en su TeX, `EVERPRESENT_vs_DESI.md` §3d).

**Con ese fondo de referencia:**
- α_th(Ω_m, h·r_d) = [(h·r_d)_f/(h·r_d)] · D̂_V(z; Ω_m)/D̂_V(z; Ω_m,f), ΛCDM plano, y α_th = 1 en
  la referencia.
- **(a) Fondo fijo** (α_th = 1): las funciones del script de DESI.
- **(b) Ω_m y h·r_d perfilados** en cada (σ_X, θ): h·r_d analíticamente, Ω_m por Brent.
- **(c) Ω_m y h·r_d marginalizados,** con prior plano, linealizados alrededor del mejor ajuste con
  σ_X = 0, como en `analisis_cmb/sn/perfil_sn.py`.
- Si un perfil no cruza −1.92 dentro de la malla (≤ 5 × 10⁻²), se extiende hasta σ_X = 1.

## Resultados

**Cota de σ_X al 95 %:**

| tratamiento | θ libre | θ = 0.001 | θ = 0.1 | θ = 1 | θ = 10 |
|---|---|---|---|---|---|
| (a) fondo fijo | **0.02633** (θ = 10⁻³, borde; idéntico bit a bit a C2) | 0.02633 | 0.0210 | 0.0144 | 0.0129 |
| (b) Ω_m, h·r_d perfilados | **0.02633** (θ = 10⁻³, borde) | 0.02633 | 0.0213 | 0.0150 | 0.0133 |
| (c) Ω_m, h·r_d marginalizados | **sin cota** (ΔlnL = −1.32 en σ_X = 1) | sin cota | 0.140 | 0.046 | 0.022 |

En los tres casos, el máximo del perfil está en σ_X → 0.

**Fondo:**

| magnitud | valor |
|---|---|
| coseno entre S y ∂ln α/∂Ω_m (h·r_d fijo, en la referencia) | **0.99996** |
| Ω_m de mejor ajuste con h·r_d libre y σ_X = 0 (7 α_DV) | 0.2978 ± 0.0124 (Fisher), h·r_d = 101.51 Mpc |
| DESI DR2 BAO solo (13 puntos), arXiv:2503.14738 | 0.2975 ± 0.0086, h·r_d = 101.54 ± 0.73 Mpc |
| σ(Ω_m) con h·r_d fijo | 0.0038 |

Con solo los 7 α_DV isótropos, σ(Ω_m) es mayor que la del artículo (0.0124 frente a 0.0086),
porque se pierde la información anisótropa D_M/D_H. El valor central coincide a 0.0003.

**El modo θ → 0:**

| magnitud | valor |
|---|---|
| σ_v = 1/√(SᵀC⁻¹S), error del modo con todo lo demás fijo | 0.00392 |
| amplitud estimada X̂ del modo a partir de los α, en unidades de σ_v | **+0.02** |
| cota de un término ln det puro, √(e^{3.84} − 1) σ_v | **0.02644** |
| supremo de C2 con θ → 0 (`analisis_cmb/reloj/RELOJ.md`) | 0.02644 |

## Por qué

1. **El modo totalmente correlacionado es la dirección de Ω_m** (coseno 0.99996), igual que en SN.
   S se definió como −∂ln D_V/∂Ω_Λ con H0 fijo; en universo plano es ∂ln D_V/∂Ω_m.
2. **La cota de C2 es un factor de Occam, no una medida de estocasticidad.**
   - Los α están referidos al mejor ajuste de DESI, así que el residuo a lo largo del modo es
     prácticamente cero: X̂ = 0.02 σ_v.
   - Entonces ΔlnL(σ) = −½ ln(1 + σ²/σ_v²), y −1.92 se alcanza en √(e^{3.84} − 1) σ_v = 6.75 ×
     0.00392 = 0.02644. Es exactamente el supremo θ → 0 de C2.
   - 2.633e-2 es ese número evaluado en el borde θ = 10⁻³ en vez de en θ = 0.
3. **Perfilar Ω_m y h·r_d no cambia nada** (0.02633 en (a) y en (b)): el mejor ajuste con los 7 α
   coincide con la referencia (0.2978 frente a 0.2974). Perfilar no elimina el término ln det; solo
   lo hace marginalizar.
4. **Marginalizando, el modo θ → 0 no tiene cota**: es degenerado con Ω_m.
5. **Con θ ≳ 1**, el residuo se descorrelaciona entre bins y deja de parecerse a un cambio de Ω_m.
   Ahí DESI lo acota con información real en los tres tratamientos:
   - fijo: 0.0144 (θ = 1) y 0.0129 (θ = 10);
   - perfilado: 0.0150 y 0.0133;
   - marginalizado: 0.046 y 0.022.
   El fondo libre afloja la cota de θ = 10 un factor 1.7, y la de θ = 1 un factor 3.2.

## Predicción

- **(c) marginalizado, θ → 0 sin cota: cumplida.** ΔlnL = −1.32 incluso en σ_X = 1.
- **(b) perfilado, θ → 0 entre 0.04 y 0.07: no se cumplió.** Da 0.02633, igual que el fondo fijo.
  - El σ relevante no es σ(Ω_m) con h·r_d libre (0.0124, que daría 6.7 × 0.0124 = 0.083).
  - Es el error del modo con todo lo demás fijo, σ_v = 0.0039, igual a σ(Ω_m) con h·r_d fijo
    (0.0038). Perfilar conserva el término ln det de esa dirección.
- **θ = 10 entre 1.5 y 4 × 10⁻² y cambio < factor 2 entre fijo y libre: a medias.**
  - El cambio sí es menor que un factor 2: 0.0129 → 0.0133 perfilado, y 0.022 marginalizado
    (factor 1.7).
  - El rango no: el fondo fijo (0.0129) y el perfilado (0.0133) quedan algo por debajo de
    1.5 × 10⁻². Solo el marginalizado (0.022) cae dentro.

## Qué parte de C2 es qué

- **La cota de C2 (σ_X < 2.6 × 10⁻²) está fijada por el modo θ → 0.** Ese modo es un cambio de Ω_m
  de fondo, y su número es el factor de Occam 6.75 σ_v frente al **mejor ajuste ΛCDM de DESI**. Es
  consistencia con el fondo de referencia, y ese fondo es el de DESI, no el de Planck. Con el fondo
  marginalizado no queda cota.
- **La cota de estocasticidad propiamente dicha, la de un residuo que varía entre redshifts
  (θ ≳ 1), es otra:**
  - 1.4–1.5 × 10⁻² con θ = 1 y 1.3 × 10⁻² con θ = 10, con el fondo fijo o perfilado;
  - 4.6 × 10⁻² y 2.2 × 10⁻² si el fondo se marginaliza.
