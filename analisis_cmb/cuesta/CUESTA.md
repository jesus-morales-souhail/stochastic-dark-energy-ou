# ¿Cuesta o temblor? Falsos positivos de w₀–wₐ producidos por un residuo OU, BAO solo (exploración)

**Predicción, escrita antes de correr `cuesta.py`:**
- θ ≥ 1 al techo: menos del 5 % de exceso sobre el control.
- θ = 0.1–0.3 al techo marginalizado: entre el 10 y el 30 %.

Exploración: no cambia `results/`, `figures/` ni ningún paper. Script `cuesta.py`; salidas
`cuesta.{json,txt}` y `cuesta_universos.npz` (Δχ², w₀, wₐ, Ω_m, h·r_d y χ² de cada ajuste de cada universo).
No se hace todavía el paso con supernovas y CMB.

## Montaje

**Universo base:** ΛCDM plano con el mejor ajuste de DESI DR2, Ω_m = 0.29740 y h·r_d = 101.538 Mpc
(`FONDO.md`), que predice los 13 números de DESI DR2 (D_V/r_d de BGS; D_M/r_d y D_H/r_d de los otros seis
bins; `desi_dr2_data.load_gaussian_bao_vector`).

**Universos simulados.** Para cada celda (σ_X, θ), 2000 universos
obs_i = pred_i · (1 + S_i X(z_i)) + n, con:
- X(z) una realización OU estacionaria en x = ln(1+z): X ~ N(0, σ_X² e^{−θ|Δx|}), la misma receta que el
  núcleo de `perfil_fondo.py`; las filas D_M y D_H del mismo bin comparten X.
- S por tipo (D_V, D_M, D_H), S = −∂ln D_tipo/∂Ω_Λ, como en `perfil_13.py`.
- n ~ N(0, C) con la covarianza oficial 13 × 13.
- **Los mismos 2000 vectores de ruido n en todas las celdas** (números aleatorios comunes); las
  realizaciones OU salen de una semilla distinta por celda. Así "exceso sobre el control" es una
  diferencia emparejada, y su error es el de esa diferencia.
- **Control σ_X = 0:** solo ruido. No depende de θ, así que es una sola celda al lado de todas.

**Malla:** θ ∈ {0.1, 0.3, 1, 3} × σ_X ∈ {0, cota al 95 % con fondo fijo en ese θ, cota al 95 % con el fondo
marginalizado en ese θ}.
- `perfil_fondo.json` solo tiene θ = 0.001, 0.1, 1, 10, y son las cotas con los 7 α_DV. Como los universos
  usan los 13 números, las cotas se recalculan aquí con la maquinaria de `perfil_13.py` (tratamientos (a)
  y (c), mismos 13 números, misma S, mismo núcleo) en los cuatro θ. En θ = 0.1 y 1 coinciden con la fila
  "13 números" de `TRECE.md` (1.833 × 10⁻², 1.280 × 10⁻²; 0.112, 3.774 × 10⁻²). Las de los 7 α_DV son un
  factor 1.1–1.25 más flojas, con lo que habrían dado algo más de falsos positivos.

| θ | σ_X, fondo fijo (a) | σ_X, fondo marginalizado (c) |
|---|---|---|
| 0.1 | 0.01833 | 0.1121 |
| 0.3 | 0.01548 | 0.0656 |
| 1 | 0.01280 | 0.0377 |
| 3 | 0.01155 | 0.0249 |

**Ajustes.** A cada universo se le ajusta ΛCDM plano (Ω_m, h·r_d) y w₀wₐCDM plano (Ω_m, h·r_d, w₀, wₐ)
por mínimo χ² con la covarianza 13 × 13:
- h·r_d analíticamente (la predicción es lineal en 1/(h·r_d)); Ω_m por Brent acotado (ΛCDM);
  (Ω_m, w₀, wₐ) por Nelder–Mead con cinco arranques, uno de ellos el mejor ΛCDM, de modo que Δχ² ≥ 0.
- Priores de DESI DR2 (sec. 6 del artículo): w₀ ∈ [−3, 1], wₐ ∈ [−3, 2], w₀ + wₐ < 0; Ω_m ∈ [0.01, 0.99].
- `scripts/eos_efectiva.py` (el de `results/eos_cpl_desi_dr2/`) fija Ω_m = 0.315 y H0 = 67.4 y ajusta solo
  (w₀, wₐ) a los 7 α; aquí se reutilizan sus ingredientes (E(z) de CPL, χ² con la covarianza oficial,
  Nelder–Mead multiarranque), pero con los 13 números y Ω_m, h·r_d libres, que es lo que hace DESI.
- Distancias por Simpson en una malla fija (128 subintervalos por tramo): coinciden con `quad` a 6 × 10⁻¹².

**Comprobaciones con los datos reales (13 números):**

| magnitud | aquí | DESI DR2 (TeX de arXiv:2503.14738) |
|---|---|---|
| ΛCDM: Ω_m, h·r_d | 0.2975, 101.54 Mpc | 0.2975 ± 0.0086, 101.54 ± 0.73 Mpc |
| χ²(ΛCDM), 11 gdl | 10.27 | — |
| Δχ² = χ²(ΛCDM) − χ²(w₀wₐ) | **4.65** | **4.7** (`tab:dchisq`, "DESI", Δχ²_MAP = −4.7, 1.7σ, ΔDIC = −0.8) |
| w₀, wₐ de mejor ajuste | −0.18, −2.72 (Ω_m = 0.387) | medias marginales: w₀ = −0.48 ^{+0.35}_{−0.17}, wₐ < −1.34 (68 %), "cortadas por los priores" (`eq:w0wa_DESI`) |

DESI no publica el mejor ajuste (MAP) de w₀, wₐ con DESI solo, solo las medias marginales; el nuestro está
sobre la misma degeneración, en el cuadrante w₀ > −1, wₐ < 0, con wₐ cerca del borde del prior.

**Falso positivo:** universo con Δχ² ≥ 4.7 **y** w₀ > −1 **y** wₐ < 0.

## Resultados

**Fracción de falsos positivos por celda (2000 universos; control σ_X = 0: FP = 4.75 ± 0.48 %):**

| θ | fondo | σ_X | FP | control | **exceso** | ± | P(Δχ² ≥ 4.7) | cuadrante | mediana Δχ² | ⟨χ²(ΛCDM)⟩ | FP con χ²(ΛCDM) ≤ 19.7 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.1 | fijo | 0.0183 | 6.3 % | 4.75 % | **+1.6 %** | 0.4 | 10.5 % | 0.49 | 1.49 | 11.5 | 4.8 % |
| 0.3 | fijo | 0.0155 | 8.5 % | 4.75 % | **+3.7 %** | 0.6 | 13.7 % | 0.49 | 1.69 | 12.1 | 5.9 % |
| 1 | fijo | 0.0128 | 12.2 % | 4.75 % | **+7.5 %** | 0.7 | 22.1 % | 0.49 | 2.13 | 13.3 | 7.6 % |
| 3 | fijo | 0.0116 | 17.6 % | 4.75 % | **+12.8 %** | 0.8 | 29.2 % | 0.50 | 2.58 | 15.3 | 9.0 % |
| 0.1 | marginalizado | 0.1121 | 32.3 % | 4.75 % | **+27.5 %** | 1.1 | 64.9 % | 0.46 | 8.14 | 33.8 | 3.8 % |
| 0.3 | marginalizado | 0.0656 | 34.1 % | 4.75 % | **+29.3 %** | 1.1 | 62.5 % | 0.49 | 7.43 | 31.5 | 4.9 % |
| 1 | marginalizado | 0.0377 | 34.7 % | 4.75 % | **+29.9 %** | 1.1 | 62.8 % | 0.51 | 7.21 | 30.8 | 5.1 % |
| 3 | marginalizado | 0.0249 | 34.5 % | 4.75 % | **+29.8 %** | 1.1 | 62.2 % | 0.50 | 6.91 | 30.8 | 5.0 % |

- "cuadrante": fracción de universos cuyo mejor ajuste cae en w₀ > −1, wₐ < 0.
- ⟨χ²(ΛCDM)⟩: media del χ² del ajuste ΛCDM (11 grados de libertad; con solo ruido, 11.0; DESI real, 10.3).
- Última columna: falsos positivos que además tienen un ajuste ΛCDM aceptable (χ² ≤ 19.7, el 95 % de χ²₁₁);
  en el control, 3.75 %.

**Control.** P(Δχ² ≥ 4.7) = 7.4 %, algo menos que el 9.5 % de una χ²₂ pura: la caja de priores (wₐ ≥ −3,
w₀ + wₐ < 0) recorta la cola; el 3.5 % de los ajustes acaba en wₐ = −3. El cuadrante de DESI es algo más
propenso que el contrario incluso sin OU: P(Δχ² ≥ 4.7 | cuadrante) = 9.7 % frente a 5.1 % fuera, porque la
degeneración w₀–wₐ permite excursiones más largas hacia wₐ → −3. Por eso FP (4.75 %) es más de la mitad
de P(Δχ² ≥ 4.7).

## Qué θ engaña más y por qué

1. **Al techo con fondo fijo, engaña más el θ grande: +1.6 % (θ = 0.1) → +3.7 % → +7.5 % → +12.8 %
   (θ = 3), aunque la amplitud permitida baja de 0.018 a 0.012.**
   - Con θ pequeño, el residuo es casi el modo S, y ese modo es la dirección de Ω_m (coseno 0.99996,
     `FONDO.md`). ΛCDM lo absorbe con Ω_m y h·r_d, y a w₀wₐ no le queda nada: ⟨χ²(ΛCDM)⟩ solo sube de 11.0
     a 11.5. La cota con fondo fijo a θ → 0 era un factor de Occam, no una señal; aquí se ve que tampoco es
     una cuesta que w₀wₐ pueda confundir con energía oscura dinámica.
   - Con θ grande, el residuo se descorrelaciona entre bins y ya no se parece a ningún cambio de fondo. Es
     ruido extra que la covarianza oficial no contiene: ⟨χ²(ΛCDM)⟩ sube a 15.3 (+4.3) y el χ² de w₀wₐ solo
     a 11.3 (+2.1). Dos parámetros más absorben más de ese ruido extra, y Δχ² se infla (mediana 1.4 → 2.6).
     Es exactamente lo que hace un temblor: no empuja hacia ningún sitio, ensancha.
2. **Al techo marginalizado, los cuatro θ dan lo mismo: +27.5 a +29.9 %, P(Δχ² ≥ 4.7) ≈ 62–65 %, mediana
   Δχ² ≈ 7–8.**
   - La cota marginalizada es, por construcción, la amplitud a la que la parte del residuo que (Ω_m, h·r_d)
     no puede absorber llega a ΔlnL = −1.92. Esa parte es justo la que ve w₀wₐ, así que las cuatro celdas
     están "al mismo nivel de detección" aunque σ_X vaya de 0.11 a 0.025.
   - Pero son universos con ⟨χ²(ΛCDM)⟩ ≈ 31–34 para 11 grados de libertad: solo el 31–33 % tiene un ajuste
     ΛCDM aceptable (χ² ≤ 19.7), frente al 95 % del control. DESI tiene χ² = 10.3. Exigiendo χ²(ΛCDM) ≤ 19.7,
     el falso positivo marginalizado cae al 3.8–5.1 %, igual que el control (3.75 %). **El techo marginalizado
     no es un techo realista para DESI: un residuo así se habría visto en el χ² de ΛCDM antes que en w₀wₐ.**
3. **Nunca es una cuesta.** En todas las celdas, la fracción de mejores ajustes en el cuadrante w₀ > −1,
   wₐ < 0 se queda en 0.46–0.51, como en el control (0.49), y el cociente FP/P(Δχ² ≥ 4.7) se queda en
   0.50–0.60 (control 0.65). El residuo OU, con cualquier θ, no tiene signo preferido: no produce la
   dirección de DESI más veces de lo que la produce el ruido. Lo que produce es un Δχ² más ancho, y el
   cuadrante lo pone el azar. Los falsos positivos que salen se parecen al ajuste real (mediana w₀ ≈ −0.3 a
   −0.5, wₐ en el borde −3), porque es donde la degeneración lleva a cualquier exceso.
4. **Con fondo fijo y χ²(ΛCDM) aceptable, el exceso real queda en +1 % (θ = 0.1) a +5 % (θ = 3)** sobre el
   3.75 % del control. Es el número que importa para el paso siguiente.

## Predicción

- **θ ≥ 1 al techo, menos del 5 % de exceso: no se cumplió.** Con el fondo fijo, +7.5 % (θ = 1) y +12.8 %
  (θ = 3); con el fondo marginalizado, +30 %. Solo se cumple si además se exige un χ²(ΛCDM) aceptable
  (+3.9 % y +5.2 %, el segundo en el borde). Me equivoqué en el mecanismo: pensé que un residuo
  descorrelacionado no engaña a un modelo suave de dos parámetros, y no engaña en dirección, pero sí infla
  Δχ² porque el ajuste lo hace con una covarianza que no lo incluye.
- **θ = 0.1–0.3 al techo marginalizado, entre el 10 y el 30 %: cumplida.** +27.5 % y +29.3 %. Pero el
  acierto vale poco: ese techo lleva χ²(ΛCDM) ≈ 32 para 11 grados de libertad, y DESI tiene 10.3.

## Para decidir el paso con supernovas y CMB

- Lo que hay que llevar es el techo con fondo fijo (o perfilado), no el marginalizado, y mirar a la vez el
  χ² de ΛCDM: en BAO solo, el exceso útil es del 1–5 % según θ, y crece con θ.
- La pregunta abierta es si SN y CMB, que fijan Ω_m y rompen la degeneración w₀–wₐ, dejan que un temblor
  de θ ≈ 1–3 con σ_X ≈ 0.012 produzca el Δχ² = −12.5 de DESI+CMB. En BAO solo no produce el signo; habría
  que ver si lo produce la combinación.
