# El reloj de discretización del Modelo 1 (exploración)

**Predicción, escrita antes de correr `reloj_m1.py`:**
- Con n = 10⁵, los tres relojes coinciden en supervivencia a ±3 puntos (≈ 25 %).
- Con n = 10⁴, el de Das (uniforme en a) se separa más de 3 puntos de los otros dos.
- Ningún χ² mínimo baja de 29.84 + 4.

Exploración: no cambia ningún resultado ni ningún paper. Script `reloj_m1.py`, salida `reloj_m1.json`.

## Pregunta

`ahmed_sim.run_das` da pasos uniformes en a ([D1] Sec. 4.2) desde a = 10⁻⁴. La caminata
S_i = S_{i−1} + α √ΔV_i ξ_i es un proceso de Wiener en V, así que en el límite continuo la fracción
que sobrevive (H² > 0 hasta hoy) no debería depender de cómo se reparten los pasos. ¿Es así a la
resolución que se usa?

## Montaje

**Común:**
- α = 0.01, arranque en a = 10⁻⁴, S_0 = 0, K = 8π/3, ω_m = 0.14301, ω_r de `ahmed_sim`.
- 2000 realizaciones. La realización i usa `default_rng(20261016 + i)`, y su paso k usa su normal k
  en los tres relojes.
- Λ fija dentro de cada paso y comprobación de H² > 0 en cada subpaso RK4, como en `run_das`.

**Los tres relojes:**
- **(a) uniforme en a:** la malla de `run_das`, `arange(i_ini, n+1)/n`.
- **(b) uniforme en ln a:** el mismo número de pasos.
- **(c) uniforme en V,** el 4-volumen de la propia simulación:
  - ΔV = (V_ref(1) − V(10⁻⁴))/n, con V_ref el 4-volumen del mismo universo sin energía oscura.
  - Cada realización avanza desde su propio estado con la ley de potencias local V ∝ a^p.
  - El ΔV real por paso es 1.000 veces el objetivo (percentiles 5–95: 0.99999–1.00001).
  - Los pasos se recortan para caer en los z de DESI y en a = 1.

`run_das` no admite otra malla, así que su bucle se copia generalizado en `reloj_m1.py` (no se toca
`ahmed_sim.py`).

**BAO:**
- χ² con el r_d de cada run (acumulado dentro del run para los 25 ω_m de `rd_real.py`) y ω_m a
  ±1.6 % de Planck, con `rd_real.minimise`: el criterio (c) sin ventana de 3e.
- En el reloj V, r_d se integra por subpasos, porque sus primeros pasos son largos.

**Comprobaciones** (`reloj_m1.json` → `check_*`):
- **§1 con sus propias semillas, con `run_das` sin tocar:** 263/2000 = 13.15 % (n = 10⁴) y
  275/1000 = 27.5 % (n = 10⁵). Se reproduce exactamente.
  - Ojo: §1 usó **una** semilla para todas las columnas, `default_rng(20261016)`, y ω_m = 0.14; no
    "20261016 + índice" ni 0.14301. Las corridas de abajo usan lo pedido, así que el reloj (a) no
    da exactamente los números de §1. Da 14.3 % y 26.7 %, a ~1σ de aquellos.
- **El bucle generalizado con la malla y el generador de `run_das` reproduce `run_das` bit a bit**
  (n = 10⁴, 2000 columnas): mismos supervivientes, H0, tiempo conforme en los z de DESI y r_d
  acumulado, todo con diferencia 0.

## Resultados (α = 0.01, 2000 semillas, las mismas en todo)

| reloj | n | sobreviven | pasos antes de a_eq = 2.93 × 10⁻⁴ | H0, p16 / p50 / p84 (σ) | Ω_Λ(0), p16 / p50 / p84 (σ) | χ²_c mínimo (ω_m) | χ²_c < 33.84 |
|---|---|---|---|---|---|---|---|
| (a) a | 10⁴ | **286 (14.3 ± 0.8 %)** | 1 | 31.6 / 40.4 / 50.0 (9.2) | −0.43 / 0.12 / 0.43 (0.51) | 130.8 (0.1439) | 0 |
| (b) ln a | 10⁴ | **525 (26.3 ± 1.0 %)** | 1165 | 31.9 / 39.6 / 49.0 (9.0) | −0.41 / 0.09 / 0.40 (0.53) | 139.3 (0.1407) | 0 |
| (c) V | 10⁴ | **334 (16.7 ± 0.8 %)** | 0 | 32.5 / 41.0 / 53.6 (13.9) | −0.35 / 0.15 / 0.50 (0.53) | 628.2 (0.1428) | 0 |
| (a) a | 10⁵ | **534 (26.7 ± 1.0 %)** | 19 | 32.7 / 40.3 / 48.9 (8.8) | −0.33 / 0.12 / 0.40 (0.52) | 282.2 (0.1453) | 0 |
| (b) ln a | 10⁵ | **387 (19.4 ± 0.9 %)** | 11652 | 32.7 / 40.1 / 49.9 (10.4) | −0.34 / 0.11 / 0.43 (0.52) | 155.9 (0.1407) | 0 |
| (c) V | 10⁵ | **357 (17.9 ± 0.9 %)** | 0 | 31.6 / 40.2 / 51.4 (12.2) | −0.43 / 0.12 / 0.46 (0.67) | 469.2 (0.1453) | 0 |
| (b) ln a | 3 × 10⁵ | **359 (18.0 ± 0.9 %)** | 34957 | 31.5 / 40.0 / 49.0 (9.4) | −0.44 / 0.11 / 0.40 (0.76) | 244.5 (0.1407) | 0 |
| (b) ln a | 10⁶ | **346 (17.3 ± 0.8 %)** | 116525 | 31.8 / 39.8 / 48.9 (9.4) | −0.42 / 0.10 / 0.40 (0.65) | 264.5 (0.1453) | 0 |

- Los errores son binomiales, ±√(p(1−p)/2000).
- Las filas de ln a con n = 3 × 10⁵ y 10⁶ salen de `reloj_m1_conv.py` (sección "Convergencia en ln a").
- Los pasos del reloj V varían por realización: mediana en las supervivientes 8639 (n = 10⁴) y 89638
  (n = 10⁵).
- **Referencia de §1** (otras semillas, ω_m = 0.14, uniforme en a): n = 10⁶ da 69/300 = 23.0 ± 2.4 %.
- **Supervivientes en común**, con las mismas semillas (`overlap_survivors`): pocas. Por ejemplo,
  ln a con n = 10⁴ y con n = 10⁵ comparten 99 de 525 y 387. Sobrevivir depende del camino, no de la
  semilla.

## Por qué difieren

La distribución de S a un V dado no depende del reloj; la supervivencia sí. Depende de dos efectos
de la discretización, de signo opuesto, que solo desaparecen cuando el paso resuelve el tiempo de
Hubble en todo el recorrido:

1. **Λ fija en un paso largo mientras ρ decae.**
   - Uniforme en a con n = 10⁴, el primer paso duplica a (Δln a = 0.69). En ese paso la radiación
     cae un factor 16 y Λ no cambia, así que un Λ negativo inocuo al principio del paso se vuelve
     letal dentro de él.
   - Por eso (a) con n = 10⁴ mata de más (14.3 %), y con n = 10⁵ (Δln a = 0.1 al principio) sube a
     26.7 %.
2. **Vigilancia discreta de un proceso continuo.** Con pasos finos en todo el recorrido, como (b),
   afinar detecta más cruces del umbral: la supervivencia baja de 26.3 % a 19.4 % al pasar de
   Δln a = 9.2 × 10⁻⁴ a 9.2 × 10⁻⁵. Es lo esperable con un error ∝ √Δ. Con solo dos puntos, una
   extrapolación ∝ √Δ daría ≈ 16 % en el continuo; es una estimación, no un resultado.
3. **El reloj V no resuelve el universo temprano.**
   - V(10⁻⁴)/V(1) = 3 × 10⁻²⁶. Con ΔV = V/n, el primer paso va de a = 10⁻⁴ a a = 0.216 (n = 10⁴,
     z = 3.6) o 0.147 (n = 10⁵, z = 5.8).
   - Durante ese paso, Λ = 0 (S_0 = 0) y no hay ningún paso antes de la igualdad: los colapsos
     tempranos no pueden ocurrir.
   - Al final del primer paso, Λ aparece de golpe con toda la fluctuación acumulada,
     KαS/V ∼ Kα ξ/√V; la mayoría de las que mueren lo hacen ahí.
   - Su supervivencia, 16.7 % y 17.9 %, es estable con n, pero mide otra cosa: solo los colapsos a
     z ≲ 6. Su dispersión en H0 es mayor (σ = 12–14 frente a 9–10).
4. **Las supervivientes se parecen en los tres relojes:** H0 con mediana 39.6–41.0 y Ω_Λ(0) con
   mediana 0.09–0.15.
5. **BAO no cambia:** χ²_c mínimo entre 131 y 628, y ninguna por debajo de 33.84 en ningún reloj.

**Qué implica para lo publicado** (sin cambiarlo):
- La fracción de runs que completan (27.6 % con α = 0.01, uniforme en a, n = 10⁵) no es un número
  convergido. A n = 10⁵ va de 17.9 % a 26.7 % según el reloj, y con pasos uniformes en ln a todavía
  baja al refinar. Lo mismo vale para el 32.3 % de [D1], también uniforme en a.
- Los resultados de BAO de la sección 3 cuentan sobre runs completados y dan 0, y las propiedades de
  las supervivientes apenas dependen del reloj. Con 2000 semillas, ningún reloj produce una
  realización que pase B.

## Predicción

- **Con n = 10⁵, los tres coinciden a ±3 puntos: no se cumplió.** Dan 26.7, 19.4 y 17.9 %, con 8.8
  puntos entre el mayor y el menor.
- **Con n = 10⁴, el de Das se separa más de 3 puntos de los otros dos: no se cumplió.** Se separa
  12.0 puntos de ln a (14.3 frente a 26.3 %), pero solo 2.4 de V (16.7 %).
- **Ningún χ² mínimo baja de 29.84 + 4: se cumplió.** El mínimo es 130.8.

## Convergencia en ln a

**Predicción, escrita antes de correr `reloj_m1_conv.py`:**
- ln a con n = 10⁶ da 17.2 ± 1.0 %.
- p∞ = 16 ± 1 %.
- Las dos estimaciones del continuo coinciden a ±1.5 puntos.

Script `reloj_m1_conv.py`, que importa `reloj_m1.py`; salida `reloj_m1_conv.json`. Mismas 2000
semillas y α = 0.01.

### Supervivencia con pasos uniformes en ln a

| n | Δln a | sobreviven | H0, p16 / p50 / p84 | Ω_Λ(0), p16 / p50 / p84 | χ²_c mínimo | tiempo |
|---|---|---|---|---|---|---|
| 10⁴ | 9.2 × 10⁻⁴ | 525 (26.3 ± 1.0 %) | 31.9 / 39.6 / 49.0 | −0.41 / 0.09 / 0.40 | 139.3 | 5 s |
| 10⁵ | 9.2 × 10⁻⁵ | 387 (19.4 ± 0.9 %) | 32.7 / 40.1 / 49.9 | −0.34 / 0.11 / 0.43 | 155.9 | 48 s |
| 3 × 10⁵ | 3.1 × 10⁻⁵ | 359 (18.0 ± 0.9 %) | 31.5 / 40.0 / 49.0 | −0.44 / 0.11 / 0.40 | 244.5 | 139 s |
| 10⁶ | 9.2 × 10⁻⁶ | **346 (17.3 ± 0.8 %)** | 31.8 / 39.8 / 48.9 | −0.42 / 0.10 / 0.40 | 264.5 | 468 s |

La supervivencia baja de forma monótona al refinar, cada vez más despacio. H0 y Ω_Λ(0) de las
supervivientes no cambian, y ningún χ²_c baja de 33.84.

### Dos estimaciones del límite continuo

**(i) Ajuste p(Δ) = p∞ + c·√Δ** a los cuatro n de ln a, por mínimos cuadrados con pesos binomiales:
- **p∞ = 16.2 ± 0.7 %**, con c = 3.30.
- χ² = 0.02 para 2 grados de libertad; los residuos son < 0.1 puntos.
- El χ² tan pequeño no indica un ajuste mejor de lo esperable. Las cuatro corridas comparten
  semillas, así que sus errores están muy correlacionados, y el ajuste los trata como
  independientes. El ±0.7 es el error binomial de un solo n, que es la incertidumbre honesta de
  p∞ con 2000 semillas.

**(ii) Corrección de continuidad de Broadie–Glasserman–Kou.**
- Para una barrera vigilada a tiempos discretos, el run discreto con la barrera acercada al proceso
  en β·σ_paso, con β = −ζ(1/2)/√(2π) = 0.5826, aproxima la vigilancia continua.
- Aquí el proceso vigilado es u = KαS/V, con S de Wiener en V, así que σ_paso = Kα√ΔV/V: es la
  σ_paso = α√ΔV del encargo, pasada a unidades de u.
- La condición de muerte al final de cada paso pasa a ser ω_m a⁻³ + ω_r a⁻⁴ + u ≤ 0.5826·σ_paso.
- **Es aproximada:** la barrera −(ω_m a⁻³ + ω_r a⁻⁴) se mueve con el tiempo, u no es un movimiento
  browniano puro (es S/V) y Λ queda fija entre pasos.

| run con BGK | sobreviven | sin corrección |
|---|---|---|
| ln a, n = 10⁵ | **314 (15.7 ± 0.8 %)** | 19.4 % |
| ln a, n = 10⁴ (comprobación) | 324 (16.2 ± 0.8 %) | 26.3 % |

La corrección BGK lleva n = 10⁴ y n = 10⁵ casi al mismo valor (16.2 y 15.7 %), aunque sin corregir
difieren 7 puntos. Las dos quedan junto a la extrapolación √Δ.

### Resultado

- **Límite continuo de la supervivencia con α = 0.01: ≈ 16 %**: 16.2 ± 0.7 % por el ajuste √Δ y
  15.7 ± 0.8 % por BGK con n = 10⁵. Difieren 0.5 puntos.
- Confirma la estimación de dos puntos del apartado "Por qué difieren" (≈ 16 %).
- Sin cambiar nada publicado: la fracción de producción (27.6 %, uniforme en a, n = 10⁵) y la de
  [D1] (32.3 %) quedan unos 11 y 16 puntos por encima del continuo. Esa diferencia es artefacto de
  la discretización (efecto 1 de "Por qué difieren"), no física del modelo.
- Los 0 de BAO de la sección 3 cuentan sobre runs completados y no dependen de esa fracción.

### Predicción

- **ln a con n = 10⁶ da 17.2 ± 1.0 %: cumplida.** Da 17.3 ± 0.8 %.
- **p∞ = 16 ± 1 %: cumplida.** Da 16.2 ± 0.7 %.
- **Las dos estimaciones del continuo coinciden a ±1.5 puntos: cumplida.** Difieren 0.5 puntos
  (16.2 frente a 15.7 %).
