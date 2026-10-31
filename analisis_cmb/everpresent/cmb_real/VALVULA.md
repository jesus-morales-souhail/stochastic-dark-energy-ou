# La válvula tanh del Modelo 2 (exploración)

**Predicción, escrita antes de correr `valvula.py`:**
- 8163 se descarta.
- 1540, 11072 y 16616 empeoran, con un factor en H entre 0.58 y 0.70.
- 12909 queda casi igual.
- Se descartan más del 20 % de las 20000.

Exploración: no cambia ningún resultado ni ningún paper. Script `valvula.py`, salida `valvula.json`.

## Qué se mide

`bao_sim.py` y `rd_real.py` usan Ω_de = tanh(Ω̂). Así 1 − Ω_de ∈ (0, 2) y
H² ≥ ½·100²[ω_m(1+z)³ + ω_r(1+z)⁴]: el Modelo 2 nunca llega a H² ≤ 0, mientras que el Modelo 1 sí y
"crashea".

**Sin válvula:**
- Ω_de = Ω̂.
- Una realización se descarta ("crash") si Ω̂ ≥ 1 en algún punto guardado. Son 458 puntos: los 401
  de `bao_sim.py` (z ∈ [0, 3] y z_rec) y los 57 de `rd_real.py` (z ∈ [10³, 10⁵]).
- El criterio (c) se recalcula con el propio `rd_real.model2()`, con tanh sustituida por la
  identidad y su npz redirigido a un directorio temporal. No se escribe nada en `resultados/`.

**Comprobaciones:**
- Ω̂ regenerado (misma semilla, mismo Cholesky) reproduce `bao_m2.json` y `bao_m2_bestfit.npz`:
  Δχ² máximo = 0.0; A = 1, B = 1 y proxy = 4569, iguales a los guardados.
- `rd_real.model2()` pasado por el mismo envoltorio, con la tanh real, reproduce
  `resultados/rd_real_m2.npz` bit a bit (las 8 matrices).

## (a) Cuánto trabaja la válvula (20000 realizaciones, α = 0.8824)

| condición | n | fracción |
|---|---|---|
| \|Ω̂\| > 1 en λ_rec | 5202 | 26.0 % (gaussiana con σ = 0.8824: 25.7 %) |
| Ω̂ ≥ 1 en λ_rec (crash) / Ω̂ ≤ −1 | 2616 / 2586 | 13.1 % / 12.9 % |
| \|Ω̂\| > 1 en algún λ ∈ [−ln 4, 0] | 9619 | 48.1 % |
| Ω̂ ≥ 1 en algún λ ∈ [−ln 4, 0] | 4860 | 24.3 % |
| Ω̂ ≥ 1 en algún punto de los 401 (z ≤ 3 y z_rec) | 6821 | 34.1 % |
| \|Ω̂\| > 1 en algún punto de los 458 | 17496 | 87.5 % |
| **Ω̂ ≥ 1 en algún punto de los 458 (se descartan sin válvula)** | **11738** | **58.7 %** |

La válvula actúa en casi todas las realizaciones: en el 87.5 %, |Ω̂| supera 1 en algún punto.

## (b) Supervivientes en z_rec

Las 65 supervivientes del Modelo 2 de `supervivientes.json`:
- 5 tendrían crash en z_rec.
- 22 lo tendrían en algún punto guardado.

La fila de cada una está en `valvula.json` → `b_survivors_model2`. Las 6 del Modelo 1 no tienen
válvula. El factor en H es frente a H sin energía oscura, 1/√(1 − Ω_de).

| índice | Ω̂(z_rec) | Ω_de(z_rec) = tanh | factor en H con válvula | factor en H sin válvula | Ω̂ máximo (458 puntos) | sin válvula |
|---|---|---|---|---|---|---|
| 1540 | −1.957 | −0.961 | 0.714 | **0.582** | 0.56 | sigue; empeora |
| 8163 | +1.270 | +0.854 | 2.616 | **crash** | 1.34 | se descarta (en z_rec) |
| 11072 | −1.321 | −0.867 | 0.732 | **0.656** | 1.57 | se descarta: Ω̂ = 1.57 en otro punto |
| 12909 | −0.261 | −0.255 | 0.893 | **0.890** | 0.69 | sigue; casi igual en z_rec |
| 16616 | −1.025 | −0.772 | 0.751 | **0.703** | 0.66 | sigue; empeora |

## (c) Criterio (c) sin válvula

| | con válvula (`rd_real.json`) | sin válvula |
|---|---|---|
| realizaciones descartadas (crash) | 0 | **11738 (58.7 %)** |
| quedan | 20000 | 8262 |
| pasan la ventana | 444 | 246 |
| **(c) con ventana: A / B** | **0 / 0** | **0 / 0** (mejor χ² = 101.2) |
| (c) sin ventana: A / B | 6 / 6 | 6 / 6 |

Qué pasa con las 5 del criterio (c), con χ²_c (ω_m a ±1.6 %, r_d propio) con y sin válvula:

| índice | χ²_c con válvula | χ²_c sin válvula |
|---|---|---|
| 1540 | 27.9 | 1209.1 |
| 8163 | 24.7 | descartada |
| 11072 | 23.9 | descartada |
| 12909 | 27.2 | 63.2 (ya no pasa B) |
| 16616 | 9.0 | 500.3 |

**Las 6 que pasan (c) sin ventana y sin válvula son otras.** Ninguna está entre las 6 que pasan con
válvula. Sin válvula mejoran mucho en BAO:

| índice | χ²_c con válvula | χ²_c sin válvula | Ω̂(z_rec) | ventana | ¿entre las 65? |
|---|---|---|---|---|---|
| 14006 | 133.7 | 24.7 | −1.44 | no | sí |
| 14159 | 84.1 | 16.4 | −0.34 | no | sí |
| 16309 | 289.7 | 13.9 | −0.18 | no | sí |
| 16725 | 86.6 | 26.9 | −0.41 | no | no |
| 16756 | 281.7 | 19.4 | −1.13 | no | no |
| 17590 | 93.2 | 13.9 | −0.55 | no | sí |

Las 6 tienen Ω̂(z_rec) < 0. Sin la válvula, la energía oscura negativa en la recombinación es más
negativa, H allí es menor y r_d mayor, y eso las ayuda en BAO con ω_m de Planck. Ninguna pasa la
ventana, así que el criterio (c) completo sigue en 0.

## Conclusión

**La válvula no rescata a nadie; sin ella el Modelo 2 queda más lejos de Planck.**
- De las 5 del criterio (c), 2 se descartan (8163, 11072). Las otras 3 empeoran en BAO: 1540 y
  16616 mucho, y 12909 deja de pasar B. Además, 1540 y 16616 se alejan más de ΛCDM en H(z_rec).
- El criterio (c) con ventana sigue en 0 con o sin válvula.
- Matiz: sin válvula, 6 realizaciones distintas pasan (c) sin ventana, con χ²_c de 84–290 a 14–27.
  Lo hacen precisamente con Ω_de más negativa en la recombinación, ninguna pasa la ventana, y el
  recuento total sin ventana es el mismo (6).

**Predicción:**
- 8163 se descarta: **cumplida**.
- 1540, 11072 y 16616 empeoran con factor en H entre 0.58 y 0.70: **cumplida en el factor** (0.582,
  0.656, 0.703; el último, en el borde por 0.003). Pero 11072 no solo empeora: se descarta, porque
  Ω̂ llega a 1.57 en otro punto guardado.
- 12909 casi igual: **cumplida en H(z_rec)** (0.893 → 0.890). En BAO no: su χ²_c sube de 27.2 a 63.2.
- Más del 20 % descartadas: **cumplida**. Son el 58.7 % con los 458 puntos, el 34.1 % con los 401 de
  BAO y el 24.3 % solo con z ≤ 3.
