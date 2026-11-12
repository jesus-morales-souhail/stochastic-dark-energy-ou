# ¿Descarta el CMB completo (ISW y lente) los caminos que sobreviven en los techos de EXACTO.md? (exploración)

**Predicción, escrita antes de correr `cmb_techos.py`:**
- Supervivientes con ρ_DE ≤ 0 en algún punto: 30–60 % en cada techo.
- Entre los que sí se pueden correr: mediana del Δχ² total < 1 en los tres techos;
  menos del 20 % con Δχ² total > 4.
- Distancia en picos, p95 < 1 en los tres techos.
- El máximo está en TT ℓ < 10 (ISW tardío), no en la lente.
- Lectura: el CMB completo no mueve las cotas de EXACTO.

Exploración: no cambia `results/`, `figures/` ni ningún paper. Script `cmb_techos.py` (etapas `prep`, `run`, `check_edge`);
salidas `cmb_techos_prep.json`, `cmb_techos.json`, `cmb_techos_cell{i}.npz`, `cmb_techos_progress.log`. Tiempo: 5 s por techo
de preparación; 3 min de CAMB (107 realizaciones, 8 procesos). **El resultado está condicionado a los caminos que CAMB
representa bien**; las fracciones excluidas se cuentan abajo, techo por techo.

## Montaje

- **Celdas:** los techos exactos de `EXACTO.md` (BAO + θ* + ω_m, 5 %): (θ, σ_X) = (1, 0.43), (3, 0.29), (10, 0.19). No son
  celdas de `exacto.py` sino cruces interpolados; como los caminos son lineales en σ_X a normales fijas, cada techo se genera con
  la semilla de la celda de `exacto.py` del mismo θ y σ_X más próximo (0.5, 0.3, 0.2), en la misma malla (3506 nodos en
  u = ln(1+z) hasta z*, con los z de DESI exactos).
- **Supervivientes:** χ²_min ≤ 17.238 con el ajuste de (h, ω_m) a las 15 filas (`exacto.fit`).
- **Suavizado** como `rapido.py`: gaussiana de desviación Δ en ln a (0.05 para θ = 1 y 3; 0.0125 para θ = 10), núcleo
  renormalizado en la malla. La malla de `exacto.py` empieza en a = 1, así que cada camino se prolonga a ln a ∈ (0, 0.4] con
  la misma recursión AR(1) (la OU estacionaria es reversible; normales nuevos, semilla [SEED, celda, 1]) para que a = 1 sea
  interior al núcleo, como el relleno de `ou_cmb.py`. Varianza conservada: 0.91 / 0.83 / 0.87 (analítica 0.95 / 0.85 / 0.87).
  χ²_min se recalcula con el camino suavizado (mediana del cambio −0.9 / −1.4 / −0.3; 16–84 %: −3.3 a +1.9 en θ = 1).
- **ρ_DE ≤ 0:** ρ_DE/ρ_crit,0 = (h² − ω_m − ω_r)/h² + X_s(z) − X_s(0) con el (h, ω_m) ajustado al camino suavizado; si es ≤ 0
  en algún punto de la tabla (10⁻³ ≤ a ≤ 1), no entra en CAMB (DarkEnergyPPF integra ρ a partir de w). Es un sesgo hacia los
  caminos mansos.
- **CAMB:** DarkEnergyPPF con la tabla w(a) = −1 − ⅓ d ln ρ_DE/d ln a (derivada en la malla rellena, tabla en 10⁻³ ≤ a ≤ 1),
  lente lineal, H0 = 100h y ω_m de cada realización, ω_b = 0.02237 (`rdcalc.OMB`, el de `exacto.py`), ω_c = ω_m − ω_b − ω_ν
  (m_ν = 0.06 eV), τ, n_s, A_s de Planck 2018, ℓ_max = 2500. **Referencia** por realización: ΛCDM con los mismos ω_b, ω_c,
  m_ν y H0 resuelto para que el 100θ* de CAMB sea el de la realización (`ou_cmb.solve_H0_for_thetastar`).
- **Estadísticos**, realización menos referencia: Δχ² frente a Planck 2018 en TT (2 ≤ ℓ < 30) y en φφ (bandas), como
  `RAPIDO.md` (`estadisticos.chi2`, `planck_datos`); distancia en picos TT/TE/EE, ℓ ∈ [30, 2500], de `picos.distance`
  (conjunta). Total = suma. Dónde cae el máximo de |Δ/σ|: el mayor entre el máximo por bin de TT ℓ < 30, el máximo por banda
  de φφ y la raíz del máximo por ℓ de la distancia conjunta en picos (calculado aquí, `joint_per_l`).
- **Comprobación por dos caminos, en todas las realizaciones:** ρ_DE que integra CAMB (`get_dark_energy_rho_w`) frente a
  la del camino. **Fiel** si máx_a |ρ_DE(CAMB) − ρ_DE(camino)|/ρ_total(a) < 10⁻³. El plan inicial comprobaba solo la tabla
  (10⁻³ ≤ a ≤ 1). **La comprobación por debajo de la tabla** (10⁻⁶ ≤ a < 10⁻³, donde CAMB extrapola w = w(10⁻³) y el camino
  vale Ω_Λ + X(z) − X(0) con X de la malla hasta z* y constante más allá) **se añadió después de ver el artefacto de θ = 10**
  (etapa `check_edge`, tras la corrida); no formaba parte del plan. El error relativo de ρ_DE se da en otra columna. Nada se
  tira en silencio: fieles y no fieles por separado; la conclusión solo usa las fieles. No se ha cambiado la tabla de w ni el
  suavizado para rescatar caminos.
- **Lanzador:** un proceso por realización con límite de tiempo (una primera corrida con `Pool` se colgó: un `STOP` de Fortran
  dentro de CAMB mata al worker y la pool espera para siempre). Dos realizaciones de θ = 10 murieron por segfault (exit −11)
  y cuentan como fallos, no como resultados.

## Resultados

**Etapas 1–3:**

| techo | supervivientes (χ²_min ≤ 17.24) | pasan tras suavizar | ρ_DE ≤ 0 (de los que pasan) | corribles | corridas | fallos CAMB |
|---|---|---|---|---|---|---|
| θ = 1, σ_X = 0.43 | 91 (4.6 %) | 82 | 55 (67 %) | 27 | 27 | 0 |
| θ = 3, σ_X = 0.29 | 93 (4.7 %) | 83 | 23 (28 %) | 60 | 40 | 0 |
| θ = 10, σ_X = 0.19 | 100 (5.0 %) | 90 | 4 (4.4 %) | 86 | 40 | 2 (segfault) |

**Fidelidad (en H²):**

| techo | corridas con resultado | fieles | no fieles | error en H² de las fieles, mediana / máx | error relativo de ρ_DE, mediana / máx | \|w + 1\| máx, mediana / p95 / máx |
|---|---|---|---|---|---|---|
| θ = 1 | 27 | **27** | 0 | 2.2 × 10⁻⁵ / 5.1 × 10⁻⁵ | 1.2 × 10⁻⁴ / 0.29 | 2.4 / 8.3 / 141 |
| θ = 3 | 40 | **40** | 0 | 2.6 × 10⁻⁵ / 8.3 × 10⁻⁵ | 1.7 × 10⁻⁴ / 1.4 | 3.3 / 7.4 / 349 |
| θ = 10 | 38 | **31** | 7 | 3.6 × 10⁻⁴ / 7.8 × 10⁻⁴ | 1.5 × 10⁻³ / 0.88 | 8.3 / 18.5 / 234 |
| θ = 10, fallos de ejecución | 2 (segfault, exit −11) | — | — | sin resultado | sin resultado | — |

Los dos segfault de θ = 10 van en su propia fila: un fallo de ejecución no dice nada de esa historia, ni a favor ni en
contra, y no entran en ninguna de las dos clases ni en los sesgos (se restan del número de corridas: 40 lanzadas, 38 con
resultado).

- Dentro de la tabla todas las 105 son fieles en H² (máximo 7.8 × 10⁻⁴, en θ = 10, cerca del umbral). El error relativo de
  ρ_DE llega a 0.29–1.4 en caminos con ρ_DE ≈ 0 sin cruzarlo (|w + 1| de 141–349 en unos pocos nodos), pero H² está bien
  ahí: era la razón de clasificar en H², no en relativo.
- **Las 7 no fieles de θ = 10 lo son por debajo de la tabla**, y se detectaron con la comprobación añadida a posteriori: con
  la del plan (solo la tabla) las 38 eran fieles y el artefacto habría entrado en los estadísticos. Tienen w(a = 10⁻³) =
  1.1–2.1 en el borde; CAMB extrapola
  w constante y ρ_DE ∝ a^{−3(1+w)} crece hacia el pasado hasta ser el 3.5 × 10⁻⁴–3 × 10⁻³ de ρ_total en a = 10⁻⁵ en cinco
  de ellas y **11 y 210 veces ρ_total** en las otras dos: energía oscura dominante antes de la recombinación, que no es
  el camino que pasó el ajuste. Sus Δχ² (36 a 65 523, todo en los picos) son el artefacto del borde de la tabla, no física
  de energía oscura tardía, y por eso van aparte. Entre las fieles, w en el borde tiene mediana −1.00 / −1.12 / −1.56 y
  máximo −0.001 / −0.09 / +0.24.

**Estadísticos, realización menos referencia (solo fieles):**

| techo | n | Δχ² total p50 / p95 / máx | Δχ²_TT p50 / p95 | Δχ²_φφ p50 / p95 | picos p50 / p95 | total > 1 | total > 4 |
|---|---|---|---|---|---|---|---|
| θ = 1 | 27 | −0.03 / 0.54 / 0.90 | 0.09 / 0.51 | −0.20 / 0.27 | 0.007 / 0.16 | 0 % | 0 % |
| θ = 3 | 40 | 0.00 / 0.78 / 1.36 | 0.20 / 0.50 | −0.24 / 0.35 | 0.016 / 0.10 | 5 % | 0 % |
| θ = 10 | 31 | −0.10 / 0.91 / 0.98 | 0.12 / 0.50 | −0.32 / 0.44 | 0.010 / 0.039 | 0 % | 0 % |

No fieles (θ = 10, 7): total mediana 231, mínimo 36, máximo 65 523; TT hasta 5.1, φφ hasta 0.28, picos de 36 a 65 525.

**Dónde cae el máximo de |Δ/σ| (fieles):**

| techo | en TT ℓ < 30 | en φφ | en picos | ℓ del máximo de TT (p16 / p50 / p84) | máx \|Δ/σ\| TT p50 / p95 | φφ p50 / p95 | picos p50 / p95 |
|---|---|---|---|---|---|---|---|
| θ = 1 | 44 % | 56 % | 0 % | 2 / 2 / 2 | 0.08 / 0.14 | 0.08 / 0.37 | 0.004 / 0.018 |
| θ = 3 | 35 % | 65 % | 0 % | 2 / 2 / 2 | 0.09 / 0.16 | 0.12 / 0.29 | 0.006 / 0.014 |
| θ = 10 | 29 % | 71 % | 0 % | 2 / 2 / 2.8 | 0.05 / 0.12 | 0.09 / 0.18 | 0.005 / 0.009 |

Cuando el máximo está en TT es siempre ℓ = 2 (ℓ < 10 en el 100 %), pero más veces está en una banda de φφ. Todos los
|Δ/σ| son ≤ 0.4: son dos señales al nivel del ruido y cuál gana es cuestión de azar.

**Artefactos.** Correlación (Spearman) del Δχ² total de las fieles con |w + 1| máx: +0.05 / −0.26 / −0.33 (θ = 1 / 3 / 10);
con la amplitud del camino (desviación típica de X_s − X_s(0) en la tabla): −0.14 / +0.13 / +0.02; con su máximo: −0.05 /
+0.19 / +0.42 (Pearson +0.40 en θ = 10). El Δχ² no sigue a |w + 1| (las correlaciones son nulas o negativas); solo en θ = 10
sigue algo al máximo del camino, débilmente. Los picos de w entran en el resultado únicamente a través del borde de la
tabla, en el grupo no fiel.

**H0.** Las realizaciones ajustadas tienen H0 = 56.7–72.0 / 56.5–77.8 / 56.7–75.3 (medianas 63.8 / 63.8 / 66.1) y su
referencia ΛCDM con el mismo θ* difiere de ellas en −4.5 a +10.9 km/s/Mpc: a estas amplitudes el residuo mueve H0 mucho, y
θ* lo compensa en la referencia. No se usa ningún prior en H0.

**Sesgos, techo por techo** (los dos grupos son disjuntos por construcción: los de ρ_DE ≤ 0 no entran en CAMB; intersección
0 en los tres techos):

| techo | no pasan tras suavizar | ρ_DE ≤ 0 | no fieles (de las corridas) | unión de los dos sesgos / supervivientes | excluidas en total / supervivientes (estimado) | fieles disponibles |
|---|---|---|---|---|---|---|
| θ = 1 | 9 / 91 | 55 | 0 / 27 | 55 / 91 = **60 %** | 64 / 91 = **70 %** | 27 de 82 que pasan |
| θ = 3 | 10 / 93 | 23 | 0 / 40 | 23 / 93 = **25 %** | 33 / 93 = **36 %** | 40 corridas de 60 corribles |
| θ = 10 | 10 / 100 | 4 | 7 / 38 (18 %) | 11 / 100 = **11 %** (conjuntos reales) | (10 + 4 + 15.8) / 100 = **30 %** (18 % de las 86 corribles) | 31 corridas de 86 corribles (más 2 segfault aparte) |

**Para θ = 1 este test no concluye, salga lo que salga:** solo 27 de los 82 caminos que pasan se pueden correr (el 67 %
tiene ρ_DE ≤ 0 en la tabla), y los 27 son, por selección, los que menos se alejan de Λ. Que sus Δχ² sean ≤ 0.9 dice que
el CMB no descarta esa tercera parte mansa, y nada de los otros dos tercios.

## Las corridas viejas no tienen el artefacto del borde

`realizaciones.py` (§3.5 y la malla de σ_X), `picos.py` y `rapido.py` usan la misma tabla (10⁻³ ≤ a ≤ 1) y la misma
extrapolación de w por debajo de a = 10⁻³. Con w(10⁻³) < 1/3, ρ_DE ∝ a^{−3(1+w)} crece hacia el pasado más despacio que la
radiación (a⁻⁴) y la energía oscura no puede ganarle antes de la recombinación: no hay artefacto. Regenerados sus caminos con
sus semillas (`w_borde_viejos.py` → `w_borde_viejos.json`, 41 casos: 28 de `realizaciones.py` con 200 realizaciones, 3 de
`picos.py` y 10 de `rapido.py` con 60), **el máximo de w en el borde a = 10⁻³ sobre todos los casos y realizaciones es
−0.705** (rapido, θ = 10, Δ = 0.0125, σ_X = 2.204 × 10⁻²); los demás extremos: −0.827 (rapido, θ = 10, Δ = 0.0125, 1.290 ×
10⁻²), −0.850 (§3.5, θ = 20, Δ = 0.05), −0.870 (rapido, θ = 10, Δ = 0.025, 2.204 × 10⁻²), −0.930 (picos, θ = 1, techo
marginalizado 4.611 × 10⁻²), −0.950 (§3.5, θ = 1, 2.633 × 10⁻²), y ≤ −0.998 en todos los casos con θ = 0.001. Todos muy
por debajo de 1/3 (y de −1/3, con lo que ρ_DE extrapolada crece más despacio que la materia). Ningún caso viejo cambia.

## Predicción

- **Supervivientes con ρ_DE ≤ 0, 30–60 % en cada techo: no se cumplió.** Salió 67 / 28 / 4.4 %: depende fuertemente de θ
  (con θ = 10 el camino suavizado casi nunca baja ρ_DE a cero; con θ = 1 lo hace dos de cada tres veces).
- **Mediana del Δχ² total < 1 en los tres techos: cumplida** (−0.03 / 0.00 / −0.10 entre las fieles). **Menos del 20 % con
  total > 4: cumplida** (0 % entre las fieles; 18 % de las corridas de θ = 10 si se contaran las 7 del artefacto del borde).
- **Picos p95 < 1 en los tres techos: cumplida entre las fieles** (0.16 / 0.10 / 0.04); no si se cuentan las 7 no fieles de
  θ = 10 (p95 = 6176).
- **El máximo en TT ℓ < 10, no en la lente: no se cumplió.** Cuando está en TT es ℓ = 2, pero más veces está en φφ
  (56 / 65 / 71 %); ambos ≤ 0.4σ.
- **Lectura, el CMB completo no mueve las cotas de EXACTO: cumplida para θ = 3 y θ = 10**, condicionada a los caminos que
  CAMB representa bien (64 % y 70 % de los supervivientes); **para θ = 1 el test no concluye** (30 % representable).

**Conclusión.** Entre los caminos que CAMB integra fielmente, el CMB completo (ISW tardío, lente, picos) añade Δχ² ≤ 1 a los
que ya pasan BAO + θ* + ω_m en los techos de θ = 3 y 10, y no los descarta; en θ = 1 dos de cada tres supervivientes tienen
ρ_DE ≤ 0 y quedan fuera del alcance de este montaje, así que las cotas de `EXACTO.md` se quedan como están, con esa
advertencia escrita en el mismo renglón.
