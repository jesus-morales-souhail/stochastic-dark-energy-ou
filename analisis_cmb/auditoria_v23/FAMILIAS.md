# FAMILIAS: ¿es el mismo el cambio permitido de la energía oscura entre hoy y z ≈ 0.5 para distintas familias de fluctuaciones?

Exploración (no toca PREPRINT ni CLAIMS). Script `familias.py`; salidas `familias.json`, `familias.txt`, `familias_chi2.npz`
(por celda: descarte, ρ_DE ≤ 0, χ²_min y Δρ(z) de cada camino en los 12 z). Tiempo: 7.6 min con 16 procesos (celdas 5.8 min,
lineal 0.4 min, regeneración en las cotas 1 min, exacta guiada 0.4 min; la guiada se repitió aparte con nodos en el máximo).

Pregunta: ¿el cambio permitido de la energía oscura entre hoy y z ≈ 0.5 es el mismo para distintas familias de fluctuaciones, aunque sus σ_X salgan distintas?

## Predicción (escrita antes de correr nada)

1. Las cotas de σ_X difieren entre familias más de un factor 1.5 en ℓ ≥ 3 (las
   familias suaves permiten σ_X mayores).
2. La raíz cuadrática media de Δρ(z) en la cota, en el z donde las familias mejor
   coinciden, coincide entre las cuatro familias y ℓ ∈ {0.3, 1, 3} dentro de un factor
   1.5 (regla del 5 %); ese z está entre 0.3 y 0.8.
3. Incluyendo ℓ = 0.1 y ℓ = 10, dentro de un factor 2.
4. Valor de esa cota de Δρ: entre 0.08 y 0.16 con la regla del 5 %; entre 0.13 y
   0.26 con la caída de 1.92 de la marginal lineal.
5. A igual ℓ ≤ 0.3, Matérn 3/2 y gaussiano difieren menos de un 25 % en σ_X.
6. Escalones: dividir el paso de la malla entre dos cambia la cota menos de un 5 %.
7. La marginal lineal reproduce la exacta guiada a menos de un 10 % en ℓ = 1, en
   las cuatro familias.

La predicción 2 no es a ciegas: el patrón ya se había visto en el OU.

Criterio fijado ahora: si en la predicción 2 o la 3 el factor pasa de 2, la nota dice
"no hay cifra universal de Δρ" y no se busca otra magnitud que coincida.

## Montaje

El de `exacto_sn.py` sin cambios: convenio B (ρ_DE/ρ_crit,0 = Ω_Λ + X(z) − X(0)), historia exacta sobre la malla de
`exacto.py` (3506 nodos, du = 2 × 10⁻³, los z de DESI como nodos exactos), BAO + θ* + ω_m + Pantheon+ (1595 filas), mismos datos,
covarianzas y tratamiento de h, ω_m y M (M explícita en el ajuste; descarte si H² ≤ 0 en el mejor ajuste ΛCDM de las 15 filas;
regla P(χ²_min ≤ χ²_obs) = 5 %, χ²_obs = 1406.327, interpolando en ln σ_X). Solo cambia el proceso X(u), u = ln(1+z),
estacionario de media cero y varianza σ_X²:

- **OU**: K = exp(−|Δu|/ℓ), θ = 1/ℓ, `exacto_sn.ou_paths` (AR(1) exacto) con las semillas de `exacto_sn.py` en las celdas
  que existen allí (ℓ ∈ {0.1, 1, 10} ↔ θ ∈ {10, 1, 0.1}; σ_X ∈ {0.02, 0.03, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7}).
- **Matérn 3/2**: K = (1 + √3 r/ℓ) exp(−√3 r/ℓ). Proceso de Markov de dos estados (X, X′) con transición exacta
  Φ(d) = e^{−λd} [[1 + λd, d], [−λ²d, 1 − λd]], λ = √3/ℓ, innovación Q = P∞ − Φ P∞ Φᵀ, P∞ = diag(1, λ²) σ_X². Exacto en malla no
  uniforme (Φ(d₁)Φ(d₂) = Φ(d₁ + d₂), comprobado a 10⁻¹⁶).
- **Gaussiano**: K = exp(−r²/(2ℓ²)). Cholesky de la matriz de correlación con el menor jitter de {10⁻¹², 10⁻¹¹, …, 10⁻⁶} que
  funciona. **Jitter usado: 10⁻¹² σ_X² en los cinco ℓ** (muy por debajo del 10⁻⁸ pedido; no se dio la condición de parada).
- **Escalones**: valor independiente N(0, σ_X²) en tramos de anchura ℓ con desfase φ ~ U(0, ℓ) por camino; K = (1 − |Δu|/ℓ)₊.
  Se sortean primero φ y los valores (independientes de la malla) y después se evalúa el camino en los nodos, así que la malla
  solo entra por la integración. El diseño de escalones con desfase uniforme al azar por camino venía del prompt; no es una
  elección de esta nota.

Celdas: ℓ ∈ {0.1, 0.3, 1, 3, 10} × σ_X ∈ {0.02, 0.03, 0.05, 0.07, 0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 1.0}, 2000 caminos por celda,
semilla por celda y familia ([20261116, familia, ℓ, σ] salvo las celdas OU compartidas). Ajuste por camino: solo la combinación
BAO + θ* + ω_m + SN. Extras: escalones en la malla de paso mitad (7011 + nodos DESI) en ℓ = 0.3 y 1 con los mismos φ y valores
(χ²_obs en esa malla 1406.331); y "escalones gaussianizados" (StepsGauss: proceso gaussiano con el núcleo triangular, por
Cholesky, jitter 10⁻¹²) por la regla del 5 % en ℓ = 0.3 y 1, para separar lo que hace el núcleo de lo que hace la no gaussianidad.

Marginal lineal: la analítica del convenio B con la respuesta completa A de `marginal_guiada.build_A` (misma malla, M como
plantilla lineal) y el núcleo de cada familia, perfil en la malla de σ extendida de `maximos_globales.py` (hasta 1.0), caída de
1.92 y 1.35 **desde el máximo global** (brentq), con la clase del máximo. Exacta guiada: `marginal_guiada.cell` con el muestreador
del prior de cada familia (para escalones, el gaussiano del núcleo triangular, que es el prior que supone la marginal lineal),
N = 5000, nodos en 0.5/1/1.5 × σ_max de la lineal y 0.6/0.8/1/1.25/1.5 × la cota lineal −1.92.

### Escalas efectivas (en unidades de ℓ; iguales en los cinco ℓ por construcción)

| familia | separación con K = e⁻¹ | escala integral ∫₀^∞ K | valores dados |
|---|---|---|---|
| OU | 1.000 ℓ | 1.000 ℓ | ℓ, ℓ |
| Matérn 3/2 | 1.239 ℓ | 1.155 ℓ | 1.24 ℓ, 1.15 ℓ |
| gaussiano | 1.414 ℓ | 1.253 ℓ | 1.41 ℓ, 1.25 ℓ |
| escalones | 0.632 ℓ | 0.500 ℓ | 0.63 ℓ, 0.5 ℓ |

Los cuatro pares coinciden con los valores dados. En z: ℓ = 0.1 son Δz ≈ 0.1 en z = 0; ℓ = 0.3, Δz ≈ 0.35; ℓ = 1, el rango
0 ≤ z ≤ 1.7; ℓ ≥ 3 cubren todo el rango de BAO y SN (u ≤ 1.2).

### Comprobaciones

- **Reproducción de EXACTO_SN**: las 24 celdas OU con semilla compartida dan χ²_min idéntico camino a camino (máx |Δχ²| = 0.0)
  y los mismos descartes.
- **Muestreadores frente a K**: varianza y correlación a los retardos 0.05, 0.1, 0.3, 1, 3 (4000 caminos, σ_X = 1): la desviación
  máxima es 0.5 errores Monte Carlo (OU, gaussiano), 0.5 (Matérn), 0.9 (escalones), 0.8 (StepsGauss). Nada que separar.
- **σ_X pequeña**: P(σ_X = 0.02) está en 0.30–0.50 en todas las familias para ℓ ≤ 3 (0.502 en escalones ℓ = 3) (0.30–0.34 en ℓ = 0.1, donde 0.02 ya
  muerde), como en el OU de EXACTO_SN. Escalones ℓ = 10: 0.64, porque el 30 % de los caminos no tiene ningún escalón en
  0 ≤ u ≤ 7 (Δρ ≡ 0, χ²_min = χ²_obs exacto, contado como ≤): es el proceso, no un fallo (891–1069 caminos por celda con
  χ² = χ²_obs). Ninguna familia da algo distinto de lo esperado con σ_X pequeña; no se dio la condición de parada.
- rms de Δρ(z) sobre el prior por Monte Carlo frente a la analítica σ_X √(2(1 − K(u_z))): coinciden al 2–3 % (2000 caminos).

## (a) Cotas de σ_X

| familia | ℓ | 5 % (bootstrap 16–84) | lineal −1.92 | lineal −1.35 | máximo lineal |
|---|---|---|---|---|---|
| OU | 0.1 | **0.070** (0.068–0.073) | 0.118 | 0.100 | 0.032, interior |
| OU | 0.3 | **0.111** (0.108–0.115) | 0.167 | 0.145 | 0.058, interior |
| OU | 1 | **0.169** (0.164–0.173) | 0.277 | 0.242 | 0.108, interior |
| OU | 3 | **0.276** (0.270–0.283) | 0.468 | 0.408 | 0.185, interior |
| OU | 10 | **0.476** (0.466–0.488) | 0.848 | 0.740 | 0.316, interior |
| Matérn 3/2 | 0.1 | **0.079** (0.076–0.081) | 0.126 | 0.108 | 0.040, interior |
| Matérn 3/2 | 0.3 | **0.147** (0.142–0.153) | 0.248 | 0.214 | 0.086, interior |
| Matérn 3/2 | 1 | **0.394** (0.379–0.409) | 0.737 | 0.607 | 0.199, interior |
| Matérn 3/2 | 3 | sin cota | sin cota | sin cota | 0.681, interior |
| Matérn 3/2 | 10 | sin cota | sin cota | sin cota | σ → 0 |
| gaussiano | 0.1 | **0.086** (0.084–0.088) | 0.154 | 0.133 | 0.050, interior |
| gaussiano | 0.3 | **0.199** (0.194–0.211) | 0.335 | 0.269 | 0.079, interior |
| gaussiano | 1 | **0.828** (0.800–0.856) | sin cota | sin cota | 0.316, interior |
| gaussiano | 3 | sin cota | sin cota | sin cota | σ → 0 |
| gaussiano | 10 | sin cota | sin cota | sin cota | σ → 0 |
| escalones | 0.1 | **0.075** (0.072–0.078) | 0.097 | 0.081 | 0.020, interior |
| escalones | 0.3 | **0.160** (0.151–0.167) | 0.155 | 0.134 | 0.054, interior |
| escalones | 1 | sin cota (P = 0.064 en σ_X = 1) | 0.268 | 0.234 | 0.100, interior |
| escalones | 3 | sin cota (P = 0.34 en 1) | 0.463 | 0.404 | 0.185, interior |
| escalones | 10 | sin cota (P = 0.60 en 1) | 0.845 | 0.738 | 0.316, interior |
| StepsGauss | 0.3 | 0.093 (0.091–0.096) | = escalones | = escalones | — |
| StepsGauss | 1 | 0.163 (0.156–0.169) | = escalones | = escalones | — |

"Sin cota" = P > 5 % en σ_X = 1.0. Ningún máximo lineal está en el borde superior (1.0): todos son interiores o σ → 0 (en los
dos σ → 0 la lineal no se separa nunca 1.92 del valor en 0 hasta σ_X = 1). P por σ y fracciones de H² ≤ 0 en `familias.txt`.

Lo que dicen las cotas:

- Para ℓ ≤ 0.3 las cuatro familias acotan parecido (0.070–0.086 en ℓ = 0.1; 0.111–0.199 en ℓ = 0.3). A ℓ = 1 se abren:
  0.169 (OU), 0.394 (Matérn), 0.828 (gaussiano), sin cota (escalones). A ℓ ≥ 3 solo el OU tiene cota.
- Las familias suaves (Matérn, gaussiano) permiten σ_X mayores a igual ℓ porque su camino en el rango de los datos (u ≤ 1.2) es
  casi una deriva lineal, que el fondo libre absorbe en gran parte; el OU, con su "ruido" a todas las escalas, se delata antes.
- Los escalones no acotan en ℓ ≥ 1 por otra razón: el primer escalón está en u = ℓ − φ, es decir, en z = e^{ℓ−φ} − 1, con φ
  uniforme en (0, ℓ); con ℓ = 1 un 8 % de los caminos no cambia hasta z ≈ 1.5 (P = 0.064 en σ_X = 1), con ℓ = 3 el 60 % no cambia
  hasta z > 2.3 y el 33 % hasta z > 6.4 (P = 0.34), y con ℓ = 10 el 88 % no cambia hasta z > 2.3 y el 30 % no cambia nunca en
  0 ≤ u ≤ 7 (P = 0.60). Un cambio ρ_DE ± ρ_crit,0 en z > 2.3 no lo ve ningún dato, ni θ* (con h y ω_m libres): esos caminos dan
  P = 0.50 ± 0.02 en todo σ_X ≤ 1, como ΛCDM (comprobación numérica en la sección "Suelo de P en los escalones"). Con el mismo núcleo pero gaussiano (StepsGauss) la cota vuelve: 0.093 (ℓ = 0.3) y 0.163 (ℓ = 1), parecida al OU
  (0.111, 0.169), que tiene casi el mismo núcleo a retardos cortos. La no gaussianidad del prior de escalones afloja la cota ×1.7
  en ℓ = 0.3 y la quita en ℓ = 1.
- **Aviso sobre EXACTO_SN.md**: con σ_X = 0.07 y 0.15 en la malla, las celdas OU (idénticas) dan la cota 5 % de
  BAO + θ* + ω_m + SN en θ = 10 como **0.070**, no 0.081 (−13 %), y en θ = 1 como 0.169, no 0.173 (−2 %); θ = 0.1 queda en 0.476.
  La diferencia es solo la interpolación en ln σ_X entre 0.05 (P = 0.107) y 0.1 (0.024), que con 0.07 (P = 0.049) en medio se
  corta antes. CONVERGENCIA_SN.md ya lo tenía con las celdas 0.07 y 0.15 y N = 20000: 0.170 / 0.102 / 0.073. No se corrige nada;
  las cifras 0.173 / 0.113 / 0.081 de EXACTO_SN.md y ERRATAS_V3.md llevan ahora una nota de malla junto a cada cita
  (BARRIDO_MAXIMOS.md no las cita). θ = 3 no se ha repetido aquí (ℓ = 0.3 es θ = 3.33).

### Suelo de P en los escalones (comprobación numérica, `familias_suelo.py`, `familias_suelo_escalones.json`)

Sin ajustes nuevos: los χ² son los de `familias_chi2.npz` y los desfases φ se regeneran de las semillas. Por ℓ, fracción del
prior cuyo primer tramo cubre todos los datos de BAO y SN (primer escalón en u₁ = ℓ − φ > 1.203, z > 2.33; las SN llegan a
z = 2.26), fracción sin ningún escalón en 0 ≤ u ≤ 7 (cubre también θ*), y el suelo de P cuando σ_X crece:

| ℓ | cubre BAO + SN | sin escalón en el rango | P(σ_X = 0.5 / 0.7 / 1.0) | P dado que cubre | P dado primer escalón en z ≤ 2.33 | supervivientes que cubren (σ_X = 1) |
|---|---|---|---|---|---|---|
| 0.1 | 0 | 0 | 0 / 0 / 0 | — | 0 | — |
| 0.3 | 0 | 0 | 0.005 / 0.003 / 0 | — | 0.005 / 0.003 / 0 | — |
| 1 | 0 | 0 | 0.122 / 0.080 / 0.064 | — | 0.122 / 0.080 / 0.064 | 0 |
| 3 | 0.60–0.62 | 0 | 0.367 / 0.373 / 0.339 | 0.50 / 0.53 / 0.50 | 0.17 / 0.12 / 0.08 | 0.82 / 0.88 / 0.92 |
| 10 | 0.88–0.90 | 0.28–0.31 | 0.626 / 0.605 / 0.602 | 0.68 / 0.67 / 0.67 (1.00 los sin escalón) | 0.17 / 0.13 / 0.08 | 0.97 / 0.98 / 0.98 |

- **ℓ = 3 y 10: coinciden.** Los caminos cuyo primer tramo cubre BAO y SN dan P = 0.50 ± 0.02 en todo σ_X (como ΛCDM: χ²_min
  queda por encima o por debajo de χ²_obs a cara o cruz; el cambio en z > 2.3 no lo ve ni θ*), y los que no tienen ningún escalón
  dan χ²_min = χ²_obs exacto, contado como ≤. El suelo es por tanto ½ × (masa que cubre) + ½ × (masa sin escalón) + el resto:
  ℓ = 3: 0.5 × 0.60 + 0.08 × 0.40 = 0.33 (medido 0.34); ℓ = 10: 0.5 × (0.88 − 0.31) + 1.0 × 0.31 + 0.08 × 0.12 = 0.60 (medido
  0.60). El 92–98 % de los supervivientes en σ_X = 1 son caminos que cubren. **El "sin cota" de los escalones en ℓ ≥ 3 lo
  produce el prior (masa sin cambio sobre los datos del 60 % y el 88 %, mayor que el 5 %), no los datos.**
- **ℓ = 1: no coinciden.** Ningún camino cubre los datos (u₁ ≤ 1 < 1.2) y el suelo es 0.064. Por z del primer escalón, en
  σ_X = 1: z₁ ≤ 0.5 (40 % del prior) P = 0.014; 0.5–1 (29 %) 0.037; 1–1.5 (22 %) 0.126; 1.5–2.33 (9 %) 0.177. El 75 % de los
  supervivientes tiene el primer escalón en z > 1. Aquí el suelo lo sostienen los datos de z > 1 (QSO, Lyα y pocas SN), que
  ven el cambio pero flojo, no una masa del prior sin cambio; con σ_X > 1 seguiría bajando (0.122 → 0.080 → 0.064 de 0.5 a 1).
- ℓ ≤ 0.3: todo camino cambia dentro de los datos y P → 0.

## (b) y (c) Δρ(z) = X(z) − X(0) en la cota

rms sobre el prior en la cota: σ_b √(2(1 − K(u_z))). El z_* se busca entre los 12 z como el que minimiza el cociente máx/mín de
esa rms entre las familias × ℓ que tienen cota (regla del 5 %, ℓ ∈ {0.3, 1, 3}); las entradas "sin cota" no entran en el cociente.

| regla | conjunto | n con cota (faltan) | z_* | máx/mín en z_* | rms en z_* | máx/mín en z = 0.5 | rms en z = 0.5 |
|---|---|---|---|---|---|---|---|
| 5 % | ℓ ∈ {0.3, 1, 3} | 8 de 12 (Matérn 3, gauss 3, escalones 1 y 3) | **0.2** | **1.86** | 0.095–0.177 (mediana 0.111) | **2.43** | 0.136–0.329 |
| 5 % | todos los ℓ | 13 de 20 (+ Matérn 10, gauss 10, escalones 10) | 0.2 | **1.96** | 0.090–0.177 (mediana 0.107) | **3.37** | 0.098–0.329 |
| −1.92 | ℓ ∈ {0.3, 1, 3} | 9 de 12 (Matérn 3, gauss 1 y 3) | 0.1 | 1.18 | 0.105–0.124 (mediana 0.117) | 2.02 | 0.204–0.412 |
| −1.92 | todos los ℓ | 15 de 20 (+ Matérn 10, gauss 10) | 0.1 | 1.27 | 0.105–0.134 (mediana 0.117) | 3.01 | 0.137–0.412 |
| −1.35 | ℓ ∈ {0.3, 1, 3} | 9 de 12 | 0.2 | 1.25 | 0.139–0.173 (mediana 0.141) | 1.92 | 0.177–0.340 |
| −1.35 | todos los ℓ | 15 de 20 | 0.1 | 1.34 | 0.084–0.113 (mediana 0.102) | 2.98 | 0.114–0.340 |

El cociente crece monótonamente con z en todas las reglas (regla 5 %, conjunto de la predicción 2: 2.08, 1.86, 1.91, 2.16, 2.43,
2.70, 2.96, 3.19, 3.41, 3.62, 4.46, 5.38 en z = 0.1 … 2.33). Quien abre el abanico son las familias suaves a ℓ = 1 (rms(0.5) =
0.22 Matérn, 0.33 gaussiano frente a 0.14 del OU) y, en z_*, los escalones a ℓ = 0.3 (0.177; sin escalones el conjunto de la
predicción 2 da 1.34 en z_* = 0.1, pero 2.43 en z = 0.5 de todas formas).

Tabla de la rms del prior en z_* = 0.2 (regla 5 % / −1.92 / −1.35), p95 de |Δρ(0.2)| y de |Δρ(0.5)| entre los supervivientes en la
cota 5 % (caminos regenerados en σ_b con la semilla de la celda más próxima; 93–103 supervivientes), y (d) fracción de
supervivientes con ρ_DE ≤ 0 en algún nodo:

| familia | ℓ | rms(0.2) 5 % | −1.92 | −1.35 | p95 \|Δρ(0.2)\| surv | p95 \|Δρ(0.5)\| surv | ρ_DE ≤ 0 surv (5 %) | ρ_DE ≤ 0 surv (−1.92) |
|---|---|---|---|---|---|---|---|---|
| OU | 0.1 | 0.090 | 0.153 | 0.130 | 0.109 | 0.185 | 0.00 | 0.00 |
| OU | 0.3 | 0.106 | 0.160 | 0.139 | 0.089 | 0.198 | 0.00 | 0.06 |
| OU | 1 | 0.097 | 0.160 | 0.140 | 0.089 | 0.182 | 0.06 | 0.12 |
| OU | 3 | 0.095 | 0.161 | 0.140 | 0.103 | 0.176 | 0.08 | 0.41 |
| OU | 10 | 0.091 | 0.161 | 0.141 | 0.095 | 0.183 | 0.19 | 0.50 |
| Matérn | 0.1 | 0.101 | 0.162 | 0.139 | 0.083 | 0.174 | 0.00 | 0.00 |
| Matérn | 0.3 | 0.110 | 0.186 | 0.161 | 0.073 | 0.179 | 0.00 | 0.24 |
| Matérn | 1 | 0.112 | 0.210 | 0.173 | 0.070 | 0.149 | 0.44 | 0.67 |
| gaussiano | 0.1 | 0.110 | 0.197 | 0.169 | 0.072 | 0.181 | 0.00 | 0.00 |
| gaussiano | 0.3 | 0.116 | 0.194 | 0.156 | 0.062 | 0.137 | 0.02 | 0.55 |
| gaussiano | 1 | 0.150 | — | — | 0.069 | 0.110 | 0.68 | — |
| escalones | 0.1 | 0.107 | 0.137 | 0.114 | 0.097 | 0.177 | 0.00 | 0.00 |
| escalones | 0.3 | 0.177 | 0.171 | 0.148 | 0.105 | 0.231 | 0.01 | 0.00 |
| escalones | 1, 3, 10 | — | 0.162, 0.161, 0.161 | 0.141 | — | — | — | 0.17, 0.38, 0.25 |

(Para las entradas sin cota, la rms que tendrían en σ_X = 1: Matérn ℓ = 3: 0.10 (z = 0.2) / 0.22 (z = 0.5); ℓ = 10: 0.03 / 0.07;
gaussiano ℓ = 3: 0.06 / 0.14; ℓ = 10: 0.02 / 0.04; escalones ℓ = 1: 0.60 / 0.90; ℓ = 3: 0.35 / 0.52; ℓ = 10: 0.19 / 0.29. Las
suaves de ℓ ≥ 3 están por debajo de la banda en z = 0.2 incluso en σ_X = 1, los escalones muy por encima.)

Lectura: la rms del prior en z = 0.5 no es la misma entre familias (×2.4–3.4). Entre los supervivientes la dispersión es menor
(p95 |Δρ(0.5)| 0.11–0.23, ×2.1; rms de supervivientes 0.064–0.117, ×1.8) y en z = 0.2 0.062–0.109 (×1.75): el dato restringe
el cambio real de los caminos que sobreviven más parecido de lo que la rms del prior en la cota sugiere, pero tampoco dentro de
×1.5. Las cotas grandes de las familias suaves a ℓ = 1 son cotas sobre caminos de los que el 44–68 % tiene ρ_DE ≤ 0 en algún z
(a z alto, donde la materia lo tapa), y a la caída −1.92 el 41–67 % en Matérn ℓ = 1, gaussiano ℓ = 0.3 y OU ℓ ≥ 3.

## (e) Exacta guiada en ℓ = 1 (N = 5000)

| familia | lineal −1.92 (analítica) | exacta −1.92 (bootstrap 16–84) | exacta, solo nodos con ESS ≥ 1 % N | −1.35 lineal / exacta | ln(Z_ex/Z_lin) en la cota lineal | ESS en la cota (mínimo entre nodos) |
|---|---|---|---|---|---|---|
| OU | 0.277 | 0.299 (0.290–0.309) | 0.289 | 0.242 / 0.250 | +0.23 | 134 (12, en 1.25 × la cota) |
| Matérn 3/2 | 0.737 | 0.750 (0.747–0.753) | 0.750 | 0.607 / 0.614 | +0.05 | 2656 (947) |
| gaussiano | sin cota (máximo en 0.316, caída 1.1 en σ_X = 1) | no se corre: no hay cota lineal que comparar | — | — | — | — |
| escalones (prior gaussiano del núcleo triangular) | 0.268 | ESS = 1 en 0.21–0.34 (un solo peso; corrección +38) | sin cota (quedan 2 nodos útiles) | 0.234 / — | +38 (un peso) | 1 (1) |

- OU: exacta +8 % sobre la lineal (+4 % con los nodos fiables); MARGINAL_GUIADA.md daba +5.6 % (0.287 frente a 0.272) con otros
  nodos y semilla. Matérn: +2 %. En los dos la lineal reproduce la exacta dentro del 10 %.
- Escalones: el muestreo guiado se desploma (ESS 1–39 en σ_X ≥ 0.1; 1 en 0.21–0.34). Diagnóstico: los pesos grandes son muestras con
  δh ≈ −0.034 y Δρ ≈ −0.9 en z > 3, donde la predicción lineal de θ* se desvía 2.4σ_θ de la exacta y χ²_exacto queda 8–12 por
  debajo de χ²_lineal; la misma población existe en el OU (ln w hasta 6 en 3000 muestras) pero allí es cola, no un peso
  dominante. No se fuerza nada: para esta familia la comparación lineal/exacta queda sin cifra.

## (f) Escalones en la malla de paso mitad

Mismos φ y valores, 7011 nodos + DESI, χ²_obs = 1406.331 (frente a 1406.327): ℓ = 0.3: 0.1602 → 0.1604 (+0.1 %); ℓ = 1: sin cota en
las dos mallas, P por σ igual a ±0.0005. La discretización no cuenta.

## Predicción, punto por punto

1. **Cumplida.** En ℓ = 3 y 10 solo el OU tiene cota (0.276, 0.476); Matérn, gaussiano y escalones no la tienen por debajo de
   1.0 (factor > 3.6 y > 2.1). Las suaves permiten σ_X mayores, sí; los escalones también, pero porque en una fracción de los
   caminos el primer cambio cae más allá de los datos, no por suavidad. El abanico se abre ya en ℓ = 1 (0.169 / 0.394 / 0.828 / sin cota).
2. **No cumplida.** z_* = 0.2 (fuera de 0.3–0.8) y el factor es 1.86 entre los 8 de 12 que tienen cota (en z = 0.5, 2.43). No
   pasa de 2 en z_*, así que el criterio fijado no salta por la letra; pero pasa de 2 en el z de la pregunta (z ≈ 0.5), y 4 de
   los 12 no tienen cota que comparar. Conclusión: **no hay cifra universal de Δρ entre hoy y z ≈ 0.5**; en z ≈ 0.2 la hay dentro
   de ×1.9 entre las familias acotadas. No se busca otra magnitud.
3. **Cumplida por la letra, no en sustancia.** Con todos los ℓ el factor en z_* = 0.2 es 1.96 (< 2) entre 13 de 20; 7 no tienen
   cota y en z = 0.5 el factor es 3.37.
4. **En parte.** Regla 5 %: 0.090–0.177 en z_* = 0.2 (mediana 0.107): dentro de 0.08–0.16 salvo los escalones ℓ = 0.3 (0.177).
   Caída −1.92: 0.105–0.134 en z_* = 0.1 (mediana 0.117), por debajo del 0.13–0.26 previsto; en z = 0.2, 0.137–0.210.
5. **En parte.** ℓ = 0.1: Matérn 0.079 frente a gaussiano 0.086 (+9 %, cumplida); ℓ = 0.3: 0.147 frente a 0.199 (+35 %, no).
   Lo mismo con la lineal (+22 % y +35 %).
6. **Cumplida.** +0.1 % en ℓ = 0.3; sin cambio en ℓ = 1.
7. **Cumplida en dos familias, indecidible en dos.** OU +8 % (+4 % con nodos fiables), Matérn +2 %; gaussiano sin cota lineal
   que comparar; escalones con ESS = 1.

## Común a todas las familias

- Para ℓ ≤ 0.3 (cambios en Δz ≲ 0.35) la cota de σ_X está en 0.07–0.09 (ℓ = 0.1) y 0.11–0.20 (ℓ = 0.3) en las cuatro familias,
  y la rms de Δρ en z = 0.2 en 0.09–0.18.
- Con la marginal lineal (−1.92) la rms de Δρ(0.1) en la cota es 0.105–0.134 en las 15 entradas con cota, y Δρ(0.2) en
  0.137–0.210: a z bajo la lineal sí da casi una cifra universal, 0.12 ± 0.01 en z = 0.1.
- Entre los supervivientes de la cota 5 %, p95 |Δρ(0.2)| = 0.06–0.11 y p95 |Δρ(0.5)| = 0.11–0.23 en las 13 entradas con cota.
- La exacta no lineal (5 %) es más estrecha que la lineal −1.92 por ×1.5–1.9 en OU, Matérn y gaussiano en todos los ℓ con cota;
  en los escalones ×1.3 (ℓ = 0.1) y ×1.0 (ℓ = 0.3). La lineal −1.35 queda en medio.
- Ningún máximo lineal está en el borde; la discretización no cuenta (escalones, malla mitad); los muestreadores reproducen K.

## Depende de la familia

- La cota de σ_X a ℓ ≥ 1: ×2.3 entre OU y Matérn en ℓ = 1, ×4.9 entre OU y gaussiano; a ℓ ≥ 3 solo el OU acota por debajo de 1.
- La rms de Δρ en z = 0.5 en la cota: ×2.4 (ℓ ∈ {0.3, 1, 3}) y ×3.4 (todos) entre familias; el cociente crece con z.
- Qué caminos sobreviven: en las suaves a ℓ = 1, la mayoría de los supervivientes en la cota tiene ρ_DE ≤ 0 en algún z
  (44 % Matérn, 68 % gaussiano); en el OU, 6 %.
- Los escalones: su prior no gaussiano afloja la cota ×1.7 en ℓ = 0.3 respecto al proceso gaussiano del mismo núcleo y la quita
  en ℓ ≥ 1: en ℓ ≥ 3 porque el 60–88 % del prior no cambia sobre los datos (suelo del prior); en ℓ = 1 porque el 31 % del prior
  cambia solo en z > 1, donde los datos lo ven flojo (suelo de los datos).
- El acuerdo lineal/exacta: +2 % (Matérn) a +8 % (OU) donde se puede medir.

## No se puede decir con estos datos

- Si existe una cifra universal de Δρ a z ≈ 0.5: no la hay entre las familias acotadas, y 7 de las 20 (familia × ℓ) no tienen
  cota por debajo de σ_X = 1, con rms(0.5) en σ_X = 1 de 0.04–0.22 (suaves) y 0.29–0.90 (escalones).
- La cota de Matérn y gaussiano a ℓ ≥ 3 (más allá de σ_X = 1, donde ρ_DE ≤ 0 es la norma); la de los escalones a ℓ ≥ 3 no
  existe con estos datos (suelo del prior) y a ℓ = 1 está por encima de 1.
- Si la marginal lineal reproduce la exacta en el gaussiano (sin cota lineal) y en los escalones (ESS = 1 con la guía lineal;
  haría falta otra propuesta, y eso ya no es `marginal_guiada.py`).
- La cota exacta de los escalones no gaussianos por marginal: la guiada supone prior gaussiano.
- El aviso sobre EXACTO_SN (θ = 10: 0.081 → 0.070 por la malla de σ) está fuera de esta pregunta; las cifras históricas no se
  han tocado, solo se han anotado.
