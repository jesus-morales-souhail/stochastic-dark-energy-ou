# MUESTRAS_SN: ¿sobreviven las cotas por familia al cambiar de muestra de supernovas?

Exploración (no toca PREPRINT ni CLAIMS). Script `muestras_sn.py`; salidas `muestras_sn.{json,txt}`.

Pregunta principal: ¿sobreviven las cotas por familia al cambiar de muestra de supernovas? Pregunta secundaria (prerregistrada,
era una observación a posteriori en FAMILIAS): ¿se mantiene la coincidencia de Δρ(0.1)?

La coincidencia de Δρ(0.1) se vio en FAMILIAS después de mirar los datos; aquí se
pone a prueba por primera vez con muestras no usadas.

Reglas fijadas: si una familia deja de tener una cota estable con otra muestra, no se
salva cambiando de estadístico ni de criterio después de ver los resultados. Las
muestras no se combinan entre sí (comparten supernovas). Una mejora del ajuste con
σ_X > 0 se lee como tensión entre datos absorbida, no como detección.

## Predicciones (escritas antes de descargar ni calcular nada)

1. Con DES-SN5YR y con Union3 las cotas (5 % y lineal −1.92) son más flojas que con
   Pantheon+, entre ×1.1 y ×1.6, en todas las celdas con cota.
2. Dovekie queda entre Pantheon+ y DES-SN5YR en cada celda.
3. Con DES-SN5YR el máximo de ln L está en σ_X > 0 con Δln L entre +1 y +4.
4. El cociente Matérn/OU y gaussiano/OU en ℓ = 1 cambia menos de un ±30 % entre
   muestras.
5. Secundaria: con la lineal −1.92, en cada muestra la dispersión de Δρ(0.1) en la
   cota entre familias es menor que un factor 1.3, y su valor está entre 0.08 y 0.18.
6. Ninguna celda de OU con θ ≥ 1 pierde su cota.

## Datos (entrada no fiable; cada muestra en `datos/<muestra>/` con `FUENTE.json`: URL oficial, fecha, sha256 y convenio; se leen solo como números, script corrido con `python -I`)

| muestra | fuente oficial | ficheros | n | convenio |
|---|---|---|---|---|
| Pantheon+ | `analisis_cmb/sn/datos` (referencia, Brout et al. 2022) | Pantheon+SH0ES.dat, STAT+SYS | 1580 (z_HD > 0.01, sin calibradores) | como `exacto_sn.py` |
| DES-SN5YR | github.com/des-science/DES-SN5YR, etiqueta 1.3, `4_DISTANCES_COVMAT` (Abbott et al. 2024, arXiv:2401.02929) | DES-SN5YR_HD.csv, STAT+SYS.txt.gz | 1829 (1635 DES + 194 de z bajo) | MU (H0 = 70); C = diag(MUERR_FINAL²) + sistemáticos (STATONLY es cero por construcción, README) |
| Dovekie | mismo repositorio, rama main (Popovic et al. 2026, arXiv:2511.07517) | DES-Dovekie_HD.csv (formato SNANA), STAT+SYS.npz | 1820 (1623 DES + 197 de z bajo) | MU (H0 = 70); el fichero es la covarianza **inversa** total (stat + sys), triángulo superior en float32; C = su inversa, como hace la verosimilitud oficial |
| Union3 | github.com/rubind/union3_release (Rubin et al. 2023, arXiv:2311.12098) | mu_mat_union3_cosmo=2_mu.fits | 22 nodos en z = 0.05 … 2.26 | FITS 23 × 23: fila 0 = z, columna 0 = μ, resto = covarianza inversa 22 × 22; se usan en tramos tal cual, con z_hel = z |

Descargas el 2026-11-12 (sha256 en `datos/FUENTES.json`). Dovekie sí tiene distancias y covarianza públicas (las del repositorio de DES en
su rama main, que ha sustituido a las de 2024; las de 2024 siguen en la etiqueta 1.3). En todas, el desplazamiento absoluto (M) se
marginaliza: parámetro libre en el ajuste no lineal (equivalente a la marginalización analítica) y plantilla lineal en la marginal.
Supernovas compartidas (por CID de los ficheros oficiales): Pantheon+ (usadas) y DES-SN5YR, **335** (190 de z bajo por nombre y 145 de
DES-3YR por número); Pantheon+ y Dovekie, 339; DES-SN5YR y Dovekie, 1718. Union3 no publica identificadores. Las muestras no se combinan.

## Comprobaciones antes de seguir

(b) ΛCDM plano con solo las supernovas de cada muestra (Ω_m libre, h fijo para Ω_r, M marginalizada), frente al valor publicado:

| muestra | n | Ω_m aquí | publicado | diferencia | de dónde se copió |
|---|---|---|---|---|---|
| Pantheon+ | 1580 | 0.3324 | 0.334 ± 0.018 | -0.0016 | abstract de arXiv:2202.04077, flat LCDM, SN solas |
| DES-SN5YR | 1829 | 0.3518 | 0.352 ± 0.017 | -0.0002 | abstract de arXiv:2401.02929 (Abbott et al.), flat LCDM, SN solas |
| Dovekie | 1820 | 0.3302 | 0.33 ± 0.015 | +0.0002 | abstract de arXiv:2511.07517 (Popovic et al.), flat LCDM |
| Union3 | 22 | 0.3558 | 0.356 ± +0.028/-0.026 | -0.0002 | Tabla 5 de arXiv:2311.12098 (Rubin et al.), fila 'SNe', flat LCDM |

Todas a menos de 0.002 (Union3 reproduce además el χ² de su Tabla 5: 23.96 frente a 24.0 con 20 grados de libertad). No saltó la parada.

(a) Con Pantheon+ las celdas de FAMILIAS se reproducen exactamente (mismas semillas, mismos caminos): OU|theta_1: Δ = 0; OU|theta_10: Δ = 0; Matern32|l_0.3: Δ = 0; Matern32|l_1: Δ = 0; Gauss|l_0.3: Δ = 0; Gauss|l_1: Δ = 0; OU|theta_0.3: máx |ΔP| = 0.0; OU|theta_3: máx |ΔP| = 0.0. Las celdas OU
con θ = 0.3 y 3 comparten semilla con EXACTO_SN y dan la misma P.

## Montaje

El de `familias.py` sin cambios (convenio B, historia exacta, BAO + θ* + ω_m + SN, descarte por H² ≤ 0 en el mejor ajuste ΛCDM de las
15 filas, regla P(χ²_min ≤ χ²_obs) = 5 % interpolando en ln σ_X, bootstrap de 1000), cambiando solo la muestra de SN. **Los 2000 caminos
de cada celda son los mismos para las cuatro muestras** (semilla por celda como en `familias.py` y `exacto_sn.py`): la comparación entre
muestras es emparejada. Celdas: OU θ ∈ {0.3, 1, 3, 10} (ℓ = 1/θ), Matérn 3/2 y gaussiano ℓ ∈ {0.3, 1}; σ_X ∈ {0.03, 0.05, 0.07, 0.1, 0.15,
0.2, 0.3, 0.5, 0.7, 1.0}. χ²_obs (ΛCDM, X = 0): Pantheon+ 1406.33, DES-SN5YR 1666.05, Dovekie 1651.88, Union3 45.34. Lineal: la
marginal analítica del convenio B con la respuesta A rehecha para cada muestra (lógica de `marginal_guiada.build_A` con la muestra como
argumento; filas de A suman cero a 10⁻¹⁵), perfil en la malla extendida hasta σ_X = 1 y caída de 1.92 / 1.35 desde el máximo global;
todos los máximos resultaron interiores. Δρ(0.1) en la cota: σ_b √(2(1 − K(ln 1.1))). Tiempo: 446 s con 16 procesos.

## Tabla final

Cota del 5 % (bootstrap 16–84), lineal −1.92 y −1.35, máximo de ln L relativo a σ_X = 0 (y su σ_X), cociente frente a Pantheon+ (5 % y
−1.92) y Δρ(0.1) en la cota (5 % y −1.92):

| celda | muestra | 5 % | −1.92 | −1.35 | máx ln L (σ_X) | /P+ 5 % | /P+ −1.92 | Δρ(0.1) 5 % | Δρ(0.1) −1.92 |
|---|---|---|---|---|---|---|---|---|---|
| OU theta_0.3 | Pantheon+ | 0.307 (0.294–0.325) | 0.492 | 0.430 | +1.17 (0.185) | — | — | 0.073 | 0.117 |
| OU theta_0.3 | DES-SN5YR | 0.470 (0.461–0.478) | 0.516 | 0.448 | +4.86 (0.199) | 1.53 | 1.05 | 0.112 | 0.122 |
| OU theta_0.3 | Dovekie | 0.371 (0.358–0.384) | 0.491 | 0.423 | +1.84 (0.158) | 1.21 | 1.00 | 0.088 | 0.116 |
| OU theta_0.3 | Union3 | 0.500 (0.489–0.524) | 0.687 | 0.602 | +3.59 (0.293) | 1.63 | 1.40 | 0.119 | 0.163 |
| OU theta_1 | Pantheon+ | 0.169 (0.164–0.173) | 0.277 | 0.242 | +1.14 (0.108) | — | — | 0.072 | 0.118 |
| OU theta_1 | DES-SN5YR | 0.270 (0.265–0.276) | 0.290 | 0.251 | +4.86 (0.108) | 1.60 | 1.04 | 0.115 | 0.123 |
| OU theta_1 | Dovekie | 0.189 (0.186–0.193) | 0.275 | 0.237 | +1.83 (0.0856) | 1.12 | 0.99 | 0.081 | 0.117 |
| OU theta_1 | Union3 | 0.274 (0.269–0.279) | 0.387 | 0.339 | +3.61 (0.158) | 1.62 | 1.39 | 0.117 | 0.165 |
| OU theta_3 | Pantheon+ | 0.109 (0.105–0.114) | 0.174 | 0.151 | +0.98 (0.0583) | — | — | 0.077 | 0.123 |
| OU theta_3 | DES-SN5YR | 0.176 (0.172–0.180) | 0.182 | 0.157 | +4.76 (0.068) | 1.61 | 1.04 | 0.124 | 0.128 |
| OU theta_3 | Dovekie | 0.131 (0.129–0.134) | 0.171 | 0.147 | +1.77 (0.05) | 1.20 | 0.98 | 0.093 | 0.121 |
| OU theta_3 | Union3 | 0.183 (0.180–0.187) | 0.245 | 0.216 | +3.55 (0.0998) | 1.67 | 1.41 | 0.129 | 0.173 |
| OU theta_10 | Pantheon+ | 0.070 (0.068–0.074) | 0.118 | 0.100 | +0.59 (0.0316) | — | — | 0.077 | 0.131 |
| OU theta_10 | DES-SN5YR | 0.134 (0.132–0.137) | 0.125 | 0.108 | +4.59 (0.05) | 1.93 | 1.06 | 0.149 | 0.138 |
| OU theta_10 | Dovekie | 0.098 (0.095–0.100) | 0.114 | 0.096 | +1.74 (0.0316) | 1.40 | 0.96 | 0.108 | 0.126 |
| OU theta_10 | Union3 | 0.136 (0.133–0.138) | 0.172 | 0.152 | +3.32 (0.0734) | 1.95 | 1.46 | 0.150 | 0.191 |
| Matern32 l_0.3 | Pantheon+ | 0.147 (0.142–0.153) | 0.248 | 0.214 | +1.40 (0.0856) | — | — | 0.067 | 0.114 |
| Matern32 l_0.3 | DES-SN5YR | 0.210 (0.199–0.220) | 0.254 | 0.217 | +4.76 (0.0793) | 1.43 | 1.03 | 0.096 | 0.117 |
| Matern32 l_0.3 | Dovekie | 0.163 (0.155–0.169) | 0.242 | 0.205 | +1.70 (0.0583) | 1.11 | 0.98 | 0.075 | 0.111 |
| Matern32 l_0.3 | Union3 | 0.240 (0.232–0.247) | 0.343 | 0.290 | +3.40 (0.116) | 1.63 | 1.39 | 0.110 | 0.158 |
| Matern32 l_1 | Pantheon+ | 0.394 (0.378–0.408) | 0.737 | 0.607 | +1.50 (0.199) | — | — | 0.062 | 0.115 |
| Matern32 l_1 | DES-SN5YR | 0.599 (0.583–0.612) | 0.695 | 0.566 | +5.01 (0.199) | 1.52 | 0.94 | 0.094 | 0.109 |
| Matern32 l_1 | Dovekie | 0.439 (0.425–0.451) | 0.645 | 0.516 | +1.95 (0.158) | 1.11 | 0.87 | 0.069 | 0.101 |
| Matern32 l_1 | Union3 | 0.643 (0.629–0.656) | 0.828 | 0.662 | +3.41 (0.232) | 1.63 | 1.12 | 0.100 | 0.129 |
| Gauss l_0.3 | Pantheon+ | 0.199 (0.193–0.211) | 0.335 | 0.269 | +1.38 (0.0793) | — | — | 0.063 | 0.105 |
| Gauss l_0.3 | DES-SN5YR | 0.296 (0.289–0.311) | 0.320 | 0.252 | +4.82 (0.0793) | 1.49 | 0.96 | 0.093 | 0.100 |
| Gauss l_0.3 | Dovekie | 0.203 (0.196–0.214) | 0.284 | 0.220 | +1.78 (0.0583) | 1.02 | 0.85 | 0.064 | 0.089 |
| Gauss l_0.3 | Union3 | 0.342 (0.327–0.355) | 0.323 | 0.263 | +3.30 (0.0998) | 1.71 | 0.97 | 0.107 | 0.101 |
| Gauss l_1 | Pantheon+ | 0.828 (0.797–0.856) | sin cota | sin cota | +1.33 (0.316) | — | — | 0.079 | — |
| Gauss l_1 | DES-SN5YR | sin cota (0.989–0.997) | sin cota | sin cota | +5.35 (0.43) | — | — | — | — |
| Gauss l_1 | Dovekie | 0.889 (0.866–0.914) | sin cota | sin cota | +2.22 (0.368) | 1.07 | — | 0.085 | — |
| Gauss l_1 | Union3 | sin cota | sin cota | sin cota | +3.55 (0.501) | — | — | — | — |

Cocientes entre familias en ℓ = 1 (Matérn/OU y gaussiano/OU; OU en θ = 1):

| muestra | Matérn/OU 5 % | Matérn/OU −1.92 | gauss/OU 5 % | gauss/OU −1.92 |
|---|---|---|---|---|
| Pantheon+ | 2.33 | 2.66 | 4.91 | — |
| DES-SN5YR | 2.22 | 2.40 | — | — |
| Dovekie | 2.32 | 2.34 | 4.69 | — |
| Union3 | 2.35 | 2.14 | — | — |

Dispersión de Δρ(0.1) en la cota entre las celdas con cota de cada muestra:

| muestra | n (5 %) | rango 5 % | máx/mín 5 % | n (−1.92) | rango −1.92 | máx/mín −1.92 |
|---|---|---|---|---|---|---|
| Pantheon+ | 8 | 0.062–0.079 | 1.28 | 7 | 0.105–0.131 | 1.24 |
| DES-SN5YR | 7 | 0.093–0.149 | 1.60 | 7 | 0.100–0.138 | 1.38 |
| Dovekie | 8 | 0.064–0.108 | 1.70 | 7 | 0.089–0.126 | 1.42 |
| Union3 | 7 | 0.100–0.150 | 1.50 | 7 | 0.101–0.191 | 1.88 |

## Lo que se ve

- **Las cotas del 5 % se aflojan con DES-SN5YR (×1.4–1.9) y Union3 (×1.6–2.0) y apenas con Dovekie (×1.0–1.4); las lineales −1.92 casi no
  se mueven con DES-SN5YR (×0.94–1.06) ni con Dovekie (×0.85–1.00) y se aflojan ×1.0–1.5 con Union3.** Las dos reglas responden a cosas
  distintas: la del 5 % compara χ²_min con el χ²_obs de ΛCDM, y en DES-SN5YR y Union3 ΛCDM ajusta peor en relación a los caminos (su
  Ω_m de SN solas, 0.352 y 0.356, está más lejos del de BAO + CMB que el de Pantheon+, 0.332, o Dovekie, 0.330), así que más caminos
  igualan a ΛCDM; la lineal mide la caída desde el máximo, y el máximo sube con la tensión (ln L máximo +4.6 a +5.4 con
  DES-SN5YR, +3.3 a +3.6 con Union3, +1.7 a +2.2 con Dovekie, +0.6 a +1.5 con Pantheon+) sin que la anchura cambie mucho. Por la
  regla fijada, la subida del máximo se lee como tensión entre datos absorbida por σ_X > 0, no como detección.
- **El gaussiano en ℓ = 1 pierde la cota del 5 % con DES-SN5YR y con Union3** (P = 0.061 y 0.076 en σ_X = 1; el 98 % y el 100 % de los
  remuestreos sin cota), la conserva con Dovekie (0.889) y ya no tenía lineal con ninguna. Por la regla fijada no se salva cambiando
  de estadístico: la cota del gaussiano a ℓ = 1 no es estable frente a la muestra.
- Los cocientes entre familias en ℓ = 1 son estables: Matérn/OU 2.2–2.35 (5 %) y 2.1–2.7 (−1.92); gaussiano/OU 4.7–4.9 donde existe.
- Δρ(0.1) en la cota −1.92: 0.105–0.131 (Pantheon+, ×1.24), 0.100–0.138 (DES-SN5YR, ×1.38), 0.089–0.126 (Dovekie, ×1.42), 0.101–0.191
  (Union3, ×1.88). Con el 5 %: ×1.28, ×1.60, ×1.70, ×1.50.

## Predicción, punto por punto

1. **No cumplida.** DES-SN5YR: 5 % ×1.43–1.93 (4 de 7 celdas en 1.1–1.6), −1.92 ×0.94–1.06 (0 de 7). Union3: 5 % ×1.62–1.95
   (0 de 7), −1.92 ×0.97–1.46 (6 de 7). Las del 5 % se aflojan más de lo previsto y las lineales menos (DES) o justo (Union3).
2. **Cumplida con el 5 %, no con la lineal.** Con el 5 % Dovekie queda entre Pantheon+ y DES-SN5YR en las 7 celdas con cota en las tres
   (y en gauss ℓ = 1, donde DES no tiene cota, queda por encima de Pantheon+). Con −1.92 queda por debajo de las dos en las 7 celdas
   (×0.85–1.00 de Pantheon+).
3. **No cumplida por poco.** Con DES-SN5YR el máximo de ln L está en σ_X > 0 en las 8 celdas (σ_X = 0.05–0.43), con Δln L +4.59 a +5.35:
   por encima del +4 previsto.
4. **Cumplida donde se puede medir.** Matérn/OU en ℓ = 1 cambia −5 % a +1 % (5 %) y −10 % a −19 % (−1.92) entre muestras; gaussiano/OU
   −4 % (5 %, solo Dovekie) y no medible con DES-SN5YR ni Union3 (sin cota del gaussiano) ni con −1.92 (sin cota en ninguna).
5. **No cumplida (secundaria).** La dispersión de Δρ(0.1) en la cota −1.92 entre familias es ×1.24 en Pantheon+ pero ×1.38, ×1.42 y ×1.88
   en DES-SN5YR, Dovekie y Union3; los valores quedan en 0.089–0.191 (en 0.08–0.18 salvo Union3 OU θ = 10, 0.191). La coincidencia
   vista en FAMILIAS no se mantiene al nivel de ×1.3 en muestras no usadas.
6. **Cumplida.** Las 12 celdas OU con θ ≥ 1 (3 θ × 4 muestras) conservan cota del 5 % y −1.92.

## Respuesta

Las cotas de OU y Matérn sobreviven al cambio de muestra: la escala de ℓ y los cocientes entre familias se mantienen, y lo que cambia
(×1.0–2.0 según regla y muestra) sigue el grado de tensión de cada muestra con BAO + CMB, no la familia. La cota del gaussiano a
ℓ = 1 no sobrevive: existe con Pantheon+ y Dovekie y no con DES-SN5YR ni Union3. La coincidencia de Δρ(0.1) se degrada de ×1.24 a
×1.4–1.9: no es una cifra universal al nivel de la predicción.

## Cierre (sin nuevas corridas de cotas)

### b. Por qué se separan las dos reglas: predicción escrita antes de calcular

"El aflojamiento del 5 % ordena las muestras igual que el exceso de χ² de ΛCDM; la lineal no."

### 1. Licencias de las muestras (comprobado el 2026-11-12; detalle en `datos/FUENTES.json` y en cada `FUENTE.json`)

| muestra | licencia encontrada | qué se hizo |
|---|---|---|
| Union3 | MIT (`LICENSE` en rubind/union3_release; mismo release en Zenodo, doi:10.5281/zenodo.15825522) | se conservan los ficheros; copia de `LICENSE` en `datos/union3/` |
| DES-SN5YR | el repositorio GitHub no tiene LICENSE; el data release en Zenodo (doi:10.5281/zenodo.12720778, DES-SN5YR-1.2.zip) es CC-BY-4.0 y pide citar DES Collaboration 2024 y Sánchez et al. 2024 | se conservan `DES-SN5YR_HD.csv` y `STAT+SYS.txt.gz`: su blob git es el mismo en las etiquetas v1.2 y 1.3 (entre ambas solo cambian README y ficheros de pipeline), así que son los del release 1.2 (no se ha abierto el zip de 1.5 GB: se asume que es la etiqueta v1.2). El `README.md` de 1.3 difiere del de v1.2 y no está cubierto: sacado del repositorio |
| Dovekie | ninguna: rama main sin LICENSE, README sin condiciones de uso, sin registro Zenodo (el de Zenodo es el release de 2024) | los cuatro ficheros sacados del repositorio (`git rm`, sin reescribir la historia); los baja `descargar_datos.py` |
| Pantheon+ | ninguna: PantheonPlusSH0ES/DataRelease sin LICENSE, READMEs sin condiciones de uso, sin Zenodo | `Pantheon+SH0ES.dat` y `Pantheon+SH0ES_STAT+SYS.cov` sacados del repositorio (`git rm --cached`, copias locales conservadas, en `.gitignore`); sha256 locales idénticos a los de GitHub el 2026-11-12; los baja `descargar_datos.py` a `analisis_cmb/sn/datos/` |

`descargar_datos.py` (solo biblioteca estándar, `python3 -I`) baja lo que falte y comprueba el sha256 de los 12 ficheros listados en
`datos/FUENTES.json`, incluidos los que sí están en el repositorio; `--verify` solo comprueba. Probado: borrado y descarga del README de
Dovekie, sha256 idéntico.

### 2a. ¿Queda σ_X = 0 fuera del intervalo de la caída de 1.92? (máximo de ln L relativo a σ_X = 0, `linear.*.lnZ_max_rel` del json)

| muestra | Δln L máximo (8 celdas) | celdas con Δln L > 1.92 (σ_X = 0 fuera) |
|---|---|---|
| Pantheon+ | +0.59 … +1.50 | 0 de 8 |
| DES-SN5YR | +4.59 … +5.35 | 8 de 8 |
| Dovekie | +1.70 … +2.22 | 2 de 8 (Matérn ℓ = 1, +1.95; gaussiano ℓ = 1, +2.22) |
| Union3 | +3.30 … +3.61 | 8 de 8 |

Por la regla fijada arriba, esto se lee como tensión entre datos absorbida por σ_X > 0, no como detección. Con las mismas supernovas
recalibradas (Dovekie frente a DES-SN5YR, 1718 en común) el máximo baja de +4.6 … +5.4 a +1.7 … +2.2: la preferencia sigue a la
calibración.

### 2b. Resultado (script `muestras_sn_lcdm.py`, salida `muestras_sn_lcdm.json`; evaluaciones de χ² en los ajustes ya guardados, ninguna cota nueva)

ΛCDM plano (X = 0) con cada muestra más BAO + θ* + ω_m: 3 parámetros (h, ω_m, M), filas = 15 + n_SN. Se dan el χ²_min (coincide con el
`chi2_obs` del json), su p-valor, el reparto del χ²_min entre las 15 filas y el bloque SN, y la tensión Δχ² = χ²_min(conjunto) −
χ²_min(15 filas solas, 17.24 con 13 gdl, p = 0.19) − χ²_min(SN solas, M marginalizada), que con un parámetro compartido (Ω_m) se
compara con χ² de 1 gdl. Aflojamiento: cociente frente a Pantheon+ de la tabla final (mediana y rango sobre las celdas con cota).

| muestra | n_SN | gdl | χ²_min | χ²_min − gdl | p (χ²_min) | 15 filas / SN en el ajuste conjunto | χ²_min SN solas | Δχ² tensión | p (1 gdl) | 5 % /P+ mediana (rango) | −1.92 /P+ mediana (rango) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Pantheon+ | 1580 | 1592 | 1406.33 | −185.7 | 0.9997 | 17.34 / 1388.99 | 1386.40 | 2.68 | 0.101 | 1 | 1 |
| Dovekie | 1820 | 1832 | 1651.88 | −180.1 | 0.9989 | 17.41 / 1634.46 | 1631.42 | 3.22 | 0.073 | 1.12 (1.02–1.40) | 0.98 (0.85–1.00) |
| DES-SN5YR | 1829 | 1841 | 1666.05 | −175.0 | 0.9985 | 17.65 / 1648.40 | 1640.08 | 8.73 | 0.003 | 1.53 (1.43–1.93) | 1.04 (0.94–1.06) |
| Union3 | 22 | 34 | 45.34 | +11.3 | 0.092 | 17.32 / 28.02 | 23.96 | 4.15 | 0.042 | 1.63 (1.62–1.95) | 1.39 (0.97–1.46) |

¿Se cumple? La predicción no fijó qué es "exceso de χ²", así que se dan las dos lecturas:

- Con el χ²_min y su p-valor (lo pedido): el p-valor ordena Pantheon+ > Dovekie > DES-SN5YR > Union3 y la mediana del aflojamiento del
  5 % ordena igual (1 < 1.12 < 1.53 < 1.63); la lineal no (0.98 < 1 < 1.04 < 1.39: Dovekie queda por debajo de Pantheon+). **Cumplida**,
  con una salvedad: en las tres muestras grandes el χ²_min queda unos 3σ por debajo de sus gdl (p = 0.9985–0.9997, como es habitual con
  estas covarianzas) y el orden entre ellas lo marcan diferencias de 0.001 en p, no un exceso.
- Con la tensión Δχ² (la parte del χ²_min que añade compartir Ω_m con BAO + CMB): ordena Pantheon+ (2.7) < Dovekie (3.2) < Union3 (4.2) <
  DES-SN5YR (8.7), y el 5 % cambia DES-SN5YR y Union3 de sitio (medianas 1.53 y 1.63, rangos 1.43–1.93 y 1.62–1.95, solapados); la lineal
  tampoco sigue ese orden. **No cumplida del todo** con esta lectura: el 5 % sigue la tensión salvo entre DES-SN5YR y Union3, que no se
  distinguen.

Por qué se separan: la regla del 5 % compara el χ²_min de cada camino con el χ²_obs de ΛCDM, y cuanto mayor es la tensión más caminos con
σ_X > 0 bajan de ese χ²_obs (el máximo de ln L de 2a sube con la tensión: +0.6 … +1.5, +1.7 … +2.2, +3.3 … +3.6, +4.6 … +5.4), así que la
cota se afloja; la lineal mide la caída de 1.92 desde el máximo, es decir la anchura del perfil, que casi no cambia con la posición del
máximo (DES-SN5YR ×0.94–1.06 pese a Δχ² = 8.7).

### 2c. Dovekie no es una muestra independiente

DES-SN5YR y Dovekie comparten 1718 supernovas (por CID; `shared_sne.des5yr_vs_dovekie`), el 94 % de las 1829 y de las 1820: Dovekie es
una prueba de calibración (mismas supernovas, nueva calibración y covarianza), no una muestra independiente. Lo que cambia entre ambas
(5 % ×1.53 → ×1.12, máximo de ln L +4.6 … +5.4 → +1.7 … +2.2) mide el efecto de la calibración, no una nueva realización de los datos.

### 2d. Dos decisiones de lectura

- **Union3 en 22 tramos con z_hel = z.** El release oficial no publica supernovas sino distancias comprimidas: su README dice de
  `mu_mat_union3_cosmo=2_mu.fits` "fits file containing compressed distances for Union3+UNITY1.5; the first row is redshift, the first
  column is distance modulus, the rest of the matrix is the inverse covariance matrix" (22 nodos en z = 0.05 … 2.26). No hay tabla por
  supernova ni redshift heliocéntrico: el factor (1 + z_hel) de cada supernova ya está dentro del ajuste UNITY que produce los nodos, así
  que el modelo se evalúa en el z de cada nodo con z_hel = z (el término 5 log₁₀[(1 + z_hel)/(1 + z_HD)] no está definido para un nodo).
  Con esto se reproducen Ω_m = 0.356 y χ² = 23.96 frente a 24.0 de su Tabla 5.
- **Dovekie: se invierte el fichero de covarianza.** `DES-Dovekie-SN_Likelihood.py` (rama main, sha256 en `FUENTE.json`) dice en la línea 73
  "# Covtot_inv is the inverse of stat+sys because it makes no sense to invert covsys", reconstruye la matriz del triángulo superior
  (líneas 79–84) y en las líneas 85–86 "# Unfortunately, to make this work in cosmosis we need to invert the inverse covariance and
  return that ... sigh ..." / `C = np.linalg.inv(inv_cov)`. `load_dovekie` hace lo mismo: C = inversa de la matriz publicada (simetrizada),
  y el ajuste usa el Cholesky de C, como en las otras muestras.

### 2e. Lectura posterior (añadida el 2026-11-12; solo nota, nada recalculado)

- El p-valor del χ²_min conjunto (≈ 0.999 en Pantheon+, Dovekie y DES-SN5YR) no mide tensión: con ~1600 puntos lo domina la
  normalización de la covarianza de SN (χ²_min queda ~3σ por debajo de los gdl en las tres), no el desacuerdo con BAO + CMB. La medida
  de tensión es la columna Δχ² de tensión de la tabla de 2b (2.7, 3.2, 8.7, 4.2), y con ella la predicción de 2b se cumple a medias: el
  5 % sigue la tensión salvo que DES-SN5YR y Union3 se cruzan (Δχ² 8.7 > 4.2 pero 5 % ×1.53 < ×1.63), y la lineal no la sigue.
- Lectura de Claude, como hipótesis (formulada después de ver estos resultados; se pone a prueba en `INYECCION.md`): "la marginal sigue
  la información de los datos y la del 5 % mezcla información y tensión: DES-SN5YR tiene la misma información que Pantheon+ y más
  tensión (solo se mueve el 5 %); Union3 tiene menos información (solo se mueve todo)". Como medida de información, el error de Ω_m con
  SN solas de cada muestra, copiado de su artículo (`datos/FUENTES.json`):

| muestra | σ(Ω_m) publicado, SN solas | fuente | Δχ² tensión | 5 % /P+ | −1.92 /P+ |
|---|---|---|---|---|---|
| Pantheon+ | ± 0.018 | abstract de arXiv:2202.04077 | 2.7 | 1 | 1 |
| DES-SN5YR | ± 0.017 | abstract de arXiv:2401.02929 | 8.7 | 1.53 | 1.04 |
| Dovekie | ± 0.015 | abstract de arXiv:2511.07517 | 3.2 | 1.12 | 0.98 |
| Union3 | +0.028 / −0.026 | Tabla 5 de arXiv:2311.12098 | 4.2 | 1.63 | 1.39 |

  Con esa medida: DES-SN5YR y Dovekie tienen tanta o más información que Pantheon+ (0.017 y 0.015 frente a 0.018) y solo se mueve el 5 %,
  según su tensión; Union3 tiene menos (0.027) y se mueven las dos reglas. Es consistente con la hipótesis, pero son cuatro muestras y la
  hipótesis se escribió después de verlas.
