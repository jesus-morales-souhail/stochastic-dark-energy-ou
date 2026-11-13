# Procedencia de la medida BAO-only de SDSS/eBOSS DR16 LRG (z_eff = 0.698) usada en `analisis_cmb/auditoria_v23/EBOSS_LRG.md`

Comprobada el 2026-11-13. Los dos ficheros se leen como números por `analisis_cmb/auditoria_v23/eboss_lrg_analisis.py`. En el paso previo a la
comparación (modo `potencia`, commit del paso 4 de EBOSS_LRG.md) el script lee **solo** la covarianza y los metadatos del vector de datos (z_eff y el
tipo de cada fila: columnas 0 y 2); los valores observados (columna 1) se leen solo en el modo `resultado`, después de ese commit. Esta nota se escribió
sin los valores; la sección "Valores" se completó en el paso 5.

## Ficheros (clon local de `CobayaSampler/bao_data`, el mismo del que salen los 13 números de DESI DR2; véase `data/PROCEDENCIA_DESI_DR2.md` §2)

| campo | valor |
|---|---|
| repositorio | https://github.com/CobayaSampler/bao_data (rama `master`; sin fichero LICENSE) |
| commit del clon local | `bb0c1c9009dc76d1391300e169e8df38fd1096db` (2025-06-26); último commit de estos dos ficheros `48dca631298d34151727cedcddd075a3c19eb809` (2022-02-15) |
| clonado (fecha de los ficheros) | 2026-07-06 |
| ruta local | `data/desi_dr2_local/dr2_data/dr2-bao-zenodo/cosmology_chains/bao_data/` |

| fichero | tamaño | sha256 | contenido |
|---|---|---|---|
| `sdss_DR16_LRG_BAO_DMDH.dat` | 78 B | `b3317e7590799fad71a9a707023d0743c14d87399d6bb4129965d6a5732d91be` | 2 filas: `z_eff  valor  tipo`, tipos `DM_over_rs` y `DH_over_rs`, z_eff = 0.698 |
| `sdss_DR16_LRG_BAO_DMDH_covtot.txt` | 81 B | `1a45e106f8e2bbf8742a6c3d4a9c11bdc288801fc6824e0db8cfbab4290f6160` | covarianza 2 × 2 de (D_M/r_d, D_H/r_d): diagonal 0.1076634, 0.2838176 (σ = 0.328, 0.533), covarianza −0.0583182 (coeficiente −0.333) |

Comprobaciones de la URL (2026-11-13):

- `https://raw.githubusercontent.com/CobayaSampler/bao_data/master/sdss_DR16_LRG_BAO_DMDH.dat` y `..._covtot.txt`: sha256 **idénticos** a los dos ficheros locales.
- Publicación original de SDSS: `https://svn.sdss.org/public/data/eboss/DR16cosmo/tags/v1_0_1/likelihoods/BAO-only/`. Su `README.txt` (leído el 2026-11-13) describe los ficheros
  `*_DMDH` como "Gaussian approximations to the joint likelihood of D_M(zeff)/r_d and D_H(zeff)/r_d", con el vector de datos y la covarianza en ficheros separados,
  r_d = 147.8 Mpc en la cosmología fiducial de eBOSS, y dice que los resultados LRG de DR16 sustituyen a los de DR12 en 0.5 < z < 0.75. En esa carpeta,
  `sdss_DR16_LRG_BAO_DMDH_covtot.txt` tiene sha256 **idéntico** al local; la URL de `sdss_DR16_LRG_BAO_DMDH.dat` devolvió 404 el 2026-11-13 (no se pudo comparar
  el vector de datos con el servidor de SDSS; se usa la copia de Cobaya, idéntica a la de GitHub).

## Qué producto es

- **BAO-only** (carpeta `BAO-only` de la publicación de SDSS; nombre `*_BAO_DMDH`), no BAO+RSD (`*_FSBAO_DMDHfs8`, que también está en el clon y no se usa) ni el
  consenso final. Valores y covarianza del mismo producto (los dos ficheros van juntos en la publicación).
- Según Bautista et al. 2021 (arXiv:2007.08993, resumen leído el 2026-11-13): la muestra combina 174 816 LRG de eBOSS con 202 642 galaxias CMASS de BOSS, en
  0.6 < z < 1, con z_eff = 0.698; los resultados finales combinan el análisis en espacio de configuración (Bautista et al.) con el de espacio de Fourier
  (Gil-Marín et al. 2020, arXiv:2007.08994). El fichero BAO-only de la publicación es el **consenso BAO-only de los dos espacios** (combinación de configuración y
  Fourier; Alam et al. 2021, arXiv:2007.08991, Tabla 3), no el de un solo análisis. Ni el `README.txt` de SDSS ni el clon de Cobaya lo dicen fichero a fichero:
  esta identificación viene de la literatura y del nombre y carpeta del fichero.
- Sin covarianza cruzada publicada entre DESI DR2 LRG2 y eBOSS DR16 LRG (no hay simulaciones comunes publicadas con esa matriz); EBOSS_LRG.md usa una
  familia hipotética ρ_x.

## Valores (completado en el paso 5, después del commit del paso 4: `b8d48a5`)

Leídos como números por `eboss_lrg_analisis.py resultado` el 2026-11-13 (sha256 del fichero sin cambio): z_eff = 0.698, D_M/r_d = 17.85824 ± 0.3281,
D_H/r_d = 19.32575 ± 0.5327, coeficiente de correlación -0.3336. Coinciden con el consenso BAO-only de DR16 LRG citado en la literatura
(D_M/r_d = 17.86 ± 0.33, D_H/r_d = 19.33 ± 0.53), distinto del consenso BAO+RSD del resumen de Bautista et al. (17.65 ± 0.30, 19.77 ± 0.47), que no se usa.
