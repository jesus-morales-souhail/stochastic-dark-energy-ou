# Procedencia de `data/desi_dr2_local/` (comprobada el 2026-11-12; los datos no se suben al repositorio)

`data/` está en `.gitignore` (salvo este fichero). Todo lo que leen los scripts del repositorio sale de dos fuentes públicas; aquí se
da la URL del registro oficial, la versión, la fecha y el sha256 de cada fichero que un script lee de verdad. Los sha256 se calcularon
sobre los ficheros presentes en disco; el md5 del zip de Zenodo se comparó con el que publica el registro.

## 1. Fuente A: paquete suplementario de DESI DR2 Results II (Zenodo)

| campo | valor |
|---|---|
| registro | https://doi.org/10.5281/zenodo.16644577 ("DESI DR2 Results II: Supplementary data", DESI Collaboration) |
| versión / fecha de publicación | v 1.0, 2025-08-20 |
| licencia (según el registro) | CC-BY-4.0 |
| fichero | `dr2-bao-zenodo.zip`, 1 306 284 840 B |
| md5 publicado por Zenodo | `425623648b4b12515fa1943a5857fd46` |
| md5 del fichero local `data/desi_dr2_local/dr2-bao-zenodo.zip` | `425623648b4b12515fa1943a5857fd46` (idéntico) |
| sha256 del fichero local | `fde0740f705c55ff01768f1b8c5f97d102775c1e3e0b25a0afe000704076d807` |
| descargado (fecha del fichero) | 2026-07-06 |
| extraído en | `data/desi_dr2_local/dr2_data/dr2-bao-zenodo/` (figure1–6 y 12, `Readme.md`, `cosmology_chains.zip` y su extracción `cosmology_chains/`) |

Fichero de este paquete que leen los scripts (`scripts/desi_dr2_data.py::load_alpha_dv`, y su copia `analisis_cmb/everpresent/desi_dr2_data.py`;
lo usan `scripts/profile_sigma_x_desi.py`, `ou_bao_likelihood.py`, `joint_w0wa_sigma_desi.py`, `eos_efectiva.py`, `tachyonic_rank1_mle.py`,
`desi_dr2_real_bao_test.py`, `analisis_cmb/auditoria_v23/verificacion.py`, `analisis_cmb/everpresent/{ahmed_sim,profile_desi_everpresent}.py`):

| fichero (ruta bajo `dr2_data/dr2-bao-zenodo/`) | tamaño | sha256 |
|---|---|---|
| `figure6/DESI_DR2_alpha_DV_over_rs.txt` (7 α_DV isótropos, Fig. 6 del paper) | 545 B | `6f9febced44a49620e7fd42c16c40ebdbd61188e128bfe734156f9202d8a9ac1` |
| `Readme.md` del paquete | 7 406 B | `a475e6eb1e1ac81113d9b6927ac6b35a3bb9fa32e9241d6a2ecb931aea3193fe` |
| `cosmology_chains.zip` (dentro del zip principal) | 1 295 421 472 B | `b679ce2749618c6e641a74bf12bbbf56a86ac72ba2cdda7dce0b6d51578ce4c6` |

`figure5/` (multipolos) solo se consulta por `multipole_dir()` para saber si existe; ningún script lee sus ficheros.

## 2. Fuente B: vector gaussiano BAO de 13 números con su covarianza 13 × 13 (repositorio `CobayaSampler/bao_data`)

El paquete de Zenodo **no** contiene `bao_data/desi_bao_dr2/` (comprobado con `unzip -l` sobre `cosmology_chains.zip`: ninguna entrada
`bao_data`). Esa carpeta es un clon de git del repositorio de datos de Cobaya, colocado a mano dentro de la extracción:

| campo | valor |
|---|---|
| repositorio | https://github.com/CobayaSampler/bao_data (rama `master`; sin fichero LICENSE; su README cita arXiv:2503.14738 y 2503.14739 como fuente de los números DESI DR2) |
| commit del clon local | `bb0c1c9009dc76d1391300e169e8df38fd1096db` (2025-06-26) |
| último commit de los dos ficheros usados | `b7b8a36e9bccb063081f811f323cada21ab5fbdd` (2025-03-20) |
| clonado (fecha de los ficheros) | 2026-07-06 |
| ruta local | `data/desi_dr2_local/dr2_data/dr2-bao-zenodo/cosmology_chains/bao_data/` |
| comprobación 2026-11-12 | los dos ficheros bajados de `raw.githubusercontent.com/CobayaSampler/bao_data/master/desi_bao_dr2/` tienen el mismo sha256 que los locales |

Ficheros leídos (`scripts/desi_dr2_data.py::load_gaussian_bao_vector`, y la copia de everpresent; los usan `analisis_cmb/fondo/{perfil_fondo,perfil_13,fisher_diseno}.py`,
`analisis_cmb/cuesta/cuesta.py`, `analisis_cmb/oscilador/oscilador.py`, `analisis_cmb/everpresent/{bao_sim,rd_planck}.py`, y `measurement_cov(use_full_cov=True)`
para proyectar la 13 × 13 a los 7 α):

| fichero (ruta bajo `bao_data/desi_bao_dr2/`) | tamaño | sha256 |
|---|---|---|
| `desi_gaussian_bao_ALL_GCcomb_mean.txt` | 472 B | `9ac154ab583ce759c0f7eef3c978c7c70a6ead2d18774caceadf1a350a640585` |
| `desi_gaussian_bao_ALL_GCcomb_cov.txt` | 2 547 B | `252a143274c8a07c78694c119617d36594f6d7965d00319ca611c6ffb886e509` |

El mismo vector y covarianza son los de la verosimilitud `bao.desi_dr2` de Cobaya (DESI DR2 Results II, Tabla IV; compresión gaussiana oficial).

## 3. Ficheros derivados (generados por scripts del repositorio, no son datos oficiales)

| fichero | lo escribe | sha256 (en disco) | lo lee alguien |
|---|---|---|---|
| `bao_data/desi_bao_dr2/desi_cov_alpha_iso.npy` (fecha 2026-11-12) | `scripts/desi_dr2_data.py::load_alpha_dv` (proyección 7 × 7) | `fdd1d37e9a880b4431171f272142cbed1d5cc20513d5afaaeb7290338cfa05ad` | no (copia de trabajo) |
| `bao_data/desi_bao_dr2/desi_cov_alpha_iso.txt` (2026-08-10) y `extraer_cov_alpha.py` | copia local sin seguir dentro del clon | `a26ec193b0730251c2274139ff7171b5c6b6972d3974a50857340535afda8035` | no |
| `data/desi_dr2_local/desi_cov_alpha_iso.npy` (2026-07-06) | `scripts/extraer_cov_cobaya.py` | `ac8f6b4a8b7b56ba87105b64503357d631ffd9d43832845925bf4d68dad63130` | no |
| `results/desi_cov/desi_gaussian_bao_ALL_GCcomb_cov.txt` (seguido por git) | `load_gaussian_bao_vector` (copia simetrizada de la 13 × 13) | `c0965b686536c50bb39a6589b1c386d234e1d7f1a9b475595c89dc50e3b0377a` | pipelines que no dependen de `data/` |
| `results/desi_cov/desi_cov_alpha_iso_7x7.txt` (seguido por git) | `load_alpha_dv` (proyección a los 7 α) | `0efd0ff9743e6ecde31998c7f19979770b1c9a76d3aeaa8730c1caabda25bae0` | ídem |

`data/desi_dr2_local/cov/` está vacía. `archive/scripts/cross_correlation_DESI.py` (retirado) apunta a `data/LRG_clustering_*.fits`,
`data/pantheon/`, `data/desi/LRG_NGC_clustering_DR1.dat.fits`: ninguno existe en `data/`; esos ficheros DR1 viven fuera del repositorio
(`~/Proyectos/01_Fisica_y_Cosmologia/proyecto_unificacion/data/desi/`), ver `manuscript/ERRATAS_V3.md` fila 128.

## 4. Para volver a montar `data/` desde cero

1. Bajar `dr2-bao-zenodo.zip` de https://zenodo.org/records/16644577, comprobar `md5 = 425623648b4b12515fa1943a5857fd46`, extraer en
   `data/desi_dr2_local/dr2_data/` y extraer dentro `cosmology_chains.zip` (queda `dr2-bao-zenodo/cosmology_chains/`).
2. `git clone https://github.com/CobayaSampler/bao_data data/desi_dr2_local/dr2_data/dr2-bao-zenodo/cosmology_chains/bao_data` y
   comprobar los dos sha256 de la sección 2 (`sha256sum bao_data/desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_{mean,cov}.txt`).
3. Comprobar `sha256sum dr2_data/dr2-bao-zenodo/figure6/DESI_DR2_alpha_DV_over_rs.txt` contra la sección 1.
