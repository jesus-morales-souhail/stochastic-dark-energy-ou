# REPRODUCCION: preparación de la reproducción a ciegas (2026-11-12)

Qué se preparó: una carpeta fuera de este repositorio, `~/Proyectos/01_Fisica_y_Cosmologia/reproduccion_ciega/`, para que una
persona sin acceso a nuestro código calcule la cota de la marginal lineal −1.92 sobre σ_X del OU con θ = 1 con BAO + θ* + ω_m +
Pantheon+. La carpeta no está en git ni se añadirá; este fichero solo deja constancia de qué contiene y del sha256 de la
especificación, sin el resultado esperado ni ninguna cifra nuestra.

Contenido de la carpeta:

- `ESPECIFICACION.md`: especificación autosuficiente (datos y columnas, covarianzas, modelo exacto en convenio B, distancias, θ* y
  su calibración, selección de supernovas y marginalización de M, priors, punto fiducial, malla en u, núcleo OU y cómo entra X(0),
  respuesta lineal con fórmulas y comprobación por diferencias finitas, verosimilitud marginal completa, malla de σ_X, máximo global,
  regla de la caída de 1.92, entregables y comprobaciones). Sin cifras de resultados.
- `datos/pantheon_plus/`: copias de `Pantheon+SH0ES.dat` y `Pantheon+SH0ES_STAT+SYS.cov` (release público de Pantheon+SH0ES).
- `datos/desi_dr2/`: copias de `desi_gaussian_bao_ALL_GCcomb_mean.txt` y `_cov.txt` (compresión gaussiana oficial de DESI DR2).
- `tablas/rd_rs_zstar_vs_omega_m.csv`: r_d(ω_m), r_s(z*; ω_m) y z*(ω_m) en 361 nodos, exportados de nuestras tablas; y
  `tablas/constantes.json` (c, Ω_r h², z* de referencia, los dos números de Planck con sus errores, la calibración de θ*, la malla).
- `MANIFIESTO_SHA256.txt`: sha256 de todos los ficheros anteriores, comprobado contra los ficheros presentes al crearlo.

No contiene ningún script de este proyecto, ninguna ruta a `analisis_cmb/` ni a `auditoria_v23/`, ni ninguna cota, ln Z o máximo
nuestros (comprobado con grep en los ficheros de texto; en los de datos no se busca).

sha256 de `ESPECIFICACION.md` (versión 1): `960e97881478cd79f47167835f4b76a507d1e80f33bd925ca9ac60ebb2b82596`

Regla fijada: la especificación se escribió antes de que el segundo investigador haga el cálculo y sin darle acceso a nuestros
scripts; cualquier duda sobre el método se responde por escrito y sin cifras. La reproducción no ha empezado.
